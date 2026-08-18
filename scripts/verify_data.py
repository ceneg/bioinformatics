#!/usr/bin/env python3
"""Verify — and, where possible, repair — the course input datasets.

Every file under `data/session*/` is tracked in git, so a student who deletes or
corrupts an input can get it back byte-exact without any download. This script
detects missing or truncated inputs and restores them from the git index.

Exit codes: 0 = all data present and valid, 1 = validation failed.
"""
import gzip
import pathlib
import subprocess
import sys

REPO_ROOT = pathlib.Path(__file__).resolve().parent.parent
DATA_DIR = REPO_ROOT / "data"

# Session 5 samples 9 and 10 are deliberately shallow (~500-800 reads); every
# other FASTQ should carry at least 1000 records.
LOW_DEPTH_SAMPLES = ("Sample9", "Sample10")
MIN_RECORDS_DEFAULT = 1000
MIN_RECORDS_LOW_DEPTH = 500

EXPECTED_VCF_VARIANTS = 500
EXPECTED_VCF_SAMPLES = 20


def is_git_worktree() -> bool:
    """True if the repository is a usable git working tree."""
    try:
        result = subprocess.run(
            ["git", "rev-parse", "--is-inside-work-tree"],
            cwd=REPO_ROOT, capture_output=True, text=True, timeout=15,
        )
        return result.returncode == 0 and result.stdout.strip() == "true"
    except (OSError, subprocess.SubprocessError):
        return False


def tracked_data_files() -> list[pathlib.Path]:
    """All input files under data/ that git knows about."""
    try:
        result = subprocess.run(
            ["git", "ls-files", "--", "data"],
            cwd=REPO_ROOT, capture_output=True, text=True, timeout=30,
        )
        if result.returncode != 0:
            return []
        return [REPO_ROOT / line for line in result.stdout.splitlines() if line.strip()]
    except (OSError, subprocess.SubprocessError):
        return []


def restore_from_git(paths: list[pathlib.Path]) -> list[pathlib.Path]:
    """Restore the given paths from HEAD. Returns the ones that came back."""
    rel = [str(p.relative_to(REPO_ROOT)) for p in paths]
    # `git restore` exists from git 2.23; fall back to `git checkout` on older versions.
    for cmd in (["git", "restore", "--source=HEAD", "--"], ["git", "checkout", "HEAD", "--"]):
        try:
            result = subprocess.run(
                cmd + rel, cwd=REPO_ROOT, capture_output=True, text=True, timeout=300
            )
        except (OSError, subprocess.SubprocessError):
            continue
        if result.returncode == 0:
            break
    return [p for p in paths if p.exists() and p.stat().st_size > 0]


def modified_data_files() -> list[pathlib.Path]:
    """Tracked input files whose content differs from HEAD.

    Catches corruption and accidental edits, which a missing/empty check alone
    would miss — a file overwritten with junk is still non-empty.
    """
    try:
        result = subprocess.run(
            ["git", "diff", "--name-only", "HEAD", "--", "data"],
            cwd=REPO_ROOT, capture_output=True, text=True, timeout=60,
        )
        if result.returncode != 0:
            return []
        return [REPO_ROOT / line for line in result.stdout.splitlines() if line.strip()]
    except (OSError, subprocess.SubprocessError):
        return []


def repair_missing_inputs() -> bool:
    """Restore any tracked input that is missing, empty or modified.

    Returns False only if something could not be restored.
    """
    tracked = tracked_data_files()
    if not tracked:
        return True  # nothing to check against (not a git checkout)

    absent = [p for p in tracked if not p.exists() or p.stat().st_size == 0]
    changed = [p for p in modified_data_files() if p not in absent]
    broken = absent + changed
    if not broken:
        return True

    print(f"\n=== Repairing {len(broken)} missing or damaged input file(s) ===")
    for p in absent:
        print(f"[MISSING] {p.relative_to(REPO_ROOT)}")
    for p in changed:
        print(f"[MODIFIED] {p.relative_to(REPO_ROOT)} (content differs from the original)")

    # Restoring over a write-protected file works because the directory is
    # writable, but clear the flag first so git never has to fight the mode.
    for p in broken:
        if p.exists():
            try:
                p.chmod(p.stat().st_mode | 0o200)
            except OSError:
                pass

    restored = restore_from_git(broken)
    for p in restored:
        print(f"[RESTORED] {p.relative_to(REPO_ROOT)}")

    still_broken = [p for p in broken if p not in restored]
    for p in still_broken:
        print(f"[FAIL] Could not restore {p.relative_to(REPO_ROOT)}")
    print()
    return not still_broken


def count_fastq_records(path: pathlib.Path) -> int:
    opener = gzip.open if path.name.endswith(".gz") else open
    with opener(path, "rt") as handle:
        return sum(1 for _ in handle) // 4


def verify_fastq_files() -> bool:
    patterns = ("session*/**/*.fastq", "session*/**/*.fastq.gz",
                "session*/*.fastq", "session*/*.fastq.gz")
    paths = sorted({p for pattern in patterns for p in DATA_DIR.glob(pattern)})

    ok = True
    for path in paths:
        record_count = count_fastq_records(path)
        size_mb = path.stat().st_size / (1024 * 1024)
        expected_min = (MIN_RECORDS_LOW_DEPTH
                        if any(s in path.name for s in LOW_DEPTH_SAMPLES)
                        else MIN_RECORDS_DEFAULT)

        if record_count < expected_min:
            print(f"[FAIL] {path.relative_to(REPO_ROOT)}: only {record_count} records "
                  f"(Expected >= {expected_min})")
            ok = False
        else:
            print(f"[OK] {path.relative_to(REPO_ROOT)}: {record_count} records ({size_mb:.2f} MB)")
    return ok


def verify_vcf() -> bool:
    vcf_path = DATA_DIR / "session4" / "population_structure.vcf"
    if not vcf_path.exists():
        print("[FAIL] VCF file does not exist and could not be restored.")
        return False

    lines = vcf_path.read_text().splitlines()
    header_lines = [line for line in lines if line.startswith("#")]
    variant_lines = [line for line in lines if not line.startswith("#")]
    samples = header_lines[-1].split("\t")[9:] if header_lines else []

    print("\n=== VCF Validation ===")
    print(f"Header lines count: {len(header_lines)}")
    print(f"Variant records count: {len(variant_lines)}")
    print(f"Samples found in VCF ({len(samples)}): {', '.join(samples)}")

    ok = True
    if len(variant_lines) != EXPECTED_VCF_VARIANTS:
        print(f"[FAIL] VCF variant count is {len(variant_lines)}, "
              f"expected {EXPECTED_VCF_VARIANTS}.")
        ok = False
    if len(samples) != EXPECTED_VCF_SAMPLES:
        print(f"[FAIL] VCF sample count is {len(samples)}, expected {EXPECTED_VCF_SAMPLES}.")
        ok = False
    return ok


def main() -> None:
    print("=== Verification of Simulated Data ===")

    if is_git_worktree():
        repairable = repair_missing_inputs()
    else:
        repairable = True
        print("[INFO] Not a git working tree — automatic repair is unavailable.")
        print("       If input files are missing, re-download the repository from")
        print("       https://github.com/ceneg/bioinformatics\n")

    fastq_ok = verify_fastq_files()
    vcf_ok = verify_vcf()

    if repairable and fastq_ok and vcf_ok:
        print("\n[SUCCESS] All datasets successfully validated.")
        sys.exit(0)

    print("\n[ERROR] Validation failed.")
    if not is_git_worktree():
        print("Re-clone the repository to restore the input data:")
        print("  git clone https://github.com/ceneg/bioinformatics")
    else:
        print("Try restoring the input data manually with:")
        print("  git restore --source=HEAD -- data/")
    sys.exit(1)


if __name__ == "__main__":
    main()
