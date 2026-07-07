# Tema 1: Sestavljanje genoma *de novo* in napovedovanje genov (Genomika)

## Cilji vaje
- Razumeti principe sestavljanja kratkih odčitkov sekvenciranja v daljše fragmente (*soseske* ali *contigs*).
- Zagnati in parametrizirati program SPAdes za sestavljanje genoma *Saccharomyces cerevisiae*.
- Uporabiti program Augustus za napovedovanje genov *ab initio* za proteinske kodirajoče regije.
- Ovrednotiti kakovost sestave genoma in napovedanih genov na podlagi bioloških podatkov.

## Uvod
Sekvenciranje naslednje generacije (NGS) generira milijone kratkih nukleotidnih zaporedij (odčitkov), ki sami po sebi ne dajejo celostne informacije o zgradbi genoma. Da bi rekonstruirali celoten genom organizma, moramo te kratke odčitke združiti. V primerih, ko nimamo na voljo referenčnega genoma, uporabimo pristop **sestavljanja genoma *de novo* (de novo assembly)**. 

Algoritmi za sestavljanje genoma *de novo* običajno temeljijo na De Bruijnovih grafih, kjer odčitke razdelimo na krajše prekrivajoče se podsekvence dolžine *k* (t.i. *k-mer-e*). S povezovanjem teh vozlišč v grafu program rekonstruira daljša zaporedja, imenovana **soseske (contigs)**, in jih nato poveže v **ogrodja (scaffolds)**. V tej vaji bomo uporabili program **SPAdes**, ki velja za zlati standard pri sestavljanju genoma *de novo* bakterijskih in manjših evkariontskih genomov.

Ko pridobimo sestavljen genom, je naslednji korak **anotacija genov**. Ker evkariontski geni vsebujejo introne in eksone, je iskanje odprtih bralnih okvirjev (ORF) zapleteno. Orodje **Augustus** uporablja skrite modele Markova (HMM) za statistično napovedovanje začetnih kodonov, intronskih izrezovalnih mest in stop kodonov neposredno iz nukleotidnega zaporedja na podlagi predhodno naučenih modelov za specifične organizme (v našem primeru za kvasovko *S. cerevisiae*).

## Bioinformatska orodja
| Orodje | Namen | Ključna lastnost |
|:---|:---|:---|
| **fastp** | Kontrola kakovosti odčitkov | Hiter in celovit pregled ter filtriranje NGS odčitkov |
| **SPAdes** | Sestavljanje genoma *de novo* | Uporablja več različnih velikosti k-merov za optimalno sestavljanje |
| **QUAST** | Ocena sestave genoma | Izračuna metrike kakovosti (npr. N50) za sestavljena ogrodja |
| **Augustus** | Napovedovanje genov *ab initio* | Statistično napovedovanje ekson-intronske strukture z uporabo HMM |

## Potek vaje

### Korak 0: Prenos gradiva in inicializacija okolja
**Opis:** Če tega še niste storili, morate najprej prenesti gradivo na svoj računalnik in inicializirati bioinformatsko okolje.

Za prenos gradiva iz repozitorija in zagon samodejne namestitve zaženite v terminalu naslednje ukaze:
```bash
curl -fsSL https://pixi.sh/install.sh | bash
source ~/.bashrc
git clone https://github.com/ceneg/bioinformatics
cd bioinformatics
pixi run test
```
*Opomba: Ukaz `pixi run test` bo samodejno namestil vsa potrebna bioinformatska orodja in prenesel začetne podatke za vse vaje.*

> **Kaj storiti v primeru izgube podatkov?**
> Če med reševanjem vaje pomotoma izbrišete vhodne datoteke ali pa se vaše okolje znajde v nepredvidenem stanju, lahko kadarkoli ponovno zaženete ukaz `pixi run test` v korenski mapi projekta (`bioinformatics`). Ta ukaz ne bo izbrisal vaših lastnih rezultatov in analiz, temveč bo zgolj ponovno prenesel manjkajoče začetne datoteke in popravil morebitne napake v namestitvi okolja.

---

### Korak 1: Priprava okolja in pregled vhodnih podatkov
**Opis:** Pred začetkom dela aktivirajte ustrezno okolje Conda/Pixi, ki vsebuje zahtevana orodja, ter preverite kakovost in prisotnost vhodnih datotek v mapi `data/session1/`.

Zaženite naslednji ukaz za aktivacijo okolja:
```bash
pixi shell -e genomika
```

⚠️ Preverite število odčitkov v vhodnih datotekah FASTQ s pomočjo ukazne vrstice:
```bash
zgrep -c "^@" data/session1/yeast_chrom1_R1.fastq.gz
```

> **Pričakovani rezultati:**
> Izpisati se mora število vrstic, ki se začnejo z znakom `@`, kar predstavlja skupno število odčitkov (približno 23.000).

---

### Korak 2: Kontrola kakovosti odčitkov s programom fastp
**Opis:** Pred sestavljanjem genoma je priporočljivo preveriti kakovost vhodnih odčitkov. Uporabili bomo orodje `fastp`.

> ⚠️ **Naloga:** Zadenite ukaz za hitro analizo odčitkov.
>
> ```bash
> fastp -i data/session1/yeast_chrom1_R1.fastq.gz -I data/session1/yeast_chrom1_R2.fastq.gz -h fastp_report.html
> ```

> **Pričakovani rezultati:**
> V trenutni mapi se bo ustvarila datoteka `fastp_report.html`. Odprite jo v spletnem brskalniku in preglejte osnovno statistiko vaših odčitkov (npr. dolžino, kakovost baz in vsebnost GC).

---

### Korak 3: Sestavljanje genoma *de novo* s programom SPAdes
**Opis:** Uporabili bomo program `spades.py` za sestavljanje naših sparjenih odčitkov v sklope. Ker gre za izobraževalni primer, bomo z parametrom `--only-assembler` preskočili fazo popravljanja napak v odčitkih, kar bo bistveno pospešilo izvedbo.

> ⚠️ **Naloga:** Napišite ukaz za zagon SPAdes sestavljanja s sparjenimi odčitki (R1 in R2) ter izhodnim direktorijem `spades_out`.
>
> ```bash
> # TODO: Napišite ukaz za zagon spades.py z -1, -2 in -o
> # Namig: Poglejte dokumentacijo orodja z `spades.py --help`
> ```

> **Pričakovani rezultati:**
> V mapi `spades_out/` bo nastalo več datotek, ključna za nas pa je `scaffolds.fasta`, ki vsebuje sestavljena ogrodja genoma.

---

### Korak 4: Ocena kakovosti sestave genoma z orodjem QUAST
**Opis:** Da bi ugotovili, kako dobro je orodje SPAdes sestavilo genom, bomo analizirali izhodno datoteko z orodjem QUAST.

> ⚠️ **Naloga:** Zaženite QUAST na sestavljenih ogrodjih.
>
> ```bash
> quast spades_out/scaffolds.fasta -o quast_out
> ```

> **Pričakovani rezultati:**
> V mapi `quast_out/` boste našli datoteko `report.html` (in `report.txt`). Preglejte metrike, kot so N50, število sosesk (contigs) in skupna dolžina sestave.

---

### Korak 5: Napovedovanje genov *ab initio* s programom Augustus
**Opis:** Iz sestavljenega genoma bomo sedaj napovedali mesta genov s programom Augustus. Parameter `--species` nastavite na `saccharomyces`.

> ⚠️ **Naloga:** Napišite ukaz za zagon programa Augustus na datoteki `spades_out/scaffolds.fasta` z izvozom v formatu GFF3 v datoteko `genes.gff3`.
>
> ```bash
> # TODO: Napišite ukaz za napovedovanje genov z augustus
> # Namig: Uporabite možnost --gff3=on in preusmeritev izhoda > genes.gff3
> ```

> **Pričakovani rezultati:**
> Datoteka `genes.gff3` bo vsebovala koordinate napovedanih genov, eksonov in intronov ter aminokislinska zaporedja napovedanih proteinov.

---

## Vprašanja za razmislek
1. Koliko ogrodij (scaffolds) je program SPAdes uspel sestaviti iz naših podatkov? (Namig: Preštejte znak `>` v datoteki `scaffolds.fasta`).
2. Zakaj se dolžine sestavljenih odrov razlikujejo od dolžine celotnega kromosoma kvasovke? Kateri deli kromosoma so najtežji za sestavljanje?
3. Kako prisotnost intronov v evkariontskih genomih otežuje napovedovanje genov v primerjavi s prokarionti?
4. Kakšna je razlika med napovedovanjem genov *ab initio* (program Augustus) in anotacijo na podlagi poravnave znanih zaporedij proteinov/mRNA?

## Slovar pojmov
| Pojem | Razlaga |
|:---|:---|
| **Soseska (Contig)** | Neprekinjeno sestavljeno zaporedje DNK, dobljeno z ujemanjem prekrivanj med odčitki |
| **Ogrodje (Scaffold)** | Zaporedje sosesk z znano orientacijo in vrstnim redom, ločenih s prazninami (N-ji) |
| **K-mer** | Podsekvenca dolžine k, pridobljena iz odčitkov za gradnjo De Bruijnovega grafa |
| **Anotacija** | Proces označevanja lokacij genov in določanja njihove biološke funkcije v genomu |
