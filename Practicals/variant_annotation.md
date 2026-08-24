
# Clinical Genomics Practical: Part 1 Variant Annotation



#### By Evelyn Collen
(Prac repurposed from the excellent Julien Soubrier)

## **1. Introduction**

### Quick refresher on clinical genomics and finding the diagnostic variant

As we touched on in the recordings, clinical genomics can be separated into two groups based on their ultimate goals: diagnostics and research.

Many labs in major hospitals, pathology centres, independent research institutes and universities do a mixture of both work.

Research has the goal of adding to human knowledge and also informing clinical practice (_i.e.,_ can we further inform a patient's genetic diagnosis? Can we help provide diagnoses at the cutting edge of research, where a standard clinical workflow couldn't?). Research asks: what unknowns are discoverable, what hypotheses can be made, what workflows and treatments can be applied and tested? 

The diagnostic mindset takes an inverse approach to these research questions. The goal here is to identify _already known_ patterns in our genomics data, with an emphasis of using well established and highly tested and validated "best practices". Diagnostic tests are also about identifying known issues in the quickest and most reliable way possible, because downstream clinical decisions are reliant on the information.

A critical thing to consider about high-throughput sequencing is the enormous amount of information that is obtained. Using the earring analogy from our lecture video, finding an earring (~1mm or basepair) on your roadtrip from Melbourne to Perth (3 billiom mm, or basepairs) is an enormous task!

Additionally, a lot of the information that we obtain from sequencing data is noisy and difficult to interpret. This is the reason why, in a lot of genetic pathology testing, we have primarily focussed on the "protein-coding" regions (~1% of the total sequence), where we know there is a good chance that a genetic variation may bring about a phenotypic change.
Non-coding and intronic regions of the genome, or the other 99%, are much more difficult to interpret, although many projects such as the [Epigenomics Roadmap](http://www.roadmapepigenomics.org/) or [Encyclopedia of DNA elements (ENCODE)](https://www.encodeproject.org/) are trying to change that. Slowly but surely, clinical is chasing after research,  by identifying functional regions of the non-coding genome that impact gene expression.

In this tutorial, we are firstly going to look at ways in which we can give context and functionality to variants, by doing a process called "variant annotation".

### 1.1 Reminder about virtual Machines

As usual we will be connecting the virtual machines: 

**Please [go here](../../Course_materials/vm_login_instructions.md) for instructions on connecting to your VM.**

### 1.2 Learning Outcomes


1. Get comfy manipulating annotation info in a vcf file
2. Understand the monumental task that is trying to find diagnostic variants in the human genome
3. Refresher on vcf file specifications
4. Learn about common variant annotations and annotation databases
5. Get an idea of common annotations and the logic behind why we include them



## **2. Adding sequential annotations**

### 2.1 This week's data

This week's tutorial uses three samples from a [gemini database tutorial](https://s3.amazonaws.com/gemini-tutorials/Gemini-Recessive-Tutorial.pdf) written by Aaron Quinlan (University of Utah).
Aaron's group has written some really helpful pieces of software including `bedtools`, `giggle` and `vcfanno`, which are all becoming standard tools in the toolkit of clinical researchers.

We'll be using variants that have already been called using the GenomeAnalysisToolkit (GATK), that was developed by the Broad Institute (Cambridge, USA).
The data has been pre-generated for you:
- 1 gzipped compressed Variant Call Format (VCF) = _trio.trim.vep.vcf.gz_
- 2 pedigree files that contain two separate examples looking at dominant and recessive disorders = _recessive.ped / dominant.ped_
- 1 SNP ID annotation file from dbSNP = _hg19.dbSNP.vcf.gz_

We will mostly use **BCFtools** for this tutorial:
https://samtools.github.io/bcftools/bcftools.html

As our VCF file already has previous annotations attached. Let's start by stripping off that information so we can start the process at the start.

Think of annotations as a map of landmarks we can use to find that pesky earring on our roadtrip!

First let's make a working directory, copy the data over, and activate a standalone conda environment that contains BCFtools:

```bash
mkdir -p ~/clinical_genomics && cd $_
cp /shared/data/clinical_genomics/*  ~/clinical_genomics/
source activate bioinf
```
It's always a good idea to stare at your files before and after you change them! Let's have a look into our pre-annotated vcf file. If you scroll down a little, you should see a whole of information in the FILTER and INFO columns:

```bash
zless trio.trim.vep.vcf.gz
```

Now let's use BCFtools to remove fields but keep GT (genotypes), and then index the output:

```bash
bcftools annotate -x FILTER,INFO,^FORMAT/GT trio.trim.vep.vcf.gz -Oz -o trio.trim.vcf.gz
bcftools index -t trio.trim.vcf.gz
zless trio.trim.vcf.gz
```

_NOTE ON INDEXING: In order to subset or retrieve data from a tab-delimited file (or any other delimited file for that matter), it is helpful to use an index. File indexes are a bit like the contact list in your phone, sorted alphabetically, so you can find your friends' phone numbers_.
_A file index, commonly created by the program [`tabix`](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3042176/), can be created on most standard bioinformatics files (BED, VCF etc)._
_It is good practice to create an index every time you make a new VCF file. A number of variant toolkit's (`gatk`, `picard`, `sambamb`, etc...) will often create an index automatically for you. You can either use the `tabix` program using the VCF prefix (`tabix -p vcf`) or use the `bcftools` sub-command `bcftools index`_

**NOTE:** If you get an _[E::hts_idx_push] Unsorted positions on sequence #1_ error on the index command, you can quickly sort the file again using the `bcftools sort` command and reindexing.

List the files in the directory and see what is produced.


### 2.2 Reminder on VCF files and genotypes

Firstly, lets review the Variant Call Format (VCF) file. 

This is the standard file for listing variants that are 'called' via a variant calling tools such as `BCFtools`, `GATK`, `freebayes` or `DeepVariant`.
The basic premise for calling variants is to identify both alleles (in a diploid genome).
To do this, sequenced DNA fragments ("reads") are aligned to the reference sequence and the base at each position is determined, counted, and filtered to determine whether the site is different from the reference sequence.
If the site is polymorphic, we count the bases at the position and determine whether it is a homozygous or a heterozygous variant.
In this file, a "0" denotes the reference base and "1" as the alternate base at that position.
So a heterozygous variant is "0/1" and a homozygous alternate variant is "1/1".
Additionally, it is possible to have multi-allelic sites, so additional alternate alleles are coded as greater than 1 (e.g. 2/2 or 0/2).

### 2.3 VCF headers refresher

Let's zoom into our vcf - like all vcfs, there is the header, with the full VCF information.
The header is denoted by lines that start with two # (i.e. ^##).
The name of the fields for the rest of the file (that contain the actual results) is denoted by lines that start with only one # (i.e. ^#)

```bash
zcat trio.trim.vcf.gz | head
```

Headers have heaps of metadata regarding the aligned reference genome, and gives you all the steps that were run to make the file, as well as the definitions of the specific fields and tags within the file.

```bash
# View ONLY the header using BCFtools
bcftools view -h trio.trim.vcf.gz
```
No need to spend a lot of time on this - but if you need, you can refer to the formal vcf specs:
[VCF file format spec](https://samtools.github.io/hts-specs/VCFv4.2.pdf)

---
**Questions:**
---

1. What was the name of the reference genome file was used to make the VCF file?
2. What program was used to call variants?
3. What are the names of the 3 individuals that are sampled in the VCF file?
4. Tricky question: the VCF file is subsetted to only include variants from specific chromosomes. Can you find out what are the names of those chromosomes, and how many individual variants are found in each, using some commandline tricks?

---
<details>
<summary>Answers</summary>
1. Homo_sapiens_assembly19.fasta 
2. UnifiedGenotyper
3. 1805, 1847, 4805 
4. You could run:

```bash
zgrep -v '^#' trio.trim.vcf.gz | cut -f1 | sort -V | uniq -c
```
to get this info, or if you're a bcftools whiz: 

```bash
bcftools query -f '%CHROM\n' trio.trim.vcf.gz | sort -V | uniq -c | awk '{print "chr"$2, $1}'
```

That should give you the following variants per chromosome:

chr2 7081
chr15 2922
chr17 4789
chr22 2445

</details>

### 2.4 Adding a FILTER tag

Not all variants are created equal.
Similar to the genotype and alignment quality metrics (base and mapping quality), the VCF file contains a QUALITY field that is also phred scaled.

QUAL is the phred-scaled quality score for the assertion made in ALT. _i.e.,_ -10log_10 prob(call in ALT is _wrong_). If ALT is ”.” (which means no variant), then this is -10log_10 probability(not no variant), and if ALT is, for example, A, this is -10log_10probability(not A). 

High QUAL scores indicate high confidence calls. Although traditionally people use integer phred scores, this field is permitted to be a floating point, to enable higher resolution for low confidence calls if desired. If unknown, the missing value should be specified.

We can add a tag (_i.e.,_ a bit of text) in the FILTER field to indicate that our variant is potentially poor quality. This is really important to know if we are going to report on the variant, which again can later inform treatment - we have to know if it is a true variant call, or maybe just noise. 

Put the text "LowQual" to the FILTER tag when QUAL<30, and print the result to the screen:

```bash
bcftools filter -m x -s LowQual -e 'QUAL<30' trio.trim.vcf.gz
```

Do that again but look at just the low quality variants that we just annotated
```bash
bcftools filter -m x -s LowQual -e 'QUAL<30' trio.trim.vcf.gz | grep 'LowQual'
```

---
**Question:**
---

1. Have a look at the options we just used (-m x, -s, -e) and their explanations by running
```bash
bcftools filter --help
```

2. Variants that are located close to indels can also indicate poor quality calls, so:
Use the `bcftools filter` sub-command to tag variants that are within 10 base-pairs of an InDel, in addition to the LowQual ones we just tagged with QUAL<30 (hint: check out the SnpGap option from the help menu to help you build your command)

---

<details>
<summary>Answers</summary>

```bash
bcftools filter -m x -g 10 -s LowQual -e 'QUAL<30' trio.trim.vcf.gz 
```
You can pipe into grep to look at the SnpGap locations:

```bash
bcftools filter -m x -g 10 -s LowQual -e 'QUAL<30' trio.trim.vcf.gz | grep SnpGap
```
You can further pipe that into wc -l to count up how many variants there are like this:

```bash
bcftools filter -m x -g 10 -s LowQual -e 'QUAL<30' trio.trim.vcf.gz | grep SnpGap | wc -l
```


</details>


### 2.5 Adding a variant ID

As you can probably imagine, each variant within the current reference genome has been extensively studied through the continual sampling of patients and individuals from around the world.
Due to this, each new variant identified is given an rsID in the [NCBI dbSNP database](https://www.ncbi.nlm.nih.gov/snp/).
These rsIds are helpful because it provides an extensive list of information about each particular variant.

(This is actually not _technically_ correct anymore, as we have sampled so many genomes recently that the rsID numbers couldn't keep up! But rsIDs can still be useful )

There are now databases such as the [genome aggregation database (gnomAD)](https://gnomad.broadinstitute.org/) that samples over 730,947 exomes and 76,215 genomes (gnomad v4, on hg38).

---
**Questions:**
---

Gnomad is helpful in knowing variant frequencies, and in finding predictive loss-of-function (pLoF) variants.
1. Search gnomAD for the gene SATB1. What are the numbers of observed and expected predicted loss of function (pLoF) variants in SATB1? How could these numbers help us determine whether this gene has a very important functional role?
2. Look up the variant `(GRCh38) 14-82099033-G-C` (rs75115269) and identify its allele frequency in African/African American population.
3. Can you find its equivalent genomics coordinates and allele frequency for the older hg19 reference genome?

---
<details>
<summary>Answers</summary>
1. 	Expected: 86.9	Observed: 11. Since the observed number of plof variants is a lot lower than expected, it likely means this gene is under pretty strong purifying/negative selection - i.e., plof variants don't want to accumulate here, because the gene is really important to survival. 
2. 0.03135
3. SNV: 14-82565377 G-C (GRCh37)

</details>

OK, back to annotating IDs.
If we have a database of known IDs, we can easily compare to the VCF file and add text to the ID field in the VCF (3rd field after CHROM and POS).
For this we can use any type of tab-delimited file, but for this week I have provided the dbSNP reference VCF which is perfect for this task.

First, let's index this dbSNP reference VCF, to help the lookup go faster:
```bash
bcftools index -t hg19.dbSNP.vcf.gz
```

Now we can add rsIDs using the `bcftools annotation` sub-command and output a new file with our IDs attached:

```bash
bcftools annotate -c CHROM,FROM,ID,REF,ALT \
    -a hg19.dbSNP.vcf.gz \
    -Oz -o trio.trim.dbSNP.vcf.gz trio.trim.vcf.gz
```

---
**Questions:**
---

1. Is the variant we explored in gnomAD (rs75115269) present in the hg19.dbSNP.vcf.gz?
2. In the output, trio.trim.dbSNP.vcf.gz, can you find the variant rs191680234, get its genomic coordinates, and identify for which individual it is variable?

---

<details>
<summary>Answers</summary>
1. Yes, we can see it at chr14 82099033, as expected. You can use either of the following to find it:

```bash
zgrep rs75115269 hg19.dbSNP.vcf.gz

bcftools view -i 'ID=="rs75115269"' hg19.dbSNP.vcf.gz
```


2. zgrep rs191680234 trio.trim.dbSNP.vcf.gz shows us the following line:
2	71797762	rs191680234	G	A	440.24	.	.	GT	0/0	0/1	0/0
We can see the second sample, 1847, has a heterozygous genotype 0/1. 

</details>


## 3. **Full variant annotation**

At the moment we have been sequentially adding layers of information in order to give context to each of the variants that we identified in the variant calling process.

To save time and space, we downloaded a fully annotated version of our VCF at the start of the tutorial which contains annotations of each variant as executed through Variant Effect Predictor (VEP).
VEP and [snpEff](http://snpeff.sourceforge.net/) are very popular variant annotation programs and use very similar effect standards.

As the name suggests, VEP uses large variant annotation databases to assign an effect and/or consequence to the particular change, which take the form of Sequence Ontology (SO) terms.

![VEP Variant Types](https://m.ensembl.org/info/genome/variation/prediction/consequences.svg)

As you can see from the table contained in the following link, each variant type is assigned a [level of variant _impact_ that is associated to the change](https://m.ensembl.org/info/genome/variation/prediction/predicted_data.html).
Some variants have a **HIGH** impact, including splice site changes, stop codon introductions and frameshift variants.
These we can interpret as potentially having a loss of function and therefore could influence phenotypes in a negative way.
**MODERATE** impact variants include missense variants and inframe deletions or insertions, that can have a significant change but can also lead to harmless changes.
**LOW** impact variants such as synonymous variants are likely to have a very low impact on overall the function of the protein.
Additional to protein-coding changes, non-coding or regulatory variant sequence ontologies are coded as **MODIFIERS**, although some non-coding regulatory sequence changes can produce a more significant impact.

So how are these included in the actual VCF file? Let's look:

```bash
# View the VCF that contains full annotation (without the header _i.e.,_ -H param)
bcftools view -H trio.trim.vep.vcf.gz | head

2	41647	.	A	G	4495.41	PASS	CSQ=intron_variant&non_coding_transcript_variant|||ENSG00000184731|FAM110C|ENST00000460464|||||processed_transcript|||||||||,intron_variant&non_coding_transcript_variant|||ENSG00000184731|FAM110C|ENST00000461026|||||processed_transcript|||||||||,intron_variant|||ENSG00000184731|FAM110C|ENST00000327669||||-/321|protein_coding|YES|CCDS42645.1|||||||	GT:AD:DP:GQ:PL	0/0:56,0:56:99:0,169,2183	0/1:33,35:68:99:1139,0,1044	0/1:119,117:237:99:3356,0,3283
2	45895	.	A	G	463.75	PASS	CSQ=missense_variant|aTc/aCc|I/T|ENSG00000184731|FAM110C|ENST00000327669|1/2|benign(0)|tolerated(0.62)|164/321|protein_coding|YES|CCDS42645.1|||||||,upstream_gene_variant|||ENSG00000184731|FAM110C|ENST00000460464|||||processed_transcript|||||||||,intron_variant&non_coding_transcript_variant|||ENSG00000184731|FAM110C|ENST00000461026|||||processed_transcript|||||||||	GT:AD:DP:GQ:PL	1/1:0,6:6:18.05:207,18,0	1/1:0,9:9:24.07:292,24,0	./.:.:.:.:.
2	224970	.	C	T	4241.64	PASS	CSQ=intron_variant|||ENSG00000035115|SH3YL1|ENST00000415006||||-/246|protein_coding||CCDS62842.1|||||||,intron_variant|||ENSG00000035115|SH3YL1|ENST00000403657||||-/227|protein_coding||CCDS62841.1|||||||,intron_variant|||ENSG00000035115|SH3YL1|ENST00000403658||||-/227|protein_coding||CCDS62841.1|||||||,intron_variant|||ENSG00000035115|SH3YL1|ENST00000405430||||-/342|protein_coding|||||||||,intron_variant&non_coding_transcript_variant|||ENSG00000035115|SH3YL1|ENST00000473104|||||processed_transcript|||||||||,intron_variant|||ENSG00000035115|SH3YL1|ENST00000451005||||-/255|protein_coding|||||||||,intron_variant|||ENSG00000035115|SH3YL1|ENST00000356150||||-/342|protein_coding|YES|CCDS42646.2|||||||,intron_variant&NMD_transcript_variant|||ENSG00000035115|SH3YL1|ENST00000479739||||-/155|nonsense_mediated_decay|||||||||,intron_variant&non_coding_transcript_variant|||ENSG00000035115|SH3YL1|ENST00000463865|||||processed_transcript|||||||||,intron_variant&non_coding_transcript_variant|||ENSG00000035115|SH3YL1|ENST00000472012|||||processed_transcript|||||||||,downstream_gene_variant|||ENSG00000035115|SH3YL1|ENST00000431160||||-/230|protein_coding|||||||||,intron_variant|||ENSG00000035115|SH3YL1|ENST00000403712||||-/323|protein_coding||CCDS54332.1|||||||,intron_variant&non_coding_transcript_variant|||ENSG00000035115|SH3YL1|ENST00000468321|||||processed_transcript|||||||||GT:AD:DP:GQ:PL	0/1:40,26:66:99:789,0,1374	0/1:47,41:88:99:1247,0,1555	0/1:93,80:175:99:2205,0,2918
...
```

My eyes!
So.....much......text......
As you can see, there is a mass of information in the INFO field, all of which starts with a consequence tag (CSQ=).
This field has a lot of information separated by pipes (|) and it is also possible to get multiple annotations per variant.
If you look at the header you can get the header information for each of these fields that are separated by |.

```bash
zgrep "^##INFO=<ID=CSQ" trio.trim.vep.vcf.gz

##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: Consequence|Codons|Amino_acids|Gene|SYMBOL|Feature|EXON|PolyPhen|SIFT|Protein_position|BIOTYPE|CANONICAL|CCDS|RadialSVM_score|RadialSVM_pred|LR_score|LR_pred|CADD_raw|CADD_phred|Reliability_index">
```

---
**Questions:**
---

For these questions, you may need to refer to the [bcftools man page:](https://samtools.github.io/bcftools/bcftools.html#expressions)



1. How many variants have a reported "missense_variant"? (TIP: You need to build a 'bcftools view' command - you will need to add the '--include' parameter; -H to get rid of the header; and the expression you want to test will be 'INFO/CSQ~"missense_variant"'). Your command should look something like this: bcftools view --include {'expression'} -H {vcf}
2. How many variants that are greater than quality 30 have a reported "missense_variant"? (Tip: the expression you need is 'QUAL>30', and you can string it together with your previous command like so: bcftools view --include '{expression} & {expression}' -H {vcf})
3. How many variants are annotated to have gained a stop codon? (Tip: the expression you need is 'INFO/CSQ~"stop_gained"')

---


<details>
<summary>Answers</summary>
1. The command to run is below, then piped into wc -l to count lines:

```bash
 bcftools view --include 'INFO/CSQ~"missense_variant"' -H trio.trim.vep.vcf.gz | wc -l 
 ```

The number of missense variants is 3287
2. 

```bash
bcftools view --include 'INFO/CSQ~"missense_variant" & QUAL>30' -H trio.trim.vep.vcf.gz  | wc -l
```
The number of missense variants with >30 quality is 3192

3. 
```bash
bcftools view --include 'INFO/CSQ~"stop_gained"' -H trio.trim.vep.vcf.gz  | wc -l
```

Number of stop gained is 47 

</details>

### 3.1 Variant scores

Of course, the impact of the variant may not always tell you the definitive pathogenicity of that variant.
LOW impact variants may also have a high impact within some systems, so a number of additional metrics have been established to add additional interpretation power to variant annotation.
These include:

- GERP
- CADD
- SIFT
- PolyPhen
- REVEL

[All the additional fields are here](https://m.ensembl.org/info/genome/variation/prediction/protein_function.html)

We don't have time to go into all of these tools but to give you a flavour, the Combined Annotation Dependent Depletion (CADD) tool scores the predicted deleteriousness of single nucleotide variants and insertion/deletions variants in the human genome by integrating multiple annotations including conservation and functional information into one metric.
Phred-style CADD raw scores are displayed and variants with higher scores are more likely to be deleterious.
CADD provides a ranking rather than a prediction or default cut-off, with higher scores more likely to be deleterious.
For example, scores above 30 are 'likely deleterious' and scores below as 'likely benign'.
Variants with scores over 30 are predicted to be the 0.1% most deleterious possible substitutions in the human genome.

### Conclusion - why annotate?


The whole purpose of variant annotation is to throw every piece of evidence we have that a variant *does something* functionally on a molecular, cellular, or even physiolological scale into our vcf. Back to our earring analogy - it's like building a map of all the places we could have lost our earring on our roadtrip, like rest areas, petrol stations, site-seeing areas. This is basically what we've covered today. 

The next process after this is to then filter out variants based on all those annotations down to a handful of candidates - and then, if the evidence is strong enough for one of them, we report that as our diagnostic variant. This process of elimination, using bits of logic to infer whether a variant could or could not be causing the observed phenotype - is called variant classification and curation. We will look at this more next prac!


