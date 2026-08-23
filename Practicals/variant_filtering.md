# Clinical Genomics Practical: Part 2 Variant Filtering and Pedigree Analyses



#### By Evelyn Collen
(prac repurposed from the excellent Julien Soubrier)

## **1. Some background**

Last prac, we annotated a bunch of information into a vcf and learnt about the annotating process. Variant filtering (also called variant priorisation) then gets us down to a handful of candidates - and then, if the evidence is strong enough for one of them, we can report that as our diagnostic variant in the variant classification and curation step.

Some single gene and well-known Mendelian genetic disorders, such as sickle-cell anemia, Tay–Sachs disease and cystic fibrosis, can be relatively straightfoward to diagnose. However, much of the time, it's a lot more complex than a single gene with a well-characterised mutation, and the diagnostic variant we are searching for has be to sifted out 3 billion base pairs' worth of possible variation. 


Today we will be looking at a variation filtering and curation method using family inheritance patterns. Modes of inheritance can provide some good logic to help us match the patterns of variant inheritance to the patterns of phenotype inheritance that we observe in the patient's family. For example, does mum or dad are carry the condition? Is it a _de novo_ mutation, that may have arisen independently of the parents?


As you may have noticed, the three samples sequenced in our data are related to each other, and form a "trio" (mother-father-son). Trios are a very standard approach in clinical genetics, especially common for prenatal testing.

### 1.1 Reminder about virtual Machines

As usual we will be connecting the virtual machines: 

**Please [go here](../../Course_materials/vm_login_instructions.md) for instructions on connecting to your VM.**

### 1.3 Learning Outcomes

CHANGEME
1. Get comfy manipulating annotation info in a vcf file
2. Understand what a monumental task it is trying to find diagnostic variants
3. Refresher on vcf file specifications
4. Learn about common variant annotations and annotation databases
5. Get an idea of common annotations and the logic behind why we include them




### 1.3 This week's tutorial

This week's tutorial is liberally taken from two tutorials written by Aaron Quinlan & his group at University of Utah, and introduces us the variant priorisation software 'Gemini'.

- [Identifying dominant gene candidates with GEMINI](https://s3.amazonaws.com/gemini-tutorials/Gemini-Dominant-Tutorial.pdf)
- [Identifying recessive gene candidates with GEMINI](https://s3.amazonaws.com/gemini-tutorials/Gemini-Recessive-Tutorial.pdf)

There are a lot of variant prioritisation programs similar to this one including Variantgrid, Emedgene, Franklin, Nostos and many others. 
[Gemini](https://gemini.readthedocs.io/en/latest/) is a database system that can read in VCF information and family/pedigree information, to enable database querying and clinical genetics analyses.
Information in gemini is stored in database system called SQL.
SQL stands for Structured Query Language, and is a SUPER popular database system in many industries.
It comes in many flavours that you might have heard before, including `MySQL`, `SQLite` and `PostgreSQL`.

Let's create and activate a conda environment with Gemini installed:
```bash
# First go to the working directory we created yesterday:
cd ~/clinical_genomics
# Activate a conda environment with Gemini installed:
source activate geminiEnv
# Let's make sure we can now use Gemini:
gemini -h
```

### 1.3 Cohort databases

Lets make some databases!
Gemini can take the VCF file and sample information in the form of a ped file (short for pedigree).
The ped file is actually a standard metadata information file that was developed for in the age of population genetics and GWAS analyses.


Unfortunately the database loading command also adds a lot of annotation information with memory requirements, so let's use a pre-generated database: ___trio.trim.vep.dominant.db___

The command you would have run is here for your info:

```bash
# Loading VCF files require all the annotation databases
## Cmd to create the gemini db for dominant study
#gemini load --cores 4 -v trio.trim.vep.vcf.gz -t VEP \
#        --skip-gene-tables -p dominant.ped trio.trim.vep.dominant.db
```

We'll use this database for our querying and `autosomal_dominant` analysis.

### 1.4 Querying a SQL database in Gemini

First - here's some quick points on SQL and examples on how to use SQL queries. 

Information is stored in tables, much like a sheet within a Microsoft Excel Spreadsheet file.
You can ask the database to give you specific, tabular information and also put restrictions on the type of information that is needed (such as a conditional like "the sex of the person must be female").

If you are interested, you can find a lot more info on [standard SQL commands](https://www.codecademy.com/articles/sql-commands).

For this prac, we only need to know how to query the database using the `SELECT` and `FROM` command.
We will also need some conditionals to  like `WHERE`, `IS`, `IS NOT` etc. 

Let's have a look at what we have in the database file:

```bash
gemini db_info trio.trim.vep.dominant.db
```

As you can see, we have a number of tables. Within those tables we have columns.
And the data within those columns has a specific data type, all of which we can query.

Let's look at some examples by using the `gemini query` sub-command.

```bash
# Get everything from the samples table with a wildcard
gemini query -q "SELECT * FROM samples" trio.trim.vep.dominant.db

# Get the names of the individuals from the samples table
gemini query -q "SELECT name FROM samples" trio.trim.vep.dominant.db

# Get the names of the individuals from the samples table that are male
gemini query -q "SELECT name FROM samples WHERE sex IS 1" trio.trim.vep.dominant.db

# Get the names of the individuals from the samples table that are male
gemini query -q "SELECT name FROM samples WHERE sex == 1" trio.trim.vep.dominant.db

# Get the names of the individuals from the samples table that are male
gemini query -q "SELECT name FROM samples WHERE sex IS NOT 2" trio.trim.vep.dominant.db

# Get the names of the individuals from the samples table that are male
gemini query -q "SELECT name, sex FROM samples WHERE sex IS NOT 2" trio.trim.vep.dominant.db
```

Note: Depending on the data type, you may need to surround character info in ''. 


---
___>>> TASK <<<___

---

#### Build your query
Now that we known the query structure and tables that we have in our database, construct some more sophisticated queries.
1. Extract the chromosome and position of all variants in the database ("FROM variants") that have a 1000 genome allele frequency in Europeans (aaf_1kg_eur) less than 0.5. How many are there?
2. Extract all the variants within the genes MAPK12 that have a variant quality > 200. How many are there? How many are also QUAL > 500?
3. How many variants that were marked as "PASS" quality were concordant?

---

CHANGEME
<details>
<summary>Answers</summary>
1. 	Expected: 86.9	Observed: 11
2. 0.03135
3. SNV: 14-82565377 G-C (GRCh37)

</details>


## **2. Autosomal Dominant disorders**

Autosomal dominant disorders are genetic disorders that do not involve the sex chromosomes and are passed down through families in a vertical transmission pattern. 
Incomplete penetrance can occur within the family, meaning that disorder may not always show phenotypically.
`Penetrance` here means "the extent to which a particular gene or set of genes is expressed in the phenotypes of individuals carrying it, measured by the proportion of carriers showing the characteristic phenotype."

![https://www.mayoclinic.org/autosomal-dominant-inheritance-pattern/img-20006210](https://www.mayoclinic.org/-/media/kcms/gbs/patient-consumer/images/2013/11/15/17/37/r7_autosomaldominantthu_jpg.jpg)

Examples of these disorders include Huntington disease, neurofibromatosis, and polycystic kidney disease.

The trio that we will be looking at today is affected by a rare disorder called [_hypobetalipoproteinemia_](https://ghr.nlm.nih.gov/condition/familial-hypobetalipoproteinemia), a disorder that impairs the body's ability to absorb and transport fats.
It occurs in about in 1 in 1,000 to 3,000 individuals (1 in ~1,000 in Europeans).

![Trio pattern](images/family.png)

Both the mother (`1805`) and the son (`4805`) are affected with the disorder. Although they have normal plasma HDL, they suffer from fat malabsorption and are in the bottom 5% for plasma cholesterol and triglycerides.

We want to be able to see these sort of relationships in the PED file so lets have a quick look:

```bash
cat dominant.ped | column -t

#family_id  sample_id  paternal_id  maternal_id  sex  phenotype  ancestry
family1     1805       -9           -9           2    2          CEU
family1     1847       -9           -9           1    1          CEU
family1     4805       1847         1805         1    2          CEU
```

As you can see, all of the relationships between the individuals are recorded in the PED file, including the sex of the individuals and the prevalence of the phenotype. -9 means NA here.
When sampling larger families, or even populations, all of the unique inheritance patterns can be recorded here and used to inform the clinical model when it comes to identifying candidate genes or variants for the disorder.

### 2.1 Sample genotype queries

Given that we will be comparing the pattern of inheritance of these variants, its easy to filter variants so that you pick up specific relationships between individuals.
For example, what if we wanted to identify variants where both 1805 (mum) and 4805 (son) have a non-reference allele?
(After all, 1805 and 4805 are both affected, so this seems likely!)

```bash
# Show all info
gemini query -q "SELECT * from variants" \
            --gt-filter "(gt_types.1805 <> HOM_REF AND gt_types.4805 <> HOM_REF)" \
            --header \
            trio.trim.vep.dominant.db

# Print just the genotypes to compare
gemini query -q "SELECT gts.1805, gts.4805 from variants" \
            --gt-filter "(gt_types.1805 <> HOM_REF AND gt_types.4805 <> HOM_REF)" \
            --header \
            trio.trim.vep.dominant.db
```

Or how about using a wildcard with the `--gt-filter` to identify all heterozygous variants.
The syntax for wildcards in `--gt-filter` follows a slightly different format:

```--gt-filter (COLUMN).(SAMPLE_WILDCARD).(SAMPLE_WILDCARD_RULE).(RULE_ENFORCEMENT)```

```bash
# Print heterozygous variants in all
gemini query -q "SELECT chrom, start, end, ref, alt, gene, impact, (gts).(*) \
                 FROM variants" \
            --gt-filter "(gt_types).(*).(==HET).(all)" \
            --header \
            trio.trim.vep.dominant.db

# Print variants where all the females variants are reference homozygous
gemini query -q "SELECT chrom, start, end, ref, alt, gene, impact, (gts).(*) \
                 FROM variants" \
            --gt-filter "(gt_types).(sex==2).(==HOM_REF).(all)" \
            --header \
            trio.trim.vep.dominant.db
```

Here you can add quality filters for each of the genotypes. 
You can look at variant depth and quality in each genotype by using the `gt_depth` and `gt_quals`

### 2.2 Using the `autosomal_dominant` tool

Gemini already comes with a number of tools that allow you to assess particular types of clinical genetic patterns, including one for autosomal dominant disorders.
This tool has the known characteristics of this type of genetic disorder hard-coded into the function.
The genotype requirements for an autosomal dominant disorder are:

1. All affected individuals must be a heterozygous variant. 
   - After all one copy comes from mum and the other from dad
   - One of which is defective
2. No unaffacted can be heterozygous or homozygous alternate
   - However they can be unknown
3. Can't be `de-novo` mutations (i.e., can't be any variants not seen in mum and/or dad)
4. Affected individual must have an affected parent
5. All affecteds must have parents where the phenotype is known
   - Otherwise throw a warning

Let's have a look at it:

```bash
gemini autosomal_dominant \
    --columns "chrom, start, end, ref, alt, gene, impact, cadd_raw" \
    trio.trim.vep.dominant.db | head | column -t
```

Here we are running the `autosomal_dominant` tool, and extracting specific columns of information from the database, printing only the first few lines and separating them out into tab delimited columns so we can see whats going on.
We want to include the important information like whether the variant is in a gene, what the impact of the variant is, and whats the pathogenicity of that variant (using the raw CADD score annotation, if you remember back from last prac).

From here we can continue the variant filtering process to a few candidates, then 'curate' them to the most likely diagnostic causal one.

---
___>>> TASK <<<___

---
#### Find candidate genes for Hypobetalipoproteinemia

1. Generate a list of the variants that have 'HIGH' and 'MODERATE' impact. How many do you have?
2. Build up some additional filters 
  - For example, the variant is likely to be rare in the European population, so it might be good to look at allele-frequency's < 0.01 in Europeans.
    - This can be either from gnomAD, ExAC or 1000genomes allele-frequencies (gnomAD generally has the most accurate frequencies)
3. Generate a list of candidate genes based on your data 
---

---
___>>> EXTRA TASK <<<___

---

#### Recessive example

1. Do the same exercise with the recessive database: ___trio.trim.vep.recessive.db___
2. What `gemini` functions can be employed to identify recessive variants?

---



```bash
gemini autosomal_recessive \
    --columns "chrom, start, end, ref, alt, gene, impact, cadd_raw" \
    trio.trim.vep.recessive.db | head | column -t
```

gemini query -q "SELECT * FROM variants WHERE (gts).0 == 1 AND (gts).1 == 1"