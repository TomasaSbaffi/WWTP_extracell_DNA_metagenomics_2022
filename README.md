# WWTP_extracell_DNA_metagenomics_2022

Bioinformatics workflow used in the article:

> Periyasamy S, Sabbatino R, Sbaffi T,... Corno G, Di Casare A. 2022. Title. ref. doi: xxxx.


## Contacts

**Tomasa Sbaffi**  
Postdoctoral Researcher  
[E-mail](mailto:tomasa.sbaffi@gmail.com)

**Andrea Di Cesare**  
Principal Investigator  
[E-mail](mailto:andrea.dicesare@cnr.it)


## Table of contents

1. [Before starting](#before-starting)
2. [Download genomes and import to anvi'o](#download-genomes-and-import-to-anvio)
3. [Annotation](#annotation)
4. [Phylogeny and phylogenomics](#phylogeny-and-phylogenomics)
5. [Abundance analysis](#abundance-analysis)


## Before starting

### You will need to have these softwares installed and added to your path

* SRA Toolkit v2.11.3: https://github.com/ncbi/sra-tools
* fasterq-dump v2.10.8: https://github.com/glarue/fasterq_dump
* GNU parallel: https://www.gnu.org/software/parallel
* BLAST v2.10.1: https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+
* seqtk v1.3: https://github.com/lh3/seqtk
* CoverM v0.6.1: https://github.com/wwood/CoverM/
* fastQC v0.11.9: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
* multiQC v1.8: https://multiqc.info/
* Cutadapt v1.16: https://cutadapt.readthedocs.io/
* MEGAHIT v1.1.2.9: https://github.com/voutcn/megahit/
* bowtie v2.3.5: http://bowtie-bio.sourceforge.net/bowtie2/
* SAMtools v1.9: http://www.htslib.org/
* GTDB-Tk v2.1.0 and GTDB release 05-RS95: https://gtdb.ecogenomic.org/
* MetaPhlan v3.0.14
* Deeparg v1.0.2
* Trimmomatic v0.39
* vsearch v2.17.1
* Trimgalore v0.6.5
* bwa-mem v0.7.17
* MetaBAT v2.15.0
* RefineM v0.1.2
* CheckM v1.2.0
* Magpurify v2.1.2
* prodigal v2.6.3
 
It sghould be possible running the analysis with any UNIX-based OS.

### Define number of threads to use

You should change below to the number of cores available in your system:

```bash
NTHREADS=40
```

# Pre-processing of raw metagenomic data

We have a total of 16 samples sequenced with Illumina.


### Check raw data with fastQC and multiQC

```bash
mkdir RAW_READS/FASTQC

fastqc RAW_READS/*.fastq.gz \
       --outdir RAW_READS/FASTQC \
       --threads $NTHREADS

multiqc RAW_READS/FASTQC \
        --outdir RAW_READS/MULTIQC \
        --interactive
```


### Run the DeepARG SS pipeline

This pipeline includes quality-filtering and PE reads merging and outputs a .clean fasta file prior annotating for Antimicrobial Resistance Genes (ARGs).

```bash
mkdir deeparg_SS_out
mkdir LOGS

# Soft links to the raw data
ln -s RAW_READS/*fastq.gz .

db_dARG_path=DB_PATH
samples=`ls *fastq.gz | awk '{split($0,x,".r"); print x[1]}' | sort | uniq`

for sample in ${samples} ; do
deeparg short_reads_pipeline \
    --forward_pe_file ${sample}.r1.fastq.gz \
    --reverse_pe_file ${sample}.r2.fastq.gz \
    --output_file deeparg_SS_out/${sample} \
    -d ${db_dARG_path} \
    --deeparg_identity 90 \
    --deeparg_evalue 1e-10 \
    --deeparg_probability 0.8 &> LOGS/${sample}.log.txt
done        
```

Log files include info about quality check single steps outcomes and 16S rRNA gene hit counts (Greengenes db) that the pipeline uses to normalise ARGs abundancies.



