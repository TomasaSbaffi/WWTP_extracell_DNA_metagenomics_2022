# WWTP Intra_extra_cellular DNA metagenomics 2022

Bioinformatics workflow used in the article:

> Periyasamy S, Sabbatino R, Sbaffi T,... Corno G, Di Casare A. 2022. Title. ref. doi: xxxx.

For the biostatistics workflow of the same paper, [Click here](Intra_Extra_DNA_script_statistical_analysis.R)


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

You should change below to the number of cores available in your system and set you project name:

```bash
NTHREADS=40
project=THIS_PROJECT
```

# Read based analysis

## Pre-processing of raw metagenomic data and deepARG

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

# Retrieving ARGs normalised pivot table
pivots=`ls *.merged.quant.subtype | ls *merged.quant.subtype | awk '{printf "%s%s", sep, $0; sep=","} END{print ""}'`
samplenames=$(basename -s .clean.deeparg.mapping.ARG.merged.quant.subtype -a *.merged.quant.subtype | awk '{printf "%s%s", sep, $0; sep=","} END{print ""}')

echo "`date +"%Y%m%d %H:%M:%S"`: merging $pivots with genetools"

genetools deeparg-table \
  --deeparg-files "$pivots" \
  --sample-names "$samplenames" \
  --output-file "${project}_pivot.tsv"
```

Log files include info about quality check single steps outcomes and 16S rRNA gene hit counts (Greengenes db) that the pipeline uses to normalise ARGs abundancies.


## Taxonomic assessment
 
```bash
# Soft links to the cleaned data
ln -s deeparg_SS_out/*.clean .

samples=`ls *.clean | awk '{split($0,x,".clean"); print x[1]}' | sort | uniq`

for sample in $samples; do
metaphlan ${sample}.clean \
   --input_type fasta \
   --bowtie2out ${sample}.bowtie2.bz2 \
   --nproc $NTHREADS \
   --stat_q 0.1 > ${sample}_profile.txt  # i deleted --ignore_eukaryotes --ignore_archaea 
done

# Retrieving percentage normalised taxonomic table, i.e. genus level

samplenames=$(basename -s _profile.txt -a *_profile.txt | awk '{printf "%s%s", sep, $0; sep="_"} END{print ""}')
merge_metaphlan_tables.py *_profile.txt > metaphlan3_${project}_pivot.txt

grep -E "g__|clade" metaphlan3_${project}_pivot.txt | sed 's/^.*g__//g' | grep -v "|s__" | cut -f1,3-50 | sed -e 's/clade_name/genus_name/g' | sed -e 's/_profile//g' > metaphlan3_${project}_pivot_genera.txt

grep "g__" metaphlan3_${project}_pivot.txt | grep -v "|s__" | cut -f1,1 | wc -l > ${project}_number_genera.txt
```


# Assembly based analysis

## Pre-processing 

```bash
# Soft links to the raw data
ln -s RAW_READS/*fastq.gz .

mkdir TRIMGALORE
mkdir FLASH

samples=`ls *fastq.gz | awk '{split($0,x,".r"); print x[1]}' | sort | uniq`

for sample in ${samples} ; do
mkdir TRIMGALORE/${sample}
trim_galore --paired ${sample}.r1.fastq.gz ${sample}.r2.fastq.gz \
  -q 20 \
  -o TRIMGALORE/${sample}/
done

for sample in ${samples} ; do
flash2 TRIMGALORE/${sample}/${sample}.r1_val_1.fq.gz TRIMGALORE/${sample}/${sample}.r2_val_2.fq.gz \
  -d FLASH/${sample}/ \
  -M 170 \
  -z \
  -t $NTHREADS
done
```

### Assemble reads with MEGAHIT
Based on the taxonomic profiles and samples origins we opted towards eight co-assemblies. Here the code for one co-assembly is provided.
Small co-assemblies were produced, the first includes the "IVB_pre" two samples.

```bash
mkdir ASSEMBLIES
        
megahit --presets meta-sensitive \
        -1 FLASH/IVB_pre_day1/out.notCombined_1.fastq.gz,FLASH/IVB_pre_day2/out.notCombined_1.fastq.gz \
        -2 FLASH/IVB_pre_day1/out.notCombined_2.fastq.gz,FLASH/IVB_pre_day2/out.notCombined_2.fastq.gz \
        -r FLASH/IVB_pre_day1/out.extendedFrags.fastq.gz,FLASH/IVB_pre_day2/out.extendedFrags.fastq.gz \
        --min-contig-len 1000\
        -t $NTHREADS \
        --out-dir IVB_pre_m1000    
```


```bash
# Indexing 
mkdir index
bwa index -p index/IVB_pre ../IVB_pre_m1000/final.contigs.fa

# Mapping contigs for all the samples included in the coassembly
mkdir SAMFILES/IVB_pre/
bwa mem index/IVB_pre /FLASH/IVB_pre_day1/out.extendedFrags.fastq.gz \
  -t $NTHREADS > SAMFILES/IVB_pre/IVB_pre_day1_se.sam
bwa mem index/IVB_pre /FLASH/IVB_pre_day1/out.notCombined_1.fastq.gz /FLASH/IVB_pre_day1/out.notCombined_2.fastq.gz \
  -t $NTHREADS > SAMFILES/IVB_pre/IVB_pre_day1_pe.sam
  
bwa mem index/IVB_pre /FLASH/IVB_pre_day2/out.extendedFrags.fastq.gz \
  -t $NTHREADS > SAMFILES/IVB_pre/IVB_pre_day2_se.sam
bwa mem index/IVB_pre /FLASH/IVB_pre_day2/out.notCombined_1.fastq.gz /FLASH/IVB_pre_day2/out.notCombined_2.fastq.gz \
  -t $NTHREADS > SAMFILES/IVB_pre/IVB_pre_day2_pe.sam

# sam2bam conversion
mkdir BAMFILES
samtools view -@ 40 -b -h SAMFILES/IVB_pre/IVB_pre_day1_pe.sam -o BAMFILES/IVB_pre/IVB_pre_day1_pe.bam
samtools view -@ 40 -b -h SAMFILES/IVB_pre/IVB_pre_day1_se.sam -o BAMFILES/IVB_pre/IVB_pre_day1_se.bam
samtools view -@ 40 -b -h SAMFILES/IVB_pre/IVB_pre_day2_pe.sam -o BAMFILES/IVB_pre/IVB_pre_day2_pe.bam
samtools view -@ 40 -b -h SAMFILES/IVB_pre/IVB_pre_day2_se.sam -o BAMFILES/IVB_pre/IVB_pre_day2_se.bam

#samtools sort
samtools sort -@ 40 -o BAMFILES/IVB_pre/sorted/IVB_pre_day1_pe_sorted.bam BAMFILES/IVB_pre/IVB_pre_day1_pe.bam
samtools sort -@ 40 -o BAMFILES/IVB_pre/sorted/IVB_pre_day1_se_sorted.bam BAMFILES/IVB_pre/IVB_pre_day1_se.bam
samtools sort -@ 40 -o BAMFILES/IVB_pre/sorted/IVB_pre_day2_pe_sorted.bam BAMFILES/IVB_pre/IVB_pre_day2_pe.bam
samtools sort -@ 40 -o BAMFILES/IVB_pre/sorted/IVB_pre_day2_se_sorted.bam BAMFILES/IVB_pre/IVB_pre_day2_se.bam

# Indexing samtools files
samtools index BAMFILES/IVB_pre/sorted/IVB_pre_day1_pe_sorted.bam
samtools index BAMFILES/IVB_pre/sorted/IVB_pre_day1_se_sorted.bam
samtools index BAMFILES/IVB_pre/sorted/IVB_pre_day2_pe_sorted.bam
samtools index BAMFILES/IVB_pre/sorted/IVB_pre_day2_se_sorted.bam
```

