# WWTP Intra_extra_cellular DNA metagenomics 2022

Bioinformatics workflow used in the article:

> Periyasamy S, Sabbatino R, Sbaffi T, Fontaneto D, Corno G, Di Casare A. 2022.  Extracellular DNA includes an important fraction of high-risk antibiotic resistance genes in treated wastewaters. Environmental Pollution, under review. doi: xxxx.

For the biostatistics workflow of the same paper, [Click here](Intra_Extra_DNA_script_statistical_analysis.R)


## Contacts

**Tomasa Sbaffi**  
Postdoctoral Researcher, bioinformatics  
[E-mail](mailto:tomasa.sbaffi@gmail.com)

**Andrea Di Cesare**  
Principal Investigator  
[E-mail](mailto:andrea.dicesare@cnr.it)


## Table of contents

1. [Before starting](#before-starting)
2. [Read based analysis](#read-based-analysis)
3. [Assembly based analysis](#assembly-based-analysis)


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
 
It should be possible running the analysis with any UNIX-based OS.

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

### Binning contigs with Metabat2

```bash
# Making the contig depth file (Make the output dir/ before running)
jgi_summarize_bam_contig_depths --outputDepth IVB_pre_m1000/depth.txt /BAMFILES/IVB_pre/sorted/*bam

# metabat using -m 1500
metabat2 --saveCls -m 1500 -v -t $NTHREADS -i /IVB_pre_m1000/final.contigs.fa -o IVB_pre_m1000/metabat2_1500/bin -a IVB_pre_m1000/depth.txt

# metabat using -m 2000
metabat2 --saveCls -m 2000 -v -t $NTHREADS -i /IVB_pre_m1000/final.contigs.fa -o IVB_pre_m1000/metabat2_2000/bin -a IVB_pre_m1000/depth.txt

# metabat using -m 2500
metabat2 --saveCls -m 2500 -v -t $NTHREADS -i /IVB_pre_m1000/final.contigs.fa -o IVB_pre_m1000/metabat2_2500/bin -a IVB_pre_m1000/depth.txt
```

### Refining bins with refineM ...
Here I show the refining and following analysis for the bins recovered with metabat -m 1500 (min lenght 1500), however we repeated the procedure with -m 2000 and -m 2500 output'ed bins and subsequently annotated the pool recovering the highest amount high quality MAGs.

```bash

DATABASES=YOUR_DATABASES_PATH

# REFINE 1500, PART I (refineM)

# calculate contigs metrics in each bin # ident. contigs with divergent genomic properties # removing contigs
refinem scaffold_stats -r -x fa -c $NTHREADS \
  IVB_pre_metasens_m1000/final.contigs.fa  \
  IVB_pre/metabat2_1500/ \
  IVB_pre/metabat2_1500/scaffold_stats/ \
  BAMFILES/IVB_pre/sorted/*.bam

refinem outliers IVB_pre/metabat2_1500/scaffold_stats/scaffold_stats.tsv \
  IVB_pre/metabat2_1500/outliers/
  
refinem filter_bins -x fa IVB_pre/metabat2_1500/ \
  IVB_pre/metabat2_1500/outliers/outliers.tsv \
  IVB_pre/metabat2_1500/filtered_bins/

## identify contigs with divergent tax assignments # predicting genes in contigs # DIAMOND # identifying the contigs # filtering bins
refinem call_genes -c $NTHREADS -x fa \
  IVB_pre/metabat2_1500/filtered_bins/ \
  IVB_pre/metabat2_1500/called_genes/
  
refinem taxon_profile -c $NTHREADS \
  IVB_pre/metabat2_1500/called_genes/ \
  IVB_pre/metabat2_1500/scaffold_stats/scaffold_stats.tsv \
  /$DATABASES/refinem/diamond_proteinDB/gtdb_r95_protein_db.2020-07-30.faa.dmnd \
  /$DATABASES/refinem/taxonomy/gtdb_r95_taxonomy.2020-07-30.tsv \
  IVB_pre/metabat2_1500/taxon_profile/

refinem taxon_filter -c $NTHREADS \
  IVB_pre/metabat2_1500/taxon_profile/ \
  IVB_pre/metabat2_1500/taxon_profile/taxon_filter.tsv

refinem filter_bins -x fa \
  IVB_pre/metabat2_1500/filtered_bins/ \
  IVB_pre/metabat2_1500/taxon_profile/taxon_filter.tsv \
  IVB_pre/metabat2_1500/filtered_bins_part2/
```

### ... and MAGpurify

```bash
# REFINE 1500, PART II (MAGpurify)

cd IVB_pre/metabat2_1500/filtered_bins_part2/
mkdir filtered_bins_final
for bin in *.fa
do
        SAMPLE=$(echo ${bin} | sed "s/.filtered.filtered.fa//")
        echo $SAMPLE
        magpurify phylo-markers $bin output/${SAMPLE} --threads $NTHREADS --db /$DATABASES/MAGpurify-db-v1.0/
        magpurify clade-markers $bin output/${SAMPLE} --db /$DATABASES/MAGpurify-db-v1.0/
        magpurify tetra-freq $bin output/${SAMPLE}
        magpurify gc-content $bin output/${SAMPLE}
        magpurify known-contam $bin output/${SAMPLE} --threads $NTHREADS --db /$DATABASES/MAGpurify-db-v1.0/
        magpurify clean-bin $bin output/${SAMPLE} filtered_bins_final/${SAMPLE}_cleaned.fna
done
#
```

### Quality checking refined bins

```bash
## using checkM to discover lineage for refined bins
checkm lineage_wf -f IVB_pre/metabat2_1500/filtered_bins_part2/filtered_bins_final/CheckM.txt \
  -t $NTHREADS \
  -x fna IVB_pre/metabat2_1500/filtered_bins_part2/filtered_bins_final/ \
  IVB_pre/metabat2_1500/filtered_bins_part2/filtered_bins_final/SCG/

## analyze final bins with checkM
checkm analyze --ali --nt -t $NTHREADS -x fna \
  IVB_pre/metabat2_1500/filtered_bins_part2/filtered_bins_final/SCG/lineage.ms \
  IVB_pre/metabat2_1500/filtered_bins_part2/filtered_bins_final/ \
  IVB_pre/metabat2_1500/filtered_bins_part2/filtered_bins_final/analyze_out/

## full qa with checkM
checkm qa -t $NTHREADS \
  IVB_pre/metabat2_1500/filtered_bins_part2/filtered_bins_final/SCG/lineage.ms \
  IVB_pre/metabat2_1500/filtered_bins_part2/filtered_bins_final/analyze_out/ \
  --tab_table \
  -f IVB_pre/metabat2_1500/filtered_bins_part2/filtered_bins_final/qa_finalbins.txt \
  -o 2

## placing bins into tree with checkM
checkm tree --ali --nt -t $NTHREADS \
  -x fna IVB_pre/metabat2_1500/filtered_bins_part2/filtered_bins_final/ \
  IVB_pre/metabat2_1500/filtered_bins_part2/filtered_bins_final/tree/

## QA on the tree information with checkM
checkm tree_qa IVB_pre/metabat2_1500/filtered_bins_part2/filtered_bins_final/tree/ \
  -f IVB_pre/metabat2_1500/filtered_bins_part2/filtered_bins_final/tree_qa.txt \
  --tab_table

### that should be it, now you have refined bins!
echo "`date +"%Y%m%d %H:%M:%S"`: FINISHED_1500!"
```
