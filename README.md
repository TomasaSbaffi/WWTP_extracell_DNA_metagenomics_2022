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

* Entrez Direct v16.2: https://www.ncbi.nlm.nih.gov/books/NBK179288/
* SRA Toolkit v2.11.3: https://github.com/ncbi/sra-tools
* fasterq-dump v2.10.8: https://github.com/glarue/fasterq_dump
* GNU parallel: https://www.gnu.org/software/parallel
* anvi’o v7.0: https://merenlab.org/software/anvio
* BLAST v2.10.1: https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+
* seqtk v1.3: https://github.com/lh3/seqtk
* MUSCLE v3.8.1551: http://www.drive5.com/muscle
* IQ-TREE v2.1.4: http://www.iqtree.org
* CoverM v0.6.1: https://github.com/wwood/CoverM/

We have used an Atos BullSequana X400 system running the Red Hat Enterprise Linux Server 7.7 (Maipo).  
You should be able to run the analysis with any UNIX-based OS.


