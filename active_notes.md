# Genome Assembly with PacBio

Need to update everything I have done so far: 
- kmer analysis of raw reads: Quality check of the raw reads
- Hifiasm to assemble the genome (purging?)
- Busco analysis to check for completeness of the genome. Quality check of the assembly
- *Merqury (I have not done this yet). Quality check of the assembly*

# Annotation
- Earl Grey: soft masking of the genome. This is a needed step for genome annotation (actually I believe this is the first step of genome annotation, but I have to check and confirm this statement). 

## Rna Seq data 

I have obtained RNA seq data from Illumina of two tissue samples: liver and kidney. Unfortunately, these samples are not from the same individual I used for the PacBio genome. However, it is an individual from the same species (*Helicops angulatus*) and same population (Orinoquian region - Meta, Colombia).

Samples are: DGC-R-7-kidney_R1_001.fastq.gz, DGC-R-7-kidney_R2_001.fastq.gz, DGC-R-7-liver_R1_001.fastq.gz, DGC-R-7-liver_R2_001.fastq.gz and they each have a size of ~ 1.5GB

I am following the steps mentioned in: https://github.com/harvardinformatics/TranscriptomeAssemblyTools/tree/master and https://www.notion.so/RNA-seq-data-processing-cc1fcd681bf040bcb181ee9f68744d9a which are the following:

1. run **fastqc** on raw fastq reads to identify potential issues with Illumina sequencing libraries such as retained adapter, over-represented (often low complexity) sequences, greater than expected decrease in base quality with read cycle.
2. 

### 1. Fastqc

I run these analyses in the Huxley cluster of AMNH using a PBS script. Personal reminder: if needed this will need to be transferred to AMNH Mendel cluster (SLURM) which is where I have my PacBio genome assembly. 

```
#!/bin/bash

#PBS -N fastqc_rna
#PBS -q batch
#PBS -m abe
#PBS -M dgarcia@amnh.org
#PBS -e /home/dgarcia/nas5/rna/fastqc
#PBS -o /home/dgarcia/nas5/rna/fastqc
#PBS -l nodes=1:ppn=28
#PBS -l walltime=40:00:00

module load fastqc-0.11.9

fastqc /home/dgarcia/nas5/rna/DGC-R-7-liver_R1_001.fastq.gz \
/home/dgarcia/nas5/rna/DGC-R-7-liver_R2_001.fastq.gz \
-o /home/dgarcia/nas5/rna/fastqc
```
And now the script for the kidney samples. For some reason, I was not able to run the fastqc of the 2 samples (liver and kidney) in the same script. This is some troubleshooting I must due to not have to run the same script mutiple times. 

```
#PBS -N fastqc_rna
#PBS -q batch
#PBS -m abe
#PBS -M dgarcia@amnh.org
#PBS -e /home/dgarcia/nas5/rna/fastqc
#PBS -o /home/dgarcia/nas5/rna/fastqc
#PBS -l nodes=1:ppn=28
#PBS -l walltime=40:00:00

module load fastqc-0.11.9

fastqc /home/dgarcia/nas5/rna/DGC-R-7-kidney_R1_001.fastq.gz \
/home/dgarcia/nas5/rna/DGC-R-7-kidney_R2_001.fastq.gz \
-o /home/dgarcia/nas5/rna/fastqc

```
