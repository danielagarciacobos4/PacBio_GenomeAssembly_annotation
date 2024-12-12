# Genome Assembly and RNA-seq Annotation

The following pipeline will show the scripts and results obtained from a PacBio (SMRT CELL) genome assembly and RNA seq annotation (still in process) for the South American aquatic snake *Helicops angulatus*. Depending on the resources and cluster capacity needed, some scripts are in the HUXLEY cluster (PBS scripting) and others are in the MENDEL cluster (SLURM- BATCH). Here is a nice picture of the species I chose to develop a high-quality annotated genome:

![IMG_5229](https://github.com/user-attachments/assets/f335c39f-083c-406d-aab0-746dcfcf26fe)


# 1. Genome assembly

I sequenced the genome of *Helicops angulatus* from Orinoquia, Colombia. This corresponds to a liver sample collected in 2022 and preserved in ethanol at 96% (IAvH-CT, Instituto Alexander von Humboldt in Colombia). DNA extractions were made with a kit for high molecular weight samples. The overall steps I did to check for the quality of the reads and genome assembly are as follows: 

1) kmer analysis of raw reads: Quality check of the raw reads
2) Hifiasm to assemble the genome (purging?)
3) Busco analysis to check for completeness of the genome. Quality check of the assembly
4) *Merqury (I have not done this yet). Quality check of the assembly*

## 1.1 kmer analysis using jellyfish for raw reads 
- [see full instructions](https://github.com/gmarcais/Jellyfish) 
- path to the analysis: /home/dgarcia/nas5/PacBio
- D6C18_Helicops_angulatus.hifireads.fastq.gz corresponds to the file with the raw reads.
- Jellyfish is a tool for fast, memory-efficient counting of k-mers in DNA. A k-mer is a substring of length k, and counting the occurrences of all such substrings is a central step in many analyses of DNA sequences.
- It outputs its k-mer counts in a binary format, which can be translated into a human-readable text format using the "jellyfish dump" command, or queried for specific k-mers with "jellyfish query".
- The command also includes a "jellyfish histo" which creates an histogram of k-mer occurrences that can be read in [GenomeScope](http://genomescope.org/)

### 1.1 Script kmer analysis using jellyfish

```
#!/bin/bash
#PBS -N kmer_dgc
#PBS -q batch
#PBS -S /bin/bash
#PBS -l nodes=1:ppn=1
#PBC -l mem=30GB
#PBS -l ncpus=40
#PBS -l walltime=50:00:00

module load jellyfish-2.3.0
cd /home/dgarcia/nas5/

jellyfish count -C -m 21 -s 1000000000 -t 38 -o D6C18_Helicops_angulatus.hifireads.fastq.jf <(zcat D6C18_Helicops_angulatus.hifireads.fastq.gz)
jellyfish histo D6C18_Helicops_angulatus.hifireads.fastq.jf -t 38 > D6C18_Helicops_angulatus.hifireads.fastq.histo
```
### 1.1 Results kmer analysis using jellyfish


![Screenshot 2024-11-27 at 4 09 35 PM](https://github.com/user-attachments/assets/33bfc30c-dfe9-4f3b-ab17-f41de558ec50)


# 2. Annotation: 

Earl Grey: soft masking of the genome. This is a needed step for genome annotation (actually I believe this is the first step of genome annotation, but I have to check and confirm this statement). 

I have obtained RNA seq data from Illumina of two tissue samples: liver and kidney. Unfortunately, these samples are not from the same individual I used for the PacBio genome. However, it is an individual from the same species (*Helicops angulatus*) and same population (Orinoquian region - Meta, Colombia).

Samples are: 
- DGC-R-7-kidney_R1_001.fastq.gz,
- DGC-R-7-kidney_R2_001.fastq.gz,
- DGC-R-7-liver_R1_001.fastq.gz,
- DGC-R-7-liver_R2_001.fastq.gz  
- They each has a size of **~ 1.5GB**

I am following the steps mentioned in: https://github.com/harvardinformatics/TranscriptomeAssemblyTools/tree/master and https://www.notion.so/RNA-seq-data-processing-cc1fcd681bf040bcb181ee9f68744d9a which are the following:

1. run **fastqc** on raw fastq reads to identify potential issues with Illumina sequencing libraries such as retained adapter, over-represented (often low complexity) sequences, greater than expected decrease in base quality with read cycle.
2. perform kmer-based read corrections with with [rCorrector](https://github.com/mourisl/Rcorrector/tree/master), see [Song and Florea 2015, Gigascience](https://gigascience.biomedcentral.com/articles/10.1186/s13742-015-0089-y)
3. remove read pairs where at least one read has been flagged by rCorrector as containing an erroneous kmer, and where it was not possible to correct the errors computationally
4. remove read pairs where at least read contains an over-represented sequence
5. perform high quality trimming with trimmomatic
6. assemble reads with [trinity](https://github.com/trinityrnaseq/trinityrnaseq)

I am running one script with PBS_ARRAY_INDEX that allows to run a large collection of PBS runs to be executed in one script. This is specified in the PBS -J 1-2 and with the ${PBS_ARRAY_INDEX}. This first script runs from step 1-5 in the above description (meaning trnity is not ran yet)


## 2.1 Script: Steps 1-5 (fastqc, Rcorrector and trimmomatic)

I run these analyses in the Huxley cluster of AMNH using a PBS script. Personal reminder: if needed this will need to be transferred to AMNH Mendel cluster (SLURM) which is where I have my PacBio genome assembly. 

```
#!/bin/bash
#PBS -q batch 
#PBS -l mem=20gb
#PBS -m abe
#PBS -M dgarcia@amnh.org
#PBS -l ncpus=8
#PBS -q batch
#PBS -N rna_trinity
#PBS -l walltime=100:00:00
#PBS -J 1-2

#READ PROCESSING
source ~/.bash_profile
export PATH=/home/dgarcia/nas4/bin/miniconda3/bin:$PATH
export OMP_NUM_THREADS=8 #we are using 8 threads, 8 cpus
conda activate prepRNA #upload trimmomatic and rcorrector since the rest of the packages are available in the Huxley cluster of AMNH.
module load fastqc-0.11.9 perl-5.26.0 jre-1.8.0_251 Trinity-2.12.0 jellyfish-2.3.0 bowtie2-2.3.5.1

wd=/home/dgarcia/nas5/rna
cd ${wd}
# The following file text has the names of the folder where I have the R1 and R2 fastq files of raw reads
directories_file=directories.txt
#PBS_ARRAY_INDEX allows a large collection of PBS runs to be executed in one script. This is specified in the PBS -J 1-2
dir=$(sed -n "${PBS_ARRAY_INDEX}p" "$directories_file")

echo ${PBS_ARRAY_INDEX}
echo $dir
R1=/home/dgarcia/nas5/rna/raw_reads/${dir}/*R1*
R2=/home/dgarcia/nas5/rna/raw_reads/${dir}/*R2*

echo $OMP_NUM_THREADS
##1. fastqc; make a fastqc folder in my wd before running this command

fastqc -o fastqc $R1 $R2 -t ${OMP_NUM_THREADS}

#2. kmer read corrections; make a folder name rcorrected before running this command. This perl script is purging and if possible fixing bad quality reads.

perl Rcorrector/run_rcorrector.pl -1 $R1 -2 $R2 -od /home/dgarcia/nas5/rna/raw_reads/${dir}/rcorrected -t ${OMP_NUM_THREADS}

#3 correct other reads, python file from the github mentioned above

if ls /home/dgarcia/nas5/rna/raw_reads/${dir}/rcorrected/*R1* 1> /dev/null 2>&1; then
	mkdir ${wd}/raw_reads/${dir}/filtered
	cd ${wd}/raw_reads/${dir}/filtered
	python /home/dgarcia/nas5/rna/FilterPEfastq.py -1 ${wd}/raw_reads/${dir}/rcorrected/*R1*.fq.gz -2 ${wd}/raw_reads/${dir}/rcorrected/*R2*.fq.gz
# -s ${wd}/transfer/${dir}/filtered
	mv ${wd}/raw_reads/${dir}/rcorrected/unfixrm_* ${wd}/raw_reads/${dir}/filtered/
fi

#4. Trim reads

if  ls /home/dgarcia/nas5/rna/raw_reads/${dir}/filtered/*_for_paired.fq.gz 1> /dev/null 2>&1; then
	echo "ran trimming"
elif ls /home/dgarcia/nas5/rna/raw_reads/${dir}/filtered/*R1* 1> /dev/null 2>&1; then
	cd ${wd}/raw_reads/${dir}/filtered
	name=$(echo *R1*| awk -F'unfixrm_|_R1' '{print $2}') # get a bash variable that is the name of the sample, taken from the directory, and retaining text prior to _
	echo $name
	trimmomatic  PE -threads ${OMP_NUM_THREADS} *R1*.fq.gz *R2*.fq.gz ${name}_for_paired.fq.gz ${name}_for_unpaired.fq.gz ${name}_rev_paired.fq.gz ${name}_rev_unpaired.fq.gz ILLUMINACLIP:/home/ddebaun/nas5/pseudo-it/TruSeq3-PE-2.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:36
fi
```

## Results: Steps 1-5 (fastqc, Rcorrector and trimmomatic)
- approximate time to obtain results of 2 tissues ~ 1.5 GB (bad coverage) = 3.5 hours

**1) Fast qc results:**
- Sequence Duplication Levels seem to be high. However, I think this makes sense for RNA sequence data since there are going to be some more overexpressed genes and thus have duplicated sequences compared to others.
<img width="899" alt="Screenshot 2024-11-26 at 10 48 07 AM" src="https://github.com/user-attachments/assets/8b23d1ab-b353-4b61-9942-c577269d0a7e">

-  Adaptor content: ~20% of the Illumina Universal Adapter.
  
**2) kmer read corrections:**
- Table summarizing the reads purged and corrected after running rcorrector:
- Around 15% of the reads were purged after correcting or eliminating bad quality reads (as a means of comparison, Dylan's RNA seq data purged ~7% of the reads).
  
![image](https://github.com/user-attachments/assets/b0069a03-c538-4072-9e60-545ab6bc062f)

**3) correct other reads, python file from the github mentioned above:**
- This step further corrects RNA reads. The results of this correction is a two file starting with the word unfixrm (e.g. unfixrm_kidney_Helicops_angulatus_IAvH-CT-36861_R2_001.cor.fq.gz). One for R1 and R2.
  
**4) Trim reads**
- This step trims the adaptors from library prep. After trimming adaptors I am left with paired reads of ~ 1 GB. This is very low-coverage reads, but this is the outcome of the initial sequencing resources requested (remember that the initial file is around ~1.5 GB - before running rcorrector).
- After running this step I am left with reverse paired, reverse unpaired, foward paired and foward unpaired file.

 
## 2.2 Script: Step 6 [Trinity](https://github.com/trinityrnaseq/trinityrnaseq)

- Trinity assembles transcript sequences from Illumina RNA-Seq data. Take a look at the [wiki](https://github.com/trinityrnaseq/trinityrnaseq/wiki) page for full step by step of how trinity works.
- path of the script: /home/dgarcia/nas5/rna
- name of the script: trinity.sh

```
#!/bin/bash
#PBS -q batch
#PBS -l mem=100gb
#PBS -m abe
#PBS -M dgarcia@amnh.org
#PBS -l ncpus=10
#PBS -q batch
#PBS -N trinityRNA
#PBS -l walltime=100:00:00
#PBS -J 1-2

cd /home/dgarcia/nas5/rna
module load perl-5.26.0 jre-1.8.0_251 Trinity-2.12.0 jellyfish-2.3.0 bowtie2-2.3.5.1

directories_file=directories.txt
dir=$(sed -n "${PBS_ARRAY_INDEX}p" "$directories_file")
echo ${PBS_ARRAY_INDEX}
echo $dir
R1=/home/dgarcia/nas5/rna/raw_reads/${dir}/filtered/*for_paired*
R2=/home/dgarcia/nas5/rna/raw_reads/${dir}/filtered/*rev_paired*
mkdir trimmed_trinity_${dir}
cd trimmed_trinity_${dir}
#make sure max mem matches mem allocation up top
Trinity --full_cleanup --seqType fq --left $R1 --right $R2 --CPU ${OMP_NUM_THREADS} --max_memory 100G --output trimmed_trinity_${dir}
```
## Results: Step 6 (Trinity)
- approximate time to obtain results of 2 tissues ~ 1.5 GB (bad coverage)
- time to run analysis for 2 tissues: 5 hours.
- obtained two files per tissue as results: 1) name_Trinity.fasta and 2) name_Trinity.fasta.gene_trans_map
- each Trinity.fasta file weights ~ 100m (kidney) and ~56m (liver)

## 2.3 Script: Step 7 [BUSCO](https://busco.ezlab.org/)
- BUSCO estimates the completeness and redundancy of processed genomic data based on universal single-copy orthologs. Read complete paper [here](https://doi.org/10.1093/molbev/msab199)
- script path and name: /nas5/dgarcia/rna/busco.sh
- Run time: started at 12:51 pm (dec 1, 2023)..

```
#!/bin/bash
#PBS -q batch 
#PBS -l mem=30gb
#PBS -m abe
#PBS -M dgarcia@amnh.org
#PBS -l ncpus=5
#PBS -q batch
#PBS -N busco_scores
#PBS -l walltime=10:00:00
#PBS -J 1-2

module load busco-5.5.0
export AUGUSTUS_CONFIG_PATH=/nas5/dgarcia/augustus-3.4.0/config
cd /nas5/dgarcia/

directories_file=rna/directories.txt
dir=$(sed -n "${PBS_ARRAY_INDEX}p" "$directories_file")
echo ${PBS_ARRAY_INDEX}
echo $dir

busco -c $OMP_NUM_THREADS -f -i /nas5/dgarcia/rna/trimmed_trinity_${dir}/trimmed_trinity_${dir}.Trinity.fasta -o /nas5/dgarcia/rna_${dir}.busco -l sauropsida_odb10 -m transcriptome --offline --download_path /home/dgarcia/nas5/busco_downloads
```
## 2.3 Results: Step 7 [BUSCO](https://busco.ezlab.org/)

- Results for kidney:
<img width="561" alt="Screenshot 2024-12-01 at 4 16 58 PM" src="https://github.com/user-attachments/assets/ea330f01-b73a-4c88-8d7c-8998ad2b2b91">

- Results for liver:
<img width="576" alt="Screenshot 2024-12-01 at 4 17 59 PM" src="https://github.com/user-attachments/assets/2711a62f-db5d-41b5-81a8-c59698046e58">



  














