# Genome Assembly and RNA-seq Annotation

The following pipeline will show the scripts and results obtained from a PacBio (SMRT CELL) genome assembly and RNA seq annotation (still in process) for the South American aquatic snake *Helicops angulatus*. Depending on the resources and cluster capacity needed, some scripts are in the HUXLEY cluster (PBS scripting) and others are in the MENDEL cluster (SLURM- BATCH).

# 1. Genome assembly

I sequenced the genome of *Helicops angulatus* from Orinoquia, Colombia. This corresponds to a liver sample collected in 2022 and preserved in ethanol at 96% (IAvH-CT, Instituto Alexander von Humboldt in Colombia). DNA extractions were made with a kit for high molecular weight samples. The overall steps I did to check for the quality of the reads and genome assembly are as follows: 

1) kmer analysis of raw reads: Quality check of the raw reads
2) Genome assembly using Hifiasm
3) Busco analysis to check for completeness of the genome
4) Genome statistics with QUAST
5) *Merqury (I have not done this yet but it is another quality check of the assembly)*

## 1.1 kmer analysis using jellyfish for raw reads 
- [see full instructions](https://github.com/gmarcais/Jellyfish)
- Analysis performed in HUXLEY cluster 
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

Interpretation of the results: 

- Blue bars represent the distribution of k-mer counts (frequencies) across the sequencing dataset.
- Sequencing error seems to be low, as shown by the orange line.
- Approximate genome length of 1.3 Gbps.
- Low heterozygosity
- ~ 17X Coverage

Overall the Kmer analysis and visualization with genome scope show that we have high-quality sequencing data with minimal errors, meaning we have good data to do a de novo genome assembly. 

## 1.2 Genome assembly using Hifiasm

- See full documentation of [Hifiasm](https://hifiasm.readthedocs.io/en/latest/#)
- Analysis performed in the MENDEL cluster
- Raw files obtained from the sequencing facility Azenta are ~ 27 Gb
- Path to raw files: /home/dgarcia/mendel-nas1/PacBio/Helicops_angulatus_Aug2024/fastqc_raw/D6C18_Helicops_angulatus.hifireads.fastq.gz
- Hifiasm is a fast haplotype-resolved de novo assembler for PacBio HiFi reads.
- Hifiasm produces primary/alternate assemblies or partially phased assemblies only with HiFi reads.
- Since MENDEL does not have the Hifiasm installed, I created a conda environment called "assembly" that contained the Hifiasm package. 

### 1.2 Script genome assembly using Hifiasm

Important flags in the script: 
- -o : output name of the assembled genome
- -t : number of threads I want the cluster to use
- path to the raw reads the file
- running time: ~48 hours

```
#!/bin/sh
#SBATCH --job-name=hifiasm_NP_DGC_2
#SBATCH --nodes=1
#SBATCH --mem=170gb
#SBATCH --cpus-per-task=30
#SBATCH --time=40:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dgarcia@amnh.org
#SBATCH --output=hifiasm_assembly_NP4-%j-%x.out


#conda init

source ~/.bash_profile
conda activate assembly

hifiasm -o Helicops_angulatus_NP4.asm -t 30 /home/dgarcia/mendel-nas1/PacBio/Helicops_angulatus_Aug2024/fastqc_raw/D6C18_Helicops_angulatus.hifireads.fastq.gz

```

### 1.2 Results genome assembly using Hifiasm

These are the files obtained after the genome assembly: 

<img width="982" alt="Screenshot 2024-12-12 at 9 18 50 AM" src="https://github.com/user-attachments/assets/11fb4315-7d67-42b8-9f0b-ce936b5cda0c" />

Check out what each file means on the following [link](https://hifiasm.readthedocs.io/en/latest/interpreting-output.html). For the following steps we will use the Helicops_angulatus_NP4.asm.bp.p_ctg.gfa file which is the primary genome assembly. 

## 1.3 [Busco](https://busco.ezlab.org/) analysis of assembled genome

- BUSCO estimates the completeness and redundancy of processed genomic data based on universal single-copy orthologs. Read complete paper [here](https://doi.org/10.1093/molbev/msab199)
- Analysis performed in the MENDEL cluster
- path for BUSCO analysis: /home/dgarcia/mendel-nas1/PacBio/Helicops_angulatus_Aug2024/busco
- I had to install busco in a conda environment named BUSCO

### 1.3 Script for Busco analysis of assembled genome

Important flags in the script: 
- -m genome: specifies we are doing a busco analysis on a genome (busco can also be calculated for transcriptome or protein)
- -o BUSCO_H.angulatus: output name
- -i $Helicops_ang: specifies the input file, which is the variable Helicops_ang pointing to the genome assembly
- -l sauropsida_odb10: Specifies the lineage dataset to use. sauropsida_odb10 is the dataset for Sauropsida (reptiles and birds), indicating that the genome is being analyzed for completeness using a reference dataset relevant to reptiles
- -f: forces overwritting of existing output
- --metaeuk: Uses MetaEuk (MetaGeneAnnotator) as the gene predictor for BUSCO analysis
- --offline: Runs BUSCO in offline mode. The clusters in AMNH DO NOT allow to connect to the internet. For this reason it is importang to previously download the sauropsida_odb10 database.

```
#!/bin/sh
#SBATCH --job-name busco_helicops
#SBATCH --nodes=10
#SBATCH --mem=100gb
#SBATCH --time=144:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dgarcia@amnh.org

source ~/.bash_profile
conda activate BUSCO

Helicops_ang="/home/dgarcia/mendel-nas1/PacBio/Helicops_angulatus_Aug2024/assembly/assemblystats_NP4/Helicops_angulatus_NP4.asm.bp.p_ctg.fa"
busco -m genome -i $Helicops_ang -o BUSCO_H.angulatus -l sauropsida_odb10 -f --metaeuk --offline --download_path /home/dgarcia/mendel-nas1/PacBio/Helicops_angulatus_Aug2024/busco

```

### 1.3 Results for Busco analysis of assembled genome

- Good Busco scores!
- Completness score could be better (ideally we are looking for higher than 96%)
- Very low percentage of duplicates 1.6%. This is good because I dont have to clean up duplicates

<img width="484" alt="Screenshot 2024-12-12 at 10 52 03 AM" src="https://github.com/user-attachments/assets/b146de7e-80d4-4db9-8ff1-1d0d890c20bd" />

## 1.4 Genome statistics with QUAST

- QUAST stands for QUality ASsessment Tool. It evaluates genome/metagenome assemblies by computing various metrics.
- Analysis performed in the HUXLEY cluster
- see full manual [here](https://quast.sourceforge.net/docs/manual.html)
 
### 1.4 Script Genome statistics with QUAST

Important flags in the script:
- -o: output directory
- -r: reference genome path

```
python quast.py -o /home/dgarcia/nas4/phd/QUAST -r /home/dgarcia/nas4/phd/QUAST/Helicops_angulatus_NP4.asm.bp.p_ctg.fa
```

### 1.5 Results Genome statistics with QUAST

- total genome size: ~ 2.1 GB
- N50: 44.6 Mbp
- L50: 9

<img width="429" alt="image" src="https://github.com/user-attachments/assets/2b91e51b-6d7a-42e3-8697-a16fcb8d5f90" />


# 2. Annotation: 

For now, this section will include the steps two main steps that are needed and will build up to the genome annotation with RNA seq: 

1) Soft mask the genome with Earl Grey
2) Assemble RNA seq for two tissues (liver and kidney)
3) Busco scores for assembled RNA

## 2.1 Earl Grey: soft masking the genome

- Given an input genome, Earl Grey will run through numerous steps to identify, curate, and annotate transposable elements (TEs).
- see full documentation of Earl Grey [here](https://github.com/TobyBaril/EarlGrey)
- This analysis was performed in the MENDEL cluster of AMNH
- path to the folder containing scripts and results: /home/dgarcia/mendel-nas1/PacBio/Helicops_angulatus_Aug2024/earl_grey

### 2.1 Script Earl Grey: soft masking the genome

Important flags of the script: 
- -g: Specifies the input genome file: /home/dgarcia/mendel-nas1/PacBio/Helicops_angulatus_Aug2024/assembly_final_18Oct2024/Helicops_angulatus_NP4.asm.bp.p_ctg.fa
- -s: Specifies the session name or identifier for this analysis
- -t: Specifies the number of threads (parallel processes) to use.
- -0: Specifies the output directory

```
#!/bin/sh

#SBATCH --job-name earlgrey2
#SBATCH --nodes=1
#SBATCH --tasks-per-node=48
#SBATCH --mem=150gb
#SBATCH --time=300:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dgarcia@amnh.org

#activate earlgrey with mamba before submitting

source ~/.bash_profile
export PATH=/home/dgarcia/mendel-nas1/miniforge3/bin:$PATH
mamba init
mamba activate earlgrey2

genome="/home/dgarcia/mendel-nas1/PacBio/Helicops_angulatus_Aug2024/assembly_final_18Oct2024/Helicops_angulatus_NP4.asm.bp.p_ctg.fa"

earlGrey -g $genome -s Helicops_angulatus_18Oct2024 -t $SLURM_NTASKS_PER_NODE -d yes -o /home/dgarcia/mendel-nas1/PacBio/Helicops_angulatus_Aug2024/earl_grey

##This in modeling my repeats, masking them, and then analysis compared to known paired database
```

### 2.1 Results Earl Grey: soft masking the genome

- the result of this analysis is a soft masked genome.
- It also gives this pie chart to have an idea of the repeat elements detected and identified in the genome:
- In this case, we ca see that this is a highly repetitive genome.

<img width="1000" alt="Screenshot 2024-12-12 at 12 38 52 PM" src="https://github.com/user-attachments/assets/e29bb427-c4a7-4736-bc82-70d4a4649c05" />


## 2.2 Assemble RNA seq with Trinity (Need to update this with new reads sequenced on june 2025)

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


## 2.2.1 Script: Steps 1-5 (fastqc, Rcorrector and trimmomatic)

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

 
## 2.2.2 Script: Step 6 [Trinity](https://github.com/trinityrnaseq/trinityrnaseq)

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
## 2.2.2 Results: Step 6 (Trinity)
- approximate time to obtain results of 2 tissues ~ 1.5 GB (bad coverage)
- time to run analysis for 2 tissues: 5 hours.
- obtained two files per tissue as results: 1) name_Trinity.fasta and 2) name_Trinity.fasta.gene_trans_map
- each Trinity.fasta file weights ~ 100m (kidney) and ~56m (liver)

## 2.2.3 Script: Step 7 [BUSCO](https://busco.ezlab.org/)
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
## 2.2.3 Results: Step 7 [BUSCO](https://busco.ezlab.org/)

- Results for kidney:
<img width="561" alt="Screenshot 2024-12-01 at 4 16 58 PM" src="https://github.com/user-attachments/assets/ea330f01-b73a-4c88-8d7c-8998ad2b2b91">

- Results for liver:
<img width="576" alt="Screenshot 2024-12-01 at 4 17 59 PM" src="https://github.com/user-attachments/assets/2711a62f-db5d-41b5-81a8-c59698046e58">



# Funannotate

### Installation: 

1) Create a conda environment named funannotate. I installed with mamba because it is faster
```
mamba create -n funannotate -c conda-forge -c bioconda funannotate

```
Keep in mind: 

- <img width="1074" height="690" alt="Screenshot 2025-08-31 at 3 30 15 PM" src="https://github.com/user-attachments/assets/0a631334-81ea-42eb-81a3-57ed656db891" />
- Path of folder: /home/dgarcia/mendel-nas1/PacBio/Helicops_angulatus_Aug2024/funannotate/
- Structure of folder: 

<img width="730" height="591" alt="Screenshot 2025-08-31 at 4 00 56 PM" src="https://github.com/user-attachments/assets/39e001b3-7298-4ea4-890e-1d97a2b6e000" />

2) Download the Funannotate databases:
```
mkdir -p $CONDA_PREFIX/etc/conda/activate.d $CONDA_PREFIX/etc/conda/deactivate.d
echo 'export FUNANNOTATE_DB=/home/dgarcia/mendel-nas1/PacBio/Helicops_angulatus_Aug2024/funannotate/02_db/funannotate_db' > $CONDA_PREFIX/etc/conda/activate.d/funannotate.sh
echo 'unset FUNANNOTATE_DB' > $CONDA_PREFIX/etc/conda/deactivate.d/funannotate.sh
```
3) Make sure we downloaded everything okay, and run a test to make sure the installation worked okay

```
funannotate database --show

funannotate test -t all --cpus 8
```

If the test is successful it shows something like this: 
<img width="1596" height="857" alt="Screenshot 2025-08-31 at 5 03 31 PM" src="https://github.com/user-attachments/assets/14703f38-6b80-4c41-818d-01ebaf8eb6f6" />
<img width="839" height="215" alt="Screenshot 2025-08-31 at 5 03 59 PM" src="https://github.com/user-attachments/assets/9d6edd31-dd56-4aba-a3cd-4944350390e5" />




  
# References: 

- Baril, T., Galbraith, J.G., and Hayward, A., Earl Grey: A Fully Automated User-Friendly Transposable Element Annotation and Analysis Pipeline, Molecular Biology and Evolution, Volume 41, Issue 4, April 2024, msae068 doi:10.1093/molbev/msae068
- Haoyu Cheng, Gregory T. Concepcion, Xiaowen Feng, Haowen Zhang & Heng Li. Haplotype-resolved de novo assembly using phased assembly graphs with hifiasm. Nature Methods. (2021).
- Gurevich, A., Saveliev, V., Vyahhi, N., & Tesler, G. (2013). QUAST: quality assessment tool for genome assemblies. Bioinformatics, 29(8), 1072-1075.
- Grabherr MG, Haas BJ, Yassour M, Levin JZ, Thompson DA, Amit I, Adiconis X, Fan L, Raychowdhury R, Zeng Q, Chen Z, Mauceli E, Hacohen N, Gnirke A, Rhind N, di Palma F, Birren BW, Nusbaum C, Lindblad-Toh K, Friedman N, Regev A. Full-length transcriptome assembly from RNA-seq data without a reference genome. Nat Biotechnol. 2011 May 15;29(7):644-52. doi: 10.1038/nbt.1883. PubMed PMID: 21572440.
- Mosè Manni, Matthew R Berkeley, Mathieu Seppey, Felipe A Simão, Evgeny M Zdobnov, BUSCO Update: Novel and Streamlined Workflows along with Broader and Deeper Phylogenetic Coverage for Scoring of Eukaryotic, Prokaryotic, and Viral Genomes, Molecular Biology and Evolution, Volume 38, Issue 10, Pages 4647–4654, https://doi.org/10.1093/molbev/msab199. (2021)













