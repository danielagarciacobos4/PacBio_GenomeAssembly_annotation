# PacBio Genome Assembly and RNA annotation
This is an ongoing project to assemble a long-read PacBio genome and annotate it with RNA sequencing data. The genome belongs to the South American aquatic snake *Helicops angulatus.* This pipeline can largely be divided into two sections: 

1) PacBio Genome assembly
2) RNA seq annotation.

Within the first section involving the genome annotation, I will address scripts and results that a) check the quality of the raw reads obtained from PacBio, b) assemble the reference genome using the Hifiasm package, and c) evaluate the genome assembly quality using BUSCO and QUAST. 

For the second section which involves the genome annotation pipeline, I will address the scripts and results that a) perform genome masking with the earlgrey pipeline, b) Rna-seq corrections and assembly using Trinity for two tissues, and c) quality check for the assembled transcriptomes. 

Here is a nice picture of the species I chose to develop a high-quality annotated genome: 


![IMG_5229](https://github.com/user-attachments/assets/27b3c5ed-6c79-4a0b-8a62-a3f32be5e566) ![IMG_6603](https://github.com/user-attachments/assets/c7e1136d-0745-4735-9386-2c26631e741f)

