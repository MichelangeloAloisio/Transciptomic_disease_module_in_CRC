Dataset Completeness Pipeline
This pipeline takes part of the scientifca paper: 'Uncovering Transcriptomic Modules in Colorectal Cancer: A Causal Inference Framework for Network-Based Analysis' 

is designed to generate bar plots for four gene categories: protein-coding, ribosomal RNA, mitochondrial RNA, and other non-coding RNA, using raw count datasets produced from BAM files (e.g., using HTSeq count or Subread featureCounts).

Pipeline Overview
The primary aim of this pipeline is to assess the completeness of raw count datasets, typically downloaded from public repositories like GEO. Since raw counts are essential for differential expression analysis using tools like DESeq2 and limma, this pipeline is useful for performing quality control of raw datasets before proceeding with differential analysis.

Key Features:
Gene Category Classification: The pipeline calculates the percentage of reads assigned to four key gene categories typically found in datasets generated using a Total RNA protocol.
Bar Plot Visualization: After classification, the percentage of counts in each category is calculated and visualized in a bar plot, making it easier to assess dataset composition.
Aim
The pipeline calculates the percentage of reads assigned to four gene categories that should be present in a dataset generated using a Total RNA protocol. For each gene in the raw RNA-Seq count tables:

The pipeline checks its classification into one of the four categories.
The percentage of counts in each category is calculated and visualized in a bar plot.
This allows researchers to assess whether the raw count data from their RNA-Seq experiments is appropriately distributed across these categories.

Pre-requisites
Before running the pipeline, you need to download the necessary input files and install the required dependencies:

GTF File
Gene annotations (specifically, the gene type) are extracted from the GTF file downloaded from ENCODE.
Download the GTF file from:
https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_35/gencode.v35.basic.annotation.gtf.gz

The Korean Dataset
The dataset (56609x342) can be downloaded from:
https://doi.org/10.5281/zenodo.8333650

The Swedish Dataset
The dataset (19765x1183) is available in the ArrayExpress database under accession number E-MTAB-12862.

Running the Pipeline
Follow the steps below to run the pipeline:

Step 1: Clone the Repository
Clone the GitHub repository to your local machine:


git clone https://github.com/MichelangeloAloisio/Transcriptomic_disease_module_in_CRC.git
cd Transcriptomic_disease_module_in_CRC/Dataset_completness/00_input_files_gtf_raw_counts
Step 2: Download Input Files
Download the following files into the 00_input_files_gtf_raw_counts folder:

The GTF file from the Gencode v35 GTF file.
The raw count datasets (Korean and Swedish datasets) from their respective links.
Step 3: Run the Pipeline
Once the input files are in place, navigate back to the main directory and execute the pipeline:


cd ../
python3 dataset_completness.py
This will process the raw count data, classify the genes into categories, and generate the corresponding bar plot.

Gene Categories
The pipeline classifies genes into the following categories based on their gene type:

1. Protein-Coding Genes:
protein_coding
TR_C_gene
TR_J_gene
TR_V_gene
TR_D_gene
IG_C_gene
IG_J_gene
IG_V_gene
IG_D_gene
2. Pseudogenes (Non-Coding):
transcribed_unitary_pseudogene
TR_J_pseudogene
translated_unprocessed_pseudogene
unprocessed_pseudogene
processed_pseudogene
IG_J_pseudogene
IG_C_pseudogene
IG_V_pseudogene
translated_processed_pseudogene
polymorphic_pseudogene
unitary_pseudogene
pseudogene
IG_pseudogene
TR_V_pseudogene
3. Ribosomal RNA (Non-Coding):
rRNA
rRNA_pseudogene
4. Mitochondrial RNA (Non-Coding):
Mt_rRNA
Mt_tRNA
5. Other Non-Coding RNA:
lncRNA
TEC
scaRNA
misc_RNA
scRNA
sRNA
transcribed_unprocessed_pseudogene
snoRNA
snRNA
miRNA
vault_RNA
ribozyme
