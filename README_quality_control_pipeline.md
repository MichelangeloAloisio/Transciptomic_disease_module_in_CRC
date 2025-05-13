RNA-Seq Quality Control Pipeline


############################################################################################################################################################################
############################################              PIPELINE OVERVIEW

This RNA-Seq Quality Control pipeline processes raw RNA-Seq count datasets to extract essential information and assess the quality of the RNA-Seq count matrix.

The pipeline performs the following tasks:

1) Sample Coverage Calculation (ONLY THIS PART OF PIPELINE IS DESCRIBED IN THE PAPER):

Computes the coverage for each sample (saved in 01_COVERAGE_PER_SAMPLE_TABLE.csv), allowing for the manual inspection of low-coverage samples.

2) Distribution of Counts:

Visualizes the distribution of counts to determine whether the dataset contains RAW counts or NORMALIZED counts.


3) Gene Category Distribution:

Calculates and visualizes the percentage of counts per gene category. This step helps to assess the completeness of the dataset. The categories are defined in gene_type tag present in GTF file (GTF file example: gencode.v35.basic.annotation.gtf.gz) 
This step is foundamental for RNA-Seq libraries prepared with ribosomal depletion. Indeed in this case, we expect high counts in non-coding genes and a lack of such counts indicates incomplete data, which is critical for the normalization step useful for differential expressiona analysis. On the contrary, for poly-A enriched libraries, non-coding genes should have negligible counts.

A total of 5 categories are evaluated by the script: 
Category 1) Protein-Coding Genes: Includes canonical protein-coding genes and immune-related genes. The list of gene_type attribute in the GTF file for this category is: 'protein_coding', 'TR_C_gene', 'TR_J_gene', 'TR_V_gene', 'TR_D_gene', 'IG_C_gene', 'IG_J_gene', 'IG_V_gene', 'IG_D_gene'


Category 2) Pseudogenes (Non-Coding), that are Genes that do not produce functional proteins.The list of gene_type attribute in the GTF file for this category is: 'transcribed_unitary_pseudogene', 'TR_J_pseudogene', 'translated_unprocessed_pseudogene', 'unprocessed_pseudogene', 'processed_pseudogene', 'IG_J_pseudogene', 'IG_C_pseudogene', 'IG_V_pseudogene', 'translated_processed_pseudogene', 'polymorphic_pseudogene', 'unitary_pseudogene', 'pseudogene', 'IG_pseudogene', 'TR_V_pseudogene'



Category 3) Ribosomal RNA (Non-Coding). The list of gene_type attribute in the GTF file for this category is: 'rRNA', 'rRNA_pseudogene'



Category 4) Mitochondrial RNA (Non-Coding).The list of gene_type attribute in the GTF file for this category is: 'Mt_rRNA', 'Mt_tRNA'


Category 5) Other Non-Coding RNA. The list of gene_type attribute in the GTF file for this category is: 'lncRNA', 'TEC', 'scaRNA', 'misc_RNA', 'scRNA', 'sRNA', 'transcribed_unprocessed_pseudogene', 'snoRNA', 'snRNA', 'miRNA', 'vault_RNA', 'ribozyme'


###############################################################################################################################################################################################
#######################################                          INPUT FILES

INPUT 1) GTF File: Contains gene annotations (gene_id, gene_name, gene_type).
Example: gencode.v35.basic.annotation.gtf.gz.
Download from Gencode.


INPUT 2) Raw Count Matrix: A tab-delimited file with gene expression counts across samples. The first column should contain gene identifiers, followed by columns for each sample.
Example: gene_counts.txt.

!!! IMPORTANT !!!
Ensure the GTF file and raw counts matrix file are in the 00_input_files_gtf_raw_counts directory.

