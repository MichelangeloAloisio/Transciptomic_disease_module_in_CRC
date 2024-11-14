RNA-Seq Quality Control Pipeline
Pipeline Overview
This RNA-Seq Quality Control pipeline processes raw RNA-Seq count datasets to extract essential information and assess the quality of the RNA-Seq count matrix.

The pipeline performs the following tasks:

Sample Coverage Calculation:

Computes the coverage for each sample (saved in 01_COVERAGE_PER_SAMPLE_TABLE.csv), allowing for the manual inspection of low-coverage samples.
Distribution of Counts:

Visualizes the distribution of counts to determine whether the dataset contains RAW counts or NORMALIZED counts.
Gene Category Distribution:

Calculates and visualizes the percentage of counts per gene category. This helps assess the completeness of the dataset.
Example: For RNA-Seq libraries prepared with ribosomal depletion, we expect high counts in non-coding genes. A lack of such counts indicates incomplete data, which is critical for the normalization step.
Example: For poly-A enriched libraries, non-coding genes should have negligible counts.
rRNA Percentage:

Calculates the percentage of rRNA to assess gene depletion in total RNA-Seq libraries (e.g., ribosomal RNA depletion).
Gene Category Summary:

Calculates the percentage of genes in five gene categories, grouping gene types by their biological function (category).
Unlike point 3, this step calculates the percentage of genes in each category based on the symbol column, not sequencing counts. Combined with point 3, this data helps understand whether rows have been removed from the table, which could affect eventual normalization.
Gene Categories
Protein-Coding Genes: Includes canonical protein-coding genes and immune-related genes. The list of gene_type in this section is:

'protein_coding', 'TR_C_gene', 'TR_J_gene', 'TR_V_gene', 'TR_D_gene', 'IG_C_gene', 'IG_J_gene', 'IG_V_gene', 'IG_D_gene'
Pseudogenes (Non-Coding): Genes that do not produce functional proteins:

'transcribed_unitary_pseudogene', 'TR_J_pseudogene', 'translated_unprocessed_pseudogene', 'unprocessed_pseudogene', 'processed_pseudogene', 'IG_J_pseudogene', 'IG_C_pseudogene', 'IG_V_pseudogene', 'translated_processed_pseudogene', 'polymorphic_pseudogene', 'unitary_pseudogene', 'pseudogene', 'IG_pseudogene', 'TR_V_pseudogene'
Ribosomal RNA (Non-Coding): Genes for ribosomal RNA (rRNA):

'rRNA', 'rRNA_pseudogene'
Mitochondrial RNA (Non-Coding): Genes encoding mitochondrial RNA:

'Mt_rRNA', 'Mt_tRNA'
Other Non-Coding RNA: Various non-coding RNAs with diverse functions:

'lncRNA', 'TEC', 'scaRNA', 'misc_RNA', 'scRNA', 'sRNA', 'transcribed_unprocessed_pseudogene', 'snoRNA', 'snRNA', 'miRNA', 'vault_RNA', 'ribozyme'
Pipeline Workflow
Extract Gene Annotations:

Processes the GTF file to extract gene annotations, including the gene_type attribute. (GTF file example: gencode.v35.basic.annotation.gtf.gz)
Process Raw Count Matrix:

Matches each gene in the count matrix to its corresponding gene type. Duplicates are checked and flagged.
Compute Sample Coverage:

Summarizes total coverage for each sample to assess sequencing depth.
Visualize Count Distribution:

Generates a boxplot of log-transformed counts across samples.
Gene-Level Summarization:

Summarizes gene expression counts across samples.
Categorize Genes:

Classifies each gene into predefined categories based on its gene_type and stores binary indicators for each category.
Merge Gene Data with Coverage:

Merges the gene categorization data with sample coverage and checks for missing values.
Generate Category Distribution Bar Plot:

Visualizes the percentage of counts in each gene category.
Calculate and Display Percentages:

Computes the percentage of genes in each category and displays it in a final bar plot.
Input Files
GTF File:

Contains gene annotations (gene_id, gene_name, gene_type).
Example: gencode.v35.basic.annotation.gtf.gz.
Download from Gencode.
Raw Count Matrix:

A tab-delimited file with gene expression counts across samples. The first column should contain gene identifiers, followed by columns for each sample.
Example: gene_counts.txt.
Required Libraries
pandas: Data manipulation
numpy: Numerical operations
matplotlib: Plotting
seaborn: Advanced visualizations
Pipeline Execution
!!! IMPORTANT !!!
Ensure the GTF file and raw counts matrix file are in the 00_input_files_gtf_raw_counts directory.

Run the pipeline in bash shell using the following command: 
python 00_RNAseq_RAW_COUNTS_Quality_Control.py


The pipeline will generate the following output files or folders:

01_RESULTS_COVERAGE_BARPLOTS: A folder that contains barplots for each column (sample) of the input matrix.
01_COVERAGE_PER_SAMPLE_TABLE_CRC.SW.mRNA.symbol.count.txt.csv: A CSV file containing the coverage (read counts) for each sample. This file is useful for recognizing low-coverage samples.
02_RESULTS_COUNT_DISTRIBUTION_BAXPLOTS: A folder containing the count distribution barplot for each sample. This is useful for verifying if the input matrix is raw or normalized.
03_PIE_CHART_COUNTS_PERCENTAGE_CRC.SW.mRNA.symbol.count.txt.png: A PNG file containing the percentage of counts for each category. This is useful for assessing if the dataset is complete or if rows (genes) have been deleted.
04_BAR_PLOT_PERCENTAGE_OF_GENE_TYPE_IN_RAW_DATA_CRC.SW.mRNA.symbol.count.txt.tiff: A TIFF file visualizing the gene names in the symbol column, divided per category.
