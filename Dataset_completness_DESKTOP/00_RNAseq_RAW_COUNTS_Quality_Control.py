#!/usr/bin/env python3
# Import required libraries
import pandas as pd
import gzip
import os, shutil
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Definitions of gene categories
# These categories represent different gene types, like coding and non-coding genes.
CODING_GENES = [
    'protein_coding', 'TR_C_gene', 'TR_J_gene', 'TR_V_gene',
    'TR_D_gene', 'IG_C_gene', 'IG_J_gene', 'IG_V_gene', 'IG_D_gene'
]

PSEUDOGENES_NON_CODING = [
    'transcribed_unitary_pseudogene', 'TR_J_pseudogene', 'translated_unprocessed_pseudogene',
    'unprocessed_pseudogene', 'processed_pseudogene', 'IG_J_pseudogene', 'IG_C_pseudogene',
    'IG_V_pseudogene', 'translated_processed_pseudogene', 'polymorphic_pseudogene',
    'unitary_pseudogene', 'pseudogene', 'IG_pseudogene', 'TR_V_pseudogene'
]

RIBOSOMAL_RNA_NON_CODING = ['rRNA', 'rRNA_pseudogene']
MITOCHONDRIAL_RNA_NON_CODING = ['Mt_rRNA', 'Mt_tRNA']
OTHER_NON_CODING_RNA = [
    'lncRNA', 'TEC', 'scaRNA', 'misc_RNA', 'scRNA',
    'sRNA', 'transcribed_unprocessed_pseudogene', 'snoRNA',
    'snRNA', 'miRNA', 'vault_RNA', 'ribozyme'
]

#########################

# Function to extract gene characteristics from a GTF file
def extract_gene_characteristics(gtf_file):
    deduplicated_gene_characteristic = set()
    with gzip.open(gtf_file, 'rt') as file:
        for line in file:
            if line.startswith('##'):  # Skip header lines
                continue
            gene_type = line.split('gene_type')[1].split(';')[0].split('"')[1]
            gene_id = line.split('gene_id')[1].split(';')[0].split('.')[0].split('"')[1]
            gene_name = line.split('gene_name')[1].split(';')[0].split('"')[1]
            deduplicated_gene_characteristic.add(f"{gene_name}\t{gene_type}\t{gene_id}")

    return sorted(list(deduplicated_gene_characteristic))

# Function to process gene expression data from a file
def process_gene_file(file_path):
    df = pd.read_csv(file_path, sep='\t')  # Read the gene expression data
    df.columns = df.columns.str.upper()  # Standardize column names to uppercase
    symbol_list = df['SYMBOL'].tolist()  # List of gene symbols
    # Remove duplicates from the list of gene symbols
    dataset_gene_deduplicated = set(map(str.strip, symbol_list))
    print(len(dataset_gene_deduplicated))
    print(f"Total genes from the file: {os.path.basename(file_path)}: {len(symbol_list)}")
    print(f"Unique genes from the file: {os.path.basename(file_path)}: {len(dataset_gene_deduplicated)}")
    if len(symbol_list) == len(dataset_gene_deduplicated):
        print('There are no duplicated genes in the DataFrame.')
    else:
        print('There are duplicated genes in the DataFrame.')
    return list(dataset_gene_deduplicated), df

# Function to create a bar plot showing the coverage per sample
def Coverage_per_sample_barplot(dataset, output_dir, file_name):
    
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)  # Remove the directory if it exists
    os.makedirs(output_dir)  # Create the output directory
    
    sum_df = pd.DataFrame(dataset.drop('SYMBOL', axis=1).sum()).reset_index()  # Sum values for each column
    sum_df.columns = ['Column', 'Coverage']
    
    for i in range(0, len(sum_df), 20):  # Generate plots in blocks of 20 columns
        subset = sum_df.iloc[i:i+20]
        plt.figure(figsize=(10, 6))
        plt.bar(subset['Column'], subset['Coverage'], color='skyblue')
        plt.title(f'Coverage per Sample (Plot {i//20 + 1})')
        plt.xlabel('Column')
        plt.ylabel('Coverage')
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
        plt.savefig(f"{output_dir}/Coverage_plot_{file_name}_{i//20 + 1}.png")
        plt.close()  # Close the plot to avoid overlapping with subsequent plots
    
    return sum_df

# Function to generate box plots of log2-transformed counts per sample
def boxplot_log_distribution_per_sample(dataset, output_dir, file_name):
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)  # Remove the directory if it exists
    os.makedirs(output_dir)  # Create the output directory
    
    dataset_no_symbol = dataset.drop('SYMBOL', axis=1)
    
    dataset_log = np.log10(dataset_no_symbol + 1)  # Log2 transformation (adding 1 to avoid log(0))
    
    for i in range(0, dataset_log.shape[1], 20):  # Generate boxplots in blocks of 20 columns
        subset = dataset_log.iloc[:, i:i+20]
        plt.figure(figsize=(10, 6))
        subset.boxplot()
        plt.title(f'Log2 RAW Count Distribution (Plot {i//20 + 1})')
        plt.xlabel('Samples')
        plt.ylabel('Log2(Counts + 1)')
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
        plt.savefig(f"{output_dir}/Boxplot_Log2_Distribution_{i//20 + 1}_{file_name}.png")
        plt.close()

# Function to calculate row-wise mean, median, and standard deviation for a DataFrame
def calculate_row_mean_median_SD(df):
    if 'SYMBOL' not in df.columns:
        raise ValueError("'SYMBOL' column is missing in the DataFrame")
    
    df_numeric = df.drop('SYMBOL', axis=1)  # Drop the SYMBOL column
    means = df_numeric.mean(axis=1)
    medians = df_numeric.median(axis=1)
    std_devs = df_numeric.std(axis=1)
    
    result_df = pd.DataFrame({
        'SYMBOL': df['SYMBOL'],
        'Mean': means,
        'Median': medians,
        'StdDev': std_devs
    })
    
    return result_df

# Function to calculate row-wise sum of gene expression data
def calculate_row_sum(df):
    if 'SYMBOL' not in df.columns:
        raise ValueError("'SYMBOL' column is missing in the DataFrame")
    
    df_numeric = df.drop('SYMBOL', axis=1)
    df_numeric = df_numeric.apply(pd.to_numeric, errors='coerce')
    
    row_sums = df_numeric.sum(axis=1, skipna=True)
    
    result_df = pd.DataFrame({
        'SYMBOL': df['SYMBOL'],
        'RowSum': row_sums
    })
    
    total_sum = row_sums.sum()
    
    return result_df, total_sum

# Function to process gene data and categorize genes based on gene type
def process_gene_data(dataset_gene_list, gene_name_type_ID):
    list_for_table = []
    for gene_name in dataset_gene_list:
        gene_features = {
            'SYMBOL': gene_name.strip(),
            'Protein_coding_genes': 0,
            'Pseudogenes': 0,
            'Ribosomal_RNA': 0,
            'Mitochondrial_RNA': 0,
            'Other_non_coding_RNA': 0
        }
        list_for_table.append(gene_features)

        for entry in gene_name_type_ID:
            name, gene_type, _ = entry.split('\t')
            if name == gene_name.strip():
                if gene_type in CODING_GENES:
                    list_for_table[-1]['Protein_coding_genes'] += 1
                elif gene_type in PSEUDOGENES_NON_CODING:
                    list_for_table[-1]['Pseudogenes'] += 1
                elif gene_type in RIBOSOMAL_RNA_NON_CODING:
                    list_for_table[-1]['Ribosomal_RNA'] += 1
                elif gene_type in MITOCHONDRIAL_RNA_NON_CODING:
                    list_for_table[-1]['Mitochondrial_RNA'] += 1
                elif gene_type in OTHER_NON_CODING_RNA:
                    list_for_table[-1]['Other_non_coding_RNA'] += 1
    df = pd.DataFrame(list_for_table)
    return df

# Function to merge DataFrames based on the 'SYMBOL' column and check for NaN values
def merge_dataframes(mean_median_SD_, gene_data_):
    # Merge the two DataFrames on the 'SYMBOL' column using a left join
    merged_df = pd.merge(mean_median_SD_, gene_data_, on='SYMBOL', how='left')
    
    # Check if there are any NaN values in the merged DataFrame
    if merged_df.isna().any().any():
        print("Warning: There are NaN values in the merged DataFrame!")
        # Optionally, display the rows containing NaN values
        print(merged_df[merged_df.isna().any(axis=1)])  # Print rows with NaN values
    
    # Return the merged DataFrame
    return merged_df


def create_pie_chart(df, output_file):
    """
    Creates a pie chart showing the percentage distribution of gene categories (e.g., Protein Coding, Pseudogenes).
    Additionally, returns a DataFrame with the total row counts per category.

    Parameters:
    df (pd.DataFrame): DataFrame containing row sums and gene categories.
    output_file (str): Path to save the output pie chart image.

    Returns:
    pd.DataFrame: A DataFrame with total counts for each category.
    
    Raises:
    ValueError: If the 'SYMBOL' column is missing from the DataFrame.
    """
    # Check if 'SYMBOL' column is present in the DataFrame
    if 'SYMBOL' not in df.columns:
        raise ValueError("'SYMBOL' column is missing in the DataFrame")
    
    # Categories corresponding to the columns in the DataFrame
    categories = ['Protein_coding_genes', 'Pseudogenes', 'Ribosomal_RNA', 'Mitochondrial_RNA', 'Other_non_coding_RNA']
    category_names = ['Protein Coding Genes', 'Pseudogenes', 'Ribosomal RNA', 'Mitochondrial RNA', 'Other Non-Coding RNA']
    
    # Initialize category sums to zero for each category
    category_sums = {category: 0 for category in categories}
    total_sum = 0  # Variable to track the total sum of expression counts
    
    # Loop through each row to calculate sums for each category
    for idx, row in df.iterrows():
        row_sum = row['RowSum']  # Get the sum of expression values for the current gene
        
        categories_count = 0  # Counter to track how many categories this gene belongs to
        
        # Check each category and add the row sum to the respective category if the gene belongs to it
        for category in categories:
            if row[category] == 1:  # If the gene belongs to this category
                category_sums[category] += row_sum  # Add the RowSum to the corresponding category
                categories_count += 1
        
        # If the gene belongs to at least one category, add its contribution to the total sum
        if categories_count > 0:
            total_sum += row_sum * categories_count  # Multiply by the number of categories the gene belongs to
    
    # Calculate the percentage distribution for each category
    category_percentages = {
        category: (value / total_sum) * 100 if total_sum > 0 else 0
        for category, value in category_sums.items()
    }

    # Create the pie chart using the calculated sums
    sizes = list(category_sums.values())  # Sizes of the pie chart slices
    
    # Create the pie chart with reduced radius (using the 'radius' argument)
    plt.figure(figsize=(6, 6))  # Keep the figure size compact
    wedges, texts = plt.pie(
        sizes, 
        startangle=140,  # Start the chart at a specific angle for aesthetics
        colors=['#66b3ff', '#99ff99', '#ffcc99', '#ff6666', '#ffb3e6'],  # Colors for each category
        wedgeprops={'edgecolor': 'black'},  # Add black edges to the wedges
        labeldistance=1.1,  # Reduce label distance to make the chart more compact
        radius=0.75  # Reduce the radius of the pie chart (default is 1)
    )

    # Create a legend with category names, counts, and percentages
    legend_labels = [
        f"{name}: ({category_percentages[cat]:.1f}%)" 
        for name, cat in zip(category_names, categories)
    ]
    
    # Position the legend to the left of the chart and adjust font size
    plt.legend(
        wedges, 
        legend_labels, 
        title="Categories", 
        loc="upper left", 
        bbox_to_anchor=(-0.4, 1),  # Move the legend outside to the left
        fontsize=9  # Smaller font for the legend
    )
    
    # Set the title of the chart
    plt.title("Read counts percentage per category", fontsize=12)  # Smaller font for the title
    plt.axis('equal')  # Ensure the pie chart is a circle

    # Save the pie chart to a file
    plt.savefig(output_file, bbox_inches='tight')  # Save with tight bounding box for better spacing
    plt.close()  # Close the plot to avoid overlapping with other plots

    # Create a DataFrame to show the total counts for each category
    category_totals_df = pd.DataFrame(
        {'Category': category_names, 'Total Count': list(category_sums.values())}
    )

    # Return the DataFrame with total counts per category
    return category_totals_df


def calculate_percentages(df):
    # Calculate the total sum of expression counts across all categories
    total_sum = df[['Protein_coding_genes', 'Pseudogenes', 'Ribosomal_RNA', 'Mitochondrial_RNA', 'Other_non_coding_RNA']].sum().sum()
    
    # Calculate the percentage distribution for each category
    percentages = (df[['Protein_coding_genes', 'Pseudogenes', 'Ribosomal_RNA', 'Mitochondrial_RNA', 'Other_non_coding_RNA']].sum() / total_sum) * 100
    
    # Create a DataFrame to store the percentages for each category
    df_percentage = pd.DataFrame(percentages, columns=['Percentage'])
    
    # Return the DataFrame with percentage values
    return df_percentage

def plot_percentage_barplot(df_percentage, title='Barplot of Percentages'):
    # Reset index and rename the columns for easier plotting
    df_plot = df_percentage.reset_index()
    df_plot.columns = ['Column', 'Percentage']

    # Define a fixed color palette for the categories
    color_palette = {
        'Protein_coding_genes': '#1f77b4',  # Blue
        'Pseudogenes': '#ff7f0e',           # Orange
        'Ribosomal_RNA': '#2ca02c',         # Green
        'Mitochondrial_RNA': '#d62728',     # Red
        'Other_non_coding_RNA': '#9467bd'   # Purple
    }

    # Map the colors to each category in the DataFrame
    df_plot['Color'] = df_plot['Column'].map(color_palette)

    # Create a bar plot using seaborn
    plt.figure(figsize=(12, 8))

    # Create the bar plot, using the fixed color palette
    barplot = sns.barplot(x='Column', y='Percentage', data=df_plot, palette=color_palette)

    # Annotate the bars with their respective percentage values
    for p in barplot.patches:
        barplot.annotate(f'{p.get_height():.1f}%',  # Annotate the height of each bar
                         (p.get_x() + p.get_width() / 2., p.get_height()),  # Positioning of the annotation
                         ha='center', va='bottom',  # Center the text horizontally and align it vertically
                         fontsize=12, color='black', weight='bold')  # Style the annotation text

    # Set title and axis labels
    plt.title(title, fontsize=16, weight='bold')
    plt.xticks(rotation=45, fontsize=12)  # Rotate x-axis labels for readability
    plt.ylabel('Percentage (%)', fontsize=14)
    plt.xlabel('Columns', fontsize=14)

    # Apply tight layout to make sure everything fits within the figure area
    plt.tight_layout()

    # Save the bar plot as a .tiff image with high resolution (300 DPI)
    plt.savefig(f"{os.getcwd()}/{title.replace(' ', '_')}.tiff", dpi=300)



#########################

# Define paths for input files
current_directory = os.getcwd()

# Define Raw count matrix file name
## Important note: the first column of the matrix must be named 'SYMBOL' or 'symbol'
raw_matrix_file = 'CMCBSN_expectedcount_342.txt'  # Korean RNASeq matrix downlaoded from https://doi.org/10.5281/zenodo.8333650 (Lee J. et al. (2024))

########################################################################
## PATH INPUT FILES
# Path to the GTF file (Gene Transfer Format file)
gtf_file_path = os.path.join(current_directory, '00_input_files_gtf_raw_counts', 'gencode.v35.basic.annotation.gtf.gz')
# Path to the raw count matrix
dataframe = os.path.join(current_directory, '00_input_files_gtf_raw_counts', raw_matrix_file)
# dataframe_swedish = os.path.join(current_directory, '00_input_files_gtf_raw_counts', 'CMCBSN_expectedcount_342.txt')  # Alternative path



########################################################################
print('STEP_1. Extract gene_type attributes from the GTF file')
# 1. Extract gene_type attributes from the GTF file
gene_name_type_ID = extract_gene_characteristics(gtf_file_path)
print(f'The number of features in the GTF file is: {len(gene_name_type_ID)}')

#########################################################################
print('STEP_2. Process the raw count matrix file, extract the gene list and the expression dataframe')

# 2. Process the raw count matrix file, extract the gene list and the expression dataframe
dataset_gene_list, dataset = process_gene_file(dataframe)
#dataset.to_csv('raw_matrix_file_output.csv', index=False)  # Save the processed dataset as a CSV

########################################################################
print('STEP_3. Compute Coverage per sample by summing counts for each sample (each column in the raw matrix)')
## 3. Compute Coverage per sample by summing counts for each sample (each column in the raw matrix)
sum_df = Coverage_per_sample_barplot(dataset, "01_RESULTS_COVERAGE_BARPLOTS",raw_matrix_file )
sum_df.to_csv(f'01_COVERAGE_PER_SAMPLE_TABLE_{raw_matrix_file}.csv', index=False)  # Save coverage per sample
print(f'The 01_COVERAGE_PER_SAMPLE_TABLE_{raw_matrix_file}.csv is saved in the current directory. Check the file to manually verify if the coverage per sample is sufficient for the analysis.')

print('STEP_4. Compute the distribution of counts per sample and create a boxplot of the log-transformed counts')
## 4. Compute the distribution of counts per sample and create a boxplot of the log-transformed counts
boxplot_log_distribution_per_sample(dataset, '02_RESULTS_COUNT_DISTRIBUTION_BAXPLOTS', raw_matrix_file)

print('STEP_5. Sum the counts in each row (summarize the data at the gene level)')
## 5. Sum the counts in each row (summarize the data at the gene level)
row_coverage_df, total_coverage = calculate_row_sum(dataset)

print('STEP_6. Generate a dataframe with gene categories (5 categories), marking each gene as 1 if it belongs to the category, or 0 if it does not')
# 6. Generate a dataframe with gene categories (5 categories), marking each gene as 1 if it belongs to the category, or 0 if it does not
gene_data = process_gene_data(dataset_gene_list, gene_name_type_ID)
#print(gene_data)
# gene_data.to_csv(f'{raw_matrix_file}.csv', index=False)  # Optionally save the gene data to a CSV file
print(f'Length of gene data: {len(gene_data)}')

print('STEP_7. Merge the gene data dataframe with the row coverage dataframe and check for NaN values')

# 7. Merge the gene data dataframe with the row coverage dataframe and check for NaN values
merged_df = merge_dataframes(row_coverage_df, gene_data)
# merged_df.to_csv(f'MERGED_DF_{raw_matrix_file}_merged_output.csv', index=False)  # Optionally save the merged dataframe

print('STEP_8. Create a pie chart to show the distribution of counts per category')
# 8. Create a pie chart to show the distribution of counts per category
create_pie_chart(merged_df, f'03_PIE_CHART_COUNTS_PERCENTAGE_{raw_matrix_file}.png', )

print('STEP_9. Calculate the percentage of genes for each category in the raw counts dataframe')
# 9. Calculate the percentage of genes for each category in the raw counts dataframe
df_percentage = calculate_percentages(gene_data)

# Create a bar plot to show the percentage distribution of gene types in the raw counts data
plot_percentage_barplot(df_percentage, f'04_BAR_PLOT_PERCENTAGE_OF_GENE_TYPE_IN_RAW_DATA_{os.path.basename(dataframe)}')

