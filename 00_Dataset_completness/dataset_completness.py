#!/usr/bin/env python3

import pandas as pd
import gzip
import os
import matplotlib.pyplot as plt
import seaborn as sns

class GeneDatasetCompleteness:
    def __init__(self, gtf_file):
        self.coding_genes = [
            'protein_coding', 'TR_C_gene', 'TR_J_gene', 'TR_V_gene',
            'TR_D_gene', 'IG_C_gene', 'IG_J_gene', 'IG_V_gene', 'IG_D_gene'
        ]

        self.pseudogenes_non_coding = [
            'transcribed_unitary_pseudogene', 'TR_J_pseudogene', 'translated_unprocessed_pseudogene',
            'unprocessed_pseudogene', 'processed_pseudogene', 'IG_J_pseudogene', 'IG_C_pseudogene',
            'IG_V_pseudogene', 'translated_processed_pseudogene', 'polymorphic_pseudogene',
            'unitary_pseudogene', 'pseudogene', 'IG_pseudogene', 'TR_V_pseudogene'
        ]

        self.ribosomal_RNA_non_coding = ['rRNA', 'rRNA_pseudogene']
        self.mitochondrial_RNA_non_coding = ['Mt_rRNA', 'Mt_tRNA']
        self.other_non_coding_RNA = [
            'lncRNA', 'TEC', 'scaRNA', 'misc_RNA', 'scRNA',
            'sRNA', 'transcribed_unprocessed_pseudogene', 'snoRNA',
            'snRNA', 'miRNA', 'vault_RNA', 'ribozyme'
        ]

        self.gtf_file = gtf_file
        self.gene_name_type_ID = self.extract_gene_characteristics()

    def process_gene_file(self, file_path):
        df = pd.read_csv(file_path, sep='\t')
        df.columns = df.columns.str.upper()  # Standardize column names
        symbol_list = df['SYMBOL'].tolist()
        dataset_gene_deduplicated = set(map(str.strip, symbol_list))

        print(f"Total genes from the file: {os.path.basename(file_path)}: {len(symbol_list)}")
        print(f"Unique genes from the file: {os.path.basename(file_path)}: {len(dataset_gene_deduplicated)}")

        if len(symbol_list) == len(dataset_gene_deduplicated):
            print('There are no duplicated genes in the DataFrame.')
            return list(dataset_gene_deduplicated)
        else:
            print('There are duplicated genes in the DataFrame.')
            exit()

    def extract_gene_characteristics(self):
        deduplicated_gene_characteristic = set()
        with gzip.open(self.gtf_file, 'rt') as file:
            for line in file:
                if line.startswith('##'):  # Skip header lines
                    continue
                gene_type = line.split('gene_type')[1].split(';')[0].split('"')[1]
                gene_id = line.split('gene_id')[1].split(';')[0].split('.')[0].split('"')[1]
                gene_name = line.split('gene_name')[1].split(';')[0].split('"')[1]
                deduplicated_gene_characteristic.add(f"{gene_name}\t{gene_type}\t{gene_id}")

        return sorted(list(deduplicated_gene_characteristic))

    def process_gene_data(self, dataset_gene_list, output_filename):
        list_for_table = []
        for gene_name in dataset_gene_list:
            gene_features = {
                'Gene_Name': gene_name.strip(),
                'Protein_coding_genes': 0,
                'Pseudogenes': 0,
                'Ribosomal_RNA': 0,
                'Mitochondrial_RNA': 0,
                'Other_non_coding_RNA': 0
            }
            list_for_table.append(gene_features)

            for entry in self.gene_name_type_ID:
                name, gene_type, _ = entry.split('\t')
                if name == gene_name.strip():
                    if gene_type in self.coding_genes:
                        list_for_table[-1]['Protein_coding_genes'] += 1
                    elif gene_type in self.pseudogenes_non_coding:
                        list_for_table[-1]['Pseudogenes'] += 1
                    elif gene_type in self.ribosomal_RNA_non_coding:
                        list_for_table[-1]['Ribosomal_RNA'] += 1
                    elif gene_type in self.mitochondrial_RNA_non_coding:
                        list_for_table[-1]['Mitochondrial_RNA'] += 1
                    elif gene_type in self.other_non_coding_RNA:
                        list_for_table[-1]['Other_non_coding_RNA'] += 1

        df = pd.DataFrame(list_for_table)
        df.to_csv(output_filename, index=False)
        return df

    def calculate_percentages(self, df):
        total_sum = df[['Protein_coding_genes', 'Pseudogenes', 'Ribosomal_RNA', 'Mitochondrial_RNA', 'Other_non_coding_RNA']].sum().sum()
        percentages = (df[['Protein_coding_genes', 'Pseudogenes', 'Ribosomal_RNA', 'Mitochondrial_RNA', 'Other_non_coding_RNA']].sum() / total_sum) * 100
        df_percentage = pd.DataFrame(percentages, columns=['Percentage'])
        return df_percentage

    def plot_percentage_barplot(self, df_percentage, title='Barplot of Percentages'):
        df_plot = df_percentage.reset_index()
        df_plot.columns = ['Column', 'Percentage']

        # Fixed colors for each category
        color_palette = {
            'Protein_coding_genes': '#1f77b4',  # Blue
            'Pseudogenes': '#ff7f0e',           # Orange
            'Ribosomal_RNA': '#2ca02c',         # Green
            'Mitochondrial_RNA': '#d62728',     # Red
            'Other_non_coding_RNA': '#9467bd'   # Purple
        }

        # Map the colors based on the 'Column' name (this should match the categories in df_plot)
        df_plot['Color'] = df_plot['Column'].map(color_palette)

        plt.figure(figsize=(12, 8))

        # Use the 'hue' parameter to differentiate colors, this is the correct usage
        barplot = sns.barplot(x='Column', y='Percentage', data=df_plot, palette=color_palette)

        # Annotate bars with their respective percentages
        for p in barplot.patches:
            barplot.annotate(f'{p.get_height():.1f}%',
                             (p.get_x() + p.get_width() / 2., p.get_height()),
                             ha='center', va='bottom',
                             fontsize=12, color='black', weight='bold')

        plt.title(title, fontsize=16, weight='bold')
        plt.xticks(rotation=45, fontsize=12)
        plt.ylabel('Percentage (%)', fontsize=14)
        plt.xlabel('Columns', fontsize=14)

        plt.tight_layout()
        plt.savefig(f"{os.getcwd()}/{title.replace(' ', '_')}.tiff", dpi=300)

    def run_analysis(self, dataset_path, output_filename):
        dataset_gene_list = self.process_gene_file(dataset_path)
        print(f'The number of features in the GTF file is: {len(self.gene_name_type_ID)}')
        gene_data = self.process_gene_data(dataset_gene_list, output_filename)
        df_percentage = self.calculate_percentages(gene_data)
        print(f'{os.path.basename(output_filename)} DataFrame:')
        print(df_percentage)
        self.plot_percentage_barplot(df_percentage, title=f'Percentages of Genes in {os.path.basename(dataset_path)}')


# Set the current working directory
current_directory = os.getcwd()
os.chdir(current_directory)

# Initialize the class with the path to the GTF file in the current directory
gtf_file_path = os.path.join(current_directory, '00_input_files_gtf_raw_counts', 'gencode.v35.basic.annotation.gtf.gz')
gene_analysis = GeneDatasetCompleteness(gtf_file_path)

# Define the paths for the Swedish and Korean datasets in the current directory
dataframe_swedish = os.path.join(current_directory, '00_input_files_gtf_raw_counts', 'CRC.SW.mRNA.symbol.count.txt')  # insert the raw count files
dataframe_korean = os.path.join(current_directory, '00_input_files_gtf_raw_counts', 'CMCBSN_expectedcount_342.txt')  # insert the raw count files

# Run the analysis for both datasets
gene_analysis.run_analysis(dataframe_swedish, 'Swedish_completeness.csv')
gene_analysis.run_analysis(dataframe_korean, 'Korean_completeness.csv')

