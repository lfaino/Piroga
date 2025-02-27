
vcf_file_path = '/data/apolloData/genomes/rucola/vcfFiles/merge.all.vcf'

import pandas as pd
import numpy as np
import re
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import seaborn as sns
from multiprocessing import Pool, cpu_count


# Function to parse a single line of genotype data
def extract_genotype_line(line):
    headers = line.strip().split('\t')[:9]  # Keep the first 9 columns for reference
    genotype_data = line.strip().split('\t')[9:]
    extracted_data = [extract_genotype(genotype_str) for genotype_str in genotype_data]
    return headers + extracted_data


# Function to extract genotype from a string
def extract_genotype(genotype_str):
    match = re.match(r'^[01][\/\|][01]', genotype_str)
    if match:
        return match.group(0)
    return np.nan  # Return NaN if no match found


# Function to parse VCF file and extract genotype matrix
def parse_vcf(file_path):
    with open(file_path, 'r') as file:
        lines = [line for line in file if not line.startswith('##')]
    headers = lines[0].strip().split('\t')

    # Use multiprocessing to process genotype data in parallel
    with Pool(cpu_count()) as pool:
        genotype_data = pool.map(extract_genotype_line, lines[1:])

    # Create DataFrame
    genotype_df = pd.DataFrame(genotype_data, columns=headers[:9] + headers[9:])
    sample_names = headers[9:]  # Extract sample names
    genotype_df = genotype_df.iloc[:, 9:]  # Drop the first 9 columns to keep only genotype data

    # Replace genotypes with numerical values
    genotype_map = {'0/0': 0, '0/1': 1, '1/1': 2, '0|0': 0, '0|1': 1, '1|1': 2, './.': 0, '.|.': 0}
    genotype_df.replace(genotype_map, inplace=True)
    genotype_df = genotype_df.apply(pd.to_numeric)
    return genotype_df, sample_names


# Parse VCF file
genotype_df, sample_names = parse_vcf(vcf_file_path)

# Handle missing values (impute or remove)
genotype_df.fillna(genotype_df.mean(), inplace=True)

# Perform PCA for dimensionality reduction (optional but recommended)
pca = PCA(n_components=2)
genotype_pca = pca.fit_transform(genotype_df.T)

# Perform K-means clustering
# kmeans = KMeans(n_clusters=3, random_state=42)
# clusters = kmeans.fit_predict(genotype_pca)
#
# # Add cluster labels to the DataFrame
# genotype_df['Cluster'] = clusters

# Visualize the clustering
plt.figure(figsize=(10, 8))
scatter = sns.scatterplot(x=genotype_pca[:, 0], y=genotype_pca[:, 1], palette='viridis')
plt.title('Clustering of Samples')
plt.xlabel('PCA Component 1')
plt.ylabel('PCA Component 2')

# Annotate each point with the corresponding sample name
for i, txt in enumerate(sample_names):
    plt.annotate(txt, (genotype_pca[i, 0], genotype_pca[i, 1]), fontsize=9, alpha=0.7)


plt.savefig('/data/apolloData/genomes/rucola/vcfFiles/merge.all.vcf.pdf')
