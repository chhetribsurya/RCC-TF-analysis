#-------------------------------------------------------------------------------------------------------------------------------------------

# SNP-GENE Association for single TF

#-------------------------------------------------------------------------------------------------------------------------------------------

import pandas as pd
import pybedtools

# Define the GWAS p-value threshold for significance
p_value_threshold = 5e-8
#p_value_threshold = 1e-6

# Step 1: Load TF peak data
tf_file = "tf_both_epigenomesignif/KLF10_ID_293.bed"
tf_df = pd.read_csv(tf_file, sep="\t", header=None, usecols=[0, 1, 2, 9], names=["Chrom", "Start", "End", "TF_CellType"])

# Split 'TF_CellType' to extract TF name and CellType, then replace "293" with "HEK293"
tf_df["TF_Name"] = tf_df["TF_CellType"].str.split("_ID_").str[0]  # Extract TF name
tf_df["CellType_TF"] = tf_df["TF_CellType"].str.split("_ID_").str[1]  # Extract CellType
tf_df["CellType_TF"] = tf_df["CellType_TF"].replace("293", "HEK293")

# Drop the original 'TF_CellType' column
tf_df.drop(columns=["TF_CellType"], inplace=True)

# Step 2: Load SNP data, filter by p-value, and convert to BED format
snp_file = "gwas_rcc/MultiAnc_RCC.tsv"
snp_df = pd.read_csv(snp_file, sep="\t")

# Filter SNPs based on the p-value threshold
significant_snp_df = snp_df[snp_df["p"] < p_value_threshold]

# Create BED format with 0-based start position
snp_df_bed = significant_snp_df[["CHR", "BP", "SNP", "A1", "A2", "p"]]
snp_df_bed["Start"] = snp_df_bed["BP"] - 1  # Convert to 0-based start position for BED format
snp_df_bed["End"] = snp_df_bed["BP"]  # Use BP as the end position
snp_df_bed = snp_df_bed[["CHR", "Start", "End", "SNP", "A1", "A2", "p"]]
snp_df_bed.columns = ["Chrom", "Start", "End", "SNP", "A1", "A2", "p"]
snp_df_bed["Chrom"] = "chr" + snp_df_bed["Chrom"].astype(str)

# Step 3: Load ABC E2G data with specific columns
e2g_file = "e2g_model/RCC_ENCFF951UIR.bed"
e2g_df = pd.read_csv(e2g_file, sep="\t", usecols=["#chr", "start", "end", "TargetGene", "TargetGeneEnsemblID", "TargetGeneTSS", "CellType"],
                     dtype={"#chr": str, "start": int, "end": int, "TargetGene": str, "TargetGeneEnsemblID": str, "TargetGeneTSS": int, "CellType": str})
e2g_df.columns = ["Chrom", "Start", "End", "TargetGene", "TargetGeneEnsemblID", "TargetGeneTSS", "CellType_E2G"]


# Step 4: Convert DataFrames to BEDTools objects
tf_bed = pybedtools.BedTool.from_dataframe(tf_df)
snp_bed = pybedtools.BedTool.from_dataframe(snp_df_bed)
e2g_bed = pybedtools.BedTool.from_dataframe(e2g_df)

# Step 5: Intersect TF peaks with E2G data
tf_e2g_intersect = tf_bed.intersect(e2g_bed, wa=True, wb=True)
tf_e2g_df = tf_e2g_intersect.to_dataframe(names=["Chrom_TF", "Start_TF", "End_TF", "TF_Name", "CellType_TF",
                                                 "Chrom_E2G", "Start_E2G", "End_E2G", "TargetGene",
                                                 "TargetGeneEnsemblID", "TargetGeneTSS", "CellType_E2G"])

# Step 6: Intersect the resulting TF-E2G data with SNP data
tf_e2g_snp_intersect = tf_e2g_intersect.intersect(snp_bed, wa=True, wb=True)
tf_e2g_snp_df = tf_e2g_snp_intersect.to_dataframe(names=["Chrom_TF", "Start_TF", "End_TF", "TF_Name", "CellType_TF",
                                                         "Chrom_E2G", "Start_E2G", "End_E2G", "TargetGene",
                                                         "TargetGeneEnsemblID", "TargetGeneTSS", "CellType_E2G",
                                                         "Chrom_SNP", "Start_SNP", "End_SNP", "SNP",
                                                         "A1", "A2", "p"])

# Step 7: Create a final DataFrame with the required columns
final_df = tf_e2g_snp_df[["TF_Name", "SNP", "TargetGene", "p", "CellType_TF", "CellType_E2G"]]
final_df.columns = ["TF", "SNP", "Gene", "GWAS_p-value", "CellType_TF", "CellType_E2G"]

# Step 8: Display the final DataFrame
print(final_df)

# Save the final DataFrame to a file
final_df.to_csv("tf_snp_gene_association.csv", index=False)


#----------------------------------------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------

# SNP-GENE Association for multiple TFs

#-------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------


import pandas as pd
import pybedtools
import glob
import os
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import BoundaryNorm

# Define the GWAS p-value threshold for significance
p_value_threshold = 5e-8

# Output directories
output_dir = "output_results"
output_dir = "output_results_eithercase"
plots_dir = os.path.join(output_dir, "plots")

# Create directories if they do not exist
os.makedirs(output_dir, exist_ok=True)
os.makedirs(plots_dir, exist_ok=True)

# Step 1: Load TF peak data
tf_files = glob.glob("tf_both_epigenomesignif/*.bed")
tf_files = glob.glob("tf_either_epigenomesignif/*.bed")
tf_dfs = []

print("Processing TF files...")
for i, tf_file in enumerate(tf_files, start=1):
    print(f"Processing TF file {i}/{len(tf_files)}: {tf_file}")
    tf_df = pd.read_csv(tf_file, sep="\t", header=None, usecols=[0, 1, 2, 9], names=["Chrom", "Start", "End", "TF_CellType"])
    # Split 'TF_CellType' to extract TF name and CellType
    tf_df["TF_Name"] = tf_df["TF_CellType"].str.split("_ID_").str[0]  # Extract TF name
    tf_df["CellType_TF"] = tf_df["TF_CellType"].str.split("_ID_").str[1]  # Extract CellType
    tf_df["CellType_TF"] = tf_df["CellType_TF"].replace("293", "HEK293")
    # Drop the original 'TF_CellType' column
    tf_df.drop(columns=["TF_CellType"], inplace=True)
    tf_dfs.append(tf_df)

# Combine all TF data into a single DataFrame
tf_combined_df = pd.concat(tf_dfs, ignore_index=True)

# Step 2: Load SNP data, filter by p-value, and convert to BED format
snp_file = "gwas_rcc/MultiAnc_RCC.tsv"
snp_df = pd.read_csv(snp_file, sep="\t")

# Filter SNPs based on the p-value threshold
significant_snp_df = snp_df[snp_df["p"] < p_value_threshold]

# Create BED format with 0-based start position
snp_df_bed = significant_snp_df[["CHR", "BP", "SNP", "A1", "A2", "p"]]
snp_df_bed["Start"] = snp_df_bed["BP"] - 1  # Convert to 0-based start position for BED format
snp_df_bed["End"] = snp_df_bed["BP"]  # Use BP as the end position
snp_df_bed = snp_df_bed[["CHR", "Start", "End", "SNP", "A1", "A2", "p"]]
snp_df_bed.columns = ["Chrom", "Start", "End", "SNP", "A1", "A2", "p"]
snp_df_bed["Chrom"] = "chr" + snp_df_bed["Chrom"].astype(str)

# Step 3: Load ABC E2G data with specific columns
e2g_files = glob.glob("e2g_model/*.bed")
e2g_dfs = []

print("\nProcessing E2G files...")
for i, e2g_file in enumerate(e2g_files, start=1):
    print(f"Processing E2G file {i}/{len(e2g_files)}: {e2g_file}")
    e2g_df = pd.read_csv(e2g_file, sep="\t", usecols=["#chr", "start", "end", "TargetGene", "TargetGeneEnsemblID", "TargetGeneTSS", "CellType"],
                         dtype={"#chr": str, "start": int, "end": int, "TargetGene": str, "TargetGeneEnsemblID": str, "TargetGeneTSS": int, "CellType": str})
    e2g_df.columns = ["Chrom", "Start", "End", "TargetGene", "TargetGeneEnsemblID", "TargetGeneTSS", "CellType_E2G"]
    e2g_dfs.append(e2g_df)

# Combine all E2G data into a single DataFrame
e2g_combined_df = pd.concat(e2g_dfs, ignore_index=True)

# Step 4: Convert DataFrames to BEDTools objects
tf_bed = pybedtools.BedTool.from_dataframe(tf_combined_df)
snp_bed = pybedtools.BedTool.from_dataframe(snp_df_bed)
e2g_bed = pybedtools.BedTool.from_dataframe(e2g_combined_df)

# Step 5: Intersect TF peaks with E2G data
tf_e2g_intersect = tf_bed.intersect(e2g_bed, wa=True, wb=True)
tf_e2g_df = tf_e2g_intersect.to_dataframe(names=["Chrom_TF", "Start_TF", "End_TF", "TF_Name", "CellType_TF",
                                                 "Chrom_E2G", "Start_E2G", "End_E2G", "TargetGene",
                                                 "TargetGeneEnsemblID", "TargetGeneTSS", "CellType_E2G"])

# Step 6: Intersect the resulting TF-E2G data with SNP data
tf_e2g_snp_intersect = tf_e2g_intersect.intersect(snp_bed, wa=True, wb=True)
tf_e2g_snp_df = tf_e2g_snp_intersect.to_dataframe(names=["Chrom_TF", "Start_TF", "End_TF", "TF_Name", "CellType_TF",
                                                         "Chrom_E2G", "Start_E2G", "End_E2G", "TargetGene",
                                                         "TargetGeneEnsemblID", "TargetGeneTSS", "CellType_E2G",
                                                         "Chrom_SNP", "Start_SNP", "End_SNP", "SNP",
                                                         "A1", "A2", "p"])

# Step 7: Create a final DataFrame with the required columns
final_df = tf_e2g_snp_df[["TF_Name", "SNP", "TargetGene", "p", "CellType_TF", "CellType_E2G"]]
final_df.columns = ["TF", "SNP", "Gene", "GWAS_p-value", "CellType_TF", "CellType_E2G"]

# Step 8: Display the final DataFrame
print("\nFinal DataFrame:")
print(final_df)

# Save the final DataFrame to a file
final_df_path = os.path.join(output_dir, "all_tf_snp_gene_association.tsv")
final_df.to_csv(final_df_path, sep="\t", index=False)


###############
# For local loading of file:
final_df = pd.read_csv("/Users/suryachhetri/Dropbox/for_diptavo/ABC_contact/all_tf_snp_gene_association.tsv", sep="\t") #both epigenome signif case
final_df = pd.read_csv("/Users/suryachhetri/Dropbox/for_diptavo/ABC_contact/output_results_eithercase/all_tf_snp_gene_association.tsv", sep="\t") #either epigenome signif case


# Step 9: Plotting and Insights

# # Plot the distribution of GWAS p-values for significant SNPs
# plt.figure(figsize=(10, 6))
# sns.histplot(final_df["GWAS_p-value"], bins=50, kde=True)
# plt.title("Distribution of GWAS p-values for Significant SNPs")
# plt.xlabel("GWAS p-value")
# plt.ylabel("Frequency")
# snp_gwas_histogram = os.path.join(plots_dir, "snp_gwas_histogram.pdf")
# plt.savefig(snp_gwas_histogram)
# plt.xscale("log")

# Convert GWAS p-values to -log10 scale
import numpy as np
final_df_hist = final_df.copy()
final_df_hist["-log10(GWAS_p-value)"] = -np.log10(final_df_hist["GWAS_p-value"])

# Plot the distribution of -log10(GWAS p-values) for significant SNPs
plt.figure(figsize=(10, 6))
sns.histplot(final_df_hist["-log10(GWAS_p-value)"], bins=50, kde=True)
plt.title("Distribution of -log10(GWAS p-values) for Significant SNPs")
plt.xlabel("-log10(GWAS p-value)")
plt.ylabel("Frequency")

# Save the histogram
snp_gwas_histogram = os.path.join(plots_dir, "snp_gwas_histogram1.pdf")
plt.savefig(snp_gwas_histogram)
plt.close()

#plt.show()

# Count the number of unique TFs, SNPs, and Genes
num_tfs = final_df["TF"].nunique()
num_snps = final_df["SNP"].nunique()
num_genes = final_df["Gene"].nunique()

print(f"Number of unique TFs: {num_tfs}")
print(f"Number of unique SNPs: {num_snps}")
print(f"Number of unique Genes: {num_genes}")

# Count the number of associations per TF
tf_association_counts = final_df["TF"].value_counts()

plt.figure(figsize=(12, 8))
sns.barplot(x=tf_association_counts.index, y=tf_association_counts.values)
plt.title("Number of Associations per TF")
plt.xlabel("Transcription Factor")
plt.ylabel("Number of Associations")
plt.xticks(rotation=90)
tf_association_counts_path = os.path.join(plots_dir, "tf_association_counts.pdf")
plt.savefig(tf_association_counts_path)
#plt.show()

# Biological insights
# - The number of unique TFs shows the diversity of regulatory elements interacting with significant SNPs and target genes.
# - The number of SNPs reflects the genetic variations potentially impacting RCC risk through these regulatory mechanisms.
# - The number of target genes indicates the breadth of regulatory influence exerted by these TFs and SNPs.

# Correlation analysis between TFs and genes
tf_gene_correlation = final_df.groupby(["TF", "Gene"]).size().unstack(fill_value=0)

plt.figure(figsize=(14, 10))
sns.heatmap(tf_gene_correlation, cmap="YlGnBu")
plt.title("Heatmap of TF-Gene Associations")
plt.xlabel("Gene")
plt.ylabel("Transcription Factor")
heatmap_path = os.path.join(plots_dir, "tf_gene_heatmap.pdf")
plt.savefig(heatmap_path)

# Debug: Inspect DataFrame
print("TF-Gene Correlation Matrix Shape:", tf_gene_correlation.shape)
print("Top TF-Gene Associations:")
print(tf_gene_correlation.head())

# Plot the heatmap of TF-Gene Associations with annotations
plt.figure(figsize=(14, 10))
sns.heatmap(tf_gene_correlation, cmap="YlGnBu", cbar_kws={"label": "Number of Associations"})
plt.title("Heatmap of TF-Gene Associations")
plt.xlabel("Gene")
plt.ylabel("Transcription Factor")
# plt.xticks(rotation=90)
# plt.yticks(rotation=0)

# Set all labels on x and y axes
plt.xticks(ticks=range(tf_gene_correlation.shape[1]), labels=tf_gene_correlation.columns, rotation=90, fontsize=3)
plt.yticks(ticks=range(tf_gene_correlation.shape[0]), labels=tf_gene_correlation.index, rotation=0, fontsize=7)

heatmap_path = os.path.join(plots_dir, "tf_gene_heatmap_ALL1.pdf")
plt.savefig(heatmap_path)
#plt.close()
#plt.show()


# Clustering might help visualize patterns better in sparse data
plt.figure(figsize=(14, 10))
sns.clustermap(tf_gene_correlation, cmap="YlGnBu", figsize=(14, 10), cbar_kws={"label": "Number of Associations"})

clustermap = os.path.join(plots_dir, "tf_gene_heatmap_clustermap_ALL.pdf")
plt.savefig(clustermap)

# Create the clustermap
g = sns.clustermap(
    tf_gene_correlation,
    cmap="YlGnBu",
    figsize=(14, 10),
    cbar_kws={"label": "Number of Associations"}
)

# Set all labels on x and y axes based on the reordered indices
g.ax_heatmap.set_xticks(range(tf_gene_correlation.shape[1]))
g.ax_heatmap.set_xticklabels(tf_gene_correlation.columns[g.dendrogram_col.reordered_ind], rotation=90, fontsize=3)

g.ax_heatmap.set_yticks(range(tf_gene_correlation.shape[0]))
g.ax_heatmap.set_yticklabels(tf_gene_correlation.index[g.dendrogram_row.reordered_ind], rotation=0, fontsize=7)

# Save the clustermap
clustermap_path = os.path.join(plots_dir, "tf_gene_heatmap_clustermap_ALL1.pdf")
g.savefig(clustermap_path)
plt.close()


# Determine custom breakpoints at intervals of 50
min_value = tf_gene_correlation.min().min()
max_value = tf_gene_correlation.max().max()

# Create breakpoints with interval of 50
breakpoints = np.arange(min_value, max_value + 25, 50)
norm = BoundaryNorm(boundaries=breakpoints, ncolors=256)

# Create the clustermap with more color breakpoints
g = sns.clustermap(
    tf_gene_correlation,
    cmap="YlGnBu",
    norm=norm,
    figsize=(14, 10),
    cbar_kws={"label": "Number of Associations"}
)

# Set all labels on x and y axes based on the reordered indices
g.ax_heatmap.set_xticks(range(tf_gene_correlation.shape[1]))
g.ax_heatmap.set_xticklabels(tf_gene_correlation.columns[g.dendrogram_col.reordered_ind], rotation=90, fontsize=3)

g.ax_heatmap.set_yticks(range(tf_gene_correlation.shape[0]))
g.ax_heatmap.set_yticklabels(tf_gene_correlation.index[g.dendrogram_row.reordered_ind], rotation=0, fontsize=7)

# Save the clustermap
clustermap_path = os.path.join(plots_dir, "tf_gene_heatmap_clustermap_ALL_breaks.pdf")
g.savefig(clustermap_path)
plt.close()
