#!/usr/bin/env python
import pandas as pd
import os

# Define input directories and file paths
input_dir = "/Users/chhetribsurya/sc1238/datasets/projects/rcc_tf_project"  # Update this to the actual directory containing the input files
output_dir = "TF_RSID_Files"

# Ensure output directory exists
os.makedirs(output_dir, exist_ok=True)

reference_file = os.path.join(input_dir, "MultiAnc_RCC.tsv")  # RSID to genomic coordinates mapping
tf_matrix_file = os.path.join(input_dir, "293_output_matrix.txt")  # RSID by TF presence/absence matrix

# Load reference file (RSID to genomic coordinates and alleles)
print("Loading reference file...")
reference_df = pd.read_csv(reference_file, sep="\t", usecols=["SNP", "CHR", "BP", "A1", "A2"])
reference_df.columns = ["RSID", "CHR", "POS", "A1", "A2"]  # Rename columns

# Add "chr" prefix to chromosome number
reference_df["CHR"] = "chr" + reference_df["CHR"].astype(str)

# Load TF matrix file (RSID by TF presence/absence matrix)
print("Loading TF matrix file...")
tf_matrix_df = pd.read_csv(tf_matrix_file, sep="\t")

# Merge the TF matrix with reference file to get genomic coordinates and alleles
print("Merging TF matrix with reference file...")
merged_df = tf_matrix_df.merge(reference_df, left_on="snpID", right_on="RSID", how="left")

# Drop redundant columns and reorder
merged_df = merged_df.drop(columns=["RSID"])
column_order = ["CHR", "POS", "snpID", "A1", "A2"] + [col for col in merged_df.columns if col not in ["CHR", "POS", "snpID", "A1", "A2"]]
merged_df = merged_df[column_order]

# Create individual TF files with RSIDs per TF
tf_list = [col for col in merged_df.columns if col not in ["CHR", "POS", "snpID", "A1", "A2"]]

print(f"Generating TF-specific RSID files in '{output_dir}' directory...")
for i, tf in enumerate(tf_list, start=1):
    print(f"Processing TF {i}/{len(tf_list)}: {tf}...")  # Show progress
    tf_rsid_df = merged_df[merged_df[tf] == 1][["CHR", "POS", "snpID", "A1", "A2"]]  # Keep only RSIDs for which TF is present
    if not tf_rsid_df.empty:
        tf_rsid_df.columns = ["CHROM", "POS", "RSID", "A1", "A2"]  # Ensure correct column names
        tf_rsid_df.to_csv(f"{output_dir}/{tf}.txt", sep="\t", index=False, header=True)

print("\n" + "TF-specific RSID file generation completed!" + "\n")
