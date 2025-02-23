import argparse
import os
import pysam
import numpy as np
import pandas as pd
import scipy.stats as stats
import matplotlib.pyplot as plt
from collections import defaultdict
from scipy.special import betaln, gammaln


# Function to count allele-specific reads from BAM
def count_allelic_reads(bam_file, snp_file, min_reads=5, apply_filter=True, effect_allele=False):
    """
    Count reads supporting each allele at given SNP locations and filter low-confidence sites.

    Args:
    - bam_file (str): Path to the BAM file.
    - snp_file (str): Path to the SNP coordinate text file.
    - min_reads (int): Minimum total reads required for a SNP to be considered.
    - apply_filter (bool): Whether to apply the read depth filter.
    - effect_allele (bool): If True, use A1 as the effect allele and A2 as the reference.

    Returns:
    - DataFrame containing allele-specific read counts at each SNP, including A1 and A2 if specified.
    """
    print("Loading SNPs from file...")
    bam = pysam.AlignmentFile(bam_file, "rb")
    snp_cols = ["CHROM", "POS"]
    
    # Check if A1 and A2 are in the SNP file
    snp_df = pd.read_csv(snp_file, sep='\t')
    if effect_allele:
        if "A1" not in snp_df.columns or "A2" not in snp_df.columns:
            raise ValueError("Error: A1 and A2 columns are required in the SNP file when using --effect-allele.\n"
                             "Please label the SNP file with 'A1' and 'A2' columns and retry.")
        snp_cols.extend(["A1", "A2"])
    
    snps = snp_df[snp_cols]
    total_snps = len(snps)
    print(f"Total SNPs to process: {total_snps}")

    allele_counts = defaultdict(lambda: {'A': 0, 'T': 0, 'C': 0, 'G': 0})

    for i, snp in snps.iterrows():
        chrom, pos = snp['CHROM'], int(snp['POS'])
        print(f"Processing SNP: {chrom}:{pos}")

        for pileupcolumn in bam.pileup(chrom, pos-1, pos, min_base_quality=0, min_mapping_quality=0):
            if pileupcolumn.pos == pos-1:
                for pileupread in pileupcolumn.pileups:
                    base = pileupread.alignment.query_sequence[pileupread.query_position]
                    if base in allele_counts[(chrom, pos)]:
                        allele_counts[(chrom, pos)][base] += 1

    df = pd.DataFrame([{"CHROM": k[0], "POS": k[1], **v} for k, v in allele_counts.items()])
    df['total_reads'] = df.iloc[:, 2:].sum(axis=1)

    if apply_filter:
        df = df[df['total_reads'] >= min_reads]  # Apply filtering if enabled

    # Add A1 and A2 columns if effect allele mode is enabled
    if effect_allele:
        df = df.merge(snps, on=["CHROM", "POS"], how="left")
        df["VAF (A1:A2)"] = df.apply(lambda row: f"{row['A1']}:{row['A2']} - {round(row[row['A1']] / (row[row['A1']] + row[row['A2']]), 2) if (row[row['A1']] + row[row['A2']]) > 0 else 'NA'}", axis=1)
    
    #if effect_allele:
    #    df = df.merge(snps, on=["CHROM", "POS"], how="left")
    #    df["VAF (A1::A2)"] = df.apply(lambda row: f"{row['A1']}::{row['A2']} - {row[row['A1']] / (row[row['A1']] + row[row['A2']]) if (row[row['A1']] + row[row['A2']]) > 0 else 'NA'}", axis=1)

    print(f"Processed {len(df)}/{total_snps} SNPs after filtering.")
    return df

## Function to perform binomial test
#def perform_binomial_test(df):
#    """
#    Perform binomial tests for allelic imbalance.
#
#    Args:
#    - df (DataFrame): DataFrame containing allele-specific read counts.
#
#    Returns:
#    - DataFrame with binomial p-values.
#    """
#    print("\nPerforming binomial test...")
#    # Ensure that only allele columns (A, T, C, G) are used for max_allele computation
#    allele_columns = ['A', 'T', 'C', 'G']
#    df['max_allele'] = df[allele_columns].idxmax(axis=1)
#
#    #df['max_allele'] = df.iloc[:, 2:].idxmax(axis=1)
#    df['max_allele_count'] = df.apply(lambda row: row[row['max_allele']], axis=1)
#
#    df['p_binomial'] = df.apply(lambda row: stats.binomtest(
#        row['max_allele_count'], row['total_reads'], p=0.5, alternative='two-sided').pvalue, axis=1)
#
#    df['significant_binomial'] = df['p_binomial'] < 0.05
#    print("Binomial test complete.")
#    return df


def perform_binomial_test(df, effect_allele=False):
    """
    Perform binomial tests for allelic imbalance.

    Args:
    - df (DataFrame): DataFrame containing allele-specific read counts.
    - effect_allele (bool): Whether to use A1 vs A2 for the binomial test.

    Returns:
    - DataFrame with binomial p-values.
    """
    print("\nPerforming binomial test...")

    if effect_allele:
        # Ensure A1 and A2 are present
        if 'A1' not in df.columns or 'A2' not in df.columns:
            raise ValueError("Error: --effect-allele is enabled, but A1 and A2 columns are missing in the dataset.")

        # Perform binomial test using A1 vs A2
        df['p_binomial'] = df.apply(lambda row: stats.binomtest(
            row[row['A1']], row[row['A1']] + row[row['A2']], p=0.5, alternative='two-sided').pvalue, axis=1)

        print("Binomial test (A1 vs A2) complete.")

    else:
        # Default: Use max_allele vs total_reads
        allele_columns = ['A', 'T', 'C', 'G']
        df['max_allele'] = df[allele_columns].idxmax(axis=1)
        df['max_allele_count'] = df.apply(lambda row: row[row['max_allele']], axis=1)

        df['p_binomial'] = df.apply(lambda row: stats.binomtest(
            row['max_allele_count'], row['total_reads'], p=0.5, alternative='two-sided').pvalue, axis=1)

        print("Binomial test (max allele vs total reads) complete.")

    df['significant_binomial'] = df['p_binomial'] < 0.05
    return df



def beta_binomial_pmf(k, n, alpha, beta):
    """
    Calculate the probability mass function of the beta-binomial distribution.
    
    Args:
    - k: number of successes
    - n: number of trials
    - alpha: first shape parameter of the beta distribution
    - beta: second shape parameter of the beta distribution
    
    Returns:
    - probability mass at k
    """
    try:
        return np.exp(
            gammaln(n + 1) - gammaln(k + 1) - gammaln(n - k + 1) +
            betaln(k + alpha, n - k + beta) - betaln(alpha, beta)
        )
    except Exception:
        return 0.0

def fit_beta_binomial(successes, trials, method='moment'):
    """
    Fit the parameters of a beta-binomial distribution.
    
    Args:
    - successes: array of success counts
    - trials: array of trial counts
    - method: fitting method ('moment' or 'mle')
    
    Returns:
    - tuple of (alpha, beta) parameters
    """
    if method == 'moment':
        # Method of moments estimation
        p = np.mean(successes / trials)
        var = np.var(successes / trials, ddof=1)
        
        if var < p * (1 - p) / np.mean(trials):
            # If variance is too small, fall back to binomial
            return 1000 * p, 1000 * (1 - p)
        
        # Calculate alpha and beta using method of moments
        moment = (p * (1 - p) - var) / (var - p * (1 - p) / np.mean(trials))
        alpha = p * moment
        beta = (1 - p) * moment
        
        return max(alpha, 0.01), max(beta, 0.01)
    
    else:
        raise ValueError("Only moment estimation is currently implemented")


def perform_beta_binomial_test(df, effect_allele=False):
    """
    Perform beta-binomial tests for allelic imbalance.
    
    Args:
    - df (DataFrame): DataFrame containing allele-specific read counts
    - effect_allele (bool): Whether to use A1 vs A2 for the test
    
    Returns:
    - DataFrame with beta-binomial p-values
    """
    print("\nPerforming beta-binomial test...")
    
    if effect_allele:
        if 'A1' not in df.columns or 'A2' not in df.columns:
            raise ValueError("Error: --effect-allele is enabled, but A1 and A2 columns are missing in the dataset.")
        

        # Get counts for A1 and A2 alleles
        # Extract counts by using the allele names to look up their counts
        successes = df.apply(lambda row: row[row['A1']], axis=1).values # equiv to a1_counts
        a2_counts = df.apply(lambda row: row[row['A2']], axis=1).values
        trials = successes + a2_counts
 
    else:
        # Use max allele vs total reads
        allele_columns = ['A', 'T', 'C', 'G']
        df['max_allele'] = df[allele_columns].idxmax(axis=1)
        df['max_allele_count'] = df.apply(lambda row: row[row['max_allele']], axis=1)
        
        successes = df['max_allele_count'].values
        trials = df['total_reads'].values
    
    # Fit beta-binomial parameters
    alpha, beta = fit_beta_binomial(successes, trials)
    
    # Calculate p-values
    def calc_pvalue(success, trial):
        # Calculate probabilities for all possible outcomes
        probs = np.array([beta_binomial_pmf(k, trial, alpha, beta) 
                         for k in range(trial + 1)])
        
        # Find outcomes with probability less than or equal to observed
        prob_observed = beta_binomial_pmf(success, trial, alpha, beta)
        extreme_outcomes = probs <= prob_observed
        
        # Two-sided p-value
        return np.sum(probs[extreme_outcomes])
    
    df['p_beta_binomial'] = [calc_pvalue(s, t) for s, t in zip(successes, trials)]
    df['significant_beta_binomial'] = df['p_beta_binomial'] < 0.05
    
    print(f"Beta-binomial test complete (fitted parameters: alpha={alpha:.3f}, beta={beta:.3f})")
    return df




## Function to perform beta-binomial test
#def beta_binomial_test(n, k, alpha=1, beta=1):
#    """
#    Perform a Beta-Binomial test for overdispersion in allelic imbalance.
#
#    Args:
#    - n (int): Total read depth.
#    - k (int): Count of the major allele.
#    - alpha (float): Prior beta distribution parameter.
#    - beta (float): Prior beta distribution parameter.
#
#    Returns:
#    - p-value from beta-binomial test.
#    """
#    log_p = betaln(k + alpha, n - k + beta) - betaln(alpha, beta)
#    return 1 - pow(10, log_p)

## Function to apply beta-binomial test to dataframe
#def perform_beta_binomial_test(df):
#    print("Performing beta-binomial test...")
#    df['p_beta_binomial'] = df.apply(lambda row: beta_binomial_test(row['total_reads'], row['max_allele_count']), axis=1)
#    df['significant_beta_binomial'] = df['p_beta_binomial'] < 0.05
#    print("Beta-binomial test complete.")
#    return df

## Function to perform beta-binomial test
#def perform_beta_binomial_test(df, effect_allele=False):
#    """
#    Perform beta-binomial test for allelic imbalance.
#
#    Args:
#    - df (DataFrame): DataFrame containing allele-specific read counts.
#    - effect_allele (bool): Whether to use A1 vs A2 for the beta-binomial test.
#
#    Returns:
#    - DataFrame with beta-binomial p-values.
#    """
#    print("\nPerforming beta-binomial test...")
#
#    if effect_allele:
#        # Ensure A1 and A2 are present
#        if 'A1' not in df.columns or 'A2' not in df.columns:
#            raise ValueError("Error: --effect-allele is enabled, but A1 and A2 columns are missing in the dataset.")
#
#        # Beta-binomial test using A1 vs A2
#        df['p_beta_binomial'] = df.apply(lambda row: beta_binomial_test(
#            row[row['A1']] + row[row['A2']], row[row['A1']]), axis=1)
#
#        print("Beta-binomial test (A1 vs A2) complete.")
#
#    else:
#        # Beta-binomial test using max allele vs total reads
#        if 'max_allele_count' not in df.columns:
#            raise ValueError("Error: 'max_allele_count' is missing. Ensure binomial test was run first.")
#
#        df['p_beta_binomial'] = df.apply(lambda row: beta_binomial_test(
#            row['total_reads'], row['max_allele_count']), axis=1)
#
#        print("Beta-binomial test (max allele vs total reads) complete.")
#
#    df['significant_beta_binomial'] = df['p_beta_binomial'] < 0.05
#    return df



# Function to save results in the specified output directory with a prefix
def save_results(df, outdir, outprefix):
    """
    Save results to separate files in the specified output directory with the given prefix.
    """
    os.makedirs(outdir, exist_ok=True)  # Ensure output directory exists

    full_output_path = os.path.join(outdir, f"{outprefix}_allelic_imbalance_results_full.tsv")
    binomial_only_output_path = os.path.join(outdir, f"{outprefix}_allelic_imbalance_results_binomial_only.tsv")

    df.to_csv(full_output_path, sep='\t', index=False)
    df.drop(columns=['p_beta_binomial', 'significant_beta_binomial']).to_csv(binomial_only_output_path, sep='\t', index=False)

    print(f"\nResults saved in {outdir}:")
    print(f"- {full_output_path}")
    print(f"- {binomial_only_output_path}")


# Function to generate and save the allelic imbalance plot
#def plot_allelic_imbalance(df, outdir, outprefix):
#    """
#    Generate and save a scatter plot showing allelic imbalance, with total read depth on the X-axis.
#
#    Args:
#    - df (pd.DataFrame): DataFrame containing allele-specific read counts and statistical results.
#    - outdir (str): Directory where the output plot should be saved.
#    - outprefix (str): Prefix for the output plot filename.
#
#    Saves:
#    - '{outdir}/{outprefix}_allelic_imbalance_plot.png': Scatter plot of allelic imbalance analysis.
#    """
#    print("Generating allelic imbalance plot...")
#
#    plt.figure(figsize=(10, 6))
#    plt.scatter(df['total_reads'], df['max_allele_count'] / df['total_reads'],
#                c=df['significant_binomial'].map({True: 'red', False: 'blue'}), alpha=0.6, label='Binomial Test')
#    plt.scatter(df['total_reads'], df['max_allele_count'] / df['total_reads'],
#                c=df['significant_beta_binomial'].map({True: 'green', False: 'blue'}), alpha=0.6, marker='x', label='Beta-Binomial Test')
#
#    plt.axhline(0.5, linestyle='--', color='black', linewidth=1)
#    plt.xlabel('Total Reads at SNP')
#    plt.ylabel('Major Allele Frequency')
#    plt.title('Allelic Imbalance Analysis by Read Depth')
#    plt.legend(title='Significance', labels=['Non-significant', 'Binomial Test Significant', 'Beta-Binomial Test Significant'])
#
#    # Save the plot
#    os.makedirs(outdir, exist_ok=True)  # Ensure output directory exists
#    plot_path = os.path.join(outdir, f"{outprefix}_allelic_imbalance_plot.png")
#    plt.savefig(plot_path, dpi=300)
#    plt.close()
#
#    print(f"Plot saved: {plot_path}")

# Function to generate and save the allelic imbalance plot
def plot_allelic_imbalance(df, outdir, outprefix, effect_allele=False):
    """
    Generate and save a scatter plot showing allelic imbalance.

    Args:
    - df (pd.DataFrame): DataFrame containing allele-specific read counts.
    - outdir (str): Directory where the output plot should be saved.
    - outprefix (str): Prefix for the output plot filename.
    - effect_allele (bool): Whether to use A1 vs A2 for the plot.

    Saves:
    - Scatter plot PNG file.
    """
    print("Generating allelic imbalance plot...")

    plt.figure(figsize=(10, 6))

    if effect_allele:
        # Ensure A1 and A2 exist
        if 'A1' not in df.columns or 'A2' not in df.columns:
            raise ValueError("Error: --effect-allele is enabled, but A1 and A2 columns are missing in the dataset.")

        # Extract numeric read counts for A1 and A2
        df['A1_count'] = df.apply(lambda row: row[row['A1']], axis=1)
        df['A2_count'] = df.apply(lambda row: row[row['A2']], axis=1)

        # Compute A1 allele frequency
        df['A1_freq'] = df['A1_count'] / (df['A1_count'] + df['A2_count'])

        # Plot A1 allele frequency
        plt.scatter(df['total_reads'], df['A1_freq'],
                    c=df['significant_binomial'].map({True: 'red', False: 'blue'}), alpha=0.6, label='Binomial Test')

        plt.ylabel('A1 Allele Frequency')

    else:
        # Default plot using max allele
        if 'max_allele_count' not in df.columns:
            raise ValueError("Error: 'max_allele_count' is missing. Ensure binomial test was run first.")

        plt.scatter(df['total_reads'], df['max_allele_count'] / df['total_reads'],
                    c=df['significant_binomial'].map({True: 'red', False: 'blue'}), alpha=0.6, label='Binomial Test')

        plt.ylabel('Major Allele Frequency')

    plt.axhline(0.5, linestyle='--', color='black', linewidth=1)
    plt.xlabel('Total Reads at SNP')
    plt.title('Allelic Imbalance Analysis')
    plt.legend(title='Significance', labels=['Non-significant', 'Binomial Test Significant'])

    os.makedirs(outdir, exist_ok=True)
    plot_path = os.path.join(outdir, f"{outprefix}_allelic_imbalance_plot.png")
    plt.savefig(plot_path, dpi=300)
    plt.close()

    print(f"Plot saved: {plot_path}")



# Main function to run the analysis
def main(bam_file, snp_file, outdir, outprefix, min_reads, no_depth_filter, effect_allele):
    print("\nStarting allelic imbalance analysis...")
    
    #allele_counts_df = count_allelic_reads(bam_file, snp_file, min_reads, apply_filter=not no_depth_filter)
    allele_counts_df = count_allelic_reads(bam_file, snp_file, min_reads, apply_filter=not no_depth_filter, effect_allele=effect_allele)
    print(allele_counts_df)

    results_df = perform_binomial_test(allele_counts_df, effect_allele)
    print(results_df)
    results_df = perform_beta_binomial_test(results_df, effect_allele)
    print(results_df)

    save_results(results_df, outdir, outprefix)
    plot_allelic_imbalance(results_df, outdir, outprefix, effect_allele)
    
    print("\nAnalysis complete.\n")


# Command-line argument parsing
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Allelic Imbalance Analysis for TF Binding")
    parser.add_argument("--bam", required=True, help="Input BAM file")
    parser.add_argument("--snp", required=True, help="Input SNP coordinate file")
    parser.add_argument("--outdir", required=True, help="Directory to save output files")
    parser.add_argument("--outprefix", required=True, help="Prefix for output file names")
    parser.add_argument("--min-reads", type=int, default=5, help="Minimum read depth filter (default: 5)")
    parser.add_argument("--no-depth-filter", action="store_true", help="Disable read depth filtering")
    parser.add_argument("--effect-allele", action="store_true", help="Use A1 as effect allele and A2 as reference allele")

    args = parser.parse_args()
    main(args.bam, args.snp, args.outdir, args.outprefix, args.min_reads, args.no_depth_filter, args.effect_allele)

