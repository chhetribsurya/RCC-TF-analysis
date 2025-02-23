#!/usr/bin/env python

import pandas as pd
import pybedtools
import seaborn as sns
import matplotlib.pyplot as plt
import argparse
import os

def split_chromatin_states(filepath, output_dir):
    """
    Splits the chromatin state data into separate BED files for each state.

    Args:
        filepath (str): Path to the chromatin state BED file.
        output_dir (str): Directory to save the split chromatin state BED files.
    """
    chromatin_data = pd.read_csv(filepath, sep='\t', header=None, usecols=[0, 1, 2, 3], names=['chr', 'start', 'end', 'state'])
    chromatin_data = chromatin_data.dropna(subset=['state'])
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    for state_id, group in chromatin_data.groupby('state'):
        safe_state_id = state_id.replace('/', '_')
        state_bed_file = os.path.join(output_dir, f'state_{safe_state_id}.bed')
        group[['chr', 'start', 'end']].astype({'start': 'int', 'end': 'int'}).to_csv(state_bed_file, sep='\t', header=False, index=False)

def read_tf_data(filepath, clean_dir):
    """
    Reads the transcription factor data from a BED file and retains only necessary columns.

    Args:
        filepath (str): Path to the transcription factor BED file.
        clean_dir (str): Directory to save the clean BED files.

    Returns:
        pybedtools.BedTool: BedTool object containing the transcription factor data.
    """
    if not os.path.exists(clean_dir):
        os.makedirs(clean_dir)
    
    tf_data = pd.read_csv(filepath, sep='\t', header=None, usecols=[0, 1, 2], names=['chr', 'start', 'end'])
    tf_bed_file = os.path.join(clean_dir, os.path.basename(filepath).replace('.bed', '_clean.bed'))
    tf_data.to_csv(tf_bed_file, sep='\t', header=False, index=False)
    return pybedtools.BedTool(tf_bed_file)

def quantify_tf_preference(tf_bed, chromatin_state_beds):
    """
    Quantifies the preference of transcription factors for different chromatin states.

    Args:
        tf_bed (pybedtools.BedTool): BedTool object containing the transcription factor data.
        chromatin_state_beds (dict): Dictionary of chromatin state BedTool objects.

    Returns:
        pd.Series: Series of fractions of TF binding sites associated with each chromatin state.
    """
    total_sites = len(tf_bed)
    state_proportions = {}
    
    for state, state_bed in chromatin_state_beds.items():
        if state == 'nan':
            continue
        intersection = tf_bed.intersect(state_bed, wa=True)
        if intersection.count() == 0:
            state_proportions[state] = 0
            continue
        unique_sites = intersection.to_dataframe().drop_duplicates(subset=['chrom', 'start', 'end'])
        state_proportions[state] = len(unique_sites) / total_sites
    
    return pd.Series(state_proportions)

def plot_heatmap(data, output_file, cmap, cluster_cols=True):
    """
    Plots a heatmap of TFs by chromatin state preferences.

    Args:
        data (pd.DataFrame): Matrix of TFs by chromatin state preferences.
        output_file (str): Path to save the heatmap.
        cmap (str): Color map for the heatmap.
        cluster_cols (bool): Whether to cluster the columns.
    """
    sns.clustermap(data, cmap=cmap, row_cluster=True, col_cluster=cluster_cols, linewidths=.5)
    plt.savefig(output_file)
    plt.close()

def main(args):
    """
    Main function to execute the script.

    Args:
        args (argparse.Namespace): Command line arguments.
    """
    chromatin_states_dir = os.path.join(args.output_dir, 'chromatin_states')
    split_chromatin_states(args.chromatin_file, chromatin_states_dir)
    
    chromatin_state_files = [os.path.join(chromatin_states_dir, f) for f in os.listdir(chromatin_states_dir) if f.endswith('.bed')]
    chromatin_state_beds = {os.path.basename(f).replace('state_', '').replace('.bed', ''): pybedtools.BedTool(f) for f in chromatin_state_files}
    
    clean_tf_dir = os.path.join(args.output_dir, 'clean_tf_files')
    tf_files = [os.path.join(args.tf_dir, f) for f in os.listdir(args.tf_dir) if f.endswith('.bed')]
    
    tf_preference_list = []
    max_preference_data = []

    for tf_file in tf_files:
        tf_bed = read_tf_data(tf_file, clean_tf_dir)
        tf_name = os.path.basename(tf_file).replace('.bed', '')
        print(f"Processing TF file: {tf_name}")
        tf_preference = quantify_tf_preference(tf_bed, chromatin_state_beds)
        tf_preference.name = tf_name
        tf_preference_list.append(tf_preference)
        max_state = tf_preference.idxmax()
        max_value = tf_preference.max()
        max_preference_data.append([tf_name, max_state, max_value])

    all_tf_preference = pd.concat(tf_preference_list, axis=1).fillna(0).T

    # Arrange columns in numeric alpha character order
    ordered_columns = sorted(all_tf_preference.columns, key=lambda x: (int(x.split('_')[0]), x))
    all_tf_preference = all_tf_preference[ordered_columns]

    # Save all_tf_preference data to file
    all_tf_preference.to_csv(os.path.join(args.output_dir, 'tf_chromatin_preference.csv'), index=True)
    print(f"Saved TF chromatin preference data to {os.path.join(args.output_dir, 'tf_chromatin_preference.csv')}")

    # Plot and save heatmap
    plot_heatmap(all_tf_preference, os.path.join(args.output_dir, 'tf_chromatin_preference_heatmap.png'), cmap='YlOrRd', cluster_cols=False)

    # Z-score normalization by columns
    zscored_data = all_tf_preference.apply(lambda x: (x - x.mean()) / x.std(), axis=0)

    # Plot and save Z-scored heatmap with diverging palette
    diverging_cmap = sns.diverging_palette(220, 20, as_cmap=True)
    #plot_heatmap(zscored_data, os.path.join(args.output_dir, 'tf_chromatin_preference_heatmap_zscored.png'), cmap=diverging_cmap, cluster_cols=False)

    # Save max preference data to file
    max_preference_df = pd.DataFrame(max_preference_data, columns=['TF', 'Chromatin State', 'Proportion'])
    max_preference_df.to_csv(os.path.join(args.output_dir, 'max_tf_chromatin_preference.csv'), index=False)
    print(f"Saved max TF chromatin preference data to {os.path.join(args.output_dir, 'max_tf_chromatin_preference.csv')}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Quantify TF preference for chromatin states and plot heatmap.")
    parser.add_argument('--chromatin_file', required=True, help="Path to the chromatin state BED file.")
    parser.add_argument('--tf_dir', required=True, help="Directory containing TF BED files.")
    parser.add_argument('--output_dir', required=True, help="Directory to save the split chromatin state BED files and output data.")
    parser.add_argument('--output_file_dir', required=True, help="Directory path to save the heatmap.")
    
    args = parser.parse_args()
    main(args)


