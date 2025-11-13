#!/usr/bin/env python3
"""
Gene × Antibiotic Correlation Analysis

Correlates individual gene abundance (RPM) with individual antibiotic exposures
for all genes and all antibiotics.

Key insight: Antibiotic exposure is CUMULATIVE
- Week 1 samples: exposure = w1 antibiotics only
- Week 3 samples: exposure = w1 + w2 antibiotics (cumulative)

Output: Gene-by-antibiotic correlation table filtered to p < 0.05 (uncorrected)
"""

import pandas as pd
import numpy as np
from pathlib import Path
from scipy.stats import spearmanr
import warnings
warnings.filterwarnings('ignore')

# Define paths
PROJECT_ROOT = Path("/home/david/projects/benchmark_biogpu")
DATA_DIR = PROJECT_ROOT / "data"
RESULTS_DIR = PROJECT_ROOT / "results" / "analysis" / "gene_antibiotic_correlations"

# Create directories
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# Significance threshold
P_THRESHOLD = 0.05

def load_data():
    """Load RPM matrix and antibiotic exposure data"""
    print("Loading data...")

    # Load RPM matrix (samples × genes)
    rpm_matrix = pd.read_csv(
        PROJECT_ROOT / "results" / "analysis" / "diamond_resistome" / "amr_rpm_matrix.tsv",
        sep="\t", index_col=0
    )
    print(f"  RPM matrix: {rpm_matrix.shape[0]} samples × {rpm_matrix.shape[1]} genes")

    # Load metadata with antibiotic exposure (fixed metadata)
    metadata = pd.read_csv(DATA_DIR / "master_metadata_fixed.tsv", sep="\t")
    print(f"  Metadata: {len(metadata)} samples")

    # Get antibiotic columns
    abx_w1_cols = [col for col in metadata.columns if col.endswith('_w1') and col != 'None_w1']
    abx_w2_cols = [col for col in metadata.columns if col.endswith('_w2') and col != 'None_w2']

    # Map w2 columns to w1 columns (same antibiotic)
    abx_pairs = {}
    for w1_col in abx_w1_cols:
        abx_name = w1_col.replace('_w1', '')
        w2_col = abx_name + '_w2'
        if w2_col in abx_w2_cols:
            abx_pairs[abx_name] = (w1_col, w2_col)

    print(f"  Found {len(abx_pairs)} antibiotics with both w1 and w2 data")

    # Keep only sample_name, Location, SampleCollectionWeek, and antibiotic columns
    keep_cols = ['sample_name', 'Location', 'SampleCollectionWeek', 'SampleType'] + abx_w1_cols + abx_w2_cols
    keep_cols = [col for col in keep_cols if col in metadata.columns]
    metadata = metadata[keep_cols].copy()

    # Fill NaN with 0 for antibiotic columns
    for col in abx_w1_cols + abx_w2_cols:
        if col in metadata.columns:
            metadata[col] = metadata[col].fillna(0)

    return rpm_matrix, metadata, abx_pairs

def calculate_cumulative_exposure(metadata, abx_pairs):
    """
    Calculate cumulative antibiotic exposure
    Week 1 samples: exposure = w1 only
    Week 3 samples: exposure = w1 + w2 (cumulative)
    """
    print("\nCalculating cumulative antibiotic exposure...")

    # Calculate cumulative exposure for each antibiotic
    for abx_name, (w1_col, w2_col) in abx_pairs.items():
        cumulative_col = f'{abx_name}_cumulative'

        # Week 1 samples: only w1 exposure
        # Week 3 samples: w1 + w2 (cumulative)
        metadata[cumulative_col] = metadata.apply(
            lambda row: row[w1_col] if row['SampleCollectionWeek'] == 'Week.1'
                        else row[w1_col] + row[w2_col] if row['SampleCollectionWeek'] == 'Week.3'
                        else 0,
            axis=1
        )

    cumulative_cols = [f'{abx}_cumulative' for abx in abx_pairs.keys()]
    print(f"  Created {len(cumulative_cols)} cumulative exposure variables")

    return metadata, cumulative_cols

def merge_data(rpm_matrix, metadata):
    """Merge RPM matrix with antibiotic exposure data"""
    print("\nMerging RPM matrix with antibiotic data...")

    # Fix ZCH sample names (add ZJH_ prefix)
    metadata_fixed = metadata.copy()
    zch_mask = metadata_fixed['Location'] == 'ZCH'
    metadata_fixed.loc[zch_mask, 'sample_name'] = 'ZJH_' + metadata_fixed.loc[zch_mask, 'sample_name']

    # Get common samples
    common_samples = rpm_matrix.index.intersection(metadata_fixed['sample_name'])
    print(f"  Samples in RPM matrix: {len(rpm_matrix)}")
    print(f"  Samples in metadata: {len(metadata_fixed)}")
    print(f"  Common samples: {len(common_samples)}")

    # Filter to common samples
    rpm_matrix_filtered = rpm_matrix.loc[common_samples]
    metadata_filtered = metadata_fixed[metadata_fixed['sample_name'].isin(common_samples)].copy()
    metadata_filtered = metadata_filtered.set_index('sample_name').loc[common_samples]

    return rpm_matrix_filtered, metadata_filtered

def correlate_genes_with_antibiotics(rpm_matrix, metadata, cumulative_cols):
    """
    Correlate each gene with each antibiotic
    Returns table with significant correlations (p < 0.05)
    """
    print("\n" + "="*60)
    print("GENE × ANTIBIOTIC CORRELATION ANALYSIS")
    print("="*60)
    print(f"\nAnalyzing {rpm_matrix.shape[1]} genes × {len(cumulative_cols)} antibiotics")
    print(f"Total comparisons: {rpm_matrix.shape[1] * len(cumulative_cols):,}")

    results = []
    total_comparisons = len(cumulative_cols) * rpm_matrix.shape[1]
    completed = 0

    for abx_col in cumulative_cols:
        abx_name = abx_col.replace('_cumulative', '')

        # Get antibiotic exposure vector
        abx_exposure = metadata[abx_col]

        # Skip if no variation
        if abx_exposure.std() == 0:
            print(f"  Skipping {abx_name} - no variation in exposure")
            completed += rpm_matrix.shape[1]
            continue

        # Number of exposed samples
        n_exposed = (abx_exposure > 0).sum()

        if n_exposed < 5:  # Need minimum samples
            print(f"  Skipping {abx_name} - only {n_exposed} exposed samples")
            completed += rpm_matrix.shape[1]
            continue

        # Correlate with each gene
        for gene in rpm_matrix.columns:
            gene_abundance = rpm_matrix[gene]

            # Skip if no variation
            if gene_abundance.std() == 0:
                completed += 1
                continue

            # Spearman correlation
            try:
                rho, pval = spearmanr(abx_exposure, gene_abundance)
            except:
                completed += 1
                continue

            # Only keep if p < threshold
            if pval < P_THRESHOLD:
                results.append({
                    'antibiotic': abx_name,
                    'gene': gene,
                    'n_exposed': n_exposed,
                    'n_samples': len(abx_exposure),
                    'spearman_rho': rho,
                    'pvalue': pval
                })

            completed += 1

        # Progress update
        if completed % 10000 == 0:
            print(f"  Completed {completed:,} / {total_comparisons:,} comparisons ({completed/total_comparisons*100:.1f}%)")

    print(f"\n✓ Completed all {total_comparisons:,} comparisons")

    # Create results dataframe
    results_df = pd.DataFrame(results)

    if len(results_df) > 0:
        # Sort by absolute correlation strength
        results_df['abs_rho'] = results_df['spearman_rho'].abs()
        results_df = results_df.sort_values('abs_rho', ascending=False)
        results_df = results_df.drop(columns=['abs_rho'])

        print(f"\nSignificant correlations (p < {P_THRESHOLD}): {len(results_df):,}")
        print(f"  Positive correlations (rho > 0): {(results_df['spearman_rho'] > 0).sum():,}")
        print(f"  Negative correlations (rho < 0): {(results_df['spearman_rho'] < 0).sum():,}")

    return results_df

def main():
    print("="*60)
    print("GENE × ANTIBIOTIC CORRELATION ANALYSIS")
    print("="*60)
    print("\nAnalyzing correlations between individual genes and antibiotics")
    print("Using cumulative antibiotic exposure (Week 1 = w1, Week 3 = w1 + w2)")
    print(f"Filtering to p < {P_THRESHOLD} (uncorrected)")

    # Load data
    rpm_matrix, metadata, abx_pairs = load_data()

    # Calculate cumulative exposure
    metadata, cumulative_cols = calculate_cumulative_exposure(metadata, abx_pairs)

    # Merge data
    rpm_matrix, metadata = merge_data(rpm_matrix, metadata)

    # Correlate genes with antibiotics
    results_df = correlate_genes_with_antibiotics(rpm_matrix, metadata, cumulative_cols)

    # Save results
    if len(results_df) > 0:
        output_file = RESULTS_DIR / f"gene_antibiotic_correlations_p{P_THRESHOLD}.tsv"
        results_df.to_csv(output_file, sep="\t", index=False)
        print(f"\n✓ Saved: {output_file}")

        # Summary statistics
        print("\n" + "="*60)
        print("SUMMARY")
        print("="*60)

        # Top positive correlations
        print("\nTop 10 positive correlations (antibiotic → higher gene abundance):")
        top_pos = results_df[results_df['spearman_rho'] > 0].head(10)
        for i, row in enumerate(top_pos.iterrows(), 1):
            _, r = row
            print(f"  {i}. {r['antibiotic']} × {r['gene']}: rho = {r['spearman_rho']:.3f}, p = {r['pvalue']:.2e}")

        # Top negative correlations
        print("\nTop 10 negative correlations (antibiotic → lower gene abundance):")
        top_neg = results_df[results_df['spearman_rho'] < 0].head(10)
        for i, row in enumerate(top_neg.iterrows(), 1):
            _, r = row
            print(f"  {i}. {r['antibiotic']} × {r['gene']}: rho = {r['spearman_rho']:.3f}, p = {r['pvalue']:.2e}")

        # By antibiotic
        print("\nSignificant correlations by antibiotic:")
        abx_summary = results_df.groupby('antibiotic').size().sort_values(ascending=False)
        for abx, count in abx_summary.head(10).items():
            print(f"  {abx}: {count} genes")

    else:
        print("\nNo significant correlations found!")

    print("\n" + "="*60)
    print("ANALYSIS COMPLETE")
    print("="*60)

if __name__ == "__main__":
    main()
