#!/usr/bin/env python3
"""
Differential Abundance Analysis for DIAMOND Resistome Data

Performs statistical comparisons of gene abundance between:
1. UCMC vs ZCH
2. Body sites (Axilla vs Groin vs Stool)
3. Week 1 vs Week 3

Uses Mann-Whitney U test (for 2 groups) and Kruskal-Wallis H test (for 3+ groups)
with FDR correction for multiple testing.
"""

import pandas as pd
import numpy as np
from pathlib import Path
from scipy import stats
from statsmodels.stats.multitest import multipletests
import warnings
warnings.filterwarnings('ignore')

# Define paths
PROJECT_ROOT = Path("/home/david/projects/benchmark_biogpu")
DATA_DIR = PROJECT_ROOT / "data"
RESULTS_DIR = PROJECT_ROOT / "results" / "analysis" / "differential_abundance"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# Analysis parameters
MIN_SAMPLES = 5  # Minimum number of samples a gene must be present in
MIN_RPM = 1.0    # Minimum RPM threshold for considering a gene "present"
FDR_THRESHOLD = 0.05  # FDR threshold for significance

def load_data():
    """Load gene matrix and metadata"""
    print("Loading data...")

    # Load RPM matrix
    rpm_matrix = pd.read_csv(RESULTS_DIR.parent / "diamond_resistome" / "amr_rpm_matrix.tsv",
                             sep="\t", index_col=0)

    # Load updated metadata (from corrected ZCH_UCMC metadata)
    metadata = pd.read_csv(DATA_DIR / "diamond_metadata_updated.tsv", sep="\t")
    metadata = metadata.set_index('sample_name')

    # Filter to samples that have metadata
    common_samples = rpm_matrix.index.intersection(metadata.index)
    rpm_matrix = rpm_matrix.loc[common_samples]
    metadata = metadata.loc[common_samples]

    print(f"Loaded {rpm_matrix.shape[0]} samples × {rpm_matrix.shape[1]} genes")
    print(f"  UCMC: {(metadata['Location']=='UCMC').sum()}")
    print(f"  ZCH: {(metadata['Location']=='ZCH').sum()}")

    return rpm_matrix, metadata

def filter_genes(rpm_matrix, min_samples=MIN_SAMPLES, min_rpm=MIN_RPM):
    """Filter genes present in at least min_samples samples"""
    print(f"\nFiltering genes...")
    print(f"  Keeping genes present (RPM > {min_rpm}) in ≥ {min_samples} samples")

    # Count samples where gene RPM > threshold
    presence_counts = (rpm_matrix > min_rpm).sum(axis=0)
    genes_to_keep = presence_counts[presence_counts >= min_samples].index

    filtered_matrix = rpm_matrix[genes_to_keep]

    print(f"  Genes before filtering: {rpm_matrix.shape[1]}")
    print(f"  Genes after filtering: {filtered_matrix.shape[1]}")
    print(f"  Genes removed: {rpm_matrix.shape[1] - filtered_matrix.shape[1]}")

    return filtered_matrix

def mann_whitney_test(rpm_matrix, metadata, group_col, group1, group2):
    """
    Perform Mann-Whitney U test for each gene between two groups
    """
    print(f"\nTesting {group1} vs {group2}...")

    # Get sample indices for each group
    group1_samples = metadata[metadata[group_col] == group1].index
    group2_samples = metadata[metadata[group_col] == group2].index

    print(f"  {group1}: n={len(group1_samples)}")
    print(f"  {group2}: n={len(group2_samples)}")

    results = []

    for gene in rpm_matrix.columns:
        group1_values = rpm_matrix.loc[group1_samples, gene]
        group2_values = rpm_matrix.loc[group2_samples, gene]

        # Calculate statistics
        mean1 = group1_values.mean()
        mean2 = group2_values.mean()
        median1 = group1_values.median()
        median2 = group2_values.median()

        # Presence/prevalence in each group
        prev1 = (group1_values > MIN_RPM).sum()
        prev2 = (group2_values > MIN_RPM).sum()

        # Mann-Whitney U test
        try:
            stat, pval = stats.mannwhitneyu(group1_values, group2_values, alternative='two-sided')
        except:
            pval = 1.0
            stat = np.nan

        # Log2 fold change (using means with pseudocount)
        log2fc = np.log2((mean2 + 0.01) / (mean1 + 0.01))

        results.append({
            'gene': gene,
            f'{group1}_mean': mean1,
            f'{group2}_mean': mean2,
            f'{group1}_median': median1,
            f'{group2}_median': median2,
            f'{group1}_prevalence': prev1,
            f'{group2}_prevalence': prev2,
            'log2_fold_change': log2fc,
            'statistic': stat,
            'p_value': pval
        })

    results_df = pd.DataFrame(results)

    # FDR correction
    results_df['fdr'] = multipletests(results_df['p_value'], method='fdr_bh')[1]

    # Sort by p-value
    results_df = results_df.sort_values('p_value')

    # Count significant genes
    n_sig = (results_df['fdr'] < FDR_THRESHOLD).sum()
    print(f"  Significant genes (FDR < {FDR_THRESHOLD}): {n_sig}")

    return results_df

def kruskal_wallis_test(rpm_matrix, metadata, group_col):
    """
    Perform Kruskal-Wallis H test for each gene across multiple groups
    """
    groups = metadata[group_col].unique()
    print(f"\nTesting across {len(groups)} groups: {', '.join(groups)}...")

    # Print sample sizes
    for group in groups:
        n = (metadata[group_col] == group).sum()
        print(f"  {group}: n={n}")

    results = []

    for gene in rpm_matrix.columns:
        group_values = [rpm_matrix.loc[metadata[metadata[group_col] == g].index, gene]
                       for g in groups]

        # Calculate statistics for each group
        group_stats = {}
        for g, vals in zip(groups, group_values):
            group_stats[f'{g}_mean'] = vals.mean()
            group_stats[f'{g}_median'] = vals.median()
            group_stats[f'{g}_prevalence'] = (vals > MIN_RPM).sum()

        # Kruskal-Wallis H test
        try:
            stat, pval = stats.kruskal(*group_values)
        except:
            pval = 1.0
            stat = np.nan

        result = {
            'gene': gene,
            'statistic': stat,
            'p_value': pval
        }
        result.update(group_stats)
        results.append(result)

    results_df = pd.DataFrame(results)

    # FDR correction
    results_df['fdr'] = multipletests(results_df['p_value'], method='fdr_bh')[1]

    # Sort by p-value
    results_df = results_df.sort_values('p_value')

    # Count significant genes
    n_sig = (results_df['fdr'] < FDR_THRESHOLD).sum()
    print(f"  Significant genes (FDR < {FDR_THRESHOLD}): {n_sig}")

    return results_df

def main():
    print("="*60)
    print("DIFFERENTIAL ABUNDANCE ANALYSIS")
    print("="*60)

    # Load data
    rpm_matrix, metadata = load_data()

    # Filter genes
    rpm_matrix_filtered = filter_genes(rpm_matrix)

    # Analysis 1: UCMC vs ZCH
    print("\n" + "="*60)
    print("1. LOCATION: UCMC vs ZCH")
    print("="*60)
    location_results = mann_whitney_test(
        rpm_matrix_filtered, metadata,
        'Location', 'UCMC', 'ZCH'
    )
    output_file = RESULTS_DIR / "location_ucmc_vs_zch.tsv"
    location_results.to_csv(output_file, sep="\t", index=False)
    print(f"✓ Saved: {output_file}")

    # Analysis 2: Body sites
    print("\n" + "="*60)
    print("2. BODY SITE: Axilla vs Groin vs Stool")
    print("="*60)
    bodysite_results = kruskal_wallis_test(
        rpm_matrix_filtered, metadata,
        'SampleType'
    )
    output_file = RESULTS_DIR / "bodysite_comparison.tsv"
    bodysite_results.to_csv(output_file, sep="\t", index=False)
    print(f"✓ Saved: {output_file}")

    # Pairwise comparisons for body sites
    print("\n  Pairwise comparisons:")

    for site1, site2 in [('Axilla', 'Groin'), ('Axilla', 'Stool'), ('Groin', 'Stool')]:
        results = mann_whitney_test(
            rpm_matrix_filtered, metadata,
            'SampleType', site1, site2
        )
        output_file = RESULTS_DIR / f"bodysite_{site1.lower()}_vs_{site2.lower()}.tsv"
        results.to_csv(output_file, sep="\t", index=False)
        print(f"  ✓ Saved: {output_file}")

    # Analysis 3: Week 1 vs Week 3
    print("\n" + "="*60)
    print("3. POSTNATAL AGE: Week 1 vs Week 3")
    print("="*60)
    week_results = mann_whitney_test(
        rpm_matrix_filtered, metadata,
        'SampleCollectionWeek', 'Week.1', 'Week.3'
    )
    output_file = RESULTS_DIR / "week_week1_vs_week3.tsv"
    week_results.to_csv(output_file, sep="\t", index=False)
    print(f"✓ Saved: {output_file}")

    # Summary
    print("\n" + "="*60)
    print("DIFFERENTIAL ABUNDANCE ANALYSIS COMPLETE")
    print("="*60)
    print(f"\nResults saved to: {RESULTS_DIR}")
    print("\nComparisons performed:")
    print("  1. Location: UCMC vs ZCH")
    print("  2. Body Site: Axilla vs Groin vs Stool (+ pairwise)")
    print("  3. Postnatal Age: Week 1 vs Week 3")
    print(f"\nSignificance threshold: FDR < {FDR_THRESHOLD}")
    print(f"Gene filtering: Present in ≥ {MIN_SAMPLES} samples (RPM > {MIN_RPM})")

if __name__ == "__main__":
    main()
