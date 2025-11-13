#!/usr/bin/env python3
"""
Antibiotic Exposure - AMR Correlation Analysis for DIAMOND Data

Key insight: Antibiotic exposure is CUMULATIVE
- Week 1 samples: exposure = w1 antibiotics only
- Week 3 samples: exposure = w1 + w2 antibiotics (cumulative)
  Note: w2 columns represent cumulative week 3 exposure in the antibiotic table

Analyses:
1. Correlate total antibiotic days with total AMR burden
2. Correlate specific antibiotics with relevant AMR genes
3. Stratify by body site
4. Test for dose-response relationships

Adapted from ZCH_UCMC_Manuscript/nicu_resistome_analysis for DIAMOND alignment results
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy.stats import spearmanr, mannwhitneyu
from statsmodels.stats.multitest import multipletests
import warnings
warnings.filterwarnings('ignore')

# Define paths
PROJECT_ROOT = Path("/home/david/projects/benchmark_biogpu")
DATA_DIR = PROJECT_ROOT / "data"
RESULTS_DIR = PROJECT_ROOT / "results" / "analysis" / "antibiotic_correlations"
FIGURES_DIR = PROJECT_ROOT / "figures" / "antibiotic_correlations"

# ZCH_UCMC metadata with antibiotic exposure data
ZCH_UCMC_ROOT = Path("/home/david/projects/ZCH_UCMC_Manuscript")
# Use fixed metadata with all sample_name fields populated
ZCH_UCMC_METADATA = PROJECT_ROOT / "data" / "master_metadata_fixed.tsv"

# Create directories
RESULTS_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

# Set plotting style
sns.set_style("whitegrid")
sns.set_context("talk")

def load_data():
    """Load DIAMOND data and merge with antibiotic exposure data"""
    print("Loading data...")

    # Load DIAMOND metadata
    diamond_metadata = pd.read_csv(DATA_DIR / "diamond_metadata.tsv", sep="\t")
    print(f"DIAMOND metadata: {len(diamond_metadata)} samples")

    # Load DIAMOND RPM matrix
    rpm_matrix = pd.read_csv(
        PROJECT_ROOT / "results" / "analysis" / "diamond_resistome" / "amr_rpm_matrix.tsv",
        sep="\t", index_col=0
    )
    print(f"RPM matrix: {rpm_matrix.shape[0]} samples × {rpm_matrix.shape[1]} genes")

    # Calculate total AMR RPM for each sample
    total_amr_rpm = rpm_matrix.sum(axis=1).reset_index()
    total_amr_rpm.columns = ['sample_name', 'total_amr_rpm']

    # Merge with DIAMOND metadata
    diamond_metadata = diamond_metadata.merge(total_amr_rpm, on='sample_name', how='left')

    # Load ZCH_UCMC metadata with antibiotic exposure data
    zch_ucmc_metadata = pd.read_csv(ZCH_UCMC_METADATA, sep="\t")
    print(f"ZCH_UCMC metadata: {len(zch_ucmc_metadata)} samples")

    # Get antibiotic columns
    abx_cols = [col for col in zch_ucmc_metadata.columns if
                (col.endswith('_w1') or col.endswith('_w2')) and
                col not in ['None_w1', 'None_w2', 'Other_w3']]

    # Keep only sample_name, Location and antibiotic columns
    abx_data = zch_ucmc_metadata[['sample_name', 'Location'] + abx_cols].copy()

    # Fix naming: ZCH samples in metadata don't have ZJH_ prefix
    # but DIAMOND data does, so add it for ZCH samples
    abx_data['sample_name_original'] = abx_data['sample_name']
    abx_data.loc[abx_data['Location'] == 'ZCH', 'sample_name'] = 'ZJH_' + abx_data.loc[abx_data['Location'] == 'ZCH', 'sample_name']

    # Drop Location column (already in diamond_metadata)
    abx_data = abx_data.drop(columns=['Location'])

    print(f"  UCMC samples in abx data: {(abx_data['sample_name_original'].str.startswith('N')).sum()}")
    print(f"  ZCH samples in abx data (with ZJH_ prefix added): {(abx_data['sample_name'].str.startswith('ZJH_')).sum()}")

    # Merge antibiotic data with DIAMOND metadata
    metadata = diamond_metadata.merge(abx_data.drop(columns=['sample_name_original']), on='sample_name', how='left')

    # Flag samples with antibiotic data
    metadata['has_abx_data'] = metadata[abx_cols].notna().any(axis=1)

    print(f"\nMerged data:")
    print(f"  Total samples: {len(metadata)}")
    print(f"  Samples with antibiotic data: {metadata['has_abx_data'].sum()}")
    print(f"  By location:")
    print(f"    UCMC with abx data: {((metadata['Location']=='UCMC') & metadata['has_abx_data']).sum()}")
    print(f"    ZCH with abx data: {((metadata['Location']=='ZCH') & metadata['has_abx_data']).sum()}")

    return metadata

def calculate_cumulative_exposure(metadata):
    """
    Calculate cumulative antibiotic exposure
    Week 1 samples: exposure = w1 only
    Week 3 samples: exposure = w1 + w2 (cumulative)
    """
    print("\nCalculating cumulative antibiotic exposure...")

    # Get all antibiotic columns (exclude Subject, Location, None, Other)
    abx_w1_cols = [col for col in metadata.columns if col.endswith('_w1') and col != 'None_w1']
    abx_w2_cols = [col for col in metadata.columns if col.endswith('_w2') and col != 'None_w2']

    # Map w2 columns to w1 columns (same antibiotic)
    abx_pairs = {}
    for w1_col in abx_w1_cols:
        abx_name = w1_col.replace('_w1', '')
        w2_col = abx_name + '_w2'
        if w2_col in abx_w2_cols:
            abx_pairs[abx_name] = (w1_col, w2_col)

    print(f"Found {len(abx_pairs)} antibiotics with data for both weeks")

    # Fill NaN with 0 for antibiotic columns
    for col in abx_w1_cols + abx_w2_cols:
        metadata[col] = metadata[col].fillna(0)

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

    # Total cumulative antibiotic days
    cumulative_cols = [col for col in metadata.columns if col.endswith('_cumulative')]
    metadata['total_abx_cumulative'] = metadata[cumulative_cols].sum(axis=1)

    # Binary: any antibiotic exposure
    metadata['any_abx_cumulative'] = (metadata['total_abx_cumulative'] > 0).astype(int)

    print(f"Created cumulative exposure variables for {len(cumulative_cols)} antibiotics")
    print(f"\nExposure summary:")
    print(f"  Samples with any antibiotic: {metadata['any_abx_cumulative'].sum()}")
    print(f"  Samples without antibiotics: {(metadata['any_abx_cumulative']==0).sum()}")

    return metadata, cumulative_cols, list(abx_pairs.keys())

def correlate_total_abx_with_amr(metadata):
    """Correlate total antibiotic exposure with total AMR burden"""
    print("\n" + "="*60)
    print("TOTAL ANTIBIOTIC EXPOSURE vs TOTAL AMR BURDEN")
    print("="*60)

    # Filter to samples with antibiotic data
    data = metadata[metadata['has_abx_data']].copy()

    # Overall correlation
    rho, pval = spearmanr(data['total_abx_cumulative'], data['total_amr_rpm'])
    print(f"\nOverall (all samples):")
    print(f"  Spearman rho = {rho:.3f}, p = {pval:.4f}")

    # By location
    print(f"\nBy Location:")
    for loc in ['UCMC', 'ZCH']:
        loc_data = data[data['Location'] == loc]
        if len(loc_data) > 10:
            rho, pval = spearmanr(loc_data['total_abx_cumulative'], loc_data['total_amr_rpm'])
            print(f"  {loc}: rho = {rho:.3f}, p = {pval:.4f} (n={len(loc_data)})")

    # By body site
    print(f"\nBy Body Site:")
    results_by_site = []
    for site in ['Axilla', 'Groin', 'Stool']:
        site_data = data[data['SampleType'] == site]
        if len(site_data) > 10:
            rho, pval = spearmanr(site_data['total_abx_cumulative'], site_data['total_amr_rpm'])
            print(f"  {site}: rho = {rho:.3f}, p = {pval:.4f} (n={len(site_data)})")
            results_by_site.append({
                'body_site': site,
                'n_samples': len(site_data),
                'spearman_rho': rho,
                'pvalue': pval
            })

    # By week
    print(f"\nBy Week:")
    for week in ['Week.1', 'Week.3']:
        week_data = data[data['SampleCollectionWeek'] == week]
        if len(week_data) > 10:
            rho, pval = spearmanr(week_data['total_abx_cumulative'], week_data['total_amr_rpm'])
            print(f"  {week}: rho = {rho:.3f}, p = {pval:.4f} (n={len(week_data)})")

    # Save results
    results_df = pd.DataFrame(results_by_site)
    results_df.to_csv(RESULTS_DIR / "total_abx_amr_correlation_by_site.tsv", sep="\t", index=False)

    return data

def correlate_specific_antibiotics(metadata, abx_names):
    """Correlate specific antibiotics with AMR burden"""
    print("\n" + "="*60)
    print("SPECIFIC ANTIBIOTIC CORRELATIONS")
    print("="*60)

    data = metadata[metadata['has_abx_data']].copy()

    results = []

    for abx in abx_names:
        cumulative_col = f'{abx}_cumulative'

        if cumulative_col not in data.columns:
            continue

        # Filter to subjects who received this antibiotic
        exposed = data[data[cumulative_col] > 0]

        if len(exposed) < 5:  # Need minimum samples
            continue

        # Overall correlation
        rho, pval = spearmanr(data[cumulative_col], data['total_amr_rpm'])

        # Compare exposed vs unexposed
        unexposed = data[data[cumulative_col] == 0]
        if len(unexposed) > 0:
            stat, mw_pval = mannwhitneyu(
                exposed['total_amr_rpm'],
                unexposed['total_amr_rpm'],
                alternative='two-sided'
            )
        else:
            mw_pval = 1.0

        results.append({
            'antibiotic': abx,
            'n_exposed': len(exposed),
            'n_unexposed': len(unexposed),
            'mean_amr_exposed': exposed['total_amr_rpm'].mean(),
            'mean_amr_unexposed': unexposed['total_amr_rpm'].mean() if len(unexposed) > 0 else np.nan,
            'spearman_rho': rho,
            'spearman_pvalue': pval,
            'mannwhitney_pvalue': mw_pval
        })

    # Create results dataframe
    results_df = pd.DataFrame(results)

    if len(results_df) > 0:
        # FDR correction
        results_df['spearman_fdr'] = multipletests(results_df['spearman_pvalue'], method='fdr_bh')[1]
        results_df['mannwhitney_fdr'] = multipletests(results_df['mannwhitney_pvalue'], method='fdr_bh')[1]

        # Sort by correlation strength
        results_df = results_df.sort_values('spearman_rho', ascending=False)

        # Save
        results_df.to_csv(RESULTS_DIR / "specific_antibiotic_correlations.tsv", sep="\t", index=False)
        print(f"\n✓ Saved correlations for {len(results_df)} antibiotics")

        # Show top correlations
        print(f"\nTop positive correlations (antibiotics → higher AMR):")
        top_pos = results_df[results_df['spearman_rho'] > 0].head(5)
        for _, row in top_pos.iterrows():
            print(f"  {row['antibiotic']}: rho = {row['spearman_rho']:.3f}, p = {row['spearman_pvalue']:.4f}, n = {row['n_exposed']}")

        print(f"\nTop negative correlations (antibiotics → lower AMR):")
        top_neg = results_df[results_df['spearman_rho'] < 0].head(5)
        for _, row in top_neg.iterrows():
            print(f"  {row['antibiotic']}: rho = {row['spearman_rho']:.3f}, p = {row['spearman_pvalue']:.4f}, n = {row['n_exposed']}")

    return results_df

def plot_abx_amr_correlation(data):
    """Plot total antibiotic exposure vs AMR burden"""
    print("\nCreating correlation plots...")

    # Figure 1: Overall scatter plot
    fig, axes = plt.subplots(2, 2, figsize=(16, 14))

    # Overall
    ax = axes[0, 0]
    ax.scatter(data['total_abx_cumulative'], data['total_amr_rpm'],
               alpha=0.5, s=50, edgecolors='black', linewidths=0.5)
    rho, pval = spearmanr(data['total_abx_cumulative'], data['total_amr_rpm'])
    ax.set_xlabel('Total Cumulative Antibiotic Days')
    ax.set_ylabel('Total AMR RPM')
    ax.set_title(f'All Samples (n={len(data)})\nSpearman rho = {rho:.3f}, p = {pval:.4f}')
    ax.grid(True, alpha=0.3)

    # By location
    ax = axes[0, 1]
    colors = {'UCMC': '#1f77b4', 'ZCH': '#ff7f0e'}
    for loc in ['UCMC', 'ZCH']:
        loc_data = data[data['Location'] == loc]
        if len(loc_data) > 0:
            ax.scatter(loc_data['total_abx_cumulative'], loc_data['total_amr_rpm'],
                       c=colors[loc], label=loc, alpha=0.6, s=50, edgecolors='black', linewidths=0.5)
    ax.set_xlabel('Total Cumulative Antibiotic Days')
    ax.set_ylabel('Total AMR RPM')
    ax.set_title('By Location')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # By body site
    ax = axes[1, 0]
    site_colors = {'Axilla': '#1f77b4', 'Groin': '#ff7f0e', 'Stool': '#2ca02c'}
    for site in ['Axilla', 'Groin', 'Stool']:
        site_data = data[data['SampleType'] == site]
        if len(site_data) > 0:
            ax.scatter(site_data['total_abx_cumulative'], site_data['total_amr_rpm'],
                       c=site_colors[site], label=site, alpha=0.6, s=50, edgecolors='black', linewidths=0.5)
    ax.set_xlabel('Total Cumulative Antibiotic Days')
    ax.set_ylabel('Total AMR RPM')
    ax.set_title('By Body Site')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # By week
    ax = axes[1, 1]
    week_colors = {'Week.1': '#2ca02c', 'Week.3': '#d62728'}
    week_markers = {'Week.1': 'o', 'Week.3': 's'}
    for week in ['Week.1', 'Week.3']:
        week_data = data[data['SampleCollectionWeek'] == week]
        if len(week_data) > 0:
            ax.scatter(week_data['total_abx_cumulative'], week_data['total_amr_rpm'],
                       c=week_colors[week], marker=week_markers[week],
                       label=week, alpha=0.6, s=50, edgecolors='black', linewidths=0.5)
    ax.set_xlabel('Total Cumulative Antibiotic Days')
    ax.set_ylabel('Total AMR RPM')
    ax.set_title('By Week')
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "antibiotic_amr_correlation.pdf", dpi=300, bbox_inches='tight')
    plt.savefig(FIGURES_DIR / "antibiotic_amr_correlation.png", dpi=300, bbox_inches='tight')
    print(f"✓ Saved: {FIGURES_DIR / 'antibiotic_amr_correlation.pdf'}")
    plt.close()

    # Figure 2: Boxplot - Any antibiotics vs no antibiotics
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))

    for i, site in enumerate(['Axilla', 'Groin', 'Stool']):
        ax = axes[i]
        site_data = data[data['SampleType'] == site]

        # Group by any antibiotic exposure
        no_abx = site_data[site_data['any_abx_cumulative'] == 0]['total_amr_rpm']
        yes_abx = site_data[site_data['any_abx_cumulative'] == 1]['total_amr_rpm']

        if len(no_abx) > 0 or len(yes_abx) > 0:
            ax.boxplot([no_abx, yes_abx], labels=['No Antibiotics', 'Any Antibiotics'])
            ax.set_ylabel('Total AMR RPM')
            ax.set_title(f'{site}\n(No abx: n={len(no_abx)}, Abx: n={len(yes_abx)})')

            # Mann-Whitney test
            if len(no_abx) > 0 and len(yes_abx) > 0:
                stat, pval = mannwhitneyu(no_abx, yes_abx, alternative='two-sided')
                ax.text(0.5, 0.95, f'p = {pval:.4f}',
                        transform=ax.transAxes, ha='center', va='top', fontsize=10)

            ax.grid(True, alpha=0.3, axis='y')

    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "antibiotic_exposure_amr_boxplot.pdf", dpi=300, bbox_inches='tight')
    plt.savefig(FIGURES_DIR / "antibiotic_exposure_amr_boxplot.png", dpi=300, bbox_inches='tight')
    print(f"✓ Saved: {FIGURES_DIR / 'antibiotic_exposure_amr_boxplot.pdf'}")
    plt.close()

def main():
    print("="*60)
    print("ANTIBIOTIC-AMR CORRELATION ANALYSIS")
    print("DIAMOND Data")
    print("="*60)
    print("\nNote: Using CUMULATIVE antibiotic exposure")
    print("  - Week 1 samples: w1 antibiotics only")
    print("  - Week 3 samples: w1 + w2 antibiotics (cumulative)")

    # Load data
    metadata = load_data()

    # Calculate cumulative exposure
    metadata, cumulative_cols, abx_names = calculate_cumulative_exposure(metadata)

    # Correlate total antibiotics with total AMR
    data_with_abx = correlate_total_abx_with_amr(metadata)

    # Correlate specific antibiotics
    specific_results = correlate_specific_antibiotics(metadata, abx_names)

    # Plot
    plot_abx_amr_correlation(data_with_abx)

    print("\n" + "="*60)
    print("ANTIBIOTIC CORRELATION ANALYSIS COMPLETE")
    print("="*60)
    print("\nOutputs:")
    print(f"  - {RESULTS_DIR / 'total_abx_amr_correlation_by_site.tsv'}")
    print(f"  - {RESULTS_DIR / 'specific_antibiotic_correlations.tsv'}")
    print("\nFigures:")
    print(f"  - {FIGURES_DIR / 'antibiotic_amr_correlation.pdf'}")
    print(f"  - {FIGURES_DIR / 'antibiotic_exposure_amr_boxplot.pdf'}")

if __name__ == "__main__":
    main()
