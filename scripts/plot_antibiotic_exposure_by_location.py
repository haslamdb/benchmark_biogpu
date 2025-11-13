#!/usr/bin/env python3
"""
Antibiotic Exposure Comparison: ZCH vs UCMC

Creates boxplots comparing antibiotic exposure between locations:
- Week 1: w1 exposure only
- Week 3: cumulative (w1 + w2) exposure

One plot per antibiotic, faceted by week
Statistical tests: Mann-Whitney U test for each comparison
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy.stats import mannwhitneyu
import warnings
warnings.filterwarnings('ignore')

# Define paths
PROJECT_ROOT = Path("/home/david/projects/benchmark_biogpu")
DATA_DIR = PROJECT_ROOT / "data"
RESULTS_DIR = PROJECT_ROOT / "results" / "analysis" / "antibiotic_exposure_comparison"
FIGURES_DIR = PROJECT_ROOT / "figures" / "antibiotic_exposure_comparison"

# Create directories
RESULTS_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

# Set plotting style
sns.set_style("whitegrid")
sns.set_context("talk")

def load_metadata():
    """Load metadata with antibiotic exposure"""
    print("Loading metadata...")

    # Load DIAMOND metadata (updated version)
    diamond_metadata = pd.read_csv(DATA_DIR / "diamond_metadata_updated.tsv", sep="\t")

    # Load full metadata with antibiotic exposure
    full_metadata = pd.read_csv(DATA_DIR / "master_metadata_fixed.tsv", sep="\t")

    # Get antibiotic columns
    abx_w1_cols = [col for col in full_metadata.columns if col.endswith('_w1') and col != 'None_w1']
    abx_w2_cols = [col for col in full_metadata.columns if col.endswith('_w2') and col != 'None_w2']

    # Keep only relevant columns from full metadata
    keep_cols = ['sample_name', 'PatientID', 'Location'] + abx_w1_cols + abx_w2_cols
    keep_cols = [col for col in keep_cols if col in full_metadata.columns]
    full_metadata = full_metadata[keep_cols].copy()

    # Fix ZCH naming (add ZJH_ prefix)
    zch_mask = full_metadata['Location'] == 'ZCH'
    full_metadata.loc[zch_mask, 'sample_name'] = 'ZJH_' + full_metadata.loc[zch_mask, 'sample_name']

    # Merge with DIAMOND metadata to get SampleCollectionWeek
    metadata = diamond_metadata[['sample_name', 'SampleCollectionWeek', 'SampleType']].merge(
        full_metadata, on='sample_name', how='inner'
    )

    print(f"  Total samples: {len(metadata)}")
    print(f"  UCMC: {(metadata['Location']=='UCMC').sum()}")
    print(f"  ZCH: {(metadata['Location']=='ZCH').sum()}")
    print(f"  Week 1: {(metadata['SampleCollectionWeek']=='Week.1').sum()}")
    print(f"  Week 3: {(metadata['SampleCollectionWeek']=='Week.3').sum()}")

    # Fill NaN with 0
    for col in abx_w1_cols + abx_w2_cols:
        if col in metadata.columns:
            metadata[col] = metadata[col].fillna(0)

    return metadata, abx_w1_cols, abx_w2_cols

def prepare_data_for_plotting(metadata, abx_w1_cols, abx_w2_cols):
    """
    Prepare data for plotting:
    - Week 1 samples: use w1 exposure
    - Week 3 samples: use cumulative (w1 + w2) exposure
    """
    print("\nPreparing data for plotting...")

    # Map w2 columns to w1 columns
    abx_pairs = {}
    for w1_col in abx_w1_cols:
        abx_name = w1_col.replace('_w1', '')
        w2_col = abx_name + '_w2'
        if w2_col in abx_w2_cols:
            abx_pairs[abx_name] = (w1_col, w2_col)

    print(f"  Found {len(abx_pairs)} antibiotics")

    # Create long-format data for plotting
    plot_data = []

    for idx, row in metadata.iterrows():
        sample_name = row['sample_name']
        location = row['Location']
        week = row['SampleCollectionWeek']
        patient_id = row['PatientID']

        for abx_name, (w1_col, w2_col) in abx_pairs.items():
            # Week 1: w1 exposure only
            # Week 3: cumulative (w1 + w2)
            if week == 'Week.1':
                exposure = row[w1_col]
            elif week == 'Week.3':
                exposure = row[w1_col] + row[w2_col]
            else:
                continue

            plot_data.append({
                'sample_name': sample_name,
                'patient_id': patient_id,
                'antibiotic': abx_name,
                'location': location,
                'week': week,
                'exposure_days': exposure
            })

    plot_df = pd.DataFrame(plot_data)

    print(f"  Created plot data: {len(plot_df)} rows")

    # Remove antibiotics with no exposure
    abx_with_exposure = plot_df.groupby('antibiotic')['exposure_days'].sum()
    abx_to_keep = abx_with_exposure[abx_with_exposure > 0].index.tolist()
    plot_df = plot_df[plot_df['antibiotic'].isin(abx_to_keep)]

    print(f"  Antibiotics with exposure: {len(abx_to_keep)}")

    return plot_df, abx_to_keep

def calculate_statistics(df, antibiotic, week, location1='UCMC', location2='ZCH'):
    """Calculate summary stats and statistical test"""

    # Filter data
    data = df[(df['antibiotic'] == antibiotic) & (df['week'] == week)]

    loc1_data = data[data['location'] == location1]['exposure_days']
    loc2_data = data[data['location'] == location2]['exposure_days']

    # Summary statistics
    stats = {
        f'{location1}_n': len(loc1_data),
        f'{location1}_mean': loc1_data.mean(),
        f'{location1}_median': loc1_data.median(),
        f'{location1}_q25': loc1_data.quantile(0.25),
        f'{location1}_q75': loc1_data.quantile(0.75),
        f'{location1}_iqr': loc1_data.quantile(0.75) - loc1_data.quantile(0.25),
        f'{location2}_n': len(loc2_data),
        f'{location2}_mean': loc2_data.mean(),
        f'{location2}_median': loc2_data.median(),
        f'{location2}_q25': loc2_data.quantile(0.25),
        f'{location2}_q75': loc2_data.quantile(0.75),
        f'{location2}_iqr': loc2_data.quantile(0.75) - loc2_data.quantile(0.25),
    }

    # Mann-Whitney U test
    if len(loc1_data) > 0 and len(loc2_data) > 0:
        try:
            stat, pval = mannwhitneyu(loc1_data, loc2_data, alternative='two-sided')
            stats['mann_whitney_stat'] = stat
            stats['mann_whitney_pval'] = pval
        except:
            stats['mann_whitney_stat'] = np.nan
            stats['mann_whitney_pval'] = np.nan
    else:
        stats['mann_whitney_stat'] = np.nan
        stats['mann_whitney_pval'] = np.nan

    return stats

def plot_antibiotic_exposure(plot_df, antibiotic, output_dir):
    """Create boxplot for one antibiotic, faceted by week"""

    # Filter to this antibiotic
    data = plot_df[plot_df['antibiotic'] == antibiotic].copy()

    if len(data) == 0:
        return None

    # Create figure with two subplots (one per week)
    fig, axes = plt.subplots(1, 2, figsize=(14, 6), sharey=True)

    weeks = ['Week.1', 'Week.3']
    week_labels = ['Week 1', 'Week 3 (cumulative)']
    colors = {'UCMC': '#1f77b4', 'ZCH': '#ff7f0e'}

    for i, (week, week_label) in enumerate(zip(weeks, week_labels)):
        ax = axes[i]
        week_data = data[data['week'] == week]

        if len(week_data) == 0:
            ax.text(0.5, 0.5, 'No data', ha='center', va='center', transform=ax.transAxes)
            ax.set_title(f'{week_label}\nNo data')
            continue

        # Create boxplot
        positions = [0, 1]
        box_data = [
            week_data[week_data['location'] == 'UCMC']['exposure_days'],
            week_data[week_data['location'] == 'ZCH']['exposure_days']
        ]

        bp = ax.boxplot(box_data, positions=positions, widths=0.6,
                        patch_artist=True, showfliers=False)

        # Color boxes
        for patch, loc in zip(bp['boxes'], ['UCMC', 'ZCH']):
            patch.set_facecolor(colors[loc])
            patch.set_alpha(0.6)

        # Overlay individual points with jitter
        for j, loc in enumerate(['UCMC', 'ZCH']):
            loc_data = week_data[week_data['location'] == loc]
            y_vals = loc_data['exposure_days']
            x_vals = np.random.normal(j, 0.04, size=len(y_vals))
            ax.scatter(x_vals, y_vals, alpha=0.4, s=30, c=colors[loc], edgecolors='black', linewidths=0.3)

        ax.set_xticks([0, 1])
        ax.set_xticklabels(['UCMC', 'ZCH'])
        ax.set_ylabel('Exposure (days)')

        # Calculate statistics
        stats = calculate_statistics(plot_df, antibiotic, week)

        # Add statistics text
        pval = stats['mann_whitney_pval']
        if not np.isnan(pval):
            if pval < 0.001:
                pval_str = 'p < 0.001'
            else:
                pval_str = f'p = {pval:.3f}'

            # Add p-value at top of plot area (within margins, below title)
            # Use axes transform to place at consistent position
            ax.text(0.5, 0.95, pval_str, ha='center', va='top',
                   transform=ax.transAxes, fontsize=11, fontweight='bold',
                   bbox=dict(boxstyle='round', facecolor='white', alpha=0.8, edgecolor='black'))

        # Add mean/median/IQR text below x-axis
        ucmc_text = (f"n={stats['UCMC_n']}\n"
                     f"Mean: {stats['UCMC_mean']:.1f}\n"
                     f"Median: {stats['UCMC_median']:.1f}\n"
                     f"IQR: {stats['UCMC_iqr']:.1f}")
        zch_text = (f"n={stats['ZCH_n']}\n"
                    f"Mean: {stats['ZCH_mean']:.1f}\n"
                    f"Median: {stats['ZCH_median']:.1f}\n"
                    f"IQR: {stats['ZCH_iqr']:.1f}")

        ax.text(0, -0.25, ucmc_text, ha='center', va='top', transform=ax.get_xaxis_transform(),
                fontsize=9, bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        ax.text(1, -0.25, zch_text, ha='center', va='top', transform=ax.get_xaxis_transform(),
                fontsize=9, bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

        ax.set_title(week_label, fontsize=13, fontweight='bold')
        ax.grid(True, alpha=0.3, axis='y')

    plt.suptitle(f'{antibiotic} Exposure: UCMC vs ZCH', fontsize=15, fontweight='bold', y=1.02)
    plt.tight_layout(rect=[0, 0.05, 1, 0.98])

    # Save
    safe_name = antibiotic.replace('/', '_').replace(' ', '_')
    output_file = output_dir / f"{safe_name}_exposure_comparison.pdf"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.savefig(output_file.with_suffix('.png'), dpi=300, bbox_inches='tight')
    plt.close()

    return output_file

def create_summary_table(plot_df, antibiotics):
    """Create summary statistics table for all antibiotics"""
    print("\nCreating summary statistics table...")

    results = []

    for antibiotic in antibiotics:
        for week in ['Week.1', 'Week.3']:
            stats = calculate_statistics(plot_df, antibiotic, week)
            stats['antibiotic'] = antibiotic
            stats['week'] = week
            results.append(stats)

    results_df = pd.DataFrame(results)

    # Reorder columns
    col_order = ['antibiotic', 'week',
                 'UCMC_n', 'UCMC_mean', 'UCMC_median', 'UCMC_q25', 'UCMC_q75', 'UCMC_iqr',
                 'ZCH_n', 'ZCH_mean', 'ZCH_median', 'ZCH_q25', 'ZCH_q75', 'ZCH_iqr',
                 'mann_whitney_stat', 'mann_whitney_pval']
    results_df = results_df[col_order]

    # Save
    output_file = RESULTS_DIR / "antibiotic_exposure_comparison_statistics.tsv"
    results_df.to_csv(output_file, sep="\t", index=False)
    print(f"✓ Saved: {output_file}")

    # Show significant differences
    sig_results = results_df[results_df['mann_whitney_pval'] < 0.05].copy()
    sig_results = sig_results.sort_values('mann_whitney_pval')

    print(f"\nSignificant differences (p < 0.05): {len(sig_results)}")
    if len(sig_results) > 0:
        print("\nTop 10 most significant differences:")
        for i, (_, row) in enumerate(sig_results.head(10).iterrows(), 1):
            ucmc_median = row['UCMC_median']
            zch_median = row['ZCH_median']
            higher = 'UCMC' if ucmc_median > zch_median else 'ZCH'
            print(f"  {i}. {row['antibiotic']} ({row['week']}): "
                  f"UCMC median={ucmc_median:.1f}, ZCH median={zch_median:.1f}, "
                  f"p={row['mann_whitney_pval']:.2e} (higher in {higher})")

    return results_df

def main():
    print("="*60)
    print("ANTIBIOTIC EXPOSURE COMPARISON: UCMC vs ZCH")
    print("="*60)

    # Load metadata
    metadata, abx_w1_cols, abx_w2_cols = load_metadata()

    # Prepare data for plotting
    plot_df, antibiotics = prepare_data_for_plotting(metadata, abx_w1_cols, abx_w2_cols)

    # Create summary statistics table
    stats_df = create_summary_table(plot_df, antibiotics)

    # Create plots for each antibiotic
    print(f"\nCreating plots for {len(antibiotics)} antibiotics...")
    for i, antibiotic in enumerate(antibiotics, 1):
        output_file = plot_antibiotic_exposure(plot_df, antibiotic, FIGURES_DIR)
        if output_file:
            print(f"  {i}/{len(antibiotics)}: {antibiotic} → {output_file.name}")

    print("\n" + "="*60)
    print("ANALYSIS COMPLETE")
    print("="*60)
    print(f"\nOutputs:")
    print(f"  Statistics: {RESULTS_DIR / 'antibiotic_exposure_comparison_statistics.tsv'}")
    print(f"  Figures: {FIGURES_DIR}/ ({len(antibiotics)} plots)")

if __name__ == "__main__":
    main()
