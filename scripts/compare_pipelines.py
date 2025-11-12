#!/usr/bin/env python3
"""
Compare DIAMOND and BioGPU pipeline results
============================================

This script compares gene detection and abundance between:
- DIAMOND blastx (traditional translated search)
- BioGPU (GPU-accelerated translated search)

Both pipelines use identical parameters (85% identity, 50% coverage)
"""

import pandas as pd
import numpy as np
from scipy.stats import spearmanr, pearsonr
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import sys
import argparse
from typing import Tuple, Dict, List

# Set style
sns.set_style("whitegrid")
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300

class PipelineComparison:
    """Compare DIAMOND and BioGPU pipeline results"""

    def __init__(self, base_dir: Path):
        self.base_dir = base_dir
        self.results_dir = base_dir / "results"
        self.diamond_dir = self.results_dir / "traditional"
        self.biogpu_dir = self.results_dir / "biogpu"
        self.output_dir = self.results_dir / "comparison"
        self.output_dir.mkdir(exist_ok=True)

        self.samples = []
        self.diamond_data = None
        self.biogpu_data = None
        self.merged_data = None
        self.timing_comparison = None

    def load_diamond_results(self) -> pd.DataFrame:
        """Load all DIAMOND pipeline results"""
        print("Loading DIAMOND results...")
        all_data = []

        for sample_dir in sorted(self.diamond_dir.glob("*/")):
            sample_name = sample_dir.name
            abundance_file = sample_dir / f"{sample_name}_abundance.tsv"

            if not abundance_file.exists():
                print(f"  Warning: Missing {abundance_file}")
                continue

            df = pd.read_csv(abundance_file, sep='\t')
            df['sample'] = sample_name
            df['pipeline'] = 'DIAMOND'
            all_data.append(df)

        if not all_data:
            raise ValueError("No DIAMOND results found!")

        combined = pd.concat(all_data, ignore_index=True)
        print(f"  Loaded {len(all_data)} samples, {len(combined)} total records")
        return combined

    def load_biogpu_results(self) -> pd.DataFrame:
        """Load all BioGPU pipeline results"""
        print("Loading BioGPU results...")
        all_data = []

        for sample_dir in sorted(self.biogpu_dir.glob("*/")):
            sample_name = sample_dir.name
            # BioGPU has nested directory structure
            # Use amr_abundance.tsv which has rpm column
            abundance_file = sample_dir / sample_name / f"{sample_name}_amr_abundance.tsv"

            if not abundance_file.exists():
                print(f"  Warning: Missing {abundance_file}")
                continue

            try:
                df = pd.read_csv(abundance_file, sep='\t', encoding='utf-8')
            except UnicodeDecodeError:
                # Try with latin-1 encoding if UTF-8 fails
                print(f"  Warning: UTF-8 decode error, trying latin-1 for {sample_name}")
                df = pd.read_csv(abundance_file, sep='\t', encoding='latin-1')

            df['sample'] = sample_name
            df['pipeline'] = 'BioGPU'
            all_data.append(df)

        if not all_data:
            raise ValueError("No BioGPU results found!")

        combined = pd.concat(all_data, ignore_index=True)
        print(f"  Loaded {len(all_data)} samples, {len(combined)} total records")
        return combined

    def load_timing_data(self) -> pd.DataFrame:
        """Load and compare timing metrics"""
        print("Loading timing data...")

        # Load DIAMOND timing
        diamond_timing = []
        for sample_dir in sorted(self.diamond_dir.glob("*/")):
            sample_name = sample_dir.name
            timing_file = sample_dir / f"{sample_name}_timing.tsv"
            if timing_file.exists():
                df = pd.read_csv(timing_file, sep='\t')
                df = df[df['step'] == 'TOTAL']
                df['sample'] = sample_name
                df['pipeline'] = 'DIAMOND'
                diamond_timing.append(df)

        # Load BioGPU timing
        biogpu_timing = []
        for sample_dir in sorted(self.biogpu_dir.glob("*/")):
            sample_name = sample_dir.name
            timing_file = sample_dir / f"{sample_name}_timing.tsv"
            if timing_file.exists():
                df = pd.read_csv(timing_file, sep='\t')
                df = df[df['step'] == 'TOTAL']
                df['sample'] = sample_name
                df['pipeline'] = 'BioGPU'
                biogpu_timing.append(df)

        if diamond_timing and biogpu_timing:
            timing = pd.concat(diamond_timing + biogpu_timing, ignore_index=True)
            print(f"  Loaded timing for {len(diamond_timing)} DIAMOND and {len(biogpu_timing)} BioGPU samples")
            return timing
        else:
            print("  Warning: No timing data found")
            return None

    def merge_results(self) -> pd.DataFrame:
        """Merge DIAMOND and BioGPU results for comparison"""
        print("Merging results...")

        # Prepare DIAMOND data - FILTER to only detected genes (read_count > 0)
        diamond = self.diamond_data[self.diamond_data['read_count'] > 0].copy()
        diamond = diamond[['sample', 'gene_name', 'read_count', 'rpm', 'tpm',
                          'mean_depth', 'percent_coverage', 'gene_length']].copy()
        diamond.columns = ['sample', 'gene_name', 'diamond_read_count', 'diamond_rpm',
                          'diamond_tpm', 'diamond_depth', 'diamond_coverage', 'gene_length']

        # Prepare BioGPU data - already filtered to detected genes
        # Note: using percent_coverage and mean_depth from amr_abundance.tsv
        biogpu = self.biogpu_data[['sample', 'gene_name', 'read_count', 'rpm',
                                    'tpkm', 'percent_coverage', 'mean_depth']].copy()
        biogpu.columns = ['sample', 'gene_name', 'biogpu_read_count', 'biogpu_rpm',
                         'biogpu_tpkm', 'biogpu_coverage', 'biogpu_depth']

        print(f"  DIAMOND detected genes (read_count > 0): {len(diamond)}")
        print(f"  BioGPU detected genes: {len(biogpu)}")

        # Merge on sample and gene_name
        merged = pd.merge(diamond, biogpu, on=['sample', 'gene_name'], how='outer', indicator=True)

        # Add detection flags
        merged['detected_diamond'] = merged['_merge'].isin(['both', 'left_only'])
        merged['detected_biogpu'] = merged['_merge'].isin(['both', 'right_only'])
        merged['detected_both'] = merged['_merge'] == 'both'

        # Fill NaN values with 0 for counts/abundances
        for col in merged.columns:
            if 'count' in col or 'rpm' in col or 'tpm' in col or 'rpkm' in col:
                merged[col] = merged[col].fillna(0)

        print(f"  Merged data: {len(merged)} total gene-sample combinations")
        print(f"  Detected by DIAMOND only: {(merged['_merge'] == 'left_only').sum()}")
        print(f"  Detected by BioGPU only: {(merged['_merge'] == 'right_only').sum()}")
        print(f"  Detected by both: {(merged['_merge'] == 'both').sum()}")

        return merged

    def calculate_correlations(self) -> Dict[str, Tuple[float, float]]:
        """Calculate correlation statistics"""
        print("\nCalculating correlations...")

        # Filter to genes detected by both methods with non-zero counts
        both = self.merged_data[
            (self.merged_data['detected_both']) &
            (self.merged_data['diamond_rpm'] > 0) &
            (self.merged_data['biogpu_rpm'] > 0)
        ].copy()

        print(f"  Genes detected by both methods (non-zero): {len(both)}")

        results = {}

        if len(both) > 3:
            # RPM correlation
            pearson_rpm = pearsonr(both['diamond_rpm'], both['biogpu_rpm'])
            spearman_rpm = spearmanr(both['diamond_rpm'], both['biogpu_rpm'])
            results['rpm'] = {
                'pearson': pearson_rpm,
                'spearman': spearman_rpm,
                'n': len(both)
            }
            print(f"  RPM - Pearson: {pearson_rpm[0]:.4f} (p={pearson_rpm[1]:.2e})")
            print(f"  RPM - Spearman: {spearman_rpm[0]:.4f} (p={spearman_rpm[1]:.2e})")

            # Read count correlation
            pearson_count = pearsonr(both['diamond_read_count'], both['biogpu_read_count'])
            spearman_count = spearmanr(both['diamond_read_count'], both['biogpu_read_count'])
            results['count'] = {
                'pearson': pearson_count,
                'spearman': spearman_count,
                'n': len(both)
            }
            print(f"  Count - Pearson: {pearson_count[0]:.4f} (p={pearson_count[1]:.2e})")
            print(f"  Count - Spearman: {spearman_count[0]:.4f} (p={spearman_count[1]:.2e})")

            # Log-transformed correlation (for better visualization)
            both['log_diamond_rpm'] = np.log10(both['diamond_rpm'] + 1)
            both['log_biogpu_rpm'] = np.log10(both['biogpu_rpm'] + 1)
            pearson_log = pearsonr(both['log_diamond_rpm'], both['log_biogpu_rpm'])
            spearman_log = spearmanr(both['log_diamond_rpm'], both['log_biogpu_rpm'])
            results['log_rpm'] = {
                'pearson': pearson_log,
                'spearman': spearman_log,
                'n': len(both)
            }
            print(f"  Log RPM - Pearson: {pearson_log[0]:.4f} (p={pearson_log[1]:.2e})")
            print(f"  Log RPM - Spearman: {spearman_log[0]:.4f} (p={spearman_log[1]:.2e})")
        else:
            print("  Warning: Not enough shared genes for correlation analysis")

        return results

    def plot_scatter_comparison(self, output_file: Path):
        """Create scatter plots comparing DIAMOND vs BioGPU"""
        print(f"\nCreating scatter plots...")

        # Filter to genes detected by both
        both = self.merged_data[
            (self.merged_data['detected_both']) &
            (self.merged_data['diamond_rpm'] > 0) &
            (self.merged_data['biogpu_rpm'] > 0)
        ].copy()

        fig, axes = plt.subplots(1, 2, figsize=(14, 6))

        # Plot 1: Linear scale
        ax = axes[0]
        ax.scatter(both['diamond_rpm'], both['biogpu_rpm'], alpha=0.3, s=10)

        # Add diagonal line
        max_val = max(both['diamond_rpm'].max(), both['biogpu_rpm'].max())
        ax.plot([0, max_val], [0, max_val], 'r--', linewidth=1, alpha=0.7, label='y=x')

        ax.set_xlabel('DIAMOND RPM', fontsize=12)
        ax.set_ylabel('BioGPU RPM', fontsize=12)
        ax.set_title('RPM Comparison (Linear Scale)', fontsize=13, fontweight='bold')
        ax.legend()

        # Add correlation text
        if hasattr(self, 'correlations') and 'rpm' in self.correlations:
            corr = self.correlations['rpm']
            text = f"n = {corr['n']:,}\n"
            text += f"Pearson r = {corr['pearson'][0]:.3f}\n"
            text += f"Spearman ρ = {corr['spearman'][0]:.3f}"
            ax.text(0.05, 0.95, text, transform=ax.transAxes,
                   verticalalignment='top', fontsize=10,
                   bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

        # Plot 2: Log scale
        ax = axes[1]
        both['log_diamond'] = np.log10(both['diamond_rpm'] + 1)
        both['log_biogpu'] = np.log10(both['biogpu_rpm'] + 1)

        ax.scatter(both['log_diamond'], both['log_biogpu'], alpha=0.3, s=10)

        # Add diagonal line
        max_log = max(both['log_diamond'].max(), both['log_biogpu'].max())
        ax.plot([0, max_log], [0, max_log], 'r--', linewidth=1, alpha=0.7, label='y=x')

        ax.set_xlabel('DIAMOND log10(RPM + 1)', fontsize=12)
        ax.set_ylabel('BioGPU log10(RPM + 1)', fontsize=12)
        ax.set_title('RPM Comparison (Log Scale)', fontsize=13, fontweight='bold')
        ax.legend()

        # Add correlation text
        if hasattr(self, 'correlations') and 'log_rpm' in self.correlations:
            corr = self.correlations['log_rpm']
            text = f"n = {corr['n']:,}\n"
            text += f"Pearson r = {corr['pearson'][0]:.3f}\n"
            text += f"Spearman ρ = {corr['spearman'][0]:.3f}"
            ax.text(0.05, 0.95, text, transform=ax.transAxes,
                   verticalalignment='top', fontsize=10,
                   bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"  Saved: {output_file}")
        plt.close()

    def plot_venn_diagram(self, output_file: Path):
        """Create Venn diagram of gene detection overlap"""
        print(f"\nCreating Venn diagram...")

        # Count unique genes per sample
        diamond_only = self.merged_data[self.merged_data['_merge'] == 'left_only'].groupby('sample').size()
        biogpu_only = self.merged_data[self.merged_data['_merge'] == 'right_only'].groupby('sample').size()
        both = self.merged_data[self.merged_data['_merge'] == 'both'].groupby('sample').size()

        fig, ax = plt.subplots(figsize=(10, 6))

        # Simple bar plot showing overlap
        samples = sorted(set(diamond_only.index) | set(biogpu_only.index) | set(both.index))
        x = np.arange(len(samples))
        width = 0.25

        diamond_vals = [diamond_only.get(s, 0) for s in samples]
        biogpu_vals = [biogpu_only.get(s, 0) for s in samples]
        both_vals = [both.get(s, 0) for s in samples]

        ax.bar(x - width, diamond_vals, width, label='DIAMOND only', alpha=0.8)
        ax.bar(x, both_vals, width, label='Both', alpha=0.8)
        ax.bar(x + width, biogpu_vals, width, label='BioGPU only', alpha=0.8)

        ax.set_xlabel('Sample', fontsize=12)
        ax.set_ylabel('Number of Genes', fontsize=12)
        ax.set_title('Gene Detection by Pipeline', fontsize=13, fontweight='bold')
        ax.set_xticks(x)
        ax.set_xticklabels(samples, rotation=45, ha='right', fontsize=8)
        ax.legend()

        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"  Saved: {output_file}")
        plt.close()

    def plot_timing_comparison(self, output_file: Path):
        """Compare pipeline timing"""
        if self.timing_comparison is None:
            print("  Skipping timing plot (no data)")
            return

        print(f"\nCreating timing comparison...")

        fig, axes = plt.subplots(1, 3, figsize=(15, 5))

        # Plot 1: Wall time comparison
        ax = axes[0]
        diamond_time = self.timing_comparison[self.timing_comparison['pipeline'] == 'DIAMOND']['wall_time_sec']
        biogpu_time = self.timing_comparison[self.timing_comparison['pipeline'] == 'BioGPU']['wall_time_sec']

        bp = ax.boxplot([diamond_time, biogpu_time], labels=['DIAMOND', 'BioGPU'])
        ax.set_ylabel('Wall Time (seconds)', fontsize=12)
        ax.set_title('Pipeline Runtime', fontsize=13, fontweight='bold')

        # Add mean values
        ax.text(1, diamond_time.mean(), f'μ={diamond_time.mean():.1f}s',
               ha='center', va='bottom', fontsize=9)
        ax.text(2, biogpu_time.mean(), f'μ={biogpu_time.mean():.1f}s',
               ha='center', va='bottom', fontsize=9)

        # Plot 2: Memory comparison
        ax = axes[1]
        diamond_mem = self.timing_comparison[self.timing_comparison['pipeline'] == 'DIAMOND']['memory_gb']
        biogpu_mem = self.timing_comparison[self.timing_comparison['pipeline'] == 'BioGPU']['memory_gb']

        ax.boxplot([diamond_mem, biogpu_mem], labels=['DIAMOND', 'BioGPU'])
        ax.set_ylabel('Peak Memory (GB)', fontsize=12)
        ax.set_title('Memory Usage', fontsize=13, fontweight='bold')

        # Plot 3: Speedup distribution
        ax = axes[2]
        # Merge timing by sample
        timing_merged = self.timing_comparison.pivot_table(
            index='sample', columns='pipeline', values='wall_time_sec'
        )
        timing_merged['speedup'] = timing_merged['DIAMOND'] / timing_merged['BioGPU']

        ax.hist(timing_merged['speedup'], bins=20, edgecolor='black', alpha=0.7)
        ax.axvline(timing_merged['speedup'].mean(), color='red', linestyle='--',
                  linewidth=2, label=f'Mean: {timing_merged["speedup"].mean():.2f}x')
        ax.axvline(1.0, color='black', linestyle=':', linewidth=1, label='Equal time')
        ax.set_xlabel('Speedup (DIAMOND time / BioGPU time)', fontsize=12)
        ax.set_ylabel('Number of Samples', fontsize=12)
        ax.set_title('Performance Ratio Distribution', fontsize=13, fontweight='bold')
        ax.legend()

        # Add text summary
        text = f"Mean speedup: {timing_merged['speedup'].mean():.2f}x\n"
        text += f"Median: {timing_merged['speedup'].median():.2f}x\n"
        if timing_merged['speedup'].mean() > 1:
            text += "DIAMOND is faster"
        else:
            text += "BioGPU is faster"
        ax.text(0.95, 0.95, text, transform=ax.transAxes,
               verticalalignment='top', horizontalalignment='right',
               fontsize=10, bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"  Saved: {output_file}")
        plt.close()

    def generate_summary_report(self, output_file: Path):
        """Generate comprehensive summary report"""
        print(f"\nGenerating summary report...")

        with open(output_file, 'w') as f:
            f.write("=" * 80 + "\n")
            f.write("DIAMOND vs BioGPU Pipeline Comparison\n")
            f.write("=" * 80 + "\n\n")

            # Overview
            f.write("## Overview\n\n")
            f.write(f"- Samples analyzed: {len(self.merged_data['sample'].unique())}\n")
            f.write(f"- Total gene-sample combinations: {len(self.merged_data):,}\n\n")

            # Detection statistics
            f.write("## Gene Detection\n\n")
            diamond_only = (self.merged_data['_merge'] == 'left_only').sum()
            biogpu_only = (self.merged_data['_merge'] == 'right_only').sum()
            both = (self.merged_data['_merge'] == 'both').sum()
            total = len(self.merged_data)

            f.write(f"- Detected by DIAMOND only: {diamond_only:,} ({100*diamond_only/total:.1f}%)\n")
            f.write(f"- Detected by BioGPU only: {biogpu_only:,} ({100*biogpu_only/total:.1f}%)\n")
            f.write(f"- Detected by both: {both:,} ({100*both/total:.1f}%)\n\n")

            # Correlation statistics
            if hasattr(self, 'correlations'):
                f.write("## Correlation Statistics\n\n")

                if 'rpm' in self.correlations:
                    corr = self.correlations['rpm']
                    f.write("### RPM (Reads Per Million)\n")
                    f.write(f"- N = {corr['n']:,} genes detected by both methods\n")
                    f.write(f"- Pearson r = {corr['pearson'][0]:.4f} (p = {corr['pearson'][1]:.2e})\n")
                    f.write(f"- Spearman ρ = {corr['spearman'][0]:.4f} (p = {corr['spearman'][1]:.2e})\n\n")

                if 'log_rpm' in self.correlations:
                    corr = self.correlations['log_rpm']
                    f.write("### Log-transformed RPM\n")
                    f.write(f"- Pearson r = {corr['pearson'][0]:.4f} (p = {corr['pearson'][1]:.2e})\n")
                    f.write(f"- Spearman ρ = {corr['spearman'][0]:.4f} (p = {corr['spearman'][1]:.2e})\n\n")

            # Timing comparison
            if self.timing_comparison is not None:
                f.write("## Performance Comparison\n\n")

                diamond_time = self.timing_comparison[self.timing_comparison['pipeline'] == 'DIAMOND']
                biogpu_time = self.timing_comparison[self.timing_comparison['pipeline'] == 'BioGPU']

                f.write("### Runtime (Wall Time)\n")
                f.write(f"- DIAMOND: {diamond_time['wall_time_sec'].mean():.1f} ± {diamond_time['wall_time_sec'].std():.1f} seconds\n")
                f.write(f"- BioGPU: {biogpu_time['wall_time_sec'].mean():.1f} ± {biogpu_time['wall_time_sec'].std():.1f} seconds\n")

                timing_merged = self.timing_comparison.pivot_table(
                    index='sample', columns='pipeline', values='wall_time_sec'
                )
                timing_merged['speedup'] = timing_merged['DIAMOND'] / timing_merged['BioGPU']
                f.write(f"- Mean speedup: {timing_merged['speedup'].mean():.2f}x ")
                if timing_merged['speedup'].mean() > 1:
                    f.write("(DIAMOND faster)\n\n")
                else:
                    f.write("(BioGPU faster)\n\n")

                f.write("### Memory Usage\n")
                f.write(f"- DIAMOND: {diamond_time['memory_gb'].mean():.2f} ± {diamond_time['memory_gb'].std():.2f} GB\n")
                f.write(f"- BioGPU: {biogpu_time['memory_gb'].mean():.2f} ± {biogpu_time['memory_gb'].std():.2f} GB\n\n")

            # Top discordant genes
            f.write("## Top Discordant Genes\n\n")
            both_detected = self.merged_data[self.merged_data['detected_both']].copy()
            both_detected['rpm_diff'] = abs(both_detected['diamond_rpm'] - both_detected['biogpu_rpm'])
            both_detected['rpm_ratio'] = both_detected[['diamond_rpm', 'biogpu_rpm']].max(axis=1) / \
                                        (both_detected[['diamond_rpm', 'biogpu_rpm']].min(axis=1) + 0.01)

            top_discordant = both_detected.nlargest(20, 'rpm_ratio')[
                ['sample', 'gene_name', 'diamond_rpm', 'biogpu_rpm', 'rpm_ratio']
            ]

            f.write("Top 20 genes with largest RPM ratio (may indicate systematic differences):\n\n")
            f.write(top_discordant.to_string(index=False))
            f.write("\n\n")

            # Summary
            f.write("## Summary\n\n")
            f.write("Both DIAMOND and BioGPU use translated search (protein space) with identical\n")
            f.write("parameters (85% identity, 50% coverage). This comparison validates:\n\n")
            f.write("1. Gene detection consistency between CPU and GPU implementations\n")
            f.write("2. Abundance quantification correlation\n")
            f.write("3. Performance differences (speed and memory)\n\n")

            if hasattr(self, 'correlations') and 'rpm' in self.correlations:
                corr_val = self.correlations['rpm']['spearman'][0]
                if corr_val > 0.9:
                    f.write(f"High correlation (ρ={corr_val:.3f}) validates BioGPU implementation.\n")
                elif corr_val > 0.7:
                    f.write(f"Good correlation (ρ={corr_val:.3f}) with some systematic differences to investigate.\n")
                else:
                    f.write(f"Moderate correlation (ρ={corr_val:.3f}) - significant differences warrant investigation.\n")

        print(f"  Saved: {output_file}")

    def save_merged_data(self, output_file: Path):
        """Save merged comparison data"""
        print(f"\nSaving merged data...")
        self.merged_data.to_csv(output_file, sep='\t', index=False)
        print(f"  Saved: {output_file}")

    def run_comparison(self):
        """Run complete comparison analysis"""
        print("\n" + "=" * 80)
        print("DIAMOND vs BioGPU Pipeline Comparison")
        print("=" * 80 + "\n")

        # Load data
        self.diamond_data = self.load_diamond_results()
        self.biogpu_data = self.load_biogpu_results()
        self.timing_comparison = self.load_timing_data()

        # Merge and analyze
        self.merged_data = self.merge_results()
        self.correlations = self.calculate_correlations()

        # Generate outputs
        self.plot_scatter_comparison(self.output_dir / "scatter_rpm_comparison.png")
        self.plot_venn_diagram(self.output_dir / "gene_detection_overlap.png")
        self.plot_timing_comparison(self.output_dir / "timing_comparison.png")
        self.generate_summary_report(self.output_dir / "comparison_summary.txt")
        self.save_merged_data(self.output_dir / "merged_comparison_data.tsv")

        print("\n" + "=" * 80)
        print("Comparison complete!")
        print(f"Results saved to: {self.output_dir}")
        print("=" * 80 + "\n")


def main():
    parser = argparse.ArgumentParser(
        description='Compare DIAMOND and BioGPU pipeline results'
    )
    parser.add_argument(
        '--base-dir',
        type=Path,
        default=Path.cwd(),
        help='Base directory of benchmark project (default: current directory)'
    )

    args = parser.parse_args()

    if not args.base_dir.exists():
        print(f"Error: Directory not found: {args.base_dir}")
        sys.exit(1)

    comparison = PipelineComparison(args.base_dir)
    comparison.run_comparison()


if __name__ == '__main__':
    main()
