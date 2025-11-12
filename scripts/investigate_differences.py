#!/usr/bin/env python3
"""
Investigate differences between DIAMOND and BioGPU results
"""

import pandas as pd
import numpy as np
from pathlib import Path

# Load merged data
merged = pd.read_csv('results/comparison/merged_comparison_data.tsv', sep='\t')

print("=" * 80)
print("Investigation of DIAMOND vs BioGPU Differences")
print("=" * 80)

# 1. Detection patterns
print("\n## 1. Gene Detection Patterns\n")

# Genes detected only by DIAMOND
diamond_only = merged[merged['_merge'] == 'left_only']
biogpu_only = merged[merged['_merge'] == 'right_only']
both = merged[merged['_merge'] == 'both']

print(f"Total unique genes across all samples:")
print(f"  DIAMOND only: {len(diamond_only):,} ({100*len(diamond_only)/len(merged):.1f}%)")
print(f"  BioGPU only: {len(biogpu_only):,} ({100*len(biogpu_only)/len(merged):.1f}%)")
print(f"  Both: {len(both):,} ({100*len(both)/len(merged):.1f}%)")

# Per-sample detection
print(f"\nPer-sample detection statistics:")
samples = merged['sample'].unique()
for sample in sorted(samples)[:5]:  # Show first 5
    sample_data = merged[merged['sample'] == sample]
    d_only = (sample_data['_merge'] == 'left_only').sum()
    b_only = (sample_data['_merge'] == 'right_only').sum()
    both_count = (sample_data['_merge'] == 'both').sum()
    print(f"  {sample}: DIAMOND={d_only}, BioGPU={b_only}, Both={both_count}")

# 2. RPM distribution comparison
print("\n## 2. RPM Distribution for Shared Genes\n")

# Filter to genes detected by both with non-zero counts
shared = both[(both['diamond_rpm'] > 0) & (both['biogpu_rpm'] > 0)]

print(f"Genes with non-zero RPM in both methods: {len(shared):,}")
print(f"\nRPM Summary Statistics:")
print(f"\nDIAMOND:")
print(f"  Mean: {shared['diamond_rpm'].mean():.2f}")
print(f"  Median: {shared['diamond_rpm'].median():.2f}")
print(f"  Min: {shared['diamond_rpm'].min():.4f}")
print(f"  Max: {shared['diamond_rpm'].max():.2f}")

print(f"\nBioGPU:")
print(f"  Mean: {shared['biogpu_rpm'].mean():.2f}")
print(f"  Median: {shared['biogpu_rpm'].median():.2f}")
print(f"  Min: {shared['biogpu_rpm'].min():.4f}")
print(f"  Max: {shared['biogpu_rpm'].max():.2f}")

# 3. Coverage comparison
print("\n## 3. Coverage Comparison for Shared Genes\n")

print(f"Coverage Summary (for genes detected by both):")
print(f"\nDIAMOND:")
print(f"  Mean coverage: {shared['diamond_coverage'].mean():.2f}%")
print(f"  Median coverage: {shared['diamond_coverage'].median():.2f}%")

print(f"\nBioGPU:")
print(f"  Mean coverage: {shared['biogpu_coverage'].mean():.2f}%")
print(f"  Median coverage: {shared['biogpu_coverage'].median():.2f}%")

# 4. Top genes by method
print("\n## 4. Top 10 Genes by RPM (Each Method)\n")

print("\nTop 10 DIAMOND-only genes (highest RPM):")
top_diamond = diamond_only.nlargest(10, 'diamond_rpm')[['sample', 'gene_name', 'diamond_rpm', 'diamond_coverage']]
print(top_diamond.to_string(index=False))

print("\nTop 10 BioGPU-only genes (highest RPM):")
top_biogpu = biogpu_only.nlargest(10, 'biogpu_rpm')[['sample', 'gene_name', 'biogpu_rpm', 'biogpu_coverage']]
print(top_biogpu.to_string(index=False))

# 5. Agreement by abundance category
print("\n## 5. Agreement by Abundance Category\n")

# Categorize shared genes by abundance
shared['abundance_category'] = pd.cut(
    shared['diamond_rpm'],
    bins=[0, 1, 10, 100, np.inf],
    labels=['Very Low (0-1)', 'Low (1-10)', 'Medium (10-100)', 'High (>100)']
)

print("Correlation by abundance category:")
for category in ['Very Low (0-1)', 'Low (1-10)', 'Medium (10-100)', 'High (>100)']:
    cat_data = shared[shared['abundance_category'] == category]
    if len(cat_data) > 3:
        from scipy.stats import spearmanr
        rho, pval = spearmanr(cat_data['diamond_rpm'], cat_data['biogpu_rpm'])
        print(f"  {category}: n={len(cat_data)}, ρ={rho:.3f}, p={pval:.2e}")

# 6. Gene families comparison
print("\n## 6. Most Discordant Gene Families\n")

# Calculate discordance
shared['ratio'] = shared[['diamond_rpm', 'biogpu_rpm']].max(axis=1) / \
                 (shared[['diamond_rpm', 'biogpu_rpm']].min(axis=1) + 0.01)
shared['abs_diff'] = abs(shared['diamond_rpm'] - shared['biogpu_rpm'])

# Group by gene_name to find consistently discordant genes
gene_stats = shared.groupby('gene_name').agg({
    'ratio': 'mean',
    'abs_diff': 'mean',
    'sample': 'count'
}).rename(columns={'sample': 'n_samples'})

print("\nGenes with highest average RPM ratio (detected in multiple samples):")
top_discordant = gene_stats[gene_stats['n_samples'] >= 3].nlargest(20, 'ratio')
print(top_discordant.to_string())

# 7. Check if DIAMOND is detecting things BioGPU filters
print("\n## 7. Low-Abundance Detection Thresholds\n")

print("\nDIAMOND genes with very low RPM (< 0.1):")
low_diamond = diamond_only[diamond_only['diamond_rpm'] < 0.1]
print(f"  Count: {len(low_diamond):,} ({100*len(low_diamond)/len(diamond_only):.1f}% of DIAMOND-only)")

print("\nBioGPU genes with very low RPM (< 0.1):")
low_biogpu = biogpu_only[biogpu_only['biogpu_rpm'] < 0.1]
print(f"  Count: {len(low_biogpu):,} ({100*len(low_biogpu)/len(biogpu_only):.1f}% of BioGPU-only)")

# 8. Coverage threshold analysis
print("\n## 8. Coverage Threshold Analysis\n")

print("\nDIAMOND-only genes by coverage:")
cov_bins = [0, 10, 25, 50, 75, 100]
for i in range(len(cov_bins)-1):
    count = ((diamond_only['diamond_coverage'] >= cov_bins[i]) &
             (diamond_only['diamond_coverage'] < cov_bins[i+1])).sum()
    print(f"  {cov_bins[i]}-{cov_bins[i+1]}%: {count:,} genes")

print("\nBioGPU-only genes by coverage:")
for i in range(len(cov_bins)-1):
    count = ((biogpu_only['biogpu_coverage'] >= cov_bins[i]) &
             (biogpu_only['biogpu_coverage'] < cov_bins[i+1])).sum()
    print(f"  {cov_bins[i]}-{cov_bins[i+1]}%: {count:,} genes")

print("\n" + "=" * 80)
print("Investigation complete!")
print("=" * 80)
