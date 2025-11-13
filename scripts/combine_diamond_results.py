#!/usr/bin/env python3
"""
Combine DIAMOND abundance results with gene annotations
Creates a unified dataset similar to biogpu format for downstream analysis
"""

import pandas as pd
import numpy as np
from pathlib import Path
import sys

def main():
    # Paths
    diamond_results = Path("/home/david/projects/benchmark_biogpu/results/traditional")
    annotation_file = Path("/home/david/projects/benchmark_biogpu/data/gene_annotations.tsv")
    output_dir = Path("/home/david/projects/benchmark_biogpu/data")
    output_file = output_dir / "diamond_amr_combined.tsv"

    print("Combining DIAMOND results with gene annotations...")
    print(f"DIAMOND results: {diamond_results}")
    print(f"Annotation file: {annotation_file}")

    # Load gene annotations
    annotations = pd.read_csv(annotation_file, sep="\t")
    print(f"\nLoaded {len(annotations)} gene annotations")

    # Find all DIAMOND abundance files
    abundance_files = list(diamond_results.glob("*/*_abundance.tsv"))
    print(f"Found {len(abundance_files)} DIAMOND abundance files")

    if len(abundance_files) == 0:
        print("ERROR: No DIAMOND abundance files found!")
        sys.exit(1)

    # Collect all data
    all_data = []

    for file in abundance_files:
        # Extract sample name from path (e.g., N01_1_2 from .../N01_1_2/N01_1_2_abundance.tsv)
        sample_name = file.parent.name

        try:
            df = pd.read_csv(file, sep="\t")

            # Filter to genes with read_count > 0
            df = df[df['read_count'] > 0].copy()

            # Add sample_name column
            df['sample_name'] = sample_name

            # Select and reorder columns to match biogpu format
            df = df[['sample_name', 'gene_name', 'read_count', 'rpm', 'mean_depth', 'percent_coverage']]

            all_data.append(df)

        except Exception as e:
            print(f"Warning: Could not read {file}: {e}")
            continue

    # Combine all samples
    combined = pd.concat(all_data, ignore_index=True)
    print(f"\nCombined data: {len(combined)} entries from {len(abundance_files)} samples")

    # Merge with gene annotations
    combined = combined.merge(
        annotations[['gene_name', 'gene_family', 'drug_class']],
        on='gene_name',
        how='left'
    )

    # Reorder columns to match biogpu format
    combined = combined[['sample_name', 'gene_name', 'gene_family', 'drug_class',
                         'read_count', 'percent_coverage', 'mean_depth', 'rpm']]

    # Sort by sample_name and gene_name
    combined = combined.sort_values(['sample_name', 'gene_name'])

    # Save
    output_dir.mkdir(parents=True, exist_ok=True)
    combined.to_csv(output_file, sep="\t", index=False)

    print(f"\nCombined DIAMOND results saved:")
    print(f"  Total entries: {len(combined)}")
    print(f"  Unique samples: {combined['sample_name'].nunique()}")
    print(f"  Unique genes: {combined['gene_name'].nunique()}")
    print(f"  Entries with drug class: {combined['drug_class'].notna().sum()}")
    print(f"  Output: {output_file}")

    # Summary statistics
    print("\nSummary by drug class:")
    drug_class_summary = combined.groupby('drug_class', dropna=False).agg({
        'gene_name': 'nunique',
        'read_count': 'sum',
        'rpm': 'sum'
    }).sort_values('rpm', ascending=False)
    print(drug_class_summary.head(15))

if __name__ == "__main__":
    main()
