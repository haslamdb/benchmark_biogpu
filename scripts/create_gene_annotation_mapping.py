#!/usr/bin/env python3
"""
Create gene annotation mapping from biogpu results
Extracts gene_name, gene_family, and drug_class for all genes
"""

import pandas as pd
from pathlib import Path
import sys

def main():
    # Paths
    biogpu_results = Path("/home/david/projects/benchmark_biogpu/results/biogpu")
    output_file = Path("/home/david/projects/benchmark_biogpu/data/gene_annotations.tsv")

    print("Creating gene annotation mapping from biogpu results...")

    # Find all biogpu AMR abundance files
    amr_files = list(biogpu_results.rglob("*_amr_abundance.tsv"))
    print(f"Found {len(amr_files)} biogpu AMR abundance files")

    if len(amr_files) == 0:
        print("ERROR: No biogpu AMR abundance files found!")
        sys.exit(1)

    # Collect all unique gene annotations
    all_genes = []

    for file in amr_files:
        try:
            df = pd.read_csv(file, sep="\t")
            # Select only annotation columns
            gene_info = df[['gene_name', 'gene_family', 'drug_class']].drop_duplicates()
            all_genes.append(gene_info)
        except Exception as e:
            print(f"Warning: Could not read {file}: {e}")
            continue

    # Combine and deduplicate
    gene_annotations = pd.concat(all_genes, ignore_index=True)
    gene_annotations = gene_annotations.drop_duplicates(subset=['gene_name'])

    # Sort by gene_name
    gene_annotations = gene_annotations.sort_values('gene_name')

    # Create output directory if needed
    output_file.parent.mkdir(parents=True, exist_ok=True)

    # Save
    gene_annotations.to_csv(output_file, sep="\t", index=False)

    print(f"\nGene annotation mapping created:")
    print(f"  Total unique genes: {len(gene_annotations)}")
    print(f"  Genes with drug class: {(gene_annotations['drug_class'] != '').sum()}")
    print(f"  Output: {output_file}")

    # Show sample
    print("\nSample entries:")
    print(gene_annotations.head(10).to_string(index=False))

if __name__ == "__main__":
    main()
