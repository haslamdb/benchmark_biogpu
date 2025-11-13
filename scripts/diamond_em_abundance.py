#!/usr/bin/env python3
"""
Apply EM Algorithm to DIAMOND Output

Resolves multi-mapping reads using Expectation-Maximization,
similar to Kallisto and BioGPU's approach.

This allows fair comparison: DIAMOND alignment + BioGPU-style EM quantification

Usage:
    python diamond_em_abundance.py \
        --diamond-output sample_diamond.tsv \
        --output sample_em_abundance.tsv \
        --min-identity 85 \
        --min-coverage 50
"""

import pandas as pd
import numpy as np
from collections import defaultdict
from pathlib import Path
import argparse
import sys

class DiamondEMQuantifier:
    """EM-based abundance quantification for DIAMOND alignments"""

    def __init__(self, min_identity=85, min_coverage=50,
                 max_iterations=100, convergence_threshold=0.001):
        self.min_identity = min_identity
        self.min_coverage = min_coverage
        self.max_iterations = max_iterations
        self.convergence_threshold = convergence_threshold

        # Data structures
        self.read_assignments = {}  # read_id -> [(gene, score), ...]
        self.gene_lengths = {}      # gene -> length (amino acids)
        self.gene_abundances = {}   # gene -> abundance (RPKM-like)
        self.gene_read_counts = {}  # gene -> fractional read count

    def load_diamond_alignments(self, diamond_file):
        """
        Load DIAMOND blastx output
        Expected format: qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore qcovhsp
        """
        print(f"Loading DIAMOND alignments from {diamond_file}...")

        df = pd.read_csv(diamond_file, sep='\t', header=None,
                         names=['qseqid', 'sseqid', 'pident', 'length',
                               'qlen', 'slen', 'qstart', 'qend',
                               'sstart', 'send', 'evalue', 'bitscore', 'qcovhsp'])

        print(f"  Loaded {len(df)} alignments")

        # Filter by identity and coverage
        df_filtered = df[
            (df['pident'] >= self.min_identity) &
            (df['qcovhsp'] >= self.min_coverage)
        ].copy()

        print(f"  After filtering (≥{self.min_identity}% ID, ≥{self.min_coverage}% cov): {len(df_filtered)} alignments")

        if len(df_filtered) == 0:
            print("ERROR: No alignments passed filters!")
            return False

        # Extract gene name from sseqid (format: 0|protein_id|...|gene_name|...)
        df_filtered['gene_name'] = df_filtered['sseqid'].str.split('|').str[4]

        # Build read assignments
        for _, row in df_filtered.iterrows():
            read_id = row['qseqid']
            gene = row['gene_name']

            # Calculate assignment score (identity * coverage * alignment length)
            score = row['pident'] * row['qcovhsp'] * row['length']

            # Store gene length (protein length in amino acids)
            if gene not in self.gene_lengths:
                self.gene_lengths[gene] = row['slen']

            # Add to read assignments
            if read_id not in self.read_assignments:
                self.read_assignments[read_id] = []
            self.read_assignments[read_id].append((gene, score))

        # Get unique reads
        unique_reads = len(self.read_assignments)
        multi_mapping = sum(1 for assignments in self.read_assignments.values() if len(assignments) > 1)

        print(f"  Unique reads: {unique_reads}")
        print(f"  Multi-mapping reads: {multi_mapping} ({100*multi_mapping/unique_reads:.1f}%)")
        print(f"  Unique genes: {len(self.gene_lengths)}")

        return True

    def initialize_abundances(self):
        """Initialize gene abundances uniformly"""
        print("\nInitializing gene abundances...")

        # Start with uniform abundance across all genes
        for gene in self.gene_lengths:
            self.gene_abundances[gene] = 1.0
            self.gene_read_counts[gene] = 0.0

        # Count reads that uniquely map to each gene (for better initialization)
        unique_counts = defaultdict(int)
        for read_id, assignments in self.read_assignments.items():
            if len(assignments) == 1:
                gene, _ = assignments[0]
                unique_counts[gene] += 1

        # Use unique counts as initial abundances (with pseudocount)
        for gene in self.gene_lengths:
            self.gene_abundances[gene] = unique_counts.get(gene, 0) + 0.1

        print(f"  Initialized {len(self.gene_abundances)} genes")
        genes_with_unique = sum(1 for count in unique_counts.values() if count > 0)
        print(f"  Genes with unique reads: {genes_with_unique}")

    def expectation_step(self):
        """E-step: Update read assignment probabilities based on current abundances"""

        # Reset gene read counts
        for gene in self.gene_read_counts:
            self.gene_read_counts[gene] = 0.0

        # For each read, assign fractionally to candidate genes
        for read_id, assignments in self.read_assignments.items():
            if len(assignments) == 1:
                # Unique mapping - assign full read
                gene, _ = assignments[0]
                self.gene_read_counts[gene] += 1.0
            else:
                # Multi-mapping - assign fractionally based on abundance and score
                genes, scores = zip(*assignments)

                # Calculate weights: score * abundance / gene_length
                weights = []
                for gene, score in assignments:
                    # Abundance normalized by gene length (longer genes produce more fragments)
                    gene_length_kb = self.gene_lengths[gene] / 1000.0
                    weight = score * (self.gene_abundances[gene] / gene_length_kb)
                    weights.append(weight)

                # Normalize to probabilities
                total_weight = sum(weights)
                if total_weight > 0:
                    probabilities = [w / total_weight for w in weights]
                else:
                    # Fallback to uniform
                    probabilities = [1.0 / len(assignments)] * len(assignments)

                # Assign fractional counts
                for gene, prob in zip(genes, probabilities):
                    self.gene_read_counts[gene] += prob

    def maximization_step(self):
        """M-step: Update gene abundances based on fractional read assignments"""

        # Calculate total assigned reads (for normalization)
        total_reads = sum(self.gene_read_counts.values())

        if total_reads == 0:
            print("WARNING: No reads assigned in M-step!")
            return

        # Update abundances: RPKM-like calculation
        # RPKM = (reads * 1e6) / (gene_length_kb * total_reads)
        for gene in self.gene_abundances:
            read_count = self.gene_read_counts[gene]
            gene_length_kb = self.gene_lengths[gene] / 1000.0

            # RPKM with pseudocount
            rpkm = ((read_count + 0.01) * 1e6) / (gene_length_kb * (total_reads + 1.0))
            self.gene_abundances[gene] = rpkm

    def run_em(self):
        """Run EM algorithm until convergence"""
        print("\nRunning EM algorithm...")

        prev_abundances = None

        for iteration in range(self.max_iterations):
            # Store previous abundances for convergence check
            prev_abundances = dict(self.gene_abundances)

            # E-step
            self.expectation_step()

            # M-step
            self.maximization_step()

            # Check convergence
            max_change = 0.0
            for gene in self.gene_abundances:
                if self.gene_read_counts[gene] > 0.1:  # Only check genes with reads
                    change = abs(self.gene_abundances[gene] - prev_abundances[gene])
                    rel_change = change / (prev_abundances[gene] + 1e-10)
                    max_change = max(max_change, rel_change)

            if iteration % 10 == 0 or iteration < 5:
                print(f"  Iteration {iteration}: max relative change = {max_change:.6f}")

            # Convergence check
            if iteration > 5 and max_change < self.convergence_threshold:
                print(f"  ✓ Converged after {iteration + 1} iterations")
                break
        else:
            print(f"  ! Reached maximum iterations ({self.max_iterations})")

    def export_abundance_table(self, output_file, total_reads_processed):
        """Export final abundance table"""
        print(f"\nExporting abundance table to {output_file}...")

        # Calculate additional metrics
        results = []

        for gene in sorted(self.gene_abundances.keys()):
            read_count = self.gene_read_counts[gene]

            # Only export genes with reads
            if read_count > 0:
                gene_length = self.gene_lengths[gene]
                gene_length_kb = gene_length / 1000.0

                # RPM (reads per million)
                rpm = (read_count * 1e6) / (total_reads_processed + 1.0)

                # TPM calculation (normalized by gene length first, then by total)
                rpk = read_count / gene_length_kb

                # RPKM
                rpkm = (read_count * 1e9) / (gene_length * 3 * total_reads_processed)

                results.append({
                    'gene_name': gene,
                    'read_count': round(read_count, 2),
                    'rpm': round(rpm, 4),
                    'rpkm': round(rpkm, 4),
                    'gene_length': gene_length,
                    'em_abundance': round(self.gene_abundances[gene], 4),
                    'evidence': 'EM'
                })

        # Calculate TPM (sum of RPK values, then normalize)
        total_rpk = sum((r['read_count'] / (r['gene_length'] / 1000.0)) for r in results)
        for r in results:
            rpk = r['read_count'] / (r['gene_length'] / 1000.0)
            r['tpm'] = round((rpk / total_rpk * 1e6) if total_rpk > 0 else 0, 4)

        # Create DataFrame and save
        df = pd.DataFrame(results)
        df = df.sort_values('read_count', ascending=False)

        df.to_csv(output_file, sep='\t', index=False)

        print(f"  ✓ Exported {len(df)} genes with reads")
        print(f"    Total fractional reads assigned: {sum(df['read_count']):.1f}")
        print(f"    Mean reads per gene: {df['read_count'].mean():.2f}")
        print(f"    Median reads per gene: {df['read_count'].median():.2f}")

        return df


def main():
    parser = argparse.ArgumentParser(
        description='Apply EM algorithm to DIAMOND output for abundance quantification',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example:
    %(prog)s \\
        --diamond-output sample_diamond.tsv \\
        --output sample_em_abundance.tsv \\
        --min-identity 85 \\
        --min-coverage 50 \\
        --max-iterations 100

DIAMOND output format should be:
    qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore qcovhsp

Generated using:
    diamond blastx --outfmt 6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore qcovhsp
        """
    )

    parser.add_argument('--diamond-output', required=True, type=Path,
                       help='DIAMOND blastx output file (tab-separated)')
    parser.add_argument('--output', required=True, type=Path,
                       help='Output abundance table')
    parser.add_argument('--min-identity', type=float, default=85,
                       help='Minimum percent identity (default: 85)')
    parser.add_argument('--min-coverage', type=float, default=50,
                       help='Minimum query coverage (default: 50)')
    parser.add_argument('--max-iterations', type=int, default=100,
                       help='Maximum EM iterations (default: 100)')
    parser.add_argument('--convergence', type=float, default=0.001,
                       help='Convergence threshold (default: 0.001)')
    parser.add_argument('--total-reads', type=int, default=None,
                       help='Total reads processed (for RPM calculation, default: auto-detect)')

    args = parser.parse_args()

    # Check input file exists
    if not args.diamond_output.exists():
        print(f"ERROR: Input file not found: {args.diamond_output}")
        sys.exit(1)

    print("="*60)
    print("DIAMOND EM Abundance Quantification")
    print("="*60)
    print(f"Input: {args.diamond_output}")
    print(f"Output: {args.output}")
    print(f"Filters: ≥{args.min_identity}% identity, ≥{args.min_coverage}% coverage")
    print(f"EM parameters: max_iter={args.max_iterations}, convergence={args.convergence}")
    print()

    # Initialize quantifier
    quantifier = DiamondEMQuantifier(
        min_identity=args.min_identity,
        min_coverage=args.min_coverage,
        max_iterations=args.max_iterations,
        convergence_threshold=args.convergence
    )

    # Load alignments
    if not quantifier.load_diamond_alignments(args.diamond_output):
        sys.exit(1)

    # Determine total reads
    if args.total_reads:
        total_reads = args.total_reads
    else:
        # Count unique read IDs
        total_reads = len(quantifier.read_assignments)

    print(f"\nTotal reads for normalization: {total_reads}")

    # Initialize and run EM
    quantifier.initialize_abundances()
    quantifier.run_em()

    # Export results
    quantifier.export_abundance_table(args.output, total_reads)

    print("\n" + "="*60)
    print("EM quantification complete!")
    print("="*60)


if __name__ == "__main__":
    main()
