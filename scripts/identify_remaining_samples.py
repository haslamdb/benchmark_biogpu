#!/usr/bin/env python3
"""
Identify remaining NICU samples that haven't been processed through DIAMOND yet
Creates a CSV of remaining samples in the same format as test_samples.csv
"""

import pandas as pd
from pathlib import Path
import argparse

def main():
    parser = argparse.ArgumentParser(description='Identify remaining unprocessed NICU samples')
    parser.add_argument('--full-list', type=Path,
                       default='/home/david/projects/biogpu/nicu_sample_key_fixed.csv',
                       help='Full sample list CSV')
    parser.add_argument('--processed-list', type=Path,
                       default='data/test_samples.csv',
                       help='Already processed samples CSV')
    parser.add_argument('--results-dir', type=Path,
                       default='results/traditional',
                       help='Results directory to check for processed samples')
    parser.add_argument('--output', type=Path,
                       default='data/remaining_samples.csv',
                       help='Output CSV for remaining samples')
    parser.add_argument('--verify-files', action='store_true',
                       help='Verify FASTQ files exist before including')

    args = parser.parse_args()

    print("=" * 80)
    print("Identifying Remaining NICU Samples for DIAMOND Processing")
    print("=" * 80)

    # Load full sample list
    print(f"\nLoading full sample list: {args.full_list}")
    full_df = pd.read_csv(args.full_list)
    print(f"  Total samples in full list: {len(full_df)}")

    # Create full paths
    full_df['r1_path'] = full_df['fastq_path'] + full_df['R1_file']
    full_df['r2_path'] = full_df['fastq_path'] + full_df['R2_file']
    full_df = full_df[['sample_name', 'r1_path', 'r2_path']]

    # Get already processed samples from the CSV
    if args.processed_list.exists():
        print(f"\nLoading processed sample list: {args.processed_list}")
        processed_df = pd.read_csv(args.processed_list)
        processed_samples = set(processed_df['sample_name'].tolist())
        print(f"  Samples in processed list: {len(processed_samples)}")
    else:
        processed_samples = set()
        print(f"\nNo processed list found at: {args.processed_list}")

    # Also check results directory for completed samples
    if args.results_dir.exists():
        print(f"\nChecking results directory: {args.results_dir}")
        result_dirs = [d.name for d in args.results_dir.iterdir() if d.is_dir()]
        processed_samples.update(result_dirs)
        print(f"  Samples with results: {len(result_dirs)}")
        print(f"  Total unique processed samples: {len(processed_samples)}")
    else:
        print(f"\nResults directory not found: {args.results_dir}")

    # Find remaining samples
    remaining_df = full_df[~full_df['sample_name'].isin(processed_samples)].copy()
    print(f"\n{'='*80}")
    print(f"REMAINING SAMPLES: {len(remaining_df)}")
    print(f"{'='*80}")

    if len(remaining_df) == 0:
        print("\n✓ All samples have been processed!")
        return 0

    # Verify files exist if requested
    if args.verify_files:
        print("\nVerifying FASTQ files exist...")
        valid_samples = []
        missing_count = 0

        for idx, row in remaining_df.iterrows():
            r1_path = Path(row['r1_path'])
            r2_path = Path(row['r2_path'])

            if r1_path.exists() and r2_path.exists():
                valid_samples.append(idx)
            else:
                missing_count += 1
                if missing_count <= 5:  # Show first 5 missing
                    print(f"  Missing: {row['sample_name']}")
                    if not r1_path.exists():
                        print(f"    R1: {r1_path}")
                    if not r2_path.exists():
                        print(f"    R2: {r2_path}")

        remaining_df = remaining_df.loc[valid_samples]

        if missing_count > 5:
            print(f"  ... and {missing_count - 5} more samples with missing files")

        print(f"\n  ✓ Valid samples (files exist): {len(remaining_df)}")
        print(f"  ✗ Invalid samples (files missing): {missing_count}")

    # Save remaining samples
    remaining_df.to_csv(args.output, index=False)
    print(f"\n✓ Saved remaining samples to: {args.output}")

    # Show sample distribution
    print("\nSample breakdown:")
    ucmc_samples = remaining_df['sample_name'].str.startswith('N').sum()
    zch_samples = remaining_df['sample_name'].str.startswith('ZJH').sum()
    print(f"  UCMC (N*): {ucmc_samples}")
    print(f"  ZCH (ZJH_*): {zch_samples}")

    # Show first 10 samples
    print("\nFirst 10 remaining samples:")
    for sample in remaining_df['sample_name'].head(10):
        print(f"  {sample}")

    if len(remaining_df) > 10:
        print(f"  ... and {len(remaining_df) - 10} more")

    return 0

if __name__ == '__main__':
    exit(main())
