#!/usr/bin/env python3
"""
Randomly select samples from NICU sample key for benchmarking
"""

import pandas as pd
import numpy as np
import argparse
from pathlib import Path

def main():
    parser = argparse.ArgumentParser(description='Randomly select NICU samples for benchmarking')
    parser.add_argument('--input', type=Path, required=True, help='Input CSV with all samples')
    parser.add_argument('--output', type=Path, required=True, help='Output CSV with selected samples')
    parser.add_argument('--n', type=int, default=50, help='Number of samples to select (default: 50)')
    parser.add_argument('--seed', type=int, default=42, help='Random seed for reproducibility (default: 42)')
    parser.add_argument('--verify-files', action='store_true', help='Verify FASTQ files exist')

    args = parser.parse_args()

    print(f"Loading samples from: {args.input}")
    df = pd.read_csv(args.input)

    print(f"Total samples in file: {len(df)}")

    # Pre-filter to only samples where files exist
    print("\nChecking which samples have existing FASTQ files...")
    valid_samples = []
    missing_count = 0

    for idx, row in df.iterrows():
        r1_path = Path(row['fastq_path']) / row['R1_file']
        r2_path = Path(row['fastq_path']) / row['R2_file']

        if r1_path.exists() and r2_path.exists():
            valid_samples.append(idx)
        else:
            missing_count += 1

    df_valid = df.loc[valid_samples]
    print(f"  ✓ Valid samples (files exist): {len(df_valid)}")
    print(f"  ✗ Invalid samples (files missing): {missing_count}")

    # Set random seed for reproducibility
    np.random.seed(args.seed)

    # Randomly select n samples from valid samples
    if args.n > len(df_valid):
        print(f"\nWarning: Requested {args.n} samples but only {len(df_valid)} valid. Using all valid samples.")
        selected = df_valid
    else:
        selected = df_valid.sample(n=args.n, random_state=args.seed)

    print(f"\nRandomly selected {len(selected)} samples (seed={args.seed})")

    # Verify files exist if requested (should always pass now)
    if args.verify_files:
        print("\nVerifying FASTQ files exist...")
        missing = []

        for idx, row in selected.iterrows():
            r1_path = Path(row['fastq_path']) / row['R1_file']
            r2_path = Path(row['fastq_path']) / row['R2_file']

            if not r1_path.exists():
                missing.append(f"R1: {r1_path}")
            if not r2_path.exists():
                missing.append(f"R2: {r2_path}")

        if missing:
            print(f"ERROR: {len(missing)} files not found:")
            for f in missing[:10]:  # Show first 10
                print(f"  {f}")
            if len(missing) > 10:
                print(f"  ... and {len(missing)-10} more")
            return 1
        else:
            print(f"✓ All {len(selected)*2} FASTQ files verified")

    # Create output with proper formatting
    # Use /bulkpool/sequence_data/mss_data/ explicitly for clarity
    output_df = selected.copy()

    # Format for benchmark pipeline: sample_name,r1_path,r2_path
    output_df['r1_path'] = output_df['fastq_path'] + output_df['R1_file']
    output_df['r2_path'] = output_df['fastq_path'] + output_df['R2_file']

    output_df = output_df[['sample_name', 'r1_path', 'r2_path']]

    # Save
    output_df.to_csv(args.output, index=False)
    print(f"\n✓ Saved selected samples to: {args.output}")

    # Print summary statistics
    print("\nSample summary:")
    print(f"  Sample name format: {output_df['sample_name'].iloc[0]}")
    print(f"  First R1 path: {output_df['r1_path'].iloc[0]}")
    print(f"  First R2 path: {output_df['r2_path'].iloc[0]}")

    # Show sample distribution
    print("\nSelected samples:")
    for sample in sorted(output_df['sample_name'].tolist()):
        print(f"  {sample}")

    return 0

if __name__ == '__main__':
    exit(main())
