#!/usr/bin/env python3
"""
Create GFF3 annotation file from biogpu AMR+stress FASTA headers
for use with htseq-count
"""

import sys
import re
import argparse
from pathlib import Path

def parse_biogpu_header(header):
    """
    Parse biogpu FASTA header format:
    >index|protein_id|nucleotide_id|?|?|gene_name|gene_name2|description coordinates

    Example:
    >0|AAA16360.1|L11078.1|1|1|stxA2b|stxA2b|Shiga_toxin_Stx2b_subunit_A L11078.1:177-1136
    """
    # Remove '>'
    header = header.strip().lstrip('>')

    # Split by pipe
    parts = header.split('|')

    if len(parts) < 7:
        return None

    index = parts[0]
    protein_id = parts[1]
    nucleotide_id = parts[2]
    gene_name = parts[5]

    # Extract gene_name2 and description (after 6th pipe)
    remainder = parts[6]
    gene_name2 = remainder.split()[0] if remainder else gene_name

    # Extract coordinates if present (format: accession:start-end)
    coords_match = re.search(r'(\S+):(\d+)-(\d+)', remainder)
    if coords_match:
        seqid = coords_match.group(1)
        start = int(coords_match.group(2))
        end = int(coords_match.group(3))
        length = end - start + 1
    else:
        seqid = nucleotide_id
        start = 1
        end = 1000  # Placeholder
        length = 1000

    return {
        'index': index,
        'protein_id': protein_id,
        'nucleotide_id': nucleotide_id,
        'gene_name': gene_name,
        'gene_name2': gene_name2,
        'seqid': seqid,
        'start': start,
        'end': end,
        'length': length
    }


def create_gff3(fasta_file, output_file):
    """
    Create GFF3 file from FASTA file

    The key insight: htseq-count matches the "seqname" column in GFF3
    to the reference name in the SAM file. Since bowtie2 will use the
    FASTA header as the reference name, we need to use a unique identifier
    that matches what bowtie2 will produce.

    Bowtie2 uses the first word of the FASTA header (up to whitespace) as
    the reference name, so we'll use that same format.
    """

    with open(output_file, 'w') as gff:
        # Write GFF3 header
        gff.write("##gff-version 3\n")

        num_genes = 0

        with open(fasta_file, 'r') as fasta:
            for line in fasta:
                if line.startswith('>'):
                    parsed = parse_biogpu_header(line)

                    if parsed is None:
                        print(f"Warning: Could not parse header: {line.strip()}", file=sys.stderr)
                        continue

                    # Create unique ID that matches bowtie2 SAM output
                    # Bowtie2 uses everything up to first whitespace
                    # We need to extract the full pipe-delimited string before the space
                    # From the original line, get everything after '>' and before first space
                    full_id = line.strip().lstrip('>').split()[0]

                    # Use the full ID for both seqname and gene_id to match bowtie2
                    seqname = full_id
                    gene_id = full_id

                    # Write GFF3 line
                    # Format: seqname source feature start end score strand frame attributes
                    # Since we don't have real coordinates, use 1 to length
                    gff.write(
                        f"{seqname}\t"
                        f"biogpu\t"
                        f"gene\t"
                        f"1\t"
                        f"{parsed['length']}\t"
                        f".\t"
                        f"+\t"
                        f".\t"
                        f"ID={gene_id};"
                        f"Name={parsed['gene_name']};"
                        f"gene={parsed['gene_name']};"
                        f"locus_tag={parsed['protein_id']}\n"
                    )

                    num_genes += 1

    return num_genes


def main():
    parser = argparse.ArgumentParser(
        description='Create GFF3 annotation file from biogpu FASTA headers'
    )
    parser.add_argument(
        'fasta_file',
        type=Path,
        help='Input FASTA file (e.g., amr_stress_dna.fasta)'
    )
    parser.add_argument(
        'output_file',
        type=Path,
        help='Output GFF3 file (e.g., amr_stress_genes.gff3)'
    )

    args = parser.parse_args()

    if not args.fasta_file.exists():
        print(f"Error: FASTA file not found: {args.fasta_file}", file=sys.stderr)
        sys.exit(1)

    print(f"Creating GFF3 from: {args.fasta_file}")
    num_genes = create_gff3(args.fasta_file, args.output_file)
    print(f"✓ Created GFF3 file: {args.output_file}")
    print(f"  Number of genes: {num_genes}")


if __name__ == '__main__':
    main()
