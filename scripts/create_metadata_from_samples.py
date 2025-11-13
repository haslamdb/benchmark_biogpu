#!/usr/bin/env python3
"""
Create metadata file for DIAMOND samples
Parses sample names to extract Location, SampleType, and Week
"""

import pandas as pd
from pathlib import Path
import re

def parse_sample_name(sample_name):
    """
    Parse sample name to extract metadata
    Format: N##_WEEK_SITE or ZJH_N##_WEEK_SITE
    Week: 1=Week.1, 2=Week.3
    Site: 2=Axilla, 3=Groin, 4=Stool
    """
    # Check if ZJH (ZCH hospital)
    if sample_name.startswith("ZJH_"):
        location = "ZCH"
        # Remove ZJH_ prefix
        sample_core = sample_name[4:]
    else:
        location = "UCMC"
        sample_core = sample_name

    # Parse the sample ID: N##_WEEK_SITE
    match = re.match(r'N(\d+)_(\d+)_(\d+)', sample_core)
    if not match:
        return None

    subject_num = match.group(1)
    week_num = int(match.group(2))
    site_num = int(match.group(3))

    # Map site number to body site
    site_map = {
        2: "Axilla",
        3: "Groin",
        4: "Stool"
    }

    # Map week number to collection week
    week_map = {
        1: "Week.1",
        2: "Week.3"
    }

    subject_id = f"N{subject_num}"
    if location == "ZCH":
        subject_id = f"ZJH_{subject_id}"

    return {
        'sample_name': sample_name,
        'SubjectID': subject_id,
        'Location': location,
        'SampleType': site_map.get(site_num, "Unknown"),
        'SampleCollectionWeek': week_map.get(week_num, "Unknown")
    }

def main():
    # Get list of all samples
    diamond_data = Path("/home/david/projects/benchmark_biogpu/data/diamond_amr_combined.tsv")
    output_file = Path("/home/david/projects/benchmark_biogpu/data/diamond_metadata.tsv")

    print("Creating metadata from DIAMOND sample names...")

    # Read sample names
    df = pd.read_csv(diamond_data, sep="\t")
    sample_names = df['sample_name'].unique()

    print(f"Found {len(sample_names)} unique samples")

    # Parse all sample names
    metadata_list = []
    for sample in sample_names:
        parsed = parse_sample_name(sample)
        if parsed:
            metadata_list.append(parsed)
        else:
            print(f"Warning: Could not parse {sample}")

    # Create metadata dataframe
    metadata = pd.DataFrame(metadata_list)

    # Sort by location, sample type, and sample name
    metadata = metadata.sort_values(['Location', 'SampleType', 'sample_name'])

    # Save
    metadata.to_csv(output_file, sep="\t", index=False)

    print(f"\nMetadata created:")
    print(f"  Total samples: {len(metadata)}")
    print(f"  Output: {output_file}")

    # Summary
    print("\nSummary by Location:")
    print(metadata.groupby('Location').size())

    print("\nSummary by Body Site:")
    print(metadata.groupby(['Location', 'SampleType']).size())

    print("\nSummary by Week:")
    print(metadata.groupby(['Location', 'SampleCollectionWeek']).size())

    # Show sample
    print("\nSample entries:")
    print(metadata.head(10).to_string(index=False))

if __name__ == "__main__":
    main()
