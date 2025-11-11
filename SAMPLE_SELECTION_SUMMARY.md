# Sample Selection Summary

## Overview

Successfully selected **50 random NICU samples** for comprehensive benchmarking analysis.

## Source Data

**Original dataset**: `/home/david/projects/biogpu/nicu_sample_key_fixed.csv`
- Total samples in CSV: **771**
- Valid samples (files exist): **644** (83.5%)
- Invalid samples (files missing): **127** (16.5%)

**FASTQ file location**: `/bulkpool/sequence_data/mss_data/`
- Symlinked from: `/home/david/sequence_data/mss_data/`

## Sample Selection

**Method**: Random sampling with fixed seed for reproducibility
- Seed: **42**
- Samples selected: **50**
- Selection from: 644 valid samples only

**Output file**: `/home/david/projects/benchmark_biogpu/data/test_samples.csv`

## Selected Samples

### Distribution

**By location**:
- UCMC (Cincinnati): 28 samples (56%)
  - Format: `N##_#_#` (e.g., N01_1_4, N02_2_4)
- ZCH (Hangzhou): 22 samples (44%)
  - Format: `ZJH_N##_#_#` (e.g., ZJH_N57_1_3, ZJH_N50_1_3)

### UCMC Samples (28)
```
N01_1_4, N02_2_4, N07_1_3, N10_2_4, N11_1_2, N12_1_4, N13_2_2, N14_2_2,
N15_1_3, N15_2_3, N20_1_3, N20_2_3, N22_1_2, N29_1_3, N33_1_4, N34_2_2,
N35_1_4, N40_1_2, N40_2_3, N41_1_2, N46_2_4, N52_1_2, N55_2_3, N56_1_3,
N56_1_4, N59_1_2, N62_1_2, N68_1_4
```

### ZCH Samples (22)
```
ZJH_N01_1_3, ZJH_N02_1_2, ZJH_N03_2_2, ZJH_N07_1_4, ZJH_N08_2_4, ZJH_N16_1_3,
ZJH_N16_2_4, ZJH_N17_1_4, ZJH_N26_1_4, ZJH_N31_1_2, ZJH_N35_1_2, ZJH_N35_1_4,
ZJH_N39_2_4, ZJH_N42_2_2, ZJH_N46_2_3, ZJH_N50_1_3, ZJH_N51_1_2, ZJH_N51_1_4,
ZJH_N55_1_2, ZJH_N56_2_2, ZJH_N57_1_3, ZJH_N59_2_3
```

## Sample Format

The naming convention appears to be:
- `[Location_Prefix_]N[Patient#]_[Week#]_[BodySite#]`

Example interpretation:
- `N56_1_3`: Patient 56, Week 1, Body site 3
- `ZJH_N57_1_3`: ZCH hospital, Patient 57, Week 1, Body site 3

Body sites (based on NICU analysis):
- Site 2: Likely Axilla
- Site 3: Likely Groin
- Site 4: Likely Stool

## Output File Format

CSV with columns: `sample_name,r1_path,r2_path`

Example:
```csv
sample_name,r1_path,r2_path
N01_1_4,/home/david/sequence_data/mss_data/N01_1_4_R1.fastq.gz,/home/david/sequence_data/mss_data/N01_1_4_R2.fastq.gz
```

## Next Steps

### 1. Run Database Setup (if not done)
```bash
conda activate benchmark_biogpu
cd /home/david/projects/benchmark_biogpu
./scripts/01_setup_databases.sh
```

**Creates**:
- Bowtie2 index for AMR+stress database
- GFF3 annotation file for htseq-count
- Takes ~10-15 minutes

### 2. Run Traditional Pipeline on All 50 Samples
```bash
# Using batch script
./scripts/03_batch_traditional_pipeline.sh
```

**Input**: `data/test_samples.csv` (automatically detected)
**Output**: `results/traditional/{sample}/` for each sample
**Time estimate**: ~15-30 minutes per sample = 12-25 hours total

**Per sample output**:
- `{sample}_abundance.tsv` - Gene counts, RPM, TPM, coverage
- `{sample}_timing.tsv` - Performance metrics
- `{sample}_stats.txt` - Human-readable summary

### 3. Run BioGPU Pipeline on Same 50 Samples
```bash
# Process each sample (can parallelize with multiple GPUs)
while IFS=, read -r sample_name r1_path r2_path; do
    if [[ "$sample_name" == "sample_name" ]]; then continue; fi

    ./scripts/04_run_biogpu_with_timing.sh \
        --sample "${sample_name}" \
        --r1 "${r1_path}" \
        --r2 "${r2_path}" \
        --gpu-id 0
done < data/test_samples.csv
```

**Output**: `results/biogpu/{sample}/` for each sample
**Time estimate**: ~5-10 minutes per sample = 4-8 hours total

**Per sample output**:
- `{sample}_amr_abundance.tsv` - AMR gene abundances
- `{sample}_mutations.tsv` - FQ resistance mutations
- `{sample}_timing.tsv` - Performance metrics
- `{sample}_stats.txt` - Human-readable summary

### 4. Compare Results

Create comparison analysis:
```python
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, pearsonr

# Load traditional results
trad_files = glob.glob('results/traditional/*/*.abundance.tsv')
trad_data = {}
for f in trad_files:
    sample = Path(f).parent.name
    df = pd.read_csv(f, sep='\t')
    trad_data[sample] = df

# Load biogpu results
biogpu_files = glob.glob('results/biogpu/*/*.amr_abundance.tsv')
biogpu_data = {}
for f in biogpu_files:
    sample = Path(f).parent.name
    df = pd.read_csv(f, sep='\t')
    biogpu_data[sample] = df

# Compare RPM values for shared genes
# ... correlation analysis, venn diagrams, etc.
```

### 5. Compare Timing/Performance

```python
# Load timing data
trad_timing = pd.read_csv('results/traditional/*/timing.tsv', sep='\t')
biogpu_timing = pd.read_csv('results/biogpu/*/timing.tsv', sep='\t')

# Calculate speedup
speedup = trad_timing['total_time'] / biogpu_timing['total_time']
print(f"Average speedup: {speedup.mean():.2f}x")
```

## Expected Results

### Gene Detection
- **High overlap**: >80% of genes detected by both methods
- **BioGPU unique**: Divergent sequences (70-90% nucleotide identity)
- **Traditional unique**: Should be minimal

### Abundance Correlation
- **Expected**: Spearman r > 0.90 for genes detected by both
- **Validates**: Both methods give similar quantification

### Performance
- **Expected speedup**: 3-8x faster with BioGPU
- **Memory**: Similar (8-10GB both)
- **Scalability**: BioGPU advantage increases with file size

## Reproducibility

**Random seed**: 42
- To select different samples: change `--seed` value
- To select same samples: use `--seed 42`

**Script location**: `/home/david/projects/benchmark_biogpu/scripts/select_random_samples.py`

**Usage**:
```bash
python3 scripts/select_random_samples.py \
    --input /home/david/projects/biogpu/nicu_sample_key_fixed.csv \
    --output data/test_samples.csv \
    --n 50 \
    --seed 42
```

## Files Created

- ✓ `data/test_samples.csv` - 50 randomly selected samples
- ✓ `scripts/select_random_samples.py` - Sample selection script
- ✓ `SAMPLE_SELECTION_SUMMARY.md` - This document

## Ready to Run!

All infrastructure is in place:
- ✓ Sample selection complete (50 samples)
- ✓ Scripts ready with timing instrumentation
- ✓ Database setup script ready
- ✓ Documentation complete

**Estimated total time**:
- Database setup: 15 min
- Traditional pipeline: 12-25 hours (can run overnight)
- BioGPU pipeline: 4-8 hours
- Analysis: 2-4 hours

**Total**: ~1-2 days of compute time
