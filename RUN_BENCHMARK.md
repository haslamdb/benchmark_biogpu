# Run Benchmark - Complete Execution Guide

## Overview

This guide walks through executing the complete benchmark analysis of **50 randomly selected NICU samples** comparing biogpu vs traditional pipelines.

---

## Prerequisites Checklist

- [ ] Conda/mamba installed
- [ ] FASTQ files accessible at `/bulkpool/sequence_data/mss_data/`
- [ ] GPU available for biogpu pipeline
- [ ] Sufficient disk space (~50GB for results)

---

## Step-by-Step Execution

### Step 1: Create Conda Environment (5 minutes)

```bash
cd /home/david/projects/benchmark_biogpu

# Create environment
mamba env create -f environment.yml

# Activate
conda activate benchmark_biogpu

# Verify tools installed
bowtie2 --version
htseq-count --version
samtools --version
python --version
```

**Expected output**: All tools should report versions

---

### Step 2: Setup Databases (10-15 minutes)

```bash
./scripts/01_setup_databases.sh
```

**What this does**:
- Links biogpu AMR+stress databases
- Copies DNA FASTA (14,280 genes)
- Creates GFF3 annotation file for htseq-count
- Builds bowtie2 index (~5-10 min)

**Output**:
```
databases/traditional/
  ├── amr_stress_dna.fasta          # 16 MB
  ├── amr_stress_genes.gff3         # Annotations
  └── amr_stress_bt2.*.bt2          # Bowtie2 index files
```

**Check**: Verify completion
```bash
ls -lh databases/traditional/
# Should see FASTA, GFF3, and *.bt2 files
```

---

### Step 3: Verify Sample Selection (done!)

Sample selection is already complete:
```bash
wc -l data/test_samples.csv
# Should show: 51 (50 samples + header)

head -5 data/test_samples.csv
```

**Details**:
- 50 random samples selected (seed=42)
- 28 UCMC samples (56%)
- 22 ZCH samples (44%)
- All files verified to exist

---

### Step 4: Run Traditional Pipeline (12-25 hours)

#### Option A: Run in screen/tmux (Recommended)
```bash
# Start screen session
screen -S benchmark_traditional

# Activate environment
conda activate benchmark_biogpu
cd /home/david/projects/benchmark_biogpu

# Run batch processing
./scripts/03_batch_traditional_pipeline.sh

# Detach: Ctrl+A, then D
# Reattach later: screen -r benchmark_traditional
```

#### Option B: Run in background
```bash
nohup ./scripts/03_batch_traditional_pipeline.sh > traditional_batch.log 2>&1 &

# Check progress
tail -f traditional_batch.log
```

**Time estimate**: 15-30 minutes per sample
- 50 samples × 20 min avg = ~17 hours

**Monitoring progress**:
```bash
# Count completed samples
ls results/traditional/*/`*_timing.tsv | wc -l

# Check most recent log
tail -f logs/batch_traditional_*.log
```

**Output per sample**:
```
results/traditional/{sample}/
  ├── {sample}_abundance.tsv      # Counts, RPM, TPM, coverage
  ├── {sample}_timing.tsv         # Performance metrics
  ├── {sample}_stats.txt          # Human-readable summary
  └── {sample}_*_time.txt         # Detailed timing
```

---

### Step 5: Run BioGPU Pipeline (4-8 hours)

After traditional pipeline completes (or in parallel on different machine):

#### Option A: Run in screen/tmux (Recommended)
```bash
# Start screen session
screen -S benchmark_biogpu

# Activate environment (if needed)
conda activate benchmark_biogpu
cd /home/david/projects/benchmark_biogpu

# Run batch processing
./scripts/05_batch_biogpu_with_timing.sh --gpu-id 0

# Detach: Ctrl+A, then D
# Reattach later: screen -r benchmark_biogpu
```

#### Option B: Run with GPU profiling (slower but detailed)
```bash
./scripts/05_batch_biogpu_with_timing.sh --gpu-id 0 --profile
```

**Time estimate**: 5-10 minutes per sample
- 50 samples × 7 min avg = ~6 hours

**Monitoring progress**:
```bash
# Count completed samples
ls results/biogpu/*/`*_timing.tsv | wc -l

# Check most recent log
tail -f logs/batch_biogpu_*.log
```

**Output per sample**:
```
results/biogpu/{sample}/
  ├── {sample}_amr_abundance.tsv      # AMR + stress genes
  ├── {sample}_mutations.tsv          # FQ resistance
  ├── {sample}_timing.tsv             # Performance metrics
  ├── {sample}_stats.txt              # Human-readable summary
  └── {sample}_*_gpu_profile.log      # GPU profiling (if --profile)
```

---

### Step 6: Aggregate Results (5 minutes)

```bash
# Create aggregated timing files
cat results/traditional/*/`*_timing.tsv | head -1 > results/traditional_timing_all.tsv
cat results/traditional/*/`*_timing.tsv | grep -v "^sample_name" | grep "TOTAL" >> results/traditional_timing_all.tsv

cat results/biogpu/*/`*_timing.tsv | head -1 > results/biogpu_timing_all.tsv
cat results/biogpu/*/`*_timing.tsv | grep -v "^sample_name" | grep "TOTAL" >> results/biogpu_timing_all.tsv

# Count successful runs
echo "Traditional pipeline: $(cat results/traditional_timing_all.tsv | grep -v "sample_name" | wc -l) samples"
echo "BioGPU pipeline: $(cat results/biogpu_timing_all.tsv | grep -v "sample_name" | wc -l) samples"
```

---

### Step 7: Compare and Analyze (2-4 hours)

Create comparison analysis script or use Python/R:

```python
import pandas as pd
import numpy as np
from scipy.stats import spearmanr, pearsonr
import matplotlib.pyplot as plt
import seaborn as sns

# Load timing data
trad_timing = pd.read_csv('results/traditional_timing_all.tsv', sep='\t')
biogpu_timing = pd.read_csv('results/biogpu_timing_all.tsv', sep='\t')

# Calculate speedup
speedup = trad_timing['wall_time_sec'] / biogpu_timing['wall_time_sec']
print(f"Speedup statistics:")
print(f"  Mean: {speedup.mean():.2f}x")
print(f"  Median: {speedup.median():.2f}x")
print(f"  Min: {speedup.min():.2f}x")
print(f"  Max: {speedup.max():.2f}x")

# Load abundance data for correlation
# ... (detailed analysis of gene detection and quantification)

# Generate figures
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# 1. Wall time comparison
axes[0,0].bar(['Traditional', 'BioGPU'],
              [trad_timing['wall_time_sec'].mean(),
               biogpu_timing['wall_time_sec'].mean()])
axes[0,0].set_ylabel('Time (seconds)')
axes[0,0].set_title('Average Pipeline Time')

# 2. Speedup distribution
axes[0,1].hist(speedup, bins=20, edgecolor='black')
axes[0,1].set_xlabel('Speedup Factor')
axes[0,1].set_ylabel('Count')
axes[0,1].set_title(f'Speedup Distribution (mean={speedup.mean():.1f}x)')
axes[0,1].axvline(speedup.mean(), color='red', linestyle='--', label='Mean')
axes[0,1].legend()

# 3. Memory comparison
axes[1,0].bar(['Traditional', 'BioGPU'],
              [trad_timing['memory_gb'].mean(),
               biogpu_timing['memory_gb'].mean()])
axes[1,0].set_ylabel('Memory (GB)')
axes[1,0].set_title('Peak Memory Usage')

# 4. Sample-by-sample comparison
axes[1,1].scatter(trad_timing['wall_time_sec'], biogpu_timing['wall_time_sec'], alpha=0.6)
axes[1,1].plot([0, trad_timing['wall_time_sec'].max()],
               [0, trad_timing['wall_time_sec'].max()],
               'r--', label='Equal time')
axes[1,1].set_xlabel('Traditional time (s)')
axes[1,1].set_ylabel('BioGPU time (s)')
axes[1,1].set_title('Per-sample Time Comparison')
axes[1,1].legend()

plt.tight_layout()
plt.savefig('results/performance_comparison.pdf', dpi=300)
print("✓ Saved: results/performance_comparison.pdf")
```

---

## Expected Timeline

| Step | Time | Can Run Overnight? |
|------|------|--------------------|
| 1. Create environment | 5 min | No |
| 2. Setup databases | 15 min | No |
| 3. Verify samples | 1 min | No |
| 4. Traditional pipeline | 12-25 hrs | **Yes** |
| 5. BioGPU pipeline | 4-8 hrs | **Yes** |
| 6. Aggregate results | 5 min | No |
| 7. Analysis | 2-4 hrs | No |
| **Total** | **~1-2 days** | |

---

## Troubleshooting

### Traditional pipeline failing
```bash
# Check specific sample log
tail -100 logs/{sample}_traditional_*.log

# Check for common issues
# - Insufficient memory: increase available RAM
# - Disk space: check `df -h`
# - Missing files: verify FASTQ paths
```

### BioGPU pipeline failing
```bash
# Check GPU availability
nvidia-smi

# Check specific sample log
tail -100 logs/{sample}_biogpu_*.log

# Common issues:
# - GPU out of memory: reduce batch size or use different GPU
# - CUDA errors: check CUDA_VISIBLE_DEVICES
```

### Slow performance
```bash
# Check system load
htop

# Check I/O wait
iostat -x 5

# If I/O bound: files on network storage?
# If CPU bound: reduce thread count
# If Memory bound: process fewer samples in parallel
```

---

## Parallelization Options

### Multiple GPUs for BioGPU
```bash
# GPU 0
./scripts/05_batch_biogpu_with_timing.sh --gpu-id 0 &

# GPU 1 (if available)
./scripts/05_batch_biogpu_with_timing.sh --gpu-id 1 &
```

Note: Need to split test_samples.csv first

### Multiple Nodes for Traditional
Process different samples on different machines (shared filesystem required)

---

## Results Structure

After completion:
```
benchmark_biogpu/
├── data/
│   └── test_samples.csv                    # 50 samples
├── results/
│   ├── traditional/
│   │   ├── N01_1_4/                        # Per-sample results
│   │   ├── N02_2_4/
│   │   └── ...                             # 50 directories
│   ├── biogpu/
│   │   ├── N01_1_4/                        # Per-sample results
│   │   ├── N02_2_4/
│   │   └── ...                             # 50 directories
│   ├── traditional_timing_all.tsv          # Aggregated traditional timing
│   └── biogpu_timing_all.tsv               # Aggregated biogpu timing
└── logs/
    ├── batch_traditional_*.log
    └── batch_biogpu_*.log
```

---

## Success Criteria

✅ **50/50 samples** processed successfully through both pipelines
✅ **Timing data** captured for all samples
✅ **Gene abundance** data generated for all samples
✅ **Correlation > 0.90** for genes detected by both methods
✅ **Speedup documented** for reviewer response

---

## For Reviewer Response

After analysis, you'll have:
- **Performance comparison**: "BioGPU was X.Xx faster"
- **Validation**: "Correlation r=0.XX for shared genes"
- **Sensitivity**: "BioGPU detected Y additional genes"
- **Figures**: Performance and correlation plots
- **Methods**: Detailed benchmarking methods

---

## Quick Reference

**Start traditional pipeline:**
```bash
screen -S trad
conda activate benchmark_biogpu
cd /home/david/projects/benchmark_biogpu
./scripts/03_batch_traditional_pipeline.sh
# Ctrl+A, D to detach
```

**Start biogpu pipeline:**
```bash
screen -S biogpu
conda activate benchmark_biogpu
cd /home/david/projects/benchmark_biogpu
./scripts/05_batch_biogpu_with_timing.sh --gpu-id 0
# Ctrl+A, D to detach
```

**Check progress:**
```bash
screen -r trad     # Reattach traditional
screen -r biogpu   # Reattach biogpu
screen -ls         # List all sessions
```

**Quick status:**
```bash
echo "Traditional: $(ls results/traditional/*/`*_timing.tsv 2>/dev/null | wc -l) / 50"
echo "BioGPU: $(ls results/biogpu/*/`*_timing.tsv 2>/dev/null | wc -l) / 50"
```

---

## Ready to Run!

Everything is set up and ready:
- ✅ 50 samples selected
- ✅ Scripts with timing instrumentation
- ✅ Documentation complete
- ✅ Batch processing scripts ready

**Let's benchmark!** 🚀
