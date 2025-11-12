# Processing All Remaining NICU Samples Through DIAMOND

## Summary

- **Total NICU samples**: 771
- **Already processed**: 51 (50 benchmark + 1 extra)
- **Remaining samples**: 594 with valid FASTQ files
- **Missing files**: 127 samples (no FASTQ files found)

## Resource Configuration

**DIAMOND uses 32 threads per sample**

With **128 cores available**, you can run:
- **4 samples in parallel** (128 / 32 = 4)
- **Optimal configuration**: MAX_JOBS=4

## Estimated Runtime

Based on benchmarking (avg 327 seconds per sample):

| Configuration | Time to Complete |
|---------------|------------------|
| Serial (1 at a time) | ~54 hours |
| Parallel (2 jobs) | ~27 hours |
| **Parallel (4 jobs)** | **~13.5 hours** ⭐ RECOMMENDED |
| Parallel (8 jobs) | ~7 hours (may overload system) |

## How to Run

### Option 1: Parallel Processing (RECOMMENDED)

**Best for 128-core system:**

```bash
cd /home/david/projects/benchmark_biogpu

# Run with 4 parallel jobs (optimal for 128 cores)
MAX_JOBS=4 ./scripts/process_remaining_samples_parallel.sh
```

**Features:**
- ✅ Progress bar showing completion
- ✅ Automatic skipping of already-processed samples
- ✅ Individual logs per sample
- ✅ Job log for tracking failures
- ✅ ~13.5 hours estimated runtime

**Adjust parallelism:**
```bash
# More conservative (2 parallel)
MAX_JOBS=2 ./scripts/process_remaining_samples_parallel.sh

# More aggressive (8 parallel) - uses only 16 threads/sample
MAX_JOBS=8 ./scripts/process_remaining_samples_parallel.sh
```

### Option 2: Serial Processing

**If you don't have GNU parallel installed or prefer simple sequential processing:**

```bash
./scripts/process_remaining_samples.sh
```

**Features:**
- ✅ Simple, reliable
- ✅ Automatic skipping of already-processed samples
- ❌ Slower (~54 hours)

## Running in Screen (RECOMMENDED)

Since this will take ~13 hours, run in a screen session:

```bash
# Start screen
screen -S diamond_all

# Activate environment
conda activate benchmark_biogpu
cd /home/david/projects/benchmark_biogpu

# Run processing
MAX_JOBS=4 ./scripts/process_remaining_samples_parallel.sh

# Detach: Ctrl+A, then D
```

**Check progress later:**
```bash
# Reattach to screen
screen -r diamond_all

# Or check progress without attaching
tail -f logs/batch_remaining_parallel_*.log

# Count completed samples
ls -d results/traditional/*/ | wc -l
# Should show: 51 (current) + processed samples
```

## Monitoring Progress

### Real-time monitoring:
```bash
# Watch the batch log
tail -f logs/batch_remaining_parallel_*.log

# Count completed samples
watch -n 60 'ls -d results/traditional/*/ | wc -l'

# Check system load
htop
```

### Check specific sample:
```bash
# Check if sample completed
ls results/traditional/N01_1_2/N01_1_2_abundance.tsv

# View sample log
tail -50 logs/N01_1_2_diamond_*.log
```

### Summary statistics:
```bash
# Total samples processed
ls -d results/traditional/*/ | wc -l

# Failed samples (if any)
grep FAILED logs/batch_remaining_parallel_*.log

# Average processing time
grep "SUCCESS:" logs/batch_remaining_parallel_*.log | \
    sed 's/.*(//;s/s)//' | awk '{sum+=$1; n++} END {print sum/n "s average"}'
```

## After Completion

Once processing is complete, all results will be in the same directory structure:

```
results/traditional/
├── N01_1_2/
│   ├── N01_1_2_abundance.tsv       # Gene counts, RPM, TPM, coverage
│   ├── N01_1_2_diamond.tsv         # Raw DIAMOND output
│   ├── N01_1_2_timing.tsv          # Performance metrics
│   └── N01_1_2_stats.txt           # Human-readable summary
├── N01_1_3/
├── ... (all 645 samples)
└── ZJH_N59_2_3/
```

## Aggregate All Results

After processing completes, aggregate into analysis-ready files:

```bash
# Create aggregated abundance matrix
python scripts/aggregate_diamond_results.py

# This will create:
# - results/diamond_abundance_matrix.tsv
# - results/diamond_rpm_matrix.tsv
# - results/diamond_tpm_matrix.tsv
# - results/diamond_timing_summary.tsv
```

(I can create this aggregation script if you need it)

## Troubleshooting

### If processing fails:

**Check failed samples:**
```bash
grep FAILED logs/batch_remaining_parallel_*.log
```

**Re-run only failed samples:**
The script automatically skips completed samples, so just re-run:
```bash
MAX_JOBS=4 ./scripts/process_remaining_samples_parallel.sh
```

### If you need to stop and resume:

The scripts check for existing `*_abundance.tsv` files and skip already-processed samples automatically. You can safely:
1. Stop the script (Ctrl+C)
2. Resume later - it will continue from where it left off

### If memory issues occur:

Reduce parallel jobs:
```bash
MAX_JOBS=2 ./scripts/process_remaining_samples_parallel.sh
```

## System Requirements

- **CPU**: 128 cores (4 samples × 32 threads/sample)
- **Memory**: ~48 GB (4 samples × ~12 GB/sample)
- **Disk**: ~100 GB for all results
- **Time**: ~13.5 hours for 594 samples

## Quick Start (TL;DR)

```bash
# 1. Start screen session
screen -S diamond_all

# 2. Activate environment
conda activate benchmark_biogpu
cd /home/david/projects/benchmark_biogpu

# 3. Run parallel processing
MAX_JOBS=4 ./scripts/process_remaining_samples_parallel.sh

# 4. Detach (Ctrl+A, then D)

# 5. Check progress later
screen -r diamond_all
# or
tail -f logs/batch_remaining_parallel_*.log
```

## Notes

- Script automatically skips already-processed samples
- Each sample takes ~5-10 minutes (avg 327s from benchmarking)
- Results are saved incrementally (no data loss if interrupted)
- All logs saved to `logs/` directory for troubleshooting
- Progress bar shows real-time completion status
