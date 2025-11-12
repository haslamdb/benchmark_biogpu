# Ready to Process All 594 Remaining NICU Samples

## What's Been Created

✅ **Sample identification script** - Finds remaining unprocessed samples
✅ **Parallel processing script** - Processes multiple samples simultaneously
✅ **Serial processing script** - Backup option for sequential processing
✅ **Sample list** - `data/remaining_samples.csv` (594 samples with valid FASTQs)

## Your System

- **CPU**: 128 cores available
- **DIAMOND per sample**: 32 threads
- **Optimal parallel jobs**: 4 samples at once (128 / 32 = 4)
- **GNU parallel**: ✅ Already installed

## Estimated Time

| Method | Jobs | Runtime |
|--------|------|---------|
| Serial | 1 | ~54 hours |
| Parallel | 2 | ~27 hours |
| **Parallel** | **4** | **~13.5 hours** ⭐ |

## Quick Start

```bash
# Start a screen session (so it keeps running)
screen -S diamond_all

# Navigate to project
conda activate benchmark_biogpu
cd /home/david/projects/benchmark_biogpu

# Run with 4 parallel jobs (optimal for your 128 cores)
MAX_JOBS=4 ./scripts/process_remaining_samples_parallel.sh

# Detach from screen: Ctrl+A, then D
```

## Monitor Progress

```bash
# Reattach to screen
screen -r diamond_all

# Or watch log without attaching
tail -f logs/batch_remaining_parallel_*.log

# Count completed samples
ls -d results/traditional/*/ | wc -l
# Currently: 51
# Target: 645 (51 + 594)
```

## What Happens

1. Script reads `data/remaining_samples.csv` (594 samples)
2. Skips any samples already in `results/traditional/`
3. Runs DIAMOND on 4 samples at a time
4. Each sample outputs to: `results/traditional/{sample_name}/`
   - `{sample}_abundance.tsv` - Gene-level results (RPM, TPM, coverage)
   - `{sample}_diamond.tsv` - Raw DIAMOND output
   - `{sample}_timing.tsv` - Performance metrics
   - `{sample}_stats.txt` - Summary

## Files Created

```
scripts/
├── identify_remaining_samples.py           # Find unprocessed samples
├── process_remaining_samples.sh            # Serial processing
└── process_remaining_samples_parallel.sh   # Parallel processing ⭐

data/
└── remaining_samples.csv                   # 594 samples to process

docs/
├── PROCESS_ALL_SAMPLES.md                  # Detailed guide
└── READY_TO_PROCESS.md                     # This file
```

## Need to Adjust?

### Run fewer parallel jobs (more conservative):
```bash
MAX_JOBS=2 ./scripts/process_remaining_samples_parallel.sh
```

### Run more parallel jobs (aggressive - may overload):
```bash
# This will use only 16 threads per sample
MAX_JOBS=8 ./scripts/process_remaining_samples_parallel.sh
```

### Run without parallel (simple sequential):
```bash
./scripts/process_remaining_samples.sh
```

## The Command You Want

Based on your 128 cores and the fact that DIAMOND uses 32 threads per sample:

```bash
screen -S diamond_all
conda activate benchmark_biogpu
cd /home/david/projects/benchmark_biogpu
MAX_JOBS=4 ./scripts/process_remaining_samples_parallel.sh
# Ctrl+A, D to detach
```

This will:
- Run 4 samples simultaneously
- Take ~13.5 hours
- Show progress bar
- Resume automatically if interrupted
- Save all results to `results/traditional/`

Ready to go! 🚀
