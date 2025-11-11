# Timing and Performance Metrics

## Overview

Both pipelines now capture comprehensive timing and performance metrics for direct comparison:

### Metrics Captured

1. **Wall Clock Time**: Total elapsed time for each step and pipeline
2. **CPU Time**: User + system CPU time (from `/usr/bin/time`)
3. **GPU Time**: GPU execution time (biogpu only, requires nvprof)
4. **Memory Usage**: Peak resident set size (RAM)
5. **Threads**: Number of threads used per step

## Traditional Pipeline Timing

### Steps Measured

| Step | Tool | Description | Typical % of Total |
|------|------|-------------|-------------------|
| 1 | bowtie2 | Nucleotide alignment | 40-60% |
| 2 | samtools | SAM→BAM conversion + sorting | 10-20% |
| 3 | htseq-count | Read quantification | 20-30% |
| 4 | bedtools | Coverage calculation | 5-10% |
| 5 | python | RPM/TPM calculation | <5% |

### Output Files

**Human-readable**: `{sample}_stats.txt`
```
Performance Metrics (Traditional Pipeline):
  Total pipeline time: 1234s (20.57m)

  Step-by-step timing (wall clock time):
    1. Bowtie2 alignment:    567s (45.9%)
    2. SAM/BAM processing:   234s (19.0%)
    3. htseq-count:          345s (27.9%)
    4. bedtools coverage:    78s (6.3%)
    5. RPM/TPM calculation:  10s (0.8%)

  CPU time (user):
    - Bowtie2:     5432s
    - samtools:    1234s
    - htseq-count: 345s
    - bedtools:    78s

  Peak memory usage:
    - Bowtie2:     4.5GB
    - samtools:    2.1GB
    - htseq-count: 1.2GB
    - bedtools:    0.8GB

  Threads used: 12
```

**Machine-readable**: `{sample}_timing.tsv`
```
sample_name  step        wall_time_sec  cpu_time_sec  memory_kb  memory_gb  threads  percent_of_total
N01_1_2      bowtie2     567           5432          4718592    4.500      12       45.94
N01_1_2      samtools    234           1234          2202009    2.100      12       18.96
N01_1_2      htseq-count 345           345           1258291    1.200      1        27.96
N01_1_2      bedtools    78            78            838861     0.800      1        6.32
N01_1_2      python      10            NA            NA         NA         1        0.81
N01_1_2      TOTAL       1234          NA            NA         NA         12       100.00
```

### Detailed Timing Files

Each step also generates a detailed `/usr/bin/time` output:
- `{sample}_bowtie2_time.txt`
- `{sample}_samtools_time.txt`
- `{sample}_htseq_time.txt`
- `{sample}_bedtools_time.txt`

These contain additional metrics like:
- Page faults
- Context switches
- I/O operations
- System calls

## BioGPU Pipeline Timing

### Steps Measured

| Step | Tool | Description | Typical % of Total |
|------|------|-------------|-------------------|
| 1 | amr_detection | AMR + stress gene detection | 60-70% |
| 2 | clean_resistance_pipeline | FQ mutation detection | 30-40% |

### Output Files

**Human-readable**: `{sample}_stats.txt`
```
Performance Metrics (BioGPU Pipeline):
  Total pipeline time: 456s (7.60m)

  Step-by-step timing (wall clock time):
    1. AMR + Stress detection: 312s (68.4%)
    2. FQ resistance detection: 144s (31.6%)

  CPU time (user + system):
    - AMR:  234s + 12s = 246s
    - FQ:   98s + 5s = 103s

  GPU time:
    - Total: 145s + 67s

  Peak memory usage:
    - AMR:  8.5GB
    - FQ:   6.2GB

Pipeline Parameters:
  - Min identity: 0.85 (85%)
  - Min coverage: 0.50 (50%)
  - Search method: Translated (6-frame)
  - GPU accelerated: Yes
```

**Machine-readable**: `{sample}_timing.tsv`
```
sample_name  step          wall_time_sec  cpu_time_sec  gpu_time  memory_kb  memory_gb  device  percent_of_total
N01_1_2      amr_detection 312           234           145s      8912896    8.500      GPU0    68.42
N01_1_2      fq_resistance 144           98            67s       6501171    6.200      GPU0    31.58
N01_1_2      TOTAL         456           NA            212s      NA         NA         GPU0    100.00
```

### GPU Profiling (Optional)

Enable with `--profile` flag to get detailed GPU metrics:
```bash
./scripts/04_run_biogpu_with_timing.sh --sample N01_1_2 --r1 ... --r2 ... --profile
```

Generates:
- `{sample}_amr_gpu_profile.log` - AMR detection GPU profile
- `{sample}_fq_gpu_profile.log` - FQ resistance GPU profile

Contains metrics like:
- Kernel execution time
- Memory transfers (host↔device)
- GPU utilization
- Memory bandwidth

## Comparing Pipelines

### Performance Comparison Script

Create a script to aggregate and compare timing across all samples:

```python
import pandas as pd
import matplotlib.pyplot as plt

# Load timing data
traditional = pd.read_csv('results/traditional_timing_summary.tsv', sep='\t')
biogpu = pd.read_csv('results/biogpu_timing_summary.tsv', sep='\t')

# Calculate speedup
speedup = traditional['total_time'] / biogpu['total_time']

# Plot comparison
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Wall time comparison
axes[0,0].bar(['Traditional', 'BioGPU'],
              [traditional['total_time'].mean(), biogpu['total_time'].mean()])
axes[0,0].set_ylabel('Time (seconds)')
axes[0,0].set_title('Average Pipeline Time')

# Speedup distribution
axes[0,1].hist(speedup, bins=20)
axes[0,1].set_xlabel('Speedup Factor')
axes[0,1].set_title(f'BioGPU Speedup (mean: {speedup.mean():.1f}x)')

# Memory comparison
axes[1,0].bar(['Traditional', 'BioGPU'],
              [traditional['max_memory_gb'].mean(), biogpu['max_memory_gb'].mean()])
axes[1,0].set_ylabel('Memory (GB)')
axes[1,0].set_title('Peak Memory Usage')

# Time per step (traditional)
traditional_steps = traditional[['bowtie2', 'samtools', 'htseq', 'bedtools']].mean()
axes[1,1].pie(traditional_steps, labels=traditional_steps.index, autopct='%1.1f%%')
axes[1,1].set_title('Traditional Pipeline Time Breakdown')

plt.tight_layout()
plt.savefig('figures/performance_comparison.pdf')
```

### Expected Results

**Hypothesis**: BioGPU should be significantly faster due to:
1. GPU acceleration
2. Translated search (more efficient than nucleotide alignment for divergent sequences)
3. Optimized custom pipeline

**Typical speedups**:
- Best case: 5-10x faster (large files, many AMR genes)
- Worst case: 2-3x faster (small files, overhead dominates)
- Average: 3-5x faster

**Trade-offs**:
- BioGPU: Higher peak memory (8-10GB), requires GPU
- Traditional: Lower memory (4-5GB), CPU-only

## Analysis Considerations

### Fair Comparison Requirements

1. **Same input data**: Identical FASTQ files
2. **Same database**: Same gene sequences (different formats)
3. **Similar coverage**: Apply 50% filter to traditional results
4. **Thread count**: Use all available cores for traditional (biogpu uses GPU)

### Factors Affecting Timing

**File size**:
- Larger files → more alignment time
- BioGPU advantage increases with file size

**AMR gene abundance**:
- More resistance genes → more alignment time
- Affects both pipelines similarly

**System factors**:
- CPU: More cores → faster traditional pipeline
- GPU: Better GPU → faster biogpu
- Memory: Insufficient RAM → swapping → much slower
- Disk I/O: Slow disk → slower both pipelines

### Reporting Performance

For reviewer response, report:

1. **Average speedup**: "BioGPU was X.Xx faster on average"
2. **Time savings**: "Reduced analysis time from Y hours to Z minutes"
3. **Scalability**: "Speedup increases with dataset size"
4. **Resource usage**: "Similar memory footprint, uses GPU acceleration"

Example text:
> "To assess computational efficiency, we compared the GPU-accelerated biogpu pipeline against a traditional CPU-based approach (bowtie2 + htseq-count) on 20 representative samples. The biogpu pipeline was 4.2× faster on average (mean time: 7.6min vs 32.1min), with speedup increasing for larger files (up to 8.5× for high-coverage samples). Both pipelines required similar peak memory (8-10GB), though biogpu utilized GPU acceleration while the traditional pipeline used 12 CPU cores."

## Timing Data Files

### Per-Sample Files

**Traditional Pipeline**:
- `results/traditional/{sample}/{sample}_timing.tsv`
- `results/traditional/{sample}/{sample}_stats.txt`
- `results/traditional/{sample}/{sample}_*_time.txt` (detailed)

**BioGPU Pipeline**:
- `results/biogpu/{sample}/{sample}_timing.tsv`
- `results/biogpu/{sample}/{sample}_stats.txt`
- `results/biogpu/{sample}/{sample}_*_gpu_profile.log` (if profiled)

### Aggregated Summary Files

Create these by combining per-sample timing files:

```bash
# Traditional pipeline summary
cat results/traditional/*/\*_timing.tsv | grep -v "^sample_name" | \
    sort -u > results/traditional_timing_summary.tsv

# BioGPU pipeline summary
cat results/biogpu/*/\*_timing.tsv | grep -v "^sample_name" | \
    sort -u > results/biogpu_timing_summary.tsv
```

## Next Steps

1. ✓ Scripts created with timing instrumentation
2. ☐ Run traditional pipeline on test samples
3. ☐ Re-run biogpu on same samples with timing
4. ☐ Create comparison script
5. ☐ Generate performance figures
6. ☐ Add to reviewer response

See `scripts/02_run_traditional_pipeline.sh` and `scripts/04_run_biogpu_with_timing.sh` for implementation.
