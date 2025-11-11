# BioGPU Pipeline Benchmarking Project

## Purpose

Compare the custom GPU-accelerated biogpu pipeline (translated search) against traditional alignment methods for antibiotic resistance gene detection. This benchmarking is critical for validating results used in the ZCH_UCMC manuscript reviewer response.

## Background

### BioGPU Pipeline (Novel Method)
- **Location**: `/home/david/projects/biogpu`
- **Method**: Translated search (6-frame translation) using GPU-accelerated alignment
- **Tools**: Custom `amr_detection` and `clean_resistance_pipeline` executables
- **Databases**:
  - AMR + Stress genes: `data/amr_stress_combined_db/` (DNA + protein)
  - FQ resistance mutations: `data/integrated_clean_db/` (nucleotide + protein)
- **Parameters**:
  - Min identity: 0.85 (85%)
  - Min coverage: 0.50 (50%)
- **Output**: Gene-level abundance in RPM (reads per million)

### Traditional Pipeline (Reference Method)
- **Method**: Nucleotide alignment + feature counting
- **Tools**:
  - bwa-mem OR bowtie2 (for read alignment)
  - featureCounts OR custom counting script (for quantification)
- **Databases**: Same reference sequences converted to nucleotide-only format
- **Output**: Gene-level read counts normalized to RPM

## Directory Structure

```
benchmark_biogpu/
├── scripts/           # Processing scripts for both pipelines
├── data/             # Test dataset (subset of NICU samples)
├── results/
│   ├── biogpu/       # Results from biogpu pipeline
│   └── traditional/  # Results from bwa-mem/bowtie2 + featureCounts
├── logs/             # Processing logs
├── databases/        # Prepared reference databases
│   ├── biogpu/      # Original biogpu format
│   └── traditional/ # BWA/bowtie2 indexed references
├── docs/             # Documentation and analysis reports
└── README.md

```

## Workflow

### Phase 1: Setup
1. ✓ Create directory structure
2. ⃞ Identify test samples (subset of NICU data for speed)
3. ⃞ Link/copy biogpu reference databases
4. ⃞ Prepare traditional alignment databases (index with bwa/bowtie2)
5. ⃞ Set up conda/mamba environment with required tools

### Phase 2: Run Pipelines
1. ⃞ Extract biogpu results for test samples (already processed)
2. ⃞ Run traditional alignment pipeline on same samples
3. ⃞ Normalize both outputs to same scale (RPM)

### Phase 3: Comparison
1. ⃞ Compare gene detection sensitivity (genes detected by each method)
2. ⃞ Compare abundance correlations (Spearman/Pearson)
3. ⃞ Identify discrepancies and investigate causes
4. ⃞ Runtime and resource usage comparison
5. ⃞ Generate summary report and figures

## Expected Outcomes

- **High correlation** (r > 0.90): Validates biogpu method
- **Sensitivity comparison**: Which method detects more genes?
- **False positive/negative analysis**: Investigate major discrepancies
- **Performance metrics**: Speed and resource usage

## Test Dataset

TBD - Will use subset of NICU samples (e.g., 10-20 samples representing different body sites and timepoints)

## Environment Setup

Tools needed:
- bwa-mem (or bowtie2)
- samtools
- featureCounts (from Subread package) OR bedtools
- Python with pandas, numpy, scipy, matplotlib, seaborn

## Key Questions to Answer

1. Do both methods identify similar sets of resistance genes?
2. Are abundance estimates correlated?
3. Which method is more sensitive for low-abundance genes?
4. Are there systematic biases (e.g., GC content, gene length)?
5. What is the computational cost difference?

## References

- BioGPU pipeline script: `/home/david/projects/biogpu/batch_process_nicu_samples_stress.sh`
- Original NICU analysis: `/home/david/projects/ZCH_UCMC_Manuscript/nicu_resistome_analysis/`
- Output location (biogpu): `/fastpool/analysis/nicu_amr_stress/`
