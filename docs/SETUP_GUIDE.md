# Benchmark Project Setup Guide

## Step 1: Create Conda Environment

```bash
cd /home/david/projects/benchmark_biogpu

# Create environment from YAML
mamba env create -f environment.yml

# Or if using conda
conda env create -f environment.yml

# Activate
conda activate benchmark_biogpu
```

## Step 2: Identify Test Samples

We need a representative subset of NICU samples for benchmarking. Criteria:
- Include all body sites (axilla, groin, stool)
- Include both timepoints (Week 1, Week 3)
- Include both locations (UCMC, ZCH)
- Include samples with varying AMR burden (low, medium, high)
- Total: ~10-20 samples (enough for statistics, fast enough to iterate)

Sample selection strategy:
1. Review NICU analysis metadata: `/home/david/projects/ZCH_UCMC_Manuscript/nicu_resistome_analysis/results/qc/master_metadata_with_qc.tsv`
2. Select stratified random samples or specific representative samples
3. Create sample list CSV

## Step 3: Link Reference Databases

### BioGPU Databases (already exist)
```bash
cd databases/biogpu

# Link AMR + Stress database
ln -s /home/david/projects/biogpu/data/amr_stress_combined_db ./amr_stress_combined_db

# Link FQ resistance database
ln -s /home/david/projects/biogpu/data/integrated_clean_db ./integrated_clean_db
```

### Traditional Pipeline Databases

**Primary Method: DIAMOND** (matching biogpu's translated search):

```bash
cd databases/traditional

# Copy protein sequences for DIAMOND
cp /home/david/projects/biogpu/data/amr_stress_combined_db/protein.fasta amr_stress_protein.fasta

# Build DIAMOND index
diamond makedb --in amr_stress_protein.fasta --db amr_stress_protein
# Creates: amr_stress_protein.dmnd
```

**Alternative Method: Bowtie2** (nucleotide alignment - optional):

```bash
# Copy DNA sequences for bowtie2
cp /home/david/projects/biogpu/data/amr_stress_combined_db/dna.fasta amr_stress_dna.fasta

# Build Bowtie2 index
bowtie2-build amr_stress_dna.fasta amr_stress_bt2
# Creates: amr_stress_bt2.*.bt2 files
```

**Note**: We now use DIAMOND blastx as the primary method because it uses translated search (protein space) just like biogpu, enabling fair comparison. Bowtie2 (nucleotide alignment) is retained as an optional alternative but is less sensitive for divergent genes.

## Step 4: Extract BioGPU Results for Test Samples

BioGPU results are in: `/fastpool/analysis/nicu_amr_stress/`

```bash
# Create links to results
cd results/biogpu

# Link to AMR+Stress results
ln -s /fastpool/analysis/nicu_amr_stress/amr_stress ./amr_stress_results

# Link to FQ resistance results
ln -s /fastpool/analysis/nicu_amr_stress/resistance ./fq_resistance_results
```

For each test sample, extract:
- `{sample_name}_amr_abundance.tsv` - Gene abundances in RPM
- Gene classifications (AMR vs Stress)

## Step 5: Locate FASTQ Files

FASTQ files location will be in the biogpu sample CSV:
```bash
cat /home/david/projects/biogpu/nicu_sample_key_fixed.csv
```

This CSV has columns: `sample_name,fastq_path,r1_file,r2_file`

## Step 6: Run DIAMOND Pipeline (Primary Method)

Script: `scripts/02_run_diamond_pipeline.sh` (already created)

Pipeline steps:
1. Run DIAMOND blastx on R1 and R2 separately (translated search)
2. Combine results and count reads per gene
3. Calculate abundance metrics (read count, RPM, TPM, coverage %)
4. Capture comprehensive timing data

Parameters:
- Identity: 85% (matching biogpu)
- Query coverage: 50% (matching biogpu)
- Mode: --sensitive

## Step 6b: Run Bowtie2 Pipeline (Optional Alternative)

Script: `scripts/02_run_traditional_pipeline.sh` (legacy method)

Pipeline steps:
1. Align reads to reference with bowtie2 (nucleotide alignment)
2. Convert SAM to BAM, sort, index
3. Count reads per gene with htseq-count
4. Calculate coverage with bedtools
5. Normalize to RPM and TPM

## Step 7: Analysis and Comparison

Python script to create: `scripts/compare_pipelines.py`

Analyses:
1. **Gene Detection Comparison**
   - Venn diagram: genes detected by each method
   - Sensitivity/specificity analysis

2. **Abundance Correlation**
   - Scatter plots: biogpu RPM vs traditional RPM
   - Spearman and Pearson correlations
   - Bland-Altman plots for agreement

3. **Discrepancy Analysis**
   - Genes detected only by biogpu (translation advantage?)
   - Genes detected only by traditional (false positives in biogpu?)
   - Large abundance differences (>2-fold) - investigate why

4. **Performance Metrics**
   - Runtime comparison
   - Memory usage
   - Disk I/O

5. **Summary Report**
   - Tables and figures for manuscript supplement
   - Recommendations for method validation

## Expected Challenges

1. **Search space matching** (SOLVED):
   - **Problem**: BioGPU uses protein space (6-frame translation), bowtie2 uses nucleotide space
   - **Solution**: Use DIAMOND blastx (translated search) with same parameters as biogpu
   - Now both methods search in protein space with 85% identity and 50% coverage

2. **Gene length normalization**:
   - Both DIAMOND and biogpu calculate RPM and TPM
   - Need to verify RPM denominator (reads vs read pairs)

3. **Multi-mapping reads**: Both methods use similar strategies
   - BioGPU: best hit per read
   - DIAMOND: top 1 hit per read (default)
   - Should produce comparable results

4. **Performance comparison**:
   - DIAMOND: CPU-based translated search (~4-5 min/sample)
   - BioGPU: GPU-accelerated translated search (~5-10 min/sample)
   - Capture comprehensive timing for both

## Success Criteria

- **Primary**: Correlation > 0.90 for genes detected by both methods
- **Secondary**: Similar gene detection (>80% overlap)
- **Acceptable**: Can explain major discrepancies with biological/technical rationale

## Next Steps After User Provides Script

Once you provide your existing traditional pipeline script:
1. Adapt it for this benchmark project
2. Select test samples
3. Run both pipelines
4. Generate comparison analysis
5. Create summary report for reviewer response
