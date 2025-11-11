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

### Traditional Pipeline Databases (need to create)

We'll use the same reference sequences but in nucleotide format for bwa-mem/bowtie2:

```bash
cd databases/traditional

# Copy protein sequences (will need to extract DNA sequences from biogpu DBs)
# These are FASTA files that need indexing

# For AMR + Stress genes
# Extract DNA sequences from: biogpu/data/amr_stress_combined_db/dna.fasta
cp /home/david/projects/biogpu/data/amr_stress_combined_db/dna.fasta amr_stress_dna.fasta

# Build BWA index
bwa index amr_stress_dna.fasta

# Build Bowtie2 index (optional - for comparison)
bowtie2-build amr_stress_dna.fasta amr_stress_dna
```

**Note**: The biogpu pipeline uses translated search (protein space), so it can detect divergent sequences. Traditional nucleotide alignment may be less sensitive for divergent genes.

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

## Step 6: Run Traditional Pipeline

Script to create: `scripts/run_traditional_pipeline.sh`

Pipeline steps:
1. Align reads to reference with bwa-mem (or bowtie2)
2. Convert SAM to BAM, sort, index
3. Count reads per gene with featureCounts
4. Normalize to RPM
5. Compare to biogpu results

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

1. **Different search spaces**: BioGPU uses protein space (6-frame translation), traditional uses nucleotide space
   - BioGPU may detect more divergent sequences
   - Fair comparison requires understanding this fundamental difference

2. **Gene length normalization**: featureCounts normalizes differently than biogpu
   - Need to ensure comparable RPM calculations

3. **Multi-mapping reads**: How each method handles reads mapping to multiple genes
   - BioGPU: best hit strategy
   - featureCounts: various options (unique, multi-overlap)

4. **Coverage thresholds**: BioGPU uses 50% coverage minimum
   - May need to apply similar filter to traditional results

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
