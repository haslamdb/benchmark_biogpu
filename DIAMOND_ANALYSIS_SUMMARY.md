# DIAMOND Resistome Analysis Summary

**Project:** Benchmark BioGPU - DIAMOND Traditional Pipeline Analysis
**Date:** November 2025
**Samples:** 641 NICU infant samples (360 UCMC, 281 ZCH)
**Database:** AMR + Stress genes (14,280 genes)

---

## Table of Contents
1. [Pipeline Overview](#pipeline-overview)
2. [Scripts Used](#scripts-used)
3. [Data Processing](#data-processing)
4. [Downstream Analyses](#downstream-analyses)
5. [Output Directories](#output-directories)
6. [Key Findings](#key-findings)

---

## Pipeline Overview

### DIAMOND Alignment Pipeline
- **Tool:** DIAMOND blastx (protein alignment)
- **Database:** `/home/david/projects/benchmark_biogpu/databases/traditional/amr_stress_protein.dmnd`
- **Input:** Paired-end FASTQ files from NICU infant stool, axilla, and groin samples
- **Samples processed:** 642 (641 with complete metadata)
  - UCMC: 360 samples
  - ZCH: 281 samples (note: original FASTQs named `ZJH_N1_*` to `ZJH_N9_*` without zero-padding, output directories have zero-padding `ZJH_N01_*` to `ZJH_N09_*`)

---

## Scripts Used

### 1. Database Setup
**Script:** `scripts/01_setup_databases.sh`
- Created DIAMOND protein database from BioGPU AMR+stress combined database
- Database: `databases/traditional/amr_stress_protein.dmnd`
- Contains ~14,280 genes (9,257 AMR genes + 5,023 stress genes)

### 2. DIAMOND Alignment
**Script:** `scripts/02_run_diamond_pipeline.sh`
- Runs DIAMOND blastx alignment for individual samples
- Calculates abundance metrics (RPM, TPM, coverage)
- Output: Per-sample abundance tables

**Batch Script:** `scripts/03_batch_traditional_pipeline.sh`
- Processes multiple samples from CSV list
- Sample list: `data/test_samples.csv`

### 3. Data Preparation & Integration

#### a. Gene Annotation Mapping
**Script:** `scripts/create_gene_annotation_mapping.py`
- Extracted gene_family and drug_class annotations from BioGPU results
- Created unified annotation file
- **Output:** `data/gene_annotations.tsv` (351 genes with annotations)

#### b. Combine DIAMOND Results
**Script:** `scripts/combine_diamond_results.py`
- Merged all DIAMOND abundance files (642 samples)
- Added gene annotations (gene_family, drug_class)
- Filtered to genes with read counts > 0
- **Output:** `data/diamond_amr_combined.tsv` (1.68M entries, 6,408 unique genes)

#### c. Metadata Creation and Correction
**Script:** `scripts/create_metadata_from_samples.py`
- Parsed sample names to extract metadata (Location, SampleType, Week)
- Fixed naming: ZCH samples in metadata lack `ZJH_` prefix (added during processing)
- **Output:** `data/diamond_metadata.tsv`

**Metadata Correction:**
- Fixed 36 ZCH samples (+ 27 UCMC) with empty `sample_name` fields in original metadata
- Reconstructed sample names from PatientID, SampleType, SampleCollectionWeek
- **Output:** `data/master_metadata_fixed.tsv` (all 669 rows now have sample_name)
- **Output:** `data/diamond_metadata_updated.tsv` (641 samples with complete metadata)

---

## Downstream Analyses

### 1. Exploratory Analysis
**Script:** `scripts/analyze_diamond_resistome.py`

**Analyses performed:**
- Created sample × gene RPM matrix (642 samples × 6,408 genes)
- Calculated diversity metrics (Shannon, Simpson, richness)
- Principal Component Analysis (PCA) stratified by body site and location
- Hierarchical clustering (Ward linkage)

**Key outputs:**
- `results/analysis/diamond_resistome/amr_rpm_matrix.tsv`
- `results/analysis/diamond_resistome/diversity_metrics.tsv`
- `results/analysis/diamond_resistome/pca_coordinates.tsv`
- `figures/diamond_resistome/pca_by_bodysite.pdf`
- `figures/diamond_resistome/diversity_metrics.pdf`

### 2. Differential Abundance Analysis
**Script:** `scripts/differential_abundance_analysis.py`

**Comparisons:**
1. **Location:** UCMC vs ZCH
2. **Body Site:** Axilla vs Groin vs Stool (+ pairwise)
3. **Postnatal Age:** Week 1 vs Week 3

**Statistical methods:**
- Mann-Whitney U test (2 groups)
- Kruskal-Wallis H test (3+ groups)
- FDR correction (Benjamini-Hochberg)
- Filtering: Genes present (RPM > 1.0) in ≥ 5 samples

**Key outputs:**
- `results/analysis/differential_abundance/location_ucmc_vs_zch.tsv` (2,666 significant genes)
- `results/analysis/differential_abundance/bodysite_*.tsv` (3,399 significant genes overall)
- `results/analysis/differential_abundance/week_week1_vs_week3.tsv` (3,262 significant genes)

### 3. Antibiotic-Resistance Correlations
**Script:** `scripts/antibiotic_resistance_correlations.py`

**Analysis:**
- Correlated antibiotic exposure with resistance gene abundance
- **Cumulative exposure model:**
  - Week 1 samples: `_w1` antibiotics only
  - Week 3 samples: `_w1 + _w2` (cumulative exposure)
- Merged with clinical metadata (antibiotic exposure data)
- **Samples analyzed:** 641 (360 UCMC + 281 ZCH, all with antibiotic data)

**Key outputs:**
- `results/analysis/antibiotic_correlations/total_abx_amr_correlation_by_site.tsv`
- `results/analysis/antibiotic_correlations/specific_antibiotic_correlations.tsv`
- `figures/antibiotic_correlations/antibiotic_amr_correlation.pdf`
- `figures/antibiotic_correlations/antibiotic_exposure_amr_boxplot.pdf`

---

## Output Directories

### Raw DIAMOND Results
```
results/traditional/
├── N01_1_2/
│   ├── N01_1_2_abundance.tsv          # Gene abundance table
│   ├── N01_1_2_diamond.tsv            # Raw DIAMOND output
│   ├── N01_1_2_stats.txt              # Alignment statistics
│   └── N01_1_2_timing.tsv             # Timing information
├── ZJH_N01_1_2/                       # ZCH samples (with ZJH_ prefix)
│   └── ...
└── [640 more sample directories]
```

### Processed Data
```
data/
├── diamond_amr_combined.tsv           # All samples combined (1.68M entries)
├── diamond_metadata.tsv               # Original metadata (642 samples)
├── diamond_metadata_updated.tsv       # Corrected metadata (641 samples)
├── master_metadata_fixed.tsv          # Fixed ZCH_UCMC metadata (669 samples)
├── gene_annotations.tsv               # Gene family & drug class annotations
└── originally_missing_samples_antibiotic_status.tsv
```

### Analysis Results
```
results/analysis/
├── diamond_resistome/
│   ├── amr_rpm_matrix.tsv             # 642 × 6,408 gene matrix
│   ├── diversity_metrics.tsv          # Shannon, Simpson, richness per sample
│   ├── pca_coordinates.tsv            # First 10 PCs
│   └── pca_variance_explained.tsv
├── differential_abundance/
│   ├── location_ucmc_vs_zch.tsv       # 2,666 significant genes
│   ├── bodysite_comparison.tsv        # 3,399 significant genes
│   ├── bodysite_axilla_vs_groin.tsv   # 1,755 significant genes
│   ├── bodysite_axilla_vs_stool.tsv   # 3,418 significant genes
│   ├── bodysite_groin_vs_stool.tsv    # 2,667 significant genes
│   └── week_week1_vs_week3.tsv        # 3,262 significant genes
└── antibiotic_correlations/
    ├── total_abx_amr_correlation_by_site.tsv
    └── specific_antibiotic_correlations.tsv
```

### Figures
```
figures/
├── diamond_resistome/
│   ├── pca_by_bodysite.pdf            # PCA stratified by body site
│   ├── pca_by_location.pdf            # UCMC vs ZCH
│   ├── pca_scree_plot.pdf             # Variance explained
│   ├── diversity_metrics.pdf          # Shannon/Simpson/richness boxplots
│   └── hierarchical_clustering.pdf
└── antibiotic_correlations/
    ├── antibiotic_amr_correlation.pdf # Scatter plots
    └── antibiotic_exposure_amr_boxplot.pdf
```

---

## Key Findings

### 1. Dataset Overview
- **Total samples analyzed:** 641 (642 processed, 1 excluded due to naming issue: ZJH_N57_1_3_diamond)
- **Samples by location:**
  - UCMC: 360 samples
  - ZCH: 281 samples
- **Samples by body site:**
  - Axilla: 215 samples
  - Groin: 207 samples
  - Stool: 219 samples
- **Samples by time point:**
  - Week 1: 337 samples
  - Week 3: 304 samples
- **Resistance genes detected:** 6,408 unique genes (4,438 after filtering)
- **Genes with drug class annotations:** 87,849 entries (5.2% of total)

### 2. Diversity Analysis

**By Location:**
- **UCMC:** Higher diversity (mean richness: 875 genes, Shannon: 5.12)
- **ZCH:** Lower diversity (mean richness: 783 genes, Shannon: 4.76)

**By Body Site:**
- **Axilla:** Highest diversity (mean richness: 1,057 genes, Shannon: 5.39)
- **Groin:** Moderate diversity
- **Stool:** High diversity

**PCA Results:**
- PC1: 12.8% variance
- PC2: 12.0% variance
- PC3: 10.7% variance
- Body site shows stronger clustering than location

### 3. Differential Abundance: UCMC vs ZCH

**Significant genes:** 2,666 (60% of tested genes)

**Top genes enriched in UCMC:**
- **Tetracycline resistance genes dominate:**
  - `tet(W)`, `tet(O)`, `tet(W/32/O)` variants
  - Nearly absent in ZCH (median = 0)
  - Highly abundant in UCMC (log2FC = -4.4 to -7.2)
- **Arsenic resistance:**
  - `arsB_Lm` (log2FC = -7.18, p < 1e-65)
- **Pattern:** UCMC infants carry substantially more tetracycline and arsenic resistance genes

### 4. Differential Abundance: Body Sites

**Overall:** 3,399 genes significantly different across sites (77% of tested genes)

**Pairwise comparisons:**
- **Axilla vs Stool:** 3,418 genes (most differences)
- **Groin vs Stool:** 2,667 genes
- **Axilla vs Groin:** 1,755 genes (least differences)

**Top genes enriched in Axilla (skin):**
- **Macrolide resistance:** `erm(X)`, `ermC`, `erm(C)` (log2FC: -2.5 to -5.1)
- **Aminoglycoside resistance:** `aph(3')-Ia`, `aac(3)-XI`, `aph(6)-Id`
- **Beta-lactam resistance:** `pbp2m` (log2FC = -7.84)
- **Pattern:** Skin microbiome harbors distinct resistance profile, particularly macrolides and beta-lactams

**Top genes enriched in Stool:**
- Typical gut-associated resistance genes

### 5. Differential Abundance: Temporal Changes (Week 1 vs Week 3)

**Significant genes:** 3,262 (74% of tested genes)

**Top genes INCREASING from Week 1 to Week 3:**
- **Quinolone resistance genes (`oqxB` variants):**
  - Dramatic increase (log2FC ≈ 1.6, ~3-fold increase)
  - All top 10 genes are `oqxB` variants
  - Week 1 median: 0.13 RPM → Week 3 median: 14.4 RPM
- **Toxin-antitoxin system:** `relE` (log2FC = 1.13)

**Pattern:** Quinolone resistance genes accumulate substantially during first 3 weeks of life

### 6. Antibiotic-Resistance Correlations

**⭐ SURPRISING KEY FINDING: NEGATIVE CORRELATION ⭐**

**Overall correlation:**
- **Total antibiotic exposure vs Total AMR burden:** rho = **-0.150**, p = **0.0001** (highly significant)
- **Interpretation:** More antibiotic exposure is associated with LOWER resistance gene abundance

**By location:**
- UCMC: rho = 0.011, p = 0.84 (not significant)
- ZCH: rho = 0.056, p = 0.35 (not significant)

**By body site (all significantly negative):**
- Axilla: rho = -0.149, p = 0.029 ⭐
- Groin: rho = -0.170, p = 0.014 ⭐
- Stool: rho = -0.136, p = 0.045 ⭐

**By week (CRITICAL FINDING):**
- **Week 1:** rho = **-0.322**, p < **0.0001** ⭐⭐⭐ (very strong negative correlation)
- Week 3: rho = -0.050, p = 0.39 (no correlation)

**Individual antibiotics:**
- **Gentamicin:** rho = 0.088, p = 0.026 ⭐ (only significant positive correlation)
- **Vancomycin:** rho = 0.061, p = 0.11 (trend toward positive)
- Most others: weak or negative correlations

**Interpretation:**
1. **Early-life effect (Week 1):** Strong negative correlation suggests antibiotics reduce overall bacterial load (including resistant bacteria) in the first week
2. **Effect disappears by Week 3:** No correlation, possibly as resistance develops or microbiome stabilizes
3. **Site-specific:** Effect seen across all body sites, not location-specific
4. **Gentamicin exception:** Only antibiotic with positive correlation (more exposure → more resistance), as expected

### 7. Metadata Quality Issues (Resolved)

**Issue identified:** 36 ZCH samples had empty `sample_name` fields in original metadata
- These samples were excluded from original BioGPU analysis (likely QC failures)
- BUT they had complete antibiotic exposure data

**Resolution:**
- Reconstructed `sample_name` from PatientID + SampleType + SampleCollectionWeek
- All 36 samples now included in analysis
- 30/36 had antibiotic exposure, 6/36 had no antibiotics (legitimate zeros)

**ZCH sample naming:**
- Original FASTQs: `ZJH_N1_*` to `ZJH_N9_*` (no zero-padding)
- Metadata: `N01_*` to `N09_*` (zero-padded, no ZJH_ prefix)
- DIAMOND output directories: `ZJH_N01_*` to `ZJH_N09_*` (zero-padded with prefix)
- **Resolution:** Added `ZJH_` prefix to ZCH metadata sample names during merging

---

## Technical Notes

### Gene Annotations
- **5.2% of gene-sample entries** have drug_class annotations (87,849/1,680,747)
- Annotations derived from BioGPU database
- Most genes lack functional annotations (future work: expand annotation coverage)

### Statistical Approach
- **Filtering:** Genes present (RPM > 1.0) in ≥ 5 samples
- **Tests:** Non-parametric (Mann-Whitney U, Kruskal-Wallis H)
- **Multiple testing correction:** Benjamini-Hochberg FDR
- **Significance threshold:** FDR < 0.05

### Sample Exclusions
- 1 sample excluded: `ZJH_N57_1_3_diamond` (malformed name, no metadata match)
- 2 ZCH samples with metadata but no DIAMOND results: `ZJH_N01_1_3`, `ZJH_N01_2_2`

---

## Future Directions

1. **Expand gene annotations:**
   - Currently only 5.2% of entries have drug class information
   - Integrate additional resistance gene databases (CARD, ResFinder, etc.)

2. **Functional class analysis:**
   - Group genes by resistance mechanism (efflux, enzyme, target modification)
   - Analyze beta-lactam, aminoglycoside, macrolide classes separately

3. **Longitudinal analysis:**
   - Track individual subjects across time points
   - Mixed-effects models accounting for subject-level variation

4. **Microbiome integration:**
   - Correlate resistance genes with taxonomic composition
   - Identify species-specific resistance patterns

5. **Clinical outcomes:**
   - Associate resistance burden with infections, NEC, other outcomes

---

## Contact & Acknowledgments

**Analysis performed:** November 2025
**Tools:** DIAMOND v2.x, Python 3.12, pandas, scipy, scikit-learn, seaborn
**Database:** BioGPU AMR+Stress combined database
**Metadata source:** ZCH_UCMC_Manuscript NICU resistome study

---

*Last updated: November 13, 2025*
