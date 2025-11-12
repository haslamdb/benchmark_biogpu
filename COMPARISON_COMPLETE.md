# DIAMOND vs BioGPU Comparison - Complete Summary

## What We Did

You ran a comprehensive benchmarking comparison of **DIAMOND blastx** (CPU-based translated search) vs **BioGPU** (GPU-accelerated translated search) on **50 NICU samples**.

**Configuration:**
- Both pipelines use protein-space alignment (6-frame translation)
- Identical parameters: 85% identity, 50% coverage
- Same reference database source: `/home/david/projects/biogpu/data/amr_stress_combined_db/`
- Outputs: Gene-level read counts, RPM, TPM, coverage

## Results Overview

### 1. Gene Detection (MAJOR DISCREPANCY)

|Method|Genes Detected|Per Sample Average|
|------|--------------|------------------|
|DIAMOND|129,540|2,591 genes/sample|
|BioGPU|5,826|117 genes/sample|
|Overlap|4,191 (3.1%)|84 genes/sample|

**DIAMOND detects 22x more genes than BioGPU.**

### 2. Correlation for Shared Genes (CONCERNING)

For the 4,191 genes detected by both methods:
- **Spearman ρ = 0.38**
- **Pearson r = 0.51**

This is **poor correlation** - we would expect ρ > 0.95 for the same algorithm.

### 3. Performance

- **DIAMOND**: 327 ± 159 seconds
- **BioGPU**: 520 ± 240 seconds
- **Winner**: DIAMOND is ~1.6x faster (unexpected!)

### 4. Highly Discordant Genes

Examples of genes with huge differences:
- **gyrA**: DIAMOND 1,046 RPM → BioGPU 0.004 RPM (75,000x difference)
- **erm(C)**: DIAMOND 1,136 RPM → BioGPU 0.023 RPM (35,000x difference)
- **bimA_Bm**: DIAMOND 0 RPM → BioGPU 22,543 RPM (detected only by BioGPU!)

## Key Findings

### Finding 1: Different Gene Sets
**DIAMOND-only top genes** (high abundance):
- Housekeeping genes: rpoB, rpoC, fusA
- Resistance genes: mph(A), erm(C), ermC
- Pathogenicity: clbB_Ec

**BioGPU-only top gene**:
- bimA_Bm (Burkholderia motility protein) - very high abundance in 10 samples

### Finding 2: Coverage Distribution Differs
- **DIAMOND-only genes**: 43% have < 25% coverage
- **BioGPU-only genes**: 74% have < 25% coverage
- Suggests BioGPU filters more aggressively despite same 50% threshold parameter

### Finding 3: Poor Correlation Across All Abundance Levels
| Abundance | Pearson r | Spearman ρ | n genes |
|-----------|-----------|------------|---------|
| Very Low (0-1 RPM) | - | 0.21 | 2,088 |
| Low (1-10 RPM) | - | 0.13 | 978 |
| Medium (10-100 RPM) | - | 0.28 | 766 |
| High (>100 RPM) | - | 0.23 | 359 |

Even high-abundance genes show poor agreement.

## Why Are Results So Different?

### Hypothesis 1: Database Differences (MOST LIKELY)
**Evidence:**
- bimA_Bm detected at very high levels by BioGPU but completely absent from DIAMOND results
- This strongly suggests the databases are NOT identical

**Next Step:** Compare the actual FASTA sequences used:
```bash
# Count sequences
grep -c "^>" databases/biogpu/amr_stress_combined_db/protein.fasta
# BioGPU: 17,099 proteins

# Check if bimA_Bm is in both databases
grep "bimA" databases/biogpu/amr_stress_combined_db/protein.fasta
# Found in BioGPU

# Check DIAMOND database source
# The .dmnd file exists but we need to verify it was built from the same FASTA
```

### Hypothesis 2: Additional Filtering in BioGPU
**Evidence:**
- BioGPU detects 22x fewer genes
- BioGPU output includes "confidence" levels (LOW/MODERATE/HIGH)
- BioGPU amr_abundance.tsv only includes genes passing some threshold

**Possible filters:**
- Minimum read count threshold
- Minimum depth/quality scores
- Coverage calculation differs
- Multi-mapping read handling

### Hypothesis 3: Algorithm Implementation Differences
**Evidence:**
- Poor correlation (ρ=0.38) even for shared genes
- DIAMOND uses `--top 1` (best hit per read)
- BioGPU alignment algorithm details unknown

### Hypothesis 4: Different Gene ID Formats
**DIAMOND format:**
```
0|AAA16360.1|1|1|stxA2b|stxA2b||1|stxA2b|STX2|Shiga_toxin...
```

**BioGPU format:**
```
Gene name: arsB_R773, bimA_Bm, etc.
```

**But:** Our comparison script uses gene_name in both, so this shouldn't cause issues unless there are name mismatches.

## Generated Files

All results in: `/home/david/projects/benchmark_biogpu/results/comparison/`

1. **comparison_summary.txt** - Statistical summary
2. **merged_comparison_data.tsv** - Full dataset (133,271 gene-sample combinations)
3. **scatter_rpm_comparison.png** - RPM correlation plots
4. **gene_detection_overlap.png** - Detection overlap by sample
5. **timing_comparison.png** - Performance comparison
6. **FINDINGS_SUMMARY.md** - Detailed analysis with recommendations

## Recommendations

### Immediate Actions (Priority Order)

**1. Verify Database Identity** ⚠️ CRITICAL
- The .dmnd file exists but we need to confirm it was built from the same protein.fasta
- Check if DIAMOND database includes all 17,099 proteins from BioGPU
- Specifically check if bimA_Bm is in DIAMOND database

**2. Compare Raw BioGPU Output**
- BioGPU likely produces more comprehensive output files
- The amr_abundance.tsv we used might be filtered
- Check if there's a raw hits file with all alignments

**3. Review BioGPU Filtering**
- Examine BioGPU source code for filtering logic
- Identify any undocumented thresholds
- Understand the "confidence" scoring system

**4. Investigate Specific Genes**
- Why is bimA_Bm absent from DIAMOND results?
- Why is gyrA/erm(C) so different between methods?
- Manually inspect alignments for a test sample

### For Your Manuscript

**Current Status:** ⚠️ **Cannot validate BioGPU results with current comparison**

The poor correlation (ρ=0.38) and 22x difference in detection indicates these are fundamentally different analytical approaches, not just CPU vs GPU implementations.

**Options:**

**Option A:** Resolve discrepancies first (RECOMMENDED)
- Timeline: 1-2 weeks
- Identify cause of low correlation
- Re-run comparison after fixes
- Validate BioGPU results properly

**Option B:** Use DIAMOND results instead
- Timeline: 2-3 weeks for full re-analysis
- More conservative (detects more genes)
- Well-established method
- Cite BioGPU as preliminary/exploratory

**Option C:** Report both methods
- Acknowledge discrepancies in methods section
- Report genes detected by both methods separately
- Note that validation is ongoing
- Weakens manuscript but honest

## Scripts Created

1. **`scripts/compare_pipelines.py`** - Main comparison script
   - Loads DIAMOND and BioGPU results
   - Merges by sample and gene_name
   - Calculates correlations
   - Generates visualizations

2. **`scripts/investigate_differences.py`** - Detailed analysis
   - Gene detection patterns
   - Abundance distributions
   - Coverage comparisons
   - Identifies discordant genes

## How to Use These Results

### View the comparison:
```bash
cd /home/david/projects/benchmark_biogpu

# View summary
cat results/comparison/comparison_summary.txt

# View detailed findings
less results/comparison/FINDINGS_SUMMARY.md

# View plots
xdg-open results/comparison/scatter_rpm_comparison.png
xdg-open results/comparison/timing_comparison.png
```

### Re-run comparison (if you fix issues):
```bash
python scripts/compare_pipelines.py
```

### Investigate specific genes:
```bash
# Find a gene in the merged data
grep "bimA_Bm" results/comparison/merged_comparison_data.tsv | head

# Check if gene exists in database
grep -i "bimA" databases/biogpu/amr_stress_combined_db/protein.fasta
```

## Next Steps

Before proceeding with the manuscript reviewer response, you need to:

1. ✅ Understand why results differ (this analysis complete)
2. ⬜ Verify database identity (check DIAMOND .dmnd contents)
3. ⬜ Investigate bimA_Bm discrepancy (why absent from DIAMOND?)
4. ⬜ Review BioGPU filtering logic (source code review)
5. ⬜ Decide which method to use for manuscript
6. ⬜ Re-run analysis if needed

## Questions to Answer

1. **Is the DIAMOND .dmnd file built from the exact same protein.fasta as BioGPU?**
   - Location: `databases/traditional/amr_stress_protein.dmnd`
   - Source should be: `databases/biogpu/amr_stress_combined_db/protein.fasta`

2. **Does BioGPU output a raw, unfiltered hits file?**
   - Check: `results/biogpu/*/ZJH_N57_1_3/ZJH_N57_1_3_hits.tsv`
   - Compare to: `amr_abundance.tsv`

3. **What is BioGPU's confidence scoring system?**
   - LOW/MODERATE/HIGH categories
   - Are LOW confidence genes excluded from downstream analysis?

4. **How does BioGPU handle multi-mapping reads?**
   - Best hit? Random? Fractional?
   - Compare to DIAMOND's `--top 1`

---

**Summary:** You successfully ran a comprehensive benchmarking comparison, but the results reveal significant methodological differences that need investigation before the methods can be considered validated. The analysis scripts and visualizations are complete and ready to use once you resolve the underlying issues.
