# Benchmark Project Status

**Created**: 2025-11-11
**Status**: Initial setup complete, awaiting user's traditional pipeline script

## Completed

- ✓ Directory structure created
- ✓ README with project overview
- ✓ Environment YAML with required tools (bowtie2, htseq-count, samtools, etc.)
- ✓ Setup guide documentation
- ✓ Identified biogpu pipeline parameters:
  - Min identity: 85%
  - Min coverage: 50%
  - Translated search (6-frame)
  - Output: RPM normalized abundance
- ✓ Located biogpu databases:
  - AMR+Stress: `/home/david/projects/biogpu/data/amr_stress_combined_db/`
    - DNA: 14,280 sequences (9,257 AMR + 5,023 stress)
    - Protein: 17,099 sequences (9,538 AMR + 7,561 stress)
  - FQ resistance: `/home/david/projects/biogpu/data/integrated_clean_db/`
- ✓ Located biogpu results:
  - `/fastpool/analysis/nicu_amr_stress/amr_stress/`
  - `/fastpool/analysis/nicu_amr_stress/resistance/`
- ✓ Found user's existing bowtie2 pipeline scripts (use bowtie2 + htseq-count)
- ✓ **Created database setup script**: `scripts/01_setup_databases.sh`
  - Links biogpu databases
  - Copies DNA FASTA for traditional pipeline
  - Creates GFF3 annotation file from FASTA headers
  - Builds bowtie2 index
- ✓ **Created GFF3 generator**: `scripts/create_gff3_from_fasta.py`
  - Parses biogpu FASTA headers
  - Creates htseq-count compatible GFF3
- ✓ **Created traditional pipeline with timing**: `scripts/02_run_traditional_pipeline.sh`
  - Runs bowtie2 → htseq-count → bedtools pipeline
  - Outputs read counts, RPM, TPM, coverage %
  - Captures comprehensive timing metrics (wall time, CPU time, memory)
- ✓ **Created batch processing script**: `scripts/03_batch_traditional_pipeline.sh`
  - Process multiple samples from CSV file
  - Progress tracking and error handling
- ✓ **Created biogpu timing script**: `scripts/04_run_biogpu_with_timing.sh`
  - Re-runs biogpu pipeline with timing instrumentation
  - Captures wall time, CPU time, GPU time, memory
  - Optional GPU profiling with nvprof
- ✓ **Created comprehensive timing documentation**:
  - `docs/TIMING_METRICS.md` - Full timing documentation
  - `TIMING_SUMMARY.md` - Quick reference guide
  - `IMPORTANT_NOTES.md` - Updated with timing section

## Next Steps

1. **Run database setup** (ready to execute):
   ```bash
   conda activate benchmark_biogpu
   cd /home/david/projects/benchmark_biogpu
   ./scripts/01_setup_databases.sh
   ```

2. **Select test samples**:
   - Create sample list (10-20 samples)
   - Stratified by body site, location, timepoint
   - Use samples with varying AMR burden

3. **Create traditional pipeline script**:
   - Adapt user's bowtie2 + htseq-count workflow
   - Process selected test samples
   - Output gene counts

4. **Create comparison script**:
   - Compare biogpu vs traditional results
   - Correlation analysis
   - Venn diagrams
   - Discrepancy investigation

5. **Run benchmark and generate report**

## Critical Considerations

### Fair Comparison Challenges

The biogpu and traditional methods are fundamentally different:

| Aspect | BioGPU | Traditional |
|--------|--------|-------------|
| Search space | Protein (6-frame translation) | Nucleotide |
| Sensitivity | Higher for divergent sequences | Limited to close matches |
| Speed | GPU-accelerated | CPU-based |
| Min identity | 85% (protein space) | ~95% typical (nucleotide) |
| Coverage | 50% minimum | Varies |

**Implications**:
- BioGPU may detect genes that traditional misses (especially divergent ones)
- This is a **feature**, not a bug, for resistance gene detection
- Fair comparison requires acknowledging different use cases

### What We're Actually Testing

1. **Agreement for high-identity matches**: Do both methods agree on abundant, conserved genes?
2. **Sensitivity gain**: How many additional genes does translated search detect?
3. **False positive rate**: Are biogpu-specific detections real or spurious?
4. **Quantification accuracy**: For genes detected by both, do abundances correlate?

### Validation Strategy

For genes detected only by biogpu:
- Check alignment quality (identity %, coverage %)
- Blast against NCBI to verify sequence identity
- Check if they're divergent homologs (expected)

For genes detected only by traditional:
- Check if they fail biogpu thresholds (50% coverage, 85% identity)
- Verify they're real hits (not mapping artifacts)

## Deliverables for Reviewer Response

1. **Correlation plot**: BioGPU vs Traditional RPM (for shared genes)
2. **Venn diagram**: Gene detection overlap
3. **Summary table**:
   - Total genes detected (each method)
   - Shared genes
   - Method-specific genes
   - Correlation statistics
4. **Brief text**: "We validated our novel translated-search pipeline against traditional nucleotide alignment methods and observed [X correlation] for conserved genes, with enhanced sensitivity for divergent resistance genes as expected from protein-space alignment."

## Files Ready

```
/home/david/projects/benchmark_biogpu/
├── README.md                    ✓ Created
├── environment.yml              ✓ Created
├── docs/
│   ├── SETUP_GUIDE.md          ✓ Created
│   └── PROJECT_STATUS.md       ✓ Created
├── scripts/                     [Empty - awaiting user script]
├── data/                        [Empty - will populate with test samples]
├── results/
│   ├── biogpu/                  [Ready to link]
│   └── traditional/             [Ready for pipeline output]
├── logs/                        [Empty - for pipeline logs]
└── databases/
    ├── biogpu/                  [Ready to link]
    └── traditional/             [Ready for indexing]
```

## References

- BioGPU pipeline: `/home/david/projects/biogpu/batch_process_nicu_samples_stress.sh`
- NICU results: `/home/david/projects/ZCH_UCMC_Manuscript/nicu_resistome_analysis/`
- Sample key: `/home/david/projects/biogpu/nicu_sample_key_fixed.csv`
