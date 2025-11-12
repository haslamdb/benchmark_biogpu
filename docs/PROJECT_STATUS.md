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
- ✓ **Switched to DIAMOND as primary method** (translated search matching biogpu)
- ✓ **Created database setup script**: `scripts/01_setup_databases.sh`
  - Links biogpu databases
  - Copies protein FASTA for DIAMOND pipeline
  - Builds DIAMOND index (.dmnd file)
  - Copies DNA FASTA for bowtie2 (optional)
  - Creates GFF3 annotation file from FASTA headers
  - Builds bowtie2 index (optional)
- ✓ **Created GFF3 generator**: `scripts/create_gff3_from_fasta.py`
  - Parses biogpu FASTA headers
  - Creates htseq-count compatible GFF3
- ✓ **Created DIAMOND pipeline with timing**: `scripts/02_run_diamond_pipeline.sh`
  - Runs DIAMOND blastx (translated search)
  - Parameters: 85% identity, 50% coverage (matching biogpu)
  - Outputs read counts, RPM, TPM, coverage %
  - Captures comprehensive timing metrics (wall time, CPU time, memory)
- ✓ **Created traditional bowtie2 pipeline**: `scripts/02_run_traditional_pipeline.sh`
  - Runs bowtie2 → htseq-count → bedtools pipeline
  - Retained as optional alternative for nucleotide-space comparison
- ✓ **Created batch processing script**: `scripts/03_batch_traditional_pipeline.sh`
  - Process multiple samples from CSV file (runs DIAMOND pipeline)
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

1. **Run DIAMOND pipeline** (~3-4 hours):
   ```bash
   conda activate benchmark_biogpu
   cd /home/david/projects/benchmark_biogpu
   ./scripts/03_batch_traditional_pipeline.sh
   ```

2. **Run biogpu pipeline with timing** (~4-8 hours):
   ```bash
   ./scripts/05_batch_biogpu_with_timing.sh --gpu-id 0
   ```

3. **Aggregate results**:
   - Combine timing data from all samples
   - Collect abundance data for comparison

4. **Create comparison script**:
   - Compare biogpu vs DIAMOND results
   - Correlation analysis (Spearman, Pearson)
   - Venn diagrams (gene detection overlap)
   - Performance metrics (speedup, memory usage)

5. **Generate figures and report** for reviewer response

## Critical Considerations

### Fair Comparison Strategy

**Key Decision**: Switch from bowtie2 to DIAMOND for fair comparison

| Aspect | BioGPU | DIAMOND (Primary) | Bowtie2 (Optional) |
|--------|--------|-------------------|-------------------|
| Search space | Protein (6-frame) | Protein (6-frame) | Nucleotide |
| Sensitivity | High (divergent genes) | High (divergent genes) | Lower (close matches) |
| Speed | GPU-accelerated | CPU-based | CPU-based |
| Min identity | 85% (protein) | 85% (protein) | ~95% typical (DNA) |
| Coverage | 50% minimum | 50% minimum | Varies |

**Why DIAMOND?**
- Both methods now search in protein space with identical parameters
- Fair, apples-to-apples comparison
- Can both detect divergent resistance genes
- Differences will reflect implementation, not fundamental methodology

**Bowtie2 retained as optional alternative**:
- Shows sensitivity difference between nucleotide and protein-space search
- Demonstrates why translated search is important for AMR detection

### What We're Actually Testing

1. **Method agreement**: Do biogpu and DIAMOND (both translated search) detect similar genes?
2. **Quantification correlation**: For shared genes, do RPM values correlate (target: r > 0.90)?
3. **Performance comparison**: GPU vs CPU for translated search - speedup and memory usage
4. **Implementation validation**: Are there systematic differences despite same methodology?

### Validation Strategy

Since both methods now use translated search with same parameters, discrepancies indicate:
- Implementation differences (algorithmic choices)
- Edge cases in alignment or counting
- Potential bugs or issues to investigate

For genes with different abundances:
- Check alignment quality metrics
- Verify read counting logic
- Check for edge cases (gene boundaries, multi-mapping)

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
├── README.md                                  ✓ Updated (DIAMOND)
├── QUICKSTART.md                              ✓ Updated (DIAMOND)
├── RUN_BENCHMARK.md                           ✓ Updated (DIAMOND)
├── IMPORTANT_NOTES.md                         ✓ Updated (DIAMOND)
├── environment.yml                            ✓ Created
├── docs/
│   ├── SETUP_GUIDE.md                        ✓ Updated (DIAMOND)
│   ├── PROJECT_STATUS.md                     ✓ Updated (DIAMOND)
│   ├── TIMING_METRICS.md                     ✓ Created
│   └── RPM_NORMALIZATION.md                  ✓ Created
├── scripts/
│   ├── 01_setup_databases.sh                 ✓ Created (DIAMOND + bowtie2)
│   ├── 02_run_diamond_pipeline.sh            ✓ Created (PRIMARY METHOD)
│   ├── 02_run_traditional_pipeline.sh        ✓ Created (bowtie2 - optional)
│   ├── 03_batch_traditional_pipeline.sh      ✓ Created (runs DIAMOND)
│   ├── 04_run_biogpu_with_timing.sh          ✓ Created
│   ├── 05_batch_biogpu_with_timing.sh        ✓ Created
│   └── create_gff3_from_fasta.py             ✓ Created
├── data/
│   └── test_samples.csv                      ✓ Created (50 samples)
├── results/
│   ├── biogpu/                               [Ready for results]
│   └── traditional/                          [Ready for DIAMOND results]
├── logs/                                      [Ready for logs]
└── databases/
    ├── biogpu/                               [Ready to link]
    └── traditional/                          [Ready for DIAMOND + bowtie2]
```

## References

- BioGPU pipeline: `/home/david/projects/biogpu/batch_process_nicu_samples_stress.sh`
- NICU results: `/home/david/projects/ZCH_UCMC_Manuscript/nicu_resistome_analysis/`
- Sample key: `/home/david/projects/biogpu/nicu_sample_key_fixed.csv`
