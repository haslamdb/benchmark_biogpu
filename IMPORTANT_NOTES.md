# Important Notes for Benchmarking

## 🔴 CRITICAL: RPM Normalization Must Match!

Before running the benchmark, we **MUST** determine exactly how biogpu calculates RPM:

### The Question:
```
RPM = (gene_read_count / DENOMINATOR) × 1,000,000
```

**What is the DENOMINATOR?**
- **Option A**: Total individual reads (R1 + R2)?
- **Option B**: Total read pairs (count of pairs)?

### Why This Matters:
If biogpu uses Option A and we use Option B (or vice versa), our RPMs will be off by 2x!

### How to Find Out:

**Method 1**: Check a FASTQ file
```bash
# Find sample N01_1_2 in sample key
grep "N01_1_2" /home/david/projects/biogpu/nicu_sample_key_fixed.csv

# Count reads
zcat R1_file.fastq.gz | wc -l  # Divide by 4
zcat R2_file.fastq.gz | wc -l  # Divide by 4

# Compare to ~107.77 million (calculated from biogpu output)
```

**Method 2**: Check biogpu source or documentation

**Method 3**: Ask someone who ran the original analysis

### Current Best Guess:
Based on reverse calculation from N01_1_2 output:
- Denominator appears to be **~107.77 million**
- Need to verify if this is:
  - ~107.77M individual reads (so ~53.9M pairs), OR
  - ~107.77M read pairs (so ~215.5M individual reads)

### Action Required:
☐ Determine biogpu RPM denominator
☐ Document it clearly
☐ Implement identical calculation in traditional pipeline
☐ Add verification step in comparison script

See `docs/RPM_NORMALIZATION.md` for detailed analysis.

---

## ⏱️ Performance Benchmarking (NEW!)

Both pipelines now capture comprehensive timing metrics:

### What's Measured:
- **Wall clock time**: Total elapsed time per step
- **CPU time**: User + system CPU time
- **GPU time**: GPU execution time (biogpu only, with nvprof)
- **Memory usage**: Peak RAM usage per step
- **Thread utilization**: Number of threads/cores used

### Output Files:
- `{sample}_timing.tsv` - Machine-readable timing data
- `{sample}_stats.txt` - Human-readable summary with timing
- `{sample}_*_time.txt` - Detailed /usr/bin/time output

### Why This Matters:
- Compare computational efficiency
- Quantify GPU speedup
- Identify bottlenecks
- Report performance for manuscript

See `docs/TIMING_METRICS.md` for full documentation.

---

## Other Key Considerations

### 1. Coverage Thresholds
**BioGPU**:
- Min identity: 85%
- Min coverage: 50%

**Traditional Pipeline**:
- bowtie2: default settings (~95% identity typical)
- No explicit coverage filter (can add with bedtools)

**Solution**: Document genes detected by each method and investigate discrepancies

### 2. Read Counting Differences

**BioGPU**: Likely counts best hit per read
**htseq-count**: Has multiple modes:
- `union` (default)
- `intersection-strict`
- `intersection-nonempty`

**Solution**: Use same logic or document differences

### 3. Multi-mapping Reads

**BioGPU**: How are multi-mapping reads handled?
**bowtie2 + htseq-count**: Default behavior needs to match

**Solution**: Check both pipelines' handling of multi-mappers

### 4. GFF3 Coordinate System

Our GFF3 uses placeholder coordinates (1 to gene_length) because:
- We're aligning to individual gene sequences
- Each gene is its own "chromosome" in bowtie2
- htseq-count just needs seqname to match between SAM and GFF3

**This should work correctly.**

### 5. Gene Length Normalization

**BioGPU outputs**:
- `read_count`: Raw counts
- `rpm`: Reads per million (normalized by library size)
- `tpkm`: Likely TPM or RPKM (normalized by gene length too)

**For comparison**: Use RPM (not TPKM) since we want library-size normalized only

### 6. Paired-End Read Counting

**Question**: Does htseq-count count:
- Each read in a pair separately (2 counts for 1 pair), OR
- Each pair as one (1 count for 1 pair)

**Default htseq-count**: Counts each read separately (so pairs = 2 counts)

**Need to verify**: How does biogpu count paired reads?

---

## Benchmark Success Criteria

### Primary Goal:
Correlation > 0.90 for genes detected by both methods

### Secondary Goals:
1. Understand which genes are detected by only one method (and why)
2. Quantify sensitivity gain from translated search
3. Validate that abundant, conserved genes show high agreement

### Deliverables for Reviewer:
1. Scatter plot: biogpu RPM vs traditional RPM
2. Venn diagram: Gene detection overlap
3. Correlation statistics (Spearman, Pearson)
4. 1-2 supplementary figures
5. Brief methods text: "We validated our translated-search pipeline..."

---

## Next Steps (In Order)

1. ☐ **Clarify RPM calculation** (CRITICAL!)
2. ☐ Create conda environment
3. ☐ Run database setup script
4. ☐ Select 10-20 test samples
5. ☐ Extract biogpu results for test samples
6. ☐ Create traditional pipeline script (with correct RPM calculation!)
7. ☐ Run traditional pipeline
8. ☐ Compare results
9. ☐ Generate figures and report
10. ☐ Add to reviewer response

---

## Questions to Resolve

- [ ] Is RPM per million reads or per million read pairs?
- [ ] How does biogpu handle multi-mapping reads?
- [ ] How does biogpu count paired-end reads (separately or as pairs)?
- [ ] Should we filter by coverage in traditional pipeline?
- [ ] What htseq-count mode should we use?

**Document all decisions for reproducibility!**
