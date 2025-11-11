# Metrics Explained - BioGPU vs Traditional Pipeline

## Output Metrics (Matching BioGPU)

Both pipelines will output the following metrics for each gene:

### 1. Read Count (Raw)
**What**: Number of reads assigned to this gene
**BioGPU**: From translated alignment
**Traditional**: From htseq-count

**Notes**:
- This is the RAW count before any normalization
- For paired-end data, need to clarify if counting reads or pairs
- htseq-count default: counts each read in a pair separately (so 1 pair = 2 counts)

### 2. RPM (Reads Per Million)
**What**: Read count normalized by library size
**Formula**:
```
RPM = (gene_read_count / total_reads) × 1,000,000
```

**CRITICAL QUESTION**: What is "total_reads"?
- **Option A**: Total individual reads (R1 + R2)
- **Option B**: Total read pairs

**Current Status**: ⚠️ MUST VERIFY how biogpu calculates this!

**Why it matters**: Using different denominators will cause 2x difference in RPM values!

### 3. TPM (Transcripts Per Million)
**What**: Read count normalized by BOTH gene length and library size
**Formula**:
```
Step 1: RPK = read_count / (gene_length_kb)
Step 2: sum_RPK = sum of all RPK values
Step 3: TPM = (RPK / sum_RPK) × 1,000,000
```

**Purpose**:
- Removes gene length bias (longer genes get more reads)
- Allows comparison between genes of different lengths
- Sum of all TPM values in a sample = 1,000,000

**Traditional pipeline**: Calculates using standard TPM formula

### 4. Mean Depth
**What**: Average number of reads covering each base of the gene
**Formula**:
```
mean_depth = total_bases_covered / gene_length
```

**BioGPU**: Outputs as `mean_depth`
**Traditional**: Calculated from bedtools coverage

**Interpretation**:
- Depth = 1: Each base covered by 1 read on average
- Depth = 10: Each base covered by 10 reads on average
- Higher depth = more confident detection

### 5. Percent Coverage
**What**: Percentage of gene bases covered by at least one read
**Formula**:
```
percent_coverage = (bases_with_coverage / gene_length) × 100
```

**BioGPU**: Outputs as `percent_coverage`
**Traditional**: From bedtools coverage

**Interpretation**:
- 100% = entire gene covered
- 50% = half of gene bases have reads
- Fragmented coverage may indicate poor alignment or chimeric sequences

**BioGPU threshold**: 50% minimum coverage filter

### 6. Gene Length
**What**: Length of the gene sequence in base pairs
**Source**: From FASTA file
**Purpose**: Reference for coverage and TPM calculations

## Comparison Strategy

### For Benchmarking, Compare:

1. **Gene Detection**
   - Which genes are detected by each method?
   - Venn diagram overlap
   - Sensitivity comparison

2. **Abundance Correlation**
   - **Primary metric**: RPM correlation (Spearman, Pearson)
   - **Secondary**: TPM correlation
   - Scatter plots biogpu vs traditional

3. **Coverage Comparison**
   - Do both methods agree on gene coverage?
   - Are there genes with high counts but low coverage? (suspicious)

4. **Discrepancy Investigation**
   - Genes detected only by biogpu: likely divergent sequences
   - Genes detected only by traditional: investigate why biogpu missed them
   - Large abundance differences (>2-fold): investigate causes

## Example Output Format

### BioGPU (_amr_abundance.tsv)
```
gene_name  gene_family  drug_class  read_count  percent_coverage  mean_depth  tpkm  rpm  evidence  confidence
arsB_pI258  arsB_pI258  ARSENIC     1409       93.26            114.23      60986  13.07  Coverage  HIGH
```

### Traditional (_abundance.tsv)
```
gene_id        gene_name  read_count  rpm     tpm      mean_depth  percent_coverage  gene_length
0|AAA21094.1  arsB_pI258  1409        13.07   12543.2  114.23      93.26            1290
```

## Key Differences to Note

### Search Space
- **BioGPU**: Protein space (6-frame translation)
  - More sensitive to divergent sequences
  - Lower nucleotide identity threshold (85%)
- **Traditional**: Nucleotide space
  - Requires high nucleotide identity (~95%)
  - May miss divergent homologs

### Coverage Filtering
- **BioGPU**: 50% minimum coverage threshold applied
- **Traditional**: No default threshold (report all)
  - Can filter post-hoc for fair comparison

### Multi-mapping
- **BioGPU**: Best hit strategy (presumably)
- **Traditional**: htseq-count default behavior
  - Need to verify handling is similar

## Expected Results

### High Agreement Zone
Genes that should agree between methods:
- High abundance genes (>100 reads)
- Conserved sequences (>95% identity)
- Full coverage (>80%)
- Well-known resistance genes

**Expected correlation**: r > 0.95

### Method-Specific Detections

**BioGPU-only genes**:
- Divergent homologs (70-90% nucleotide identity)
- Frameshift-containing genes
- Genes from distantly related species

**Traditional-only genes**:
- Possibly spurious low-count genes
- Genes with chimeric alignments
- Should be rare if methods agree

## Action Items Before Running

- [ ] **CRITICAL**: Determine biogpu RPM denominator (reads vs pairs)
- [ ] Verify biogpu multi-mapping strategy
- [ ] Verify biogpu paired-end counting method
- [ ] Document all assumptions
- [ ] Test on 1-2 samples first before batch processing

## Files to Compare

### From BioGPU
```
/fastpool/analysis/nicu_amr_stress/amr_stress/{sample}_amr_abundance.tsv
```
Columns: gene_name, gene_family, drug_class, read_count, percent_coverage, mean_depth, tpkm, rpm

### From Traditional
```
/home/david/projects/benchmark_biogpu/results/traditional/{sample}/{sample}_abundance.tsv
```
Columns: gene_id, gene_name, read_count, rpm, tpm, mean_depth, percent_coverage, gene_length

## Normalization Best Practices

### For Fair Comparison
1. ✓ Use same gene database
2. ✓ Use same RPM calculation method
3. ✓ Document any filtering differences
4. ✓ Report correlation for genes detected by both methods
5. ✓ Investigate and explain method-specific detections

### For Future Work (Bacterial Read Normalization)
Instead of total reads, use bacterial reads as denominator:
```
AMR_RPM = (AMR_reads / total_bacterial_reads) × 1,000,000
```

This requires:
1. Run pathogen profiler or similar tool
2. Quantify total bacterial reads
3. Use as denominator instead of raw reads

**Advantage**: Removes bias from host/non-bacterial contamination
