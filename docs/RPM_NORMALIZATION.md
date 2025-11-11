# RPM Normalization - Critical for Fair Comparison

## Current Situation

### BioGPU Pipeline
- Outputs both `read_count` (raw counts) and `rpm` (reads per million)
- Used RPM for all NICU resistome analysis
- **Question**: Is RPM calculated per million **reads** or per million **read pairs**?

### Formula
```
RPM = (gene_read_count / total_reads) × 1,000,000
```

## Need to Determine

### Example Calculation (Sample N01_1_2)

From `/fastpool/analysis/nicu_amr_stress/amr_stress/N01_1_2_amr_abundance.tsv`:

| Gene | read_count | rpm |
|------|-----------|-----|
| arsB_pI258 | 1409 | 13.0737 |
| arsB_R773 | 413 | 3.83209 |
| bimA_Bm | 3116 | 32.8029 |

**Reverse calculation**:
```
total = (read_count / rpm) × 1,000,000

arsB_pI258: 1409 / 13.0737 × 1,000,000 = 107,782,688
arsB_R773:   413 / 3.83209 × 1,000,000 = 107,767,598
bimA_Bm:    3116 / 32.8029 × 1,000,000 =  95,011,583
```

**Average**: ~107.77 million for first two (third is different, maybe different sample?)

## Critical Question for Fair Comparison

**Is this 107.77 million:**
1. **Total individual reads** (R1 + R2 combined)?
   - Meaning ~53.9M read pairs
   - RPM = (gene_reads / (R1_reads + R2_reads)) × 1M

2. **Total read pairs**?
   - Meaning ~215.5M individual reads
   - RPM = (gene_reads / num_read_pairs) × 1M

## Why This Matters

For benchmarking, we need to calculate RPM the **exact same way** in the traditional pipeline:

### If using total individual reads:
```bash
# Count total reads in FASTQ files
total_r1=$(zcat sample_R1.fastq.gz | wc -l | awk '{print $1/4}')
total_r2=$(zcat sample_R2.fastq.gz | wc -l | awk '{print $1/4}')
total_reads=$((total_r1 + total_r2))

# Calculate RPM for each gene
rpm = (htseq_count / total_reads) × 1,000,000
```

### If using read pairs:
```bash
# Count read pairs (just count R1)
total_pairs=$(zcat sample_R1.fastq.gz | wc -l | awk '{print $1/4}')

# Calculate RPM for each gene
rpm = (htseq_count / total_pairs) × 1,000,000
```

## How to Determine

### Option 1: Check FASTQ file
```bash
# For sample N01_1_2, find FASTQ files
grep "N01_1_2" /home/david/projects/biogpu/nicu_sample_key_fixed.csv

# Count reads in R1
zcat [path_to_R1].fastq.gz | wc -l | awk '{print $1/4}'

# Count reads in R2
zcat [path_to_R2].fastq.gz | wc -l | awk '{print $1/4}'

# Compare to 107.77M
```

### Option 2: Check biogpu source code
Look for where RPM is calculated in the C++ source:
```bash
grep -r "rpm\|RPM\|normali" /home/david/projects/biogpu/src/
```

### Option 3: Check biogpu documentation
Look for normalization documentation:
```bash
grep -r "reads per million\|RPM\|normalization" /home/david/projects/biogpu/docs/
```

## Recommendation

**IMPORTANT**: We need to clarify this with you before running the benchmark!

Once we know the denominator, we'll document it clearly and ensure both pipelines use identical normalization.

## Future Normalization Strategy

You mentioned wanting to normalize by bacterial reads instead:

```
AMR RPM = (AMR_gene_reads / total_bacterial_reads) × 1,000,000
```

This makes more biological sense because:
- Removes bias from host DNA contamination
- Removes bias from non-bacterial reads (fungi, viruses, etc.)
- More accurately reflects AMR burden in the bacterial community

**Implementation**:
1. Use pathogen profiler to quantify bacterial reads
2. Sum all reads assigned to bacterial taxa
3. Use this as denominator instead of total raw reads

This would require re-running the analysis with bacterial read quantification.

## Action Items

- [ ] Determine if biogpu uses total reads or read pairs for RPM
- [ ] Document the exact denominator used
- [ ] Implement identical normalization in traditional pipeline
- [ ] Add this information to benchmark comparison script

## Notes

For the current benchmark, we'll match whatever biogpu is using.
For future work, bacterial read normalization would be ideal but requires additional analysis.
