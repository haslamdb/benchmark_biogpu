#!/bin/bash
################################################################################
# Traditional Pipeline: bowtie2 + htseq-count + bedtools
# Outputs: read counts, RPM, coverage %, TPM (matching biogpu output)
################################################################################

set -euo pipefail

BENCHMARK_DIR="/home/david/projects/benchmark_biogpu"
DB_DIR="${BENCHMARK_DIR}/databases/traditional"
OUTPUT_DIR="${BENCHMARK_DIR}/results/traditional"
LOG_DIR="${BENCHMARK_DIR}/logs"
TEMP_DIR="${BENCHMARK_DIR}/temp"

mkdir -p "${OUTPUT_DIR}" "${LOG_DIR}" "${TEMP_DIR}"

# Database files
BT2_INDEX="${DB_DIR}/amr_stress_bt2"
GFF3_FILE="${DB_DIR}/amr_stress_genes.gff3"
FASTA_FILE="${DB_DIR}/amr_stress_dna.fasta"

# Processing parameters
THREADS=32

################################################################################
# Usage
################################################################################

usage() {
    cat << EOF
Usage: $0 --sample SAMPLE_NAME --r1 R1_FILE --r2 R2_FILE [options]

Required:
  --sample NAME    Sample name (e.g., N01_1_2)
  --r1 FILE        Path to R1 FASTQ file (can be .gz)
  --r2 FILE        Path to R2 FASTQ file (can be .gz)

Optional:
  --threads N      Number of threads (default: 12)
  --keep-sam       Keep SAM file (default: delete after processing)
  --keep-bam       Keep BAM files (default: delete after processing)

Output files created in ${OUTPUT_DIR}/SAMPLE_NAME/:
  - SAMPLE_NAME_counts.tsv          # htseq-count output (raw counts)
  - SAMPLE_NAME_abundance.tsv       # Read counts + RPM + TPM + Coverage
  - SAMPLE_NAME_coverage.bed        # bedtools coverage output (detailed)
  - SAMPLE_NAME_stats.txt           # Alignment statistics

Example:
  $0 --sample N01_1_2 \\
     --r1 /path/to/sample_R1.fastq.gz \\
     --r2 /path/to/sample_R2.fastq.gz \\
     --threads 16

EOF
    exit 1
}

################################################################################
# Parse arguments
################################################################################

SAMPLE_NAME=""
R1_FILE=""
R2_FILE=""
KEEP_SAM=false
KEEP_BAM=false

while [[ $# -gt 0 ]]; do
    case $1 in
        --sample)
            SAMPLE_NAME="$2"
            shift 2
            ;;
        --r1)
            R1_FILE="$2"
            shift 2
            ;;
        --r2)
            R2_FILE="$2"
            shift 2
            ;;
        --threads)
            THREADS="$2"
            shift 2
            ;;
        --keep-sam)
            KEEP_SAM=true
            shift
            ;;
        --keep-bam)
            KEEP_BAM=true
            shift
            ;;
        -h|--help)
            usage
            ;;
        *)
            echo "Unknown option: $1"
            usage
            ;;
    esac
done

# Check required arguments
if [[ -z "${SAMPLE_NAME}" ]] || [[ -z "${R1_FILE}" ]] || [[ -z "${R2_FILE}" ]]; then
    echo "ERROR: Missing required arguments"
    usage
fi

# Check files exist
if [[ ! -f "${R1_FILE}" ]]; then
    echo "ERROR: R1 file not found: ${R1_FILE}"
    exit 1
fi

if [[ ! -f "${R2_FILE}" ]]; then
    echo "ERROR: R2 file not found: ${R2_FILE}"
    exit 1
fi

################################################################################
# Setup
################################################################################

SAMPLE_DIR="${OUTPUT_DIR}/${SAMPLE_NAME}"
mkdir -p "${SAMPLE_DIR}"

LOGFILE="${LOG_DIR}/${SAMPLE_NAME}_traditional_$(date +%Y%m%d_%H%M%S).log"

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "${LOGFILE}"
}

log "========================================="
log "Traditional Pipeline: ${SAMPLE_NAME}"
log "========================================="
log "R1: ${R1_FILE}"
log "R2: ${R2_FILE}"
log "Output: ${SAMPLE_DIR}"
log "Threads: ${THREADS}"
log ""

################################################################################
# Step 1: Bowtie2 alignment
################################################################################

log "Step 1: Running bowtie2 alignment..."

SAM_FILE="${SAMPLE_DIR}/${SAMPLE_NAME}.sam"

# Track both wall time and CPU time
STEP_START=$(date +%s)
/usr/bin/time -v -o "${SAMPLE_DIR}/${SAMPLE_NAME}_bowtie2_time.txt" \
    bowtie2 \
        --no-unal \
        -p ${THREADS} \
        -q \
        -x "${BT2_INDEX}" \
        -1 "${R1_FILE}" \
        -2 "${R2_FILE}" \
        -S "${SAM_FILE}" \
        2> "${SAMPLE_DIR}/${SAMPLE_NAME}_bowtie2.log"

STEP_END=$(date +%s)
STEP_DURATION=$((STEP_END - STEP_START))

# Extract CPU time from /usr/bin/time output
CPU_TIME=$(grep "User time" "${SAMPLE_DIR}/${SAMPLE_NAME}_bowtie2_time.txt" | awk '{print $4}')
SYS_TIME=$(grep "System time" "${SAMPLE_DIR}/${SAMPLE_NAME}_bowtie2_time.txt" | awk '{print $4}')
MAX_MEM=$(grep "Maximum resident set size" "${SAMPLE_DIR}/${SAMPLE_NAME}_bowtie2_time.txt" | awk '{print $6}')

log "✓ Bowtie2 completed"
log "  Wall time: ${STEP_DURATION}s"
log "  CPU time: ${CPU_TIME}s user + ${SYS_TIME}s system"
log "  Max memory: ${MAX_MEM} KB"

# Store timing info
BT2_WALL_TIME=${STEP_DURATION}
BT2_CPU_TIME=${CPU_TIME}
BT2_MEM=${MAX_MEM}

# Extract alignment stats
TOTAL_PAIRS=$(grep "reads; of these:" "${SAMPLE_DIR}/${SAMPLE_NAME}_bowtie2.log" | awk '{print $1}')
ALIGNED=$(grep "aligned concordantly exactly 1 time" "${SAMPLE_DIR}/${SAMPLE_NAME}_bowtie2.log" | awk '{print $1}')
log "  Total pairs: ${TOTAL_PAIRS}"
log "  Aligned pairs: ${ALIGNED}"

################################################################################
# Step 2: Convert to BAM and sort
################################################################################

log ""
log "Step 2: Converting SAM to sorted BAM..."

BAM_FILE="${SAMPLE_DIR}/${SAMPLE_NAME}.bam"
SORTED_BAM="${SAMPLE_DIR}/${SAMPLE_NAME}.sorted.bam"

STEP_START=$(date +%s)
/usr/bin/time -v -o "${SAMPLE_DIR}/${SAMPLE_NAME}_samtools_time.txt" bash -c "
    samtools view -@ ${THREADS} -bS '${SAM_FILE}' > '${BAM_FILE}' && \
    samtools sort -@ ${THREADS} '${BAM_FILE}' -o '${SORTED_BAM}' && \
    samtools index '${SORTED_BAM}'
"
STEP_END=$(date +%s)
STEP_DURATION=$((STEP_END - STEP_START))

SAMTOOLS_CPU_TIME=$(grep "User time" "${SAMPLE_DIR}/${SAMPLE_NAME}_samtools_time.txt" | awk '{print $4}')
SAMTOOLS_MEM=$(grep "Maximum resident set size" "${SAMPLE_DIR}/${SAMPLE_NAME}_samtools_time.txt" | awk '{print $6}')

log "✓ BAM file created and sorted"
log "  Wall time: ${STEP_DURATION}s"
log "  CPU time: ${SAMTOOLS_CPU_TIME}s"
log "  Max memory: ${SAMTOOLS_MEM} KB"

SAMTOOLS_WALL_TIME=${STEP_DURATION}

################################################################################
# Step 3: htseq-count for read quantification
################################################################################

log ""
log "Step 3: Running htseq-count..."

COUNTS_FILE="${SAMPLE_DIR}/${SAMPLE_NAME}_counts.tsv"

STEP_START=$(date +%s)
/usr/bin/time -v -o "${SAMPLE_DIR}/${SAMPLE_NAME}_htseq_time.txt" \
    htseq-count \
        "${SAM_FILE}" \
        "${GFF3_FILE}" \
        -f sam \
        -a 1 \
        --idattr ID \
        --type gene \
        -s no \
        --additional-attr=Name \
        > "${COUNTS_FILE}" 2> "${SAMPLE_DIR}/${SAMPLE_NAME}_htseq.log"

STEP_END=$(date +%s)
STEP_DURATION=$((STEP_END - STEP_START))

HTSEQ_CPU_TIME=$(grep "User time" "${SAMPLE_DIR}/${SAMPLE_NAME}_htseq_time.txt" | awk '{print $4}')
HTSEQ_MEM=$(grep "Maximum resident set size" "${SAMPLE_DIR}/${SAMPLE_NAME}_htseq_time.txt" | awk '{print $6}')

# Check for errors in htseq-count
if grep -q "Error\|Warning" "${SAMPLE_DIR}/${SAMPLE_NAME}_htseq.log"; then
    log "⚠ htseq-count produced warnings/errors - check ${SAMPLE_NAME}_htseq.log"
fi

# Count genes with reads
NUM_GENES=$(grep -v "^__" "${COUNTS_FILE}" | awk '$2 > 0' | wc -l)
log "✓ htseq-count completed: ${NUM_GENES} genes with reads"
log "  Wall time: ${STEP_DURATION}s"
log "  CPU time: ${HTSEQ_CPU_TIME}s"
log "  Max memory: ${HTSEQ_MEM} KB"

HTSEQ_WALL_TIME=${STEP_DURATION}

################################################################################
# Step 4: bedtools coverage for coverage statistics
################################################################################

log ""
log "Step 4: Calculating gene coverage with bedtools..."

COVERAGE_FILE="${SAMPLE_DIR}/${SAMPLE_NAME}_coverage.bed"

STEP_START=$(date +%s)
/usr/bin/time -v -o "${SAMPLE_DIR}/${SAMPLE_NAME}_bedtools_time.txt" \
    bedtools coverage \
        -a "${GFF3_FILE}" \
        -b "${SORTED_BAM}" \
        > "${COVERAGE_FILE}"

STEP_END=$(date +%s)
STEP_DURATION=$((STEP_END - STEP_START))

BEDTOOLS_CPU_TIME=$(grep "User time" "${SAMPLE_DIR}/${SAMPLE_NAME}_bedtools_time.txt" | awk '{print $4}')
BEDTOOLS_MEM=$(grep "Maximum resident set size" "${SAMPLE_DIR}/${SAMPLE_NAME}_bedtools_time.txt" | awk '{print $6}')

log "✓ Coverage calculation completed"
log "  Wall time: ${STEP_DURATION}s"
log "  CPU time: ${BEDTOOLS_CPU_TIME}s"
log "  Max memory: ${BEDTOOLS_MEM} KB"

BEDTOOLS_WALL_TIME=${STEP_DURATION}

################################################################################
# Step 5: Calculate RPM, TPM, and combine all metrics
################################################################################

log ""
log "Step 5: Calculating RPM, TPM, and combining metrics..."

# Use Python to combine all metrics
STEP_START=$(date +%s)

python3 << PYTHON_SCRIPT > "${SAMPLE_DIR}/${SAMPLE_NAME}_abundance.tsv"
import sys
import re

# Parse command line args passed via environment
sample_name = "${SAMPLE_NAME}"
counts_file = "${COUNTS_FILE}"
coverage_file = "${COVERAGE_FILE}"
fasta_file = "${FASTA_FILE}"
total_pairs = int("${TOTAL_PAIRS}")

# TODO: CRITICAL - Determine if biogpu uses total_reads or total_pairs!
# For now, we'll calculate both and include notes
# Assumption: htseq-count counts individual reads, so a pair = 2 counts
# Therefore, if biogpu uses pairs, we need total_pairs as denominator
# If biogpu uses reads, we need total_pairs * 2 as denominator

# **PLACEHOLDER**: Set this based on biogpu behavior
USE_PAIRS = True  # If True, RPM = (count / total_pairs) * 1M
                   # If False, RPM = (count / (total_pairs*2)) * 1M

if USE_PAIRS:
    total_for_rpm = total_pairs
    rpm_note = "per million pairs"
else:
    total_for_rpm = total_pairs * 2
    rpm_note = "per million reads"

# Read gene lengths from FASTA
gene_lengths = {}
current_gene = None

with open(fasta_file, 'r') as f:
    for line in f:
        if line.startswith('>'):
            # Parse header: >index|protein_id|...
            header = line.strip()[1:]
            parts = header.split('|')
            if len(parts) >= 2:
                gene_id = f"{parts[0]}|{parts[1]}"
                current_gene = gene_id
                gene_lengths[gene_id] = 0
        else:
            if current_gene:
                gene_lengths[current_gene] += len(line.strip())

# Read htseq-count output
# Format: gene_id \t gene_name \t count (due to --additional-attr=Name)
counts = {}
gene_names = {}

with open(counts_file, 'r') as f:
    for line in f:
        if line.startswith('__'):
            continue
        parts = line.strip().split('\t')
        if len(parts) >= 3:
            # 3-column format: gene_id, gene_name, count
            gene_id = parts[0]
            gene_name = parts[1]
            count = int(parts[2])
            counts[gene_id] = count
            gene_names[gene_id] = gene_name
        elif len(parts) >= 2:
            # 2-column format: gene_id, count (fallback)
            gene_id = parts[0]
            count = int(parts[1])
            counts[gene_id] = count
            gene_names[gene_id] = gene_id

# Read bedtools coverage output
# Format: chr, start, end, feature, score, strand, frame, attributes, coverage_depth, bases_covered, gene_length, fraction_covered
coverage_data = {}

with open(coverage_file, 'r') as f:
    for line in f:
        if line.startswith('#'):
            continue
        parts = line.strip().split('\t')
        if len(parts) >= 10:
            gene_id = parts[0]  # First column is the seqname (gene_id)
            # bedtools coverage last 4 columns:
            # - number of features/bases covered
            # - length of feature
            # - fraction of feature covered
            # We want the last column (fraction covered)
            bases_covered = int(parts[-3])
            gene_length = int(parts[-2])
            fraction_covered = float(parts[-1])

            # Calculate mean depth (total coverage / gene length)
            # Actually, bedtools gives us covered bases, not depth
            # For depth, we need column 9 (number of overlaps)
            # Let me recalculate...
            # Column layout for bedtools coverage -a GFF -b BAM:
            # 1-9: GFF fields
            # 10: Number of features in B that overlapped A
            # 11: Number of bases in A that had coverage
            # 12: Length of A
            # 13: Fraction of A covered

            num_overlaps = int(parts[9]) if len(parts) > 9 else 0
            mean_depth = num_overlaps / gene_length if gene_length > 0 else 0
            percent_coverage = fraction_covered * 100

            coverage_data[gene_id] = {
                'mean_depth': mean_depth,
                'percent_coverage': percent_coverage,
                'bases_covered': bases_covered
            }

# Calculate TPM
# TPM formula:
# 1. Normalize by gene length: RPK = count / (gene_length / 1000)
# 2. Sum all RPKs
# 3. TPM = (RPK / sum_RPK) * 1,000,000

rpk_values = {}
for gene_id, count in counts.items():
    if count > 0 and gene_id in gene_lengths:
        length_kb = gene_lengths[gene_id] / 1000.0
        rpk_values[gene_id] = count / length_kb if length_kb > 0 else 0
    else:
        rpk_values[gene_id] = 0

sum_rpk = sum(rpk_values.values())

# Output combined table
print("gene_id\tgene_name\tread_count\trpm\ttpm\tmean_depth\tpercent_coverage\tgene_length")

for gene_id in sorted(counts.keys()):
    count = counts[gene_id]
    gene_name = gene_names.get(gene_id, gene_id)
    gene_length = gene_lengths.get(gene_id, 0)

    # RPM
    rpm = (count / total_for_rpm) * 1_000_000 if total_for_rpm > 0 else 0

    # TPM
    rpk = rpk_values.get(gene_id, 0)
    tpm = (rpk / sum_rpk) * 1_000_000 if sum_rpk > 0 else 0

    # Coverage
    cov = coverage_data.get(gene_id, {})
    mean_depth = cov.get('mean_depth', 0)
    percent_coverage = cov.get('percent_coverage', 0)

    print(f"{gene_id}\t{gene_name}\t{count}\t{rpm:.4f}\t{tpm:.4f}\t{mean_depth:.4f}\t{percent_coverage:.2f}\t{gene_length}")

# Print metadata to stderr
print(f"# Total read pairs: {total_pairs}", file=sys.stderr)
print(f"# RPM calculated {rpm_note}", file=sys.stderr)
print(f"# TPM calculated using standard formula", file=sys.stderr)

PYTHON_SCRIPT

STEP_END=$(date +%s)
STEP_DURATION=$((STEP_END - STEP_START))

log "✓ Combined abundance table created"
log "  Output: ${SAMPLE_NAME}_abundance.tsv"
log "  Columns: gene_id, gene_name, read_count, rpm, tpm, mean_depth, percent_coverage, gene_length"
log "  Wall time: ${STEP_DURATION}s"

PYTHON_WALL_TIME=${STEP_DURATION}

################################################################################
# Step 6: Create summary statistics with timing
################################################################################

log ""
log "Step 6: Creating summary statistics with timing..."

# Calculate total time
PIPELINE_TOTAL=$((BT2_WALL_TIME + SAMTOOLS_WALL_TIME + HTSEQ_WALL_TIME + BEDTOOLS_WALL_TIME + PYTHON_WALL_TIME))

cat > "${SAMPLE_DIR}/${SAMPLE_NAME}_stats.txt" << EOF
Sample: ${SAMPLE_NAME}
Date: $(date)

Alignment Statistics:
  Total read pairs: ${TOTAL_PAIRS}
  Aligned pairs: ${ALIGNED}
  Alignment rate: $(awk "BEGIN {printf \"%.2f\", (${ALIGNED}/${TOTAL_PAIRS})*100}")%

Gene Detection:
  Genes with reads: ${NUM_GENES}

Performance Metrics (Traditional Pipeline):
  Total pipeline time: ${PIPELINE_TOTAL}s ($(awk "BEGIN {printf \"%.2f\", ${PIPELINE_TOTAL}/60}")m)

  Step-by-step timing (wall clock time):
    1. Bowtie2 alignment:    ${BT2_WALL_TIME}s ($(awk "BEGIN {printf \"%.1f\", (${BT2_WALL_TIME}/${PIPELINE_TOTAL})*100}")%)
    2. SAM/BAM processing:   ${SAMTOOLS_WALL_TIME}s ($(awk "BEGIN {printf \"%.1f\", (${SAMTOOLS_WALL_TIME}/${PIPELINE_TOTAL})*100}")%)
    3. htseq-count:          ${HTSEQ_WALL_TIME}s ($(awk "BEGIN {printf \"%.1f\", (${HTSEQ_WALL_TIME}/${PIPELINE_TOTAL})*100}")%)
    4. bedtools coverage:    ${BEDTOOLS_WALL_TIME}s ($(awk "BEGIN {printf \"%.1f\", (${BEDTOOLS_WALL_TIME}/${PIPELINE_TOTAL})*100}")%)
    5. RPM/TPM calculation:  ${PYTHON_WALL_TIME}s ($(awk "BEGIN {printf \"%.1f\", (${PYTHON_WALL_TIME}/${PIPELINE_TOTAL})*100}")%)

  CPU time (user):
    - Bowtie2:     ${BT2_CPU_TIME}s
    - samtools:    ${SAMTOOLS_CPU_TIME}s
    - htseq-count: ${HTSEQ_CPU_TIME}s
    - bedtools:    ${BEDTOOLS_CPU_TIME}s

  Peak memory usage:
    - Bowtie2:     ${BT2_MEM} KB ($(awk "BEGIN {printf \"%.2f\", ${BT2_MEM}/1024/1024}")GB)
    - samtools:    ${SAMTOOLS_MEM} KB ($(awk "BEGIN {printf \"%.2f\", ${SAMTOOLS_MEM}/1024/1024}")GB)
    - htseq-count: ${HTSEQ_MEM} KB ($(awk "BEGIN {printf \"%.2f\", ${HTSEQ_MEM}/1024/1024}")GB)
    - bedtools:    ${BEDTOOLS_MEM} KB ($(awk "BEGIN {printf \"%.2f\", ${BEDTOOLS_MEM}/1024/1024}")GB)

  Threads used: ${THREADS}

Files Created:
  - ${SAMPLE_NAME}_counts.tsv         # Raw htseq-count output
  - ${SAMPLE_NAME}_abundance.tsv      # Combined metrics (counts, RPM, TPM, coverage)
  - ${SAMPLE_NAME}_coverage.bed       # Detailed bedtools coverage
  - ${SAMPLE_NAME}_stats.txt          # This file
  - ${SAMPLE_NAME}_timing.tsv         # Machine-readable timing data
  - ${SAMPLE_NAME}_*_time.txt         # Detailed /usr/bin/time output per step

Normalization:
  RPM: TODO - Verify if using pairs or reads as denominator!
  TPM: Standard formula (normalized by gene length and library size)
  Coverage: Fraction of gene bases covered by at least 1 read

Pipeline: bowtie2 -> htseq-count -> bedtools -> RPM/TPM calculation
EOF

# Create machine-readable timing file
cat > "${SAMPLE_DIR}/${SAMPLE_NAME}_timing.tsv" << EOF
sample_name	step	wall_time_sec	cpu_time_sec	memory_kb	memory_gb	threads	percent_of_total
${SAMPLE_NAME}	bowtie2	${BT2_WALL_TIME}	${BT2_CPU_TIME}	${BT2_MEM}	$(awk "BEGIN {printf \"%.3f\", ${BT2_MEM}/1024/1024}")	${THREADS}	$(awk "BEGIN {printf \"%.2f\", (${BT2_WALL_TIME}/${PIPELINE_TOTAL})*100}")
${SAMPLE_NAME}	samtools	${SAMTOOLS_WALL_TIME}	${SAMTOOLS_CPU_TIME}	${SAMTOOLS_MEM}	$(awk "BEGIN {printf \"%.3f\", ${SAMTOOLS_MEM}/1024/1024}")	${THREADS}	$(awk "BEGIN {printf \"%.2f\", (${SAMTOOLS_WALL_TIME}/${PIPELINE_TOTAL})*100}")
${SAMPLE_NAME}	htseq-count	${HTSEQ_WALL_TIME}	${HTSEQ_CPU_TIME}	${HTSEQ_MEM}	$(awk "BEGIN {printf \"%.3f\", ${HTSEQ_MEM}/1024/1024}")	1	$(awk "BEGIN {printf \"%.2f\", (${HTSEQ_WALL_TIME}/${PIPELINE_TOTAL})*100}")
${SAMPLE_NAME}	bedtools	${BEDTOOLS_WALL_TIME}	${BEDTOOLS_CPU_TIME}	${BEDTOOLS_MEM}	$(awk "BEGIN {printf \"%.3f\", ${BEDTOOLS_MEM}/1024/1024}")	1	$(awk "BEGIN {printf \"%.2f\", (${BEDTOOLS_WALL_TIME}/${PIPELINE_TOTAL})*100}")
${SAMPLE_NAME}	python	${PYTHON_WALL_TIME}	NA	NA	NA	1	$(awk "BEGIN {printf \"%.2f\", (${PYTHON_WALL_TIME}/${PIPELINE_TOTAL})*100}")
${SAMPLE_NAME}	TOTAL	${PIPELINE_TOTAL}	NA	NA	NA	${THREADS}	100.00
EOF

log "✓ Summary statistics created with timing metrics"
log "  - ${SAMPLE_NAME}_stats.txt (human-readable)"
log "  - ${SAMPLE_NAME}_timing.tsv (machine-readable)"

################################################################################
# Step 7: Cleanup
################################################################################

log ""
log "Step 7: Cleanup..."

if [[ "${KEEP_SAM}" == false ]]; then
    rm -f "${SAM_FILE}"
    log "✓ Removed SAM file"
fi

if [[ "${KEEP_BAM}" == false ]]; then
    rm -f "${BAM_FILE}" "${SORTED_BAM}" "${SORTED_BAM}.bai"
    log "✓ Removed BAM files"
else
    log "✓ Kept BAM files"
fi

################################################################################
# Done
################################################################################

log ""
log "========================================="
log "Sample ${SAMPLE_NAME} completed successfully"
log "========================================="
log "Output directory: ${SAMPLE_DIR}"
log "Main output: ${SAMPLE_NAME}_abundance.tsv"
log "Log file: ${LOGFILE}"
log ""

exit 0
