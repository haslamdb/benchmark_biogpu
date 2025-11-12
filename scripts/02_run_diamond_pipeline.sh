#!/bin/bash
################################################################################
# Traditional Pipeline: DIAMOND blastx + read quantification
# Outputs: read counts, RPM, coverage %, TPM (matching biogpu output)
# Uses same thresholds as biogpu: 85% identity, 50% coverage
################################################################################

set -euo pipefail

################################################################################
# Configuration
################################################################################

BENCHMARK_DIR="/home/david/projects/benchmark_biogpu"
DB_DIR="${BENCHMARK_DIR}/databases/traditional"
RESULTS_DIR="${BENCHMARK_DIR}/results/traditional"
LOG_DIR="${BENCHMARK_DIR}/logs"

# Databases
DIAMOND_DB="${DB_DIR}/amr_stress_protein"
PROTEIN_FASTA="/home/david/projects/biogpu/data/amr_stress_combined_db/protein.fasta"

# Executables (full paths to ensure they work in all environments)
DIAMOND="/home/david/miniforge3/envs/benchmark_biogpu/bin/diamond"
PYTHON="/home/david/miniforge3/envs/benchmark_biogpu/bin/python3"

# Processing parameters (matching biogpu)
THREADS=32
MIN_IDENTITY=85  # 85% identity
MIN_COVERAGE=50  # 50% query coverage

################################################################################
# Usage
################################################################################

usage() {
    cat << EOF
Usage: $0 --sample SAMPLE_NAME --r1 R1_FILE --r2 R2_FILE [options]

Required:
  --sample NAME    Sample name (e.g., N01_1_2)
  --r1 FILE        Path to R1 FASTQ file (can be gzipped)
  --r2 FILE        Path to R2 FASTQ file (can be gzipped)

Optional:
  --threads N      Number of threads (default: 32)

Pipeline: DIAMOND blastx (translated search)
Thresholds: ${MIN_IDENTITY}% identity, ${MIN_COVERAGE}% coverage
Output matching biogpu format

Example:
  $0 --sample N01_1_2 \\
     --r1 /path/to/sample_R1.fastq.gz \\
     --r2 /path/to/sample_R2.fastq.gz

EOF
    exit 1
}

################################################################################
# Parse arguments
################################################################################

SAMPLE_NAME=""
R1_FILE=""
R2_FILE=""

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
        -h|--help)
            usage
            ;;
        *)
            echo "Unknown option: $1"
            usage
            ;;
    esac
done

# Validate required arguments
if [[ -z "${SAMPLE_NAME}" ]] || [[ -z "${R1_FILE}" ]] || [[ -z "${R2_FILE}" ]]; then
    echo "Error: Missing required arguments"
    usage
fi

# Validate input files
if [[ ! -f "${R1_FILE}" ]]; then
    echo "Error: R1 file not found: ${R1_FILE}"
    exit 1
fi

if [[ ! -f "${R2_FILE}" ]]; then
    echo "Error: R2 file not found: ${R2_FILE}"
    exit 1
fi

################################################################################
# Setup
################################################################################

SAMPLE_DIR="${RESULTS_DIR}/${SAMPLE_NAME}"
mkdir -p "${SAMPLE_DIR}" "${LOG_DIR}"

LOG_FILE="${LOG_DIR}/${SAMPLE_NAME}_diamond_$(date +%Y%m%d_%H%M%S).log"

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "${LOG_FILE}"
}

log "========================================="
log "DIAMOND Pipeline: ${SAMPLE_NAME}"
log "========================================="
log "R1: ${R1_FILE}"
log "R2: ${R2_FILE}"
log "Output: ${SAMPLE_DIR}"
log "Threads: ${THREADS}"
log ""

PIPELINE_START=$(date +%s)

################################################################################
# Step 1: Count total reads
################################################################################

log "Step 1: Counting total reads..."

STEP_START=$(date +%s)

# Count reads from R1 file
if [[ "${R1_FILE}" == *.gz ]]; then
    TOTAL_READS=$(zcat "${R1_FILE}" | wc -l)
else
    TOTAL_READS=$(cat "${R1_FILE}" | wc -l)
fi

TOTAL_READS=$((TOTAL_READS / 4))  # FASTQ has 4 lines per read
TOTAL_PAIRS=${TOTAL_READS}

STEP_END=$(date +%s)
STEP_DURATION=$((STEP_END - STEP_START))

log "✓ Read counting completed"
log "  Total read pairs: ${TOTAL_PAIRS}"
log "  Wall time: ${STEP_DURATION}s"
log ""

################################################################################
# Step 2: Run DIAMOND blastx on R1 and R2 separately
################################################################################

log "Step 2: Running DIAMOND blastx (translated search)..."
log "  Parameters: ${MIN_IDENTITY}% identity, ${MIN_COVERAGE}% query coverage"

DIAMOND_OUTPUT="${SAMPLE_DIR}/${SAMPLE_NAME}_diamond.tsv"
DIAMOND_R1_OUTPUT="${SAMPLE_DIR}/${SAMPLE_NAME}_diamond_R1.tsv"
DIAMOND_R2_OUTPUT="${SAMPLE_DIR}/${SAMPLE_NAME}_diamond_R2.tsv"

STEP_START=$(date +%s)

# Run DIAMOND on R1
log "  Processing R1..."
/usr/bin/time -v -o "${SAMPLE_DIR}/${SAMPLE_NAME}_diamond_R1_time.txt" \
    "${DIAMOND}" blastx \
        --query "${R1_FILE}" \
        --db "${DIAMOND_DB}" \
        --out "${DIAMOND_R1_OUTPUT}" \
        --outfmt 6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore qcovhsp \
        --id ${MIN_IDENTITY} \
        --query-cover ${MIN_COVERAGE} \
        --threads ${THREADS} \
        --sensitive \
        --top 1 \
        > "${SAMPLE_DIR}/${SAMPLE_NAME}_diamond_R1.log" 2>&1

# Run DIAMOND on R2
log "  Processing R2..."
/usr/bin/time -v -o "${SAMPLE_DIR}/${SAMPLE_NAME}_diamond_R2_time.txt" \
    "${DIAMOND}" blastx \
        --query "${R2_FILE}" \
        --db "${DIAMOND_DB}" \
        --out "${DIAMOND_R2_OUTPUT}" \
        --outfmt 6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore qcovhsp \
        --id ${MIN_IDENTITY} \
        --query-cover ${MIN_COVERAGE} \
        --threads ${THREADS} \
        --sensitive \
        --top 1 \
        > "${SAMPLE_DIR}/${SAMPLE_NAME}_diamond_R2.log" 2>&1

# Combine R1 and R2 outputs
cat "${DIAMOND_R1_OUTPUT}" "${DIAMOND_R2_OUTPUT}" > "${DIAMOND_OUTPUT}"

STEP_END=$(date +%s)
STEP_DURATION=$((STEP_END - STEP_START))

# Extract timing from R1 (R2 should be similar)
DIAMOND_CPU_TIME=$(grep "User time" "${SAMPLE_DIR}/${SAMPLE_NAME}_diamond_R1_time.txt" | awk '{print $4}')
DIAMOND_SYS_TIME=$(grep "System time" "${SAMPLE_DIR}/${SAMPLE_NAME}_diamond_R1_time.txt" | awk '{print $4}')
DIAMOND_MEM=$(grep "Maximum resident set size" "${SAMPLE_DIR}/${SAMPLE_NAME}_diamond_R1_time.txt" | awk '{print $6}')

DIAMOND_R2_CPU_TIME=$(grep "User time" "${SAMPLE_DIR}/${SAMPLE_NAME}_diamond_R2_time.txt" | awk '{print $4}')
DIAMOND_R2_SYS_TIME=$(grep "System time" "${SAMPLE_DIR}/${SAMPLE_NAME}_diamond_R2_time.txt" | awk '{print $4}')

TOTAL_DIAMOND_CPU=$(awk "BEGIN {printf \"%.2f\", ${DIAMOND_CPU_TIME} + ${DIAMOND_R2_CPU_TIME}}")
TOTAL_DIAMOND_SYS=$(awk "BEGIN {printf \"%.2f\", ${DIAMOND_SYS_TIME} + ${DIAMOND_R2_SYS_TIME}}")

# Count alignments
TOTAL_ALIGNMENTS=$(wc -l < "${DIAMOND_OUTPUT}")

log "✓ DIAMOND blastx completed"
log "  Wall time: ${STEP_DURATION}s"
log "  CPU time: ${TOTAL_DIAMOND_CPU}s user + ${TOTAL_DIAMOND_SYS}s system"
log "  Max memory: ${DIAMOND_MEM} KB ($(awk "BEGIN {printf \"%.2f\", ${DIAMOND_MEM}/1024/1024}")GB)"
log "  Total alignments: ${TOTAL_ALIGNMENTS}"
log ""

# Store timing
DIAMOND_WALL_TIME=${STEP_DURATION}
DIAMOND_CPU_COMBINED=${TOTAL_DIAMOND_CPU}

################################################################################
# Step 3: Process DIAMOND output to gene-level counts
################################################################################

log "Step 3: Calculating gene-level counts, RPM, TPM, and coverage..."

STEP_START=$(date +%s)

# Use Python to process DIAMOND output
"${PYTHON}" << PYTHON_SCRIPT > "${SAMPLE_DIR}/${SAMPLE_NAME}_abundance.tsv"
import sys
from collections import defaultdict

# Parse variables from bash
sample_name = "${SAMPLE_NAME}"
diamond_output = "${DIAMOND_OUTPUT}"
protein_fasta = "${PROTEIN_FASTA}"
total_pairs = int("${TOTAL_PAIRS}")

# Read gene lengths from protein FASTA
gene_info = {}
current_gene = None
current_seq = []

with open(protein_fasta, 'r') as f:
    for line in f:
        if line.startswith('>'):
            # Save previous gene
            if current_gene:
                seq = ''.join(current_seq)
                gene_info[current_gene]['protein_length'] = len(seq)
                gene_info[current_gene]['dna_length'] = len(seq) * 3  # Codon = 3 nucleotides

            # Parse new gene
            # Format: >0|AAA16360.1|1|1|stxA2b|stxA2b||1|stxA2b|STX2|Shiga_toxin_Stx2b_subunit_A
            header = line.strip()[1:]
            parts = header.split('|')

            if len(parts) >= 5:
                gene_id = parts[0] + '|' + parts[1]  # e.g., "0|AAA16360.1"
                gene_name = parts[4] if parts[4] else parts[1]

                gene_info[gene_id] = {
                    'gene_name': gene_name,
                    'full_id': gene_id,
                    'header': header
                }
                current_gene = gene_id
                current_seq = []
        else:
            if current_gene:
                current_seq.append(line.strip())

# Save last gene
if current_gene:
    seq = ''.join(current_seq)
    gene_info[current_gene]['protein_length'] = len(seq)
    gene_info[current_gene]['dna_length'] = len(seq) * 3

# Count reads per gene and collect coverage info
gene_counts = defaultdict(int)
gene_positions = defaultdict(set)  # Track which positions are covered

with open(diamond_output, 'r') as f:
    for line in f:
        if not line.strip():
            continue

        parts = line.strip().split('\t')
        if len(parts) < 13:
            continue

        read_id = parts[0]
        subject_id = parts[1]  # Full gene ID from database
        pident = float(parts[2])
        qcovhsp = float(parts[12])

        # Extract gene_id (first two pipe-delimited fields)
        subject_parts = subject_id.split('|')
        if len(subject_parts) >= 2:
            gene_id = subject_parts[0] + '|' + subject_parts[1]
        else:
            gene_id = subject_id

        gene_counts[gene_id] += 1

        # Track covered positions (subject start to end)
        if len(parts) >= 10:
            sstart = int(parts[8])
            send = int(parts[9])
            for pos in range(min(sstart, send), max(sstart, send) + 1):
                gene_positions[gene_id].add(pos)

# Calculate RPM and TPM
# RPM = (reads / total_pairs) * 1,000,000
# TPM = (RPM / gene_length_kb) / sum(RPM / gene_length_kb) * 1,000,000

rpm_per_length = {}
for gene_id in gene_info:
    count = gene_counts.get(gene_id, 0)
    rpm = (count / total_pairs) * 1e6
    length_kb = gene_info[gene_id]['dna_length'] / 1000.0
    rpm_per_length[gene_id] = rpm / length_kb if length_kb > 0 else 0

sum_rpm_per_length = sum(rpm_per_length.values())

# Output table
print("gene_id\tgene_name\tread_count\trpm\ttpm\tmean_depth\tpercent_coverage\tgene_length")

for gene_id in sorted(gene_info.keys()):
    info = gene_info[gene_id]
    count = gene_counts.get(gene_id, 0)
    rpm = (count / total_pairs) * 1e6

    length_kb = info['dna_length'] / 1000.0
    if sum_rpm_per_length > 0:
        tpm = (rpm_per_length[gene_id] / sum_rpm_per_length) * 1e6
    else:
        tpm = 0.0

    # Calculate coverage
    protein_length = info['protein_length']
    if gene_id in gene_positions and protein_length > 0:
        covered_positions = len(gene_positions[gene_id])
        percent_coverage = (covered_positions / protein_length) * 100
        mean_depth = count / protein_length if protein_length > 0 else 0
    else:
        percent_coverage = 0.0
        mean_depth = 0.0

    print(f"{info['full_id']}\t{info['gene_name']}\t{count}\t{rpm:.4f}\t{tpm:.4f}\t{mean_depth:.4f}\t{percent_coverage:.2f}\t{info['dna_length']}")

PYTHON_SCRIPT

STEP_END=$(date +%s)
STEP_DURATION=$((STEP_END - STEP_START))

log "✓ Gene-level metrics calculated"
log "  Output: ${SAMPLE_NAME}_abundance.tsv"
log "  Wall time: ${STEP_DURATION}s"
log ""

PYTHON_WALL_TIME=${STEP_DURATION}

################################################################################
# Step 4: Create summary statistics
################################################################################

log "Step 4: Creating summary statistics with timing..."

PIPELINE_END=$(date +%s)
PIPELINE_TOTAL=$((PIPELINE_END - PIPELINE_START))

# Count genes with reads
GENES_WITH_READS=$(awk -F'\t' 'NR>1 && $3>0 {count++} END {print count}' "${SAMPLE_DIR}/${SAMPLE_NAME}_abundance.tsv")

# Create human-readable summary
cat > "${SAMPLE_DIR}/${SAMPLE_NAME}_stats.txt" << EOF
Sample: ${SAMPLE_NAME}
Date: $(date)

Alignment Statistics:
  Total read pairs: ${TOTAL_PAIRS}
  DIAMOND alignments: ${TOTAL_ALIGNMENTS}
  Alignment rate: $(awk "BEGIN {printf \"%.2f\", (${TOTAL_ALIGNMENTS}/(${TOTAL_PAIRS}*2))*100}")%

Gene Detection:
  Genes with reads: ${GENES_WITH_READS}

Performance Metrics (DIAMOND Pipeline):
  Total pipeline time: ${PIPELINE_TOTAL}s ($(awk "BEGIN {printf \"%.2f\", ${PIPELINE_TOTAL}/60}")m)

  Step-by-step timing (wall clock time):
    1. DIAMOND blastx (R1+R2): ${DIAMOND_WALL_TIME}s ($(awk "BEGIN {printf \"%.1f\", (${DIAMOND_WALL_TIME}/${PIPELINE_TOTAL})*100}")%)
    2. Gene quantification:     ${PYTHON_WALL_TIME}s ($(awk "BEGIN {printf \"%.1f\", (${PYTHON_WALL_TIME}/${PIPELINE_TOTAL})*100}")%)

  CPU time (user):
    - DIAMOND: ${DIAMOND_CPU_COMBINED}s

  Peak memory usage:
    - DIAMOND: ${DIAMOND_MEM} KB ($(awk "BEGIN {printf \"%.2f\", ${DIAMOND_MEM}/1024/1024}")GB)

  Threads used: ${THREADS}

Search Parameters:
  Min identity: ${MIN_IDENTITY}%
  Min query coverage: ${MIN_COVERAGE}%
  Search mode: Translated (blastx)

Files Created:
  - ${SAMPLE_NAME}_abundance.tsv       # Gene counts, RPM, TPM, coverage
  - ${SAMPLE_NAME}_diamond.tsv         # Raw DIAMOND output
  - ${SAMPLE_NAME}_stats.txt           # This file
  - ${SAMPLE_NAME}_timing.tsv          # Machine-readable timing data

Pipeline: DIAMOND blastx (translated search, matching biogpu parameters)
EOF

log "✓ Summary statistics created"
log ""

# Create machine-readable timing file
cat > "${SAMPLE_DIR}/${SAMPLE_NAME}_timing.tsv" << EOF
sample_name	step	wall_time_sec	cpu_time_sec	memory_kb	memory_gb	threads	percent_of_total
${SAMPLE_NAME}	diamond	${DIAMOND_WALL_TIME}	${DIAMOND_CPU_COMBINED}	${DIAMOND_MEM}	$(awk "BEGIN {printf \"%.3f\", ${DIAMOND_MEM}/1024/1024}")	${THREADS}	$(awk "BEGIN {printf \"%.2f\", (${DIAMOND_WALL_TIME}/${PIPELINE_TOTAL})*100}")
${SAMPLE_NAME}	python	${PYTHON_WALL_TIME}	NA	NA	NA	1	$(awk "BEGIN {printf \"%.2f\", (${PYTHON_WALL_TIME}/${PIPELINE_TOTAL})*100}")
${SAMPLE_NAME}	TOTAL	${PIPELINE_TOTAL}	NA	NA	NA	${THREADS}	100.00
EOF

################################################################################
# Step 5: Cleanup
################################################################################

log "Step 5: Cleanup..."

# Keep main outputs, remove intermediate files
rm -f "${DIAMOND_R1_OUTPUT}" "${DIAMOND_R2_OUTPUT}"

log "✓ Intermediate files removed"
log ""

################################################################################
# Done
################################################################################

log "========================================="
log "Sample ${SAMPLE_NAME} completed successfully"
log "========================================="
log "Output directory: ${SAMPLE_DIR}"
log "Main output: ${SAMPLE_NAME}_abundance.tsv"
log "Log file: ${LOG_FILE}"
log ""

cat << EOF
# Total read pairs: ${TOTAL_PAIRS}
# RPM calculated per million pairs
# TPM calculated using standard formula
# Search parameters: ${MIN_IDENTITY}% identity, ${MIN_COVERAGE}% coverage
EOF

exit 0
