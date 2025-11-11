#!/bin/bash
################################################################################
# Run BioGPU Pipeline with Comprehensive Timing
# Re-runs selected samples to capture CPU, GPU, and wall clock time
################################################################################

set -euo pipefail

BENCHMARK_DIR="/home/david/projects/benchmark_biogpu"
BIOGPU_DIR="/home/david/projects/biogpu"
OUTPUT_DIR="${BENCHMARK_DIR}/results/biogpu"
LOG_DIR="${BENCHMARK_DIR}/logs"
TEMP_DIR="${BENCHMARK_DIR}/temp"

mkdir -p "${OUTPUT_DIR}" "${LOG_DIR}" "${TEMP_DIR}"

# BioGPU executables
AMR_EXE="${BIOGPU_DIR}/build/runtime/kernels/genes/amr_detection"
FQ_EXE="${BIOGPU_DIR}/build/runtime/kernels/resistance/clean_resistance_pipeline"

# BioGPU databases
AMR_DB_DNA="${BIOGPU_DIR}/data/amr_stress_combined_db/dna.fasta"
AMR_DB_PROTEIN="${BIOGPU_DIR}/data/amr_stress_combined_db/protein.fasta"
AMR_PROTEIN_DB="${BIOGPU_DIR}/data/amr_stress_combined_protein_db"
FQ_DB_NUC="${BIOGPU_DIR}/data/integrated_clean_db/nucleotide"
FQ_DB_PROTEIN="${BIOGPU_DIR}/data/integrated_clean_db/protein"

# Parameters (matching original pipeline)
MIN_IDENTITY="0.85"
MIN_COVERAGE="0.50"

################################################################################
# Usage
################################################################################

usage() {
    cat << EOF
Usage: $0 --sample SAMPLE_NAME --r1 R1_FILE --r2 R2_FILE [options]

Required:
  --sample NAME    Sample name (e.g., N01_1_2)
  --r1 FILE        Path to R1 FASTQ file
  --r2 FILE        Path to R2 FASTQ file

Optional:
  --gpu-id N       GPU device ID (default: 0)
  --profile        Use nvprof for GPU profiling (requires CUDA toolkit)

Output files created in ${OUTPUT_DIR}/SAMPLE_NAME/:
  - SAMPLE_NAME_amr_abundance.tsv      # AMR gene abundances
  - SAMPLE_NAME_fq_mutations.tsv       # FQ resistance mutations
  - SAMPLE_NAME_timing.tsv             # Comprehensive timing data
  - SAMPLE_NAME_stats.txt              # Human-readable summary

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
GPU_ID=0
USE_NVPROF=false

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
        --gpu-id)
            GPU_ID="$2"
            shift 2
            ;;
        --profile)
            USE_NVPROF=true
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

################################################################################
# Setup
################################################################################

SAMPLE_DIR="${OUTPUT_DIR}/${SAMPLE_NAME}"
mkdir -p "${SAMPLE_DIR}"

LOGFILE="${LOG_DIR}/${SAMPLE_NAME}_biogpu_$(date +%Y%m%d_%H%M%S).log"

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "${LOGFILE}"
}

log "========================================="
log "BioGPU Pipeline with Timing: ${SAMPLE_NAME}"
log "========================================="
log "R1: ${R1_FILE}"
log "R2: ${R2_FILE}"
log "GPU ID: ${GPU_ID}"
log "Output: ${SAMPLE_DIR}"
log ""

# Create temporary CSV for this sample
TMP_CSV="${TEMP_DIR}/${SAMPLE_NAME}.csv"
echo "SampleName,FilePath,R1 file,R2 file" > "${TMP_CSV}"
echo "${SAMPLE_NAME},$(dirname ${R1_FILE}),$(basename ${R1_FILE}),$(basename ${R2_FILE})" >> "${TMP_CSV}"

# Verify files exist
if [[ ! -f "${R1_FILE}" ]]; then
    log "ERROR: R1 file not found: ${R1_FILE}"
    exit 1
fi

if [[ ! -f "${R2_FILE}" ]]; then
    log "ERROR: R2 file not found: ${R2_FILE}"
    exit 1
fi

################################################################################
# Step 1: AMR + Stress Gene Detection
################################################################################

log "Step 1: Running AMR + Stress gene detection..."

export CUDA_VISIBLE_DEVICES=${GPU_ID}

AMR_START=$(date +%s)

if [[ "${USE_NVPROF}" == true ]]; then
    log "  GPU profiling enabled (nvprof)"
    /usr/bin/time -v -o "${SAMPLE_DIR}/${SAMPLE_NAME}_amr_time.txt" \
        nvprof --log-file "${SAMPLE_DIR}/${SAMPLE_NAME}_amr_gpu_profile.log" \
        "${AMR_EXE}" \
            "${AMR_DB_DNA},${AMR_DB_PROTEIN}" \
            "${TMP_CSV}" \
            "${SAMPLE_DIR}" \
            --protein-db "${AMR_PROTEIN_DB}" \
            --no-merge \
            --min-identity "${MIN_IDENTITY}" \
            --min-coverage "${MIN_COVERAGE}" \
            > "${SAMPLE_DIR}/${SAMPLE_NAME}_amr.log" 2>&1
else
    /usr/bin/time -v -o "${SAMPLE_DIR}/${SAMPLE_NAME}_amr_time.txt" \
        "${AMR_EXE}" \
            "${AMR_DB_DNA},${AMR_DB_PROTEIN}" \
            "${TMP_CSV}" \
            "${SAMPLE_DIR}" \
            --protein-db "${AMR_PROTEIN_DB}" \
            --no-merge \
            --min-identity "${MIN_IDENTITY}" \
            --min-coverage "${MIN_COVERAGE}" \
            > "${SAMPLE_DIR}/${SAMPLE_NAME}_amr.log" 2>&1
fi

AMR_END=$(date +%s)
AMR_WALL_TIME=$((AMR_END - AMR_START))

# Extract timing from /usr/bin/time
AMR_CPU_TIME=$(grep "User time" "${SAMPLE_DIR}/${SAMPLE_NAME}_amr_time.txt" | awk '{print $4}')
AMR_SYS_TIME=$(grep "System time" "${SAMPLE_DIR}/${SAMPLE_NAME}_amr_time.txt" | awk '{print $4}')
AMR_MEM=$(grep "Maximum resident set size" "${SAMPLE_DIR}/${SAMPLE_NAME}_amr_time.txt" | awk '{print $6}')

log "✓ AMR + Stress detection completed"
log "  Wall time: ${AMR_WALL_TIME}s"
log "  CPU time: ${AMR_CPU_TIME}s user + ${AMR_SYS_TIME}s system"
log "  Max memory: ${AMR_MEM} KB ($(awk "BEGIN {printf \"%.2f\", ${AMR_MEM}/1024/1024}")GB)"

# Count detected genes
if [[ -f "${SAMPLE_DIR}/${SAMPLE_NAME}_amr_abundance.tsv" ]]; then
    AMR_GENES=$(tail -n +2 "${SAMPLE_DIR}/${SAMPLE_NAME}_amr_abundance.tsv" | wc -l)
    STRESS_GENES=$(tail -n +2 "${SAMPLE_DIR}/${SAMPLE_NAME}_amr_abundance.tsv" | grep "STRESS" | wc -l)
    log "  Detected: $((AMR_GENES - STRESS_GENES)) AMR genes, ${STRESS_GENES} stress genes"
fi

################################################################################
# Step 2: FQ Resistance Mutation Detection
################################################################################

log ""
log "Step 2: Running FQ resistance mutation detection..."

FQ_START=$(date +%s)

if [[ "${USE_NVPROF}" == true ]]; then
    log "  GPU profiling enabled (nvprof)"
    /usr/bin/time -v -o "${SAMPLE_DIR}/${SAMPLE_NAME}_fq_time.txt" \
        nvprof --log-file "${SAMPLE_DIR}/${SAMPLE_NAME}_fq_gpu_profile.log" \
        "${FQ_EXE}" \
            "${FQ_DB_NUC}" \
            "${FQ_DB_PROTEIN}" \
            --csv "${TMP_CSV}" \
            --no-bloom \
            --output-dir "${SAMPLE_DIR}" \
            > "${SAMPLE_DIR}/${SAMPLE_NAME}_fq.log" 2>&1
else
    /usr/bin/time -v -o "${SAMPLE_DIR}/${SAMPLE_NAME}_fq_time.txt" \
        "${FQ_EXE}" \
            "${FQ_DB_NUC}" \
            "${FQ_DB_PROTEIN}" \
            --csv "${TMP_CSV}" \
            --no-bloom \
            --output-dir "${SAMPLE_DIR}" \
            > "${SAMPLE_DIR}/${SAMPLE_NAME}_fq.log" 2>&1
fi

FQ_END=$(date +%s)
FQ_WALL_TIME=$((FQ_END - FQ_START))

# Extract timing
FQ_CPU_TIME=$(grep "User time" "${SAMPLE_DIR}/${SAMPLE_NAME}_fq_time.txt" | awk '{print $4}')
FQ_SYS_TIME=$(grep "System time" "${SAMPLE_DIR}/${SAMPLE_NAME}_fq_time.txt" | awk '{print $4}')
FQ_MEM=$(grep "Maximum resident set size" "${SAMPLE_DIR}/${SAMPLE_NAME}_fq_time.txt" | awk '{print $6}')

log "✓ FQ resistance detection completed"
log "  Wall time: ${FQ_WALL_TIME}s"
log "  CPU time: ${FQ_CPU_TIME}s user + ${FQ_SYS_TIME}s system"
log "  Max memory: ${FQ_MEM} KB ($(awk "BEGIN {printf \"%.2f\", ${FQ_MEM}/1024/1024}")GB)"

# Count mutations
if [[ -f "${SAMPLE_DIR}/${SAMPLE_NAME}/${SAMPLE_NAME}_mutations.tsv" ]]; then
    FQ_MUTATIONS=$(tail -n +2 "${SAMPLE_DIR}/${SAMPLE_NAME}/${SAMPLE_NAME}_mutations.tsv" | wc -l)
    log "  Detected: ${FQ_MUTATIONS} resistance mutations"
fi

################################################################################
# Step 3: Parse GPU timing (if nvprof was used)
################################################################################

GPU_TIME="NA"
GPU_MEM="NA"

if [[ "${USE_NVPROF}" == true ]]; then
    log ""
    log "Step 3: Parsing GPU profiling data..."

    # Extract GPU time from nvprof logs
    if [[ -f "${SAMPLE_DIR}/${SAMPLE_NAME}_amr_gpu_profile.log" ]]; then
        AMR_GPU_TIME=$(grep "GPU activities:" -A 1 "${SAMPLE_DIR}/${SAMPLE_NAME}_amr_gpu_profile.log" | tail -1 | awk '{print $2}')
        log "  AMR GPU time: ${AMR_GPU_TIME}"
    fi

    if [[ -f "${SAMPLE_DIR}/${SAMPLE_NAME}_fq_gpu_profile.log" ]]; then
        FQ_GPU_TIME=$(grep "GPU activities:" -A 1 "${SAMPLE_DIR}/${SAMPLE_NAME}_fq_gpu_profile.log" | tail -1 | awk '{print $2}')
        log "  FQ GPU time: ${FQ_GPU_TIME}"
    fi

    GPU_TIME="${AMR_GPU_TIME:-NA} + ${FQ_GPU_TIME:-NA}"
fi

################################################################################
# Step 4: Create Timing Summary
################################################################################

log ""
log "Step 4: Creating timing summary..."

PIPELINE_TOTAL=$((AMR_WALL_TIME + FQ_WALL_TIME))

cat > "${SAMPLE_DIR}/${SAMPLE_NAME}_stats.txt" << EOF
Sample: ${SAMPLE_NAME}
Date: $(date)
GPU Device: ${GPU_ID}

Detection Results:
  AMR genes: $((${AMR_GENES:-0} - ${STRESS_GENES:-0}))
  Stress genes: ${STRESS_GENES:-0}
  FQ mutations: ${FQ_MUTATIONS:-0}

Performance Metrics (BioGPU Pipeline):
  Total pipeline time: ${PIPELINE_TOTAL}s ($(awk "BEGIN {printf \"%.2f\", ${PIPELINE_TOTAL}/60}")m)

  Step-by-step timing (wall clock time):
    1. AMR + Stress detection: ${AMR_WALL_TIME}s ($(awk "BEGIN {printf \"%.1f\", (${AMR_WALL_TIME}/${PIPELINE_TOTAL})*100}")%)
    2. FQ resistance detection: ${FQ_WALL_TIME}s ($(awk "BEGIN {printf \"%.1f\", (${FQ_WALL_TIME}/${PIPELINE_TOTAL})*100}")%)

  CPU time (user + system):
    - AMR:  ${AMR_CPU_TIME}s + ${AMR_SYS_TIME}s = $(awk "BEGIN {print ${AMR_CPU_TIME} + ${AMR_SYS_TIME}}")s
    - FQ:   ${FQ_CPU_TIME}s + ${FQ_SYS_TIME}s = $(awk "BEGIN {print ${FQ_CPU_TIME} + ${FQ_SYS_TIME}}")s

  GPU time:
    - Total: ${GPU_TIME}

  Peak memory usage:
    - AMR:  ${AMR_MEM} KB ($(awk "BEGIN {printf \"%.2f\", ${AMR_MEM}/1024/1024}")GB)
    - FQ:   ${FQ_MEM} KB ($(awk "BEGIN {printf \"%.2f\", ${FQ_MEM}/1024/1024}")GB)

Pipeline Parameters:
  - Min identity: ${MIN_IDENTITY} (85%)
  - Min coverage: ${MIN_COVERAGE} (50%)
  - Search method: Translated (6-frame)
  - GPU accelerated: Yes

Files Created:
  - ${SAMPLE_NAME}_amr_abundance.tsv    # AMR + stress gene abundances
  - ${SAMPLE_NAME}_mutations.tsv        # FQ resistance mutations
  - ${SAMPLE_NAME}_stats.txt            # This file
  - ${SAMPLE_NAME}_timing.tsv           # Machine-readable timing
  - ${SAMPLE_NAME}_*_time.txt           # Detailed /usr/bin/time output
EOF

if [[ "${USE_NVPROF}" == true ]]; then
    echo "  - ${SAMPLE_NAME}_*_gpu_profile.log    # GPU profiling data" >> "${SAMPLE_DIR}/${SAMPLE_NAME}_stats.txt"
fi

# Create machine-readable timing file
cat > "${SAMPLE_DIR}/${SAMPLE_NAME}_timing.tsv" << EOF
sample_name	step	wall_time_sec	cpu_time_sec	gpu_time	memory_kb	memory_gb	device	percent_of_total
${SAMPLE_NAME}	amr_detection	${AMR_WALL_TIME}	${AMR_CPU_TIME}	${AMR_GPU_TIME:-NA}	${AMR_MEM}	$(awk "BEGIN {printf \"%.3f\", ${AMR_MEM}/1024/1024}")	GPU${GPU_ID}	$(awk "BEGIN {printf \"%.2f\", (${AMR_WALL_TIME}/${PIPELINE_TOTAL})*100}")
${SAMPLE_NAME}	fq_resistance	${FQ_WALL_TIME}	${FQ_CPU_TIME}	${FQ_GPU_TIME:-NA}	${FQ_MEM}	$(awk "BEGIN {printf \"%.3f\", ${FQ_MEM}/1024/1024}")	GPU${GPU_ID}	$(awk "BEGIN {printf \"%.2f\", (${FQ_WALL_TIME}/${PIPELINE_TOTAL})*100}")
${SAMPLE_NAME}	TOTAL	${PIPELINE_TOTAL}	NA	${GPU_TIME}	NA	NA	GPU${GPU_ID}	100.00
EOF

log "✓ Timing summary created"
log "  - ${SAMPLE_NAME}_stats.txt (human-readable)"
log "  - ${SAMPLE_NAME}_timing.tsv (machine-readable)"

################################################################################
# Cleanup
################################################################################

log ""
log "Cleanup..."
rm -f "${TMP_CSV}"
log "✓ Removed temporary files"

################################################################################
# Done
################################################################################

log ""
log "========================================="
log "Sample ${SAMPLE_NAME} completed successfully"
log "========================================="
log "Output directory: ${SAMPLE_DIR}"
log "Total time: ${PIPELINE_TOTAL}s ($(awk "BEGIN {printf \"%.2f\", ${PIPELINE_TOTAL}/60}")m)"
log "Log file: ${LOGFILE}"
log ""

exit 0
