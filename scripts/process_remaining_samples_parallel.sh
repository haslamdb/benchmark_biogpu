#!/bin/bash
################################################################################
# Process All Remaining NICU Samples Through DIAMOND Pipeline (PARALLEL)
# Uses GNU parallel for efficient batch processing
################################################################################

set -euo pipefail

################################################################################
# Configuration
################################################################################

BENCHMARK_DIR="/home/david/projects/benchmark_biogpu"
SCRIPT_DIR="${BENCHMARK_DIR}/scripts"
SAMPLE_LIST="${BENCHMARK_DIR}/data/remaining_samples.csv"
LOG_DIR="${BENCHMARK_DIR}/logs"
RESULTS_DIR="${BENCHMARK_DIR}/results/traditional"

# Pipeline script
DIAMOND_PIPELINE="${SCRIPT_DIR}/02_run_diamond_pipeline.sh"

# Parallelization settings
# DIAMOND uses 32 threads per sample, so limit concurrent samples based on CPU count
# Recommended: total_cores / 32 = max parallel jobs
# Example: 128 cores / 32 = 4 parallel samples
MAX_JOBS=${MAX_JOBS:-4}  # Number of samples to run in parallel

################################################################################
# Setup
################################################################################

mkdir -p "${LOG_DIR}" "${RESULTS_DIR}"

BATCH_LOG="${LOG_DIR}/batch_remaining_parallel_$(date +%Y%m%d_%H%M%S).log"

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "${BATCH_LOG}"
}

log "========================================="
log "Batch Processing Remaining NICU Samples (PARALLEL)"
log "========================================="
log "Sample list: ${SAMPLE_LIST}"
log "Max parallel jobs: ${MAX_JOBS}"
log "Log file: ${BATCH_LOG}"
log ""

################################################################################
# Validate inputs
################################################################################

if [ ! -f "${SAMPLE_LIST}" ]; then
    log "ERROR: Sample list not found: ${SAMPLE_LIST}"
    log "Run: python scripts/identify_remaining_samples.py --verify-files"
    exit 1
fi

if [ ! -f "${DIAMOND_PIPELINE}" ]; then
    log "ERROR: DIAMOND pipeline script not found: ${DIAMOND_PIPELINE}"
    exit 1
fi

# Check if GNU parallel is available
if ! command -v parallel &> /dev/null; then
    log "ERROR: GNU parallel not found. Install with: conda install -c conda-forge parallel"
    log "Alternatively, use the non-parallel version: scripts/process_remaining_samples.sh"
    exit 1
fi

################################################################################
# Count samples
################################################################################

TOTAL_SAMPLES=$(tail -n +2 "${SAMPLE_LIST}" | wc -l)
log "Total samples to process: ${TOTAL_SAMPLES}"

# Count already processed
ALREADY_PROCESSED=0
while IFS=, read -r SAMPLE_NAME R1_PATH R2_PATH; do
    ABUNDANCE_FILE="${RESULTS_DIR}/${SAMPLE_NAME}/${SAMPLE_NAME}_abundance.tsv"
    if [ -f "${ABUNDANCE_FILE}" ]; then
        ALREADY_PROCESSED=$((ALREADY_PROCESSED + 1))
    fi
done < <(tail -n +2 "${SAMPLE_LIST}")

REMAINING=$((TOTAL_SAMPLES - ALREADY_PROCESSED))
log "Already processed: ${ALREADY_PROCESSED}"
log "Remaining to process: ${REMAINING}"
log ""

if [ ${REMAINING} -eq 0 ]; then
    log "✓ All samples already processed!"
    exit 0
fi

################################################################################
# Estimate runtime
################################################################################

AVG_TIME_PER_SAMPLE=327  # seconds (from benchmarking)
TOTAL_TIME_SERIAL=$((REMAINING * AVG_TIME_PER_SAMPLE))
TOTAL_TIME_PARALLEL=$((TOTAL_TIME_SERIAL / MAX_JOBS))

HOURS_PARALLEL=$((TOTAL_TIME_PARALLEL / 3600))
MINS_PARALLEL=$(((TOTAL_TIME_PARALLEL % 3600) / 60))

log "Estimated time:"
log "  Serial (1 at a time): $(((TOTAL_TIME_SERIAL / 3600)))h $(((TOTAL_TIME_SERIAL % 3600) / 60))m"
log "  Parallel (${MAX_JOBS} at a time): ${HOURS_PARALLEL}h ${MINS_PARALLEL}m"
log ""

################################################################################
# Process samples in parallel
################################################################################

log "Starting parallel processing..."
log "Press Ctrl+C to cancel"
log ""

START_TIME=$(date +%s)

# Create a function to process one sample
process_sample() {
    local SAMPLE_NAME=$1
    local R1_PATH=$2
    local R2_PATH=$3
    local RESULTS_DIR=$4
    local DIAMOND_PIPELINE=$5
    local LOG_DIR=$6

    # Check if already processed
    ABUNDANCE_FILE="${RESULTS_DIR}/${SAMPLE_NAME}/${SAMPLE_NAME}_abundance.tsv"
    if [ -f "${ABUNDANCE_FILE}" ]; then
        echo "SKIP: ${SAMPLE_NAME} (already processed)"
        return 0
    fi

    # Run DIAMOND pipeline
    SAMPLE_START=$(date +%s)

    if bash "${DIAMOND_PIPELINE}" \
        --sample "${SAMPLE_NAME}" \
        --r1 "${R1_PATH}" \
        --r2 "${R2_PATH}" \
        >> "${LOG_DIR}/${SAMPLE_NAME}_parallel.log" 2>&1; then

        SAMPLE_END=$(date +%s)
        SAMPLE_DURATION=$((SAMPLE_END - SAMPLE_START))
        echo "SUCCESS: ${SAMPLE_NAME} (${SAMPLE_DURATION}s)"
        return 0
    else
        SAMPLE_END=$(date +%s)
        SAMPLE_DURATION=$((SAMPLE_END - SAMPLE_START))
        echo "FAILED: ${SAMPLE_NAME} (${SAMPLE_DURATION}s)"
        return 1
    fi
}

# Export function and variables for parallel
export -f process_sample
export RESULTS_DIR
export DIAMOND_PIPELINE
export LOG_DIR

# Run parallel processing
# --jobs: number of parallel jobs
# --colsep: CSV delimiter
# --bar: show progress bar
# --results: save individual outputs
# --joblog: log all jobs
tail -n +2 "${SAMPLE_LIST}" | \
    parallel \
        --jobs ${MAX_JOBS} \
        --colsep ',' \
        --bar \
        --joblog "${LOG_DIR}/parallel_joblog_$(date +%Y%m%d_%H%M%S).txt" \
        --results "${LOG_DIR}/parallel_results" \
        process_sample {1} {2} {3} "${RESULTS_DIR}" "${DIAMOND_PIPELINE}" "${LOG_DIR}" \
    | tee -a "${BATCH_LOG}"

################################################################################
# Summary
################################################################################

END_TIME=$(date +%s)
TOTAL_DURATION=$((END_TIME - START_TIME))
HOURS=$((TOTAL_DURATION / 3600))
MINS=$(((TOTAL_DURATION % 3600) / 60))
SECS=$((TOTAL_DURATION % 60))

log ""
log "========================================="
log "Batch Processing Complete"
log "========================================="

# Count results
COMPLETED=0
for SAMPLE_DIR in "${RESULTS_DIR}"/*/; do
    SAMPLE=$(basename "${SAMPLE_DIR}")
    if [ -f "${SAMPLE_DIR}/${SAMPLE}_abundance.tsv" ]; then
        COMPLETED=$((COMPLETED + 1))
    fi
done

log "Total samples in list: ${TOTAL_SAMPLES}"
log "Successfully completed (total): ${COMPLETED}"
log "Processing time: ${HOURS}h ${MINS}m ${SECS}s"
log "Average time per sample: $((TOTAL_DURATION / REMAINING))s"
log ""
log "Results directory: ${RESULTS_DIR}"
log "Batch log: ${BATCH_LOG}"
log "Job log: ${LOG_DIR}/parallel_joblog_*.txt"
log "========================================="

# Check for failures
FAILED=$(grep -c "FAILED:" "${BATCH_LOG}" || true)
if [ ${FAILED} -gt 0 ]; then
    log ""
    log "WARNING: ${FAILED} samples failed. Check logs for details:"
    log "  grep FAILED ${BATCH_LOG}"
    exit 1
fi

log ""
log "✓ All samples processed successfully!"

exit 0
