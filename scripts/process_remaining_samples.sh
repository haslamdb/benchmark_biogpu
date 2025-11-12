#!/bin/bash
################################################################################
# Process All Remaining NICU Samples Through DIAMOND Pipeline
# Reads from data/remaining_samples.csv and processes each sample
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
MAX_PARALLEL=${MAX_PARALLEL:-1}  # Number of samples to run in parallel (default: 1)

################################################################################
# Setup
################################################################################

mkdir -p "${LOG_DIR}" "${RESULTS_DIR}"

BATCH_LOG="${LOG_DIR}/batch_remaining_$(date +%Y%m%d_%H%M%S).log"

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "${BATCH_LOG}"
}

log "========================================="
log "Batch Processing Remaining NICU Samples"
log "========================================="
log "Sample list: ${SAMPLE_LIST}"
log "Max parallel: ${MAX_PARALLEL}"
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

################################################################################
# Count samples
################################################################################

TOTAL_SAMPLES=$(tail -n +2 "${SAMPLE_LIST}" | wc -l)
log "Total samples to process: ${TOTAL_SAMPLES}"
log ""

################################################################################
# Process each sample
################################################################################

PROCESSED=0
FAILED=0
SKIPPED=0

START_TIME=$(date +%s)

# Read sample list (skip header)
tail -n +2 "${SAMPLE_LIST}" | while IFS=, read -r SAMPLE_NAME R1_PATH R2_PATH; do
    PROCESSED=$((PROCESSED + 1))

    log "----------------------------------------"
    log "Processing sample ${PROCESSED}/${TOTAL_SAMPLES}: ${SAMPLE_NAME}"
    log "----------------------------------------"

    # Check if already processed
    SAMPLE_DIR="${RESULTS_DIR}/${SAMPLE_NAME}"
    ABUNDANCE_FILE="${SAMPLE_DIR}/${SAMPLE_NAME}_abundance.tsv"

    if [ -f "${ABUNDANCE_FILE}" ]; then
        log "✓ SKIPPED: ${SAMPLE_NAME} (already processed)"
        SKIPPED=$((SKIPPED + 1))
        continue
    fi

    # Run DIAMOND pipeline
    SAMPLE_START=$(date +%s)

    if bash "${DIAMOND_PIPELINE}" \
        --sample "${SAMPLE_NAME}" \
        --r1 "${R1_PATH}" \
        --r2 "${R2_PATH}" \
        >> "${BATCH_LOG}" 2>&1; then

        SAMPLE_END=$(date +%s)
        SAMPLE_DURATION=$((SAMPLE_END - SAMPLE_START))

        log "✓ SUCCESS: ${SAMPLE_NAME} (${SAMPLE_DURATION}s)"

        # Estimate time remaining
        ELAPSED=$((SAMPLE_END - START_TIME))
        REMAINING=$((TOTAL_SAMPLES - PROCESSED))
        if [ ${PROCESSED} -gt 0 ]; then
            AVG_TIME=$((ELAPSED / PROCESSED))
            EST_REMAINING=$((AVG_TIME * REMAINING))
            EST_HOURS=$((EST_REMAINING / 3600))
            EST_MINS=$(((EST_REMAINING % 3600) / 60))
            log "  Progress: ${PROCESSED}/${TOTAL_SAMPLES} (${EST_HOURS}h ${EST_MINS}m remaining)"
        fi
    else
        SAMPLE_END=$(date +%s)
        SAMPLE_DURATION=$((SAMPLE_END - SAMPLE_START))

        log "✗ FAILED: ${SAMPLE_NAME} (${SAMPLE_DURATION}s)"
        log "  Check log: ${LOG_DIR}/${SAMPLE_NAME}_diamond_*.log"
        FAILED=$((FAILED + 1))

        # Don't exit on failure, continue with next sample
    fi

    log ""
done

################################################################################
# Summary
################################################################################

END_TIME=$(date +%s)
TOTAL_DURATION=$((END_TIME - START_TIME))
HOURS=$((TOTAL_DURATION / 3600))
MINS=$(((TOTAL_DURATION % 3600) / 60))
SECS=$((TOTAL_DURATION % 60))

log "========================================="
log "Batch Processing Complete"
log "========================================="
log "Total samples: ${TOTAL_SAMPLES}"
log "Successfully processed: ${PROCESSED}"
log "Skipped (already done): ${SKIPPED}"
log "Failed: ${FAILED}"
log "Total time: ${HOURS}h ${MINS}m ${SECS}s"
log ""
log "Results directory: ${RESULTS_DIR}"
log "Batch log: ${BATCH_LOG}"
log "========================================="

if [ ${FAILED} -gt 0 ]; then
    log ""
    log "WARNING: ${FAILED} samples failed. Check logs for details."
    exit 1
fi

exit 0
