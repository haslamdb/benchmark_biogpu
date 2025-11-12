#!/bin/bash
################################################################################
# Batch Process Multiple Samples Through Traditional Pipeline
################################################################################

set -euo pipefail

BENCHMARK_DIR="/home/david/projects/benchmark_biogpu"
SAMPLE_LIST="${BENCHMARK_DIR}/data/test_samples.csv"
LOG_DIR="${BENCHMARK_DIR}/logs"

mkdir -p "${LOG_DIR}"

BATCH_LOG="${LOG_DIR}/batch_traditional_$(date +%Y%m%d_%H%M%S).log"

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "${BATCH_LOG}"
}

log "========================================="
log "Batch Traditional Pipeline Processing"
log "========================================="
log "Sample list: ${SAMPLE_LIST}"
log ""

if [[ ! -f "${SAMPLE_LIST}" ]]; then
    log "ERROR: Sample list not found: ${SAMPLE_LIST}"
    log "Please create a CSV file with columns: sample_name,r1_path,r2_path"
    exit 1
fi

# Count samples
TOTAL_SAMPLES=$(tail -n +2 "${SAMPLE_LIST}" | wc -l)
log "Total samples to process: ${TOTAL_SAMPLES}"
log ""

# Process each sample
CURRENT=0
SUCCESS=0
FAILED=0
START_TIME=$(date +%s)

while IFS=, read -r sample_name r1_path r2_path; do
    # Skip header
    if [[ "${sample_name}" == "sample_name" ]] || [[ -z "${sample_name}" ]]; then
        continue
    fi

    CURRENT=$((CURRENT + 1))

    log "========================================="
    log "Processing sample ${CURRENT}/${TOTAL_SAMPLES}: ${sample_name}"
    log "========================================="

    # Run DIAMOND pipeline
    if "${BENCHMARK_DIR}/scripts/02_run_diamond_pipeline.sh" \
        --sample "${sample_name}" \
        --r1 "${r1_path}" \
        --r2 "${r2_path}" \
        --threads 32; then

        SUCCESS=$((SUCCESS + 1))
        log "✓ Sample ${sample_name} completed"
    else
        FAILED=$((FAILED + 1))
        log "✗ Sample ${sample_name} failed"
    fi

    log ""

    # Progress update every 5 samples
    if (( CURRENT % 5 == 0 )); then
        ELAPSED=$(($(date +%s) - START_TIME))
        AVG_TIME=$((ELAPSED / CURRENT))
        REMAINING=$((TOTAL_SAMPLES - CURRENT))
        EST_REMAINING=$((AVG_TIME * REMAINING))

        log "=== PROGRESS SUMMARY ==="
        log "Completed: ${CURRENT}/${TOTAL_SAMPLES} ($(( CURRENT * 100 / TOTAL_SAMPLES ))%)"
        log "Successful: ${SUCCESS}"
        log "Failed: ${FAILED}"
        log "Avg time per sample: ${AVG_TIME}s"
        log "Est. time remaining: $((EST_REMAINING / 3600))h $((EST_REMAINING % 3600 / 60))m"
        log "======================="
        log ""
    fi

done < "${SAMPLE_LIST}"

# Final summary
END_TIME=$(date +%s)
TOTAL_DURATION=$((END_TIME - START_TIME))

log ""
log "========================================="
log "BATCH PROCESSING COMPLETE"
log "========================================="
log "Total samples: ${TOTAL_SAMPLES}"
log "Successful: ${SUCCESS}"
log "Failed: ${FAILED}"
log "Total duration: $((TOTAL_DURATION / 3600))h $((TOTAL_DURATION % 3600 / 60))m $((TOTAL_DURATION % 60))s"
log "Average time per sample: $((TOTAL_DURATION / TOTAL_SAMPLES))s"
log ""
log "Results saved to: ${BENCHMARK_DIR}/results/traditional/"
log "Log file: ${BATCH_LOG}"
log ""

if [[ ${FAILED} -gt 0 ]]; then
    log "⚠ ${FAILED} samples failed. Check individual logs for details."
    exit 1
fi

log "✓ All samples processed successfully!"
exit 0
