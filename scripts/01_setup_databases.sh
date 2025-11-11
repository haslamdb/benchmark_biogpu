#!/bin/bash
################################################################################
# Setup Databases for BioGPU Benchmarking
# Creates bowtie2 index and GFF3 file from biogpu AMR+stress database
################################################################################

set -euo pipefail

BENCHMARK_DIR="/home/david/projects/benchmark_biogpu"
BIOGPU_DATA="/home/david/projects/biogpu/data"

# Output directories
DB_DIR="${BENCHMARK_DIR}/databases"
BIOGPU_DB_DIR="${DB_DIR}/biogpu"
TRAD_DB_DIR="${DB_DIR}/traditional"
LOG_DIR="${BENCHMARK_DIR}/logs"

mkdir -p "${BIOGPU_DB_DIR}" "${TRAD_DB_DIR}" "${LOG_DIR}"

LOGFILE="${LOG_DIR}/database_setup_$(date +%Y%m%d_%H%M%S).log"

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "${LOGFILE}"
}

log "========================================="
log "BioGPU Benchmark Database Setup"
log "========================================="

################################################################################
# Step 1: Link BioGPU databases
################################################################################

log "Step 1: Linking BioGPU databases..."

if [ ! -L "${BIOGPU_DB_DIR}/amr_stress_combined_db" ]; then
    ln -s "${BIOGPU_DATA}/amr_stress_combined_db" "${BIOGPU_DB_DIR}/amr_stress_combined_db"
    log "✓ Linked AMR+Stress combined database"
else
    log "✓ AMR+Stress database already linked"
fi

if [ ! -L "${BIOGPU_DB_DIR}/integrated_clean_db" ]; then
    ln -s "${BIOGPU_DATA}/integrated_clean_db" "${BIOGPU_DB_DIR}/integrated_clean_db"
    log "✓ Linked FQ resistance database"
else
    log "✓ FQ resistance database already linked"
fi

################################################################################
# Step 2: Copy and prepare DNA FASTA for bowtie2
################################################################################

log ""
log "Step 2: Preparing DNA FASTA for bowtie2..."

DNA_FASTA="${BIOGPU_DATA}/amr_stress_combined_db/dna.fasta"
TRAD_FASTA="${TRAD_DB_DIR}/amr_stress_dna.fasta"

if [ ! -f "${TRAD_FASTA}" ]; then
    log "Copying DNA FASTA to traditional database directory..."
    cp "${DNA_FASTA}" "${TRAD_FASTA}"
    log "✓ DNA FASTA copied ($(grep -c '^>' ${TRAD_FASTA}) sequences)"
else
    log "✓ DNA FASTA already exists"
fi

################################################################################
# Step 3: Create GFF3 file for htseq-count
################################################################################

log ""
log "Step 3: Creating GFF3 annotation file for htseq-count..."

GFF3_FILE="${TRAD_DB_DIR}/amr_stress_genes.gff3"

if [ ! -f "${GFF3_FILE}" ]; then
    log "Generating GFF3 from FASTA headers..."

    # Use Python script to create GFF3
    python3 "${BENCHMARK_DIR}/scripts/create_gff3_from_fasta.py" \
        "${TRAD_FASTA}" \
        "${GFF3_FILE}" \
        >> "${LOGFILE}" 2>&1

    # Count entries
    NUM_GENES=$(grep -c "^[^#]" "${GFF3_FILE}")
    log "✓ GFF3 file created with ${NUM_GENES} gene entries"
else
    log "✓ GFF3 file already exists"
fi

################################################################################
# Step 4: Build bowtie2 index
################################################################################

log ""
log "Step 4: Building bowtie2 index..."

BT2_INDEX="${TRAD_DB_DIR}/amr_stress_bt2"

if [ ! -f "${BT2_INDEX}.1.bt2" ]; then
    log "Building bowtie2 index (this may take several minutes)..."

    START_TIME=$(date +%s)
    bowtie2-build \
        --threads 12 \
        "${TRAD_FASTA}" \
        "${BT2_INDEX}" \
        >> "${LOGFILE}" 2>&1

    END_TIME=$(date +%s)
    DURATION=$((END_TIME - START_TIME))

    log "✓ Bowtie2 index built in ${DURATION}s"
    log "  Index files: ${BT2_INDEX}.*.bt2"
else
    log "✓ Bowtie2 index already exists"
fi

################################################################################
# Step 5: Verify setup
################################################################################

log ""
log "Step 5: Verifying database setup..."

# Check all required files exist
ERRORS=0

if [ ! -f "${TRAD_FASTA}" ]; then
    log "ERROR: DNA FASTA not found: ${TRAD_FASTA}"
    ERRORS=$((ERRORS + 1))
fi

if [ ! -f "${GFF3_FILE}" ]; then
    log "ERROR: GFF3 file not found: ${GFF3_FILE}"
    ERRORS=$((ERRORS + 1))
fi

if [ ! -f "${BT2_INDEX}.1.bt2" ]; then
    log "ERROR: Bowtie2 index not found: ${BT2_INDEX}"
    ERRORS=$((ERRORS + 1))
fi

if [ ${ERRORS} -eq 0 ]; then
    log "✓ All database files verified"
else
    log "ERROR: ${ERRORS} issues found. Check log: ${LOGFILE}"
    exit 1
fi

################################################################################
# Summary
################################################################################

log ""
log "========================================="
log "Database Setup Complete"
log "========================================="
log "DNA FASTA: ${TRAD_FASTA}"
log "GFF3 file: ${GFF3_FILE}"
log "Bowtie2 index: ${BT2_INDEX}"
log ""
log "Gene count: $(grep -c '^>' ${TRAD_FASTA})"
log "  AMR genes: ~9,257"
log "  Stress genes: ~5,023"
log "  Total: ~14,280"
log ""
log "Log file: ${LOGFILE}"
log "========================================="

exit 0
