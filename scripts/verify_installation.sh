#!/bin/bash
################################################################################
# Verify Installation of All Required Tools
################################################################################

set -euo pipefail

echo "========================================="
echo "Verifying benchmark_biogpu Environment"
echo "========================================="
echo ""

ERRORS=0

# Function to check tool
check_tool() {
    local tool=$1
    local version_flag=${2:---version}

    if command -v ${tool} &> /dev/null; then
        echo "✓ ${tool} found"
        ${tool} ${version_flag} 2>&1 | head -1 | sed 's/^/  /'
    else
        echo "✗ ${tool} NOT FOUND"
        ERRORS=$((ERRORS + 1))
    fi
    echo ""
}

echo "Checking bioinformatics tools..."
echo "---------------------------------"
check_tool bowtie2 --version

# bwa doesn't have --version, just check if it exists
if command -v bwa &> /dev/null; then
    echo "✓ bwa found"
    bwa 2>&1 | grep "Version" | head -1 | sed 's/^/  /' || echo "  bwa is installed"
else
    echo "✗ bwa NOT FOUND"
    ERRORS=$((ERRORS + 1))
fi
echo ""

check_tool samtools --version
check_tool htseq-count --version
check_tool bedtools --version
check_tool featureCounts -v

echo "Checking Python and packages..."
echo "--------------------------------"
check_tool python --version

# Check Python packages
python << 'PYTHON_CHECK'
import sys
packages = [
    'numpy',
    'pandas',
    'scipy',
    'matplotlib',
    'seaborn',
    'Bio',  # biopython
    'pysam',
    'HTSeq'
]

print("Python packages:")
all_ok = True
for pkg in packages:
    try:
        if pkg == 'Bio':
            import Bio
            version = Bio.__version__
            actual_name = 'biopython'
        elif pkg == 'HTSeq':
            import HTSeq
            version = HTSeq.__version__
            actual_name = 'HTSeq'
        else:
            mod = __import__(pkg)
            version = getattr(mod, '__version__', 'unknown')
            actual_name = pkg

        print(f"  ✓ {actual_name}: {version}")
    except ImportError:
        print(f"  ✗ {pkg}: NOT FOUND")
        all_ok = False

sys.exit(0 if all_ok else 1)
PYTHON_CHECK

if [ $? -ne 0 ]; then
    ERRORS=$((ERRORS + 1))
fi

echo ""
echo "Checking system utilities..."
echo "----------------------------"
check_tool pigz --version
check_tool parallel --version

echo "==========================================="
if [ ${ERRORS} -eq 0 ]; then
    echo "✓ All dependencies installed successfully!"
    echo "==========================================="
    echo ""
    echo "Environment is ready to use."
    echo "Activate with: conda activate benchmark_biogpu"
    echo ""
    echo "Next steps:"
    echo "  1. Run database setup: ./scripts/01_setup_databases.sh"
    echo "  2. Process samples: ./scripts/03_batch_traditional_pipeline.sh"
    exit 0
else
    echo "✗ ${ERRORS} errors found. Please fix before proceeding."
    echo "==========================================="
    exit 1
fi
