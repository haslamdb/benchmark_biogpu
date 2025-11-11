# Quick Start Guide - BioGPU Benchmarking

## What's Ready

The benchmark project is set up and ready to go! Here's what we've created:

### Project Structure
```
benchmark_biogpu/
├── README.md                    # Project overview
├── QUICKSTART.md               # This file
├── environment.yml             # Conda environment with tools
├── scripts/
│   ├── 01_setup_databases.sh   # Database setup (READY TO RUN)
│   └── create_gff3_from_fasta.py  # GFF3 generator
├── docs/
│   ├── SETUP_GUIDE.md          # Detailed setup instructions
│   └── PROJECT_STATUS.md       # Current status
├── data/                       # Test samples (to be added)
├── results/
│   ├── biogpu/                 # BioGPU results
│   └── traditional/            # Traditional pipeline results
├── logs/                       # Processing logs
└── databases/
    ├── biogpu/                 # Links to biogpu DBs
    └── traditional/            # Bowtie2 indexed DBs
```

## Step 1: Create Environment (5 minutes)

```bash
cd /home/david/projects/benchmark_biogpu

# Create conda environment with all required tools
mamba env create -f environment.yml

# Or if using conda
conda env create -f environment.yml

# Activate
conda activate benchmark_biogpu
```

**Tools included**:
- bowtie2 (alignment)
- htseq-count (read counting - same as your existing scripts)
- samtools (SAM/BAM manipulation)
- Python with pandas, numpy, scipy, matplotlib, seaborn
- bedtools, parallel, and other utilities

## Step 2: Setup Databases (10-15 minutes)

This script will:
1. Link biogpu AMR+stress databases (14,280 genes)
2. Copy DNA FASTA for traditional alignment
3. Create GFF3 annotation file for htseq-count
4. Build bowtie2 index

```bash
./scripts/01_setup_databases.sh
```

**What it creates**:
- `databases/traditional/amr_stress_dna.fasta` - DNA sequences (16 MB)
- `databases/traditional/amr_stress_genes.gff3` - Annotations for htseq-count
- `databases/traditional/amr_stress_bt2.*.bt2` - Bowtie2 index files

**Log file**: Check `logs/database_setup_YYYYMMDD_HHMMSS.log` for details

## Step 3: What's Next

Now you're ready to:

1. **Select test samples** (~10-20 NICU samples)
   - Stratified by body site, location, timepoint
   - Include varying AMR burden levels

2. **Adapt your bowtie2 + htseq-count workflow**
   - We found your existing scripts in `/home/david/projects/ShellScripts/`
   - They use: `bowtie2` → `htseq-count` → `bedtools coverage`
   - Need to adapt for our test samples and new database

3. **Run both pipelines**
   - Extract biogpu results for test samples
   - Run traditional pipeline on same samples

4. **Compare results**
   - Correlation analysis
   - Gene detection overlap
   - Investigate discrepancies

## Key Database Info

**BioGPU AMR+Stress Database**:
- Source: `/home/david/projects/biogpu/data/amr_stress_combined_db/`
- Total genes: 14,280 (DNA), 17,099 (protein)
- AMR genes: ~9,257 (DNA), ~9,538 (protein)
- Stress genes: ~5,023 (DNA), ~7,561 (stress)
- Categories:
  - Antibiotic resistance (CARD, ResFinder)
  - Oxidative stress (catalases, SODs, peroxiredoxins)
  - Heat shock (chaperones, proteases)
  - DNA repair (SOS response, recombination)
  - Efflux systems
  - And more...

**FASTA Header Format**:
```
>index|protein_id|nucleotide_id|?|?|gene_name|gene_name2|description coordinates
```

Example:
```
>0|AAA16360.1|L11078.1|1|1|stxA2b|stxA2b|Shiga_toxin_Stx2b_subunit_A L11078.1:177-1136
```

**GFF3 Format** (for htseq-count):
```
0|AAA16360.1  biogpu  gene  1  960  .  +  .  ID=0|AAA16360.1;Name=stxA2b;gene=stxA2b;locus_tag=AAA16360.1
```

## Comparison Strategy

### BioGPU Pipeline
- **Method**: Translated search (6-frame) in protein space
- **Parameters**: 85% identity, 50% coverage
- **Sensitivity**: Detects divergent sequences
- **Already run**: Results in `/fastpool/analysis/nicu_amr_stress/`

### Traditional Pipeline
- **Method**: Nucleotide alignment with bowtie2
- **Parameters**: Default bowtie2 settings
- **Counting**: htseq-count (same as your existing workflow)
- **To run**: On same samples as biogpu

### What We Expect
- **High correlation** for conserved genes (r > 0.90)
- **More genes detected by biogpu** (divergent sequences)
- **Similar quantification** for abundant, conserved genes

## Troubleshooting

**Issue**: `conda: command not found`
- **Solution**: Load conda first: `module load conda` or check `~/.bashrc`

**Issue**: `bowtie2-build` fails
- **Solution**: Ensure enough memory (~4-8 GB) and disk space (~1 GB)

**Issue**: GFF3 doesn't match SAM alignments
- **Solution**: Check that bowtie2 uses same sequence IDs as GFF3
- Our script handles this automatically by parsing FASTA headers

**Issue**: htseq-count reports "no features"
- **Solution**: Verify GFF3 seqnames match SAM reference names
- Use `samtools view sample.bam | head` to check SAM reference names
- Use `grep "^[^#]" genes.gff3 | cut -f1 | sort -u` to check GFF3 seqnames

## Questions?

- **Project README**: `README.md` - Overview and background
- **Setup Guide**: `docs/SETUP_GUIDE.md` - Detailed instructions
- **Project Status**: `docs/PROJECT_STATUS.md` - Current progress
- **BioGPU Docs**: `/home/david/projects/biogpu/docs/STRESS_GENE_INTEGRATION_2025-11-01.md`

## Ready to Go!

Everything is set up. Just need to:
1. Create conda environment
2. Run database setup script
3. Select test samples
4. You're ready to benchmark!

Let me know when you're ready for the next steps!
