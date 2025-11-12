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
- DIAMOND (translated search - primary method matching biogpu)
- bowtie2 (nucleotide alignment - optional alternative)
- htseq-count (read counting for bowtie2 pipeline)
- samtools (SAM/BAM manipulation)
- Python with pandas, numpy, scipy, matplotlib, seaborn
- bedtools, parallel, and other utilities

## Step 2: Setup Databases (10-15 minutes)

This script will:
1. Link biogpu AMR+stress databases (14,280 DNA genes, 17,099 protein genes)
2. Copy DNA FASTA for bowtie2 pipeline (optional)
3. Copy protein FASTA for DIAMOND pipeline (primary method)
4. Create GFF3 annotation file for htseq-count (bowtie2 pipeline)
5. Build DIAMOND index (.dmnd file)
6. Build bowtie2 index (optional)

```bash
./scripts/01_setup_databases.sh
```

**What it creates**:
- `databases/traditional/amr_stress_protein.fasta` - Protein sequences for DIAMOND
- `databases/traditional/amr_stress_protein.dmnd` - DIAMOND index
- `databases/traditional/amr_stress_dna.fasta` - DNA sequences for bowtie2 (optional)
- `databases/traditional/amr_stress_genes.gff3` - Annotations for htseq-count (bowtie2)
- `databases/traditional/amr_stress_bt2.*.bt2` - Bowtie2 index files (optional)

**Log file**: Check `logs/database_setup_YYYYMMDD_HHMMSS.log` for details

## Step 3: What's Next

Now you're ready to:

1. **Select test samples** (~10-20 NICU samples)
   - Stratified by body site, location, timepoint
   - Include varying AMR burden levels

2. **Run DIAMOND pipeline** (primary method)
   - Script: `scripts/02_run_diamond_pipeline.sh`
   - Uses DIAMOND blastx (translated search matching biogpu)
   - Same parameters: 85% identity, 50% coverage
   - Single sample: `./scripts/02_run_diamond_pipeline.sh --sample NAME --r1 R1.fq.gz --r2 R2.fq.gz`
   - Batch: `./scripts/03_batch_traditional_pipeline.sh`

3. **Run biogpu pipeline** with timing
   - Script: `scripts/04_run_biogpu_with_timing.sh`
   - Single sample: `./scripts/04_run_biogpu_with_timing.sh --sample NAME --r1 R1.fq.gz --r2 R2.fq.gz --gpu-id 0`
   - Batch: `./scripts/05_batch_biogpu_with_timing.sh --gpu-id 0`

4. **(Optional) Run bowtie2 pipeline** for nucleotide comparison
   - Uses bowtie2 → htseq-count → bedtools coverage
   - Script: `scripts/02_run_traditional_pipeline.sh`

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

### Traditional Pipeline (DIAMOND - Primary)
- **Method**: Translated search with DIAMOND blastx (matching biogpu's protein space)
- **Parameters**: 85% identity, 50% coverage (matching biogpu)
- **Mode**: --sensitive
- **To run**: `./scripts/02_run_diamond_pipeline.sh` or batch with `03_batch_traditional_pipeline.sh`

### Alternative: Bowtie2 Pipeline (Optional)
- **Method**: Nucleotide alignment with bowtie2
- **Parameters**: Default bowtie2 settings
- **Counting**: htseq-count
- **To run**: `./scripts/02_run_traditional_pipeline.sh`
- **Note**: Nucleotide-only alignment, less sensitive than translated search

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
