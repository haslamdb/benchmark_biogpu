#!/usr/bin/env python3
"""
DIAMOND Resistome Analysis
Adapted from ZCH_UCMC_Manuscript resistome analysis for DIAMOND alignment results

Performs:
1. Create AMR gene matrix (samples × genes with RPM values)
2. PCA analysis stratified by body site
3. Hierarchical clustering
4. Diversity metrics (Shannon, Simpson)
5. Generate visualization figures
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import pdist, squareform
from scipy.stats import entropy
import warnings
warnings.filterwarnings('ignore')

# Define paths
PROJECT_ROOT = Path("/home/david/projects/benchmark_biogpu")
DATA_DIR = PROJECT_ROOT / "data"
RESULTS_DIR = PROJECT_ROOT / "results" / "analysis" / "diamond_resistome"
FIGURES_DIR = PROJECT_ROOT / "figures" / "diamond_resistome"

# Create directories
RESULTS_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

# Set plotting style
sns.set_style("whitegrid")
sns.set_context("talk")

def load_data():
    """Load metadata and AMR data"""
    print("Loading data...")

    # Load updated metadata (from corrected ZCH_UCMC metadata)
    metadata = pd.read_csv(DATA_DIR / "diamond_metadata_updated.tsv", sep="\t")

    # Load AMR data in chunks to save memory
    print("Loading AMR/resistance gene data...")
    chunks = []
    chunksize = 100000
    for chunk in pd.read_csv(DATA_DIR / "diamond_amr_combined.tsv", sep="\t", chunksize=chunksize):
        chunks.append(chunk)

    amr_data = pd.concat(chunks, ignore_index=True)
    del chunks  # Free memory

    print(f"Metadata: {len(metadata)} samples")
    print(f"Resistance gene data: {len(amr_data)} entries")
    print(f"  Unique samples: {amr_data['sample_name'].nunique()}")
    print(f"  Unique genes: {amr_data['gene_name'].nunique()}")
    print(f"  Genes with drug_class annotation: {amr_data['drug_class'].notna().sum()} ({amr_data['drug_class'].notna().sum()/len(amr_data)*100:.1f}%)")

    return metadata, amr_data

def create_amr_matrix(metadata, amr_data):
    """Create sample × gene RPM matrix"""
    print("\nCreating AMR gene matrix...")

    # Create matrix: samples (rows) × genes (columns) with RPM values
    amr_matrix = amr_data.pivot_table(
        index='sample_name',
        columns='gene_name',
        values='rpm',
        fill_value=0
    )

    print(f"AMR matrix shape: {amr_matrix.shape}")
    print(f"  Samples: {amr_matrix.shape[0]}")
    print(f"  Genes: {amr_matrix.shape[1]}")

    # Save matrix
    amr_matrix.to_csv(RESULTS_DIR / "amr_rpm_matrix.tsv", sep="\t")
    print(f"✓ Saved: {RESULTS_DIR / 'amr_rpm_matrix.tsv'}")

    return amr_matrix

def calculate_diversity(amr_matrix, metadata):
    """Calculate diversity metrics"""
    print("\nCalculating diversity metrics...")

    # Shannon diversity
    def shannon_diversity(row):
        # Remove zeros
        abundances = row[row > 0]
        if len(abundances) == 0:
            return 0
        # Normalize to proportions
        props = abundances / abundances.sum()
        return entropy(props, base=np.e)

    # Simpson diversity (1 - D)
    def simpson_diversity(row):
        abundances = row[row > 0]
        if len(abundances) == 0:
            return 0
        props = abundances / abundances.sum()
        return 1 - np.sum(props ** 2)

    # Calculate for each sample
    diversity = pd.DataFrame({
        'richness': (amr_matrix > 0).sum(axis=1),
        'shannon': amr_matrix.apply(shannon_diversity, axis=1),
        'simpson': amr_matrix.apply(simpson_diversity, axis=1),
        'total_rpm': amr_matrix.sum(axis=1)
    }, index=amr_matrix.index)

    diversity = diversity.reset_index().rename(columns={'index': 'sample_name'})

    # Merge with metadata
    diversity = diversity.merge(
        metadata[['sample_name', 'SubjectID', 'Location', 'SampleType', 'SampleCollectionWeek']],
        on='sample_name',
        how='left'
    )

    # Save
    diversity.to_csv(RESULTS_DIR / "diversity_metrics.tsv", sep="\t", index=False)
    print(f"✓ Saved: {RESULTS_DIR / 'diversity_metrics.tsv'}")

    # Summary statistics
    print("\nDiversity by Location:")
    print(diversity.groupby('Location')[['richness', 'shannon', 'simpson']].agg(['mean', 'std']))

    print("\nDiversity by Body Site:")
    print(diversity.groupby('SampleType')[['richness', 'shannon', 'simpson']].agg(['mean', 'std']))

    return diversity

def run_pca_analysis(amr_matrix, metadata):
    """Run PCA on AMR matrix, stratified by body site"""
    print("\nRunning PCA analysis...")

    # Filter to samples that have metadata
    samples_with_metadata = metadata[metadata['sample_name'].isin(amr_matrix.index)].copy()
    samples_with_metadata = samples_with_metadata.set_index('sample_name')

    # Get common samples
    common_samples = amr_matrix.index.intersection(samples_with_metadata.index)
    print(f"  Samples in matrix: {len(amr_matrix)}")
    print(f"  Samples with metadata: {len(samples_with_metadata)}")
    print(f"  Common samples: {len(common_samples)}")

    # Filter both to common samples
    amr_matrix = amr_matrix.loc[common_samples]
    samples_with_metadata = samples_with_metadata.loc[common_samples]

    # Log-transform (log10(RPM + 1))
    amr_log = np.log10(amr_matrix + 1)

    # Standardize
    scaler = StandardScaler()
    amr_scaled = scaler.fit_transform(amr_log)

    # PCA
    pca = PCA()
    pca_coords = pca.fit_transform(amr_scaled)

    # Create results dataframe
    pca_df = pd.DataFrame(
        pca_coords[:, :10],  # First 10 PCs
        columns=[f'PC{i+1}' for i in range(10)],
        index=amr_matrix.index
    )

    # Add metadata
    pca_df = pca_df.join(samples_with_metadata[['SubjectID', 'Location', 'SampleType', 'SampleCollectionWeek']])

    # Save PCA results
    pca_df.to_csv(RESULTS_DIR / "pca_coordinates.tsv", sep="\t")

    # Save variance explained
    var_exp = pd.DataFrame({
        'PC': [f'PC{i+1}' for i in range(len(pca.explained_variance_ratio_))],
        'variance_explained': pca.explained_variance_ratio_,
        'cumulative_variance': np.cumsum(pca.explained_variance_ratio_)
    })
    var_exp.to_csv(RESULTS_DIR / "pca_variance_explained.tsv", sep="\t", index=False)

    print(f"  PC1: {pca.explained_variance_ratio_[0]*100:.1f}% variance")
    print(f"  PC2: {pca.explained_variance_ratio_[1]*100:.1f}% variance")
    print(f"  PC3: {pca.explained_variance_ratio_[2]*100:.1f}% variance")
    print(f"✓ Saved PCA results")

    return pca_df, pca, var_exp

def plot_pca(pca_df, var_exp):
    """Create PCA visualization"""
    print("\nCreating PCA plots...")

    # Variance explained by first 10 PCs
    pc1_var = var_exp.loc[0, 'variance_explained'] * 100
    pc2_var = var_exp.loc[1, 'variance_explained'] * 100
    pc3_var = var_exp.loc[2, 'variance_explained'] * 100

    # Figure 1: PC1 vs PC2, separate panels by body site
    fig, axes = plt.subplots(1, 3, figsize=(24, 7))

    sites = ['Axilla', 'Groin', 'Stool']
    colors = {'UCMC': '#1f77b4', 'ZCH': '#ff7f0e'}

    for i, site in enumerate(sites):
        ax = axes[i]
        site_data = pca_df[pca_df['SampleType'] == site]

        for loc in ['UCMC', 'ZCH']:
            loc_data = site_data[site_data['Location'] == loc]
            if len(loc_data) > 0:
                ax.scatter(
                    loc_data['PC1'], loc_data['PC2'],
                    c=colors[loc], label=loc, s=80, alpha=0.7, edgecolors='black', linewidths=0.5
                )

        ax.set_xlabel(f'PC1 ({pc1_var:.1f}%)')
        ax.set_ylabel(f'PC2 ({pc2_var:.1f}%)')
        ax.set_title(f'{site} (n={len(site_data)})')
        ax.legend()
        ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "pca_by_bodysite.pdf", dpi=300, bbox_inches='tight')
    plt.savefig(FIGURES_DIR / "pca_by_bodysite.png", dpi=300, bbox_inches='tight')
    print(f"✓ Saved: {FIGURES_DIR / 'pca_by_bodysite.pdf'}")
    plt.close()

    # Figure 2: PC1 vs PC2, colored by location
    fig, ax = plt.subplots(1, 1, figsize=(10, 8))

    for loc in ['UCMC', 'ZCH']:
        loc_data = pca_df[pca_df['Location'] == loc]
        ax.scatter(
            loc_data['PC1'], loc_data['PC2'],
            c=colors[loc], label=loc, s=80, alpha=0.7, edgecolors='black', linewidths=0.5
        )

    ax.set_xlabel(f'PC1 ({pc1_var:.1f}%)')
    ax.set_ylabel(f'PC2 ({pc2_var:.1f}%)')
    ax.set_title('PCA by Location (UCMC vs ZCH)')
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "pca_by_location.pdf", dpi=300, bbox_inches='tight')
    plt.savefig(FIGURES_DIR / "pca_by_location.png", dpi=300, bbox_inches='tight')
    print(f"✓ Saved: {FIGURES_DIR / 'pca_by_location.pdf'}")
    plt.close()

    # Figure 3: Variance explained scree plot
    fig, ax = plt.subplots(1, 1, figsize=(10, 6))

    pcs = range(1, min(11, len(var_exp)+1))
    ax.bar(pcs, var_exp.iloc[:10]['variance_explained'] * 100, color='steelblue', edgecolor='black')
    ax.set_xlabel('Principal Component')
    ax.set_ylabel('Variance Explained (%)')
    ax.set_title('PCA Scree Plot')
    ax.set_xticks(pcs)
    ax.grid(True, alpha=0.3, axis='y')

    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "pca_scree_plot.pdf", dpi=300, bbox_inches='tight')
    plt.savefig(FIGURES_DIR / "pca_scree_plot.png", dpi=300, bbox_inches='tight')
    print(f"✓ Saved: {FIGURES_DIR / 'pca_scree_plot.pdf'}")
    plt.close()

def plot_diversity(diversity):
    """Plot diversity metrics"""
    print("\nCreating diversity plots...")

    # Filter out samples without body site information
    diversity_clean = diversity[diversity['SampleType'].isin(['Axilla', 'Groin', 'Stool'])].copy()

    fig, axes = plt.subplots(2, 3, figsize=(20, 12))

    # Row 1: By Location
    # Shannon
    sns.boxplot(data=diversity_clean, x='Location', y='shannon', hue='SampleType', ax=axes[0, 0])
    axes[0, 0].set_title('Shannon Diversity by Location & Body Site')
    axes[0, 0].set_ylabel('Shannon Index')

    # Simpson
    sns.boxplot(data=diversity_clean, x='Location', y='simpson', hue='SampleType', ax=axes[0, 1])
    axes[0, 1].set_title('Simpson Diversity by Location & Body Site')
    axes[0, 1].set_ylabel('Simpson Index')

    # Richness
    sns.boxplot(data=diversity_clean, x='Location', y='richness', hue='SampleType', ax=axes[0, 2])
    axes[0, 2].set_title('Gene Richness by Location & Body Site')
    axes[0, 2].set_ylabel('Number of Genes')

    # Row 2: By Body Site
    # Shannon
    sns.boxplot(data=diversity_clean, x='SampleType', y='shannon', hue='Location', ax=axes[1, 0])
    axes[1, 0].set_title('Shannon Diversity by Body Site & Location')
    axes[1, 0].set_ylabel('Shannon Index')
    axes[1, 0].set_xlabel('Body Site')

    # Simpson
    sns.boxplot(data=diversity_clean, x='SampleType', y='simpson', hue='Location', ax=axes[1, 1])
    axes[1, 1].set_title('Simpson Diversity by Body Site & Location')
    axes[1, 1].set_ylabel('Simpson Index')
    axes[1, 1].set_xlabel('Body Site')

    # Richness
    sns.boxplot(data=diversity_clean, x='SampleType', y='richness', hue='Location', ax=axes[1, 2])
    axes[1, 2].set_title('Gene Richness by Body Site & Location')
    axes[1, 2].set_ylabel('Number of Genes')
    axes[1, 2].set_xlabel('Body Site')

    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "diversity_metrics.pdf", dpi=300, bbox_inches='tight')
    plt.savefig(FIGURES_DIR / "diversity_metrics.png", dpi=300, bbox_inches='tight')
    print(f"✓ Saved: {FIGURES_DIR / 'diversity_metrics.pdf'}")
    plt.close()

def run_clustering(amr_matrix, metadata):
    """Hierarchical clustering"""
    print("\nRunning hierarchical clustering...")

    # Filter to samples that have metadata
    samples_with_metadata = metadata[metadata['sample_name'].isin(amr_matrix.index)].copy()
    samples_with_metadata = samples_with_metadata.set_index('sample_name')

    # Get common samples
    common_samples = amr_matrix.index.intersection(samples_with_metadata.index)

    # Filter both to common samples
    amr_matrix = amr_matrix.loc[common_samples]
    samples_with_metadata = samples_with_metadata.loc[common_samples]

    # Log-transform
    amr_log = np.log10(amr_matrix + 1)

    # Compute distance matrix (Euclidean on log-transformed data)
    distances = pdist(amr_log, metric='euclidean')

    # Hierarchical clustering
    linkage_matrix = linkage(distances, method='ward')

    # Create dendrogram with colored labels by location
    fig, ax = plt.subplots(1, 1, figsize=(20, 10))

    # Color mapping
    location_colors = {'UCMC': '#1f77b4', 'ZCH': '#ff7f0e'}

    # No labels for large dendrograms (too many samples)
    dendrogram(
        linkage_matrix,
        no_labels=True,
        ax=ax,
        color_threshold=0
    )

    ax.set_title('Hierarchical Clustering of Samples by AMR Gene Profile (Ward Linkage)')
    ax.set_xlabel('Sample')
    ax.set_ylabel('Distance')

    # Legend
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor=location_colors[loc], label=loc) for loc in ['UCMC', 'ZCH']]
    ax.legend(handles=legend_elements, loc='upper right')

    plt.tight_layout()
    plt.savefig(FIGURES_DIR / "hierarchical_clustering.pdf", dpi=300, bbox_inches='tight')
    plt.savefig(FIGURES_DIR / "hierarchical_clustering.png", dpi=300, bbox_inches='tight')
    print(f"✓ Saved: {FIGURES_DIR / 'hierarchical_clustering.pdf'}")
    plt.close()

def main():
    print("="*60)
    print("DIAMOND RESISTOME ANALYSIS")
    print("="*60)

    # Load data
    metadata, amr_data = load_data()

    # Create AMR matrix
    amr_matrix = create_amr_matrix(metadata, amr_data)

    # Calculate diversity
    diversity = calculate_diversity(amr_matrix, metadata)

    # PCA
    pca_df, pca, var_exp = run_pca_analysis(amr_matrix, metadata)

    # Plot PCA
    plot_pca(pca_df, var_exp)

    # Plot diversity
    plot_diversity(diversity)

    # Clustering
    run_clustering(amr_matrix, metadata)

    print("\n" + "="*60)
    print("DIAMOND RESISTOME ANALYSIS COMPLETE")
    print("="*60)
    print("\nOutputs:")
    print(f"  - AMR RPM matrix: {RESULTS_DIR / 'amr_rpm_matrix.tsv'}")
    print(f"  - Diversity metrics: {RESULTS_DIR / 'diversity_metrics.tsv'}")
    print(f"  - PCA coordinates: {RESULTS_DIR / 'pca_coordinates.tsv'}")
    print(f"  - PCA variance: {RESULTS_DIR / 'pca_variance_explained.tsv'}")
    print(f"\nFigures:")
    print(f"  - PCA by body site: {FIGURES_DIR / 'pca_by_bodysite.pdf'}")
    print(f"  - PCA by location: {FIGURES_DIR / 'pca_by_location.pdf'}")
    print(f"  - PCA scree plot: {FIGURES_DIR / 'pca_scree_plot.pdf'}")
    print(f"  - Diversity metrics: {FIGURES_DIR / 'diversity_metrics.pdf'}")
    print(f"  - Hierarchical clustering: {FIGURES_DIR / 'hierarchical_clustering.pdf'}")

if __name__ == "__main__":
    main()
