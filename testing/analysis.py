import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import StandardScaler
from scipy import stats
import os

# Create output directory for results
os.makedirs('analysis_results', exist_ok=True)

# Set random seed for reproducibility
np.random.seed(42)

def preprocess_data(adata):
    """Preprocess the data for trajectory analysis"""
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    sc.pp.pca(adata)
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    sc.tl.diffmap(adata)
    sc.tl.dpt(adata, n_dcs=10)
    return adata

def analyze_cell_type_bias(adata):
    """Analyze potential biases in cell type distributions and trajectories"""
    # Calculate distribution statistics with renamed column
    cell_type_stats = pd.DataFrame({
        'cell_count': adata.obs['paul15_clusters'].value_counts(),
        'percentage': adata.obs['paul15_clusters'].value_counts(normalize=True) * 100
    })
    
    # Calculate pseudotime statistics with explicit column names
    pseudotime_stats = adata.obs.groupby('paul15_clusters', observed=True)['dpt_pseudotime'].agg([
        ('pseudotime_mean', 'mean'),
        ('pseudotime_std', 'std'),
        ('pseudotime_count', 'count')
    ])
    
    # Merge the dataframes
    combined_stats = pd.concat([cell_type_stats, pseudotime_stats], axis=1)
    
    # Save statistics to CSV
    combined_stats.to_csv('analysis_results/cell_type_statistics.csv')
    return combined_stats

def save_visualizations(adata, stats):
    """Create and save visualizations"""
    # 1. Cell type distribution
    plt.figure(figsize=(10, 6))
    stats['percentage'].plot(kind='bar')
    plt.title('Cell Type Distribution')
    plt.xticks(rotation=45)
    plt.ylabel('Percentage')
    plt.savefig('analysis_results/cell_type_distribution.png')
    plt.close()
    
    # 2. UMAP by cell type
    plt.figure(figsize=(30, 15))
    sc.pl.umap(adata, color='paul15_clusters', show=False)
    plt.title('UMAP by Cell Type')
    plt.subplots_adjust(right=0.62)  # Increase space for legend
    plt.savefig('analysis_results/umap_cell_types.png')
    plt.close()    

    # 3. Pseudotime distribution
    plt.figure(figsize=(12, 6))
    sns.boxplot(data=adata.obs, x='paul15_clusters', y='dpt_pseudotime')
    plt.xticks(rotation=45)
    plt.title('Pseudotime Distribution by Cell Type')
    plt.savefig('analysis_results/pseudotime_distribution.png')
    plt.close()
    
    # 4. Gene expression patterns
    plt.figure(figsize=(15, 6))
    key_genes = ['Mpo', 'Elane', 'Prtn3']
    sc.pl.violin(adata, key_genes, groupby='paul15_clusters', show=False)
    plt.title('Key Gene Expression by Cell Type')
    plt.savefig('analysis_results/gene_expression_patterns.png')
    plt.close()
    
    # 5. Trajectory density
    plt.figure(figsize=(10, 6))
    sns.kdeplot(data=adata.obs, x='dpt_pseudotime', hue='paul15_clusters')
    plt.title('Trajectory Density by Cell Type')
    plt.savefig('analysis_results/trajectory_density.png')
    plt.close()

    #6. PAGA Graph
    plt.figure(figsize=(20, 15))
    sc.tl.paga(adata, groups='paul15_clusters')
    sc.pl.paga(adata, color='paul15_clusters', node_size_scale=2, show=False)
    plt.title('PAGA Graph')
    plt.savefig('analysis_results/paga_graph.png')
    plt.close()

    #7. PAGA Graph with pseudotime
    plt.figure(figsize=(20, 15))
    sc.pl.paga(adata, color='dpt_pseudotime', node_size_scale=2, show=False)
    plt.title('PAGA Graph with Pseudotime')
    plt.savefig('analysis_results/paga_graph_pseudotime.png')
    plt.close()
    
def calculate_fairness_metrics(adata):
    """Calculate fairness metrics across cell types"""
    fairness_metrics = {}
    total_cells = len(adata)
    expected_proportion = 1.0 / len(adata.obs['paul15_clusters'].unique())
    
    for cell_type in adata.obs['paul15_clusters'].unique():
        cell_count = sum(adata.obs['paul15_clusters'] == cell_type)
        actual_proportion = cell_count / total_cells
        disparity_ratio = actual_proportion / expected_proportion
        
        cell_type_data = adata[adata.obs['paul15_clusters'] == cell_type]
        mean_pseudotime = np.mean(cell_type_data.obs['dpt_pseudotime'])
        std_pseudotime = np.std(cell_type_data.obs['dpt_pseudotime'])
        
        fairness_metrics[cell_type] = {
            'representation_ratio': disparity_ratio,
            'mean_pseudotime': mean_pseudotime,
            'std_pseudotime': std_pseudotime,
            'sample_size': cell_count
        }
    
    # Save fairness metrics to CSV
    metrics_df = pd.DataFrame(fairness_metrics).T
    metrics_df.to_csv('analysis_results/fairness_metrics.csv')
    return fairness_metrics

def main():
    print("Loading Paul15 dataset...")
    adata = sc.datasets.paul15()
    
    print("Preprocessing data...")
    adata = preprocess_data(adata)
    
    print("Running cell type bias analysis...")
    cell_type_stats = analyze_cell_type_bias(adata)
    
    print("Calculating fairness metrics...")
    fairness_metrics = calculate_fairness_metrics(adata)
    
    print("Generating and saving visualizations...")
    save_visualizations(adata, cell_type_stats)
    
    # Save summary statistics to a text file
    with open('analysis_results/analysis_summary.txt', 'w') as f:
        f.write("Cell Type Statistics:\n")
        f.write(str(cell_type_stats))
        f.write("\n\nFairness Metrics:\n")
        for cell_type, metrics in fairness_metrics.items():
            f.write(f"\n{cell_type}:\n")
            for metric, value in metrics.items():
                f.write(f"  {metric}: {value:.3f}\n")
        
        # Add ANOVA test results
        f_stat, p_val = stats.f_oneway(*[
            adata[adata.obs['paul15_clusters'] == ct].obs['dpt_pseudotime']
            for ct in adata.obs['paul15_clusters'].unique()
        ])
        f.write(f"\nANOVA test for pseudotime differences:\n")
        f.write(f"F-statistic: {f_stat:.3f}\n")
        f.write(f"p-value: {p_val:.3e}\n")
    
    print("\nAnalysis complete! Results have been saved to the 'analysis_results' directory.")
    print("Files generated:")
    print("- cell_type_statistics.csv")
    print("- fairness_metrics.csv")
    print("- analysis_summary.txt")
    print("- cell_type_distribution.png")
    print("- umap_cell_types.png")
    print("- pseudotime_distribution.png")
    print("- gene_expression_patterns.png")
    print("- trajectory_density.png")
    print("- paga_graph.png")
    print("- paga_graph_pseudotime.png")

if __name__ == "__main__":
    main()
