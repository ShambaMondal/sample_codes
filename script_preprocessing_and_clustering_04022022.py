'''
This sample script was written for anlayzing scRNA-seq data (24hpf during zebrafish heart development)
'''

import numpy as np
import pandas as pd
import scanpy as sc
from matplotlib import pyplot as plt
import seaborn as sns

sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.autosave= True
sc.settings.figdir='./'
#sc.settings.set_figure_params(dpi=80, facecolor='white') # original
sc.settings.set_figure_params(dpi_save=150, facecolor='white', format='png')

results_file = 'write/zebrafish_nkx2.5_heart_24hpf.h5ad'  # the file that will store the analysis results

# Read in the count matrix into an AnnData object, which holds many slots for annotations and different representations of the data. It also comes with its own HDF5-based file format: .h5ad.
mtx_data_dir = '/mnt/hdd/sam/zdg_research/Eugeniusz_Tralle_single_cell_project_data_30122021/cellranger_pipeline_19012022/run_count_Eugeniusz_Tralle_scRNA_first_sample_19012022_GFP/outs/filtered_feature_bc_matrix/'

adata = sc.read_10x_mtx(
    mtx_data_dir,  # the directory with the `.mtx` file
    var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
    cache=True)                              # write a cache file for faster subsequent reading


adata.var_names_make_unique()  # this is unnecessary if using `var_names='gene_ids'` in `sc.read_10x_mtx`
print('1', adata) # n_obs × n_vars = 9434 × 26573

### Preprocessing

# Show those genes that yield the highest fraction of counts in each single cell, across all cells.
#'''
sc.pl.highest_expr_genes(adata, n_top=20, )
#plt.show()
#'''

# Basic filtering:
sc.pp.filter_cells(adata, min_genes=200) # filtered out 29 cells that have less than 200 genes expressed
sc.pp.filter_genes(adata, min_cells=3) # filtered out 4166 genes that are detected in less than 3 cells

print('2', adata) # n_obs × n_vars = 9405 × 22407

# Let’s assemble some information about mitochondrial genes, which are important for quality control.
#Citing from “Simple Single Cell” workflows (Lun, McCarthy & Marioni, 2017):
#    High proportions are indicative of poor-quality cells (Islam et al. 2014; Ilicic et al. 2016), possibly because of loss of cytoplasmic RNA from perforated cells. The reasoning is that mitochondria are larger than individual transcript molecules and less likely to escape through tears in the cell membrane.
#With pp.calculate_qc_metrics, we can compute many metrics very efficiently.

adata.var['mt'] = adata.var_names.str.startswith('mt-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True) # originally in the preprocessing pipeline

print('3', adata)

# visualize qcmetrics outcome (https://scanpy.readthedocs.io/en/stable/generated/scanpy.pp.calculate_qc_metrics.html#scanpy.pp.calculate_qc_metrics)
#'''
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=True, inplace=True) # modified for seaborn plots
sns.jointplot(
    data=adata.obs,
    x="log1p_total_counts",
    y="log1p_n_genes_by_counts",
    kind="hex",
)
plt.savefig('Figure_2_sns_jointplot_log1p_total_counts_vs_log1p_n_genes_by_counts_06022022.png', dpi=600)
plt.clf()
sns.histplot(adata.obs["pct_counts_mt"])
plt.savefig('Figure_2_sns_histplot_pct_counts_mt_06022022.png', dpi=600)
#'''

# A violin plot of some of the computed quality measures:
#   - the number of genes expressed in the count matrix
#   - the total counts per cell
#   - the percentage of counts in mitochondrial genes
#'''
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True, save= '_Figure_2_quality_measures_num_genes_expressed_in_count_matrix_and_tot_counts_per_cell_and_pct_of_counts_in_mito_genes_05022022.png')
#'''

# Remove cells that have too many mitochondrial genes expressed or too many total counts:
#'''
sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt', save='_Figure_3a_total_counts_vs_pct_counts_mt_05022022.png')

sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', save= '_Figure_3b_total_counts_vs_n_genes_by_counts_05022022.png')
#'''

# Actually do the filtering by slicing the AnnData object.
print('4', adata[adata.obs.pct_counts_mt < 5, :]) # n_obs × n_vars = 7911 × 22407
print(adata[adata.obs.pct_counts_mt < 10, :]) # n_obs × n_vars = 8854 × 22407
print(adata[adata.obs.pct_counts_mt < 20, :]) # n_obs × n_vars = 9299 × 22407
adata = adata[adata.obs.n_genes_by_counts < 6000, :]
adata = adata[adata.obs.pct_counts_mt < 10, :]

# Total-count normalize (library-size correct) the data matrix X to 10,000 reads per cell, so that counts become comparable among cells.
# We decided to normalize it to 55,000 reads per cell (mean reads per cell in our data)
sc.pp.normalize_total(adata, target_sum=55000)

# Logarithmize the data:
sc.pp.log1p(adata)

# Identify highly-variable genes
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pl.highly_variable_genes(adata, save='_Figure_4_highly_variable_genes_06022022.png')

# Set the .raw attribute of the AnnData object to the normalized and logarithmized raw gene expression for later use in differential testing and visualizations of gene expression. This simply freezes the state of the AnnData object.
#NOTE: You can get back an AnnData of the object in .raw by calling .raw.to_adata().
adata.raw = adata

#If you don’t proceed below with correcting the data with sc.pp.regress_out and scaling it via sc.pp.scale, you can also get away without using .raw at all.
#The result of the previous highly-variable-genes detection is stored as an annotation in .var.highly_variable and auto-detected by PCA and hence, sc.pp.neighbors and subsequent manifold/graph tools. In that case, the step actually do the filtering below is unnecessary, too.

# Actually do the filtering
#adata = adata[:, adata.var.highly_variable]

# Regress out effects of total counts per cell and the percentage of mitochondrial genes expressed. Scale the data to unit variance.
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])

# Scale each gene to unit variance. Clip values exceeding standard deviation 10.
sc.pp.scale(adata, max_value=10)

### Principal component analysis
# Reduce the dimensionality of the data by running principal component analysis (PCA), which reveals the main axes of variation and denoises the data.
sc.tl.pca(adata, svd_solver='arpack')

# We can make a scatter plot in the PCA coordinates, but we will not use that later on.
sc.pl.pca(adata, color='nkx2.5', save= '_Figure_5_scatter_plot_PCA_coordinates_06022022.png')

# Let us inspect the contribution of single PCs to the total variance in the data. This gives us information about how many PCs we should consider in order to compute the neighborhood relations of cells, e.g. used in the clustering function sc.tl.louvain() or tSNE sc.tl.tsne(). In our experience, often a rough estimate of the number of PCs does fine.
sc.pl.pca_variance_ratio(adata, log=True)

# Save the result.
adata.write(results_file)

print(adata)

### Computing the neighborhood graph

# Let us compute the neighborhood graph of cells using the PCA representation of the data matrix. You might simply use default values here. For the sake of reproducing Seurat’s results, let’s take the following values.
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)

### Embedding the neighborhood graph
## We suggest embedding the graph in two dimensions using UMAP (McInnes et al., 2018), see below. It is potentially more faithful to the global connectivity of the manifold than tSNE, i.e., it better preserves trajectories. In some ocassions, you might still observe disconnected clusters and similar connectivity violations. They can usually be remedied by running:
## tl.paga(adata)
## pl.paga(adata, plot=False)  # remove `plot=False` if you want to see the coarse-grained graph
## tl.umap(adata, init_pos='paga')

sc.tl.umap(adata)
sc.pl.umap(adata, save='_raw_True_06022022.png')

# As we set the .raw attribute of adata, the previous plots showed the “raw” (normalized, logarithmized, but uncorrected) gene expression. You can also plot the scaled and corrected gene expression by explicitly stating that you don’t want to use .raw.
sc.pl.umap(adata, use_raw=False, save='_raw_False_06022022.png')

### Clustering the neighborhood graph
# As with Seurat and many other frameworks, we recommend the Leiden graph-clustering method (community detection based on optimizing modularity) by Traag *et al.* (2018). Note that Leiden clustering directly clusters the neighborhood graph of cells, which we already computed in the previous section.
sc.tl.leiden(adata)

# Plot the clusters
##sc.pl.umap(adata, color=['leiden', 'nkx2.5'], save='_leiden_clusters_06022022.png')
sc.pl.umap(adata, color='leiden', save='_leiden_clusters_06022022.png')
sc.pl.umap(adata, color='leiden', legend_loc='on data', save='_leiden_clusters_marked_on_data_06022022.png')
sc.pl.umap(adata, color='nkx2.5', save='_leiden_clusters_colored_by_nkx2.5_06022022.png')
sc.pl.umap(adata, color=['nkx2.5', 'GFP'], save='_leiden_clusters_colored_by_nkx2.5_and_GFP_06022022.png')
sc.pl.umap(adata, color=['nkx2.7', 'GFP'], save='_leiden_clusters_colored_by_nkx2.7_and_GFP_06022022.png')
sc.pl.umap(adata, color=['nkx2.5', 'nkx2.7', 'GFP'], save='_leiden_clusters_colored_by_nkx2.5_nkx2.7_and_GFP_06022022.png')

# Save the result.
adata.write(results_file)

### Finding marker genes
#Let us compute a ranking for the highly differential genes in each cluster. For this, by default, the .raw attribute of AnnData is used in case it has been initialized before. The simplest and fastest method to do so is the t-test.
sc.tl.rank_genes_groups(adata, 'leiden', method='t-test')
#sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False, save='_marker_genes_t_test_06022022.png')
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False, save='_marker_genes_t_test_06022022.pdf')


sc.settings.verbosity = 2  # reduce the verbosity

# The result of a Wilcoxon rank-sum (Mann-Whitney-U) test is very similar. We recommend using the latter in publications, see e.g., Sonison & Robinson (2018). You might also consider much more powerful differential testing packages like MAST, limma, DESeq2 and, for python, the recent diffxpy.

sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
#sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False, save='_marker_genes_Wilcoxon_test_06022022.png')
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False, save='_marker_genes_Wilcoxon_test_06022022.pdf')

# Save the result.
adata.write(results_file)

# As an alternative, let us rank genes using logistic regression. For instance, this has been suggested by Natranos et al. (2018). The essential difference is that here, we use a multi-variate appraoch whereas conventional differential tests are uni-variate. Clark et al. (2014) has more details.
sc.tl.rank_genes_groups(adata, 'leiden', method='logreg')
#sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False, save='_marker_genes_logistic_regression_06022022.png')
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False, save='_marker_genes_logistic_regression_06022022.pdf')

