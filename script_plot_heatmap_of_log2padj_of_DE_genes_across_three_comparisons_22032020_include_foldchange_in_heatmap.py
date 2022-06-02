import pandas as pd
import seaborn as sns
import numpy as np
from matplotlib import pyplot as plt
#from matplotlib import cm
from matplotlib import patches as mpatches
from matplotlib import gridspec as mgridspec


df = pd.read_table('dataframe_of_DE_genes_across_3_comparisons_with_log2foldchange_labels_rowcolors.csv')
print(df)

#'''
df = pd.read_table('final_dataframe_of_DE_genes_across_3_comparisons.csv')
print(df)
# shorten names of columns
df = df.rename(columns={'GCM6h_vs_UNTR': '6h', 'GCM24h_vs_UNTR': '24h', 'GCM48h_vs_UNTR': '48h'})
print(df)

### label rows with significant comparison
#df.loc[:, 'label'] = ['']*len(df)

#df[df.GCM6h_vs_UNTR < 0.05)]
#print(df < 0.05)
cond = df < 0.05
print(cond)
label = []
for val in cond.values:#[:5]:
    label.append(', '.join(list(cond.columns[val].values)))

df.loc[:, 'label'] = label
print(df)
print(df.label.unique())

#lut = dict(zip(df.label.unique(), "rbgcmyk"))
#lut = dict(zip(df.label.unique(), cm.Paired))
#print(lut)
#row_colors = df.label.map(lut)
#print(row_colors)

#lut = {'GCM24h_vs_UNTR': 'g', 'GCM6h_vs_UNTR': 'r', 'GCM24h_vs_UNTR, GCM48h_vs_UNTR': 'b', 'GCM48h_vs_UNTR': 'm', 'GCM6h_vs_UNTR, GCM24h_vs_UNTR, GCM48h_vs_UNTR': 'k', 'GCM6h_vs_UNTR, GCM48h_vs_UNTR': 'y', 'GCM6h_vs_UNTR, GCM24h_vs_UNTR': 'c'}
lut = {'24h': 'g', '6h': 'r', '24h, 48h': 'b', '48h': 'm', '6h, 24h, 48h': 'k', '6h, 48h': 'y', '6h, 24h': 'c'}
row_colors = df.label.map(lut)
print(row_colors)
for col in ['6h', '24h', '48h']: #['GCM6h_vs_UNTR', 'GCM24h_vs_UNTR', 'GCM48h_vs_UNTR']:
    df.loc[:, col] = -np.log10(df.loc[:, col])
#g = sns.clustermap(iris, row_colors=row_colors)
#print(df)

### add row_colors column to the dataframe with labels
df.loc[:, 'row_colors'] = row_colors
#'''

### the following block was used to create combinatorial clustermaps (method and metric combo)
'''
for method in ['single', 'complete', 'average', 'weighted', 'centroid', 'median', 'ward']:
    for metric in ['euclidean', 'minkowski', 'cityblock', 'seuclidean', 'sqeuclidean', 'cosine', 'correlation', 'chebyshev', 'canberra', 'braycurtis', 'mahalanobis', 'wminkowski']:
        try:
            g = sns.clustermap(df.loc[:, ['GCM6h_vs_UNTR', 'GCM24h_vs_UNTR', 'GCM48h_vs_UNTR']], method=method, metric=metric, figsize=(12, 20), col_cluster=False, cmap='Blues', vmax=20, yticklabels='', row_colors=row_colors)
            plt.title('method=%s,\nmetric=%s'%(method, metric))
            plt.savefig('clustermap_%s_%s.png'%(method, metric))
            plt.close('all')
        except ValueError:
            pass
'''

#'''
#sns.clustermap(df.loc[:, ['GCM6h_vs_UNTR', 'GCM24h_vs_UNTR', 'GCM48h_vs_UNTR']], method='centroid', metric='euclidean', col_cluster=False, cmap='Blues', vmax=20, yticklabels='', row_colors=row_colors)
#plt.legend()
g = sns.clustermap(cond, row_colors=row_colors, yticklabels='')
#plt.show()
print(dir(g))
#['__class__', '__delattr__', '__dict__', '__dir__', '__doc__', '__eq__', '__format__', '__ge__', '__getattribute__', '__gt__', '__hash__', '__init__', '__init_subclass__', '__le__', '__lt__', '__module__', '__ne__', '__new__', '__reduce__', '__reduce_ex__', '__repr__', '__setattr__', '__sizeof__', '__str__', '__subclasshook__', '__weakref__', '_clean_axis', '_get_palette', '_legend_out', '_margin_titles', '_preprocess_colors', '_update_legend_data', 'add_legend', 'ax_col_colors', 'ax_col_dendrogram', 'ax_heatmap', 'ax_row_colors', 'ax_row_dendrogram', 'cax', 'col_color_labels', 'col_colors', 'color_list_to_matrix_and_cmap', 'data', 'data2d', 'dendrogram_col', 'dendrogram_row', 'dim_ratios', 'fig', 'format_data', 'gs', 'mask', 'plot', 'plot_colors', 'plot_dendrograms', 'plot_matrix', 'row_color_labels', 'row_colors', 'savefig', 'set', 'standard_scale', 'z_score']
#print(g.dendrogram_row.reordered_ind)
reordered_index = np.array(g.dendrogram_row.reordered_ind)

### now reset index of df
df = df.reset_index()
#df = df.rename(columns={'index': 'gene'})
print(df)
to_plot_df = df.loc[reordered_index, :]
print(to_plot_df)
to_plot_df = to_plot_df.set_index('index')
to_plot_df.index.name = None
print(to_plot_df)
### save to_plot_df to a file
to_plot_df.to_csv('dataframe_of_DE_genes_across_3_comparisons_with_-log10padj_labels_rowcolors.csv', index=True, sep='\t')

sns.set(font_scale=1.0, palette='muted') # this sets the general font size, also changes the color palette to "deep". See: https://seaborn.pydata.org/generated/seaborn.set.html

### adding legends for row_colors: step 1
lut_reordered = {'6h': 'r', '24h': 'g', '6h, 24h': 'c', '6h, 24h, 48h': 'k','24h, 48h': 'b', '6h, 48h': 'y', '48h': 'm'}
lut_reordered_for_legend = {'6h_only': 'r', '24h_only': 'g', '6h, 24h': 'c', '6h, 24h, 48h': 'k','24h, 48h': 'b', '6h, 48h': 'y', '48h_only': 'm'}
#legend_rowcolors = [mpatches.Patch(color=c, label=l) for c,l in zip(lut_reordered.values(), lut_reordered.keys())] #df[['tissue type','label']].drop_duplicates().values]
legend_rowcolors = [mpatches.Patch(color=c, label=l) for c,l in zip(lut_reordered_for_legend.values(), lut_reordered_for_legend.keys())]

### plot clustermap with reordered index, with row_cluster=False
#g2 = sns.clustermap(to_plot_df.loc[:, ['GCM6h_vs_UNTR', 'GCM24h_vs_UNTR', 'GCM48h_vs_UNTR']], row_cluster=False, col_cluster=False, cmap='Blues', vmax=20, yticklabels='', row_colors=row_colors)
g2 = sns.clustermap(to_plot_df.loc[:, ['6h', '24h', '48h']], row_cluster=False, col_cluster=False, cmap='Blues', vmax=20, yticklabels='', row_colors=row_colors)

### adding -log10pval near the colorbar
g2.cax.set_title('-log10(adj. pvalue)', fontsize=14)


### adding legends for row_colors: step 2
#l2=g2.ax_heatmap.legend(loc='center left',bbox_to_anchor=(1.01,0.85),handles=legend_rowcolors,frameon=True)
l2=g2.ax_heatmap.legend(loc='upper center', bbox_to_anchor=(0.46, 1.34), handles=legend_rowcolors,frameon=False)
l2.set_title(title='Subclasses of DE genes', prop={'size':14})

g2.ax_row_colors.set_xticklabels(['subclass'])

### adding -log10pval near the colorbar
#plt.title('-log10(adjusted pvalue)')
#g2.cax.set_title('-log10(adjusted pvalue)', fontsize=14)

#g2.ax_row_colors.set_xlabel('DE group')
#g2.ax_row_colors.set_label('DE group')
#g2.ax_row_colors.set_xticklabels(['DE class'])
#g2.ax_row_colors.set_ylabel('test_row_color')
#g2.ax_heatmap.annotate('annotate_axh', (2.0, 2.0), fontsize=20)
#g2.cax.annotate('annotate_cax', (0.2, -3.0), fontsize=20)
#plt.annotate('plt_annot', (5.0, 5.0), fontsize=20) # this defaults to pvalue color_bar axis
#print(dir(g2.ax_row_colors))
#print(dir(g2.ax_row_dendogram))
#g2.ax_heatmap.scatter([1.0, 2.0], [2.0, 4.0])
#g2.ax_heatmap.set_label('test')
#g2.ax_heatmap.set_ylabel('test')
#print(dir(g2.col_color_labels))
#print(dir(g2.ax_heatmap))
#print(dir(g2.ax_heatmap.set_xticklabels))
#plt.title('test, test')
#plt.show()


# set the gridspec to only cover half of the figure
g2.gs.update(left=0.05, right=0.45)

sns.set(style="whitegrid", palette='muted')
#create new gridspec for the right part
#gs2 = mgridspec.GridSpec(1,1, left=0.55)
gs2 = mgridspec.GridSpec(2,1, left=0.55, height_ratios=[0.25, 0.75])
# create axes within this new gridspec
ax2 = g2.fig.add_subplot(gs2[0])
ax3 = g2.fig.add_subplot(gs2[1])
# plot boxplot in the new axes
#sns.boxplot(data=to_plot_df.loc[:, ['GCM6h_vs_UNTR', 'GCM24h_vs_UNTR', 'GCM48h_vs_UNTR']], orient="h", palette="Set2", ax = ax2)
#sns.boxplot(data=to_plot_df.loc[:, ['6h', '24h', '48h']], orient="h", palette="Set2", ax = ax2)

### violinplot of significant genes in comparisons: 6h, 24h, 48h (vs UNTR)
fc_df = pd.read_table('dataframe_of_DE_genes_across_3_comparisons_with_log2foldchange_labels_rowcolors.csv')
#sns.boxplot(data=fc_df.loc[:, ['6h', '24h', '48h']], orient="h", palette=['r', 'g', 'm'], ax = ax2)
#sns.violinplot(data=fc_df.loc[:, ['6h', '24h', '48h']], orient="h", palette=['r', 'g', 'm'], ax = ax2)

### make a longform dataframe and plot violinplot for 7 subgroups
fc_df = pd.read_table('dataframe_of_DE_genes_across_3_comparisons_with_log2foldchange_labels_rowcolors.csv')
print(fc_df)
#sns.boxplot(data=fc_df.loc[:, ['6h', '24h', '48h']], orient="h", palette=['r', 'g', 'm'], ax = ax2)
#sns.violinplot(data=fc_df.loc[:, ['6h', '24h', '48h']], orient="h", palette=['r', 'g', 'm'], ax = ax2)

### Add a kdeplot of log2foldchange of 6h/24h/48h
##fc_df = pd.read_table('dataframe_of_DE_genes_across_3_comparisons_with_log2foldchange_labels_rowcolors.csv')
##sns.set(style='whitegrid', palette='muted')
#sns.kdeplot(fc_df.loc[:, '6h'].dropna(), shade=True, color='r', ax=ax2)#, palette={'6h': 'r', '24h': 'g', '48h': 'm'})
#sns.kdeplot(fc_df.loc[:, '24h'].dropna(), shade=True, color='g', ax=ax2)
#sns.kdeplot(fc_df.loc[:, '48h'].dropna(), shade=True, color='m', ax=ax2)
sns.kdeplot(fc_df.loc[:, '6h'].dropna(), color='k', linestyle='solid', ax=ax2)
sns.kdeplot(fc_df.loc[:, '24h'].dropna(), color='k', linestyle='dotted', ax=ax2)
sns.kdeplot(fc_df.loc[:, '48h'].dropna(), color='k', linestyle='dashed', ax=ax2)
ax2.set_xlabel('log2FoldChange')
ax2.set_ylabel('density')
ax2.legend(title='All DE Genes at', loc='upper right', frameon=False, bbox_to_anchor=(1.04, 1.0))#, prop={'size':14})

### make a longform dataframe and plot violinplot for 7 subgroups
fc_df = fc_df.reset_index() # this transforms the older index into a column named 'index'
fc_df = fc_df.melt(id_vars=['label', 'row_colors'], value_vars=['6h', '24h', '48h'], var_name='time', value_name='log2FoldChange')
#fc_df = fc_df.rename(columns={'value': 'log2FoldChange'}) # the default column name for ex-columns is 'variable'(6h, 24h, 48h are values now), and the values (of 6h,24h,48h) is 'value'
print(fc_df)
#lut_reordered = {'6h': 'r', '24h': 'g', '6h, 24h': 'c', '6h, 24h, 48h': 'k','24h, 48h': 'b', '6h, 48h': 'y', '48h': 'm'} #already defined above
#sns.violinplot(data=fc_df, x='variable', y='value', hue='label', palette=lut_reordered) #orient='h')
violin = sns.violinplot(data=fc_df, y='time', x='log2FoldChange', hue='label', palette=lut_reordered, orient='h', ax= ax3,)
violin.legend_.remove()


g2.savefig('z_test.png')
#g2.savefig('z_test_hd.png', dpi=1000)
#'''
