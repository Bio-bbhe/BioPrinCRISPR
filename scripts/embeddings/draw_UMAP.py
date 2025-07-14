import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset

csv_path = 'umap_embeddings_with_clusters.csv'

umap_df = pd.read_csv(csv_path, sep='\t', index_col=0)
if 'Sequence_ID' not in umap_df.columns:
    umap_df.reset_index(inplace=True)

num_unique_clusters = umap_df['Cluster_ID'].nunique()

def highlight_type(seq_id):
    if seq_id.startswith('Cas13an'):
        return 'Cas13an'
    elif seq_id.startswith('MEK') or seq_id.startswith('MEQ'):
        return 'MEK/MEQ'
    else:
        return 'Other'

umap_df['Highlight_Type'] = umap_df['Sequence_ID'].apply(highlight_type)

palette = sns.color_palette("tab20", n_colors=num_unique_clusters)

fig, ax = plt.subplots(figsize=(12, 10))

sns.scatterplot(
    x='UMAP_X', y='UMAP_Y',
    hue='Cluster_ID',
    palette=palette,
    data=umap_df[umap_df['Highlight_Type']=='Other'],
    legend=False,
    alpha=0.7,
    s=25,
    linewidth=0,
    ax=ax
)

cas13an_points = umap_df[umap_df['Highlight_Type'] == 'Cas13an']
sns.scatterplot(
    x='UMAP_X', y='UMAP_Y',
    hue='Cluster_ID',
    palette=palette,
    data=cas13an_points,
    legend=False,
    edgecolor='black',
    linewidth=2,
    alpha=1,
    s=80,
    ax=ax,
    label='Known Cas13an'
)

mekmeq_color = '#D32F2F' 
mekmeq_points = umap_df[umap_df['Highlight_Type'] == 'MEK/MEQ']
sns.scatterplot(
    x='UMAP_X', y='UMAP_Y',
    hue='Cluster_ID',
    palette=palette,
    data=mekmeq_points,
    legend=False,
    edgecolor=mekmeq_color,
    linewidth=2,
    alpha=1,
    s=80,
    ax=ax,
    label='PbCas13an\nSpCas13an'
)
#ax.set_title(f'UMAP Projection with K-Means Clusters (k={num_unique_clusters})',
#             fontsize=18, weight='bold')
ax.set_xlabel('')
ax.set_ylabel('')
ax.tick_params(axis='both', labelsize=14)

for spine in ax.spines.values():
    spine.set_visible(True)
    spine.set_linewidth(2)
    spine.set_color('black')

ax.grid(True, linestyle='--', alpha=0.5)

ax.legend(loc='lower right', fontsize=16)

highlight_points = umap_df[umap_df['Highlight_Type'] != 'Other']
margin = 0.3

x_min, x_max = highlight_points['UMAP_X'].min() - margin, highlight_points['UMAP_X'].max() + margin
y_min, y_max = highlight_points['UMAP_Y'].min() - margin, highlight_points['UMAP_Y'].max() + margin

min_range = 0.8
if (x_max - x_min) < min_range:
    mid_x = (x_max + x_min) / 2
    x_min, x_max = mid_x - min_range / 2, mid_x + min_range / 2
if (y_max - y_min) < min_range:
    mid_y = (y_max + y_min) / 2
    y_min, y_max = mid_y - min_range / 2, mid_y + min_range / 2

axins = inset_axes(ax, width=2.5, height=2.5, loc='upper right', borderpad=2)

sns.scatterplot(
    x='UMAP_X', y='UMAP_Y',
    hue='Cluster_ID',
    palette=palette,
    data=umap_df[umap_df['Highlight_Type']=='Other'],
    legend=False,
    alpha=0.7,
    s=40,
    linewidth=0,
    ax=axins
)

sns.scatterplot(
    x='UMAP_X', y='UMAP_Y',
    hue='Cluster_ID',
    palette=palette,
    data=cas13an_points,
    legend=False,
    edgecolor='black',
    linewidth=2,
    alpha=1,
    s=100,
    ax=axins
)

sns.scatterplot(
    x='UMAP_X', y='UMAP_Y',
    hue='Cluster_ID',
    palette=palette,
    data=mekmeq_points,
    legend=False,
    edgecolor=mekmeq_color,
    linewidth=2,
    alpha=1,
    s=100,
    ax=axins
)

axins.set_xlim(x_min, x_max)
axins.set_ylim(y_min, y_max)
axins.set_xticks([])
axins.set_yticks([])

ax.set_xlabel('UMAP_X', fontsize=16)
ax.set_ylabel('UMAP_Y', fontsize=16)


for spine in axins.spines.values():
    spine.set_linewidth(1.5)
    spine.set_color('black')


mark_inset(ax, axins, loc1=2, loc2=4, fc="none", ec="0.5", lw=2)

plt.show()
