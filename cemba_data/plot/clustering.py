import numpy as np
import pandas as pd
import seaborn as sns

from .cell import categorical_scatter


def _plot_kmode_overall_single_r(ax, overall_df, pureness_ylim=None, completeness_ylim=None):
    overall_df = overall_df.copy()
    overall_df['k'] = overall_df['k'].astype(int)

    sns.lineplot(data=overall_df, color='orange',
                 x='k', y='pureness', ax=ax, label='Pureness',
                 markers=True, dashes=True, legend=None)
    sns.scatterplot(data=overall_df, color='orange',
                    x='k', y='pureness', ax=ax, legend=None)
    ax_twin = ax.twinx()
    sns.lineplot(data=overall_df, color='steelblue',
                 x='k', y='completeness', ax=ax_twin, label='Completeness',
                 markers=True, dashes=True, legend=None, )
    sns.scatterplot(data=overall_df, color='steelblue',
                    x='k', y='completeness', ax=ax_twin, legend=None)
    ax.grid()
    if pureness_ylim is None:
        pureness_ylim = tuple(np.quantile(overall_df['pureness'], q=(0.01, 0.99)))
    ax.set(ylim=pureness_ylim)

    if completeness_ylim is None:
        completeness_ylim = tuple(np.quantile(overall_df['completeness'], q=(0.01, 0.99)))
    ax_twin.set(ylim=completeness_ylim)
    return ax


def plot_kmode_overall(axes, overall_df, pureness_ylim=None, completeness_ylim=None):
    for ax, (resolution, sub_df) in zip(axes.flat, overall_df.groupby('resolution')):
        _plot_kmode_overall_single_r(ax, sub_df, pureness_ylim=pureness_ylim, completeness_ylim=completeness_ylim)
        ax.set_title(f'Resolution = {resolution}')
    return axes


def plot_kmode_stats(axes, coord_data, cell_profile, coord_base='umap'):
    plot_data = coord_data.copy()
    plot_data['cell_ambiguity'] = cell_profile['cell_ambiguity']
    plot_data['cluster'] = cell_profile['cluster']
    plot_data['cluster_pureness'] = cell_profile['cluster_pureness']
    plot_data['cluster_completeness'] = cell_profile['cluster_completeness']

    if axes.size != 4:
        raise ValueError('Number of axes is not 4.')
    ax1, ax2, ax3, ax4 = axes.flatten()
    categorical_scatter(data=plot_data, hue='cluster', ax=ax1,
                        coord_base=coord_base, palette='tab20')
    if ax1.is_first_row():
        ax1.set_title('Cluster')
    sns.scatterplot(data=plot_data, x=f'{coord_base}_0', y=f'{coord_base}_1',
                    hue=1 - plot_data['cell_ambiguity'], ax=ax2, hue_norm=(0.9, 1),
                    palette='viridis', s=3, linewidth=0, legend=None)
    if ax2.is_first_row():
        ax2.set_title('1 - Cell Ambiguity')
    sns.scatterplot(data=plot_data, x=f'{coord_base}_0', y=f'{coord_base}_1',
                    hue='cluster_pureness', ax=ax3, hue_norm=(0, 1),
                    palette='viridis', s=3, linewidth=0, legend=None)
    if ax3.is_first_row():
        ax3.set_title('Cluster Pureness')
    sns.scatterplot(data=plot_data, x=f'{coord_base}_0', y=f'{coord_base}_1',
                    hue='cluster_completeness', ax=ax4, hue_norm=(0, 1),
                    palette='viridis', s=3, linewidth=0, legend=None)
    if ax4.is_first_row():
        ax4.set_title('Cluster Completeness')
    for ax in axes.flat:
        ax.axis('off')
    return axes


def auto_order(confusion_matrix, row_order):
    confusion_matrix = confusion_matrix.reindex(row_order).copy()
    records = []
    for col, row in confusion_matrix.idxmax().iteritems():
        value = confusion_matrix.loc[row, col]
        records.append([col, row, value])
    df = pd.DataFrame(records)
    order_dict = {}
    for row, sub_df in df.sort_values(2).groupby(1):
        order_dict[row] = sub_df[0].tolist()
    col_order = []
    for row in row_order:
        try:
            col_order += order_dict[row]
        except KeyError:
            pass
    return col_order


def plot_confusion_matrix(matrix, ax, row_order, col_order='auto', cmap='viridis',
                          reverse_diagnal=False, **heatmap_kws):
    if col_order == 'auto':
        col_order = auto_order(matrix, row_order)
    if reverse_diagnal:
        col_order = col_order[::-1]
    nrows, ncols = matrix.shape
    ordered_matrix = matrix.loc[row_order[::-1], col_order]
    sns.heatmap(ordered_matrix, ax=ax, cmap=cmap, **heatmap_kws)
    ax.set_ylim(-0.5, nrows + 0.5)
    return ordered_matrix
