"""
Group visualization

Input:
- Tidy Data: row is cell/group, groupby columns (x, y), groupby orders (x, y)
- Matrix Data: row is x|y, col is y|x. groupby and groupby orders provide separately
- Dendrogram: provided separately

Plot Type:
- Dotplot/Matrixplot: tidy data, single axes, multi visual map (size, color, marker)
- Heatmap: matrix data, single ax, single visual map
- Stack-violin: matrix data, multi axes, single visual map
- TracePlot: matrix data, multi axes, single visual map

Appendix:
- dendrogram
- category color blocks: ax, orientation, order
- category brackets
"""
from .tree import dendrogram
from .color import level_one_palette
import seaborn as sns
import numpy as np


def fancy_heatmap(ax, tidy_data, row_col, col_col,
                  size_col, sig_cutoff, color_col, palette='viridis',
                  row_order=None, col_order=None,
                  heatmap_scatter_kws=None, sig_scatter_kws=None):
    data = tidy_data.copy()
    data[row_col] = data[row_col].astype('category')
    data[col_col] = data[col_col].astype('category')

    if row_order is None:
        row_dendro = dendrogram(data, groupby=row_col)
        row_order = row_dendro['dendrogram_info']['ivl']
    data['row_i'] = data[row_col].map(lambda i: row_order.index(i))

    if col_order is None:
        col_dendro = dendrogram(data, groupby=col_col, cor_method='spearman')
        col_order = col_dendro['dendrogram_info']['ivl']
    data['col_i'] = data[col_col].map(lambda i: col_order.index(i))
    data['sig_marker'] = data[size_col] > sig_cutoff

    # plot color blocks
    _heatmap_scatter_kws = dict(hue_norm=(-1, 1), size_norm=(sig_cutoff / 2, sig_cutoff * 2), sizes=(20, 800),
                                legend=None, marker='s')
    if heatmap_scatter_kws is not None:
        _heatmap_scatter_kws.update(heatmap_scatter_kws)
    sns.scatterplot(x='row_i', y='col_i',
                    data=data, palette=palette,
                    hue=color_col, size=size_col, **_heatmap_scatter_kws, ax=ax)
    # plot sig marker
    _sig_scatter_kws = dict(color='white', linewidth=1, s=100, marker='+')
    if sig_scatter_kws is not None:
        _sig_scatter_kws.update(sig_scatter_kws)
    sns.scatterplot(x='row_i', y='col_i',
                    data=data[data['sig_marker']], **_sig_scatter_kws, ax=ax)

    ax.set(xlim=(-0.5, data[row_col].unique().size - 0.5),
           yticks=range(data[col_col].unique().size),
           yticklabels=col_order,
           xticks=range(data[row_col].unique().size),
           xticklabels=row_order)
    for label in ax.xaxis.get_ticklabels():
        label.set(rotation=45, rotation_mode='anchor', ha='right')
    sns.despine(ax=ax, left=True, bottom=True)
    ax.set(xlabel=col_col, ylabel=row_col)
    return ax


def stack_bar_plot(ax, data, group_col, hue, palette=None, orient='h'):
    """
    Plot stacked barplot between group_col and count_col. E.g. Sample portion in each cluster
    Parameters
    ----------
    ax
    data
    group_col
    hue
    palette
    orient

    Returns
    -------

    """
    raw_count = data.groupby(group_col)[hue].value_counts().unstack().fillna(0)
    norm_by_row_sum = raw_count.divide(raw_count.sum(axis=1), axis=0)
    row_cum_sum = np.cumsum(norm_by_row_sum, axis=1)

    if palette is None:
        palette = level_one_palette(row_cum_sum.columns)
    elif isinstance(palette, str):
        palette = level_one_palette(row_cum_sum.columns, palette=palette)
    else:
        pass

    for col_name, data in list(row_cum_sum.iteritems())[::-1]:
        if orient == 'h':
            sns.barplot(y=row_cum_sum.index, x=data,
                        label=col_name, color=palette[col_name],
                        orient=orient, ax=ax)
        else:
            sns.barplot(x=row_cum_sum.index, y=data,
                        label=col_name, color=palette[col_name],
                        orient=orient, ax=ax)
    ax.set(ylabel=group_col, xlabel='Portion')
    return ax
