import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np


def plot_on_plate(data, value_col, groupby, ncols=4,
                  plate_base=384, figsize=(5, 3.5),
                  vmin=0, vmax=1, heatmap_kws=None, aggregation_func=None):
    """
    Plot metadata into 384 or 96 plate view (heatmap)

    Parameters
    ----------
    data
        dataframe contain all kinds of metadata
    value_col
        value to be plotted on plate view
    groupby
        groupby column, typically groupby plate id column(s) to plot each plate separately
    ncols
        number of column for axes, nrows will be calculated accordingly
    plate_base
        {384, 96} size of the plate view
    figsize
        matplotlib.Figure figsize
    vmin
        cmap vmin
    vmax
        cmap vmax
    heatmap_kws
        kws pass to sns.heatmap
    aggregation_func
        apply to reduce rows after groupby if the row is not unique
    """

    heatmap_data_list = []
    heatmap_names = []
    for plate, sub_df in data.groupby(groupby):
        # check if plate base are duplicated
        duplicated = sub_df[[f'{plate_base}_col', f'{plate_base}_row']].duplicated().sum() != 0
        if duplicated:
            if aggregation_func is None:
                raise ValueError('Row after groupby is not unique, aggregation_func can not be None')
            heatmap_data = sub_df.groupby([f'{plate_base}_col', f'{plate_base}_row'])[value_col]\
                .apply(aggregation_func).unstack().fillna(-999)
        else:
            heatmap_data = sub_df.set_index([f'{plate_base}_col', f'{plate_base}_row'])[value_col]\
                .unstack().fillna(-999)
        heatmap_data_list.append(heatmap_data)
        if isinstance(plate, str):
            heatmap_names.append(plate)
        else:
            heatmap_names.append('\n'.join(plate))

    nrows = round(len(heatmap_data_list) / ncols)
    fig, axes = plt.subplots(figsize=(figsize[0] * ncols, figsize[1] * nrows),
                             ncols=ncols,
                             nrows=nrows)
    fig.suptitle(f'{value_col} on {len(heatmap_data_list)}*{plate_base} plates \n Color Range [{vmin}, {vmax}]',
                 fontsize=16)

    cmap = plt.cm.viridis
    cmap.set_under(color='#FFFFFF')
    for heatmap_data, heatmap_name, ax in zip(heatmap_data_list, heatmap_names, np.ravel(axes)):
        sns.heatmap(heatmap_data, vmin=vmin, vmax=vmax,
                    cmap=cmap, ax=ax, **heatmap_kws)
        ax.set_title(heatmap_name)
    fig.tight_layout(rect=[0, 0.05, 1, 0.95])
    return fig, axes


def simple_violin(data, x, y, rotate_x=True):
    fig, ax = plt.subplots(figsize=(5, 3))
    sns.violinplot(x=x, y=y, data=data, ax=ax)
    if rotate_x:
        ax.tick_params(axis='x', rotation=90)
    return fig, ax


def cutoff_vs_cell_remain(data, cutoff_num=1000,
                          xlim_quantile=(0.001, 0.999), distribution_ylim=None,
                          bins=100, kde=False, count_percent=True):
    xlim = tuple(np.quantile(data, xlim_quantile))
    x = np.linspace(xlim[0], xlim[1], cutoff_num)
    count_list = np.array([(data > i).sum() for i in x])
    original_total_data = data.size
    if count_percent:
        count_list = count_list / original_total_data
    data = data[(data < xlim[1]) & (data > xlim[0])]

    fig, ax1 = plt.subplots()
    ax1.grid()
    ax1 = sns.distplot(data, bins=bins, kde=kde,
                       ax=ax1)
    ax1.set_xlim(xlim)
    ax1.set_xlabel('Data Point Count')
    if distribution_ylim is not None:
        ax1.set_ylim(*distribution_ylim)

    ax2 = ax1.twinx()
    sns.scatterplot(x=x, y=count_list, linewidth=0, legend=None,
                    s=5, hue=x, palette='viridis', ax=ax2)
    ax2.set_ylabel('% of Data Pass Filter')
    return fig, (ax1, ax2)


def success_vs_fail(data, filter_col, filter_cutoff, x, y):
    use_data = data.copy()
    use_data['filter'] = (data[filter_col] > filter_cutoff).apply(lambda i: 'Success' if i else 'Fail')
    fig, ax = plt.subplots()
    sns.violinplot(x=x, y=y, data=use_data, hue='filter',
                   ax=ax)
    ax.tick_params(axis='x', rotation=90)
    return fig, ax
