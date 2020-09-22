import copy

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.colors import Normalize


def plot_on_plate(data,
                  hue,
                  groupby,
                  ncols=4,
                  plate_base=384,
                  figsize_scale=1,
                  row='Row384',
                  col='Col384',
                  vmin=0,
                  vmax=1,
                  aggregation_func=None):
    """
    Plot metadata into 384 or 96 plate view (heatmap)
    Parameters
    ----------
    data
        dataframe contain plate postion and metric used for color
    hue
        int/float column name used as hue
    groupby
        groupby column, typically groupby plate id column(s) to plot each plate separately
    ncols
        number of column for axes, nrows will be calculated accordingly
    plate_base
        {384, 96} size of the plate view
    figsize_scale
        scale of figure size
    row
        column name for rows
    col
        column name for columns
    vmin
        cmap vmin
    vmax
        cmap vmax
    aggregation_func
        apply to reduce rows after groupby if the row is not unique
    """

    if plate_base == 384:
        plate_nrows, plate_ncols = 16, 24

    elif plate_base == 96:
        plate_nrows, plate_ncols = 8, 12
    else:
        raise ValueError(f'Plate base {plate_base} unknown')

    plot_data_list = []
    plate_names = []
    for plate, sub_df in data.groupby(groupby):
        # check if plate base are duplicated
        duplicated = sub_df[[row, col]].duplicated().sum() != 0
        if duplicated:
            if aggregation_func is None:
                raise ValueError(
                    'Row after groupby is not unique, aggregation_func can not be None'
                )
            plot_data = sub_df.groupby([row,
                                        col])[[hue]].apply(aggregation_func)
        else:
            plot_data = sub_df.set_index([row, col])[[hue]]
        # reindex, missing value will keep as NA
        full_index = pd.MultiIndex.from_tuples([(i, j)
                                                for i in range(plate_nrows)
                                                for j in range(plate_ncols)],
                                               names=[row, col])
        plot_data = plot_data.reindex(full_index).reset_index()
        plot_data_list.append(plot_data)
        if isinstance(plate, str):
            plate_names.append(plate)
        else:
            plate_names.append('\n'.join(plate))

    ncols = min(len(plot_data_list), ncols)
    nrows = int(np.ceil(len(plot_data_list) / ncols))
    cbar_frac = 0.06

    fig = plt.figure(figsize=((6.2 * ncols) * (1 + cbar_frac) * figsize_scale,
                              4 * nrows * figsize_scale))
    gs = fig.add_gridspec(nrows, ncols, wspace=0.1)
    cmap = copy.copy(mpl.cm.get_cmap("viridis"))
    cmap.set_under(color='#440154')
    cmap.set_over(color='#FDE725')
    cmap.set_bad(color='#FFFFFF')
    cnorm = Normalize(vmin, vmax)

    for i, (name, data) in enumerate(zip(plate_names, plot_data_list)):
        ax_row = i // ncols
        ax_col = i % ncols

        ax = fig.add_subplot(gs[ax_row, ax_col])
        ax.scatter(
            x=data[col],
            y=data[row],
            # have to do this, otherwise NaN is skipped.
            c=[cmap(cnorm(v)) for v in data[hue]],
            s=100,
            linewidth=1,
            edgecolor='lightgray')
        ax.set(title=name,
               ylabel='',
               ylim=(plate_nrows, -1),
               yticks=list(range(16)),
               yticklabels=[chr(i + 65) for i in range(0, 16)],
               xlabel='',
               xticks=range(24),
               xticklabels=range(1, 25))
        ax.xaxis.set_tick_params(labelsize=8)
        ax.yaxis.set_tick_params(labelsize=8)
        ax.xaxis.tick_top()
    fig.colorbar(mpl.cm.ScalarMappable(norm=cnorm, cmap=cmap),
                 ax=fig.axes,
                 shrink=0.6,
                 fraction=cbar_frac,
                 label=hue)
    return fig, plate_names, plot_data_list


def cutoff_vs_cell_remain(data,
                          xlim_quantile=(0.01, 0.99),
                          distribution_ylim=None,
                          bins=100,
                          kde=False):
    xlim = tuple(np.quantile(data, xlim_quantile))
    x = np.linspace(xlim[0], xlim[1], 500)
    count_list = np.array([(data > i).sum() for i in x])
    original_total_data = data.size
    count_list = count_list / original_total_data * 100
    data = data[(data < xlim[1]) & (data > xlim[0])]

    fig, ax1 = plt.subplots(figsize=(6, 3))
    ax1 = sns.distplot(data, bins=bins, kde=kde, ax=ax1)
    ax1.set_xlim(xlim)
    ax1.set_xlabel(data.name)
    if distribution_ylim is not None:
        ax1.set_ylim(*distribution_ylim)

    ax2 = ax1.twinx()
    ax2.plot(x, count_list, linewidth=1, c='grey')
    ax2.set_ylabel('% of Cell Remained')
    ax2.set(ylim=(0, 100), yticks=range(0, 101, 10))
    ax2.grid(c='lightgray', linestyle='--', linewidth=0.5)
    return fig, xlim
