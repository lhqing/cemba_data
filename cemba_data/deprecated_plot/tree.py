import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def generate_tree_dict(data):
    """
    Helper function for echarts tree structures

    Parameters
    ----------
    data

    Returns
    -------

    """
    children = []
    if data.shape[1] == 1:
        cell_counts = data.iloc[:, 0].value_counts()
        for child, count in cell_counts.items():
            children.append({
                'name': child,
                'value': count
            })
    else:
        for child, sub_df in data.groupby(data.columns[0]):
            count = sub_df.shape[0]
            children.append({
                'name': child,
                'value': count,
                'children': generate_tree_dict(sub_df.iloc[:, 1:])
            })
    return children


def dendrogram(data_df, groupby,
               cor_method='pearson',
               linkage_method='complete',
               linkage_metric='euclidean',
               dendogram_kws=None):
    """
    Computes a hierarchical clustering for the given `groupby` categories on data_df.

    Average values of columns in data_df are used to compute a correlation matrix.

    .. note::
        The computation of the hierarchical clustering is based on predefined groups and not
        per cell. The correlation matrix is computed using by default pearson but other methods
        are available.

    Parameters
    ----------
    data_df
    groupby
    cor_method
    linkage_method
    linkage_metric
    dendogram_kws

    Returns
    -------
    dendogram_dict
    """
    data_df = data_df.copy()

    if isinstance(groupby, pd.Series):
        groupby = groupby.astype('category')
    elif isinstance(groupby, str):
        groupby = data_df.pop(groupby)
    categories = groupby.cat.categories

    # aggregate values within categories using 'mean'
    mean_df = data_df.groupby(groupby).mean()

    import scipy.cluster.hierarchy as sch

    corr_matrix = mean_df.T.corr(method=cor_method)
    z_var = sch.linkage(corr_matrix,
                        method=linkage_method,
                        metric=linkage_metric,
                        optimal_ordering=False)

    if dendogram_kws is None:
        dendogram_kws = {}
    dendro_info = sch.dendrogram(z_var, labels=categories, no_plot=True, **dendogram_kws)

    # order of groupby categories
    categories_idx_ordered = dendro_info['leaves']

    dendogram_dict = {'groupby': groupby,  # groupby series
                      'cor_method': cor_method,  # correlation method
                      'linkage_method': linkage_method,  # linkage method
                      'correlation_matrix': corr_matrix.values,  # correlation matrix used for linkage
                      'linkage': z_var,  # linkage data
                      'dendrogram_info': dendro_info,  # core dendogram data
                      'categories_idx_ordered': categories_idx_ordered,  # ordered categories
                      }
    return dendogram_dict


def artificial_dendrogram():
    """
    This function take known hierarchical groups and build artificial dendrogram.
    The leafs will be

    Returns
    -------

    """
    raise NotImplementedError


def _translate_pos(pos_list, new_ticks, old_ticks):
    """
    transforms the dendrogram coordinates to a given new position.
    The xlabel_pos and orig_ticks should be of the same
    length.

    This is mostly done for the heatmap case, where the position of the
    dendrogram leaves needs to be adjusted dependening on the size of the category.

    Parameters
    ----------
    pos_list
        list of dendrogram positions that should be translated
    new_ticks
        sorted list of goal tick positions (e.g. [0,1,2,3] )
    old_ticks
        sorted list of original tick positions (e.g. [5, 15, 25, 35]), This list is
        usually the default position used by scipy.cluster.hierarchy.dendrogram`

    Returns
    -------
    translated list of positions
    """
    # of given coordinates.
    new_xs = []
    for x_val in pos_list:
        if x_val in old_ticks:
            new_x_val = new_ticks[old_ticks.index(x_val)]
        else:
            # find smaller and bigger indices
            idx_next = np.searchsorted(old_ticks, x_val, side="left")
            idx_prev = idx_next - 1
            old_min = old_ticks[idx_prev]
            old_max = old_ticks[idx_next]
            new_min = new_ticks[idx_prev]
            new_max = new_ticks[idx_next]
            new_x_val = ((x_val - old_min) / (old_max - old_min)) * (new_max - new_min) + new_min
        new_xs.append(new_x_val)
    return new_xs


def _plot_dendrogram(ax, dendogram_dict, orientation='right', remove_labels=True,
                     ticks=None):
    """
    Plots a dendrogram on the given ax using the precomputed dendrogram information
    stored in .uns[dendrogram_key]
    """

    dendro_info = dendogram_dict['dendrogram_info']
    leaves = dendro_info["ivl"]
    icoord = np.array(dendro_info['icoord'])
    dcoord = np.array(dendro_info['dcoord'])

    orig_ticks = np.arange(5, len(leaves) * 10 + 5, 10).astype(float)
    # check that ticks has the same length as orig_ticks
    if ticks is not None and len(orig_ticks) != len(ticks):
        print("ticks argument does not have the same size as orig_ticks. The argument will be ignored")
        ticks = None

    for xs, ys in zip(icoord, dcoord):
        if ticks is not None:
            xs = _translate_pos(xs, ticks, orig_ticks)
        if orientation in ['right', 'left']:
            ax.plot(ys, xs, color='#555555')
        else:
            ax.plot(xs, ys, color='#555555')

    ax.tick_params(bottom=False, top=False, left=False, right=False)
    ticks = ticks if ticks is not None else orig_ticks
    if orientation in ['right', 'left']:
        ax.set_yticks(ticks)
        ax.set_yticklabels(leaves, fontsize='small', rotation=0)
        ax.tick_params(labelbottom=False, labeltop=False)
        if orientation == 'left':
            xmin, xmax = ax.get_xlim()
            ax.set_xlim(xmax, xmin)
            ax.tick_params(labelleft=False, labelright=True)
    else:
        ax.set_xticks(ticks)
        ax.set_xticklabels(leaves, fontsize='small', rotation=90)
        ax.tick_params(labelleft=False, labelright=False)
        if orientation == 'bottom':
            ymin, ymax = ax.get_ylim()
            ax.set_ylim(ymax, ymin)
            ax.tick_params(labeltop=True, labelbottom=False)

    if remove_labels:
        ax.tick_params(labelbottom=False, labeltop=False, labelleft=False, labelright=False)

    ax.grid(False)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    return ax


def _plot_categories_as_color_blocks(ax, groupby, order, colors=None, orientation='left', cmap_name='tab20'):
    """
    Plots categories as colored blocks. If orientation is 'left', the categories are plotted vertically, otherwise
    they are plotted horizontally.

    Parameters
    ----------
    ax : matplotlib ax
    colors : list of valid color names optional (default: `None`)
        Color to use for each category.
    orientation : `str`, optional (default: `left`)
    cmap_name : `str`
        Name of colormap to use, in case colors is None

    Returns
    -------
    ticks position, labels, colormap
    """
    groupby = pd.concat([groupby[groupby == i] for i in order])

    groupby = groupby.astype('category')
    if order is None:
        order = groupby.cat.categories
    else:
        order = pd.Index(order)
    group_label = groupby.name

    from matplotlib.colors import ListedColormap, BoundaryNorm
    if colors is None:
        groupby_cmap = plt.get_cmap(cmap_name)
    else:
        groupby_cmap = ListedColormap(colors, group_label + '_cmap')
    norm = BoundaryNorm(np.arange(groupby_cmap.N + 1) - .5, groupby_cmap.N)

    # determine group_label positions such that they appear
    # centered next/below to the color code rectangle assigned to the category
    value_sum = 0
    ticks = []  # list of centered position of the labels
    labels = []
    label2code = {}  # dictionary of numerical values asigned to each label
    for code, (label, value) in enumerate(groupby.value_counts(sort=False)
                                                  .reindex(order)
                                                  .iteritems()):
        ticks.append(value_sum + (value / 2))
        labels.append(label)
        value_sum += value
        label2code[label] = code

    if orientation == 'left':
        ax.imshow(np.matrix(groupby.map(label2code).values).T,
                  aspect='auto', cmap=groupby_cmap, norm=norm)
        if len(labels) > 1:
            ax.set_yticks(ticks)
            ax.set_yticklabels(labels)
        # remove y ticks
        ax.tick_params(axis='y', left=False, labelsize='small')
        # remove x ticks and labels
        ax.tick_params(axis='x', bottom=False, labelbottom=False)
        ax.set_ylabel(group_label)
    else:
        ax.imshow(np.matrix(groupby.map(label2code).values),
                  aspect='auto', cmap=groupby_cmap, norm=norm)
        if len(labels) > 1:
            ax.set_xticks(ticks)
            if max([len(x) for x in labels]) < 3:
                # if the labels are small do not rotate them
                rotation = 0
            else:
                rotation = 90
            ax.set_xticklabels(labels, rotation=rotation)
        # remove x ticks
        ax.tick_params(axis='x', bottom=False, labelsize='small')
        # remove y ticks and labels
        ax.tick_params(axis='y', left=False, labelleft=False)
        ax.set_xlabel(group_label)

    ax.grid(False)
    # remove surrounding lines
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)

    return ticks, labels, groupby_cmap, norm


def _plot_groups_brackets(ax, group_positions, group_labels,
                          left_adjustment=-0.3, right_adjustment=0.3,
                          rotation=None, orientation='top'):
    """
    Draws brackets that represent groups of genes on the give axis.
    For best results, this axis is located on top of an image whose
    x axis contains gene names.

    The ax should share the x axis with the main ax.

    Eg: ax = fig.add_subplot(axs[0, 0], sharex=dot_ax)

    This function is used by dotplot, heatmap etc.

    Parameters
    ----------
    ax : matplotlib axis
        In this axis the gene marks are drawn
    group_positions : list of `tuples`
        Each item in the list, should contain the start and end position that the
        bracket should cover.
        Eg. [(0, 4), (5, 8)] means that there are two brackets, one for the var_names (eg genes)
        in positions 0-4 and other for positions 5-8
    group_labels
        List of group labels
    left_adjustment : `float`
        adjustment to deprecated_plot the bracket start slightly before or after the first gene position.
        If the value is negative the start is moved before.
    right_adjustment : `float`
        adjustment to deprecated_plot the bracket end slightly before or after the last gene position
        If the value is negative the start is moved before.
    rotation : `float` (default None)
        rotation degrees for the labels. If not given, small labels (<4 characters) are not
        rotated, otherwise, they are rotated 90 degrees
    orientation : `str` (default `top`)
        location of the brackets. Either `top` or `right`
    Returns
    -------
    None
    """
    import matplotlib.patches as patches
    from matplotlib.path import Path

    # get the 'brackets' coordinates as lists of start and end positions
    left = [x[0] + left_adjustment for x in group_positions]
    right = [x[1] + right_adjustment for x in group_positions]

    # verts and codes are used by PathPatch to make the brackets
    verts = []
    codes = []
    if orientation == 'top':
        # rotate labels if any of them is longer than 4 characters
        if rotation is None and group_labels is not None and len(group_labels) > 0:
            if max([len(x) for x in group_labels]) > 4:
                rotation = 90
            else:
                rotation = 0
        for idx in range(len(left)):
            verts.append((left[idx], 0))  # lower-left
            verts.append((left[idx], 0.6))  # upper-left
            verts.append((right[idx], 0.6))  # upper-right
            verts.append((right[idx], 0))  # lower-right
            codes.append(Path.MOVETO)
            codes.append(Path.LINETO)
            codes.append(Path.LINETO)
            codes.append(Path.LINETO)

            group_x_center = left[idx] + float(right[idx] - left[idx]) / 2
            ax.text(group_x_center, 1.1, group_labels[idx], ha='center',
                    va='bottom', rotation=rotation)
    else:
        top = left
        bottom = right
        for idx in range(len(top)):
            verts.append((0, top[idx]))  # upper-left
            verts.append((0.15, top[idx]))  # upper-right
            verts.append((0.15, bottom[idx]))  # lower-right
            verts.append((0, bottom[idx]))  # lower-left
            codes.append(Path.MOVETO)
            codes.append(Path.LINETO)
            codes.append(Path.LINETO)
            codes.append(Path.LINETO)

            diff = bottom[idx] - top[idx]
            group_y_center = top[idx] + float(diff) / 2
            if diff * 2 < len(group_labels[idx]):
                # cut label to fit available space
                group_labels[idx] = group_labels[idx][:int(diff * 2)] + "."
            ax.text(0.6, group_y_center, group_labels[idx], ha='right',
                    va='center', rotation=270, fontsize='small')

    path = Path(verts, codes)
    patch = patches.PathPatch(path, facecolor='none', lw=1.5)
    ax.add_patch(patch)

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.grid(False)
    # remove y ticks
    ax.tick_params(axis='y', left=False, labelleft=False)
    # remove x ticks and labels
    ax.tick_params(axis='x', bottom=False, labelbottom=False, labeltop=False)
    return ax
