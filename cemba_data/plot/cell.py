import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib as mpl
from sklearn.neighbors import LocalOutlierFactor
from matplotlib.cm import get_cmap
from matplotlib.colors import Normalize, LogNorm
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from .color import level_one_palette
from .utilities import tight_hue_range, plot_colorbar


def _robust_scale_0_1(coord, expand_border_scale=0.1, border_quantile=0.01):
    """linear scale any coord range to [0, 1]"""
    robust_min = coord.quantile(border_quantile)
    robust_max = coord.quantile(1 - border_quantile)
    expand_border = (robust_max - robust_min) * expand_border_scale
    true_min = robust_min - expand_border
    true_max = robust_max + expand_border
    true_delta = true_max - true_min
    true_coord = (coord - true_min) / true_delta
    return true_coord


def _make_tiny_axis_lable(ax, coord_name, arrow_kws=None, fontsize=5):
    """this function assume coord is [0, 1]"""
    # clean ax axises
    ax.set(xticks=[], yticks=[], xlabel=None, ylabel=None)
    sns.despine(ax=ax, left=True, bottom=True)

    _arrow_kws = dict(width=0.003, linewidth=0, color='black')
    if arrow_kws is not None:
        _arrow_kws.update(arrow_kws)

    ax.arrow(0.06, 0.06, 0, 0.06, **_arrow_kws)
    ax.arrow(0.06, 0.06, 0.06, 0, **_arrow_kws)
    ax.text(0.09, 0.03, f'{coord_name} 1',
            fontdict=dict(fontsize=fontsize,
                          horizontalalignment='center',
                          verticalalignment='center'))
    ax.text(0.03, 0.09, f'{coord_name} 2',
            fontdict=dict(fontsize=fontsize, rotation=90,
                          horizontalalignment='center',
                          verticalalignment='center'))
    return


def _dodge(point_array, max_movement=0.05):
    """dodge point to prevent overlap (very naive)"""
    # add some noise to prevent exact same points
    noise = (np.random.random(size=point_array.shape) - 0.5) / 10000
    point_array += noise

    force_list = []
    for point in point_array:
        delta = point - point_array
        delta = delta[delta.sum(axis=1) != 0]
        # scale y axis to prefer movement on y
        delta[:, 1] = delta[:, 1] / 2
        distance = np.sqrt(np.sum(delta ** 2, axis=1))
        force = (delta / (distance ** 5)[:, None]).sum(axis=0)  # exaggerate distance**4
        force_list.append(force)
    force_array = np.array(force_list)
    movement_array = force_array / force_array.max() * max_movement  # max movement 0.05
    adj_points = point_array + movement_array
    return adj_points


def _text_anno_scatter(data, ax, dodge, anno_col='text_anno', text_anno_kws=None):
    """Add text annotation to a scatter plot"""
    _text_anno_kws = dict(fontsize=6,
                          fontweight='black',
                          horizontalalignment='center',
                          verticalalignment='center')
    if text_anno_kws is not None:
        _text_anno_kws.update(text_anno_kws)

    if dodge is None:
        # text annotation
        for text, sub_df in data.groupby(anno_col):
            if (not isinstance(text, str)) or (text == ''):
                continue
            x, y, *_ = sub_df.median()
            ax.text(x, y, text,
                    fontdict=_text_anno_kws,
                    bbox=dict(boxstyle="round",
                              ec=(0.5, 0.5, 0.5, 0.2),
                              fc=(0.8, 0.8, 0.8, 0.2)))
    else:
        # text annotation
        points = []
        texts = []
        for text, sub_df in data.groupby(anno_col):
            x, y = sub_df.median()
            points.append((x, y))
            texts.append(text)
        points = np.array(points)

        # force layout to make label separate
        adj_points = _dodge(points, max_movement=dodge)
        for text, (x, y), (adj_x, adj_y) in zip(texts, points, adj_points):
            ax.text(adj_x, adj_y, text,
                    fontdict=_text_anno_kws,
                    bbox=dict(boxstyle="round",
                              ec=(0.5, 0.5, 0.5, 0.4),
                              fc=(0.9, 0.9, 0.9, 0.6)))
            adj_distance = np.sqrt((adj_x - x) ** 2 + (adj_y - y) ** 2)
            if adj_distance > 0.05:
                ax.arrow(adj_x, adj_y, x - adj_x, y - adj_y,
                         color=(0.5, 0.5, 0.5, 0.8),
                         width=0.003, linewidth=0)
    return


def _sizebar(ax, color=(0.5, 0.5, 0.5), lw=0.5):
    """plot a triangle sizebar"""
    from matplotlib.patches import PathPatch
    from matplotlib.path import Path

    # make triangle
    verts = [(1, 0), (1, 1), (0, 1), (1, 0)]
    codes = [Path.MOVETO] + (len(verts) - 1) * [Path.LINETO]
    tri = Path(verts, codes)
    triangle = PathPatch(tri, facecolor=color, lw=lw)
    ax.add_patch(triangle)

    sns.despine(ax=ax, bottom=True, left=True, right=True)
    ax.yaxis.tick_right()
    ax.yaxis.set_label_position('right')

    ax.set(ylim=(0, 1), xlim=(-0.1, 1), xticks=[], yticks=[])
    return ax


def density_based_sample(data, coords, portion=None, size=None, seed=None):
    clf = LocalOutlierFactor(n_neighbors=20, algorithm='auto',
                             leaf_size=30, metric='minkowski',
                             p=2, metric_params=None)
    data_coords = data[coords]
    clf.fit(data_coords)
    # original score is negative, the larger the denser
    density_score = clf.negative_outlier_factor_
    delta = density_score.max() - density_score.min()
    # density score to probability: the denser the less probability to be picked up
    probability_score = 1 - (density_score - density_score.min()) / delta
    probability_score = np.sqrt(probability_score)
    probability_score = probability_score / probability_score.sum()

    if size is not None:
        pass
    elif portion is not None:
        size = int(data_coords.index.size * portion)
    else:
        raise ValueError('Either portion or size should be provided.')
    if seed is not None:
        np.random.seed(seed)
    selected_cell_index = np.random.choice(data_coords.index,
                                           size=size,
                                           replace=False,
                                           p=probability_score)
    return data.reindex(selected_cell_index)


def categorical_scatter(data, ax, coord_base='tsne', hue=None, palette='tab10',
                        expand_border_scale=0.1, border_quantile=0.01, text_anno=None, dodge=None,
                        scatter_kws=None, text_anno_kws=None, axis_format='tiny',
                        show_legend=False, legend_kws=None, max_points=None):
    # down sample plot data if needed.
    if max_points is not None:
        if data.shape[0] > max_points:
            data = density_based_sample(data, seed=1, size=max_points,
                                        coords=[f'{coord_base}_0',
                                                f'{coord_base}_1'])
    _scatter_kws = dict(linewidth=0, s=7, legend=None, palette=palette)
    if scatter_kws is not None:
        _scatter_kws.update(scatter_kws)

    _coord_base = coord_base
    coord_base = coord_base.lower()
    if coord_base == 'tsne':
        real_coord_name = 'tSNE'
    elif coord_base == 'umap':
        real_coord_name = 'UMAP'
    elif coord_base in ['pca', 'pc']:
        real_coord_name = 'PC'
    else:
        real_coord_name = _coord_base

    scaled_x = _robust_scale_0_1(data[f'{coord_base}_0'], expand_border_scale, border_quantile)
    scaled_y = _robust_scale_0_1(data[f'{coord_base}_1'], expand_border_scale, border_quantile)
    _data = pd.DataFrame({'x': scaled_x,
                          'y': scaled_y})
    if hue is not None:
        if isinstance(hue, str):
            _data[hue] = data[hue]
        else:
            _data['hue'] = hue
            hue = 'hue'

        # deal with color palette
        palette = _scatter_kws['palette']
        if isinstance(palette, str) or isinstance(palette, list):
            palette_dict = level_one_palette(_data[hue], order=None, palette=palette)
        elif isinstance(palette, dict):
            palette_dict = palette
        else:
            raise TypeError(f'Palette can only be str, list or dict, '
                            f'got {type(palette)}')
        _scatter_kws['palette'] = palette_dict

    if text_anno is not None:
        if isinstance(text_anno, str):
            _data[text_anno] = data[text_anno]
        else:
            _data['text_anno'] = text_anno
            text_anno = 'text_anno'

    sns.scatterplot(x='x', y='y', hue=hue,
                    data=_data, ax=ax, **_scatter_kws)

    # clean axis
    if axis_format == 'tiny':
        _make_tiny_axis_lable(ax, real_coord_name)
    elif axis_format == 'despine':
        sns.despine(ax=ax)
    elif (axis_format == 'empty') or axis_format is None:
        ax.set(xticks=[], yticks=[], xlabel=None, ylabel=None)
        sns.despine(ax=ax, left=True, bottom=True)
    elif axis_format == 'normal':
        pass
    else:
        raise ValueError('axis_format need to be one of these: tiny, despine, empty.')

    if text_anno:
        _text_anno_scatter(_data[['x', 'y', text_anno]], ax=ax, dodge=dodge,
                           anno_col=text_anno, text_anno_kws=text_anno_kws)

    if show_legend:
        hue_length = len(_data[hue])
        _legend_kws = dict(ncol=(1 if hue_length <= 14 else 2 if hue_length <= 30 else 3))
        if legend_kws is not None:
            _legend_kws.update(legend_kws)
        ax.legend(**_legend_kws)
    return ax


def continuous_scatter(data, ax, coord_base='tsne',
                       sample_portion=None, sample_size=None, seed=0,
                       expand_border_scale=0.1, border_quantile=0.01, text_anno=None, dodge=None,
                       hue=None, hue_portion=0.8, cmap='viridis', colorbar=True,
                       size=None, size_portion=0.8, label_fontsize=5, sizebar=True,
                       scatter_kws=None, text_anno_kws=None, axis_format='tiny', max_points=None):
    # down sample plot data if needed.
    if max_points is not None:
        if data.shape[0] > max_points:
            data = density_based_sample(data, seed=1, size=max_points,
                                        coords=[f'{coord_base}_0',
                                                f'{coord_base}_1'])

    if (sample_portion is not None) or (sample_size is not None):
        data = density_based_sample(data.copy(),
                                    coords=[f'{coord_base}_0', f'{coord_base}_1'],
                                    portion=sample_portion, size=sample_size, seed=seed)

    _coord_base = coord_base
    coord_base = coord_base.lower()
    if coord_base == 'tsne':
        real_coord_name = 'tSNE'
    elif coord_base == 'umap':
        real_coord_name = 'UMAP'
    elif coord_base in ['pca', 'pc']:
        real_coord_name = 'PC'
    else:
        real_coord_name = _coord_base

    _scatter_kws = dict(linewidth=0, s=7, legend=None)
    if scatter_kws is not None:
        _scatter_kws.update(scatter_kws)

    scaled_x = _robust_scale_0_1(data[f'{coord_base}_0'], expand_border_scale, border_quantile)
    scaled_y = _robust_scale_0_1(data[f'{coord_base}_1'], expand_border_scale, border_quantile)
    _data = pd.DataFrame({'x': scaled_x,
                          'y': scaled_y})
    if hue is not None:
        if isinstance(hue, str):
            _data[hue] = data[hue]
        else:
            _data['hue'] = hue
            hue = 'hue'
        # get the smallest range that include "hue_portion" of data
        hue_norm = tight_hue_range(_data[hue], hue_portion)
        # from here, cmap become colormap object
        cmap = mpl.cm.get_cmap(cmap)
        # cnorm is the normalizer for color
        cnorm = mpl.colors.Normalize(vmin=hue_norm[0],
                                     vmax=hue_norm[1])
        cmap.set_bad(color=(0.5, 0.5, 0.5, 0.5))
    else:
        hue_norm = None
        cnorm = None

    if size is not None:
        if isinstance(size, str):
            _data[size] = data[size]
        else:
            _data['size'] = size
            size = 'size'
        # get the smallest range that include "size_portion" of data
        size_norm = tight_hue_range(_data[hue], size_portion)
        # snorm is the normalizer for size
        # for size, LogNorm is more suitable
        snorm = LogNorm(vmin=size_norm[0],
                        vmax=size_norm[1])
        # replace s with sizes
        s = _scatter_kws.pop('s')
        sizes = (1, s)
    else:
        size_norm = None
        sizes = None
        snorm = None

    if text_anno is not None:
        if isinstance(text_anno, str):
            _data[text_anno] = data[text_anno]
        else:
            _data['text_anno'] = text_anno
            text_anno = 'text_anno'

    sns.scatterplot(x='x', y='y', data=_data,
                    hue=hue, palette=cmap, hue_norm=cnorm,
                    size=size, sizes=sizes, size_norm=snorm,
                    ax=ax, **_scatter_kws)

    # clean axis
    if axis_format == 'tiny':
        _make_tiny_axis_lable(ax, real_coord_name)
    elif axis_format == 'despine':
        sns.despine(ax=ax)
    elif (axis_format == 'empty') or axis_format is None:
        ax.set(xticks=[], yticks=[], xlabel=None, ylabel=None)
        sns.despine(ax=ax, left=True, bottom=True)
    elif axis_format == 'normal':
        pass
    else:
        raise ValueError('axis_format need to be one of these: tiny, despine, empty.')

    return_axes = [ax]

    if colorbar and (hue is not None):
        # small ax for colorbar
        # inset_axes add a special ax in figure level: mpl_toolkits.axes_grid1.parasite_axes.AxesHostAxes
        # so there is not way to access cax within ax, we should return cax as well
        # In matplotlib 3.0, ax.inset_axes() can access its own inset_axes, but it is experimental
        cax = inset_axes(ax, width="3%", height="25%",
                         loc='lower right', borderpad=0)

        cax = plot_colorbar(cax, cmap=cmap, cnorm=cnorm, hue_norm=hue_norm,
                            label_kws=dict(label='Normalized mCH%', labelpad=8,
                                           rotation=270, fontsize=label_fontsize))
        return_axes.append(cax)

    if sizebar and (size is not None):
        # small ax for sizebar
        sax = inset_axes(ax, width="3%", height="25%",
                         loc='upper right', borderpad=0)
        # same as cax, we should return this
        return_axes.append(sax)

        sax = _sizebar(sax)
        sax.set_yticks([0, 0.5, 1])  # the y axis for sizebar is fixed
        ticklabels = list(map(lambda i: f'{i:.0f}', [size_norm[0], sum(size_norm) / 2, size_norm[1]]))
        sax.set_yticklabels(ticklabels)
        sax.set_ylabel('cov', labelpad=10, rotation=270, fontsize=label_fontsize)
        sax.yaxis.set_tick_params(labelsize=label_fontsize, width=0.5, pad=1)

    if text_anno is not None:
        _text_anno_scatter(_data[['x', 'y', text_anno]], ax=ax, dodge=dodge,
                           anno_col=text_anno, text_anno_kws=text_anno_kws)
        # TODO adjust label color, turn white when background is dark
    return tuple(return_axes)
