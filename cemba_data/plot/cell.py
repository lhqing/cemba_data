import matplotlib as mpl
import numpy as np
import pandas as pd
import seaborn as sns
from adjustText import adjust_text
from matplotlib.cm import get_cmap
from matplotlib.colors import Normalize, LogNorm
from matplotlib.lines import Line2D
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from sklearn.neighbors import LocalOutlierFactor

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


def _make_tiny_axis_lable(ax, coord_name, arrow_kws=None, fontsize=6):
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


def _text_anno_scatter(data, ax, edge_color=(0.5, 0.5, 0.5, 0.2), face_color=(0.8, 0.8, 0.8, 0.2), palette=None,
                       dodge_text=False, anno_col='text_anno', text_anno_kws=None, dodge_kws=None):
    """Add text annotation to a scatter plot"""
    _text_anno_kws = dict(fontsize=10,
                          fontweight='black',
                          horizontalalignment='center',
                          verticalalignment='center')
    if text_anno_kws is not None:
        _text_anno_kws.update(text_anno_kws)
    text_list = []
    for text, sub_df in data.groupby(anno_col):
        text = str(text)
        if text in ['', 'nan']:
            continue
        x, y, *_ = sub_df.median()

        if palette is not None:
            _fc = palette[text]
        else:
            _fc = face_color
        text = ax.text(x, y, text,
                       fontdict=_text_anno_kws,
                       bbox=dict(boxstyle="round",
                                 ec=edge_color,
                                 fc=_fc))
        text_list.append(text)

    if dodge_text:
        _dodge_parms = dict(force_points=(0.02, 0.05),
                            arrowprops=dict(arrowstyle="fancy",
                                            fc=edge_color,
                                            ec="none",
                                            connectionstyle="angle3,angleA=0,angleB=-90"),
                            autoalign='xy')
        if dodge_kws is not None:
            _dodge_parms.update(dodge_kws)
        adjust_text(text_list, x=data['x'], y=data['y'], **_dodge_parms)
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
                             p=2, metric_params=None, contamination=0.1)
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


def _translate_coord_base(coord_base):
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
    return real_coord_name


def categorical_scatter(data, ax, coord_base='umap', scatter_kws=None,  # about basic scatter
                        expand_border_scale=0.1, border_quantile=0.01,  # about xlim and ylim
                        hue=None, palette='tab10',  # about hue
                        text_anno=None, dodge_text=False, dodge_kws=None,  # about text anno
                        text_anno_kws=None,  text_anno_palette=None,  # about text anno
                        show_legend=False, legend_kws=None,  # about legend
                        axis_format='tiny', max_points=5000):  # other adjustment
    data = data.copy()
    # down sample plot data if needed.
    if max_points is not None:
        if data.shape[0] > max_points:
            data = density_based_sample(data, seed=1, size=max_points,
                                        coords=[f'{coord_base}_0',
                                                f'{coord_base}_1'])
    _scatter_kws = dict(linewidth=0, s=7, legend=None, palette=palette)
    if scatter_kws is not None:
        _scatter_kws.update(scatter_kws)
    real_coord_name = _translate_coord_base(coord_base)

    scaled_x = _robust_scale_0_1(data[f'{coord_base}_0'], expand_border_scale, border_quantile)
    scaled_y = _robust_scale_0_1(data[f'{coord_base}_1'], expand_border_scale, border_quantile)
    _data = pd.DataFrame({'x': scaled_x,
                          'y': scaled_y})

    palette_dict = None
    if hue is not None:
        if isinstance(hue, str):
            _data[hue] = data[hue].astype('category')
        else:
            _data['hue'] = hue.astype('category')
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
        _text_anno_scatter(_data[['x', 'y', text_anno]], ax=ax,
                           dodge_text=dodge_text, dodge_kws=dodge_kws, palette=text_anno_palette,
                           anno_col=text_anno, text_anno_kws=text_anno_kws)

    if show_legend and (hue is not None):
        n_hue = len(palette_dict)
        _legend_kws = dict(ncol=(1 if n_hue <= 14 else 2 if n_hue <= 30 else 3), fontsize=10)
        if legend_kws is not None:
            _legend_kws.update(legend_kws)

        handles = []
        labels = []
        exist_hues = _data[hue].unique()
        for hue_name, color in palette_dict.items():
            if hue_name not in exist_hues:
                # skip hue_name that do not appear in the plot
                continue
            handle = Line2D([0], [0], marker='o', color='w',
                            markerfacecolor=color, markersize=_legend_kws['fontsize'])
            handles.append(handle)
            labels.append(hue_name)
        _legend_kws['handles'] = handles
        _legend_kws['labels'] = labels
        ax.legend(**_legend_kws)
    return ax, _data


def continuous_scatter(data, ax, coord_base='umap', scatter_kws=None,
                       expand_border_scale=0.1, border_quantile=0.01,
                       hue=None, hue_norm=None, hue_portion=0.8, cmap='viridis',
                       colorbar=True, colorbar_label_kws=None,
                       size=None, size_norm=None, size_portion=0.8, sizes=None,
                       sizebar=True, sizebar_label_kws=None,
                       text_anno=None, dodge_text=False, dodge_params=None,
                       text_anno_kws=None, text_anno_palette=None,
                       axis_format='tiny', max_points=5000):
    data = data.copy()
    # down sample plot data if needed.
    if max_points is not None:
        if data.shape[0] > max_points:
            data = density_based_sample(data, seed=1, size=max_points,
                                        coords=[f'{coord_base}_0',
                                                f'{coord_base}_1'])

    _scatter_kws = dict(linewidth=0, s=7, legend=None)
    if scatter_kws is not None:
        _scatter_kws.update(scatter_kws)
    real_coord_name = _translate_coord_base(coord_base)

    scaled_x = _robust_scale_0_1(data[f'{coord_base}_0'], expand_border_scale, border_quantile)
    scaled_y = _robust_scale_0_1(data[f'{coord_base}_1'], expand_border_scale, border_quantile)
    _data = pd.DataFrame({'x': scaled_x,
                          'y': scaled_y})
    if hue is not None:
        if isinstance(hue, str):
            _data[hue] = data[hue].astype(float)
        else:
            _data['hue'] = hue.astype(float)
            hue = 'hue'
        if hue_norm is None:
            # get the smallest range that include "hue_portion" of data
            hue_norm = tight_hue_range(_data[hue], hue_portion)
        if isinstance(cmap, str):
            # from here, cmap become colormap object
            cmap = mpl.cm.get_cmap(cmap)
            cmap.set_bad(color=(0.5, 0.5, 0.5, 0.5))
        else:
            if not isinstance(cmap, mpl.cm.ScalarMappable):
                raise TypeError(f'cmap can only be str or ScalarMappable, got {type(cmap)}')
        # cnorm is the normalizer for color
        cnorm = mpl.colors.Normalize(vmin=hue_norm[0],
                                     vmax=hue_norm[1])
    else:
        hue_norm = None
        cnorm = None

    if size is not None:
        if isinstance(size, str):
            _data[size] = data[size].astype(float)
        else:
            _data['size'] = size.astype(float)
            size = 'size'
        # make sure size >= 1  why?
        _data[size] = _data[size] - _data[size].min() + 1

        if size_norm is None:
            # get the smallest range that include "size_portion" of data
            size_norm = tight_hue_range(_data[size], size_portion)

            # snorm is the normalizer for size
            # for size, LogNorm is more suitable
            size_norm = LogNorm(vmin=size_norm[0],
                                vmax=size_norm[1])
        # replace s with sizes
        s = _scatter_kws.pop('s')
        if sizes is None:
            sizes = (1, s)
    else:
        size_norm = None
        sizes = None

    if text_anno is not None:
        if isinstance(text_anno, str):
            _data[text_anno] = data[text_anno]
        else:
            _data['text_anno'] = text_anno
            text_anno = 'text_anno'

    sns.scatterplot(x='x', y='y', data=_data,
                    hue=hue, palette=cmap, hue_norm=cnorm,
                    size=size, sizes=sizes, size_norm=size_norm,
                    ax=ax, **_scatter_kws)

    if text_anno is not None:
        _text_anno_scatter(_data[['x', 'y', text_anno]], ax=ax, dodge_kws=dodge_params, palette=text_anno_palette,
                           anno_col=text_anno, text_anno_kws=text_anno_kws, dodge_text=dodge_text)

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
        _colorbar_label_kws = dict(fontsize=10, label=hue, labelpad=10, rotation=270)
        if colorbar_label_kws is not None:
            _colorbar_label_kws.update(colorbar_label_kws)

        # small ax for colorbar
        # inset_axes add a special ax in figure level: mpl_toolkits.axes_grid1.parasite_axes.AxesHostAxes
        # so there is not way to access cax within ax, we should return cax as well
        # In matplotlib 3.0, ax.inset_axes() can access its own inset_axes, but it is experimental
        cax = inset_axes(ax, width="3%", height="25%",
                         loc='lower right', borderpad=0)
        cax = plot_colorbar(cax, cmap=cmap, cnorm=cnorm, hue_norm=hue_norm,
                            label_kws=_colorbar_label_kws)
        return_axes.append(cax)

    if sizebar and (size is not None):
        _sizebar_label_kws = dict(fontsize=10, ylabel=size, labelpad=20, rotation=270)
        if sizebar_label_kws is not None:
            _sizebar_label_kws.update(sizebar_label_kws)

        # small ax for sizebar
        sax = inset_axes(ax, width="3%", height="25%",
                         loc='upper right', borderpad=0)
        # same as cax, we should return this
        return_axes.append(sax)

        sax = _sizebar(sax)
        # set simple tick label, let the user do custom works
        ticklabels = ['low', 'high']
        sax.yaxis.set(ticks=[0, 1], ticklabels=ticklabels)
        sax.set_ylabel(**_sizebar_label_kws)

    return tuple(return_axes), _data


def scatter_density(ax, data, groupby, coord_base='umap', contour_levels=(-0.5,),
                    color=None, palette=None, contour_kws=None, lof_kws=None):
    data = data.copy()

    if isinstance(groupby, str):
        groupby = data[groupby]
    else:
        data['groupby'] = groupby
        groupby = data['groupby']

    _contour_kws = dict(linewidths=1, levels=contour_levels, linestyles='dashed')
    if contour_kws is not None:
        _contour_kws.update(contour_kws)
    _lof_kws = dict(n_neighbors=20, novelty=True, contamination='auto')
    if lof_kws is not None:
        _lof_kws.update(lof_kws)

    for group, sub_data in data[[f'{coord_base}_0', f'{coord_base}_1']].groupby(groupby):
        xx, yy = np.meshgrid(np.linspace(*ax.get_xlim(), 500), np.linspace(*ax.get_ylim(), 500))
        clf = LocalOutlierFactor(**_lof_kws)
        clf.fit(sub_data)
        z = clf.decision_function(np.c_[xx.ravel(), yy.ravel()])
        z = z.reshape(xx.shape)

        if palette is None:
            if color is None:
                _color = 'lightgray'
            else:
                _color = color
        else:
            _color = palette[group]

        # plot contour line(s)
        ax.contour(xx, yy, z, colors=_color, **_contour_kws)
    return ax
