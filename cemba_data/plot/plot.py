import holoviews as hv
import seaborn as sns
import matplotlib.pyplot as plt
from math import log2


def sample_data(df, sample=None):
    if sample is not None and sample < df.shape[0]:
        print(
            f'Sampling {sample} data from total. '
            f'Set sample=None if don\'t want to sample, '
            f'which may cause interactive plots huge and slow.')
        new_df = df.sample(sample)
        # TODO: fancy sampling, remain/exclude extreme value, keep data structure
        return new_df
    else:
        return df


def pair_plot(data, hue=None, interactive=False, hue_overlay=False,
              sample=500, alpha=0.8, size=3, height=2.5, aspect=1,
              plot_space=0.1, sns_kws=None, palette=None):
    # hue is in the data col, all columns are kdims
    if interactive:
        # sampling to speed up and prevent over-plotting
        data = sample_data(data, sample)

        # only remain int/float columns + hul_col
        number_cols = []
        removed_cols = []
        for col_name, col_data in data.iteritems():
            if 'int' in col_data.dtype.name or 'float' in col_data.dtype.name:
                number_cols.append(col_name)
            elif col_name == hue:
                number_cols.append(col_name)
            else:
                removed_cols.append(col_name)
        data = data[number_cols]
        if len(removed_cols) != 0:
            print('Nonnumerical columns are not displayed in interactive mode:', removed_cols)

        table = hv.Table(data)
        if hue_overlay:
            plot_data = table.groupby(hue, container_type=hv.NdOverlay)
        else:
            plot_data = table
        g = hv.operation.gridmatrix(plot_data,
                                    chart_type=hv.Points)
        g = g.opts({'Points': {'plot': {'tools': ['box_select', 'lasso_select', 'hover']}}})
    else:
        # do not sample data for non-interactive plots
        if sns_kws is None:
            sns_kws = {}
        g = sns.PairGrid(data=data, hue=hue, palette=palette, height=height, aspect=aspect, **sns_kws)
        g = g.map_offdiag(plt.scatter, edgecolor='w', s=size * 10, alpha=alpha)
        g = g.map_diag(sns.kdeplot, lw=2, legend=False, shade=True)
        g.fig.subplots_adjust(wspace=plot_space, hspace=plot_space)
        g.add_legend()
    return g


def large_categorical_scatter(x, y, color, data=None, sample=3000, size=5,
                              height=600, interactive=True, show_legend=False):
    if interactive:
        # sampling to speed up and prevent over-plotting
        data = sample_data(data, sample)
        key_dimensions = [x, color]
        value_dimensions = [y]
        data_table = hv.Table(data, key_dimensions, value_dimensions)

        opts_dict = {'Scatter': {'plot': dict(color_index=color,
                                              width=height,
                                              height=height,
                                              show_legend=show_legend),
                                 'style': dict(size=size)}}
        _scatter = data_table.to.scatter(x, y)
        g = _scatter.overlay(color).opts(opts_dict)
    else:
        # TODO
        g = None
    return g


def large_continuous_scatter(x, y, color, data=None, sample=3000, size=5,
                             height=600, interactive=True, add_hist=True):
    if interactive:
        # sampling to speed up and prevent over-plotting
        data = sample_data(data, sample)
        data = data[[x, y, color]]
        opts_dict = {'Points': {'plot': dict(color_index=color,
                                             scaling_factor=50,
                                             tools=['lasso_select'],
                                             height=height,
                                             width=int(height*1.2),
                                             colorbar=True),
                                'style': dict(cmap='viridis',
                                              size=size)}}
        # can't set alpha here, conflict with point selection
        g = hv.Points(data, vdims=[color]).opts(opts_dict)
        if add_hist:
            g = g.hist()
    else:
        # TODO
        g = None
    return g


def reduced_scatter_plot(data, x, y, category_col,
                         name_col='auto_gen', color_col=None, size_col='item_count',
                         return_point=True, return_text=True,
                         font_size=1.5, height=400, alpha=0.8):
    category_df = data.groupby(category_col).median()
    category_df['item_count'] = data[category_col].value_counts()

    if name_col == 'auto_gen':
        name_col = 'name'
        category_df['name'] = category_df.index.map(lambda n: f'{category_col}_{n}')

    key_dimensions = [x]
    value_dimensions = [y, size_col, name_col]
    if color_col is not None:
        value_dimensions.append(color_col)

    if not return_point:
        # still plot points, but invisible
        alpha = 0
    opts_dict = {'Scatter': {'plot': dict(color_index=color_col,
                                          size_index=size_col,
                                          tools=['hover'],
                                          height=height,
                                          width=height,
                                          show_legend=False),
                             'style': dict(alpha=alpha,
                                           cmap='viridis')}}
    ds = hv.Dataset(category_df, key_dimensions, value_dimensions)
    point_g = ds.to(hv.Scatter, [x, y]).opts(opts_dict)

    if return_text:
        text_annos = [hv.Text(row[x], row[y], row[name_col]).options(
            text_font_size=f'{log2(row[size_col])*font_size}pt')
            for i, row in category_df.iterrows()]
        text_g = hv.Overlay(text_annos).options(tools=['hover'])
        return point_g * text_g
    else:
        return point_g
