import pandas as pd
import numpy as np
from collections import OrderedDict


def sunbrust(pie_data, ax, inner_radius=0.25, outer_radius=1,
             anno_col=None, text_anno='text', anno_layer_size=0.05,
             col_color_dict=None, startangle=0, anno_ang_min=5,
             anno_border=1.2, text_expend=1.2,
             uniform_section=False):
    if uniform_section:
        dedup_groups = pie_data.reset_index().set_index(pie_data.columns.tolist()).index.drop_duplicates()
        dedup_groups = pd.DataFrame(dedup_groups.tolist(), columns=pie_data.columns)
        pie_data = dedup_groups

    _col_color_dict = {col: None for col in pie_data.columns}
    if col_color_dict is not None:
        _col_color_dict.update(col_color_dict)

    ncols = pie_data.columns.size
    if anno_col is None:
        anno_layer_size = 0
    outer_radius = outer_radius - anno_layer_size
    layer_size = (outer_radius - inner_radius) / ncols

    previous_order = pd.Series([])
    anno_wedges = []
    anno_names = []
    for col, col_name in enumerate(pie_data.columns):
        cur_radius = inner_radius + (col + 1) * layer_size
        col_pie_data = pie_data[col_name].value_counts()

        # manage order
        if col == 0:
            _ordered_data = col_pie_data
            previous_order = _ordered_data
        else:
            records = []
            for section in previous_order.index:
                section_subs = pie_data[pie_data.iloc[:, col - 1] == section][col_name].unique()
                records.append(col_pie_data.reindex(pd.Index(section_subs)).sort_values(ascending=False))
            _ordered_data = pd.concat(records)
            previous_order = _ordered_data

        pie_color = _col_color_dict[col_name]
        if isinstance(pie_color, dict):
            pie_color = [pie_color[i] for i in _ordered_data.index]
        ax.pie(_ordered_data, radius=cur_radius,
               colors=pie_color, startangle=startangle,
               wedgeprops=dict(width=layer_size, edgecolor='w'))

        if anno_col == col:
            wedges, texts = ax.pie(_ordered_data, radius=outer_radius + anno_layer_size,
                                   colors=pie_color, startangle=startangle,
                                   wedgeprops=dict(width=anno_layer_size, edgecolor='w'))
            if text_anno:
                anno_wedges = wedges
                anno_names = _ordered_data.index.tolist()

    # wedges annotation
    bbox_props = dict(boxstyle="round,pad=0.2", fc="#FFFFFF88", ec="#00000022", lw=0.72)
    kw = dict(xycoords='data', textcoords='data', arrowprops=dict(arrowstyle="-", color='#00000055'),
              bbox=bbox_props, zorder=0, va="center")

    # separate all y
    y_niche = np.arange(-anno_border, anno_border, anno_border / 10)
    allow_y_niche = OrderedDict(enumerate(y_niche))

    # annotate wedges
    for i, p in enumerate(anno_wedges):
        delta_ang = p.theta2 - p.theta1
        if delta_ang < anno_ang_min:
            continue
        ang = (p.theta2 - p.theta1) / 2. + p.theta1
        y = np.sin(np.deg2rad(ang))
        x = np.cos(np.deg2rad(ang))

        if text_anno == 'anno_box':
            horizontalalignment = {-1: "right", 1: "left"}[int(np.sign(x))]
            connectionstyle = f"angle,angleA=0,angleB={ang},rad=5"
            kw["arrowprops"].update({"connectionstyle": connectionstyle})

            suitable_niche = np.abs(np.array(list(allow_y_niche.values())) - anno_border * y).argmin()
            suitable_niche = list(allow_y_niche.keys())[suitable_niche]
            y_niche = allow_y_niche.pop(suitable_niche)
            ax.annotate(anno_names[i], xy=(x, y), xytext=(anno_border * np.sign(x), y_niche),
                        horizontalalignment=horizontalalignment, **kw)
        elif text_anno == 'text':
            if x < 0:
                ang += 180
            if ang > 180:
                ang = ang - 360
            ax.text(x * text_expend, y * text_expend, anno_names[i],
                    fontdict=None, withdash=False,
                    rotation=ang, va='center', ha='center')
        elif text_anno is None:
            pass
        else:
            raise ValueError(f'text_anno can only be "text", "anno_box" or None, got {text_anno}')
