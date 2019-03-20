import numpy as np
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle


class Sankey:
    def __init__(self, ax, data, col_order=None,
                 left_pad=0.05, right_pad=0.05,
                 label_width=0.05, gap_height=0.02,
                 column_height_scale=1., label_orders=None, strip_min_n=1):
        self.ax = ax
        if col_order is not None:
            for col in col_order:
                if col not in data.columns:
                    raise KeyError(f'column {col} not in data, '
                                   f'data columns: {data.columns}, '
                                   f'col_order: {col_order}')
            self.data = data[col_order].copy()
        else:
            self.data = data.copy()
        if label_orders is None:
            label_orders = {}
        self.label_orders = {col: label_orders.pop(col, None) for col in self.data.columns}

        self.label_width = label_width
        self.gap_height = gap_height
        if isinstance(column_height_scale, float):
            self.column_height_scale = {col: column_height_scale for col in self.data.columns}
        else:
            self.column_height_scale = column_height_scale
        self.ncols = self.data.shape[1]
        if self.ncols < 2:
            raise ValueError(f'At least 2 columns needed for making Sankey plot, '
                             f'got {self.ncols} only (after apply col_order if provided).')
        self.col_centers = np.linspace(left_pad, 1 - right_pad, self.ncols)

        # add columns
        self.columns = []
        self._add_cols()

        # add strips
        self.strip_min_n = strip_min_n
        self.strips = []
        self._add_strips()

    def _add_cols(self):
        for col_id, (col_name, col) in enumerate(self.data.iteritems()):
            label_n_dict = col.value_counts().to_dict()
            _label_order = self.label_orders[col_name]
            _column_height_scale = self.column_height_scale[col_name]
            self.columns.append(Column(col_id=col_id,
                                       name=col_name,
                                       col_center=self.col_centers[col_id],
                                       label_width=self.label_width,
                                       label_n_dict=label_n_dict,
                                       label_order=_label_order,
                                       gap_height=self.gap_height,
                                       column_height_scale=_column_height_scale))
        return

    def _add_strips(self):
        # add strips
        for i in range(self.ncols - 1):
            left_col = self.columns[i]
            right_col = self.columns[i + 1]

            pair_data = self.data.iloc[:, [i, i + 1]]
            pair_count = pair_data.groupby(pair_data.columns.tolist()) \
                .apply(lambda sub: sub.shape[0])

            for left_label in left_col.labels:
                for right_label in right_col.labels:
                    try:
                        n = pair_count[(left_label.name, right_label.name)]
                    except KeyError:
                        continue
                    if n >= self.strip_min_n:
                        self.strips.append(Strip(left_label=left_label,
                                                 right_label=right_label,
                                                 n=n))
        return

    def draw(self, curve_value=50, fill_kws=None,
             col_label_facecolor='steelblue', strip_color='left'):
        # draw rectangles
        rectangles = []
        facecolors = []
        for col in self.columns:
            for label in col.labels:
                if isinstance(col_label_facecolor, dict):
                    color = col_label_facecolor[(col.name, label.name)]
                else:
                    color = col_label_facecolor
                facecolors.append(color)
                rectangles.append(label.get_rectangle(facecolor=color))
        # Create patch collection with specified colour/alpha
        pc = PatchCollection(rectangles, facecolors=facecolors)
        self.ax.add_collection(pc)

        # draw strips
        for strip in self.strips:
            if strip_color == 'left':
                facecolor = strip.left_label.color
            elif strip_color == 'right':
                facecolor = strip.right_label.color
            else:
                # assume this is a color that can be understood in ax.fill_between
                facecolor = strip_color
            _fill_kws = {'facecolor': facecolor}
            if fill_kws is not None:
                _fill_kws.update(fill_kws)
            strip.plot_strip(ax=self.ax, curve_value=curve_value, **_fill_kws)
        return self.ax


class Column:
    def __init__(self, col_id, name, col_center, label_width, label_n_dict,
                 label_order=None, gap_height=0.02, column_height_scale=1.):
        if column_height_scale <= 0:
            raise ValueError(f'column_height_scale must > 0, got {column_height_scale}')
        self.col_id = col_id
        self.name = name
        self.col_center = col_center
        self.n_label = len(label_n_dict)

        cur_bottom = 0.5 - column_height_scale / 2  # col is centered in 0.5
        # each gap is same
        total_gap_height = (self.n_label - 1) * gap_height
        unit_height = column_height_scale * (1 - total_gap_height) / sum(label_n_dict.values())
        gap_height = gap_height * column_height_scale

        self.labels = []
        if label_order is None:
            label_order = label_n_dict.keys()

        label_left = col_center - label_width / 2
        for label in label_order:
            label_n = label_n_dict[label]
            self.labels.append(Label(name=label, left=label_left, n=label_n,
                                     width=label_width, unit_height=unit_height,
                                     bottom=cur_bottom))
            cur_bottom = cur_bottom + label_n * unit_height + gap_height


class Label:
    def __init__(self, name, left, n, width, unit_height, bottom):
        self.name = name
        self.left = left
        self.bottom = bottom
        self.n = n
        self.height = unit_height * n
        self.unit_height = unit_height
        self.width = width

        self.left_strip_n_cum = 0
        self.right_strip_n_cum = 0
        self.color = 'steelblue'

    def get_rectangle(self, **rect_kws):
        if 'facecolor' in rect_kws:
            # update label color
            self.color = rect_kws['facecolor']
        # note: the facecolor will not work here...
        # the real color parameter will be collected in PatchCollection in the Sankey.draw method
        rect = Rectangle((self.left, self.bottom),
                         self.width, self.height, **rect_kws)
        return rect

    def link_left_strip(self, strip_n):
        self.left_strip_n_cum += strip_n

    def link_right_strip(self, strip_n):
        self.right_strip_n_cum += strip_n


class Strip:
    def __init__(self, left_label, right_label, n):
        self.left_label: Label = left_label
        self.right_label: Label = right_label
        self.n = n

        self.left = left_label.left + left_label.width  # right side of left label
        self.right = right_label.left  # left side of right label

        # calculate the real bottom and height based on label
        # on the right of left label
        self.left_bottom = left_label.bottom + left_label.unit_height * left_label.right_strip_n_cum
        self.left_height = left_label.unit_height * n
        # on the left of right label
        self.right_bottom = right_label.bottom + right_label.unit_height * right_label.left_strip_n_cum
        self.right_height = right_label.unit_height * n

        # call link_*_strip will change Label.*strip_cum,
        # which is used to determine real strip position for next strip
        left_label.link_right_strip(self.n)
        right_label.link_left_strip(self.n)
        return

    def plot_strip(self, ax, curve_value=50, **fill_kws):
        """
        This function taken from here:
        https://github.com/anazalea/pySankey
        """
        if curve_value < 40:
            curve_value = 40
            print(f'Curve value >= 40, got {curve_value}, set to 40 instead.')
        ys_d = np.array(curve_value * [self.left_bottom] + \
                        curve_value * [self.right_bottom])
        ys_d = np.convolve(ys_d, 0.05 * np.ones(20), mode='valid')
        ys_d = np.convolve(ys_d, 0.05 * np.ones(20), mode='valid')

        ys_u = np.array(curve_value * [self.left_bottom + self.left_height] + \
                        curve_value * [self.right_bottom + self.right_height])
        ys_u = np.convolve(ys_u, 0.05 * np.ones(20), mode='valid')
        ys_u = np.convolve(ys_u, 0.05 * np.ones(20), mode='valid')

        _fill_kws = {'alpha': 0.65}
        if fill_kws is not None:
            _fill_kws.update(fill_kws)

        ax.fill_between(
            np.linspace(self.left, self.right, len(ys_d)),
            ys_d, ys_u, **_fill_kws)
        return ax


""" test

farms = np.random.randint(1, 5, size=300)
stations = np.random.randint(1, 5, size=300)
customers = np.random.randint(1, 5, size=300)
gender = np.random.randint(1, 3, size=300)

sankey_data = pd.DataFrame({'farms': farms, 
                            'stations': stations,
                            'customers': customers,
                            'gender': gender})

color_dict = {('farms', 1): 'red', 
             ('farms', 2): 'green', 
             ('farms', 3): 'blue', 
             ('farms', 4): 'orange', 
             ('customers', 1): 'red', 
             ('customers', 2): 'green', 
             ('customers', 3): 'blue', 
             ('customers', 4): 'orange',  
             ('stations', 1): 'red', 
             ('stations', 2): 'green', 
             ('stations', 3): 'blue', 
             ('stations', 4): 'orange', 
             ('gender', 1): 'blue', 
             ('gender', 2): 'orange'}

fig, ax = plt.subplots(figsize=(8, 8))
sankey = Sankey(ax=ax, data=sankey_data, gap_height=0.005, strip_min_n=0)
sankey.draw(curve_value=50, col_label_facecolor=color_dict, strip_color='left')
"""
