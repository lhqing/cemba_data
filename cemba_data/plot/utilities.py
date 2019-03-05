import numpy as np
from matplotlib.colorbar import ColorbarBase


def tight_hue_range(hue_data, portion):
    """Automatic select a SMALLEST data range that covers [portion] of the data"""
    hue_quantiles = hue_data.quantile(q=np.arange(0, 1, 0.01))
    min_window_right = hue_quantiles.rolling(window=int(portion * 100)) \
        .apply(lambda i: i.max() - i.min(), raw=True) \
        .idxmin()
    min_window_left = max(0, min_window_right - portion)
    vmin, vmax = tuple(hue_data.quantile(q=[min_window_left,
                                            min_window_right]))
    if np.isnan(vmin):
        vmin = hue_data.min()
    else:
        vmin = max(hue_data.min(), vmin)
    if np.isnan(vmax):
        vmax = hue_data.max()
    else:
        vmax = min(hue_data.max(), vmax)
    return vmin, vmax


def plot_colorbar(cax, cmap, cnorm, hue_norm, label_kws):
    _label_kws = {'fontsize': 3, 'label': ''}
    _label_kws.update(label_kws)

    colorbar = ColorbarBase(cax, cmap=cmap, norm=cnorm,
                            orientation='vertical', extend='both')
    colorbar.set_label(**_label_kws)
    colorbar_ticks = [hue_norm[0], sum(hue_norm) / 2, hue_norm[1]]
    # TODO automatic ticklabel format, auto sci-format, float trim etc
    colorbar_ticklabels = list(map(lambda i: f'{i:.2f}', colorbar_ticks))
    colorbar.set_ticks(colorbar_ticks)
    colorbar.set_ticklabels(colorbar_ticklabels)
    colorbar.outline.set_linewidth(0.5)
    colorbar.ax.tick_params(size=_label_kws['fontsize'], labelsize=_label_kws['fontsize'],
                            pad=1, width=0.5)
    return cax
