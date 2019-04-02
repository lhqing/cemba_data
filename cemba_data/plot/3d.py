from .color import level_one_palette
from matplotlib import animation


def clear_axes3d_background(ax):
    # Get rid of the panes
    ax.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.w_yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.w_zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))

    # Get rid of the spines
    ax.w_xaxis.line.set_color((1.0, 1.0, 1.0, 0.0))
    ax.w_yaxis.line.set_color((1.0, 1.0, 1.0, 0.0))
    ax.w_zaxis.line.set_color((1.0, 1.0, 1.0, 0.0))

    # Get rid of the ticks
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_zticks([])
    return ax


def plot_cat_3d_scatter(ax_3d, x, y, z, hue, data, hue_palette=None, scatter_kws=None,
                        remove_coords=False):
    data = data.copy()

    if isinstance(hue, str):
        hue = hue
    else:
        data['hue'] = hue
        hue = 'hue'
    if hue_palette is None:
        hue_palette = level_one_palette(hue)

    xs, ys, zs = data[x], data[y], data[z]
    colors = data[hue].map(hue_palette).tolist()
    _scatter_kws = dict(alpha=0.8, linewidth=0, s=3)
    if scatter_kws is not None:
        _scatter_kws.update(scatter_kws)
    ax_3d.scatter(xs, ys, zs, c=colors, **_scatter_kws)
    if remove_coords:
        ax_3d = clear_axes3d_background(ax_3d)
    return ax_3d


def make_animation(fig, ax, out_path, frames=360, interval=10):
    def update(i):
        ax.view_init(10, i)
    ani = animation.FuncAnimation(
        fig, func=update, frames=frames, interval=interval, blit=False)
    ani.save(out_path, writer='imagemagic')
    return ani
