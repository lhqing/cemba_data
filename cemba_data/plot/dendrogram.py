import matplotlib.patches as mpatches
import matplotlib.path as mpath
from scipy.cluster.hierarchy import dendrogram

Path = mpath.Path


def straight_branch(ax, a, b, plot_kws):
    """Draw link line between a and b"""
    a_x, ay = a
    bx, by = b
    branch_x = [a_x, bx, bx]
    branch_y = [ay, ay, by]
    if plot_kws is None:
        plot_kws = {}
    return ax.plot(branch_x, branch_y, **plot_kws)


def curve_branch(ax, a, b, patch_kws):
    a_x, ay = a
    bx, by = b
    mid_x = (bx + a_x) / 2
    mid_y = (by + ay) / 2
    branch_x = [a_x, a_x, mid_x, bx, bx]
    branch_y = [ay, mid_y, mid_y, mid_y, by]

    if patch_kws is None:
        patch_kws = {}

    pp1 = mpatches.PathPatch(Path(
        [(x, y) for x, y in zip(branch_x, branch_y)],
        [Path.MOVETO, Path.CURVE3, Path.CURVE3, Path.CURVE3, Path.CURVE3]),
        fc="none", **patch_kws)
    return ax.add_patch(pp1)


def plot_dendrogram(linkage_df,
                    labels_list,
                    dendro_kws=None,
                    ax=None,
                    branch_type='straight', plot_kws=None):
    if plot_kws is None:
        plot_kws = {}

    _dendro_kws = dict(no_plot=True)
    if dendro_kws is not None:
        _dendro_kws.update(dendro_kws)
    # all we need is the leaves order from dendrogram,
    # bellow we recalculate the node position to match the node id,
    # so we can control each node
    dendro = dendrogram(linkage_df, labels=labels_list, **_dendro_kws)

    node_pos = {}  # all node including singleton and non-singleton
    for leaf_x, leaf in enumerate(dendro['leaves']):
        node_pos[int(leaf)] = (leaf_x, 0)

    direct_link_map = {}  # node linkage, keys only contain non-singleton
    for i, (left, right, height, _) in linkage_df.iterrows():
        node_id = int(i + linkage_df.shape[0] + 1)
        left = int(left)
        right = int(right)
        node_x = (node_pos[left][0] + node_pos[right][0]) / 2
        node_pos[node_id] = [node_x, height]
        direct_link_map[node_id] = [int(left), int(right)]

    if ax is not None:
        for node_id, (node_x, node_y) in node_pos.items():
            # plot node
            ax.text(node_x, node_y, node_id, fontsize=4, ha='center', va='center')

            # plot branch
            # only non-singleton node has branch:
            if node_id in direct_link_map:
                # get child
                left_child, right_child = direct_link_map[node_id]

                # plot left branch
                if branch_type == 'straight':
                    # left_branch_lines
                    _ = straight_branch(ax, (node_x, node_y),
                                        node_pos[left_child], plot_kws)
                else:
                    # left_branch_lines
                    _ = curve_branch(ax, (node_x, node_y),
                                     node_pos[left_child], plot_kws)

                # plot right branch
                if branch_type == 'straight':
                    # right_branch_lines
                    _ = straight_branch(ax, (node_x, node_y),
                                        node_pos[right_child], plot_kws)
                else:
                    # right_branch_lines
                    _ = curve_branch(ax, (node_x, node_y),
                                     node_pos[right_child], plot_kws)
    return node_pos, direct_link_map, dendro
