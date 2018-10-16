import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects
import scanpy.api as sc
from ..data.hdf5 import Study
from .zoo import *


def get_feature_dispersion(ann, min_mean=0.2, max_mean=1.25, min_disp=0.3):
    filter_result = sc.pp.filter_genes_dispersion(ann.X, min_mean=min_mean,
                                                  max_mean=max_mean, min_disp=min_disp)
    fig, axes = plt.subplots(nrows=1, ncols=2)
    fig.set_size_inches(10, 3)
    axes[0].set_title('Feature Raw Dispersions')
    axes[1].set_title('Feature Normalized Dispersions')
    axes[0].set_ylabel('log Raw Dispersion')
    axes[1].set_ylabel('log Norm Dispersion')
    axes[0].set_xlabel('log1p Region mean mC%')
    axes[1].set_xlabel('log1p Region mean mC%')
    axes[0].scatter(x=filter_result.means,
                    y=filter_result.dispersions,
                    c=['black' if i else 'grey' for i in filter_result.gene_subset], s=1)
    axes[1].scatter(x=filter_result.means,
                    y=filter_result.dispersions_norm,
                    c=['black' if i else 'grey' for i in filter_result.gene_subset], s=1)
    fig.tight_layout()
    print(f'{filter_result.gene_subset.sum()} features been selected.')
    return filter_result


def ann_to_study(ann):
    col_dict = ann.var.to_dict('series')
    col_dict['col_names'] = ann.var.index
    row_dict = ann.obs.to_dict('series')
    row_dict['row_names'] = ann.obs.index
    uns_dict = ann.uns
    study = Study(ann.X, col_dict=col_dict, row_dict=row_dict, uns_dict=uns_dict, study_name=None)
    return study


def plot_feature_set(ann):
    features = ['region', 'replicate', 'louvain', 'Mapped reads', 'overall_mch']
    if 'manual_anno' in ann.obs_keys():
        features.append('manual_anno')
    colors = plt.get_cmap('Paired').colors
    random.shuffle(list(colors))

    fig, axes = plt.subplots(nrows=2, ncols=len(features))
    base_l = ['UMAP'] * len(features) + ['tSNE'] * len(features)

    x_l = [ann.obsm[f'X_{base.lower()}'][:, 0] for base in base_l]
    y_l = [ann.obsm[f'X_{base.lower()}'][:, 1] for base in base_l]

    for x, y, base, ax, feature in zip(x_l, y_l, base_l, np.ravel(axes), features * 2):
        feature_value = ann.obs[feature]
        if feature_value.dtype.name == 'object':
            feature_value = ann.obs[feature].astype('category')
        if feature_value.dtype.name != 'category':
            if feature == 'overall_mch':
                im = ax.scatter(x=x, y=y, alpha=0.5, s=10, c=feature_value.tolist(), cmap='viridis', vmin=0, vmax=0.03)
            elif feature == 'Mapped reads':
                im = ax.scatter(x=x, y=y, alpha=0.5, s=10, c=feature_value.tolist(), cmap='viridis', vmin=5e5, vmax=3e6)
            else:
                print(feature)
                raise ValueError('Unknown feature name', feature)
            fig.colorbar(im, ax=ax, orientation='vertical', format='%.0e')
        else:
            if feature in ['louvain', 'manual_anno']:
                legend_loc = 'on data'
            else:
                legend_loc = 'right margin'
            if base == 'tSNE':
                sc.pl.tsne(ann, color=[feature], legend_loc=legend_loc, ax=ax, show=False, alpha=0.7)
            else:
                sc.pl.umap(ann, color=[feature], legend_loc=legend_loc, ax=ax, show=False, alpha=0.7)
        ax.set_title(feature)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xlabel(f'{base} 1')
        ax.set_ylabel(f'{base} 2')
    fig.set_size_inches(5 * len(features), 8)
    fig.tight_layout()
    return fig, axes


def plot_marker_violin(ann, groups, n_gene=6):
    nplot = len(groups)
    fig, axes = plt.subplots(nrows=1, ncols=nplot)
    if nplot == 1:
        axes = [axes]
    cat_set = set()
    for ax, group in zip(axes, groups):
        sc.pl.rank_genes_groups_violin(ann, groups=group, n_genes=n_gene, ax=ax, show=False)
        gene_list = ann.uns['rank_genes_groups']['names'][group][:n_gene]
        cat_list = [cat.get_gene_name(i) for i in gene_list]
        [cat_set.add(i) for i in cat_list]
        ax.set_xticklabels(cat_list)
        ax.set_ylabel('Normalized mCH%')
        # ax.set_ylim(0, 5)
        ax.set_title(f'{group} vs. rest, top {n_gene} marker')
        fig.set_size_inches(5 * nplot, 4)
    return fig, axes, list(cat_set)


def plot_genes(gene_list, ann, s=1, vmin=0.01, vmax=0.99, species='mouse'):
    marker_list = [cat.get_gene_name(i, species=species) for i in gene_list]
    plot_number = len(marker_list)

    fig, axes = plt.subplots(2, ncols=plot_number)
    base_l = ['UMAP'] * len(gene_list) + ['tSNE'] * len(gene_list)

    x_l = [ann.obsm[f'X_{base.lower()}'][:, 0] for base in base_l]
    y_l = [ann.obsm[f'X_{base.lower()}'][:, 1] for base in base_l]

    for x, y, base, ax, marker, title in zip(x_l, y_l, base_l, np.ravel(axes), marker_list * 2, gene_list * 2):
        ax.scatter(x=x, y=y, s=s, c=list(np.ravel(ann[:, marker].X)), cmap='viridis', vmin=vmin, vmax=vmax)
        ax.set_title(title)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xlabel(f'{base} 1')
        ax.set_ylabel(f'{base} 2')
    fig.set_size_inches(3 * len(marker_list), 6)
    fig.tight_layout()
    return fig, axes


def anno_gene_ann(gene_ann, compute_ann, cell_meta):
    for k, v in compute_ann.uns.items():
        gene_ann.uns[k] = v

    gene_ann = gene_ann[compute_ann.obs_names, :]
    for var_name, var_data in compute_ann.var.iteritems():
        gene_ann.var[var_name] = var_data

    for obs_name, obs_data in compute_ann.obs.iteritems():
        gene_ann.obs[obs_name] = obs_data

    for k in compute_ann.obsm.keys():
        gene_ann.obsm[k] = compute_ann.obsm[k]

    for col, data in cell_meta.iteritems():
        if col in ['region', 'replicate', 'Total reads', 'Mapped reads', 'mCG/CG', 'mCH/CH', '% Genome covered']:
            gene_ann.obs[col] = data
    return gene_ann


def category_color_transfer(categories, cmap=None, shuffle=False):
    category_set = set(categories)
    if cmap is None:
        if len(category_set) < 15:
            cmap = 'tab10'
        else:
            cmap = 'tab20'
    colors = list(plt.get_cmap(cmap).colors)
    if shuffle:
        random.shuffle(colors)
    index_dict = {cat: i % len(colors) for i, cat in enumerate(category_set)}
    return [colors[index_dict[i]] for i in categories]


def plot_cat_scatter(ax, coord_df, colors, coord_name, title, alpha=0.8,
                     xlim=None, ylim=None, text_size='small', linewidth=1,
                     auto_border_quantile=0.99, auto_expand=0.1, cmap=None, s=10):
    ax.scatter(x=coord_df.iloc[:, 0], y=coord_df.iloc[:, 1],
               c=category_color_transfer(colors, cmap=cmap), s=s, alpha=alpha)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel(coord_name + ' 1')
    ax.set_ylabel(coord_name + ' 2')
    ax.set_title(title)

    if xlim is not None:
        ax.set_xlim(*xlim)
    else:
        amin = coord_df.iloc[:, 0].quantile(1 - auto_border_quantile)
        amax = coord_df.iloc[:, 0].quantile(auto_border_quantile)
        ax.set_xlim(amin - (amax - amin) * auto_expand, amax + (amax - amin) * auto_expand)
    if ylim is not None:
        ax.set_ylim(*ylim)
    else:
        amin = coord_df.iloc[:, 1].quantile(1 - auto_border_quantile)
        amax = coord_df.iloc[:, 1].quantile(auto_border_quantile)
        ax.set_ylim(amin - (amax - amin) * auto_expand, amax + (amax - amin) * auto_expand)
    plot_cluster_text(coord_df, colors, ax, fontsize=text_size, linewidth=linewidth)
    return ax


def plot_cluster_text(coord_df, cluster, ax,
                      fontsize='small', edgecolor='white', linewidth=2.5, facecolor='black',
                      fix_outsider=True):
    cur_xlim = ax.get_xlim()
    cur_ylim = ax.get_ylim()
    for c in cluster.groupby(cluster):
        c_id = str(c[0])
        cells = c[1].index
        cell_cord = coord_df.loc[cells].median().tolist()
        if fix_outsider:
            cell_cord[0] = min(cell_cord[0], cur_xlim[1])
            cell_cord[0] = max(cell_cord[0], cur_xlim[0])
            cell_cord[1] = min(cell_cord[1], cur_ylim[1])
            cell_cord[1] = max(cell_cord[1], cur_ylim[0])

        t = ax.text(*cell_cord, c_id, fontsize=fontsize, fontweight='semibold')
        t.set_path_effects([path_effects.PathPatchEffect(edgecolor=edgecolor,
                                                         linewidth=linewidth,
                                                         facecolor=facecolor)])
    return ax
