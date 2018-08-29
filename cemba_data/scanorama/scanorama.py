from annoy import AnnoyIndex
from intervaltree import IntervalTree
from itertools import cycle, islice
import numpy as np
import operator
import random
from sklearn.manifold import TSNE
from sklearn.metrics.pairwise import rbf_kernel
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import normalize
import sys

from .process import merge_datasets
from .t_sne_approx import TSNEApprox
from .utils import plt, dispersion, reduce_dimensionality
from .utils import visualize_cluster, visualize_expr

np.random.seed(0)
random.seed(0)

# Default parameters.
ALPHA = 0.10
APPROX = True
DIMRED = 100
HVG = 10000
KNN = 20
N_ITER = 500
PERPLEXITY = 1200
SIGMA = 150
VERBOSE = 2


# Do batch correction on the data.
def correct(datasets_full, genes_list, hvg=HVG, verbose=VERBOSE,
            sigma=SIGMA, ds_names=None):
    datasets, genes = merge_datasets(datasets_full, genes_list)
    datasets_dimred, genes = process_data(datasets, genes, hvg=hvg)

    assemble(
        datasets_dimred,  # Assemble in low dimensional space.
        expr_datasets=datasets,  # Modified in place.
        verbose=verbose, knn=KNN, sigma=sigma, approx=APPROX,
        ds_names=ds_names
    )

    return datasets, genes


# Randomized SVD.
def dimensionality_reduce(datasets, dimred=DIMRED):
    X = np.concatenate(datasets)
    X = reduce_dimensionality(X, dim_red_k=dimred)
    datasets_dimred = []
    base = 0
    for ds in datasets:
        datasets_dimred.append(X[base:(base + ds.shape[0]), :])
        base += ds.shape[0]
    return datasets_dimred


# Normalize and reduce dimensionality.
def process_data(datasets, genes, hvg=HVG, dimred=DIMRED):
    # Only keep highly variable genes
    if hvg > 0:
        X = np.concatenate(datasets)
        disp = dispersion(X)
        top_genes = set(genes[
                            list(reversed(np.argsort(disp)))[:HVG]
                        ])
        for i in range(len(datasets)):
            gene_idx = [idx for idx, g_i in enumerate(genes)
                        if g_i in top_genes]
            datasets[i] = datasets[i][:, gene_idx]
        genes = np.array(sorted(top_genes))

    # Normalize.
    datasets = [normalize(ds, axis=1) for ds in datasets]

    # Compute compressed embedding.
    if dimred > 0:
        datasets_dimred = dimensionality_reduce(datasets)
        return datasets_dimred, genes

    return datasets, genes


# Exact nearest neighbors search.
def nn(ds1, ds2, knn=KNN, metric_p=2):
    # Find nearest neighbors of first dataset.
    nn_ = NearestNeighbors(knn, p=metric_p)
    nn_.fit(ds2)
    ind = nn_.kneighbors(ds1, return_distance=False)

    match = set()
    for a, b in zip(range(ds1.shape[0]), ind):
        for b_i in b:
            match.add((a, b_i))

    return match


# Approximate nearest neighbors using locality sensitive hashing.
def nn_approx(ds1, ds2, knn=KNN, metric='angular', n_trees=10):
    # Build index.
    a = AnnoyIndex(ds2.shape[1], metric=metric)
    for i in range(ds2.shape[0]):
        a.add_item(i, ds2[i, :])
    a.build(n_trees)

    # Search index.
    ind = []
    for i in range(ds1.shape[0]):
        ind.append(a.get_nns_by_vector(ds1[i, :], knn, search_k=-1))
    ind = np.array(ind)

    # Match.
    match = set()
    for a, b in zip(range(ds1.shape[0]), ind):
        for b_i in b:
            match.add((a, b_i))

    return match


# Find mutual nearest neighbors.
def mnn(ds1, ds2, knn=KNN, approx=APPROX):
    # Find nearest neighbors in first direction.
    if approx:
        match1 = nn_approx(ds1, ds2, knn=knn)
    else:
        match1 = nn(ds1, ds2, knn=knn)

    # Find nearest neighbors in second direction.
    if approx:
        match2 = nn_approx(ds2, ds1, knn=knn)
    else:
        match2 = nn(ds2, ds1, knn=knn)

    # Compute mutual nearest neighbors.
    mutual = match1 & set([(b, a) for a, b in match2])

    return mutual


# Visualize alignment between two datasets.
def plot_mapping(curr_ds, curr_ref, ds_ind, ref_ind):
    tsne = TSNE(n_iter=400, verbose=VERBOSE, random_state=69)

    tsne.fit(curr_ds)
    plt.figure()
    coords_ds = tsne.embedding_[:, :]
    coords_ds[:, 1] += 100
    plt.scatter(coords_ds[:, 0], coords_ds[:, 1])

    tsne.fit(curr_ref)
    coords_ref = tsne.embedding_[:, :]
    plt.scatter(coords_ref[:, 0], coords_ref[:, 1])

    x_list, y_list = [], []
    for r_i, c_i in zip(ds_ind, ref_ind):
        x_list.append(coords_ds[r_i, 0])
        x_list.append(coords_ref[c_i, 0])
        x_list.append(None)
        y_list.append(coords_ds[r_i, 1])
        y_list.append(coords_ref[c_i, 1])
        y_list.append(None)
    plt.plot(x_list, y_list, 'b-', alpha=0.3)
    plt.show()


# Populate a table (in place) that stores mutual nearest neighbors
# between datasets.
def fill_table(table, i, curr_ds, datasets, base_ds=0,
               knn=KNN, approx=APPROX, metric='angular'):
    curr_ref = np.concatenate(datasets)
    if approx:
        match = nn_approx(curr_ds, curr_ref, knn=knn, metric=metric)
    else:
        match = nn(curr_ds, curr_ref, knn=knn, metric_p=1)

    # Build interval tree.
    itree_ds_idx = IntervalTree()
    itree_pos_base = IntervalTree()
    pos = 0
    for j in range(len(datasets)):
        n_cells = datasets[j].shape[0]
        itree_ds_idx[pos:(pos + n_cells)] = base_ds + j
        itree_pos_base[pos:(pos + n_cells)] = pos
        pos += n_cells

    # Store all mutual nearest neighbors between datasets.
    for d, r in match:
        interval = itree_ds_idx[r]
        assert (len(interval) == 1)
        j = interval.pop().data
        interval = itree_pos_base[r]
        assert (len(interval) == 1)
        base = interval.pop().data
        if not (i, j) in table:
            table[(i, j)] = set()
        table[(i, j)].add((d, r - base))
        assert (r - base >= 0)


# Fill table of alignment scores.
def find_alignments_table(datasets, knn=KNN, approx=APPROX,
                          verbose=VERBOSE, prenormalized=False, metric='angular'):
    if not prenormalized:
        datasets = [normalize(ds, axis=1) for ds in datasets]

    table = {}
    for i in range(len(datasets)):
        if len(datasets[:i]) > 0:
            fill_table(table, i, datasets[i], datasets[:i], knn=knn,
                       approx=approx, metric=metric)
        if len(datasets[i + 1:]) > 0:
            fill_table(table, i, datasets[i], datasets[i + 1:],
                       knn=knn, base_ds=i + 1, approx=approx, metric=metric)

    # Count all mutual nearest neighbors between datasets.
    matches = {}
    table1 = {}
    table_print = np.zeros((len(datasets), len(datasets)))
    for i in range(len(datasets)):
        for j in range(len(datasets)):
            if i >= j:
                continue
            if not (i, j) in table or not (j, i) in table:
                continue
            match_ij = table[(i, j)]
            match_ji = set([(b, a) for a, b in table[(j, i)]])
            matches[(i, j)] = match_ij & match_ji
            table1[(i, j)] = (max(
                float(len(set([idx for idx, _ in matches[(i, j)]]))) /
                datasets[i].shape[0],
                float(len(set([idx for _, idx in matches[(i, j)]]))) /
                datasets[j].shape[0]
            ))
            table_print[i, j] += table1[(i, j)]
    if verbose > 1:
        print(table_print)

    return table1, table_print, matches


# Find the matching pairs of cells between datasets.
def find_alignments(datasets, knn=KNN, approx=APPROX, verbose=VERBOSE,
                    prenormalized=False, metric='angular'):
    table1, _, matches = find_alignments_table(
        datasets, knn=knn, approx=approx, verbose=verbose,
        prenormalized=prenormalized, metric=metric
    )

    alignments = [(i, j) for (i, j), val in reversed(
        sorted(table1.items(), key=operator.itemgetter(1))
    ) if val > ALPHA]

    return alignments, matches


# Find connections between datasets to identify panoramas.
def connect(datasets, knn=KNN, approx=APPROX, verbose=VERBOSE, metric='angular'):
    # Find alignments.
    alignments, _ = find_alignments(
        datasets, knn=knn, approx=approx, verbose=verbose, metric=metric
    )
    if verbose:
        print(alignments)

    panoramas = []
    connected = set()
    for i, j in alignments:
        # See if datasets are involved in any current panoramas.
        panoramas_i = [panoramas[p] for p in range(len(panoramas))
                       if i in panoramas[p]]
        assert (len(panoramas_i) <= 1)
        panoramas_j = [panoramas[p] for p in range(len(panoramas))
                       if j in panoramas[p]]
        assert (len(panoramas_j) <= 1)

        if len(panoramas_i) == 0 and len(panoramas_j) == 0:
            panoramas.append([i])
            panoramas_i = [panoramas[-1]]

        if len(panoramas_i) == 0:
            panoramas_j[0].append(i)
        elif len(panoramas_j) == 0:
            panoramas_i[0].append(j)
        elif panoramas_i[0] != panoramas_j[0]:
            panoramas_i[0] += panoramas_j[0]
            panoramas.remove(panoramas_j[0])

        connected.add(i)
        connected.add(j)

    for i in range(len(datasets)):
        if not i in connected:
            panoramas.append([i])

    return panoramas


# Compute nonlinear translation vectors between dataset
# and a reference.
def transform(curr_ds, curr_ref, ds_ind, ref_ind, sigma):
    # Compute the matching.
    match_ds = curr_ds[ds_ind, :]
    match_ref = curr_ref[ref_ind, :]
    bias = match_ref - match_ds
    weights = rbf_kernel(curr_ds, match_ds, gamma=0.5 * sigma)
    avg_bias = np.dot(weights, bias) / \
               np.tile(np.sum(weights, axis=1),
                       (curr_ds.shape[1], 1)).T
    return avg_bias


# Finds alignments between datasets and uses them to construct
# panoramas. "Merges" datasets by correcting gene expression
# values.
def assemble(datasets, verbose=VERBOSE, view_match=False, knn=KNN,
             sigma=SIGMA, approx=APPROX, expr_datasets=None,
             ds_names=None, metric='angular'):
    if len(datasets) == 1:
        return datasets

    alignments, matches = find_alignments(datasets, knn=knn, approx=approx,
                                          verbose=verbose, metric=metric)

    ds_assembled = {}
    panoramas = []
    for i, j in alignments:
        if verbose:
            if ds_names is None:
                print('Processing datasets {}'.format((i, j)))
            else:
                print('Processing datasets {} <=> {}'.
                      format(ds_names[i], ds_names[j]))

        # Only consider a dataset a fixed amount of times.
        if not i in ds_assembled:
            ds_assembled[i] = 0
        ds_assembled[i] += 1
        if not j in ds_assembled:
            ds_assembled[j] = 0
        ds_assembled[j] += 1
        if ds_assembled[i] > 3 and ds_assembled[j] > 3:
            continue

        # See if datasets are involved in any current panoramas.
        panoramas_i = [panoramas[p] for p in range(len(panoramas))
                       if i in panoramas[p]]
        assert (len(panoramas_i) <= 1)
        panoramas_j = [panoramas[p] for p in range(len(panoramas))
                       if j in panoramas[p]]
        assert (len(panoramas_j) <= 1)

        if len(panoramas_i) == 0 and len(panoramas_j) == 0:
            if datasets[i].shape[0] < datasets[j].shape[0]:
                i, j = j, i
            panoramas.append([i])
            panoramas_i = [panoramas[-1]]

        # Map dataset i to panorama j.
        if len(panoramas_i) == 0:
            curr_ds = datasets[i]
            curr_ref = np.concatenate([datasets[p] for p in panoramas_j[0]])

            match = []
            base = 0
            for p in panoramas_j[0]:
                if i < p and (i, p) in matches:
                    match.extend([(a, b + base) for a, b in matches[(i, p)]])
                elif i > p and (p, i) in matches:
                    match.extend([(b, a + base) for a, b in matches[(p, i)]])
                base += datasets[p].shape[0]

            ds_ind = [a for a, _ in match]
            ref_ind = [b for _, b in match]

            bias = transform(curr_ds, curr_ref, ds_ind, ref_ind, sigma)
            datasets[i] = curr_ds + bias

            if expr_datasets:
                curr_ds = expr_datasets[i]
                curr_ref = np.concatenate([expr_datasets[p]
                                           for p in panoramas_j[0]])
                bias = transform(curr_ds, curr_ref, ds_ind, ref_ind, sigma)
                expr_datasets[i] = curr_ds + bias

            panoramas_j[0].append(i)

        # Map dataset j to panorama i.
        elif len(panoramas_j) == 0:
            curr_ds = datasets[j]
            curr_ref = np.concatenate([datasets[p] for p in panoramas_i[0]])

            match = []
            base = 0
            for p in panoramas_i[0]:
                if j < p and (j, p) in matches:
                    match.extend([(a, b + base) for a, b in matches[(j, p)]])
                elif j > p and (p, j) in matches:
                    match.extend([(b, a + base) for a, b in matches[(p, j)]])
                base += datasets[p].shape[0]

            ds_ind = [a for a, _ in match]
            ref_ind = [b for _, b in match]

            bias = transform(curr_ds, curr_ref, ds_ind, ref_ind, sigma)
            datasets[j] = curr_ds + bias

            if expr_datasets:
                curr_ds = expr_datasets[j]
                curr_ref = np.concatenate([expr_datasets[p]
                                           for p in panoramas_i[0]])
                bias = transform(curr_ds, curr_ref, ds_ind, ref_ind, sigma)
                expr_datasets[j] = curr_ds + bias

            panoramas_i[0].append(j)

        # Merge two panoramas together.
        else:
            curr_ds = np.concatenate([datasets[p] for p in panoramas_i[0]])
            curr_ref = np.concatenate([datasets[p] for p in panoramas_j[0]])

            # Find base indices into each panorama.
            base_i = 0
            for p in panoramas_i[0]:
                if p == i: break
                base_i += datasets[p].shape[0]
            base_j = 0
            for p in panoramas_j[0]:
                if p == j: break
                base_j += datasets[p].shape[0]

            # Find matching indices.
            match = []
            base = 0
            for p in panoramas_i[0]:
                if p == i and j < p and (j, p) in matches:
                    match.extend([(b + base, a + base_j)
                                  for a, b in matches[(j, p)]])
                elif p == i and j > p and (p, j) in matches:
                    match.extend([(a + base, b + base_j)
                                  for a, b in matches[(p, j)]])
                base += datasets[p].shape[0]
            base = 0
            for p in panoramas_j[0]:
                if p == j and i < p and (i, p) in matches:
                    match.extend([(a + base_i, b + base)
                                  for a, b in matches[(i, p)]])
                elif p == j and i > p and (p, i) in matches:
                    match.extend([(b + base_i, a + base)
                                  for a, b in matches[(p, i)]])
                base += datasets[p].shape[0]

            ds_ind = [a for a, _ in match]
            ref_ind = [b for _, b in match]

            # Apply transformation to entire panorama.
            bias = transform(curr_ds, curr_ref, ds_ind, ref_ind, sigma)
            curr_ds += bias
            base = 0
            for p in panoramas_i[0]:
                n_cells = datasets[p].shape[0]
                datasets[p] = curr_ds[base:(base + n_cells), :]
                base += n_cells

            if expr_datasets:
                curr_ds = np.concatenate([expr_datasets[p]
                                          for p in panoramas_i[0]])
                curr_ref = np.concatenate([expr_datasets[p]
                                           for p in panoramas_j[0]])
                bias = transform(curr_ds, curr_ref, ds_ind, ref_ind, sigma)
                curr_ds += bias
                base = 0
                for p in panoramas_i[0]:
                    n_cells = expr_datasets[p].shape[0]
                    expr_datasets[p] = curr_ds[base:(base + n_cells), :]
                    base += n_cells

            # Merge panoramas i and j and delete one.
            if panoramas_i[0] != panoramas_j[0]:
                panoramas_i[0] += panoramas_j[0]
                panoramas.remove(panoramas_j[0])

        # Visualize.
        if view_match:
            plot_mapping(curr_ds, curr_ref, ds_ind, ref_ind)

    return datasets


# Non-optimal dataset assembly. Simply accumulate datasets into a
# reference.
def assemble_accum(datasets, verbose=VERBOSE, knn=KNN, sigma=SIGMA,
                   approx=APPROX):
    if len(datasets) == 1:
        return datasets

    for i in range(len(datasets) - 1):
        j = i + 1

        if verbose:
            print('Processing datasets {}'.format((i, j)))

        ds1 = datasets[j]
        ds2 = np.concatenate(datasets[:i + 1])
        match = mnn(ds1, ds2, knn=knn, approx=approx)

        ds_ind = [a for a, _ in match]
        ref_ind = [b for _, b in match]

        bias = transform(ds1, ds2, ds_ind, ref_ind, sigma)
        datasets[j] = ds1 + bias

    return datasets


def interpret_alignments(datasets, expr_datasets, genes,
                         verbose=VERBOSE, knn=KNN, approx=APPROX,
                         n_permutations=None, metric='angular'):
    if n_permutations is None:
        n_permutations = float(len(genes) * 30)

    alignments, matches = find_alignments(
        datasets, knn=knn, approx=approx, verbose=verbose, metric=metric
    )

    for i, j in alignments:
        # Compute average bias vector that aligns two datasets together.
        ds_i = expr_datasets[i]
        ds_j = expr_datasets[j]
        if i < j:
            match = matches[(i, j)]
        else:
            match = matches[(j, i)]
        i_ind = [a for a, _ in match]
        j_ind = [b for _, b in match]
        avg_bias = np.absolute(
            np.mean(ds_j[j_ind, :] - ds_i[i_ind, :], axis=0)
        )

        # Construct null distribution and compute p-value.
        null_bias = (
                ds_j[np.random.randint(ds_j.shape[0], size=n_permutations), :] -
                ds_i[np.random.randint(ds_i.shape[0], size=n_permutations), :]
        )
        p = ((np.sum(np.greater_equal(
            np.absolute(np.tile(avg_bias, (n_permutations, 1))),
            np.absolute(null_bias)
        ), axis=0, dtype=float) + 1) / (n_permutations + 1))

        print('>>>> Stats for alignment {}'.format((i, j)))
        for k in range(len(p)):
            print('{}\t{}'.format(genes[k], p[k]))


if __name__ == '__main__':
    pass
