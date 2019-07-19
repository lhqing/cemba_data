"""
Scrublet

See licence here:
https://github.com/AllonKleinLab/scrublet

TODO Read their original code, and do some cleaning. Integrate with Dask

"""

from scrublet.helper_functions import *
from sklearn.decomposition import PCA
import pandas as pd
import logging


def calculate_posterior_mc_rate(mc_da, cov_da, normalize_per_cell=True, clip_norm_value=10):
    # TODO calculate cell_a, cell_b separately
    # so we can do post_rate only in a very small set of gene to prevent memory issue
    raw_rate = mc_da / cov_da
    cell_rate_mean = np.nanmean(raw_rate, axis=1)[:, None]  # this skip na
    cell_rate_var = np.nanvar(raw_rate, axis=1)[:, None]  # this skip na
    # based on beta distribution mean, var
    # a / (a + b) = cell_rate_mean
    # a * b / ((a + b) ^ 2 * (a + b + 1)) = cell_rate_var
    # calculate alpha beta value for each cell
    cell_a = (1 - cell_rate_mean) * (cell_rate_mean ** 2) / cell_rate_var - cell_rate_mean
    cell_b = cell_a * (1 / cell_rate_mean - 1)
    # cell specific posterior rate
    post_rate = (mc_da + cell_a) / (cov_da + cell_a + cell_b)
    if normalize_per_cell:
        # there are two ways of normalizing per cell, by posterior or prior mean:
        # prior_mean = cell_a / (cell_a + cell_b)
        # posterior_mean = post_rate.mean(dim=var_dim)
        # Here I choose to use prior_mean to normalize cell,
        # therefore all cov == 0 features will have normalized rate == 1 in all cells.
        # i.e. 0 cov feature will provide no info
        prior_mean = cell_a / (cell_a + cell_b)
        post_rate = post_rate / prior_mean
        if clip_norm_value is not None:
            post_rate[post_rate > clip_norm_value] = clip_norm_value
    return post_rate


def highly_variable_methylation_feature(X, feature_mean_cov, var_names,
                                        min_disp=0.5, max_disp=None, min_mean=0, max_mean=5,
                                        n_top_feature=None, bin_min_features=5, mean_binsize=0.05, cov_binsize=100):
    """
    Adapted from Scanpy, see license above
    The main difference is that, this function normalize dispersion based on both mean and cov bins.
    """
    # RNA is only scaled by mean, but methylation is scaled by both mean and cov
    log = logging.getLogger()
    log.info('extracting highly variable features')
    if n_top_feature is not None:
        log.info('If you pass `n_top_feature`, all cutoffs are ignored.')
    # warning for extremely low cov
    low_cov_portion = (feature_mean_cov < 30).sum() / feature_mean_cov.size
    if low_cov_portion > 0.2:
        log.warning(f'{int(low_cov_portion * 100)}% feature with < 10 mean cov, '
                    f'consider filter by cov before find highly variable feature. '
                    f'Otherwise some low coverage feature may be elevated after normalization.')
    cov = feature_mean_cov
    mean = np.mean(X, axis=0)
    var = np.var(X, axis=0, ddof=1)
    # now actually compute the dispersion
    if 0 in mean:
        mean[mean == 0] = 1e-12  # set entries equal to zero to small value
    # raw dispersion is the variance normalized by mean
    dispersion = var / mean
    if 0 in dispersion:
        dispersion[dispersion == 0] = np.nan
    dispersion = np.log(dispersion)
    # all of the following quantities are "per-feature" here
    df = pd.DataFrame(index=var_names)
    df['mean'] = mean
    df['dispersion'] = dispersion
    df['cov'] = cov
    # instead of n_bins, use bin_size, because cov and mc are in different scale
    df['mean_bin'] = (df['mean'] / mean_binsize).astype(int)
    df['cov_bin'] = (df['cov'] / cov_binsize).astype(int)
    # save bin_count df, gather bins with more than bin_min_features features
    bin_count = df.groupby(['mean_bin', 'cov_bin']).apply(lambda i: i.shape[0]).reset_index().sort_values(0, ascending=False)
    bin_count.head()
    bin_more_than = bin_count[bin_count[0] > bin_min_features]
    if bin_more_than.shape[0] == 0:
        raise ValueError(f'No bin have more than {bin_min_features} features, uss larger bin size.')
    # for those bin have too less features, merge them with closest bin in manhattan distance
    # this usually don't cause much difference (a few hundred features), but the scatter plot will look more nature
    index_map = {}
    for _, (mean_id, cov_id, count) in bin_count.iterrows():
        if count > 1:
            index_map[(mean_id, cov_id)] = (mean_id, cov_id)
        manhattan_dist = (bin_more_than['mean_bin'] - mean_id).abs() + (bin_more_than['cov_bin'] - cov_id).abs()
        closest_more_than = manhattan_dist.sort_values().index[0]
        closest = bin_more_than.loc[closest_more_than]
        index_map[(mean_id, cov_id)] = tuple(closest.tolist()[:2])
    # apply index_map to original df
    raw_bin = df[['mean_bin', 'cov_bin']].set_index(['mean_bin', 'cov_bin'])
    raw_bin['use_mean'] = pd.Series(index_map).apply(lambda i: i[0])
    raw_bin['use_cov'] = pd.Series(index_map).apply(lambda i: i[1])
    df['mean_bin'] = raw_bin['use_mean'].values
    df['cov_bin'] = raw_bin['use_cov'].values
    # calculate bin mean and std, now disp_std_bin shouldn't have NAs
    disp_grouped = df.groupby(['mean_bin', 'cov_bin'])['dispersion']
    disp_mean_bin = disp_grouped.mean()
    disp_std_bin = disp_grouped.std(ddof=1)
    # actually do the normalization
    _mean_norm = disp_mean_bin.loc[list(zip(df['mean_bin'], df['cov_bin']))]
    _std_norm = disp_std_bin.loc[list(zip(df['mean_bin'], df['cov_bin']))]
    df['dispersion_norm'] = (df['dispersion'].values - _mean_norm.values) / _std_norm.values
    dispersion_norm = df['dispersion_norm'].values.astype('float32')
    log.info('    finished')
    return dispersion_norm

class Scrublet():
    def __init__(self, mc, cov, bins, sim_doublet_ratio=2.0, n_neighbors=None, expected_doublet_rate=0.1, stdev_doublet_rate=0.02):
        ''' Initialize Scrublet object with counts matrix and doublet prediction parameters

        Parameters
        ----------
        counts_matrix : scipy sparse matrix or ndarray, shape (n_cells, n_genes)
            Matrix containing raw (unnormalized) UMI-based transcript counts.
            Converted into a scipy.sparse.csc_matrix.

        total_counts : ndarray, shape (n_cells,), optional (default: None)
            Array of total UMI counts per cell. If `None`, this is calculated
            as the row sums of `counts_matrix`.

        sim_doublet_ratio : float, optional (default: 2.0)
            Number of doublets to simulate relative to the number of observed
            transcriptomes.

        n_neighbors : int, optional (default: None)
            Number of neighbors used to construct the KNN graph of observed
            transcriptomes and simulated doublets. If `None`, this is
            set to round(0.5 * sqrt(n_cells))

        expected_doublet_rate : float, optional (default: 0.1)
            The estimated doublet rate for the experiment.

        stdev_doublet_rate : float, optional (default: 0.02)
            Uncertainty in the expected doublet rate.

        Attributes
        ----------
        predicted_doublets_ : ndarray, shape (n_cells,)
            Boolean mask of predicted doublets in the observed
            transcriptomes.

        doublet_scores_obs_ : ndarray, shape (n_cells,)
            Doublet scores for observed transcriptomes.

        doublet_scores_sim_ : ndarray, shape (n_doublets,)
            Doublet scores for simulated doublets.

        doublet_errors_obs_ : ndarray, shape (n_cells,)
            Standard error in the doublet scores for observed
            transcriptomes.

        doublet_errors_sim_ : ndarray, shape (n_doublets,)
            Standard error in the doublet scores for simulated
            doublets.

        threshold_: float
            Doublet score threshold for calling a transcriptome
            a doublet.

        z_scores_ : ndarray, shape (n_cells,)
            Z-score conveying confidence in doublet calls.
            Z = `(doublet_score_obs_ - threhsold_) / doublet_errors_obs_`

        detected_doublet_rate_: float
            Fraction of observed transcriptomes that have been called
            doublets.

        detectable_doublet_fraction_: float
            Estimated fraction of doublets that are detectable, i.e.,
            fraction of simulated doublets with doublet scores above
            `threshold_`

        overall_doublet_rate_: float
            Estimated overall doublet rate,
            `detected_doublet_rate_ / detectable_doublet_fraction_`.
            Should agree (roughly) with `expected_doublet_rate`.

        manifold_obs_: ndarray, shape (n_cells, n_features)
            The single-cell "manifold" coordinates (e.g., PCA coordinates)
            for observed transcriptomes. Nearest neighbors are found using
            the union of `manifold_obs_` and `manifold_sim_` (see below).

        manifold_sim_: ndarray, shape (n_doublets, n_features)
            The single-cell "manifold" coordinates (e.g., PCA coordinates)
            for simulated doublets. Nearest neighbors are found using
            the union of `manifold_obs_` (see above) and `manifold_sim_`.

        doublet_parents_ : ndarray, shape (n_doublets, 2)
            Indices of the observed transcriptomes used to generate the
            simulated doublets.

        doublet_neighbor_parents_ : list, length n_cells
            A list of arrays of the indices of the doublet neighbors of
            each observed transcriptome (the ith entry is an array of
            the doublet neighbors of transcriptome i).
        '''

        # initialize counts matrices
        self._M_obs = mc
        self._T_obs = cov
        rateb = calculate_posterior_mc_rate(mc, tc)
        disp = highly_variable_methylation_feature(rateb, np.mean(tc, axis=0), bins)
        idx = np.argsort(disp)[::-1]
        self._hvg_filter = idx[:2000]
        self._E_obs = rateb[:, self._hvg_filter]
        self._E_sim = None
        self._embeddings = {}
        self.sim_doublet_ratio = sim_doublet_ratio
        self.n_neighbors = n_neighbors
        self.expected_doublet_rate = expected_doublet_rate
        self.stdev_doublet_rate = stdev_doublet_rate
        if self.n_neighbors is None:
            self.n_neighbors = int(round(0.5*np.sqrt(self._E_obs.shape[0])))

    ######## Core Scrublet functions ########

    def scrub_doublets(self, synthetic_doublet_umi_subsampling=1.0, use_approx_neighbors=True,
                       distance_metric='euclidean', get_doublet_neighbor_parents=False,
                       min_counts=3, min_cells=3, min_gene_variability_pctl=85, log_transform=False,
                       mean_center=True, normalize_variance=True, n_prin_comps=30, verbose=True):
        ''' Standard pipeline for preprocessing, doublet simulation, and doublet prediction

        Automatically sets a threshold for calling doublets, but it's best to check
        this by running plot_histogram() afterwards and adjusting threshold
        with call_doublets(threshold=new_threshold) if necessary.

        Arguments
        ---------
        synthetic_doublet_umi_subsampling : float, optional (defuault: 1.0)
            Rate for sampling UMIs when creating synthetic doublets. If 1.0,
            each doublet is created by simply adding the UMIs from two randomly
            sampled observed transcriptomes. For values less than 1, the
            UMI counts are added and then randomly sampled at the specified
            rate.

        use_approx_neighbors : bool, optional (default: True)
            Use approximate nearest neighbor method (annoy) for the KNN
            classifier.

        distance_metric : str, optional (default: 'euclidean')
            Distance metric used when finding nearest neighbors. For list of
            valid values, see the documentation for annoy (if `use_approx_neighbors`
            is True) or sklearn.neighbors.NearestNeighbors (if `use_approx_neighbors`
            is False).

        get_doublet_neighbor_parents : bool, optional (default: False)
            If True, return the parent transcriptomes that generated the
            doublet neighbors of each observed transcriptome. This information can
            be used to infer the cell states that generated a given
            doublet state.

        min_counts : float, optional (default: 3)
            Used for gene filtering prior to PCA. Genes expressed at fewer than
            `min_counts` in fewer than `min_cells` (see below) are excluded.

        min_cells : int, optional (default: 3)
            Used for gene filtering prior to PCA. Genes expressed at fewer than
            `min_counts` (see above) in fewer than `min_cells` are excluded.

        min_gene_variability_pctl : float, optional (default: 85.0)
            Used for gene filtering prior to PCA. Keep the most highly variable genes
            (in the top min_gene_variability_pctl percentile), as measured by
            the v-statistic [Klein et al., Cell 2015].

        log_transform : bool, optional (default: False)
            If True, log-transform the counts matrix (log10(1+TPM)).
            `sklearn.decomposition.TruncatedSVD` will be used for dimensionality
            reduction, unless `mean_center` is True.

        mean_center : bool, optional (default: True)
            If True, center the data such that each gene has a mean of 0.
            `sklearn.decomposition.PCA` will be used for dimensionality
            reduction.

        normalize_variance : bool, optional (default: True)
            If True, normalize the data such that each gene has a variance of 1.
            `sklearn.decomposition.TruncatedSVD` will be used for dimensionality
            reduction, unless `mean_center` is True.

        n_prin_comps : int, optional (default: 30)
            Number of principal components used to embed the transcriptomes prior
            to k-nearest-neighbor graph construction.

        verbose : bool, optional (default: True)
            If True, print progress updates.

        Sets
        ----
        doublet_scores_obs_, doublet_errors_obs_,
        doublet_scores_sim_, doublet_errors_sim_,
        predicted_doublets_, z_scores_
        threshold_, detected_doublet_rate_,
        detectable_doublet_fraction_, overall_doublet_rate_,
        doublet_parents_, doublet_neighbor_parents_
        '''
        t0 = time.time()

        print_optional('Simulating doublets...', verbose)
        self.simulate_doublets(sim_doublet_ratio=self.sim_doublet_ratio,
                               synthetic_doublet_umi_subsampling=synthetic_doublet_umi_subsampling)

        print_optional('Embedding transcriptomes using PCA...', verbose)
        self.pipeline_pca(n_prin_comps=n_prin_comps)

        print_optional('Calculating doublet scores...', verbose)
        self.calculate_doublet_scores(use_approx_neighbors=use_approx_neighbors, distance_metric=distance_metric,
                                      get_doublet_neighbor_parents=get_doublet_neighbor_parents)
        self.call_doublets(verbose=verbose)

        t1=time.time()
        print_optional('Elapsed time: {:.1f} seconds'.format(t1 - t0), verbose)
        return self.doublet_scores_obs_, self.predicted_doublets_, self.doublet_scores_sim_

    def simulate_doublets(self, sim_doublet_ratio=None, synthetic_doublet_umi_subsampling=1.0):
        ''' Simulate doublets by adding the counts of random observed transcriptome pairs.

        Arguments
        ---------
        sim_doublet_ratio : float, optional (default: None)
            Number of doublets to simulate relative to the number of observed
            transcriptomes. If `None`, self.sim_doublet_ratio is used.

        synthetic_doublet_umi_subsampling : float, optional (defuault: 1.0)
            Rate for sampling UMIs when creating synthetic doublets. If 1.0,
            each doublet is created by simply adding the UMIs from two randomly
            sampled observed transcriptomes. For values less than 1, the
            UMI counts are added and then randomly sampled at the specified
            rate.

        Sets
        ----
        doublet_parents_
        '''

        if sim_doublet_ratio is None:
            sim_doublet_ratio = self.sim_doublet_ratio
        else:
            self.sim_doublet_ratio = sim_doublet_ratio

        n_obs = self._E_obs.shape[0]
        n_sim = int(n_obs * sim_doublet_ratio)
        pair_ix = np.random.randint(0, n_obs, size=(n_sim, 2))

        M1 = self._M_obs[pair_ix[:,0],:]
        M2 = self._M_obs[pair_ix[:,1],:]
        T1 = self._T_obs[pair_ix[:,0],:]
        T2 = self._T_obs[pair_ix[:,1],:]

        self._E_sim = calculate_posterior_mc_rate(M1+M2, T1+T2)[:, self._hvg_filter]
        self.doublet_parents_ = pair_ix
        return

    def pipeline_pca(self, n_prin_comps=50):
        X_obs = self._E_obs
        X_sim = self._E_sim
        pca = PCA(n_components=n_prin_comps).fit(X_obs)
        self.manifold_obs_ = pca.transform(X_obs)
        self.manifold_sim_ = pca.transform(X_sim)
        return

    def calculate_doublet_scores(self, use_approx_neighbors=True, distance_metric='euclidean', get_doublet_neighbor_parents=False):
        ''' Calculate doublet scores for observed transcriptomes and simulated doublets

        Requires that manifold_obs_ and manifold_sim_ have already been set.

        Arguments
        ---------
        use_approx_neighbors : bool, optional (default: True)
            Use approximate nearest neighbor method (annoy) for the KNN
            classifier.

        distance_metric : str, optional (default: 'euclidean')
            Distance metric used when finding nearest neighbors. For list of
            valid values, see the documentation for annoy (if `use_approx_neighbors`
            is True) or sklearn.neighbors.NearestNeighbors (if `use_approx_neighbors`
            is False).

        get_doublet_neighbor_parents : bool, optional (default: False)
            If True, return the parent transcriptomes that generated the
            doublet neighbors of each observed transcriptome. This information can
            be used to infer the cell states that generated a given
            doublet state.

        Sets
        ----
        doublet_scores_obs_, doublet_scores_sim_,
        doublet_errors_obs_, doublet_errors_sim_,
        doublet_neighbor_parents_

        '''

        self._nearest_neighbor_classifier(
            k=self.n_neighbors,
            exp_doub_rate=self.expected_doublet_rate,
            stdev_doub_rate=self.stdev_doublet_rate,
            use_approx_nn=use_approx_neighbors,
            distance_metric=distance_metric,
            get_neighbor_parents=get_doublet_neighbor_parents
            )
        return self.doublet_scores_obs_

    def _nearest_neighbor_classifier(self, k=40, use_approx_nn=True, distance_metric='euclidean',
                                     exp_doub_rate=0.1, stdev_doub_rate=0.03, get_neighbor_parents=False):
        manifold = np.vstack((self.manifold_obs_, self.manifold_sim_))
        doub_labels = np.concatenate((np.zeros(self.manifold_obs_.shape[0], dtype=int),
                                      np.ones(self.manifold_sim_.shape[0], dtype=int)))

        n_obs = np.sum(doub_labels == 0)
        n_sim = np.sum(doub_labels == 1)

        # Adjust k (number of nearest neighbors) based on the ratio of simulated to observed cells
        k_adj = int(round(k * (1+n_sim/float(n_obs))))

        # Find k_adj nearest neighbors
        neighbors = get_knn_graph(manifold, k=k_adj, dist_metric=distance_metric, approx=use_approx_nn, return_edges=False)

        # Calculate doublet score based on ratio of simulated cell neighbors vs. observed cell neighbors
        doub_neigh_mask = doub_labels[neighbors] == 1
        n_sim_neigh = doub_neigh_mask.sum(1)
        n_obs_neigh = doub_neigh_mask.shape[1] - n_sim_neigh

        rho = exp_doub_rate
        r = n_sim / float(n_obs)
        nd = n_sim_neigh.astype(float)
        ns = n_obs_neigh.astype(float)
        N = float(k_adj)

        # Bayesian
        q=(nd+1)/(N+2)
        Ld = q*rho/r/(1-rho-q*(1-rho-rho/r))

        se_q = np.sqrt(q*(1-q)/(N+3))
        se_rho = stdev_doub_rate

        se_Ld = q*rho/r / (1-rho-q*(1-rho-rho/r))**2 * np.sqrt((se_q/q*(1-rho))**2 + (se_rho/rho*(1-q))**2)

        self.doublet_scores_obs_ = Ld[doub_labels == 0]
        self.doublet_scores_sim_ = Ld[doub_labels == 1]
        self.doublet_errors_obs_ = se_Ld[doub_labels==0]
        self.doublet_errors_sim_ = se_Ld[doub_labels==1]

        # get parents of doublet neighbors, if requested
        neighbor_parents = None
        if get_neighbor_parents:
            parent_cells = self.doublet_parents_
            neighbors = neighbors - n_obs
            neighbor_parents = []
            for iCell in range(n_obs):
                this_doub_neigh = neighbors[iCell,:][neighbors[iCell,:] > -1]
                if len(this_doub_neigh) > 0:
                    this_doub_neigh_parents = np.unique(parent_cells[this_doub_neigh,:].flatten())
                    neighbor_parents.append(this_doub_neigh_parents)
                else:
                    neighbor_parents.append([])
            self.doublet_neighbor_parents_ = np.array(neighbor_parents)
        return

    def call_doublets(self, threshold=None, verbose=True):
        ''' Call trancriptomes as doublets or singlets

        Arguments
        ---------
        threshold : float, optional (default: None)
            Doublet score threshold for calling a transcriptome
            a doublet. If `None`, this is set automatically by looking
            for the minimum between the two modes of the `doublet_scores_sim_`
            histogram. It is best practice to check the threshold visually
            using the `doublet_scores_sim_` histogram and/or based on
            co-localization of predicted doublets in a 2-D embedding.

        verbose : bool, optional (default: True)
            If True, print summary statistics.

        Sets
        ----
        predicted_doublets_, z_scores_, threshold_,
        detected_doublet_rate_, detectable_doublet_fraction,
        overall_doublet_rate_
        '''

        if threshold is None:
            # automatic threshold detection
            # http://scikit-image.org/docs/dev/api/skimage.filters.html
            from skimage.filters import threshold_minimum
            try:
                threshold = threshold_minimum(self.doublet_scores_sim_)
                if verbose:
                    print("Automatically set threshold at doublet score = {:.2f}".format(threshold))
            except:
                self.predicted_doublets_ = None
                if verbose:
                    print("Warning: failed to automatically identify doublet score threshold. Run `call_doublets` with user-specified threshold.")
                return self.predicted_doublets_

        Ld_obs = self.doublet_scores_obs_
        Ld_sim = self.doublet_scores_sim_
        se_obs = self.doublet_errors_obs_
        Z = (Ld_obs - threshold) / se_obs
        self.predicted_doublets_ = Ld_obs > threshold
        self.z_scores_ = Z
        self.threshold_ = threshold
        self.detected_doublet_rate_ = (Ld_obs>threshold).sum() / float(len(Ld_obs))
        self.detectable_doublet_fraction_ = (Ld_sim>threshold).sum() / float(len(Ld_sim))
        self.overall_doublet_rate_ = self.detected_doublet_rate_ / self.detectable_doublet_fraction_

        if verbose:
            print('Detected doublet rate = {:.1f}%'.format(100*self.detected_doublet_rate_))
            print('Estimated detectable doublet fraction = {:.1f}%'.format(100*self.detectable_doublet_fraction_))
            print('Overall doublet rate:')
            print('\tExpected   = {:.1f}%'.format(100*self.expected_doublet_rate))
            print('\tEstimated  = {:.1f}%'.format(100*self.overall_doublet_rate_))

        return self.predicted_doublets_


