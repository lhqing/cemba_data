import glob
import pathlib
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import pandas as pd
from pybedtools import BedTool
from scipy import stats


def _calc_robust_mean(data, low=1, high=1):
    """Calculate robust mean using scipy.stats.sigmaclip,
    fall back to some simple mean if the data size is too small."""
    clipped, lower, upper = stats.sigmaclip(data, low=low, high=high)
    if clipped.size == 0:
        portion = 0.68
        skip_n = int(min(data.size / 2, np.round(data.size * (1 - portion) / 2)))
        if skip_n != 0:
            return np.sort(data)[skip_n:-skip_n].mean()
        else:
            return data.mean()
    return clipped.mean()


def _filter_mc_rate_table(dms_table, frac_cols, out_path, low_delta, high_delta, low_sd, high_sd):
    """Calculate robust mean, and then filter by robust mean and delta,
    keep all DMS that have at least one column pass filter"""
    mc_rate = dms_table.loc[:, frac_cols].astype(float)
    # remember methylpy results use -1 as nan...
    mc_rate_robust_mean = mc_rate.apply(lambda i: _calc_robust_mean(i[i > -1],
                                                                    low=low_sd,
                                                                    high=high_sd),
                                        axis=1)
    judge_df = (mc_rate.values > (mc_rate_robust_mean + high_delta).values[:, None]) | \
               (mc_rate.values < (mc_rate_robust_mean - low_delta).values[:, None])
    judge_df = ~mc_rate.where(judge_df).isna() & (mc_rate.values != -1)
    dms_table = dms_table[judge_df.sum(axis=1) > 0]
    dms_table.to_msgpack(out_path)
    return


def filter_dms(dms_path, output_path, cpu=1, pvalue_cutoff=0.01,
               low_mc_delta=0.3, high_mc_delta=0.3, low_mc_sd=1, high_mc_sd=1):
    """
    Pre-filter methylpy DMS results. Filter by both:
    1. P-value Cutoff
    2. Calculate significant hypo or hyper DMS for each column based on robust mean and delta value,
    keep DMS that have more than 1 column significant.

    Save all the output in chunks into pandas HDFStore

    Parameters
    ----------
    dms_path
        methylpy output "*_rms_results.tsv.gz"
    output_path
        Output path of the filtered results
    cpu
        number of CPU to use
    pvalue_cutoff
        P-value of DMS cutoff
    low_mc_delta
        Lower delta to the robust mean for column specific hypo-DMS
    high_mc_delta
        Upper delta to the robust mean for column specific hypo-DMS
    low_mc_sd
        Lower SD for column to be included in robust mean
    high_mc_sd
        Upper SD for column to be included in robust mean

    Returns
    -------

    """

    cpu = min(cpu, 30)
    chunksize = 100000
    with ProcessPoolExecutor(cpu) as executor:
        dms_table = pd.read_csv(dms_path, nrows=1, sep='\t', header=0, index_col=None)
        frac_cols = dms_table.columns[dms_table.columns.str.startswith('frac')]

        futures = {}
        for chunk_id, dms_table in enumerate(pd.read_csv(dms_path,
                                                         sep='\t', header=0,
                                                         index_col=None,
                                                         chunksize=chunksize)):
            dms_table = dms_table[dms_table['pvalue'] < pvalue_cutoff].copy()
            if dms_table.shape[0] < 1:
                # nothing remain after pvalue filter
                continue
            chunk_out_path = output_path + f'.chunk.{chunk_id}.msg'
            future = executor.submit(_filter_mc_rate_table, dms_table=dms_table,
                                     frac_cols=frac_cols, out_path=chunk_out_path,
                                     low_delta=low_mc_delta, high_delta=high_mc_delta,
                                     low_sd=low_mc_sd, high_sd=high_mc_sd)
            futures[future] = chunk_id
        for future in as_completed(futures):
            chunk_id = futures[future]
            try:
                future.result()
            except Exception as e:
                print(chunk_id, 'error')
                raise e

    output_path = pathlib.Path(output_path)
    chunk_paths = pd.Series({int(i.name.split('.')[-2]): str(i)
                             for i in output_path.parent.glob(f'{output_path.name}*chunk*msg')}).sort_index()
    with pd.HDFStore(output_path, complevel=5, complib='bzip2') as store:
        for chunk, path in chunk_paths.iteritems():
            store[f'chunk_{chunk}'] = pd.read_msgpack(path)
            subprocess.run(['rm', '-f', path])


def _get_specific_dms(dms_chunk, out_prefix, low_sd, high_sd, low_delta, high_delta):
    """Filter and save hyper/hypo DMS for each chunk and each column"""
    frac_df = dms_chunk.loc[:, dms_chunk.columns.str.startswith('frac')]
    robust_mean = frac_df.apply(lambda i: _calc_robust_mean(i[i > -1].astype(float),
                                                            low=low_sd,
                                                            high=high_sd),
                                axis=1)
    if robust_mean.isna().sum() != 0:
        print(dms_chunk.loc[robust_mean.isna()])
        raise
    hyper_judge_df = (frac_df.values > (robust_mean + high_delta).values[:, None])
    hypo_judge_df = (frac_df.values < (robust_mean - low_delta).values[:, None])
    hyper_judge_df = ~frac_df.where(hyper_judge_df).isna() & (frac_df.values != -1)
    hypo_judge_df = ~frac_df.where(hypo_judge_df).isna() & (frac_df.values != -1)

    for col_name, col in hyper_judge_df.iteritems():
        data = dms_chunk.loc[col[col].index, ['chr', 'pos']]
        data.to_msgpack(out_prefix + f'.hyper_dms.{col_name.lstrip("frac_")}.msg')

    for col_name, col in hypo_judge_df.iteritems():
        data = dms_chunk.loc[col[col].index, ['chr', 'pos']]
        data.to_msgpack(out_prefix + f'.hypo_dms.{col_name.lstrip("frac_")}.msg')
    return


def dms_to_dmr(dms_df, min_dms_distance=500, min_dms_number=2):
    """
    Merge DMS to DMR using bedtools merge

    Parameters
    ----------
    dms_df
        dms_dataframe have 2 columns ["chr", "pos"]
    min_dms_distance
    min_dms_number

    Returns
    -------

    """
    dms_df.columns = ['chrom', 'start']
    dms_df['end'] = dms_df['start']

    # merge dms to dmr
    dms_bed = BedTool.from_dataframe(dms_df)
    dmr_bed = dms_bed.merge(d=min_dms_distance, c=1, o='count')

    # filter dmr by n_dms
    dmr_df = dmr_bed.to_dataframe()
    dmr_df.columns = ['chrom', 'start', 'end', 'n_dms']
    dmr_df = dmr_df[dmr_df['n_dms'] >= min_dms_number].copy()
    return dmr_df


def extract_dmr(store_path, out_prefix,
                low_mc_sd=1, high_mc_sd=1,
                low_mc_delta=0.3, high_mc_delta=0.2,
                min_dms_distance=500, min_dms_number=2,
                save_dms=True, cpu=1):
    """
    Extract DMR from filtered DMS HDFStore

    Parameters
    ----------
    store_path
        Path of DMS pandas HDFStore
    out_prefix
        Output prefix of DMR files
    low_mc_sd
        Lower SD for column to be included in robust mean
    high_mc_sd
        Upper SD for column to be included in robust mean
    low_mc_delta
        Lower delta to the robust mean for column specific hypo-DMS
    high_mc_delta
        Upper delta to the robust mean for column specific hypo-DMS
    min_dms_distance
        Minimum distance between DMSs to merge them
    min_dms_number
        Minimum DMS number to keep DMR
    save_dms
        Whether save the DMS file too
    cpu
        Number of CPU to use

    Returns
    -------

    """
    out_prefix = out_prefix.rstrip('.')
    with ProcessPoolExecutor(cpu) as executor:
        future_dict = {}
        with pd.HDFStore(store_path, 'r') as store:
            key_order = pd.Series({int(key.split('_')[1]): key for key in store.keys()}).sort_index()
            merge_chunk_step = 10
            for chunk_id, chunk_start in enumerate(range(0, key_order.size, merge_chunk_step)):
                keys = key_order.iloc[chunk_start:chunk_start + merge_chunk_step]
                dms_chunk = pd.concat([store[key] for key in keys])
                _out_prefix = out_prefix + f'.{chunk_id}'
                future = executor.submit(_get_specific_dms,
                                         dms_chunk=dms_chunk, out_prefix=_out_prefix,
                                         low_sd=low_mc_sd, high_sd=high_mc_sd,
                                         low_delta=low_mc_delta, high_delta=high_mc_delta)
                future_dict[future] = chunk_id
        for future in as_completed(future_dict):
            chunk_id = future_dict[future]
            try:
                future.result()
            except Exception as e:
                print(f'{chunk_id} failed')
                raise e

    for dms_type in ['hypo', 'hyper']:
        chunk_paths = list(glob.glob(out_prefix + f'*{dms_type}_dms*.msg'))
        records = []
        for chunk_path in chunk_paths:
            *_, chunk_id, _, col_name, _ = pathlib.Path(chunk_path).name.split('.')
            chunk_id = int(chunk_id)
            records.append([chunk_id, col_name, chunk_path])
        chunk_path_df = pd.DataFrame(records, columns=['chunk_id', 'col', 'path']).set_index('chunk_id')

        for col, sub_df in chunk_path_df.groupby('col'):
            sub_df = sub_df.sort_index()
            chunks = []
            for path in sub_df['path']:
                chunks.append(pd.read_msgpack(path))
                subprocess.run(['rm', '-f', path])
            total_col_data = pd.concat(chunks)
            if save_dms:
                total_col_data.to_msgpack(out_prefix + f'.{dms_type}_dms.{col}.msg')
            dmr_df = dms_to_dmr(total_col_data,
                                min_dms_distance=min_dms_distance,
                                min_dms_number=min_dms_number)
            dmr_df.to_msgpack(out_prefix + f'.{dms_type}_dmr.{col}.msg')
    return
