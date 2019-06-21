import pathlib
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import pandas as pd
from scipy import stats


def calc_robust_mean(data, low=1, high=1):
    clipped, lower, upper = stats.sigmaclip(data, low=low, high=high)
    if clipped.size == 0:
        portion = 0.68
        skip_n = int(min(data.size / 2, np.round(data.size * (1 - portion) / 2)))
        return np.sort(data)[skip_n:-skip_n].mean()
    return clipped.mean()


def filter_mc_rate_table(dms_table, frac_cols, out_path, low_delta, high_delta, low_sd, high_sd):
    mc_rate = dms_table.loc[:, frac_cols]
    mc_rate_robust_mean = mc_rate.apply(lambda i: calc_robust_mean(i[i > -1].astype(float),
                                                                   low=low_sd,
                                                                   high=high_sd),
                                        axis=1)
    if mc_rate_robust_mean.isna().sum() != 0:
        print(dms_table.loc[mc_rate_robust_mean.isna(), :])
        raise
    judge_df = (mc_rate.values > (mc_rate_robust_mean + high_delta).values[:, None]) | \
               (mc_rate.values < (mc_rate_robust_mean - low_delta).values[:, None])
    judge_df = ~mc_rate.where(judge_df).isna() & (mc_rate.values != -1)
    dms_table = dms_table[judge_df.sum(axis=1) > 0]
    dms_table.to_msgpack(out_path)
    print(out_path, 'saved')
    return


def filter_dms(dms_path, out_path, cpu=1, pvalue_cutoff=0.01,
               low_mc_delta=0.3, high_mc_delta=0.3, low_mc_sd=1, high_mc_sd=1):
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
            chunk_out_path = out_path + f'.chunk.{chunk_id}.msg'
            print(chunk_id, 'submitted')
            future = executor.submit(filter_mc_rate_table, dms_table=dms_table,
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

    out_path = pathlib.Path(out_path)
    chunk_paths = pd.Series({int(i.name.split('.')[-2]): str(i)
                             for i in out_path.parent.glob(f'{out_path.name}*chunk*msg')}).sort_index()
    print(f'Aggregate all chunks into {out_path}')
    with pd.HDFStore(out_path, complevel=5, complib='bzip2') as store:
        for chunk, path in chunk_paths.iteritems():
            store[f'chunk_{chunk}'] = pd.read_msgpack(path)
            subprocess.run(['rm', '-f', path])
