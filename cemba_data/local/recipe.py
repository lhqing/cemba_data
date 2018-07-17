from .prepare_study import prepare_study
from ..data.hdf5 import Study
import os


def region_preprocess_recipe(project_name, study_name, cell_meta,
                             cell_na_max=0.3, region_na_max=0.05,
                             chrom100k_ch_cov_cutoff=100,
                             chrom100k_cg_cov_cutoff=10,
                             gene_ch_cov_cutoff=10,
                             gene_cg_cov_cutoff=5,
                             save=True,
                             out_dir=''):

    if os.path.exists(out_dir + f'mch_chrom100k_{study_name}.h5ad'):
        mch_100k = Study.from_file(out_dir + f'mch_chrom100k_{study_name}.h5ad')
    else:
        # prepare chrom100k mCH study
        mch_100k = prepare_study(project_name=project_name, study_name=study_name, cell_list=cell_meta.index,
                                 region='chrom100k', region_context='CHN',
                                 coverage_cutoff=chrom100k_ch_cov_cutoff, out_dir=None)
        mch_100k.filter_cell(cell_na_max)  # cell no more than 50% na
        mch_100k.filter_region(region_na_max)  # region no more than 5% na
        mch_100k.add_row_feature('overall_mch', cell_meta['mCH/CH'])
        mch_100k.normalize_by_row_feature(feature_name='overall_mch')
        mch_100k.imputation(strategy='mean', axis=0)

    if os.path.exists(out_dir + f'mcg_chrom100k_{study_name}.h5ad'):
        mcg_100k = Study.from_file(out_dir + f'mcg_chrom100k_{study_name}.h5ad')
    else:
        # prepare chrom100k mCG study
        mcg_100k = prepare_study(project_name='CEMBA_RS1_METHY', study_name=study_name, cell_list=cell_meta.index,
                                 region='chrom100k', region_context='CGN',
                                 coverage_cutoff=chrom100k_cg_cov_cutoff, out_dir=None)
        mcg_100k.filter_cell(cell_na_max)
        mcg_100k.filter_region(region_na_max)
        mcg_100k.add_row_feature('overall_mcg', cell_meta['mCG/CG'])
        mcg_100k.normalize_by_row_feature(feature_name='overall_mcg')
        mcg_100k.imputation(strategy='mean', axis=0)

    if os.path.exists(out_dir + f'mch_gene_{study_name}.h5ad'):
        mch_gene = Study.from_file(out_dir + f'mch_gene_{study_name}.h5ad')
    else:
        # prepare gene mCH study
        mch_gene = prepare_study(project_name='CEMBA_RS1_METHY', study_name=study_name, cell_list=cell_meta.index,
                                 region='gene', region_context='CHN',
                                 coverage_cutoff=gene_ch_cov_cutoff, out_dir=None)
        mch_gene.add_row_feature('overall_mch', cell_meta['mCH/CH'])
        mch_gene.normalize_by_row_feature(feature_name='overall_mch')

    if os.path.exists(out_dir + f'mcg_gene_{study_name}.h5ad'):
        mcg_gene = Study.from_file(out_dir + f'mcg_gene_{study_name}.h5ad')
    else:
        # prepare gene mCG study
        mcg_gene = prepare_study(project_name='CEMBA_RS1_METHY', study_name=study_name, cell_list=cell_meta.index,
                                 region='gene', region_context='CGN',
                                 coverage_cutoff=gene_cg_cov_cutoff, out_dir=None)
        mcg_gene.add_row_feature('overall_mcg', cell_meta['mCG/CG'])
        mcg_gene.normalize_by_row_feature(feature_name='overall_mcg')

    if save:
        mch_100k.save(out_dir + f'mch_chrom100k_{study_name}.h5ad', compression='gzip', compression_opts=None)
        mcg_100k.save(out_dir + f'mcg_chrom100k_{study_name}.h5ad', compression='gzip', compression_opts=None)
        mch_gene.save(out_dir + f'mch_gene_{study_name}.h5ad', compression='gzip', compression_opts=None)
        mcg_gene.save(out_dir + f'mcg_gene_{study_name}.h5ad', compression='gzip', compression_opts=None)

    return mch_100k.to_ann(), mcg_100k.to_ann(), mch_gene.to_ann(), mcg_gene.to_ann()


