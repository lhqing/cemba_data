import pandas as pd
from pathlib import Path
import shutil
from .mc_bulk_multigroup_template import MERGE_TEMPLATE, MERGE_EXTRACT_TEMPLATE


'''
the group_file is in csv format with the first column as allcpath and different group the others.
the group_file needs a header.
'''

def merge_bulk_multigroup(group_path, output_path, chrom_size_path, 
                          n_cpu=10, elem_snakegroup_num = 50,
                          cate_snakegroup_num = 10, ):
    

    outdir = Path(output_path)
    outdir.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(group_path, outdir/'GROUP.csv')

    df = pd.read_csv(group_path)
    df = df.rename(columns={df.columns[0]:'_path'})
    sample_cates = df.columns[1:]

    df['_elem'] = pd.factorize(df[sample_cates].astype(str).apply('-'.join, axis=1))[0]
    countdict = df['_elem'].value_counts().to_dict()
    df['_elem'] = df['_elem'].apply(lambda x: f'{x}_{countdict[x]}')

    elem_grp_df = df.groupby('_elem')['_path'].apply(lambda x: x.unique()).to_frame()
    elem_grp_df.index.name = '_sample'
    elem_grp_df['_cate'] = '_elem'
    elem_grp_df = elem_grp_df.reset_index()[['_cate','_sample','_path']]

    df = df[df.columns[1:]].drop_duplicates()

    df['_path'] = output_path+'/_elem/'+df['_elem']+'.allc.tsv.gz'

    cate_grp_df = []
    for cate in sample_cates:
        catedf = df[[cate,'_path']].groupby(cate)['_path'].apply(lambda x: x.unique()).to_frame()
        catedf['_cate'] = cate
        catedf.index.name = '_sample'
        catedf = catedf.reset_index()
        cate_grp_df.append(catedf)
    cate_grp_df = pd.concat(cate_grp_df).reset_index(drop=True)[['_cate','_sample','_path']]   


    def prepare_snakefiles(grp_df, output_path, tag, n_per_snake=None, template=MERGE_TEMPLATE):
        outdir = Path(output_path)
        snkdir = outdir/'snakefiles'
        snkdir.mkdir(exist_ok=True)

        for cate in grp_df['_cate'].unique():
            catedir = outdir/cate
            catedir.mkdir(exist_ok=True)

        for _,(cate,sample,paths) in grp_df.iterrows():
            catedir = outdir/cate
            with open(catedir/f'{sample}.pathlist','w') as f:
                f.write('\n'.join(paths))

        if n_per_snake is None:
            n_per_snake = len(grp_df)

        snk_ids = []
        for i, snkdf in grp_df.groupby(grp_df.index%n_per_snake):
            snk_id = f'{tag}_{i}'

            tocp_df = snkdf[snkdf['_path'].apply(len)==1]
            tomg_df = snkdf[snkdf['_path'].apply(len)>1]

            with open(snkdir/f'{snk_id}.snakefile', 'w') as f:
                f.write(
f'''merge_allc_cpu = {n_cpu}
mcg_context = 'CGN'
chrom_size_path = '{chrom_size_path}'
merge_sample_prefixes = [{','.join("'"+tomg_df['_cate']+'/'+tomg_df['_sample']+"'")}]
copy_sample_prefixes = [{','.join("'"+tocp_df['_cate']+'/'+tocp_df['_sample']+"'")}]
group = "{snk_id}"
'''
                )
                f.write(template)
            snk_ids.append(snk_id)

        return snk_ids

    elem_snk_ids = prepare_snakefiles(elem_grp_df, output_path, 'elem',elem_snakegroup_num, template=MERGE_TEMPLATE)
    cate_snk_ids = prepare_snakefiles(cate_grp_df, output_path, 'cate',cate_snakegroup_num, template=MERGE_EXTRACT_TEMPLATE)

    def prepare_commands(snake_ids):
        cmds = [f'snakemake  -d {outdir.resolve()} --snakefile {outdir.resolve()}/snakefiles/{snkid}.snakefile '
                f'-j {n_cpu} --default-resources mem_mb=100 --resources mem_mb=1000 --rerun-incomplete' \
                for snkid in snake_ids]
        return cmds

    
    
    with open(outdir/'run_snakemake_cmds_1.txt', 'w') as f:
        f.write('\n'.join(prepare_commands(elem_snk_ids)))
    with open(outdir/'run_snakemake_cmds_2.txt', 'w') as f:
        f.write('\n'.join(prepare_commands(cate_snk_ids)))    
