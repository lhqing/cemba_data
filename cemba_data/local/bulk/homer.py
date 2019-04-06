import pandas as pd
import re


def parse_homer_known_results(result_path, anno_path):
    anno_df = pd.read_csv(anno_path, index_col='gene_name')
    df = pd.read_csv(result_path, sep='\t')
    records = []
    p = re.compile(r'(?P<gene>.+)\((?P<DBD>.+)\)')
    for _, name in df.iloc[:, 0].iteritems():
        name, detail, source = name.split('/')
        m = p.search(name)
        try:
            gene = m.group('gene')
            dbd = m.group('DBD')
        except AttributeError:
            if name == 'ZNF652':
                records.append({'gene_name': 'ZNF652',
                                'DBD': None,
                                'gene_detail': None,
                                'source': None})
            else:
                print(name)
            continue
        records.append({'gene_name': gene,
                        'DBD': dbd,
                        'gene_detail': detail,
                        'source': source})
    name_df = pd.DataFrame(records)
    name_df = pd.concat([name_df, df], axis=1)
    name_df['gene_name'] = name_df['gene_name'].map(anno_df.reset_index().set_index('raw_name')['gene_name'].to_dict())
    name_df = name_df.groupby('gene_name').apply(lambda i: i.sort_values('Log P-value').iloc[0, :])
    name_df['gene_id'] = name_df['gene_name'].map(anno_df['gene_id'].to_dict())
    target_portion = name_df['% of Target Sequences with Motif'].str[:-1].astype(float)
    bg_portion = name_df['% of Background Sequences with Motif'].str[:-1].astype(float)
    name_df['fold_change'] = target_portion / bg_portion
    name_df['-lgp'] = -name_df['Log P-value']
    return name_df
