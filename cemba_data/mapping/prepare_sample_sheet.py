import pathlib
from collections import OrderedDict

import pandas as pd

import cemba_data

PACKAGE_DIR = pathlib.Path(cemba_data.__path__[0])

with open(PACKAGE_DIR / 'local/mapping/single_cell_info/plate_info_template.txt') as f:
    PLATE_INFO_TEMPLATE = f.read()

with open(PACKAGE_DIR / 'local/mapping/single_cell_info/sample_sheet_header.txt') as f:
    SAMPLESHEET_DEFAULT_HEADER = f.read()
SECTIONS = ['[CriticalInfo]', '[LibraryInfo]', '[PlateInfo]']
LIMITED_CHOICES = {
    'n_random_index': [8, 384],
    'input_plate_size': [384, 1536]}
CRITICAL_INFO_KEYS = ['n_random_index', 'input_plate_size',
                      'pool_id', 'tube_label', 'email']
# key (n_random_index, input_plate_size)
BARCODE_TABLE = {
    ('8', '384'): '/gale/netapp/home/hanliu/ref/inhouse/IMPORTANT_cemba_i7_i5_index.tsv'
}


def get_kv_pair(line):
    try:
        k, v = line.split('=')
        return k, v
    except ValueError:
        raise ValueError(f'Each key=value line must contain a "=" to separate key and value. Got {line}')


def read_plate_info(plateinfo_path):
    cur_section = ''
    cur_section_id = -1

    critical_info = {}
    library_info = OrderedDict()
    plate_header = True
    plate_info = []

    with open(plateinfo_path) as f:
        for line in f:
            line = line.strip('\n')
            if line == '' or line.startswith('#'):
                continue
            # print(line)

            # determine section
            if line.startswith('['):
                cur_section_id += 1
                if line == SECTIONS[cur_section_id]:
                    cur_section = line
                else:
                    raise ValueError(
                        f'Section name and order must be [CriticalInfo] [LibraryInfo] [PlateInfo], '
                        f'got {line} at {cur_section_id + 1} section.')
            elif cur_section == '[CriticalInfo]':
                k, v = get_kv_pair(line)
                if k not in CRITICAL_INFO_KEYS:
                    raise ValueError(f'Unknown key {k} in [CriticalInfo]')
                else:
                    critical_info[k] = v
            elif cur_section == '[LibraryInfo]':
                k, v = get_kv_pair(line)
                if (k in critical_info.keys()) or (k in library_info.keys()):
                    raise ValueError(f'Found duplicated key {k}')
                else:
                    library_info[k] = v
            elif cur_section == '[PlateInfo]':
                ll = line.split('\t')
                if plate_header:
                    plate_header = False
                    if (ll[0] != 'plate_id') or (ll[1] != 'primer_quarter'):
                        raise ValueError(
                            'The first line of [PlateInfo] must be header line, '
                            'and starts with "plate_id\tprimer_quarter"')
                plate_info.append(ll)
            else:
                raise ValueError(f'Got a malformed line: {line}')

    for k in CRITICAL_INFO_KEYS:
        if k not in critical_info:
            raise ValueError(f'[CriticalInfo] missing key-value pair "{k}"')

    header = plate_info[0]
    plate_info = pd.DataFrame(plate_info[1:], columns=plate_info[0])
    for k, v in library_info.items():
        if k in plate_info.columns:
            raise ValueError(f'Found duplicated key {k}')
        plate_info[k] = v
    col_order = header[:2] + list(library_info.keys()) + header[2:]
    plate_info = plate_info[col_order].copy()
    plate_info['sample_id_prefix'] = plate_info.apply(lambda i: '-'.join(i[2:].astype(str).tolist()), axis=1)

    # after getting sample_id_prefix, add critical info into plate_info too
    for k, v in critical_info.items():
        if k in plate_info.columns:
            raise ValueError(f'Found duplicated key {k}')
        plate_info[k] = v

    return critical_info, plate_info


def plate_384_random_index_8(plate_info, barcode_table):
    records = []

    # check plate_info primer compatibility
    for primer_quarter, n_plate in plate_info['primer_quarter'].value_counts().iteritems():
        if n_plate < 2:
            raise ValueError(f'{primer_quarter} only have 1 plate in the table, are you really sure?')
        elif n_plate == 2:
            pass
        else:
            raise ValueError(f'{primer_quarter} have {n_plate} plates in the table, that is impossible.')

    for primer_quarter, plate_pair in plate_info.groupby('primer_quarter'):
        plate1, plate2 = plate_pair['plate_id']
        # check plate pair info consistance
        for col_name, col in plate_pair.iteritems():
            if col.unique().size != 1:
                if col_name != 'plate_id':
                    print(f'{col_name} contains different information between {plate1} and {plate2}, '
                          f'Will put {plate1} prefix into sample_id. This should not happen normally.')
        for col in 'ABCDEFGH':
            for row in range(1, 13):
                plate_pos = f'{col}{row}'
                cur_row = barcode_table.loc[(primer_quarter, plate_pos)]
                i5_barcode = cur_row['i5_index_sequence']
                i7_barcode = cur_row['i7_index_sequence']
                sample_id_prefix = plate_pair['sample_id_prefix'].iloc[0]
                sample_id = f'{sample_id_prefix}-{plate1}-{plate2}-{plate_pos}'

                # THIS IS BASED ON FORMAT BCL2FASTQ NEEDS
                records.append({'Sample_ID': sample_id,
                                'index': i5_barcode,
                                'index2': i7_barcode,
                                'Sample_Project': plate_pair['tube_label'].iloc[0],
                                'Description': plate_pair['email'].iloc[0]})
    # THIS IS BASED ON FORMAT BCL2FASTQ NEEDS
    sample_sheet = pd.DataFrame(records)
    sample_sheet['Sample_Name'] = ''
    sample_sheet['Sample_Well'] = ''
    sample_sheet['I7_Index_ID'] = ''
    sample_sheet['I5_Index_ID'] = ''
    sample_sheet['I7_Index_ID'] = ''
    sample_sheet['Sample_Plate'] = 'Plate'

    miseq_sample_sheet = sample_sheet[['Sample_ID', 'Sample_Name', 'Sample_Plate',
                                       'Sample_Well', 'I7_Index_ID', 'index',
                                       'I5_Index_ID', 'index2', 'Sample_Project',
                                       'Description']].copy()

    lanes = []
    for i in range(1, 5):
        lane_df = miseq_sample_sheet.copy()
        lane_df['Lane'] = i
        lanes.append(lane_df)
    nova_sample_sheet = pd.concat(lanes)
    nova_sample_sheet = nova_sample_sheet[['Lane', 'Sample_ID', 'Sample_Name', 'Sample_Plate',
                                           'Sample_Well', 'I7_Index_ID', 'index',
                                           'I5_Index_ID', 'index2', 'Sample_Project',
                                           'Description']].copy()

    return miseq_sample_sheet, nova_sample_sheet


def make_sample_sheet(plate_info_paths, output_prefix, header_path=None):
    if isinstance(plate_info_paths, str):
        plate_info_paths = [plate_info_paths]
    critical_infos = []
    plate_infos = []

    miseq_sample_sheets = []
    nova_sample_sheets = []

    for plate_info_path in plate_info_paths:
        critical_info, plate_info = read_plate_info(plate_info_path)
        critical_infos.append(critical_info)
        plate_infos.append(plate_info)

        barcode_table_path = BARCODE_TABLE[critical_info['n_random_index'], critical_info['input_plate_size']]
        barcode_table = pd.read_csv(barcode_table_path, sep='\t')
        barcode_table['primer_quarter'] = barcode_table['Index_set'] + "_" + barcode_table['Index_quarter']
        barcode_table.set_index(['primer_quarter', 'plate_pos'], inplace=True)

        miseq_sample_sheet, nova_sample_sheet = plate_384_random_index_8(plate_info, barcode_table)
        miseq_sample_sheets.append(miseq_sample_sheet)
        nova_sample_sheets.append(nova_sample_sheet)
    miseq_sample_sheet = pd.concat(miseq_sample_sheets)
    nova_sample_sheet = pd.concat(nova_sample_sheets)

    with open(output_prefix + '.miseq.sample_sheet.csv', 'w') as output_f:
        if header_path is None:
            output_f.write(SAMPLESHEET_DEFAULT_HEADER)
        else:
            with open(header_path) as hf:
                output_f.write(hf.read())
        output_f.write(miseq_sample_sheet.to_csv(index=None))

    with open(output_prefix + '.novaseq.sample_sheet.csv', 'w') as output_f:
        if header_path is None:
            output_f.write(SAMPLESHEET_DEFAULT_HEADER)
        else:
            with open(header_path) as hf:
                output_f.write(hf.read())
        output_f.write(nova_sample_sheet.to_csv(index=None))
    return
