import pandas as pd
import numpy as np
import holoviews as hv
from allensdk.core.reference_space import ReferenceSpace
from allensdk.api.queries.ontologies_api import OntologiesApi
from allensdk.core.structure_tree import StructureTree
from allensdk.api.queries.synchronization_api import SynchronizationApi
from allensdk.api.queries.image_download_api import ImageDownloadApi
from IPython.core.display import display, HTML
import json
import nrrd
import os
import logging

logging.getLogger().setLevel(logging.WARNING)

SLICE_REF_MICRON = {1: 2695, 2: 3290, 3: 3885, 4: 4480, 5: 5075, 6: 5670,
                    7: 6265, 8: 6860, 9: 7455, 10: 8050, 11: 8645, 12: 9240,
                    13: 9835, 14: 10430, 15: 11025, 16: 11620, 17: 12215, 18: 12810}
ALLEN_DIR = '/gale/netapp/scratch/hanliu/allen/'
# ALLEN_DIR = '/Users/hq/Documents/src/learn/allen_api'

ACRONYM_TO_FINEST = json.load(open(f'{ALLEN_DIR}/StructureTree/acronym_to_finest_nodes.json'))
GENE_SECTIONDATASET = pd.read_csv(f'{ALLEN_DIR}/SectionDataset/All_Mouse_Brain_ISH_Exp.csv')
ALLEN_GENE = pd.read_csv(f'{ALLEN_DIR}/SectionDataset/All_Mouse_Gene.csv',
                         index_col='gene_symbol', dtype={'id': str, 'gene_name': str,
                                                         'entrez_gene_id': str,
                                                         'homologene_group_id': str})
RESOLUTION = 25
ANNOTATION, META = nrrd.read(f'{ALLEN_DIR}/AnnotationGrid/annotation_{RESOLUTION}.nrrd')
MOUSE_CORONAL_REF_SPACE = 9
MOUSE_SAGITTAL_REF_SPACE = 10
MOUSE_CORONAL_ATLAS_ID = 602630314
MOUSE_CORONAL_HALF_ATLAS_ID = 1
MOUSE_SAGITTAL_ATLAS_ID = 2


def get_ref_space():
    oapi = OntologiesApi()
    structure_graph = oapi.get_structures_with_sets([1])  # 1 is the id of the adult mouse structure graph
    # This removes some unused fields returned by the query
    structure_graph = StructureTree.clean_structures(structure_graph)
    tree = StructureTree(structure_graph)
    rsp = ReferenceSpace(tree, ANNOTATION, [RESOLUTION, RESOLUTION, RESOLUTION])
    return rsp


def get_ref_coord_mask(acronym, coronal_slice, use_left=True):
    """
    get the ReferenceSpace coordinates and masks
    :param acronym:
    :param coronal_slice:
    :param use_left:
    :return:
    """
    if isinstance(acronym, str):
        acronym = [acronym]

    # get finest descendants
    total_nodes = []
    for acro in acronym:
        total_nodes += ACRONYM_TO_FINEST[acro]
    total_nodes = list(set(total_nodes))
    # get coronal pos and region mask
    coronal_micron = SLICE_REF_MICRON[coronal_slice]
    coronal_anno_pos = int(coronal_micron / RESOLUTION)
    coronal_slice_data = ANNOTATION[coronal_anno_pos, :, :]
    coronal_mask = np.isin(coronal_slice_data, total_nodes)
    if coronal_mask.sum() == 0:
        print(acronym, 'none of these regions show in the slice', coronal_slice)
        return None, None

    # get sagittal and horizontal pos and mask
    nonzero_x = np.nonzero(coronal_mask.sum(axis=0))
    if use_left:
        sagittal_anno_pos = int(np.percentile(nonzero_x, 75))
        sagittal_micron = sagittal_anno_pos * RESOLUTION
    else:
        sagittal_anno_pos = int(np.percentile(nonzero_x, 25))
        sagittal_micron = sagittal_anno_pos * RESOLUTION
    sagittal_slice_data = ANNOTATION[:, :, sagittal_anno_pos]
    sagittal_mask = np.isin(sagittal_slice_data, total_nodes)

    horizontal_anno_pos = int(np.percentile(np.nonzero(coronal_mask.sum(axis=1)), 50))
    horizontal_micron = horizontal_anno_pos * RESOLUTION

    ref_coord = (coronal_micron, sagittal_micron, horizontal_micron)

    return ref_coord, coronal_mask, sagittal_mask


def get_section_dataset_id(gene, slice_axis='all', remain_one=True):
    """
    get the section dataset id for a gene
    :param gene:
    :param slice_axis:
    :param remain_one:
    :return:
    """
    if slice_axis == 'all':
        slice_axis = ['coronal', 'sagittal']
    else:
        slice_axis = [slice_axis]
    gene = gene.lower()
    datasets = GENE_SECTIONDATASET[GENE_SECTIONDATASET['gene'].apply(lambda i: i.lower() == gene) &
                                   (GENE_SECTIONDATASET['plane'].isin(slice_axis))]
    if remain_one:
        datasets = datasets.drop_duplicates(subset='plane')
    return datasets.set_index('section_data_set_id')


def sync_ref_to_section_image(ref_coord, section_data_set_id, section_data_set_plane):
    """
    synchronize the ref coord to a image in the section dataset
    :param ref_coord:
    :param section_data_set_id:
    :param section_data_set_plane:
    :return:
    """
    sync_api = SynchronizationApi()
    x, y, z = ref_coord
    if section_data_set_plane in ['s', 'sagittal']:
        reference_space_id = MOUSE_SAGITTAL_REF_SPACE
    elif section_data_set_plane in ['c', 'coronal']:
        reference_space_id = MOUSE_CORONAL_REF_SPACE
    else:
        raise ValueError('Unknown section_data_set_plane value:',
                         section_data_set_plane, 'Supported are s|sagittal, c|coronal')
    result = sync_api.get_reference_to_image(reference_space_id,
                                             x, y, z,
                                             [section_data_set_id])
    return pd.Series(result[0]['image_sync'])


def sync_image_to_atlas(section_image_id, x, y, section_data_set_plane,
                        coronal_atlas=MOUSE_CORONAL_ATLAS_ID,
                        sagittal_atlas=MOUSE_SAGITTAL_ATLAS_ID):
    """
    synchronize the image in the section dataset to an annotated ref atlas
    :param section_image_id:
    :param x:
    :param y:
    :param section_data_set_plane:
    :param coronal_atlas:
    :param sagittal_atlas:
    :return:
    """
    if section_data_set_plane in ['c', 'coronal']:
        atlas_id = coronal_atlas
    elif section_data_set_plane in ['s', 'sagittal']:
        atlas_id = sagittal_atlas
    else:
        raise ValueError('Unknown section_data_set_plane value:',
                         section_data_set_plane, 'Supported are s|sagittal, c|coronal')
    sync_api = SynchronizationApi()
    result = sync_api.get_image_to_atlas(section_image_id, int(x), int(y), atlas_id)
    return pd.Series(result['image_sync']).fillna(atlas_id).astype(int)


def download_image(dataset_id, image_id, downsample=2, view=None, quality=80):
    """
    Download image from SectionDataset, cache all downloaded images
    :param dataset_id:
    :param image_id:
    :param downsample:
    :param view:
    :param quality:
    :return:
    """
    image_api = ImageDownloadApi()
    file_path = f'{ALLEN_DIR}/ImageDownload/{dataset_id}/' \
                f'{image_id}_view-{view}_downsample-{downsample}_quality-{quality}.jpg'
    try:
        os.mkdir(f'{ALLEN_DIR}/ImageDownload/{dataset_id}')
    except FileExistsError:
        pass
    if os.path.exists(file_path):
        return file_path
    kwargs = {'downsample': downsample, 'quality': quality}
    if view is None:
        view = 'ISH'
    if view != 'ISH':
        kwargs['view'] = view
    image_api.download_image(image_id, file_path=file_path, endpoint=None, **kwargs)
    return file_path


def download_atlas(dataset_id, image_id, downsample=1, annotation=True, quality=80):
    """
    Download annotated atlas image
    :param dataset_id:
    :param image_id:
    :param downsample:
    :param annotation:
    :param quality:
    :return:
    """
    image_api = ImageDownloadApi()
    file_path = f'{ALLEN_DIR}/AtlasDownload/{dataset_id}/' \
                f'{image_id}annotation-{annotation}_downsample-{downsample}_quality-{quality}.jpg'
    try:
        os.mkdir(f'{ALLEN_DIR}/AtlasDownload/{dataset_id}')
    except FileExistsError:
        pass
    if os.path.exists(file_path):
        return file_path
    kwargs = {'downsample': downsample,
              'quality': quality,
              'annotation': annotation,
              'atlas': dataset_id}
    image_api.download_atlas_image(image_id, file_path=file_path, **kwargs)
    return file_path


def allen_image_viewer_url(dataset_id, image_id,
                           base='http://mouse.brain-map.org/experiment/siv?id={dataset_id}&imageId={image_id}'):
    return base.format(dataset_id=dataset_id, image_id=image_id)


def allen_gene_url(gene, base='http://mouse.brain-map.org/gene/show/{gene_id}'):
    return base.format(gene_id=ALLEN_GENE.loc[gene]['id'])


def allen_3d_grid_url(dataset_id, base='aibe://mouse.brain-map.org/grid_data/v1/visualize/{dataset_id}?atlas=310'):
    return base.format(dataset_id=dataset_id)


def query_allen(gene, coronal_slice, region_acronym, remain_one=True, use_left=True,
                image_downsample=3, atlas_downsample=2, quality=80):
    # 1. In ReferenceSpace, get the xyz coord and mask (coronal)
    ref_coord, coronal_mask, sagittal_mask = get_ref_coord_mask(region_acronym, coronal_slice, use_left=use_left)
    # 2. Get SectionDataSet ID for the gene
    dataset_df = get_section_dataset_id(gene, slice_axis='all', remain_one=remain_one)
    # 3. Sync ref coord to SectionImage for each SectionDataSet
    image_df = dataset_df.apply(
        lambda i: sync_ref_to_section_image(ref_coord=ref_coord,
                                            section_data_set_id=i.name,
                                            section_data_set_plane=i['plane']), axis=1)
    # 4. Download SectionImage for each SectionDataSet
    image_df['ish_file_path'] = image_df.apply(
        lambda i: download_image(dataset_id=i['section_data_set_id'],
                                 image_id=i['section_image_id'],
                                 downsample=image_downsample,
                                 quality=quality), axis=1)
    image_df['expression_file_path'] = image_df.apply(
        lambda i: download_image(dataset_id=i['section_data_set_id'],
                                 image_id=i['section_image_id'],
                                 downsample=image_downsample, view='expression',
                                 quality=quality), axis=1)
    # 5. Sync AtlasImage for each SectionImage
    atlas_df = image_df.apply(
        lambda i: sync_image_to_atlas(
            section_image_id=i['section_image_id'], x=i['x'], y=i['y'],
            section_data_set_plane=dataset_df.loc[int(i['section_data_set_id'])]['plane']), axis=1)
    atlas_df.rename(columns={'section_data_set_id': 'atlas_id',
                             'section_image_id': 'atlas_image_id'}, inplace=True)
    # 6. Download AtlasImage for each SectionImage
    atlas_df['ref_file_path'] = atlas_df.apply(
        lambda i: download_atlas(dataset_id=i['atlas_id'],
                                 image_id=i['atlas_image_id'],
                                 downsample=atlas_downsample, annotation=True), axis=1)
    # 7. Merge all the results into single df
    total_df = pd.concat([dataset_df, image_df, atlas_df], axis=1)
    # 8. add some convenient urls
    total_df['allen_viewer_url'] = total_df.apply(
        lambda i: allen_image_viewer_url(i['section_data_set_id'], i['section_image_id']), axis=1)
    total_df['allen_gene_url'] = total_df.apply(
        lambda i: allen_gene_url(i['gene']), axis=1)
    total_df['allen_3d_grid_url'] = total_df.apply(
        lambda i: allen_3d_grid_url(i['section_data_set_id']), axis=1)
    total_df.reset_index(drop=True, inplace=True)

    return total_df


def plot_query_result(query_df, height=300, aspect=1.3):
    width = int(height * aspect)
    tmp_df = query_df.rename(columns={'ish_file_path': 'ISH',
                                      'expression_file_path': 'Expression',
                                      'ref_file_path': 'Atlas'})
    # plain table
    df = tmp_df.melt(id_vars=['gene', 'plane', 'section_data_set_id', 'section_image_id'],
                     value_vars=['ISH', 'Expression', 'Atlas'])
    curve_dict = {}
    for _, row in df.iterrows():
        curve_dict[(row['gene'], row['plane'], row['variable'])] = \
            hv.RGB.load_image(row['value']).options(width=width, height=height)
    kdims = [hv.Dimension('Gene'),
             hv.Dimension('Plane'),
             hv.Dimension('Image Type')]
    hm = hv.HoloMap(curve_dict, kdims=kdims)
    return hm


def print_allen_url(df):
    """
    Print some Allen url for convenience, only use in notebook
    :param df:
    :return:
    """
    display(HTML('<h4>Some url for your convenience</h4>'))
    display(HTML('<p>For Gene summary:</p>'))
    for gene in df['gene'].unique():
        url = df[df['gene'] == gene]['allen_gene_url'].iloc[0]
        text = f'Go to Allen\'s {gene} summary page.'
        display(HTML(f'<a href="{url}">{text}</a>'))
    display(HTML('<p>For high resolution ISH images viewer:</p>'))
    for _, row in df.iterrows():
        plane, gene, dataset_id, url = row[['plane', 'gene', 'section_data_set_id', 'allen_viewer_url']]
        text = f'Go to Allen\'s Viewer for {plane} {gene} experiment (id={dataset_id}).'
        display(HTML(f'<a href="{url}">{text}</a>'))
    display(HTML('<p>For 3D expression viewer '
                 '(Need to install <a href="http://mouse.brain-map.org/static/brainexplorer">'
                 'Allen Brain Explorer</a>):</p>'))
    for section_data_set_id in df['section_data_set_id'].unique():
        url = df[df['section_data_set_id'] == section_data_set_id]['allen_3d_grid_url'].iloc[0]
        gene = df[df['section_data_set_id'] == section_data_set_id]['gene'].iloc[0]
        text = f'Go to Brain Explorer for {gene} 3D expressions (id={section_data_set_id}).'
        display(HTML(f'<a href="{url}">{text}</a>'))
    return



