"""
Using strategy described here:
http://xarray.pydata.org/en/stable/io.html#combining-multiple-files

Input:
1. cell list, consistent to all datasets and filtered by some criteria.
   if None, use all cells from the datasets.
2. datasets

Step:
1. split the cell list into chunks
2. for each cell list chunk, search all datasets.
   Because the search is only performed on dataset index,
   and the index is cached, this is actually fast.
3. return the loaded dataset for each chunk and concat into one study.

Output:
A study, actually have the same structure as dataset.

"""

import xarray as xr
import glob
import pathlib
import functools
import logging
from ...tools.hdf5.netndf import MCDS

log = logging.getLogger()


@functools.lru_cache(maxsize=500)
def _get_ds_indexes(ds_path):
    """
    Get full dataset indexes.
    """
    with xr.open_dataset(ds_path) as ds:
        return ds.indexes


def _process_one_path_rna(ds_path,
                          cell_list,
                          region_select_dict,
                          dataarray_select,
                          mc_type,
                          count_type):
        # TODO write function for RNA da
        return


def _process_one_path_mc(ds_path,
                         cell_list,
                         region_select_dict,
                         dataarray_select,
                         mc_type,
                         count_type):
    """
    query one dataset with a list of cell_id.
    if no cell in this dataset, return None,
    else return loaded subset of the dataset
    """

    log.debug(f'Process dataset:{ds_path}')
    ds_indexes = _get_ds_indexes(ds_path)
    selection_dict = {}

    # select cells
    ds_cell_index = ds_indexes['cell']
    ds_select_cell_list = [cell for cell in cell_list if cell in ds_cell_index]
    if len(ds_select_cell_list) == 0:
        # no cell matched for this dataset
        return None
    selection_dict['cell'] = ds_select_cell_list

    # select mc_type
    if mc_type is not None:
        ds_mc_type_index = ds_indexes['mc_type']
        for _mc_type in mc_type:
            if _mc_type not in ds_mc_type_index:
                raise ValueError(_mc_type, 'not in dataset mc_type coords.')
        selection_dict['mc_type'] = mc_type

    # select count_type
    if count_type is not None:
        ds_count_type_index = ds_indexes['count_type']
        for _count_type in count_type:
            if _count_type not in ds_count_type_index:
                raise ValueError(_count_type, 'not in dataset count_type coords.')
        selection_dict['count_type'] = count_type

    # regions
    if region_select_dict is not None:
        for region, region_list in region_select_dict.items():
            if region not in ds_indexes:
                raise ValueError(region, 'not in dataset dims.')
            selection_dict[region] = region_list

    with xr.open_dataset(ds_path) as ds:
        # apply selection dict
        selected_ds = ds.sel(**selection_dict)
        # select dataarray
        if dataarray_select is not None:
            for da in dataarray_select:
                if da not in selected_ds.data_vars:
                    raise ValueError(f'{da} not in dataset.data_vars, '
                                     f'the acceptable data_vars are {selected_ds.data_vars}')
            selected_ds = selected_ds[dataarray_select]
        # really load selected data
        selected_ds.load()
        return selected_ds


def _process_cell_list(netcdf_files,
                       cell_list,
                       region_select_dict,
                       dataarray_select,
                       mc_type,
                       count_type):
    """
    For a given cell list and datasets, search and combine dataset subsets for this cell list.
    """
    datasets = []
    for file_path in netcdf_files:
        # TODO add functions for different data type
        dataset = _process_one_path_mc(file_path,
                                       cell_list,
                                       region_select_dict,
                                       dataarray_select,
                                       mc_type,
                                       count_type)
        if dataset is not None:
            datasets.append(dataset)

    if len(datasets) != 0:
        combined = xr.concat(datasets, dim='cell')
        return combined
    else:
        # print('No data remained.')
        return None


def _prepare_study_generator(
        netcdf_files,
        cell_list,
        region_select_dict,
        dataarray_select,
        mc_type,
        count_type,
        iter_chunksize):
    """
    generator for iterating through all netcdf_files, only call this from prepare study
    """
    cell_list_chunks = (cell_list[i:i + iter_chunksize]
                        for i in range(0, len(cell_list), iter_chunksize))
    for cell_chunk in cell_list_chunks:
        combined_data = _process_cell_list(netcdf_files,
                                           cell_chunk,
                                           region_select_dict,
                                           dataarray_select,
                                           mc_type,
                                           count_type)
        if combined_data is None:
            continue
        yield MCDS(combined_data)


def prepare_study(netcdf_files,
                  cell_list=None,
                  region_select_dict=None,  # dict if multiple da
                  dataarray_select=None,  # list if multiple da
                  mc_type=None,
                  count_type=None,
                  save_path=None,
                  compression=False,
                  compression_level=1,
                  iter_chunksize=None):
    """
    Prepare study based on a list of cell id and a list of netcdf_files

    Parameters
    ----------
    netcdf_files
        List of MCDS paths for checking cell_list
    cell_list
        List of cell_ids used for filtering MCDS cell dimension
    region_select_dict
        List of region/feature ids use for filtering MCDS corresponding dimension
    dataarray_select
        Str/List of dataarray_name
    mc_type
        Str/List of mc_types
    count_type
        Str/List of count_types
    save_path
        Output path if you want to save the study
    compression
        Compression the MCDS output
    compression_level
        Compression level
    iter_chunksize
        If not None, this function will return a iterator with this chunksize

    Returns
    -------
    dataset if iter_chunksize is None
        study is still a dataset whose cells combined from different datasets by certain purpose.
    iterator of dataset if iter_chunksize is not None
    """
    # netcdf_files
    if isinstance(netcdf_files, str):
        netcdf_files = sorted(glob.glob(netcdf_files))
    else:
        if not isinstance(netcdf_files, list):
            raise TypeError(f'netcdf_files should be either str or list, '
                            f'provided {type(netcdf_files)}')

    if len(netcdf_files) == 0:
        print('No valid path provided.')
        return None
    for file_path in netcdf_files:
        if not pathlib.Path(file_path).exists():
            raise ValueError(f'{file_path} do not exist.')

    if cell_list is None:
        # no cell list, select all cells
        cell_list = []
        for file_path in netcdf_files:
            ds_index = _get_ds_indexes(file_path)
            cell_list += ds_index['cell'].tolist()
        cell_list = list(set(cell_list))

    if isinstance(dataarray_select, str):
        dataarray_select = [dataarray_select]
    if isinstance(mc_type, str):
        mc_type = [mc_type]
    if isinstance(count_type, str):
        count_type = [count_type]

    if iter_chunksize is None:
        combined_data = _process_cell_list(netcdf_files,
                                           cell_list,
                                           region_select_dict,
                                           dataarray_select,
                                           mc_type,
                                           count_type)
        if combined_data is None:
            return None
        if save_path is not None:
            encoding_dict = {}
            for da in combined_data.data_vars:
                encoding_dict[da] = {'zlib': compression,
                                     'complevel': compression_level}
            combined_data.to_netcdf(save_path, encoding=encoding_dict)
            return None
        else:
            return MCDS(combined_data)
    else:
        # generator mode
        return _prepare_study_generator(
            netcdf_files=netcdf_files,
            cell_list=cell_list,
            region_select_dict=region_select_dict,  # dict if multiple da
            dataarray_select=dataarray_select,  # list if multiple da
            mc_type=mc_type,
            count_type=count_type,
            iter_chunksize=iter_chunksize
        )
