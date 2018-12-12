import xarray as xr
import glob
import pathlib
import functools


@functools.lru_cache(maxsize=500)
def _get_ds_indexes(ds_path):
    with xr.open_dataset(ds_path) as ds:
        return ds.indexes


def _process_one_path(ds_path,
                      cell_list,
                      region_select_dict,
                      dataarray_select,
                      mc_type,
                      count_type):
    """
    :param ds_path:
    :param cell_list:
    :param region_select_dict: dict of region select or None
    :param dataarray_select: list of da names or None
    :param mc_type:
    :param count_type:
    :return:
    """
    # print('Process dataset:', ds_path)
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
                raise ValueError(_count_type, 'not in dataset mc_type coords.')
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
                    raise ValueError(da, 'not in dataset.data_vars')
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
    datasets = []
    for file_path in netcdf_files:
        dataset = _process_one_path(file_path,
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
    # netcdf_files
    if isinstance(netcdf_files, str):
        netcdf_files = sorted(glob.glob(netcdf_files))
    else:
        if not isinstance(netcdf_files, list):
            raise TypeError('netcdf_files should be either str or list, provided',
                            type(netcdf_files))
    if len(netcdf_files) == 0:
        print('No valid path provided.')
        return None
    for file_path in netcdf_files:
        if not pathlib.Path(file_path).exists():
            raise ValueError(file_path, 'do not exist.')

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
            return combined_data
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
        yield combined_data
