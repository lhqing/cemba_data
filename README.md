# CEMBA Data Processing (0.1.2)
[](http://www.network-science.de/ascii/)
<pre>
      ___           ___           ___                         ___
     /\__\         /\__\         /\  \         _____         /\  \
    /:/  /        /:/ _/_       |::\  \       /::\  \       /::\  \
   /:/  /        /:/ /\__\      |:|:\  \     /:/\:\  \     /:/\:\  \
  /:/  /  ___   /:/ /:/ _/_   __|:|\:\  \   /:/ /::\__\   /:/ /::\  \
 /:/__/  /\__\ /:/_/:/ /\__\ /::::|_\:\__\ /:/_/:/\:|__| /:/_/:/\:\__\
 \:\  \ /:/  / \:\/:/ /:/  / \:\~~\  \/__/ \:\/:/ /:/  / \:\/:/  \/__/
  \:\  /:/  /   \::/_/:/  /   \:\  \        \::/_/:/  /   \::/__/
   \:\/:/  /     \:\/:/  /     \:\  \        \:\/:/  /     \:\  \
    \::/  /       \::/  /       \:\__\        \::/  /       \:\__\
     \/__/         \/__/         \/__/         \/__/         \/__/
</pre>

## Giant workflow
![](/doc/image/pipeline.svg)
- **Stage 1** Mapping and QC. Prepare single cell base level methylation raw data.
- **Stage 2** Single cell level analysis. Prepare dataset (xarray based) and do clustering analysis together with other packages such as scanpy.
- **Stage 3** Cluster level analysis. Merge single cell data into cluster and identify characteristic of clusters.

## Documentation
[TODO] [Here](https://cemba-data.readthedocs.io/en/latest)

## Install
Before install this package, please install [Anaconda](https://www.anaconda.com/download/) to get most of the required python packages. **Python 3.6 or above is needed**, old version won't work.

```bash
# install
git clone https://github.com/lhqing/cemba_data.git
cd cemba_data
pip install .
```

## Quick Start
```bash
# to use the command line interface
yap -h
```

```python
# to use cemba_data as a python module
import cemba_data as cd
```

For more information, see the [documentation](https://cemba-data.readthedocs.io/en/latest).

## Important Reference
### Single cell:
- [scanpy](https://github.com/theislab/scanpy)
- [AnnData](https://github.com/theislab/anndata)
### Methylation and Epigenomics:
- [Methylpy](https://github.com/yupenghe/methylpy)
- [REPTILE](https://github.com/yupenghe/REPTILE)

