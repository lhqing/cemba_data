# CEMBA Data Processing (0.1.2)
- **Data mapping**
- **Data integration using xarray**
- **Data preprocessing, computation and visualization**
- **Regulatory region analysis**

## Giant workflow
![](/doc/image/pipeline.svg)
- **Stage 1** Mapping and QC. Prepare single cell base level methylation raw data.
- **Stage 2** Single cell level analysis. Prepare dataset (xarray based) and preprocessing for clustering analysis.
- **Stage 3** Cluster level analysis. Merge single cell data into cluster and identify the features of clusters.

## Documentation
[TODO] [Here](https://cemba-data.readthedocs.io/en/latest)

## Get Start
Before install this package, please install [Anaconda](https://www.anaconda.com/download/) to get most of the required python packages. **Python 3.6 or above is needed**, old version won't work.


```bash
# install
git clone https://github.com/lhqing/cemba_data.git
cd cemba_data
pip install .

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

