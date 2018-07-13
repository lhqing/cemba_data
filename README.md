NOTE: Currently, this pipeline is still under development, I can't promise the API won't change. I will remove this note once I fill confident about the general structure.

# CEMBA Data Processing (0.1.1)
- **Data integration using HDF5**
- **Data preprocessing, computation and visualization**
- **Regulatory region analysis**

## Giant workflow
[TODO]

## Documentation
[TODO] [Here](https://cemba-data.readthedocs.io/en/latest)

## Get Start
Before install this package, please install [Anaconda](https://www.anaconda.com/download/) to get most of the required python packages. **Python 3.6 or above is needed**, old version won't work.


```bash
# install
git clone https://github.com/lhqing/cemba_data.git
cd cemba_data
python setup.py install

# to use the command line interface
yap -h
```
```python
# to use cemba_data as a python module
import cemba_data as cd
from cemba_data.data.dataset import Dataset, Study

# read a dataset
dataset = Dataset(path_to_dataset)
# read a study
study = Study(path_to_study)
```
For more information, see the [documentation]((https://cemba-data.readthedocs.io/en/latest)).

