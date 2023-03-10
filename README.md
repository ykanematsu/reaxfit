# ReaxFit
Parameter-fitting module for lammps-reaxff with differential_evolution of scipy.
## Requirements
- lammps (library)
- numpy
- scipy
## Install
- use anaconda
```sh
conda install -c conda-forge lammps
conda install scipy numpy
pip3 install git+https://github.com/ykanematsu/reaxfit.git
```
- [option] install lammps from source code
```sh
mkdir build
cd build
cmake ../cmake -DLAMMPS_EXCEPTIONS=yes -DBUILD_SHARED_LIBS=yes -DMLIAP_ENABLE_PYTHON=yes -DPKG_PYTHON=yes -DPKG_MANYBODY=yes -DPKG_REAXFF=yes -DPYTHON_EXECUTABLE:FILEPATH=`which python3`
make -j4
make install-python
```
## Useage 
- See [A sample on colab](https://colab.research.google.com/github/ykanematsu/reaxfit/blob/master/reaxfit_sample.ipynb)
