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
pip3 install reaxfit
```
For Windows, lammps from conda forge is not available. You can alternatively download the lammps binary with the python library from [lammps.org](https://packages.lammps.org/windows.html).
- [option] install jupyterlab and ase
It will be convenient to use reaxff with jupyterlab and ase.
```sh
conda install -c conda-forge jupyterlab ase
```
## Useage 
- See [A sample on colab](https://colab.research.google.com/github/ykanematsu/reaxfit/blob/main/reaxfit_sample.ipynb)
