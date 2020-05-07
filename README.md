# Exo_k

Author: Jeremy Leconte (CNRS/LAB/Univ. Bordeaux)

The goal of `exo_k` is to provide a library to:

* Interpolate efficiently and easily in correlated-k and cross section tables
* Convert easily correlated-k and cross section tables from one format to another (pickle, hdf5, LMDZ GCM format).
* Adapt precomputed correlated-k tables to your need:
  * by changing the resolution
  * by changing the pressure/temperature grid
* To create tables for a mix of gases using tables for individual gases.
* Test various assumptions in an integrated radiative transfer framework.

## Installation

Exo_k can be installed using pip:
```
pip install -e .
```
Or by running the [setup.py](./setup.py) script:
```
python setup.py install
```
To generate the documentation, you will need to install the following packages:
```
pip install nbsphinx sphinx-autoapi sphinx_rtd_theme
conda install sphinx # installs more (required) dependencies than pip
```
You can then generate the documentation by running:
```
python setup.py doc
```
(or by simply running `make` in the `doc/` folder). The documentation will be generated in the doc/html folder (you can open the [index.html](doc/html/index.html) file to check it out).

## Usage

To learn how to use `exo_k`, you can follow the [tutorial jupyter notebook](doc/tutorial-exo_k.ipynb) or read the [documentation](doc/html/index.html) (see `Installation` above to generate it).

Have fun!