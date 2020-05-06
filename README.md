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

To learn how to use `exo_k`, you can generate the documentation by running `make -C doc` (or by simply running `make` in the `doc/` folder). The documentation will be generated in the doc/html folder (you can open the index.html file to check it out).
You can also simply follow the tutorial jupyter notebook `doc/tutorial-exo_k.ipynb`.

Have fun!