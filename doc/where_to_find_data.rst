
Where to find data to use with `Exo_k`?
=======================================

`Exo_k` is a software library to handle radiative data for atmospheric applications --
data that you may directly use within the atmospheric models of `exo_k`, or convert to
the format needed by your favorite code. 

However, if you are not yourself producing such data (correlated-k tables, cross section tables,
or high-resolution spectra), you might ask yourself where to find some.

Fortunately, several projects have released publicly available datasets. Hereafter, I have
listed some of these public sources of radiative data. If they are useful
to your work, please cite the relevant articles listed below or
refer to the corresponding websites to know the proper way to acknowledge
their contribution (in addition to acknowledging `exo_k` if you use the library itself). 

Cross sections and k-coefficient tables from the Exomol project
---------------------------------------------------------------

While known for its linelists, the `Exomol <http://exomol.com/>`_ project
has recently put a lot of efforts in providing the community with
cross-section and k-coefficient tables computed for a wide range of temperatures
and pressures (`Chubb et al. 2020 <https://ui.adsabs.harvard.edu/abs/2020arXiv200900687C/abstract>`_).

These tables can be downloaded on their `website <http://exomol.com/data/data-types/opacity/>`_

k-coefficient tables from `ExoREM`
----------------------------------

`ExoREM <https://ui.adsabs.harvard.edu/abs/2015A%26A...582A..83B/abstract>`_
is a 1D model of the atmosphere of ultracool dwarfs and hot giant planets,
whose radiative transfer module is based on the correlated-k approach.

Recently we have made the k-coefficient tables used in ExoREM
publicly available: https://lesia.obspm.fr/exorem/

These tables use our extended HDF5 format that is completely compatible with ExoMOL data. 


k-coefficient tables from the LMDZ Generic GCM database
-------------------------------------------------------

The LMDZ Generic is a Generic climate model whose radiative transfer is based on correlated-k tables.
This climate model can be used to simulate planets with very different atmospheres.
For this, one of the main ingredients to change are the opacities used by the code. 

The LMDZ developpers team has gathered a sample of k-tables that have been used with the GCM along the years.
The specificities of each dataset is listed below. The data themselves can be found at:
https://www.lmd.jussieu.fr/~lmdz/planets/LMDZ.GENERIC/datagcm/corrk_data/

 * CO2_H2Ovar --> correlated-k coefficients initially used in
   `Wordsworth et al. (2013) <https://ui.adsabs.harvard.edu/abs/2013Icar..222....1W/abstract>`_
   to simulate the atmosphere of early Mars. Designed to simulate a CO2-dominated atmosphere,
   with H2O as a variable species (H2O mixing ratio < 10%,
   total pressure < 100bar, temperature < 400K).
   (H2O self and foreign far wing absorption --> not included)

 * Earth_110-710K --> initially used in 
   `Leconte et al. (2013) <https://ui.adsabs.harvard.edu/abs/2013Natur.504..268L/abstract>`_
   to simulate the atmosphere
   of Earth entering in runaway greenhouse.
   Designed to simulate the atmosphere of the Earth (N2-dominated + 376ppm of CO2)
   with H2O as a variable species (H2O mixing ratio up to 100%, temperature up to 710K).
   (H2O self and foreign far wing absorption --> need to be included separately with MT CKD)

 * megaCO2 --> designed to simulate a pure CO2 atmosphere.

 * early_earth_CO2_XXX_CH4_YYY --> initially used in
   `Charnay et al. (2013) <https://ui.adsabs.harvard.edu/abs/2013JGRD..11810414C/abstract>`_.
   Designed to simulate a N2-dominated atmosphere with CO2 and CH4 mixing ratio fixed
   (as in the name of the tables), with H2O as a variable species.
   Corresponds to case A, B, and C in Charnay et al. 2014 JGR.
   (H2O self and foreign far wing absorption --> not included)

 * present_day_Venus --> designed to simulate the atmosphere of present-day Venus.
   No variable gas. These coefficients were built by Martin Turbet and Sebastien Lebonnois
   to simulate present-day Venus with the generic model.
   Note that the coefficients are valid only for the temperature/pressure
   conditions along present-day Venus thermal profile.

 * Venus_H2Ovar_extreme_16g_with_cont --> designed to simulate the atmosphere
   of Venus with variable H2O. Source: Martin Turbet (not been validated yet).
   All continua are included in the spectra.

 * CO2_H2Ovar_extreme_16g_with_cont --> designed to simulate a thick CO2+H2O atmosphere 
   were both gases can become dominant.
   Source: Martin Turbet. Used in
   `Pluriel et al. 2019 <https://ui.adsabs.harvard.edu/abs/2019Icar..317..583P/abstract>`_
   (magma ocean planets),
   as well as `Turbet et al. 2019 <https://ui.adsabs.harvard.edu/abs/2020Icar..33513419T/abstract>`_
   (to simulate the effect of very large impacts on early Mars).
   All continua are included in the spectra.

 * CO2_SO2step_H2Ovar --> Designed to simulate a thick CO2 atmosphere with various concentrations
   of SO2. H2O is treated as a variable species. Used by
   `Kerber et al. (2015) <https://ui.adsabs.harvard.edu/abs/2015Icar..261..133K/abstract>`_
   to simulate volcanism on early Mars.

High-resolution cross sections from the petitRADTRANS database
--------------------------------------------------------------

PetitRADTRANS (`Molliere et al. (2019) <https://ui.adsabs.harvard.edu/abs/2019A%26A...627A..67M/abstract>`_)
is a 1D model calculating transmission and emission spectra of exoplanets.

A full documentation can be found there:
https://petitradtrans.readthedocs.io/en/latest/content/installation.html

Along with the code, one can download a set of high-resolution cross-section files:
https://petitradtrans.readthedocs.io/en/latest/content/installation.html#download-opacity-data

These binary files can be used by `Exo_k` to produce cross-section and k-coefficient tables
(see the :ref:`'neat examples'<neat_examples>` section).


Collision Induced Absorption (CIA)
----------------------------------

HITRAN has published a quite extensive list of CIA tables that can be found
`here <https://hitran.org/cia/>`_. They can be read directly by `exo_k`

Other more specific sets of CIA absorptions that have been used with LMDG
can be found `there <https://www.lmd.jussieu.fr/~lmdz/planets/LMDZ.GENERIC/datagcm/continuum_data/>`_.
The `README <https://www.lmd.jussieu.fr/~lmdz/planets/LMDZ.GENERIC/datagcm/continuum_data/README_continuum_files>`_
file contains a lot of information (but be sure to specify the right units!).

For a subset of very common background species (H2, He, N2, H2O), we have started to convert
some of the data above to an hdf5 format which contains unit information so that you can use
it more safely!

See: https://mycore.core-cloud.net/index.php/s/04jH1uXC9ZnGFm4