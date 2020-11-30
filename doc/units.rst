Notes on units and formats
==========================

Keeping track of units
----------------------

It is fair to say that units have often been a source of error, and that there are many unit systems 
in use when dealing with opacities and radiative data!


To avoid as much confusion as possible, `exo_k` keeps track of the units
for the pressure and cross sections along with the data themselves (as attributes of the various
classes in the library). To do so, the library needs to
know the units used in the input files read. 

Available formats
-----------------

The currently supported formats are the following
(for the moment, the format is recognized by the library using the extension of the file).

.. list-table::
   :widths: 5 5 5 5 5
   :header-rows: 1

   * - Name 
     - Data type 
     - Ext.
     - k unit
     - P unit
   * - `ExoMol <http://exomol.com/data/data-types/opacity/>`_
     - k-table/x-sec
     - .h5
     - cm^2/molec
     - bar
   * - LMDZ
     - k-tables
     - 
     - cm^2/molec
     - mbar
   * - Nemesis 
     - k-tables
     - .kta
     - 10^-20 cm^2/molec
     - bar
   * - ARCIS
     - k-tables
     - .fits 
     - cm^2/molec
     - bar
   * - `ExoREM <https://lesia.obspm.fr/exorem/>`_
     - k-tables
     - .h5
     - cm^2/molec
     - bar
   * - `Exo_transmit <https://github.com/elizakempton/Exo_Transmit>`_
     - x-sec
     - .dat
     - m^2/molec
     - Pa 
   * - kspectrum
     - high-resolution
     - ASCII
     - m^-1
     - 
   * - `Helios-k <https://helios-k.readthedocs.io/en/latest/>`_
     - high-resolution
     - ASCII
     - cm^2/molec
     - 
   * - HITRAN
     - CIA coefficients
     - .cia
     - cm^5/molec
     - 
   * - `petit RADTRANS <https://www.dropbox.com/sh/w7sa20v8qp19b4d/AABKF0GsjghsYLJMUJXDgrHma?dl=0>`_
     - high-resolution
     - binary
     -
     -     

.. important::
    Self defining formats, like HDF5, make it possible to
    specify directly the units of the various physical quantities
    in the file as attributes.

    `Exo_k` takes full advantage of this feature and is able to read/write
    such hdf5 files with any combination of units. In addition, while being completely
    compatible with the Exomol format, our hdf5 format adds several variables and attributes
    to render the files more self-sufficient.

    Data output in one of the other formats discussed above
    is automatically converted back to the right units before writing. 



For high-resolution spectra, any file provided in ascii format with at least a column for wavenumber and one for opacity will work.


If you would like to see your favorite format up there, please go see this :ref:`section<contribute_format>`

In any case, the code can be forced to assume that the input file is using different units
by using the `file_kdata_unit` and `file_p_unit` keywords in the initialization methods for
cross sections and correlated-k tables (see the Getting Started section).

Then, the data can be converted to another unit system by simply calling the
`convert_kdata_unit` and `convert_p_unit` methods with the new desired units
(see the Getting Started section). 

Spectral units
--------------

The basic quantity used for the spectral dimension in the code is wavenumber.
As much as it pains me to write this, the community of infrared spectroscopists seems
very attached to inverse cemtimeters (cm^-1), so this is the unit used in `exo_k`. 

While all the computations in the library use wavenumbers, some functions and methods can be
given wavelengths as arguments.
Then, `wl` is usually in the name (as in `wl_range`). For the moment, whenever 
wavelengths are involved, microns are used.

Spectral information can be plotted using both wavenumbers or wavelength.

Units in the forward model
--------------------------

For numerical efficiency, the forward model integrated to `exo_k`
is written with a specific unit system (to avoid conversions).
We chose International (SI/MKS) units for cross-sections (m^2/molecule), pressures (Pa),
and CIA coefficients (m^-5).
We use cm^-1 for wavenumbers.
If the model is given radiative data in another unit system, it will complain. 

For this reason, if you want to use the forward model provided (or simply if you want to use a consistent set of units
without too much thinking about it), we recommend the use of
MKS units that can be enforced throughout the code with a single line:

>>> exo_k.Settings().set_mks(True)

After this line, all input data will automatically be converted to the MKS system when loaded. 


