Notes on units and formats
==========================

It is fair to say that units have often been a source of error, and that there are many unit systems 
in use when dealing with opacities and radiative data!

Keeping track of units
----------------------

To avoid as much confusion as possible, `exo_k` keeps track of the units
for the pressure and cross sections along with the data themselves (as attributes of the various
classes in the library). To do so, the library needs to
know the units used in the input files read. There are two possible cases:
  1. Self defining formats can directly specify the units in the file. 
     This is the case of all the hdf5 an pickle files created with `exo_k` and for new generation
     Exomol hdf5 files. 
  2. Some specific formats use fixed, known units. The ones we handle for the moment are:

     * Old pickle Exomol files (cm^2/molecule, bar)
     * LMDZ correlated-k tables (cm^2/molecule, mbar)
     * Exo_transmit (m^2/molecule, Pa). See https://github.com/elizakempton/Exo_Transmit or Kempton et al. (2016) for details
     * kspectrum (m^2/molecule)
     * HITRAN cia tables (cia coefficient in cm^5/molecule)

In any case, the code can be forced to assume that the input file is using different units
by using the `old_kdata_unit` and `old_p_unit` keywords in the initialization methods for
cross sections and correlated-k tables (see the Getting Started section).

Then, the data can be converted to another unit system by simply calling the
`convert_kdata_unit` and `convert_p_unit` methods with the new desired units
(see the Getting Started section). 

Data output by the code in a self defining format save the units used as well so they can be read
and used when this file is used elsewhere. Data output in one of the formats discussed above
is automatically converted back to the right units before writing. 

Units in the forward model
--------------------------

For numerical efficiency, however, the forward model integrated to `exo_k`
is written with a specific unit system (to avoid conversions).
We chose International (SI/MKS) units (m^2/molecule and Pa).
If the model is given radiative data in another unit system, it will complain. 

For this reason, you want to use the forward model provided, we recommand the use of
MKS units that can be enforced throughout the code with a single line::
    exo_k.Settings().set_mks(True)

After this line, all input data will automatically be converted to the MKS system. 


