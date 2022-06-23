Exo_k
*****

About `Exo_k`
=============

`Exo_k` is a Python 3 based library to handle radiative opacities from various sources for atmospheric applications.
It enables you to:

* Interpolate efficiently and easily in correlated-k and cross section tables.
* Convert easily correlated-k and cross section tables from one format to another
  (hdf5, LMDZ GCM, Exomol, Nemesis, PetitCode, TauREx, ExoREM, ARCIS, etc.).
* Adapt precomputed correlated-k tables to your needs by changing:

  * the spectral and quadrature (g) grids,
  * the pressure/temperature grid.
* Create tables for a mix of gases using tables for individual gases.
* Create your own tables from high-resolution spectra (for example from K-spectrum, Helios-K, etc.).
* Use your data in an integrated radiative transfer framework to simulate planetary atmospheres.
  
On this website, you'll find a 'Getting Started' section that will show you how to do all that
with concrete examples that you can run on your own machine through the 
`tutorial jupyter notebook <https://forge.oasu.u-bordeaux.fr/jleconte/exo_k-public/-/blob/public/tutorial-exo_k.ipynb>`_
provided in the repository. Many important concepts and options are
presented along the way.

The API reference section also lists systematically all the classes
and methods available in the library and details the necessary inputs
and available options for each of them. This documentation is
searchable with the search bar in the top left corner. 

If you find this library useful in your reasearch, please acknowledge it by
citing `Leconte (2021) <https://ui.adsabs.harvard.edu/abs/2021A%26A...645A..20L/abstract>`_:

  * Spectral binning of precomputed correlated-k coefficients. **Astronomy and Astrophysics** 645. Leconte, J. 2021. doi:10.1051/0004-6361/202039040 

Enjoy!

`Jeremy Leconte <http://perso.astrophy.u-bordeaux.fr/~jleconte/>`_.

Recent releases
===============

v1.1.0 (August 2021): New scheme for the computation of atmospheric emission/transmission
to ensure an improved numerical accuracy. The variable names to instantiate atm objects have
changed accordingly (see tutorial). 

v1.0.2 (June 2021): Adds a few missing dependencies. Enables computation of thermal
emission spectra with scattering through the two-stream method (full documentation pending). 
Enables creating Xtables for a mix of gases (CIA can be added as well). Solves some issues
with the 2018 Hitran CIA format.

v1.0.1 (Jan 2021): Solves a binary/string conversion issue introduced by version 3 of h5py.
Enables linear interpolation in pressure (default is log). Enables creation of
empty tables to be filled later and extension of the spectral range of existing tables. 

v1.0.0 (Dec 2020): Finally our first official version. Creation of a
:ref:`'Neat examples'<neat_examples>` section with fully worked out use cases for `Exo_k`. 

v0.0.5 (Oct 2020): Ensures compatibility with latest Exomol correlated-k
and cross-section tables
(`Chubb et al. 2020 <https://ui.adsabs.harvard.edu/abs/2020arXiv200900687C/abstract>`_).


Doc Contents
============

.. toctree::
   :maxdepth: 2

   installation
   basic_principles
   units
   getting_started
   where_to_find_data
   practical_guide_to_atm
   exo_k_evol
   neat_examples
   API reference <autoapi/index>
   genindex
   how_to_contribute
   dev_and_doc

Other Links
===========

* Project homepage: http://perso.astrophy.u-bordeaux.fr/~jleconte/
* Code repository: https://forge.oasu.u-bordeaux.fr/jleconte/exo_k-public
* Documentation: http://perso.astrophy.u-bordeaux.fr/~jleconte/exo_k-doc/index.html
* Contact: jeremy.leconte at u-bordeaux.fr

Acknowledgements
================

This project has received funding from the European Research Council (ERC)
under the European Union's Horizon 2020 research and innovation programme
(grant agreement n° 679030/WHIPLASH).

The framework for this documentation has been developped by Aurelien Falco using Sphinx. 

