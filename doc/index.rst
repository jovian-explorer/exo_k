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

The API reference section also lists systematically all the classes and methods available in the library and details the necessary inputs
and available options for each of them. This documentation is
searchable with the search bar in the top left corner. 

Enjoy!

J. Leconte

Recent releases
===============

v1.0.1 (Jan 2021): Solves a binary/string conversion issue on some platforms.
Enables linear interpolation in pressure (default is log). Enables creation of
empty tables to be filled later and spectral extension of existing tables. 

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
   neat_examples
   API reference <autoapi/index>
   genindex
   how_to_contribute
   build_doc

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

