Exo_k
*****

About `Exo_k`
=============

`Exo_k` is a Python 3 based library to handle radiative opacities from various sources for atmospheric applications.
It enables you to:

* Interpolate efficiently and easily in correlated-k and cross section tables.
* Convert easily correlated-k and cross section tables from one format to another
  (pickle, hdf5, LMDZ GCM, Nemesis, PetitCode, TauREx, etc.).
* Adapt precomputed correlated-k tables to your needs by changing:

  * the resolution and quadrature (g) grid,
  * the pressure/temperature grid.
* Create tables for a mix of gases using tables for individual gases.
* Create your own tables from high-resolution spectra (for example from K-spectrum, Helios-K, etc.).
* Use your data in an integrated radiative transfer framework to simulate planetary atmospheres.
  
On this website, you'll find a 'Getting Started' section that will show you how to do all that
with concrete examples that you can run on your own machine through the tutorial-exo_k.ipynb
jupyter notebook provided in the repository. Many important concepts and options are
presented along the way.

The API reference section also lists systematically all the classes and methods available in the library and details the necessary inputs
and available options for each of them. This documentation is
searchable with the search bar in the top left corner. 

Enjoy!

J. Leconte


Doc Contents
============

.. toctree::
   :maxdepth: 2

   installation
   units
   getting_started
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

