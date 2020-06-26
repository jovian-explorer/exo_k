How to contribute
=================

Suggest improvements
--------------------

`Exo_k` is still evolving.
Whether you have noticed a bug that needs fixing, or thought of a feature that you think should be added,
the best way to make it happen for the moment is to send me an email at jeremy dot leconte at u-bordeaux.fr.

In order to speed up the process, please follow these steps:

 * In case of a bug, describe the problem and try to provide a minimal (not) working example
   exhibiting the unintended behavior. As much as possible, please provide a link to a minimal
   dataset if needed to reproduce the bug. 
 * Please provide the number or commit SHA of the version you have been using.
 * Although it might sometimes seem obvious, please describe the intended behavior.
 * If you have come up with a solution, please provide it so it can be tested and implemented. 

.. _contribute_format:

What about your data format?
----------------------------

For the specific case where you would like to see a particular format implemented (e.g. linked to a specific code
generating or using radiative data), please provide a complete description of the format in question, including
how to read/write the following attributes (some of these may be relevant to some specific data types only):

 * The arrays containing the temperature and pressure grids.
 * The arrays containing the wavenumber (or wavelength) of the centers (and edges for k-coefficients) of the spectral bins.
 * The array containing the opacity data
   (with a description of the opacity type: cross section per molecule or per unit mass, absorption coefficient, etc.).
 * The units used for all these parameters.
 * The name or identifier of the molecule described.
 * The quadrature weights (and possibly abcissa).

As nobody understands your data better than you, the easiest way to do this
and make sure your data will be read/written properly is to provide a Python functions with the following
format:

>>> def read_my_format(filename):
>>>     'Do what you gotta do'
>>>     return pgrid, tgrid, wn_grid, ggrid, weights, opacity_data, p_unit, wn_unit, opacity_unit


>>> def write_my_format(filename, pgrid, tgrid, wn_grid, ggrid, weights, opacity_data):
>>>     'Creates the relevant file'
>>>     return None
