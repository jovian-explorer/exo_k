# -*- coding: utf-8 -*-
"""
@author: jeremy leconte

Library of useful functions for users (only).

Functions here CANNOT be called into the library itself: importing this module in others would
lead to recursion problems.
"""
import numpy as np
from exo_k.ktable import Ktable
from exo_k.ktable5d import Ktable5d
from exo_k.xtable import Xtable
from exo_k.util.hires_spectrum import Hires_spectrum

def hires_to_ktable(filename_grid=None, xgrid=None, **kwargs):
    """Emulates :func:`exo_k.ktable.Ktable.hires_to_ktable`
    and :func:`exo_k.ktable5d.Ktable5d.hires_to_ktable` as functions and not methods
    so that the user can call them without creating a Ktable first.

    See those methods for details on the available arguments and options.
    """
    if filename_grid is not None:
        if len(np.array(filename_grid).shape) == 3 and xgrid is None:
            raise RuntimeError("""From the shape of filename_grid
                it seems you want have an xgrid dimension, but xgrid is None!
                """)    
    if xgrid is not None: # creates a Ktable5d
        if filename_grid is not None:
            if len(np.array(filename_grid).shape) != 3:
                raise RuntimeError("""You provided a xgrid of vmrs, 
                    but the shape of filename_grid is not 
                    of the type (Np, Nt, Nx)!
                    """)
        res=Ktable5d()
        res.hires_to_ktable(filename_grid=filename_grid, xgrid=xgrid, **kwargs)
    else:
        res=Ktable()
        res.hires_to_ktable(filename_grid=filename_grid, **kwargs)
    return res

def hires_to_xtable(**kwargs):
    """Emulates :func:`exo_k.xtable.Xtable.hires_to_xtable`
    as a function and not a method
    so that the user can call it without creating an Xtable first.

    See :func:`~exo_k.xtable.Xtable.hires_to_xtable` for details
    on the available arguments and options.
    """
    res=Xtable()
    res.hires_to_xtable(**kwargs)
    return res

def convert_kspectrum_to_hdf5(file_in, file_out=None, **kwargs):
    """Converts kspectrum like spectra to hdf5 format for speed and space.
    Helper function. Real work done in :class:`~exo_k.util.hires_spectrum.Hires_spectrum`
    __init__ funtion.

    Go there to see all the available arguments and options.

    Parameters
    ----------
        file_in: str
            Initial kspectrum filename.
        file_out: str
            Name of the final hdf5 file to be created. If not provided,
            'file_in.h5' will be used. 
    """
    tmp=Hires_spectrum(file_in, **kwargs)
    if file_out is None: file_out=file_in
    tmp.write_hdf5(file_out)

