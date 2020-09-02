# -*- coding: utf-8 -*-
"""
@author: jeremy leconte

Library of useful functions for users (only).

Functions here CANNOT be called into the library itself: importing this module in others would
lead to recursive import problems.
"""
import os
import numpy as np
from exo_k.ktable import Ktable
from exo_k.ktable5d import Ktable5d
from exo_k.xtable import Xtable
from exo_k.hires_spectrum import Hires_spectrum
from exo_k.util.cst import PI

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

def create_fname_grid(base_string, logpgrid=None, tgrid=None, xgrid=None,
        p_kw=None, t_kw=None, x_kw=None):
    """Creates a grid of filenames from an array of pressures, temperatures (
    and vmr if there is a variable gas).

    Parameters
    ----------
        base_string: str
            Generic name of the spectra files with specific keywords to be replaced 
            by the relevant numerical values
        logpgrid: Array
            Grid in log(pressure/Pa) of the input
        tgrid: Array
            Grid in temperature of the input
        xgrid: Array
            Input grid in vmr of the variable gas
        
    .. warning::
        The result of this function is much more predictable 
        if the values in the above arrays are given as integers. 
        If you want to use floats anyway, good luck. 

    Parameters
    ----------
        p_kw: str
        t_kw: str
        x_kw: str
            The pattern string that will be recognized as keywords between
            {} in base_string (See examples).

    Examples
    --------

        >>> logpgrid=[1,2]
        >>> tgrid=np.array([100.,150.,200.])
        >>> file_grid=exo_k.create_fname_grid('spectrum_CO2_1e{logp}Pa_{t}K.hdf5',
                  logpgrid=logpgrid,tgrid=tgrid,p_kw='logp',t_kw='t')                  
        array([['spectrum_CO2_1e1Pa_100K.hdf5', 'spectrum_CO2_1e1Pa_150K.hdf5',
        'spectrum_CO2_1e1Pa_200K.hdf5'],
        ['spectrum_CO2_1e2Pa_100K.hdf5', 'spectrum_CO2_1e2Pa_150K.hdf5',
        'spectrum_CO2_1e2Pa_200K.hdf5']], dtype='<U28')

    """
    logpgrid=np.array(logpgrid)
    tgrid=np.array(tgrid)
    res=[]
    if xgrid is None:
        for iP in range(logpgrid.size):
            for iT in range(tgrid.size):
                dict_opt={p_kw:str(logpgrid[iP]),t_kw:str(tgrid[iT])}
                fname=base_string.format(**dict_opt)
                res.append(fname)
        return np.array(res).reshape((logpgrid.size,tgrid.size))
    else:
        xgrid=np.array(xgrid)
        for iP in range(logpgrid.size):
            for iT in range(tgrid.size):
                for iX in range(xgrid.size):
                    dict_opt={p_kw:str(logpgrid[iP]), \
                        t_kw:str(tgrid[iT]),x_kw:str(xgrid[iX])}
                    fname=base_string.format(**dict_opt)
                    res.append(fname)
        return np.array(res).reshape((logpgrid.size,tgrid.size,xgrid.size))

def finalize_LMDZ_dir(corrkname, IRsize, VIsize):
    """Creates the right links for a LMDZ type directory to be read by the LMDZ generic GCM.

    You will need to create a proper Q.dat before using with the LMDZ GCM. 

    .. important::
        You must have run :func:`exo_k.ktable.Ktable.write_LMDZ` or
        :func:`exo_k.ktable5d.Ktable5d.write_LMDZ` for both of your IR and VI channels
        beforehand.

    Parameters
    ----------
        corrkname: str
            Path to the directory with the LMDZ ktable to finalize
        IRsize: int
            Number of IR spectral bins
        VIsize: int
            Number of VI spectral bins
    """
    newdir=os.path.join(corrkname,str(IRsize)+'x'+str(VIsize))
    try:
        os.mkdir(newdir)
    except FileExistsError:
        os.system('rm -rf '+newdir)
        os.mkdir(newdir)
    os.symlink('../IR'+str(IRsize)+'/corrk_gcm_IR.dat',os.path.join(newdir,'corrk_gcm_IR.dat'))
    os.symlink('../IR'+str(IRsize)+'/narrowbands_IR.in',os.path.join(newdir,'narrowbands_IR.in'))
    os.symlink('../VI'+str(VIsize)+'/corrk_gcm_VI.dat',os.path.join(newdir,'corrk_gcm_VI.dat'))
    os.symlink('../VI'+str(VIsize)+'/narrowbands_VI.in',os.path.join(newdir,'narrowbands_VI.in'))
    print('Everything went ok. Your ktable is in:',newdir)
    print("You'll probably need to add Q.dat before using it though!")

def mmr_to_number_density(mmr, gas_density, r_eff, condensate_density):
    """Converts a mass mixing ratio (mmr or q) in a number density of particles
    (in number per unit volume)

    Parameters
    ----------
        mmr: float or array
            Mass mixing ratio (in kg per kg of air)
        gas_density: float or array
            Density of the gas (in kg/m^3)
        r_eff: float or array
            Effective radius of the particles
        condensate_density: float or array
            Density of the constituent of the condensed particles (in kg/m^3)
    """
    particle_mass=4.*PI*r_eff**3*condensate_density/3.
    return mmr*gas_density/particle_mass
    