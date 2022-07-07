# -*- coding: utf-8 -*-
"""
@author: jeremy leconte

Library of useful functions for handling filenames
"""
import numpy as np
import h5py

class EndOfFile(Exception):
    """Error for an end of file
    """

def select_kwargs(kwargs, filter_keys_list=[]):
    """Function to select only some keyword arguments from a
    kwargs dict to pass to a function afterward.
    
    Parameters
    ----------
        kwargs: dict
            Dictionary of keyword arguments and values.
        filter_keys_list: list of str
            Names of the keys to select.

    Returns
    -------
        filtered_kwargs: dict
            A dictionary with only the selected keys (if they were
            present in kwargs).
        
    Examples
    --------
        >>> def func(allo=None):
        >>>     print(allo)
        >>>
        >>> def bad_global_func(**kwargs):
        >>>     print(kwargs)
        >>>     func(**kwargs)
        >>>
        >>> def good_global_func(**kwargs):
        >>>     print(kwargs)
        >>>     func(**select_kwargs(kwargs,['allo']))

        >>> bad_global_func(allo=3.,yeah=None)        
        {'allo': 3.0, 'yeah': None}
        TypeError: func() got an unexpected keyword argument 'yeah'

        >>> good_global_func(allo=3.,yeah=None)        
        {'allo': 3.0, 'yeah': None}
        3.0

    """
    filtered_kwargs=dict()
    for key in filter_keys_list:
        if key in kwargs.keys():
            filtered_kwargs[key]=kwargs[key]
    return filtered_kwargs

def create_fname_grid_Kspectrum_LMDZ(Np, Nt, Nx=None, suffix='', nb_digit=3):
    """Creates a grid of filenames consistent with Kspectrum/LMDZ
    from the number of pressure, temperatures (, and vmr) points (respectively) in the grid.

    Parameters
    ----------
        Np, Nt: int
            Number of Pressure and temperature points.
        Nx: int, optional
            Number of vmr points if there is a variable gas
        suffix: str, optional
            String to add behind the 'k001'
        nb_digit: int, optional
            Total number of digit including the zeros ('k001' is 3)
    
    Returns
    -------
        list of str
            List of the files in the right order and formating to be 
            given to :func:`exo_k.ktable.Ktable.hires_to_ktable` or 
            :func:`exo_k.xtable.Xtable.hires_to_xtable`.

    Examples
    --------

        >>> exo_k.create_fname_grid_Kspectrum_LMDZ(2,3)                  
        array([['k001', 'k002', 'k003'],
        ['k004', 'k005', 'k006']], dtype='<U4')

        >>> exo_k.create_fname_grid_Kspectrum_LMDZ(2,3,suffix='.h5')                  
        array([['k001.h5', 'k002.h5', 'k003.h5'],
        ['k004.h5', 'k005.h5', 'k006.h5']], dtype='<U7')

    """    
    res=[]
    ii=1
    if Nx is None:
        for _ in range(Np):
            for _ in range(Nt):
                res.append('k'+str(ii).zfill(nb_digit)+suffix)
                ii+=1
        return np.array(res).reshape((Np,Nt))
    else:
        for _ in range(Nx):
            for _ in range(Np):
                for _ in range(Nt):
                    res.append('k'+str(ii).zfill(nb_digit)+suffix)
                    ii+=1
        res=np.array(res).reshape((Nx,Np,Nt))
        return np.transpose(res,(2,1,0))

def convert_exo_transmit_to_hdf5(file_in, file_out, mol='unspecified'):
    """Converts exo_transmit like spectra to hdf5 format for speed and space.

    Parameters
    ----------
        file_in: str
            Initial exo_transmit filename.
        file_out: str
            Name of the final hdf5 file to be created.
    """
    tmp_wlgrid=[]
    tmp_kdata=[]
    with open(file_in, 'r') as file:
        tmp = file.readline().split()
        tgrid=np.array([float(ii) for ii in tmp])
        tmp = file.readline().split()
        pgrid=np.array([float(ii) for ii in tmp])
        Np=pgrid.size
        Nt=tgrid.size
        while True:
            line=file.readline()
            if line is None or line=='': break
            tmp_wlgrid.append(float(line.split()[0])) # wavelength in m
            tmp_kdata.append([])
            for _ in range(Np):
                tmp_kdata[-1].append(file.readline().split()[1:])
    tmp_wlgrid=np.array(tmp_wlgrid)
    Nw=tmp_wlgrid.size
    kdata=np.zeros((Np,Nt,Nw))
    for iP in range(Np):
        for iT in range(Nt):
            for iW in range(Nw):
                kdata[iP,iT,iW]=tmp_kdata[Nw-iW-1][iP][iT]
    print(kdata.shape)
    wns=0.01/tmp_wlgrid[::-1]
    
    if not file_out.lower().endswith(('.hdf5', '.h5')):
        fullfilename=file_out+'.h5'
    else:
        fullfilename=file_out
    compression="gzip"
    with h5py.File(fullfilename, 'w') as f:
        f.attrs["mol_name"] = mol
        f.create_dataset("p", data=pgrid,compression=compression)
        f["p"].attrs["units"] = 'Pa'
        f.create_dataset("t", data=tgrid,compression=compression)
        f.create_dataset("xsecarr", data=kdata,compression=compression)
        f["xsecarr"].attrs["units"] = 'm^2'
        f.create_dataset("bin_edges", data=wns,compression=compression)

def _read_array(file, Nvalue, N_per_line=5, Nline=None, revert=False):
    """Reads an array in a .dat filestream. 
    Assumes that the arrays are arranged N_per_line values per line.

    Parameters
    ----------
        file: file stream
            File to be read.
        Nvalue: int
            Number of values to be read. 
        N_per_line: int
            Number of values assumed per line.
        revert: boolean (optional)
            Whether or not to revert the array (because some arrays are
            in decreasing order).
        Nline: int (optional)
            Slightly speeds up the process if the user specifies
            the number of lines that will need to be read. 

    Returns
    -------
        Array
            A numpy array with the values.
    """
    
    if Nline is None:
        Nline=Nvalue//N_per_line
        if Nvalue%N_per_line != 0:
            Nline+=1
    new_array=[]
    for _ in range(Nline):
        line=[float(i) for i in file.readline().split()]
        new_array+=line
    if revert :
        new_array=np.array(new_array)[::-1]
    else:
        new_array=np.array(new_array)
    #print(new_array)
    return new_array

def _read_exorem_k_array(file, Nline):
    """Reads an array in a .dat filestream. 
    Assumes that the arrays are arranged N_per_line values per line.

    Parameters
    ----------
        file: file stream
            File to be read.
        Nline: int
            Number of values/lines to be read. 

    Returns
    -------
        Array
            A numpy array with the values.
    """
    
    new_array=[]
    for _ in range(Nline):
        srt=file.readline()
        try:
            new_array.append(float(srt))
        except ValueError:
            i = srt.index('-')
            srt = srt[:i] + 'e' + srt[i:]
            new_array.append(float(srt))
    new_array=np.array(new_array)
    return new_array



####### OLD FUNCTIONS ###################
def convert_old_kspectrum_to_hdf5(file_in, file_out, skiprows=0):
    """Deprecated.
    
    Converts kspectrum like spectra to hdf5 format for speed and space

    Parameters
    ----------
        file_in: str
            Initial kspectrum filename.
        file_out: str
            Name of the final hdf5 file to be created.
    """
    wn_hr,k_hr=np.loadtxt(file_in,skiprows=skiprows,unpack=True) 
    with h5py.File(file_out, 'w') as f:
        f.create_dataset("wns", data=wn_hr,compression="gzip")
        f.create_dataset("k", data=k_hr,compression="gzip")
