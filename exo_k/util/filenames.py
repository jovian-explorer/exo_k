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

def create_fname_grid(base_string,logpgrid=None, tgrid=None, xgrid=None,
        p_kw=None, t_kw=None, x_kw=None):
    """Creates a grid of filenames from an array of pressures, temperatures,
    and vmr if there is a variable gas.
    Parameters:
        base_string: str
            Generic name of the spectra files with specific keywords to be replaced 
            by the relevant numerical values
        logpgrid: Array
            Grid in log(pressure/Pa) of the input
        tgrid: Array
            Grid in temperature of the input
        xgrid: Array
            Input grid in vmr of the variable gas
        The values in the above arrays must be convertible in integers
        p_kw: str
        t_kw: str
        x_kw: str
            Input grid in vmr of the variable gas

    Example:
        ```
        logpgrid=[1,2]
        tgrid=np.array([100.,150.,200.])
        print(logpgrid,tgrid)
        file_grid=exo_k.create_fname_grid('spectrum_CO2_1e{logp}Pa_{t}K.hdf5',
                  logpgrid=logpgrid,tgrid=tgrid,p_kw='logp',t_kw='t')
        ```
        Results in 
        ```
        array([['spectrum_CO2_1e1Pa_100K.hdf5', 'spectrum_CO2_1e1Pa_150K.hdf5',
        'spectrum_CO2_1e1Pa_200K.hdf5'],
       ['spectrum_CO2_1e2Pa_100K.hdf5', 'spectrum_CO2_1e2Pa_150K.hdf5',
        'spectrum_CO2_1e2Pa_200K.hdf5']], dtype='<U28')
        ```
    """
    logpgrid=np.array(logpgrid)
    tgrid=np.array(tgrid)
    res=[]
    if xgrid is None:
        for iP in range(logpgrid.size):
            for iT in range(tgrid.size):
                dict_opt={p_kw:str(int(logpgrid[iP])),t_kw:str(int(tgrid[iT]))}
                fname=base_string.format(**dict_opt)
                res.append(fname)
        return np.array(res).reshape((logpgrid.size,tgrid.size))
    else:
        xgrid=np.array(xgrid)
        for iP in range(logpgrid.size):
            for iT in range(tgrid.size):
                for iX in range(xgrid.size):
                    dict_opt={p_kw:str(int(logpgrid[iP])), \
                        t_kw:str(int(tgrid[iT])),x_kw:str(int(xgrid[iX]))}
                    fname=base_string.format(**dict_opt)
                    res.append(fname)
        return np.array(res).reshape((logpgrid.size,tgrid.size,xgrid.size))

def finalize_LMDZ_dir(corrkname,IRsize,VIsize):
    """Creates the right links for a LMDZ type directory to be read by the LMDZ generic GCM.
      => you must have run write_LMDZcorrk for your IR and VI channels.
    """
    import os
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
    print("You'll probably need to add Q.dat before using it thouh!")

def convert_kspectrum_to_hdf5(file_in,file_out,skiprows=0):
    """Converts kspectrum like spectra to hdf5 format for speed and space
    """
    wn_hr,k_hr=np.loadtxt(file_in,skiprows=skiprows,unpack=True) 
    f = h5py.File(file_out, 'w')
    f.create_dataset("wns", data=wn_hr,compression="gzip")
    f.create_dataset("k", data=k_hr,compression="gzip")
    f.close()    

def convert_exo_transmit_to_hdf5(file_in,file_out,mol='unspecified'):
    """Converts exo_transmit like spectra to hdf5 format for speed and space
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
        fullfilename=file_out+'.hdf5'
    else:
        fullfilename=file_out
    compression="gzip"
    f = h5py.File(fullfilename, 'w')
    f.attrs["mol_name"] = mol
    f.create_dataset("p", data=pgrid,compression=compression)
    f["p"].attrs["units"] = 'Pa'
    f.create_dataset("t", data=tgrid,compression=compression)
    f.create_dataset("xsecarr", data=kdata,compression=compression)
    f["xsecarr"].attrs["units"] = 'm^2'
    f.create_dataset("bin_edges", data=wns,compression=compression)
    f.close()    

