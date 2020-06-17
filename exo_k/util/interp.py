# -*- coding: utf-8 -*-
"""
@author: jeremy leconte
Library of useful functions for interpolation
"""
import numpy as np
from numpy.polynomial.legendre import leggauss
import numba
from numba.typed import List
import astropy.units as u

@numba.njit(nogil=True,fastmath=True)
def bilinear_interpolation(z00, z10, z01, z11, x, y):
    """
    2D interpolation
    Applies linear interpolation across x and y between xmin,xmax and ymin,ymax

    Parameters
    ----------
        z00: array
            Array corresponding to xmin,ymin
        z10: array
            Array corresponding to xmax,ymin
        z01: array
            Array corresponding to xmin,ymax
        z11: array
            Array corresponding to xmax,ymax
        x: float
            weight on x coord
        y: float
            weight on y coord
    """
    xy = x*y
    res = np.zeros_like(z00)
    for i in range(z00.shape[0]):
        res[i] = (z11[i]-z01[i]+z00[i]-z10[i])*xy +(z01[i]-z00[i])*y +(z10[i]-z00[i])*x +z00[i]
    return res

@numba.njit(nogil=True,fastmath=True)
def linear_interpolation(z00, z10, x):
    """1D interpolation.

    Applies linear interpolation in x between xmin, xmax.

    Parameters
    ----------
        z00: array
            Array corresponding to xmin
        z10: array
            Array corresponding to xmax
        x: float
            weight on x coord
    """
    res = np.zeros_like(z00)
    for i in range(z00.shape[0]):
#        res[i] = z10[i]*x+z00[i]*x2
        tmp1=z00[i]
        tmp2=z10[i]
        tmp2-=tmp1
        tmp1+=tmp2*x
        res[i] = tmp1
    return res


def interp_ind_weights(x_to_interp,x_grid):
    """Finds the indices and weights to interpolate from a x_grid to a x_to_interp.
    """
    xmin=x_grid.min()
    xmax=x_grid.max()
    used_x=np.where(x_to_interp>xmax,xmax,x_to_interp)
    used_x=np.where(used_x<xmin,xmin,used_x)
    indices=np.searchsorted(x_grid,used_x)
    indices=np.where(indices==0,1,indices)
    return indices,(used_x-x_grid[indices-1])/(x_grid[indices]-x_grid[indices-1])

def rebin_ind_weights(old_bin_grid, new_bin_grid):
    """Computes the indices and weights to be used when
    rebinning an array of values f of size N in the bins delimited by old_bin_grid (of size N+1)
    onto new bins delimited by new_bin_grid (of size Nnew+1)
    Returns indices, weights to be used as follows to proceed with the rebinning

    >>> f_rebinned=[np.dot(f[indicestosum[ii]-1:indicestosum[ii+1]],final_weights[ii])
    for ii in range(Nnew)]
    """
#    new_bin_grid_used=np.where(new_bin_grid<old_bin_grid[0],old_bin_grid[0],new_bin_grid)
    indicestosum=np.searchsorted(old_bin_grid, new_bin_grid, side='right')
    if indicestosum[-1]==old_bin_grid.size : indicestosum[-1]-=1
    dg=old_bin_grid[1:]-old_bin_grid[:-1]
    final_weights=List()
    for ig in range(new_bin_grid.size-1):
        if indicestosum[ig+1]>indicestosum[ig]:
            weights=np.concatenate((old_bin_grid[[indicestosum[ig]]]-new_bin_grid[[ig]],\
                dg[indicestosum[ig]:indicestosum[ig+1]-1],  \
                new_bin_grid[[ig+1]]-old_bin_grid[[indicestosum[ig+1]-1]]))
        else:
            weights=np.array([new_bin_grid[ig+1]-new_bin_grid[ig]])
        weights=weights/np.sum(weights)
        final_weights.append(weights)
    return indicestosum, final_weights

@numba.njit()
def kdata_conv_loop(kdata1,kdata2,kdataconv,shape):
    """Computes the convolution of two kdata tables.

    Parameters
    ----------
        kdata1,kdata2 : arrays
            The two ktable.kdata tables to convolve.
        shape : array
            shape of the ktabke.kdata tables to convolve (last size is Ng).
        kdataconv : array
            Result table where the last dimension as a length equal to Ng^2.
    """
    Ng=shape[-1]
    for i in range(shape[0]):
     for j in range(shape[1]):
      for k in range(shape[2]):
       for l in range(Ng):
        for m in range(Ng):
            kdataconv[i,j,k,l*Ng+m]=kdata1[i,j,k,m]+kdata2[i,j,k,l]

@numba.njit(nogil=True,fastmath=True)
def kdata_conv_loop_profile(kdata1,kdata2,kdataconv,Nlay,Nw,Ng):
    """Computes the convolution of two atmospheric kdata profiles.

    Nothing is returned. kdataconv is changed in place. 

    Parameters
    ----------
        kdata1,kdata2 : arrays
            The two ktable.kdata tables to convolve.
        kdataconv : array
            Result table where the last dimension as a length equal to Ng^2.
        Nlay : int
            Number of atmospheric layers
        Nw : int
            Number of wavenumber points
        Ng : int
            Number of g-points
    """
    for i in range(Nlay):
        for j in range(Nw):
            for l in range(Ng):
                for m in range(Ng):
                    kdataconv[i,j,l*Ng+m]=kdata1[i,j,m]+kdata2[i,j,l]

@numba.njit(fastmath=True)
def rebin(f_fine,fine_grid,coarse_grid):
    """Computes the binned version of a function on a coarser grid.
    The two grids do not need to have the same boundaries.

    Parameters
    ----------
        f_fine: array
            Function to be rebinned, given on the fine_grid.
        fine_grid: array
            The high resolution grid edges we start with.
        coarse_grid: array
            The coarser resolution grid edges inside which we want to bin the function f.
       """
    indicestosum=np.searchsorted(fine_grid,coarse_grid,side='right')
    Nfine=fine_grid.size
    if indicestosum[-1]==Nfine : 
        indicestosum=np.where(indicestosum<Nfine,indicestosum,Nfine-1)
    dgrid=np.diff(fine_grid)
    Nn=coarse_grid.size-1
    f_coarse=np.zeros(Nn)
    ifine=indicestosum[0]-1
    for ii in range(Nn):
      if indicestosum[ii+1]>indicestosum[ii]:
        tmp_w=fine_grid[ifine+1]-coarse_grid[ii]
        tmp_f=f_fine[ifine]*tmp_w
        for iifine in range(indicestosum[ii],indicestosum[ii+1]-1):
          tmp_f+=dgrid[iifine]*f_fine[iifine]
          tmp_w+=dgrid[iifine]
        ifine=indicestosum[ii+1]-1
        tmp_dw=(coarse_grid[ii+1]-fine_grid[ifine])
        #print(ifine,tmp_dw,f_fine[ifine])
        tmp_f+=f_fine[ifine]*tmp_dw
        tmp_w+=tmp_dw
        f_coarse[ii]=tmp_f/tmp_w
      else:
        f_coarse[ii]=f_fine[ifine]
    return f_coarse

#@numba.njit(nogil=True,fastmath=True)
def RandOverlap_2_kdata_prof(Nlay, Nw, Ng, kdata1, kdata2, weights, ggrid):
    """Function to randomely mix the opacities of 2 species in an atmospheric profile.

    Parameters
    ----------
        Nlay, Nw, Ng: int
            Number of layers, spectral bins, and gauss points. 
        kdata1, kdata2: arrays of size (Nlay, Nw, Ng)
            vmr weighted cross-sections for the two species.
        weights: array
            gauss weights.
        ggrid: array
            g-points.
    
    Returns
    -------
        array
            k-coefficient array of the mix over the atmospheric column.
    """
    kdataconv=np.zeros((Nlay,Nw,Ng**2))
    weightsconv=np.zeros(Ng**2)
    newkdata=np.zeros((Nlay,Nw,Ng))

    for jj in range(Ng):
        for ii in range(Ng):
            weightsconv[ii*Ng+jj]=weights[jj]*weights[ii]

    kdata_conv_loop_profile(kdata1,kdata2,kdataconv,Nlay,Nw,Ng)

    for ilay in range(Nlay):
        for iW in range(Nw):
            tmp=kdataconv[ilay,iW,:]
            indices=np.argsort(tmp)
            kdatasort=tmp[indices]
            weightssort=weightsconv[indices]
            newggrid=np.cumsum(weightssort)
            #ind=np.searchsorted(newggrid,ggrid,side='left')
            #newkdata[ilay,iW,:]=kdatasort[ind]
            newkdata[ilay,iW,:]=np.interp(ggrid,newggrid,kdatasort)
    return newkdata

def unit_convert(quantity,unit_file='unspecified',unit_in='unspecified',unit_out='unspecified'):
    """Chooses the final unit to use and the conversion factor to apply.

    Parameters
    ----------
        quantity: str
            The name of the pysical quantity handled for potential error messages
        unit_file/unit_in/unit_out: str

            Respectively:

              * String with the unit found in (or assumed from the format of) the initial data,
              * The unit we think the initial data are in if unit_file is 'unspecified'
                or (we believe) wrong,
              * The unit we want to convert to. If unspecified, we do not convert.
    Returns
    -------
        unit_to_write: str
            Resulting unit.
        conversion_factor: float
            A multiplicating factor for the data to proceed to the conversion.
    """
    if unit_in!='unspecified':
      if((unit_file!='unspecified') and (unit_file!=unit_in)):
        print('Be careful, you are assuming that '+quantity+' is '+unit_in)
        print('but the input file says that it is '+unit_file+'. The former will be used')
      starting_unit=unit_in
    else:
      if unit_file=='unspecified':
        raise NotImplementedError("""I could not find the {quantity} used in the file.
            So you should tell me what unit you think is used by specifying 
            the keyword argument: file_{quantity}= a unit recognized by the astropy.units library.
            If you want to convert to another unit, add the {quantity} keyword
            with the desired unit.""".format(quantity=quantity))
      else:
        starting_unit=unit_file
    if unit_out=='unspecified':
        unit_to_write=starting_unit
    else:
        unit_to_write=unit_out
    return unit_to_write,u.Unit(starting_unit).to(u.Unit(unit_to_write))

def rm_molec(unit_name):
    """Removes "/molecule" or "/molec" from a unit string.

    Parameters
    ----------
        unit_name: str
            String to be changed.
        
    Returns
    -------
        str
            The unit name without the ending "/molecule" or "/molec"
    """
    return unit_name.replace('/molecule','').replace('/molec','')

@numba.njit
def is_sorted(a):
    """Finds out if an array is sorted. Returns True if it is.
    """
    for i in range(a.size-1):
         if a[i+1] < a[i] :
               return False
    return True


def gauss_legendre(order):
    """Computes the weights and abscissa for a Gauss Legendre quadrature of order `order`

    Parameters
    ----------
        order: int
            Order of the quadrature wanted.

    Returns
    -------
        weights: array(order)
            Weights to be used in the quadrature.
        ggrid: array(order)
            Abscissa to be used for the quadrature.
        gedges: array(order+1)
            Cumulative sum of the weights. Goes from 0 to 1.
       """ 
    ggrid,weights=leggauss(order)
    weights=weights/2.
    ggrid=(ggrid+1.)/2.
    gedges=np.insert(np.cumsum(weights),0,0.)
    return weights,ggrid,gedges

@numba.njit
def spectrum_to_kdist(k_hr,wn_hr,dwn_hr,wnedges,ggrid):
    """Creates the k-distribution from a single high resolution spectrum
    
    Parameters
    ----------
        k_hr : array
            Spectrum
        wn_hr : array
            Wavenumber grid
        dwn_hr : array
            Width of the high resolution wavenumber bins.
        wnedges : array
            Lower resolution wavenumber bin edges inside which the k-dist will be computed.
        ggrid : array
            Grid of g-point abscissas

    Returns
    -------
        kdata : array
            k coefficients for each (Wn, g) bin.
    """
    pos=np.searchsorted(wn_hr,wnedges)
    #print(pos)
    kdata=np.zeros((wnedges.size-1,ggrid.size))
    for ib in range(wnedges.size-1):
        if pos[ib+1]==pos[ib]:
            continue
            # JL20 why did I write that ??
            # for the moment, this leaves zeros but we should probably make the code stop
        tmp_k=k_hr[pos[ib]:pos[ib+1]]
        totwg=np.sum(dwn_hr[pos[ib]:pos[ib+1]])
        wg=dwn_hr[pos[ib]:pos[ib+1]]/totwg
        sort_ind=np.argsort(tmp_k)
        newgfine=np.cumsum(wg[sort_ind])
        # choosing the line below or the two after actually makes a significant difference
        kdata[ib]=np.interp(ggrid,newgfine,tmp_k[sort_ind])
        #g_pos=np.searchsorted(newgfine,ggrid,side='right')
        #kdata[ib]=tmp_k[sort_ind[g_pos]] 
        # an interpolation between tmp_k[sort_ind[g_pos]] and tmp_k[sort_ind[g_pos-1]] 
        # would probably be cleaner.
    return kdata



@numba.njit()
def bin_down_corrk_numba(newshape, kdata, old_ggrid, new_ggrid, gedges, indicestosum, \
        wngrid_filter, wn_weigths, num, use_rebin):
    """bins down a kcoefficient table (see :func:`~exo_k.ktable.Ktable.bin_down` for details)
    
    Parameters
    ----------
        newshape : array
            Shape of kdata.
        kdata : array
            table to bin down.
        old_ggrid : array
            Grid of old g-point abscissas for kdata.
        new_ggrid : array
            New g-points for the binned-down k-coefficients.
        gedges : array
            Cumulative sum of the weights. Goes from 0 to 1.
            Used only if use_rebin=True
        indicestosum: list of lists
            Indices of wavenumber bins to be used for the averaging
        wngrid_filter: array
            Indices of the new table where there will actually be data (zero elsewhere)
        wn_weigths: list of lists
            Weights to be used for the averaging (same size as indicestosum)
        num: int
            Number of points to fine sample the g function in log-k space
        use_rebin: boolean
            Whether to use rebin or interp method. 
    """
    ggrid0to1=np.copy(old_ggrid)
    ggrid0to1[0]=0.
    ggrid0to1[-1]=1.
    newkdata=np.zeros(newshape, dtype=np.float64)
    for iW,newiW in enumerate(wngrid_filter[0:-1]):
        tmp_dwn=wn_weigths[iW]
        for iP in range(newshape[0]):
            for iT in range(newshape[1]):            
                tmp_logk=np.log(kdata[iP,iT,indicestosum[iW]-1:indicestosum[iW+1]])
                #if iP==0 and iT==0: print(tmp_logk)
                logk_min=np.amin(tmp_logk[:,0])
                logk_max=np.amax(tmp_logk[:,-1])
                if logk_min==logk_max:
                    newkdata[iP,iT,newiW,:]=np.ones(newshape[-1])*np.exp(logk_max)
                else:
                    logk_max=logk_max+(logk_max-logk_min)/(num-3.)
                    logk_min=logk_min-(logk_max-logk_min)/(num-3.)
                    logkgrid=np.linspace(logk_min,logk_max,num)
                    newg=np.zeros(logkgrid.size)
                    for ii in range(tmp_logk.shape[0]):
                        newg+=np.interp(logkgrid,tmp_logk[ii],ggrid0to1)*tmp_dwn[ii]
                        #newg+=np.interp(logkgrid,tmp_logk[ii],ggrid0to1,
                        #    left=0., right=1.)*tmp_dwn[ii]
                    if use_rebin:
                        newkdata[iP,iT,newiW,:]=rebin(np.exp(logkgrid), newg, gedges)
                    else:
                        newkdata[iP,iT,newiW,:]=np.interp(new_ggrid, newg, np.exp(logkgrid))
    return newkdata

## OLD FUNCTIONS

@numba.njit()
def kdata_conv(kdata1,kdata2,kdataconv,Ng):
    """Deprecated.
    
    Performs the convolution of the kdata over g-space.
    """
    for jj in range(Ng):
        for ii in range(Ng):
            kdataconv[:,:,:,ii*Ng+jj]=kdata1[:,:,:,jj]+kdata2[:,:,:,ii]

@numba.njit()
def kdata_conv_loop_bad(kdata1,kdata2,kdataconv,shape):
    """Deprecated. 
    
    Performs the convolution of the kdata over g-space.
    But the loop is in bad order. Keep for test. 
    """
    Ng=shape[-1]
    for jj in range(Ng):
     for ii in range(Ng):
      for k in range(shape[2]):
       for l in range(shape[1]):
        for m in range(shape[0]):
            kdataconv[m,l,k,ii*Ng+jj]=kdata1[m,l,k,jj]+kdata2[m,l,k,ii]


#@numba.njit(fastmath=True)
def g_interp(logk,Nw,Ng,Num,ggrid,wl_weights,indices):
    """Interpolates logk on a new g-grid. 
    """
    #write = False
    #if write : print(logk);print(Nw,Ng,Num);print(ggrid);print(wl_weights);print(indices)
    lkmin=logk[0,0]
    for iW in range(1,Nw):
        if logk[iW,0] < lkmin : lkmin = logk[iW,0]
    lkmax=logk[0,-1]
    for iW in range(1,Nw):
        if logk[iW,-1] > lkmax : lkmax = logk[iW,-1]
    dlk=(lkmax-lkmin)/(Num-3)
    logkgrid=np.zeros(Num)
    Ovgedges=1./np.diff(ggrid)
    newg=np.zeros(Num)
    lk=lkmin-dlk
    logkgrid[0]=lk
    for iN in range(1,Num):
        lk+=dlk
        logkgrid[iN]=lk
        tmpg=0.
        for iW in range(Nw):
            ind=indices[iW]
            #if write : print(ind,iW,lk,logk[iW,ind])
            if lk < logk[iW,ind]:
                #if write : print('lk<logk')
                if ind!=0 :
                    #if write: print('ind!=0')
                    #if write: print( \
                    # (ggrid[ind]*(lk-logk[iW,ind-1])+ggrid[ind-1]*(logk[iW,ind]-lk))\
                    #   *Ovgedges[ind-1]*wl_weights[iW])
                    tmpg+=(ggrid[ind]*(lk-logk[iW,ind-1])+ggrid[ind-1]*(logk[iW,ind]-lk))  \
                        *Ovgedges[ind-1]*wl_weights[iW]
            else:
                #if write : print('lk>logk')
                if ind!=Ng-1:
                    #if write : print('ind!=Ng-1')
                    ind+=1
                    indices[iW]=ind
                    #if write : print(  \
                    # (ggrid[ind]*(lk-logk[iW,ind-1])+ggrid[ind-1]*(logk[iW,ind]-lk)) \
                    # *Ovgedges[ind-1]*wl_weights[iW])
                    tmpg+=(ggrid[ind]*(lk-logk[iW,ind-1])+ggrid[ind-1]*(logk[iW,ind]-lk)) \
                        *Ovgedges[ind-1]*wl_weights[iW]
                else:
                    #if write : print('ind=Ng-1')
                    tmpg+=wl_weights[iW]
        newg[iN]=tmpg
    #if write : print(newg);print(logkgrid)

    return newg,logkgrid

