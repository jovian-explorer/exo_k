# -*- coding: utf-8 -*-
"""
Created on Jun 20 2019

@author: jeremy leconte

Library of useful functions for radiative transfer calculations
"""
import numpy as np
import numba

#def solve_2stream_nu(NG, source_nu, dtau_nu, omega0_nu=0., g_asym_nu=0.,
#                mu0=0.5, flux_top_dw_nu=0., alb_surf=0., verbose=False):
#    """Deals with the spectral axis
#    """
#    pass

@numba.jit(nopython=True,fastmath=True)
def solve_2stream_nu(NG, source_nu,dtau_nu,omega0_nu,g_asym_nu,mu0=0.5,
        flux_top_dw_nu=0., alb_surf=0.):
#      
#      use radinc_h
#      use radcommon_h, only: planckir, tlimit,sigma, gweight
#      use comcstfi_mod, only: pi
#      
#      implicit none
#      #      
#      real*8 TLEV(L_LEVELS), PLEV(L_LEVELS)
#      real*8 TAUCUMI(L_LEVELS,L_NSPECTI,L_NGAUSS)
#      real*8 FMNETI(L_NLAYRAD)
#      real*8 WNOI(L_NSPECTI), DWNI(L_NSPECTI)
#      real*8 dtau_nu(L_NLAYRAD,L_NSPECTI,L_NGAUSS)
#      real*8 FMUPI(L_NLEVRAD), FMDI(L_NLEVRAD)
#      real*8 COSBI(L_NLAYRAD,L_NSPECTI,L_NGAUSS)
#      real*8 omega0_nu(L_NLAYRAD,L_NSPECTI,L_NGAUSS)
#      real*8 NFLUXTOPI
#      real*8 NFLUXTOPI_nu(L_NSPECTI)
#      real*8 fluxupi_nu(L_NLAYRAD,L_NSPECTI)
#      real*8 FTOPUP
#      
#      real*8 mu0, RSFI, TSURF, BSURF, TTOP, BTOP, TAUTOP
#      real*8 PLANCK, PLTOP
#      real*8 fluxupi(L_NLAYRAD), fluxdni(L_NLAYRAD)
#      
#      real*8 PLANCKSUM,PLANCKREF
#      
## AB : variables for interpolation
#      REAL*8 C1
#      REAL*8 C2
#      REAL*8 P1
#      
##======================================================================#
      
    Nlay,Nw,Ng   = dtau_nu.shape
    FLUXUPI_nu   = np.zeros((Nlay+1, Nw, Ng))
    FLUXDWI_nu   = np.zeros((Nlay+1, Nw, Ng))
    FNETI_nu     = np.zeros((Nlay+1, Nw, Ng))
              
# WE NOW ENTER A MAJOR LOOP OVER SPECTRAL INTERVALS IN THE INFRARED
# TO CALCULATE THE NET FLUX IN EACH SPECTRAL INTERVAL
            
    for NW in range(Nw):         
        for iG in range(Ng):            
            
# SET UP THE UPPER AND LOWER BOUNDARY CONDITIONS ON THE IR
# CALCULATE THE DOWNWELLING RADIATION AT THE TOP OF THE MODEL
# OR THE TOP LAYER WILL COOL TO SPACE UNPHYSICALLY
            
#            TAUTOP = dtau_nu[1,NW,iG]*PLEV[2]/(PLEV[4]-PLEV[2])
#JL21 simple boundary condition for now to compare with toon
#            TAUTOP = TAUCUMI[0,NW,iG]
#            BTOP   = (1.e0-np.exp(-TAUTOP/mu0))*PLTOP
            
# WE CAN NOW SOLVE FOR THE COEFFICIENTS OF THE TWO STREAM
# CALL A SUBROUTINE THAT SOLVES  FOR THE FLUX TERMS
# WITHIN EACH INTERVAL AT THE MIDPOINT WAVENUMBER 
            
            FMUPI, FMDI, FNETI = solve_2stream(source_nu[:,NW], dtau_nu[:,NW,iG],
                    omega0_nu[:,NW,iG], g_asym_nu[:,NW,iG], mu0=mu0,
                    flux_top_dw=flux_top_dw_nu, alb_surf=alb_surf)
            FLUXUPI_nu[:,NW,iG] = FMUPI
            FLUXDWI_nu[:,NW,iG] = FMDI
            FNETI_nu[:,NW,iG] = FNETI
         
## NOW CALCULATE THE CUMULATIVE IR NET FLUX
#            NFLUXTOPI += FTOPUP*GWEIGHT[iG]
#            
## and same thing by spectral band... (RDW)
#            NFLUXTOPI_nu[NW] += FTOPUP * GWEIGHT[iG] 
#            
#            for L in range(L_NLEVRAD-1):
##           CORRECT FOR THE WAVENUMBER INTERVALS
#                FMNETI[L]  += (FMUPI[L]-FMDI[L]) * GWEIGHT[iG]
#                FLUXUPI[L] += FMUPI[L]*GWEIGHT[iG]
#                FLUXDNI[L] += FMDI[L]*GWEIGHT[iG]
##         and same thing by spectral band... (RW)
#                FLUXUPI_nu[L,NW] += FMUPI[L] * GWEIGHT[iG] 
#      #End NGAUSS LOOP
         
# *** END OF MAJOR SPECTRAL INTERVAL LOOP IN THE INFRARED****
      
    return FLUXUPI_nu, FLUXDWI_nu, FNETI_nu

@numba.jit(nopython=True,fastmath=True)
def solve_2stream(source, dtau, omega0, g_asym, mu0=0.5, flux_top_dw=0., 
                      alb_surf=0.):
    """Based on gfluxi in LMDZ code      
##-----------------------------------------------------------------------
##  THIS SUBROUTINE TAKES THE OPTICAL CONSTANTS AND BOUNDARY CONDITIONS
##  FOR THE INFRARED FLUX AT ONE WAVELENGTH AND SOLVES FOR THE FLUXES AT
##  THE LEVELS.  THIS VERSION IS SET UP TO WORK WITH LAYER OPTICAL DEPTHS
##  MEASURED FROM THE TOP OF EACH LAYER.  THE TOP OF EACH LAYER HAS  
##  OPTICAL DEPTH ZERO.  IN THIS SUB LEVEL N IS ABOVE LAYER N. THAT IS LAYER N
##  HAS LEVEL N ON TOP AND LEVEL N+1 ON BOTTOM. OPTICAL DEPTH INCREASES
##  FROM TOP TO BOTTOM.  SEE C.P. MCKAY, TGM NOTES.
##  THE TRI-DIAGONAL MATRIX SOLVER IS DSOLVER AND IS DOUBLE PRECISION SO MANY 
##  VARIABLES ARE PASSED AS SINGLE THEN BECOME DOUBLE IN DSOLVER
##
## NLL            = NUMBER OF LEVELS (NLAYERS + 1) MUST BE LESS THAT NL (101)
## TLEV(L_LEVELS) = ARRAY OF TEMPERATURES AT GCM LEVELS
## WAVEN          = WAVELENGTH FOR THE COMPUTATION
## DW             = WAVENUMBER INTERVAL
## dtau(NLAYER)   = ARRAY OPTICAL DEPTH OF THE LAYERS
## omega0(NLEVEL)     = SINGLE SCATTERING ALBEDO
## g_asym(NLEVEL) = ASYMMETRY FACTORS, 0=ISOTROPIC
## mu0          = AVERAGE ANGLE, MUST BE EQUAL TO 0.5 IN IR
## alb_surf            = SURFACE REFLECTANCE
## BTOP           = UPPER BOUNDARY CONDITION ON IR INTENSITY (NOT FLUX)
## BSURF          = SURFACE EMISSION = (1-RSFI)*PLANCK, INTENSITY (NOT FLUX)
## FP(NLEVEL)     = UPWARD FLUX AT LEVELS
## FM(NLEVEL)     = DOWNWARD FLUX AT LEVELS
## FMIDP(NLAYER)  = UPWARD FLUX AT LAYER MIDPOINTS
## FMIDM(NLAYER)  = DOWNWARD FLUX AT LAYER MIDPOINTS
##-----------------------------------------------------------------------
"""
#      
#      INTEGER NLL, NLAYER, L, NW, NT, NT2
#      REAL*8  TERM, CPMID, CMMID
#      REAL*8  PLANCK
#      REAL*8  EM,EP
#      REAL*8  g_asym(L_NLAYRAD), omega0(L_NLAYRAD), dtau(L_NLAYRAD)
#      REAL*8  TAUCUM(L_LEVELS), DTAUK
#      REAL*8  TLEV(L_LEVELS)
#      REAL*8  WAVEN, DW, mu0, alb_surf
#      REAL*8  BTOP, BSURF, FMIDP(L_NLAYRAD), FMIDM(L_NLAYRAD)
#      REAL*8  B0(L_NLAYRAD)
#      REAL*8  B1(L_NLAYRAD)
#      REAL*8  ALPHA(L_NLAYRAD)
#      REAL*8  LAMDA(L_NLAYRAD),XK1(L_NLAYRAD),XK2(L_NLAYRAD)
#      REAL*8  GAMA(L_NLAYRAD),CP(L_NLAYRAD),CM(L_NLAYRAD)
#      REAL*8  CPM1(L_NLAYRAD),CMM1(L_NLAYRAD),E1(L_NLAYRAD)
#      REAL*8  E2(L_NLAYRAD)
#      REAL*8  E3(L_NLAYRAD)
#      REAL*8  E4(L_NLAYRAD)
#      REAL*8  FTOPUP, FLUXUP, FLUXDN
#      REAL*8 :: TAUMAX = L_TAUMAX
#
## AB : variables for interpolation
#      REAL*8 C1
#      REAL*8 C2
#      REAL*8 P1
#      REAL*8 P2
#      
##=======================================================================
##     WE GO WITH THE HEMISPHERIC CONSTANT APPROACH IN THE INFRARED
#      
    TAUMAX=8.
    Nlay=dtau.size
    FMIDP=np.zeros((Nlay+1))
    FMIDM=np.zeros((Nlay+1))
    ALPHA=np.zeros((Nlay))
    LAMDA=np.zeros((Nlay))
    B0=np.zeros((Nlay))
    B1=np.zeros((Nlay))
    GAMA=np.zeros((Nlay))
    CPM1=np.zeros((Nlay))
    CMM1=np.zeros((Nlay))
    CP=np.zeros((Nlay))
    CM=np.zeros((Nlay))
    E1=np.zeros((Nlay))
    E2=np.zeros((Nlay))
    E3=np.zeros((Nlay))
    E4=np.zeros((Nlay))

    for L in range(Nlay):

        if omega0[L]==1.e0:
           omega0[L] = 0.99999e0
        
        ALPHA[L] = np.sqrt((1.e0-omega0[L])/(1.e0-omega0[L]*g_asym[L]) )
        LAMDA[L] = ALPHA[L]*(1.e0-omega0[L]*g_asym[L])/mu0
        
        B0[L] = source[L]
        B1[L] = (source[L+1] - B0[L]) / dtau[L]

        GAMA[L] = (1.e0-ALPHA[L])/(1.e0+ALPHA[L])
        TERM    = mu0/(1.e0-omega0[L]*g_asym[L])
         
# CPM1 AND CMM1 ARE THE CPLUS AND CMINUS TERMS EVALUATED
# AT THE TOP OF THE LAYER, THAT IS ZERO OPTICAL DEPTH
         
        CPM1[L] = B0[L]+B1[L]*TERM
        CMM1[L] = B0[L]-B1[L]*TERM
         
# CP AND CM ARE THE CPLUS AND CMINUS TERMS EVALUATED AT THE
# BOTTOM OF THE LAYER.  THAT IS AT dtau OPTICAL DEPTH.
# JL18 put CP and CM after the calculation of CPM1 and CMM1 to avoid unecessary calculations. 
         
        CP[L] = CPM1[L] +B1[L]*dtau[L] 
        CM[L] = CMM1[L] +B1[L]*dtau[L] 
      
# NOW CALCULATE THE EXPONENTIAL TERMS NEEDED
# FOR THE TRIDIAGONAL ROTATED LAYERED METHOD
# WARNING IF dtau[J) IS GREATER THAN ABOUT 35 (VAX)
# WE CLIP IT TO AVOID OVERFLOW.
      
    for L in range(Nlay):
        EP    = np.exp( min((LAMDA[L]*dtau[L]),TAUMAX)) # CLIPPED EXPONENTIAL
        EM    = 1.e0/EP
        E1[L] = EP+GAMA[L]*EM
        E2[L] = EP-GAMA[L]*EM
        E3[L] = GAMA[L]*EP+EM
        E4[L] = GAMA[L]*EP-EM
      
# DOUBLE PRECISION TRIDIAGONAL SOLVER
    BTOP=flux_top_dw
    BSURF=(1.-alb_surf)*source[-1]
 
    XK1,XK2 = DSOLVER(Nlay,GAMA,CP,CM,CPM1,CMM1,E1,E2,E3,E4,BTOP,
                BSURF,alb_surf)
      
# NOW WE CALCULATE THE FLUXES
      
    for L in range(Nlay):
        DTAUK = 0.
        EP    = np.exp(min(LAMDA[L]*DTAUK,TAUMAX)) # CLIPPED EXPONENTIAL 
        EM    = 1.e0/EP
        TERM  = mu0/(1.e0-omega0[L]*g_asym[L])
         
# CP AND CM ARE THE CPLUS AND CMINUS TERMS EVALUATED AT THE
# TOP OF THE LAYER.  THAT IS AT 0  OPTICAL DEPTH
         
        CPMID    = B0[L]+B1[L]*DTAUK +B1[L]*TERM
        CMMID    = B0[L]+B1[L]*DTAUK -B1[L]*TERM
        FMIDP[L] = XK1[L]*EP + GAMA[L]*XK2[L]*EM + CPMID
        FMIDM[L] = XK1[L]*EP*GAMA[L] + XK2[L]*EM + CMMID
        
# And now, for the special bottom layer

    L    = Nlay-1
    DTAUK=dtau[L]
    #DTAUK= 0.
    EP   = np.exp(min((LAMDA[L]*DTAUK),TAUMAX)) # CLIPPED EXPONENTIAL 
    EM   = 1.e0/EP
    TERM = mu0/(1.e0-omega0[L]*g_asym[L])

# CP AND CM ARE THE CPLUS AND CMINUS TERMS EVALUATED AT THE
# BOTTOM OF THE LAYER.  THAT IS AT dtau  OPTICAL DEPTH

    CPMID    = B0[L]+B1[L]*DTAUK +B1[L]*TERM
    CMMID    = B0[L]+B1[L]*DTAUK -B1[L]*TERM
    FMIDP[L+1] = XK1[L]*EP + GAMA[L]*XK2[L]*EM + CPMID
    FMIDM[L+1] = XK1[L]*EP*GAMA[L] + XK2[L]*EM + CMMID
 
    return   FMIDP, FMIDM, FMIDP-FMIDM

@numba.jit(nopython=True,fastmath=True)
def DSOLVER(NL,GAMA,CP,CM,CPM1,CMM1,E1,E2,E3,E4,BTOP,
                       BSURF,alb_surf):
    """
#  GCM2.0  Feb 2003
#
# DOUBLE PRECISION VERSION OF SOLVER

#*********************************************************
#* THIS SUBROUTINE SOLVES FOR THE COEFFICIENTS OF THE    *
#* TWO STREAM SOLUTION FOR GENERAL BOUNDARY CONDITIONS   *
#* NO ASSUMPTION OF THE DEPENDENCE ON OPTICAL DEPTH OF   *
#* C-PLUS OR C-MINUS HAS BEEN MADE.                      *
#* NL     = NUMBER OF LAYERS IN THE MODEL                *
#* CP     = C-PLUS EVALUATED AT TAO=0 (TOP)              *
#* CM     = C-MINUS EVALUATED AT TAO=0 (TOP)             *
#* CPM1   = C-PLUS  EVALUATED AT TAOSTAR (BOTTOM)        *
#* CMM1   = C-MINUS EVALUATED AT TAOSTAR (BOTTOM)        *
#* EP     = EXP(LAMDA*dtau)                              *
#* EM     = 1/EP                                         *
#* E1     = EP + GAMA *EM                                *
#* E2     = EP - GAMA *EM                                *
#* E3     = GAMA*EP + EM                                 *
#* E4     = GAMA*EP - EM                                 *
#* BTOP   = THE DIFFUSE RADIATION INTO THE MODEL AT TOP  *
#* BSURF  = THE DIFFUSE RADIATION INTO THE MODEL AT      *
#*          THE BOTTOM: INCLUDES EMMISION AND REFLECTION *
#*          OF THE UNATTENUATED PORTION OF THE DIRECT    *
#*          BEAM. BSTAR+alb_surf*FO*EXP(-TAOSTAR/U0)          *
#* alb_surf    = REFLECTIVITY OF THE SURFACE                  *
#* XK1    = COEFFICIENT OF THE POSITIVE EXP TERM         *
#* XK2    = COEFFICIENT OF THE NEGATIVE EXP TERM         *
#*********************************************************
"""
###      PARAMETER (NMAX=201)
#      IMPLICIT REAL*8  (A-H,O-Z)
#      DIMENSION GAMA(NL),CP(NL),CM(NL),CPM1(NL),CMM1(NL),XK1(NL),
#     *          XK2(NL),E1(NL),E2(NL),E3(NL),E4(NL)
#      DIMENSION AF(2*NL),BF(2*NL),CF(2*NL),DF(2*NL),XK(2*NL)
#======================================================================#

    L=2*NL

#    ************MIXED COEFFICENTS**********
#    THIS VERSION AVOIDS SINGULARITIES ASSOC.
#    WITH omega0=0 BY SOLVING FOR XK1+XK2, AND XK1-XK2.
    AF=np.empty((L))
    BF=np.empty((L))
    CF=np.empty((L))
    DF=np.empty((L))
    XK=np.empty((L))
    XK1=np.empty((NL))    
    XK2=np.empty((NL))    

    AF[0] = 0.0
    BF[0] = GAMA[0]+1.
    CF[0] = GAMA[0]-1.
    DF[0] = BTOP-CMM1[0]
    N     = 0
    LM2   = L-2

#     EVEN TERMS
 
    for I in range(1,LM2,2):
    #  DO I=2,LM2,2
        #print(I,N)
        AF[I] = (E1[N]+E3[N])*(GAMA[N+1]-1.)       
        BF[I] = (E2[N]+E4[N])*(GAMA[N+1]-1.)
        CF[I] = 2.0*(1.-GAMA[N+1]**2)
        DF[I] = (GAMA[N+1]-1.) * (CPM1[N+1] - CP[N]) + \
               (1.-GAMA[N+1])* (CM[N]-CMM1[N+1])
        N     = N+1
 
    N   = 0
    LM1 = L-1
    for I in range(2,LM1,2):
        AF[I] = 2.0*(1.-GAMA[N]**2)
        BF[I] = (E1[N]-E3[N])*(1.+GAMA[N+1])
        CF[I] = (E1[N]+E3[N])*(GAMA[N+1]-1.)
        DF[I] = E3[N]*(CPM1[N+1] - CP[N]) + E1[N]*(CM[N] - CMM1[N+1])
        N     = N+1
 
    AF[-1] = E1[-1]-alb_surf*E3[-1]
    BF[-1] = E2[-1]-alb_surf*E4[-1]
    CF[-1] = 0.0
    DF[-1] = BSURF-CP[-1]+alb_surf*CM[-1]

    XK = DTRIDGL(L,AF,BF,CF,DF)
 
#     ***UNMIX THE COEFFICIENTS****

    for N in range(NL):
        XK1[N] = XK[2*N]+XK[2*N+1]
        XK2[N] = XK[2*N]-XK[2*N+1]

#       NOW TEST TO SEE IF XK2 IS REALLY ZERO TO THE LIMIT OF THE
#       MACHINE ACCURACY  = 1 .E -30
#       XK2 IS THE COEFFICEINT OF THE GROWING EXPONENTIAL AND MUST
#       BE TREATED CAREFULLY

        if XK2[N] == 0.0: continue
#        IF (ABS (XK2(N)/XK(2*N-1)) .LT. 1.E-30) XK2(N)=0.0

        if (np.abs(XK2[N]/(XK[2*N-1]+1.e-20)) <= 1.E-30): XK2[N]=0.0   # For debug only (with -Ktrap=fp option)

    return XK1, XK2

@numba.jit(nopython=True,fastmath=True)
def DTRIDGL(L,AF,BF,CF,DF):
    """
    !  GCM2.0  Feb 2003

    !     DOUBLE PRECISION VERSION OF TRIDGL

          DIMENSION AF(L),BF(L),CF(L),DF(L),XK(L)
          DIMENSION AS(2*L),DS(2*L)

    !*    THIS SUBROUTINE SOLVES A SYSTEM OF TRIDIAGIONAL MATRIX
    !*    EQUATIONS. THE FORM OF THE EQUATIONS ARE:
    !*    A(I)*X(I-1) + B(I)*X(I) + C(I)*X(I+1) = D(I)

    !======================================================================!
    """
    AS=np.empty_like(AF)
    DS=np.empty_like(AF)
    XK=np.empty_like(AF)
    AS[-1] = AF[-1]/BF[-1]
    DS[-1] = DF[-1]/BF[-1]

    for I in range(1,L):
        X         = 1./(BF[L+1-I-2] - CF[L+1-I-2]*AS[L+2-I-2])
        AS[L+1-I-2] = AF[L+1-I-2]*X
        DS[L+1-I-2] = (DF[L+1-I-2]-CF[L+1-I-2]*DS[L+2-I-2])*X
 
    XK[0]=DS[0]
    for I in range(1,L):
        XKB   = XK[I-1]
        XK[I] = DS[I]-AS[I]*XKB
    return XK
    