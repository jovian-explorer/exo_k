# -*- coding: utf-8 -*-
"""
Created in Jan 2021

@author: jeremy leconte
"""
import numpy as np
import numba

jit=False

@numba.jit(nopython=True,fastmath=True)
def solve_2stream_nu_xsec(source_nu, dtau_nu, omega0_nu, g_asym_nu,
                flux_top_dw_nu, mu0=2./3., alb_surf=0., verbose=False, mid_layer=False):
    """Deals with the spectral axis
    """
    NLEV, NW = source_nu.shape
    flux_up=np.zeros((NLEV, NW))
    flux_dw=np.zeros((NLEV, NW))
    flux_net=np.zeros((NLEV, NW))
    kernel=np.zeros((NLEV, NLEV, NW))
    for iW in range(NW):
        flux_up[:,iW], flux_dw[:,iW], flux_net[:,iW], kernel[:,:,iW] = \
            solve_2stream(source_nu[:,iW], dtau_nu[:,iW],
                omega0_nu[:,iW], g_asym_nu[:,iW],
                mu0=mu0, flux_top_dw=flux_top_dw_nu[iW],
                alb_surf=alb_surf, verbose=verbose)
    return flux_up, flux_dw, flux_up-flux_dw, kernel
    
@numba.jit(nopython=True,fastmath=True)
def solve_2stream_nu_corrk(source_nu, dtau_nu, omega0_nu, g_asym_nu,
                flux_top_dw_nu, mu0=2./3., alb_surf=0., verbose=False, mid_layer=False):
    """Deals with the spectral axis
    """
    NLEV, NW = source_nu.shape
    NG = dtau_nu.shape[-1]
    flux_up=np.zeros((NLEV, NW, NG))
    flux_dw=np.zeros((NLEV, NW, NG))
    flux_net=np.zeros((NLEV, NW, NG))
    kernel=np.zeros((NLEV, NLEV, NW, NG))
    for iW in range(NW):
        for iG in range(NG):
            flux_up[:,iW, iG], flux_dw[:,iW,iG], flux_net[:,iW,iG], kernel[:,:,iW,iG] = \
                solve_2stream(source_nu[:,iW], dtau_nu[:,iW,iG], 
                    omega0_nu[:,iW,iG], g_asym_nu[:,iW,iG],
                    mu0=mu0, flux_top_dw=flux_top_dw_nu[iW],
                    alb_surf=alb_surf, verbose=verbose)
    return flux_up, flux_dw, flux_up-flux_dw, kernel


@numba.jit(nopython=True,fastmath=True)
def solve_2stream(source, dtau, omega0, g_asym, mu0=2./3., flux_top_dw=0.,
            alb_surf=0., verbose=False):
    """
    Inherited from ExoRem (Blain et al. 2020):
    https://gitlab.obspm.fr/dblain/exorem/-/blob/master/src/fortran/exorem/radiative_transfer.f90

    PURPOSE:                                                        

    THIS subroutine COMPUTES THE UPWARD, DOWNWARD AND NET THERMAL   
    FLUX IN AN INHOMOGENEOUS ABSORBING, SCATTERING ATMOSPHERE.      
    THE TWO_STREAM APPROXIMATION (QUADRATURE WITH COS OF            
    AVERAGE ANGLE = mu0) IS USED TO FIND THE DIFFUSE   	        
    REFLECTIVITY AND TRANSMISSIVITY AND THE TOTAL UPWARD AND        
    DOWNWARD FLUXES FOR EACH OF THE NLAY HOMOGENEOUS LAYERS.        
    THE ADDING METHOD IS then USED TO COMBINE THESE LAYERS.  IF     
    ANY LAYER IS THICKER THAN DTAU = *EMAX1*, IT IS ASSUMED TO BE   
    SEMI-INFINITE.  LAYERS THICKER THAN DTAU = *EMAX2*, ARE TREATED 
    AS INFINITE LAYERS.                                             

    NOTE: TO ACCOUNT FOR A DIFFUSE FLUX AT THE TOP OF THE ATMOS-    
          PHERE, THE USER MUST SET flux_dw(1) EQUAL TO THAT VALUE.     
          FOR NO DOWNWARD DIFFUSE FLUX AT THE TOP OF THE ATMOSPHERE,
          THE USER MUST INITIALIZE flux_dw(1) TO ZERO IN THE CALLING   
          PROGRAM.                                                  


    Parameters
    ----------
        source: array
            PLANCK function AT EACH LEVEL FROM TOP OF THE          
            ATMOSPHERE (L=0) TO THE SURFACE (L=NLAY) (NLAY+1 VALUES)  
            JL: actually seems to be pi times the planck function
            as defined in the rest of the code. 
        dtau: array
            ARRAY OF NORMAL INCIDENCE OPTICAL DEPTHS IN EACH       
            HOMOGENEOUS MODEL LAYER. (NLAY VALUES)                 
        omega0: array
            ARRAY OF SINGLE SCATTERING ALBEDOS FOR EACH HOMO-      
            GENEOUS MODEL LAYER. (NLAY VALUES)                     
        g_asym: array
            ARRAY OF ASSYMETRY parameterS FOR EACH HOMOGENEOUS     
            MODEL LAYER. (NLAY VALUES)                             
    
    Returns
    -------                                                        

        flux_up: array
            UPWARD FLUX AT NLAY+1 LAYER BOUNDARIES.                 
            (flux_up(L) REFERS TO THE UPWARD FLUX AT THE TOP          
            OF LAYER L)                                           
        flux_dw: array
            DOWNWARD FLUX AT NLAY+1 LAYER BOUNDARIES.               
            (flux_dw(L) REFERS TO THE DOWNWARD FLUX AT THE BOTTOM     
            OF LAYER L-1) 
        flux_up-flux_dw: array
            Net flux at the same levels

    Other Parameters
    ----------------
        mu0: float
            Cos of quadrature angle (Original mu0 was 2/3.)
        flux_top_dw: float
            Diffuse downward flux at upper interface.
    """
#            se atmosphere, only : n_levels
#
#            implicit none
#
#            
#            ********** COMMON BLOCKS USED IN DELTA-EDDINGTON ROUTINES.
#            
#            integer, intent(in) :: nlay
#            doubleprecision, intent(in) :: dtau(n_levels - 1),
#                   g_asym(n_levels - 1), source(n_levels)
#
#            doubleprecision, intent(inout) :: omega0(n_levels - 1)
#
#            doubleprecision, intent(out) :: flux_up(n_levels), flux_dw(n_levels),
#                 DKERNEL(n_levels, n_levels)
#
#            doubleprecision :: source2(n_levels), OMP(n_levels - 1), dtauP(n_levels - 1), 
#                RU(n_levels), &
#                DD(n_levels), DFLUX(n_levels - 1), G1(n_levels - 1), G2(n_levels - 1), &
#                DU(n_levels - 1), UFL(n_levels), DFL(n_levels), UFLUX(n_levels), 
#                EKT(n_levels - 1), &
#                SK(n_levels - 1), RL(n_levels), temperatures_layers(n_levels), RS(n_levels)
#
#            integer :: j, l, np2, nlev
#            doubleprecision :: skt, prec, emax1, emax2, alb, mu0, denom, e2ktm, emis,
#                flux_top_dw
    PREC = 1.e-10
    EMAX1 = 8.
    EMAX2 = 24.
    #EMAX1, EMAX2: optical depth at which layers are treated as semi infinite
    #        and infinite (resp). 
    NLAY = dtau.size
    NLEV = NLAY + 1
    NP2 = NLAY + 2

    dtauP=np.zeros_like(dtau)
    OMP=np.zeros_like(dtau)
    G1=np.zeros_like(dtau)
    G2=np.zeros_like(dtau)
    SK=np.zeros_like(dtau)
    DU=np.zeros_like(dtau)
    EKT=np.zeros_like(dtau)
    DFLUX=np.zeros_like(dtau)

    RL=np.zeros_like(source)
    RS=np.zeros_like(source)
    RU=np.zeros_like(source)
    DD=np.zeros_like(source)
    DFL=np.zeros_like(source)
    UFLUX=np.zeros_like(source)
    UFL=np.zeros_like(source)
    source2=np.zeros_like(source)
    temperatures_layers=np.zeros_like(source)
    flux_up=np.zeros_like(source)
    flux_dw=np.zeros_like(source)

    DKERNEL=np.zeros((NLEV,NLEV))
    #
    #****   SCALE THE OPTICAL DEPTHS, SINGLE SCATTERING ALBEDOS AND THE
    #       SCATTERING ASSYMETRY FACTORS FOR USE IN THE TWO-STREAM
    #       APPROXIMATION.  INITIALIZE OTHER QUANTITIES.  USE WISCOMBE'S
    #       TRICK TO SUBTRACT A SMALL VALUE FROM THE SINGLE SCATTERING
    #       FOR THE CASE OF A CONSERVATIVE ATMOSPHERE
    #
    for L in range(NLAY):
        if (1. - omega0[L])  <  PREC: omega0[L] = 1. - PREC
        dtauP[L] = (1. - omega0[L] * g_asym[L]) * dtau[L]
        OMP[L] = (1. - g_asym[L]) * omega0[L] / (1. - omega0[L] * g_asym[L])
        G1[L] = (1. - 0.5 * OMP[L]) / mu0
        G2[L] = 0.5 * OMP[L] / mu0
        SK[L] = np.sqrt(G1[L] * G1[L] - G2[L] * G2[L])
    #
    #**** DIFFUSE TRANSMITTANCE AND REFLECTANCE OF EACH HOMOGENEOUS LAYER
    #
    #     THE EQUATIONS FOR RL AND temperatures_layers WERE DERIVED BY SOLVING THE HOMOGENEOUS
    #     PART OF THE EQUATION OF TRANSFER (IE. NO SOURCE TERM)
    #
    for L in range(NLAY):
        SKT = SK[L] * dtauP[L]
        #print(L)
        if SKT  > EMAX2 :
        #
        #     INFINITE LAYERS
        #
            #print('EMAX2', L, SKT, G1[L], SK[L], (G1[L] + SK[L]))
            RL[L] = G2[L] / (G1[L] + SK[L])
            temperatures_layers[L] = 0.e0
            continue
        EKT[L] = np.exp(SKT)
        if SKT  >  EMAX1 :
        #
        #     SEMI-INFINITE LAYERS
        #
            #print('EMAX1', L, SKT, (G1[L] + SK[L]),(EKT[L] * (G1[L] + SK[L])))
            RL[L] = G2[L] / (G1[L] + SK[L])
            temperatures_layers[L] = 2.e0 * SK[L] / (EKT[L] * (G1[L] + SK[L]))
            continue
        E2KTM = EKT[L] * EKT[L] - 1.
        DENOM = G1[L] * E2KTM + SK[L] * (E2KTM + 2.e0)

        RL[L] = G2[L] * E2KTM / DENOM
        temperatures_layers[L] = 2.e0 * SK[L] * EKT[L] / DENOM
    #
    #****   SET THE "REFLECTIVITY", "TRANSMISSIVITY" AND "EMISSIVITY" OF
    #       THE "SURFACE" ASSUMING SEMI-INFINITE LAYER WITH SAME OMP
    #       AS LAYER NLAY. FOR "EMISSIVITY", ASSUME SAME dB/dTau AS LAYER NLAY
    #ALB = G2[NLAY-1] / (G1[NLAY-1] + SK[NLAY-1])
    #RL[NLEV-1] = (1. - ALB)
    #EMIS = 1. - ALB + (1. + ALB) * mu0 * \
    #    (1. - source[NLEV - 2] / source[NLEV-1]) / dtauP[NLAY-1]
    ##JL21 below, I put a new boundary condition at the bottom to mimick a dark surface.
    ALB=alb_surf
    RL[NLEV-1]=1.
    EMIS=1.-ALB
    temperatures_layers[NLEV-1] = 0.e0
    #
    #****   USE ADDING METHOD TO FIND THE REFLECTANCE AND TRANSMITTANCE
    #       OF COMBINED LAYERS.  ADD DOWNWARD FROM THE TOP AND UPWARD
    #       FROM THE BOTTOM AT THE SAME TIME.
    #
    RS[0] = RL[0]
    RU[NLEV-1] = ALB

    for L in range(NLAY):
        DD[L] = 1. / (1. - RS[L] * RL[L + 1])
        RS[L + 1] = RL[L + 1] + temperatures_layers[L + 1] * \
            temperatures_layers[L + 1] * RS[L] * DD[L]
        DU[NLEV - L - 2] = 1. / (1. - RL[NLEV - L - 2] * RU[NP2 - L - 2])
        RU[NLEV - L - 2] = RL[NLEV - L - 2] + temperatures_layers[NLEV - L - 2] * \
            temperatures_layers[NLEV - L - 2] * RU[NP2 - L - 2] * DU[NLEV - L - 2]
    #
    #****   COMPUTE THE UPWARD AND DOWNWARD FLUX FOR EACH HOMOGENEOUS LAYER
    #
    #**** LOOP OVER J FOR KERNEL
    for J in range(NLEV + 1):
        DFL[NLEV-1] = 0.e0
        if J <= NLEV-1:
            flux_dw[0] = 0.e0
            for L in range(NLEV):
                DKERNEL[J, L] = 0.e0
                source2[L] = 0.e0
            source2[J] = 1.
        else:
            flux_dw[0] = flux_top_dw
            for L in range(NLEV):
                source2[L] = source[L]
        UFL[NLEV-1] = EMIS * source2[NLEV-1]
        for L in range(NLAY):
            UFL[L], DFL[L]=_DEDIR1(source2[L], source2[L + 1], omega0[L],
                                g_asym[L], dtau[L], mu0)
        #
        #****   USE ADDING METHOD TO FIND UPWARD AND DOWNWARD FLUXES
        #       FOR COMBINED LAYERS.  START AT TOP
        #
        DFLUX[0] = DFL[0] + temperatures_layers[0] * flux_dw[0]

        for L in range(NLAY-1):
            DFLUX[L + 1] = temperatures_layers[L + 1] * \
                (RS[L] * (RL[L + 1] * DFLUX[L] + UFL[L + 1]) * DD[L] + DFLUX[L]) + DFL[L + 1]
        flux_dw[NLEV-1] = (DFLUX[NLAY-1] + RS[NLAY-1] * EMIS * source2[NLEV-1]) / \
            (1. - RS[NLAY-1] * ALB)
        #
        #****   USE ADDING METHOD TO FIND UPWARD AND DOWNWARD FLUXES
        #       FOR COMBINED LAYERS.  START AT BOTTOM.
        #
        UFLUX[NLEV-1] = UFL[NLEV-1]
        for L in range(NLAY):
            UFLUX[NLEV - L - 2] = temperatures_layers[NLEV - L - 2] * (UFLUX[NP2 - L - 2] + \
                    RU[NP2 - L - 2] * DFL[NLEV - L - 2]) * DU[NLEV - L - 2] + \
                    UFL[NLEV - L - 2]
        #
        #****   FIND THE TOTAL UPWARD AND DOWNWARD FLUXES AT INTERFACES
        #       BETWEEN INHOMOGENEOUS LAYERS.
        #
        flux_up[0] = UFLUX[0] + RU[0] * flux_dw[0]
        flux_up[NLEV-1] = ALB * flux_dw[NLEV-1] + EMIS * source2[NLEV-1]
        for L in range(NLAY-1):
            flux_up[NLEV - L - 2] = (UFLUX[NLEV - L - 2] + RU[NLEV - L - 2] * \
                    DFLUX[NLAY - L - 2]) / (1. - RU[NLEV - L - 2] * \
                    RS[NLAY - L - 2])
        for L in range(NLAY-1):
            flux_dw[L + 1] = DFLUX[L] + RS[L] * flux_up[L + 1]
        #
        if J <= NLEV-1:
            for L in range(NLEV):
                DKERNEL[J, L] = flux_up[L] - flux_dw[L]

    return flux_up, flux_dw, flux_up-flux_dw, DKERNEL

@numba.jit(nopython=True,fastmath=True)
def _DEDIR1(BATM1, BATM2, AIR, GIR, TAU, mu0):
    """            
            CCCCCCCCCCCCCCCCCCCCCCCCCCC  D E D I R 1  CCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            C                                                                    CC
            C    PURPOSE :                                                       CC
            C                                                                    CC
            C    THIS subroutine USES THE TWO-STREAM ROUTINE (QUADRATURE WITH    CC
            C    COS OF AVERAGE ANGLE = mu0) TO FIND THE UPWARD AND DOWNWARD    CC
            C    THERMAL FLUXES EMITTED FROM A HOMOGENEOUS LAYER.                CC
            C                                                                    CC
            C    AUTHORS:  DAVID PAIGE AND DAVID CRISP                           CC
            C                                                                    CC
            C    INPUT:                                                     CC
            C                                                                    CC
            C    BTATM1 IS ATMOSPHERIC EMISSION at level L (ARBITRARY UNITS)     CC
            C    BTATM2 IS ATMOSPHERIC EMISSION at level L+1 (ARBITRARY UNITS)   CC
            C    AIR IS IR SINGLE SCATTERING ALBEDO                              CC
            C    GIR IS IR ASYMMETRY parameter                                   CC
            C    TAU IS OPTICAL DEPTH (EXTINCTION = ABSORPTION + SCATTERING)     CC
            C                                                                    CC
            C    OUTPUT :                                                   CC
            C                                                                    CC
            C    FUPTOP IS UPWARD FLUX AT TOP OF ATM (SAME UNITS AS BTATM1 or 2) CC
            C    FDNBOT IS DOWNWARD FLUX AT SURFACE                              CC
            C    FUPBOT IS UPWARD FLUX AT SURFACE                                CC
            C                                                                    CC
            CCCCCCCCCCCCCCCCCCCCCCCCCCC  D E D I R 1  CCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    """        
#            doubleprecision, intent(in) :: BATM1, BATM2, air, gir, tau
#            doubleprecision, intent(out) :: fuptop, fdnbot
#
#            doubleprecision :: mu0, c1, c2, cap, dair, db, dtauIR, emkt, \
#                          epkt, omttp, opttp, v1, v2, v3, v4, v5, v6
    DAIR = AIR * (1. - GIR) / (1. - GIR * AIR)
    CAP = np.sqrt(1. - DAIR) / mu0
    OPTTP = np.sqrt(1. - 0.5 * DAIR + np.sqrt(1. - DAIR))
    OMTTP = np.sqrt(1. - 0.5 * DAIR - np.sqrt(1. - DAIR))
    dtauIR = (1. - GIR * AIR) * TAU
    db = (BATM2 - BATM1) * mu0 / dtauIR

    if CAP * dtauIR < 24.:
        #
        #****     THE LAYER THICKNESS IS FINITE
        #
        EPKT = np.exp(+CAP * dtauIR)
        EMKT = 1. / EPKT
        #
        #****     SET TOP B.C.
        #
        V1 = OPTTP
        V2 = OMTTP
        V3 = -(BATM1 - db)
        #
        #****     SET LOWER B.C.
        #
        V4 = EMKT * OMTTP
        V5 = EPKT * OPTTP
        V6 = -(BATM2 + db)
        #
        #****     SOLVE SYSTEM OF EQUATIONS
        #
        C1 = (V3 * V5 - V6 * V2) / (V1 * V5 - V4 * V2)
        C2 = (V1 * V6 - V3 * V4) / (V1 * V5 - V4 * V2)

        FUPTOP = C1 * OMTTP + C2 * OPTTP + BATM1 + db
        FDNBOT = C1 * EMKT * OPTTP + C2 * EPKT * OMTTP + BATM2 - db

    else:
        #
        #****     ASSUME LAYER IS SEMI-INFINITE.
        #
        FUPTOP = BATM1 * (1. - OMTTP / OPTTP) + db * (1. + OMTTP / OPTTP)
        FDNBOT = BATM2 * (1. - OMTTP / OPTTP) - db * (1. + OMTTP / OPTTP)
    return FUPTOP, FDNBOT
