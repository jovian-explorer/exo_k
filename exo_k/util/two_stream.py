# -*- coding: utf-8 -*-
"""
Created on Jun 20 2019

@author: jeremy leconte

Library of useful functions for radiative transfer calculations
"""
import numpy as np
import numba


#@numba.jit(nopython=True,fastmath=True)
def two_stream_quad_nu(NG, SOURCE_NU, DTAU_NU, OM_NU, G_NU,
                AMU1=2./3., FLX_TOP_DOWN_NU=0., PREC = 1.e-10,
                EMAX1 = 8., EMAX2 = 24.):
    """Deals with the spectral axis
    """
    NLEV=SOURCE_NU.shape[0]
    NW = DTAU_NU.shape[1]
    if NG is None:
        FLXU=np.zeros((NLEV, NW))
        FLXD=np.zeros((NLEV, NW))
        for iW in range(NW):
            FLXU[:,iW], FLXD[:,iW] = two_stream_quad(SOURCE_NU[:,iW], DTAU_NU[:,iW],
                OM_NU[:,iW], G_NU[:,iW],
                AMU1=AMU1, FLX_TOP_DOWN=FLX_TOP_DOWN_NU,
                PREC = PREC, EMAX1 = EMAX1, EMAX2 = EMAX2)
    else:
        FLXU=np.zeros((NLEV, NW, NG))
        FLXD=np.zeros((NLEV, NW, NG))
        #SOURCE_NU2=np.zeros((NLEV, NW, NG), dtype=float)
        #print(FLXU.shape, FLXD.shape, SOURCE_NU.shape, DTAU_NU.shape, OM_NU.shape, G_NU.shape)
        for iW in range(NW):
            for iG in range(NG):
                #print(iG,NG)
                #SOURCE_NU2[:,iW,iG]=SOURCE_NU[:,iW]
                FLXU[:,iW, iG],FLXD[:,iW, iG]=two_stream_quad(SOURCE_NU[:,iW], DTAU_NU[:,iW,iG], 
                    OM_NU[:,iW,iG], G_NU[:,iW,iG],
                    AMU1=AMU1, FLX_TOP_DOWN=FLX_TOP_DOWN_NU,
                    PREC = PREC, EMAX1 = EMAX1, EMAX2 = EMAX2)
    return FLXU, FLXD
    



@numba.jit(nopython=True,fastmath=True)
def two_stream_quad(SOURCE, DTAU, OM, G, AMU1=2./3., FLX_TOP_DOWN=0., PREC = 1.e-10,
                EMAX1 = 8., EMAX2 = 24.):
    """
    Inherited from ExoRem (Blain et al. 2020):
    https://gitlab.obspm.fr/dblain/exorem/-/blob/master/src/fortran/exorem/radiative_transfer.f90

    PURPOSE:                                                        

    THIS subroutine COMPUTES THE UPWARD, DOWNWARD AND NET THERMAL   
    FLUX IN AN INHOMOGENEOUS ABSORBING, SCATTERING ATMOSPHERE.      
    THE TWO_STREAM APPROXIMATION (QUADRATURE WITH COS OF            
    AVERAGE ANGLE = AMU1) IS USED TO FIND THE DIFFUSE   	        
    REFLECTIVITY AND TRANSMISSIVITY AND THE TOTAL UPWARD AND        
    DOWNWARD FLUXES FOR EACH OF THE NLAY HOMOGENEOUS LAYERS.        
    THE ADDING METHOD IS then USED TO COMBINE THESE LAYERS.  IF     
    ANY LAYER IS THICKER THAN DTAU = *EMAX1*, IT IS ASSUMED TO BE   
    SEMI-INFINITE.  LAYERS THICKER THAN DTAU = *EMAX2*, ARE TREATED 
    AS INFINITE LAYERS.                                             

    NOTE: TO ACCOUNT FOR A DIFFUSE FLUX AT THE TOP OF THE ATMOS-    
          PHERE, THE USER MUST SET FLXD(1) EQUAL TO THAT VALUE.     
          FOR NO DOWNWARD DIFFUSE FLUX AT THE TOP OF THE ATMOSPHERE,
          THE USER MUST INITIALIZE FLXD(1) TO ZERO IN THE CALLING   
          PROGRAM.                                                  


    Parameters
    ----------
        DTAU: array
            ARRAY OF NORMAL INCIDENCE OPTICAL DEPTHS IN EACH       
            HOMOGENEOUS MODEL LAYER. (NLAY VALUES)                 
        OM: array
            ARRAY OF SINGLE SCATTERING ALBEDOS FOR EACH HOMO-      
            GENEOUS MODEL LAYER. (NLAY VALUES)                     
        G: array
            ARRAY OF ASSYMETRY parameterS FOR EACH HOMOGENEOUS     
            MODEL LAYER. (NLAY VALUES)                             
        SOURCE: array
            PLANCK function AT EACH LEVEL FROM TOP OF THE          
            ATMOSPHERE (L=0) TO THE SURFACE (L=NLAY) (NLAY+1 VALUES)             
    
    Returns
    -------                                                        

        FLXU: array
            UPWARD FLUX AT NLAY+1 LAYER BOUNDARIES.                 
            (FLXU(L) REFERS TO THE UPWARD FLUX AT THE TOP          
            OF LAYER L)                                           
        FLXD: array
            DOWNWARD FLUX AT NLAY+1 LAYER BOUNDARIES.               
            (FLXD(L) REFERS TO THE DOWNWARD FLUX AT THE BOTTOM     
            OF LAYER L-1)                                         

    Other Parameters
    ----------------
        AMU1: float
            Cos of quadrature angle (Original AMU1 was 2/3.)
        FLX_TOP_DOWN: float
            Diffsue downward flux at upper interface.
        PREC: float
            precision
        EMAX1, EMAX2: float
            optical depth at which layers are treated as semi infinite
            and infinite (resp). 
    """
#            se atmosphere, only : n_levels
#
#            implicit none
#
#            
#            ********** COMMON BLOCKS USED IN DELTA-EDDINGTON ROUTINES.
#            
#            integer, intent(in) :: nlay
#            doubleprecision, intent(in) :: DTAU(n_levels - 1), G(n_levels - 1), SOURCE(n_levels)
#
#            doubleprecision, intent(inout) :: OM(n_levels - 1)
#
#            doubleprecision, intent(out) :: FLXU(n_levels), FLXD(n_levels),
#                 DKERNEL(n_levels, n_levels)
#
#            doubleprecision :: SOURCE2(n_levels), OMP(n_levels - 1), DTAUP(n_levels - 1), 
#                RU(n_levels), &
#                DD(n_levels), DFLUX(n_levels - 1), G1(n_levels - 1), G2(n_levels - 1), &
#                DU(n_levels - 1), UFL(n_levels), DFL(n_levels), UFLUX(n_levels), 
#                EKT(n_levels - 1), &
#                SK(n_levels - 1), RL(n_levels), temperatures_layers(n_levels), RS(n_levels)
#
#            integer :: j, l, np2, nlev
#            doubleprecision :: skt, prec, emax1, emax2, alb, AMU1, denom, e2ktm, emis,
#                FLX_TOP_DOWN

    NLAY = DTAU.size
    NLEV = NLAY + 1
    NP2 = NLAY + 2

    DTAUP=np.zeros_like(DTAU)
    OMP=np.zeros_like(DTAU)
    G1=np.zeros_like(DTAU)
    G2=np.zeros_like(DTAU)
    SK=np.zeros_like(DTAU)
    DU=np.zeros_like(DTAU)
    EKT=np.zeros_like(DTAU)
    DFLUX=np.zeros_like(DTAU)

    RL=np.zeros_like(SOURCE)
    RS=np.zeros_like(SOURCE)
    RU=np.zeros_like(SOURCE)
    DD=np.zeros_like(SOURCE)
    DFL=np.zeros_like(SOURCE)
    UFLUX=np.zeros_like(SOURCE)
    UFL=np.zeros_like(SOURCE)
    SOURCE2=np.zeros_like(SOURCE)
    temperatures_layers=np.zeros_like(SOURCE)
    FLXU=np.zeros_like(SOURCE)
    FLXD=np.zeros_like(SOURCE)

    DKERNEL=np.zeros((NLEV,NLEV))
    #
    #****   SCALE THE OPTICAL DEPTHS, SINGLE SCATTERING ALBEDOS AND THE
    #       SCATTERING ASSYMETRY FACTORS FOR USE IN THE TWO-STREAM
    #       APPROXIMATION.  INITIALIZE OTHER QUANTITIES.  USE WISCOMBE'S
    #       TRICK TO SUBTRACT A SMALL VALUE FROM THE SINGLE SCATTERING
    #       FOR THE CASE OF A CONSERVATIVE ATMOSPHERE
    #
    for L in range(NLAY):
        if (1. - OM[L])  <  PREC: OM[L] = 1. - PREC
        DTAUP[L] = (1. - OM[L] * G[L]) * DTAU[L]
        OMP[L] = (1. - G[L]) * OM[L] / (1. - OM[L] * G[L])
        G1[L] = (1. - 0.5 * OMP[L]) / AMU1
        G2[L] = 0.5 * OMP[L] / AMU1
        SK[L] = np.sqrt(G1[L] * G1[L] - G2[L] * G2[L])
    #
    #**** DIFFUSE TRANSMITTANCE AND REFLECTANCE OF EACH HOMOGENEOUS LAYER
    #
    #     THE EQUATIONS FOR RL AND temperatures_layers WERE DERIVED BY SOLVING THE HOMOGENEOUS
    #     PART OF THE EQUATION OF TRANSFER (IE. NO SOURCE TERM)
    #
    for L in range(NLAY):
        SKT = SK[L] * DTAUP[L]
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
    #JL21 below, I put a new boundary condition at the bottom to mimick a dark surface.
    #JL ALB = G2[NLAY-1] / (G1[NLAY-1] + SK[NLAY-1])
    ALB=0.
    #JL RL[NLEV-1] = (1. - ALB)
    RL[NLEV-1]=1.
    temperatures_layers[NLEV-1] = 0.e0
    #JL EMIS = 1. - ALB + (1. + ALB) * AMU1 * (1. - SOURCE[NLEV - 2] / SOURCE[NLEV-1]) / DTAUP[NLAY-1]
    EMIS=1.
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
            FLXD[0] = 0.e0
            for L in range(NLEV):
                DKERNEL[L, J] = 0.e0
                SOURCE2[L] = 0.e0
            SOURCE2[J] = 1.
        else:
            FLXD[0] = FLX_TOP_DOWN
            for L in range(NLEV):
                SOURCE2[L] = SOURCE[L]
        UFL[NLEV-1] = EMIS * SOURCE2[NLEV-1]
        for L in range(NLAY):
            UFL[L], DFL[L]=DEDIR1(SOURCE2[L], SOURCE2[L + 1], OM[L], G[L], DTAU[L])
        #
        #****   USE ADDING METHOD TO FIND UPWARD AND DOWNWARD FLUXES
        #       FOR COMBINED LAYERS.  START AT TOP
        #
        DFLUX[0] = DFL[0] + temperatures_layers[0] * FLXD[0]

        for L in range(NLAY-1):
            DFLUX[L + 1] = temperatures_layers[L + 1] * \
                (RS[L] * (RL[L + 1] * DFLUX[L] + UFL[L + 1]) * DD[L] + DFLUX[L]) + DFL[L + 1]
        FLXD[NLEV-1] = (DFLUX[NLAY-1] + RS[NLAY-1] * EMIS * SOURCE2[NLEV-1]) / \
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
        FLXU[0] = UFLUX[0] + RU[0] * FLXD[0]
        FLXU[NLEV-1] = ALB * FLXD[NLEV-1] + EMIS * SOURCE2[NLEV-1]
        for L in range(NLAY-1):
            FLXU[NLEV - L - 2] = (UFLUX[NLEV - L - 2] + RU[NLEV - L - 2] * \
                    DFLUX[NLAY - L - 2]) / (1. - RU[NLEV - L - 2] * \
                    RS[NLAY - L - 2])
        for L in range(NLAY-1):
            FLXD[L + 1] = DFLUX[L] + RS[L] * FLXU[L + 1]
        #
        if J <= NLEV-1:
            for L in range(NLEV):
                DKERNEL[L, J] = FLXU[L] - FLXD[L]

    return FLXU, FLXD #, DKERNEL

@numba.jit(nopython=True,fastmath=True)
def DEDIR1(BATM1, BATM2, AIR, GIR, TAU):
    """            
            CCCCCCCCCCCCCCCCCCCCCCCCCCC  D E D I R 1  CCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            C                                                                    CC
            C    PURPOSE :                                                       CC
            C                                                                    CC
            C    THIS subroutine USES THE TWO-STREAM ROUTINE (QUADRATURE WITH    CC
            C    COS OF AVERAGE ANGLE = AMU1) TO FIND THE UPWARD AND DOWNWARD    CC
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
#            doubleprecision :: AMU1, c1, c2, cap, dair, db, DTAUIR, emkt, \
#                          epkt, omttp, opttp, v1, v2, v3, v4, v5, v6
    AMU1 = 2. / 3.
    DAIR = AIR * (1. - GIR) / (1. - GIR * AIR)
    CAP = np.sqrt(1. - DAIR) / AMU1
    OPTTP = np.sqrt(1. - 0.5 * DAIR + np.sqrt(1. - DAIR))
    OMTTP = np.sqrt(1. - 0.5 * DAIR - np.sqrt(1. - DAIR))
    DTAUIR = (1. - GIR * AIR) * TAU
    db = (BATM2 - BATM1) * AMU1 / DTAUIR

    if CAP * DTAUIR < 24.:
        #
        #****     THE LAYER THICKNESS IS FINITE
        #
        EPKT = np.exp(+CAP * DTAUIR)
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
