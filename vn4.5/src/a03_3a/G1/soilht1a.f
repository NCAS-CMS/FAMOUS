C ******************************COPYRIGHT******************************
C (c) CROWN COPYRIGHT 1995, METEOROLOGICAL OFFICE, All Rights Reserved.
C
C Use, duplication or disclosure of this code is subject to the
C restrictions as set forth in the contract.
C
C                Meteorological Office
C                London Road
C                BRACKNELL
C                Berkshire UK
C                RG12 2SZ
C 
C If no contract has been raised with this copy of the code, the use,
C duplication or disclosure of it is strictly prohibited.  Permission
C to do so must first be obtained in writing from the Head of Numerical
C Modelling at the above address.
C ******************************COPYRIGHT******************************
C
C*LL  SUBROUTINE SOIL_HTF-----------------------------------------------
CLL
CLL  Purpose:  For land points, calculates the heat flux between the
CLL            Top two soil layers, and updates the deep soil
CLL            temperatures.  For other points, or if there is only
CLL            one soil layer, the flux is simply set to zero.
CLL
CLL            Note the terminology used: the scheme has n soil layers,
CLL            where n is given by the argument NSOIL ( = DS_LEVELS+1).
CLL            Layer 1 is the top soil layer, assumed to be at
CLL            temperature TSTAR.  The layers beneath this one are
CLL            termed "deep" soil layers, of which there are obviously
CLL            n-1.  The first layer below the top soil layer is
CLL            referred to indifferently as "soil layer 2" or "deep soil
CLL            layer 1".
CLL       Before March 1990 this routine was specific to a 4 soil layer
CLL            scheme.  To use any other number of layers, replace the
CLL            definition of PARAMETER PSOIL, which is used to dimension
CLL            a couple of small work arrays.
CLL
CLL
CLL F.Hewer     <- programmer of some or all of previous code or changes
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  Date
CLL   3.4   06/06/94  DEF TIMER replaced by LOGICAL LTIMER
CLL                   Argument LTIMER added
CLL                                                 S.J.Swarbrick
CLL   4.1   08/05/96  decks A03_2C and A03_3B removed
CLL                                     S D Jackson
CLL   4.2   Oct. 96   T3E migration - *DEF CRAY removed
CLL                                     S J Swarbrick                   
!LL   4.3   06/03/97  Dimension T_DEEP_SOIL by LAND_FIELD
!LL                   for non-MPP and MPP runs. D. Robinson.
!LL   4.5   18/06/98  Changed Timer calls to indicate non-barrier
!LL                                                   P.Burton
CLL   4.5    Jul. 98  Kill the IBM specific lines (JCThil)
CLL
CLL  Programming standard: Unified Model Documentation Paper No.4
CLL                        version no. 2, dated 18/1/90.
CLL
CLL  System component covered: P242.
CLL
CLL  Documentation: Unified Model Documentation Paper No 24.
C*
C*L ARGUMENTS:----------------------------------------------------------
      SUBROUTINE SOIL_HTF(
     + HCAP,HCON,LAYER_DEPTH,LYING_SNOW,TSTAR,LAND_MASK,TIMESTEP,
     + P_FIELD,LAND_FIELD,P_POINTS,P1,
     + LAND_PTS,LAND1,LAND_INDEX,
     + NSOIL,T_DEEP_SOIL,ASOIL_1,SOIL_HT_FLUX
     +,LTIMER)
      IMPLICIT NONE
      LOGICAL LTIMER
      INTEGER
     + P_FIELD               ! IN Number of gridpoints in field.
     +,LAND_FIELD            ! IN Number of land points in field.
     +,P_POINTS              ! IN Number of gridpoints to be processed.
     +,P1                    ! IN First point processed within field.
     +,NSOIL                 ! IN Number of soil layers (N_S in doc).
     +,LAND_PTS              ! IN Number of land points to be processed.
     +,LAND1                 ! IN First land point to be processed.
     +,LAND_INDEX(P_FIELD)   ! IN Index of land points on the P-grid.
C                            !    The ith element contains the position
C                            !    in whole grid of the ith land point.
C                            !    (Must match parameter PSOIL.)
      REAL
     + HCAP(LAND_FIELD)      ! IN Soil heat capacity in J/K/m**3
C                            !    (C_S in documentation).
     +,HCON(LAND_FIELD)      ! IN Soil thermal conductivity in W/m/K
C                            !    (LAMBDA_S in documentation).
     +,LAYER_DEPTH(NSOIL)    ! IN Soil layer depths as multiples of
C                            !    depth of top layer (ZETA_r in
C                            !    documentation).  The following were
C                            !    used in the GCM (5th Ann Cyc) :-
C
C            LAYER_DEPTH /1.000,3.908,14.05,44.65/  (see eqns P242.15)
C
C           (LAYER_DEPTH(1) must of course be 1.0, by definition.)
C
     +,LYING_SNOW(P_FIELD)   ! IN Lying snow amount (kg per sq m).
     +,TSTAR(P_FIELD)        ! IN Surface (i.e. soil layer 1) temp (K).
      LOGICAL
     + LAND_MASK(P_FIELD)    ! IN Land mask (T if land, else F).
      REAL TIMESTEP          ! IN Timestep (s).
      REAL
     + T_DEEP_SOIL(LAND_FIELD,NSOIL-1)
C                            ! INOUT Deep soil temperatures (K).
      REAL
     + ASOIL_1(P_FIELD)      ! OUT Soil coefficient used elsewhere in
C                            !     P24 (sq m K per Joule * timestep).
     +,SOIL_HT_FLUX(P_FIELD) ! OUT Heat flux from soil layer 1 to
C                            !     soil layer 2, i.e. +ve values for
C                            !     downward flux (Watts per sq m).
      INTEGER ERROR          ! OUT Error flag: 0 if AOK;
     +                       !     1 if NSOIL & PSOIL are not the same.
C*
C  Local physical constants --------------------------------------------
C  (Must be here because contains PSOIL, used for later dimensioning.)
C RHO_WATER removed to avoid clash with declaration in C_DENSTY
C J.Smith 28/06/95
      REAL OMEGA1,RHO_SNOW,DEFF_SNOW,SNOW_HCON,SNOW_HCAP
      INTEGER PSOIL
      PARAMETER (
     + PSOIL=4                  ! No. of soil layers (must = NSOIL).
     +,OMEGA1=3.55088E-4        ! Tunable characteristic freq (rad/s).
     +,RHO_SNOW=250.0           ! Density of lying snow (kg per m**3).
     +,DEFF_SNOW=0.1            ! Depth of `effective' snow surface
C                               ! layer (m).
     +,SNOW_HCON=0.265          ! Thermal conductivity of lying snow
C                               ! (Watts per m per K).
     +,SNOW_HCAP=0.63E6         ! Thermal capacity of lying snow
C                               ! (J/K/m3)
     +)
C*L  Workspace usage ---------------------------------------------------
C  Two blocks of real space, plus a few odd bits, are needed.
      REAL
     + DELT1(LAND_FIELD) ! Temperature difference across "this" soil lay
     +,DELT2(LAND_FIELD) ! Temperature difference across "next" soil lay
      REAL ASOIL(NSOIL),BSOIL(NSOIL)
C
      EXTERNAL TIMER
C*
C  Local variables -----------------------------------------------------
      INTEGER
     + I           ! Loop counter; horizontal land point field index.
     +,J           ! Loop counter; horizontal land and sea field index.
     +,THIS_LEVEL  ! Loop counter for loop through (deep) soil levels.
     +,NEXT_LEVEL  ! "Next" level during loops through soil levels.
      REAL
     + DS_RATIO   ! 2 * (Ratio of actual snowdepth to depth of top
C                 !      soil layer).
     +,SIFACT     ! Snow Insulation FACTor, GAMMA_SNOW in documentation.
     +,Z2_PLUS_1  ! Combined depth of top 2 soil layers.
      IF (LTIMER) THEN
      CALL TIMER('SOILHTF ',103)
      ENDIF
      Z2_PLUS_1 = LAYER_DEPTH(1) + LAYER_DEPTH(2)
      IF (NSOIL.GT.1) THEN
        DO 1 THIS_LEVEL = 1,NSOIL-1
          NEXT_LEVEL = THIS_LEVEL + 1
C-----------------------------------------------------------------------
CL 1. Set conductivity coefficients ASOIL and BSOIL, absorbing TIMESTEP
CL    from the LHS of eqns P242.3 and P242.4, where they are used.
CL    ASOIL/TIMESTEP is ASr, and BSOIL/TIMESTEP is BSr, in eqns
CL    P242.12, P242.13.  ASOIL(1), i.e. A1, is a function of soil type
CL    and therefore varies from point to point.  Therefore it is set in
CL    loop round points, and not here.
C-----------------------------------------------------------------------
          ASOIL(NEXT_LEVEL) = -OMEGA1 * TIMESTEP /
     +      ( LAYER_DEPTH(NEXT_LEVEL) *
     +        (LAYER_DEPTH(THIS_LEVEL) + LAYER_DEPTH(NEXT_LEVEL)) )
          BSOIL(THIS_LEVEL) = OMEGA1 * TIMESTEP /
     +      ( LAYER_DEPTH(THIS_LEVEL) *
     +        (LAYER_DEPTH(THIS_LEVEL) + LAYER_DEPTH(NEXT_LEVEL)) )
    1   CONTINUE
        BSOIL(NSOIL)=0.0
CMIC$ DO ALL VECTOR SHARED(P_POINTS, P1, P_FIELD, LAND_MASK,
CMIC$1   SOIL_HT_FLUX) PRIVATE(J)
CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
        DO 2 J=P1,P1+P_POINTS-1
C-----------------------------------------------------------------------
CL 2.  Set soil heat flux to zero over sea points (including sea-ice).
C-----------------------------------------------------------------------
          IF (.NOT.LAND_MASK(J)) THEN
            SOIL_HT_FLUX(J)=0.0
          ENDIF                      !  (If a land point.)
    2   CONTINUE
CMIC$ DO ALL VECTOR SHARED(LAND_PTS, LAND1,TIMESTEP, Z2_PLUS_1, P_FIELD,
CMIC$1   LAND_FIELD, SOIL_HT_FLUX, HCON, HCAP, ASOIL_1, LYING_SNOW,
CMIC$2   T_DEEP_SOIL, TSTAR, DELT1, BSOIL, DELT2, ASOIL, LAND_INDEX)
CMIC$3   PRIVATE(SIFACT, I, J, DS_RATIO)
CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
        DO 21 I=LAND1,LAND1+LAND_PTS-1
            J = LAND_INDEX(I)
C-----------------------------------------------------------------------
CL 3.   Land point calculations multi-soil-layer model.
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
CL 3.1  Set A1 - see eqns P242.11, P242.14  Units: J-1 m2 K * timestep.
CL      (The variable here is set to A1 * timestep.)
C-----------------------------------------------------------------------
            ASOIL_1(J) = SQRT( OMEGA1 / (2.0 * HCON(I) * HCAP(I)) )
     +                   * TIMESTEP
C-----------------------------------------------------------------------
CL 3.2 Initialise factor used to modify conductivity coefficients in
CL     the presence of lying snow (Snow Insulation FACTor).
C-----------------------------------------------------------------------
            SIFACT=1.0
C-----------------------------------------------------------------------
CL 3.3 Calculate factor which modifies the thermal conductivity between
CL    the top two subsurface layers in order to represent the insulating
CL    effects of lying snow.  See eqn P242.19.
C     Note that Z2_PLUS_1 is in fact (1+ZETA2), as the depth of the top
C     soil layer is by definition (of the units) 1.0.
C     Also note, LYING_SNOW is held on the whole grid, not just land pts
C-----------------------------------------------------------------------
            IF (LYING_SNOW(J).GT.0.0) THEN
              DS_RATIO = 2.0 *
     +          LYING_SNOW(J) / (RHO_SNOW)  ! Actual depth of lying snow
C                                           ! in metres.
     +          *
     +        ( ASOIL_1(J) * HCAP(I)  ! Reciprocal of thickness of top
     +          / TIMESTEP )          ! soil layer in metres. Equivalent
C                                     ! formula for thickness :-
C                                     ! SQRT(2*HCON/(HCAP*OMEGA1)).
C
              IF (DS_RATIO.LE.1.0) THEN
                SIFACT = 1.0 / (1.0 + DS_RATIO/Z2_PLUS_1)
              ELSEIF (LYING_SNOW(J) .LT. 5.0E3) THEN  ! See below ***
                SIFACT = Z2_PLUS_1 / (
     +                   (HCON(I)/SNOW_HCON) * (DS_RATIO-1.0)
     +                   + 1.0 + Z2_PLUS_1
     +                   )
              ELSE     ! *** See final paragraph of P242 documentation.
                SIFACT=1.0
              ENDIF    ! DSRATIO <= 0
            ENDIF      ! Lying snow
C-----------------------------------------------------------------------
CL 3.4 Calculate heat flux from top to next-to-top soil layer (+ve
CL     downwards), in Watts per square metre.  See eqn P242.8 (middle
CL     line), as modified according to P242.21.
C-----------------------------------------------------------------------
            DELT1(I) = T_DEEP_SOIL(I,1) - TSTAR(J)
            SOIL_HT_FLUX(J) = -SIFACT
     +                        * ( BSOIL(1)/(ASOIL_1(J)) )
     +                        * DELT1(I)
C-----------------------------------------------------------------------
CL 3.5 Update deep soil temperatures.
C-----------------------------------------------------------------------
CL 3.6 Update first deep soil layer (i.e. soil layer 2).  This is done
CL     separately because of the need to modify ASOIL(2) to account for
CL     snow insulation.  Eqn P242.3, modified according to P242.20.
C-----------------------------------------------------------------------
            DELT2(I) = T_DEEP_SOIL(I,2) - T_DEEP_SOIL(I,1)
            T_DEEP_SOIL(I,1) = T_DEEP_SOIL(I,1) +
     +        SIFACT * ASOIL(2) * DELT1(I) +
     +          BSOIL(2)*DELT2(I)
            DELT1(I)=DELT2(I)
   21   CONTINUE
C-----------------------------------------------------------------------
CL 3.7 Update deep soil layers 2 to (NSOIL-2), i.e. soil layers 3 to
CL     NSOIL-1, if these layers exist in the current model.  Eqn P242.3.
C      Loop counter (THIS_LEVEL) ranges over deep soil layers.
C-----------------------------------------------------------------------
            IF (NSOIL.GT.3) THEN
              DO 3 THIS_LEVEL = 2,NSOIL-2
                NEXT_LEVEL = THIS_LEVEL + 1
CMIC$ DO ALL VECTOR SHARED(LAND_PTS, LAND1, NEXT_LEVEL, THIS_LEVEL
CMIC$1   , LAND_FIELD, LAND_MASK, NSOIL, T_DEEP_SOIL, DELT2, ASOIL
CMIC$2   , DELT1, BSOIL) PRIVATE(I)
CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
                DO 31 I=LAND1,LAND1+LAND_PTS-1
                    DELT2(I) = T_DEEP_SOIL(I,NEXT_LEVEL)
     +                             - T_DEEP_SOIL(I,THIS_LEVEL)
                    T_DEEP_SOIL(I,THIS_LEVEL)=T_DEEP_SOIL(I,THIS_LEVEL)
     +                + ASOIL(NEXT_LEVEL) * DELT1(I)
     +                  + BSOIL(NEXT_LEVEL) * DELT2(I)
                    DELT1(I)=DELT2(I)
   31           CONTINUE
    3         CONTINUE
            ENDIF
C-----------------------------------------------------------------------
CL 3.8 Update deepest soil layer temperature, if this hasn't been done
CL     already.  This is done separately simply to avoid an IF in the
CL     DO 2 loop (there is obviously no next-layer temperature for the
CL     bottom layer).  See eqn P242.4.
C-----------------------------------------------------------------------
            IF (NSOIL.GT.2) THEN
CMIC$ DO ALL VECTOR SHARED(LAND_PTS, LAND1, NSOIL, LAND_FIELD,
CMIC$1   LAND_MASK, T_DEEP_SOIL, ASOIL, DELT1) PRIVATE(I)
CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
              DO 4 I=LAND1,LAND1+LAND_PTS-1
                  T_DEEP_SOIL(I,NSOIL-1) = T_DEEP_SOIL(I,NSOIL-1)
     +              + ASOIL(NSOIL) * DELT1(I)
    4         CONTINUE
            ENDIF
C-----------------------------------------------------------------------
CL 4. Finally, tidy up for the remaining case.
C-----------------------------------------------------------------------
      ELSE                           !  (If 1-layer soil scheme).
CMIC$ DO ALL VECTOR SHARED(P_POINTS,P1,P_FIELD,SOIL_HT_FLUX)PRIVATE(J)
CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
        DO 5 J=P1,P1+P_POINTS-1
          SOIL_HT_FLUX(J)=0.0
    5   CONTINUE
      ENDIF                          !  ENDIF for > 1 soil layer.
    6 CONTINUE                       !  Branch for error exit.
      IF (LTIMER) THEN
      CALL TIMER('SOILHTF ',104)
      ENDIF
      RETURN
      END
