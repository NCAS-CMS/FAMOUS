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
C*LL  SUBROUTINE TR_MIX ------------------------------------------------
CLL
CLL  Purpose: Calculate tracer flux and pass through to IMP_MIX to solve
CLL
CLL  Suitable for single column use; activate *IF definition IBM.
CLL
CLL  SDJ  <- Programmers of some or all of previous code or changes
CLL
CLL  Model           Modification history:
CLL version  Date
CLL
CLL   3.4  18/10/94   *DECK inserted into UM version 3.4. S Jackson
CLL   4.1  02/05/96  Surface emissions and dry deposition coefficients
CLL                  added as input arguments; surface deposition flux
CLL                  added as output argument.         M.Woodage
CLL   4.2   Oct. 96   T3E migration - *DEF CRAY removed
CLL                                     S J Swarbrick
CLL   4.5    Jul. 98  Kill the IBM specific lines (JCThil)
CLL
CLL  Programming standard: UM Documentation Paper No 4, Version 2,
CLL                        dated 18/1/90
CLL
CLL  System component covered: P244
CLL
CLL  Project task: P24
CLL
CLL  Documentation: UM Documentation Paper No 24.
CLL
C*----------------------------------------------------------------------
C*L  Arguments :-
      SUBROUTINE TR_MIX (
     & P_FIELD,BL_LEVELS,FIRST_ROW,ROW_LENGTH,N_ROWS
     &,DELTA_AK,DELTA_BK
     &,GAMMA_RHOKH_RDZ,RHOKH_1
     &,PSTAR,TIMESTEP
     &,F_FIELD,FIELD
     &,SURF_EM,RES_FACTOR,SURF_DEP_FLUX
     &,NRML
     &,ERROR,LTIMER
     &)

      IMPLICIT NONE
      INTEGER
     & P_FIELD                     ! IN No. of points in P-grid.
     &,BL_LEVELS                   ! IN No. of atmospheric levels for
C                                  !    which boundary layer fluxes are
C                                  !    calculated.
     &,FIRST_ROW                   ! IN First row of data to be treated,
C                                  !    referred to P-grid (must be > 1
C                                  !    since "polar" rows are never
C                                  !    treated).
     &,ROW_LENGTH                  ! IN No. of points in latitude row.
     &,N_ROWS                      ! IN No. of rows of data to be
C                                  !    treated, referred to P-grid.
C                                  !    FIRST_ROW+N_ROWS-1 must be less
C                                  !    than P_ROWS, since "polar" rows
C                                  !    are never treated.
      REAL
     & DELTA_AK(BL_LEVELS)         ! IN Difference of hybrid 'A' across
C                                  !    boundary layers (K-1/2 to K+1/2)
C                                  !    (upper minus lower).
     &,DELTA_BK(BL_LEVELS)         ! IN Difference of hybrid 'B' across
C                                  !    boundary layers (K-1/2 to K+1/2)
C                                  !    (upper minus lower).
     &,GAMMA_RHOKH_RDZ(P_FIELD,2:BL_LEVELS)
C                                  ! IN Mixing coeff. above surface
C                                  !    = GAMMA(K)*RHOKH(,K)*RDZ(K)
C                                  !    for K>=2 (from KMKH).
     &,RHOKH_1(P_FIELD)            ! IN  Surface exchange coeff.
C                                  !     from P243 (SF_EXCH)
     &,PSTAR(P_FIELD)              ! IN Surface pressure (Pa).
     &,TIMESTEP                    ! IN Timestep in seconds.
     &,SURF_EM(P_FIELD)            ! IN, Surface emissions in kg/m2/s
     &,RES_FACTOR(P_FIELD)         ! IN, dry dep coeff=Ra/(Ra+Rb+Rc)
C
      REAL
     & F_FIELD(P_FIELD,BL_LEVELS)  ! OUT Flux of tracer in kg/m2/s.
     &,FIELD(P_FIELD,BL_LEVELS)    ! INOUT Tracer amount in kg/kg.
     &,SURF_DEP_FLUX(P_FIELD)      ! OUT, surface deposn flux (kg/m2/s)
C
      INTEGER
     & NRML(P_FIELD)               ! IN The number of model layers
C                                  !    in the unstable rapidly mixing
C                                  !    layer. Zero if surface layer
C                                  !    is stable.
     &,ERROR                       ! OUT 1 if bad arguments, else 0.
C
      LOGICAL
     & LTIMER                      ! IN Logical switch for TIMER
C                                  !    diagnostics
C*
C*L  External references :-
      EXTERNAL IMP_MIX,TIMER
C*
C*L  Local and other symbolic constants :-
C*L------------------COMDECK C_G----------------------------------------
C G IS MEAN ACCEL DUE TO GRAVITY AT EARTH'S SURFACE
      REAL G

      PARAMETER(G=9.80665)
C*----------------------------------------------------------------------
C*L------------------COMDECK C_GAMMA------------------------------------
C GAMMA holds the implicit weighting coeff. for up to 30 B.L. levels.
C It is only required for the the number of B.L. levels actually used,
C so it does not need to be set up to 30 when less BL levels are used.
      REAL GAMMA(30)       ! Max of 30 Boundary Layer levels assumed.
C
      DATA GAMMA / 2 * 2.0 , 1.5 , 27 * 1.0 /
C*----------------------------------------------------------------------
C*
C*L Workspace :-
C
      REAL
     & RHOK_DEP(P_FIELD)     ! Surface deposition coeficient
C*
C  Local scalars :-
      INTEGER
     & P_POINTS ! Number of points on P-grid
     &,P1       ! First point in P-grid
     &,P_ROWS   ! Number of rows in P-grid
     &,I        ! Loop counter (horizontal field index).
     &,K        ! Loop counter (vertical index).
     &,KM1      ! K minus 1.
     &,KP1      ! K plus 1.
C
C-----------------------------------------------------------------------
CL  0.  Check that the scalars input to define the grid are consistent.
C       See comments to routine SF_EXCH for details.
C-----------------------------------------------------------------------
C
      IF (LTIMER) THEN
        CALL TIMER('TRMIX   ',3)
      ENDIF
C-----------------------------------------------------------------------
CL    Set pointers, etc.
C-----------------------------------------------------------------------
C
      P_POINTS = N_ROWS * ROW_LENGTH
      P1 = 1 + (FIRST_ROW-1)*ROW_LENGTH
      P_ROWS = N_ROWS
C

      ERROR=0
C
C
C-----------------------------------------------------------------------
CL 1.  Calculate flux of tracer:
C-----------------------------------------------------------------------
CL 1.1 Above the surface
C-----------------------------------------------------------------------
C
      DO 1 K=2,BL_LEVELS
        DO 11 I=P1,P1+P_POINTS-1
          F_FIELD(I,K) = - (GAMMA_RHOKH_RDZ(I,K) / GAMMA(K)) *
     &                                   (FIELD(I,K) - FIELD(I,K-1))
   11   CONTINUE
   1  CONTINUE
C
C-----------------------------------------------------------------------
CL 1.2 At the surface: (i) set surface flux equal to input emissions
CL                   (should be passed in as ancillary file, else ZERO)
CL                     (ii) Use input resistance factors to calculate
CL                   surface deposition (if ZERO then no dry deposition)
C-----------------------------------------------------------------------
C
        DO 21 I=P1,P1+P_POINTS-1
          F_FIELD(I,1) = SURF_EM(I)       ! Inject surface emissions
          RHOK_DEP(I) = RES_FACTOR(I) * RHOKH_1(I)
          SURF_DEP_FLUX(I) = -RHOK_DEP(I) * FIELD(I,1)
          RHOK_DEP(I) = GAMMA(1) * RHOK_DEP(I)
   21   CONTINUE
C
C-----------------------------------------------------------------------
CL 2.  Call routine IMPL_CAL to calculate incrememnts to tracer field
CL     and suface deposition flux for output
C-----------------------------------------------------------------------
C
      CALL IMP_MIX (
     & P_FIELD,P1,P_POINTS,BL_LEVELS,DELTA_AK,DELTA_BK
     &,GAMMA_RHOKH_RDZ,RHOK_DEP
     &,PSTAR,TIMESTEP
     &,F_FIELD,SURF_DEP_FLUX,FIELD
     &,NRML
     &,ERROR,LTIMER
     & )

!

      IF (LTIMER) THEN
        CALL TIMER('TRMIX   ',4)
      ENDIF

      RETURN
      END
