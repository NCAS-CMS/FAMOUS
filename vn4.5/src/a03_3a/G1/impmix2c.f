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
C*LL  SUBROUTINE IMP_MIX -----------------------------------------------
CLL
CLL  Purpose: Calculate turbulent mixing increments for a passive tracer
CLL           using an implicit numerical scheme.  The tridiagonal
CLL           matices are inverted using simple Gaussian elimination.
CLL
CLL
CLL  Model           Modification history :
CLL version  Date
CLL   3.4  18/10/94   *DECK inserted into UM version 3.4. S Jackson
CLL   4.2   Oct. 96   T3E migration - *DEF CRAY removed
CLL                                     S J Swarbrick
CLL  4.5    Jul. 98  Kill the IBM specific lines (JCThil)
CLL
CLL  SDJ  <- Programmers of some or all of previous code or changes
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
      SUBROUTINE IMP_MIX (
     & P_FIELD,P1,P_POINTS,BL_LEVELS,DELTA_AK,DELTA_BK
     &,GAMMA_RHOKH_RDZ,GAMMA_RHOK_DEP
     &,PSTAR,TIMESTEP
     &,F_FIELD,SURF_DEP_FLUX,FIELD
     &,NRML
     &,ERROR,LTIMER
     & )
      IMPLICIT NONE
C
C  Inputs :-
C
      INTEGER
     & P_FIELD                     ! IN No. of points in P-grid.
     &,P1                          ! IN First point to be processed in
C                                  !    P-grid.
     &,P_POINTS                    ! IN Number of P-grid points to be
C                                  !    processed.
     &,BL_LEVELS                   ! IN No. of atmospheric levels for
C                                  !    which boundary layer fluxes are
C                                  !    calculated.
      REAL
     & DELTA_AK(BL_LEVELS)         ! IN Difference of hybrid 'A' across
C                                  !    boundary layers (K-1/2 to K+1/2)
C                                  !    (upper minus lower).
     &,DELTA_BK(BL_LEVELS)         ! IN Difference of hybrid 'B' across
C                                  !    boundary layers (K-1/2 to K+1/2)
C                                  !    (upper minus lower).
     &,GAMMA_RHOKH_RDZ(P_FIELD,2:BL_LEVELS)
C                                  ! IN Turbulent mixing coefs. above
C                                  !    surface, =GAMMA(K)*RHOKH(,K)
C                                  !    *RDZ(K) for K>=2 (from KMKH).
     &,GAMMA_RHOK_DEP(P_FIELD)     ! IN Surface exchange coefficient
C                                  !    for surface deposition*GAMMA(1)
     &,PSTAR(P_FIELD)              ! IN Surface pressure (Pa).
     &,TIMESTEP                    ! IN Timestep in seconds.
C
C  Next 2 arrays are IN as "explicit" fluxes and OUT as "implicit"
C  fluxes.
C
      REAL
     & F_FIELD(P_FIELD,BL_LEVELS)  ! INOUT Flux of tracer
     &,SURF_DEP_FLUX(P_FIELD)      ! INOUT surface deposition flux
     &,FIELD(P_FIELD,BL_LEVELS)    ! INOUT Amount of tracer

      INTEGER
     & NRML(P_FIELD)               ! IN The number of model layers
C                                  !    in the unstable rapidly mixing
C                                  !    layer. Zero if surface layer
C                                  !    is stable.
     &,ERROR                       ! OUT 1 if bad arguments, else 0.

      LOGICAL
     & LTIMER                      ! IN Logical switch for TIMER diagnos

C*
C*L  External references :-

      EXTERNAL TIMER

C*
C*L  Local and other symbolic constants :-
C*L------------------COMDECK C_G----------------------------------------
C G IS MEAN ACCEL DUE TO GRAVITY AT EARTH'S SURFACE
      REAL G

      PARAMETER(G=9.80665)
C*----------------------------------------------------------------------
C*
C*L Workspace :-
C   4*BL_LEVELS + 4 blocks of real workspace are required.
      REAL
     & AF(P_FIELD,BL_LEVELS)       ! Elements in rows in matrix
C                                  ! equation (modified during
C                                  ! Gaussian elimination calculations).
     &,AF_RML(P_FIELD)             ! Matrix element for field in
C                                  ! rapidly mixing layer.
     &,DELTAP(P_FIELD,BL_LEVELS)   ! Vertical pressure difference across
C                                  ! hybrid layers (upper minus lower)
C                                  ! (Pa).
     &,DELTAP_RML(P_FIELD)         ! Vertical pressure difference across
C                                  ! the rapidly mixing layer (Pa).
     &,D_FIELD(P_FIELD,BL_LEVELS)  ! Delta FIELD (tracer field)
C                                  ! elements of vector on RHS, then
C                                  ! LHS, of eqn P244.79.
     &,DFLD_RML(P_FIELD)           ! Delta FIELD for rapidly
C                                  ! mixing layer.
     &,DTRDZ(P_FIELD,BL_LEVELS)    ! -g.dt/dp for bottom BL_LEVELS
C                                  ! model layers (needed in P245).
     &,DTRDZ_RML(P_FIELD)          ! -g.dt/dp for the rapidly
C                                  ! mixing layer (needed in P245).
C*
C  Local scalars :-
      REAL
     & CF       ! Matrix element for local increments to rml
     &,CF_RML   ! As above but for rapidly mixing layer increment.
     &,RBF      ! Reciprocal of B for local increments to rml
     &,RBF_RML  ! As above but for the rapidly mixing layer increment.
     &,DTIMEG   ! TIMESTEP * G (used in several places).
      INTEGER
     & BLM1     ! BL_LEVELS minus 1.
     &,NRMLP1   ! NRML plus 1.
     &,I        ! Loop counter (horizontal field index).
     &,J        ! Offset version of I.
     &,K        ! Loop counter (vertical index).
     &,KM1      ! K minus 1.
     &,KP1      ! K plus 1.
C
      IF (LTIMER) THEN
        CALL TIMER('IMPMIX  ',3)
      ENDIF

      ERROR=0

      DTIMEG = TIMESTEP * G
      BLM1 = BL_LEVELS-1

C
CL----------------------------------------------------------------------
CL (A) Calculations on P-grid.
CL----------------------------------------------------------------------
C
C-----------------------------------------------------------------------
CL 1.  Calculate pressure across layer (in hybrid coordinates), DELTAP,
CL     and then -gdt/dP = dt/rho*dz for use throughout section (A)
C-----------------------------------------------------------------------
C
      DO 1 K=1,BL_LEVELS
        DO 11 I=P1,P1+P_POINTS-1
          DELTAP(I,K) = DELTA_AK(K) + PSTAR(I)*DELTA_BK(K)
          DTRDZ(I,K) = -DTIMEG / DELTAP(I,K)
   11   CONTINUE
   1  CONTINUE
C
C-----------------------------------------------------------------------
CL 2.  Calculate implicit FIELD increments due to local mixing within
CL     the rapidly mixing layer (where it exists).
CL     The surface fluxes F_FIELD(I,1) and SURF_DEP_FLUX(I) are used for
CL     calculating the rapidly mixing layer increments so not here.
CL     Therefore the matrix equation we must solve to find the implicit
CL     FIELD increments due to local mixing within the rml does not
CL     have a "surface" row.
C-----------------------------------------------------------------------
C
C-----------------------------------------------------------------------
CL 2.1 Start 'upward sweep' with lowest model layer, which will be the
CL     bottom of the rapidly mixing layer (rml) if it exists.
C-----------------------------------------------------------------------
C
      DO 21 I=P1,P1+P_POINTS-1
        IF (NRML(I) .GE. 2) THEN
C
C  "Explicit" increments due to local mixing within the rml.
C
          D_FIELD(I,1) = -DTRDZ(I,1) * F_FIELD(I,2)
C
C  Define matrix elements (CF always zero for this case).
C
          AF(I,1) = -DTRDZ(I,1) * GAMMA_RHOKH_RDZ(I,2)
          RBF = 1.0 / ( 1.0 - AF(I,1) )
C
C  Now start Gaussian elimination
C
          D_FIELD(I,1) = RBF * D_FIELD(I,1)
          AF(I,1) = RBF * AF(I,1)
C
C  Start calculating DELTAP_RML. Mid-level depths added in 2.2 below.
C
          DELTAP_RML(I) = DELTAP(I,1)
        ELSE ! No rapidly mixing layer calculations
          DTRDZ_RML(I) = 1.E30
          DFLD_RML(I) = 1.E30
          AF_RML(I) = 1.E30
          DELTAP_RML(I) = 1.E30
        ENDIF
   21 CONTINUE
C
C-----------------------------------------------------------------------
CL 2.2 Continue upward sweep through middle of the rapidly mixing layer
CL     (if it exists) and to its top. NB NRML is always < or = BLM1.
C-----------------------------------------------------------------------
C
      DO 22 K=2,BLM1
        KP1 = K+1
        KM1 = K-1
        DO 221 I=P1,P1+P_POINTS-1
C
C   If in the top rapidly mixing layer then do not include flux at its
C   top in the calculation, i.e. F_FIELD(I,NRML+1) is not
C   included here; it is taken account of in the non-local mixing
C   through the "rapidly mixing layer".
C
          IF ( K .EQ. NRML(I) ) THEN
C
C   Add final DELTAP contribution to DELTAP_RML and then calculate
C   DTRDZ_RML.  Lower level contributions added in 2.1 and below.
C
            DELTAP_RML(I) = DELTAP_RML(I) + DELTAP(I,K)
            DTRDZ_RML(I) = -DTIMEG / DELTAP_RML(I)
C
C  "Explicit" flux divergence across layer giving explicit
C  increment due to the local mixing at the top of rml.
C
            D_FIELD(I,K) = DTRDZ(I,K) * F_FIELD(I,K)
C
C  Define matrix elements (AF always zero for this case).
C
            CF = -DTRDZ(I,K) * GAMMA_RHOKH_RDZ(I,K)
            RBF = 1.0 / ( 1.0 - CF*( 1.0 + AF(I,KM1) ) )
C
C  Now start Gaussian elimination
C
            D_FIELD(I,K) = RBF * ( D_FIELD(I,K)
     &                                   - CF*D_FIELD(I,KM1) )
          ELSEIF (K .LT. NRML(I)) THEN
C
C  Add layer depths to form total rml depth.
C
            DELTAP_RML(I) = DELTAP_RML(I) + DELTAP(I,K)
C
C  "Explicit" flux divergence across layer giving explicit
C  increment due to the local mixing.
C
            D_FIELD(I,K) = -DTRDZ(I,K) * (F_FIELD(I,KP1) - F_FIELD(I,K))
C
C  Define matrix elements.
C
            AF(I,K) = -DTRDZ(I,K) * GAMMA_RHOKH_RDZ(I,KP1)
            CF = -DTRDZ(I,K) * GAMMA_RHOKH_RDZ(I,K)
            RBF = 1.0 / ( 1.0 - AF(I,K) - CF * ( 1.0 + AF(I,KM1) ) )
C
C  Now start Gaussian elimination
C
            D_FIELD(I,K) = RBF * ( D_FIELD(I,K) - CF*D_FIELD(I,KM1) )
            AF(I,K) = RBF * AF(I,K)
          ENDIF
  221   CONTINUE
   22 CONTINUE
C
C-----------------------------------------------------------------------
CL 2.3 Downward sweep through matrix. Add implicit increments due to
CL     local mixing within the rapidly mixing layer.  Update surface
CL     deposition flux and the top-of-rml tracer flux using
CL     local mixing increments for layers 1 and NRML respectively.
C-----------------------------------------------------------------------
C
      DO 23 K=BLM1,1,-1
        KP1 = K + 1
        DO 231 I=P1,P1+P_POINTS-1
          IF ((NRML(I) .GE. 2) .AND. (K .EQ. NRML(I))) THEN
            FIELD(I,K) = FIELD(I,K) + D_FIELD(I,K)
            F_FIELD(I,KP1) = F_FIELD(I,KP1)
     &                          + GAMMA_RHOKH_RDZ(I,KP1)*D_FIELD(I,K)
          ELSEIF ((NRML(I) .GE. 2) .AND. (K .LT. NRML(I))) THEN
            D_FIELD(I,K) = D_FIELD(I,K)
     &                           - AF(I,K)*D_FIELD(I,KP1)
            FIELD(I,K) = FIELD(I,K) + D_FIELD(I,K)
          ENDIF
          IF ((NRML(I) .GE. 2) .AND. (K .EQ. 1)) THEN
            SURF_DEP_FLUX(I) = SURF_DEP_FLUX(I)
     &                           - GAMMA_RHOK_DEP(I) * D_FIELD(I,1)
          ENDIF
  231   CONTINUE
   23 CONTINUE
C
C-----------------------------------------------------------------------
CL 4.  Calculate those matrix and vector elements on the LHS of eqn
CL     which are to do with implicit solution of the tracer
CL     transport problem at the surface, above the rml (if it exists)
CL     and between all levels if it does not.
CL     Begin with "upward sweep" through lower half of matrix).
C-----------------------------------------------------------------------
C
      DO 31 I=P1,P1+P_POINTS-1
C-----------------------------------------------------------------------
CL 4.2 Lowest atmospheric layer FIELD row of matrix.
C-----------------------------------------------------------------------
        IF (NRML(I) .GE. 2) THEN
          NRMLP1 = NRML(I) + 1
C
C  "Explicit" rapidly mixing layer increment for FIELD.
C
          DFLD_RML(I) = -DTRDZ_RML(I) *
     &                 (F_FIELD(I,NRMLP1) - F_FIELD(I,1)
     &                                    - SURF_DEP_FLUX(I))
          AF_RML(I) = -DTRDZ_RML(I) * GAMMA_RHOKH_RDZ(I,NRMLP1)
          CF_RML = -DTRDZ_RML(I) * GAMMA_RHOK_DEP(I)
          RBF_RML = 1.0 / ( 1.0 - AF_RML(I) - CF_RML )
          DFLD_RML(I) = RBF_RML * DFLD_RML(I)
          AF_RML(I) = RBF_RML * AF_RML(I)
        ELSE
C
C "Explicit" increment to FIELD(1) when there is no rapidly mixing layer
C  or it does not extend beyond the bottom model layer.
C
          D_FIELD(I,1) = -DTRDZ(I,1) *
     &                   ( F_FIELD(I,2) - F_FIELD(I,1)
     &                                  - SURF_DEP_FLUX(I) )

          CF = -DTRDZ(I,1) * GAMMA_RHOK_DEP(I)
          AF(I,1) = -DTRDZ(I,1) * GAMMA_RHOKH_RDZ(I,2)
          RBF = 1.0 / ( 1.0 - AF(I,1) - CF )
          D_FIELD(I,1) = RBF * D_FIELD(I,1)
          AF(I,1) = RBF * AF(I,1)
        ENDIF
   31 CONTINUE
C-----------------------------------------------------------------------
CL 4.3 Rows of matrix applying to FIELD transport into model layers in
CL     the "middle" of the "boundary" layer, i.e. all but the bottom
CL     layer and the top "boundary" layer.
C-----------------------------------------------------------------------
      DO 43 K=2,BLM1
        KP1 = K+1
        KM1 = K-1
        DO 431 I=P1,P1+P_POINTS-1
C
C   "Explicit" flux divergence across layer giving explicit FIELD
C   increment due to mixing above rml if it exists or for all levels if
C   it does not.
C
          NRMLP1 = NRML(I) + 1
          IF (K .GT. NRML(I)) THEN
            D_FIELD(I,K) = -DTRDZ(I,K) * (F_FIELD(I,KP1) - F_FIELD(I,K))
            AF(I,K) = -DTRDZ(I,K) * GAMMA_RHOKH_RDZ(I,KP1)
            CF = -DTRDZ(I,K) * GAMMA_RHOKH_RDZ(I,K)
            IF ((NRML(I) .GE. 2) .AND. (K .EQ. NRMLP1)) THEN
              RBF = 1.0 / ( 1.0 - AF(I,K)
     &                            - CF * ( 1.0 + AF_RML(I) ) )
              D_FIELD(I,K) = RBF * ( D_FIELD(I,K)
     &                                    - CF*DFLD_RML(I) )
            ELSE
              RBF = 1.0 / ( 1.0 - AF(I,K)
     &                            - CF * ( 1.0 + AF(I,KM1) ) )
              D_FIELD(I,K) = RBF * ( D_FIELD(I,K)
     &                                    - CF*D_FIELD(I,KM1) )
            ENDIF
            AF(I,K) = RBF * AF(I,K)
          ENDIF
  431   CONTINUE
   43 CONTINUE
C-----------------------------------------------------------------------
CL 4.4 Top "boundary" layer FIELD row of matrix. FIELD for this layer
CL     can then be, and is, updated.
C-----------------------------------------------------------------------
      DO 44 I=P1,P1+P_POINTS-1
        D_FIELD(I,BL_LEVELS) = DTRDZ(I,BL_LEVELS) * F_FIELD(I,BL_LEVELS)
C
        CF = -DTRDZ(I,BL_LEVELS) * GAMMA_RHOKH_RDZ(I,BL_LEVELS)
        IF (NRML(I) .EQ. BLM1) THEN
          RBF = 1.0 / ( 1.0 - CF*( 1.0 + AF_RML(I) ) )
          D_FIELD(I,BL_LEVELS) = RBF * ( D_FIELD(I,BL_LEVELS)
     &                                        - CF*DFLD_RML(I) )
        ELSE
          RBF = 1.0 / ( 1.0 - CF*( 1.0 + AF(I,BLM1) ) )
          D_FIELD(I,BL_LEVELS) = RBF * ( D_FIELD(I,BL_LEVELS)
     &                                        - CF*D_FIELD(I,BLM1) )
        ENDIF
        FIELD(I,BL_LEVELS) = FIELD(I,BL_LEVELS) + D_FIELD(I,BL_LEVELS)
   44 CONTINUE
C
C-----------------------------------------------------------------------
CL 5.  "Downward sweep" through whole matrix.  FIELD is updated when
CL     the final implicit increments have been calculated.
C-----------------------------------------------------------------------
CL 5.1 Remaining FIELD rows of matrix and add implicit increments above
CL     the rml or at all levels if it is less than two layers thick.
C-----------------------------------------------------------------------
C
      DO 51 K=BLM1,1,-1
        DO 511 I=P1,P1+P_POINTS-1
          IF ( (K .GT. NRML(I)) .OR. (NRML(I) .LT. 2) ) THEN
            D_FIELD(I,K) = D_FIELD(I,K) - AF(I,K)*D_FIELD(I,K+1)
            FIELD(I,K) = FIELD(I,K) + D_FIELD(I,K)
          ENDIF
  511   CONTINUE
   51 CONTINUE
C
C-----------------------------------------------------------------------
CL 5.2 Rapidly mixing layer increments
C-----------------------------------------------------------------------
      DO 52 I=P1,P1+P_POINTS-1
        IF ( NRML(I) .GE. 2 ) THEN
          NRMLP1 = NRML(I) + 1
          DFLD_RML(I) = DFLD_RML(I) - AF_RML(I) * D_FIELD(I,NRMLP1)
          FIELD(I,1) = FIELD(I,1) + DFLD_RML(I)
        ENDIF
   52 CONTINUE
C
C-----------------------------------------------------------------------
CL 5.4 Add implicit increments due to rapid mixing if in rapid mixing
CL     layer.
C-----------------------------------------------------------------------
      DO 54 K=2,BL_LEVELS
        DO 541 I=P1,P1+P_POINTS-1
          IF (K .LE. NRML(I)) THEN
C
C  Add the increments due to rapid mixing if in the rapidly mixing layer
C
            FIELD(I,K) = FIELD(I,K) + DFLD_RML(I)
          ENDIF
  541   CONTINUE
   54 CONTINUE
C
C-----------------------------------------------------------------------
CL 6.  Calculate final implicit flux of tracer.
C-----------------------------------------------------------------------
CL 6.1 Surface fluxes.
C-----------------------------------------------------------------------
C
      DO 61 I=P1,P1+P_POINTS-1
        IF ( NRML(I) .GE. 2 ) THEN
          SURF_DEP_FLUX(I) = SURF_DEP_FLUX(I) - GAMMA_RHOK_DEP(I) *
     &                                                       DFLD_RML(I)
        ELSE
          SURF_DEP_FLUX(I) = SURF_DEP_FLUX(I) - GAMMA_RHOK_DEP(I) *
     &                                                      D_FIELD(I,1)
        ENDIF
   61 CONTINUE
C
C-----------------------------------------------------------------------
CL 6.2 Fluxes at layer interfaces above the surface.
C-----------------------------------------------------------------------
      DO 62 K=2,BL_LEVELS
        KM1 = K-1
        DO 621 I=P1,P1+P_POINTS-1
C
C  Calculate and store fluxes due to local mixing.
C  F_FIELD(local mixing) stored in array AF.
C
          NRMLP1 = NRML(I) + 1
          IF ((NRML(I) .GE. 2) .AND. (K .EQ. NRMLP1)) THEN
            AF(I,K) = F_FIELD(I,K) - GAMMA_RHOKH_RDZ(I,K)
     &                              * ( D_FIELD(I,K) - DFLD_RML(I) )
          ELSE
            AF(I,K) = F_FIELD(I,K) - GAMMA_RHOKH_RDZ(I,K)
     &                              * ( D_FIELD(I,K) - D_FIELD(I,KM1) )
          ENDIF
C
C  Now calculate the implicit fluxes including both local mixing and
C  if appropriate also the fluxes due to rapid mixing through layers.
C
          IF ( K .EQ. 2 ) THEN
            IF ( NRML(I) .GE. 2 ) THEN
              F_FIELD(I,K) = AF(I,K)
     &                         + F_FIELD(I,KM1) + SURF_DEP_FLUX(I)
     &                         - DFLD_RML(I) / DTRDZ(I,KM1)
            ELSE
              F_FIELD(I,K) = AF(I,K)
            ENDIF
          ELSEIF ( K .LE. NRML(I) ) THEN
            F_FIELD(I,K) = AF(I,K) - AF(I,KM1)
     &                     + F_FIELD(I,KM1) - DFLD_RML(I) / DTRDZ(I,KM1)
          ELSE
            F_FIELD(I,K) = AF(I,K)
          ENDIF
  621   CONTINUE
   62 CONTINUE
C
  999 CONTINUE   ! Branch for error exit.

      IF (LTIMER) THEN
        CALL TIMER('IMPMIX  ',4)
      ENDIF

      RETURN
      END
