C *****************************COPYRIGHT******************************
C (c) CROWN COPYRIGHT 1997, METEOROLOGICAL OFFICE, All Rights Reserved.
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
!!!  SUBROUTINE IM_CAL_UV ---------------------------------------------
!!!
!!!  Purpose: Calculate increments for U or V in the boundary layer,
!!!           using an implicit numerical scheme.  The tridiagonal
!!!           matrices are inverted using simple Gaussian elimination.
!!!
!!!
!!!  Model           Modification history
!!! version  Date
CLL  4.5    Jul. 98  Kill the IBM specific lines (JCThil)
!!!
!!!
!!! JJ - Programmers of some or all of previous code or changes
!!!
!!!
!!!  Programming standard: UM Documentation Paper No 4, Version 2,
!!!                        dated 18/1/90
!!!
!!!  System component covered: P244
!!!
!!!  Project task: P24
!!!
!!!  Documentation: UM Documentation Paper No 24.
!!!
!!!---------------------------------------------------------------------

!  Arguments :-
      SUBROUTINE IM_CAL_UV (
     & U_V_FIELD,U1_V1
     &,U_V_POINTS,BL_LEVELS,ROW_LENGTH
     &,GAMMA
     &,RHOKM_U_V
     &,U_V,U0_V0,TIMESTEP
     &,RHOKM_1_U_V,DU_NT_DV_NT,DU_DV
     &,DTRDZ_U_V,RDZ_U_V,TAU_X_Y
     &,LTIMER
     &)

      IMPLICIT NONE
      LOGICAL LTIMER

      INTEGER
     & U_V_FIELD                   ! IN No. of points in U_V-grid.
     &,U1_V1                       ! IN First point to be processed in
!                                       U_V-grid.
     &,U_V_POINTS                  ! IN Number of U_V-grid points.
     &,BL_LEVELS                   ! IN No. of atmospheric levels for
!                                       which boundary layer fluxes are
!                                       calculated.
     &,ROW_LENGTH                  ! IN No. of points in latitude row.

      REAL
     & DTRDZ_U_V(U_V_FIELD,BL_LEVELS)
!                                    IN -g.dt/dp for model wind layers
     &,DU_NT_DV_NT(U_V_FIELD,BL_LEVELS)
!                                    IN u_v non-turbulent increments.
     &,GAMMA(BL_LEVELS)            ! IN Implicit weighting coef.
     &,RDZ_U_V(U_V_FIELD,2:BL_LEVELS)
!                                    IN Reciprocal of the vertical
!                                       distance from level K-1 to
!                                       level K. (K > 1) on wind levels
     &,RHOKM_1_U_V(U_V_FIELD)      ! IN Level 1 exchange coefficient for
!                                       momentum
     &,RHOKM_U_V(U_V_FIELD,2:BL_LEVELS)
!                                    IN Exchange coefficients for
!                                       momentum, on UV-grid with
!                                       first and last rows ignored.
!                                       for K>=2 (from KMKH).
     &,U0_V0(U_V_FIELD)            ! IN Westerly_Southerly component of
!                                       surface current
!                                       (m/s; 0 over land) UVG.
     &,U_V(U_V_FIELD,BL_LEVELS)    ! IN Westerly_Southerly component of
!                                       wind.
     &,TIMESTEP                    ! IN Timestep in seconds.

! INOUT
      REAL
     & TAU_X_Y(U_V_FIELD,BL_LEVELS)! INOUT x_y-component of turbulent
!                                        stress at levels k-1/2;
!                                        eg. TAUX(,1) is surface stress.
!                                        UV-grid, 1st and last rows set
!                                        to "missing data". (N/sq m)
!                                        IN as "explicit" fluxes from
!                                        ex_flux_uv, OUT as "implicit

!OUT
      REAL
     & DU_DV(U_V_FIELD,BL_LEVELS)  ! OUT delta (U or V) elements of
!                                        vector on RHS, then LHS, of
!                                        eqn P244.80.


!  External references :-
      EXTERNAL TIMER


! Workspace :-
!   6*BL_LEVELS + 2 blocks of real workspace are required.
      REAL
     & AQ_AM(U_V_FIELD,BL_LEVELS)  ! As AQ: "Q" elements of matrix in
!                                    eqn P244.79 (modified during
!                                    Gaussian elimination process).
!                                    As AM: elements of matrix in eqn
!                                    P244.80 (also get modified).

!  Local scalars :-
      REAL
     & CM       ! Matrix element in eqn P244.80.
     &,RBM      ! Reciprocal of BM(') (eqns P244.81, 85, 89).

      INTEGER
     & BLM1     ! BL_LEVELS minus 1.
     &,I        ! Loop counter (horizontal field index).
     &,K        ! Loop counter (vertical index).

!-----------------------------------------------------------------------
!!  0.  Check that the scalars input to define the grid are consistent.
!       See comments to routine SF_EXCH for details.
!-----------------------------------------------------------------------

      IF (LTIMER) THEN
        CALL TIMER('IMCALUV ',3)
      ENDIF

      BLM1 = BL_LEVELS-1

!-----------------------------------------------------------------------
!! 1.  Solve matrix equation P244.80 for implicit increments to U or V.
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!! 1.1 Initial calculations and "upward sweep".
!! (a) "Surface" fluxes.
!-----------------------------------------------------------------------

        DO I=U1_V1+ROW_LENGTH,U1_V1+U_V_POINTS-ROW_LENGTH-1

!  "Explicit" increments to U(1) and V(1) when there is no rapidly
!  mixing layer or it does not extend beyond the bottom model layer.

        DU_DV(I,1) = DTRDZ_U_V(I,1) * ( TAU_X_Y(I,2) - TAU_X_Y(I,1) )
!                                                            ... P244.67
! cjj addition of non-turbulent increments
        DU_DV(I,1) = DU_DV(I,1) + DU_NT_DV_NT(I,1)


        CM = -DTRDZ_U_V(I,1) * GAMMA(1)*RHOKM_1_U_V(I)         ! P244.66
        AQ_AM(I,1) = -DTRDZ_U_V(I,1) * GAMMA(2)*RHOKM_U_V(I,2)
     &               *RDZ_U_V(I,2)                             ! P244.64
        RBM = 1.0 / ( 1.0 - AQ_AM(I,1) - CM )    ! Reciprocal of P244.81
        DU_DV(I,1) = RBM * DU_DV(I,1)                          ! P244.82
        AQ_AM(I,1) = RBM * AQ_AM(I,1)                          ! P244.84

      ENDDO ! loop over U_V_POINTS

!-----------------------------------------------------------------------
!! (b) Fluxes at (or rows representing) layer interfaces above the
!!     surface but below the top of the boundary layer.
!-----------------------------------------------------------------------

      DO K=2,BLM1

        DO I=U1_V1+ROW_LENGTH,U1_V1+U_V_POINTS-ROW_LENGTH-1


          DU_DV(I,K) = DTRDZ_U_V(I,K) *
     &                   ( TAU_X_Y(I,K+1) - TAU_X_Y(I,K) )     ! P244.74
! cjj addition of non-turbulent increments
          DU_DV(I,K) = DU_DV(I,K) + DU_NT_DV_NT(I,K)

          AQ_AM(I,K) = -DTRDZ_U_V(I,K) * GAMMA(K+1)*RHOKM_U_V(I,K+1)*
     &                    RDZ_U_V(I,K+1)                       ! P244.71
          CM = -DTRDZ_U_V(I,K) * GAMMA(K)*RHOKM_U_V(I,K)*
     &          RDZ_U_V(I,K)                                   ! P244.72
          RBM = 1.0 / ( 1.0 - AQ_AM(I,K) -CM*( 1.0 + AQ_AM(I,K-1) ) )
!                     1                      2                    2 1
!                                              ... Reciprocal of P244.85

          DU_DV(I,K) = RBM * ( DU_DV(I,K) - CM*DU_DV(I,K-1) )
!                                                            ... P244.86

          AQ_AM(I,K) = RBM * AQ_AM(I,K)                        ! P244.88
        ENDDO !loop over u_v_points
      ENDDO ! loop over 2,BLM1


!-----------------------------------------------------------------------
!! (c) Top "boundary" layer; also increment U and V here, as implicit
!!     flux for this layer is got from "upward sweep" alone.
!-----------------------------------------------------------------------

        DO I=U1_V1+ROW_LENGTH,U1_V1+U_V_POINTS-ROW_LENGTH-1


        DU_DV(I,BL_LEVELS) = -DTRDZ_U_V(I,BL_LEVELS) *
     &  TAU_X_Y(I,BL_LEVELS)
!                                                            ... P244.78
! cjj addition of non-turbulent increments
        DU_DV(I,BL_LEVELS) = DU_DV(I,BL_LEVELS)
     &                       + DU_NT_DV_NT(I,BL_LEVELS)

        CM = -DTRDZ_U_V(I,BL_LEVELS) * GAMMA(BL_LEVELS)*
     &        RHOKM_U_V(I,BL_LEVELS)*RDZ_U_V(I,BL_LEVELS)
!                                                            ... P244.76
        RBM = 1.0 / ( 1.0 - CM*( 1.0 + AQ_AM(I,BLM1) ) )

!                                              ... Reciprocal of P244.89

        DU_DV(I,BL_LEVELS) = RBM * ( DU_DV(I,BL_LEVELS)    ! P244.90
     &                                    - CM*DU_DV(I,BLM1) )
      ENDDO

!-----------------------------------------------------------------------
!! 1.2 Complete solution of matrix equation by performing "downward
!!     sweep", then update U and V.
!-----------------------------------------------------------------------

      DO K=BLM1,1,-1

        DO I=U1_V1+ROW_LENGTH,U1_V1+U_V_POINTS-ROW_LENGTH-1


          DU_DV(I,K) = DU_DV(I,K) - AQ_AM(I,K)*DU_DV(I,K+1) ! P244.92
        ENDDO
      ENDDO

!-----------------------------------------------------------------------
!! 2.  Essentially diagnostic calculations, though some values (i.e. the
!!     surface wind stresses) are required by the coupled version of the
!!     model.
!-----------------------------------------------------------------------
!! 2.1 Surface wind stress components.
!!     Pass out the value of RHOKM(,1) in GAMMA(*)_RHOKM_1 to satisfy
!!     STASH. GAMMA(*)_RHOKM_RDZ will contain precisely that on output.
!-----------------------------------------------------------------------

        DO I=U1_V1+ROW_LENGTH,U1_V1+U_V_POINTS-ROW_LENGTH-1


        TAU_X_Y(I,1) = TAU_X_Y(I,1) +
     &                 GAMMA(1)*RHOKM_1_U_V(I)*DU_DV(I,1)   !... P244.61

      ENDDO  !u_v_points

!-----------------------------------------------------------------------
!! 2.2 Wind stress components at layer interfaces above the surface.
!-----------------------------------------------------------------------

      DO K=2,BL_LEVELS
        DO I=U1_V1+ROW_LENGTH,U1_V1+U_V_POINTS-ROW_LENGTH-1


          AQ_AM(I,K) = TAU_X_Y(I,K) +
     &    GAMMA(K) * RHOKM_U_V(I,K) * RDZ_U_V(I,K)    ! P244.61
     &                        *( DU_DV(I,K) - DU_DV(I,K-1) )
          TAU_X_Y(I,K) = AQ_AM(I,K)

        ENDDO !u_v_points
      ENDDO ! bl_levels

      IF (LTIMER) THEN
        CALL TIMER('IMCALUV ',4)
      ENDIF

      RETURN
      END
