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
!!!  Purpose: Calculate explicit flux of momentum in u or v direction
!!!
!!!
!!!  Model           Modification history
!!! version  Date
!!!  4.3    23/2/97   New deck
!!!  4.5    Jul. 98  Replace IBM with SCMA  (JCThil)
!!!
!!!  JJ, SDJ <- Programmers of some or all of previous code or changes
!!!
!!!
!!!  Programming standard: UMDP3
!!!
!!!  System component covered: P244
!!!
!!!  Project task: P24
!!!
!!!  Documentation: UM Documentation Paper No 24.
!!!
!!!---------------------------------------------------------------------

! SUBROUTINE EX_FLUX_UV

      SUBROUTINE EX_FLUX_UV (
     &  U_V_POINTS
     &, U_V_FIELD
     &, ROW_LENGTH
     &, BL_LEVELS
     &, U1_V1
     &, U_V, U_0_V_0
     &, RDZ_U_V
     &, RHOKM_U_V
     &, TAU_X_Y
     &, LTIMER
     &  )

      IMPLICIT NONE

! ARGUMENTS WITH INTENT IN. IE: INPUT VARIABLES.

      LOGICAL
     &  LTIMER                 ! IN Flag for TIMER diagnostics

      INTEGER
     &  U_V_POINTS             ! IN Number of U_V-grid points.
     &, U_V_FIELD              ! IN No. of points in U_V-grid.
     &, ROW_LENGTH             ! IN No. of points in latitude row.
     &, BL_LEVELS              ! IN No. of atmospheric levels for
!                                  which boundary layer fluxes are
!                                  calculated.

     &, U1_V1                  ! IN First point to be processed in
!                                  U_V-grid.


      REAL
     &  RDZ_U_V (U_V_FIELD, 2:BL_LEVELS)
!                                IN Reciprocal of the vertical
!                                   distance from level K-1 to
!                                   level K. (K > 1) on wind levels
     &, RHOKM_U_V (U_V_FIELD, BL_LEVELS)
!                                IN Exchange coefficients for
!                                   momentum, on UV-grid with
!                                   first and last rows ignored.
!                                   for K>=2 (from KMKH).
     &, U_V (U_V_FIELD, BL_LEVELS)
!                                IN Westerly_Southerly component of
!                                   wind.

     &, U_0_V_0 (U_V_FIELD)    ! IN Westerly_Southerly component of
!                                   surface current
!                                   (m/s; 0 over land) UVG.


! ARGUMENTS WITH INTENT OUT. IE: INPUT VARIABLES CHANGED ON OUTPUT.

      REAL
     &  TAU_X_Y (U_V_FIELD, BL_LEVELS)
!                                OUT explicit x_y-component of
!                                    turbulent stress at levels
!                                    k-1/2; eg. TAUX(,1) is surface
!                                    stress. UV-grid, 1st and last rows
!                                    set to "missing data". (N/sq m)

! LOCAL VARIABLES.

      INTEGER
     &  I
     &, J
     &, K
     &, ERROR

! External routines
      EXTERNAL TIMER


      IF (LTIMER) THEN
        CALL TIMER('EXFLUXUV',3)
      ENDIF

!-----------------------------------------------------------------------
!!  0.  Check that the scalars input to define the grid are consistent.
!-----------------------------------------------------------------------

      ERROR = 0

!-----------------------------------------------------------------------
!!  1.  Calculate "explicit" surface fluxes of momentum
!-----------------------------------------------------------------------

! Level 1. Formerly calculated at end of SFEXCH

        DO I=U1_V1+ROW_LENGTH,U1_V1+U_V_POINTS-ROW_LENGTH-1

           TAU_X_Y(I,1) = RHOKM_U_V(I,1) * ( U_V(I,1) - U_0_V_0(I) )

        ENDDO

! Other Levels. Formerly calculated at end of KMKH
      DO K = 2,BL_LEVELS
        DO I=U1_V1+ROW_LENGTH,U1_V1+U_V_POINTS-ROW_LENGTH-1
          TAU_X_Y(I,K) = RHOKM_U_V(I,K) * ( U_V(I,K) - U_V(I,K-1) )
     &    *RDZ_U_V(I,K)

        END DO
      END DO

    6 CONTINUE ! exit error point

      IF (LTIMER) THEN
        CALL TIMER('EXFLUXUV',4)
      ENDIF
      RETURN

      END
