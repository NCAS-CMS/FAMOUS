C ******************************COPYRIGHT******************************
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
C
!!!   SUBROUTINE PHI_M_H ----------------------------------------------
!!!
!!!  Purpose: Calculate the integrated froms of the Monin-Obukhov
!!!           stability functions for surface exchanges.
!!!
!!!
!!!  Model            Modification history:
!!! version  Date
!!!
!!!   4.4  12/05/97   *DECK created by R.N.B.Smith
!!!
!!!  Programming standard:
!!!
!!!  System component covered: Part of P243.
!!!
!!!  Documentation: UM Documentation Paper No 24.
!!!
!*L  Arguments:---------------------------------------------------------
      SUBROUTINE PHI_M_H(
     & P_POINTS,P_FIELD,P1,L_LAND,LAND_MASK,
     & RECIP_L_MO,Z_UV,Z_TQ,Z0M,Z0H,PHI_M,PHI_H,LTIMER
     &)
      IMPLICIT NONE

      INTEGER
     & P_POINTS           ! IN Number of gridpoints treated.
     &,P_FIELD            ! IN Size of field on p-grid.
     &,P1                 ! IN First p-point to be treated.

      LOGICAL
     & LTIMER
     &,L_LAND             ! IN If true, treat land points only.
     &,LAND_MASK(P_FIELD) ! IN Land mask


      REAL
     & RECIP_L_MO(P_FIELD)
!                    ! IN Reciprocal of the Monin-Obukhov length (m^-1).
     &,Z_UV(P_FIELD) ! IN Height of wind level above roughness height (m
     &,Z_TQ(P_FIELD) ! IN Height of temperature, moisture and scalar lev
!                    !    above the roughness height (m).
     &,Z0M(P_FIELD)  ! IN Roughness length for momentum (m).
     &,Z0H(P_FIELD)  ! IN Roughness length for heat/moisture/scalars (m)
!
      REAL
     & PHI_M(P_FIELD)! OUT Stability function for momentum.
     &,PHI_H(P_FIELD)! OUT Stability function for heat/moisture/scalars.
!
!*L  Workspace usage----------------------------------------------------
!    No work areas are required.
!
!*----------------------------------------------------------------------
!*L  External subprograms called:

      EXTERNAL TIMER

!*----------------------------------------------------------------------
!  Common and local physical constants.
!
!  None.
!
!  Define local variables.
!
      INTEGER I,L     ! Loop counter; horizontal field index.
!
      REAL
     & PHI_MN         ! Neutral value of stability function for momentum
     &,PHI_HN         ! Neutral value of stability function for scalars.
     &,ZETA_UV        ! Temporary in calculation of PHI_M.
     &,ZETA_0M        ! Temporary in calculation of PHI_M.
     &,ZETA_TQ        ! Temporary in calculation of PHI_H.
     &,ZETA_0H        ! Temporary in calculation of PHI_H.
     &,X_UV_SQ        ! Temporary in calculation of PHI_M.
     &,X_0M_SQ        ! Temporary in calculation of PHI_M.
     &,X_UV           ! Temporary in calculation of PHI_M.
     &,X_0M           ! Temporary in calculation of PHI_M.
     &,Y_TQ           ! Temporary in calculation of PHI_H.
     &,Y_0H           ! Temporary in calculation of PHI_H.

      IF (LTIMER) THEN
        CALL TIMER('PHI_M_H ',3)
      ENDIF

      I=0  ! Reset loop counter

      DO L=P1,P1+P_POINTS-1
!       IF (L_LAND.AND.LAND_MASK(L)) THEN
!         I=I+1
!       ELSEIF (.NOT.L_LAND) THEN
          I=L
!       ENDIF

!       IF(.NOT.L_LAND.OR.(L_LAND.AND.LAND_MASK(L))) THEN

!
!-----------------------------------------------------------------------
!! 1. Calculate neutral values of PHI_M and PHI_H.
!-----------------------------------------------------------------------
!
          PHI_MN = LOG( (Z_UV(I) + Z0M(I)) / Z0M(I) )
          PHI_HN = LOG( (Z_TQ(I) + Z0M(I)) / Z0H(I) )
!
!-----------------------------------------------------------------------
!! 2. Calculate stability parameters.
!-----------------------------------------------------------------------
!
          ZETA_UV = (Z_UV(I) + Z0M(I)) * RECIP_L_MO(I)
          ZETA_TQ = (Z_TQ(I) + Z0M(I)) * RECIP_L_MO(I)
          ZETA_0M = Z0M(I) * RECIP_L_MO(I)
          ZETA_0H = Z0H(I) * RECIP_L_MO(I)
!
!-----------------------------------------------------------------------
!! 3. Calculate PHI_M and PHI_H for neutral and stable conditions.
!-----------------------------------------------------------------------
!
          IF (RECIP_L_MO(I) .GE. 0.0) THEN
            PHI_M(I) = PHI_MN + 4.0 * (ZETA_UV - ZETA_0M)
            PHI_H(I) = PHI_HN +
     &                 (1.0 + 2.0*ZETA_TQ) * (1.0 + 2.0*ZETA_TQ) -
     &                 (1.0 + 2.0*ZETA_0H) * (1.0 + 2.0*ZETA_0H)
!
!-----------------------------------------------------------------------
!! 4. Calculate PHI_M and PHI_H for unstable conditions.
!-----------------------------------------------------------------------
!
          ELSE

            X_UV_SQ = SQRT(1.0 - 16.0*ZETA_UV)
            X_0M_SQ = SQRT(1.0 - 16.0*ZETA_0M)
            X_UV = SQRT(X_UV_SQ)
            X_0M = SQRT(X_0M_SQ)
            PHI_M(I) = PHI_MN - 2.0*LOG( (1.0+X_UV) / (1.0+X_0M) )
     &                      - LOG( (1.0+X_UV_SQ) / (1.0+X_0M_SQ) )
     &                      + 2.0*( ATAN(X_UV) - ATAN(X_0M) )

            Y_TQ = SQRT(1.0 - 16.0*ZETA_TQ)
            Y_0H = SQRT(1.0 - 16.0*ZETA_0H)
            PHI_H(I) = PHI_HN - 2.0*LOG( (1.0+Y_TQ) / (1.0+Y_0H) )

          ENDIF

!       ENDIF ! L_LAND and LAND_MASK

      ENDDO

      IF (LTIMER) THEN
        CALL TIMER('PHI_M_H ',4)
      ENDIF

      RETURN
      END
