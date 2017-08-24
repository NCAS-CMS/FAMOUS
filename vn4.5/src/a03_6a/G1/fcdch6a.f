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
!!!   SUBROUTINE FCDCH--------------------------------------------------
!!!
!!!  Purpose: Calculate surface transfer coefficients at one or more
!!!           gridpoints.
!!!
!!!  Model            Modification history:
!!! version  Date
!!!
!!!   4.4   08/05/97  New version of subroutine using Monin-Obukhov
!!!                   stability functions created.
!!!                                                 R.N.B.Smith
!!!
!!!  Programming standard:
!!!
!!!  System component covered: Part of P243.
!!!
!!!  Documentation: UM Documentation Paper No 24, section P243.
!!!
!!  Arguments:---------------------------------------------------------
      SUBROUTINE FCDCH(
     & P_POINTS,P_FIELD,P1,L_LAND,LAND_MASK,
     & DB,VSHR,Z0M,Z0H,ZH,Z1_UV,Z1_TQ,WIND_PROFILE_FACTOR,
     & CDV,CHV,CDV_STD,V_S,V_S_STD,RECIP_L_MO,LTIMER
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
     & DB(P_FIELD)   ! IN Buoyancy difference between surface and lowest
!                    !    temperature and humidity level in the
!                    !    atmosphere (m/s^2).
     &,VSHR(P_FIELD) ! IN Wind speed difference between the surface and
!                    !    the lowest wind level in the atmosphere (m/s).
     +,Z0M(P_FIELD)  ! IN Roughness length for momentum transport (m).
     +,Z0H(P_FIELD)  ! IN Roughness length for heat and moisture (m).
     +,ZH(P_FIELD)   ! IN Depth of boundary layer (m).
     +,Z1_UV(P_FIELD)! IN Height of lowest wind level (m).
     +,Z1_TQ(P_FIELD)! IN Height of lowest temperature and humidity
!                    !    level (m).
     &,WIND_PROFILE_FACTOR(P_FIELD)
!                    ! IN for adjusting the surface transfer
!                    !    coefficients to remove form drag effects.

      REAL
     & CDV(P_FIELD)  ! OUT Surface transfer coefficient for momentum
!                    !     including orographic form drag (m/s).
     +,CHV(P_FIELD)  ! OUT Surface transfer coefficient for
!                    !     heat, moisture & other scalars (m/s).
     &,CDV_STD(P_FIELD)
!                    ! OUT Surface transfer coefficient for momentum
!                    !     excluding orographic form drag (m/s).
     &,V_S(P_FIELD)  ! OUT Surface layer scaling velocity
!                    !     including orographic form drag (m/s).
     &,V_S_STD(P_FIELD)
!                    ! OUT Surface layer scaling velocity
!                    !     excluding orographic form drag (m/s).
     &,RECIP_L_MO(P_FIELD)
!                    ! OUT Reciprocal of the Monin-Obukhov length
!                    !     (m^-1).

!*L  Workspace usage----------------------------------------------------
!
!     Local work arrays.
!
      REAL
     & PHI_M(P_FIELD) ! Monin-Obukhov stability function for momentum
!                     ! integrated to the model's lowest wind level.
     &,PHI_H(P_FIELD) ! Monin-Obukhov stability function for scalars
!                     ! integrated to the model's lowest temperature
!                     ! and humidity level.
!
!*----------------------------------------------------------------------

      EXTERNAL TIMER,PHI_M_H

!*----------------------------------------------------------------------
!  Common and local constants.
C*L-----------------COMDECK C_VKMAN-------------------------------------
C VKMAN IS VON KARMAN'S CONSTANT
      REAL VKMAN

      PARAMETER(VKMAN=0.4)
C*----------------------------------------------------------------------

      REAL BETA,THIRD
      PARAMETER (
     & BETA=0.08,   ! Tunable parameter in the surface layer scaling
!                   ! velocity formula (multiplying the turbulent
!                   ! convective scaling velocity).
     + THIRD=1./3.  ! One third.
     +)
      INTEGER N_ITS ! Number of iterations for Monin-Obukhov length
!                   ! and stability functions.
      PARAMETER (
     & N_ITS=5
     &)
!
!  Define local variables
!
      INTEGER I,L   ! Loop counter; horizontal field index.
      INTEGER IT    ! Iteration loop counter.

      REAL
     & B_FLUX       ! Surface bouyancy flux over air density.
     &,U_S          ! Surface friction velocity (effective value).
     &,U_S_STD      ! Surface friction velocity (standard value).
     &,W_S          ! Surface turbulent convective scaling velocity.

      IF (LTIMER) THEN
        CALL TIMER('FCDCH   ',3)
      ENDIF

!
!-----------------------------------------------------------------------
!! 1. Set initial values for the iteration.
!-----------------------------------------------------------------------
!
      I=0  ! Reset loop counter

      DO L=P1,P1+P_POINTS-1
!       IF (L_LAND.AND.LAND_MASK(L)) THEN
!         I=I+1
!       ELSEIF (.NOT.L_LAND) THEN
          I=L
!       ENDIF

!       IF(.NOT.L_LAND.OR.(L_LAND.AND.LAND_MASK(L))) THEN

          IF (DB(I) .LT. 0.0 .AND. VSHR(I) .LT. 2.0) THEN
!-----------------------------------------------------------------------
!           Start the iteration from the convective limit.
!-----------------------------------------------------------------------
            RECIP_L_MO(I) = -VKMAN/(BETA*BETA*BETA*ZH(I))
          ELSE
!-----------------------------------------------------------------------
!           Start the iteration from neutral values.
!-----------------------------------------------------------------------
            RECIP_L_MO(I) = 0.0
          ENDIF
!       ENDIF ! L_LAND and LAND_MASK
      ENDDO
      CALL PHI_M_H (P_POINTS,P_FIELD,P1,L_LAND,LAND_MASK,
     &              RECIP_L_MO,Z1_UV,Z1_TQ,Z0M,Z0H,
     &              PHI_M,PHI_H,
     &              LTIMER)
!
      I=0  ! Reset loop counter

      DO L=P1,P1+P_POINTS-1
!       IF (L_LAND.AND.LAND_MASK(L)) THEN
!         I=I+1
!       ELSEIF (.NOT.L_LAND) THEN
          I=L
!       ENDIF

!       IF(.NOT.L_LAND.OR.(L_LAND.AND.LAND_MASK(L))) THEN

          IF (DB(I) .LT. 0.0 .AND. VSHR(I) .LT. 2.0) THEN
!-----------------------------------------------------------------------
!           Start the iteration from the convective limit.
!-----------------------------------------------------------------------
            V_S_STD(I) = BETA *
     &          SQRT( BETA * ( VKMAN / PHI_H(I) ) * ZH(I) * (-DB(I)) )
            V_S(I) = V_S_STD(I)
          ELSE
!-----------------------------------------------------------------------
!           Start the iteration from neutral values.
!-----------------------------------------------------------------------
            V_S(I) = ( VKMAN / PHI_M(I) ) * VSHR(I)
            V_S_STD(I) = V_S(I) * WIND_PROFILE_FACTOR(I)
          ENDIF
          CHV(I) = ( VKMAN / PHI_H(I) ) * V_S_STD(I)
          CDV(I) = ( VKMAN / PHI_M(I) ) * V_S(I)
          CDV_STD(I) = CDV(I) * ( V_S_STD(I) / V_S(I) ) *
     &                          WIND_PROFILE_FACTOR(I)
!       ENDIF ! L_LAND and LAND_MASK
      ENDDO
!-----------------------------------------------------------------------
!! 2. Iterate to obtain sucessively better approximations for CD & CH.
!-----------------------------------------------------------------------
      DO IT = 1,N_ITS
!
        I=0  ! Reset loop counter

        DO L=P1,P1+P_POINTS-1
!         IF (L_LAND.AND.LAND_MASK(L)) THEN
!           I=I+1
!         ELSEIF (.NOT.L_LAND) THEN
            I=L
!         ENDIF

!         IF(.NOT.L_LAND.OR.(L_LAND.AND.LAND_MASK(L))) THEN

            B_FLUX = -CHV(I) * DB(I)
            U_S = SQRT( CDV(I) * VSHR(I) )
            U_S_STD = SQRT( CDV_STD(I) * VSHR(I) )
            IF (DB(I) .LT. 0.0) THEN
              W_S = (ZH(I) * B_FLUX)**THIRD
              V_S(I) = SQRT(U_S*U_S + BETA*BETA*W_S*W_S)
              V_S_STD(I) = SQRT(U_S_STD*U_S_STD + BETA*BETA*W_S*W_S)
            ELSE
              V_S(I) = U_S
              V_S_STD(I) = U_S_STD
            ENDIF
            RECIP_L_MO(I) = -VKMAN * B_FLUX /
     &                       (V_S(I)*V_S(I)*V_S(I))
!         ENDIF ! L_LAND and LAND_MASK
        ENDDO
        CALL PHI_M_H (P_POINTS,P_FIELD,P1,L_LAND,LAND_MASK,
     &                RECIP_L_MO,Z1_UV,Z1_TQ,Z0M,Z0H,
     &                PHI_M,PHI_H,
     &                LTIMER)
!
        I=0  ! Reset loop counter

        DO L=P1,P1+P_POINTS-1
!         IF (L_LAND.AND.LAND_MASK(L)) THEN
!           I=I+1
!         ELSEIF (.NOT.L_LAND) THEN
            I=L
!         ENDIF

!         IF(.NOT.L_LAND.OR.(L_LAND.AND.LAND_MASK(L))) THEN

            CHV(I) = ( VKMAN / PHI_H(I) ) * V_S_STD(I)
            CDV(I) = ( VKMAN / PHI_M(I) ) * V_S(I)
            CDV_STD(I) = CDV(I) * ( V_S_STD(I) / V_S(I) ) *
     &                            WIND_PROFILE_FACTOR(I)
!         ENDIF ! L_LAND and LAND_MASK
        ENDDO
      ENDDO ! Iteration loop

      IF (LTIMER) THEN
        CALL TIMER('FCDCH   ',4)
      ENDIF

      RETURN
      END
