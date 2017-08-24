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
!***********************************************************************
! Calculates the soil respiration based on a simplified version of the
! model of Raich et al. (1991).
!***********************************************************************
      SUBROUTINE MICROBE (LAND_FIELD,LAND_PTS,LAND1
     &,                   CS,STH_SOIL,V_SAT,V_WILT,TSOIL,RESP_S)

      IMPLICIT NONE

      INTEGER
     & LAND_FIELD                 ! IN Total number of land points.
     &,LAND_PTS                   ! IN Number of land points to be
!                                 !    processed.
     &,LAND1                      ! IN First land point to be
!                                 !    processed.

      REAL
     & CS(LAND_FIELD)             ! IN Soil carbon (kg C/m2).
     &,STH_SOIL(LAND_FIELD)       ! IN Top layer soil moisture as a
!                                 !    fraction of saturation (m3/m3).
     &,V_SAT(LAND_FIELD)          ! IN Volumetric soil moisture
!                                 !    concentration at saturation
!                                 !    (m3 H2O/m3 soil).
     &,V_WILT(LAND_FIELD)         ! IN Volumetric soil moisture
!                                 !    concentration below which
!                                 !    stomata close (m3 H2O/m3 soil).
                                  !    as a fraction of saturation.
     &,TSOIL(LAND_FIELD)          ! IN Soil temperature (K).
     &,RESP_S(LAND_FIELD)         ! OUT Soil respiration (kg C/m2/s).
     &,FSTH,FTEMP                 ! WORK Factors describing the
!                                 !      influence of soil moisture and
!                                 !      soil temperature respectively
!                                 !      on the soil respiration.
     &,STH_OPT                    ! WORK Fractional soil moisture at
!                                 !      which respiration is maximum.
     &,STH_WILT                   ! WORK Wilting soil moisture as a
!                                 !      fraction of saturation.
      INTEGER
     & L                          ! Loop counter

!-----------------------------------------------------------------------
! Local parameters
!-----------------------------------------------------------------------
      REAL
     & KAPS                       ! Specific soil respiration rate
!                                 ! at 25 deg ! and optimum soil
!                                 ! moisture (/s).
     &,Q10                        ! Q10 factor for soil respiration.
      PARAMETER (KAPS = 0.5E-8, Q10 = 2.0)


      DO L=LAND1,LAND1+LAND_PTS-1

        IF (V_SAT(L) .GT. 0.0) THEN

          STH_WILT = V_WILT(L) / V_SAT(L)
          STH_OPT = 0.5 * (1 + STH_WILT)

          IF (STH_SOIL(L) .LE. STH_WILT) THEN
            FSTH = 0.2
          ELSEIF (STH_SOIL(L) .GT. STH_WILT .AND.
     &            STH_SOIL(L) .LE. STH_OPT) THEN
            FSTH = 0.2 + 0.8 * ((STH_SOIL(L) - STH_WILT)
     &                        / (STH_OPT - STH_WILT))
          ELSEIF (STH_SOIL(L) .GT. STH_OPT) THEN
            FSTH = 1 - 0.8 * (STH_SOIL(L) - STH_OPT)
          ENDIF

          FTEMP = Q10 ** (0.1 * (TSOIL(L) - 298.15))

          RESP_S(L) = KAPS * CS(L) * FSTH * FTEMP

        ELSE

          RESP_S(L) = 0.0

        ENDIF

      ENDDO

      RETURN

      END
