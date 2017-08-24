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
!**********************************************************************
! Routine to calculate the aerodynamic resistance
!
! Written by Peter Cox (June 1997)
!**********************************************************************
      SUBROUTINE RAERO (LAND_FIELD,LAND_INDEX,P_FIELD
     &,                 VEG_PTS,VEG_INDEX
     &,                 RIB,WIND,Z0H,Z0M,Z1,RA)

      IMPLICIT NONE

      INTEGER
     & LAND_FIELD                 ! IN Total number of land points.
     &,LAND_INDEX(LAND_FIELD)     ! IN Index of land points on the
!                                 !    P-grid.
     &,P_FIELD                    ! IN Total number of P-gridpoints.
     &,VEG_PTS                    ! IN Number of vegetated points.
     &,VEG_INDEX(LAND_FIELD)      ! IN Index of vegetated points
!                                 !    on the land grid.
     &,I,J,L                      ! WORK Loop counters.

      REAL
     & RIB(P_FIELD)               ! IN Bulk Richardson number.
     &,WIND(P_FIELD)              ! IN Windspeed (m/s).
     &,Z0H(LAND_FIELD)            ! IN Roughness length for heat (m).
     &,Z0M(LAND_FIELD)            ! IN Roughness length for momentum (m)
     &,Z1(P_FIELD)                ! IN Reference level (m).
     &,RA(LAND_FIELD)             ! OUT Aerodynamic resistance (s/m).
     &,BH                         ! WORK Stability coefficient.
     &,CHN(LAND_FIELD)            ! WORK Neutral drag coefficient.
     &,FH(LAND_FIELD)             ! WORK Stability factor.
     &,ZETAH,ZETAM                ! WORK Tempories in calculation of
!                                 !      CHN.
!-----------------------------------------------------------------------
! Local parameters
!-----------------------------------------------------------------------
      REAL
     & AH,CZ                      ! Stability parameters.
     &,KARMAN                     ! Von Karman's constant.
      PARAMETER (AH = 10.0, CZ = 4.0, KARMAN = 0.4)

!-----------------------------------------------------------------------
! Calculate the neutral bulk tranfer coefficient.
!-----------------------------------------------------------------------
      DO J=1,VEG_PTS
        L = VEG_INDEX(J)
        I = LAND_INDEX(L)
        ZETAM = LOG((Z1(I) + Z0M(L)) / Z0M(L))
        ZETAH = LOG((Z1(I) + Z0M(L)) / Z0H(L))
        CHN(L) = (KARMAN * KARMAN) / (ZETAH * ZETAM)
      ENDDO

!-----------------------------------------------------------------------
! Calculate the stability factor.
!-----------------------------------------------------------------------
      DO J=1,VEG_PTS
        L = VEG_INDEX(J)
        I = LAND_INDEX(L)
        BH = AH * CHN(L) * CZ * SQRT (Z1(I) / Z0H(L))
        IF (RIB(I) .GE. 0.0) THEN
          FH(L) = 1.0 / (1 + AH * RIB(I))
        ELSE
          FH(L) = 1 - AH * RIB(I) / (1 + BH * SQRT(-RIB(I)))
        ENDIF
      ENDDO

!-----------------------------------------------------------------------
! Calculate the aerodynamic resistance.
!-----------------------------------------------------------------------
      DO J=1,VEG_PTS
        L = VEG_INDEX(J)
        I = LAND_INDEX(L)
        RA(L) = 1.0 / (FH(L) * CHN(L) * WIND(I))
      ENDDO

      RETURN
      END
