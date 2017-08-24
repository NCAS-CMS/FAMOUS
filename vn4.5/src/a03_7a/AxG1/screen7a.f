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
!!!   SUBROUTINE SCREEN_TQ----------------------------------------------
!!!
!!!  Purpose: Diagnose temperature and/or specific humidity at screen
!!!           height (1.5 metres), as requested via the STASH flags.
!!!           This version outputs gridbox-average diagnostics over
!!!           land surface tiles, but diagnostics for individual tiles
!!!           could be made available as well.
!!!
!!!
!!!  Model            Modification history:
!!! version  Date
!!!
!!!   4.4    8/97     New deck for MOSES II (R. Essery)
!!!
!!!---------------------------------------------------------------------
      SUBROUTINE SCREEN_TQ (
     & P_POINTS,P_FIELD,P1,LAND1,LAND_PTS,LAND_FIELD,NTYPE,
     & LAND_INDEX,TILE_INDEX,TILE_PTS,LAND_MASK,
     & SQ1P5,ST1P5,CHR1P5M,CHR1P5M_SICE,PSTAR,QW_1,RESFT,
     & TILE_FRAC,TL_1,TSTAR,TSTAR_TILE,
     & Z0H,Z0H_TILE,Z0M,Z0M_TILE,Z1,
     & Q1P5M,T1P5M
     & )

      IMPLICIT NONE

      INTEGER
     & P_POINTS             ! IN Number of P-grid points to be
!                           !  processed.
     &,P_FIELD              ! IN Total number of P-grid points.
     &,P1                   ! IN First P-point to be processed.
     &,LAND1                ! IN First land point to be processed.
     &,LAND_PTS             ! IN Number of land points to be processed.
     &,LAND_FIELD           ! IN Total number of land points.
     &,NTYPE                ! IN Number of tiles per land point.
     &,LAND_INDEX(P_FIELD)  ! IN Index of land points.
     &,TILE_INDEX(LAND_FIELD,NTYPE)
!                           ! IN Index of tile points.
     &,TILE_PTS(NTYPE)      ! IN Number of tile points.

      LOGICAL
     & LAND_MASK(P_FIELD)   ! IN T for land points, F otherwise.
     &,SQ1P5                ! IN STASH flag for 1.5-metre sp humidity.
     &,ST1P5                ! IN STASH flag for 1.5-metre temperature.

      REAL
     & CHR1P5M(LAND_FIELD,NTYPE)
!                           ! IN Ratio of coefficients for
!                           !    calculation of 1.5 m T.
     &,CHR1P5M_SICE(P_FIELD)! IN Ratio of coefficients for
!                           !    calculation of 1.5 m T.
     &,PSTAR(P_FIELD)       ! IN Surface pressure (Pa).
     &,QW_1(P_FIELD)        ! IN Total water content of lowest
!                                atmospheric layer (kg per kg air).
     &,RESFT(LAND_FIELD,NTYPE)
!                           ! IN Surface resistance factor.
     &,TILE_FRAC(LAND_FIELD,NTYPE)
!                           ! IN Tile fractions.
     &,TL_1(P_FIELD)        ! IN Liquid/frozen water temperature for
!                                lowest atmospheric layer (K).
     &,TSTAR(P_FIELD)       ! IN Gridbox mean surface temperature (K).
     &,TSTAR_TILE(LAND_FIELD,NTYPE)
!                           ! IN Tile surface temperatures (K).
     &,Z0H(P_FIELD)         ! IN Roughness length for heat and
!                           !    moisture (m).
     &,Z0H_TILE(LAND_FIELD,NTYPE)
!                           ! IN Tile roughness lengths for heat and
!                           !    moisture (m).
     &,Z0M(P_FIELD)         ! IN Roughness length for momentum (m).
     &,Z0M_TILE(LAND_FIELD,NTYPE)
!                           ! IN Tile roughness lengths for momentum (m)
     &,Z1(P_FIELD)          ! IN Height of lowest atmospheric level (m).

      REAL
     & Q1P5M(P_FIELD)       ! OUT Specific humidity at screen height of
!                           !     1.5 metres (kg water per kg air).
     &,T1P5M(P_FIELD)       ! OUT Temperature at screen height of
!                           !     1.5 metres (K).

      REAL
     & CER1P5M              ! Ratio of coefficients reqd for
!                           ! calculation of 1.5 m Q.
     &,PSTAR_LAND(LAND_FIELD)! Surface pressure for land points.
     &,QS(P_FIELD)          ! Surface saturated sp humidity.
     &,QS_TILE(LAND_FIELD)  ! Surface saturated sp humidity.
     &,Q                    ! Local Q at 1.5 m.
     &,T                    ! Local T at 1.5 m.

       INTEGER
     & I           ! Loop counter (horizontal field index).
     &,J           ! Loop counter (tile point index).
     &,L           ! Loop counter (land point field index).
     &,N           ! Loop counter (tile index).

! Local and other symbolic constants used :-

C*L------------------COMDECK C_G----------------------------------------
C G IS MEAN ACCEL DUE TO GRAVITY AT EARTH'S SURFACE
      REAL G

      PARAMETER(G=9.80665)
C*----------------------------------------------------------------------
C*L-----------COMDECK C_HT_M FOR SUBROUTINE SF_EXCH----------
C Z10M  = height of 10m level for diagnostic calculations (m).
C Z1P5M = height of 1.5m level for diagnostic calculations (m).
      REAL Z10M,Z1P5M

      PARAMETER(Z10M  = 10.0,
     &          Z1P5M = 1.5)
C*----------------------------------------------------------------------
C*L------------------COMDECK C_R_CP-------------------------------------
C R IS GAS CONSTANT FOR DRY AIR
C CP IS SPECIFIC HEAT OF DRY AIR AT CONSTANT PRESSURE
C PREF IS REFERENCE SURFACE PRESSURE
      REAL R,CP,KAPPA,PREF  ! PREFTOKI

      PARAMETER(R=287.05,
     &          CP=1005.,
     &          KAPPA=R/CP,
     &          PREF=100000.)
C*----------------------------------------------------------------------

      REAL GRCP
      PARAMETER ( GRCP = G / CP )

!-----------------------------------------------------------------------
! Diagnose local and GBM temperatures at 1.5 m if requested via ST1P5
!-----------------------------------------------------------------------
      IF (ST1P5) THEN

        DO I=P1,P1+P_POINTS-1
          T1P5M(I) = 0.
          IF ( .NOT.LAND_MASK(I) ) THEN
            T1P5M(I) = TSTAR(I) - GRCP*Z1P5M + CHR1P5M_SICE(I) *
     &                 (TL_1(I) - TSTAR(I) + GRCP*(Z1(I)+Z0M(I)-Z0H(I)))
          ENDIF
        ENDDO

        DO N=1,NTYPE
          DO J=1,TILE_PTS(N)
            L = TILE_INDEX(J,N)
            I = LAND_INDEX(L)
            T = TSTAR_TILE(L,N) - GRCP*Z1P5M + CHR1P5M(L,N) *
     &                      ( TL_1(I) - TSTAR_TILE(L,N) +
     &                        GRCP*(Z1(I)+Z0M_TILE(L,N)-Z0H_TILE(L,N)) )
            T1P5M(I) = T1P5M(I) + TILE_FRAC(L,N)*T
          ENDDO
        ENDDO

      ENDIF

!-----------------------------------------------------------------------
! Diagnose local and GBM humidities at 1.5 m if requested via SQ1P5
!-----------------------------------------------------------------------
      IF (SQ1P5) THEN

        CALL QSAT(QS(P1),TSTAR(P1),PSTAR(P1),P_POINTS)
        DO I=P1,P1+P_POINTS-1
          Q1P5M(I) = 0.
          IF ( .NOT.LAND_MASK(I) ) THEN
            CER1P5M = CHR1P5M_SICE(I) - 1.
            Q1P5M(I) = QW_1(I) + CER1P5M*( QW_1(I) - QS(I) )
          ENDIF
        ENDDO

        DO L=LAND1,LAND1+LAND_PTS-1
          I = LAND_INDEX(L)
          PSTAR_LAND(L) = PSTAR(I)
        ENDDO

        DO N=1,NTYPE
          CALL QSAT(QS_TILE(LAND1),TSTAR_TILE(LAND1,N),
     &              PSTAR_LAND(LAND1),LAND_PTS)
          DO J=1,TILE_PTS(N)
            L = TILE_INDEX(J,N)
            I = LAND_INDEX(L)
            CER1P5M = RESFT(L,N)*(CHR1P5M(L,N) - 1.)
            Q = QW_1(I) + CER1P5M*( QW_1(I) - QS_TILE(L) )
            Q1P5M(I) = Q1P5M(I) + TILE_FRAC(L,N)*Q
          ENDDO
        ENDDO

      ENDIF

      RETURN
      END
