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
CLL  SUBROUTINE SURF_HYD-----------------------------------------------
CLL
CLL  PURPOSE : TO CARRY OUT CANOPY AND SURFACE HYDROLOGY CALCULATIONS
CLL
CLL            CANOPY WATER CONTENT IS DEPRECIATED BY EVAPORATION
CLL
CLL            SNOWMELT IS RUNOFF THE SURFACE WITHOUT INTERACTING
CLL            WITH THE CANOPY
CLL
CLL            THE CANOPY INTERCEPTION AND SURFACE RUNOFF OF
CLL            LARGE-SCALE RAIN IS CALCUALTED
CLL
CLL            THE CANOPY INTERCEPTION AND SURFACE RUNOFF OF
CLL            CONVECTIVE RAIN IS CALCUALTED
CLL
CLL
CLL  SUITABLE FOR SINGLE COLUMN MODEL USE
CLL
CLL  WRITTEN FOR CRAY-YMP BY S.ALLEN-BETT AND D.GREGORY
CLL
CLL  MODEL            MODIFICATION HISTORY FROM MODEL VERSION 3.0:
CLL VERSION  DATE
CLL
CLL  PROGRAMMING STANDARDS : UNIFIED MODEL DOCUMENTATION PAPER NO. 4
CLL  VERSION NO. 1 18/1/90
CLL
CLL  SYSTEM TASK : P252
CLL
CLL  DOCUMENTATION : UNIFIED MODEL DOCUMENTATION PAPER NO 25
CLL
CLLEND-----------------------------------------------------------------
C
C*L  ARGUMENTS---------------------------------------------------------
C
      SUBROUTINE SURF_HYD (NPNTS,NTYPE,TILE_PTS,TILE_INDEX,
     &                     LICE_PTS,LICE_INDEX,
     &                     CAN_CPY,E_CANOPY,FRAC,INFIL,CON_RAIN,
     &                     LS_RAIN,SNOW_FRAC,SNOW_MELT,TIMESTEP,
     &                     CAN_WCNT,
     &                     CAN_WCNT_GB,DSMC_DT,SURF_ROFF,TOT_TFALL)

      IMPLICIT NONE

      INTEGER
     & NPNTS                ! IN Total number of land points.
     &,NTYPE                ! IN Number of tiles.
     &,TILE_PTS(NTYPE)      ! IN Number of tile points.
     &,TILE_INDEX(NPNTS,NTYPE)
!                           ! IN Index of tile points.
     &,LICE_INDEX(NPNTS)    ! IN Array of land ice points.
     &,LICE_PTS             ! IN Number of land ice points.

      REAL
     & CAN_CPY(NPNTS,NTYPE-1)! IN Canopy capacity for snow-free
!                           !     land tiles (kg/m2).
     &,E_CANOPY(NPNTS,NTYPE-1)
!                           ! IN Snow-free canopy evaporation (kg/m2/s).
     &,FRAC(NPNTS,NTYPE)    ! IN Tile fractions.
     &,INFIL(NPNTS,NTYPE)   ! IN Infiltration rate (kg/m2/s).
     &,CON_RAIN(NPNTS)      ! IN Convective rain (kg/m2/s).
     &,LS_RAIN(NPNTS)       ! IN Large-scale rain (kg/m2/s).
     &,SNOW_FRAC(NPNTS)     ! IN Fraction of gridbox with snow cover.
     &,SNOW_MELT(NPNTS)     ! IN Snow melt (kg/m2/s).
     &,TIMESTEP             ! IN Timestep (s).

      REAL
     & CAN_WCNT(NPNTS,NTYPE-1)
!                           ! INOUT Snow-free tile canopy water contents
!                           !       (kg/m2).

      REAL
     & CAN_WCNT_GB(NPNTS)   ! OUT Gridbox canopy water content (kg/m2).
     &,DSMC_DT(NPNTS)       ! OUT Rate of change of soil moisture
                            !     content (kg/m2/s).
     &,SURF_ROFF(NPNTS)     ! OUT Cumulative surface runoff (kg/m2/s).
     &,TOT_TFALL(NPNTS)     ! OUT Cumulative canopy throughfall
                            !     (kg/m2/s).

!  Workspace -----------------------------------------------------------
      REAL
     & CAN_COND(NPNTS)      ! Canopy condensation (kg/m2/s).

      INTEGER
     & I                    ! Loop counter (land points).
     &,J                    ! Loop counter (tile points).
     &,N                    ! Loop counter (tiles).

      EXTERNAL FRUNOFF,SIEVE

! Zero cumulative stores
      DO I=1,NPNTS
        CAN_WCNT_GB(I) = 0.0
        TOT_TFALL(I) = 0.0
        SURF_ROFF(I) = 0.0
        DSMC_DT(I)   = 0.0
      ENDDO

! Land-ice points
      DO J=1,LICE_PTS
        I = LICE_INDEX(J)
        TOT_TFALL(I) = LS_RAIN(I) + CON_RAIN(I)
        SURF_ROFF(I) = TOT_TFALL(I) + SNOW_MELT(I)
      ENDDO

! Surface runoff of snowmelt
      DO N=1,NTYPE-1
        DO J=1,TILE_PTS(N)
          I = TILE_INDEX(J,N)
          IF ( SNOW_MELT(I) .GT. 0. )
     &      SURF_ROFF(I) = SURF_ROFF(I) + FRAC(I,N)*SNOW_MELT(I) *
     &                     EXP( - SNOW_FRAC(I)*INFIL(I,N)/SNOW_MELT(I) )
!                                                            ... P252.14
        ENDDO
      ENDDO

! Reduce canopy water content by evaporation
      DO N=1,NTYPE-1
        DO J=1,TILE_PTS(N)
        I = TILE_INDEX(J,N)
          IF (E_CANOPY(I,N) .GT. 0.0)
     &      CAN_WCNT(I,N) =
     &                 MAX( CAN_WCNT(I,N) - E_CANOPY(I,N)*TIMESTEP, 0. )
        ENDDO
      ENDDO

      DO N=1,NTYPE-1

! Define canopy condensation when evaporation is negative
        DO J=1,TILE_PTS(N)
        I = TILE_INDEX(J,N)
          IF ( E_CANOPY(I,N) .LT. 0. ) THEN
           CAN_COND(I) = - E_CANOPY(I,N)
          ELSE
           CAN_COND(I) = 0.
          ENDIF
        ENDDO

! Canopy interception, throughfall and surface runoff for condensation,
! assumed to cover 100% of gridbox
        CALL SIEVE (NPNTS,TILE_PTS(N),TILE_INDEX(1,N),1.,
     &              CAN_CPY(1,N),CAN_COND,FRAC(1,N),TIMESTEP,
     &              CAN_WCNT(1,N),TOT_TFALL)
        CALL FRUNOFF (NPNTS,TILE_PTS(N),TILE_INDEX(1,N),1.,
     &                CAN_CPY(1,N),CAN_WCNT(1,N),INFIL(1,N),CAN_COND,
     &                FRAC(1,N),TIMESTEP,
     &                SURF_ROFF)

! Canopy interception, throughfall and surface runoff for large-scale
! rain, assumed to cover 100% of gridbox
        CALL SIEVE (NPNTS,TILE_PTS(N),TILE_INDEX(1,N),1.,
     &              CAN_CPY(1,N),LS_RAIN,FRAC(1,N),TIMESTEP,
     &              CAN_WCNT(1,N),TOT_TFALL)
        CALL FRUNOFF (NPNTS,TILE_PTS(N),TILE_INDEX(1,N),1.,
     &                CAN_CPY(1,N),CAN_WCNT(1,N),INFIL(1,N),LS_RAIN,
     &                FRAC(1,N),TIMESTEP,
     &                SURF_ROFF)

! Canopy interception, throughfall and surface runoff for convective
! rain, assumed to cover 30% of gridbox
        CALL SIEVE (NPNTS,TILE_PTS(N),TILE_INDEX(1,N),0.3,
     &              CAN_CPY(1,N),CON_RAIN,FRAC(1,N),TIMESTEP,
     &              CAN_WCNT(1,N),TOT_TFALL)
        CALL FRUNOFF (NPNTS,TILE_PTS(N),TILE_INDEX(1,N),0.3,
     &                CAN_CPY(1,N),CAN_WCNT(1,N),INFIL(1,N),CON_RAIN,
     &                FRAC(1,N),TIMESTEP,
     &                SURF_ROFF)

        DO I=1,NPNTS
          CAN_WCNT_GB(I) = CAN_WCNT_GB(I) + FRAC(I,N)*CAN_WCNT(I,N)
        ENDDO

      ENDDO

      DO I=1,NPNTS
        DSMC_DT(I) = TOT_TFALL(I) + SNOW_MELT(I) - SURF_ROFF(I)
      ENDDO

      RETURN
      END
