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
!    SUBROUTINE SMC_EXT-----------------------------------------------
!
! Subroutine Interface:
      SUBROUTINE SMC_EXT (NPNTS,NSHYD,TILE_PTS,TILE_INDEX
     &,                   F_ROOT,FRAC,STHU,V_CRIT,V_SAT,V_WILT
     &,                   WT_EXT,FSMC)

      IMPLICIT NONE
!
! Description:
!     Calculates the soil moisture availability factor and
!     the fraction of the transpiration which is extracted from each
!     soil layer.
!
!
! Documentation : UM Documentation Paper 25
!
! Current Code Owner : David Gregory
!
! History:
! Version   Date     Comment
! -------   ----     -------
!  4.4               New deck   Peter Cox
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!
! System component covered: P25
! System Task: P25
!

! Subroutine arguments:
!   Scalar arguments with intent(IN) :
      INTEGER
     & NPNTS                ! IN Number of gridpoints.
     &,NSHYD                ! IN Number of soil moisture layers.
     &,TILE_PTS             ! IN Number of points containing the
C                           !    given surface type.
     &,TILE_INDEX(NPNTS)    ! IN Indices on the land grid of the
C                           !    points containing the given
C                           !    surface type.


!   Array arguments with intent(IN) :
      REAL
     & F_ROOT(NSHYD)        ! IN Fraction of roots in each soil
!                           !    layer.
     &,FRAC(NPNTS)          ! IN Tile fraction.
     &,STHU(NPNTS,NSHYD)    ! IN Unfrozen soil moisture content of
!                           !    each layer as a fraction of
!                           !    saturation.
     &,V_CRIT(NPNTS)        ! IN Volumetric soil moisture
!                           !    concentration above which
!                           !    evapotranspiration is not sensitive
!                           !    to soil water (m3 H2O/m3 soil).
     &,V_SAT(NPNTS)         ! IN Volumetric soil moisture
!                           !    concentration at saturation
!                           !    (m3 H2O/m3 soil).
     &,V_WILT(NPNTS)        ! IN Volumetric soil moisture
!                           !    concentration below which
!                           !    stomata close (m3 H2O/m3 soil).

!   Array arguments with intent(INOUT) :
      REAL
     & WT_EXT(NPNTS,NSHYD)  ! OUT Cummulative fraction of transpiration
!                           !     extracted from each soil layer
!                           !     (kg/m2/s).


!   Array arguments with intent(OUT) :
      REAL
     & FSMC(NPNTS)          ! OUT Soil moisture availability
!                           !     factor.

! Local scalars:
      INTEGER
     & I,J,N                ! WORK Loop counters

! Local arrays:
      REAL
     & FSMC_L(NPNTS,NSHYD)  ! WORK Soil moisture availability
!                           !      factor for each soil layer.

!----------------------------------------------------------------------
! Initialisations
!----------------------------------------------------------------------
      DO I=1,NPNTS
        FSMC(I)=0.0
      ENDDO

!----------------------------------------------------------------------
! Calculate the soil moisture availability factor for each layer and
! weight with the root fraction to calculate the total availability
! factor.
!----------------------------------------------------------------------
      DO N=1,NSHYD
        DO J=1,TILE_PTS
          I=TILE_INDEX(J)

          FSMC_L(I,N)=(STHU(I,N)*V_SAT(I)-V_WILT(I))
     &               /(V_CRIT(I)-V_WILT(I))
          FSMC_L(I,N)=MAX(FSMC_L(I,N),0.0)
          FSMC_L(I,N)=MIN(FSMC_L(I,N),1.0)

          FSMC(I)=FSMC(I)+F_ROOT(N)*FSMC_L(I,N)

        ENDDO
      ENDDO

!----------------------------------------------------------------------
! Calculate the fraction of the tranpiration which is extracted from
! each soil layer.
!----------------------------------------------------------------------
      DO N=1,NSHYD
        DO J=1,TILE_PTS
          I=TILE_INDEX(J)
          IF (FSMC(I) .GT. 0.0)
     &      WT_EXT(I,N) = WT_EXT(I,N) +
     &                            FRAC(I)*F_ROOT(N)*FSMC_L(I,N)/FSMC(I)
        ENDDO
      ENDDO

      RETURN
      END
