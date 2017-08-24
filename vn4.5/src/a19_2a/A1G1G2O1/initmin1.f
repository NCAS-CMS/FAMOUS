C *****************************COPYRIGHT******************************
C (c) CROWN COPYRIGHT 1998, METEOROLOGICAL OFFICE, All Rights Reserved.
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
! Ensures that PFT fractions are greater than non-zero minimum fraction
!
! Subroutine Interface:
      SUBROUTINE INIT_MIN(LAND_FIELD,LAND1,LAND_PTS,FRAC,CS)

      IMPLICIT NONE
!
! Description:
!   If fractions of any PFTs are less than a non-zero minimum fraction
!   on land points that are not entirely (or mostly) covered by ice,
!   water or urban, initialise the PFT fractions to the minimum fraction
!   and take the excess proportionally from other PFTs and soil to
!   ensure that the total fractional cover of all PFTs + soil remains
!   unchanged.
!
! Method:
!   For PFTs with fraction < minimum fraction, reset fraction to minimum
!   fraction and find the total increase for all PFTs.  For all PFTS,
!   define "available fraction" as the difference between fraction
!   and minimum fraction, and find "available fraction" from sum
!   of all PFT available fractions plus fraction of soil (this is the
!   "available fraction" for soil; the minimum fraction for soil is
!   zero).  Reduce fractions of all PFTs and soil by amounts weighted
!   by "available fraction" / "total available fraction" such that the
!   sum of the reductions equals the total increase made earlier.  On
!   points with insufficent veg or soil to do this, take no action as
!   vegetation will not be modelled on these points.
!
!
! Current Code Owner: Richard Betts
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   4.5    8/5/98    Original code.  Richard Betts
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!
      INTEGER
     + NNVG                       ! Number of non-vegetation surface
C                                 ! types.
     +,NPFT                       ! Number of plant functional types.
     +,NTYPE                      ! Number of surface types.
     +,SOIL                       ! Index of the surface type 'Soil'
      PARAMETER (NNVG=4, NPFT=5, NTYPE=9, SOIL=8)
C                                 ! Land surface types :
C                                 !     1 - Broadleaf Tree
C                                 !     2 - Needleleaf Tree
C                                 !     3 - C3 Grass
C                                 !     4 - C4 Grass
C                                 !     5 - Shrub
C                                 !     6 - Urban
C                                 !     7 - Water
C                                 !     8 - Soil
C                                 !     9 - Ice

! Subroutine arguments
!   Scalar arguments with intent(in):
      INTEGER
     & LAND_FIELD          ! IN Total number of land points.
     &,LAND1               ! IN First land point to be processed.
     &,LAND_PTS            ! IN Number of land point to be processed.

      REAL
     & CS(LAND_FIELD)         ! INOUT Soil carbon content (kg C/m2).
     &,FRAC(LAND_FIELD,NTYPE) ! INOUT Fractions of surface types.

      INTEGER
     & L                   ! Loop counter for land points
     &,N                   ! Loop counter for surface types

      REAL
     & FRAC_AVAIL(LAND_FIELD,NTYPE)! LOCAL The part of FRAC that is
!                                  !       available for "donation"
     &,TOT_FRAC_NEED(LAND_FIELD)   ! LOCAL Total fraction needed to make
!                                  !       PFT fractions up to minimum
     &,TOT_FRAC_AVAIL(LAND_FIELD)  ! LOCAL Total fractional area
!                                  !       available to give to PFTs
!                                  !       with less than minimum frac.

!----------------------------------------------------------------------
! Local parameters
!----------------------------------------------------------------------
      REAL
     + CS_MIN                     ! Minimum soil carbon (kg C/m2).
      PARAMETER(CS_MIN = 1.0E-6)
      REAL
     + FRAC_MIN                   ! Minimum ("seed") areal fraction.
      PARAMETER(FRAC_MIN = 0.01)

      DO L=LAND1,LAND1+LAND_PTS-1
        TOT_FRAC_NEED(L) = 0.0
        TOT_FRAC_AVAIL(L) = 0.0

!-----------------------------------------------------------------------
! Find total fraction available for donation to PFTs with less than
! the minimum coverage
!-----------------------------------------------------------------------
        DO N=1,NPFT
          IF (FRAC(L,N).LT.FRAC_MIN) THEN
              TOT_FRAC_NEED(L) = TOT_FRAC_NEED(L) +
     &        (FRAC_MIN - FRAC(L,N))
          ELSE IF (FRAC(L,N).GE.FRAC_MIN) THEN
            FRAC_AVAIL(L,N) = FRAC(L,N)-FRAC_MIN
            TOT_FRAC_AVAIL(L) = TOT_FRAC_AVAIL(L) + FRAC_AVAIL(L,N)
          ENDIF
        ENDDO
        N=SOIL
        FRAC_AVAIL(L,N) = FRAC(L,N)
        TOT_FRAC_AVAIL(L) = TOT_FRAC_AVAIL(L) + FRAC(L,N)

!-----------------------------------------------------------------------
! If sufficient total fraction is available, modify fractions of veg and
! soil and also modify soil carbon.  If insufficient fraction available,
! do neither of these as TRIFFID will not operate on such points.
!-----------------------------------------------------------------------
        IF (TOT_FRAC_AVAIL(L).GE.TOT_FRAC_NEED(L)) THEN

!-----------------------------------------------------------------------
! i)  If PFT fraction is less than the minimum fraction, increase it
!     to the minimum fraction.
!-----------------------------------------------------------------------
          DO N=1,NPFT
            IF (FRAC(L,N).LT.FRAC_MIN) THEN
              FRAC(L,N) = FRAC_MIN
              FRAC_AVAIL(L,N) = 0.0
            ELSEIF (FRAC(L,N).EQ.FRAC_MIN) THEN
              FRAC_AVAIL(L,N) = 0.0
            ENDIF
          ENDDO

!-----------------------------------------------------------------------
! ii) Scale other PFTs and soil to keep total coverage of veg+soil
!     unchanged.  The relative proportions of the soil fraction and the
!     PFT fractions greater than the minimum fraction remain constant.
!-----------------------------------------------------------------------
          DO N=1,NPFT
            FRAC(L,N) = FRAC(L,N) -
     &      ( (FRAC_AVAIL(L,N)/TOT_FRAC_AVAIL(L)) * TOT_FRAC_NEED(L) )
          ENDDO

          N=SOIL
          FRAC(L,N) = FRAC(L,N) -
     &    ( (FRAC_AVAIL(L,N)/TOT_FRAC_AVAIL(L)) * TOT_FRAC_NEED(L) )

!-----------------------------------------------------------------------
! iii) If soil carbon content is less than minimum allowed, increase
!      it to the minimum.
!-----------------------------------------------------------------------
          IF (CS(L).LT.CS_MIN) THEN
            CS(L) = CS_MIN
          ENDIF

        ENDIF

      ENDDO

      RETURN
      END
