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
!
! Counts the number of points containing each surface type and creates
! a TILE_INDEX array specifying the location of these points on the land
! grid.
!
! Subroutine Interface:

      SUBROUTINE TILEPTS(P_FIELD,LAND_FIELD,LAND1,LAND_PTS,
     &                   FRAC,TILE_PTS,TILE_INDEX
     &                   )


      IMPLICIT NONE
!
! Description:
!
! Method:
!
! Current Code Owner:  Richard Betts
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   4.4    16/10/97   Original code.  Peter Cox
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.

      INTEGER
     & P_FIELD               ! IN Number of P-points in whole grid
     &,LAND_FIELD            ! IN Number of land points in whole grid.
     &,LAND1                 ! IN First land point to be processed.
     &,LAND_PTS              ! IN Number of land points to be processed.

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

      REAL
     & FRAC(LAND_FIELD,NTYPE)       ! IN Fractions of surface types.

      INTEGER
     & TILE_PTS(NTYPE)              ! OUT Number of land points which
C                                   !     include the nth surface type.
     &,TILE_INDEX(LAND_FIELD,NTYPE) ! OUT Indices of land points which
C                                   !     include the nth surface type.
     &,L,N                          ! WORK Loop counters.
C-----------------------------------------------------------------------
C Local parameters
C-----------------------------------------------------------------------
      REAL
     + FRAC_MIN                   ! Minimum ("seed") areal fraction.
      PARAMETER(FRAC_MIN = 0.01)

C-----------------------------------------------------------------------
C Create the TILE_INDEX array of land points with each surface type
C-----------------------------------------------------------------------
      DO N=1,NTYPE
        TILE_PTS(N) = 0
        DO L=LAND1,LAND1+LAND_PTS-1
          TILE_INDEX(L,N)=0
          IF (FRAC(L,N).GT.0.0) THEN
            TILE_PTS(N) = TILE_PTS(N) + 1
            TILE_INDEX(TILE_PTS(N),N) = L
          ENDIF
        ENDDO
      ENDDO

      RETURN
      END
