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
C*LL  SUBROUTINE SF_RESIST----------------------------------------------
CLL
CLL  Purpose: Calculate surface moisture flux resistance factors.
CLL
CLL
CLLEND-----------------------------------------------------------------
C*
C*L  Arguments --------------------------------------------------------
      SUBROUTINE SF_RESIST (
     & P_FIELD,LAND_FIELD,TILE_PTS,LAND_INDEX,TILE_INDEX,
     & CANOPY,CATCH,CH,DQ,EPDT,GC,VSHR,
     & FRACA,RESFS,RESFT,LTIMER
     & )

      IMPLICIT NONE

      INTEGER
     & P_FIELD             ! IN Total number of P-grid points.
     &,LAND_FIELD          ! IN Total number of land points.
     &,TILE_PTS            ! IN Number of tile points.
     &,LAND_INDEX(P_FIELD )! IN Index of land points.
     &,TILE_INDEX(LAND_FIELD)
!                          ! IN Index of tile points.

      LOGICAL
     & LTIMER              ! IN Logical switch for TIMER diags

      REAL
     & CANOPY(LAND_FIELD)  ! IN Surface water (kg per sq metre).  F642.
     &,CATCH(LAND_FIELD)   ! IN Surface capacity (max. surface water)
!                          !    (kg per sq metre).  F6416.
     &,CH(LAND_FIELD)      ! IN Transport coefficient for heat and
!                          !    moisture transport
     &,DQ(LAND_FIELD)      ! IN Sp humidity difference between surface
!                          !    and lowest atmospheric level (Q1 - Q*).
     &,EPDT(LAND_FIELD)    ! IN "Potential" Evaporation * Timestep.
!                          !    Dummy variable for first call to routine
     &,GC(LAND_FIELD)      ! IN Interactive canopy conductance
!                          !    to evaporation (m/s)
     &,VSHR(P_FIELD)       ! IN Magnitude of surface-to-lowest-level
!                          !    windshear

      REAL
     & FRACA(LAND_FIELD)   ! OUT Fraction of surface moisture flux with
!                          !     only aerodynamic resistance.
     &,RESFS(LAND_FIELD)   ! OUT Combined soil, stomatal and aerodynamic
!                          !     resistance factor for fraction 1-FRACA.
     &,RESFT(LAND_FIELD)   ! OUT Total resistance factor
!                          !     FRACA+(1-FRACA)*RESFS.

! Workspace -----------------------------------------------------------
      INTEGER
     & I           ! Horizontal field index.
     &,J           ! Tile field index.
     &,L           ! Land field index.

      IF (LTIMER) THEN
        CALL TIMER('SFRESIST',3)
      ENDIF

!-----------------------------------------------------------------------
!     Evaporation over land surfaces without snow is limited by
!     soil moisture availability and stomatal resistance.
!     Set FRACA (= fA in the documentation) according to P243.68,
!     and RESFS (= fS) according to P243.75 and P243.61.
!-----------------------------------------------------------------------
      DO J=1,TILE_PTS
        L = TILE_INDEX(J)
        I = LAND_INDEX(L)

!-----------------------------------------------------------------------
! Calculate the fraction of the flux with only aerodynamic resistance
! (canopy evaporation).
! Set to 1 for negative moisture flux (no surface/stomatal resistance to
! condensation).
!-----------------------------------------------------------------------
        FRACA(L) = 1.0
        IF ( DQ(L).LT.0.0 ) FRACA(L) = 0.0
        IF ( DQ(L).LT.0.0 .AND. CATCH(L).GT.0.0 )
     &    FRACA(L) = CANOPY(L) / ( EPDT(L) + CATCH(L) )
        FRACA(L) = MIN(FRACA(L),1.0)

!-----------------------------------------------------------------------
! Calculate resistance factors for transpiration from vegetation tiles
! and bare soil evaporation from soil tiles.
!-----------------------------------------------------------------------
        RESFS(L) = GC(L) / ( GC(L) + CH(L)*VSHR(I) )
        RESFT(L) = FRACA(L) + (1.0 - FRACA(L)) * RESFS(L)

      ENDDO

      IF (LTIMER) THEN
        CALL TIMER('SFRESIST',4)
      ENDIF

      RETURN
      END
