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
!!!   SUBROUTINE SF_RESIST----------------------------------------------
!!!
!!!  Purpose: Calculate surface resistances for surface layer
!!!
!!! SDJ         <- programmer of some or all of previous code changes
!!!
C Modification History:
C Version Date     Change
C  4.5    Jul. 98  Kill the IBM specific lines (JCThil)
!!!--------------------------------------------------------------------

!  Arguaments -------------------------------------------------------
      SUBROUTINE SF_RESIST (
     & P_POINTS,LAND_PTS,P_FIELD,LAND_FIELD,LAND_MASK,INT_STOM,P1,LAND1,
     & LAND_INDEX,
     & ROOTD,SMVCCL,SMVCWT,SMC,V_SOIL,VFRAC,CANOPY,CATCH,DQ,
     & EPDT,LYING_SNOW,GC,RESIST,CHV,PSIS,FRACA,
     & RESFS,F_SE,RESFT,LTIMER
     & )
      IMPLICIT NONE


      INTEGER              !    Variables defining grid.
     & P_POINTS            ! IN Number of P-grid points to be processed
     &,LAND_PTS            ! IN Number of land points to be processed.
     &,P1                  ! IN First P-point to be processed.
     &,LAND1               ! IN First land point to be processed.
     &,P_FIELD             ! IN Total number of P-grid points.
     &,LAND_FIELD          ! IN Total number of land points.

     &,LAND_INDEX(LAND_FIELD)! IN Index for compressed land point array;
!                               i'th element holds position in the FULL
!                               field of the ith land pt to be
!                               processed
      LOGICAL
     & INT_STOM            ! IN T for interactive stomatal resistance.
     &,LTIMER

      REAL
     & CANOPY(LAND_FIELD)  ! IN Surface water (kg per sq metre).  F642.
     &,CATCH(LAND_FIELD)   ! IN Surface capacity (max. surface water)
!                               (kg per sq metre).  F6416.
     &,CHV(P_FIELD)        ! IN Transport coefficient for heat and
!                               moisture transport
     &,DQ(P_FIELD)         ! IN Sp humidity difference between surface
!                               and lowest atmospheric level (Q1 - Q*).
!                               Holds value over sea-ice where
!                               ICE_FRACT>0 i.e. Leads contribution not
!                               included.
     &,EPDT(P_FIELD)       ! IN "Potential" Evaporation * Timestep.
!                               Dummy variable for first call to routine
     &,GC(LAND_FIELD)      ! IN Interactive canopy conductance
     &,LYING_SNOW(P_FIELD) ! IN Lying snow amount (kg per sq metre).
     &,RESIST(LAND_FIELD)  ! IN "Stomatal" resistance to evaporation
!                               (seconds per metre).  F6415.
     &,ROOTD(LAND_FIELD)   ! IN "Root depth" (metres).  F6412.
     &,SMC(LAND_FIELD)     ! IN Soil moisture content (kg per sq m).
!                               F621.
     &,SMVCCL(LAND_FIELD)  ! IN Critical volumetric SMC (cubic metres
!                               per cubic metre of soil).  F6232.
     &,SMVCWT(LAND_FIELD)  ! IN Volumetric wilting point (cubic m of
!                               water per cubic m of soil).  F6231.
!                               Note: (SMVC!! - SMVCWT) is the critical
!                               volumetric available soil
     &,V_SOIL(LAND_FIELD)  ! IN Volumetric soil moisture concentration
!                               in the top soil layer (m3 H2O/m3 soil).
     &,VFRAC(LAND_FIELD)   ! IN Vegetated fraction.

      LOGICAL
     & LAND_MASK(P_FIELD)  ! IN .TRUE. for land; .FALSE. elsewhere. F60.

!  Output variables.

      REAL
     & FRACA(P_FIELD)      ! OUT Fraction of surface moisture flux with
!                                only aerodynamic resistance.
     &,PSIS(P_FIELD)       !     Soil moisture availability factor.
     &,RESFS(P_FIELD)      ! OUT Combined soil, stomatal and aerodynamic
!                                resistance factor = PSIS/(1+RS/RA) for
!                                fraction (1-FRACA)
     &,F_SE(P_FIELD)       ! OUT Fraction of the evapotranspiration
!                                which is bare soil evaporation.
     &,RESFT(P_FIELD)      ! OUT Total resistance factor
!                                FRACA+(1-FRACA)*RESFS.

!   Define local storage.

!   (a) Workspace.

C*L-----------COMDECK C_DENSTY FOR SUBROUTINE SF_EXCH----------
C RHOSEA    = density of sea water (kg/m3)
C RHO_WATER = density of pure water (kg/m3)
      REAL RHOSEA,RHO_WATER

      PARAMETER(RHOSEA = 1000.0,
     &          RHO_WATER = 1000.0)
C*----------------------------------------------------------------------
      REAL
     + DZSOIL(4)               ! Soil layer thicknesses (m).
      DATA DZSOIL /0.100, 0.250, 0.650, 2.000 /
C-----------------------------------------------------------------------
! History:
! Version  Date  Comment
!  3.4   18/5/94 Add PP missing data indicator. J F Thomson
C*L------------------COMDECK C_MDI-------------------------------------
      REAL    RMDI_PP   ! PP missing data indicator (-1.0E+30)
      REAL    RMDI_OLD  ! Old real missing data indicator (-32768.0)
      REAL    RMDI      ! New real missing data indicator (-2**30)
      INTEGER IMDI      ! Integer missing data indicator

      PARAMETER(
     & RMDI_PP  =-1.0E+30,
     & RMDI_OLD =-32768.0,
     & RMDI =-32768.0*32768.0,
     & IMDI =-32768)
C*----------------------------------------------------------------------


!  Workspace --------------------------------------------------------
      INTEGER
     & I           ! Loop counter (horizontal field index).
     &,L           ! Loop counter (land field index).

      REAL
     & FSMC        ! Soil moisture factor for bare soil evaporation.
     &,SMCRIT      ! "Critical" available SMC in kg per sq m.


      EXTERNAL TIMER

      IF (LTIMER) THEN
        CALL TIMER('SFRESIST',3)
      ENDIF

!-----------------------------------------------------------------------
!!  1 Evaporation over land surfaces without snow is limited by
!!    soil moisture availability and stomatal resistance.
!!    Set FRACA (= fA in the documentation) according to P243.68,
!!    PSIS according to P243.65, and RESFS (= fS) according to P243.75
!!    and P243.61, using neutral-stability value of CH (as explained
!!    in section (v) of the P243 documentation).
!-----------------------------------------------------------------------

      DO I=P1,P1+P_POINTS-1

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!  1.1 Set parameters (workspace) to values relevant for non-land
!!      points.  Only aerodynamic resistance applies.
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        FRACA(I) = 1.0
        PSIS(I)  = 0.0
        RESFT(I) = 1.0
        RESFS(I) = 0.0
      ENDDO

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!  1.2 Over-write workspace for other points.
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
      DO L=LAND1,LAND1+LAND_PTS-1
        I = LAND_INDEX(L)

!-----------------------------------------------------------------------
! If the interactive stomatal resistance is being used, calculate the
! soil water factor for bare soil evaporation
!-----------------------------------------------------------------------

        IF (INT_STOM) THEN

          IF (V_SOIL(L) .GT. SMVCCL(L)) THEN
            FSMC = 1.0
          ELSEIF (V_SOIL(L) .LE. SMVCWT(L)) THEN
            FSMC = 0.0
          ELSE
            FSMC = (V_SOIL(L) - SMVCWT(L))
     &           / (SMVCCL(L) - SMVCWT(L))
          ENDIF

        ELSE
!  Calculate the soil moisture availability factor, PSIS.

        SMCRIT = RHO_WATER * ROOTD(L) * (SMVCCL(L)-SMVCWT(L))
!                                                            ... P243.66

        PSIS(I) = 0.0
        IF (SMC(L).GE.SMCRIT .AND. SMCRIT.GT.0.0)
     &    PSIS(I) = 1.0
        IF (SMC(L).LT.SMCRIT .AND. SMC(L).GT.0.0)
     &    PSIS(I) = SMC(L)/SMCRIT


        ENDIF ! end of INT_STOM block

!  For snow-covered land or land points with negative moisture flux
!  set the fraction of the flux with only aerodynamic resistance to 1
!  (no surface/stomatal resistance to evap from snow or condensation).

        FRACA(I) = 1.0

!  When there is positive moisture flux from snow-free land, calculate
!  the fraction of the flux from the "canopy".

        IF (DQ(I).LT.0.0 .AND. LYING_SNOW(I).LE.0.0) FRACA(I) = 0.0
        IF (DQ(I).LT.0.0.AND.LYING_SNOW(I).LE.0.0.AND.CATCH(L).GT.0.0)
     &    FRACA(I) = CANOPY(L)/(EPDT(I) + CATCH(L))


!-----------------------------------------------------------------------
! If the interactive stomatal resistance is being used calculate
! separate resistance factors for bare soil evaporation and
! transpiration. Assume a surface resistance of 100 s/m for bare soil.
!-----------------------------------------------------------------------

        IF (INT_STOM) THEN       ! Interactive Canopy Resistance

!-----------------------------------------------------------------------
! Set resistance and moisture availability factors to zero for land ice
!-----------------------------------------------------------------------
          IF (GC(L).EQ.RMDI) THEN  ! land-ice points

            PSIS(I) = 0.0
            RESFS(I) = 0.0
            F_SE(I) = 0.0

          ELSE

!-----------------------------------------------------------------------
! If the interactive stomatal resistance is being used set the moisture
! availability factor to one, since moisture stress is already taken
! account of in SF_STOM  (Peter Cox 21/11/95).
!-----------------------------------------------------------------------

            PSIS(I) = 1.0
            RESFS(I) = VFRAC(L) * GC(L) / ( GC(L) + CHV(I))
     &       + (1 - VFRAC(L)) * FSMC / (1.0 + CHV(I)*100.0)

            F_SE(I) = 0.0

            IF (RESFS(I) .GT. 0.0) THEN
              F_SE(I) =  (1 - VFRAC(L)) * FSMC
     &                 / (RESFS(I)*(1.0 + CHV(I)*100.0))
            ENDIF

          ENDIF

        ELSE

          RESFS(I) = PSIS(I) / ( 1.0 + CHV(I)*RESIST(L) )
          F_SE(I)=0

        ENDIF

        RESFT(I) = FRACA(I) + (1.0 - FRACA(I)) * RESFS(I)


      ENDDO         ! Evaporation over land points only, section 3.4.2

      IF (LTIMER) THEN
        CALL TIMER('SFRESIST',4)
      ENDIF

      RETURN
      END
