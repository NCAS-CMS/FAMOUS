C *****************************COPYRIGHT******************************
C (c) CROWN COPYRIGHT 1996, METEOROLOGICAL OFFICE, All Rights Reserved.
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
CLL  Purpose: Calculate bulk Richardson number for surface layer
CLL
C Modification History:
C Version Date     Change
C  4.5    Jul. 98  Kill the IBM specific lines (JCThil)
CLL
CLLEND-----------------------------------------------------------------
C*
C*L  Arguaments -------------------------------------------------------
      SUBROUTINE SF_RESIST (
     & P_POINTS,LAND_PTS,LAND_MASK,INT_STOM,
     & P1,LAND_INDEX,
     & ROOTD,SMVCCL,SMVCWT,SMC,V_SOIL,VFRAC,CANOPY,CATCH,DQ,
     & EPDT,LYING_SNOW,GC,RESIST,VSHR,CH,PSIS,FRACA,
     & RESFS,F_SE,RESFT,LTIMER
     & )
      IMPLICIT NONE


      INTEGER              !    Variables defining grid.
     & P_POINTS            ! IN Number of P-grid points to be processed
     &,LAND_PTS            ! IN Number of land points to be processed.

     &,LAND_INDEX(LAND_PTS)! IN Index for compressed land point array;
C                          !    i'th element holds position in the FULL
C                          !    field of the ith land pt to be
C                          !    processed
     &,P1                  ! IN First P-point to be processed.

      LOGICAL
     & INT_STOM            ! IN T for interactive stomatal resistance.

      REAL
     & CANOPY(LAND_PTS)    ! IN Surface water (kg per sq metre).  F642.
     &,CATCH(LAND_PTS)     ! IN Surface capacity (max. surface water)
C                          !    (kg per sq metre).  F6416.
     &,CH(P_POINTS)        ! IN Transport coefficient for heat and
C                          !    moisture transport
     &,DQ(P_POINTS)        ! IN Sp humidity difference between surface
C                          !    and lowest atmospheric level (Q1 - Q*).
C                          !    Holds value over sea-ice where
C                          !    ICE_FRACT>0 i.e. Leads contribution not
C                          !    included.
     &,EPDT(P_POINTS)      ! IN "Potential" Evaporation * Timestep.
C                          !    Dummy variable for first call to routine
     &,LYING_SNOW(P_POINTS)! IN Lying snow amount (kg per sq metre).
     &,GC(LAND_PTS)        ! IN Interactive canopy conductance
C                          !    to evaporation (m/s)
     &,RESIST(LAND_PTS)    ! IN Fixed "stomatal" resistance
C                          !    to evaporation (s/m)
     &,ROOTD(LAND_PTS)     ! IN "Root depth" (metres).  F6412.
     &,SMC(LAND_PTS)       ! IN Soil moisture content (kg per sq m).
C                          !    F621.
     &,SMVCCL(LAND_PTS)    ! IN Critical volumetric SMC (cubic metres
C                          !    per cubic metre of soil).  F6232.
     &,SMVCWT(LAND_PTS)    ! IN Volumetric wilting point (cubic m of
C                          !    water per cubic m of soil).  F6231.
     &,V_SOIL(LAND_PTS)    ! IN Volumetric soil moisture concentration
C                          !    in the top soil layer (m3 H2O/m3 soil).
     &,VFRAC(LAND_PTS)     ! IN Vegetated fraction.
     &,VSHR(P_POINTS)      ! IN Magnitude of surface-to-lowest-lev. wind
C
C    Note: (SMVCCL - SMVCWT) is the critical volumetric available soil
      LOGICAL
     & LAND_MASK(P_POINTS) ! IN .TRUE. for land; .FALSE. elsewhere. F60.
C
C  Output variables.
C
      REAL
     & FRACA(P_POINTS)    ! OUT Fraction of surface moisture flux with
C                         !     only aerodynamic resistance.
     &,PSIS(P_POINTS)     ! OUT Soil moisture availability factor.
     &,RESFS(P_POINTS)    ! OUT Combined soil, stomatal and aerodynamic
C                         !     resistance factor = PSIS/(1+RS/RA) for
C                         !     fraction (1-FRACA).
     &,F_SE(P_POINTS)     ! OUT Fraction of the evapotranspiration
C                         !     which is bare soil evaporation.
     &,RESFT(P_POINTS)    ! OUT Total resistance factor
C                         !     FRACA+(1-FRACA)*RESFS.

      LOGICAL LTIMER      ! Logical switch for TIMER diags
C   Define local storage.
C
C   (a) Workspace.
C
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


C*L  Workspace --------------------------------------------------------
      INTEGER
     & I           ! Loop counter (horizontal field index).
     &,J           ! Loop counter (land field index).
     &,L           ! Loop counter (land field index).
      REAL
     & FSMC        ! Soil moisture factor for bare soil evaporation.
     &,SMCRIT      ! "Critical" available SMC in kg per sq m.

      IF (LTIMER) THEN
        CALL TIMER('SFRESIST',3)
      ENDIF
C-----------------------------------------------------------------------
CL  6 Evaporation over land surfaces without snow is limited by
CL    soil moisture availability and stomatal resistance.
C     Set FRACA (= fA in the documentation) according to P243.68,
C     PSIS according to P243.65, and RESFS (= fS) according to P243.75
C     and P243.61, using neutral-stability value of CH (as explained
C     in section (v) of the P243 documentation).
C-----------------------------------------------------------------------
      DO I=1,P_POINTS
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
CL  6.1 Set parameters (workspace) to values relevant for non-land
CL      points.  Only aerodynamic resistance applies.
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        FRACA(I) = 1.0
        PSIS(I) = 0.0
        RESFT(I) = 1.0
        RESFS(I) = 0.0
      ENDDO
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
CL  6.2 Over-write workspace for other points.
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
      DO L=1,LAND_PTS
        I = LAND_INDEX(L) - (P1-1)

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
C
C  Calculate the soil moisture availability factor, PSIS.
C
          SMCRIT = RHO_WATER * ROOTD(L) * (SMVCCL(L)-SMVCWT(L))
C                                                            ... P243.66
C
          PSIS(I) = 0.0
          IF (SMC(L).GE.SMCRIT .AND. SMCRIT.GT.0.0)
     &      PSIS(I) = 1.0
          IF (SMC(L).LT.SMCRIT .AND. SMC(L).GT.0.0)
     &      PSIS(I) = SMC(L)/SMCRIT

        ENDIF
C
C  For snow-covered land or land points with negative moisture flux
C  set the fraction of the flux with only aerodynamic resistance to 1
C  (no surface/stomatal resistance to evap from snow or condensation).
C
        FRACA(I) = 1.0
C
C  When there is positive moisture flux from snow-free land, calculate
C  the fraction of the flux from the "canopy".
C
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
            RESFS(I) = VFRAC(L) * GC(L) / ( GC(L) + CH(I)*VSHR(I))
     &       + (1 - VFRAC(L)) * FSMC / (1.0 + CH(I)*VSHR(I)*100.0)

            F_SE(I) = 0.0

            IF (RESFS(I) .GT. 0.0) THEN
              F_SE(I) =  (1 - VFRAC(L)) * FSMC
     &                 / (RESFS(I)*(1.0 + CH(I)*VSHR(I)*100.0))
            ENDIF

          ENDIF

        ELSE

          RESFS(I) = PSIS(I) / ( 1.0 + CH(I)*VSHR(I)*RESIST(L))
          F_SE(I) = 0.0

        ENDIF

        RESFT(I) = FRACA(I) + (1.0 - FRACA(I)) * RESFS(I)
      ENDDO         ! Evaporation over land points only, section 3.4.2

      IF (LTIMER) THEN
        CALL TIMER('SFRESIST',4)
      ENDIF
      RETURN
      END
