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
!    SUBROUTINE SMC_ROOT----------------------------------------------
!
! Subroutine Interface:
      SUBROUTINE SMC_ROOT (NPNTS,NSHYD,F_TYPE,DZ,ROOTD,STHU,VFRAC,    
     &               V_SAT,V_WILT,SMC,V_ROOT,V_SOIL,WT_EXT,LTIMER)      

      IMPLICIT NONE
!
! Description:
!     Calculates the volumetric soil moisture in the top soil layer,
!     the volumetric soil moisture in the rootzone, the gridbox mean
!     available soil moisture, and the fraction of the transpiration
!     which is extracted from each soil layer.             (Cox, 2/96)
!
!
! Documentation : UM Documentation Paper 25
!
! Current Code Owner : David Gregory
!
! History:
! Version   Date     Comment
! -------   ----     -------
!  4.1               New deck   Peter Cox
!LL   4.5   18/06/98  Changed Timer calls to indicate non-barrier
!LL                                                   P.Burton
!  4.5      6/98     Optional exponential root profile Peter Cox
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!
! System component covered: P25
! System Task: P25
!

! Global variables:
C*L-----------COMDECK C_DENSTY FOR SUBROUTINE SF_EXCH----------
C RHOSEA    = density of sea water (kg/m3)
C RHO_WATER = density of pure water (kg/m3)
      REAL RHOSEA,RHO_WATER

      PARAMETER(RHOSEA = 1000.0,
     &          RHO_WATER = 1000.0)
C*----------------------------------------------------------------------
C RHO_WATER removed to avoid clash with declaration in C_DENSTY
C J.Smith 28/06/95
      REAL OMEGA1,RHO_SNOW,DEFF_SNOW,SNOW_HCON,SNOW_HCAP
      INTEGER PSOIL
      PARAMETER (
     + PSOIL=4                  ! No. of soil layers (must = NSOIL).
     +,OMEGA1=3.55088E-4        ! Tunable characteristic freq (rad/s).
     +,RHO_SNOW=250.0           ! Density of lying snow (kg per m**3).
     +,DEFF_SNOW=0.1            ! Depth of `effective' snow surface
C                               ! layer (m).
     +,SNOW_HCON=0.265          ! Thermal conductivity of lying snow
C                               ! (Watts per m per K).
     +,SNOW_HCAP=0.63E6         ! Thermal capacity of lying snow
C                               ! (J/K/m3)
     +)
      INTEGER
     & CAN_MODEL            ! 1 for no thermal canopy (pre 4.5 UM).
C                           ! 2 for radiative coupling between
C                           !   vegetated surface and first soil
C                           !   temperature.
C                           ! 3 for radiative coupling between
C                           !   vegetated surface and first soil
C                           !   temperature, plus canopy thermal
C                           !   capacity.
     &,REX_MODEL            ! 1 for uniform root density profile
C                           !   (pre 4.5 UM) with MOSES I
C                           !   rootdepths.
C                           ! 2 for exponential root density
C                           !   profile with "old" rootdepths.
     &,TF_MODEL             ! 1 for uniformily distributed canopy
C                           !   water (pre 4.5 UM).
C                           ! 2 for bimodally distributed canopy
C                           !   water with random overlap.
C                           ! 3 for bimodally distributed canopy
C                           !   water with maximum overlap.
C
C For pre 4.5 UM MOSES I choose:
C     PARAMETER (CAN_MODEL=1, REX_MODEL=1, TF_MODEL=1)
C
      PARAMETER (CAN_MODEL=1, REX_MODEL=1, TF_MODEL=1)

! Subroutine arguments:
!   Scalar arguments with intent(IN) :
      INTEGER
     & NPNTS                ! IN Number of gridpoints.
     &,NSHYD                ! IN Number of soil moisture levels.
     &,F_TYPE(NPNTS)        ! IN Plant functional type:
!                           !     1 - Broadleaf tree
!                           !     2 - Needleleaf tree
!                           !     3 - C3 Grass
!                           !     4 - C4 Grass

!   Array arguments with intent(IN) :
      REAL
     & DZ(NSHYD)            ! IN Soil layer thicknesses (m).
     &,ROOTD(NSHYD)         ! IN Rootdepth (m).
     &,STHU(NPNTS,NSHYD)    ! IN Unfrozen soil moisture content of
!                           !    each layer as a fraction of
!                           !    saturation.
     &,VFRAC(NPNTS)         ! IN Vegetated fraction.
     &,V_SAT(NPNTS)         ! IN Volumetric soil moisture
!                           !    concentration at saturation
!                           !    (m3 H2O/m3 soil).
     &,V_WILT(NPNTS)        ! IN Volumetric soil moisture
!                           !    concentration below which
!                           !    stomata close (m3 H2O/m3 soil).

      LOGICAL LTIMER        ! Logical switch for TIMER diags

!   Array arguments with intent(OUT) :
      REAL
     & SMC(NPNTS)           ! OUT Available soil moisture in the
!                           !     rootzone (kg/m2).
     &,V_ROOT(NPNTS)        ! OUT Volumetric soil moisture
!                           !     concentration in the rootzone
!                           !     (m3 H2O/m3 soil).
     &,V_SOIL(NPNTS)        ! OUT Volumetric soil moisture
!                           !     concentration in the top soil
!                           !     layer (m3 H2O/m3 soil).
     &,WT_EXT(NPNTS,NSHYD)  ! OUT Fraction of transpiration extracted
!                           !     from each soil layer (kg/m2/s).
! Local scalars:
      INTEGER
     & I,J,N                ! WORK Loop counters

! Local arrays:
      REAL
     & RHO_ROOT(NPNTS,NSHYD)! WORK Density of roots in each soil layer
!                           !      (normalised).
     &,RHO_RNORM(NPNTS)     ! WORK Normalisation factor for RHO_ROOT.
     &,SMCLA(NPNTS,NSHYD)   ! WORK Available soil moisture in each
!                           !      layer (scaled with root density)
!                           !      (kg/m2)
     &,Z(0:NSHYD)           ! WORK Depths of soil layer boundaries (m).
     &,Z_ROOT(NPNTS)        ! WORK Rootdepth of the vegetated area (m).

!----------------------------------------------------------------------
! Functional type dependent parameters
!----------------------------------------------------------------------
      INTEGER
     & R_LAYERS(4)          ! Number of soil layers from which water
!                           ! can be extracted.
!----------------------------------------------------------------------
!                     BT   NT  C3G  C4G
!----------------------------------------------------------------------
      DATA R_LAYERS/   4,   4,   3,   3/


      IF (LTIMER) THEN
        CALL TIMER('SMROOT  ',103)
      ENDIF

!----------------------------------------------------------------------
! Initialisations
!----------------------------------------------------------------------
      DO I=1,NPNTS
        SMC(I)=0.0
        V_ROOT(I)=0.0
        Z_ROOT(I)=0.0
        RHO_RNORM(I)=1.0
        V_SOIL(I)=STHU(I,1)*V_SAT(I)
      ENDDO


!----------------------------------------------------------------------
! Calculate the root density in each layer, assuming either:     
!   An exponential profile
!----------------------------------------------------------------------
      IF (REX_MODEL .EQ. 2) THEN

        Z(0)=0.0
        DO N=1,NSHYD
          Z(N)=Z(N-1)+DZ(N)
        ENDDO

        DO I=1,NPNTS
!----------------------------------------------------------------------
! Assume here that the gridbox mean rootdepth includes a contribution
! of 0.1m from the non-vegetated area
!----------------------------------------------------------------------
          IF (VFRAC(I).GT.0.0) THEN
            Z_ROOT(I)=(ROOTD(I)-0.1*(1.0-VFRAC(I)))/VFRAC(I)
          ELSE
            Z_ROOT(I)=0.0
          ENDIF

          IF (Z_ROOT(I).GT.0.0) THEN
            RHO_RNORM(I)=(1-EXP(-Z(NSHYD)/Z_ROOT(I)))/Z_ROOT(I)
          ENDIF
        ENDDO

        DO N=1,NSHYD
          DO I=1,NPNTS                                                  
            IF (V_SAT(I).GT.0.0 .AND. Z_ROOT(I).GT.0.0) THEN
              RHO_ROOT(I,N)=(EXP(-Z(N-1)/Z_ROOT(I))
     &           -EXP(-Z(N)/Z_ROOT(I)))/(DZ(N)*RHO_RNORM(I))
            ENDIF                                                       
          ENDDO                                                         
        ENDDO                                                           

!----------------------------------------------------------------------
!   A uniform profile
!----------------------------------------------------------------------
      ELSE

        DO N=1,NSHYD                                                    
          DO I=1,NPNTS                                                  
            IF (R_LAYERS(F_TYPE(I)).GE.N) THEN    
              RHO_ROOT(I,N)=1.0                    
              Z_ROOT(I)=Z_ROOT(I)+DZ(N)          
            ELSE                                
              RHO_ROOT(I,N)=0.0                
            ENDIF                             
          ENDDO                              
        ENDDO

      ENDIF

!----------------------------------------------------------------------
! Calculate the volumetric soil moisture in the rootzone and the
! available moisture in each soil layer. Only do calculation for
! non land-ice points (V_SAT is set to zero for land-ice points).
!----------------------------------------------------------------------
      DO N=1,NSHYD
        DO I=1,NPNTS

          IF (V_SAT(I).GT.0.0 .AND. Z_ROOT(I).GT.0.0) THEN
            V_ROOT(I)=V_ROOT(I)+RHO_ROOT(I,N)*STHU(I,N)
     &                         *V_SAT(I)*(DZ(N)/Z_ROOT(I))
            SMCLA(I,N)=RHO_ROOT(I,N)*(STHU(I,N)*V_SAT(I)-V_WILT(I))
     &                              *RHO_WATER*DZ(N)
            SMCLA(I,N)=MAX(SMCLA(I,N),0.0)
            SMC(I)=SMC(I)+SMCLA(I,N)
          ENDIF

        ENDDO
      ENDDO

!----------------------------------------------------------------------
! Calculate the fraction of the tranpiration which is extracted from
! each soil layer.
!----------------------------------------------------------------------
      DO N=1,NSHYD
        DO I=1,NPNTS

          IF (V_SAT(I).GT.0.0) THEN

            IF (SMC(I).GT.0.0) THEN
              WT_EXT(I,N)=SMCLA(I,N)/SMC(I)
            ELSE
              WT_EXT(I,N)=0.0
            ENDIF
          ELSE
              WT_EXT(I,N)=0.0
          ENDIF

        ENDDO
      ENDDO

!---------------------------------------------------------------------
! Diagnose the gridbox mean available soil moisture
!---------------------------------------------------------------------
      DO I=1,NPNTS
        SMC(I)=VFRAC(I)*SMC(I)
     &        +(1-VFRAC(I))*RHO_WATER*DZ(1)*
     &          MAX(0.0,(V_SOIL(I)-V_WILT(I)))
      ENDDO

      IF (LTIMER) THEN
        CALL TIMER('SMROOT  ',104)
      ENDIF

      RETURN
      END
