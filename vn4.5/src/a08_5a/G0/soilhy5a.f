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
!    SUBROUTINE SOIL_HYD-----------------------------------------------
!
! Subroutine Interface:
      SUBROUTINE SOIL_HYD (
     & NPNTS,NSHYD,SOIL_PTS,SOIL_INDEX,B,DZ,EXT,FW,KS,SATHH,TIMESTEP
     &,V_SAT,SLOW_RUNOFF,SMCL,STHU,W_FLUX,STF_SLOW_RUNOFF,LTIMER
     &)

      IMPLICIT NONE
!
! Description:
!     Increments the layer soil moisture contents and calculates
!     calculates gravitational runoff. Calls the following:
!
!     HYD_CON - to calculate the hydraulic conductivity
!                                                     (Cox, 6/95)
!
!     DARCY - to calculate the Darcian fluxes between soil layers
!                                                     (Cox, 6/95)
!
!
! Documentation : UM Documentation Paper 25
!
! Current Code Owner : David Gregory
!
! History:
! Version   Date     Comment
! -------   ----     -------
!  4.1      6/96     New deck.   Peter Cox
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
      REAL
     & POND_MAX        ! Maximum amount of ponded water (kg/m2).
     &,POND_FRAC_MAX   ! Maximum ponded fraction of the gridbox.
      PARAMETER (POND_MAX = 0.0, POND_FRAC_MAX = 0.5)
C-----------------------------------------------------------------------
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

! Subroutine arguments
!   Scalar arguments with intent(IN) :
      INTEGER
     & NPNTS                ! IN Number of gridpoints.
     &,NSHYD                ! IN Number of soil moisture levels.
     &,SOIL_PTS             ! IN Number of soil points.

      REAL
     & TIMESTEP             ! IN Model timestep (s).

      LOGICAL
     & STF_SLOW_RUNOFF      ! IN Stash flag for sub-surface runoff.



!   Array arguments with intent(IN) :
      INTEGER
     & SOIL_INDEX(NPNTS)    ! IN Array of soil points.

      REAL
     & B(NPNTS)             ! IN Clapp-Hornberger exponent.
     &,DZ(NSHYD)            ! IN Thicknesses of the soil layers (m).
     &,EXT(NPNTS,NSHYD)     ! IN Extraction of water from each soil
!                           !    layer (kg/m2/s).
     &,FW(NPNTS)            ! IN Throughfall from canopy plus snowmelt
!                           !    minus surface runoff (kg/m2/s).
     &,KS(NPNTS)            ! IN Saturated hydraulic conductivity
!                           !    (kg/m2/s).
     &,SATHH(NPNTS)         ! IN Saturated soil water pressure (m).
     &,V_SAT(NPNTS)         ! IN Volumetric soil moisture
!                           !    concentration at saturation
!                           !    (m3 H2O/m3 soil).
!
      LOGICAL LTIMER        ! Logical switch for TIMER diags

!   Array arguments with intent(OUT) :
      REAL
     & SLOW_RUNOFF(NPNTS)   ! OUT Drainage from the base of the
!                           !     soil profile (kg/m2/s).
     &,SMCLSAT(NPNTS,NSHYD) ! OUT The saturation moisture content of
!                           !     each layer (kg/m2).
     &,W_FLUX(NPNTS,0:NSHYD)! OUT The fluxes of water between layers
!                           !     (kg/m2/s).

!   Array arguments with intent(INOUT) :
      REAL
     & SMCL(NPNTS,NSHYD)    ! INOUT Total soil moisture contents
!                           !       of each layer (kg/m2).
     &,STHU(NPNTS,NSHYD)    ! INOUT Unfrozen soil moisture content of ea
!                           !       layer as a fraction of saturation.


! Local scalars:
      INTEGER
     & I,J,N                ! WORK Loop counters.

! Local arrays:
      REAL
     & DSMCL(NPNTS)         ! WORK The transfer of soil
!                           !      moisture (kg/m2/timestep).
     &,EXCESS(NPNTS)        ! WORK Excess soil moisture (kg/m2).
     &,SMCLMAX(NPNTS,NSHYD) ! WORK The maximum moisture content
!                           !      of each layer (kg/m2).
     &,SMCLU(NPNTS,NSHYD)   ! WORK Unfrozen soil moisture contents
!                           !      of each layer (kg/m2).
     &,STHUK(NPNTS)         ! WORK Fractional saturation of lowest
!                           !      layer.

! Function & Subroutine calls:
      EXTERNAL
     & HYD_CON,DARCY

      IF (LTIMER) THEN
        CALL TIMER('SOILHYD ',103) 
      ENDIF
!-----------------------------------------------------------------------
! Calculate the unfrozen soil moisture contents and the saturation
! total soil moisture for each layer.
!-----------------------------------------------------------------------
! CDIR$ IVDEP here would force vectorization but changes results!
! Fujitsu vectorization directive
!OCL NOVREC
      DO N=1,NSHYD
        DO J=1,SOIL_PTS
          I=SOIL_INDEX(J)
          SMCLSAT(I,N)=RHO_WATER*DZ(N)*V_SAT(I)
          SMCLU(I,N)=STHU(I,N)*SMCLSAT(I,N)
        ENDDO
      ENDDO

!-----------------------------------------------------------------------
! Top boundary condition
!-----------------------------------------------------------------------
      DO J=1,SOIL_PTS
        I=SOIL_INDEX(J)
        W_FLUX(I,0)=FW(I)
      ENDDO

!-----------------------------------------------------------------------
! Define the maximum water content that may exist in each layer
!-----------------------------------------------------------------------
      DO N=2,NSHYD
        DO J=1,SOIL_PTS
          I=SOIL_INDEX(J)
          SMCLMAX(I,N)=SMCLSAT(I,N)
        ENDDO
      ENDDO

!-----------------------------------------------------------------------
! Allow for some ponding of water by permitting excess moisture to
! remain in the top layer.
!-----------------------------------------------------------------------
      DO J=1,SOIL_PTS
        I=SOIL_INDEX(J)
        SMCLMAX(I,1)=SMCLSAT(I,1)+POND_MAX
      ENDDO

!-----------------------------------------------------------------------
! Increment the soil moisture contents for each layer, beginning with
! the bottom.
!-----------------------------------------------------------------------
      DO J=1,SOIL_PTS
        I=SOIL_INDEX(J)
        STHUK(I)=STHU(I,NSHYD)
      ENDDO
      CALL HYD_CON (NPNTS,SOIL_PTS,SOIL_INDEX,B,KS,STHUK,
     &              W_FLUX(1,NSHYD),LTIMER)

CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
      DO J=1,SOIL_PTS
        I=SOIL_INDEX(J)
        DSMCL(I)=(W_FLUX(I,NSHYD)+EXT(I,NSHYD))*TIMESTEP
        IF (DSMCL(I).GT.SMCLU(I,NSHYD)) THEN
          DSMCL(I)=SMCLU(I,NSHYD)
          W_FLUX(I,NSHYD)=DSMCL(I)/TIMESTEP-EXT(I,NSHYD)
        ENDIF

        SMCL(I,NSHYD)=SMCL(I,NSHYD)-TIMESTEP*(W_FLUX(I,NSHYD)
     &                              +EXT(I,NSHYD))
        SMCLU(I,NSHYD)=SMCLU(I,NSHYD)-TIMESTEP*(W_FLUX(I,NSHYD)
     &                              +EXT(I,NSHYD))
        STHU(I,NSHYD)=SMCLU(I,NSHYD)/SMCLSAT(I,NSHYD)
      ENDDO

!-----------------------------------------------------------------------
! Layers NSHYD to 1
!-----------------------------------------------------------------------
      DO N=NSHYD,2,-1
        CALL DARCY (NPNTS,SOIL_PTS,SOIL_INDEX,B,KS,SATHH
     &,             STHU(1,N-1),DZ(N-1),STHU(1,N),DZ(N)
     &,             W_FLUX(1,N-1),LTIMER)

CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
        DO J=1,SOIL_PTS
          I=SOIL_INDEX(J)
          DSMCL(I)=(W_FLUX(I,N-1)+EXT(I,N-1))*TIMESTEP
          IF (DSMCL(I).GT.SMCLU(I,N-1)) THEN
            DSMCL(I)=SMCLU(I,N-1)
            W_FLUX(I,N-1)=DSMCL(I)/TIMESTEP-EXT(I,N-1)
          ELSEIF (DSMCL(I).LT.(SMCL(I,N-1)-SMCLMAX(I,N-1))) THEN
            DSMCL(I)=SMCL(I,N-1)-SMCLMAX(I,N-1)
            W_FLUX(I,N-1)=DSMCL(I)/TIMESTEP-EXT(I,N-1)
          ENDIF

          DSMCL(I)=W_FLUX(I,N-1)*TIMESTEP
          IF (DSMCL(I).GT.(SMCLMAX(I,N)-SMCL(I,N))) THEN
            DSMCL(I)=SMCLMAX(I,N)-SMCL(I,N)
            W_FLUX(I,N-1)=DSMCL(I)/TIMESTEP
          ELSEIF (DSMCL(I).LT.(-SMCLU(I,N))) THEN
            DSMCL(I)=-SMCLU(I,N)
            W_FLUX(I,N-1)=DSMCL(I)/TIMESTEP
          ENDIF

          SMCL(I,N)=SMCL(I,N)+TIMESTEP*W_FLUX(I,N-1)
          SMCLU(I,N)=SMCLU(I,N)+TIMESTEP*W_FLUX(I,N-1)
          STHU(I,N)=SMCLU(I,N)/SMCLSAT(I,N)
          SMCL(I,N-1)=SMCL(I,N-1)-TIMESTEP*(W_FLUX(I,N-1)
     &                               +EXT(I,N-1))
          SMCLU(I,N-1)=SMCLU(I,N-1)-TIMESTEP*(W_FLUX(I,N-1)
     &                               +EXT(I,N-1))
          STHU(I,N-1)=SMCLU(I,N-1)/SMCLSAT(I,N-1)

        ENDDO
      ENDDO

CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
      DO J=1,SOIL_PTS
        I=SOIL_INDEX(J)
        SMCL(I,1)=SMCL(I,1)+TIMESTEP*W_FLUX(I,0)
        SMCLU(I,1)=SMCLU(I,1)+TIMESTEP*W_FLUX(I,0)
        STHU(I,1)=SMCLU(I,1)/SMCLSAT(I,1)
      ENDDO

!-----------------------------------------------------------------------
! If a layer is supersaturated move the excess moisture to the layer
! below.
!-----------------------------------------------------------------------
      DO N=1,NSHYD-1
CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
        DO J=1,SOIL_PTS
          I=SOIL_INDEX(J)

          EXCESS(I)=MAX((SMCL(I,N)-SMCLMAX(I,N)),0.0)

          IF (EXCESS(I).GT.0.0) THEN

            W_FLUX(I,N)=W_FLUX(I,N)+EXCESS(I)/TIMESTEP

            SMCL(I,N)=SMCL(I,N)-EXCESS(I)
            SMCLU(I,N)=SMCLU(I,N)-EXCESS(I)
            STHU(I,N)=SMCLU(I,N)/SMCLSAT(I,N)

            SMCL(I,N+1)=SMCL(I,N+1)+EXCESS(I)
            SMCLU(I,N+1)=SMCLU(I,N+1)+EXCESS(I)
            STHU(I,N+1)=SMCLU(I,N+1)/SMCLSAT(I,N+1)

          ENDIF

        ENDDO
      ENDDO

!-----------------------------------------------------------------------
! If there is still excess moisture add this to the flux from the base
! of the soil profile.
!-----------------------------------------------------------------------
CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
      DO J=1,SOIL_PTS
        I=SOIL_INDEX(J)
        EXCESS(I)=MAX((SMCL(I,NSHYD)-SMCLMAX(I,NSHYD)),0.0)
        IF (EXCESS(I).GT.0.0) THEN

          W_FLUX(I,NSHYD)=W_FLUX(I,NSHYD)+EXCESS(I)/TIMESTEP

          SMCL(I,NSHYD)=SMCL(I,NSHYD)-EXCESS(I)
          SMCLU(I,NSHYD)=SMCLU(I,NSHYD)-EXCESS(I)
          STHU(I,NSHYD)=SMCLU(I,NSHYD)/SMCLSAT(I,NSHYD)

        ENDIF
      ENDDO

!-----------------------------------------------------------------------
! Output slow runoff (drainage) diagnostic.
!-----------------------------------------------------------------------
      IF (STF_SLOW_RUNOFF) THEN
        DO I=1,NPNTS
          SLOW_RUNOFF(I)=0.0
        ENDDO
        DO J=1,SOIL_PTS
          I=SOIL_INDEX(J)
          SLOW_RUNOFF(I)=W_FLUX(I,NSHYD)
        ENDDO
      ENDIF
      IF (LTIMER) THEN
        CALL TIMER('SOILHYD ',104)
      ENDIF
      RETURN
      END
