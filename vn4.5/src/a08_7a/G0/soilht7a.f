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
!    SUBROUTINE SOIL_HTC------------------------------------------------
!
! Subroutine Interface:
      SUBROUTINE SOIL_HTC (
     & NPNTS,NSHYD,SOIL_PTS,SOIL_INDEX,B,DZ,HCAP,HCON,SATHH
     &,SNOW_SOIL_HTF,SOIL_SURF_HTF,TIMESTEP,V_SAT,W_FLUX
     &,SMCL,STHU,STHF,TSOIL
     +,LTIMER
     +)

      IMPLICIT NONE
!
! Description:
!     Updates deep soil temperatures, frozen and unfrozen
!     frozen soil water content. Calls the following subroutines:
!
!     HEAT_CON - to calculate the soil thermal conductivity
!                                                     (Cox, 6/95)
!
! Documentation : UM Documentation Paper 25
!
! Current Code Owner : David Gregory
!
! History:
! Version   Date     Comment
! -------   ----     -------
!  4.1               New deck.   Peter Cox
!  4.4      7/97     MOSES II.   Richard Essery
!  4.5   26/5/98     Correction to stop overflow when SMCL very small.
!                    C.Bunton    
!LL   4.5   18/06/98  Changed Timer calls to indicate non-barrier
!LL                                                   P.Burton
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!
! System component covered: P25
! System Task: P25
!

! Global variables:
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


!   Array arguments with intent(IN) :
      INTEGER
     & SOIL_INDEX(NPNTS)    ! IN Array of soil points.

      REAL
     & B(NPNTS)             ! IN Clapp-Hornberger exponent.
     &,DZ(NSHYD)            ! IN Thicknesses of the soil layers (m).
     &,HCAP(NPNTS)          ! IN Soil heat capacity (J/K/m3).
     &,HCON(NPNTS)          ! IN Soil thermal conductivity (W/m/K).
     &,SATHH(NPNTS)         ! IN Saturated soil water pressure (m).
     &,SMCL(NPNTS,NSHYD)    ! IN Soil moisture content of each
!                           !    layer (kg/m2).
     &,SNOW_SOIL_HTF(NPNTS) ! IN Heat flux from snow to soil (W/m2).
     &,SOIL_SURF_HTF(NPNTS) ! IN Net downward surface heat flux (W/m2).
     &,V_SAT(NPNTS)         ! IN Volumetric soil moisture
!                           !    concentration at saturation
!                           !    (m3 H2O/m3 soil).
     +,W_FLUX(NPNTS,0:NSHYD)! IN The fluxes of water between layers
!                           !    (kg/m2/s).
C
      LOGICAL LTIMER        ! Logical switch for TIMER diags

!   Array arguments with intent(OUT) :

!   Array arguments with intent(INOUT) :
      REAL
     & STHF(NPNTS,NSHYD)    ! INOUT Frozen soil moisture content of
!                           !       each layer as a fraction of
!                           !       saturation.
     &,STHU(NPNTS,NSHYD)    ! INOUT Unfrozen soil moisture content of
!                           !       each layer as a fraction of
!                           !       saturation.
     &,TSOIL(NPNTS,NSHYD)   ! INOUT Sub-surface temperatures (K).

! Local scalars:
      INTEGER
     & I,J,M,N              ! WORK Loop counters.
     &,J_ITER               ! WORK Number of soil points which require
!                           !      iteration.
     &,MMAX                 ! WORK Maximum number of iterations on
!                           !      temperature.
     &,FIRST_SOIL_PT        !      First soil point

      REAL
     & FACUR                ! WORK Required flux conservation accuracy
!                           !      (W/m2).
     &,TACUR                ! WORK Required accuracy of temperature
!                           !      calculation (Celsius).
      PARAMETER( MMAX=3, FACUR=0.01, TACUR=0.00000 )

! Local arrays:
      INTEGER
     & ITER_PTS(NPNTS)      ! WORK Array of soil points which require
!                           !      iteration.

      REAL
     & CEACUR(NPNTS)        ! WORK Flux conservation accuracy of the
!                           !      calculation (W/m2)
     &,DHSL0(NPNTS,NSHYD)   ! WORK Total heat increment to the layer
!                           !      (J/m2/timestep)
     &,DHSL(NPNTS,NSHYD)    ! WORK The heat available to update the
!                           !      layer temperature (J/m2/timestep)
     &,DHSLA(NPNTS,NSHYD)   ! WORK The heat increment due to advection o
!                           !      heat by moisture fluxes (J/m2/timeste
     &,DSMCLF(NPNTS,NSHYD)  ! WORK The increment to the layer frozen
!                           !      soil moisture (kg/m2/timestep).
     &,DTHU(NPNTS,NSHYD)    ! WORK Rate of change of volumetric unfrozen
!                           !      soil moisture concentration with
!                           !      temperature (m3 liquid H2O/m3 soil/K)
     &,DTSL(NPNTS,NSHYD)    ! WORK The increment to the layer temperatur
!                           !      (K/timestep).
     &,HCAPT(NPNTS)         ! WORK The total volumetric heat capacity
!                           !      (soil+water) of the layer (J/m3/K).
     &,HC(NPNTS,NSHYD)      ! WORK The thermal conductivity of each
!                           !      layer (W/m/K).
     &,HCONS(NPNTS)         ! WORK The thermal conductivity between
!                           !      adjacent soil layers (W/m/K).
     &,H_FLUX(NPNTS,0:NSHYD) !WORK The fluxes of heat between layers
!                           !      (W/m2).
     &,SMCLF(NPNTS,NSHYD)   ! WORK Frozen moisture content of each
!                           !      soil layer (kg/m2).
     &,SMCLF0(NPNTS,NSHYD)  ! WORK Previous value of SMCLF (kg/m2).
     &,SMCLSAT(NPNTS,NSHYD) ! WORK The saturation moisture content of
!                           !      each layer (kg/m2).
     &,SMCLU(NPNTS,NSHYD)   ! WORK Unfrozen moisture content of each
!                           !      soil layer (kg/m2).
     &,SMCLU0(NPNTS,NSHYD)  ! WORK Previous value of SMCLU (kg/m2).
     &,SMCFU(NPNTS)         ! WORK Fractional saturation (unfrozen water
!                           !      at layer boundaries.
     &,SMCFF(NPNTS)         ! WORK Fractional saturation (frozen water)
!                           !      at layer boundaries.
     &,TMAX(NPNTS,NSHYD)    ! WORK Temperature above which all water is
!                           !      unfrozen (Celsius)
     &,TSL(NPNTS,0:NSHYD)   ! WORK Soil layer temperatures (Celsius)
!                           !      TSL(0) temperature of incoming water.
!                           !      TSL(1:NSHYD) sub-surface soil
!                           !      temperatures .
     &,TSL0(NPNTS,0:NSHYD)  ! WORK Previous value of TSL (Celsius).

! Function & Subroutine calls:
      EXTERNAL
     & HEAT_CON

C*L-----------COMDECK C_DENSTY FOR SUBROUTINE SF_EXCH----------
C RHOSEA    = density of sea water (kg/m3)
C RHO_WATER = density of pure water (kg/m3)
      REAL RHOSEA,RHO_WATER

      PARAMETER(RHOSEA = 1000.0,
     &          RHO_WATER = 1000.0)
C*----------------------------------------------------------------------
C*L------------------COMDECK C_LHEAT------------------------------------
C LC IS LATENT HEAT OF CONDENSATION OF WATER AT 0DEGC
C LF IS LATENT HEAT OF FUSION AT 0DEGC
      REAL LC,LF

      PARAMETER(LC=2.501E6,
     &          LF=0.334E6)
C*----------------------------------------------------------------------
      REAL HCAPW,HCAPI,RHO_ICE,DPSIDT
      PARAMETER(
     + HCAPW=4180      ! Specific heat capacity of water (J/kg/K)
     +,HCAPI=2100      ! Specific heat capacity of ice (J/kg/K)
     +,RHO_ICE=917     ! Density of ice (kg/m3)
     +,DPSIDT=114.3    ! Rate of change of ice potential with temperatur
C                      ! RHO_ICE*LF/ZERODEGC*1/(RHO_WATER*G) (m/K)
     +)
C-----------------------------------------------------------------------
C*L------------------COMDECK C_O_DG_C-----------------------------------
C ZERODEGC IS CONVERSION BETWEEN DEGREES CELSIUS AND KELVIN
C TFS IS TEMPERATURE AT WHICH SEA WATER FREEZES
C TM IS TEMPERATURE AT WHICH FRESH WATER FREEZES AND ICE MELTS
      REAL ZERODEGC,TFS,TM

      PARAMETER(ZERODEGC=273.15,
     &          TFS=271.35,
     &          TM=273.15)
C*----------------------------------------------------------------------


      IF (LTIMER) THEN
        CALL TIMER('SOILHTC ',103)
      ENDIF

!-----------------------------------------------------------------------
! Initialisations
!-----------------------------------------------------------------------
      FIRST_SOIL_PT=SOIL_INDEX(1)

      DO J=1,SOIL_PTS
        I=SOIL_INDEX(J)
        TSL(I,0)=TSOIL(I,1)-ZERODEGC
      ENDDO

      DO N=1,NSHYD
        DO J=1,SOIL_PTS
          I=SOIL_INDEX(J)
!-----------------------------------------------------------------------
! Define soil layer temperatures TSL (in celsius).
!-----------------------------------------------------------------------
          TSL(I,N)=TSOIL(I,N)-ZERODEGC
        ENDDO
      ENDDO

      DO N=1,NSHYD
! CDIR$ IVDEP here would force vectorization but changes results!
        DO J=1,SOIL_PTS
          I=SOIL_INDEX(J)
!-----------------------------------------------------------------------
! Diagnose the frozen and unfrozen water.
!-----------------------------------------------------------------------

          SMCLSAT(I,N)=RHO_WATER*DZ(N)*V_SAT(I)
          SMCLF(I,N)=SMCLSAT(I,N)*STHF(I,N)
          SMCLU(I,N)=SMCL(I,N)-SMCLF(I,N)
        ENDDO                                  !  J=1,SOIL_PTS
      ENDDO                                    ! N=1,NSHYD

!----------------------------------------------------------------------
! Calculate the heat conductivity between adjacent layers
!----------------------------------------------------------------------

      DO N=1,NSHYD-1
        DO J=1,SOIL_PTS
          I=SOIL_INDEX(J)
          SMCFU(I)=(DZ(N+1)*STHU(I,N)+DZ(N)*STHU(I,N+1))
     &              /(DZ(N+1)+DZ(N))
          SMCFF(I)=(DZ(N+1)*STHF(I,N)+DZ(N)*STHF(I,N+1))
     &              /(DZ(N+1)+DZ(N))
        ENDDO
        CALL HEAT_CON(NPNTS-FIRST_SOIL_PT+1,HCON(FIRST_SOIL_PT)
     &,              SMCFU(FIRST_SOIL_PT),SMCFF(FIRST_SOIL_PT)         
     &,              V_SAT(FIRST_SOIL_PT),HCONS(FIRST_SOIL_PT),LTIMER)
        DO J=1,SOIL_PTS
          I=SOIL_INDEX(J)
          HC(I,N)=HCONS(I)
        ENDDO
      ENDDO

!--------------------------------------------------------------------
! Calculate heat fluxes across layer boundaries
!--------------------------------------------------------------------
      DO N=1,NSHYD-1
        DO J=1,SOIL_PTS
          I=SOIL_INDEX(J)
          H_FLUX(I,N)=-HC(I,N)*2.0*(TSL(I,N+1)-TSL(I,N))
     &                            /(DZ(N+1)+DZ(N))
        ENDDO
      ENDDO
CDIR$ IVDEP
      DO J=1,SOIL_PTS
        I=SOIL_INDEX(J)
        H_FLUX(I,NSHYD)=0.0
        H_FLUX(I,0) = SOIL_SURF_HTF(I) + SNOW_SOIL_HTF(I)
      ENDDO

!--------------------------------------------------------------------
! Calculate the advection of heat by moisture fluxes
!--------------------------------------------------------------------
      DO N=2,NSHYD-1
        DO J=1,SOIL_PTS
          I=SOIL_INDEX(J)
          DHSLA(I,N)=TIMESTEP*HCAPW*DZ(N)*
     &    (W_FLUX(I,N-1)*(TSL(I,N-1)-TSL(I,N))/(DZ(N)+DZ(N-1))
     &     +W_FLUX(I,N)*(TSL(I,N)-TSL(I,N+1))/(DZ(N)+DZ(N+1)))
        ENDDO
      ENDDO

CDIR$ IVDEP
      DO J=1,SOIL_PTS
        I=SOIL_INDEX(J)
        DHSLA(I,1)=TIMESTEP*HCAPW*DZ(1)*
     &   (W_FLUX(I,0)*(TSL(I,0)-TSL(I,1))/DZ(1)
     &   +W_FLUX(I,1)*(TSL(I,1)-TSL(I,2))/(DZ(1)+DZ(2)))
        DHSLA(I,NSHYD)=TIMESTEP*HCAPW*DZ(NSHYD)*
     &               W_FLUX(I,NSHYD-1)*(TSL(I,NSHYD-1)-TSL(I,NSHYD))
     &               /(DZ(NSHYD)+DZ(NSHYD-1))
      ENDDO

      DO N=1,NSHYD
CDIR$ IVDEP
        DO J=1,SOIL_PTS
          I=SOIL_INDEX(J)
!-----------------------------------------------------------------------
! Calculate TMAX, the temperature above which all soil water is
! unfrozen. Check that (SMCLSAT/SMCL)**B will not overflow when SMCL is
! very small. The function EPSILON  gives the number of type (real) of 
! SMCL that is negligeable compared to 1.                              
!-----------------------------------------------------------------------
          IF (SMCL(I,N) .GT. SMCLSAT(I,N)*(EPSILON(SMCL(I,N)) ) )THEN 

            TMAX(I,N)=-SATHH(I)/DPSIDT
     &             *(SMCLSAT(I,N)/SMCL(I,N))**(B(I))
          ELSE
            TMAX(I,N)=-ZERODEGC
          ENDIF
          TMAX(I,N)=MAX(TMAX(I,N),-ZERODEGC)

          DHSL0(I,N)=TIMESTEP*(H_FLUX(I,N-1)-H_FLUX(I,N))
     &                 +DHSLA(I,N)


          DHSL(I,N)=DHSL0(I,N)

!-----------------------------------------------------------------------
! Calculate the soil temperature increments
!-----------------------------------------------------------------------
          TSL0(I,N)=TSL(I,N)
          SMCLF0(I,N)=SMCLF(I,N)
          SMCLU0(I,N)=SMCLU(I,N)

          IF (TSL(I,N).GT.TMAX(I,N)) THEN         ! All water unfrozen
            DTHU(I,N)=0.0
          ELSEIF (TSL(I,N).EQ.TMAX(I,N).AND.      ! Remains unfrozen
     &            DHSL(I,N).GE.0.0) THEN
            DTHU(I,N)=0.0
          ELSE                                    ! Phase changes
            DTHU(I,N)=DPSIDT*SMCLSAT(I,N)
     &               /(B(I)*SATHH(I)*RHO_WATER*DZ(N))
     &               *(-DPSIDT*TSL(I,N)/SATHH(I))**(-1.0/B(I)-1.0)
          ENDIF

          HCAPT(I)=HCAP(I)+(HCAPW-HCAPI)*SMCLU(I,N)/DZ(N)
     &            +HCAPI*SMCL(I,N)/DZ(N)
     &            +RHO_WATER*DTHU(I,N)*((HCAPW-HCAPI)*TSL(I,N)+LF)

          DTSL(I,N)=1.0/(HCAPT(I)*DZ(N))*DHSL(I,N)
          TSL(I,N)=TSL(I,N)+DTSL(I,N)

!-----------------------------------------------------------------------
! If the temperature increment is small and frozen water exists
! assume that the excess energy goes into phase change
!-----------------------------------------------------------------------
           IF (ABS(DTSL(I,N)).LT.TACUR.AND.
     &         TSL(I,N).LE.TMAX(I,N)) THEN
             DSMCLF(I,N)=-DHSL(I,N)/((HCAPW-HCAPI)*TSL(I,N)+LF)
             DSMCLF(I,N)=MAX(DSMCLF(I,N),-SMCLF(I,N))
             DSMCLF(I,N)=MIN(DSMCLF(I,N),SMCLU(I,N))
             SMCLU(I,N)=SMCLU(I,N)-DSMCLF(I,N)
             SMCLF(I,N)=SMCLF(I,N)+DSMCLF(I,N)
           ENDIF

!-----------------------------------------------------------------------
! Check to see if the discontinuity in HCAPT at TMAX has been crossed
! - if it has return to TMAX
!-----------------------------------------------------------------------
          IF ((TSL(I,N)-TMAX(I,N))*(TSL0(I,N)-TMAX(I,N))
     &        .LT.0.0) THEN
            DTSL(I,N)=DTSL(I,N)+(TMAX(I,N)-TSL(I,N))
            TSL(I,N)=TMAX(I,N)
          ENDIF

!-----------------------------------------------------------------------
! Diagnose unfrozen and frozen water contents
!-----------------------------------------------------------------------
          IF (TSL(I,N).GE.TMAX(I,N)) THEN
            SMCLU(I,N)=SMCL(I,N)
            SMCLF(I,N)=0.0
          ELSE
            SMCLU(I,N)=SMCLSAT(I,N)
     &                *(-DPSIDT*TSL(I,N)/SATHH(I))**(-1.0/B(I))
            SMCLF(I,N)=SMCL(I,N)-SMCLU(I,N)
          ENDIF

!-----------------------------------------------------------------------
! Calculate the error in heat conservation
!-----------------------------------------------------------------------
          DSMCLF(I,N)=SMCLF(I,N)-SMCLF0(I,N)
          DHSL(I,N)=DHSL(I,N)-(HCAP(I)*DZ(N)+HCAPW*SMCLU0(I,N)
     &                        +HCAPI*SMCLF0(I,N))*DTSL(I,N)
     &                   -DSMCLF(I,N)*((HCAPI-HCAPW)*TSL0(I,N)-LF)

!-----------------------------------------------------------------------
! Calculate the error in flux conservation
!-----------------------------------------------------------------------
          CEACUR(I)=ABS(DHSL(I,N))/TIMESTEP

        ENDDO

!-----------------------------------------------------------------------
! Define the array of points which fail to meet the flux criterion
!-----------------------------------------------------------------------
        J_ITER=0
C          CDIR$ IVDEP
        DO J=1,SOIL_PTS
          I=SOIL_INDEX(J)

          IF (CEACUR(I) .GT. FACUR) THEN
            J_ITER=J_ITER+1
            ITER_PTS(J_ITER)=J
          ENDIF

        ENDDO


!-----------------------------------------------------------------------
! Iterate on these points
!-----------------------------------------------------------------------
        DO M=1,MMAX-1
CDIR$ IVDEP
          DO J=1,J_ITER
            I=SOIL_INDEX(ITER_PTS(J))

            TSL0(I,N)=TSL(I,N)
            SMCLF0(I,N)=SMCLF(I,N)
            SMCLU0(I,N)=SMCLU(I,N)

            IF (TSL(I,N).GT.TMAX(I,N)) THEN         ! All water unfrozen
              DTHU(I,N)=0.0
            ELSEIF (TSL(I,N).EQ.TMAX(I,N).AND.      ! Remains unfrozen
     &              DHSL(I,N).GE.0.0) THEN
              DTHU(I,N)=0.0
            ELSE                                    ! Phase changes
              DTHU(I,N)=DPSIDT*SMCLSAT(I,N)
     &                 /(B(I)*SATHH(I)*RHO_WATER*DZ(N))
     &                 *(-DPSIDT*TSL(I,N)/SATHH(I))**(-1.0/B(I)-1.0)
            ENDIF

            HCAPT(I)=HCAP(I)+(HCAPW-HCAPI)*SMCLU(I,N)/DZ(N)
     &              +HCAPI*SMCL(I,N)/DZ(N)
     &              +RHO_WATER*DTHU(I,N)*((HCAPW-HCAPI)*TSL(I,N)+LF)

            DTSL(I,N)=1.0/(HCAPT(I)*DZ(N))*DHSL(I,N)
            TSL(I,N)=TSL(I,N)+DTSL(I,N)

!-----------------------------------------------------------------------
! If the temperature increment is small and frozen water exists
! assume that the excess energy goes into phase change
!-----------------------------------------------------------------------
            IF (ABS(DTSL(I,N)).LT.TACUR.AND.
     &          TSL(I,N).LE.TMAX(I,N)) THEN
              DSMCLF(I,N)=-DHSL(I,N)/((HCAPW-HCAPI)*TSL(I,N)+LF)
              DSMCLF(I,N)=MAX(DSMCLF(I,N),-SMCLF(I,N))
              DSMCLF(I,N)=MIN(DSMCLF(I,N),SMCLU(I,N))
              SMCLU(I,N)=SMCLU(I,N)-DSMCLF(I,N)
              SMCLF(I,N)=SMCLF(I,N)+DSMCLF(I,N)
            ENDIF

!-----------------------------------------------------------------------
! Check to see if the discontinuity in HCAPT at TMAX has been crossed
! - if it has return to TMAX
!-----------------------------------------------------------------------
            IF ((TSL(I,N)-TMAX(I,N))*(TSL0(I,N)-TMAX(I,N))
     &          .LT.0.0) THEN
              DTSL(I,N)=DTSL(I,N)+(TMAX(I,N)-TSL(I,N))
              TSL(I,N)=TMAX(I,N)
            ENDIF

!-----------------------------------------------------------------------
! Diagnose unfrozen and frozen water contents
!-----------------------------------------------------------------------
            IF (TSL(I,N).GE.TMAX(I,N)) THEN
              SMCLU(I,N)=SMCL(I,N)
              SMCLF(I,N)=0.0
            ELSE
              SMCLU(I,N)=SMCLSAT(I,N)
     &                  *(-DPSIDT*TSL(I,N)/SATHH(I))**(-1.0/B(I))
              SMCLF(I,N)=SMCL(I,N)-SMCLU(I,N)
            ENDIF

!-----------------------------------------------------------------------
! Calculate the error in heat conservation
!-----------------------------------------------------------------------
            DSMCLF(I,N)=SMCLF(I,N)-SMCLF0(I,N)
            DHSL(I,N)=DHSL(I,N)-(HCAP(I)*DZ(N)+HCAPW*SMCLU0(I,N)
     &                          +HCAPI*SMCLF0(I,N))*DTSL(I,N)
     &                     -DSMCLF(I,N)*((HCAPI-HCAPW)*TSL0(I,N)-LF)

          ENDDO
        ENDDO
      ENDDO

!-----------------------------------------------------------------------
! Diagnose soil temperatures (K) and fractional values of unfrozen and
! frozen water.
!-----------------------------------------------------------------------

      DO N=1,NSHYD
        DO J=1,SOIL_PTS
          I=SOIL_INDEX(J)
          TSOIL(I,N)=TSL(I,N)+ZERODEGC
          STHU(I,N)=SMCLU(I,N)/SMCLSAT(I,N)
          STHF(I,N)=SMCLF(I,N)/SMCLSAT(I,N)
        ENDDO
      ENDDO
      IF (LTIMER) THEN
        CALL TIMER('SOILHTC ',104)
      ENDIF

      RETURN
      END
