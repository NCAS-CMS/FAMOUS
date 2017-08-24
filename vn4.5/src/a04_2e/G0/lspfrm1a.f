C ******************************COPYRIGHT******************************
C (c) CROWN COPYRIGHT 1995, METEOROLOGICAL OFFICE, All Rights Reserved.
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
C*LL  SUBROUTINE LSP_FRMT-----------------------------------------------
CLL
CLL  Purpose: Adjust partition of cloud water between ice and liquid,
CLL           so that it is consistent with the temperature.  Then
CLL           freeze or melt precipitation falling into the layer from
CLL           above, when necessary.
CLL
CLL           In each case latent heating or cooling modifies the
CLL           temperature, and in both cases all the water is assumed
CLL           to undergo the phase change, unless this would take the
CLL           temperature the other side of freezing, in which case the
CLL           amount frozen or melted is limited to the amount needed
CLL           to take the temperature to zero C precisely.
CLL
CLL  Called by LS_PPN (P26) once for each model layer.
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  date
CLL
CLL  Programming standard: Unified Model Documentation Paper No 4,
CLL                        Version 1, dated 12/9/89.
CLL
CLL  System component covered: Part of P261.
CLL
CLL  Documentation: Unified Model Documentation Paper No 26.
C*
C*L  Arguments:---------------------------------------------------------
      SUBROUTINE LSP_FRMT
     +(RHODZ,TIMESTEP,POINTS,QCF,QCL,RAIN,SNOW,T)
      IMPLICIT NONE
      INTEGER POINTS        ! IN No. of gridpoints in batch.
      REAL
     + RHODZ(POINTS)        ! IN Mass of air in layer p.u.a. (kg/sq m).
      REAL TIMESTEP         ! IN Timestep (seconds).
      REAL
     + QCF(POINTS)          ! INOUT Cloud ice (kg water per kg air).
     +,QCL(POINTS)          ! INOUT Cloud liquid water (kg per kg air).
     +,RAIN(POINTS)         ! INOUT Rainfall rate (kg per sq m per s).
     +,SNOW(POINTS)         ! INOUT Snowfall rate (kg per sq m per s).
     +,T(POINTS)            ! INOUT Temperature (K).
C*
C*L  No workspace nor external subprograms required---------------------
C*  Common physical constants-------------------------------------------
C*L------------------COMDECK C_R_CP-------------------------------------
C R IS GAS CONSTANT FOR DRY AIR
C CP IS SPECIFIC HEAT OF DRY AIR AT CONSTANT PRESSURE
C PREF IS REFERENCE SURFACE PRESSURE
      REAL R,CP,KAPPA,PREF  ! PREFTOKI

      PARAMETER(R=287.05,
     &          CP=1005.,
     &          KAPPA=R/CP,
     &          PREF=100000.)
C*----------------------------------------------------------------------

C*L------------------COMDECK C_LHEAT------------------------------------
C LC IS LATENT HEAT OF CONDENSATION OF WATER AT 0DEGC
C LF IS LATENT HEAT OF FUSION AT 0DEGC
      REAL LC,LF

      PARAMETER(LC=2.501E6,
     &          LF=0.334E6)
C*----------------------------------------------------------------------
C*L------------------COMDECK C_O_DG_C-----------------------------------
C ZERODEGC IS CONVERSION BETWEEN DEGREES CELSIUS AND KELVIN
C TFS IS TEMPERATURE AT WHICH SEA WATER FREEZES
C TM IS TEMPERATURE AT WHICH FRESH WATER FREEZES AND ICE MELTS
      REAL ZERODEGC,TFS,TM

      PARAMETER(ZERODEGC=273.15,
     &          TFS=271.35,
     &          TM=273.15)
C*----------------------------------------------------------------------

C   Local physical constants ---------------------------------------
      REAL CPRLF,LFRCP
      PARAMETER (
     + LFRCP=LF/CP          ! Latent heat of fusion/Cp (K kg air/kg wat)
     +,CPRLF=1./LFRCP       ! Reciprocal of LFRCP.
     +)
C  Local variables------------------------------------------------------
C  (a) Real scalar effectively expanded to workspace by the Cray, using
C      vector registers.
      REAL
     + WPC                  ! LOCAL Amounts of water undergoing phase
C                           !       change.  2 different units are used.
C  (b) Other scalar.
      INTEGER I             ! Loop counter (horizontal field index).
C*
C-----------------------------------------------------------------------
CL  Loop round gridpoints.
C-----------------------------------------------------------------------
      DO 1 I=1,POINTS
C-----------------------------------------------------------------------
CL 1. Adjust cloud water and temperature to make them consistent.
CL    See equations P26.13 - P26.20.
C-----------------------------------------------------------------------
        IF(T(I).LE.TM)THEN
          WPC=MIN(QCL(I),CPRLF*(TM-T(I)))                      ! P26.13
          QCL(I)=QCL(I)-WPC                                    ! P26.15
          QCF(I)=QCF(I)+WPC                                    ! P26.16
          T(I)=T(I)+WPC*LFRCP                                  ! P26.14
        ELSE
          WPC=MIN(QCF(I),CPRLF*(T(I)-TM))                      ! P26.17
          QCL(I)=QCL(I)+WPC                                    ! P26.19
          QCF(I)=QCF(I)-WPC                                    ! P26.20
          T(I)=T(I)-WPC*LFRCP                                  ! P26.18
        ENDIF
C-----------------------------------------------------------------------
CL 2. Freeze or melt precipitation, on basis of updated temperature.
CL    See equations P26.21 - P26.28.
C-----------------------------------------------------------------------
        IF(T(I).LE.TM)THEN
          WPC=MIN(
     +            RAIN(I),
     +            CPRLF*(TM-T(I))*RHODZ(I)/TIMESTEP
     +            )                                            ! P26.21
          RAIN(I)=RAIN(I)-WPC                                  ! P26.23
          SNOW(I)=SNOW(I)+WPC                                  ! P26.24
          T(I)=T(I)+WPC*TIMESTEP*LFRCP/RHODZ(I)                ! P26.22
        ELSE
          WPC=MIN(
     +            SNOW(I),
     +            CPRLF*(T(I)-TM)*RHODZ(I)/TIMESTEP
     +            )                                            ! P26.21
          RAIN(I)=RAIN(I)+WPC                                  ! P26.23
          SNOW(I)=SNOW(I)-WPC                                  ! P26.24
          T(I)=T(I)-WPC*TIMESTEP*LFRCP/RHODZ(I)                ! P26.22
        ENDIF
    1 CONTINUE
      RETURN
      END
