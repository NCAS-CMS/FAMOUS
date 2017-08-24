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
C*LL  SUBROUTINE SOIL_HYD ----------------------------------------------
CLL
CLL  Purpose:  Firstly, the soil moisture content (SMC) is incremented
CLL            by the infiltration amount (FW) calculated at P252 and
CLL            decremented by the evaporation (E_SOIL) calculated at
CLL            P24.  Secondly, the sub-surface ("slow") runoff of soil
CLL            moisture is calculated, and the SMC decreased
CLL            appropriately.
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  Date
CLL   4.5    Jul. 98  Kill the IBM specific lines (JCThil)
CLL
CLL  Programming standard: Unified Model Documentation Paper No.4
CLL                        version no. 2, dated 18/1/90.
CLL
CLL  Logical component covered: P253.
CLL
CLL  System task:
CLL
CLL  Documentation: UM Documentation Paper No 25.
CLLEND------------------------------------------------------------------
C*
C*L  ARGUMENTS:---------------------------------------------------------
      SUBROUTINE SOIL_HYD (
     + ROOT_DEPTH,BS,C_EAG,E_SOIL,FW,TIMESTEP,POINTS,SMC
     +,SLOW_RUNOFF,STF_SLOW_RUNOFF
     +)
      IMPLICIT NONE
      INTEGER POINTS
      REAL
     + E_SOIL(POINTS)     ! IN Evaporation from soil surface (+ve
C                         !    upwards) (kg per sq m per sec).
     +,FW(POINTS)         ! IN Throughfall from canopy plus snowmelt
C                         !    minus surface runoff (kg per sq m per s).
     +,ROOT_DEPTH(POINTS) ! IN Root depth (m).
     +,BS(POINTS)         ! IN Coefficient in conductivity calculation,
C                         !    defined by eq P253.4 (kg per sq m per s).
C                         !    This is BS in the documentation.
     +,C_EAG(POINTS)      ! IN Exponent in conductivity calculation.
C                         !    This is C in the documentation.
      REAL
     + TIMESTEP           ! IN Timestep in seconds.
      REAL
     + SMC(POINTS)        ! INOUT Soil moisture content (kg per sq m).
      REAL
     + SLOW_RUNOFF(POINTS) ! OUT "Slow" runoff (kg per sq m per sec).
      LOGICAL STF_SLOW_RUNOFF ! IN stash flag for slow runoff
C-----------------------------------------------------------------------
C   workspace required.
C   one array required.
C  No EXTERNAL subprograms are called.
C-----------------------------------------------------------------------
      REAL SRO_DT(POINTS)  ! Slow runoff (kg per sq m per timestep).
C*
C  Model constants
      REAL RHO_WATER              ! DENSITY OF WATER (KG/M3)
C
      PARAMETER (RHO_WATER = 1000.0)
C

C
C-----------------------------------------------------------------------
C  Define local variables.
C-----------------------------------------------------------------------
      INTEGER I            ! Loop counter (horizontal field index).
C-----------------------------------------------------------------------
CL     No significant structure, apart from loop round points :-
C-----------------------------------------------------------------------
      DO 1 I=1,POINTS
C-----------------------------------------------------------------------
CL  1. Modify SMC as a result of infiltration and evaporation.  P253.1.
C-----------------------------------------------------------------------
        SMC(I) = SMC(I) + TIMESTEP*(FW(I) - E_SOIL(I))
        SRO_DT(I) = 0.0
        IF(SMC(I).GT.0.0 .AND. ROOT_DEPTH(I).GT.0.0)THEN
C-----------------------------------------------------------------------
CL  2. Calculate subsurface runoff integrated over timestep, in kg per
CL     sq metre (SRO_DT), eqn P253.6 in different units.  This is
CL     equivalent to calculating the actual hydraulic conductivity of
CL     the soil as a function of SMC (see DCTN 38, eqn 3.10).
CL     Then update the SMC (eqn P253.7).
C-----------------------------------------------------------------------
          SRO_DT(I) = BS(I)
     +            * (( SMC(I) /(RHO_WATER*ROOT_DEPTH(I)))**C_EAG(I))
     +            * TIMESTEP
          IF (SRO_DT(I).GT.SMC(I)) SRO_DT(I)=SMC(I)
          SMC(I) = SMC(I) - SRO_DT(I)
        ENDIF
    1 CONTINUE
C-----------------------------------------------------------------------
CL  3. Convert runoff to kg per sq metre per second for diagnostic.
C-----------------------------------------------------------------------
      IF (STF_SLOW_RUNOFF) THEN
        DO I=1,POINTS
          SLOW_RUNOFF(I) = SRO_DT(I) / TIMESTEP
        ENDDO
      ENDIF
      RETURN
      END
