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
CLL  SUBROUTINE HYDROL-------------------------------------------------
CLL
CLL  PURPOSE : CALLS SUBROUTINE SFSNOM TO CALCULATE SNOW PROCESSES
CLL            AT THE SURFACE
CLL
CLL            CALLS SUBROUTINE SURF_HYD TO CALCULATE CANOPY
CLL            AND SURFACE HYDROLOGY
CLL
CLL            CALLS SUBROUTINE SOIL_HYD TO CALCULATE SOIL HYDROLOGY
CLL
CLL  SUITABLE FOR SINGLE COLUMN MODEL USE
CLL
CLL  WRITTEN FOR CYBER/ETA-10 BY S.ALLEN AND D.GREGORY
CLL
CLL  MODEL            MODIFICATION HISTORY FROM MODEL VERSION 3.0:
CLL VERSION  DATE
CLL   4.4  29/10/97  MODIFIED FOR PROGNOSTIC SNOW ALBEDO SCHEME
CLL                                                   R. ESSERY
CLL   4.5    Jul. 98  Kill the IBM specific lines (JCThil)
CLL
CLL  PROGRAMMING STANDARDS : UNIFIED MODEL DOCUMENTATION PAPER NO. 4
CLL  VERSION NO. 1 DATED 18/1/90
CLL
CLL  LOGICAL COMPONENTS COVERED: P27
CLL
CLL  SYSTEM TASK :
CLL
CLL  DOCUMENTATION : UNIFIED MODEL DOCUMENTATION PAPER NO 25
CLL
CLLEND-----------------------------------------------------------------
C
C*L  ARGUMENTS---------------------------------------------------------
C
      SUBROUTINE HYDROL (NPNTS,E_CANOPY,LS_RAIN,CON_RAIN,
     *                    CAN_CPY,INFIL,CON_SNOW,
     *                    HCAP,HCON,LS_SNOW,SNOWSUB,
     *                    ROOTD,BS,C_EAG,SOILMC,
     *                    CAN_WCNT,RGRAIN,L_SNOW_ALBEDO,
     *                    SNOW_DEP,TSTAR,SNOW_MELT,
     *                    TOT_TFALL,SURF_ROFF,TIMESTEP,
     *                    E_SOIL
     *                   ,HF_SNOW_MELT,STF_HF_SNOW_MELT
     *                   ,SUB_SURF_ROFF,STF_SUB_SURF_ROFF
     *                    )
C
      IMPLICIT NONE
C
C
C-----------------------------------------------------------------------
C VECTOR LENGTHS
C-----------------------------------------------------------------------
C
      INTEGER NPNTS              ! IN FULL VECTOR LENGTH
C
C
C-----------------------------------------------------------------------
C VARIABLES WHICH ARE INPUTS
C
C NOTE THAT FULL REFERS TO FULL LENGTH DATA
C-----------------------------------------------------------------------
C
      REAL TIMESTEP              ! IN MODEL TIMESTEP (S)
C
      REAL E_SOIL(NPNTS)         ! IN SURFACE EVAPORATION
                                 !    (KG/M2/S)
C
      REAL E_CANOPY(NPNTS)       ! IN CANOPY EVAPORATION(KG/M2/S)
C
      REAL LS_RAIN(NPNTS)        ! IN LARGE-SCALE RAIN (KG/M2/S)
C
      REAL CON_RAIN(NPNTS)       ! IN CONVECTIVE RAIN (KG/M2/S)
C
      REAL CON_SNOW(NPNTS)       ! IN CONVECTIVE SNOWFALL (KG/M2/S)
C
      REAL LS_SNOW(NPNTS)        ! IN LARGE-SCALE SNOWFALL (KG/M2/S)
C
      REAL SNOWSUB(NPNTS)        ! IN SUBLIMATION OF LYING SNOW
                                 !    (KG/M2/S)
C
      REAL CAN_CPY(NPNTS)        ! IN CANOPY CAPACITY (KG/M2)
C
      REAL INFIL(NPNTS)          ! IN RATE AT WHICH WATER INFILTRATES
                                 !    INTO THE SOIL
C
      REAL HCAP(NPNTS)           ! IN SOIL HEAT CAPACITY (J/K/M3)
C
      REAL HCON(NPNTS)           ! IN SOIL THERMAL CONDUCTIVITY (W/M/K)
C
      REAL ROOTD(NPNTS)          ! IN ROOT DEPTH (M)
C
      REAL BS(NPNTS)             ! IN BS USED IN CALCULATION OF
                                 !    SUB-SURFACE RUNOFF (SEE EQN
                                 !    (P253.4), UM DOC PAPER NO 25)
C
      REAL C_EAG(NPNTS)          ! IN EXPONENT USED IN CALCULATION OF
                                 !    SUB-SURFACE RUNOFF USING
                                 !    EAGLESON' EMPIRICAL FORMULA
                                 !    ( SEE EQN(P253.4), UM DOC
                                 !    PAPER NO 25 )
C
      LOGICAL L_SNOW_ALBEDO      ! IN FLAG FOR PROGNOSTIC SNOW ALBEDO   
C                                                                       
      LOGICAL STF_HF_SNOW_MELT   ! IN STASH FLAG FOR SNOW MELT
                                 ! HEAT FLUX
C
      LOGICAL STF_SUB_SURF_ROFF  ! IN STASH FLAG FOR SUB-SURFACE
                                 !    RUNOFF
C
C
C-----------------------------------------------------------------------
C VARIABLES WHICH ARE INPUT AND OUTPUT
C-----------------------------------------------------------------------
C
      REAL SOILMC(NPNTS)         ! INOUT
                                 ! IN SOIL MOISTURE CONTENT
                                 !    (KG/M2)
                                 ! OUT UPDATED SOIL MOISTURE
                                 !     CONTENT (KG/M2)
C
      REAL CAN_WCNT(NPNTS)       ! INOUT
                                 ! IN CANOPY WATER CONTENT
                                 !    (KG/M2)
                                 ! OUT UPDATED CANOPY WATER
                                 !     CONTENT (KG/M2)
C
      REAL RGRAIN(NPNTS)         ! INOUT Snow grain size (microns).
C                                                                       
      REAL SNOW_DEP(NPNTS)       ! INOUT
                                 ! IN SNOW DEPTH (KG OF WATER/M2)
                                 ! OUT UPDATED SNOW DEPTH
                                 !     (KG OF WATER/M2)
C
      REAL TSTAR(NPNTS)          ! INOUT
                                 ! IN SURFACE TEMPERATURE (K)
                                 ! OUT UPDATED SURFACE TEMPERATURE (K)
C
C
C-----------------------------------------------------------------------
C VARIABLES WHICH ARE OUTPUT
C-----------------------------------------------------------------------
C
      REAL TOT_TFALL(NPNTS)      ! OUT TOTAL THROUGHFALL (KG/M2/S)
C
      REAL SURF_ROFF(NPNTS)      ! OUT SURFACE RUNOFF (KG/M2/S)
C
      REAL SNOW_MELT(NPNTS)      ! OUT SNOWMELT (KG/M2/S)
C
      REAL HF_SNOW_MELT(NPNTS)   ! OUT SNOWMELT HEAT FLUX (WATTS/M2)
C
      REAL SUB_SURF_ROFF(NPNTS)  ! OUT 'SLOW' RUNOFF (KG/M2/S)
C
C
C-----------------------------------------------------------------------
C VARIABLES WHICH ARE DEFINED LOCALLY
C-----------------------------------------------------------------------
C
      REAL DSMC_DT(NPNTS)        ! RATE OF CHANGE OF
                                 ! SOIL MOISTURE CONTENT DUE TO
                                 ! WATER FALLING ONTO THE SURFACE
                                 ! AFTER SURFACE RUNOFF
                                 ! (KG/M2/S)
C
C
C----------------------------------------------------------------------
C  EXTERNAL ROUTINE CALLS
C----------------------------------------------------------------------
C
      EXTERNAL SFSNOW,SURF_HYD,SOIL_HYD
C
C*---------------------------------------------------------------------
CL
CL----------------------------------------------------------------------
CL CALL SFSNOW
CL
CL CARRIES OUT INCREMENT TO THE SNOW DEPTH AND
CL SNOW MELT CALCULATIONS
CL
CL UM DOCUMENTATION PAPER NO 25
CL SECTION (2)
CL
CL COMPONENT P251
CL----------------------------------------------------------------------
CL
      CALL SFSNOW(CON_SNOW,HCAP,HCON,LS_SNOW,SNOWSUB,
     *             TIMESTEP,NPNTS,RGRAIN,L_SNOW_ALBEDO,
     *             SNOW_DEP,TSTAR,SNOW_MELT
     *            ,HF_SNOW_MELT,STF_HF_SNOW_MELT
     *            )
CL
CL----------------------------------------------------------------------
CL CALL SURF_HYD
CL
CL CARRIES OUT CANOPY AND SURFACE HYDROLOGY CALCULATIONS
CL
CL UM DOCUMENTATION PAPER NO 25
CL SECTION (3)
CL
CL COMPONENT P252
CL----------------------------------------------------------------------
CL
      CALL SURF_HYD(NPNTS,E_CANOPY,SNOW_MELT,LS_RAIN,
     *              CON_RAIN,DSMC_DT,SURF_ROFF,CAN_WCNT,
     *              CAN_CPY,INFIL,TOT_TFALL,TIMESTEP)
CL
CL----------------------------------------------------------------------
CL CALL SOIL_HYD
CL
CL CARRIES OUT SOIL HYDROLOGY CALCULATIONS
CL
CL UM DOCUMENTATION PAPER NO 25
CL SECTION (4)
CL
CL COMPONENT P253
CL----------------------------------------------------------------------
CL
      CALL SOIL_HYD(ROOTD,BS,C_EAG,E_SOIL,DSMC_DT,
     *              TIMESTEP,NPNTS,SOILMC
     *             ,SUB_SURF_ROFF,STF_SUB_SURF_ROFF
     *             )
      RETURN
      END
