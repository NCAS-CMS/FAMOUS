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
CLL  SUBROUTINE SURF_HYD-----------------------------------------------
CLL
CLL  PURPOSE : TO CARRY OUT CANOPY AND SURFACE HYDROLOGY CALCULATIONS
CLL
CLL            CANOPY WATER CONTENT IS DEPRECIATED BY EVAPORATION
CLL
CLL            SNOWMELT IS RUNOFF THE SURFACE WITHOUT INTERACTING
CLL            WITH THE CANOPY
CLL
CLL            THE CANOPY INTERCEPTION AND SURFACE RUNOFF OF
CLL            LARGE-SCALE RAIN IS CALCUALTED
CLL
CLL            THE CANOPY INTERCEPTION AND SURFACE RUNOFF OF
CLL            CONVECTIVE RAIN IS CALCUALTED
CLL
CLL
CLL  SUITABLE FOR SINGLE COLUMN MODEL USE
CLL
CLL  WRITTEN FOR CRAY-YMP BY S.ALLEN-BETT AND D.GREGORY
CLL
CLL  MODEL            MODIFICATION HISTORY FROM MODEL VERSION 3.0:
CLL VERSION  DATE
CLL  4.5    01/10/98  Removed old section-version defs. K Rogers
CLL   4.5    Jul. 98  Kill the IBM specific lines (JCThil)
CLL
CLL  PROGRAMMING STANDARDS : UNIFIED MODEL DOCUMENTATION PAPER NO. 4
CLL  VERSION NO. 1 18/1/90
CLL
CLL  SYSTEM TASK : P252
CLL
CLL  DOCUMENTATION : UNIFIED MODEL DOCUMENTATION PAPER NO 25
CLL
CLLEND-----------------------------------------------------------------
C
C*L  ARGUMENTS---------------------------------------------------------
C
      SUBROUTINE SURF_HYD (NPNTS,E_CANOPY,SNOW_MELT,LS_RAIN,
     *                     CON_RAIN,DSMC_DT,SURF_ROFF,CAN_WCNT,
     *                     CAN_CPY,INFIL,TOT_TFALL,TIMESTEP)
C
      IMPLICIT NONE
C
C-----------------------------------------------------------------------
C VECTOR LENGTHS
C-----------------------------------------------------------------------
C
      INTEGER NPNTS         ! IN VECTOR LENGTH
C
C
C-----------------------------------------------------------------------
C VARIABLES WHICH ARE INPUT
C-----------------------------------------------------------------------
C
      REAL TIMESTEP         ! IN MODEL TIMESTEP (S)
C
      REAL E_CANOPY(NPNTS)  ! IN CANOPY EVAPORATION (KG/M2/S)
C
      REAL SNOW_MELT(NPNTS) ! IN SNOW MELT (KG/M2/S)
C
      REAL LS_RAIN(NPNTS)   ! IN LARGE-SCALE RAIN (KG/M2/S)
C
      REAL CON_RAIN(NPNTS)  ! IN CONVECTIVE RAIN (KG/M2/S)
C
      REAL CAN_CPY(NPNTS)   ! IN CANOPY CAPACITY (KG/M2)
C
      REAL INFIL(NPNTS)     ! IN INFILTRATION RATE(KG/M2/S)
C
C
C-----------------------------------------------------------------------
C VARIABLES WHICH ARE INPUT AND OUTPUT
C-----------------------------------------------------------------------
C
      REAL CAN_WCNT(NPNTS)  ! INOUT
                            ! IN CANOPY WATER CONTENT (KG/M2)
                            ! OUT UPDATED CANOPY WATER CONTENT
                            !     (KG/M2)
C
C
C-----------------------------------------------------------------------
C VARIABLES WHICH ARE OUTPUT
C-----------------------------------------------------------------------
C
      REAL DSMC_DT(NPNTS)   ! OUT CUMULATIVE RATE OF CHANGE OF SOIL
                            !     MOISTURE CONTENT WITH TIMES
                            !     (KG/M2/S)
C
      REAL SURF_ROFF(NPNTS) ! OUT CUMULATIVE SURFACE RUNOFF (KG/M2/S)
C
      REAL TOT_TFALL(NPNTS) ! OUT TOTAL CANOPY THROUGHFALL (KG/M2/S)
C
C
C-----------------------------------------------------------------------
C VARIABLES WHICH ARE DEFINED LOCALLY
C
C ON THE CRAY ARRAYS ARE DYNAMICALLY ALLOCATED
C
C WORK SPACE USAGE : 3 * NPNTS
C
C-----------------------------------------------------------------------
C
C
      REAL CAN_COND(NPNTS)  ! CANOPY CONDENSATION (KG/M2/S)
C
      REAL RUNOFF(NPNTS)    ! SURFACE RUNOFF FROM SINGLE WATER
                            ! TYPE (KG/M2/S)
C
      REAL TFALL(NPNTS)     ! THROUGHFALL FROM SINGLE WATER
                            ! TYPE (KG/M2/S)
C
C
C-----------------------------------------------------------------------
C INTERNAL LOOP COUNTERS
C-----------------------------------------------------------------------
C
      INTEGER I             ! LOOP COUNTER
C
C
C----------------------------------------------------------------------
C EXTERNAL SUBROUTINE CALLS
C----------------------------------------------------------------------
C
      EXTERNAL FRUNOFF,SIEVE
C
C*--------------------------------------------------------------------
C
      DO 10 I=1,NPNTS
CL
CL----------------------------------------------------------------------
CL DEPRECIATE CANOPY WATER CONTENT BY EVAPORATION
CL
CL UM DOCUMENTATION PAPER P25
CL SECTION (3), EQN(P252.17)
CL----------------------------------------------------------------------
CL
        IF (E_CANOPY(I) .GT. 0.0)
     *   CAN_WCNT(I) = CAN_WCNT(I)-E_CANOPY(I)*TIMESTEP
C
C-----------------------------------------------------------------------
C NEXT IF TEST INTRODUCED TO PREVENT CANOPY WATER CONTENT BECOMING
C SLIGHTLY NEGATIVE DUE TO ROUNDING ERROR
C-----------------------------------------------------------------------
C
        IF (CAN_WCNT(I) . LT. 0.0) CAN_WCNT(I) = 0.0
C
C-----------------------------------------------------------------------
C  ZERO CUMULATIVE STORES
C-----------------------------------------------------------------------
C
        TOT_TFALL(I) = 0.0
        SURF_ROFF(I) = 0.0
        DSMC_DT(I)   = 0.0
  10  CONTINUE
CL
CL----------------------------------------------------------------------
CL CALCULATION OF SURFACE RUN0FF OF SNOWMELT
CL
CL UM DOCUMENTATION PAPER NO 25
CL SECTION (3)
CL----------------------------------------------------------------------
CL
       CALL FRUNOFF (NPNTS,SNOW_MELT,SNOW_MELT,CAN_CPY,CAN_CPY,INFIL,   
     *               1.0,RUNOFF,TIMESTEP)                               
C
      DO 20 I=1,NPNTS
        IF(SNOW_MELT(I) .GT. 0.0)THEN
C
C----------------------------------------------------------------------
C CUMULATE SURFACE RUNOFF AND RATE OF CHANGE OF
C SOIL MOISTURE CONTENT DUE TO INFILTRATION OF SURFACE WATER
C----------------------------------------------------------------------
C
          SURF_ROFF(I) = RUNOFF(I)
          DSMC_DT(I) = (SNOW_MELT(I) - RUNOFF(I))
        END IF
C
C-----------------------------------------------------------------------
C DEFINE CANOPY CONDENSATION (WHEN CANOPY EVAPORATION IS NEGATIVE)
C-----------------------------------------------------------------------
C
          IF (E_CANOPY(I) .LT. 0.0) THEN
           CAN_COND(I) = -E_CANOPY(I)
          ELSE
           CAN_COND(I) =0.0
          END IF
  20  CONTINUE
CL
CL----------------------------------------------------------------------
CL CALCULATE CANOPY INTERCEPTION, THROUGHFALL AND SURFACE RUNOFF FOR
CL CANOPY CONDENSATION
CL
CL UM DOCUMENTATION PAPER NO 25
CL SECTION (3)
CL
CL FRACTIONAL AREA OF GRID BOX WHERE
CL CANOPY CONDENSATION IS ASSUMED TO FALL = 1.0
CL----------------------------------------------------------------------
CL
      CALL SIEVE (NPNTS,CAN_COND,CAN_WCNT,CAN_CPY,1.0,TFALL,
     *           TIMESTEP)
C
      CALL FRUNOFF (NPNTS,CAN_COND,TFALL,CAN_WCNT,CAN_CPY,INFIL,1.0,   
     *              RUNOFF,TIMESTEP)                      
C
      DO 30 I=1,NPNTS
        IF (CAN_COND(I) .GT. 0.0) THEN
C
C-----------------------------------------------------------------------
C UPDATE CANOPY WATER CONTENT FOR INTERCEPTION OF CANOPY
C CONDENSATION
C
C UM DOCUMENTATION PAPER NO 25
C SECTION (3(I)), EQN(P252.10)
C-----------------------------------------------------------------------
C
          CAN_WCNT(I) = CAN_WCNT(I) + (CAN_COND(I) - TFALL(I))*TIMESTEP
C
C----------------------------------------------------------------------
C CUMULATE THROUGHFALL, SURAFCE RUNOFF AND RATE OF CHANGE OF
C SOIL MOISTURE CONTENT DUE TO INFILTRATION OF SURFACE WATER
C----------------------------------------------------------------------
C
          TOT_TFALL(I) = TFALL(I)
          SURF_ROFF(I) = SURF_ROFF(I) + RUNOFF(I)
          DSMC_DT(I) = DSMC_DT(I) + (TFALL(I) - RUNOFF(I))
C
        END IF
  30  CONTINUE
CL
CL----------------------------------------------------------------------
CL CALCULATE CANOPY INTERCEPTION, THROUGHFALL AND SURFACE RUNOFF FOR
CL LARGE-SCALE RAIN
CL
CL UM DOCUMENTATION PAPER NO 25
CL SECTION (3)
CL
CL FRACTIONAL AREA OF GRID BOX WHERE
CL LARGE-SCALE RAIN IS ASSUMED TO FALL = 0.5
CL----------------------------------------------------------------------
CL
      CALL SIEVE (NPNTS,LS_RAIN,CAN_WCNT,CAN_CPY,1.0,TFALL,
     *            TIMESTEP)
C
      CALL FRUNOFF (NPNTS,LS_RAIN,TFALL,CAN_WCNT,CAN_CPY,INFIL,1.0,   
     *              RUNOFF,TIMESTEP)                         
C
      DO 40 I=1,NPNTS
        IF (LS_RAIN(I) .GT. 0.0) THEN
C
C-----------------------------------------------------------------------
C UPDATE CANOPY WATER CONTENT FOR INTERCEPTION OF
C LARGE-SCALE RAIN
C
C UM DOCUMENTATION PAPER  NO 25
C SECTION (3(I)), EQN(P252.10)
C-----------------------------------------------------------------------
C
          CAN_WCNT(I) = CAN_WCNT(I) + (LS_RAIN(I) - TFALL(I))*TIMESTEP
C
C----------------------------------------------------------------------
C CUMULATE THROUGHFALL, SURAFCE RUNOFF AND RATE OF CHANGE OF
C SOIL MOISTURE CONTENT DUE TO INFILTRATION OF SURFACE WATER
C----------------------------------------------------------------------
C
          TOT_TFALL(I) = TOT_TFALL(I)+TFALL(I)
          SURF_ROFF(I) = SURF_ROFF(I)+RUNOFF(I)
          DSMC_DT(I) = DSMC_DT(I) + (TFALL(I) - RUNOFF(I))
C
        END IF
C
  40  CONTINUE
CL
CL----------------------------------------------------------------------
CL CALCULATE CANOPY INTERCEPTION, THROUGHFALL AND SURFACE RUNOFF FOR
CL CONVECTIVE RAIN
CL
CL UM DOCUMENTATION PAPER NO 25
CL SECTION (3)
CL
CL FRACTIONAL AREA OF GRID BOX WHERE
CL CONVECTIVE RAIN IS ASSUMED TO FALL = 0.1
CL----------------------------------------------------------------------
CL
      CALL SIEVE (NPNTS,CON_RAIN,CAN_WCNT,CAN_CPY,0.3,TFALL,
     *            TIMESTEP)
C
      CALL FRUNOFF (NPNTS,CON_RAIN,TFALL,CAN_WCNT,CAN_CPY,INFIL,0.3,   
     *              RUNOFF,TIMESTEP)          
C
      DO 50 I= 1,NPNTS
        IF (CON_RAIN(I) .GT. 0.0) THEN
C
C-----------------------------------------------------------------------
C UPDATE CANOPY WATER CONTENT FOR INTERCEPTION OF
C CONVECTIVE RAIN
C
C UM DOCUMENTATION PAPER NO 25
C SECTION (3(I)), EQN(P252.10)
C-----------------------------------------------------------------------
C
          CAN_WCNT(I) = CAN_WCNT(I) + (CON_RAIN(I) - TFALL(I))*TIMESTEP
C
C----------------------------------------------------------------------
C CUMULATE THROUGHFALL, SURAFCE RUNOFF AND RATE OF CHANGE OF
C SOIL MOISTURE CONTENT DUE TO INFILTRATION OF SURFACE WATER
C----------------------------------------------------------------------
C
          TOT_TFALL(I) = TOT_TFALL(I)+TFALL(I)
          SURF_ROFF(I) = SURF_ROFF(I)+RUNOFF(I)
          DSMC_DT(I) = DSMC_DT(I) + (TFALL(I) - RUNOFF(I))
C
        END IF
C
  50  CONTINUE
C
      RETURN
      END
