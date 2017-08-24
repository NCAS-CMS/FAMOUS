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
CLL  SUBROUTINE CON_RAD------------------------------------------------
CLL
CLL  PURPOSE : CALCULATES CONVECTIVE CLOUD TOP, BASE AND
CLL            AMOUNT
CLL
CLL  SUITABLE FOR SINGLE COLUMN MODEL USE
CLL
CLL  CODE REWORKED FOR CRAY Y-MP BY D.GREGORY AUTUMN/WINTER 1989/90
CLL
CLL  MODEL            MODIFICATION HISTORY FROM MODEL VERSION 3.0:
CLL VERSION  DATE
CLL  3.3  23/12/93 Change to cloud top because of
CLL                change to detrainment rate calculation.  D.Gregory.
CLL
CLL  3.4  21/03/94  Add lowest conv.cloud diagnostics.  R.T.H.Barnes.
CLL
CLL  4.4  26/09/97  Pass in extra cloud water variable to allow rain
CLL                 out in CLOUDW before calculation of water path
CLL                 if L_CCW is set to .TRUE. in CLOUDW.      J.M.G.
CLL
CLL  PROGRAMMING STANDARDS : UNIFIED MODEL DOCUMENTATION PAPER NO. 4
CLL  VERSION NO. 1
CLL
CLL  LOGICAL COMPONENT NUMBER: P27
CLL
CLL  SYSTEM TASK : P27
CLL
CLL  DOCUMENTATION : UNIFIED MODEL DOCUMENTATION PAPER P27
CLL                  SECTION (9)
CLL
CLLEND-----------------------------------------------------------------
C
C*L  ARGUMENTS---------------------------------------------------------
C
      SUBROUTINE CON_RAD (K,XPK,XPKP1,FLXKP1,BTERM,CCA,ICCB,ICCT,TCW,
     *              CCW,CCLWP,DELPKP1,LCCA,LCBASE,LCTOP,LCCLWP,NPNTS)
C
      IMPLICIT NONE
C
C
C----------------------------------------------------------------------
C MODEL CONSTANTS
C----------------------------------------------------------------------
C
C*L------------------COMDECK C_G----------------------------------------
C G IS MEAN ACCEL DUE TO GRAVITY AT EARTH'S SURFACE
      REAL G

      PARAMETER(G=9.80665)
C*----------------------------------------------------------------------
C
C----------------------------------------------------------------------
C VECTOR LENGTH AND LOOP VARIABLES
C----------------------------------------------------------------------
C
      INTEGER NPNTS        ! IN VECTOR LENGTH
C
      INTEGER K            ! IN PRESENT MODEL LAYER
C
      INTEGER I            ! LOOP COUNTER
C
C
C----------------------------------------------------------------------
C VARIABLES WHICH ARE INPUT
C----------------------------------------------------------------------
C
      REAL XPK(NPNTS)      ! IN PARCEL CLOUD WATER IN LAYER K (KG/KG)
C
      REAL XPKP1(NPNTS)    ! IN PARCEL CLOUD WATER IN LAYER K+1 (KG/KG)
C
      LOGICAL BTERM(NPNTS) ! IN MASK FOR POINTS WHERE CONVECTION
                           !    IS ENDING
C
      REAL FLXKP1(NPNTS)   ! IN PARCEL MASSFLUX IN LAYER K+1 (PA/S)
C
      REAL DELPKP1(NPNTS)  ! IN PRESSURE DIFFERENCE ACROSS LAYER K+1
C
      REAL CCW(NPNTS)      ! IN PARCEL CLOUD WATER AS CALCULATED BEFORE
                           !    PRECIPITATION. LAYER K+1 (KG/KG)
C
C
C----------------------------------------------------------------------
C VARIABLES WHICH ARE INPUT AND OUTPUT
C----------------------------------------------------------------------
C
      REAL TCW(NPNTS)      ! INOUT
                           ! IN  TOTAL CONDENSED WATER SUMMED TO
                           !     LAYER K (KG/M**2/S)
                           ! OUT TOTAL CONDENSED WATER SUMMED TO
                           !     LAYER K+1 OR IF CONVECTION HAS
                           !     TERMINATED ZEROED (KG/M**2/S)
C
      REAL CCLWP(NPNTS)    ! INOUT
                           ! IN  TOTAL CLOUD LIQUID WATER PATH
                           !     SUMMED TO LAYER K  (KG/M**2)
                           ! OUT TOTAL CLOUD LIQUID WATER PATH
                           !     SUMMED TO LAYER K+1 (KG/M**2)
      REAL LCCA(NPNTS)      ! INOUT LOWEST CONV.CLOUD AMOUNT (%)
C
      INTEGER LCBASE(NPNTS) ! INOUT LOWEST CONV.CLOUD BASE LEVEL
C
      INTEGER LCTOP(NPNTS)  ! INOUT LOWEST CONV.CLOUD TOP LEVEL
C
C
C----------------------------------------------------------------------
C VARIABLES WHICH ARE AND OUTPUT
C----------------------------------------------------------------------
C
      REAL CCA(NPNTS)      ! OUT CONVECTIVE CLOUD AMOUNT (%)
C
      INTEGER ICCB(NPNTS)   ! OUT CONVECTIVE CLOUD BASE LEVEL
C
      INTEGER ICCT(NPNTS)   ! OUT CONVECTIVE CLOUD TOP LEVEL
C
      REAL LCCLWP(NPNTS)    ! OUT LOWEST CONV.CLOUD LIQ.WATER PATH
C
C
C*---------------------------------------------------------------------
CL
CL---------------------------------------------------------------------
CL CALCULATE CLOUD BASE and Lowest Cloud Base
CL
CL WHEN CLOUD BASE SET ZERO TOTAL CONDENSED WATER
CL---------------------------------------------------------------------
CL
      DO  I = 1,NPNTS
        IF ( XPK(I) .LE. 0.0 .AND. CCW(I) .GT. 0 ) THEN
          ICCB(I)=K+1
          CCLWP(I)=0.0
        END IF

        IF ( XPK(I) .LE. 0.0 .AND. CCW(I) .GT. 0.0 .AND.
     &       LCBASE(I) .EQ. 0 ) THEN
          LCBASE(I)=K+1
          LCCLWP(I)=0.0
        END IF
CL
CL---------------------------------------------------------------------
CL CALCULATE CLOUD TOP and Lowest Cloud Top
CL---------------------------------------------------------------------
CL
        IF (BTERM(I) .AND.
     *      ((CCW(I).GT.0.0).OR.(XPK(I).GT.0.0)) ) ICCT(I) = K+1

        IF (BTERM(I) .AND.  LCTOP(I).EQ.0 .AND.
     *      ((CCW(I).GT.0.0).OR.(XPK(I).GT.0.0)) ) THEN
          LCTOP(I) = K+1
        END IF
C
        IF ( FLXKP1(I) .GT. 0.0) THEN
CL
CL---------------------------------------------------------------------
CL SUM TOTAL CONDENSED WATER PER SECOND - ASSUMES THAT THE INITIAL
CL CONVECTIVE LAYER IS UNSATURATED
CL---------------------------------------------------------------------
CL
          TCW(I) = TCW(I) + FLXKP1(I) * CCW(I) / G
CL
CL---------------------------------------------------------------------
CL SUM CONV CONDENSED WATER PATH - ASSUMES THAT THE INITIAL
CL CONVECTIVE LAYER IS UNSATURATED
CL---------------------------------------------------------------------
CL
          CCLWP(I) = CCLWP(I) + XPKP1(I) * DELPKP1(I) / G
CL
CL---------------------------------------------------------------------
CL SUM CONV CONDENSED WATER PATH up to lowest conv.cloud
CL ASSUMES THAT THE INITIAL CONVECTIVE LAYER IS UNSATURATED
CL---------------------------------------------------------------------
CL
          IF (LCCA(I).LE.0.0) THEN
            LCCLWP(I) = LCCLWP(I) + CCW(I) * DELPKP1(I) / G
          END IF
C
        END IF
CL
CL---------------------------------------------------------------------
CL CALCULATE CONVECTIVE CLOUD AMOUNT IF CONVECTION TERMINATES IN
CL LAYER K AND TOTAL CONDENSED WATER PATH OVER A TIME STEP
CL
CL UM DOCUMENTATION PAPER P27
CL SECTION (9), EQUATION (37)
CL---------------------------------------------------------------------
CL
        IF( BTERM(I) .AND. TCW(I).GT.0.0 ) THEN
C
          IF ( TCW(I) .LT. 2.002E-6 ) TCW(I) = 2.002E-6
C
          CCA(I) = 0.7873 + 0.06 * LOG(TCW(I))
          IF (CCA(I) .GT. 1.0) CCA(I) = 1.0
C
          IF (LCCA(I).LE.0.0) THEN
            LCCA(I) = 0.7873 + 0.06 * LOG(TCW(I))
            IF (LCCA(I) .GT. 1.0) LCCA(I) = 1.0
          END IF
C
          TCW(I) = 0.0
C
        END IF
      END DO ! I loop over NPNTS
C
      RETURN
      END
