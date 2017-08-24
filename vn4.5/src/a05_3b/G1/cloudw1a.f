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
CLL  SUBROUTINE CLOUD_W------------------------------------------------
CLL
CLL  PURPOSE : CLOUD MICROPHYSICS ROUTINE
CLL
CLL            CALCULATES PRECIPITATION PRODUCED IN LIFTING PARCEL
CLL            FROM LAYER K TO K+1
CLL
CLL            CALL CON_RAD TO CALCULATE PARAMETERS FOR RADIATION
CLL            CALCULATION
CLL
CLL  SUITABLE FOR SINGLE COLUMN MODEL USE
CLL
CLL  CODE REWORKED FOR CRAY Y-MP BY D.GREGORY AUTUMN/WINTER 1989/90
CLL
CLL  MODEL            MODIFICATION HISTORY FROM MODEL VERSION 3.0:
CLL VERSION  DATE
CLL  3.2  8/07/93 : added convective cloud condensed water diagnostic
CLL               : P Inness
CLL  3.4  21/03/94  Add lowest conv.cloud diagnostics.  R.T.H.Barnes.
CLL  4.4  26/09/97  Logical L_CCW passed in to determine if convective
CLL                 precip is included in water path (no, if .T.)
CLL
CLL  PROGRAMMING STANDARDS : UNIFIED MODEL DOCUMENTATION PAPER NO. 4
CLL  VERSION NO. 1
CLL
CLL  PROJECT TASK : P27
CLL
CLL  DOCUMENTATION : UNIFIED MODEL DOCUMENTATION PAPER P27
CLL
CLLEND-----------------------------------------------------------------
C
C*L  ARGUMENTS---------------------------------------------------------
C
      SUBROUTINE CLOUD_W (K,NPNTS,XPKP1,PREKP1,XSQKP1,BLOWST,
     *                    FLXKP1,XPK,THEKP1,QEKP1,BWKP1,BLAND,
     *                    QSEKP1,BGMKP1,BTERM,CCA,ICCB,ICCT,TCW,DEPTH,
     *                    EKP14,EKP34,DELEXKP1,CCLWP,DELPKP1,CCW,
     *                    LCCA,LCBASE,LCTOP,LCCLWP,L_SHALLOW,
     *                    L_CCW
     &                   ,MPARWTR
     &                   ,UD_FACTOR
     &                    )
C
      IMPLICIT NONE
C
C----------------------------------------------------------------------
C MODEL CONSTANTS
C----------------------------------------------------------------------
C
      REAL CRITDSEA  !  CRITICAL DEPTH OF CLOUD FOR THE FORMATION OF
                     !  CONVECTIVE PRECIPITATION OVER SEA (M)
      PARAMETER (CRITDSEA = 1.5E3)
C
      REAL CRITDLND  !  CRITICAL DEPTH OF CLOUD FOR THE FORMATION OF
                     !  CONVECTIVE PRECIPITATION OVER LAND (M)
      PARAMETER (CRITDLND = 4.0E3)
C
      REAL CRITDICE  !  CRITICAL DEPTH OF A GLACIATED CLOUD FOR THE
                     !  FORMATION OF CONVECTIVE PRECIPITATION (M)
      PARAMETER (CRITDICE = 1.0E3)
C
      REAL MPARWTR  !  MINIMUM PARCEL WATER IN GRAMS PER KILOGRAM
                    !  BEFORE PRECIPITATION IS ALLOWED (KG/KG)
C
C*L------------------COMDECK C_G----------------------------------------
C G IS MEAN ACCEL DUE TO GRAVITY AT EARTH'S SURFACE
      REAL G

      PARAMETER(G=9.80665)
C*----------------------------------------------------------------------
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

C*L------------------COMDECK C_EPSLON-----------------------------------
C EPSILON IS RATIO OF MOLECULAR WEIGHTS OF WATER AND DRY AIR
      REAL EPSILON,C_VIRTUAL

      PARAMETER(EPSILON=0.62198,
     &          C_VIRTUAL=1./EPSILON-1.)
C*----------------------------------------------------------------------

C
C----------------------------------------------------------------------
C VECTOR LENGTH AND LOOP COUNTERS
C----------------------------------------------------------------------
C
      INTEGER NPNTS          ! IN VECTOR LENGTH
C
      INTEGER K              ! IN PRESENT MODEL LAYER
C
      INTEGER I              ! LOOP COUNTER
C
C
C----------------------------------------------------------------------
C VARIABLES WHICH ARE INPUT
C----------------------------------------------------------------------
C
      REAL THEKP1(NPNTS)     ! IN POTENTIAL TEMPERATURE OF CLOUD
                             !    ENVIRONMENT IN LAYER K+1 (K)
C
      REAL QEKP1(NPNTS)      ! IN MIXING RATIO OF CLOUD ENVIRONMENT
                             !    IN LAYER K+1 (KG/KG)
C
      REAL QSEKP1(NPNTS)     ! IN SATURATION MIXING RATIO OF CLOUD
                             !    ENVIRONMENT IN LAYER K+1 (KG/KG)
C
      REAL XPK(NPNTS)        ! IN PARCEL CLOUD WATER IN LAYER K (KG/KG)
C
      LOGICAL BWKP1(NPNTS)   ! IN MASK FOR WHETHER CONDENSATE IS
                             !    LIQUID IN LAYER K+1
C
      LOGICAL BGMKP1(NPNTS)  ! IN MASK FOR PARCELS WHICH ARE
                             !    SATURATED IN LAYER K+1
C
      LOGICAL BLAND(NPNTS)   ! IN LAND/SEA MASK
C
      LOGICAL BTERM(NPNTS)   ! IN MASK FOR PARCELS WHICH TERMINATE IN
                             !    LAYER K+1
C
      LOGICAL BLOWST(NPNTS)  ! IN MASK FOR THOSE POINTS AT WHICH
                             !    STABILITY IS LOW ENOUGH FOR
                             !    CONVECTION TO OCCUR
C
      LOGICAL L_SHALLOW(NPNTS) ! IN MASK FOR POINTS WHERE SHALLOW 
                               !    CONVECTION IS LIKELY
      LOGICAL L_CCW            ! IN SWITCH FOR CLOUD WATER CHANGES:
                               !    (PRECIP NOT INC. IN WATER PATH)
C
      REAL FLXKP1(NPNTS)     ! IN PARCEL MASSFLUX IN LAYER K+1 (PA/S)
C
      REAL XSQKP1(NPNTS)     ! IN EXCESS PARCEL MIXING RATIO IN
                             !    LAYER K+1 (KG/KG)
C
      REAL EKP14(NPNTS)      ! IN ENTRAINMENT RATE AT LEVEL K+1/4
                             !    MULTIPLIED BY APPROPRIATE LAYER
                             !    THICKNESS
C
      REAL EKP34(NPNTS)      ! IN ENTRAINEMNT RATE AT LEVEL K+3/4
                             !    MULTIPLIED BY APPROPRIATE LAYER
                             !    THICKNESS
C
      REAL DELEXKP1(NPNTS)   ! IN DIFFERENCE IN EXNER RATIO ACROSS
                             !    LAYER K+1 (PA)
C
      REAL DELPKP1(NPNTS)    ! IN PRESSURE DIFFERENCE ACROSS LAYER K+1
C
C
C----------------------------------------------------------------------
C VARIABLES WHICH ARE INPUT AND OUTPUT
C----------------------------------------------------------------------
C
      REAL TCW(NPNTS)        ! INOUT
                             ! IN  TOTAL CONDENSED WATER SUMMED UPTO
                             !     LAYER K (KG/M**2/S)
                             ! OUT TOTAL CONDENSED WATER SUMMED UPTO
                             !     LAYER K+1 (KG/M**2/S)
C
      REAL DEPTH(NPNTS)      ! INOUT
                             ! IN  DEPTH OF CONVECTIVE CLOUD TO
                             !     LAYER K (M)
                             ! OUT DEPTH OF CONVECTIVE CLOUD TO
                             !     LAYER K+1 (M)
C
      REAL CCLWP(NPNTS)      ! INOUT
                             ! IN  CONDENSED WATER PATH SUMMED UPTO
                             !     LAYER K (KG/M**2)
                             ! OUT CONDENSED WATER PATH SUMMED UPTO
                             !     LAYER K+1 (KG/M**2)
C
C
C----------------------------------------------------------------------
C VARIABLES WHICH ARE OUTPUT
C----------------------------------------------------------------------
C
      REAL PREKP1(NPNTS)     ! OUT PRECIPITATION FROM PARCEL AS IT
                             !     RISES FROM LAYER K TO K+1 (KG/M**2/S)
C
      REAL XPKP1(NPNTS)      ! OUT PARCEL CLOUD WATER IN LAYER K+1
                             !     (KG/KG)
C
      REAL CCA(NPNTS)        ! OUT CONVECTIVE CLOUD AMOUNT (%)
C
      INTEGER ICCB(NPNTS)    ! OUT CONVECTIVE CLOUD BASE LEVEL
C
      INTEGER ICCT(NPNTS)    ! OUT CONVECTIVE CLOUD TOP LEVEL
C
      REAL CCW(NPNTS)        ! OUT CONVECTIVE CLOUD LIQUID WATER
                             ! (G/KG) ON MODEL LEVELS
C
      REAL LCCA(NPNTS)       ! OUT LOWEST CONV.CLOUD AMOUNT (%)
C
      INTEGER LCBASE(NPNTS)  ! OUT LOWEST CONV.CLOUD BASE LEVEL
C
      INTEGER LCTOP(NPNTS)   ! OUT LOWEST CONV.CLOUD TOP LEVEL
C
      REAL LCCLWP(NPNTS)     ! OUT LOWEST CONV.CLOUD LIQ.WATER PATH
C
C
C----------------------------------------------------------------------
C VARIABLES WHICH ARE LOCALLY DEFINED
C----------------------------------------------------------------------
C
      REAL DCRIT             ! CRITICAL DEPTH AT WHICH PRECIPITATION
                             ! MAY FORM (M)
C
      REAL XMIN              ! AMOUNT OF CLOUD WATER RETAINED BY THE
                             ! PARCEL ON PRECIPITATION (KG/KG)
C
      REAL EPSS              ! (1.0+EKP14)*(1.0+EKP34)
C
      REAL CCW_UD(NPNTS)    ! Cloud water for radiation
C                           !
      REAL UD_FACTOR        ! Factor to multiply ccw by
C----------------------------------------------------------------------
C  EXTERNAL ROUTINES CALLED
C----------------------------------------------------------------------
C
      EXTERNAL CON_RAD
C
C*---------------------------------------------------------------------
CL
CL----------------------------------------------------------------------
CL  CALCULATE CLOUD WATER BEFORE PRECIPITATION
CL
CL UM DOCUMENTATION PAPER P27
CL SECTION (2B), EQUATION (13A)
CL----------------------------------------------------------------------
CL
      DO 10 I=1,NPNTS
       EPSS = (1.+EKP14(I))*(1.+EKP34(I))
       XPKP1(I) = (XPK(I)/EPSS) + XSQKP1(I)
   10 CONTINUE
CL
CL----------------------------------------------------------------------
CL STORE CONVECTIVE CLOUD LIQUID WATER BEFORE PRECIPITATION
CL----------------------------------------------------------------------
CL
      DO I=1,NPNTS
        CCW(I) = XPKP1(I)
      END DO
      IF (.NOT. L_CCW) THEN
CL
CL----------------------------------------------------------------------
CL CALCULATE CONVECTIVE CLOUD BASE, CONVECTIVE CLOUD TOP , TOTAL
CL CONDENSED WATER/ICE AND CONVECTIVE CLOUD AMOUNT
CL
CL SUBROUTINE CON_RAD
CL
CL UM DOCUMENTATION PAPER P27
CL SECTION (9)
CL----------------------------------------------------------------------
CL
      CALL CON_RAD(K,XPK,XPKP1,FLXKP1,BTERM,CCA,ICCB,ICCT,
     *       TCW,CCW,CCLWP,DELPKP1,LCCA,LCBASE,LCTOP,LCCLWP,NPNTS)      
      ENDIF
CL
CL----------------------------------------------------------------------
CL CALCULATE CLOUD DEPTH AND ASSIGN CRITICAL CLOUD DEPTHS
CL
CL UM DOCUMENTATION PAPER P27
CL SECTION (8), EQUATION (34), (35)
CL----------------------------------------------------------------------
CL
      DO 30 I=1,NPNTS
       IF ( BLOWST(I) ) DEPTH(I) = 0.
C
       IF ( BGMKP1(I) )
     *   DEPTH(I) = DEPTH(I) + ( CP * THEKP1(I) *
     *                             (1.0+C_VIRTUAL*QEKP1(I)) *
     *                                          DELEXKP1(I)/G )
C
       IF (.NOT.BWKP1(I)) THEN
          DCRIT = CRITDICE
       ELSE IF (BLAND(I)) THEN
          DCRIT = CRITDLND
       ELSE
          DCRIT = CRITDSEA
       ENDIF
CL
CL----------------------------------------------------------------------
CL CALCULATE PRECIPITATION FROM LAYER K+1 AND ADJUST CLOUD WATER
CL
CL UM DOCUMENTATION PAPER P27
CL SECTION (8), EQUATION (36)
CL----------------------------------------------------------------------
CL
       XMIN = MIN (MPARWTR , 0.5*QSEKP1(I))
       IF (     ( (DEPTH(I) .GT. DCRIT) .OR. 
     *            ((.NOT. L_SHALLOW(I)) .AND. (L_CCW)) )
     *     .AND. (XPKP1(I) .GT. XMIN)) THEN
          PREKP1(I) = (XPKP1(I) - XMIN) * FLXKP1(I) / G
          XPKP1(I) = XMIN
          CCW_UD(I)=XPKP1(I)*UD_FACTOR
       ELSE
          PREKP1 (I) = 0.
          CCW_UD(I)=XPKP1(I)
       ENDIF
   30  CONTINUE
      IF (L_CCW) THEN
CL
CL----------------------------------------------------------------------
CL CALCULATE CONVECTIVE CLOUD BASE, CONVECTIVE CLOUD TOP , TOTAL
CL CONDENSED WATER/ICE AND CONVECTIVE CLOUD AMOUNT
CL
CL SUBROUTINE CON_RAD - MOVED TO AFTER RAIN OUT HAS OCCURRED IF L_CCW
CL IS TRUE (SET IN UMUI).
CL UM DOCUMENTATION PAPER P27
CL SECTION (9)
CL----------------------------------------------------------------------
CL
        CALL CON_RAD(K,XPK,CCW_UD,FLXKP1,BTERM,CCA,ICCB,ICCT,
     *       TCW,CCW,CCLWP,DELPKP1,LCCA,LCBASE,LCTOP,LCCLWP,NPNTS)
CL
CL----------------------------------------------------------------------
CL STORE CONVECTIVE CLOUD LIQUID WATER AFTER PRECIPITATION
CL----------------------------------------------------------------------
CL
        DO I=1,NPNTS
          CCW(I) = CCW_UD(I)
        END DO
      ENDIF
C
      RETURN
      END
