C ******************************COPYRIGHT******************************
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
C
CLL  SUBROUTINE CONVEC2------------------------------------------------
CLL
CLL  PURPOSE : COMPLETES LIFTING OF THE PARCEL FROM LAYER K TO K+1
CLL
CLL            CALL SUBROUTINE PARCEL AND ENVIRON
CLL
CLL            SUBROUTINE PARCEL CALCULATES AN INITIAL MASS FLUX,
CLL            CARRIES OUT THE DETRAINMENT CALCULATION, TESTS
CLL            TO SEE IF CONVECTION IS TERMINATING AND CALCULATES THE
CLL            PRECIPITATION RATE FROM LAYER K+1
CLL
CLL            SUBROUTINE ENVIRON CALCULATES THE EFFECT OF CONVECTION
CLL            UPON THE LARGE-SCALE ATMOSPHERE
CLL
CLL  SUITABLE FOR SINGLE COLUMN MODEL USE
CLL
CLL  MODEL            MODIFICATION HISTORY:
CLL VERSION  DATE
CLL   4.2   1/11/96   New deck version based on CONV22A with HADCM2
CLL                   specific modifications: R Jones
!LL   4.5   19/5/98   Correction of CAPE diagnostic.   Julie Gregory
!LL  4.5   23/02/98  Call comdecks. D. Robinson
CLL
CLL  PROGRAMMING STANDARDS : UNIFIED MODEL DOCUMENTATION PAPER NO. 3
CLL
CLL  LOGICAL COMPONENT NUMBER: P27
CLL
CLL  SYSTEM TASK : P27
CLL
CLL  DOCUMENTATION : UNIFIED MODEL DOCUMENTATION PAPER 27
CLL
CLLEND-----------------------------------------------------------------
C
C*L  ARGUMENTS---------------------------------------------------------
C
      SUBROUTINE CONVEC2 (NPNTS,NLEV,K,THEK,THEKP1,QEK,QEKP1,QSEKP1,
     *                   DQSKP1,PSTAR,THPK,QPK,THPKP1,QPKP1,XSQKP1,
     *                   RBUOY,QSEK,DQSK,THPI,QPI,XPK,FLXK,BWKP1,BGMKP1,
     *                   BGMK,BLOWST,BLAND,BTERM,DEPTH,PREKP1,DTHEK,
     *                   DQEK,DTHEKP1,DQEKP1,BINIT,CCA,ICCB,ICCT,
     *                   TCW,EKP14,EKP34,AMDETK,PK,PKP1,
     *                   EXK,EXKP1,DELEXKP1,DELPK,DELPKP1,
     *                   CCLWP,CCW,LCCA,LCBASE,LCTOP,LCCLWP,
     *                   CAPE)
C
      IMPLICIT NONE
C
C----------------------------------------------------------------------
C  MODEL CONSTANTS
C----------------------------------------------------------------------
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
C
C----------------------------------------------------------------------
C VECTOR LENGTHS AND LOOP COUNTERS
C----------------------------------------------------------------------
C
      INTEGER NPNTS          ! IN VECTOR LENGTH
C
      INTEGER NLEV           ! IN NUMBER OF MODEL LAYERS
C
      INTEGER I              ! LOOP COUNTER
C
      INTEGER K              ! PRESENT MODEL LAYER
C
C
C----------------------------------------------------------------------
C VARIABLES WHICH ARE INPUT
C----------------------------------------------------------------------
C
      REAL THEK(NPNTS)       ! IN POTENTIAL TEMPERATURE OF CLOUD
                             !    ENVIRONMENT IN LAYER K (K)
C
      REAL THEKP1(NPNTS)     ! IN POTENTIAL TEMPERATURE OF CLOUD
                             !    ENVIRONMENT IN LAYER K+1 (K)
C
      REAL QEK(NPNTS)        ! IN MIXING RATIO OF CLOUD
                             !    ENVIRONMENT IN LAYER K (KG/KG)
C
      REAL QEKP1(NPNTS)      ! IN MIXING RATIO OF CLOUD
                             !    ENVIRONMENT IN LAYER K+1 (KG/KG)
C
      REAL QSEKP1(NPNTS)     ! IN SATURATION MIXING RATIO OF CLOUD
                             !    ENVIRONMENT IN LAYER K+1 (KG/KG)
C
      REAL DQSKP1(NPNTS)     ! IN GRADIENT OF SATURATION MIXING RATIO
                             !    WITH POTENTIAL TEMPERATURE FOR THE
                             !    CLOUD ENVIRONMENT IN LAYER K+1
                             !    (KG/KG)
C
      REAL PSTAR(NPNTS)      ! IN SURFACE PRESSURE (PA)
C
      REAL THPKP1(NPNTS)     ! IN PARCEL POTENTIAL TEMPERATURE
                             !    IN LAYER K (K)
C
      REAL QPKP1(NPNTS)      ! IN PARCEL MIXING RATIO IN LAYER K+1
                             !    (KG/KG)
C
      REAL XSQKP1(NPNTS)     ! IN EXCESS WATER IN PARCEL AFTER LIFTING
                             !    LAYER K TO K+1 (KG/KG)
C
      REAL RBUOY(NPNTS)      ! IN PARCEL BUOYANCY IN LAYER K+1 (K)
C
      REAL QSEK(NPNTS)       ! IN SATURATION MIXING RATIO OF CLOUD
                             !    ENVIRONMENT IN LAYER K (KG/KG)
C
      REAL DQSK(NPNTS)       ! IN GRADIENT OF SATURATION MIXING RATIO
                             !    WITH POTENTIAL TEMPERATURE FOR THE
                             !    CLOUD ENVIRONMENT OF LAYER K
                             !    (KG/KG/K)
C
      REAL THPI(NPNTS)       ! IN INITIAL PARCEL POTENTIAL TEMPERATURE
                             !    (K)
C
      REAL QPI(NPNTS)        ! IN INITIAL PARCEL MIXING RATIO
                             !    (KG/KG)
C
      LOGICAL BWKP1(NPNTS)   ! IN MASK FOR WHETHER CONDENSATE IS
                             !    LIQUID IN LAYER K+1
C
      LOGICAL BGMKP1(NPNTS)  ! IN MASK FOR PARCELS WHICH ARE
                             !    SATURATED IN LAYER K+1
C
      LOGICAL BLAND(NPNTS)   ! IN LAND/SEA MASK
C
      LOGICAL BINIT(NPNTS)   ! IN MASK FOR THOSE POINTS AT WHICH
                             !    CONVECTION IS OCCURING
C
      LOGICAL BLOWST(NPNTS)  ! IN MASK FOR THOSE POINTS AT WHICH
                             !    STABILITY IS LOW ENOUGH FOR
                             !    CONVECTION TO OCCUR
C
      REAL EKP14(NPNTS)      ! IN ENTRAINMENT COEFFICIENT AT LEVEL
                             !    K+1/4 MULTIPLIED BY APPROPRIATE
                             !    LAYER THICKNESS
C
      REAL EKP34(NPNTS)      ! IN ENTRAINMENT COEFFICIENT AT LEVEL
                             !    K+1/4 MULTIPLIED BY APPROPRIATE
                             !    LAYER THICKNESS
C
      REAL AMDETK(NPNTS)     ! IN MIXING DETRAINMENT COEFFICIENT
                             !    AT LEVEL K MULTIPLIED BY APPROPRIATE
                             !    LAYER THICKNESS
C
      REAL DELPKP12(NPNTS)   ! IN PRESSURE DIFFERENCE BETWEEN
                             !    MID-POINTS OF LAYERS K AND K+1
                             !    (PA)
C
      REAL PK(NPNTS)         ! IN PRESSURE AT MID-POINT OF LAYER K
                             !    (PA)
C
      REAL PKP1(NPNTS)       ! IN PRESSURE AT MID-POINT OF LAYER K+1
                             !    (PA)
C
      REAL EXK(NPNTS)        ! IN EXNER RATIO AT MID-POINT OF LAYER K
C
      REAL EXKP1(NPNTS)      ! IN EXNER RATIO AT MID-POINT OF LAYER K+1
C
      REAL DELEXKP1(NPNTS)   ! IN DIFFERENCE IN EXNER RATIO BETWEEN
                             !    MID-POINTS OF LAYERS K AND K+1
C
      REAL DELPK(NPNTS)      ! IN DIFFERENCE IN PRESSURE ACROSS LAYER K
                             !    (PA)
C
      REAL DELPKP1(NPNTS)    ! IN DIFFERENCE IN PRESSURE ACROSS
                             !    LAYER K+1 (PA)
C
C
C----------------------------------------------------------------------
C VARIABLES WHICH ARE INPUT BUT WHICH ARE ALSO UPDATED IN THIS ROUTINE
C----------------------------------------------------------------------
C
      REAL THPK(NPNTS)       ! INOUT
                             ! IN  PARCEL POTENTIAL TEMPERATURE
                             !     IN LAYER K (K)
                             ! OUT PARCEL POTENTIAL TEMPERATURE
                             !     IN LAYER K+1 (K)
C
      REAL QPK(NPNTS)        ! INOUT
                             ! IN  PARCEL MIXING RATIO IN LAYER K
                             !     (KG/KG)
                             ! OUT PARCEL MIXING RATIO IN LAYER K+1
                             !     (KG/KG)
C
      REAL XPK(NPNTS)        ! INOUT
                             ! IN  PARCEL CLOUD WATER IN LAYER K
                             !     (KG/KG)
                             ! OUT PARCEL CLOUD WATER IN LAYER K+1
                             !     (KG/KG)
C
      REAL FLXK(NPNTS)       ! INOUT
                             ! IN  PARCEL MASSFLUX IN LAYER K (PA/S)
                             ! OUT PARCEL MASSFLUX IN LAYER K+1 (PA/S)
C
      LOGICAL BGMK(NPNTS)    ! INOUT
                             ! IN  MASK FOR PARCELS WHICH ARE
                             !     SATURATED IN LAYER K
                             ! OUT MASK FOR PARCELS WHICH ARE
                             !     SATURATED IN LAYER K+1
C
      REAL DTHEK(NPNTS)      ! INOUT
                             ! IN  INCREMENT TO MODEL POTENTIAL
                             !     TEMPERATURE IN LAYER K DUE TO
                             !     CONVECTION (MAY BE NONE ZERO
                             !     DUE TO A PREVIOUS SPLIT FINAL
                             !     DETRAINEMNT CALCULATION) (K/S)
                             ! OUT UPDATED INCREMENT TO MODEL POTENTIAL
                             !     TEMPERATURE IN LAYER K DUE TO
                             !     CONVECTION (K/S)
C
      REAL DQEK(NPNTS)       ! INOUT
                             ! IN  INCREMENT TO MODEL MIXING RATIO
                             !     IN LAYER K DUE TO CONVECTION
                             !     (MAY BE NONE ZERO DUE TO A
                             !     PREVIOUS SPLIT FINAL DETRAINEMNT
                             !     CALCULATION) (KG/KG/S)
                             ! OUT UPDATED INCREMENT TO MODEL MIXING
                             !     RATIO IN LAYER K DUE TO CONVECTION
                             !     (KG/KG/S)
C
      REAL TCW(NPNTS)        ! INOUT
                             ! IN  TOTAL CONDENSED WATER SUMMED TO
                             !     LAYER K (KG/M**2/S)
                             ! OUT UPDATED TOTAL CONDENSED WATER
                             !     SUMMED TO LAYER K+1 (KG/M**2/S)
C
      REAL DEPTH(NPNTS)      ! INOUT
                             ! IN  DEPTH OF CONVECTIVE CLOUD TO
                             !     LAYER K (M)
                             ! OUT UPDATED DEPTH OF CONVECTIVE CLOUD
                             !     TO LAYER K+1 (M)
C
      REAL CCLWP(NPNTS)      ! INOUT
                             ! IN  CONDENSED WATER PATH SUMMED TO
                             !     LAYER K (KG/M**2)
                             ! OUT UPDATED CONDENSED WATER PATH
                             !     SUMMED TO LAYER K+1 (KG/M**2)
C
      REAL CAPE(NPNTS)       ! IN  CONVECTIVE AVAILABLE POTENTIAL
                             !     ENERGY UP TO THE CURRENT CONVECTING
                             !     LAYER
                             ! OUT CONVECTIVE AVAILABLE POTENTIAL
                             !     ENERGY INCLUDING ADDITION DUE TO
                             !     THE CAPE WITHIN THE CURRENT LAYER
C
C
C
C----------------------------------------------------------------------
C VARIABLES WHICH ARE OUTPUT
C----------------------------------------------------------------------
C
      LOGICAL BTERM(NPNTS)   ! OUT MASK FOR PARCELS WHICH TERMINATE IN
                             !     LAYER K+1
C
      REAL PREKP1(NPNTS)     ! OUT PRECIPITATION FROM PARCEL AS IT
                             !     RISES FROM LAYER K TO K+1 (KG/M**2/S)
C
      REAL DTHEKP1(NPNTS)    ! OUT INCREMENT TO MODEL POTENTIAL
                             !     TEMPERATURE IN LAYER K+1 DUE TO
                             !     CONVECTION (K/S)
C
      REAL DQEKP1(NPNTS)     ! OUT INCREMENT TO MODEL MIXING RATIO
                             !     IN LAYER K+1 DUE TO CONVECTION
                             !     (KG/KG/S)
C
      REAL CCA(NPNTS)        ! OUT CONVECTIVE CLOUD AMOUNT (%)
C
      INTEGER ICCB(NPNTS)     ! OUT CONVECTIVE CLOUD BASE LEVEL
C
      INTEGER ICCT(NPNTS)     ! OUT CONVECTIVE CLOUD TOP LEVEL
C
      REAL CCW(NPNTS)        ! OUT CONVECTIVE CLOUD WATER(G/KG) ON
                             ! MODEL LEVELS
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
C VARIABLES WHICH ARE DEFINED LOCALLY
C
C ON THE IBM ARRAYS ARE ALLOCATED USING A PARAMETER STATEMENT
C
C ON THE CRAY ARRAYS ARE DYNAMICALLY ALLOCATED
C----------------------------------------------------------------------
C
      REAL THRK(NPNTS)       ! PARCEL DETRAINMENT POTENTIAL
                             ! TEMPERATURE IN LAYER K (K)
C
      REAL QRK(NPNTS)        ! PARCEL DETRAINMENT MIXING RATIO
                             ! IN LAYER K (KG/KG)
C
      REAL XPKP1(NPNTS)      ! PARCEL CLOUD WATER IN LAYER K+1 (KG/KG)
C
      REAL FLXKP1(NPNTS )    ! PARCEL MASSFLUX IN LAYER K+1 (PA/S)
C
      REAL DELTAK(NPNTS)     ! PARCEL FORCED DETRAINMENT RATE
                             ! IN LAYER K MULTIPLIED BY APPROPRIATE
                             ! LAYER THICKNESS
C
      REAL THVP,THVE,RHO     ! VIRTUAL TEMPERATURE OF PARCEL, VIRTUAL
                             ! TEMPERATURE OF ENVIRONMENT AND DENSITY,
                             ! USED IN CALCULATION OF CAPE.
C
      LOGICAL L_SHALLOW         ! \  DUMMY Variables used in parcel in
     &,       L_CCW             !  | versions 3A and 3B. Added here to
      PARAMETER (L_CCW=.FALSE.) ! /  keep argument lists consistent.
!
C----------------------------------------------------------------------
C EXTERNAL ROUTINES CALLED
C----------------------------------------------------------------------
C
      EXTERNAL PARCEL,ENVIRON
C
C*---------------------------------------------------------------------
CL
CL----------------------------------------------------------------------
CL COMPLETE LIFTING PARCELS TO LAYER K+1
CL
CL SUBROUTINE PARCEL
CL
CL UM DOCUMENTATION PAPER P27
CL SECTIONS (5),(6),(7),(8),(9)
CL----------------------------------------------------------------------
CL
       CALL PARCEL (K,NPNTS,NLEV,PSTAR,THEKP1,THEK,QEKP1,QEK,
     *              QSEK,QSEKP1,DQSK,DQSKP1,BLAND,BWKP1,
     *              DELTAK,FLXK,THPK,QPK,THRK,QRK,
     *              BTERM,THPKP1,QPKP1,PREKP1,XPK,XPKP1,FLXKP1,
     *              XSQKP1,THPI,QPI,BGMK,BGMKP1,BLOWST,RBUOY,
     *              CCA,ICCB,ICCT,TCW,DEPTH,
     *              EKP14,EKP34,AMDETK,DELPKP1,PK,PKP1,
     *              EXK,EXKP1,DELEXKP1,CCLWP,CCW,
     &              LCCA,LCBASE,LCTOP,LCCLWP,L_SHALLOW,L_CCW)
CL
CL----------------------------------------------------------------------
CL CALCULATE THE EFFECT ON THE ENVIRONMENT (EXCEPT FOR THE
CL THE EVAPORATION OF PRECIPITATION AND CHANGE OF PHASE)
CL
CL SUBROUTINE ENVIRON
CL
CL UM DOCUMENTATION PAPER P27
CL SECTION (10)
CL----------------------------------------------------------------------
CL
       CALL ENVIRON (NPNTS,DTHEK,DQEK,DTHEKP1,DQEKP1,
     *               THEK,QEK,DELTAK,FLXK,
     *               THPK,QPK,THRK,QRK,THEKP1,QEKP1,
     *               BTERM,THPKP1,QPKP1,XPK,XPKP1,BWKP1,FLXKP1,
     *               BLOWST,EKP14,EXK,EXKP1,DELPK,DELPKP1,
     *               AMDETK)
C
      DO 10 I=1,NPNTS
C
C-----------------------------------------------------------------------
C RESET BINIT WHERE CONVECTION HAS TERMINATED
C-----------------------------------------------------------------------
C
        BINIT(I) = .NOT.BTERM(I)
   10 CONTINUE
C
CL---------------------------------------------------------------------
CL CALCULATE CONTRIBUTION TO CAPE AND RATE OF CHANGE OF CAPE DUE TO
CL THE UPDRAUGHT
CL---------------------------------------------------------------------
C
      DO I=1,NPNTS
        THVP=THPK(I)*(1.0+C_VIRTUAL*QPK(I))
        THVE=THEK(I)*(1.0+C_VIRTUAL*QEK(I))
        RHO=PK(I)/(R*THEK(I)*EXK(I))
C
        CAPE(I)=CAPE(I)+(THVP-THVE)*DELPK(I)/(RHO*THVE)                 
C
       IF(BTERM(I))THEN
        THVP=THPKP1(I)*(1.0+C_VIRTUAL*QPKP1(I))
        THVE=THEKP1(I)*(1.0+C_VIRTUAL*QEKP1(I))
        RHO=PKP1(I)/(R*THEKP1(I)*EXKP1(I))
C
        CAPE(I)=CAPE(I)+(THVP-THVE)*DELPKP1(I)/(RHO*THVE)               
C
       END IF
C
      END DO
C
CL
CL---------------------------------------------------------------------
CL SWAP PARCEL VALUES READY FOR THE NEXT PART OF ASCENT
CL FROM LAYER K+1 TO K+2
CL---------------------------------------------------------------------
CL
      IF (K .LT. NLEV-1) THEN
        DO 30 I=1,NPNTS
          IF (BINIT(I)) THEN
            THPK(I) = THPKP1(I)
            QPK(I) = QPKP1(I)
            XPK(I) = XPKP1(I)
            FLXK(I) = FLXKP1(I)
            BGMK(I) = BGMKP1(I)
          END IF
 30     CONTINUE
      END IF
C
      RETURN
      END
