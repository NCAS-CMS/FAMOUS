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
CLL  SUBROUTINE ENVIRON------------------------------------------------
CLL
CLL  PURPOSE : CALCULATE THE EFFECT OF CONVECTION UPON THE
CLL            LARGE-SCALE ATMOSPHERE
CLL
CLL  SUITABLE FOR SINGLE COLUMN MODEL USE
CLL
CLL  CODE REWORKED FOR CRAY Y-MP BY D.GREGORY AUTUMN/WINTER 1989/90
CLL
CLL  MODEL            MODIFICATION HISTORY FROM MODEL VERSION 3.0:
CLL VERSION  DATE
CLL
CLL  PROGRAMMING STANDARDS : UNIFIED MODEL DOCUMENTATION PAPER NO. 4
CLL  VERSION NO. 1
CLL
CLL  LOGICAL COMPONENTS COVERED: P27
CLL
CLL  SYSTEM TASK :
CLL
CLL  DOCUMENTATION : UNIFIED MODEL DOCUMENTATION PAPER P27
CLL
CLLEND-----------------------------------------------------------------
C
C*L  ARGUMENTS---------------------------------------------------------
C
      SUBROUTINE ENVIRON (NPNTS,DTHEK,DQEK,DTHEKP1,DQEKP1,
     *                    THEK,QEK,DELTAK,FLXK,THPK,QPK,
     *                    THRK,QRK,THEKP1,QEKP1,BTERM,THPKP1,
     *                    QPKP1,XPK,XPKP1,BWKP1,FLXKP1,BLOWST,
     *                    EKP14,EXK,EXKP1,DELPK,DELPKP1,AMDETK)
C
      IMPLICIT NONE
C
C-----------------------------------------------------------------------
C MODEL CONSTANTS
C-----------------------------------------------------------------------
C
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
      REAL THPIXS, QPIXS ! INITIAL EXCESS POTENTIAL TEMPERATURE (K)
                         ! AND MIXING RATIO (KG/KG)
       PARAMETER ( THPIXS = 0.2, QPIXS = 0.0 )
C
C
C-----------------------------------------------------------------------
C VECTOR LENGTHS AND LOOP COUNTERS
C-----------------------------------------------------------------------
C
      INTEGER NPNTS          ! IN VECTOR LENGTH
C
      INTEGER I              ! LOOP COUNTER
C
C
C-----------------------------------------------------------------------
C VARIABLES THAT ARE INPUT
C-----------------------------------------------------------------------
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
      REAL THPK(NPNTS)       ! IN PARCEL POTENTIAL TEMPERATURE IN
                             !    LAYER K (K)
C
      REAL QPK(NPNTS)        ! IN PARCEL MIXING RATIO IN LAYER K (KG/KG)
C
      REAL THPKP1(NPNTS)     ! IN PARCEL POTENTIAL TEMPERATURE IN
                             !    LAYER K+1 (K)
C
      REAL QPKP1(NPNTS)      ! IN PARCEL MIXING RATIO IN LAYER K+1
                             !    (KG/KG)
C
      REAL XPK(NPNTS)        ! IN PARCEL CLOUD WATER IN LAYER K (KG/KG)
C
      REAL FLXK(NPNTS)       ! IN PARCEL MASSFLUX IN LAYER K (PA/S)
C
      LOGICAL BWKP1(NPNTS)   ! IN MASK FOR WHETHER CONDENSATE IS
                             !    LIQUID IN LAYER K+1
C
      LOGICAL BTERM(NPNTS)   ! IN MASK FOR PARCELS WHICH TERMINATE IN
                             !    LAYER K+1
C
      LOGICAL BLOWST(NPNTS)  ! IN MASK FOR THOSE POINTS AT WHICH
                             !    STABILITY IS LOW ENOUGH FOR
                             !    CONVECTION TO OCCUR
C
      REAL THRK(NPNTS)       ! IN PARCEL DETRAINMENT POTENTIAL
                             !    TEMPERATURE IN LAYER K (K)
C
      REAL QRK(NPNTS)        ! IN PARCEL DETRAINMENT MIXING RATIO
                             !    IN LAYER K (KG/KG)
C
      REAL XPKP1(NPNTS)      ! IN PARCEL CLOUD WATER IN LAYER K+1
                             !    (KG/KG)
C
      REAL FLXKP1(NPNTS)     ! IN PARCEL MASSFLUX IN LAYER K+1 (PA/S)
C
      REAL DELTAK(NPNTS)     ! IN PARCEL FORCED DETRAINMENT RATE
                             !    IN LAYER K MULTIPLIED BY APPROPRIATE
                             !    LAYER THICKNESS
C
      REAL EKP14(NPNTS)      ! IN ENTRAINMENT RATE FOR LEVEL K+1/4
                             !    MULTIPLIED BY APPROPRIATE LAYER
                             !    THICKNESS
C
      REAL EXK(NPNTS)        ! IN EXNER RATIO FOR MID-POINT OF LAYER K
C
      REAL EXKP1(NPNTS)      ! IN EXNER RATIO FOR MID-POINT OF
                             !    LAYER K+1
C
      REAL DELPK(NPNTS)      ! IN PRESSURE DIFFERENCE ACROSS LAYER K
                             !    (PA)
C
      REAL DELPKP1(NPNTS)    ! IN PRESSURE DIFFERENCE ACROSS LAYER K+1
                             !    (PA)
C
      REAL AMDETK(NPNTS)     ! IN MIXING DETRIANMENT AT LEVEL K
                             !    MULTIPLIED BY APPROPRIATE LAYER
                             !    THICKNESS
C
C
C-----------------------------------------------------------------------
C VARIABLES THAT ARE INPUT AND OUTPUT
C-----------------------------------------------------------------------
C
      REAL DTHEK(NPNTS)      ! INOUT
                             ! IN  INCREMENT TO MODEL POTENTIAL
                             !     TEMPERATURE IN LAYER K DUE TO
                             !     CONVECTION (MAY BE NONE ZERO
                             !     DUE TO A PREVIOUS SPLIT FINAL
                             !     DETRAINMENT CALCULATION) (K/S)
                             ! OUT UPDATED INCREMENT TO MODEL POTENTIAL
                             !     TEMPERATURE IN LAYER K DUE TO
                             !     CONVECTION (K/S)
C
      REAL DQEK(NPNTS)       ! INOUT
                             ! IN  INCREMENT TO MODEL MIXING RATIO
                             !     IN LAYER K DUE TO CONVECTION
                             !     (MAY BE NONE ZERO
                             !     DUE TO A PREVIOUS SPLIT FINAL
                             !     DETRAINMENT CALCULATION) (KG/KG/S)
                             ! OUT UPDATED INCREMENT TO MODEL MIXING
                             !     RATIO IN LAYER K DUE TO
                             !     CONVECTION (KG/KG/S)
C
C
C-----------------------------------------------------------------------
C VARIABLES THAT ARE OUTPUT
C-----------------------------------------------------------------------
C
      REAL DTHEKP1(NPNTS)    ! OUT INCREMENT TO MODEL POTENTIAL
                             !     TEMPERATURE IN LAYER K+1 DUE TO
                             !     CONVECTION (K/S)
C
      REAL DQEKP1(NPNTS)     ! OUT INCREMENT TO MODEL MIXING RATIO
                             !     IN LAYER K+1 DUE TO CONVECTION
                             !     (KG/KG/S)
C
C
C-----------------------------------------------------------------------
C VARIABLES THAT ARE DEFINED LOCALLY
C-----------------------------------------------------------------------
C
      REAL EL                ! LATENT HEAT OF CONDENSATION OR
                             ! (CONDENSATION + FUSION) (J/KG)
C
      REAL TEMPRY            ! TEMPORARY ARRAY
C
C*---------------------------------------------------------------------
C
      DO 10 I=1,NPNTS
C
C-----------------------------------------------------------------------
C   CREATE A VECTOR OF LATENT HEATS
C-----------------------------------------------------------------------
C
       IF (BWKP1(I)) THEN
          EL = LC
       ELSE
          EL = LC + LF
       ENDIF
C
C----------------------------------------------------------------------
C CALCULATE PARCEL MASSFLUX DIVIDED BY THE THICKNESS OF LAYER K
C THIS VALUE IS USED IN SEVERAL PLACES IN THE SUBROUTINE
C----------------------------------------------------------------------
C
       TEMPRY = FLXK(I)/DELPK(I)
C
       IF (BLOWST(I)) THEN
CL
CL----------------------------------------------------------------------
CL AT THE LOWEST CONVECTIVE LAYER, THE PARCEL MASS FLUX IS A FLUX FROM
CL THE ENVIRONMENT. IE. THE INITIAL MASS FLUX IS ENTRAINED WITH EXCESS
CL POTENTIAL TEMPERATURE AND MIXING RATIO TPIXS, QPIXS
CL
CL UM DOCUMENTATIO PAPER P27
CL SECTION (10), EQUATION (39)
CL----------------------------------------------------------------------
CL
         DTHEK(I) = DTHEK(I) - TEMPRY*THPIXS
         DQEK(I) = DQEK(I) - TEMPRY*QPIXS
       ENDIF
CL
CL---------------------------------------------------------------------
CL EFFECT OF CONVECTION UPON POTENTIAL TEMPERATURE OF LAYER K
CL
CL UM DOCUMENTATION PAPER P27
CL SECTION (10), EQUATION (38A)
CL--------------------------------------------------------------------
CL
       DTHEK(I) = DTHEK(I) + TEMPRY * (
     *
     *           (1+EKP14(I)) * (1.0-DELTAK(I)) *        ! COMPENSATING
     *           (1-AMDETK(I)) * (THEKP1(I)-THEK(I))     ! SUBSIDENCE
     *         +
     *           DELTAK(I) * (1.0-AMDETK(I)) *           ! FORCED
     *           (THRK(I)-THEK(I)-                       ! DETRAINMENT
     *                    ((EL/CP)*XPK(I)/EXK(I)))
     *         +
     *           AMDETK(I) * (THPK(I)-THEK(I)-           ! MIXING
     *                    ((EL/CP)*XPK(I)/EXK(I)))       ! DETRAINMENT
     *         )
CL
CL---------------------------------------------------------------------
CL EFFECT OF CONVECTION UPON MIXING RATIO OF LAYER K
CL
CL UM DOCUMENTATION PAPER P27
CL SECTION (10), EQUATION (38B)
CL--------------------------------------------------------------------
CL
       DQEK(I) = DQEK(I) + TEMPRY * (
     *
     *           (1+EKP14(I)) * (1.0-DELTAK(I)) *        ! COMPENSATING
     *           (1-AMDETK(I)) * (QEKP1(I)-QEK(I))       ! SUBSIDENCE
     *         +
     *           DELTAK(I) * (1.0-AMDETK(I)) *           ! FORCED
     *           (QRK(I)-QEK(I)+XPK(I))                  ! DETRAINMENT
     *         +
     *           AMDETK(I) * (QPK(I)-QEK(I)+             ! MIXING
     *                                XPK(I))            ! DETRAINMENT
     *         )
CL
CL----------------------------------------------------------------------
CL TERMINAL DETRAINMENT AND SUBSIDENCE IN TERMINAL LAYER
CL
CL UM DOCUMENTATION PAPER P27
CL SECTION (10), EQUATION (40)
CL--------------------------------------------------------------------
CL
       IF ( BTERM(I) ) THEN
          TEMPRY = FLXKP1(I)/DELPKP1(I)
          DTHEKP1(I) = DTHEKP1(I) + TEMPRY*((THPKP1(I)-THEKP1(I))
     *                                   - EL*XPKP1(I)/(EXKP1(I)*CP))
          DQEKP1(I)  = DQEKP1(I) + TEMPRY*(QPKP1(I)-QEKP1(I)
     *                                                + XPKP1(I))
       ENDIF
  10   CONTINUE
C
      RETURN
      END
