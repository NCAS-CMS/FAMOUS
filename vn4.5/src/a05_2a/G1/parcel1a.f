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
CLL  SUBROUTINE PARCEL-------------------------------------------------
CLL
CLL  PURPOSE : COMPLETES LIFTING OF THE PARCEL FROM LAYER K TO K+1
CLL
CLL            CALL SUBROUTINE DETRAIN, TERM_CON, CLOUD_W
CLL
CLL            AN INITIAL MASS FLUX IS CALCULATED
CLL
CLL            SUBROUTINE DETRAIN CARRIES OUT THE FORCED DETRAINMENT
CLL            CALCULATION
CLL
CLL            SUBROUTINE TERM_CON TESTS FOR ANY CONVECTION WHICH IS
CLL            TERMINATING IN LAYER K+1
CLL
CLL            SUBROUTINE CLOUD_W CARRIES OUT THE CLOUD MICROPHYSICS
CLL            CALCULATION
CLL
CLL  SUITABLE FOR SINGLE COLUMN MODEL USE
CLL
CLL  CODE REWORKED FOR CRAY Y-MP BY D.GREGORY AUTUMN/WINTER 1989/90
CLL                                   THE RELEVANT POINTS IN CONVECT
CLL
CLL  MODEL            MODIFICATION HISTORY FROM MODEL VERSION 3.0:
CLL VERSION  DATE
CLL  3.2  8/07/93 : added convective cloud condensed water diagnostic
CLL               : P Inness
CLL  3.4  21/03/94  Add lowest conv.cloud diagnostics.  R.T.H.Barnes.
CLL   4.2    Oct. 96  T3E migration: *DEF CRAY removed 
CLL                                    S.J.Swarbrick
CLL  4.4  26/09/97  Logical L_CCW passed in to determine if precip is
CLL                 included in water path. MPARWTR passed down as an
CLL                 argument at versions 3A and 3B.       J.M.Gregory
CLL   4.5    Jul. 98  Kill the IBM specific lines (JCThil)
CLL
CLL  PROGRAMMING STANDARDS : UNIFIED MODEL DOCUMENTATION PAPER NO. 4
CLL  VERSION NO. 1
CLL
CLL  LOGICAL COMPONENTS COVERED: P27
CLL
CLL  DOCUMENTATION : UNIFIED MODEL DOCUMENTATION PAPER P27
CLL
CLLEND-----------------------------------------------------------------
C
C*L  ARGUMENTS---------------------------------------------------------
      SUBROUTINE PARCEL (K,NPNTS,NLEV,PSTAR,THEKP1,THEK,QEKP1,QEK,
     *                   QSEK,QSEKP1,DQSK,DQSKP1,BLAND,BWKP1,
     *                   DELTAK,FLXK,THPK,QPK,THRK,QRK,
     *                   BTERM,THPKP1,QPKP1,PREKP1,XPK,XPKP1,FLXKP1,
     *                   XSQKP1,THPI,QPI,BGMK,BGMKP1,BLOWST,RBUOY,
     *                   CCA,ICCB,ICCT,TCW,DEPTH,
     *                   EKP14,EKP34,AMDETK,DELPKP1,PK,PKP1,
     *                   EXK,EXKP1,DELEXKP1,CCLWP,CCW,
     *                LCCA,LCBASE,LCTOP,LCCLWP,L_SHALLOW,L_CCW
     &               ,UD_FACTOR
     &                )
C
      IMPLICIT NONE
C
C----------------------------------------------------------------------
C MODEL CONSTANTS
C----------------------------------------------------------------------
C
      REAL XSBMIN !  MINIMUM EXCESS BUOYANCY TO CONTINUE PARCEL ASCENT
                  !  (K)
      PARAMETER (XSBMIN = 0.2)
C
      REAL C,D  !  CONSTANTS USED TO DETERMINE THE INITIAL CONVECTIVE
                !  MASS FLUX FROM PARCEL BUOYANCY
                ! MI = 1.E-3 * (D + C*BUOYANCY/DELP)
       PARAMETER ( C = 5.17E-4, D = 0.0 )
C
C
C----------------------------------------------------------------------
C VECTOR LENGTHS AND LOOP COUNTERS
C----------------------------------------------------------------------
C
      INTEGER NPNTS          ! IN VECTOR LENGTH
C
      INTEGER NLEV           ! IN NUMBER OF MODEL LEVELS
C
      INTEGER NDET           ! COMPRESSED VECTOR LENGTH FOR
                             ! FORCED DETRAINMENT CALCULATION
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
                             !    (KG/KG/K)
C
      REAL PSTAR(NPNTS)      ! IN SURFACE PRESSURE (PA)
C
      REAL THPK(NPNTS)       ! IN PARCEL POTENTIAL TEMPERATURE
                             !    IN LAYER K (KG/KG)
C
      REAL QPK(NPNTS)        ! IN PARCEL MIXING RATIO IN LAYER K (KG/KG)
C
      REAL XSQKP1(NPNTS)     ! IN EXCESS PARCEL WATER AFER LIFTING FROM
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
      REAL QPI(NPNTS)        ! IN INITIAL PARCEL MIXING RATIO (KG/KG)
C
      REAL XPK(NPNTS)        ! IN PARCEL CLOUD WATER IN LAYER K (KG/KG)
C
      LOGICAL BWKP1(NPNTS)   ! IN MASK FOR WHETHER CONDENSATE IS
                             !    LIQUID IN LAYER K+1
C
      LOGICAL BGMK(NPNTS)    ! IN MASK FOR PARCELS WHICH ARE
                             !    SATURATED IN LAYER K
C
      LOGICAL BLAND(NPNTS)   ! IN LAND/SEA MASK
C
      LOGICAL BLOWST(NPNTS)  ! IN MASK FOR THOSE POINTS AT WHICH
                             !    STABILITY IS LOW ENOUGH FOR
                             !    CONVECTION TO OCCUR
C
      LOGICAL L_SHALLOW(NPNTS) ! IN MASK FOR POINTS WHERE CONVECTION
                               !    IS EXPECTED TO BE SHALLOW
      LOGICAL L_CCW            ! IN SWITCH FOR CLOUD WATER CHANGES:
                               !    (PRECIP NOT INC. IN WATER PATH)
C
      REAL EKP14(NPNTS)      ! IN ENTRAINMENT COEFFICIENT AT LEVEL
                             !    K+1/4 MULTIPLIED BY APPROPRIATE
                             !    LAYER THICKNESS
C
      REAL EKP34(NPNTS)      ! IN ENTRAINMENT COEFFICIENT AT LEVEL
                             !   K+3/4 MULTIPLIED BY APPROPRIATE
                             !   LAYER THICKNESS
C
      REAL AMDETK(NPNTS)     ! IN MIXING DETRAINMENT COEFFICIENT
                             !    AT LEVEL K MULTIPLIED BY
                             !    APPROPORIATE LAYER THICKNESS
C
      REAL DELPKP1(NPNTS)    ! IN PRESSURE DIFFERENCE ACROSS
                             !    LAYER K+1 (PA)
C
      REAL PK(NPNTS)         ! IN PRESSURE AT LEVEL K (PA)
C
      REAL PKP1(NPNTS)       ! IN PRESSURE AT LEVEL K+1 (PA)
C
      REAL EXK(NPNTS)        ! IN EXNER FUNCTION AT LEVEL K
C
      REAL EXKP1(NPNTS)      ! IN EXNER FUNCTION AT LEVEL K+1
C
      REAL DELEXKP1(NPNTS)   ! IN DIFFERENCE IN EXNER FUNCTION ACROSS
                             !    LAYER K+1
!
      REAL UD_FACTOR         ! IN Updraught factor. Set in user 
                             !    interface to reduce cloud water
                             !    seen by radiation. Used if L_CCW
                             !    is set to true.
C
C
C---------------------------------------------------------------------
C VARAIBLES WHICH ARE BOTH INPUT AND OUTPUT
C---------------------------------------------------------------------
C
      REAL THPKP1(NPNTS)     ! INOUT
                             ! IN  ESTIMATE OF PARCEL POTENTIAL
                             !     TEMPERATURE IN LAYER K+1 AFTER
                             !     ENTRAINMENT AND LATENT HEATING (K)
                             ! OUT FINAL PARCEL POTENTIAL TEMPERATURE
                             !     IN LAYER K+1 (AFTER FORCED
                             !     DETRAINEMNT) (K)
C
      REAL QPKP1(NPNTS)      ! INOUT
                             ! IN  ESTIMATE OF PARCEL MIXING RATIO
                             !     IN LAYER K+1 AFTER ENTRAINMENT AND
                             !     LATENT HEATING (KG/KG)
                             ! OUT FINAL PARCEL MIXING RATIO
                             !     IN LAYER K+1 (AFTER FORCED
                             !     DETRAINEMNT) (KG/KG)
C
      REAL FLXK(NPNTS)       ! INOUT
                             ! IN  PARCEL MASSFLUX IN LAYER K
                             !     (NON-ZERO IF CONVECTION IS NOT
                             !     INITIATED FROM LAYER K) (PA/S)
                             ! OUT PARCEL MASSFLUX IN LAYER K
                             !     (SET IF CONVECTION IS INITIATED
                             !     IN LAYER K) (PA/S)
C
      LOGICAL BGMKP1(NPNTS)  ! INOUT
                             ! IN  MASK FOR PARCELS WHICH ARE
                             !     SATURATED IN LAYER K+1
                             !     CALCULATED ON THE BASIS OF
                             !     INPUT PARCEL POTENTIAL TEMPERATURE
                             !     AND MIXING RATIO
                             ! OUT MASK FOR PARCELS WHICH ARE
                             !     SATURATED IN LAYER K+1 CALCULATED
                             !     FORM PARCEL TEMPERATURE AND
                             !     MIXING RATIO AFTER FORCED
                             !     DETARINMENT CALCULATION
C
      REAL TCW(NPNTS)        ! INOUT
                             ! IN  TOTAL CONDENSED WATER CONTENT
                             !     SUMMED UPTO LAYER K (KG/M**2/S)
                             ! OUT UPDATED TOTAL CONDENSED WATER
                             !     CONTENT SUMMED UPTO LAYER K+1
                             !     (KG/M**2/S)
C
      REAL DEPTH(NPNTS)      ! INOUT
                             ! IN  DEPTH OF CONVECTIVE CLOUD TO
                             !     LAYER K (M)
                             ! OUT UPDATED DEPTH OF CONVECTIVE
                             !     CLOUD TO LAYER K+1 (M)
C
      REAL CCLWP(NPNTS)      ! INOUT
                             ! IN  CONDENSED WATER PATH
                             !     SUMMED UPTO LAYER K (KG/M**2)
                             ! OUT UPDATED CONDENSED WATER PATH
                             !     SUMMED UPTO LAYER K+1 (KG/M**2)
C
C
C---------------------------------------------------------------------
C VARIABLES WHICH ARE OUTPUT
C---------------------------------------------------------------------
C
      LOGICAL BTERM(NPNTS)   ! OUT MASK FOR PARCELS WHICH TERMINATE IN
                             !     LAYER K+1
C
      REAL PREKP1(NPNTS)     ! OUT PRECIPITATION FROM PARCEL AS IT
                             !     RISES FROM LAYER K TO K+1 (KG/M**2/S)
C
      REAL THRK(NPNTS)       ! OUT PARCEL DETRAINMENT POTENTIAL
                             !     TEMPERATURE IN LAYER K (K)
C
      REAL QRK(NPNTS)        ! OUT PARCEL DETRAINMENT MIXING RATIO
                             !     IN LAYER K (KG/KG)
C
      REAL XPKP1(NPNTS)      ! OUT PARCEL CLOUD WATER IN LAYER K+1
                             !     (KG/KG)
C
      REAL FLXKP1(NPNTS)     ! OUT PARCEL MASSFLUX IN LAYER K+1 (PA/S)
C
      REAL DELTAK(NPNTS)     ! OUT PARCEL FORCED DETRAINMENT
                             !     COEFFICIENT IN LAYER K
                             !     MULTIPLIED BY APPROPRIATE
                             !     LAYER THICKNESS
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
C---------------------------------------------------------------------
C VARIABLES WHICH ARE DEFINED LOCALLY
C---------------------------------------------------------------------
C
      REAL THEK_C(NPNTS)     ! COMPRESSED POTENTIAL TEMPERATURE OF
                             ! CLOUD ENVIRONMENT IN LAYER K (K)
C
      REAL THEKP1_C(NPNTS)   ! COMPRESSED POTENTIAL TEMPERATURE OF
                             ! CLOUD ENVIRONMENT IN LAYER K+1 (K)
C
      REAL QEK_C(NPNTS)      ! COMPRESSED MIXING RATIO OF CLOUD
                             ! ENVIRONMENT IN LAYER K (KG/KG)
C
      REAL QEKP1_C(NPNTS)    ! COMPRESSED MIXING RATIO OF CLOUD
                             ! ENVIRONMENT IN LAYER K+1 (KG/KG)
C
      REAL QSEK_C(NPNTS)     ! COMPRESSED SATURATION MIXING RATIO OF
                             ! CLOUD ENVIRONMENT IN LAYER K (KG/KG)
C
      REAL DQSK_C(NPNTS)     ! COMPRESSED GRADIENT OF SATURATION MIXING
                             ! RATIO WITH POTENTIAL TEMPERATURE FOR THE
                             ! CLOUD ENVIRONMENT OF LAYER K (KG/KG/K)
C
      REAL QSEKP1_C(NPNTS)   ! COMPRESSED SATURATION MIXING RATIO OF
                             ! CLOUD ENVIRONMENT IN LAYER K+1 (KG/KG)
C
      REAL DQSKP1_C(NPNTS)   ! COMPRESSED GRADIENT OF SATURATION MIXING
                             ! RATIO WITH POTENTIAL TEMPERATURE FOR
                             ! THE CLOUD ENVIRONMENT IN LAYER K+1
                             ! (KG/KG/K)
C
      REAL THPK_C(NPNTS)     ! COMPRESSED PARCEL POTENTIAL
                             ! TEMPERATURE IN LAYER K (K)
C
      REAL QPK_C(NPNTS)      ! COMPRESSED PARCEL MIXING RATIO IN
                             ! LAYER K (KG/KG)
C
      REAL THPKP1_C(NPNTS)   ! COMPRESSED PARCEL POTENTIAL
                             ! TEMPERATURE IN LAYER K+1 (K)
C
      REAL QPKP1_C(NPNTS)    ! COMPRESSED PARCEL MIXING RATIO
                             ! IN LAYER K+1 (KG/KG)
C
      REAL XSQKP1_C(NPNTS)   ! EXCESS PARCEL WATER AFER LIFTING
                             ! FROM LAYER K TO K+1 (KG/KG)
C
      REAL THRK_C(NPNTS)     ! COMPRESSED PARCEL DETRAINMENT
                             ! POTENTIAL TEMPERATURE IN LAYER K (K)
C
      REAL QRK_C(NPNTS)      ! COMPRESSED PARCEL DETRAINMENT MIXING
                             ! RATIO IN LAYER K (KG/KG)
C
      REAL DELTAK_C(NPNTS)   ! COMPRESSED PARCEL FORCED DETRAINMENT
                             ! COEFFICIENT IN LAYER K
                             ! MULTIPLIED BY APPROPRIATE
                             ! LAYER THICKNESS
C
      REAL EKP14_C(NPNTS)    ! COMPRESSED IN ENTRAINMENT COEFFICIENT AT
                             ! LEVEL K+1/4 MULTIPLIED BY APPROPRIATE
                             ! LAYER THICKNESS
C
      REAL EKP34_C(NPNTS)    ! COMPRESSED ENTRAINMENT COEFFICIENT AT
                             ! LEVEL K+3/4 MULTIPLIED BY APPROPRIATE
                             ! LAYER THICKNESS
C
      REAL PK_C(NPNTS)       ! COMPRESSED PRESSURE AT LEVEL K (PA)
C
      REAL PKP1_C(NPNTS)     ! COMPRESSED PRESSURE AT LEVEL K+1 (PA)
C
      REAL EXK_C(NPNTS)      ! COMPRESSED EXNER FUNCTION AT LEVEL K
C
      REAL EXKP1_C(NPNTS)    ! COMPRESSED EXNER FUNCTION AT LEVEL K+1
C
      LOGICAL BWKP1_C(NPNTS) ! COMPRESSED MASK FOR WHETHER CONDENSATE
                             ! IS LIQUID IN LAYER K+1
C
      LOGICAL BGMK_C(NPNTS)  ! COMPRESSED MASK FOR PARCELS WHICH ARE
                             ! SATURATED IN LAYER K
C
      LOGICAL BGMKP1_C(NPNTS) ! COMPRESSED MASK FOR PARCELS
                              ! WHICH ARESATURATED IN LAYER K+1
C
      INTEGER INDEX1(NPNTS)  ! INDEX FOR COMPRESS AND EXPAND
C
      LOGICAL BDETK(NPNTS)   ! MASK FOR POINTS UNDERGOING
                             ! FORCED DETRAINMENT
C
C
C----------------------------------------------------------------------
C EXTERNAL ROUTINES CALLED
C----------------------------------------------------------------------
C
      EXTERNAL DETRAIN,TERM_CON,CLOUD_W
C
C*--------------------------------------------------------------------
C
C
      DO 5 I=1,NPNTS
CL
CL---------------------------------------------------------------------
CL CALCULATE MASK FOR THOSE POINTS UNDERGOING FORCED DETRAINMENT
CL
CL UM DOCUMENTATION PAPER P27
CL SECTION (6), EQUATION (23)
CL---------------------------------------------------------------------
CL
       BDETK(I) = RBUOY(I) .LT. XSBMIN
C
   5  CONTINUE
CL
CL----------------------------------------------------------------------
CL  COMPRESS ALL INPUT ARRAYS FOR THE FORCED DETRAINMENT CALCULATIONS
CL----------------------------------------------------------------------
CL
      NDET = 0
      DO 10 I=1,NPNTS
       IF (BDETK(I))THEN
         NDET = NDET + 1
         INDEX1(NDET) = I
       END IF
  10  CONTINUE
C
      IF (NDET .NE. 0) THEN
        DO 35 I=1,NDET
          THEK_C(I)  = THEK(INDEX1(I))
          QEK_C(I)   = QEK(INDEX1(I))
          THPK_C(I)  = THPK(INDEX1(I))
          QPK_C(I)   = QPK(INDEX1(I))
          QSEK_C(I)  = QSEK(INDEX1(I))
          DQSK_C(I)  = DQSK(INDEX1(I))
          THEKP1_C(I)= THEKP1(INDEX1(I))
          QEKP1_C(I) = QEKP1(INDEX1(I))
          THPKP1_C(I)= THPKP1(INDEX1(I))
          QPKP1_C(I) = QPKP1(INDEX1(I))
          QSEKP1_C(I)= QSEKP1(INDEX1(I))
          DQSKP1_C(I)= DQSKP1(INDEX1(I))
          XSQKP1_C(I)= XSQKP1(INDEX1(I))
          EKP14_C(I) = EKP14(INDEX1(I))
          EKP34_C(I) = EKP34(INDEX1(I))
          PK_C(I)    = PK(INDEX1(I))
          PKP1_C(I)  = PKP1(INDEX1(I))
          EXK_C(I)   = EXK(INDEX1(I))
          EXKP1_C(I) = EXKP1(INDEX1(I))
C
          BGMK_C(I)  = BGMK(INDEX1(I))
          BGMKP1_C(I)= BGMKP1(INDEX1(I))
          BWKP1_C(I) = BWKP1(INDEX1(I))
  35    CONTINUE
CL
CL-------------------------------------------------------------------
CL DETRAINMENT CALCULATION
CL
CL SUBROUTINE DETRAIN
CL
CL UM DOCUMENTATION PAPER P27
CL SECTION (6)
CL-------------------------------------------------------------------
CL
         CALL DETRAIN (NDET,THEK_C,QEK_C,THPK_C,QPK_C,
     *                 QSEK_C,DQSK_C,BGMK_C,THEKP1_C,
     *                 QEKP1_C,THPKP1_C,QPKP1_C,QSEKP1_C,
     *                 DQSKP1_C,BGMKP1_C,BWKP1_C,
     *                 XSQKP1_C,DELTAK_C,
     *                 THRK_C,QRK_C,EKP14_C,EKP34_C,
     *                 PK_C,PKP1_C,EXK_C,EXKP1_C)
C
C-----------------------------------------------------------------------
C   DECOMPRESS/EXPAND OUTPUT ARRAYS FROM THE DETRAINMENT CALCULATIONS
C-----------------------------------------------------------------------
C
C
CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
        DO 40 I=1,NDET
          THPKP1(INDEX1(I)) = THPKP1_C(I)
          QPKP1(INDEX1(I))  = QPKP1_C(I)
          XSQKP1(INDEX1(I)) = XSQKP1_C(I)
C
          BGMKP1(INDEX1(I)) = BGMKP1_C(I)
  40    CONTINUE
      ENDIF
C
      DO 45 I=1,NPNTS
        DELTAK(I) = 0.0
        THRK(I) = 0.0
        QRK(I) = 0.0
  45  CONTINUE
C
CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
      DO 50 I=1,NDET
        DELTAK(INDEX1(I)) = DELTAK_C(I)
        THRK(INDEX1(I))   = THRK_C(I)
        QRK(INDEX1(I))    = QRK_C(I)
  50  CONTINUE
CL
CL----------------------------------------------------------------------
CL  CALCULATE MASS FLUX AT LEVEL K+1.
CL
CL  UM DOCUMENTATION PAPER P27
CL  SECTION (2B), EQUATION (10A)
CL----------------------------------------------------------------------
CL
      DO 60 I=1,NPNTS
       FLXKP1(I) = FLXK(I)*(1.+EKP14(I))*(1.+EKP34(I))*(1.-DELTAK(I))*
     *                                           (1.-AMDETK(I))
  60  CONTINUE
CL
CL---------------------------------------------------------------------
CL TEST FOR POINTS AT WHICH CONVECTION TERMINATES IN LAYER K+1
CL
CL SUBROUTINE TERM_CON
CL
CL UM DOCUMENTATION PAPER P27
CL SECTION (7)
CL---------------------------------------------------------------------
CL
      CALL TERM_CON (NPNTS,NLEV,K,BTERM,BWKP1,FLXKP1,THEKP1,QEKP1,THPI,
     *               QPI,QSEKP1,DELTAK,EXKP1,EKP14,EKP34,PSTAR)
CL
CL----------------------------------------------------------------------
CL CLOUD MICROPHYSICS CALCULATION
CL
CL SUBROUTINE CLOUD_W
CL
CL UM DOCUMENTATION PAPER P27
CL SECTION (8), (9)
CL----------------------------------------------------------------------
CL
      CALL CLOUD_W (K,NPNTS,XPKP1,PREKP1,XSQKP1,BLOWST,FLXKP1,
     *              XPK,THEKP1,QEKP1,BWKP1,BLAND,QSEKP1,BGMKP1,
     *              BTERM,CCA,ICCB,ICCT,TCW,DEPTH,EKP14,EKP34,DELEXKP1,
     *     CCLWP,DELPKP1,CCW,LCCA,LCBASE,LCTOP,LCCLWP,L_SHALLOW,
     *     L_CCW
     &    ,UD_FACTOR
     &     )
C
      RETURN
      END
