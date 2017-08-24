! ******************************COPYRIGHT******************************
! (c) CROWN COPYRIGHT 1998, METEOROLOGICAL OFFICE, All Rights Reserved.
!
! Use, duplication or disclosure of this code is subject to the
! restrictions as set forth in the contract.
!
!                Meteorological Office
!                London Road
!                BRACKNELL
!                Berkshire UK
!                RG12 2SZ
!
! If no contract has been raised with this copy of the code, the use,
! duplication or disclosure of it is strictly prohibited.  Permission
! to do so must first be obtained in writing from the Head of Numerical
! Modelling at the above address.
! ******************************COPYRIGHT******************************
!
!+ Large-scale Cloud Scheme.
! Subroutine Interface:
      SUBROUTINE LS_CLD(
!      Pressure related fields
     & AK, BK, PSTAR
!      Array dimensions
     &,LEVELS, RHCPT, POINTS, PFIELD
!      Prognostic Fields
     &,T, CF, Q, QCF, QCL
!      Liquid and frozen ice cloud fractions
     &,CFL, CFF
     &,ERROR)
!
      IMPLICIT NONE
!
! Purpose:
!   This subroutine calculates liquid and ice cloud fractional cover
!   for use with the enhanced precipitation microphysics scheme.
!
! Method:
!   Statistical cloud scheme separates input moisture into specific
!   humidity and cloud liquid water. Temperature calculated from liquid
!   water temperature. Cloud fractions calculated from statistical
!   relation between cloud fraction and cloud liquid/ice water content.
!   Critical relative humidity now specified for all grid cells.
!
! Current Owner of Code: S. Cusack
!
! History:
! Version   Date     Comment
!  4.5    14-05-98   Original Code     S. Cusack
!
! Description of Code:
!   FORTRAN 77  + common extensions also in Fortran90.
!   This code is written to UMDP3 version 6 programming standards.
!
!   System component covered: P292
!
!   Documentation: UMDP No.29
!
!  Global Variables:----------------------------------------------------
! History:
! Version  Date  Comment
!  3.4   18/5/94 Add PP missing data indicator. J F Thomson
C*L------------------COMDECK C_MDI-------------------------------------
      REAL    RMDI_PP   ! PP missing data indicator (-1.0E+30)
      REAL    RMDI_OLD  ! Old real missing data indicator (-32768.0)
      REAL    RMDI      ! New real missing data indicator (-2**30)
      INTEGER IMDI      ! Integer missing data indicator

      PARAMETER(
     & RMDI_PP  =-1.0E+30,
     & RMDI_OLD =-32768.0,
     & RMDI =-32768.0*32768.0,
     & IMDI =-32768)
C*----------------------------------------------------------------------
C*L------------------COMDECK C_PI---------------------------------------
CLL
CLL 4.0 19/09/95  New value for PI. Old value incorrect
CLL               from 12th decimal place. D. Robinson
CLL
      REAL PI,PI_OVER_180,RECIP_PI_OVER_180

      PARAMETER(
     & PI=3.14159265358979323846, ! Pi
     & PI_OVER_180 =PI/180.0,   ! Conversion factor degrees to radians
     & RECIP_PI_OVER_180 = 180.0/PI ! Conversion factor radians to
     &                              ! degrees
     & )
C*----------------------------------------------------------------------

!
!  Subroutine Arguments:------------------------------------------------
      INTEGER           !, INTENT(IN)
     & LEVELS
!       No. of levels being processed.
     &,POINTS
!       No. of gridpoints being processed.
     &,PFIELD
!       No. of points in global field (at one vertical level).
!
      REAL              !, INTENT(IN)
     & QCF(PFIELD,LEVELS)
!       Cloud ice content at processed levels (kg water per kg air).
     &,PSTAR(PFIELD)
!       Surface pressure (Pa).
     &,AK(LEVELS)
!       Hybrid "A" co-ordinate.
     &,BK(LEVELS)
!       Hybrid "B" co-ordinate.
     &,RHCPT(PFIELD,LEVELS)
!       Critical relative humidity for all points
!
      REAL              !, INTENT(INOUT)
     & Q(PFIELD,LEVELS)
!       On input : Total water content (QW) (kg per kg air).
!       On output: Specific humidity at processed levels
!                   (kg water per kg air).
     &,T(PFIELD,LEVELS)
!       On input : Liquid/frozen water temperature (TL) (K).
!       On output: Temperature at processed levels (K).
!
      REAL              !, INTENT(OUT)
     & CF(PFIELD,LEVELS)
!       Cloud fraction at processed levels (decimal fraction).
     &,QCL(PFIELD,LEVELS)
!       Cloud liquid water content at processed levels (kg per kg air).
     &,CFL(PFIELD,LEVELS)
!       Liquid cloud fraction at processed levels (decimal fraction).
     &,CFF(PFIELD,LEVELS)
!       Frozen cloud fraction at processed levels (decimal fraction).
!
!     Error Status:
      INTEGER ERROR     !, INTENT(OUT)  0 if OK; 1 if bad arguments.
!
!  Local parameters and other physical constants------------------------
      REAL ROOTWO       ! Sqrt(2.)
!
!  Local scalars--------------------------------------------------------
!
!  (a) Scalars effectively expanded to workspace by the Cray (using
!      vector registers).
      REAL
     & QCFRBS           ! qCF / bs
     &,PHIQCF           ! Arc-cosine term in Cloud ice fraction calc.
     &,COSQCF           ! Cosine term in Cloud ice fraction calc.
!
!  (b) Others.
      INTEGER K,I       ! Loop counters: K - vertical level index.
!                       !                I - horizontal field index.
      INTEGER QC_POINTS ! No. points with non-zero cloud
!
!  Local dynamic arrays-------------------------------------------------
!    7 blocks of real workspace are required.
      REAL
     & P(POINTS)
!       Pressure at successive levels (Pa).
     &,QSL(POINTS)
!       Saturated specific humidity for temp TL or T.
     &,QN(POINTS)
!       Cloud water normalised with BS.
     &,GRID_QC(POINTS,LEVELS)
!       Gridbox mean saturation excess at processed levels
!        (kg per kg air). Set to RMDI when cloud is absent.
     &,BS(POINTS,LEVELS)
!       Maximum moisture fluctuation /6*sigma at processed levels
!        (kg per kg air). Set to RMDI when cloud is absent.
      LOGICAL
     & LQC(POINTS)         ! True for points with non-zero cloud
      INTEGER
     & INDEX(POINTS)       ! Index for points with non-zero cloud
!
!  External subroutine calls: ------------------------------------------
      EXTERNAL QSAT,QSAT_WAT,LS_CLD_C
!- End of Header
! ----------------------------------------------------------------------
!  Check input arguments for potential over-writing problems.
! ----------------------------------------------------------------------
      ERROR=0
      IF (POINTS.GT.PFIELD) THEN
        ERROR=1
        GO TO 9999
      END IF
!
! ==Main Block==--------------------------------------------------------
! Subroutine structure :
! Loop round levels to be processed.
! ----------------------------------------------------------------------
! Levels_do1:
      DO K=1,LEVELS
!
! ----------------------------------------------------------------------
! 1. Calculate QSAT at liquid/ice water temperature, TL, and initialize
!    cloud water, sub-grid distribution and fraction arrays.
!    This requires a preliminary calculation of the pressure.
!    NB: On entry to the subroutine 'T' is TL and 'Q' is QW.
! ----------------------------------------------------------------------
! Points_do1:
        DO I=1, POINTS
          P(I) = AK(K) + PSTAR(I)*BK(K)
          QCL(I,K) = 0.0
          CFL(I,K) = 0.0
          GRID_QC(I,K) = RMDI
          BS(I,K) = RMDI
        END DO ! Points_do1
!
        CALL QSAT_WAT(QSL,T(1,K),P,POINTS)
!
! Points_do2:
        DO I=1, POINTS
! Rhcrit_if:
          IF (RHCPT(I,K) .LT. 1.) THEN
! ----------------------------------------------------------------------
! 2. Calculate the quantity QN = QC/BS = (QW/QSL-1)/(1-RHcrit)
!    if RHcrit is less than 1
! ----------------------------------------------------------------------
!
            QN(I) = (Q(I,K) / QSL(I) - 1.) / (1. - RHCPT(I,K))
!
! ----------------------------------------------------------------------
! 3. Set logical variable for cloud, LQC, for the case RHcrit < 1;
!    where QN > -1, i.e. qW/qSAT(TL,P) > RHcrit, there is cloud
! ----------------------------------------------------------------------
!
            LQC(I) = (QN(I) .GT. -1.)
          ELSE
! ----------------------------------------------------------------------
! 2.a Calculate QN = QW - QSL if RHcrit equals 1
! ----------------------------------------------------------------------
!
            QN(I) = Q(I,K) - QSL(I)
!
! ----------------------------------------------------------------------
! 3.a Set logical variable for cloud, LQC, for the case RHcrit = 1;
!     where QN > 0, i.e. qW > qSAT(TL,P), there is cloud
! ----------------------------------------------------------------------
!
            LQC(I) = (QN(I) .GT. 0.)
          END IF ! Rhcrit_if
        END DO ! Points_do2
!
! ----------------------------------------------------------------------
! 4. Form index of points where non-zero liquid cloud fraction
! ----------------------------------------------------------------------
!
! Points_do3:
        QC_POINTS=0
        DO I=1,POINTS
          IF (LQC(I)) THEN
            QC_POINTS = QC_POINTS + 1
            INDEX(QC_POINTS) = I
          END IF
        END DO ! Points_do3
!
! ----------------------------------------------------------------------
! 5. Call LS_CLD_C to calculate cloud water content, specific humidity,
!                  water cloud fraction and determine temperature.
! ----------------------------------------------------------------------
! Qc_points_if:
        IF (QC_POINTS .GT. 0) THEN
          CALL LS_CLD_C(P,QSL,QN,Q(1,K),T(1,K)
     &                 ,QCL(1,K),CFL(1,K),GRID_QC(1,K),BS(1,K)
     &                 ,INDEX,QC_POINTS,POINTS,RHCPT(1,K))
        END IF ! Qc_points_if
!
! ----------------------------------------------------------------------
! 6. Calculate cloud fractions for ice clouds.
!    THIS IS STILL HIGHLY EXPERIMENTAL.
!    Begin by calculating Qsat(T,P), at Temperature, for estimate of bs.
! ----------------------------------------------------------------------
        CALL QSAT_WAT(QSL,T(1,K),P,POINTS)
        ROOTWO = SQRT(2.)
!
! Points_do4:
        DO I=1, POINTS
! ----------------------------------------------------------------------
! 6a Calculate qCF/bs.
! ----------------------------------------------------------------------
          QCFRBS =  QCF(I,K) / ((1. - RHCPT(I,K)) * QSL(I))
!
! ----------------------------------------------------------------------
! 6b Calculate frozen cloud fraction from frozen cloud water content.
! ----------------------------------------------------------------------
          IF (QCFRBS .LE. 0.) THEN
            CFF(I,K) = 0.0
          ELSEIF (0. .LT. QCFRBS  .AND. (6. * QCFRBS) .LE. 1.) THEN
            CFF(I,K) = 0.5 * ((6. * QCFRBS)**(2./3.))
          ELSEIF (1. .LT. (6.*QCFRBS) .AND. QCFRBS .LT. 1.) THEN
            PHIQCF = ACOS(ROOTWO * 0.75 * (1. - QCFRBS))
            COSQCF = COS((PHIQCF + (4. * PI)) / 3.)
            CFF(I,K) = 1. - (4. * COSQCF * COSQCF)
          ELSEIF (QCFRBS .GE. 1.) THEN
            CFF(I,K) = 1.
          END IF
! ----------------------------------------------------------------------
! 6c Calculate combined cloud fraction.
! ----------------------------------------------------------------------
!         Use maximum overlap condition
!         CF(I,K) = MAX(CFL(I,K), CFF(I,K))
!
!         Use minimum overlap condition
          CF(I,K) = MIN(CFL(I,K)+CFF(I,K), 1.0)
!
        END DO ! Points_do4
!
      END DO ! Levels_do
!
 9999 CONTINUE ! Error exit
      RETURN
      END
! ======================================================================
!
!+ Large-scale Cloud Scheme Compression routine (Cloud points only).
! Subroutine Interface:
      SUBROUTINE LS_CLD_C(
     & P_F,QSL_F,QN_F,Q_F,T_F
     &,QCL_F,CF_F,GRID_QC_F,BS_F
     &,INDEX,POINTS,POINTS_F,RHCPT_F)
      IMPLICIT NONE
!
! Purpose: Calculates liquid cloud water amounts and cloud amounts,
!          temperature and specific humidity from cloud-conserved and
!          other model variables. This is done for one model level.
!
! Current Owner of Code: S. Cusack
!
! History:
! Version   Date     Comment
!  4.5    12-05-98   Original Code
!
! Description of Code:
!   FORTRAN 77  + common extensions also in Fortran90.
!   This code is written to UMDP3 version 6 programming standards.
!
!   System component covered: P292
!
!   Documentation: UMDP No.29
!
!  Global Variables:----------------------------------------------------
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

C*L------------------COMDECK C_LHEAT------------------------------------
C LC IS LATENT HEAT OF CONDENSATION OF WATER AT 0DEGC
C LF IS LATENT HEAT OF FUSION AT 0DEGC
      REAL LC,LF

      PARAMETER(LC=2.501E6,
     &          LF=0.334E6)
C*----------------------------------------------------------------------
!
!  Subroutine Arguments:------------------------------------------------
      INTEGER           !, INTENT(IN)
     & POINTS_F
!       No. of gridpoints being processed.
     &,POINTS
!       No. of gridpoints with non-zero cloud
     &,INDEX(POINTS)
!       INDEX for  points with non-zero cloud from lowest model level.
!
      REAL              !, INTENT(IN)
     & P_F(POINTS_F)
!       pressure (Pa).
     &,QSL_F(POINTS_F)
!       saturated humidity at temperature TL, and pressure P_F
     &,QN_F(POINTS_F)
!       Normalised super/subsaturation ( = QC/BS).
     &,RHCPT_F(POINTS_F)
!       Critical relative humidity in all grid-cells.
!
      REAL              !, INTENT(INOUT)
     & Q_F(POINTS_F)
!       On input : Vapour + liquid water content (QW) (kg per kg air).
!       On output: Specific humidity at processed levels
!                   (kg water per kg air).
     &,T_F(POINTS_F)
!       On input : Liquid water temperature (TL) (K).
!       On output: Temperature at processed levels (K).
!
      REAL              !, INTENT(OUT)
     & QCL_F(POINTS_F)
!       Cloud liquid water content at processed levels (kg per kg air).
     &,CF_F(POINTS_F)
!       Liquid cloud fraction at processed levels.
     &,GRID_QC_F(POINTS_F)
!       Super/subsaturation on processed levels. Input initially RMDI.
     &,BS_F(POINTS_F)
!       Value of bs at processed levels. Input initialized to RMDI.
!
!  Local parameters and other physical constants------------------------
      REAL ALPHL,LCRCP                  ! Derived parameters.
      PARAMETER (
     & ALPHL=EPSILON*LC/R               ! For liquid AlphaL calculation.
     &,LCRCP=LC/CP                      ! Lat ht of condensation/Cp.
     &)
      REAL WTN                          ! Weighting for ALPHAL iteration
      INTEGER
     & ITS                              ! Total number of iterations
      PARAMETER (ITS=5,WTN=0.75)
!
!  Local scalars--------------------------------------------------------
!
!  (a) Scalars effectively expanded to workspace by the Cray (using
!      vector registers).
      REAL
     & AL                ! LOCAL AL (see equation P292.6).
     &,ALPHAL            ! LOCAL ALPHAL (see equation P292.5).
!
!  (b) Others.
      INTEGER   I,II,N   ! Loop counters: I,II - horizontal field index.
!                                       : N - iteration counter.
!
!  Local dynamic arrays-------------------------------------------------
!    8 blocks of real workspace are required.
      REAL
     & P(POINTS)
!       Pressure  (Pa).
     &,QS(POINTS)
!       Saturated spec humidity for temp T.
     &,QCN(POINTS)
!       Cloud water normalised with BS.
     &,T(POINTS)
!       temperature.
     &,Q(POINTS)
!       specific humidity.
     &,BS(POINTS)
!       Sigmas*sqrt(6): sigmas the parametric standard deviation of
!       local cloud water content fluctuations.
     &,ALPHAL_NM1(POINTS)
!       ALPHAL at previous iteration.
!
!  External subroutine calls: ------------------------------------------
      EXTERNAL QSAT_WAT
!
!- End of Header
!
! ==Main Block==--------------------------------------------------------
! Operate on INDEXed points with non-zero cloud fraction.
! ----------------------------------------------------------------------
! Points_do1:
      DO I=1, POINTS
        II = INDEX(I)
        P(I)  = P_F(II)
        QCN(I)= QN_F(II)
! ----------------------------------------------------------------------
! 1. Calculate ALPHAL (eq P292.5) and AL (P292.6).
!    CAUTION: T_F acts as TL (input value) until update in final section
!    CAUTION: Q_F acts as QW (input value) until update in final section
! ----------------------------------------------------------------------
!
        ALPHAL = ALPHL * QSL_F(II) / (T_F(II) * T_F(II))       ! P292.5
        AL = 1.0 / (1.0 + (LCRCP * ALPHAL))                    ! P292.6
        ALPHAL_NM1(I) = ALPHAL
!
! Rhcrit_if1:
        IF (RHCPT_F(II) .LT. 1.) THEN
! ----------------------------------------------------------------------
! 2. Calculate cloud fraction CF, BS (ie. sigma*sqrt(6), where sigma is
!    as in P292.14) and normalised cloud water QCN=qc/BS, using eqs
!    P292.15 & 16 if RHcrit < 1.
! N.B. QN (input) is initially in QCN
! N.B. QN does not depend on AL and so CF and QCN can be calculated
!      outside the iteration (which is performed in LS_CLD_C).
!      QN is > -1 for all points processed so CF > 0.
! ----------------------------------------------------------------------
!
          BS(I) = (1.0 - RHCPT_F(II)) * AL * QSL_F(II)         ! P292.14
          IF (QCN(I) .LE. 0.) THEN
            CF_F(II) = 0.5 * (1. + QCN(I)) * (1. + QCN(I))
            QCN(I)= (1. + QCN(I)) * (1. + QCN(I)) * (1. + QCN(I)) / 6.
          ELSEIF (QCN(I) .LT. 1.) THEN
            CF_F(II) = 1. - 0.5 * (1. - QCN(I)) * (1. - QCN(I))
            QCN(I)=QCN(I) + (1.-QCN(I)) * (1.-QCN(I)) * (1.-QCN(I))/6.
          ELSE ! QN .GE. 1
            CF_F(II) = 1.
          END IF ! QCN_if
        ELSE ! i.e. if RHcrit = 1
! ----------------------------------------------------------------------
! 3.a If RHcrit = 1., all points processed have QN > 0 and CF = 1.
! ----------------------------------------------------------------------
          BS(I) = AL
          CF_F(II) = 1.
        END IF ! Rhcrit_if1
!
! ----------------------------------------------------------------------
! 3.1 Calculate 1st approx. to qc (store in QCL)
! ----------------------------------------------------------------------
!
        QCL_F(II) = QCN(I) * BS(I)
!
! ----------------------------------------------------------------------
! 3.2 Calculate 1st approx. specific humidity (total minus cloud water)
! ----------------------------------------------------------------------
!
        Q(I) = Q_F(II) - QCL_F(II)
!
! ----------------------------------------------------------------------
! 3.3 Calculate 1st approx. to temperature, adjusting for latent heating
! ----------------------------------------------------------------------
!
        T(I) = T_F(II) + LCRCP*QCL_F(II)
      END DO ! Points_do1
!
! ----------------------------------------------------------------------
! 4. Iteration to find better cloud water values.
! ----------------------------------------------------------------------
! Its_if:
      IF (ITS .GE. 2) THEN
! Its_do:
        DO N=2, ITS
!
          CALL QSAT_WAT(QS,T,P,POINTS)
! Points_do2:
          DO I=1, POINTS
            II = INDEX(I)
! T_if:
            IF (T(I) .GT. T_F(II)) THEN
!           NB. T > TL implies cloud fraction > 0.
              ALPHAL = (QS(I) - QSL_F(II)) / (T(I) - T_F(II))
              ALPHAL = WTN * ALPHAL + (1.0 - WTN) * ALPHAL_NM1(I)
              ALPHAL_NM1(I) = ALPHAL
              AL = 1.0 / (1.0 + (LCRCP * ALPHAL))
! Rhcrit_if2:
              IF (RHCPT_F(II) .LT. 1.) THEN
                BS(I) = (1.0 - RHCPT_F(II)) * AL * QSL_F(II)   ! P292.1
              ELSE
                BS(I) = AL
              END IF  ! Rhcrit_if2
!
! ----------------------------------------------------------------------
! 4.1 Calculate Nth approx. to qc (store in QCL).
! ----------------------------------------------------------------------
!
              QCL_F(II) = QCN(I) * BS(I)
!
! ----------------------------------------------------------------------
! 4.2 Calculate Nth approx. spec. humidity (total minus cloud water).
! ----------------------------------------------------------------------
!
              Q(I) = Q_F(II) - QCL_F(II)
!
! ----------------------------------------------------------------------
! 4.3 Calculate Nth approx. to temperature, adjusting for latent heating
! ----------------------------------------------------------------------
!
              T(I) = T_F(II) + LCRCP * QCL_F(II)
!
            END IF ! T_if
          END DO ! Points_do2
        END DO ! Its_do
      END IF ! Its_if
!
! ----------------------------------------------------------------------
! 5. Finally scatter back cloud point results to full field arrays.
!    CAUTION: T_F updated from TL (input) to T (output)
!    CAUTION: Q_F updated from QW (input) to Q (output)
! ----------------------------------------------------------------------
!
CDIR$ IVDEP
! Points_do3:
      DO I=1,POINTS
        Q_F(INDEX(I)) = Q(I)
        T_F(INDEX(I)) = T(I)
        GRID_QC_F(INDEX(I)) = BS(I) * QN_F(INDEX(I))
        BS_F(INDEX(I)) = BS(I)
      END DO ! Points_do3
!
      RETURN
      END
! ======================================================================
