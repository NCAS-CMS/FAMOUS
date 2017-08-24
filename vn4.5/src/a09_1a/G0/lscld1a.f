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
C*LL  SUBROUTINE LS_CLD and other---------------------------------------
CLL
CLL  Purpose: Calculates cloud amount, cloud water amounts (ice and
CLL           liquid), and temperature and specific humidity increments
CLL           due to cloud formation, from cloud-conserved and other
CLL           model variables.  This is done for levels 1 to LEVELS
CLL           (specified in the argument list).
CLL  NB: Throughout, levels are counted from the bottom up, i.e. the
CLL      lowest level under consideration is level 1, the next lowest
CLL      level 2, and so on.
CLL
CLL  Suitable for single-column use.
CLL
CLL C.Wilson    <- programmer of some or all of previous code or changes
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  Date
CLL   3.4   5/08/94 Remove calls to TIMER (under *DEF TIMER). R.Rawlins
CLL
!LL   4.0   9/05/95 Changed argument list to export mean cloud water
!LL                 content, QC, and bs for precipitation. A.C.Bushell
!LL
!LL   4.1   9/02/96 Set default QC and bs to RMDI in LS_CLD and removed
!LL                 zero cloud initialization in LS_CLD_C. A.C.Bushell
!LL   4.2    Oct. 96  T3E migration: *DEF CRAY removed; was used to
!LL                   switch on dynamic allocation & WHENIMD. 
!LL                                    S.J.Swarbrick
!LL
!LL   4.4   13/8/97 Several gathers and scatters removed from LSCLD_C
!LL                 in order to reduce run time.
!LL                 D. Salmond
!LL  4.5  27/04/98  Add Fujitsu vectorization directives - needed 
!LL                 because of vn4.4 optimization. RBarnes@ecmwf.int
!LL
CLL  Programming standard: Unified Model Documentation Paper No 4,
CLL                        Version 2, dated 18/1/90.
CLL
CLL  System component covered: P292
CLL
CLL  Documentation: UMDP No.29
CLL
C*L  Arguments:---------------------------------------------------------
      SUBROUTINE LS_CLD(
     + AK,BK,PSTAR,RHCRIT,LEVELS,POINTS,PFIELD,
     & T,CF,Q,QCF,QCL,
     & GRID_QC,BS,ERROR
     +)
      IMPLICIT NONE
      INTEGER
     + LEVELS              ! IN No. of levels being processed.
     +,POINTS              ! IN No. of gridpoints being processed.
     +,PFIELD              ! IN No. of points in global field (at one
C                          !    vertical level).
      REAL
     + PSTAR(PFIELD)       ! IN Surface pressure (Pa).
     +,RHCRIT(LEVELS)      ! IN Critical relative humidity.  See the
C                          !    the paragraph incorporating eqs P292.11
C                          !    to P292.14; the values need to be tuned
C                          !    for the given set of levels.
     +,AK(LEVELS)          ! IN Hybrid "A" co-ordinate.
     +,BK(LEVELS)          ! IN Hybrid "B" co-ordinate.
      REAL
     + Q(PFIELD,LEVELS)    ! INOUT On input: Total water content (QW)
C                          !       (kg per kg air).
C                          !       On output: Specific humidity at
C                          !       processed levels (kg water per kg
C                          !       air).
     +,T(PFIELD,LEVELS)    ! INOUT On input: Liquid/frozen water
C                          !       temperature (TL) (K).
C                          !       On output: Temperature at processed
C                          !       levels (K).
      REAL
     + CF(PFIELD,LEVELS)   ! OUT Cloud fraction at processed levels
C                          !     (decimal fraction).
     +,QCF(PFIELD,LEVELS)  ! OUT Cloud ice content at processed levels
C                          !     (kg per kg air).
     +,QCL(PFIELD,LEVELS)  ! OUT Cloud liquid water content at
C                          !     processed levels (kg per kg air).
     &,GRID_QC(PFIELD,LEVELS)  ! OUT Gridbox mean cloud condensate at
!                                    processed levels (kg per kg air).
!                                    Set to RMDI when cloud is absent.
     &,BS(PFIELD,LEVELS)   ! OUT Maximum moisture fluctuation /6*sigma
!                                at processed levels (kg per kg air).
!                                Set to RMDI when cloud is absent.
      INTEGER ERROR        ! OUT 0 if OK; 1 if bad arguments.
C
C*L  Workspace usage----------------------------------------------------
C    5 blocks of real workspace are required.
      REAL                 ! "Automatic" arrays on Cray.
     & P(POINTS)           ! WORK Pressure at successive levels (Pa).
     &,QSL(POINTS)         ! WORK Saturated spec humidity for temp TL.
     &,QN(POINTS)          ! WORK Cloud water normalised with BS.
      LOGICAL
     & LQC(POINTS)         ! WORK True for points with non-zero cloud
      INTEGER
     & INDEX(POINTS)       ! WORK Index for points with non-zero cloud
C*  Local and other physical constants----------------------------------
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
C*L  External subroutine called ----------------------------------------
      EXTERNAL QSAT,LS_CLD_C
C* Local, including SAVE'd, storage------------------------------------
C
C  (a) Scalars effectively expanded to workspace by the Cray (using
C      vector registers).
C     REAL - None
C
C  (b) Others.
      INTEGER K,I       ! Loop counters: K - vertical level index.
C                       !                I - horizontal field index.
      INTEGER QC_POINTS ! No. points with non-zero cloud

C-----------------------------------------------------------------------
C  Check input arguments for potential over-writing problems.
C-----------------------------------------------------------------------
      ERROR=0
      IF(POINTS.GT.PFIELD)THEN
        ERROR=1
        GOTO1000
      ENDIF
C
C-----------------------------------------------------------------------
CL Subroutine structure :
CL Loop round levels to be processed.
C-----------------------------------------------------------------------
C
      DO K=1,LEVELS
C
C-----------------------------------------------------------------------
CL 1. Calculate QSAT at liquid/ice water temperature, TL,
CL    and initialise cloud ice, water and fraction arrays.
C     This requires a preliminary calculation of the pressure.
C     NB: On entry to the subroutine 'T' is TL and 'Q' is QW.
C-----------------------------------------------------------------------
C
        DO I=1,POINTS
          P(I)=AK(K)+PSTAR(I)*BK(K)
          QCF(I,K)=0.0
          QCL(I,K)=0.0
          CF(I,K) =0.0
          GRID_QC(I,K) = RMDI
          BS(I,K) =RMDI
        ENDDO ! Loop over points
C
        CALL QSAT(QSL,T(1,K),P,POINTS)
C
        DO I=1,POINTS
          IF (RHCRIT(K) .LT. 1.) THEN
C
C-----------------------------------------------------------------------
CL 2. Calculate the quantity QN = QC/BS = (QW/QSL-1)/(1-RHcrit)
CL    if RHcrit is less than 1
C-----------------------------------------------------------------------
C
            QN(I) = (Q(I,K)/QSL(I)-1.)/(1.-RHCRIT(K))
C
C-----------------------------------------------------------------------
CL 3. Set logical variable for cloud, LQC, for the case RHcrit < 1;
C     where QN > -1, i.e. qW/qSAT(TL,P) > RHcrit, there is cloud
C-----------------------------------------------------------------------
C
            LQC(I) = (QN(I) .GT. -1.)
          ELSE
C
C-----------------------------------------------------------------------
CL 2.a Calculate QN = QW - QSL if RHcrit equals 1
C-----------------------------------------------------------------------
C
            QN(I) = Q(I,K) - QSL(I)
C
C-----------------------------------------------------------------------
CL 3.a Set logical variable for cloud, LQC, for the case RHcrit = 1;
CL     where QN > 0, i.e. qW > qSAT(TL,P), there is cloud
C-----------------------------------------------------------------------
C
            LQC(I) = (QN(I) .GT. 0.)
          ENDIF ! Test on RHCRIT
        ENDDO ! Loop over points
C
C-----------------------------------------------------------------------
CL 4. Form index of points where non-zero cloud fraction
C-----------------------------------------------------------------------
C
        QC_POINTS=0
        DO I=1,POINTS
          IF(LQC(I)) THEN
            QC_POINTS=QC_POINTS+1
            INDEX(QC_POINTS)=I
          ENDIF
        ENDDO ! Loop over points
C
C-----------------------------------------------------------------------
CL 5. Call LS_CLD_C to calculate cloud ice and water contents, cloud
CL                  fractions, spec. humidity and determine temperature
C-----------------------------------------------------------------------
C
        IF(QC_POINTS.GT.0) THEN
          CALL LS_CLD_C(P,RHCRIT(K),QSL,QN,Q(1,K),T(1,K),QCF(1,K),
     &                  QCL(1,K),CF(1,K),GRID_QC(1,K),BS(1,K),
     &                  INDEX,QC_POINTS,POINTS)
        ENDIF ! qc_points > 0
C
      ENDDO ! Loop over levels
C
 1000 CONTINUE ! Error exit
      RETURN
      END

C*LL  SUBROUTINE LS_CLD_C-----------------------------------------------
CLL
CLL  Language: FORTRAN 77;  runs under at least IBM and Cray compilers,
CLL            after going through a Cray update-like preprocessor.
CLL
CLL  Suitable for single-column use.
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  Date
CLL
CLL  Programming standard: Unified Model Documentation Paper No 4,
CLL                        Version 1, dated 07/2/91.
CLL
CLL  System component covered: P292
CLL
CLL  Purpose: Calculates cloud water amounts (ice and liquid), cloud
CLL           amounts and temperature and specific humidity
CLL           from cloud-conserved and other model variables.
CLL           This is done for one level.
CLL           Iteration is used to improve the determination of
CLL           ALPHAL, hence AL and so QCF, QCL, Q and T.
CLL
CLL  Documentation: UMDP No.29
CLL
CLL
C*L  Arguments:---------------------------------------------------------
      SUBROUTINE LS_CLD_C(
     & P_F,RHCRIT,QSL_F,QN_F,Q_F,T_F
     &,QCF_F,QCL_F,CF_F,GRID_QC_F,BS_F
     &,INDEX,POINTS,POINTS_F)
      IMPLICIT NONE
      INTEGER
     + POINTS_F            ! IN No. of gridpoints being processed.
     +,POINTS              ! IN No. of gridpoints with non-zero cloud
     +,INDEX(POINTS)       ! IN INDEX for  points with non-zero cloud
C                          !    from lowest model level.
      REAL
     + P_F(POINTS_F)       ! IN pressure (Pa).
     +,QSL_F(POINTS_F)     ! IN saturated humidity at temperature TL,
C                          !    and pressure P_F
     +,QN_F(POINTS_F)      ! IN Normalised super/subsaturation (=QC/BS).
     +,RHCRIT              ! IN Critical relative humidity.  See the
C                          !    the paragraph incorporating eqs P292.11
C                          !    to P292.14;
      REAL
     + Q_F(POINTS_F)       ! INOUT On input: Total water content (QW)
C                          !       (kg per kg air).
C                          !       On output: Specific humidity at
C                          !       processed levels (kg water per kg
C                          !       air).
     +,T_F(POINTS_F)       ! INOUT On input: Liquid/frozen water
C                          !       temperature (TL) (K).
C                          !       On output: Temperature at processed
C                          !       levels (K).
      REAL
     + QCF_F(POINTS_F)     ! OUT Cloud ice content at processed levels
C                          !     (kg per kg air).
     +,QCL_F(POINTS_F)     ! OUT Cloud liquid water content at
C                          !     processed levels (kg per kg air).
     +,CF_F(POINTS_F)      ! OUT Cloud fraction at processed levels.
C
     &,GRID_QC_F(POINTS_F) ! OUT Super/subsaturation on processed levels
!                                Input initialized to RMDI.
     &,BS_F(POINTS_F)      ! OUT Value of bs at processed levels.
!                                Input initialized to RMDI.
C*L  Workspace usage----------------------------------------------------
!    14 blocks of real workspace are required.
      REAL                 ! "Automatic" arrays on Cray.
     & P(POINTS)           ! WORK Pressure  (Pa).
     &,QS(POINTS)          ! WORK Saturated spec humidity for temp T.
     &,QCN(POINTS)         ! WORK Cloud water normalised with BS.
     &,T(POINTS)           ! WORK temperature.
     &,Q(POINTS)           ! WORK specific humidity.
     &,BS(POINTS)          ! WORK Sigmas*sqrt(6): sigmas the parametric
!                                 standard deviation of local cloud
!                                 water content fluctuations.
     &,ALPHAL_NM1(POINTS)  ! WORK ALPHAL at previous iteration.
C*  Local and other physical constants----------------------------------
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
C*L------------------COMDECK C_O_DG_C-----------------------------------
C ZERODEGC IS CONVERSION BETWEEN DEGREES CELSIUS AND KELVIN
C TFS IS TEMPERATURE AT WHICH SEA WATER FREEZES
C TM IS TEMPERATURE AT WHICH FRESH WATER FREEZES AND ICE MELTS
      REAL ZERODEGC,TFS,TM

      PARAMETER(ZERODEGC=273.15,
     &          TFS=271.35,
     &          TM=273.15)
C*----------------------------------------------------------------------

      REAL ALPHF,ALPHL,LSRCP,LCRCP      ! Derived parameters.
     +,LFRCP,CPRLF                      !
      PARAMETER (
     + ALPHF=EPSILON*(LF+LC)/R          ! For frozen AlphaL calculation.
     +,ALPHL=EPSILON*LC/R               ! For liquid AlphaL calculation.
     +,LSRCP=(LF+LC)/CP                 ! Lat ht of sublimation/Cp.
     +,LCRCP=LC/CP                      ! Lat ht of condensation/Cp.
     +,LFRCP=LF/CP                      ! Lat ht of fusion/Cp.
     +,CPRLF=CP/LF                      ! Cp/lat ht of fusion.
     +)
      REAL WTN
      INTEGER
     & ITS                              ! Total number of iterations
     &,N                                ! Iteration counter
      PARAMETER (ITS=5,WTN=0.75)
C*L  External subroutine called ----------------------------------------
      EXTERNAL QSAT
C* Local, including SAVE'd, storage------------------------------------
C
C  (a) Scalars effectively expanded to workspace by the Cray (using
C      vector registers).
      REAL
     + AL                  ! LOCAL AL (see equation P292.6).
     +,ALPHAL              ! LOCAL ALPHAL (see equation P292.5).
     +,TESTT               ! LOCAL temporary temperature for partition
C                          !       of cloud water into ice and liquid
     +,FRACF               ! Fraction of cloud water which is frozen.
C
C  (b) Others.
      INTEGER   I       ! Loop counters: I - horizontal field index.
      INTEGER II
!
!-----------------------------------------------------------------------
!L Gather points with non-zero cloud fraction.
C-----------------------------------------------------------------------
C
        DO I=1,POINTS
          P(I)=P_F(INDEX(I))
          QCN(I)=QN_F(INDEX(I))
        END DO ! Loop over points
C
C-----------------------------------------------------------------------
CL Loop over points with cloud.
C-----------------------------------------------------------------------
!
! Fujitsu vectorization directive 
! Needed because of indirect addressing introduced at vn4.4
!OCL NOVREC
        DO I=1,POINTS
        II=INDEX(I)
!-----------------------------------------------------------------------
!L 1. Calculate ALPHAL (eq P292.5) and AL (P292.6).
!-----------------------------------------------------------------------
!
          IF (T_F(II) .GT. TM) THEN
            ALPHAL = ALPHL * QSL_F(II) / (T_F(II) * T_F(II))   ! P292.5
            AL = 1.0 / (1.0 + (LCRCP * ALPHAL))         ! P292.6
          ELSE
            ALPHAL = ALPHF * QSL_F(II) / (T_F(II) * T_F(II))   ! P292.5
            AL = 1.0 / (1.0 + (LSRCP * ALPHAL))         ! P292.6
          ENDIF
          ALPHAL_NM1(I) = ALPHAL
!
          IF (RHCRIT .LT. 1.) THEN
!-----------------------------------------------------------------------
!L 2. Calculate cloud fraction C, BS (ie. sigma*sqrt(6), where sigma is
!L    as in P292.14) and normalised cloud water QCN=qc/BS, using eqs
!L    P292.15 & 16 if RHcrit < 1.
!  N.B. QN (input) is initially in QCN
!  N.B. QN does not depend on AL and so CF and QCN can be calculated
!       outside the iteration (which is performed in LS_CLD_C).
!       QN is > -1 for all points processed so CF > 0.
!-----------------------------------------------------------------------
!
            BS(I) = (1.0 - RHCRIT) * AL * QSL_F(II)    ! P292.14
            IF (QCN(I) .LE. 0.) THEN
              CF_F(II)=0.5*(1.+QCN(I))*(1.+QCN(I))
              QCN(I)=(1.+QCN(I))*(1.+QCN(I))*(1.+QCN(I))/6.
            ELSEIF (QCN(I) .LT. 1.) THEN
              CF_F(II)=1.-0.5*(1.-QCN(I))*(1.-QCN(I))
              QCN(I)=QCN(I) + (1.-QCN(I))*(1.-QCN(I))*(1.-QCN(I))/6.
            ELSE ! QN .GE. 1
              CF_F(II)=1.
            ENDIF ! Tests on QN
          ELSE ! i.e. if RHcrit =1
C-----------------------------------------------------------------------
!L 3.a Set the cloud fraction to 1 if RHcrit = 1.
C      For the case RHcrit =1, QN is > 0 for all points processed
C      so CF =1.
C-----------------------------------------------------------------------
            BS(I) = AL
            CF_F(II) = 1.
          ENDIF ! Test on RHCRIT
C
C-----------------------------------------------------------------------
CL 3.1 Calculate 1st approx. to qc (store in QCL)
C-----------------------------------------------------------------------
C
          QCL_F(II)=QCN(I)*BS(I)
C
C-----------------------------------------------------------------------
CL 3.2 Calculate 1st approx. specific humidity (total minus cloud water)
C-----------------------------------------------------------------------
C
          Q(I)=Q_F(II)-QCL_F(II)
C
C-----------------------------------------------------------------------
CL 3.3 Perform  partition of cloud water into liquid and ice
CL     components, and calculate 1st approx. to temperature,
CL     accounting for latent heating.
CL     First assume cloud water is all liquid.
C-----------------------------------------------------------------------
C
          T(I)=T_F(II)+LCRCP*QCL_F(II)
          IF(T(I) .GT. TM) THEN          ! Liquid case
            QCF_F(II)=0.0
          ELSE                           ! Frozen or mixed phase
C
C-----------------------------------------------------------------------
CL 3.4 Cloud ice present; either all cloud water is cloud ice and T<TM
CL     or a mixture of ice and liquid and T=TM
C
C      Form test temperature assuming all ice
C      N.B. total cloud water stored in QCL at this stage
C-----------------------------------------------------------------------
C
            TESTT =T(I)+LFRCP*QCL_F(II)
            IF(TESTT .LT. TM) THEN       ! Frozen case
              QCF_F(II)=QCL_F(II)
              T(I)=TESTT
            ELSE                         ! Mixed phase case
              QCF_F(II)= CPRLF*(TM-T(I))
              T(I)=TM
            ENDIF                        ! End frozen
          ENDIF                          ! End liquid
C
C-----------------------------------------------------------------------
CL 3.5 Calculate 1st approx. to cloud liquid water content.
C-----------------------------------------------------------------------
C
          QCL_F(II) = QCL_F(II) - QCF_F(II)
        ENDDO ! Loop over points
C
C-----------------------------------------------------------------------
CL 4. Iteration to find better cloud water values.
C-----------------------------------------------------------------------
C
        IF(ITS.GE.2) THEN
         DO N=2,ITS
C
          CALL QSAT(QS,T,P,POINTS)
C
! Fujitsu vectorization directive
! Needed because of indirect addressing introduced at vn4.4
!OCL NOVREC
          DO I=1,POINTS
           II=INDEX(I)
           IF(T(I).GT.T_F(II)) THEN
C           ! N.B. Cloud water > 0 implies T > TL and so the
C           !        denominator in the following statement is non-zero.
            ALPHAL=(QS(I)-QSL_F(II))/(T(I)-T_F(II))
            ALPHAL=WTN*ALPHAL+(1.0-WTN)*ALPHAL_NM1(I)
            ALPHAL_NM1(I)=ALPHAL
            FRACF=QCF_F(II)/(QCL_F(II)+QCF_F(II))
            AL=1.0/(1.0 + (LCRCP+FRACF*LFRCP)*ALPHAL)
            IF (RHCRIT .LT. 1.) THEN
              BS(I) = (1.0 - RHCRIT) * AL * QSL_F(II)         ! P292.14
            ELSE
              BS(I) = AL
            ENDIF
C
C-----------------------------------------------------------------------
CL 4.1 Calculate Nth approx. to qc (store in QCL).
C-----------------------------------------------------------------------
C
            QCL_F(II)=QCN(I)*BS(I)
C
C-----------------------------------------------------------------------
CL 4.2 Calculate Nth approx. spec. humidity (total minus cloud water).
C-----------------------------------------------------------------------
C
            Q(I)=Q_F(II)-QCL_F(II)
C
C-----------------------------------------------------------------------
CL 4.3 Perform  partition of cloud water into liquid and ice
CL     components, and calculate Nth approx. to temperature,
CL     accounting for latent heating.
CL     First assume cloud water is all liquid.
C-----------------------------------------------------------------------
C
            T(I)=T_F(II)+LCRCP*QCL_F(II)
            IF(T(I) .GT. TM) THEN        ! Liquid case
              QCF_F(II)=0.0
            ELSE                         ! Frozen or mixed phase
C
C-----------------------------------------------------------------------
CL 4.4 Cloud ice present; either all cloud water is cloud ice and T<TM
CL     or a mixture of ice and liquid and T=TM
C
C      Form test temperature assuming all ice
C      N.B. total cloud water stored in QCL at this stage
C-----------------------------------------------------------------------
C
              TESTT =T(I)+LFRCP*QCL_F(II)
              IF(TESTT .LT. TM) THEN     ! Frozen case
                QCF_F(II)=QCL_F(II)
                T(I)=TESTT
              ELSE                       ! Mixed phase case
                QCF_F(II)= CPRLF*(TM-T(I))
                T(I)=TM
              ENDIF                      ! End frozen
            ENDIF                        ! End liquid
C
C-----------------------------------------------------------------------
CL 4.5 Calculate Nth approx. to cloud liquid water content.
C-----------------------------------------------------------------------
C
            QCL_F(II) = QCL_F(II) - QCF_F(II)
           ENDIF ! T > TL
          ENDDO ! Loop over points
         ENDDO ! Loop over iterations
        ENDIF ! ITS ge 2
C
C-----------------------------------------------------------------------
CL 5. Finally scatter back results
C-----------------------------------------------------------------------
C
CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
      DO I=1,POINTS
        Q_F(INDEX(I)) = Q(I)
        T_F(INDEX(I)) = T(I)
        GRID_QC_F(INDEX(I)) = BS(I) * QN_F(INDEX(I))
        BS_F(INDEX(I)) = BS(I)
      END DO ! Loop over points
C
      RETURN
      END
