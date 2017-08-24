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
C*LL  SUBROUTINE LSP_FORM-----------------------------------------------
!LL
!LL  Purpose: Form or augment precipitation at the expense of cloud
!LL           water in one model layer.
!LL
!LL C.Senior    <- programmer of some or all of previous code or changes
!LL
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!LL 3.1      23/02/93 Fully divergent form of ice fallout including:
!LL                   1) fully implicit finite difference scheme
!LL                   2) total cloud water used in determining RCL and
!LL                      RCF instead of QL and QF used previously
!LL                   3) separate precipitation rates for ice (RCF) and
!LL                      liquid water (RCL) in mixed phase region
!LL                   4) minimum cloud for precipitation, CFMIN set
!LL                      equal to zero, in comdeck C_LSPFRM.
!LL                                             Ruth Carnell 26/02/93
!LL  Model            Modification history from model version 3.0:
!LL version  Date
!LL 3.3      22/11/94 Distribution of moisture within gridbox added:
!LL                   follows form of LSCLD1A routine / 2A variant.
!LL                                           Andrew Bushell 15/12/94
!LL
!LL 3.3      04/04/95 Possible fix to long timestep solution of ice
!LL                   fallout + explicit formulation for short ts case.
!LL                                    D Gregory & A Bushell 04/04/95
!LL
!LL 3.5      28/09/95 Changed *CALL FOCWWIL to CALL LSP_FOCWWIL to allow
!LL                   precompilation of alternative FOCWWIL options.
!LL                                           Andrew Bushell 28/09/95
!LL
!LL 4.1      14/02/96 Tidied code to remove intermittent problem with
!LL                   uninitialized LSCQC and LSCBS. Bit comparable.
!LL                                           Andrew Bushell 14/02/96
!LL
!LL 4.5      31/03/98 Code change to avoid underflow divide by zero on  
!LL                   some machines. Bit comparable. Andrew Bushell
!LL                                                                     
!LL  Programming standard: Unified Model Documentation Paper No 4,
!LL                        Version 1, dated 12/9/89.
!LL
!LL  Logical component covered: Part of P26.
!LL
!LL  System task:
!LL
!LL  Documentation: Unified Model Documentation Paper No 26.
C*
C*L  Arguments:---------------------------------------------------------
      SUBROUTINE LSP_FORM(
     & CF,P,Q,RHODZ,T,TIMESTEP,POINTS,QCF,QCL,RAIN,SNOW
     &,F_DELTA_SNOW,BLAND,CW_SEA,CW_LAND,LSCQC,LSCBS,VFALL
     &)
      IMPLICIT NONE
      INTEGER         ! Input integer scalar :-
     & POINTS         ! IN Number of points to be processed.
      REAL            ! Input real arrays :-
     & CF(POINTS)     ! IN Cloud fraction in this layer.
     &,P(POINTS)      ! IN Air pressure at this level (Pa).
     &,Q(POINTS)      ! IN Sp humidity at this level (kg wat per kg air)
     &,RHODZ(POINTS)  ! IN Air mass p.u.a. in this layer (kg per sq m).
     &,T(POINTS)      ! IN Temperature at this level (K).
     &,LSCQC(POINTS)  ! IN Large scale cloud Qc value.
     &,LSCBS(POINTS)  ! IN Large scale cloud bs value.
      REAL            ! Input real scalar :-
     & TIMESTEP       ! IN Timestep (s).
     &,CW_SEA         ! IN threshold cloud liquid water content
!                          over sea for conversion to ppn
!                          (kg water per m**3)
     &,CW_LAND        ! IN threshold cloud liquid water content
!                          over land for conversion to ppn
!                          (kg water per m**3)
      REAL            ! Updated real arrays :-
     & QCF(POINTS)    ! INOUT Cloud ice (kg water per kg air).
     &,QCL(POINTS)    ! INOUT Cloud liquid water (kg water per kg air).
     &,RAIN(POINTS)   ! INOUT Rate of rainfall entering this layer from
!                             above (on i/p), or leaving it (on o/p)
!                             (kg per sq m per s).
     &,SNOW(POINTS)   ! INOUT Rate of snowfall entering this layer from
!                             above (on i/p), or leaving it (on o/p)
!                             (kg per sq m per s).
     &,F_DELTA_SNOW(POINTS) ! INOUT snow fraction from ice falling as
!                                   water
     &,VFALL(POINTS)  ! INOUT Fall velocity of ice into layer (m/s) i/p
!                             Fall velocity of ice from layer (m/s) o/p
      LOGICAL BLAND(POINTS)  ! IN Land/sea mask
C*L  Workspace usage----------------------------------------------------
!  1 real, 0 logical and 0 integer blocks are required, as follows :-
      REAL
     & FL(POINTS)     ! Fraction of total cloud water which is liquid.
C*L   External subprogram called :-
      EXTERNAL LSP_FOCWWIL  !NB: ERF_LSP is an in-line comdeck
      REAL ERF_LSP, XERF, TERF
C*
!-----------------------------------------------------------------------
!  Common, then local, physical constants.
!-----------------------------------------------------------------------
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
C*LL--------------------------------------------------------------------
!LL  Tunable constants local to routine LSP_FORM, called by P26
!LL  LS_PPN.
!LL---------------------------------------------------------------------
      REAL CA,PCF,CFMIN,CT,FPL,VF1,VFPOWER
      PARAMETER (
     &CA=1.0,           ! Accretion constant (kg per m**2).
     &PCF=1.01086E-3,   ! Tunable in Heymsfield formula (kg/m**3).
     &CFMIN=0.0,        ! Minimum cloud amount for ppn calculations.
     &CT=1.0E-4,        ! Reciprocal time constant for conversion of
!                         cloud liquid water to ppn (Hertz).
     &VF1=1.0,          ! Tunable speed in Heymsfield formula (m/s).
     &VFPOWER=0.17      ! Local parameter, used in Heymsfield formula.
     &)
C*
!
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

      REAL      SQRT_PI   ! Defined for use by ERF
!-----------------------------------------------------------------------
!  Define local scalars.
!-----------------------------------------------------------------------
!  (a) Reals effectively expanded to workspace by the Cray (using
!      vector registers).
      REAL            ! Real workspace.  At end of DO loop, contains :-
     & QC             ! Total cloud water (liquid + ice).
     &,RCL            ! Rate of liquid water exchange
     &,RCF            ! Rate of ice exchange
     &,DELTA          ! 'ice' to be treated as water
     &,DELTA_SNOW     ! Snow formed by ice falling as water
     &,RHO            ! Density of air in the layer.
     &,VFRDZ          ! Ice-particle fallout speed (m/s), divided by
!                       layer depth (m).
     &,QCF_U          ! Cloud ice nominal update value (kg/kg).
     &,VTEMP          ! Virtual temperature as at start of loop.
     &,CW             ! Parameter for efficient conversion of
!                       liquid cloud water to ppn
     &,QCERR          ! Calculated qc from bs and Qc(lscld)
     &,QSUM           ! Qc(lscld) + bs
     &,QDIF           ! Qc(lscld) - bs
     &,NORCT          ! (Ct * Cw**3) / (qc * rho**3 * bs**2 * 2.)
     &,NERF_SUM       ! ERF_LSP(qsum * rho / Cw)
     &,NERF_DIF       ! ERF_LSP(qdif * rho / Cw)
!  (b) Others.
      INTEGER I       ! Loop counter (horizontal field index).
     &,NV             ! Index counter for 'when' routine.
!->>>-------------------------------------------------------------------
      SQRT_PI = SQRT(PI)
!
!-----------------------------------------------------------------------
!L 0. Calculate fraction of cloud water which is liquid (FL(I)),
!L    according to equation P26.50.
!-----------------------------------------------------------------------
      CALL LSP_FOCWWIL(T, POINTS, FL)
!
      DO I=1,POINTS
!-----------------------------------------------------------------------
!     Store total cloud water in QC (see eq P26.30) and ice falling as
!     water in DELTA. Also update QCF.
!-----------------------------------------------------------------------
!
        QC=QCL(I)+QCF(I)
        DELTA = (FL(I) * QC) - QCL(I)
        QCF(I) = (1. - FL(I)) * QC
!
!-----------------------------------------------------------------------
!  1. Perform error check on input LSCQC and LSCBS which should be able
!     to reproduce QC. Commented out as this will stop multitasking.
!-----------------------------------------------------------------------
!----
!       IF (LSCBS(I) .LT. 0.) THEN
!-------Should equal RMDI indicating no cloud = no cloud distribution---
!        IF (QC .GT. 0.) WRITE(6,*) '%Error @ ',I,'QC: ',QC,':',LSCBS(I)
!         QCERR = 0.0
!       ELSE
!         IF (LSCQC(I) .LE. 0.) THEN
!           QSUM = LSCQC(I) + LSCBS(I)
!           QCERR = QSUM * QSUM * QSUM / (6. * LSCBS(I) * LSCBS(I))
!         ELSEIF (0. .LT. LSCQC(I) .AND. LSCQC(I) .LT. LSCBS(I)) THEN
!           QSUM = LSCQC(I) + LSCBS(I)
!           QCERR = ((QSUM*QSUM*QSUM) - (2.*LSCQC(I)*LSCQC(I)*LSCQC(I)))
!    &              / (6. * LSCBS(I) * LSCBS(I))
!         ELSE   ! LSCQC(I) .GE. LSCBS(I)
!           QCERR = LSCQC(I)
!         ENDIF  ! LSC_QC ranges
!----
!         IF (QC .LE. 0.) THEN
!---------Might happen due to rounding errors. Large discrepancy would
!---------indicate inconsistency: cloud distribution won't be used tho'.
!           WRITE (6,*) '%QC @ ',I,' : ',QC,' /= QCERR ',QCERR
!         ELSEIF (((QCERR/QC)-1.) * ((QCERR/QC)-1.) .GT. 1.0E-18) THEN
!-------Only serious if discrepancy is above bit resolution error level-
!---------(Needs to be set for given machine)---------------------------
!           IF (QC .GT. 1.0E-18) THEN
!             WRITE (6,*) '%CAUTION @ ',I,' : QC ',QC,' /= QCERR ',QCERR
!           ENDIF
!         ENDIF ! Qc le 0
!       ENDIF  ! LSC_bs missing data
!-----------------------------------------------------------------------
!L 2. Calculate fractional rate of conversion of liquid cloud water to
!L    precipitation.
!     Store in RCL.
!-----------------------------------------------------------------------
!
!  (a) Calculate density of air, RHO, via virtual temperature.
!
        VTEMP=T(I)*(1.+C_VIRTUAL*Q(I)-QC)      ! Virtual temperature
        RHO=P(I)/(R*VTEMP)
!
!  (b) Calculate P/CA.  Store in RCL.
!
        RCL=(RAIN(I)+SNOW(I))/CA
!
!  (c) Calculate term based on moisture distribution in gridbox as long
!      as there is cloud present (QC from LS_CLD > -bs). Total qC used.
!      (Also initialise variable calculated in sect 3 below, q.v.)
!
        VFRDZ=0.0
!
        IF ((QC * LSCBS(I) * LSCBS(I)) .GT. 0.0) THEN 
!       ie. cloud water content > 0
! NB:QC LE 0 or ((-LSCBS(I)) .GE. LSCQC(I)) => leave RCL, RCF unchanged.
!
          QSUM = LSCQC(I) + LSCBS(I)
          QDIF = LSCQC(I) - LSCBS(I)
!
          IF (BLAND(I)) THEN
            CW=CW_LAND
          ELSE
            CW=CW_SEA
          ENDIF
!
          NORCT= CW / RHO
          NORCT= 0.5 * CT * NORCT * NORCT * NORCT / (QC * LSCBS(I)
     &               * LSCBS(I))
!
          XERF = (RHO * QSUM / CW)
C*LL  IN-LINE FUNCTION ERF_LSP------------------------------------------
!LL
!LL     Purpose:
!LL     Vectorizable calculation of (x**3 / 3) - x + (PI**0.5 * ERF(x)/2
!LL     needed for LSPFOR routine evaluation of autoconversion rate.
!LL     Need to declare REAL ERF_LSP, XERF, TERF and to calculate
!LL     SQRT_PI = SQRT(PI) in calling routine.
!LL
!LL  Author: A Bushell
!LL
!LL  Model            Modification history from model version 4.0:
!LL version  date
!LL
!LL  Documentation:
!LL For abs(XERF) <= 1.161, use the lowest order terms of the integrated
!LL                         exp(-XERF**2) Maclaurin expansion. For a
!LL                         mantissa precision error 1E-13, this will be
!LL                         EXACT for abs(XERF) < 0.03 with max error at
!LL                         1.161 of about 0.01%.
!LL For abs(XERF) >  1.161, Press, Teukolsky, Vettering & Flannery
!LL                         NUMERICAL RECIPES IN FORTRAN eq 6.2.6 gives
!LL                         a continued fraction representation. The 1st
!LL                         4 terms give ERF to high accuracy for
!LL                         abs(XERF) > 2 & max error < 0.02% at 1.161.
!LL
!LL NOTE: ERF = (2/SQRT(Pi)) * (XERF - ((XERF**3)/3) + ERF_LSP)
!LL   BUT I have a provisional ERF COMDECK if you are interested. ACB.
C*----------------------------------------------------------------------
!
      TERF = XERF * XERF
!
      IF (ABS(XERF) .LE. 1.161) THEN
!------Hardwired version of small x approximation-----------------------
        ERF_LSP = XERF * TERF * TERF * 0.1
     &   * (1. - ((5. / 21.) * TERF * (1. - ((7. / 36.) * TERF
     &   * (1. - ((9. / 55.) * TERF * (1. - ((11. / 78.) * TERF ))))
     &     ))))
!
      ELSEIF (ABS(XERF) .GT. 1.161  .AND.  ABS(XERF) .LT. 6.) THEN
!------Continued fraction approximation at large xerf-------------------
        ERF_LSP = SIGN((SQRT_PI * 0.5),XERF) - (XERF * (1. - (TERF / 3.)
     &   + (EXP(- TERF) / (2. * (TERF + 0.5 - (0.5 / (TERF + 2.5
     &   - (3. / (TERF + 4.5 - (7.5 / (TERF + 6.5)))))))) )))
!
      ELSE
      ERF_LSP = SIGN((SQRT_PI * 0.5),XERF) + (XERF * ((TERF / 3.) - 1.))
!
      ENDIF
C*----------------------------------------------------------------------
          NERF_SUM = ERF_LSP
!
!         Note 1st 2 conditions = .FALSE. if Bs = 0.
!
          IF((- LSCBS(I)) .LT. LSCQC(I) .AND. LSCQC(I) .LE. 0.0) THEN
            RCL=RCL + (NORCT * NERF_SUM)
!   UMDP26 Eq(P26.53) ... -bs < Qc <= 0.
          ELSEIF(0.0 .LT. LSCQC(I)  .AND.  LSCQC(I) .LT. LSCBS(I)) THEN
            XERF = (LSCQC(I) * RHO / CW)
C*LL  IN-LINE FUNCTION ERF_LSP------------------------------------------
!LL
!LL     Purpose:
!LL     Vectorizable calculation of (x**3 / 3) - x + (PI**0.5 * ERF(x)/2
!LL     needed for LSPFOR routine evaluation of autoconversion rate.
!LL     Need to declare REAL ERF_LSP, XERF, TERF and to calculate
!LL     SQRT_PI = SQRT(PI) in calling routine.
!LL
!LL  Author: A Bushell
!LL
!LL  Model            Modification history from model version 4.0:
!LL version  date
!LL
!LL  Documentation:
!LL For abs(XERF) <= 1.161, use the lowest order terms of the integrated
!LL                         exp(-XERF**2) Maclaurin expansion. For a
!LL                         mantissa precision error 1E-13, this will be
!LL                         EXACT for abs(XERF) < 0.03 with max error at
!LL                         1.161 of about 0.01%.
!LL For abs(XERF) >  1.161, Press, Teukolsky, Vettering & Flannery
!LL                         NUMERICAL RECIPES IN FORTRAN eq 6.2.6 gives
!LL                         a continued fraction representation. The 1st
!LL                         4 terms give ERF to high accuracy for
!LL                         abs(XERF) > 2 & max error < 0.02% at 1.161.
!LL
!LL NOTE: ERF = (2/SQRT(Pi)) * (XERF - ((XERF**3)/3) + ERF_LSP)
!LL   BUT I have a provisional ERF COMDECK if you are interested. ACB.
C*----------------------------------------------------------------------
!
      TERF = XERF * XERF
!
      IF (ABS(XERF) .LE. 1.161) THEN
!------Hardwired version of small x approximation-----------------------
        ERF_LSP = XERF * TERF * TERF * 0.1
     &   * (1. - ((5. / 21.) * TERF * (1. - ((7. / 36.) * TERF
     &   * (1. - ((9. / 55.) * TERF * (1. - ((11. / 78.) * TERF ))))
     &     ))))
!
      ELSEIF (ABS(XERF) .GT. 1.161  .AND.  ABS(XERF) .LT. 6.) THEN
!------Continued fraction approximation at large xerf-------------------
        ERF_LSP = SIGN((SQRT_PI * 0.5),XERF) - (XERF * (1. - (TERF / 3.)
     &   + (EXP(- TERF) / (2. * (TERF + 0.5 - (0.5 / (TERF + 2.5
     &   - (3. / (TERF + 4.5 - (7.5 / (TERF + 6.5)))))))) )))
!
      ELSE
      ERF_LSP = SIGN((SQRT_PI * 0.5),XERF) + (XERF * ((TERF / 3.) - 1.))
!
      ENDIF
C*----------------------------------------------------------------------
            RCL=RCL + (NORCT * (NERF_SUM - (2. * ERF_LSP)) )
!   UMDP26 Eq(P26.53) ... 0 < Qc < bs.
          ELSEIF(LSCBS(I) .LE. LSCQC(I)) THEN
            XERF = (RHO * QDIF / CW)
C*LL  IN-LINE FUNCTION ERF_LSP------------------------------------------
!LL
!LL     Purpose:
!LL     Vectorizable calculation of (x**3 / 3) - x + (PI**0.5 * ERF(x)/2
!LL     needed for LSPFOR routine evaluation of autoconversion rate.
!LL     Need to declare REAL ERF_LSP, XERF, TERF and to calculate
!LL     SQRT_PI = SQRT(PI) in calling routine.
!LL
!LL  Author: A Bushell
!LL
!LL  Model            Modification history from model version 4.0:
!LL version  date
!LL
!LL  Documentation:
!LL For abs(XERF) <= 1.161, use the lowest order terms of the integrated
!LL                         exp(-XERF**2) Maclaurin expansion. For a
!LL                         mantissa precision error 1E-13, this will be
!LL                         EXACT for abs(XERF) < 0.03 with max error at
!LL                         1.161 of about 0.01%.
!LL For abs(XERF) >  1.161, Press, Teukolsky, Vettering & Flannery
!LL                         NUMERICAL RECIPES IN FORTRAN eq 6.2.6 gives
!LL                         a continued fraction representation. The 1st
!LL                         4 terms give ERF to high accuracy for
!LL                         abs(XERF) > 2 & max error < 0.02% at 1.161.
!LL
!LL NOTE: ERF = (2/SQRT(Pi)) * (XERF - ((XERF**3)/3) + ERF_LSP)
!LL   BUT I have a provisional ERF COMDECK if you are interested. ACB.
C*----------------------------------------------------------------------
!
      TERF = XERF * XERF
!
      IF (ABS(XERF) .LE. 1.161) THEN
!------Hardwired version of small x approximation-----------------------
        ERF_LSP = XERF * TERF * TERF * 0.1
     &   * (1. - ((5. / 21.) * TERF * (1. - ((7. / 36.) * TERF
     &   * (1. - ((9. / 55.) * TERF * (1. - ((11. / 78.) * TERF ))))
     &     ))))
!
      ELSEIF (ABS(XERF) .GT. 1.161  .AND.  ABS(XERF) .LT. 6.) THEN
!------Continued fraction approximation at large xerf-------------------
        ERF_LSP = SIGN((SQRT_PI * 0.5),XERF) - (XERF * (1. - (TERF / 3.)
     &   + (EXP(- TERF) / (2. * (TERF + 0.5 - (0.5 / (TERF + 2.5
     &   - (3. / (TERF + 4.5 - (7.5 / (TERF + 6.5)))))))) )))
!
      ELSE
      ERF_LSP = SIGN((SQRT_PI * 0.5),XERF) + (XERF * ((TERF / 3.) - 1.))
!
      ENDIF
C*----------------------------------------------------------------------
            NERF_DIF = ERF_LSP
!
            XERF = (LSCQC(I) * RHO / CW)
C*LL  IN-LINE FUNCTION ERF_LSP------------------------------------------
!LL
!LL     Purpose:
!LL     Vectorizable calculation of (x**3 / 3) - x + (PI**0.5 * ERF(x)/2
!LL     needed for LSPFOR routine evaluation of autoconversion rate.
!LL     Need to declare REAL ERF_LSP, XERF, TERF and to calculate
!LL     SQRT_PI = SQRT(PI) in calling routine.
!LL
!LL  Author: A Bushell
!LL
!LL  Model            Modification history from model version 4.0:
!LL version  date
!LL
!LL  Documentation:
!LL For abs(XERF) <= 1.161, use the lowest order terms of the integrated
!LL                         exp(-XERF**2) Maclaurin expansion. For a
!LL                         mantissa precision error 1E-13, this will be
!LL                         EXACT for abs(XERF) < 0.03 with max error at
!LL                         1.161 of about 0.01%.
!LL For abs(XERF) >  1.161, Press, Teukolsky, Vettering & Flannery
!LL                         NUMERICAL RECIPES IN FORTRAN eq 6.2.6 gives
!LL                         a continued fraction representation. The 1st
!LL                         4 terms give ERF to high accuracy for
!LL                         abs(XERF) > 2 & max error < 0.02% at 1.161.
!LL
!LL NOTE: ERF = (2/SQRT(Pi)) * (XERF - ((XERF**3)/3) + ERF_LSP)
!LL   BUT I have a provisional ERF COMDECK if you are interested. ACB.
C*----------------------------------------------------------------------
!
      TERF = XERF * XERF
!
      IF (ABS(XERF) .LE. 1.161) THEN
!------Hardwired version of small x approximation-----------------------
        ERF_LSP = XERF * TERF * TERF * 0.1
     &   * (1. - ((5. / 21.) * TERF * (1. - ((7. / 36.) * TERF
     &   * (1. - ((9. / 55.) * TERF * (1. - ((11. / 78.) * TERF ))))
     &     ))))
!
      ELSEIF (ABS(XERF) .GT. 1.161  .AND.  ABS(XERF) .LT. 6.) THEN
!------Continued fraction approximation at large xerf-------------------
        ERF_LSP = SIGN((SQRT_PI * 0.5),XERF) - (XERF * (1. - (TERF / 3.)
     &   + (EXP(- TERF) / (2. * (TERF + 0.5 - (0.5 / (TERF + 2.5
     &   - (3. / (TERF + 4.5 - (7.5 / (TERF + 6.5)))))))) )))
!
      ELSE
      ERF_LSP = SIGN((SQRT_PI * 0.5),XERF) + (XERF * ((TERF / 3.) - 1.))
!
      ENDIF
C*----------------------------------------------------------------------
!
            RCL=RCL + (NORCT * (NERF_SUM - (2. * ERF_LSP) + NERF_DIF) )
!   UMDP26 Eq(P26.53) ... bs <= Qc.
          ENDIF ! QC conditions for rain calculation.
!
!-----------------------------------------------------------------------
!L 3. Calculate ice particle fall-out speed VFALL (eq P26.33), divided
!L    by layer depth dz (metres), again assuming moisture distribution
!L    across gridbox, with local mean LSCQC and fluctuation LSCBS.
!-----------------------------------------------------------------------
!      Store VFALL/dz in VFRDZ. Note use of RHO/RHODZ for 1/dz.
!
          VFRDZ=VF1*((RHO/PCF)**VFPOWER)/((2.+VFPOWER)*(3.+VFPOWER))
          VFRDZ=VFRDZ*RHO/(RHODZ(I)*QC*LSCBS(I)*LSCBS(I))
!      We have already said IF(LSCQC(I) .GT. (-1*LSCBS(I)) ), so
          IF(LSCQC(I) .LE. 0.0) THEN
            VFRDZ=VFRDZ*(QSUM**(3.+VFPOWER))
!   UMDP26 Eq(P26.55) ... -bs < Qc <= 0.
          ELSEIF(0.0 .LT. LSCQC(I) .AND. LSCQC(I) .LE. LSCBS(I)) THEN
        VFRDZ=VFRDZ*((QSUM**(3.+VFPOWER))-(2.*(LSCQC(I)**(3.+VFPOWER))))
!   UMDP26 Eq(P26.55) ... 0 < Qc <= bs.
          ELSE ! Lscbs(I) .LT. lscqc(I)
        VFRDZ=VFRDZ*((QSUM**(3.+VFPOWER))-(2.*(LSCQC(I)**(3.+VFPOWER)))
     &      + ((LSCQC(I)-LSCBS(I))**(3.+VFPOWER)))
!   UMDP26 Eq(P26.55) ... bs < Qc.
          ENDIF ! QC conditions for ice calculation.
        ENDIF ! Cloud water content > 0
!
!-----------------------------------------------------------------------
!L 4. Calculate fraction of snow falling like water and adjust snow flux
!L    to that of snow falling as ice.
!-----------------------------------------------------------------------
        DELTA_SNOW=F_DELTA_SNOW(I)*SNOW(I)
        SNOW(I)=(1.-F_DELTA_SNOW(I))*SNOW(I)
!
!-----------------------------------------------------------------------
!L 5. Calculate ice water content and frozen water precipitation flux.
!L    Distinguish short timestep and long timestep icefall cases on
!L    basis of VFALL(I) from layer above (at this stage of routine).
!-----------------------------------------------------------------------
!
        IF ((VFALL(I) * TIMESTEP) .LE. (RHODZ(I) / RHO)) THEN
!-----------------------------------------------------------------------
!       SHORT Timestep case (also no precip falling in from above)
!-----------------------------------------------------------------------
!
          VFALL(I) = VFRDZ * RHODZ(I) / RHO
          IF ((VFALL(I) * TIMESTEP) .LE. (RHODZ(I) / RHO)) THEN
            RCF = QCF(I) * RHO * VFALL(I)
          ELSE
!         Cannot allow more ice to leave than was already there---------
            RCF = QCF(I) * RHODZ(I) / TIMESTEP
          ENDIF
          QCF(I) = QCF(I) + (TIMESTEP * (SNOW(I) - RCF) / RHODZ(I))
!
        ELSE
!-----------------------------------------------------------------------
!       LONG Timestep case
!-----------------------------------------------------------------------
!
          QCF_U = SNOW(I) / (RHO * VFALL(I))
!
          RCF = SNOW(I) - (RHODZ(I) * (QCF_U - QCF(I)) / TIMESTEP)
!         IF (RCF .LT. 0.)  THEN  ! Cannot, due to long TS condition
          QCF(I) = QCF_U
          IF (VFALL(I) .LT. (VFRDZ * RHODZ(I) / RHO)) VFALL(I) = VFRDZ *
     &                                                    RHODZ(I) / RHO
        ENDIF  ! VFALL * Timestep =< Rhodz / Rho
!
!-----------------------------------------------------------------------
!L 6. Update cloud water components (QCL, DELTA) as per eqn P26.36.
!     Note use of TIMESTEP to integrate changes over the timestep.
!-----------------------------------------------------------------------
!
        QCL(I)=QCL(I)/(1.0+TIMESTEP*RCL)
        DELTA=DELTA/(1.0+TIMESTEP*RCL)
!
!-----------------------------------------------------------------------
!L 7. Increase rates of rain- and snowfall (in kg per sq m per sec), as
!L    per eqs P26.38, 39.
!-----------------------------------------------------------------------
!
        RAIN(I)=RAIN(I) + RHODZ(I)*RCL*QCL(I)
        DELTA_SNOW=DELTA_SNOW+RHODZ(I)*RCL*DELTA
        SNOW(I) = RCF + DELTA_SNOW
!
!-----------------------------------------------------------------------
!L 8. Remember fraction of SNOW from DELTA so that it can be allowed
!L    not to fall into the next layer.
!L    Update QCF.
!-----------------------------------------------------------------------
!
        IF (SNOW(I).GT.0.0) THEN
          F_DELTA_SNOW(I)=DELTA_SNOW/SNOW(I)
        ELSE
          F_DELTA_SNOW(I)=0.0
        ENDIF
!
        QCF(I)=QCF(I) + DELTA
!
      END DO ! Loop over points
!
      RETURN
      END
