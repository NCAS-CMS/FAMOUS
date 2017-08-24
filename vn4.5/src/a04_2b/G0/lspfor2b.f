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
C*LL  SUBROUTINE LSP_FORM-----------------------------------------------
CLL
CLL  Purpose: Form or augment precipitation at the expense of cloud
CLL           water in one model layer.
CLL
CLL C.Senior    <- programmer of some or all of previous code or changes
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  Date
CLL 3.1      23/02/93 Fully divergent form of ice fallout including:
CLL                   1) fully implicit finite difference scheme
CLL                   2) total cloud water used in determining RCL and
CLL                      RCF instead of QL and QF used previously
CLL                   3) separate precipitation rates for ice (RCF) and
CLL                      liquid water (RCL) in mixed phase region
CLL                   4) minimum cloud for precipitation, CFMIN set
CLL                      equal to zero, in comdeck C_LSPFRM.
CLL                                             Ruth Carnell 26/02/93
CLL
!LL 4.0      05/10/95 Modified to CALL LSP_FOCWWIL not *CALL FOCWWIL.
!LL                                           Andrew Bushell 05/10/95
!LL  4.2    Oct. 96   T3E migration: *DEF CRAY removed
!LL                            S.J.Swarbrick
!LL
CLL  Programming standard: Unified Model Documentation Paper No 4,
CLL                        Version 1, dated 12/9/89.
CLL
CLL  Logical component covered: Part of P26.
CLL
CLL  System task:
CLL
CLL  Documentation: Unified Model Documentation Paper No 26.
C*
C*L  Arguments:---------------------------------------------------------
      SUBROUTINE LSP_FORM(
     + CF,P,Q,RHODZ,T,TIMESTEP,POINTS,QCF,QCL,RAIN,SNOW
     +,F_DELTA_SNOW,BLAND,CW_SEA,CW_LAND
     +)
      IMPLICIT NONE
      INTEGER         ! Input integer scalar :-
     + POINTS         ! IN Number of points to be processed.
      REAL            ! Input real arrays :-
     + CF(POINTS)     ! IN Cloud fraction in this layer.
     +,P(POINTS)      ! IN Air pressure at this level (Pa).
     +,Q(POINTS)      ! IN Sp humidity at this level (kg wat per kg air)
     +,RHODZ(POINTS)  ! IN Air mass p.u.a. in this layer (kg per sq m).
     +,T(POINTS)      ! IN Temperature at this level (K).
      REAL            ! Input real scalar :-
     + TIMESTEP       ! IN Timestep (s).
     +,CW_SEA         ! IN threshold cloud liquid water content
C                     !    over sea for conversion to ppn
C                     !    (kg water per m**3)
     +,CW_LAND        ! IN threshold cloud liquid water content
C                     !    over land for conversion to ppn
C                     !    (kg water per m**3)
      REAL            ! Updated real arrays :-
     + QCF(POINTS)    ! INOUT Cloud ice (kg water per kg air).
     +,QCL(POINTS)    ! INOUT Cloud liquid water (kg water per kg air).
     +,RAIN(POINTS)   ! INOUT Rate of rainfall entering this layer from
C                     !       above (on i/p), or leaving it (on o/p)
C                     !       (kg per sq m per s).
     +,SNOW(POINTS)   ! INOUT Rate of snowfall entering this layer from
C                     !       above (on i/p), or leaving it (on o/p)
C*                    !       (kg per sq m per s).
     +,F_DELTA_SNOW(POINTS) ! INOUT snow fraction from ice falling as
C                     !             water
      LOGICAL BLAND(POINTS)  ! IN Land/sea mask
C*L  Workspace usage----------------------------------------------------
!  1 real, 0 logical and 0 integer blocks are required, as follows :-
      REAL
     & FL(POINTS)     ! Fraction of total cloud water which is liquid.
C*L   External subprogram called :-
      EXTERNAL LSP_FOCWWIL
C*
C-----------------------------------------------------------------------
C  Common, then local, physical constants.
C-----------------------------------------------------------------------
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

C*LL--------------------------------------------------------------------
CLL  Tunable constants local to routine LSP_FORM, called by P26
CLL  LS_PPN.
CLL---------------------------------------------------------------------
      REAL CA,PCF,CFMIN,CT,FPL,VF1,VFPOWER
      PARAMETER (
     +CA=1.0,           ! Accretion constant (kg per m**2).
     +PCF=1.01086E-3,   ! Tunable in Heymsfield formula (kg/m**3).
     +CFMIN=0.0,        ! Minimum cloud amount for ppn calculations.
     +CT=1.0E-4,        ! Reciprocal time constant for conversion of
C                       ! cloud liquid water to ppn (Hertz).
     +FPL=0.0,          ! Fraction of snow which falls through the
C                       ! layer, i.e. is treated as per rain.
C                       ! NB: This will eventually be a function of
C                       !     some model variables - parameter used
C                       !     only for temporary convenience.
     +VF1=1.0,          ! Tunable speed in Heymsfield formula (m/s).
     +VFPOWER=0.17      ! Local parameter, used in Heymsfield formula.
     +)
C*

C
C
C-----------------------------------------------------------------------
C  Define local scalars.
C-----------------------------------------------------------------------
C  (a) Reals effectively expanded to workspace by the Cray (using
C      vector registers).
      REAL            ! Real workspace.  At end of DO loop, contains :-
     & QC             ! Total cloud water (liquid + ice).
     +,RCL            ! Rate of liquid water exchange
     +,RCF            ! Rate of ice exchange
     +,DELTA          ! 'ice' to be treated as water
     +,DELTA_SNOW     ! Snow formed by ice falling as water
     +,RHO            ! Density of air in the layer.
     +,TEMP           ! Scalar temporary store.
     +,VFRDZ          ! Ice-particle fallout speed (m/s), divided by
C                     ! layer depth (m).
     +,VTEMP          ! Virtual temperature as at start of loop.
     +,CW               ! Parameter for efficient conversion of
C                       ! liquid cloud water to ppn
C  (b) Others.
      INTEGER I       ! Loop counter (horizontal field index).
     +,NV             ! Index counter for 'when' routine.
!
!-----------------------------------------------------------------------
!L 0. Calculate fraction of cloud water which is liquid (FL),
!L    according to equation P26.29.
!-----------------------------------------------------------------------
      CALL LSP_FOCWWIL(T, POINTS, FL)
!
      DO 1 I=1,POINTS
C-----------------------------------------------------------------------
!  1. Store total cloud water in QC (see eq P26.30) and ice falling as
!     falling as water in DELTA. Also update QCF.
C-----------------------------------------------------------------------
C
        QC=QCL(I)+QCF(I)
        DELTA = (FL(I) * QC) - QCL(I)
        QCF(I) = (1. - FL(I)) * QC
C
C-----------------------------------------------------------------------
CL 2. Calculate fractional rate of conversion of liquid cloud water to
CL    precipitation.
C     Store in RCL.
C-----------------------------------------------------------------------
C
C  (a) Calculate density of air, RHO, via virtual temperature.
C
        VTEMP=T(I)*(1.+C_VIRTUAL*Q(I)-QC)      ! Virtual temperature
        RHO=P(I)/(R*VTEMP)
C
C  (b) Calculate P/CA.  Store in RCL.
C
        RCL=(RAIN(I)+SNOW(I))/CA
C
C  (c) Calculate "term involving CT" - assume zero if cloud fraction
C      so small as to make the exp very close to 1.
C      Note that QC is now used in this expression.
C      (Also initialise variable calculated in sect 3 below, q.v.)
C
        TEMP=0.0
        VFRDZ=0.0
C
      IF (BLAND(I)) THEN
           CW=CW_LAND
      ELSE
           CW=CW_SEA
      ENDIF
C
        IF(CF(I).GT.CFMIN)THEN
          TEMP=RHO*QC/(CF(I)*CW)
          RCL=CT*(1.0-EXP(-TEMP*TEMP))
     +                    + RCL
C
C-----------------------------------------------------------------------
CL 3. Calculate ice particle fall-out speed VF (eq P26.33), divided
CL    by layer depth dz (metres).
C-----------------------------------------------------------------------
C
C      Store VF/dz in VFRDZ. Note (1) QC is now used in this
C      expression, (2) use of RHO/RHODZ for 1/dz.
C
          VFRDZ=VF1*
     +      ((RHO*QC/(CF(I)*PCF))**VFPOWER)
          VFRDZ=VFRDZ*RHO/RHODZ(I)
        ENDIF
C
C-----------------------------------------------------------------------
CL 4. Calculate fractional rate of conversion of cloud ice to
CL    precipitation.
C-----------------------------------------------------------------------
C
        RCF=VFRDZ
C
C-----------------------------------------------------------------------
CL 5. Update cloud water components (QCL, QCF) as per eqs P26.36, 37.
C     Note use of TIMESTEP to integrate changes over the timestep.
CL    Also update DELTA_SNOW, SNOW and DELTA
C-----------------------------------------------------------------------
C
        DELTA_SNOW=F_DELTA_SNOW(I)*SNOW(I)
        SNOW(I)=(1.-F_DELTA_SNOW(I))*SNOW(I)
C
        QCL(I)=QCL(I)/(1.0+TIMESTEP*RCL)
        QCF(I)=QCF(I)/(1.0+TIMESTEP*RCF) +
     +   (1.0-FPL)*TIMESTEP*SNOW(I)/(RHODZ(I)*(1.+TIMESTEP*RCF))
        DELTA=DELTA/(1.0+TIMESTEP*RCL)
C
C-----------------------------------------------------------------------
CL 6. Increase rates of rain- and snowfall (in kg per sq m per sec), as
CL    per eqs P26.38, 39.
C-----------------------------------------------------------------------
C
        RAIN(I)=RAIN(I) + RHODZ(I)*RCL*QCL(I)
        SNOW(I)=SNOW(I)*FPL + RHODZ(I)*RCF*QCF(I)
        DELTA_SNOW=DELTA_SNOW+RHODZ(I)*RCL*DELTA
        SNOW(I)=SNOW(I)+DELTA_SNOW
C
C-----------------------------------------------------------------------
CL 7. Remember fraction of SNOW from DELTA so that it can be allowed
CL    not to fall into the next layer.
CL    Update QCF.
C-----------------------------------------------------------------------
C
        IF (SNOW(I).GT.0.0) THEN
        F_DELTA_SNOW(I)=DELTA_SNOW/SNOW(I)
        ELSE
        F_DELTA_SNOW(I)=0.0
        ENDIF
C
        QCF(I)=QCF(I) + DELTA
C
    1 CONTINUE
      RETURN
      END
