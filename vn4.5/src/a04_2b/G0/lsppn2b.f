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
C*LL  SUBROUTINES LS_PPN LS_PPNC and LSP_FORM---------------------------
CLL  Purpose:
CLL          LS_PPN and LS_PPNC:
CLL           Calculate large-scale (dynamical) precipitation.
CLL           LS_PPNC is the gather/scatter routine which then
CLL           calls LSP_EVAP,LSP_FRMT,LSP_FORM.
CLL           Treatment of evaporation made implicit. C Wilson 18/09/90.
CLL  Note: in all cases, level counters (incl subscripts) run from 1
CLL        (lowest model layer) to Q_LEVELS (topmost "wet" model
CLL        layer) - it is assumed that the bottom Q_LEVELS layers are
CLL        the "wet" layers.
CLL
CLL  Put through fpp on Cray.  Activate *IF definition CRAY if running
CLL  on the Cray.  Function FOCWWIL is now a COMDECK
CLL                (This function is called by LSP_FORM.)
CLL
CLL  This routine is suitable for single-column use.
CLL
CLL C.Wilson    <- programmer of some or all of previous code or changes
CLL C.Senior    <- programmer of some or all of previous code or changes
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  Date
CLL 3.1      23/02/93 LS_PPN and LS_PPNC
CLL                   Inclusion of F_DELTA_SNOW (fraction of snow from
CLL                   ice falling as water) for use in LSP_FORM with
CLL                   fully divergent ice fallout.
CLL 3.4      15/08/94 Include layer rain and snow deltas for aerosol.
CLL
CLL                                             Ruth Carnell 26/02/93
CLL 3.4    21/09/94 Correct gather/scatter index code (only accessed
CLL                 when *DEF CRAY is not active). Note that *DEF CRAY
CLL                 is used inappropriately in this routine to
CLL                 separate code for the full model and the single
CLL                 column model.                         R. Rawlins
CLL   4.2    Oct. 96  T3E migration: *DEF CRAY removed, whenimd removed.
CLL                                    S.J.Swarbrick
CLL 4.5  01/05/98  Restrict murk aerosol calculations to aerosol
CLL                levels=boundary levels. P.Clark
CLL
CLL  Programming standard: Unified Model Documentation Paper No 4,
CLL                        Version 2, dated 18/1/90.
CLL
CLL  Logical component covered: P26.
CLL
CLL  Project task:
CLL
CLL  Documentation: UM Documentation Paper 26.
CLL
C*L  Arguments:---------------------------------------------------------
      SUBROUTINE LS_PPN(
     +AK,BK,CF,DELTA_AK,DELTA_BK,PSTAR,TIMESTEP,BLAND,
     +CW_SEA,CW_LAND,Q_LEVELS,PFIELD,
     & POINTS,K1STPT,A_LEVELS,Q,QCF,QCL,T,
     & AEROSOL,L_MURK,LSRAIN,LSSNOW,            
     & ERROR
     +)
      IMPLICIT NONE
      INTEGER
     + Q_LEVELS ! IN Number of "wet" levels in the model.
     &,A_LEVELS ! IN Number of aerosol levels
     +,PFIELD   ! IN Number of gridpoints in one field (at one level).
     +,POINTS   ! IN Number of gridpoints being processed.
     +,K1STPT   ! IN First gridpoint processed within complete field.
      REAL
     + CF(PFIELD,Q_LEVELS) ! IN Cloud fraction.
     +,PSTAR(PFIELD)       ! IN Surface pressure (Pa).
     +,AK(Q_LEVELS)        ! IN Hybrid co-ordinate for centre of layer.
     +,BK(Q_LEVELS)        ! IN Hybrid co-ordinate for centre of layer.
     +,DELTA_AK(Q_LEVELS)  ! IN Change of hybrid co-ord across layer.
C                          !    (Upper minus lower).
     +,DELTA_BK(Q_LEVELS)  ! IN Change of hybrid co-ord across layer.
C                          !    (Upper minus lower).
      REAL TIMESTEP        ! IN Timestep (sec).
      REAL CW_SEA,          ! IN threshold cloud liquid water content
C                           !    over sea for conversion to ppn
C                           !    (kg water per m**3)
     +     CW_LAND          ! IN threshold cloud liquid water content
C                           !    over land for conversion to ppn
C                           !    (kg water per m**3)
      LOGICAL BLAND(PFIELD) ! IN Land/sea mask
     & , L_MURK                  ! IN : Aerosol needs scavenging.
      REAL
     + Q(PFIELD,Q_LEVELS)   ! INOUT Specific humidity (kg water/kg air).
     +,QCF(PFIELD,Q_LEVELS) ! INOUT Cloud ice (kg per kg air).
     +,QCL(PFIELD,Q_LEVELS) ! INOUT Cloud liquid water (kg per kg air).
     +,T(PFIELD,Q_LEVELS)   ! INOUT Temperature (K).
     &,AEROSOL(PFIELD,A_LEVELS)  ! INOUT Aerosol
      REAL
     + LSRAIN(PFIELD) ! OUT Surface rainfall rate (kg per sq m per s).
     +,LSSNOW(PFIELD) ! OUT Surface snowfall rate (kg per sq m per s).
      INTEGER
     + ERROR          ! OUT Return code - 0 if OK,
C                     !                   1 if bad arguments.
C*L  Workspace usage ---------------------------------------------------
C  0 real,1 logical and 1 integer blocks are required, as follows :-
      LOGICAL
     + H(PFIELD)      ! Used as "logical" in compression.
     & ,L_SCAVENGE  ! scavenge aerosol on level.
      INTEGER
     + IX(PFIELD)     ! Index for compress/expand.
      REAL F_DELTA_SNOW(PFIELD) ! snow fraction from ice falling
C                               ! as water
C  External subroutines called -----------------------------------------
      EXTERNAL LS_PPNC
C*----------------------------------------------------------------------
C  Physical constants -------------------------------------------------
      REAL CFMIN
      PARAMETER (
     + CFMIN=1.0E-3        ! Used for LS_PPNC  compress.
     +)
C  Define local variables ----------------------------------------------
      INTEGER I,K     ! Loop counters: I - horizontal field index;
C                     !                K - vertical level index.
     +,N              ! "nval" for WHEN routine.
!
      ERROR=0
      IF((K1STPT+POINTS-1).GT.PFIELD)THEN
        ERROR=1
        GOTO20
      ENDIF
C-----------------------------------------------------------------------
CL Internal structure.
CL 1. Initialise rain and snow to zero.
C-----------------------------------------------------------------------
      DO 1 I=K1STPT,K1STPT+POINTS-1
        LSRAIN(I)=0.0
        LSSNOW(I)=0.0
        F_DELTA_SNOW(I)=0.0
    1 CONTINUE
C
C-----------------------------------------------------------------------
CL 2. Loop round levels from top down (counting bottom level as level 1,
CL    as is standard in the Unified model).
C-----------------------------------------------------------------------
C
      DO 5 K=Q_LEVELS,1,-1
C-----------------------------------------------------------------------
CL 2.5 Form INDEX IX to gather/scatter variables in LS_PPNC
C-----------------------------------------------------------------------
C
C  Set index where cloud fraction > CFMIN or where non-zero pptn
C  Note: whenimd is functionally equivalent to WHENILE (but autotasks).
C
        N=0
        DO 10 I=K1STPT,K1STPT+POINTS-1
          IF (CF(I,K).GT.CFMIN .OR. (LSRAIN(I)+LSSNOW(I)).GT.0.0) THEN
            N=N+1
            IX(N)=I-K1STPT+1
          ENDIF
   10   CONTINUE
C
        L_SCAVENGE = L_MURK .AND. (K.LE.A_LEVELS)
        IF(N.GT.0)THEN

          CALL LS_PPNC(IX,N,TIMESTEP,POINTS,PSTAR(K1STPT),
     &                 LSRAIN(K1STPT),LSSNOW(K1STPT),CF(K1STPT,K),
     &                 QCF(K1STPT,K),QCL(K1STPT,K),T(K1STPT,K),
     &                 AEROSOL(K1STPT,MIN(K,A_LEVELS)),L_SCAVENGE,
     &                Q(K1STPT,K),AK(K),BK(K),DELTA_AK(K),DELTA_BK(K),
     &                F_DELTA_SNOW(K1STPT),BLAND(K1STPT),CW_SEA,CW_LAND)
        ENDIF
C
    5 CONTINUE
   20 CONTINUE
      RETURN
      END
C*LL  SUBROUTINE LS_PPNC------------------------------------------------
C*L  Arguments:---------------------------------------------------------
      SUBROUTINE LS_PPNC(
     & IX,N,TIMESTEP,POINTS,PSTAR,LSRAIN,LSSNOW,
     & CF,QCF,QCL,T,AEROSOL,L_MURK,Q,
     + AK,BK,DELTA_AK,DELTA_BK
     +,F_DELTA_SNOW,BLAND,CW_SEA,CW_LAND
     +)
      IMPLICIT NONE
      INTEGER
     + N        ! IN Number of points where pptn non-zero from above
     +,IX(N)    ! IN gather/scatter index
     +,POINTS   ! IN Number of gridpoints being processed.
     +          !    or where CF>CFMIN
      REAL
     + PSTAR(POINTS)       ! IN Surface pressure (Pa).
     +,CF(POINTS)          ! IN Cloud fraction.
     +,AK                  ! IN Hybrid co-ordinate for centre of layer.
     +,BK                  ! IN Hybrid co-ordinate for centre of layer.
     +,DELTA_AK            ! IN Change of hybrid co-ord across layer.
C                          !    (Upper minus lower).
     +,DELTA_BK            ! IN Change of hybrid co-ord across layer.
C                          !    (Upper minus lower).
      REAL CW_SEA,          ! IN threshold cloud liquid water content
C                           !    over sea for conversion to ppn
C                           !    (kg water per m**3)
     +     CW_LAND          ! IN threshold cloud liquid water content
C                           !    over land for conversion to ppn
C                           !    (kg water per m**3)
      LOGICAL BLAND(POINTS)   ! IN Land/sea mask
     & , L_MURK                  ! IN : Aerosol needs scavenging.
      REAL TIMESTEP        ! IN Timestep (sec).
      REAL
     + Q(POINTS)            ! INOUT Specific humidity (kg water/kg air).
     +,QCF(POINTS)          ! INOUT Cloud ice (kg per kg air).
     +,QCL(POINTS)          ! INOUT Cloud liquid water (kg per kg air).
     +,T(POINTS)            ! INOUT Temperature (K).
     &,AEROSOL(POINTS)      ! INOUT Aerosol (K).
      REAL
     + LSRAIN(POINTS) !INOUT Surface rainfall rate (kg per sq m per s).
     +,LSSNOW(POINTS) !INOUT Surface snowfall rate (kg per sq m per s).
     +,F_DELTA_SNOW(POINTS) !INOUT snow fraction from ice falling as
C                           !water.
C*L  Workspace usage ---------------------------------------------------
C  10 real,0 logical and 0 integer blocks are required, as follows :-
      REAL
     + PSTAR_C(N)          ! gathered Surface pressure (Pa).
     +,CF_C(N)             ! gathered Cloud fraction.
     +,Q_C(N)            ! gathered Specific humidity (kg water/kg air).
     +,QCF_C(N)          ! gathered Cloud ice (kg per kg air).
     +,QCL_C(N)          ! gathered Cloud liquid water (kg per kg air).
     +,T_C(N)            ! gathered Temperature (K).
     &,AERO_C(N)      !  Aerosol
     +,LSRAIN_C(N) !gathered Surface rainfall rate (kg per sq m per s).
     +,LSSNOW_C(N) !gathered Surface snowfall rate (kg per sq m per s).
     +,F_DELTA_SNOW_C(N) !gathered fraction of snow
      REAL
     + RHODZ(N)       ! WORK Used for air mass p.u.a. in successive
C                     !      layers.
     +,P(N)           ! WORK Used for pressure at successive levels.
      LOGICAL BLAND_C(N)          ! gathered land/sea mask
C  External subroutines called -----------------------------------------
      EXTERNAL LSP_EVAP,LSP_FORM,LSP_FRMT,LSP_SCAV
C*----------------------------------------------------------------------
C  Physical constants -------------------------------------------------
C*L------------------COMDECK C_G----------------------------------------
C G IS MEAN ACCEL DUE TO GRAVITY AT EARTH'S SURFACE
      REAL G

      PARAMETER(G=9.80665)
C*----------------------------------------------------------------------
      REAL P1UPONG
      PARAMETER (
     + P1UPONG=1./G        ! One upon g (sq seconds per m).
     +)
C  Define local variables ----------------------------------------------
      INTEGER I       ! Loop counters: I - horizontal field index;
C
C-----------------------------------------------------------------------
CL Internal structure.
CL 1. gather variables using index
C-----------------------------------------------------------------------
      DO 1 I=1,N
        LSRAIN_C(I)=LSRAIN(IX(I))
        LSSNOW_C(I)=LSSNOW(IX(I))
        PSTAR_C(I) =PSTAR(IX(I))
        BLAND_C(I) =BLAND(IX(I))
        CF_C(I)=CF(IX(I))
        QCF_C(I)=QCF(IX(I))
        QCL_C(I)=QCL(IX(I))
        Q_C(I)=Q(IX(I))
        T_C(I)=T(IX(I))
        IF (L_MURK) AERO_C(I)=AEROSOL(IX(I))
        F_DELTA_SNOW_C(I)=F_DELTA_SNOW(IX(I))
    1 CONTINUE
C
C-----------------------------------------------------------------------
CL 2   Calculate pressure at current level, and air mass p.u.a. of
CL     current layer.
C      (Negative in RHODZ formula takes account of sign of DELTAs.)
C-----------------------------------------------------------------------
        DO 2 I=1,N
          P(I)=AK+PSTAR_C(I)*BK
          RHODZ(I)=-P1UPONG*(DELTA_AK+PSTAR_C(I)*DELTA_BK)
    2   CONTINUE
C
C-----------------------------------------------------------------------
CL 3 If there is precipitation falling from above, then :-
C-----------------------------------------------------------------------
C
CL 3.1 Evaporate from precipitation, and calculate the effect on the
CL     temperature and specific humidity.  Do this by calling LSP_EVAP.
C
          CALL LSP_EVAP(P,RHODZ,TIMESTEP,N,
     +      Q_C,LSRAIN_C,LSSNOW_C,T_C )
C
C
CL 3.2 Change phase of precipitation where necessary, and calculate
CL     the effect on the temperature and specific humidity.  Also set
CL     rain/snow indicator (after any incrementing of the temperature).
CL     All this is done by calling LSP_FRMT.
C
          CALL LSP_FRMT(RHODZ,TIMESTEP,N,QCF_C,
     +      QCL_C,LSRAIN_C,LSSNOW_C,T_C )
C-----------------------------------------------------------------------
CL 3.3 Form (or augment) precipitation at the expense of cloud water.
CL     Do this by calling LSP_FORM.
C-----------------------------------------------------------------------
C
          CALL LSP_FORM(CF_C,P,Q_C,RHODZ,T_C,TIMESTEP,N,QCF_C,QCL_C,
     +                 LSRAIN_C,LSSNOW_C,F_DELTA_SNOW_C,BLAND_C,
     +                 CW_SEA,CW_LAND)
C
        IF (L_MURK) THEN
C Lose aerosol by scavenging.
C
          CALL LSP_SCAV(TIMESTEP,N,LSRAIN_C,LSSNOW_C,AERO_C)
C
        END IF
C-----------------------------------------------------------------------
CL 4  Scatter back arrays which will have been changed.
CL
C-----------------------------------------------------------------------
C
CDIR$ IVDEP
          DO 4 I=1,N
            T(IX(I))=T_C(I)
            Q(IX(I))=Q_C(I)
            QCF(IX(I))=QCF_C(I)
            QCL(IX(I))=QCL_C(I)
            IF (L_MURK) AEROSOL(IX(I))=AERO_C(I)
            LSRAIN(IX(I))=LSRAIN_C(I)
            LSSNOW(IX(I))=LSSNOW_C(I)
            F_DELTA_SNOW(IX(I)) = F_DELTA_SNOW_C(I)
   4      CONTINUE
      RETURN
      END
