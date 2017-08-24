*****************************COPYRIGHT******************************
C (c) CROWN COPYRIGHT 1998, METEOROLOGICAL OFFICE, All Rights Reserved.
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
C*LL  SUBROUTINES LS_PPN and LS_PPNC------------------------------------
!LL  Purpose:
!LL          LS_PPN and LS_PPNC:
!LL           Calculate large-scale (dynamical) precipitation.
!LL           LS_PPNC is the gather/scatter routine which then
!LL           calls LSP_ICE.
!LL  Note: in all cases, level counters (incl subscripts) run from 1
!LL        (lowest model layer) to Q_LEVELS (topmost "wet" model
!LL        layer) - it is assumed that the bottom Q_LEVELS layers are
!LL        the "wet" layers.
!LL
!LL  Put through fpp on Cray.  Activate *IF definition CRAY if running
!LL  on the Cray.
!LL
!LL  Modification History from Version 4.4
!LL     Version    Date
!LL       4.5      March 98          New Deck        Damian Wilson
!LL
!LL  Programming standard: Unified Model Documentation Paper No 4,
!LL                        Version 2, dated 18/1/90.
!LL
!LL  Logical component covered: P26.
!LL
!LL  Project task:
!LL
!LL  Documentation: UM Documentation Paper 26.
!LL
C*L  Arguments:---------------------------------------------------------
      SUBROUTINE LS_PPN(
     &AK,BK,CF,DELTA_AK,DELTA_BK,PSTAR,TIMESTEP,BLAND,
     &CW_SEA,CW_LAND,
     &CLOUD_LIQ_FRAC,CLOUD_ICE_FRAC,
     &RHCRIT,
     &RHCPT, L_RHCPT,
     &Q_LEVELS,PFIELD,
     &POINTS,K1STPT,LSPICE_DIM1,LSPICE_DIM2,A_LEVELS,Q,QCF,QCL,T,
     &SO2,L_SULPC_SO2,
     &NH3,L_SULPC_NH3,               !MW
     &SO4_AIT,SO4_ACC,SO4_DIS,
     & AGED_SOOT, L_SOOT,
     &AEROSOL,L_MURK,
     &LSRAIN,LSSNOW,
     &LSRAIN3D,LSSNOW3D,
     &LSCAV_SO2,LSCAV_SO4AIT,LSCAV_SO4ACC,LSCAV_SO4DIS,
     &LSCAV_NH3,                      !MW
     & LSCAV_AGEDSOOT,
     &ERROR
     &)
      IMPLICIT NONE
      INTEGER
     & Q_LEVELS ! IN Number of "wet" levels in the model.
     &,PFIELD   ! IN Number of gridpoints in one field (at one level).
     &,POINTS   ! IN Number of gridpoints being processed.
     &,K1STPT   ! IN First gridpoint processed within complete field.
     &,LSPICE_DIM1  ! Dimension of arrays LSRAIN3D and LSSNOW3D.
     &,LSPICE_DIM2  ! Dimension of arrays LSRAIN3D and LSSNOW3D.
     &,A_LEVELS ! IN Number of aerosol levels.
      REAL
     & CF(PFIELD,Q_LEVELS) ! IN Cloud fraction.
     &,PSTAR(PFIELD)       ! IN Surface pressure (Pa).
     &,AK(Q_LEVELS)        ! IN Hybrid co-ordinate for centre of layer.
     &,BK(Q_LEVELS)        ! IN Hybrid co-ordinate for centre of layer.
     &,DELTA_AK(Q_LEVELS)  ! IN Change of hybrid co-ord across layer.
!                               (Upper minus lower).
     &,DELTA_BK(Q_LEVELS)  ! IN Change of hybrid co-ord across layer.
!                               (Upper minus lower).
     &,RHCRIT(Q_LEVELS)    ! IN Critical humidity for cloud formation.
     &,RHCPT(PFIELD,Q_LEVELS) ! IN: Crit. hum. for cloud formation
     &,CLOUD_LIQ_FRAC(PFIELD,Q_LEVELS) !IN Cloud liquid fraction.
     &,CLOUD_ICE_FRAC(PFIELD,Q_LEVELS) !IN Cloud ice fraction.
      REAL TIMESTEP        ! IN Timestep (sec).
     &    ,CW_SEA          ! IN threshold cloud liquid water content
!                               over sea for conversion to ppn
!                               (kg water per m**3)
     &    ,CW_LAND         ! IN threshold cloud liquid water content
!                               over land for conversion to ppn
!                               (kg water per m**3)
      LOGICAL BLAND(PFIELD) ! IN Land/sea mask
     &,       L_MURK        ! IN Aerosol needs scavenging.
      LOGICAL L_RHCPT  ! Indicates whether RHcrit parametrization is on.
      LOGICAL L_SULPC_SO2  !IN Sulphur Cycle on, tracers scavenged if T
     &,       L_SULPC_NH3  !IN indicates if NH3 present    !MW
     &,       L_SOOT       !IN indicates whether soot present
!
      REAL
     & Q(PFIELD,Q_LEVELS)   ! INOUT Specific humidity (kg water/kg air).
     &,QCF(PFIELD,Q_LEVELS) ! INOUT Cloud ice (kg per kg air).
     &,QCL(PFIELD,Q_LEVELS) ! INOUT Cloud liquid water (kg per kg air).
     &,T(PFIELD,Q_LEVELS)   ! INOUT Temperature (K).
     &,AEROSOL(PFIELD,A_LEVELS) ! INOUT Aerosol (K).
      REAL               !INOUT, Sulphur Cycle tracers (mmr kg/kg)
     &    SO2(PFIELD,Q_LEVELS)
     &   ,NH3(PFIELD,Q_LEVELS)                 !MW
     &   ,SO4_AIT(PFIELD,Q_LEVELS)
     &   ,SO4_ACC(PFIELD,Q_LEVELS)
     &   ,SO4_DIS(PFIELD,Q_LEVELS)
     &   ,AGED_SOOT(PFIELD,Q_LEVELS)
!
      REAL
     & LSRAIN(PFIELD) ! OUT Surface rainfall rate (kg per sq m per s).
     &,LSSNOW(PFIELD) ! OUT Surface snowfall rate (kg per sq m per s).
      REAL               ! OUT column totals of S Cycle tracers scavngd
     &    LSCAV_SO2(PFIELD)
     &   ,LSCAV_NH3(PFIELD)                      !MW
     &   ,LSCAV_SO4AIT(PFIELD)
     &   ,LSCAV_SO4ACC(PFIELD)
     &   ,LSCAV_SO4DIS(PFIELD)
     &   ,LSCAV_AGEDSOOT(PFIELD)
!
      REAL
     &    LSRAIN3D(LSPICE_DIM1,LSPICE_DIM2) ! OUT Rain rate out of
!                                             each model layer
     &   ,LSSNOW3D(LSPICE_DIM1,LSPICE_DIM2) ! OUT Snow rate out of
!                                             each model layer
!
      INTEGER
     & ERROR          ! OUT Return code - 0 if OK,
!                                         1 if bad arguments.
C*L  Workspace usage ---------------------------------------------------
!  0 real,1 logical and 2 integer blocks are required, as follows :-
      LOGICAL
     & H(PFIELD)      ! Used as "logical" in compression.
     &,L_SCAVENGE     ! scavenge aerosol on level.
      INTEGER
     & IX(PFIELD)     ! Index for compress/expand.
      REAL F_DELTA_SNOW(PFIELD) ! snow fraction from ice falling
!                                 as water
      REAL VFALL(PFIELD)        ! snow fall velocity (m per s).
! Allocate CX and CONSTP arrays
! Sets up the size of arrays for CX and CONSTP
      REAL CX(16),CONSTP(16)
!
!  External subroutines called -----------------------------------------
      EXTERNAL LS_PPNC,LSPCON
C*----------------------------------------------------------------------
!  Physical constants -------------------------------------------------
      REAL CFMIN
      PARAMETER (
     & CFMIN=1.0E-3        ! Used for LS_PPNC  compress.
     &)
!  Define local variables ----------------------------------------------
      INTEGER I,K     ! Loop counters: I - horizontal field index;
!                                      K - vertical level index.
     &,N              ! "nval" for WHEN routine.
!
      ERROR=0
      IF((K1STPT+POINTS-1).GT.PFIELD)THEN
        ERROR=1
        GOTO9999
      ENDIF
! Define CX and CONSTP values
      CALL LSPCON(CX,CONSTP)
!-----------------------------------------------------------------------
!L Internal structure.
!L 1. Initialise rain and snow to zero.
!   Initialise scavenged amounts of S Cycle tracers to 0 for full field
!-----------------------------------------------------------------------
      DO I=K1STPT,K1STPT+POINTS-1
        LSRAIN(I)=0.0
        LSSNOW(I)=0.0
        F_DELTA_SNOW(I)=0.0
        VFALL(I)=0.0
      END DO ! Loop over points
!
       DO I=1,PFIELD
        LSCAV_SO2(I)=0.0
        LSCAV_NH3(I)=0.0                        !MW
        LSCAV_SO4AIT(I)=0.0
        LSCAV_SO4ACC(I)=0.0
        LSCAV_SO4DIS(I)=0.0
        LSCAV_AGEDSOOT(I)=0.0
       END DO
!
!-----------------------------------------------------------------------
!L 2. Loop round levels from top down (counting bottom level as level 1,
!L    as is standard in the Unified model).
!-----------------------------------------------------------------------
!
      DO K=Q_LEVELS,1,-1
!-----------------------------------------------------------------------
!L 2.5 Form INDEX IX to gather/scatter variables in LS_PPNC
!-----------------------------------------------------------------------
!
!  Set index where cloud fraction > CFMIN or where non-zero pptn
!  Note: whenimd is functionally equivalent to WHENILE (but autotasks).
!
!
        N=0
        DO I=K1STPT,K1STPT+POINTS-1
          IF (CLOUD_LIQ_FRAC(I,K).GT.CFMIN
     &          .OR. (LSRAIN(I)+LSSNOW(I)).GT.0.0
     &          .OR. QCF(I,K).GT.0.0 ) THEN
            N=N+1
            IX(N)=I - K1STPT + 1
          ENDIF
        END DO ! Loop over points
!
        L_SCAVENGE = L_MURK .AND. (K.LE.A_LEVELS)
!
        IF(N.GT.0)THEN

          CALL LS_PPNC(IX,N,TIMESTEP,POINTS,PSTAR(K1STPT),
     &                 LSRAIN(K1STPT),LSSNOW(K1STPT),CF(K1STPT,K),
     &                 QCF(K1STPT,K),QCL(K1STPT,K),T(K1STPT,K),
     &           SO2(K1STPT,K),L_SULPC_SO2,
     &           NH3(K1STPT,K),L_SULPC_NH3,                     !MW
     &           SO4_AIT(K1STPT,K),SO4_ACC(K1STPT,K),SO4_DIS(K1STPT,K),
     &                 AGED_SOOT(K1STPT,K), L_SOOT,
! Aerosol is only defined on A_LEVELS not Q_LEVELS so limit index K
     &                 AEROSOL(K1STPT,MIN(K,A_LEVELS)),L_SCAVENGE,
     &           LSCAV_NH3(K1STPT),                             !MW
     &           LSCAV_SO2(K1STPT),LSCAV_SO4AIT(K1STPT),
     &           LSCAV_SO4ACC(K1STPT),LSCAV_SO4DIS(K1STPT),
     &                 LSCAV_AGEDSOOT(K1STPT),
     &                 Q(K1STPT,K),AK(K),BK(K),DELTA_AK(K),DELTA_BK(K),
     &                 F_DELTA_SNOW(K1STPT),BLAND(K1STPT),CW_SEA,
     &                 CW_LAND,
     &                CLOUD_LIQ_FRAC(K1STPT,K),CLOUD_ICE_FRAC(K1STPT,K),
     &                 RHCRIT(K),
     &                 RHCPT(K1STPT,K), L_RHCPT,
     &                 VFALL(K1STPT),CX,CONSTP)
        ENDIF
!
! Copy rainfall and snowfall rates to 3D fields for diagnostic output
!
        IF (LSPICE_DIM1 .EQ. PFIELD
     &      .AND. LSPICE_DIM2 .EQ. Q_LEVELS) THEN
! Only copy rain and snow to 3D fields if arrays are dimensionalized.
          DO I=K1STPT,K1STPT+POINTS-1
            LSRAIN3D(I,K)=LSRAIN(I)
            LSSNOW3D(I,K)=LSSNOW(I)
          END DO ! Loop over I for 3D diagnostics
        ENDIF
!
      END DO ! Loop over K
 9999 CONTINUE   ! Branch for error exit
      RETURN
      END
C*LL  SUBROUTINE LS_PPNC------------------------------------------------
C*L  Arguments:---------------------------------------------------------
      SUBROUTINE LS_PPNC(
     & IX,N,TIMESTEP,POINTS,PSTAR,LSRAIN,LSSNOW
     &,CF,QCF,QCL,T
     &,SO2,L_SULPC_SO2
     &,NH3,L_SULPC_NH3                           !MW
     &,SO4_AIT,SO4_ACC,SO4_DIS
     &,AGED_SOOT, L_SOOT
     &,AEROSOL,L_MURK
     &,LSCAV_NH3                                  !MW
     &,LSCAV_SO2,LSCAV_SO4AIT,LSCAV_SO4ACC,LSCAV_SO4DIS
     &,LSCAV_AGEDSOOT,Q
     &,AK,BK,DELTA_AK,DELTA_BK
     &,F_DELTA_SNOW,BLAND,CW_SEA,CW_LAND
!    &,LSC_QC,LSC_BS
     &,CLOUD_LIQ_FRAC,CLOUD_ICE_FRAC
     &,RHCRIT
     &,RHCPT, L_RHCPT
     &,VFALL,CX,CONSTP
     &)
      IMPLICIT NONE
      INTEGER
     & N        ! IN Number of points where pptn non-zero from above
!                    or where CF>CFMIN
     &,IX(N)    ! IN gather/scatter index
     &,POINTS   ! IN Number of gridpoints being processed.
      REAL
     & PSTAR(POINTS)  ! IN Surface pressure (Pa).
     &,CF(POINTS)     ! IN Cloud fraction.
     &,AK             ! IN Hybrid co-ordinate for centre of layer.
     &,BK             ! IN Hybrid co-ordinate for centre of layer.
     &,DELTA_AK       ! IN Change of hybrid co-ord across layer.
!                          (Upper minus lower).
     &,DELTA_BK       ! IN Change of hybrid co-ord across layer.
!                          (Upper minus lower).
     &,RHCRIT        ! IN Critical humidity for cloud formation.
     &,RHCPT(POINTS)     ! IN Critical humidity for cloud formation.
!    &,LSC_QC(POINTS) ! IN Large scale cloud Qc (kg/kg air).
!    &,LSC_BS(POINTS) ! IN Large scale cloud bs, moisture fluctuation.
     &,CLOUD_LIQ_FRAC(POINTS) ! IN Cloud liquid fraction.
     &,CLOUD_ICE_FRAC(POINTS) ! IN Cloud ice fraction.
     &,TIMESTEP       ! IN Timestep (sec).
     &,CW_SEA         ! IN threshold cloud liquid water content over sea
!                          for conversion to ppn (kg water per m**3).
     &,CW_LAND        ! IN threshold cloud liq. water content over land
!                          for conversion to ppn (kg water per m**3).
      LOGICAL BLAND(POINTS) ! IN Land/sea mask
     &,L_MURK         ! IN Aerosol needs scavenging.
      LOGICAL L_RHCPT ! Indicates whether RHcrit parametrization is on.
      LOGICAL L_SULPC_SO2     !IN Sulphur Cycle on, tracers scavngd if T
     &,       L_SULPC_NH3     !IN indicates if NH3 present      !MW
     &,       L_SOOT
!
      REAL
     & Q(POINTS)            ! INOUT Specific humidity (kg water/kg air).
     &,QCF(POINTS)          ! INOUT Cloud ice (kg per kg air).
     &,QCL(POINTS)          ! INOUT Cloud liquid water (kg per kg air).
     &,T(POINTS)            ! INOUT Temperature (K).
     &,AEROSOL(POINTS)      ! INOUT Aerosol (K).
     &,LSRAIN(POINTS) !INOUT Surface rainfall rate (kg per sq m per s).
     &,LSSNOW(POINTS) !INOUT Surface snowfall rate (kg per sq m per s).
     &,F_DELTA_SNOW(POINTS) ! INOUT snow fraction from ice falling as
!                                   water.
     &,VFALL(POINTS)        ! INOUT fall velocity of ice (m per s).
      REAL                    !INOUT S Cycle tracers & scavngd amounts
     &    SO2(POINTS)
     &   ,NH3(POINTS)                         !MW
     &   ,SO4_AIT(POINTS)
     &   ,SO4_ACC(POINTS)
     &   ,SO4_DIS(POINTS)
     &   ,LSCAV_SO2(POINTS)
     &   ,LSCAV_NH3(POINTS)                   !MW
     &   ,LSCAV_SO4AIT(POINTS)
     &   ,LSCAV_SO4ACC(POINTS)
     &   ,LSCAV_SO4DIS(POINTS)
     &   ,AGED_SOOT(POINTS)
     &   ,LSCAV_AGEDSOOT(POINTS)
!
C*L  Workspace usage ---------------------------------------------------
!
      REAL
     & PSTAR_C(N)        ! gathered Surface pressure (Pa).
     &,CF_C(N)           ! gathered Cloud fraction.
     &,Q_C(N)            ! gathered Specific humidity (kg water/kg air).
     &,QCF_C(N)          ! gathered Cloud ice (kg per kg air).
     &,QCL_C(N)          ! gathered Cloud liquid water (kg per kg air).
     &,T_C(N)            ! gathered Temperature (K).
     &,AERO_C(N)         ! gathered Aerosol.
     &,LSRAIN_C(N) !gathered Surface rainfall rate (kg per sq m per s).
     &,LSSNOW_C(N) !gathered Surface snowfall rate (kg per sq m per s).
     &,F_DELTA_SNOW_C(N) ! gathered fraction of snow
!    &,LSC_QC_C(N)       ! gathered Large scale cloud Qc (kg per kg air)
!    &,LSC_BS_C(N)       ! gathered Large scale cloud bs.
     &,CLF_C(N)          ! gathered Cloud liquid fraction.
     &,CIF_C(N)          ! gathered Cloud ice fraction.
     &,VFALL_C(N)        ! gathered fall velocity (m per s).
     &,RHCPT_C(N)        ! gathered Critical relative humidity
      REAL                     ! gathered S Cycle tracer arrays
     &    SO2_C(N)
     &   ,NH3_C(N)                              !MW
     &   ,SO4_AIT_C(N)
     &   ,SO4_ACC_C(N)
     &   ,SO4_DIS_C(N)
     &   ,LSCAV_SO2_C(N)
     &   ,LSCAV_NH3_C(N)                        !MW
     &   ,LSCAV_SO4AIT_C(N)
     &   ,LSCAV_SO4ACC_C(N)
     &   ,LSCAV_SO4DIS_C(N)
     &   ,AGED_SOOT_C(N)
     &   ,LSCAV_AGEDSOOT_C(N)
!
      REAL
     & RHODZ(N)       ! WORK Used for air mass p.u.a. in successive
!                            layers.
     &,P(N)           ! WORK Used for pressure at successive levels.
      LOGICAL BLAND_C(N)          ! gathered land/sea mask
!
! Call size of CX and CONSTP
! Sets up the size of arrays for CX and CONSTP
      REAL CX(16),CONSTP(16)
!
!
! Call comdecks containing large scale precipitation scavenging
! coefficients for Sulphur Cycle and soot variables
!--------------------COMDECK C_SULLSP-----------------------------
!
! LARGE SCAL PPN SCAVENGING COEFFS FOR SULPHUR CYCLE
!
      REAL                   ! parameters for S Cycle scavenging
     &     KLRAIN_SO2              ! scavenging coeff for SO2 in rain
     &    ,KLSNOW_SO2              ! scavenging coeff for SO2 in snow
     &    ,KLRAIN_NH3              ! scavenging coeff for NH3 in rain
     &    ,KLSNOW_NH3              ! scavenging coeff for NH3 in snow
     &    ,KLRAIN_SO4AIT           ! scav coeff for SO4_AITKEN in rain
     &    ,KLSNOW_SO4AIT           ! scav coeff for SO4_AITKEN in snow
     &    ,KLRAIN_SO4ACC           ! scav coeff for SO4_ACCU   in rain
     &    ,KLSNOW_SO4ACC           ! scav coeff for SO4_ACCU   in snow
     &    ,KLRAIN_SO4DIS           ! scav coeff for SO4_DISS   in rain
     &    ,KLSNOW_SO4DIS           ! scav coeff for SO4_DISS   in snow
C
      PARAMETER (                  ! Provisional values
     &     KLRAIN_SO2   = 0.5E-4   ! May be set to 0 for no scavenging
     &,KLSNOW_SO2    = 0.0
     &,KLRAIN_NH3    = 0.5E-4
     &,KLSNOW_NH3    = 0.0
     &,KLRAIN_SO4AIT = 0.0
     &,KLSNOW_SO4AIT = 0.0
     &,KLRAIN_SO4ACC = 0.0
     &,KLSNOW_SO4ACC = 0.0
     &    ,KLRAIN_SO4DIS = 0.0      ! N.B. Not same as SCONSCV value
     &    ,KLSNOW_SO4DIS = 0.0
     &    )
!
!*----------------------------------------------------------------------
!--------------------COMDECK C_ST_LSP-----------------------------
!
! LARGE SCALE PPN SCAVENGING COEFFS FOR SOOT
!
      REAL                         ! parameters for soot scavenging
     &     KLRAIN_AGEDSOOT         ! scav coeff for soot in rain
     &    ,KLSNOW_AGEDSOOT         ! scav coeff for soot in snow
C
      PARAMETER (                  ! Provisional values to turn off
!                                    washout. Motivation: at present
!                                    no LS washout of soot is used.
     &     KLRAIN_AGEDSOOT = 0.0   !
     &    ,KLSNOW_AGEDSOOT = 0.0   !
     &    )            ! Values of 0.5E-4 have been used previously.
!
!--------------------------------------------------------------------
!
!  External subroutines called -----------------------------------------
      EXTERNAL LSP_ICE,LSP_SCAV
     &        ,SLSPSCV
C*----------------------------------------------------------------------
!  Physical constants -------------------------------------------------
C*L------------------COMDECK C_G----------------------------------------
C G IS MEAN ACCEL DUE TO GRAVITY AT EARTH'S SURFACE
      REAL G

      PARAMETER(G=9.80665)
C*----------------------------------------------------------------------
      REAL P1UPONG
      PARAMETER (
     & P1UPONG=1./G        ! One upon g (sq seconds per m).
     &)
!  Define local variables ----------------------------------------------
      INTEGER I       ! Loop counters: I - horizontal field index;
!
!-----------------------------------------------------------------------
!L Internal structure.
!L 1. gather variables using index
!-----------------------------------------------------------------------
      DO I=1,N
        LSRAIN_C(I)=LSRAIN(IX(I))
        LSSNOW_C(I)=LSSNOW(IX(I))
        PSTAR_C(I) =PSTAR(IX(I))
        BLAND_C(I) =BLAND(IX(I))
        CF_C(I)=CF(IX(I))
!       LSC_QC_C(I)=LSC_QC(IX(I))
!       LSC_BS_C(I)=LSC_BS(IX(I))
        CLF_C(I)=CLOUD_LIQ_FRAC(IX(I))
        CIF_C(I)=CLOUD_ICE_FRAC(IX(I))
        QCF_C(I)=QCF(IX(I))
        QCL_C(I)=QCL(IX(I))
        Q_C(I)=Q(IX(I))
        T_C(I)=T(IX(I))
        IF (L_MURK) AERO_C(I)=AEROSOL(IX(I))
        F_DELTA_SNOW_C(I)=F_DELTA_SNOW(IX(I))
        VFALL_C(I)=VFALL(IX(I))
      END DO ! Loop over points
      IF (L_RHCPT) THEN
        DO I=1,N
          RHCPT_C(I)=RHCPT(IX(I))
        ENDDO
      ELSE
        DO I=1,N
          RHCPT_C(I)=RHCRIT
        ENDDO
      ENDIF
!
      IF (L_SULPC_SO2) THEN        ! gather S Cycle tracers
        DO I=1,N
          SO2_C(I)=SO2(IX(I))
          SO4_AIT_C(I)=SO4_AIT(IX(I))
          SO4_ACC_C(I)=SO4_ACC(IX(I))
          SO4_DIS_C(I)=SO4_DIS(IX(I))
          LSCAV_SO2_C(I)=LSCAV_SO2(IX(I))
          LSCAV_SO4AIT_C(I)=LSCAV_SO4AIT(IX(I))
          LSCAV_SO4ACC_C(I)=LSCAV_SO4ACC(IX(I))
          LSCAV_SO4DIS_C(I)=LSCAV_SO4DIS(IX(I))
        END DO
!                                                          !MW
        IF (L_SULPC_NH3) THEN                              !MW
          DO I=1,N                                         !MW
          NH3_C(I)=NH3(IX(I))                              !MW
          LSCAV_NH3_C(I)=LSCAV_NH3(IX(I))                  !MW
          END DO                                           !MW
        END IF                                             !MW
!                                                          !MW
      END IF
!
      IF (L_SOOT) THEN
        DO I=1,N
          AGED_SOOT_C(I)=AGED_SOOT(IX(I))
          LSCAV_AGEDSOOT_C(I)=LSCAV_AGEDSOOT(IX(I))
        ENDDO
      END IF
!
!-----------------------------------------------------------------------
!L 2  Calculate pressure at current level, and air mass p.u.a. of
!L    current layer.
!     (Negative in RHODZ formula takes account of sign of DELTAs.)
!-----------------------------------------------------------------------
      DO I=1,N
        P(I)=AK+PSTAR_C(I)*BK
        RHODZ(I)=-P1UPONG*(DELTA_AK+PSTAR_C(I)*DELTA_BK)
      END DO ! Loop over points
!
!-----------------------------------------------------------------------
! ICE FORMATION/EVAPORATION/MELTING
! WATER CLOUD AND RAIN FORMATION/EVAPORATION
!-----------------------------------------------------------------------
! The call to LSP_ICE replaces the calls to LSP_EVAP, LSPFRMT
! and LSP_FORM.
! CLF_C contains cloud fraction for ice
! CIF_C contains cloud fraction for water
          CALL LSP_ICE(P,RHODZ,TIMESTEP,N,
     &      RHCPT_C,
! Make available sulphur cycle parameters to influence the cloud
! microphysics as INTENT IN varables
     &      SO4_ACC_C,SO4_DIS_C,
     &      QCF_C,QCL_C,Q_C,LSRAIN_C,LSSNOW_C,VFALL_C,T_C,
     &      CLF_C,CIF_C,BLAND_C,CX,CONSTP)
!-----------------------------------------------------------------------
!L 3.4 Lose aerosol by scavenging: call LSP_SCAV
!-----------------------------------------------------------------------
!
      IF (L_MURK)  THEN
        CALL LSP_SCAV(TIMESTEP,N,LSRAIN_C,LSSNOW_C,AERO_C)
      ENDIF
!
!L  3.4.1 Scavenge Sulphur Cycle tracers: call SLSPSCV
!
       IF (L_SULPC_SO2) THEN
!
!  scavenge SO2
         IF (KLRAIN_SO2.GT.0.0 .OR. KLSNOW_SO2.GT.0.0) THEN      !MW
         CALL SLSPSCV(SO2_C,LSCAV_SO2_C,
     &                KLRAIN_SO2,KLSNOW_SO2,
     &                RHODZ,TIMESTEP,N,LSRAIN_C,LSSNOW_C)
         END IF                                                  !MW
!
!  scavenge NH3 if present                                        !MW
        IF (L_SULPC_NH3) THEN                                     !MW
!                                                                 !MW
          IF (KLRAIN_NH3.GT.0.0 .OR. KLSNOW_NH3.GT.0.0) THEN      !MW
          CALL SLSPSCV(NH3_C,LSCAV_NH3_C,                         !MW
     &                KLRAIN_NH3,KLSNOW_NH3,                      !MW
     &                RHODZ,TIMESTEP,N,LSRAIN_C,LSSNOW_C)         !MW
          END IF                                                  !MW
!                                                                 !MW
        END IF          ! end L_SULPC_NH3 condition               !MW
!                                                                 !MW
!  scavenge SO4_AIT
         IF (KLRAIN_SO4AIT.GT.0.0 .OR. KLSNOW_SO4AIT.GT.0.0) THEN  !MW
         CALL SLSPSCV(SO4_AIT_C,LSCAV_SO4AIT_C,
     &                KLRAIN_SO4AIT,KLSNOW_SO4AIT,
     &                RHODZ,TIMESTEP,N,LSRAIN_C,LSSNOW_C)
         END IF                                                    !MW
!
!  scavenge SO4_ACC
         IF (KLRAIN_SO4ACC.GT.0.0 .OR. KLSNOW_SO4ACC.GT.0.0) THEN  !MW
         CALL SLSPSCV(SO4_ACC_C,LSCAV_SO4ACC_C,
     &                KLRAIN_SO4ACC,KLSNOW_SO4ACC,
     &                RHODZ,TIMESTEP,N,LSRAIN_C,LSSNOW_C)
         END IF                                                    !MW
!
!  scavenge SO4_DIS
         IF (KLRAIN_SO4DIS.GT.0.0 .OR. KLSNOW_SO4DIS.GT.0.0) THEN
         CALL SLSPSCV(SO4_DIS_C,LSCAV_SO4DIS_C,
     &                KLRAIN_SO4DIS,KLSNOW_SO4DIS,
     &                RHODZ,TIMESTEP,N,LSRAIN_C,LSSNOW_C)
         END IF                                                    !MW
!
       END IF
!
!
      IF (L_SOOT) THEN !  Scavenge soot.
         CALL SLSPSCV(AGED_SOOT_C,LSCAV_AGEDSOOT_C,
     &                KLRAIN_AGEDSOOT,KLSNOW_AGEDSOOT,
     &                RHODZ,TIMESTEP,N,LSRAIN_C,LSSNOW_C)
      END IF
!

!-----------------------------------------------------------------------
!L 4  Scatter back arrays which will have been changed.
!L
!-----------------------------------------------------------------------
!
CDIR$ IVDEP
      DO I=1,N
        T(IX(I))=T_C(I)
        Q(IX(I))=Q_C(I)
        QCF(IX(I))=QCF_C(I)
        QCL(IX(I))=QCL_C(I)
        IF (L_MURK) AEROSOL(IX(I))=AERO_C(I)
        LSRAIN(IX(I))=LSRAIN_C(I)
        LSSNOW(IX(I))=LSSNOW_C(I)
        F_DELTA_SNOW(IX(I)) = F_DELTA_SNOW_C(I)
        VFALL(IX(I))=VFALL_C(I)
      END DO ! Loop over points
!
      IF (L_SULPC_SO2) THEN       ! scatter back S Cycle tracer arrays
        DO I=1,N
          SO2(IX(I))=SO2_C(I)
          SO4_AIT(IX(I))=SO4_AIT_C(I)
          SO4_ACC(IX(I))=SO4_ACC_C(I)
          SO4_DIS(IX(I))=SO4_DIS_C(I)
          LSCAV_SO2(IX(I))=LSCAV_SO2_C(I)
          LSCAV_SO4AIT(IX(I))=LSCAV_SO4AIT_C(I)
          LSCAV_SO4ACC(IX(I))=LSCAV_SO4ACC_C(I)
          LSCAV_SO4DIS(IX(I))=LSCAV_SO4DIS_C(I)
        END DO
!                                                           !MW
        IF (L_SULPC_NH3) THEN                               !MW
          DO I=1,N                                          !MW
          NH3(IX(I))=NH3_C(I)                               !MW
          LSCAV_NH3(IX(I))=LSCAV_NH3_C(I)                   !MW
          END DO                                            !MW
        END IF                                              !MW
!                                                           !MW
      END IF
!
      IF (L_SOOT) THEN
        DO I=1,N
          AGED_SOOT(IX(I))=AGED_SOOT_C(I)
          LSCAV_AGEDSOOT(IX(I))=LSCAV_AGEDSOOT_C(I)
        ENDDO
      END IF
!
      RETURN
      END
