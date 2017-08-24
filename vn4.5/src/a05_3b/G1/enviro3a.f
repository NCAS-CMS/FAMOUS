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
CLL
CLL  MODEL            MODIFICATION HISTORY FROM MODEL VERSION 3.0:
CLL VERSION  DATE
CCL   4.0    5/05/95   : New deck added for version 3A of convection
CLL                      scheme, based on ENVIRO2B.
CLL                      Includes updating of tracers and momentum.
CLL                      Pete Inness.
CLL
CLL   4.3    03/02/97  Allow facility to switch off, under contol of a
CLL                    logical, code which cools and dries the layer
CLL                    where convection initiates to compensate for the
CLL                    initial parcel excess temperature and moisture;
CLL                    also, separately, the code which sets the 
CLL                    parcel excess in model layer 1 to the s.d. of
CLL                    the turbulent fluctuations.
CLL                                                  R.N.B.Smith
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
      SUBROUTINE ENVIRON (K,NPNTS,NP_FULL,DTHEK,DQEK,DTHEKP1,DQEKP1,
     *                    THEK,QEK,DELTAK,FLXK,THPK,QPK,
     *                    THRK,QRK,THEKP1,QEKP1,BTERM,THPKP1,
     *                    QPKP1,XPK,XPKP1,BWKP1,FLXKP1,BLOWST,
     *                    EKP14,EXK,EXKP1,DELPK,DELPKP1,AMDETK,T1_SD,
     *                    Q1_SD,L_MOM,DUEK,DVEK,DUEKP1,DVEKP1,UEK,VEK,
     *                    UPK,VPK,UEKP1,VEKP1,UPKP1,VPKP1,EFLUX_U_UD,
     *                    EFLUX_V_UD,L_SHALLOW,
     *                    L_MID,L_TRACER,NTRA,DTRAEK,DTRAEKP1,
     *                    TRAEK,TRAPK,TRAEKP1,TRAPKP1,L_XSCOMP,L_SDXS)  
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
      REAL THPIXS_DEEP,    ! INITIAL EXCESS POTENTIAL TEMPERATURE (K)
     *     QPIXS_DEEP      ! AND MIXING RATIO (KG/KG) FOR DEEP
                           ! CONVECTION
C
      REAL THPIXS_SHALLOW, ! INITIAL EXCESS POTENTIAL TEMPERATURE (K)
     *     QPIXS_SHALLOW   ! AND MIXING RATIO (KG/KG) FOR SHALLOW
                           ! CONVECTION
C
      REAL THPIXS_MID,     ! INITIAL EXCESS POTENTIAL TEMPERATURE (K)
     *     QPIXS_MID       ! AND MIXING RATIO (KG/KG) FOR MID-LEVEL
                           ! CONVECTION
C
      PARAMETER (THPIXS_DEEP= 0.2, QPIXS_DEEP =0.0)
      PARAMETER (THPIXS_SHALLOW = 0.2, QPIXS_SHALLOW =0.0)
      PARAMETER (THPIXS_MID= 0.2, QPIXS_MID =0.0)
C
C-----------------------------------------------------------------------
C VECTOR LENGTHS AND LOOP COUNTERS
C-----------------------------------------------------------------------
C
      INTEGER NPNTS          ! IN VECTOR LENGTH
C
      INTEGER NP_FULL        ! IN FULL VECTOR LENGTH
C
      INTEGER NTRA           ! IN NUMBER OF TRACERS
C
      INTEGER I,KTRA         ! LOOP COUNTERS
C
      INTEGER K              ! IN NUMBER OF MODEL LEVELS
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
      REAL UEK(NPNTS)        ! IN ENVIRONMENT U IN LAYER K (M/S)
C
      REAL UEKP1(NPNTS)      ! IN ENVIRONMENT U IN LAYER K+1 (M/S)
C
      REAL VEK(NPNTS)        ! IN ENVIRONMENT V IN LAYER K (M/S)
C
      REAL VEKP1(NPNTS)      ! IN ENVIRONMENT V IN LAYER K+1 (M/S)
C
      REAL TRAEK(NP_FULL,    ! IN TRACER OF CLOUD ENVIRONMENT
     *           NTRA)       !    IN LAYER K (KG/KG)
C
      REAL TRAEKP1(NP_FULL,  ! IN TRACER OF CLOUD ENVIRONMENT
     *             NTRA)     !    IN LAYER K+1 (KG/KG)
C
      REAL THPK(NPNTS)       ! IN PARCEL POTENTIAL TEMPERATURE IN
                             !    LAYER K (K)
C
      REAL QPK(NPNTS)        ! IN PARCEL MIXING RATIO IN LAYER K (KG/KG)
C
      REAL UPK(NPNTS)        ! IN PARCEL U IN LAYER K (M/S)
C
      REAL VPK(NPNTS)        ! IN PARCEL V IN LAYER K (M/S)
C
      REAL TRAPK(NP_FULL,    ! IN PARCEL TRACER IN LAYER K (KG/KG)
     *           NTRA)
C
      REAL THPKP1(NPNTS)     ! IN PARCEL POTENTIAL TEMPERATURE IN
                             !    LAYER K+1 (K)
C
      REAL QPKP1(NPNTS)      ! IN PARCEL MIXING RATIO IN LAYER K+1
                             !    (KG/KG)
C
      REAL UPKP1(NPNTS)      ! IN PARCEL U IN LAYER K+1 (M/S)
C
      REAL VPKP1(NPNTS)      ! IN PARCEL V IN LAYER K+1 (M/S)
C
      REAL TRAPKP1(NP_FULL,  ! IN PARCEL TRACER IN LAYER K+1
     *             NTRA)     !    (KG/KG)
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
      LOGICAL L_SHALLOW(NPNTS), !IN SWITCHES FOR TYPE OF CONVECTION
     *          L_MID(NPNTS)    !   LIKELY TO DEVELOP
C
      LOGICAL L_TRACER       ! IN SWITCH FOR INCLUSION OF TRACERS
C
      LOGICAL L_MOM          ! IN SWITCH FOR INCLUSION OF
                             !    MOMENTUM TRANSPORTS
C
      LOGICAL L_XSCOMP       ! IN Switch for allowing compensating
                             !    cooling and drying of the environment
                             !    in initiating layer
C
      LOGICAL L_SDXS         ! IN Switch for allowing parcel excess to
                             !    be set to s.d. of turbulent
                             !    fluctuations in lowest model layer
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
      REAL T1_SD(NPNTS)      ! IN Standard deviation of turbulent
C                            !    fluctuations of layer 1
C                            !    temperature (K).
      REAL Q1_SD(NPNTS)      ! IN Standard deviation of turbulent
C                            !    fluctuations of layer 1
C                            !    humidity (kg/kg).
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
      REAL DUEK(NPNTS)       ! INOUT
                             ! IN  INCREMENT TO MODEL U DUE TO
                             !     CONVECTION (M/S)
                             ! OUT UPDATED INCREMENT TO MODEL U
                             !     DUE TO CONVECTION
C
      REAL DVEK(NPNTS)       ! INOUT
                             ! IN  INCREMENT TO MODEL V DUE TO
                             !     CONVECTION (M/S)
                             ! OUT UPDATED INCREMENT TO MODEL V
                             !     DUE TO CONVECTION
C
      REAL DTRAEK(NP_FULL,   ! INOUT
     *            NTRA)      ! IN INCREMENT TO MODEL TRACER IN
                             !    LAYER K DUE TO CONVECTION
                             !    (MAY BE NON ZERO DUE TO
                             !    A PREVIOUS SPLIT FINAL DETRAINMENT
                             !    CALCULATION (KG/KG/S)
                             ! OUT UPDATED INCREMENT TO MODEL TRACER
                             !    IN LAYER K DUE TO CONVECTION
                             !    (KG/KG/S)
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
      REAL DUEKP1(NPNTS)     ! OUT INCREMENT TO MODEL U IN LAYER K+1
                             !     DUE TO CONVECTION
C
      REAL DVEKP1(NPNTS)     ! OUT INCREMENT TO MODEL V IN LAYER K+1
                             !     DUE TO CONVECTION
C
      REAL DTRAEKP1(NP_FULL, ! OUT INCREMENT TO MODEL TRACER
     *              NTRA)    !     IN LAYER K+1 DUE TO CONVECTION
                             !     (KG/KG)
C
      REAL EFLUX_U_UD(NPNTS),   ! INOUT
     *     EFLUX_V_UD(NPNTS)    ! IN  EDDY FLUX OF MOMENTUM AT BOTTOM
                                !     OF A LAYER DUE TO UD
                                ! OUT EDDY FLUX OF MOMENTUM AT TOP
                                !     OF A LAYER DUE TO UD
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
      REAL THPIXS,QPIXS      ! PARCEL EXCESS POTENTIAL TEMP(K)
                             ! AND MOISTURE(KG/KG)
C
      REAL FLX_U_KP0P5       ! FLUX OF ZONAL MOMENTUM IN CLOUD AT TOP
                             ! OF CURRENT LAYER
C
      REAL FLX_V_KP0P5       ! FLUX OF MERIDIONAL MOM. IN CLOUD AT TOP
                             ! OF CURRENT LAYER
C
C*---------------------------------------------------------------------
C
      DO I=1,NPNTS
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
       IF (BLOWST(I) .AND. L_XSCOMP) THEN
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
         IF(L_SHALLOW(I))THEN
CL
         THPIXS=THPIXS_SHALLOW
         QPIXS=QPIXS_SHALLOW
CL
         ELSEIF(L_MID(I))THEN
         THPIXS=THPIXS_MID
         QPIXS=QPIXS_MID
CL
         ELSE
CL
         THPIXS=THPIXS_DEEP
         QPIXS=QPIXS_DEEP
CL
         ENDIF
CL
         IF ( L_SDXS .AND. K .EQ. 1 ) THEN
           DTHEK(I) = DTHEK(I) - TEMPRY*MAX(THPIXS , T1_SD(I)/EXK(I))
           DQEK(I) = DQEK(I) - TEMPRY*MAX(QPIXS , Q1_SD(I))
         ELSE
           DTHEK(I) = DTHEK(I) - TEMPRY*THPIXS
           DQEK(I) = DQEK(I) - TEMPRY*QPIXS
         ENDIF
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
C
       END IF
CL
      END DO
CL
CL---------------------------------------------------------------------
CL CALCULATE EFFECT OF CONVECTION UPON MOMENTUM OF LAYER K
CL AND DO TERMINAL DETRAINMENT OF MOMENTUM
CL
CL RATE OF CHANGE OF WIND FIELD BY CONVECTION IS ESTIMATED USING A
CL DIVERGENCE OF VERTICAL EDDY MOMENTUM FLUX ACROSS THE LAYER
CL AN UPSTREAM ASSUMPTION IS USED IN THE FINITE DIFFERENCE
CL APPROXIMATIONS
CL--------------------------------------------------------------------
CL
      IF(L_MOM)THEN
C
      DO I=1,NPNTS
C----------------------------------------------------------------------
C ESTIMATE EDDY FLUX AT TOP OF CURRENT LAYER DUE TO CONVECTION
C----------------------------------------------------------------------
       FLX_U_KP0P5 = FLXK(I) * (1.0-AMDETK(I)) * (1.0-DELTAK(I)) *
     *               (1.0+EKP14(I)) * (UPK(I)-UEKP1(I))
       FLX_V_KP0P5 = FLXK(I) * (1.0-AMDETK(I)) * (1.0-DELTAK(I)) *
     *               (1.0+EKP14(I)) * (VPK(I)-VEKP1(I))
C
       IF (BLOWST(I)) THEN
C----------------------------------------------------------------------
C INITIAL CONVECTING LAYER - NO FLUX AT BASE OF LAYER
C----------------------------------------------------------------------
       DUEK(I) = DUEK(I) - FLX_U_KP0P5 / DELPK(I)
       DVEK(I) = DVEK(I) - FLX_V_KP0P5 / DELPK(I)
C----------------------------------------------------------------------
C STORE EDDY FLUX AT TOP OF LAYER READY FOR CALCULATION OF NEXT LAYER
C----------------------------------------------------------------------
       EFLUX_U_UD(I) = FLX_U_KP0P5
       EFLUX_V_UD(I) = FLX_V_KP0P5
C
       ELSE
C----------------------------------------------------------------------
C CONVECTING LAYER - TAKE EDDY FLUX DIVERGENCE ACROSS THE LAYER
C----------------------------------------------------------------------
       DUEK(I) = DUEK(I) - ( (FLX_U_KP0P5 - EFLUX_U_UD(I)) /
     *                                                DELPK(I) )
       DVEK(I) = DVEK(I) - ( (FLX_V_KP0P5 - EFLUX_V_UD(I)) /
     *                                                DELPK(I) )
C----------------------------------------------------------------------
C STORE EDDY FLUX AT TOP OF LAYER READY FOR CALCULATION OF NEXT LAYER
C----------------------------------------------------------------------
       EFLUX_U_UD(I) = FLX_U_KP0P5
       EFLUX_V_UD(I) = FLX_V_KP0P5
C
      END IF
C
      IF(BTERM(I))THEN
C----------------------------------------------------------------------
C CONVECTION TERMINATES - CALCULATE INCREMENT DUE TO CONVECTION
C IN TOP LAYER - NO FLUX OUT OF TOP OF LAYER
C----------------------------------------------------------------------
          DUEKP1(I)  = EFLUX_U_UD(I) / DELPKP1(I)
          DVEKP1(I)  = EFLUX_V_UD(I) / DELPKP1(I)
C----------------------------------------------------------------------
C ZERO EDDY FLUX OUT OF TOP OF LAYER
C----------------------------------------------------------------------
          EFLUX_U_UD(I) = 0.0
          EFLUX_V_UD(I) = 0.0
C
      END IF
C
      END DO
C
      END IF
C
CL_____________________________________________________________________
CL
CL EFFECT OF CONVECTION ON TRACER CONTENT OF LAYER K
CL (LOOPING OVER NUMBER OF TRACER VARIABLES)
CL AND DO TERMINAL DETRAINMENT OF TRACER
CL_____________________________________________________________________
CL
      IF(L_TRACER)THEN
CL
      DO KTRA = 1,NTRA
      DO I = 1,NPNTS
CL
       TEMPRY = FLXK(I)/DELPK(I)
       DTRAEK(I,KTRA) = DTRAEK(I,KTRA) + TEMPRY * (
     *
     * (1+EKP14(I)) * (1.0-DELTAK(I)) *                 ! COMPENSATING
     * (1-AMDETK(I)) * (TRAEKP1(I,KTRA)-TRAEK(I,KTRA))  ! SUBSIDENCE
     *+
     * DELTAK(I) * (1.0-AMDETK(I)) *                    ! FORCED
     * (TRAPK(I,KTRA)-TRAEK(I,KTRA))                    ! DETRAINMENT
     *+
     * AMDETK(I) * (TRAPK(I,KTRA)-TRAEK(I,KTRA))        ! MIXING
     *           )                                      ! DETRAINMENT
CL
C
        IF(BTERM(I))THEN
          TEMPRY = FLXKP1(I)/DELPKP1(I)
          DTRAEKP1(I,KTRA) = DTRAEKP1(I,KTRA) +TEMPRY*
     *                       (TRAPKP1(I,KTRA)-TRAEKP1(I,KTRA))
        END IF
C
      END DO
C
      END DO
C
      END IF
      RETURN
      END
