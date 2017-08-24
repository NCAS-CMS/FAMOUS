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
CLL  SUBROUTINE DD_INIT------------------------------------------------
CLL
CLL  PURPOSE : ROUTINE TO INITIALISE THE DOWNDRAUGHT
CLL
CLL  SUITABLE FOR SINGLE COLUMN MODEL USE
CLL
CLL
CLL  MODEL            MODIFICATION HISTORY FROM MODEL VERSION 3.0:
CLL VERSION  DATE
CLL   4.0   5/5/95    New deck added for version 3A of convection
CLL                   scheme. Includes tracers and momentum in the
CLL                   convective parcel. removes model level
CLL                   dependence in initiation of downdraught.
CLL                   Pete Inness.
!     4.2  10/01/97   Split up IF statement to prevent use of
!                     uninitialised values in BUOY array. D. Robinson.
!     4.3  19/03/97   Split up another IF statement to prevent use of
!                     uninitialised values in BUOY array. D. Robinson.
CLL   4.5    Jul. 98  Kill the IBM specific lines (JCThil)
CLL
CLL  PROGRAMMING STANDARDS : UNIFIED MODEL DOCUMENTATION PAPER NO. 3
CLL  VERSION NO. 4  DATED 5/2/92
CLL
CLL  LOGICAL COMPONENTS COVERED:
CLL
CLL  SYSTEM TASK : P27
CLL
CLL  DOCUMENTATION : UNIFIED MODEL DOCUMENTATION PAPER P27
CLL
CLLEND-----------------------------------------------------------------
C
C*L  ARGUMENTS---------------------------------------------------------
C
      SUBROUTINE DD_INIT(NPNTS,NP_FULL,TH_UD_K,Q_UD_K,THE_K,QE_K,PK,
     &                   EXK,THDD_K,QDD_K,DELTD,DELQD,BDD_START,K,BDDI,
     &                   BDD_ON,L_MOM,U_UD_K,V_UD_K,UE_K,VE_K,UDD_K,
     &                   VDD_K,DELUD,DELVD,L_TRACER,NTRA,TRA_UD_K,
     &                   TRAE_K,TRADD_K,DELTRAD)
C
      IMPLICIT NONE
C
C-----------------------------------------------------------------------
C MODEL CONSTANTS
C-----------------------------------------------------------------------
C*L------------------COMDECK C_EPSLON-----------------------------------
C EPSILON IS RATIO OF MOLECULAR WEIGHTS OF WATER AND DRY AIR
      REAL EPSILON,C_VIRTUAL

      PARAMETER(EPSILON=0.62198,
     &          C_VIRTUAL=1./EPSILON-1.)
C*----------------------------------------------------------------------

C-----------------------------------------------------------------------
C VECTOR LENGTHS AND LOOP COUNTERS
C-----------------------------------------------------------------------
C
C
      INTEGER I,KTRA            ! LOOP COUNTERS
C
      INTEGER NPNTS             ! VECTOR LENGTH
C
      INTEGER NP_FULL           ! FULL VECTOR LENGTH
C
      INTEGER NTRA              ! NUMBER OF TRACER VARIABLES
C
      INTEGER K                 ! IN PRESENT MODEL LAYER
C
C-----------------------------------------------------------------------
C VARIABLES WHICH ARE INPUT
C-----------------------------------------------------------------------
C
      REAL THE_K(NPNTS)         ! IN POTENTIAL TEMPERATURE OF
                                !    ENVIRONMENT IN LAYER K (K)
C
      REAL TH_UD_K(NPNTS)       ! IN PARCEL POTENTIAL TEMPERATURE OF
                                !    UPDRAUGHT, LAYER K (K)
C
      REAL QE_K(NPNTS)          ! IN MIXING RATIO OF ENVIRONMENT IN
                                !    LAYER K (KG/KG)
C
      REAL Q_UD_K(NPNTS)        ! IN PARCEL MIXING RATIO OF UPDRAUGHT,
                                !    LAYER K (KG/KG)
C
      REAL UE_K(NPNTS)          ! IN U IN ENVIRONMENT IN LAYER K (M/S)
C
      REAL U_UD_K(NPNTS)        ! IN PARCEL U OF UPDRAUGHT IN LAYER K
                                !    (M/S)
C
      REAL VE_K(NPNTS)          ! IN V IN ENVIRONMENT IN LAYER K (M/S)
C
      REAL V_UD_K(NPNTS)        ! IN PARCEL V OF UPDRAUGHT IN LAYER K
                                !    (M/S)
C
      REAL TRAE_K(NP_FULL,NTRA) ! IN TRACER CONTENT OF ENVIRONMENT
                                !    IN LAYER K (KG/KG)
C
      REAL TRA_UD_K(NP_FULL,    ! IN PARCEL TRACER CONTENT OF
     *              NTRA)       !    UPDRAUGHT IN LAYER K (KG/KG)
C
      REAL EXK(NPNTS)           ! IN EXNER RATIO OF LAYER K
C
      REAL PK(NPNTS)            ! IN PRESSURE OF LAYER K (PA)
C
      LOGICAL BDDI(NPNTS)       ! IN MASK FOR THOSE POINTS WHERE
                                !    DOWNDRAUGHT MAY INITIATE
C
      LOGICAL BDD_ON(NPNTS)     ! IN MASK FOR THOSE POINTS WHERE
                                !    DOWNDRAUGHT IS ON
C
      LOGICAL L_TRACER          ! IN SWITCH FOR INCLUSION OF TRACERS
C
      LOGICAL L_MOM             ! IN SWITCH FOR INCLUSION OF
                                !    MOMENTUM TRANSPORTS
C
C-----------------------------------------------------------------------
C VARIABLES WHICH ARE INPUT AND OUTPUT
C-----------------------------------------------------------------------
C
      LOGICAL BDD_START(NPNTS)  ! INOUT
                                ! IN  MASK FOR THOSE POINT WHERE
                                !     DOWNDRAUGHT MAY START
                                ! OUT MASK FOR THOSE POINTS WHERE
                                !
C
C-----------------------------------------------------------------------
C VARIABLES WHICH ARE OUTPUT
C-----------------------------------------------------------------------
C
      REAL THDD_K(NPNTS)        ! OUT DOWNDRAUGHT POTENTIAL TEMPERATURE
                                !     OF LAYER K
C
      REAL QDD_K(NPNTS)         ! OUT DOWNDRAUGHT MIXING RATIO OF
                                !     LAYER K
C
      REAL UDD_K(NPNTS)         ! OUT DOWNDRAUGHT U IN LAYER K (M/S)
C
      REAL VDD_K(NPNTS)         ! OUT DOWNDRAUGHT V IN LAYER K (M/S)
C
      REAL TRADD_K(NP_FULL,     ! OUT DOWNDRAUGHT TRACER CONTENT OF
     *             NTRA)        !     LAYER K
C
      REAL DELTD(NPNTS)         ! OUT COOLING NECESSARY TO ACHIEVE
                                !     SATURATION
C
      REAL DELQD(NPNTS)         ! OUT MOISTENING NECESSARY TO ACHIEVE
                                !     SATURATION
C
      REAL DELUD(NPNTS)         ! OUT CHANGE TO ENVIRONMENT U DUE TO
                                !     DOWNDRAUGHT FORMATION (M/S)
C
      REAL DELVD(NPNTS)         ! OUT CHANGE TO ENVIRONMENT V DUE TO
                                !     DOWNDRAUGHT FORMATION (M/S)
C
      REAL DELTRAD(NP_FULL,NTRA)! OUT DEPLETION OF ENVIRONMENT TRACER
                                !     DUE TO FORMATION OF DOWNDRAUGHT
C
C-----------------------------------------------------------------------
C VARIABLES WHICH ARE DEFINED LOCALLY
C-----------------------------------------------------------------------
C
C
      REAL TH_MEAN(NPNTS)       ! MEAN POTENTIAL TEMPERATURE USED IN
                                ! CALCULATION OF SATURATED DOWNDRAUGHT
                                ! POTENTIAL TEMPERATURE IN LAYER K
C
      REAL Q_MEAN(NPNTS)        ! MEAN MIXING RATIO USED IN CALCULATION
                                ! OF SATURATED DOWNDRAUGHT
                                ! MIXING RATIO FOR LAYER K
C
      REAL U_MEAN(NPNTS)        ! MEAN U USED IN CALCULATION OF DELUD
                                ! FOR LAYER K
C
      REAL V_MEAN(NPNTS)        ! MEAN V USED IN CALCULATION OF DELVD
                                ! FOR LAYER K
C
      REAL TRA_MEAN(NPNTS,NTRA) ! MEAN TRACER USED AS INITIAL TRACER
                                ! CONTENT OF DOWNDRAUGHT IN LAYER K
                                ! (KG/KG)
C
      REAL T_MEAN(NPNTS)        ! MEAN TEMPERATURE USED IN CALCULATION
                                ! OF SATURATED DOWNDRAUGHT POTENTIAL
                                ! TEMPERATURE OF LAYER K (K)
C
      REAL THDDS(NPNTS)         ! SATURATED DOWNDRAUGHT POTENTIAL
                                ! TEMPERATURE IN LAYER K (K)
C
      REAL QDDS(NPNTS)          ! SATURATED DOWNDRAUGHT MIXING RATIO
                                ! IN LAYER K (KG/KG)
C
      REAL BUOY(NPNTS)          ! BUOYANCY OF PARCEL IN LAYER K
C
C
      REAL THDD_V               ! VIRTUAL POTENTIAL TEMPERATURE OF
                                ! PARCEL IN LAYER K
C
      REAL THE_V                ! VIRTUAL POTENTIAL TEMPERATURE OF
                                ! ENVIRONMENT IN LAYER K
C
C-----------------------------------------------------------------------
C EXTERNAL ROUTINES CALLED
C-----------------------------------------------------------------------
C
      EXTERNAL SATCAL
C
C-----------------------------------------------------------------------
C CALCULATE MEAN TEMPERATURE, MIXING RATIO, U, V AND TRACER
C-----------------------------------------------------------------------
C
      DO I=1,NPNTS
       TH_MEAN(I) = (THE_K(I)+TH_UD_K(I))*0.5
       Q_MEAN(I) = (QE_K(I)+Q_UD_K(I))*0.5
       T_MEAN(I) = TH_MEAN(I)*EXK(I)
      END DO
C
      IF(L_MOM)THEN
       DO I=1,NPNTS
        U_MEAN(I) = (UE_K(I)+U_UD_K(I))*0.5
        V_MEAN(I) = (VE_K(I)+V_UD_K(I))*0.5
       END DO
      END IF
      IF(L_TRACER)THEN
C
      DO KTRA=1,NTRA
        DO I=1,NPNTS
          TRA_MEAN(I,KTRA) = (TRAE_K(I,KTRA)+TRA_UD_K(I,KTRA))*0.5
        END DO
      END DO
C
      END IF
C
C
C-----------------------------------------------------------------------
C CALCULATE SATURATED DOWNDRAUGHT POTENTIAL TEMPERATURE FOR LAYER K
C-----------------------------------------------------------------------
C
      CALL SATCAL(NPNTS,T_MEAN,TH_MEAN,PK,QDDS,THDDS,K,EXK,Q_MEAN,
     *            THE_K)
C
C-----------------------------------------------------------------------
C IS SATURATED PARCEL NEGATIVELY BUOYANT COMPARED TO ENVIRONMENT
C-----------------------------------------------------------------------
C
      DO I=1,NPNTS
       IF (.NOT. BDD_ON(I) .AND. BDDI(I) ) THEN
          THDD_V = THDDS(I)*(1.0+C_VIRTUAL*QDDS(I))
          THE_V = THE_K(I)*(1.0+C_VIRTUAL*QE_K(I))
          BUOY(I) = THDD_V - THE_V
C
          IF (BUOY(I) .LT. 0.5 ) THEN
C
C-----------------------------------------------------------------------
C INITIATE DOWNDRAUGHT
C-----------------------------------------------------------------------
C
             THDD_K(I) = THDDS(I)
             QDD_K(I) = QDDS(I)
             BDD_START(I) = .TRUE.
C
C-----------------------------------------------------------------------
C CALCULATE COOLING AND MOISTENING TO ACHIEVE SATURATION
C-----------------------------------------------------------------------
C
             DELTD(I) = THDDS(I)-THE_K(I)
             DELQD(I) = QDDS(I)-QE_K(I)
          END IF
       END IF
      END DO
C
      IF(L_MOM)THEN
        DO I=1,NPNTS
          IF(.NOT.BDD_ON(I).AND.BDDI(I))THEN
            IF(BUOY(I).LT.0.5)THEN
             UDD_K(I) = U_MEAN(I)
             VDD_K(I) = V_MEAN(I)
             DELUD(I) = UDD_K(I)-UE_K(I)
             DELVD(I) = VDD_K(I)-VE_K(I)
            END IF
          END IF
        END DO
      END IF
C
C
      IF(L_TRACER)THEN
C
        DO KTRA=1,NTRA
          DO I=1,NPNTS
            IF(.NOT.BDD_ON(I).AND.BDDI(I).AND.K.GE.4)THEN
              IF(BUOY(I).LT.0.5)THEN
              TRADD_K(I,KTRA) = TRA_MEAN(I,KTRA)
              DELTRAD(I,KTRA) = TRADD_K(I,KTRA)-TRAE_K(I,KTRA)
              END IF
            END IF
          END DO
        END DO
C
      END IF
      RETURN
      END
C
