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
CLL  SUBROUTINE LIFT_PAR-----------------------------------------------
CLL
CLL  PURPOSE : LIFTS THE PARCEL FROM LAYER K TO K+1
CLL            TAKING ENTRAINEMNT AND MOIST PROCESSES INTO ACOUNT
CLL
CLL            SUBROUTINE LATENT_H CALCULATES THE MOIST PROCESSES
CLL
CLL  SUITABLE FOR SINGLE COLUMN MODEL USE
CLL
CLL
CLL  MODEL            MODIFICATION HISTORY FROM MODEL VERSION 3.0:
CLL VERSION  DATE
CLL   4.0    5/05/95  : New deck for version 3A of convection scheme.
CLL                     Includes tracers and momentum in the convective
CLL                     parcel.
CLL                     Pete Inness.
CLL   4.5    Jul. 98  Kill the IBM specific lines (JCThil)
CLL
CLL  PROGRAMMING STANDARDS : UNIFIED MODEL DOCUMENTATION PAPER NO. 4
CLL  VERSION NO. 1
CLL
CLL  LOGICAL COMPONENTS COVERED: P27
CLL
CLL  DOCUMENTATION : UNIFIED MODEL DOCUMENTATION PAPER NO. ##
CLL
CLLEND-----------------------------------------------------------------
C
C*L  ARGUMENTS---------------------------------------------------------
C
      SUBROUTINE LIFT_PAR (NPNTS,NP_FULL,THPKP1,QPKP1,XSQKP1,BGMKP1,
     *                     BWKP1,THPK,QPK,THEKP1,QEKP1,THEK,QEK,QSEKP1,
     *                     DQSKP1,PKP1,EXKP1,EKP14,EKP34,L_MOM,UPKP1,
     *                     VPKP1,UPK,VPK,UEK,UEKP1,VEK,VEKP1,L_TRACER,
     *                     NTRA,TRAPKP1,TRAPK,TRAEKP1,TRAEK,L_SHALLOW)
C
      IMPLICIT NONE
C
C----------------------------------------------------------------------
C VECTOR LENGTHS AND LOOP COUNTERS
C----------------------------------------------------------------------
C
      INTEGER NPNTS          ! IN VECTOR LENGTH
C
      INTEGER NP_FULL        ! IN FULL VECTOR LENGTH
C
      INTEGER I,KTRA         ! LOOP COUNTERS
C
      INTEGER NTRA           ! IN NUMBER OF TRACER VARIABLES
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
      REAL UEK(NPNTS)        ! IN U OF ENVIRONMENT IN LAYER K (M/S)
C
      REAL UEKP1(NPNTS)      ! IN U OF ENVIRONMENT IN LAYER K+1 (M/S)
C
      REAL VEK(NPNTS)        ! IN V OF ENVIRONMENT IN LAYER K (M/S)
C
      REAL VEKP1(NPNTS)      ! IN V OF ENVIRONMENT IN LAYER K+1 (M/S)
C
      REAL TRAEK(NP_FULL,    ! IN TRACER CONTENT OF CLOUD
     *           NTRA)       !    ENVIRONMENT IN LAYER K (KG/KG)
C
      REAL TRAEKP1(NP_FULL,  ! IN TRACER CONTENT OF CLOUD
     *             NTRA)     !    ENVIRONMENT IN LAYER K+1 (KG/KG)
C
      REAL QSEKP1(NPNTS)     ! IN SATURATION MIXING RATIO OF CLOUD
                             !    ENVIRONMENT IN LAYER K+1 (KG/KG)
C
      REAL DQSKP1(NPNTS)     ! IN GRADIENT OF SATURATION MIXING RATIO
                             !    WITH POTENTIAL TEMPERATURE FOR THE
                             !    CLOUD ENVIRONMENT IN LAYER K+1
                             !    (KG/KG/K)
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
      REAL TRAPK(NP_FULL,    ! IN PARCEL TRACER CONTENT IN LAYER K
     *           NTRA)       !    (KG/KG)
C
      LOGICAL BWKP1(NPNTS)   ! IN MASK FOR WHETHER CONDENSATE IS
                             !    LIQUID IN LAYER K+1
C
      REAL PKP1(NPNTS)       ! IN PRESSURE AT LEVEL K+1 (PA)
C
      REAL EXKP1(NPNTS)      ! IN EXNER RATIO AT MID-POINT OF LAYER K+1
C
      REAL EKP14(NPNTS)      ! IN ENTRAINMENT COEFFICIENT AT LEVEL
                             !    K+1/4 MULTIPLIED BY APPROPRIATE
                             !    LAYER THICKNESS
C
      REAL EKP34(NPNTS)      ! IN ENTRAINMENT COEFFICIENT AT LEVEL
                             !    K+3/4 MULTIPLIED BY APPROPRIATE
                             !    LAYER THICKNESS
C
      LOGICAL L_TRACER       ! IN LOGICAL SWITCH FOR INCLUSION OF
                             !    TRACERS
C
      LOGICAL L_MOM          ! IN LOGICAL SWITCH FOR INCLUSION OF
                             !    MOMENTUM TRANSPORTS
C
      LOGICAL L_SHALLOW(NPNTS) ! IN LOGICAL INDICATOR OF WHETHER
                               !    CONVECTION AT A GRIDPOINT
                               !    TERMINATES IN THE BOUNDARY LAYER.
C
C----------------------------------------------------------------------
C VARIABLES WHICH ARE OUTPUT
C----------------------------------------------------------------------
C
      REAL THPKP1(NPNTS)     ! OUT PARCEL POTENTIAL TEMPERATURE IN
                             !     LAYER K+1 AFTER ENTRAINMENT AND
                             !     LATENT HEATING (K)
C
      REAL QPKP1(NPNTS)      ! OUT PARCEL MIXING RATIO IN LAYER K+1
                             !     AFTER ENTRAINMENT AND LATENT HEATING
                             !     (KG/KG)
C
      REAL UPKP1(NPNTS)      ! OUT PARCEL U IN LAYER K+1 AFTER
                             !     ENTRAINMENT (M/S)
C
      REAL VPKP1(NPNTS)      ! OUT PARCEL V IN LAYER K+1 AFTER
                             !     ENTRAINMENT (M/S)
C
      REAL TRAPKP1(NP_FULL,  ! OUT PARCEL TRACER CONTENT IN LAYER
     *             NTRA)     !     K+1 AFTER ENTRAINMENT. (KG/KG)
C
      REAL XSQKP1(NPNTS)     ! OUT EXCESS PARCEL WATER AFTER
                             !     LIFTING FROM LAYER K TO K+1
                             !     (KG/KG)
C
      LOGICAL BGMKP1(NPNTS)  ! OUT MASK FOR PARCELS WHICH ARE
                             !     SATURATED IN LAYER K+1 AFTER
                             !     ENTRAINMENT AND LATENT HEATING
C
C
C----------------------------------------------------------------------
C VARIABLES WHICH ARE DEFINED LOCALLY
C
C ON THE IBM ARRAYS ARE ALLOCATED USING A PARAMETER STATEMENT
C
C ON THE CRAY ARRAYS ARE DYNAMICALLY ALLOCATED
C----------------------------------------------------------------------
C
      REAL THPKP1T(NPNTS)    ! INITIAL ESTIMATE OF PARCEL TEMPERATURE
                             ! IN LAYER K+1 AFTER ENTRAINMENT (K)
C
      REAL TT(NPNTS)         ! TEMPORARY TEMPERATURE USED IN CALCULATION
                             ! OF SATURATION MIXING RATIO (K)
C
      REAL QSPKP1(NPNTS)     ! SATURATION MIXING RATIO OF PARCEL
                             ! AFTER DRY ASCENT (KG/KG)
C
C
C----------------------------------------------------------------------
C EXTERNAL ROUTINES CALLED
C----------------------------------------------------------------------
C
      EXTERNAL QSAT,LATENT_H
C
C*---------------------------------------------------------------------
C
      DO I=1,NPNTS
CL
CL----------------------------------------------------------------------
CL  LIFT PARCEL MIXING RATIO, POTENTIAL TEMPERATURE, U, V AND TRACER
CL  TO THE NEXT LEVEL
CL----------------------------------------------------------------------
CL
CL----------------------------------------------------------------------
CL  INITIAL 'DRY' ASCENT
CL
CL  UM DOCUMENTATION PAPER P27
CL  SECTION (3), EQUATIONS (11B), (12B)
CL----------------------------------------------------------------------
CL
       THPKP1(I) = (  THPK(I)
     *             + EKP14(I)*THEK(I) + EKP34(I)*(1.+EKP14(I))*THEKP1(I)
     *             ) / ((1.+EKP14(I))*(1.+EKP34(I)))
C
       QPKP1(I) = (  QPK(I)
     *             + EKP14(I)*QEK(I) + EKP34(I)*(1.+EKP14(I))*QEKP1(I)
     *             ) / ((1.+EKP14(I))*(1.+EKP34(I)))
C
       END DO
C
      IF(L_MOM)THEN
C
      DO I=1,NPNTS
       UPKP1(I) = (  UPK(I)
     *             + EKP14(I)*UEK(I) + EKP34(I)*(1.+EKP14(I))*UEKP1(I)
     *             ) / ((1.+EKP14(I))*(1.+EKP34(I)))
C
       VPKP1(I) = (  VPK(I)
     *             + EKP14(I)*VEK(I) + EKP34(I)*(1.+EKP14(I))*VEKP1(I)
     *             ) / ((1.+EKP14(I))*(1.+EKP34(I)))
C----------------------------------------------------------------------
C IF CONVECTION IS DEEP OR MID-LEVEL, ADD AN IN-CLOUD PRESSURE
C GRADIENT TERM TO THE MOMENTUM INCREMENTS
C----------------------------------------------------------------------
C
       IF(.NOT.L_SHALLOW(I))THEN
        UPKP1(I) = UPKP1(I) - (0.7*(UEK(I)-UEKP1(I))/(1.0+EKP34(I)))
        VPKP1(I) = VPKP1(I) - (0.7*(VEK(I)-VEKP1(I))/(1.0+EKP34(I)))
       END IF
      END DO
C
      END IF
C
       IF(L_TRACER)THEN
C
       DO KTRA = 1,NTRA
       DO I = 1,NPNTS
C
       TRAPKP1(I,KTRA) = ( TRAPK(I,KTRA)
     * + EKP14(I)*TRAEK(I,KTRA) + EKP34(I)*(1.+EKP14(I))*TRAEKP1(I,KTRA)
     *             ) / ((1.+EKP14(I))*(1.+EKP34(I)))
C
       END DO
       END DO
C
       END IF
C
C-----------------------------------------------------------------------
C   CALCULATE WHERE THE PARCEL IS SUPERSATURATED (IE WHERE GAMMA(K+1)=1
C   SEE DCTN 29 PAGE 123)
C-----------------------------------------------------------------------
C
C
C-----------------------------------------------------------------------
C CONVERT POTENTIAL TEMPERATURE TO TEMPERATURE AND CALCULATE
C PRESSURE OF LAYER K FOR CALCULATION OF SATURATED
C MIXING RATIO
C-----------------------------------------------------------------------
C
      DO I=1,NPNTS
       TT(I) = THPKP1(I)*EXKP1(I)
      END DO
      CALL QSAT (QSPKP1,TT,PKP1,NPNTS)
C
      DO 20 I=1,NPNTS
       BGMKP1(I) = QPKP1(I) .GT. QSPKP1(I)
CL
CL----------------------------------------------------------------------
CL  CONDENSATION CALCULATION
CL
CL  SUBROUTINE LATENT_H
CL
CL  UM DOCUMENTATION PAPER P27
CL  SECTION (4)
CL----------------------------------------------------------------------
CL
       THPKP1T(I) = THPKP1(I)
   20 CONTINUE
C
      CALL LATENT_H (NPNTS,THPKP1T,QPKP1,THEKP1,QSEKP1,DQSKP1,
     *               BGMKP1,BWKP1,EXKP1)
C
C-----------------------------------------------------------------------
C   CALCULATE A MORE ACCURATE PARCEL SATURATED MIXING RATIO AND CONDENSE
C   OUT ANY EXCESS WATER VAPOUR. STORE THE EXCESS AMOUNTS IN 'XSQKP1'
C   FOR LATER.  SET PARCEL POTENTIAL TEMPERATURES TO THE PROVISIONAL
C   VALUES EXCEPT WHERE THE PARCEL IS NOT SUPERSATURATED WRT THE NEW
C   SATURATED MIXING RATIO. RECALCULATE BIT VECTOR 'BGMKP1'.
C-----------------------------------------------------------------------
C
C
C-----------------------------------------------------------------------
C CONVERT POTENTIAL TEMPERATURE TO TEMPERATURE AND CALCULATE
C PRESSURE OF LAYER K FOR CALCULATION OF SATURATED
C MIXING RATIO
C-----------------------------------------------------------------------
C
      DO 35 I = 1,NPNTS
       TT(I) = THPKP1T(I)*EXKP1(I)
   35 CONTINUE
      CALL QSAT (QSPKP1,TT,PKP1,NPNTS)
C
      DO 40 I=1,NPNTS
       XSQKP1(I) = QPKP1(I) - QSPKP1(I)
C
       IF(XSQKP1(I) .LE. 0.0) THEN
         BGMKP1(I) = .FALSE.
         XSQKP1(I) = 0.0
       ELSE
         BGMKP1(I) = .TRUE.
         THPKP1(I) = THPKP1T(I)
       END IF
C
       QPKP1(I) = QPKP1(I) - XSQKP1(I)
   40 CONTINUE
C
      RETURN
      END
