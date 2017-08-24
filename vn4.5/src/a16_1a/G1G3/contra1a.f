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
CLL  SUBROUTINE CONTRAIL------------------------------------------------
CLL
CLL PURPOSE: TO CALCULATE CONTRAIL UPPER AND LOWER HEIGHT
CLL
CLL  WRITTEN  BY C.M.ROBERTS
CLL  MODIFIED BY J.T.HEMING
CLL
CLL  MODEL            MODIFICATION HISTORY FROM MODEL VERSION 3.0:
CLL VERSION  DATE
CLL  4.4  30/07/97  Change limit for TERM6 to prevent calculations
CLL                 going out of T3E real number range. D. Robinson
!LL  4.5  20/04/98  Start-end args added to enable dupicate halo
!LL                 calculations to be avoided. S.D.Mullerworth
CLL
CLL  PROGRAMMING STANDARD: UNIFIED MODEL DOCUMENTATION PAPER NO. 4,
CLL  VERSION 2, DATED 18/01/90
CLL
CLL  LOGICAL COMPONENT NUMBER: D431
CLL
CLL  PROJECT TASK:
CLL
CLL  DOCUMENTATION:
CLL
CLLEND------------------------------------------------------------------
C
C*L ARGUMENTS:----------------------------------------------------------
      SUBROUTINE CONTRAIL(
C data in
     & P,T,PSTAR,
     & P_EXNER_HALF,THETA,
C data out
     & UPPER_CONTRAIL,LOWER_CONTRAIL,
C constants in
     & POINTS,P_LEVELS,
! Range of points to calculate
     & START,END)
C-----------------------------------------------------------------------
      IMPLICIT NONE
C-----------------------------------------------------------------------
      EXTERNAL ICAO_HT
C
      INTEGER
     * POINTS         ! IN  NO OF POINTS
     *,P_LEVELS       ! IN  NO OF MODEL LEVELS
     &,START,END      ! IN  Range of points to calculate
C
      REAL
     * T(POINTS,P_LEVELS)     ! IN  TEMPERATURE AT FULL LEVELS
     *,P(POINTS,P_LEVELS)     ! IN  PRESSURE    "   "     "
     *,PSTAR(POINTS)          ! IN  SURFACE PRESSURE
     *,P_EXNER_HALF(POINTS,P_LEVELS+1) ! IN  EXNER PRESSURE AT MODEL
     *                                 !     HALF LEVELS
     *,THETA(POINTS,P_LEVELS)          ! IN  POTENTIAL TEMPERATURE AT
     *                                 !     MODEL FULL LEVELS
     *,UPPER_CONTRAIL(POINTS) ! OUT  UPPER VALUE OF THE CONTRAIL
     *,LOWER_CONTRAIL(POINTS) ! OUT  LOWER VALUE OF THE CONTRAIL
C-----------------------------------------------------------------------
C Local Variables
C-----------------------------------------------------------------------
      INTEGER
     * I,J           !  LOOP COUNTERS

      REAL
     * DEL_EXNER_JM1                ! EXNER PRESSURE DIFF BET LEVELS
     *                              !                  J-3/2 AND J-1/2
     *,DEL_EXNER_J                  !   "     "  " " "J-1/2 AND J+1/2
     *,TERM1                        ! \
     *,TERM2                        !  } TEMPORARY STORAGE VARIABLES
     *,TERM3                        !  }
     *,TERM4                        !  } USED IN SEVERAL DIFFERENT
     *,TERM5                        !  }       CALCULATIONS
     *,TERM6                        ! /
     *,MINTRA_M_14(POINTS,P_LEVELS) !'MINTRA line minus 14 degrees' temp
     *,EXPONENT                     ! lapse rate between two levels *R/g
     *,P_M(POINTS)      ! Pressure of point of intersection
     *,P_MAX(POINTS)    ! Maximum pressure of intersection in a column
     *,P_MIN(POINTS)    ! Minimum pressure of intersection in a column
     *,TEST1            ! Test variable
     *,TEST2            ! Test variable
     *,P_LIMIT          ! Upper limit (lowest pressure) at which to
     *                  ! perform calculations
C-----------------------------------------------------------------------
C Logicals
C-----------------------------------------------------------------------
      LOGICAL
     * INTERSECTION(POINTS) ! TRUE if intersection occurs.
C-----------------------------------------------------------------------
C Constants
C-----------------------------------------------------------------------
C*
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

C*L------------------COMDECK C_G----------------------------------------
C G IS MEAN ACCEL DUE TO GRAVITY AT EARTH'S SURFACE
      REAL G

      PARAMETER(G=9.80665)
C*----------------------------------------------------------------------
C-----------------------------------------------------------------------
      REAL CP_OVER_TWO_G
     *,A0,A1         !   The 'MINTRA line minus 14 degrees' is
C                    !   approximated by Tm=A0*(P**A1)
      PARAMETER(CP_OVER_TWO_G=CP/(2.*G)
     &         ,A0=137.14816    ! This constant set for pressure in Pa
     &         ,A1=0.046822
     &         ,P_LIMIT=5000.0) ! No calculations above this level in Pa
C*----------------------------------------------------------------------
CL    Initialise logical array and max and min pressure arrays
C-----------------------------------------------------------------------
      DO I=START,END
        INTERSECTION(I)=.FALSE.
        P_MAX(I)=0.0
        P_MIN(I)=PSTAR(I)
      ENDDO
C-----------------------------------------------------------------------
CL  Calculate the approximation to the 'MINTRA-line minus 14 degrees'
CL   at full levels for each point
C-----------------------------------------------------------------------
      DO J=1,P_LEVELS
        DO I=START,END
          MINTRA_M_14(I,J)=A0*P(I,J)**A1
        ENDDO
      ENDDO
C-----------------------------------------------------------------------
CL Loop round all points and all levels
C-----------------------------------------------------------------------
      DO 2 J=1,P_LEVELS
        DO 1 I=START,END
          IF (J.EQ.1) THEN
C-----------------------------------------------------------------------
CL Check if level 1 temperature is less than MINTRA-14
C-----------------------------------------------------------------------
            IF (MINTRA_M_14(I,1).GT.T(I,1)) THEN
              P_M(I)=PSTAR(I)
              INTERSECTION(I)=.TRUE.
            ENDIF
          ELSEIF(P(I,J).GE.P_LIMIT)THEN
C-----------------------------------------------------------------------
CL  Continue if below pressure level P_LIMIT
CL  Calculate 'MINTRA line minus 14 degrees' - 'environment curve'
CL  at J-1 and J
C-----------------------------------------------------------------------
            TERM3=MINTRA_M_14(I,J)-T(I,J)
            TERM4=MINTRA_M_14(I,J-1)-T(I,J-1)
C-----------------------------------------------------------------------
CL  If TERM3 and TERM4 have different signs or either is equal to zero
CL  there is a point of intersection between levels j-1 and j ;
CL  otherwise there is not.
C-----------------------------------------------------------------------
            TEST1=TERM3*TERM4
            IF (TEST1.LE.0.0) THEN
              INTERSECTION(I)=.TRUE.
C-----------------------------------------------------------------------
CL Exner pressure difference across layers
C-----------------------------------------------------------------------
              DEL_EXNER_J=P_EXNER_HALF(I,J)-P_EXNER_HALF(I,J+1)
              DEL_EXNER_JM1=P_EXNER_HALF(I,J-1)-P_EXNER_HALF(I,J)
C-----------------------------------------------------------------------
CL Numerator
C-----------------------------------------------------------------------
              TERM1=T(I,J-1)-T(I,J)
C-----------------------------------------------------------------------
CL Denominator
C-----------------------------------------------------------------------
              TERM2=CP_OVER_TWO_G*(THETA(I,J-1)*DEL_EXNER_JM1
     *        +THETA(I,J)*DEL_EXNER_J)
C-----------------------------------------------------------------------
CL Lapse rate between level j-1 and j =TERM1/TERM2
CL Exponent=Lapse rate *(R/g)
C-----------------------------------------------------------------------
              EXPONENT=TERM1/TERM2*R/G
C-----------------------------------------------------------------------
CL  Find P_M, the pressure of intersection
C-----------------------------------------------------------------------
              TERM5=T(I,J-1)/A0*P(I,J-1)**(-EXPONENT)
              TERM6=A1-EXPONENT
              IF (ABS(TERM6).LT.1.0E-5) TERM6=1.0E-5 
              P_M(I)=TERM5**(1.0/TERM6)
C-----------------------------------------------------------------------
CL  Check that P_M lies between P(I,J) and P(I,J-1)
CL  If not set P_M to the upper level P(I,J)
C-----------------------------------------------------------------------
              TEST2=(P_M(I)-P(I,J))*(P_M(I)-P(I,J-1))
              IF(TEST2.GT.0.0) P_M(I)=P(I,J)
            ENDIF
          ENDIF
   1    CONTINUE
C-----------------------------------------------------------------------
CL Is P_M a maximum or minimum pressure of intersection?
C  i.e. a solution to the equation at this point.
C-----------------------------------------------------------------------
        DO I=START,END
          IF (INTERSECTION(I)) THEN
            IF (P_M(I).GT.P_MAX(I)) P_MAX(I)=P_M(I)
            IF (P_M(I).LT.P_MIN(I)) P_MIN(I)=P_M(I)
          ENDIF
        ENDDO
   2  CONTINUE
C-----------------------------------------------------------------------
CL If only one one intersection set P_MIN to zero
CL If no intersections set P_MAX and P_MIN to zero
C-----------------------------------------------------------------------
      DO I=START,END
        IF (INTERSECTION(I).AND.(P_MAX(I).EQ.P_MIN(I))) P_MIN(I)=0.0
        IF (.NOT.INTERSECTION(I)) THEN
          P_MAX(I)=0.0
          P_MIN(I)=0.0
        ENDIF
      ENDDO
C-----------------------------------------------------------------------
CL  Convert P_MAX and P_MIN into ICAO heights in thousands of feet
C-----------------------------------------------------------------------
      CALL ICAO_HT(P_MAX(START),END-START+1,LOWER_CONTRAIL(START))
      CALL ICAO_HT(P_MIN(START),END-START+1,UPPER_CONTRAIL(START))
C-----------------------------------------------------------------------
      RETURN
      END
C-----------------------------------------------------------------------
