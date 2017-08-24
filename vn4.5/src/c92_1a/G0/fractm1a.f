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
CLL  SUBROUTINE FRAC_TIM------------------------------------------------
CLL
CLL  Purpose:  Calculates fractional time at which DATA field changes
CLL            from zero to non-zero or vice versa. The algorithm
CLL            assumes that the changes progress in the latitudindal
CLL            direction from time T1 - T2. Used for snow depth and
CLL            ice-fraction.
CLL
CLL  Written by A. Dickinson 30/03/90
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  date
CLL
CLL   3.1  15/01/93    Correct error in logic
CLL                    Author: A. Dickinson    Reviewer: C. Jones
CLL
CLL   3.2  15/03/93    Interpret RMDI a legitimate end of
CLL                    transition sequence.
CLL                    Author: A. Dickinson    Reviewer: C. Jones
CLL
CLL  Programming standard:
CLL           Unified Model Documentation Paper No 3
CLL           Version No 1 15/1/90
CLL
CLL  System component:S192
CLL
CLL  System task: S1
CLL
CLL
CLL  Documentation:
CLL            The interpolation formulae are described in
CLL            unified model on-line documentation paper S1.
CLL
CLL  -------------------------------------------------------------------
C*L  Arguments:---------------------------------------------------------

      SUBROUTINE FRAC_TIM(DATA_T1,DATA_T2,FRAC_TIME,P_ROWS,ROW_LENGTH)

      IMPLICIT NONE

      INTEGER
     * POINTS     !IN No of points to be processed
     *,ROW_LENGTH !IN Length of row
     *,P_ROWS     !IN Number of rows

      REAL
     * DATA_T1(ROW_LENGTH,P_ROWS)    !IN Data at time T1
     *,DATA_T2(ROW_LENGTH,P_ROWS)    !IN Data at time T2 where T2>T1
     *,FRAC_TIME(ROW_LENGTH,P_ROWS)  !OUT Fractional time at which DATA
     *             !changes between zero and non-zero in this time range


C Local arrays:---------------------------------------------------------
       INTEGER TYPE(ROW_LENGTH,P_ROWS) !Latitudinal transition indicator
C ----------------------------------------------------------------------
C*L External subroutines called:----------------------------------------
C None
C*----------------------------------------------------------------------
C Local variables:------------------------------------------------------
      REAL
     * ALPHA  !Fractional time
     *,A      !Denominator in fractional time calculation

      INTEGER
     * I,J,JJ,JJJ,J1,J2 !Indices
     *,ITREND           !N-S trend
     *,JHEM             !Row of equator
     *,J1J2             !Average value of J1 & J2
C ----------------------------------------------------------------------
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

C Calculate equator row
          JHEM=P_ROWS/2

CL 1. Set transition indicators  1= zero -> non zero
CL                              -1= non zero -> zero
CL                               0= no transition

      DO 100 J=1,P_ROWS
      DO 110 I=1,ROW_LENGTH

      IF(DATA_T1(I,J).EQ.0..AND.DATA_T2(I,J).GT.0.)THEN
        TYPE(I,J)=1
      ELSEIF(DATA_T2(I,J).EQ.0..AND.DATA_T1(I,J).GT.0.)THEN
        TYPE(I,J)=-1
      ELSE
        TYPE(I,J)=0
      ENDIF

C Initialise fractional time to missing data indicator or
C 0.5 if transition
      FRAC_TIME(I,J)=RMDI
      IF(TYPE(I,J).NE.0)FRAC_TIME(I,J)=.5

110   CONTINUE
100   CONTINUE


Cl 2. Search by longitude for groups where transition occurs

      DO 200 I=1,ROW_LENGTH
      DO 210 J=1,P_ROWS-1
      IF(TYPE(I,J).EQ.0.AND.TYPE(I,J+1).NE.0)THEN
        J1=J+1
        J2=J+1
        DO 220 JJ=J1+1,P_ROWS-1
        IF(TYPE(I,JJ).EQ.TYPE(I,J1))THEN
          J2=J2+1
        ELSEIF(TYPE(I,JJ).EQ.0)THEN

          A=1./FLOAT(J2-J1+2)

C Compute transition indicators
          J1J2=(J1+J2)/2
          IF((DATA_T1(I,J1-1).EQ.0.AND.DATA_T1(I,J2+1).EQ.0)
     *    .OR.(DATA_T1(I,J1-1).GT.0.AND.DATA_T1(I,J2+1).GT.0))THEN
C No transition
            ITREND=0
          ELSEIF(TYPE(I,J1).EQ.1)THEN
C Transition is zero to non-zero
            IF(J1J2.LT.JHEM)THEN
              ITREND=-1
            ELSE
              ITREND=1
            ENDIF
          ELSEIF(TYPE(I,J1).EQ.-1)THEN
C Transition is non-zero to zero
            IF(J1J2.LT.JHEM)THEN
              ITREND=1
            ELSE
              ITREND=-1
            ENDIF
          ENDIF

C ITREND indicates how FRAC_TIME varies with latitude
          IF(ITREND.EQ.1)THEN
            DO 230 JJJ=J1,J2
            FRAC_TIME(I,JJJ)=A*(J2-JJJ+1)
230         CONTINUE
          ELSE IF(ITREND.EQ.-1) THEN
            DO 240 JJJ=J1,J2
            FRAC_TIME(I,JJJ)=A*(JJJ-J1+1)
240         CONTINUE
          ENDIF

          GOTO 210
        ENDIF
220     CONTINUE
      ENDIF

210   CONTINUE
200   CONTINUE

      RETURN
      END
