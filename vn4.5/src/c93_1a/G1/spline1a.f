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
CLL  SUBROUTINE SPLINE--------------------------------------------------
CLL
CLL  PURPOSE:   To fit a cubic spline to vertical data profiles.
CLL             Outputs the second derivative at the input data points.
CLL             To evaluate the spline and its first two derivatives at
CLL             a particular point call subroutine EVAL_SP.
CLL  Tested under compiler CFT77
CLL  Tested under OS version 5.1
CLL
CLL  Author J.T.Heming
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  date
CLL vn4.2  07/11/96:Allow this routine to be used in C93_2A (Farnon)
CLL
CLL  Logical components covered: D413
CLL
CLL  Project TASK:
CLL
CLL  Programming standard: U M DOC  Paper NO. 4,
CLL
CLL  External documentation
CLL
CLLEND------------------------------------------------------------------
C
C*L  ARGUMENTS:---------------------------------------------------------
      SUBROUTINE SPLINE(
C data in
     & X,Y,H_POINTS,V_POINTS,
C data out
     & AMU)
C*L
C-----------------------------------------------------------------------
      IMPLICIT NONE
C*----------------------------------------------------------------------
      INTEGER
     * H_POINTS ! The number of sets of data for which splines to be
     *          !   calculated
     *,V_POINTS ! The number of points at which each spline is to be
     *          !   evaluated
C-----------------------------------------------------------------------
      REAL
     * X(V_POINTS)  ! Indendent variable array - n.b. assumes that x-
     *              ! profile is the same for all H_POINTS on a V_POINT
     *,Y(H_POINTS,V_POINTS)      ! Dependent variable array
     *,AMU(H_POINTS,V_POINTS)    ! Output array of second derivative
C*
C*L
C-----------------------------------------------------------------------
C Local Variables
C-----------------------------------------------------------------------
      INTEGER
     * I,K      ! Loop counters
C-----------------------------------------------------------------------
      REAL
     * H(V_POINTS)
     *,RH(V_POINTS)
     *,A(V_POINTS)
     *,RB
     *,WORK1(H_POINTS)
     *,WORK2(H_POINTS)
C-----------------------------------------------------------------------
CL    Initialise top and bottom values of AMU
C-----------------------------------------------------------------------
      DO I=1,H_POINTS
        AMU(I,1)=0.0
        AMU(I,V_POINTS)=0.0
      ENDDO
C-----------------------------------------------------------------------
CL  Calculate the difference between adjacent linear first derivatives
C-----------------------------------------------------------------------
      DO K=1,V_POINTS-1
        H(K)=X(K+1)-X(K)
        RH(K)=1.0/H(K)
        DO I=1,H_POINTS
          WORK1(I)=(Y(I,K+1)-Y(I,K))*6.0*RH(K)
        ENDDO
        IF(K.GT.1)THEN
          DO I=1,H_POINTS
            AMU(I,K)=WORK1(I)-WORK2(I)
          ENDDO
        ENDIF
        DO I=1,H_POINTS
          WORK2(I)=WORK1(I)
        ENDDO
      ENDDO
C-----------------------------------------------------------------------
      H(V_POINTS)=0.0
      RH(V_POINTS)=0.0
      A(2)=0.5
      DO K=2,V_POINTS
        IF(K.GT.2)THEN
          A(K)=H(K-1)*RB*RH(K-2)
        ENDIF
        RB=1.0/(2.0*(H(K-1)+H(K))*RH(K-1)-A(K))
        DO I=1,H_POINTS
          AMU(I,K)=RH(K)*AMU(I,K)
        ENDDO
        DO I=1,H_POINTS
          AMU(I,K)=RB*(AMU(I,K)-AMU(I,K-1))
        ENDDO
      ENDDO
C-----------------------------------------------------------------------
CL    Calculate second derivatives at all but top level
C-----------------------------------------------------------------------
CL *** Following loop labelled to workaround fmp mistranslation
C
      DO 100 K=1,V_POINTS-1
        DO I=1,H_POINTS
          AMU(I,V_POINTS-K)=AMU(I,V_POINTS-K)-(A(V_POINTS-K+1)*
     *    AMU(I,V_POINTS-K+1))
        ENDDO
 100  CONTINUE
C=======================================================================
C     END OF SUBROUTINE SPLINE
C=======================================================================
      RETURN
      END
C=======================================================================
