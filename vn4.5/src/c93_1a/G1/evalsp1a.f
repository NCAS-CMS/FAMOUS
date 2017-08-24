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
CLL  SUBROUTINE EVAL_SP-------------------------------------------------
CLL
CLL  PURPOSE:   To evaluate a cubic spline at a given point.
CLL             Outputs spline value,1st and 2nd derivatives.
CLL             Can be called after subroutine SPLINE which outputs the
CLL             second derivative at all data points.
CLL  Tested under compiler CFT77
CLL  Tested under OS version 5.1
CLL
CLL  Author J.T.Heming
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  Date
CLL vn4.2  07/11/96:Allow this routine to be used in C93_2A (Farnon)
CLL
CLL  Logical components covered :D413
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
      SUBROUTINE EVAL_SP(
C data in
     & X,Y,AMU,H_POINTS,V_POINTS,X0,
C data out
     & Y0,DERIV1,DERIV2)
C*L
C-----------------------------------------------------------------------
      IMPLICIT NONE
C*----------------------------------------------------------------------
      INTEGER
     * H_POINTS     ! The number of splines to be evaluated
     *,V_POINTS     ! The number of points in each spline profile
C-----------------------------------------------------------------------
      REAL
     * X(V_POINTS)            ! IN Independent variable array input data
     *,Y(H_POINTS,V_POINTS)   ! IN Dependent variable array input data
     *,AMU(H_POINTS,V_POINTS) ! IN Array of 2nd derivatives
     *,X0                     ! IN X-value at which spline is evaluated
     *,Y0(H_POINTS)           ! OUT Array of spline values
     *,DERIV1(H_POINTS)       ! OUT Array of 1st derivatives
     *,DERIV2(H_POINTS)       ! OUT Array of 2nd derivatives
C*
C*L
C-----------------------------------------------------------------------
C Local Variables
C-----------------------------------------------------------------------
      INTEGER
     * I,K,J    ! Loop counters
C-----------------------------------------------------------------------
      REAL
     * XD1,XD2,XD3,RXD3,XD12,RXD3O6 ! Work variables
     *,WORK1(H_POINTS)              ! Work array
C-----------------------------------------------------------------------
      LOGICAL
     * FOUND    ! Set true once spline evaluated
C-----------------------------------------------------------------------
CL    Check that X0 lies within range of X-values
CL    N.B. This routine assumes X decreases with height
C-----------------------------------------------------------------------
      IF(X0.LE.X(1).AND.X0.GE.X(V_POINTS))THEN
C-----------------------------------------------------------------------
        FOUND=.FALSE.
C-----------------------------------------------------------------------
CL    Find adjacent X values between which X0 lies
C-----------------------------------------------------------------------
        DO J=1,V_POINTS-1
          IF(X0.GT.X(J+1).AND.(.NOT.FOUND))THEN
            FOUND=.TRUE.
            K=J
          ELSEIF(X0.EQ.X(V_POINTS).AND.(.NOT.FOUND))THEN
            FOUND=.TRUE.
            K=V_POINTS-1
          ENDIF
        ENDDO
C-----------------------------------------------------------------------
CL    Set local variables
C-----------------------------------------------------------------------
        XD1=X(K+1)-X0
        XD2=X0-X(K)
        XD3=X(K+1)-X(K)
        RXD3=1.0/XD3
        RXD3O6=RXD3/6.0
        XD12=XD1*XD2/6.0
C-----------------------------------------------------------------------
CL    Loop through horizontal field
C-----------------------------------------------------------------------
        DO I=1,H_POINTS
C-----------------------------------------------------------------------
CL    Calculate factor used in calculation of Y0 and DERIV1
C-----------------------------------------------------------------------
          WORK1(I)=(AMU(I,K+1)*(XD2+XD3)+AMU(I,K)*(XD1+XD3))
     *      *RXD3O6
C-----------------------------------------------------------------------
CL    Evaluate spline value Y0 at point X0
C-----------------------------------------------------------------------
          Y0(I)=((Y(I,K+1)*XD2)+(Y(I,K)*XD1))*RXD3
          Y0(I)=Y0(I)-WORK1(I)*XD1*XD2
C-----------------------------------------------------------------------
CL    Evaluate first derivative at point X0
C-----------------------------------------------------------------------
          DERIV1(I)=(Y(I,K+1)-Y(I,K)+(AMU(I,K+1)-AMU(I,K))*XD12)
     *      *RXD3
          DERIV1(I)=DERIV1(I)+WORK1(I)*(XD2-XD1)
C-----------------------------------------------------------------------
CL    Evaluate second deriavtive at point X0
C-----------------------------------------------------------------------
          DERIV2(I)=((AMU(I,K+1)*XD1)+(AMU(I,K)*XD2))*RXD3
        ENDDO
C-----------------------------------------------------------------------
C     Error message
C-----------------------------------------------------------------------
      ELSE
        WRITE(6,999)X0
 999    FORMAT(' Spline not evaluated - X0 out of range. X0=',F10.4)
      ENDIF
C=======================================================================
C     END OF SUBROUTINE EVAL_SP
C=======================================================================
      RETURN
      END
C=======================================================================
