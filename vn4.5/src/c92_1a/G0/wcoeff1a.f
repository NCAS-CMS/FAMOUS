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
CLL  Subroutine W_COEFF------------------------------------------------
CLL
CLL  Purpose:
CLL          Calculates coefficients used to translate u and v compo-
CLL          nents of wind between equatorial (eq) latitude-longitude
CLL          grid and standard latitude-longitude grid (or visa versa).
CLL          Input latitudes and longitudes are in degrees.
CLL
CLL  Written by A. Dickinson
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  date
CLL  4.5  2/07/98 Correction for special case when no rotation occurs ie
CLL               rotated grid is identical to original. R. Rawlins
CLL
CLL  Documentation: The transformation formulae are described in
CLL                 unified model on-line documentation paper S1.
CLL
CLL  ------------------------------------------------------------------
C
C*L  Arguments:--------------------------------------------------------
      SUBROUTINE W_COEFF
     *(COEFF1,COEFF2,LAMBDA,LAMBDA_EQ,PHI_POLE,LAMBDA_POLE,POINTS)

      IMPLICIT NONE

      INTEGER
     * POINTS            !IN  Number of points to be processed

      REAL
     * COEFF1(POINTS)    !OUT Coefficient of rotation no 1
     *,COEFF2(POINTS)    !OUT Coefficient of rotation no 2
     *,LAMBDA(POINTS)    !IN  Longitude
     *,LAMBDA_EQ(POINTS) !IN  Longitude in equatorial lat-lon coords
     *,PHI_POLE          !IN  Latitude of equatorial lat-lon pole
     *,LAMBDA_POLE       !IN  Longitude of equatorial lat-lon pole
C Workspace usage:-----------------------------------------------------
C None
C----------------------------------------------------------------------
C External subroutines called:-----------------------------------------
C None
C*---------------------------------------------------------------------
C Define local varables:-----------------------------------------------
      REAL A_LAMBDA,E_LAMBDA,SIN_E_LAMBDA,SIN_PHI_POLE
     *    ,COS_PHI_POLE,C1,C2,LAMBDA_ZERO
      INTEGER I
C----------------------------------------------------------------------
C Constants from comdecks:---------------------------------------------
C*L------------------COMDECK C_PI---------------------------------------
CLL
CLL 4.0 19/09/95  New value for PI. Old value incorrect
CLL               from 12th decimal place. D. Robinson
CLL
      REAL PI,PI_OVER_180,RECIP_PI_OVER_180

      PARAMETER(
     & PI=3.14159265358979323846, ! Pi
     & PI_OVER_180 =PI/180.0,   ! Conversion factor degrees to radians
     & RECIP_PI_OVER_180 = 180.0/PI ! Conversion factor radians to
     &                              ! degrees
     & )
C*----------------------------------------------------------------------

C----------------------------------------------------------------------

CL 0.1 Check for special case of no rotation (or complete N-S rotation)
C
      
      IF(ABS(PHI_POLE).LT.90.) THEN       ! Normal rotation

CL 1. Initialise local constants
C
C Longitude of zeroth meridian
      LAMBDA_ZERO=LAMBDA_POLE+180.

C Sine and cosine of latitude of eq pole

      SIN_PHI_POLE=SIN(PI_OVER_180*PHI_POLE)
      COS_PHI_POLE=COS(PI_OVER_180*PHI_POLE)

CL 2. Evaluate translation coefficients
C
      DO 200 I=1,POINTS

C Actual longitude converted to radians

      A_LAMBDA=PI_OVER_180*(LAMBDA(I)-LAMBDA_ZERO)

C Convert eq longitude to radians and take sine

      E_LAMBDA=LAMBDA_EQ(I)*PI_OVER_180
      SIN_E_LAMBDA=SIN(E_LAMBDA)

C Formulae used are from eqs (4.19) and (4.21)

      C1=SIN(A_LAMBDA)*SIN_E_LAMBDA*SIN_PHI_POLE
     *           +COS(A_LAMBDA)*COS(E_LAMBDA)
      COEFF1(I)=C1
      C2=SQRT(1.0-C1*C1)
      COEFF2(I)=SIGN(C2,SIN_E_LAMBDA)

200   CONTINUE

      ELSE       ! Special case: no rotation (or complete N-S rotation)

       C1=SIGN(1.0,PHI_POLE)           ! =1.0 if no rotation
       DO I=1,POINTS
        COEFF1(I)=C1
        COEFF2(I)=0.
       ENDDO

      ENDIF      ! End of test for special case of no rotation

      RETURN
      END
