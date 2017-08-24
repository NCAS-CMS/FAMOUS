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
CLL  Subroutine EQTOLL-------------------------------------------------
CLL
CLL  Purpose:  Calculates latitude and longitude on standard grid
CLL            from input arrays of latitude and longitude on
CLL            equatorial latitude-longitude (eq) grid used
CLL            in regional models. Both input and output latitudes
CLL            and longitudes are in degrees.
CLL
CLL  Written by A. Dickinson
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  Date
CLL
CLL  Documentation: The transformation formulae are described in
CLL                 unified model on-line documentation paper S1.
CLL
CLL Logical components covered : S131
CLL
CLL Project task :
CLL
CLL External documentation:
CLL
CLLEND-----------------------------------------------------------------
C
C*L  Arguments:--------------------------------------------------------
      SUBROUTINE EQTOLL
     *(PHI_EQ,LAMBDA_EQ,PHI,LAMBDA,PHI_POLE,LAMBDA_POLE,POINTS)

      IMPLICIT NONE

      INTEGER
     * POINTS            !IN  Number of points to be processed

      REAL
     * PHI(POINTS)       !OUT Latitude
     *,LAMBDA(POINTS)    !OUT Longitude (0 =< LON < 360)
     *,LAMBDA_EQ(POINTS) !IN  Longitude in equatorial lat-lon coords
     *,PHI_EQ(POINTS)    !IN  Latitude in equatorial lat-lon coords
     *,PHI_POLE          !IN  Latitude of equatorial lat-lon pole
     *,LAMBDA_POLE       !IN  Longitude of equatorial lat-lon pole

C Workspace usage:-----------------------------------------------------
C None
C----------------------------------------------------------------------
C External subroutines called:-----------------------------------------
C None
C*---------------------------------------------------------------------
C Local varables:------------------------------------------------------
      REAL E_LAMBDA,E_PHI,A_LAMBDA,ARG,A_PHI,SIN_PHI_POLE,COS_PHI_POLE
      REAL TERM1,TERM2,SMALL,LAMBDA_ZERO
      INTEGER I
      PARAMETER(SMALL=1.0E-6)
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

CL 1. Initialise local constants
C
C Latitude of zeroth meridian
      LAMBDA_ZERO=LAMBDA_POLE+180.
C Sine and cosine of latitude of eq pole
      SIN_PHI_POLE=SIN(PI_OVER_180*PHI_POLE)
      COS_PHI_POLE=COS(PI_OVER_180*PHI_POLE)

CL 2. Transform from equatorial to standard latitude-longitude

      DO 200 I= 1,POINTS

C Scale eq longitude to range -180 to +180 degs

      E_LAMBDA=LAMBDA_EQ(I)
      IF(E_LAMBDA.GT. 180.0) E_LAMBDA=E_LAMBDA-360.0
      IF(E_LAMBDA.LT.-180.0) E_LAMBDA=E_LAMBDA+360.0

C Convert eq latitude & longitude to radians

      E_LAMBDA=PI_OVER_180*E_LAMBDA
      E_PHI=PI_OVER_180*PHI_EQ(I)

C Compute latitude using equation (4.7)

      ARG=COS_PHI_POLE*COS(E_LAMBDA)*COS(E_PHI)
     *                   +SIN(E_PHI)*SIN_PHI_POLE
      ARG=MIN(ARG, 1.0)
      ARG=MAX(ARG,-1.0)
      A_PHI=ASIN(ARG)
      PHI(I)=RECIP_PI_OVER_180*A_PHI

C Compute longitude using equation (4.8)

      TERM1 =(COS(E_PHI)*COS(E_LAMBDA)*SIN_PHI_POLE
     *       -SIN(E_PHI)*COS_PHI_POLE)
      TERM2=COS(A_PHI)
      IF(TERM2.LT.SMALL) THEN
        A_LAMBDA=0.0
      ELSE
        ARG=TERM1/TERM2
        ARG=MIN(ARG, 1.0)
        ARG=MAX(ARG,-1.0)
        A_LAMBDA=RECIP_PI_OVER_180*ACOS(ARG)
        A_LAMBDA=SIGN(A_LAMBDA,E_LAMBDA)
        A_LAMBDA=A_LAMBDA+LAMBDA_ZERO
      END IF

C Scale longitude to range 0 to 360 degs

      IF(A_LAMBDA.GE.360.0) A_LAMBDA=A_LAMBDA-360.0
      IF(A_LAMBDA.LT.0.0) A_LAMBDA=A_LAMBDA+360.0
      LAMBDA(I)=A_LAMBDA

200   CONTINUE

      RETURN
      END
