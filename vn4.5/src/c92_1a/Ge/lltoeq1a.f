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
CLL  Subroutine LLTOEQ-------------------------------------------------
CLL
CLL  Purpose:  Calculates latitude and longitude on equatorial
CLL            latitude-longitude (eq) grid used in regional
CLL            models from input arrays of latitude and
CLL            longitude on standard grid. Both input and output
CLL            latitudes and longitudes are in degrees.
CLL
CLL  Written by A. Dickinson
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  date
CLL    4.5   18/08/98  Added DEF,FLDOP   (A Van der Wal)     
CLL
CLL Programming standard :
CLL
CLL Logical components covered : S132
CLL
CLL Project task :
CLL
CLL  Documentation: The transformation formulae are described in
CLL                 unified model on-line documentation paper S1.
CLLEND -----------------------------------------------------------------
C
C*L  Arguments:--------------------------------------------------------
      SUBROUTINE LLTOEQ
     *(PHI,LAMBDA,PHI_EQ,LAMBDA_EQ,PHI_POLE,LAMBDA_POLE,POINTS)

      IMPLICIT NONE

      INTEGER
     * POINTS            !IN  Number of points to be processed

      REAL
     * PHI(POINTS)       !IN  Latitude
     *,LAMBDA(POINTS)    !IN  Longitude
     *,LAMBDA_EQ(POINTS) !OUT Longitude in equatorial lat-lon coords
     *,PHI_EQ(POINTS)    !OUT Latitude in equatorial lat-lon coords
     *,PHI_POLE          !IN  Latitude of equatorial lat-lon pole
     *,LAMBDA_POLE       !IN  Longitude of equatorial lat-lon pole
C Workspace usage:-----------------------------------------------------
C None
C ---------------------------------------------------------------------
C External subroutines called:-----------------------------------------
C None
C*---------------------------------------------------------------------
C Define local varables:-----------------------------------------------
      REAL A_LAMBDA,A_PHI,E_LAMBDA,ARG,E_PHI,SIN_PHI_POLE,COS_PHI_POLE
      REAL TERM1,TERM2,SMALL,LAMBDA_ZERO
      INTEGER I
      PARAMETER(SMALL=1.0E-6)
C ---------------------------------------------------------------------
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

C ---------------------------------------------------------------------

CL 1. Initialise local constants
C
C Latitude of zeroth meridian
      LAMBDA_ZERO=LAMBDA_POLE+180.
C Sine and cosine of latitude of eq pole
      SIN_PHI_POLE=SIN(PI_OVER_180*PHI_POLE)
      COS_PHI_POLE=COS(PI_OVER_180*PHI_POLE)

CL 2. Transform from standard to equatorial latitude-longitude

      DO 200 I= 1,POINTS

C Scale longitude to range -180 to +180 degs

      A_LAMBDA=LAMBDA(I)-LAMBDA_ZERO
      IF(A_LAMBDA.GT. 180.0)A_LAMBDA=A_LAMBDA-360.
      IF(A_LAMBDA.LE.-180.0)A_LAMBDA=A_LAMBDA+360.

C Convert latitude & longitude to radians

      A_LAMBDA=PI_OVER_180*A_LAMBDA
      A_PHI=PI_OVER_180*PHI(I)

C Compute eq latitude using equation (4.4)

      ARG=-COS_PHI_POLE*COS(A_LAMBDA)*COS(A_PHI)
     *                   +SIN(A_PHI)*SIN_PHI_POLE
      ARG=MIN(ARG, 1.0)
      ARG=MAX(ARG,-1.0)
      E_PHI=ASIN(ARG)
      PHI_EQ(I)=RECIP_PI_OVER_180*E_PHI

C Compute eq longitude using equation (4.6)

      TERM1 =(COS(A_PHI)*COS(A_LAMBDA)*SIN_PHI_POLE
     *       +SIN(A_PHI)*COS_PHI_POLE)
      TERM2=COS(E_PHI)
      IF(TERM2.LT.SMALL) THEN
        E_LAMBDA=0.0
      ELSE
        ARG=TERM1/TERM2
        ARG=MIN(ARG, 1.0)
        ARG=MAX(ARG,-1.0)
        E_LAMBDA=RECIP_PI_OVER_180*ACOS(ARG)
        E_LAMBDA=SIGN(E_LAMBDA,A_LAMBDA)
      ENDIF

C Scale longitude to range 0 to 360 degs

      IF(E_LAMBDA.GE.360.0) E_LAMBDA=E_LAMBDA-360.0
      IF(E_LAMBDA.LT.0.0) E_LAMBDA=E_LAMBDA+360.0
      LAMBDA_EQ(I)=E_LAMBDA

200   CONTINUE

      RETURN
      END
