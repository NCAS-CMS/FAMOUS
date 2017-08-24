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
CLL  Subroutine W_EQTOLL----------------------------------------------
CLL
CLL  Purpose:  Calculates u and v components of wind on standard
CLL            latitude-longitude grid by rotating wind
CLL            components on equatorial latitude-longitude (eq)
CLL            grid.
CLL
CLL  Written by A. Dickinson
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  date
CLL  4.1   31/05/96     The number of v points to be processed on a
CLL                     C grid differs by row_length. u,v therefore
CLL                     treated separately.
CLL                     Author I.Edmond       Reviewer D. Goddard
CLL
CLL Programming standard :
CLL
CLL Logical components covered : S133
CLL
CLL Project task :
CLL
CLL  Documentation: The transformation formulae are described in
CLL                 unified model on-line documentation paper S1.
CLL
CLLEND -----------------------------------------------------------------
C
C*L  Arguments:--------------------------------------------------------
      SUBROUTINE W_EQTOLL(COEFF1,COEFF2,U_EQ,V_EQ,U,V,POINTS,POINTS2)

      IMPLICIT NONE

      INTEGER
     * POINTS            !IN  Number of points to be processed
     *,POINTS2    ! IN  Number of v points to be processed

      REAL
     * COEFF1(POINTS)    !IN  Coefficient of rotation no 1
     *,COEFF2(POINTS)    !IN  Coefficient of rotation no 2
     *,U_EQ(POINTS)      !IN  u component of wind on equatorial grid
     *,V_EQ(POINTS)      !IN  v component of wind on equatorial grid
     *,U(POINTS)         !OUT u component of wind on lat-lon grid
     *,V(POINTS)         !OUT v component of wind on lat-lon grid
C Workspace usage:-----------------------------------------------------
C None
C----------------------------------------------------------------------
C External subroutines called:-----------------------------------------
C None
C*---------------------------------------------------------------------
C Define local varables:-----------------------------------------------
      INTEGER I
C----------------------------------------------------------------------

CL 1. Transform wind components
C
C Formulae used are from eq (4.14)

      DO 100 I = 1,POINTS
       U(I)=COEFF1(I)*U_EQ(I)+COEFF2(I)*V_EQ(I)
100   CONTINUE
      ! On a C grid number of u,v points processed differs by row_length
      DO I = 1,POINTS2
       V(I)=COEFF1(I)*V_EQ(I)-COEFF2(I)*U_EQ(I)
      END DO

      RETURN
      END
