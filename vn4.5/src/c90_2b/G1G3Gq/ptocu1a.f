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
CLL  SUBROUTINE P_TO_CU---------------------------------------------
CLL
CLL  Purpose:  Interpolates a horizontal field from pressure to wind
CLL            points on an Arakawa C grid. This routine carries out
CLL            E-W interpolation to u-point. Under UPDATE
CLL            identifier GLOBAL the data is assumed periodic along
CLL            rows. Otherwise, the first and last value on each row are
CLL            calculated using one-sided differencing.
CLL
CLL  Not suitable for single column use.
CLL
CLL  Written 12/9/91 by A. Dickinson
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  date
CLL
CLL  Programming standard: Unified Model Documentation Paper No 3
CLL                        Version No 3 18/12/90
CLL
CLL  System component: S101
CLL
CLL  System task: S1
CLL
CLL  Documentation:  The equation used is (2.1)
CLL                  in unified model documentation paper No. S1
CLL
CLLEND-------------------------------------------------------------

C
C*L  Arguments:---------------------------------------------------
      SUBROUTINE P_TO_CU
     1  (P_DATA,U_DATA,P_FIELD,U_FIELD,ROW_LENGTH,ROWS)

      IMPLICIT NONE

      INTEGER
     *  ROWS               !IN    Number of rows to be updated.
     *, ROW_LENGTH         !IN    Number of points per row
     *, P_FIELD            !IN    Number of points in input field
     *, U_FIELD            !IN    Number of points in output field

      REAL
     * P_DATA(P_FIELD)     !INOUT Data on p points
     *,U_DATA(U_FIELD)     !  OUT Data on uv points
C*---------------------------------------------------------------------

C*L  Local arrays:-----------------------------------------------------
C    None
C*---------------------------------------------------------------------

C*L  External subroutine calls:---------------------------------------
C    None
C*---------------------------------------------------------------------

C----------------------------------------------------------------------
C    Define local variables
C----------------------------------------------------------------------
      INTEGER
     *  U_POINTS      !     Number of values at u points
     *,I,M            !     Horizontal loop indices

C---------------------------------------------------------------------
CL    1.     Initialise local constants
C---------------------------------------------------------------------

      U_POINTS      =  ROW_LENGTH * ROWS

C---------------------------------------------------------------------
CL    2.     Calculate horizontal average at u points
C---------------------------------------------------------------------

      DO 200 I=1,U_POINTS-1
       U_DATA(I)=0.5*(P_DATA(I)+P_DATA(I+1))
200   CONTINUE

C  End points

C Cyclic wrap around
      DO 201 I=ROW_LENGTH,U_POINTS,ROW_LENGTH
       U_DATA(I)=0.5*(P_DATA(I)+P_DATA(I+1-ROW_LENGTH))
201   CONTINUE

      RETURN
      END
