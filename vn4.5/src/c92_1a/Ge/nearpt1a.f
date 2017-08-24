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
CLL  SUBROUTINE NEAR_PT-------------------------------------------------
CLL
CLL  Purpose:  To produce gather indices which map each point on the
CLL            target grid onto the nearest point on the source grid.
CLL            This allows interpolation by choosing the value of the
CLL            nearest neighbour. The code uses coefficients and gather
CLL            indices calculated by subroutine H_INT_CO.
CLL
CLL  Written by A. Dickinson
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  date
CLL
CLL  Programming standard:
CLL           Unified Model Documentation Paper No 3
CLL           Version No 1 15/1/90
CLL
CLL  System component:
CLL
CLL  System task: S123
CLL
CLL  Documentation: The interpolation formulae are described in
CLL                 unified model on-line documentation paper S1.
CLL
CLL  -------------------------------------------------------------------
C*L  Arguments:---------------------------------------------------------

      SUBROUTINE NEAR_PT
     *(INDEX_B_L,INDEX_B_R,WEIGHT_T_R,WEIGHT_B_R,WEIGHT_T_L,WEIGHT_B_L
     *,POINTS,POINTS_LAMBDA_SRCE,INDEX_NEAREST)

      IMPLICIT NONE

      INTEGER
     * POINTS_LAMBDA_SRCE   !IN Number of lambda points on source grid
     *,POINTS               !IN Total number of points on target grid
     *,INDEX_B_L(POINTS)    !IN  Index of bottom lefthand corner
     *                      !    of source gridbox
     *,INDEX_B_R(POINTS)    !IN  Index of bottom righthand corner
     *                      !    of source gridbox
     *,INDEX_NEAREST(POINTS)!OUT Index of nearest source point to
     *                      ! each target point

      REAL
     * WEIGHT_T_R(POINTS) !IN  Weight applied to value at top right
     *                    !    hand corner of source gridbox
     *,WEIGHT_B_L(POINTS) !IN  Weight applied to value at bottom left
     *                    !    hand corner of source gridbox
     *,WEIGHT_B_R(POINTS) !IN  Weight applied to value at bottom right
     *                    !    hand corner of source gridbox
     *,WEIGHT_T_L(POINTS) !IN  Weight applied to value at top left
     *                    !    hand corner of source gridbox

C Local arrays:---------------------------------------------------------
      INTEGER
     * INDEX_TEMP(POINTS,4)   ! Index of 4 sourrounding source points
     *                        ! ordered by distance

      REAL
     * MAX_WEIGHT(POINTS,4)   ! Linear interpolation weights ordered by
                              ! distance

C*L External subroutines called:----------------------------------------
C None
C*----------------------------------------------------------------------
C Local variables:------------------------------------------------------
      REAL TEMP
      INTEGER I,J,K,ITEMP
C ----------------------------------------------------------------------

C 1.  Accumulate source weights and indices associated with
C     each coastal point on target grid.

      DO 100 I=1,POINTS

      MAX_WEIGHT(I,1)=WEIGHT_B_L(I)
      MAX_WEIGHT(I,2)=WEIGHT_B_R(I)
      MAX_WEIGHT(I,3)=WEIGHT_T_L(I)
      MAX_WEIGHT(I,4)=WEIGHT_T_R(I)
      INDEX_TEMP(I,1)=INDEX_B_L(I)
      INDEX_TEMP(I,2)=INDEX_B_R(I)
      INDEX_TEMP(I,3)=INDEX_B_L(I)
     *                -POINTS_LAMBDA_SRCE
      INDEX_TEMP(I,4)=INDEX_B_R(I)
     *                -POINTS_LAMBDA_SRCE
100   CONTINUE

C 2.  Sort gather indices of the 4 surrounding source
C     gridpoints according to distance from target gridpoint;
C     arranged so that nearest point comes first in list (ie K=1).

      DO 200 K=1,3
      DO 200 J=K+1,4
      DO 210 I=1,POINTS
      IF(MAX_WEIGHT(I,K).LT.MAX_WEIGHT(I,J))THEN
      TEMP=MAX_WEIGHT(I,K)
      MAX_WEIGHT(I,K)=MAX_WEIGHT(I,J)
      MAX_WEIGHT(I,J)=TEMP
      ITEMP=INDEX_TEMP(I,K)
      INDEX_TEMP(I,K)=INDEX_TEMP(I,J)
      INDEX_TEMP(I,J)=ITEMP
      ENDIF
210   CONTINUE
200   CONTINUE

C 3. Assign index of nearest source point to output array

      DO 300 I=1,POINTS
      INDEX_NEAREST(I)=INDEX_TEMP(I,1)
300   CONTINUE

      RETURN
      END
