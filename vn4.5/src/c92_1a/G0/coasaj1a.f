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
CLL  SUBROUTINE COAST_AJ-----------------------------------------------
CLL
CLL  Written by A. Dickinson
CLL
CLL  Purpose:
CLL       (i)  To produce gather indices which map each coastal point
CLL            on the target grid onto its nearest point on the source
CLL            grid. This allows correction of those surface fields
CLL            which are non-homogeneous across land/sea boundaries
CLL            after horizontal interpolation by subroutine H_INT.
CLL            The algorithm uses linear interpolation weights and
CLL            gather indices calculated by subroutine H_INT_CO.
CLL
CLL       (ii) If a land-sea mask for the target grid is not provided,
CLL            one is created. When a target land/sea mask is provided,
CLL            further index is output containing those points on the
CLL            target grid for which the 4 surrounding source points
CLL            are not of the same land/sea type as the target point.
CLL            These points will generally be new islands etc resolved b
CLL            a hires land/sea mask. They should be set to appropriate
CLL            values (eg climatology) as required.
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  Date
CLL
CLL   3.1  10/02/93    Bug removed in initialisation of variables
CLL                    Author: A. Dickinson    Reviewer: F. Rawlins
CLL   3.4  25/10/93    Output arguments INDEX_TARG_SEA_UNRES and
CLL                    INDEX_TARG_LAND_UNRES declared as arrays.
CLL                    Although previously declared as scalars, no
CLL                    detremental effects were observed from the error.
CLL                    Author: D.M.Goddard     Reviewer:
CLL
CLL  Programming standard: Unified Model Documentation Paper No 3
CLL                        Version No 1 15/1/90
CLL
CLL  Logical component number: S122
CLL
CLL  Project task: S1
CLL
CLL  Documentation: The interpolation formulae are described in
CLL                 unified model on-line documentation paper S1.
CLL
CLL  ------------------------------------------------------------------
C
C*L  Arguments:---------------------------------------------------------

      SUBROUTINE COAST_AJ
     *(INDEX_B_L,INDEX_B_R,WEIGHT_T_R,WEIGHT_B_R,WEIGHT_T_L,WEIGHT_B_L
     *,POINTS_LAMBDA_SRCE,POINTS_PHI_SRCE,POINTS,LAND_SEA_SRCE
     *,LAND_SEA_TARG,INDEX_TARG,INDEX_SRCE,COASTAL_POINTS,MASK
     *,INDEX_TARG_SEA_UNRES,SEA_POINTS_UNRES
     *,INDEX_TARG_LAND_UNRES,LAND_POINTS_UNRES)

      IMPLICIT NONE

      INTEGER
     * POINTS_LAMBDA_SRCE   !IN Number of lambda points on source grid
     *,POINTS_PHI_SRCE      !IN Number of phi points on source grid
     *,POINTS               !IN Total number of points on target grid
     *,COASTAL_POINTS       !OUT Number of coastal points on target grid
     *,LAND_POINTS_UNRES    !OUT No of unresolved land pts when MASK=T
     *,SEA_POINTS_UNRES     !OUT No of unresolved sea pts when MASK=T
     *,INDEX_B_L(POINTS)    !IN  Index of bottom lefthand corner
     *                      !    of source gridbox
     *,INDEX_B_R(POINTS)    !IN  Index of bottom righthand corner
     *                      !    of source gridbox
     *,LAND_SEA_TARG(POINTS)!INOUT Land/sea mask on target grid. If MASK
     *                      !then precalculated land/sea mask on input.
     *,LAND_SEA_SRCE(POINTS_LAMBDA_SRCE*POINTS_PHI_SRCE)
     *                      !IN Land/sea mask on source grid

      INTEGER
     * INDEX_TARG(POINTS)   !OUT Index of target coastal points
     *,INDEX_SRCE(POINTS)   !OUT Index of source points mapped onto
     *                      !target coastal points
     *,INDEX_TARG_SEA_UNRES(POINTS)
     *                      !OUT Index of sea pts on target grid which a
     *                      !unresolved when MASK=T
     *,INDEX_TARG_LAND_UNRES(POINTS)
     *                      !OUT Index of land pts on target grid which
     *                      !unresolved when MASK=T

      REAL
     * WEIGHT_T_R(POINTS) !IN  Weight applied to value at top right
     *                    !    hand corner of source gridbox
     *,WEIGHT_B_L(POINTS) !IN  Weight applied to value at bottom left
     *                    !    hand corner of source gridbox
     *,WEIGHT_B_R(POINTS) !IN  Weight applied to value at bottom right
     *                    !    hand corner of source gridbox
     *,WEIGHT_T_L(POINTS) !IN  Weight applied to value at top left
     *                    !    hand corner of source gridbox

      LOGICAL
     * MASK      !IN =F, then land/sea mask estimated
     *           !   =T, then land/sea mask input as LAND_SEA_TARG

C Local arrays:---------------------------------------------------------
      INTEGER
     * LAND_SEA_TEMP(POINTS)  ! Array used to accumulate output land/sea
     *,INDEX_TEMP(POINTS,4)   ! Index of 4 sourrounding source points
     *                        ! order by distance
     *,LAND_SEA_COAST(POINTS) ! Mask of coastal points on target grid

      REAL
     * MAX_WEIGHT(POINTS,4)   ! Linear interpolation weights ordered by

      LOGICAL
     * LOGIC_TEST(POINTS)     ! Logical string used to accumulate tests

C*L External subroutines called:----------------------------------------
C*----------------------------------------------------------------------
C Local variables:------------------------------------------------------
      REAL TEMP
      INTEGER I,J,K,START,ITEMP
C ----------------------------------------------------------------------

CL 1. Sum of land/sea mask values on source grid surrounding each point
CL   on target grid. This assumes that the mask is integer 1's and 0's.

      DO 100 I=1,POINTS

      LAND_SEA_TEMP(I)=LAND_SEA_SRCE(INDEX_B_L(I))
     *           +LAND_SEA_SRCE(INDEX_B_R(I))
     *           +LAND_SEA_SRCE(INDEX_B_L(I)-POINTS_LAMBDA_SRCE)
     *           +LAND_SEA_SRCE(INDEX_B_R(I)-POINTS_LAMBDA_SRCE)

100   CONTINUE

CL 2. Generate coastal gather indices and new land/sea mask

C 2.1 Gather index for coastal points on target grid

      DO 210 I=1,POINTS
      LOGIC_TEST(I)=LAND_SEA_TEMP(I).GT.0.AND.LAND_SEA_TEMP(I).LT.4
210   CONTINUE

      COASTAL_POINTS = 0
      DO I=1,POINTS
        IF(LOGIC_TEST(I))THEN
          COASTAL_POINTS=COASTAL_POINTS + 1
          INDEX_TARG(COASTAL_POINTS) = I
        END IF
      END DO


C Set target land-sea mask to land if all 4 surrounding
C source points are also land.

      DO 211 I=1,POINTS
      IF(LAND_SEA_TEMP(I).EQ.4)LAND_SEA_TEMP(I)=1
211   CONTINUE


C 2.2 Gather index for points which differ between input
C land/sea mask and first estimate and are not first
C estimate coastal points.

      IF(MASK)THEN
      J=COASTAL_POINTS

      DO 220 I=1,POINTS
      LAND_SEA_COAST(I)=0
220   CONTINUE
      DO 221 I=1,COASTAL_POINTS
      LAND_SEA_COAST(INDEX_TARG(I))=1
221   CONTINUE

      DO 222 I=1,POINTS
      LOGIC_TEST(I)=LAND_SEA_TEMP(I).NE.LAND_SEA_TARG(I)
     *.AND.LAND_SEA_COAST(I).EQ.0
222   CONTINUE

      START=COASTAL_POINTS
      COASTAL_POINTS = 0
      DO I=1,POINTS
        IF(LOGIC_TEST(I))THEN
          COASTAL_POINTS=COASTAL_POINTS + 1
          INDEX_TARG(START+COASTAL_POINTS) = I
        END IF
      END DO
      COASTAL_POINTS=COASTAL_POINTS+J

      ENDIF

C 2.3 Accumulate source weights and indices associated with
C     each coastal point on target grid.

      DO 230 I=1,COASTAL_POINTS

      MAX_WEIGHT(I,1)=WEIGHT_B_L(INDEX_TARG(I))
      MAX_WEIGHT(I,2)=WEIGHT_B_R(INDEX_TARG(I))
      MAX_WEIGHT(I,3)=WEIGHT_T_L(INDEX_TARG(I))
      MAX_WEIGHT(I,4)=WEIGHT_T_R(INDEX_TARG(I))
      INDEX_TEMP(I,1)=INDEX_B_L(INDEX_TARG(I))
      INDEX_TEMP(I,2)=INDEX_B_R(INDEX_TARG(I))
      INDEX_TEMP(I,3)=INDEX_B_L(INDEX_TARG(I))
     *                -POINTS_LAMBDA_SRCE
      INDEX_TEMP(I,4)=INDEX_B_R(INDEX_TARG(I))
     *                -POINTS_LAMBDA_SRCE
230   CONTINUE

C 2.4 Sort gather indices of the 4 surrounding source
C     gridpoints according to distance from target gridpoint;
C     arranged so that nearest point comes first in list (ie K=1).

      DO 240 K=1,3
      DO 240 J=K+1,4
      DO 241 I=1,COASTAL_POINTS
      IF(MAX_WEIGHT(I,K).LT.MAX_WEIGHT(I,J))THEN
      TEMP=MAX_WEIGHT(I,K)
      MAX_WEIGHT(I,K)=MAX_WEIGHT(I,J)
      MAX_WEIGHT(I,J)=TEMP
      ITEMP=INDEX_TEMP(I,K)
      INDEX_TEMP(I,K)=INDEX_TEMP(I,J)
      INDEX_TEMP(I,J)=ITEMP
      ENDIF
241   CONTINUE
240   CONTINUE


C 2.5 Initialise source gather index as nearest point on source grid
C     to target coastal point.
      DO 250  I=1,COASTAL_POINTS
      INDEX_SRCE(I)=INDEX_TEMP(I,1)
250   CONTINUE

C 2.6 Select source gather index as nearest point of same land/sea type
C     when a land/sea mask has been input on target grid
      IF(MASK) THEN

      DO 260 K=4,1,-1
      DO 261 I=1,COASTAL_POINTS
      IF(LAND_SEA_TARG(INDEX_TARG(I)).EQ.LAND_SEA_SRCE(INDEX_TEMP(I,K)))
     * INDEX_SRCE(I)=INDEX_TEMP(I,K)
261   CONTINUE
260   CONTINUE

      ENDIF

C 2.7 Update coastal values of land-sea mask

      DO 270 I=1,COASTAL_POINTS
      LAND_SEA_TEMP(INDEX_TARG(I))=LAND_SEA_SRCE(INDEX_SRCE(I))
270   CONTINUE


C 2.8 Overwrite target land/sea mask with estimate or output
C indices of target land/sea points for which none of 4 surrounding
C source points are of same type.

      IF(.NOT.MASK)THEN
        DO 280 I=1,POINTS
        LAND_SEA_TARG(I)=LAND_SEA_TEMP(I)
280     CONTINUE
      LAND_POINTS_UNRES=0
      SEA_POINTS_UNRES=0

      ELSE

      DO 281 I=1,POINTS
      LOGIC_TEST(I)=LAND_SEA_TARG(I).EQ.0.AND.
     *              (LAND_SEA_TEMP(I).NE.LAND_SEA_TARG(I))
281   CONTINUE

      SEA_POINTS_UNRES = 0
      DO I=1,POINTS
        IF(LOGIC_TEST(I))THEN
          SEA_POINTS_UNRES=SEA_POINTS_UNRES + 1
          INDEX_TARG_SEA_UNRES(SEA_POINTS_UNRES) = I
        END IF
      END DO

      DO 282 I=1,POINTS
      LOGIC_TEST(I)=LAND_SEA_TARG(I).EQ.1.AND.
     *              (LAND_SEA_TEMP(I).NE.LAND_SEA_TARG(I))
282   CONTINUE

      LAND_POINTS_UNRES = 0
      DO I=1,POINTS
        IF(LOGIC_TEST(I))THEN
          LAND_POINTS_UNRES=LAND_POINTS_UNRES + 1
          INDEX_TARG_LAND_UNRES(LAND_POINTS_UNRES) = I
        END IF
      END DO

      ENDIF

      RETURN
      END
