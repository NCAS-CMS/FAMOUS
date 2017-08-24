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
CLL  Subroutine TO_LAND_POINTS-----------------------------------------
CLL
CLL  Purpose:  Selects land points values from input horizontal field DA
CLL            and writes then as contiguous values to array
CLL            DATA_LAND_POINTS.
CLL
CLL  Written by A. Dickinson 11/9/90
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  date
CLL
CLL Programming standard :
CLL
CLL Logical components covered : S171
CLL
CLL Project task :
CLL
CLL  Documentation:  None
CLL
CLL  ------------------------------------------------------------------
C
C*L  Arguments:--------------------------------------------------------

      SUBROUTINE TO_LAND_POINTS
     * (DATA,DATA_LAND_POINTS,LAND_SEA_MASK,POINTS,LAND_POINTS)

      IMPLICIT NONE

      INTEGER
     * POINTS !IN Total no of both land and sea points to be processed
     *,LAND_POINTS !OUT No of land points

      LOGICAL
     * LAND_SEA_MASK(POINTS)    !IN Land-sea mask

      REAL
     * DATA_LAND_POINTS(POINTS) !OUT Data on land points only
     *,DATA(POINTS)             !IN Data on land and sea points

C Workspace usage:-----------------------------------------------------
      INTEGER INDEX_LAND_POINTS(POINTS) !Gather index for land points
C----------------------------------------------------------------------
C External subroutines called:-----------------------------------------
C*---------------------------------------------------------------------
C Local varables:------------------------------------------------------
      INTEGER I ! Integer index
C----------------------------------------------------------------------

CL Compute gather index for land points

      LAND_POINTS = 0
      DO I=1,POINTS
        IF(LAND_SEA_MASK(I))THEN
          LAND_POINTS=LAND_POINTS + 1
          INDEX_LAND_POINTS(LAND_POINTS) = I
        END IF
      END DO

CL Gather land points from input array DATA

      DO I=1,LAND_POINTS
        DATA_LAND_POINTS(I)=DATA(INDEX_LAND_POINTS(I))
      ENDDO

      RETURN
      END
