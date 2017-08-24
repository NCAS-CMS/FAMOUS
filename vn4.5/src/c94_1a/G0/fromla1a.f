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
CLL  Subroutine FROM_LAND_POINTS---------------------------------------
CLL
CLL  Written by A. Dickinson
CLL
CLL  Purpose:  Selects successive values from input array
CLL            DATA_LAND_POINTS and writes them to land points
CLL            in horizontal field array DATA.
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  date
CLL
CLL  Documentation: None
CLL
CLL  ------------------------------------------------------------------
C
C*L  Arguments:--------------------------------------------------------

      SUBROUTINE FROM_LAND_POINTS
     * (DATA,DATA_LAND_POINTS,LAND_SEA_MASK,POINTS,LAND_POINTS)

      IMPLICIT NONE

      INTEGER
     * POINTS !IN Total no of both land and sea points to be processed
     *,LAND_POINTS !OUT No of land points

      LOGICAL
     * LAND_SEA_MASK(POINTS)    !IN Land-sea mask

      REAL
     * DATA_LAND_POINTS(POINTS) !IN Data on land points only
     *,DATA(POINTS)             !OUT Data on land and sea points

C Workspace usage:-----------------------------------------------------
      INTEGER INDEX_LAND_POINTS(POINTS) !Scatter index for land points
C----------------------------------------------------------------------
C External subroutines called:-----------------------------------------
C*---------------------------------------------------------------------
C Local varables:------------------------------------------------------
      INTEGER I ! Integer index
C----------------------------------------------------------------------
C Constants from comdecks:---------------------------------------------
! History:
! Version  Date  Comment
!  3.4   18/5/94 Add PP missing data indicator. J F Thomson
C*L------------------COMDECK C_MDI-------------------------------------
      REAL    RMDI_PP   ! PP missing data indicator (-1.0E+30)
      REAL    RMDI_OLD  ! Old real missing data indicator (-32768.0)
      REAL    RMDI      ! New real missing data indicator (-2**30)
      INTEGER IMDI      ! Integer missing data indicator

      PARAMETER(
     & RMDI_PP  =-1.0E+30,
     & RMDI_OLD =-32768.0,
     & RMDI =-32768.0*32768.0,
     & IMDI =-32768)
C*----------------------------------------------------------------------
C----------------------------------------------------------------------

CL Initialise sea points to MDI

      DO I=1,POINTS
        DATA(I)=RMDI
      ENDDO

CL Calculate scatter index for land points

      LAND_POINTS = 0
      DO I=1,POINTS
        IF(LAND_SEA_MASK(I))THEN
          LAND_POINTS=LAND_POINTS + 1
          INDEX_LAND_POINTS(LAND_POINTS) = I
        END IF
      END DO

CL Scatter land points to array DATA
CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
      DO I=1,LAND_POINTS
        DATA(INDEX_LAND_POINTS(I))=DATA_LAND_POINTS(I)
      ENDDO

      RETURN
      END
