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
CLL  SUBROUTINE T_INT:--------------------------------------------------
CLL
CLL  Purpose:
CLL       Carries out linear interpolation in time between two fields.
CLL       If the missing data indicator is present at one of the
CLL       times, the value at the other time is used.
CLL
CLL  Written by A. Dickinson 30/03/90
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  date
CLL
CLL  Programming standard:
CLL           Unified Model Documentation Paper No 3
CLL           Version No 1 15/1/90
CLL
CLL  System component:S190
CLL
CLL  System task: S1
CLL
CLL  Documentation:
CLL       The interpolation formulae are described in
CLL       unified model on-line documentation paper S1.
CLL
CLL  -------------------------------------------------------------------
C*L  Arguments:---------------------------------------------------------

      SUBROUTINE T_INT(DATA_T1,T1,DATA_T2,T2,DATA_T3,T3,POINTS)

      IMPLICIT NONE

      INTEGER
     * POINTS  !IN No of points to be processed

      REAL
     * DATA_T1(POINTS) !IN Data at T1
     *,DATA_T2(POINTS) !IN Data at T2
     *,DATA_T3(POINTS) !OUT Data at T3
     *,T1 !IN Time of first data field
     *,T2 !IN Time of second data field
     *,T3 !IN Time at which new field is required T1<=T3<=T2


C Local arrays:---------------------------------------------------------
C None
C ----------------------------------------------------------------------
C*L External subroutines called:----------------------------------------
C None
C*----------------------------------------------------------------------
C Local variables:------------------------------------------------------
      REAL
     * ALPHA !Fractional time

      INTEGER
     * I     !Loop index
C ----------------------------------------------------------------------
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


      ALPHA=(T3-T1)/(T2-T1)
      DO 100 I=1,POINTS
      DATA_T3(I)=DATA_T2(I)*ALPHA+DATA_T1(I)*(1-ALPHA)
      IF(DATA_T1(I).EQ.RMDI)DATA_T3(I)=DATA_T2(I)
      IF(DATA_T2(I).EQ.RMDI)DATA_T3(I)=DATA_T1(I)
100   CONTINUE

      RETURN
      END

