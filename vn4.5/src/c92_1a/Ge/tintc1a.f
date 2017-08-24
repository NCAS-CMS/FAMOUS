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
CLL  SUBROUTINE T_INT_C-------------------------------------------------
CLL
CLL  Purpose:
CLL    Carries out linear interpolation in time between two fields at
CLL    times T1 and T2. If the missing data indicator is present at one
CLL    of the times, the value at the other time is used. The interpolat
CLL    is controlled by a field ZI. A prescribed value is inserted where
CLL    If ZI changes between 0 and non-zero in the period T1 - T2, then
CLL    the field is linearly interpolated between its value at the time
CLL    ZI is non-zero and the prescibed value at the time when ZI become
CLL    The fractional time at which ZI changes between 0 and non-zero in
CLL    period T1 - T2 must be provided as input.
CLL
CLL  Written by A. Dickinson 30/03/90
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  date
CLL   3.2  17/03/93 : Correct for rounding problem, ie case where alpha
CLL                   should be exactly equal to frac_time.
CLL                   Author: R.A Stratton
CLL
CLL  Programming standard: Unified Model Documentation Paper No 3
CLL                        Version No 1 15/1/90
CLL
CLL  System component:S191
CLL
CLL  System task: S1
CLL
CLL  Documentation:  The interpolation formulae are described in
CLL            unified model on-line documentation paper S1.
CLL
CLL  -------------------------------------------------------------------
C*L  Arguments:---------------------------------------------------------

      SUBROUTINE T_INT_C(DATA_T1,T1,DATA_T2,T2,DATA_T3,T3,POINTS
     *,FRAC_TIME,ZI_T1,PRES_VALUE)

      IMPLICIT NONE

      INTEGER
     * POINTS             !IN No of points to be processed

      REAL
     * DATA_T1(POINTS)    !IN Data at T1
     *,DATA_T2(POINTS)    !IN Data at T2
     *,DATA_T3(POINTS)    !OUT D_ta at T3
     *,ZI_T1(POINTS)      !IN Value of controlling fieled at T1
     *,PRES_VALUE(POINTS) !IN Prescribed value of Data when ZI=0
     *,FRAC_TIME(POINTS)  !IN Fractional time at which ZI changes betwee
     *                    !zero and non-zero in this time range
     *,T1 !IN Time of first data field
     *,T2 !In Time of second data field
     *,T3 !IN Time at which new field is required T1<=T3<=T2


C Local arrays:---------------------------------------------------------
       REAL INT_TIME(POINTS)
C ----------------------------------------------------------------------
C*L External subroutines called:----------------------------------------
      EXTERNAL T_INT
C*----------------------------------------------------------------------
C Local variables:------------------------------------------------------
      REAL
     * ALPHA !Fractional time
     *,EPSILON     ! rounding error
     *,ALPHA_PLUS  ! add rounding error to alpha
     *,ALPHA_MINUS ! alpha minus rounding error

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

C set rounding error
      EPSILON=1.0E-6

      CALL T_INT(DATA_T1,T1,DATA_T2,T2,DATA_T3,T3,POINTS)

      ALPHA=(T3-T1)/(T2-T1)
      ALPHA_PLUS=ALPHA+EPSILON
      ALPHA_MINUS=ALPHA-EPSILON

      DO 100 I=1,POINTS

      IF(FRAC_TIME(I).NE.RMDI)THEN
        IF(ZI_T1(I).EQ.0.0)THEN
           IF(ALPHA_MINUS.LT.FRAC_TIME(I))THEN
             DATA_T3(I)=PRES_VALUE(I)
           ELSE
             INT_TIME(I)=(ALPHA-FRAC_TIME(I))/(1.-FRAC_TIME(I))
             DATA_T3(I)=PRES_VALUE(I)*(1.-INT_TIME(I))
     *                 +DATA_T2(I)*INT_TIME(I)
           ENDIF
        ELSE
           IF(ALPHA_PLUS.GT.FRAC_TIME(I))THEN
             DATA_T3(I)=PRES_VALUE(I)
           ELSE
             INT_TIME(I)=(FRAC_TIME(I)-ALPHA)/(FRAC_TIME(I))
             DATA_T3(I)=PRES_VALUE(I)*(1.-INT_TIME(I))
     *                 +DATA_T1(I)*INT_TIME(I)
           ENDIF
        ENDIF
      ENDIF

100   CONTINUE




      RETURN
      END

