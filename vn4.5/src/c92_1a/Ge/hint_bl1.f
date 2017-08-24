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
!+ Performs Bi-linear horizitontal interpolation
!
! Subroutine Interface:
      SUBROUTINE H_INT_BL(ROWS_IN,ROW_LENGTH_IN,LEN_FIELD_OUT
     &,                   INDEX_B_L,INDEX_B_R,DATA_IN
     &,                   WEIGHT_B_L,WEIGHT_B_R,WEIGHT_T_L,WEIGHT_T_R
     &,                   DATA_OUT)

CLL  System component: S121
CLL
CLL  System task: S1
CLL
CLL  Purpose:
CLL
CLL  Documentation:
CLL            The interpolation formulae are described in
CLL            unified model on-line documentation paper S1.
CLL
      IMPLICIT NONE
!
! Description:
!   Carries out bi-linear horizontal interpolation using coefficients
!   and gather indices calculated in subroutine H_INT_CO
!
! Method:
!   See UMDP S1 for full desciption
!
! Current Code Owner: D.M. Goddard
!
! History:
! Version   Date     Comment
! -------   ----     -------
! 3.0      ??????   Original code. A.Dickinson
! 4.0      160395   Renamed and brought up to new standard D.M. Goddard
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!
! System component covered: S121
! System Task:              S1
!
! Declarations:
!   These are of the form:-
!     INTEGER      ExampleVariable      !Description of variable
!
! Subroutine arguments
!   Scalar arguments with intent(in):
      INTEGER  ROWS_IN              !No of P rows on source grid
      INTEGER  ROW_LENGTH_IN        !No of pts per row on source grid
      INTEGER  LEN_FIELD_OUT        !No of points on target grid

!   Array  arguments with intent(in):
      INTEGER  INDEX_B_L(LEN_FIELD_OUT)
                                     !Index of bottom lefthand corner
                                     !  of source gridbox
      INTEGER  INDEX_B_R(LEN_FIELD_OUT)
                                     !Index of bottom righthand corner
                                     !  of source gridbox
      REAL     DATA_IN(ROWS_IN*ROW_LENGTH_IN)
                                      !Data before interpolation
      REAL     WEIGHT_B_L(LEN_FIELD_OUT)
                                     !Weight applied to value at bottom
                                     !lefthand corner of source gridbox
      REAL     WEIGHT_B_R(LEN_FIELD_OUT)
                                     !Weight applied to value at bottom
                                     !righthand corner of source gridbox
      REAL     WEIGHT_T_L(LEN_FIELD_OUT)
                                     !Weight applied to value at top
                                     !lefthand corner of source gridbox
      REAL     WEIGHT_T_R(LEN_FIELD_OUT)
                                     !Weight applied to value at top
                                     !righthand corner of source gridbox

!   Array  arguments with intent(out):
      REAL     DATA_OUT(LEN_FIELD_OUT) !Data after interpolation

! Local scalars:
      INTEGER      I

! Function & Subroutine calls:
!     External None

!- End of header

!     1. Carry out horizontal interpolation using equation (2.1)

      DO I=1,LEN_FIELD_OUT

      DATA_OUT(I)=WEIGHT_B_L(I)*DATA_IN(INDEX_B_L(I))
     &           +WEIGHT_B_R(I)*DATA_IN(INDEX_B_R(I))
     &           +WEIGHT_T_L(I)*DATA_IN(INDEX_B_L(I)-ROW_LENGTH_IN)
     &           +WEIGHT_T_R(I)*DATA_IN(INDEX_B_R(I)-ROW_LENGTH_IN)

      END DO

      RETURN
      END

