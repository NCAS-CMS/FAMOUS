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
!+ Performs Area weighted horizitontal interpolation
!
! Subroutine Interface:
      SUBROUTINE H_INT_AW(ROWS_IN,ROWS_OUT
     &,                   ROW_LENGTH_IN,ROW_LENGTH_OUT,GLOBAL
     &,                   AW_INDEX_TARG_LHS,AW_INDEX_TARG_TOP
     &,                   AW_COLAT_T,AW_LONG_L,DATA_IN,DATA_OUT)

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
! 4.0     16/03/95   Original code. D.M. Goddard
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
      INTEGER      ROWS_IN          !No of rows on source grid
      INTEGER      ROWS_OUT         !No of rows on target grid
      INTEGER      ROW_LENGTH_IN    !No of pts per row on source grid
      INTEGER      ROW_LENGTH_OUT   !No of pts per row on target grid
      LOGICAL      GLOBAL           !True if global area required

!   Array  arguments with intent(in):
      INTEGER      AW_INDEX_TARG_LHS(ROW_LENGTH_OUT+1)
                                    !Index of source box overlapping
                                    !lhs of target grid-box
      INTEGER      AW_INDEX_TARG_TOP(ROWS_OUT+1)
                                    !Index of source box overlapping
                                    !top of target grid-box
      REAL         AW_COLAT_T(ROWS_OUT+1)
                                    !Colatitude of top of target grd-box
                                    ! (in units of DELTA_LAT_SRCE)
      REAL         AW_LONG_L(ROW_LENGTH_OUT+1)
                                    !Left longitude of target grid-box
                                    ! (in units of DELTA_LONG_SRCE)
      REAL         DATA_IN(ROW_LENGTH_IN*ROWS_IN)
                                    !Data before interpolation

!   Array  arguments with intent(out):
      REAL         DATA_OUT(ROW_LENGTH_OUT*ROWS_OUT)
                                    !Data after interpolation

! Local scalars:
      INTEGER      I

! Local arrays:
      REAL         BOXSUM(ROW_LENGTH_OUT,ROWS_OUT)
                                    !Sum of data on target grid

! Function & Subroutine calls:
      External BOX_SUM

!- End of header

!     1. Calculate sum of contribution from gridboxes

      CALL BOX_SUM(ROW_LENGTH_IN,ROWS_IN,ROW_LENGTH_OUT,ROWS_OUT
     &,            AW_LONG_L,AW_COLAT_T
     &,            AW_INDEX_TARG_LHS,AW_INDEX_TARG_TOP
     &,            GLOBAL,DATA_OUT,DATA_IN)


      RETURN
      END

