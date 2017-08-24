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
      SUBROUTINE HorizontalInterpLinear
     &              (LowerBound,
     &               Len1In,  Len2In,
     &               Len1Out, Len2Out,
     &               DataExt,
     &               WtLambda,
     &               WtPhi,
     &               IOut,
     &               JOut,
     &               DataOut)
!
! Description: Performs linear interpolation of a 2-d field to a 2-d
!              set of points defined by IOut, JOut and WtLambda,
!              WtPhi.
!
! Method: This is a modified version of the routine tri_linear written
!         by Mark Mawson and described in:
!
!                The proposed semi-Lagrangian advection scheme for the
!                   semi-Implicit Unified Model integration scheme.
!                         F.R. Division working paper No 162.
!                                    Mark H. Mawson
!
!
! Owner: Stuart Bell
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   4.0   6/6/95   Equiv. to VAR code as at time of build:Stuart Bell
!
! Code Description:
!   Language:           Fortran 77 plus
!   Software Standards: "UM and Met O standards".
!
!
! Declarations:

        IMPLICIT NONE

!* Subroutine arguments
! Scalar arguments with INTENT(in):
        INTEGER     LowerBound    ! lower bounds of DataExt
        INTEGER     Len1In    ! Dimension of DataIn in i direction.
        INTEGER     Len2In    ! Dimension of DataIn in j direction.
        INTEGER     Len1Out   ! Dimension of DataOut in i direction.
        INTEGER     Len2Out   ! Dimension of DataOut in j direction.

! Array  arguments with INTENT(in):
        INTEGER     IOut (Len1Out,Len2Out)   ! Point such that
        INTEGER     JOut (Len1Out,Len2Out)   ! the desired output point
!                                            ! lies between it and it+1.
        REAL   DataExt(LowerBound:Len1In+1-LowerBound,
     &          LowerBound:Len2In+1-LowerBound)  ! Data interpolated
        REAL   WtLambda (Len1Out,Len2Out)   ! A number between 0 & 1.
        REAL   WtPhi (Len1Out,Len2Out)      ! A number between 0 & 1.

! Array  arguments with INTENT(out):
        REAL    DataOut (Len1Out,Len2Out)    ! Data interpolated to
!                                            ! desired locations.

!* End of Subroutine arguments

! Local scalars:
        INTEGER             i                 ! } Loop
        INTEGER             j                 ! } indices.
!- End of header -------------------------------------------------------

!  1.0   Perform linear interpolation in i and j  simultaneously.
! ----------------------------------------------------------------------
        DO j = 1, Len2Out
         DO i = 1, Len1Out

        DataOut (i,j)= (1.0 - WtLambda(i,j))
     &               * (1.0 - WtPhi(i,j))
     &               * DataExt (IOut(i,j), JOut(i,j))
     &              +
     &                 WtLambda(i,j)
     &               * (1.0 - WtPhi(i,j))
     &               * DataExt (IOut(i,j) + 1, JOut(i,j))
     &              +
     &                 (1.0 - WtLambda(i,j))
     &               * WtPhi(i,j)
     &               * DataExt (IOut(i,j), JOut(i,j) + 1)
     &              +
     &                 WtLambda(i,j)
     &               * WtPhi(i,j)
     &               * DataExt (IOut(i,j) + 1, JOut(i,j) + 1)

         END DO      ! Close i loop
        END DO       ! Close j loop

! End of routine.
      RETURN
      END
