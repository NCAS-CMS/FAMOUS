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
      SUBROUTINE HorizontalInterpCubic
     &              (LowerBound,
     &               Len1In,  Len2In,
     &               Len1Out, Len2Out,
     &               DataExt,
     &               WtLambda,
     &               WtPhi,
     &               IOut,
     &               JOut,
     &               DataOut)

! Description: Performs cubic Lagrange interpolation of a 2-d field to a
!              2-d set of points defined by IOut, JOut, and WtLambda,
!              WtPhi.
!
! Method: This is a modified version of the routine cubic_lagrange
!         written by Mark Mawson and described in:
!
!           The proposed semi-Lagrangian advection scheme for the
!               semi-Implicit Unified Model integration scheme.
!                     F.R. Division working paper No 162.
!                              Mark H. Mawson
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
        INTEGER             i           !} Loop
        INTEGER             j           !} indices.

        REAL                Recip6             !} Useful local scalars.
        REAL                ZIminus            !}
        REAL                ZI                 !}
        REAL                ZIplus             !}
        REAL                ZIplus2            !}
        REAL                phi_3              !}
        REAL                phi_2              !}
        REAL                phi                !}
        REAL                lambda_3           !}
        REAL                lambda_2           !}
        REAL                lambda             !}
        REAL                Coeffminus         !}
        REAL                Coeffzero          !}
        REAL                Coeffplus          !}
        REAL                Coeffplus2         !}
        REAL                CoeffLminus        !}
        REAL                CoeffLzero         !}
        REAL                CoeffLplus         !}
        REAL                CoeffLplus2        !}

!- End of header -------------------------------------------------------

!-----------------------------------------------------------------------
!  1.0  Set up useful scalar and loop over position.
!-----------------------------------------------------------------------

        Recip6 = 1.0 / 6.0

        DO j = 1, Len2Out
         DO i = 1, Len1Out

! ----------------------------------------------------------------------
!  2.0   Perform cubic interpolation in j direction.
! ----------------------------------------------------------------------

       phi   = WtPhi (i,j)
       phi_2 = phi * phi
       phi_3 = phi * phi_2

       Coeffplus2 = Recip6 * ( phi_3 - phi )

       Coeffplus  = 0.5 * ( phi_3  - phi_2  - 2.0*phi )

       Coeffzero  = 0.5 * ( phi_3  - 2.0*phi_2 - phi + 2.0 )

       Coeffminus = Recip6 * ( phi_3 - 3.0*phi_2 + 2.0*phi )

       ZIminus = Coeffplus2 *  DataExt (IOut(i,j) - 1, JOut(i,j) + 2)
     &         - Coeffplus  *  DataExt (IOut(i,j) - 1, JOut(i,j) + 1)
     &         + Coeffzero  *  DataExt (IOut(i,j) - 1, JOut(i,j) )
     &         - Coeffminus *  DataExt (IOut(i,j) - 1, JOut(i,j) - 1)

       ZI      = Coeffplus2 *  DataExt (IOut(i,j), JOut(i,j) + 2)
     &         - Coeffplus  *  DataExt (IOut(i,j), JOut(i,j) + 1)
     &         + Coeffzero  *  DataExt (IOut(i,j), JOut(i,j) )
     &         - Coeffminus *  DataExt (IOut(i,j), JOut(i,j) - 1)

       ZIplus  = Coeffplus2 *  DataExt (IOut(i,j) + 1, JOut(i,j) + 2)
     &         - Coeffplus  *  DataExt (IOut(i,j) + 1, JOut(i,j) + 1)
     &         + Coeffzero  *  DataExt (IOut(i,j) + 1, JOut(i,j) )
     &         - Coeffminus *  DataExt (IOut(i,j) + 1, JOut(i,j) - 1)

       ZIplus2 = Coeffplus2 *  DataExt (IOut(i,j) + 2, JOut(i,j) + 2)
     &         - Coeffplus  *  DataExt (IOut(i,j) + 2, JOut(i,j) + 1)
     &         + Coeffzero  *  DataExt (IOut(i,j) + 2, JOut(i,j) )
     &         - Coeffminus *  DataExt (IOut(i,j) + 2, JOut(i,j) - 1)

!-----------------------------------------------------------------------
!  3.0  Interpolate in i direction and calculate final answer.
!-----------------------------------------------------------------------

       lambda   = WtLambda (i,j)
       lambda_2 = lambda * lambda
       lambda_3 = lambda * lambda_2

       CoeffLplus2 = Recip6 * ( lambda_3  - lambda  )

       CoeffLplus  = - 0.5 * ( lambda_3 - lambda_2 - 2.0*lambda )

       CoeffLzero  = 0.5 * ( lambda_3 - 2.0*lambda_2 - lambda + 2.0 )

       CoeffLminus = - Recip6 *( lambda_3 - 3.0*lambda_2 + 2.0*lambda )

       DataOut (i,j) = CoeffLplus2 *   ZIplus2
     &                     + CoeffLplus  *   ZIplus
     &               + CoeffLzero  *   ZI
     &               + CoeffLminus *   ZIminus


         END DO       ! Close i loop.
        END DO        ! Close j loop.

! End of routine.
      RETURN
      END
