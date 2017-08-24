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
      SUBROUTINE HorizontalInterpQuintic
     &              (LowerBound,
     &               Len1In,  Len2In,
     &               Len1Out, Len2Out,
     &               DataExt,
     &               WtLambda,
     &               WtPhi,
     &               IOut,
     &               JOut,
     &               DataOut)

! Description: Performs quintic Lagrange interpolation of a 2-d field to
!              a 2-d set of points defined by IOut, JOut and WtLambda,
!              WtPhi.
!
! Method: This is a modified version of the routine quintic_lagrange
!         written by Sue Coulter and described in:
!
!           The proposed semi-Lagrangian advection scheme for the
!              semi-Implicit Unified Model integration scheme.
!                     F.R. Division working paper No 162.
!                              Mark H. Mawson
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
        INTEGER     IOut (Len1Out,Len2Out)  ! Point such that
        INTEGER     JOut (Len1Out,Len2Out)  ! the desired output point
!                                           ! lies between it and it+1.
        REAL   DataExt(LowerBound:Len1In+1-LowerBound,
     &          LowerBound:Len2In+1-LowerBound)  ! Data interpolated
        REAL   WtLambda (Len1Out,Len2Out)   ! A number between 0 & 1.
        REAL   WtPhi (Len1Out,Len2Out)      ! A number between 0 & 1.

! Array  arguments with INTENT(out):
        REAL    DataOut (Len1Out,Len2Out)   ! Data interpolated to
!                                           ! desired locations.
!* End of Subroutine arguments

! Local scalars:
        INTEGER             i                ! } Loop
        INTEGER             j                ! } indices.

        REAL                Recip6           ! } A multitude of useful
        REAL                Recip12          ! } local scalars...
        REAL                Recip24          ! }
        REAL                Recip120         ! }
        REAL                ZIminus          ! }
        REAL                ZIminus2         ! }
        REAL                ZI               ! }
        REAL                ZIplus           ! }
        REAL                ZIplus2          ! }
        REAL                ZIplus3          ! }
        REAL                phi              ! }
        REAL                phi_3            ! }
        REAL                phi_2            ! }
        REAL                phi_4            ! }
        REAL                phi_5            ! }
        REAL                lambda_3         ! }
        REAL                lambda_2         ! }
        REAL                lambda           ! }
        REAL                lambda_4         ! }
        REAL                lambda_5         ! }
        REAL                Coeffminus2      ! }
        REAL                Coeffminus       ! }
        REAL                Coeffzero        ! }
        REAL                Coeffplus        ! }
        REAL                Coeffplus2       ! }
        REAL                Coeffplus3       ! }
        REAL                CoeffLminus2     ! }
        REAL                CoeffLminus      ! }
        REAL                CoeffLzero       ! }
        REAL                CoeffLplus       ! }
        REAL                CoeffLplus2      ! }
        REAL                CoeffLplus3      ! }

!- End of header -------------------------------------------------------

!-----------------------------------------------------------------------
!  1.0  Set up some useful scalars.
!-----------------------------------------------------------------------

        Recip6   = 1.0 / 6.0
        Recip12  = 1.0 / 12.0
        Recip24  = 1.0 / 24.0
        Recip120 = 1.0 / 120.0


        DO j = 1, Len2Out
         DO i = 1, Len1Out

!-----------------------------------------------------------------------
!  2.0   Perform quintic Lagrange interpolation in j direction
!-----------------------------------------------------------------------

       phi    = WtPhi (i,j)
       phi_2  = phi * phi
       phi_3  = phi_2 * phi
       phi_4  = phi_2 * phi_2
       phi_5  = phi_2 * phi_3

       Coeffplus3  = Recip120 * ( phi_5 - 5.0 * phi_3 + 4.0 * phi )

       Coeffplus2  = - Recip24 *
     &             ( phi_5 - phi_4 - 7.0 * phi_3 + phi_2 + 6.0 * phi )

       Coeffplus   = Recip12 *
     &             ( phi_5 - 2.0 * phi_4 - 7.0 * phi_3 +
     &               8.0 * phi_2  + 12.0 * phi )

       Coeffzero   = - Recip12 *
     &             ( phi_5 - 3.0 * phi_4 - 5.0 * phi_3 +
     &               15.0 * phi_2 + 4.0 * phi -12.0 )

       Coeffminus  =  Recip24 *
     &             ( phi_5 - 4.0 * phi_4 - phi_3 +
     &               16.0 * phi_2 - 12.0 * phi )

       Coeffminus2 = - Recip120 *
     &             (  phi_5 - 5.0 * phi_4 + 5.0 * phi_3 +
     &                5.0 * phi_2  - 6.0 * phi )

       ZIminus2 = Coeffplus3 *
     &                DataExt (IOut(i,j) - 2, JOut(i,j) + 3)
     &         + Coeffplus2 *
     &                DataExt (IOut(i,j) - 2, JOut(i,j) + 2)
     &         + Coeffplus *
     &                DataExt (IOut(i,j) - 2, JOut(i,j) + 1)
     &         + Coeffzero *
     &                DataExt (IOut(i,j) - 2, JOut(i,j))
     &         + Coeffminus *
     &                DataExt (IOut(i,j) - 2, JOut(i,j) - 1)
     &         + Coeffminus2 *
     &                DataExt (IOut(i,j) - 2, JOut(i,j) - 2)

       ZIminus = Coeffplus3 *
     &                DataExt (IOut(i,j) - 1, JOut(i,j) + 3)
     &       + Coeffplus2 *
     &                DataExt (IOut(i,j) - 1, JOut(i,j) + 2)
     &       + Coeffplus *
     &                DataExt (IOut(i,j) - 1, JOut(i,j) + 1)
     &       + Coeffzero *
     &                DataExt (IOut(i,j) - 1, JOut(i,j))
     &       + Coeffminus *
     &                DataExt (IOut(i,j) - 1, JOut(i,j) - 1)
     &       + Coeffminus2 *
     &                DataExt (IOut(i,j) - 1,JOut(i,j) - 2)

       ZI     = Coeffplus3 *
     &                DataExt (IOut(i,j), JOut(i,j) + 3)
     &       + Coeffplus2 *
     &                DataExt (IOut(i,j), JOut(i,j) + 2)
     &       + Coeffplus *
     &                DataExt (IOut(i,j), JOut(i,j) + 1)
     &       + Coeffzero *
     &                DataExt (IOut(i,j), JOut(i,j))
     &       + Coeffminus *
     &                DataExt (IOut(i,j), JOut(i,j) - 1)
     &       + Coeffminus2 *
     &                DataExt (IOut(i,j), JOut(i,j) - 2)

       ZIplus  = Coeffplus3 *
     &                DataExt (IOut(i,j) + 1, JOut(i,j) + 3)
     &       + Coeffplus2 *
     &                DataExt (IOut(i,j) + 1, JOut(i,j) + 2)
     &       + Coeffplus *
     &                DataExt (IOut(i,j) + 1, JOut(i,j) + 1)
     &       + Coeffzero *
     &                DataExt (IOut(i,j) + 1, JOut(i,j))
     &       + Coeffminus *
     &                DataExt (IOut(i,j) + 1, JOut(i,j) - 1)
     &       + Coeffminus2 *
     &                DataExt (IOut(i,j) + 1, JOut(i,j) - 2)

       ZIplus2 =  Coeffplus3 *
     &                DataExt (IOut(i,j) + 2, JOut(i,j) + 3)
     &        +  Coeffplus2 *
     &                DataExt (IOut(i,j) + 2, JOut(i,j) + 2)
     &        +  Coeffplus *
     &                DataExt (IOut(i,j) + 2, JOut(i,j) + 1)
     &        +  Coeffzero *
     &                DataExt (IOut(i,j) + 2,JOut(i,j))
     &        +  Coeffminus *
     &                DataExt (IOut(i,j) + 2, JOut(i,j) - 1)
     &        +  Coeffminus2 *
     &                DataExt (IOut(i,j) + 2,JOut(i,j) - 2)


       ZIplus3 =  Coeffplus3 *
     &                DataExt (IOut(i,j) + 3, JOut(i,j) + 3)
     &        +  Coeffplus2 *
     &                DataExt (IOut(i,j) + 3, JOut(i,j) + 2)
     &        +  Coeffplus *
     &                DataExt (IOut(i,j) + 3, JOut(i,j) + 1)
     &        +  Coeffzero *
     &                DataExt (IOut(i,j) + 3, JOut(i,j))
     &        +  Coeffminus *
     &                DataExt (IOut(i,j) + 3, JOut(i,j) - 1)
     &        +  Coeffminus2 *
     &                DataExt (IOut(i,j) + 3, JOut(i,j) - 2)

!-----------------------------------------------------------------------
!  3.0  Interpolate in i direction and calculate final answer.
!-----------------------------------------------------------------------

       lambda    = WtLambda (i,j)
       lambda_2  = lambda * lambda
       lambda_3  = lambda_2 * lambda
       lambda_4  = lambda_2 * lambda_2
       lambda_5  = lambda_2 * lambda_3

       CoeffLplus3 = Recip120 *
     &                   ( lambda_5 - 5.0 * lambda_3 + 4.0 * lambda )

       CoeffLplus2 = - Recip24 *
     &             ( lambda_5 - lambda_4 - 7.0 * lambda_3 +
     &               lambda_2  + 6.0 * lambda )

       CoeffLplus  = Recip12 *
     &             ( lambda_5 - 2.0 * lambda_4 - 7.0 * lambda_3 +
     &               8.0 * lambda_2 + 12.0 * lambda )

       CoeffLzero  = - Recip12 *
     &             ( lambda_5 - 3.0 * lambda_4 - 5.0 * lambda_3 +
     &               15.0 * lambda_2 + 4.0 * lambda - 12.0 )

       CoeffLminus =  Recip24 *
     &             ( lambda_5 - 4.0 * lambda_4 - lambda_3 +
     &                16.0 * lambda_2 - 12.0 * lambda )

       CoeffLminus2 = - Recip120 *
     &             ( lambda_5 - 5.0 * lambda_4 + 5.0 * lambda_3 +
     &               5.0 *lambda_2  - 6.0 * lambda )


       DataOut (i,j) =   CoeffLplus3  *  ZIplus3
     &               +  CoeffLplus2  *  ZIplus2
     &               +  CoeffLplus   *  ZIplus
     &               +  CoeffLzero   *  ZI
     &               +  CoeffLminus  *  ZIminus
     &               +  CoeffLminus2 *  ZIminus2

         END DO                               ! Close i loop.
        END DO                                ! Close j loop.


! End of routine.
      RETURN
      END
