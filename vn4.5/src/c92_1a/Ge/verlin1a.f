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
        Subroutine VerticalInterpLinear
     &      ( len1In,
     &          len2In,
     &          len3In,
     &          len3Out,
     &          KOut,
     &          RIn,
     &          ROut,
     &          DataIn,
     &          DataOut )

!   Description
!   Linearly interpolates DataIn at RIn
!   to DataOut at ROut in the vertical
!   i.e. 3rd dimension of the arrays.
!
! Method:
!   Simplified version of Mark Mawson's general interpolation routine
!   TriLinear.
!                The proposed semi-Lagrangian advection scheme for the
!                   semi-Implicit Unified Model integration scheme.
!                         F.R. Division working paper No 162.
!
! Current Code Owner: Stuart Bell
!
! History:
! Version   Date     Comment
! -------   ----     -------
! 4.0      11/7/95   Equiv. to VAR code as at time of build:Stuart Bell
!
! Code Description:
!   Language:           Fortran 77 plus
!   Software Standards: "UM and Met O standards".
!
! Declarations:

        Implicit none

! Subroutine arguments
! Scalar arguments with : intent(in)
        Integer   len1In    ! Extent of DataIn in i direction.
        Integer   len2In    ! Extent of DataIn in j direction.
        Integer   len3In    ! Extent of DataIn in k direction.
        Integer   len3Out   ! Extent of DataOut in k direction.

! Array  arguments with : intent(in)
        Integer   KOut (len1In,len2In,len3Out) ! level in RIn below ROut

        Real  RIn  (len1In,len2In,len3In)  ! Vertical coordinate input

        Real  ROut (len1In,len2In,len3Out)  ! Vertical coordinate output

        Real  DataIn (len1In,len2In,len3In) ! Data on the original
!                                           ! vertical levels.

! Array  arguments with : intent(out)
        Real     DataOut (len1In,len2In,len3Out) ! Data interpolated to
!                                                ! new vertical levels.
! Local scalars:
        Integer                    i, j, k       ! Loop indices
        Real                 rhere       ! } Temporary data stores for
        Real                 rhereplus   ! } elements of the RIn array

!- End of header

         Do k = 1, len3Out
          Do j = 1, len2In
           Do i = 1, len1In
             rhere           = RIn (i,j,KOut(i,j,k))
             rhereplus      = RIn (i,j,KOut(i,j,k) + 1)

             DataOut(i,j,k) = ( ( (ROut(i,j,k) - rhere)
     &                         * DataIn (i,j,KOut(i,j,k) + 1) )
     &                      - ( (ROut(i,j,k) - rhereplus)
     &                         * DataIn (i,j,KOut(i,j,k))     )
     &                      )  / ( rhereplus - rhere )

          End Do
         End Do
        End Do

        Return
        End
