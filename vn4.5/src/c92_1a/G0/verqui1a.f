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
        Subroutine VerticalInterpQuintic
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
!   Quintic Lagrange interpolation of DataIn at RIn
!   to DataOut at ROut in the vertical
!   i.e. 3rd dimension of the arrays.
!
! Method:
!   Simplified version of Mark Mawson's general interpolation routine
!   Quintic Lagrange
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
! Scalar arguments with :intent(in)
        Integer   len1In    ! Extent of DataIn in i direction.
        Integer   len2In    ! Extent of DataIn in j direction.
        Integer   len3In    ! Extent of DataIn in k direction.
        Integer   len3Out   ! Extent of DataOut in k direction.

! Array  arguments with :intent(in)
        Integer   KOut (len1In,len2In,len3Out) ! level in RIn below ROut

        Real  RIn  (len1In,len2In,len3In)  ! Vertical coordinate input

        Real  ROut (len1In,len2In,len3Out)  ! Vertical coordinate output

        Real  DataIn (len1In,len2In,len3In) ! Data on the original
!                                           ! vertical levels.

! Array  arguments with :intent(out)
        Real     DataOut (len1In,len2In,len3Out) ! Data interpolated to
!                                                ! new vertical levels.

! Local scalars:
        Integer  i, j, k                       ! Loop indices

        Real     rhereminus2                 ! } Temporary data stores
        Real     rhereminus                  ! } elements of the RIn
        Real     rhere                       ! } array
        Real     rhereplus                   ! }
        Real     rhereplus2                  ! }
        Real     rhereplus3                  ! }
        Real     Z1, Z2, Z3, Z4                ! } Stores for Partial
        Real     Z1A, Z2A, Z3A, Z4A, Z5A, Z6A  ! } sums
        Real     Z1B, Z2B, Z3B, Z4B, Z5B, Z6B  ! }

!- End of header

        Do k = 1, len3Out
         Do j = 1, len2In
          Do i = 1, len1In
           If (KOut(i,j,k) .eq. 1 .or.
     &       KOut(i,j,k) .eq. len3In-1) then
! use the linear method to interpolate between bottom most levels
! or top most levels.
           rhere           = RIn (i,j,KOut(i,j,k))
           rhereplus      = RIn (i,j,KOut(i,j,k) + 1)

             DataOut(i,j,k) = ( ( (ROut(i,j,k) - rhere)
     &                         * DataIn (i,j,KOut(i,j,k) + 1) )
     &                      - ( (ROut(i,j,k) - rhereplus)
     &                         * DataIn (i,j,KOut(i,j,k))     )
     &                      )  / ( rhereplus - rhere )

         Else  if ( KOut(i,j,k) .eq. 2 .or.
     &              KOut(i,j,k) .eq. len3In-2) then
! use cubic interpolation.
          rhereminus     = RIn (i,j,KOut(i,j,k) - 1)
          rhere           = RIn (i,j,KOut(i,j,k))
          rhereplus      = RIn (i,j,KOut(i,j,k) + 1)
          rhereplus2     = RIn (i,j,KOut(i,j,k) + 2)

           Z1 =  ( (ROut(i,j,k) - rhere)
     &           * (ROut(i,j,k) - rhereplus)
     &           * (ROut(i,j,k) - rhereplus2)
     &           *  DataIn (i,j,KOut(i,j,k)-1)  )
     &           /
     &           ( (rhereminus - rhere)
     &           * (rhereminus - rhereplus)
     &           * (rhereminus - rhereplus2)  )

           Z2 =  ( (ROut(i,j,k) - rhereminus)
     &           * (ROut(i,j,k) - rhereplus)
     &           * (ROut(i,j,k) - rhereplus2)
     &           * DataIn (i,j,KOut(i,j,k))     )
     &           /
     &           ( (rhere - rhereminus)
     &           * (rhere - rhereplus)
     &           * (rhere - rhereplus2)        )

           Z3 =  ( (ROut(i,j,k) - rhereminus)
     &           * (ROut(i,j,k) - rhere)
     &           * (ROut(i,j,k) - rhereplus2)
     &           * DataIn (i,j,KOut(i,j,k)+1)   )
     &            /
     &            ( (rhereplus - rhereminus)
     &           * (rhereplus - rhere)
     &           * (rhereplus - rhereplus2)   )

           Z4 =  ( (ROut(i,j,k) - rhereminus)
     &           * (ROut(i,j,k) - rhere)
     &           * (ROut(i,j,k) - rhereplus)
     &           * DataIn (i,j,KOut(i,j,k)+2)   )
     &           /
     &           ( (rhereplus2 - rhereminus)
     &           * (rhereplus2 - rhere)
     &           * (rhereplus2 - rhereplus)   )

           DataOut (i,j,k) = Z1+Z2+Z3+Z4

         Else ! use quintic interpolation
          rhereminus2    = RIn (i,j,KOut(i,j,k) - 2)
          rhereminus     = RIn (i,j,KOut(i,j,k) - 1)
          rhere           = RIn (i,j,KOut(i,j,k))
          rhereplus      = RIn (i,j,KOut(i,j,k) + 1)
          rhereplus2     = RIn (i,j,KOut(i,j,k) + 2)
          rhereplus3     = RIn (i,j,KOut(i,j,k) + 3)

          Z1A =     (ROut(i,j,k) - rhereminus)
     &            * (ROut(i,j,k) - rhere)
     &            * (ROut(i,j,k) - rhereplus)
     &            * (ROut(i,j,k) - rhereplus2)
     &            * (ROut(i,j,k) - rhereplus3)
     &            * DataIn (i,j,KOut(i,j,k)-2)

          Z1B =     (rhereminus2 - rhereminus)
     &            * (rhereminus2 - rhere)
     &            * (rhereminus2 - rhereplus)
     &            * (rhereminus2 - rhereplus2)
     &            * (rhereminus2 - rhereplus3)

          Z2A =     (ROut(i,j,k) - rhereminus2)
     &            * (ROut(i,j,k) - rhere)
     &            * (ROut(i,j,k) - rhereplus)
     &            * (ROut(i,j,k) - rhereplus2)
     &            * (ROut(i,j,k) - rhereplus3)
     &            * DataIn (i,j,KOut(i,j,k)-1)

          Z2B =     (rhereminus - rhereminus2)
     &            * (rhereminus - rhere)
     &            * (rhereminus - rhereplus)
     &            * (rhereminus - rhereplus2)
     &            * (rhereminus - rhereplus3)

          Z3A =     (ROut(i,j,k) - rhereminus2)
     &            * (ROut(i,j,k) - rhereminus)
     &            * (ROut(i,j,k) - rhereplus)
     &            * (ROut(i,j,k) - rhereplus2)
     &            * (ROut(i,j,k) - rhereplus3)
     &            *  DataIn (i,j,KOut(i,j,k))

          Z3B =     (rhere - rhereminus2)
     &            * (rhere - rhereminus)
     &            * (rhere - rhereplus)
     &            * (rhere - rhereplus2)
     &            * (rhere - rhereplus3)

          Z4A =     (ROut(i,j,k) - rhereminus2)
     &            * (ROut(i,j,k) - rhereminus)
     &            * (ROut(i,j,k) - rhere)
     &            * (ROut(i,j,k) - rhereplus2)
     &            * (ROut(i,j,k) - rhereplus3)
     &            * DataIn (i,j,KOut(i,j,k)+1)

          Z4B =     (rhereplus - rhereminus2)
     &            * (rhereplus - rhereminus)
     &            * (rhereplus - rhere)
     &            * (rhereplus - rhereplus2)
     &            * (rhereplus - rhereplus3)

          Z5A =     (ROut(i,j,k) - rhereminus2)
     &            * (ROut(i,j,k) - rhereminus)
     &            * (ROut(i,j,k) - rhere)
     &            * (ROut(i,j,k) - rhereplus)
     &            * (ROut(i,j,k) - rhereplus3)
     &            * DataIn (i,j,KOut(i,j,k)+2)

          Z5B =     (rhereplus2 - rhereminus2)
     &            * (rhereplus2 - rhereminus)
     &            * (rhereplus2 - rhere)
     &        * (rhereplus2 - rhereplus)
     &            * (rhereplus2 - rhereplus3)

          Z6A =     (ROut(i,j,k) - rhereminus2)
     &            * (ROut(i,j,k) - rhereminus)
     &            * (ROut(i,j,k) - rhere)
     &            * (ROut(i,j,k) - rhereplus)
     &            * (ROut(i,j,k) - rhereplus2)
     &            * DataIn (i,j,KOut(i,j,k)+3)

          Z6B =     (rhereplus3 - rhereminus2)
     &            * (rhereplus3 - rhereminus)
     &            * (rhereplus3 - rhere)
     &            * (rhereplus3 - rhereplus)
     &            * (rhereplus3 - rhereplus2)

          DataOut (i,j,k) = Z1A/Z1B+Z2A/Z2B+Z3A/Z3B
     &                     +Z4A/Z4B+Z5A/Z5B+Z6A/Z6B

           End If
          End Do
         End Do
        End Do

        Return
        End
