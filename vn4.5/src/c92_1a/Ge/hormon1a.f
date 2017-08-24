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
      SUBROUTINE HorizontalInterpMonotone
     &                        (LowerBound,
     &                         Len1In,  Len2In,
     &                         Len1Out, Len2Out,
     &                         DataExt,
     &                         DataHigh,
     &                         DataMono,
     &                         IOut,
     &                         JOut,
     &                         DataOut)

! Description: Ensures monotonicity if high order non-monotone
!              interpolation has been used.
!
! Method: This is a modified version of the routine mono_conserv
!         written by Mark Mawson and described in:
!
!           The proposed semi-Lagrangian advection scheme for the
!              semi-Implicit Unified Model integration scheme.
!                     F.R. Division working paper No 162.
!                              Mark H. Mawson
!
!          and is based on Priestley, 1993 (the full reference
!          for which can be found in the above documentation).
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
        REAL    DataHigh (Len1Out,Len2Out)    ! Data interpolated
!                                             ! by high order scheme
        REAL     DataMono (Len1Out,Len2Out)   ! Data interpolated
!                                             ! by monotone scheme.

! Array  arguments with INTENT(out):
        REAL    DataOut (Len1Out,Len2Out)     ! Data interpolated to
!                                             ! desired locations.

!* End of Subroutine arguments

! Local scalars:
        INTEGER             i                ! } Loop
        INTEGER             j                ! } indices.

        REAL                MaxMono
        REAL                HighlessMono

! Local arrays:
        REAL     MaxAlpha (Len1Out,Len2Out)
        REAL     MinMono (Len1Out,Len2Out)

!- End of header -------------------------------------------------------


!-----------------------------------------------------------------------
!  1.0 Find Max and min values of alpha allowed for montonicity.
!               (Equation 1.8 in Priestley 1993)
!-----------------------------------------------------------------------

        DO j = 1, Len2Out
         DO i = 1, Len1Out

! Find max and min monotone values for the point concerned.

         MaxMono = MAX ( DataExt(IOut(i,j), JOut(i,j)),
     &                 DataExt(IOut(i,j)+1, JOut(i,j)),
     &                 DataExt(IOut(i,j), JOut(i,j)+1),
     &                 DataExt(IOut(i,j)+1, JOut(i,j)+1)   )

         MinMono(i,j) = MIN ( DataExt(IOut(i,j), JOut(i,j)),
     &                      DataExt(IOut(i,j)+1, JOut(i,j)),
     &                      DataExt(IOut(i,j), JOut(i,j)+1),
     &                      DataExt(IOut(i,j)+1, JOut(i,j)+1)  )

         HighlessMono = DataHigh (i,j) - DataMono (i,j)

         MaxAlpha(i,j) = 0.0

         IF (HighlessMono .gt. 0.0) THEN
        MaxAlpha(i,j) = MAX (  0.0,
     &               (MaxMono - DataMono(i,j)) / HighlessMono )

         ELSE IF (HighlessMono .lt. 0.0) THEN
        MaxAlpha(i,j) = MAX (  0.0,
     &               (MinMono(i,j) - DataMono(i,j)) / HighlessMono )
         END IF

         MaxAlpha(i,j) = MIN (1.0, MaxAlpha(i,j))

         END DO
        END DO


!-----------------------------------------------------------------------
!  2.0  Form output data given the alpha values.
!-----------------------------------------------------------------------

        DO j = 1, Len2Out
         DO i = 1, Len1Out

         DataOut(i,j) = (1.0 -  MaxAlpha(i,j)) * DataMono(i,j)
     &              + MaxAlpha(i,j) * DataHigh(i,j)

! ...still need to check value not less then minimum because of
! rounding error problems on the Cray...

         IF (DataOut(i,j) .lt. MinMono(i,j)) DataOut(i,j) = MinMono(i,j)

         END DO
        END DO

! End of routine.
      RETURN
      END
