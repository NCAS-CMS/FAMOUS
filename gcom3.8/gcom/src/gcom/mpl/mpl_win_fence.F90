! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file Copyright
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Subroutine MPL_Win_Fence()
!
!     ******************************************************************
!     * Purpose:
!     *
!     *  Syncronises RMA calls
!     *
!     *  Output:  error
!     *
!     ******************************************************************

Subroutine MPL_Win_Fence (assert, win, error)

Use mpl, Only :                               &
    GC_Int_Kind,  GC_Log_Kind,  GC_Real_Kind, &
    MPL_Int_Kind, MPL_Log_Kind, MPL_Real_Kind

Implicit None

! Arguments and Variables at Model/GCOM precision level
Integer (KIND=GC_INT_KIND) :: &
  assert,                     & 
  win,                        &
  error                       

! Arguments and Variables at Internal/MPI Library precision level
Integer (KIND=MPL_INT_KIND) ::  &
  l_assert,                     &
  l_win,                        &
  l_error                        


!=======================================================================

l_assert = assert
l_win    = win

Call MPI_Win_Fence(l_assert, l_win, l_error)

error = l_error

Return
End Subroutine MPL_Win_Fence
