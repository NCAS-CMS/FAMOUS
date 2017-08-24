! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file Copyright
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Subroutine MPL_Win_Create()
!
!     ******************************************************************
!     * Purpose:
!     *
!     *  Create an MPI Window object for one-sided communication
!     *
!     *  Output:  win, error
!     *
!     ******************************************************************

Subroutine MPL_Win_Create (base, isize, disp_unit, info, comm, &
                           win, error)

Use mpl, Only :                                &
    GC_Int_Kind,  GC_Log_Kind,  GC_Real_Kind,  &
    MPL_Int_Kind, MPL_Log_Kind, MPL_Real_Kind, &
    MPL_ADDRESS_KIND

Implicit None

! Arguments and Variables at Model/GCOM precision level
Integer (KIND=GC_INT_KIND) :: &
  base(*),                    &
  disp_unit,                  &
  info,                       &
  comm,                       & 
  win,                        &
  error                       

! Arguments and Variables at Internal/MPI Library precision level
Integer (KIND=MPL_INT_KIND) ::  &
  l_disp_unit,                  &
  l_info,                       &
  l_comm,                       &
  l_win,                        &
  l_error                        

Integer (KIND=MPL_ADDRESS_KIND) :: &
  isize

!=======================================================================

l_comm      = comm
l_disp_unit = disp_unit
l_info      = info

Call MPI_Win_Create(base, isize, l_disp_unit, l_info, l_comm, &
                    l_win, l_error)

win   = l_win
error = l_error

Return
End Subroutine MPL_Win_Create
