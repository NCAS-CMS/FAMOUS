! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file Copyright
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Subroutine MPL_Comm_Split()
!
!     ******************************************************************
!     * Purpose:
!     *
!     *  Splits a communicator
!     *
!     *  Output:  error
!     *
!     ******************************************************************

Subroutine MPL_Comm_Split (comm, color, key, newcomm, error)

use mpl, Only :                               &
    GC_Int_Kind,  GC_Log_Kind,  GC_Real_Kind, &
    MPL_Int_Kind, MPL_Log_Kind, MPL_Real_Kind

Implicit None

! Arguments and Variables at Model/GCOM precision level
Integer (KIND=GC_INT_KIND) :: &
  comm,                       & 
  color,                      &
  key,                        &
  newcomm,                    &
  error                       
  

! Arguments and Variables at Internal/MPI Library precision level
Integer (KIND=MPL_INT_KIND) ::  &
  l_comm,                       &
  l_color,                      &
  l_key,                        &
  l_newcomm,                    &
  l_error                        

!=======================================================================

l_comm  = comm
l_color = color
l_key   = key

Call MPI_Comm_Split(l_comm, l_color, l_key, l_newcomm, l_error)

newcomm = l_newcomm
error   = l_error

Return
End Subroutine MPL_Comm_Split
