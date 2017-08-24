! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file Copyright
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Subroutine MPL_Comm_Group()
!
!     ******************************************************************
!     * Purpose:
!     *
!     *  Accesses group with given communicator
!     *
!     *  Output:  error
!     *
!     ******************************************************************

Subroutine MPL_Comm_Group (comm, group, error)

use mpl, Only :                               &
    GC_Int_Kind,  GC_Log_Kind,  GC_Real_Kind, &
    MPL_Int_Kind, MPL_Log_Kind, MPL_Real_Kind

Implicit None

! Arguments and Variables at Model/GCOM precision level
Integer (KIND=GC_INT_KIND) :: &
  comm,                       & 
  group,                      &
  error                       
  

! Arguments and Variables at Internal/MPI Library precision level
Integer (KIND=MPL_INT_KIND) ::  &
  l_comm,                       &
  l_group,                      &
  l_error                        

!=======================================================================

l_comm = comm

Call MPI_Comm_Group(l_comm,l_group,l_error)

group = l_group
error = l_error

Return
End Subroutine MPL_Comm_Group
