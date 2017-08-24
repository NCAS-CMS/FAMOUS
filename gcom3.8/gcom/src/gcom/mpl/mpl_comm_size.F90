! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file Copyright
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Subroutine MPL_Comm_Size ()
!
!     ******************************************************************
!     * Purpose:
!     *
!     *  Return the size of the comunicating group.
!     *
!     * Input:
!     *  comm - the comunicator(group size) for which the size is
!     *         requested
!     *
!     * Output:
!     *  size - the size of the group
!     *  err  - an error code
!     *
!     ******************************************************************

Subroutine MPL_Comm_Size (comm, size, ierr)

use mpl, Only :                               &
    GC_Int_Kind,  GC_Log_Kind,  GC_Real_Kind, &
    MPL_Int_Kind, MPL_Log_Kind, MPL_Real_Kind

Implicit None

! Arguments and Variables at Model/GCOM precision level
Integer (KIND=GC_INT_KIND) ::      &
  comm,                            &
  size,                            &
  ierr

! Arguments and Variables at Internal/MPI Library precision level
Integer (KIND=MPL_INT_KIND)  ::    &
  l_comm,                          &
  l_size,                          &
  l_ierr

!Cast Model precision input arguments to MPI library precision.
l_comm = comm

!Call MPI routine
Call MPI_Comm_Size(l_comm, l_size, l_ierr)

!Cast returned arguments to Model Preciison
ierr = l_ierr
size = l_size

Return
End Subroutine MPL_Comm_Size
