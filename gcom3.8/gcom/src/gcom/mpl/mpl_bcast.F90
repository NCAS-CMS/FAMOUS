! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file Copyright
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Subroutine MPL_BCast ()
!
!     ******************************************************************
!     * Purpose:
!     *
!     *  Broadcast data
!     *
!     * Output: ierr
!     *  
!     ******************************************************************

Subroutine MPL_BCast (buffer, count, datatype, root, comm, ierr )

use mpl, Only :                               &
    GC_Int_Kind,  GC_Log_Kind,  GC_Real_Kind, &
    MPL_Int_Kind, MPL_Log_Kind, MPL_Real_Kind

Implicit None

! Arguments and Variables at Model/GCOM precision level
Integer (KIND=GC_INT_KIND) ::      &
  buffer(*),                       &
  count,                           &
  datatype,                        &
  root,                            &
  comm,                            &
  ierr

! Arguments and Variables at Internal/MPI Library precision level
Integer (KIND=MPL_INT_KIND)  ::    &
  l_count,                         &
  l_datatype,                      &
  l_root,                          &
  l_comm,                          &
  l_ierr

!Cast Model precision input arguments to MPI library precision.
l_count     = count
l_datatype  = datatype
l_root      = root
l_comm      = comm

!Call MPI routine
Call MPI_Bcast(buffer, l_count, l_datatype, l_root, l_comm, l_ierr)

!Cast returned arguments to Model Preciison
ierr = l_ierr

Return
End Subroutine MPL_BCast
