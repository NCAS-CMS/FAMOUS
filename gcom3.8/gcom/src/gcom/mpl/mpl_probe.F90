! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file Copyright
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Subroutine MPL_Probe()
!
!     ******************************************************************
!     * Purpose:
!     *
!     *  Probe unreceived message to find status
!     *
!     *  Output:  error, request
!     *
!     ******************************************************************

Subroutine MPL_Probe (source, tag, comm, status, error)

use mpl, Only :                                &
    GC_Int_Kind,  GC_Log_Kind,  GC_Real_Kind,  &
    MPL_Int_Kind, MPL_Log_Kind, MPL_Real_Kind, &
    MPL_STATUS_SIZE

Implicit None

! Arguments and Variables at Model/GCOM precision level
Integer (KIND=GC_INT_KIND) :: &
  source,                     &
  tag,                        &
  comm,                       &
  status(MPL_STATUS_SIZE),    &
  error                       
  

! Arguments and Variables at Internal/MPI Library precision level
Integer (KIND=MPL_INT_KIND) ::  &
  l_source,                     &
  l_tag,                        & 
  l_comm,                       &
  l_status(MPL_STATUS_SIZE),    &
  l_error

!=======================================================================

l_source   = source
l_tag      = tag
l_comm     = comm

Call MPI_Probe(l_source, l_tag, l_comm, l_status, l_error)

status(:) = l_status(:)
error     = l_error

Return
End Subroutine MPL_Probe
