! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Subroutine MPL_Iprobe()
!
!     ******************************************************************
!     * Purpose:
!     *
!     *  Returns flag=true if there is a message that can be received
!     *  that matches the pattern specified by the arguments source,
!     *  tag and comm. 
!     *
!     *  Output:  error, request
!     *
!     ******************************************************************

Subroutine MPL_Iprobe (source, tag, comm, flag, istatus, error)

use mpl, Only :                                &
    GC_Int_Kind,  GC_Log_Kind,                 &
    MPL_Int_Kind, MPL_Log_Kind,                &
    MPL_STATUS_SIZE

Implicit None

! Arguments and Variables at Model/GCOM precision level
Integer (KIND=GC_INT_KIND) :: &
  source,                     &
  tag,                        &
  comm,                       &
  istatus(MPL_STATUS_SIZE),   &
  error                       
Logical (KIND=GC_LOG_KIND) :: flag

! Arguments and Variables at Internal/MPI Library precision level
Integer (KIND=MPL_INT_KIND) ::  &
  l_source,                     &
  l_tag,                        & 
  l_comm,                       &
  l_status(MPL_STATUS_SIZE),    &
  l_error
Logical (KIND=MPL_LOG_KIND) :: l_flag

!=======================================================================

l_source   = source
l_tag      = tag
l_comm     = comm

Call MPI_Iprobe(l_source, l_tag, l_comm, l_flag, l_status, l_error)

flag       = l_flag
istatus(:) = l_status(:)
error      = l_error

Return
End Subroutine MPL_Iprobe

