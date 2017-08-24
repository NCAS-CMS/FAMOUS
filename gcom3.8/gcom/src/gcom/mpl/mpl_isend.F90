! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file Copyright
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Subroutine MPL_Isend()
!
!     ******************************************************************
!     * Purpose:
!     *
!     *  Immediate Send
!     *
!     *  Output:  error, request
!     *
!     ******************************************************************

Subroutine MPL_Isend (buf, count, datatype, dest, tag, comm, request, error)

use mpl, Only :                               &
    GC_Int_Kind,  GC_Log_Kind,  GC_Real_Kind, &
    MPL_Int_Kind, MPL_Log_Kind, MPL_Real_Kind

Implicit None

! Arguments and Variables at Model/GCOM precision level
Integer (KIND=GC_INT_KIND) :: &
  buf(*),                     &
  count,                      &
  datatype,                   &
  dest,                       &
  tag,                        &
  comm,                       &
  request,                    &
  error                       
  

! Arguments and Variables at Internal/MPI Library precision level
Integer (KIND=MPL_INT_KIND) ::  &
  l_count,                      & 
  l_dest,                       &
  l_tag,                        & 
  l_datatype,                   &
  l_comm,                       &
  l_request,                    &
  l_error

!=======================================================================

l_count    = count
l_dest     = dest
l_tag      = tag
l_datatype = datatype
l_comm     = comm

Call MPI_Isend(buf, l_count, l_datatype, l_dest, l_tag, l_comm, &
               l_request, l_error)

request = l_request
error   = l_error

Return
End Subroutine MPL_Isend
