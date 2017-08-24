! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file Copyright
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Subroutine MPL_Send()
!
!     ******************************************************************
!     * Purpose:
!     *
!     *  Blocking Send
!     *
!     *  Output:  error, request
!     *
!     ******************************************************************

Subroutine MPL_Send (buf, count, datatype, dest, tag, comm, error)

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
  error                       
  

! Arguments and Variables at Internal/MPI Library precision level
Integer (KIND=MPL_INT_KIND) ::  &
  l_count,                      & 
  l_dest,                       &
  l_tag,                        & 
  l_datatype,                   &
  l_comm,                       &
  l_error

!=======================================================================

l_count    = count
l_dest     = dest
l_tag      = tag
l_datatype = datatype
l_comm     = comm

Call MPI_Send(buf, l_count, l_datatype, l_dest, l_tag, l_comm, &
              l_error)

error   = l_error

Return
End Subroutine MPL_Send
