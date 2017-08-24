! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file Copyright
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Subroutine MPL_Sendrecv()
!
!     ******************************************************************
!     * Purpose:
!     *
!     *  Send and recv in one call
!     *
!     *  Output:  error
!     *
!     ******************************************************************

Subroutine MPL_Sendrecv (sendbuf, sendcount, sendtype, dest,   sendtag, &
                         recvbuf, recvcount, recvtype, source, recvtag, &
                         comm, status, error)

use mpl, Only :                                &
    GC_Int_Kind,  GC_Log_Kind,  GC_Real_Kind,  &
    MPL_Int_Kind, MPL_Log_Kind, MPL_Real_Kind, &
    MPL_Status_Size

Implicit None

! Arguments and Variables at Model/GCOM precision level
Integer (KIND=GC_INT_KIND) :: &
  sendbuf(*),                 &
  sendcount,                  &
  sendtype,                   &
  dest,                       &
  sendtag,                    &
  recvbuf(*),                 &
  recvcount,                  &
  recvtype,                   &
  source,                     &
  recvtag,                    &
  comm,                       &
  status(MPL_STATUS_SIZE),    &
  error                       
  

! Arguments and Variables at Internal/MPI Library precision level
Integer (KIND=MPL_INT_KIND) ::  &
  l_sendcount,                  & 
  l_sendtype,                   &
  l_dest,                       &
  l_sendtag,                    & 
  l_recvcount,                  &
  l_recvtype,                   &
  l_source,                     &
  l_recvtag,                    &
  l_comm,                       &
  l_status(MPL_STATUS_SIZE),    &
  l_error

!=======================================================================

l_sendcount    = sendcount
l_sendtype     = sendtype
l_dest         = dest
l_sendtag      = sendtag
l_recvcount    = recvcount
l_recvtype     = recvtype
l_source       = source
l_recvtag      = recvtag
l_comm         = comm

Call MPI_Sendrecv(sendbuf, l_sendcount, l_sendtype, l_dest,   l_sendtag, &
                  recvbuf, l_recvcount, l_recvtype, l_source, l_recvtag, &
                  l_comm, l_status, l_error)

status(:) = l_status(:)
error     = l_error

Return
End Subroutine MPL_Sendrecv
