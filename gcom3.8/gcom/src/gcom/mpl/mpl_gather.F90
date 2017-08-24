! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file Copyright
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Subroutine MPL_Gather()
!
!     ******************************************************************
!     * Purpose:
!     *
!     *  Gather from many processors to 1
!     *
!     *  Output:  error, recvbuf
!     *
!     ******************************************************************

Subroutine MPL_Gather(sendbuf, sendcnt, sendtype,          &
                      recvbuf, recvcnt, recvtype,          &
                      root,    comm,    error)

USE mpl, ONLY :                               &
    GC_Int_Kind,  GC_Log_Kind,  GC_Real_Kind, &
    MPL_Int_Kind, MPL_Log_Kind, MPL_Real_Kind

Implicit None

! Arguments and Variables at Model/GCOM precision level
Integer (KIND=GC_INT_KIND) :: &
  sendbuf(*),                 &
  recvbuf(*),                 &
  recvcnt,                    &
  sendcnt,                    &
  sendtype,                   &
  recvtype,                   &
  root,                       &
  comm,                       &
  error                       

! Arguments and Variables at Internal/MPI Library precision level
Integer (KIND=MPL_INT_KIND) ::  &
  l_sendcnt,                    & 
  l_sendtype,                   &
  l_recvtype,                   &
  l_recvcnt,                    &
  l_root,                       &
  l_comm,                       &
  l_error

!=======================================================================

l_root     = root
l_comm     = comm
l_recvtype = recvtype
l_sendtype = sendtype
l_sendcnt  = sendcnt
l_recvcnt = recvcnt

Call MPI_Gather(sendbuf, l_sendcnt, l_sendtype,            &
                recvbuf, l_recvcnt, l_recvtype,            &
                l_root,  l_comm,    l_error)

error   = l_error

Return
End Subroutine MPL_Gather
