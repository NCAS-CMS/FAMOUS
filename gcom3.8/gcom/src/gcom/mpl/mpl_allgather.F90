! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file Copyright
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Subroutine MPL_Allgather()
!
!     ******************************************************************
!     * Purpose:
!     *
!     *  The block of data sent from the jth process is received by
!     *  every process and placed in the jth block of the receive buffer
!     *
!     *  Output:  error, recvbuf
!     *
!     ******************************************************************

Subroutine MPL_Allgather(sendbuf, sendcount,  sendtype,          &
                         recvbuf, recvcount,  recvtype,          &
                         comm,    error)

use mpl, Only :                               &
    GC_Int_Kind,  GC_Log_Kind,  GC_Real_Kind, &
    MPL_Int_Kind, MPL_Log_Kind, MPL_Real_Kind

Implicit None

! Arguments and Variables at Model/GCOM precision level
Integer (KIND=GC_INT_KIND) :: &
  sendbuf(*),                 &
  recvbuf(*),                 &
  sendcount,                   &
  recvcount,                  &
  sendtype,                   &
  recvtype,                   &
  comm,                       &
  error                       
  

! Arguments and Variables at Internal/MPI Library precision level
Integer (KIND=MPL_INT_KIND) ::  &
  l_sendcount,                  & 
  l_recvcount,                  &
  l_sendtype,                   &
  l_recvtype,                   &
  l_comm,                       &
  l_error

!=======================================================================

l_comm      = comm
l_recvtype  = recvtype
l_sendtype  = sendtype
l_sendcount = sendcount
l_recvcount = recvcount

Call MPI_Allgather(sendbuf, l_sendcount,  l_sendtype,  &
                   recvbuf, l_recvcount,  l_recvtype,  &
                   l_comm,  l_error)

error   = l_error

Return
End Subroutine MPL_Allgather
