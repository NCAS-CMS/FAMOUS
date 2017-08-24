! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file Copyright
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Subroutine MPL_Alltoallv()
!
!     ******************************************************************
!     * Purpose:
!     *
!     *  The block of data sent from the jth process is received by
!     *  every process and placed in the jth block of the receive buffer
!     *  Sizes/offsets may vary.
!     *
!     *
!     *  Output:  error, recvbuf
!     *
!     ******************************************************************

Subroutine MPL_Alltoallv(sendbuf, sendcnts, senddispls, sendtype, &
                         recvbuf, recvcnts, recvdispls, recvtype, &
                         comm,    error)

use mpl, Only :                               &
    GC_Int_Kind,  GC_Log_Kind,  GC_Real_Kind, &
    MPL_Int_Kind, MPL_Log_Kind, MPL_Real_Kind

Implicit None

! Arguments and Variables at Model/GCOM precision level
Integer (KIND=GC_INT_KIND) :: &
  sendbuf(*),                 &
  sendcnts(*),                &
  senddispls(*),              &
  sendtype,                   &
  recvbuf(*),                 &
  recvcnts(*),                &
  recvdispls(*),              &
  recvtype,                   &
  comm,                       &
  error                       

! Arguments and Variables at Internal/MPI Library precision level
Integer (KIND=MPL_INT_KIND) ::  &
  l_sendtype,                   &
  l_recvtype,                   &
  l_comm,                       &
  l_error,                      &
  ssize

Integer (KIND=MPL_INT_KIND), Allocatable :: &
  l_sendcnts(:),                            &
  l_senddispls(:),                          &
  l_recvcnts(:),                            &
  l_recvdispls(:)

Integer (KIND=MPL_INT_KIND), Save :: &
  ssize_save  = -1234,               &
  l_comm_save = -1234

!=======================================================================

l_comm     = comm
l_recvtype = recvtype
l_sendtype = sendtype

! Unlike most of MPL, we need here to know the number of processors
! we are dealing with to allocate some buffers for copying.
! Fortran can't tell how big they are due to implicit sizing.
If (l_comm == l_comm_save) Then
  ! Can use saved size
  ssize = ssize_save
Else
  ! Need to know how many processors we are dealing with
  Call MPI_Comm_Size(l_comm, ssize, l_error)
  ssize_save  = ssize
  l_comm_save = l_comm
End If

Allocate( l_sendcnts(ssize) )
Allocate( l_senddispls(ssize) )
Allocate( l_recvcnts(ssize) )
Allocate( l_recvdispls(ssize) )

l_sendcnts(1:ssize)   = sendcnts(1:ssize)
l_senddispls(1:ssize) = senddispls(1:ssize)
l_recvcnts(1:ssize)   = recvcnts(1:ssize)
l_recvdispls(1:ssize) = recvdispls(1:ssize)

Call MPI_Alltoallv(sendbuf, l_sendcnts, l_senddispls, l_sendtype,  &
                   recvbuf, l_recvcnts, l_recvdispls, l_recvtype,  &
                   l_comm,  l_error)

error   = l_error
Deallocate( l_sendcnts )
Deallocate( l_senddispls )
Deallocate( l_recvcnts )
Deallocate( l_recvdispls )

Return
End Subroutine MPL_Alltoallv
