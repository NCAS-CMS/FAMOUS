! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file Copyright
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Subroutine MPL_Scatterv()
!
!     ******************************************************************
!     * Purpose:
!     *
!     *  Scatter a vector from 1 processor to many
!     *
!     *  Output:  error, recvbuf
!     *
!     ******************************************************************

Subroutine MPL_Scatterv(sendbuf, sendcnts, displs, sendtype,          &
                        recvbuf, recvcnt,  recvtype, root, comm, error)

use mpl, Only :                               &
    GC_Int_Kind,  GC_Log_Kind,  GC_Real_Kind, &
    MPL_Int_Kind, MPL_Log_Kind, MPL_Real_Kind

Implicit None

! Arguments and Variables at Model/GCOM precision level
Integer (KIND=GC_INT_KIND) :: &
  sendbuf(*),                 &
  recvbuf(*),                 &
  sendcnts(*),                &
  displs(*),                  &
  recvcnt,                    &
  sendtype,                   &
  recvtype,                   &
  root,                       &
  comm,                       &
  error                       
  

! Arguments and Variables at Internal/MPI Library precision level
Integer (KIND=MPL_INT_KIND) ::  &
  l_recvcnt,                    & 
  l_sendtype,                   &
  l_recvtype,                   &
  l_root,                       &
  l_comm,                       &
  l_error,                      &
  ssize

Integer (KIND=MPL_INT_KIND), Allocatable :: &
  l_sendcnts(:),                            &
  l_displs(:)

Integer (KIND=MPL_INT_KIND), Save :: &
  ssize_save  = -1234,               &
  l_comm_save = -1234

!=======================================================================

l_root     = root
l_comm     = comm
l_recvtype = recvtype
l_sendtype = sendtype
l_recvcnt  = recvcnt

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
Allocate( l_displs(ssize) )

l_sendcnts(1:ssize) = sendcnts(1:ssize)
l_displs(1:ssize)   = displs(1:ssize)


Call MPI_Scatterv(sendbuf, l_sendcnts, l_displs, l_sendtype,  &
                  recvbuf, l_recvcnt,  l_recvtype,            &
                  l_root,  l_comm,     l_error)

error   = l_error
Deallocate( l_sendcnts )
Deallocate( l_displs )

Return
End Subroutine MPL_Scatterv
