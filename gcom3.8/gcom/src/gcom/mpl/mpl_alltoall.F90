! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file Copyright
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Subroutine MPL_Alltoall()
!
!     ******************************************************************
!     * Purpose:
!     *
!     *  Communication from all to all processors
!     *
!     *  Output:  error
!     *
!     ******************************************************************

Subroutine MPL_Alltoall (sendbuf, sendcount, sendtype, &
                         recvbuf, recvcount, recvtype, &
                         comm, error)

use mpl, Only :                               &
    GC_Int_Kind,  GC_Log_Kind,  GC_Real_Kind, &
    MPL_Int_Kind, MPL_Log_Kind, MPL_Real_Kind

Implicit None

! Arguments and Variables at Model/GCOM precision level
Integer (KIND=GC_INT_KIND) :: &
  sendbuf(*),                 &
  sendcount,                  &
  sendtype,                   &
  recvbuf(*),                 &
  recvcount,                  &
  recvtype,                   &
  comm,                       &
  error                       
  

! Arguments and Variables at Internal/MPI Library precision level
Integer (KIND=MPL_INT_KIND) ::  &
  l_sendcount,                  & 
  l_sendtype,                   &
  l_recvcount,                  & 
  l_recvtype,                   &
  l_comm,                       &
  l_error

!=======================================================================

l_sendcount    = sendcount
l_sendtype     = sendtype
l_recvcount    = recvcount
l_recvtype     = recvtype
l_comm     = comm

Call MPI_Alltoall(sendbuf, l_sendcount, l_sendtype, &
                  recvbuf, l_recvcount, l_recvtype, &
                  l_comm, l_error)

error = l_error

Return
End Subroutine MPL_Alltoall
