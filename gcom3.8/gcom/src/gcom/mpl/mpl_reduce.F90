! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file Copyright
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Subroutine MPL_Reduce()
!
!     ******************************************************************
!     * Purpose:
!     *
!     *  Reductions with answers going to root.
!     *
!     *  Input:   sendbuf, count, datatype, op, root, comm
!     *  Output:  recvbuf, error
!     *
!     ******************************************************************

Subroutine MPL_Reduce (sendbuf, recvbuf, count, datatype, op, root, comm, error)

Use mpl, Only :                               &
    GC_Int_Kind,  GC_Log_Kind,  GC_Real_Kind, &
    MPL_Int_Kind, MPL_Log_Kind, MPL_Real_Kind

Implicit None

! Arguments and Variables at Model/GCOM precision level
Integer (KIND=GC_INT_KIND) :: &
  sendbuf(*),                 &
  recvbuf(*),                 &
  count,                      &
  datatype,                   &
  op,                         &
  root,                       &
  comm,                       &
  error                       
  

! Arguments and Variables at Internal/MPI Library precision level
Integer (KIND=MPL_INT_KIND) ::  &
  l_count,                      & 
  l_op,                         & 
  l_root,                       & 
  l_datatype,                   &
  l_comm,                       &
  l_error

!=======================================================================

l_count    = count
l_op       = op
l_datatype = datatype
l_comm     = comm
l_root     = root

Call MPI_Reduce(sendbuf, recvbuf, l_count, l_datatype, l_op, l_root, &
                l_comm, l_error)

error = l_error

Return
End Subroutine MPL_Reduce
