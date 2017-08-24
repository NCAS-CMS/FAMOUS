! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Subroutine MPL_Comm_Get_Errhandler()
!
!     ******************************************************************
!     * Purpose:
!     *
!     *  Gets the error handler currently attached to comm
!     *
!     *  Output: error handler, error code
!     *
!     ******************************************************************

Subroutine MPL_Comm_Get_Errhandler (comm, errhandler, error)

use mpl, Only :                               &
    GC_Int_Kind,                              &
    MPL_Int_Kind

Implicit None

! Arguments and Variables at Model/GCOM precision level
Integer (KIND=GC_INT_KIND) :: &
  comm,                       &
  errhandler,                 &
  error                       
  

! Arguments and Variables at Internal/MPI Library precision level
Integer (KIND=MPL_INT_KIND) :: &
  l_comm,                      &
  l_errhandler,                &
  l_error                       

!=======================================================================

l_comm     = comm

Call MPI_Comm_Get_Errhandler(l_comm, l_errhandler, l_error)

errhandler = l_errhandler
error      = l_error

Return
End Subroutine MPL_Comm_Get_Errhandler
