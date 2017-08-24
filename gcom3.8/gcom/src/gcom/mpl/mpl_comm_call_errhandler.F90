! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Subroutine MPL_Comm_Call_Errhandler()
!
!     ******************************************************************
!     * Purpose:
!     *
!     *  Invokes the error handler currently attached to comm
!     *
!     *  Output: error (if error handler returns)
!     *
!     ******************************************************************

Subroutine MPL_Comm_Call_Errhandler (comm, errorcode, error)

use mpl, Only :                               &
    GC_Int_Kind,                              &
    MPL_Int_Kind

Implicit None

! Arguments and Variables at Model/GCOM precision level
Integer (KIND=GC_INT_KIND) :: &
  comm,                       &
  errorcode,                  &
  error                       
  

! Arguments and Variables at Internal/MPI Library precision level
Integer (KIND=MPL_INT_KIND) :: &
  l_comm,                      &
  l_errorcode,                 &
  l_error                       

!=======================================================================

l_comm      = comm
l_errorcode = errorcode

Call MPI_Comm_Call_Errhandler(l_comm, l_errorcode, l_error)

error      = l_error

Return
End Subroutine MPL_Comm_Call_Errhandler
