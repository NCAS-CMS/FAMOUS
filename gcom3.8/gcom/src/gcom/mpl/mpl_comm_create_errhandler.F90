! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Subroutine MPL_Comm_Create_Errhandler()
!
!     ******************************************************************
!     * Purpose:
!     *
!     *  Create an MPI Error Handler
!     *
!     *  Output:  handle, and error code
!     *
!     ******************************************************************

Subroutine MPL_Comm_Create_Errhandler (func, errhandler, error)

use mpl, Only :                               &
    GC_Int_Kind,                              &
    MPL_Int_Kind

Implicit None

External :: func

! Arguments and Variables at Model/GCOM precision level
Integer (KIND=GC_INT_KIND) :: &
  errhandler,                 &
  error                       
  

! Arguments and Variables at Internal/MPI Library precision level
Integer (KIND=MPL_INT_KIND) :: &
  l_errhandler,                &
  l_error                       

!=======================================================================

Call MPI_Comm_Create_Errhandler(func, l_errhandler, l_error)

errhandler = l_errhandler
error      = l_error

Return
End Subroutine MPL_Comm_Create_Errhandler
