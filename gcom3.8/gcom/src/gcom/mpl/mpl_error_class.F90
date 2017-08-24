! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Subroutine MPL_Error_Class()
!
!     ******************************************************************
!     * Purpose:
!     *
!     *  Map an MPI error code to an error class
!     *
!     *  Input:  errorcode
!     *  Output: errorclass that errorcode maps to
!     *          error (whether the call worked)
!     *
!     ******************************************************************

Subroutine MPL_Error_Class (errorcode, errorclass, error)

use mpl, Only :  &
    GC_Int_Kind, &
    MPL_Int_Kind

Implicit None

! Arguments and Variables at Model/GCOM precision level
Integer (KIND=GC_INT_KIND) :: &
  errorcode,                  &
  errorclass,                 &
  error                       
  

! Arguments and Variables at Internal/MPI Library precision level
Integer (KIND=MPL_INT_KIND) :: &
  l_errorcode,                 &
  l_errorclass,                &
  l_error                       

!=======================================================================

l_errorcode  = errorcode

Call MPI_Error_Class(l_errorcode, l_errorclass, l_error)

errorclass   = l_errorclass
error        = l_error

Return
End Subroutine MPL_Error_Class
