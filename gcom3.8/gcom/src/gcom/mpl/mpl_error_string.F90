! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Subroutine MPL_Error_String()
!
!     ******************************************************************
!     * Purpose:
!     *
!     *  Returns the error string associated with an error code or class.
!     *
!     *  Input:  errorcode
!     *  Output: error string and its length
!     *          error (whether the call worked)
!     *
!     ******************************************************************

Subroutine MPL_Error_String (errorcode, errorstring, resultlen, error)

use mpl, Only :          &
    GC_Int_Kind,         & 
    MPL_Int_Kind        

Implicit None

! Arguments and Variables at Model/GCOM precision level
Integer (KIND=GC_INT_KIND) :: &
  errorcode,                  &
  resultlen,                  &
  error                       

Character(Len=*) :: errorstring
  

! Arguments and Variables at Internal/MPI Library precision level
Integer (KIND=MPL_INT_KIND) :: &
  l_errorcode,                 &
  l_resultlen,                 &
  l_error                       

!=======================================================================

l_errorcode  = errorcode

Call MPI_Error_String(l_errorcode, errorstring, l_resultlen, l_error)

resultlen    = l_resultlen
error        = l_error

Return
End Subroutine MPL_Error_String
