! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file Copyright
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Subroutine MPL_Abort()
!
!     ******************************************************************
!     * Purpose:
!     *
!     *  Error report
!     *
!     *  Output:  error
!     *
!     ******************************************************************

Subroutine MPL_Abort (comm, code, error)

use mpl, Only :                               &
    GC_Int_Kind,  GC_Log_Kind,  GC_Real_Kind, &
    MPL_Int_Kind, MPL_Log_Kind, MPL_Real_Kind

Implicit None

! Arguments and Variables at Model/GCOM precision level
Integer (KIND=GC_INT_KIND) :: &
  comm,                       & 
  code,                       &
  error                       
  

! Arguments and Variables at Internal/MPI Library precision level
Integer (KIND=MPL_INT_KIND) ::  &
  l_comm,                       &
  l_code,                       &
  l_error                        

!=======================================================================

l_comm = comm
l_code = code

Call MPI_Abort(l_comm,l_code,l_error)

error = l_error

Return
End Subroutine MPL_Abort
