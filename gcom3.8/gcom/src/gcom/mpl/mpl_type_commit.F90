! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file Copyright
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Subroutine MPL_Type_Commit()
!
!     ******************************************************************
!     * Purpose:
!     *
!     *  Commit Type to MPI for use
!     *
!     *  Output:  error
!     *
!     ******************************************************************

Subroutine MPL_Type_Commit (datatype, error)

use mpl, Only :                               &
    GC_Int_Kind,  GC_Log_Kind,  GC_Real_Kind, &
    MPL_Int_Kind, MPL_Log_Kind, MPL_Real_Kind

Implicit None

! Arguments and Variables at Model/GCOM precision level
Integer (KIND=GC_INT_KIND) :: &
  datatype,                   &
  error                       
  

! Arguments and Variables at Internal/MPI Library precision level
Integer (KIND=MPL_INT_KIND) ::  &
  l_datatype,                   &
  l_error                        

!=======================================================================

l_datatype = datatype

Call MPI_Type_Commit(l_datatype,l_error)

datatype = l_datatype
error    = l_error

Return
End Subroutine MPL_Type_Commit
