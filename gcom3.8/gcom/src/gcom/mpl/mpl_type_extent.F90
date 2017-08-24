! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file Copyright
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Subroutine MPL_Type_Extent()
!
!     ******************************************************************
!     * Purpose:
!     *
!     *  Define MPI vector type
!     *
!     *  Output:  error
!     *
!     ******************************************************************

Subroutine MPL_Type_Extent (datatype, extent, error)

use mpl, Only :                               &
    GC_Int_Kind,  GC_Log_Kind,  GC_Real_Kind, &
    MPL_Int_Kind, MPL_Log_Kind, MPL_Real_Kind

Implicit None

! Arguments and Variables at Model/GCOM precision level
Integer (KIND=GC_INT_KIND) :: &
  datatype,                   &
  extent,                     &
  error                       
  

! Arguments and Variables at Internal/MPI Library precision level
Integer (KIND=MPL_INT_KIND) :: &
  l_datatype,                  &
  l_extent,                    &
  l_error                       

!=======================================================================

l_datatype    = datatype

Call MPI_Type_Extent(l_datatype, l_extent, l_error)

extent   = l_extent
error    = l_error

Return
End Subroutine MPL_Type_Extent
