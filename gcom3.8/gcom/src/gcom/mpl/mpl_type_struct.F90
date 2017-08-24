! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file Copyright
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Subroutine MPL_Type_Struct()
!
!     ******************************************************************
!     * Purpose:
!     *
!     *  Define MPI vector type
!     *
!     *  Output:  error
!     *
!     ******************************************************************

Subroutine MPL_Type_Struct (count, array_of_blocklengths,             &
                            array_of_displacements, array_of_types,   &
                            newtype, error)

use mpl, Only :                               &
    GC_Int_Kind,  GC_Log_Kind,  GC_Real_Kind, &
    MPL_Int_Kind, MPL_Log_Kind, MPL_Real_Kind

Implicit None

! Arguments and Variables at Model/GCOM precision level
Integer (KIND=GC_INT_KIND) ::    &
  count,                         &
  array_of_blocklengths(count),  &
  array_of_displacements(count), &
  array_of_types(count),         &
  newtype,                       &
  error                       
  

! Arguments and Variables at Internal/MPI Library precision level
Integer (KIND=MPL_INT_KIND) ::       &
  l_count,                           &
  l_array_of_blocklengths(count),    &
  l_array_of_displacements(count),   &
  l_array_of_types(count),           &
  l_newtype,                         &
  l_error                       

!=======================================================================

l_count                     = count
l_array_of_blocklengths(:)  = array_of_blocklengths(:)
l_array_of_displacements(:) = array_of_displacements(:)
l_array_of_types(:)         = array_of_types(:)

Call MPI_Type_Struct(l_count, l_array_of_blocklengths,             &
                     l_array_of_displacements, l_array_of_types,   &
                     l_newtype, l_error)

newtype  = l_newtype
error    = l_error

Return
End Subroutine MPL_Type_Struct
