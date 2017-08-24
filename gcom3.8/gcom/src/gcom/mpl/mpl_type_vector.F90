! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file Copyright
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Subroutine MPL_Type_Vector()
!
!     ******************************************************************
!     * Purpose:
!     *
!     *  Define MPI vector type
!     *
!     *  Output:  error
!     *
!     ******************************************************************

Subroutine MPL_Type_Vector (count, blocklength, stride, oldtype, newtype, &
                            error)

use mpl, Only :                               &
    GC_Int_Kind,  GC_Log_Kind,  GC_Real_Kind, &
    MPL_Int_Kind, MPL_Log_Kind, MPL_Real_Kind

Implicit None

! Arguments and Variables at Model/GCOM precision level
Integer (KIND=GC_INT_KIND) :: &
  count,                      &
  blocklength,                &
  stride,                     &
  oldtype,                    &
  newtype,                    &
  error                       
  

! Arguments and Variables at Internal/MPI Library precision level
Integer (KIND=MPL_INT_KIND) :: &
  l_count,                     &
  l_blocklength,               &
  l_stride,                    &
  l_oldtype,                   &
  l_newtype,                   &
  l_error                       

!=======================================================================

l_count       = count
l_blocklength = blocklength
l_stride      = stride
l_oldtype     = oldtype

Call MPI_Type_Vector(l_count, l_blocklength, l_stride, l_oldtype, &
                     l_newtype ,l_error)

newtype  = l_newtype
error    = l_error

Return
End Subroutine MPL_Type_Vector
