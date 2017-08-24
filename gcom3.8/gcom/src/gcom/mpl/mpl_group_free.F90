! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file Copyright
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Subroutine MPL_Group_Free()
!
!     ******************************************************************
!     * Purpose:
!     *
!     *  Frees a group
!     *
!     *  Output:  error
!     *
!     ******************************************************************

Subroutine MPL_Group_Free(group, error)

use mpl, Only :                               &
    GC_Int_Kind,  GC_Log_Kind,  GC_Real_Kind, &
    MPL_Int_Kind, MPL_Log_Kind, MPL_Real_Kind

Implicit None

! Arguments and Variables at Model/GCOM precision level
Integer (KIND=GC_INT_KIND) :: &
  group,                      &
  error                       
  

! Arguments and Variables at Internal/MPI Library precision level
Integer (KIND=MPL_INT_KIND) ::  &
  l_group,                      &
  l_error                        

!=======================================================================

l_group = group

Call MPI_Group_Free(l_group, l_error)

group = l_group
error = l_error

Return
End Subroutine MPL_Group_Free
