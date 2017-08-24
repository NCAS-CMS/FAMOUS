! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file Copyright
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Subroutine MPL_Group_Translate_Ranks()
!
!     ******************************************************************
!     * Purpose:
!     *
!     *  Translates ranks of processes in one group to another
!     *
!     *  Output:  error
!     *
!     ******************************************************************

Subroutine MPL_Group_Translate_Ranks (group1, n, ranks1, &
                                      group2, ranks2, error)

use mpl, Only :                               &
    GC_Int_Kind,  GC_Log_Kind,  GC_Real_Kind, &
    MPL_Int_Kind, MPL_Log_Kind, MPL_Real_Kind

Implicit None

! Arguments and Variables at Model/GCOM precision level
Integer (KIND=GC_INT_KIND) :: &
  group1,                     &
  n,                          &
  ranks1(n),                  &
  group2,                     &
  ranks2(n),                  &
  error                       
  

! Arguments and Variables at Internal/MPI Library precision level
Integer (KIND=MPL_INT_KIND) ::  &
  l_group1,                     &
  l_n,                          &
  l_ranks1(n),                  &
  l_group2,                     &
  l_ranks2(n),                  &
  l_error                        

!=======================================================================

l_group1    = group1
l_n         = n
l_group2    = group2
l_ranks1(:) = ranks1(:)

Call MPI_Group_Translate_Ranks(l_group1, l_n, l_ranks1, &
                               l_group2, l_ranks2, l_error)

ranks2(:) = l_ranks2(:)
error     = l_error

Return
End Subroutine MPL_Group_Translate_Ranks
