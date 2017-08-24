! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file Copyright
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Subroutine MPL_Group_Incl()
!
!     ******************************************************************
!     * Purpose:
!     *  Produces a group by reordering an existing group and taking
!     *  only listed members
!     *
!     *  Input:   group, n, ranks
!     *  Output:  newgroup, error
!     *
!     ***************************************

Subroutine MPL_Group_Incl (group, n, irank, newgroup, error)

Use mpl, Only :                               &
    GC_Int_Kind,  GC_Log_Kind,  GC_Real_Kind, &
    MPL_Int_Kind, MPL_Log_Kind, MPL_Real_Kind

Implicit None

! Arguments and Variables at Model/GCOM precision level
Integer (KIND=GC_INT_KIND) :: &
  group,                      &
  n,                          &
  irank(n),                   &
  newgroup,                   &
  error                       
  

! Arguments and Variables at Internal/MPI Library precision level
Integer (KIND=MPL_INT_KIND) ::  &
  l_group,                      &
  l_n,                          &
  l_irank(n),                   &
  l_newgroup,                   &
  l_error                        

!=======================================================================

l_group = group
l_n     = n
l_irank(:) = irank(:)

Call MPI_Group_Incl(l_group, l_n, l_irank, l_newgroup, l_error)

newgroup = l_newgroup
error    = l_error

Return
End Subroutine MPL_Group_Incl
