! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file Copyright
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Subroutine MPL_Comm_Get_Attr ()
!
!     ******************************************************************
!     * Purpose:
!     *
!     *  Get attributes associated with communicator
!     *
!     ******************************************************************

Subroutine MPL_Comm_Get_Attr (comm, keyval, attribute_val, flag, ierror)

use mpl, Only :                                &
    GC_Int_Kind,  GC_Log_Kind,  GC_Real_Kind,  &
    MPL_Int_Kind, MPL_Log_Kind, MPL_Real_Kind, &
    MPL_ADDRESS_KIND

Implicit None

! Arguments and Variables at Model/GCOM precision level
Integer (KIND=GC_INT_KIND) ::      &
  comm,                            &
  keyval,                          &
  ierror

Integer (KIND=MPL_ADDRESS_KIND) :: &
  attribute_val

Logical (KIND=GC_INT_KIND) ::      &
  flag

! Arguments and Variables at Internal/MPI Library precision level
Integer (KIND=MPL_INT_KIND)  ::    &
  l_comm,                          &
  l_keyval,                        &
  l_ierror

Logical (KIND=MPL_INT_KIND) ::     &
  l_flag

!Cast Model precision input arguments to MPI library precision.
l_comm   = comm
l_keyval = keyval
l_flag   = flag

!Call MPI routine
Call MPI_Comm_Get_Attr(l_comm, l_keyval, attribute_val, l_flag, l_ierror)

!Cast returned arguments to Model Preciison
flag   = l_flag
ierror = l_ierror

Return

End Subroutine MPL_Comm_Get_Attr
