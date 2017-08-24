! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file Copyright
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Subroutine MPL_Comm_Dup()
!
!     ******************************************************************
!     * Purpose:
!     *
!     *  Clones a communicator
!     *
!     * Input:
!     *  comm  - the comunicator
!     *
!     * Output:
!     *  newcomm - new communicator
!     *  error  - an error code
!     *
!     ******************************************************************
Subroutine MPL_Comm_Dup (comm, newcomm, error)

Use mpl, Only :                               &
    GC_Int_Kind,  GC_Log_Kind,  GC_Real_Kind, &
    MPL_Int_Kind, MPL_Log_Kind, MPL_Real_Kind

Implicit None

! Arguments and Variables at Model/GCOM precision level
Integer (KIND=GC_INT_KIND) :: &
  comm,                       & 
  newcomm,                    &
  error                       
  

! Arguments and Variables at Internal/MPI Library precision level
Integer (KIND=MPL_INT_KIND) ::  &
  l_comm,                       &
  l_newcomm,                    &
  l_error                        

!=======================================================================

l_comm  = comm

Call MPI_Comm_Dup(l_comm, l_newcomm, l_error)

newcomm = l_newcomm
error   = l_error

Return
End Subroutine MPL_Comm_Dup
