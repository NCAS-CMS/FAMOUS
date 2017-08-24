! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file Copyright
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Subroutine mpl_init_thread()
!
!     ******************************************************************
!     * Purpose: Initialize the (threaded) MPI execution environment 
!     * 
!     * Input:  required
!     * Output: provided, error
!     ******************************************************************


Subroutine Mpl_Init_Thread (required, provided, error)

Use mpl, Only :                               &
    GC_Int_Kind,  GC_Log_Kind,  GC_Real_Kind, &
    MPL_Int_Kind, MPL_Log_Kind, MPL_Real_Kind

Implicit None

! Arguments and Variables at Model/GCOM precision level
Integer (KIND=GC_INT_KIND) ::error   
Integer (KIND=GC_INT_KIND) ::required
Integer (KIND=GC_INT_KIND) ::provided

! Arguments and Variables at Internal/MPI Library precision level
Integer (KIND=MPL_INT_KIND) ::  l_error
Integer (KIND=MPL_INT_KIND) ::  l_provided
Integer (KIND=MPL_INT_KIND) ::  l_required

!=======================================================================

l_required=required

Call MPI_Init_Thread(l_required, l_provided, l_error)

provided=l_provided
error = l_error

Return
End Subroutine Mpl_Init_Thread
