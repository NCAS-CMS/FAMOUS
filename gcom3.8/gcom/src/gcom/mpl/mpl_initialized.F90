! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file Copyright
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Subroutine MPL_Initialized ()
!
!     ******************************************************************
!     * Purpose:
!     *
!     *  Check if MPI is initialised
!     *
!     *  Output:  error
!     *
!     ******************************************************************

Subroutine MPL_Initialized (flag, error)

use mpl, Only :                               &
    GC_Int_Kind,  GC_Log_Kind,  GC_Real_Kind, &
    MPL_Int_Kind, MPL_Log_Kind, MPL_Real_Kind

Implicit None

! Arguments and Variables at Model/GCOM precision level
Integer (KIND=GC_LOG_KIND) :: flag      ! Return status of init
Integer (KIND=GC_INT_KIND) :: error     ! Error code to return

! Arguments and Variables at Internal/MPI Library precision level
Integer (KIND=MPL_INT_KIND) ::  &
  l_error                       ! Error code returned by MPI_INIT
Integer (KIND=MPL_LOG_KIND) ::  &
  l_flag                        ! Flag returned by MPI_INIT

!=======================================================================

Call MPI_Initialized(l_flag, l_error)

flag  = l_flag
error = l_error

Return
End Subroutine MPL_Initialized
