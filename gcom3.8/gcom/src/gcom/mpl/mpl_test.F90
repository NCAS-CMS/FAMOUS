! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file Copyright
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Subroutine MPL_Test()
!
!     ******************************************************************
!     * Purpose:
!     *
!     *  Test requests have completed
!     *
!     *  Input:   request
!     *  Output:  error, flag, status, request
!     *
!     ******************************************************************

Subroutine MPL_Test (request, flag, status, error)

Use mpl, Only :                                &
    GC_Int_Kind,  GC_Log_Kind,  GC_Real_Kind,  &
    MPL_Int_Kind, MPL_Log_Kind, MPL_Real_Kind, &
    MPL_STATUS_SIZE

Implicit None

! Arguments and Variables at Model/GCOM precision level
Integer (KIND=GC_INT_KIND) ::       &
  count,                            &
  request        ,                  &
  status(MPL_STATUS_SIZE),          &
  error                       
  
Logical (KIND=GC_LOG_KIND) :: flag

! Arguments and Variables at Internal/MPI Library precision level
Integer (KIND=MPL_INT_KIND) ::       &
  l_count,                           & 
  l_request,                         &
  l_status(MPL_STATUS_SIZE),         &
  l_error

Logical (KIND=MPL_LOG_KIND) :: l_flag

!=======================================================================

l_request = request

Call MPI_Test(l_request, l_flag, l_status, l_error)

flag          = l_flag
status(:)     = l_status(:)
request       = l_request
error         = l_error

Return
End Subroutine MPL_Test
