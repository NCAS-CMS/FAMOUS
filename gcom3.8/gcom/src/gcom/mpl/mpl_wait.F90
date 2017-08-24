! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file Copyright
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Subroutine MPL_Waitall()
!
!     ******************************************************************
!     * Purpose:
!     *
!     *  Wait until request has completed
!     *
!     *  Input:   request
!     *  Output:  error, request, status
!     *
!     ******************************************************************

Subroutine MPL_Wait ( request, status, error)

Use mpl, Only :                                &
    GC_Int_Kind,  GC_Log_Kind,  GC_Real_Kind,  &
    MPL_Int_Kind, MPL_Log_Kind, MPL_Real_Kind, &
    MPL_STATUS_SIZE

Implicit None

! Arguments and Variables at Model/GCOM precision level
Integer (KIND=GC_INT_KIND) ::       &
  request,                          &
  status(MPL_STATUS_SIZE),          &
  error
  
! Arguments and Variables at Internal/MPI Library precision level
Integer (KIND=MPL_INT_KIND) ::       &
  l_request,                         &
  l_status(MPL_STATUS_SIZE),         &
  l_error

!=======================================================================

l_request = request

Call MPI_Wait( l_request, l_status, l_error)

status(:) = l_status(:)
request   = l_request
error     = l_error

Return
End Subroutine MPL_Wait
