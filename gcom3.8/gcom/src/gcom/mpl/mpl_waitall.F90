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
!     *  Wait until all requests have completed
!     *
!     *  Output:  error, flag, statuses
!     *
!     ******************************************************************

Subroutine MPL_Waitall (count, requests, statuses, error)

use mpl, Only :                                &
    GC_Int_Kind,  GC_Log_Kind,  GC_Real_Kind,  &
    MPL_Int_Kind, MPL_Log_Kind, MPL_Real_Kind, &
    MPL_STATUS_SIZE

Implicit None

! Arguments and Variables at Model/GCOM precision level
Integer (KIND=GC_INT_KIND) ::       &
  count,                            &
  requests(count),                  &
  statuses(MPL_STATUS_SIZE,count),  &
  error                       
  
! Arguments and Variables at Internal/MPI Library precision level
Integer (KIND=MPL_INT_KIND) ::       &
  l_count,                           & 
  l_requests(count),                 &
  l_statuses(MPL_STATUS_SIZE,count), &
  l_error

!=======================================================================

l_count       = count
l_requests(:) = requests(:)

Call MPI_Waitall(l_count, l_requests, l_statuses, l_error)

statuses(:,:) = l_statuses(:,:)
requests(:)   = l_requests(:)
error         = l_error

Return
End Subroutine MPL_Waitall
