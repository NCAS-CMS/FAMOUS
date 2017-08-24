! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file Copyright
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Subroutine MPL_Testall()
!
!     ******************************************************************
!     * Purpose:
!     *
!     *  Test requests have completed
!     *
!     *  Output:  error, flag, statuses
!     *
!     ******************************************************************

Subroutine MPL_Testall (count, requests, flag, statuses, error)

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
  
Logical (KIND=GC_LOG_KIND) :: flag

! Arguments and Variables at Internal/MPI Library precision level
Integer (KIND=MPL_INT_KIND) ::       &
  l_count,                           & 
  l_requests(count),                 &
  l_statuses(MPL_STATUS_SIZE,count), &
  l_error

Logical (KIND=MPL_LOG_KIND) :: l_flag

!=======================================================================

l_count       = count
l_requests(:) = requests(:)

Call MPI_Testall(l_count, l_requests, l_flag, l_statuses, l_error)

flag          = l_flag
statuses(:,:) = l_statuses(:,:)
requests(:)   = l_requests(:)
error         = l_error

Return
End Subroutine MPL_Testall
