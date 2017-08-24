! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file Copyright
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Subroutine MPL_Get()
!
!     ******************************************************************
!     * Purpose:
!     *
!     *  Gets data from remote processor
!     *
!     *  Output:  origin_addr, error
!     *
!     ******************************************************************

Subroutine MPL_Get (origin_addr, origin_count, origin_datatype,              &
                    target_rank, target_disp, target_count, target_datatype, &
                    win, error)

Use mpl, Only :                                &
    GC_Int_Kind,  GC_Log_Kind,  GC_Real_Kind,  &
    MPL_Int_Kind, MPL_Log_Kind, MPL_Real_Kind, &
    MPL_ADDRESS_KIND

Implicit None

! Arguments and Variables at Model/GCOM precision level
Integer (KIND=GC_INT_KIND) :: &
  origin_addr(*),             &
  origin_count,               &
  origin_datatype,            &
  target_rank,                & 
  target_count,               &
  target_datatype,            &
  win,                        &
  error                       

! Arguments and Variables at Internal/MPI Library precision level
Integer (KIND=MPL_INT_KIND) ::  &
  l_origin_count,               &
  l_origin_datatype,            &
  l_target_rank,                &
  l_target_count,               &
  l_target_datatype,            &
  l_win,                        &
  l_error                        

Integer (KIND=MPL_ADDRESS_KIND) :: &
  target_disp

!=======================================================================

l_origin_count    = origin_count
l_origin_datatype = origin_datatype
l_target_rank     = target_rank
l_target_count    = target_count
l_target_datatype = target_datatype
l_win             = win

Call MPI_Get(   origin_addr, l_origin_count, l_origin_datatype,              &
              l_target_rank, target_disp, l_target_count, l_target_datatype, &
              l_win, l_error)

error = l_error

Return
End Subroutine MPL_Get
