! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file Copyright
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Subroutine MPL_Type_Create_Resized()
!
!     ******************************************************************
!     * Purpose:
!     *
!     *  Define MPI Type resized
!     *
!     *  Input :  oldtype, lb, extent
!     *  Output:  newtype, error
!     *
!     ******************************************************************

Subroutine MPL_Type_Create_Resized (oldtype, lb, extent, &
                                    newtype, error)

USE mpl, ONLY :                               &
    GC_Int_Kind,  GC_Log_Kind,  GC_Real_Kind, &
    MPL_Int_Kind, MPL_Log_Kind, MPL_Real_Kind, &
    MPL_ADDRESS_KIND

Implicit None

! Arguments and Variables at Model/GCOM precision level
Integer (KIND=GC_INT_KIND) :: &
  oldtype,                    &
  lb,                         &
  extent,                     &
  newtype,                    &
  error

! Arguments and Variables at Internal/MPI Library precision level
Integer (KIND=MPL_INT_KIND) :: &
  l_oldtype,                   &
  l_newtype,                   &
  l_error

Integer (KIND=MPL_ADDRESS_KIND) :: &
  l_lb,                            &
  l_extent

!=======================================================================

l_oldtype    = oldtype
l_lb         = lb
l_extent     = extent

Call MPI_Type_Create_Resized(l_oldtype, l_lb, l_extent,             &
                             l_newtype, l_error)

newtype      = l_newtype
error        = l_error

Return
End Subroutine MPL_Type_Create_Resized
