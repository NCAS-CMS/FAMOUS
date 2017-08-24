! *****************************COPYRIGHT*******************************
! (c) CROWN COPYRIGHT, Met Office, All Rights Reserved.           
! Please refer to Copyright file in top level GCOM directory
!                 for further details
! *****************************COPYRIGHT*******************************

! Fortran header file. PLEASE use the parameter variables in user
! routines calling GC and NOT the numeric values. The latter are
! potentially subject to change without further notice.

#include "gc_kinds.h"

!     GC general options
      INTEGER (KIND=GC_INT_KIND), PARAMETER :: GC_OK         =     0
      INTEGER (KIND=GC_INT_KIND), PARAMETER :: GC_FAIL       =    -1
      INTEGER (KIND=GC_INT_KIND), PARAMETER :: GC_NONE       =     0
      INTEGER (KIND=GC_INT_KIND), PARAMETER :: GC_ANY        =    -1
      INTEGER (KIND=GC_INT_KIND), PARAMETER :: GC_DONTCARE   =    -1

!     GC groups (GCG) support
      INTEGER (KIND=GC_INT_KIND), PARAMETER :: GC_ALLGROUP = 0
      INTEGER (KIND=GC_INT_KIND), PARAMETER :: GCG_ALL = GC_ALLGROUP

!     GC reserved message tags
      INTEGER (KIND=GC_INT_KIND), PARAMETER :: GC_MTAG_LOW   = 999999901
      INTEGER (KIND=GC_INT_KIND), PARAMETER :: GC_MTAG_HIGH  = 999999999

!     GCG_RALLETOALLE index parameters
      INTEGER (KIND=GC_INT_KIND), PARAMETER :: S_DESTINATION_PE = 1
      INTEGER (KIND=GC_INT_KIND), PARAMETER :: S_BASE_ADDRESS_IN_SEND_ARRAY = 2
      INTEGER (KIND=GC_INT_KIND), PARAMETER :: S_NUMBER_OF_ELEMENTS_IN_ITEM = 3
      INTEGER (KIND=GC_INT_KIND), PARAMETER :: S_STRIDE_IN_SEND_ARRAY = 4
      INTEGER (KIND=GC_INT_KIND), PARAMETER :: S_ELEMENT_LENGTH = 5
      INTEGER (KIND=GC_INT_KIND), PARAMETER :: S_BASE_ADDRESS_IN_RECV_ARRAY = 6
      INTEGER (KIND=GC_INT_KIND), PARAMETER :: S_STRIDE_IN_RECV_ARRAY = 7

      INTEGER (KIND=GC_INT_KIND), PARAMETER :: R_SOURCE_PE = 1
      INTEGER (KIND=GC_INT_KIND), PARAMETER :: R_BASE_ADDRESS_IN_RECV_ARRAY = 2
      INTEGER (KIND=GC_INT_KIND), PARAMETER :: R_NUMBER_OF_ELEMENTS_IN_ITEM = 3
      INTEGER (KIND=GC_INT_KIND), PARAMETER :: R_STRIDE_IN_RECV_ARRAY = 4
      INTEGER (KIND=GC_INT_KIND), PARAMETER :: R_ELEMENT_LENGTH = 5
      INTEGER (KIND=GC_INT_KIND), PARAMETER :: R_BASE_ADDRESS_IN_SEND_ARRAY = 6
      INTEGER (KIND=GC_INT_KIND), PARAMETER :: R_STRIDE_IN_SEND_ARRAY = 7

!     Options
      INTEGER (KIND=GC_INT_KIND), PARAMETER :: GC_FORCE_BITREP     = 1
      INTEGER (KIND=GC_INT_KIND), PARAMETER :: GC_ALLTOALL_VERSION = 2

!     Option Status
      INTEGER (KIND=GC_INT_KIND), PARAMETER :: GC_ON  = 1
      INTEGER (KIND=GC_INT_KIND), PARAMETER :: GC_OFF = 0
      INTEGER (KIND=GC_INT_KIND), PARAMETER :: GC_ALLTOALL_ORIG    = 1
      INTEGER (KIND=GC_INT_KIND), PARAMETER :: GC_ALLTOALL_MULTI   = 2
