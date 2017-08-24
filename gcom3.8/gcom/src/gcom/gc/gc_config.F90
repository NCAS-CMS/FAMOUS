! *****************************COPYRIGHT*******************************
! (c) CROWN COPYRIGHT, Met Office, All Rights Reserved.
! Please refer to Copyright file in top level GCOM directory
!                 for further details
! *****************************COPYRIGHT*******************************

#include "gc_prolog.h"

SUBROUTINE GC_CONFIG (MXPROC, MXCOLL, MXPT2PT, INTF)
!     ******************************************************************
!     * Purpose:
!     *
!     *  Return information about the GC configuration.
!     *
!     * Output:
!     *  MXPROC    - maximum numbers of processors compiled into the
!     *              interface - deprecated
!     *  MXCOLL    - maximum number of elements for collective
!     *              operations  - deprecated
!     *  MXPT2PT   - maximum number of elements for point to point
!     *              operations
!     *  INTF      - name of interface selected at compile time
!     *
!     * NOTES:
!     *    
!     ******************************************************************

USE GCOM_MOD, ONLY :                                                   &
    GC_NONE
    
IMPLICIT NONE

#include "gc_kinds.h"

INTEGER (KIND=GC_INT_KIND) :: MXPROC, MXCOLL, MXPT2PT
CHARACTER(LEN=*)           :: INTF


MXCOLL = 0            ! Deprecated
MXPROC = 0            ! Deprecated
MXPT2PT = GC_NONE
INTF = 'GCOM Version ' //                                              &
  GC_VERSION //                                                        &
  ' built at ' //                                                      &
  GC_BUILD_DATE //                                                     &
  ' Interface: ' //                                                    &
#if defined(PREC_32B)
     GC_DESCRIP //                                                     &
  ' 32'
#else
     GC_DESCRIP //                                                     &
  ' 64'
#endif
END
