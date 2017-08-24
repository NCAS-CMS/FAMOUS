! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file Copyright
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
Module MPL

!     ******************************************************************
!     * Purpose:
!     *
!     *  A Module to select the Correct kind for Integer, Real and Log
!     *  variables expected by the MPI library.
!     *  Cast the MPI parameters used in GCOM/the model to thier MPL
!     *  equivalents.
!     *  And define MPL_INteger and MPL_Real to match MPI_INTEGER(8)
!     *  or MPI_REAL(8) as appropriate.
!     *
!     ******************************************************************

#if defined (MPI_SRC)
Use MPI
#endif
#include "gc_kinds.h"

! Select the kind to use for variables to match the MPI Library
#if defined(MPILIB_32B)
Integer, Parameter :: MPL_Int_Kind  = gc_integer32
Integer, Parameter :: MPL_Log_Kind  = gc_integer32
Integer, Parameter :: MPL_Real_Kind = gc_real32
#else
Integer, Parameter :: MPL_Int_Kind  = gc_integer64
Integer, Parameter :: MPL_Log_Kind  = gc_integer64
Integer, Parameter :: MPL_Real_Kind = gc_real64
#endif


#if defined (MPI_SRC)
! Cast the parameters used in the MPI library to subtly remnamed
! versions for the MPL interface.
Integer (Kind=GC_INT_KIND), Parameter :: MPL_COMM_WORLD   = MPI_COMM_WORLD
Integer (Kind=GC_INT_KIND), Parameter :: MPL_COMM_SELF    = MPI_COMM_SELF 
Integer (Kind=GC_INT_KIND), Parameter :: MPL_STATUS_SIZE  = MPI_STATUS_SIZE
Integer (Kind=GC_INT_KIND), Parameter :: MPL_TAG_UB       = MPI_TAG_UB     
Integer (Kind=GC_INT_KIND), Parameter :: MPL_ANY_SOURCE   = MPI_ANY_SOURCE 
Integer (Kind=GC_INT_KIND), Parameter :: MPL_ANY_TAG      = MPI_ANY_TAG
Integer (Kind=GC_INT_KIND), Parameter :: MPL_SUM          = MPI_SUM
Integer (Kind=GC_INT_KIND), Parameter :: MPL_MIN          = MPI_MIN
Integer (Kind=GC_INT_KIND), Parameter :: MPL_MAX          = MPI_MAX
Integer (Kind=GC_INT_KIND), Parameter :: MPL_ADDRESS_KIND = MPI_ADDRESS_KIND
Integer (Kind=GC_INT_KIND), Parameter :: MPL_TAG          = MPI_TAG
Integer (Kind=GC_INT_KIND), Parameter :: MPL_SOURCE       = MPI_SOURCE
Integer (Kind=GC_INT_KIND), Parameter :: MPL_ERROR        = MPI_ERROR
Integer (Kind=GC_INT_KIND), Parameter :: MPL_SUCCESS      = MPI_SUCCESS
Integer (Kind=GC_INT_KIND), Parameter :: MPL_INFO_NULL    = MPI_INFO_NULL
Integer (Kind=GC_INT_KIND), Parameter :: MPL_REQUEST_NULL = MPI_REQUEST_NULL
Integer (Kind=GC_INT_KIND), Parameter :: MPL_MAX_ERROR_STRING  = MPI_MAX_ERROR_STRING
Integer (Kind=GC_INT_KIND), Parameter :: MPL_THREAD_MULTIPLE   = MPI_THREAD_MULTIPLE
Integer (Kind=GC_INT_KIND), Parameter :: MPL_THREAD_SINGLE     = MPI_THREAD_SINGLE
Integer (Kind=GC_INT_KIND), Parameter :: MPL_THREAD_FUNNELED   = MPI_THREAD_FUNNELED
Integer (Kind=GC_INT_KIND), Parameter :: MPL_THREAD_SERIALIZED = MPI_THREAD_SERIALIZED


! Set parameters for MPI types, including some set according to
! GC precision (integer, real)
Integer (Kind=GC_INT_KIND), Parameter :: MPL_BYTE       = MPI_BYTE
Integer (Kind=GC_INT_KIND), Parameter :: MPL_CHARACTER  = MPI_CHARACTER
Integer (Kind=GC_INT_KIND), Parameter :: MPL_INTEGER4   = MPI_INTEGER4
Integer (Kind=GC_INT_KIND), Parameter :: MPL_INTEGER8   = MPI_INTEGER8
Integer (Kind=GC_INT_KIND), Parameter :: MPL_REAL4      = MPI_REAL4
Integer (Kind=GC_INT_KIND), Parameter :: MPL_REAL8      = MPI_REAL8
Integer (Kind=GC_INT_KIND), Parameter :: MPL_LOGICAL4   = MPI_INTEGER4
Integer (Kind=GC_INT_KIND), Parameter :: MPL_LOGICAL8   = MPI_INTEGER8

#if defined(PREC_32B)
Integer (Kind=GC_INT_KIND), Parameter :: MPL_LOGICAL    = MPI_LOGICAL
Integer (Kind=GC_INT_KIND), Parameter :: MPL_INTEGER    = MPI_INTEGER
Integer (Kind=GC_INT_KIND), Parameter :: MPL_REAL       = MPI_REAL
#else
! mpich2 does not have MPI_LOGICAL8, so use MPI_INTEGER8
Integer (Kind=GC_INT_KIND), Parameter :: MPL_LOGICAL    = MPI_INTEGER8
Integer (Kind=GC_INT_KIND), Parameter :: MPL_INTEGER    = MPI_INTEGER8
Integer (Kind=GC_INT_KIND), Parameter :: MPL_REAL       = MPI_REAL8
#endif

! Dummy version for NON-MPI GCOM
#else
Integer (Kind=GC_INT_KIND), Parameter :: MPL_COMM_WORLD   = -99
Integer (Kind=GC_INT_KIND), Parameter :: MPL_STATUS_SIZE  = -99
Integer (Kind=GC_INT_KIND), Parameter :: MPL_TAG_UB       = -99
Integer (Kind=GC_INT_KIND), Parameter :: MPL_ANY_SOURCE   = -99
Integer (Kind=GC_INT_KIND), Parameter :: MPL_ANY_TAG      = -99
Integer (Kind=GC_INT_KIND), Parameter :: MPL_SUM          = -99
Integer (Kind=GC_INT_KIND), Parameter :: MPL_MIN          = -99
Integer (Kind=GC_INT_KIND), Parameter :: MPL_MAX          = -99
Integer (Kind=GC_INT_KIND), Parameter :: MPL_TAG          = -99
Integer (Kind=GC_INT_KIND), Parameter :: MPL_SOURCE       = -99
Integer (Kind=GC_INT_KIND), Parameter :: MPL_ERROR        = -99
Integer (Kind=GC_INT_KIND), Parameter :: MPL_SUCCESS      = -99
Integer (Kind=GC_INT_KIND), Parameter :: MPL_INFO_NULL    = -99
Integer (Kind=GC_INT_KIND), Parameter :: MPL_REQUEST_NULL = -99
Integer (Kind=GC_INT_KIND), Parameter :: MPL_THREAD_MULTIPLE   = -99
Integer (Kind=GC_INT_KIND), Parameter :: MPL_THREAD_SINGLE     = -99
Integer (Kind=GC_INT_KIND), Parameter :: MPL_THREAD_FUNNELED   = -99
Integer (Kind=GC_INT_KIND), Parameter :: MPL_THREAD_SERIALIZED = -99

Integer (Kind=GC_INT_KIND), Parameter :: MPL_ADDRESS_KIND = GC_INT_KIND
Integer (Kind=GC_INT_KIND), Parameter :: MPL_BYTE         = -99
Integer (Kind=GC_INT_KIND), Parameter :: MPL_CHARACTER    = -99
Integer (Kind=GC_INT_KIND), Parameter :: MPL_INTEGER      = -99
Integer (Kind=GC_INT_KIND), Parameter :: MPL_REAL         = -99
Integer (Kind=GC_INT_KIND), Parameter :: MPL_LOGICAL      = -99
Integer (Kind=GC_INT_KIND), Parameter :: MPL_REAL4        = -99
Integer (Kind=GC_INT_KIND), Parameter :: MPL_REAL8        = -99
#endif


End Module MPL
