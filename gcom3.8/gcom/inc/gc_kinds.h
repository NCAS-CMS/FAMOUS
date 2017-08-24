! Description:
!   Contains parameters defining kinds for 32 and 64 integers
!   and reals.

! Parameters for 32 and 64 bit kinds

! Precision and range for 64 bit real
Integer, Parameter :: gc_prec64  = 15
Integer, Parameter :: gc_range64 = 307

! Precision and range for 32 bit real
Integer, Parameter :: gc_prec32  = 6
Integer, Parameter :: gc_range32 = 37

! Range for integers
Integer, Parameter :: gc_irange64=15
Integer, Parameter :: gc_irange32=9

! Kind for 64 bit real
Integer, Parameter :: gc_real64  =                                &
                      selected_real_kind(gc_prec64,gc_range64)
! Kind for 32 bit real
Integer, Parameter :: gc_real32  =                                &
                      selected_real_kind(gc_prec32,gc_range32)
! Kind for 64 bit integer
Integer, Parameter :: gc_integer64 =                              &
                      selected_int_kind(gc_irange64)
! Kind for 32 bit integer
Integer, Parameter :: gc_integer32 =                              &
                      selected_int_kind(gc_irange32)

#if defined(PREC_32B)
Integer, Parameter :: gc_int_kind  = gc_integer32
Integer, Parameter :: gc_log_kind  = gc_integer32
Integer, Parameter :: gc_real_kind = gc_real32
#else
Integer, Parameter :: gc_int_kind  = gc_integer64
Integer, Parameter :: gc_log_kind  = gc_integer64
Integer, Parameter :: gc_real_kind = gc_real64
#endif


! End GC_KINDS
