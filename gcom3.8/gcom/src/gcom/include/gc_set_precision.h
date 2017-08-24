#if defined(PREC_32B)
#undef PREC_64B
#define GC__ISIZE 4_GC_INT_KIND
#define GC__RSIZE 4_GC_REAL_KIND
#else

#if !defined(PREC_64B)
#define PREC_64B
#endif

#define GC__ISIZE 8_GC_INT_KIND
#define GC__RSIZE 8_GC_REAL_KIND
#endif
