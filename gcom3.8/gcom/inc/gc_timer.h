!   Variables required for the GCOM timer
#if defined(GCOM_TIMER)
#define TIMER_NTYPES     6

#define TIMET_BARRIER    1
#define TIMET_SEND       2
#define TIMET_RECV       3
#define TIMET_BCAST_SEND 4
#define TIMET_BCAST_RECV 5
#define TIMET_COLL       6

#define THIS_COUNT       1

REAL (KIND=GC_REAL_KIND)   :: GCOM_TIMER_TIME(TIMER_NTYPES)
                                                         ! timing
INTEGER (KIND=GC_INT_KIND) ::                                     &
  GCOM_TIMER_COUNT(TIMER_NTYPES) ! Number of calls                &
, GCOM_TIMER_DATA(TIMER_NTYPES)  ! Amount of data

COMMON /GCOM_TIMER_COMMON/                                        &
  GCOM_TIMER_TIME                                                 &
, GCOM_TIMER_COUNT                                                &
, GCOM_TIMER_DATA

INTEGER (KIND=GC_INT_KIND) :: start_timestamp(4)  ! Local variable

INTERFACE
  FUNCTION GCOM_TIMESTAMP() RESULT (TIMESTAMP)
    INTEGER (KIND=GC_INT_KIND) :: TIMESTAMP(4)
  END FUNCTION GCOM_TIMESTAMP
END INTERFACE

INTERFACE
  FUNCTION GCOM_TIMEDIFF(TIMESTAMP) RESULT (TIME_SECONDS)
    REAL (KIND=GC_REAL_KIND)                :: TIME_SECONDS
    INTEGER (KIND=GC_INT_KIND), INTENT (IN) :: TIMESTAMP(4)
  END FUNCTION GCOM_TIMEDIFF
END INTERFACE
#endif

