! *****************************COPYRIGHT*******************************
! (c) CROWN COPYRIGHT, Met Office, All Rights Reserved.
! Please refer to Copyright file in top level GCOM directory
!                 for further details
! *****************************COPYRIGHT*******************************

! A number of functions are provided for timing internal
! to GCOM, and an API to allow code to access GCOM timing information

#include "gc_prolog.h"
#if defined(GCOM_TIMER)
FUNCTION GCOM_TIMESTAMP() RESULT(TIMESTAMP)

IMPLICIT NONE

#include "gc_kinds.h"

INTEGER (KIND=GC_INT_KIND)       :: TIMESTAMP(4)

INTEGER (KIND=GC_INT_KIND), SAVE :: LAST_COUNT,ROLLOVER

DATA LAST_COUNT/0/
DATA ROLLOVER/0/

CALL SYSTEM_CLOCK(TIMESTAMP(1),TIMESTAMP(3),TIMESTAMP(4))

IF (TIMESTAMP(1)  <   LAST_COUNT) THEN
  ROLLOVER=ROLLOVER+1
ENDIF
LAST_COUNT=TIMESTAMP(1)

TIMESTAMP(2)=ROLLOVER

RETURN
END FUNCTION GCOM_TIMESTAMP

FUNCTION GCOM_TIMEDIFF(TIMESTAMP_START) RESULT(TIME_SECONDS)

IMPLICIT NONE

#include "gc_kinds.h"

REAL (KIND=GC_REAL_KIND)   :: TIME_SECONDS
INTEGER (KIND=GC_INT_KIND) :: TIMESTAMP_START(4),                 &
                              TIMESTAMP_NOW

REAL (KIND=GC_REAL_KIND)   :: ONE_OVER_RATE

INTERFACE
  FUNCTION GCOM_TIMESTAMP() RESULT (TIMESTAMP)
    INTEGER (KIND=GC_INT_KIND) :: TIMESTAMP(4)
  END FUNCTION GCOM_TIMESTAMP
END INTERFACE

TIMESTAMP_NOW=GCOM_TIMESTAMP()
ONE_OVER_RATE=1.0/TIMESTAMP_NOW(3)


! First work out the time for any roll-over period      
TIME_SECONDS=(TIMESTAMP_NOW(2)-                                   &
              TIMESTAMP_START(2))*                                &
             TIMESTAMP_NOW(4)*ONE_OVER_RATE

! And now the straightforward time from the count
TIME_SECONDS=TIME_SECONDS+                                        &
             (TIMESTAMP_NOW(1)-                                   &
              TIMESTAMP_START(1))*ONE_OVER_RATE

RETURN
END FUNCTION GCOM_TIMEDIFF

#endif      
FUNCTION GC_GET_TIMER_TIME(TIMER_TYPE) RESULT(TIME_SECONDS)

IMPLICIT NONE

#include "gc_kinds.h"

REAL (KIND=GC_REAL_KIND)                :: TIME_SECONDS
INTEGER (KIND=GC_INT_KIND) , INTENT(IN) :: TIMER_TYPE

#if defined(GCOM_TIMER)      
#include "gc_timer.h"

TIME_SECONDS=GCOM_TIMER_TIME(TIMER_TYPE)
#else
TIME_SECONDS=0.0
#endif      
RETURN
END FUNCTION GC_GET_TIMER_TIME

FUNCTION GC_GET_TIMER_COUNT(TIMER_TYPE) RESULT(COUNT)

IMPLICIT NONE

#include "gc_kinds.h"

INTEGER (KIND=GC_INT_KIND)              :: COUNT
INTEGER (KIND=GC_INT_KIND), INTENT(IN)  :: TIMER_TYPE

#if defined(GCOM_TIMER)           
#include "gc_timer.h"

COUNT=GCOM_TIMER_COUNT(TIMER_TYPE)
#else
COUNT=0
#endif      
     
RETURN
END FUNCTION GC_GET_TIMER_COUNT

FUNCTION GC_GET_TIMER_DATA(TIMER_TYPE) RESULT(DATA_SIZE)

IMPLICIT NONE

#include "gc_kinds.h"

INTEGER (KIND=GC_INT_KIND)              :: DATA_SIZE
INTEGER (KIND=GC_INT_KIND), INTENT(IN)  :: TIMER_TYPE

#if defined(GCOM_TIMER)                 
#include "gc_timer.h"

DATA_SIZE=GCOM_TIMER_DATA(TIMER_TYPE)
#else
DATA_SIZE=0
#endif      
RETURN
END FUNCTION GC_GET_TIMER_DATA
