C ******************************COPYRIGHT******************************
C (c) CROWN COPYRIGHT 1995, METEOROLOGICAL OFFICE, All Rights Reserved.
C
C Use, duplication or disclosure of this code is subject to the
C restrictions as set forth in the contract.
C
C                Meteorological Office
C                London Road
C                BRACKNELL
C                Berkshire UK
C                RG12 2SZ
C 
C If no contract has been raised with this copy of the code, the use,
C duplication or disclosure of it is strictly prohibited.  Permission
C to do so must first be obtained in writing from the Head of Numerical
C Modelling at the above address.
C ******************************COPYRIGHT******************************
C
! Gets the CPU time
!
! Function interface:
      REAL*8 FUNCTION SECOND()

      IMPLICIT NONE
!
! Description:
!   SECOND calls the HP system function ETIME which returns the number
!   of CPU seconds which have elapsed.
!   Function SECOND having a similarly named inbuilt cray function 'seco
!   acts as an interface to ETIME thus allowing TIMER1A to replace the H
!   version TIMER2A.
!
! Current Code Owner: Ian Edmond
!
! History:
! Version    Date     Comment
! -------    ----     -------
! vn3.4      18/11/94  Original code. Ian Edmond
!   4.5      28/08/98  Remove call to TSECND.
!                      Bob Carruthers. Cray Research.
!
! Code Description:
!   Language: FORTRAN 77 + some CRAY extensions
!   This code is written to UMDP3 v6 programming standards.
!
! System component covered:
! System Task:
!
!- End of header
      REAL*4 ETIME, dummy(2), ELTIME2
      REAL*8 ELTIME
      ELTIME2=ETIME(dummy)
      ELTIME=ELTIME2
      SECOND=ELTIME
       RETURN
       END
! Gets the elapsed time
!
! Subroutine Interface:
      SUBROUTINE TIMEF(elptime)

      IMPLICIT NONE
!
! Description:
!   TIMEF calls the HP system function SECNDS which returns the number
!   of seconds which have elapsed. Subroutine TIMEF having a similarly n
!   inbuilt cray function 'timef' acts as an interface to
!   SECNDS thus allowing TIMER1A to replace the HP version TIMER2A.
!
! Current Code Owner: Ian Edmond
!
! History:
! Version   Date     Comment
! -------   ----     -------
! vn3.4    18/11/94 Original code. Ian Edmond
! vn4.0    23/03/95 Required for 64bit precision, REAL*4 val
! vn4.3    19/03/97 Corrected scaling factor for T3E IRTC
!                   and added DEFS to allow compilation only on
!                   non CRAY PVP.                   P.Burton
!  4.5  12/06/98  Use GETTOD for elapsed time on Fujitsu.
!                                        RBarnes@ecmwf.int
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!
! System component covered:
! System Task:
!
!- End of header
!
      REAL LAST_TIME  ! value of ELPTIME last time we called
      REAL OFFSET     ! offset to add onto elptime to take account
!                     ! of clock reset at midnight
      DATA LAST_TIME,OFFSET/0.0,0.0/
      SAVE LAST_TIME,OFFSET

      REAL ONE_DAY    ! number of seconds in day
      PARAMETER(ONE_DAY=24*60*60)
       REAL ELPTIME
       REAL*4 val
      REAL*8 SECOND
       val=0.0
      ELPTIME=SECOND(val)+ OFFSET
      IF (ELPTIME .LT. LAST_TIME) THEN
        ELPTIME=ELPTIME+ONE_DAY
        OFFSET=OFFSET+ONE_DAY
      ENDIF
      LAST_TIME=ELPTIME
      ELPTIME=ELPTIME*1000.0
       RETURN
       END
