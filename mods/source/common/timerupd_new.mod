*ID TIMUPD
*/
*/ TIMEF is an SGI intrinsic function and a UM subroutine name so to avoid
*/ confusion remove all calls to TIMEF, replace with call to TIMEFN.
*/ TIMEFN is the same as TIMEF subroutine with machine specific
*/ calls removed and replaced with f90 intrinsic function system_clock.
*/
*DECLARE TIMER1A
*D PXTIMEFN.26,PXTIMEFN.28
*D PXTIMEFN.29,PXTIMEFN.32
      CALL TIMEFN(ELPSTART)
*D PXTIMEFN.33,PXTIMEFN.36
      CALL TIMEFN(ELPEND)
*D PXTIMEFN.37,PXTIMEFN.40
         CALL TIMEFN(ELPEND)
*D PXTIMEFN.41,PXTIMEFN.44
      CALL TIMEFN(ELPSTART)
*D PXTIMEFN.45,PXTIMEFN.48
      CALL TIMEFN(ELPEND)
*D PXTIMEFN.49,PXTIMEFN.52
      CALL TIMEFN(ELPSTART)
*/
*DECLARE TIMER3A
*D PXTIMEFN.53,PXTIMEFN.59
      CALL TIMEFN(temp)
*/
*/ Replace the wallclock timer with vn5.X version
*/
*DECLARE TIMEFN2A
*D PXTIMEFN.14,PXTIMEFN.25
!----------------------------------------------------------------------
!
!+ Gets the elapsed time from the system

! Function Interface:
      SUBROUTINE TIMEFN(Get_Wallclock_Time)

      IMPLICIT NONE
!
! Description:
!   The system function SYSTEM_CLOCK is used to return the numbers of
!   seconds which have elapsed.
!
! Current Code Owner: Anette Van der Wal
!
! IE, PB, Richard Barnes  <- programmer of some or all of previous
!                            code or changes
!
! History:
! Version   Date      Comment
! -------   ----      -------
! 5.3       24/04/01  Complete re-write to simplify.  A Van der Wal
!
! Code Description:
!   Language: FORTRAN 77 plus some Fortran 90 features
!   This code is written to UMDP3 v6 programming standards
!
! System component covered:
! System Task:
!
!- End of header
      REAL Get_Wallclock_Time

! Local variables
      INTEGER, SAVE :: Start_Count=-1
      INTEGER, SAVE :: Old_Count=0
      REAL, SAVE    :: Rollover=0.0
      INTEGER       :: Count, Count_Rate, Count_Max, Elapsed_Count
      REAL, SAVE    :: OneOver_Count_Rate=0.0

! Intrinsic procedures called:
      INTRINSIC SYSTEM_CLOCK

      CALL SYSTEM_CLOCK(Count=Count,Count_Rate=Count_Rate,
     &                  Count_Max=Count_Max)

      IF ((Old_Count .LT. Start_Count) .AND.
     &    ((Count .LT. Old_Count) .OR. (Count .GT. Start_Count))) THEN
        Rollover=Rollover+(REAL(Count_Max)/REAL(Count_Rate))
      ENDIF

      IF (Start_Count .EQ. -1) THEN
        Start_Count=Count
        OneOver_Count_Rate=1.0/REAL(Count_Rate)
      ENDIF

      Elapsed_Count=Count-Start_Count
      IF (Elapsed_Count .LT. 0) Elapsed_Count=Elapsed_Count+Count_Max

      Get_Wallclock_Time = Rollover+
     &                    (REAL(Elapsed_Count)*OneOver_Count_Rate)
      Old_Count=Count

      RETURN
      END
