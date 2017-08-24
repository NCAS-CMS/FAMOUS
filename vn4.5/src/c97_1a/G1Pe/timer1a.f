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
      SUBROUTINE TIMER(SUB,I)
CFPP$ NOCONCUR R
!   ..................................................................
!   SUBROUTINE TIMER
!   ----------------
!   A PROGRAM TO RECORD THE FLOW THROUGH AND TIME SPENT IN EACH
!   SUBROUTINE OF A MULTI SUBROUTINED PROGRAM.
!   CALLS TO TIMER MUST BE MADE BEFORE THE FIRST EXECUTABLE STATEMENT
!   OF EACH SUBROUTINE AND BEFORE EACH RETURN POINT.
!
!   PARAMETERS:
!   SUB - 8 BYTE CHARACTER STRING CONTAINING NAME OF CALLING
!         PROGRAM WHICH HAS BEEN LEFT ADJUSTED
!
!   I  -  I=1 - FIRST CALL FROM MAIN (OR THE HIGHEST LEVEL PROGRAM)
!         I=2 - LAST CALL FROM MAIN ( OR THE HIGHEST LEVEL PROGRAM)
!         I=3 - FIRST CALL FROM LOWER SUBROUTINE
!         I=4 CALL BEFORE RETURN STATEMENTS IN LOWER SUBROUTINE
!
!   ---NOTE :THE CRAY FACILITY PERFTRACE IS MORE APPROPRIATE FOR SINGLE
!   TASKED TIMING DIAGNOSTICS, BUT A BREAKDOWN OF ELAPSE TIME BY
!   SUBROUTINE OF MULTITASKED JOBS IS UNAVAILABLE WITHOUT THIS HOMEMADE
!   TIMER ROUTINE
!
!   REVISED TO CONFORM WITH UM SOFTWARE STANDARDS, THIS VERSION RUNS
!   ON BOTH THE CRAY & HP WORKSTATIONS. IF THERE ARE > 200 DIFFERENT
!   SUBROUTINE CALLS TO TIMER, FINAL TABLE REPLACED BY WARNING MESSAGE.
!   ..................................................................
!   ---AUTHOR OF THIS REVISION IAN EDMOND
!   ---DATE 12/04/94
!
!    Model            Modification history from model version 3.0:
!   version  Date
CLL   4.1       12/03/96  Changes to avoid clashes with TIMER2A
CLL                       P.Burton
!LL   4.5       17/04/98  Modified to cope with action (I) > 100
!LL                       intended for MPP timer.       P.Burton
!  4.5  12/06/98  Use CLOCK(...,0,2) for CPUtime on Fujitsu.
!                                        RBarnes@ecmwf.int
!
!   Programming standard :
!
!   Logical components covered : D66
!
!   Project task :
!
!   External documentation:
!
!  ---------------------------------------------------------------------
!

       IMPLICIT NONE

       INTEGER KLIMIT
       PARAMETER(KLIMIT=200)
       CHARACTER*8 SUB,SUBNAME(KLIMIT),RETURNAM(KLIMIT),SWORK
       INTEGER NENTER(KLIMIT),I,L,J,K,NSUBS,ISTOP,NCALLS,IWORK
       INTEGER ICALL
       REAL ELAPSE(KLIMIT),TOTE,ELPEND,ELPSTART,CPUSTART,TOT,AT,CPUEND
       REAL TOTELAP,AVELAP,PCENT,RWORK,TOTLAP,SPEEDUP,P,TIME(KLIMIT)
       REAL UWORK
       SAVE SUBNAME,RETURNAM,NENTER,ELAPSE,TIME
       SAVE K,J,NSUBS,ELPSTART,ISTOP,CPUSTART

       REAL SECOND
!      -----------------------------------------------------------------
       IF (I .GT. 100) THEN
         ICALL=I-100
       ELSE
         ICALL=I
       ENDIF
       IF (ICALL.EQ.1) THEN

!      First call to timer from the main program
!      Set up initial values of variables

         K         = 1
         J         = 0
         ISTOP     = 0
         NSUBS     = 1
         SUBNAME(1)= SUB

         DO L=1,KLIMIT
           ELAPSE(L) = 0.0
           TIME(L) = 0.0
           NENTER(L) = 0

         ENDDO

      NENTER(1)=1
      CALL TIMEF(ELPSTART)
         CPUSTART=SECOND(CPUSTART)

!      -----------------------------------------------------------------
       ELSEIF ((ICALL.EQ.2).AND.(ISTOP.EQ.0)) THEN

!      Last call to timer from main program
!      Print out table of results

!        Stop timer
      CPUEND=SECOND()
      CALL TIMEF(ELPEND)
         ELAPSE(1)=ELAPSE(1)+(ELPEND-ELPSTART)*1.E-3
         TIME(1)=TIME(1)+CPUEND-CPUSTART

!        Calculate total time in program
         TOTE=0.0
         TOT=0.0
         DO K=1,NSUBS

           TOTE=TOTE+ELAPSE(K)
           TOT=TOT+TIME(K)

         ENDDO

!        Sort subroutines into time order

         DO K=1,(NSUBS-1)

           DO J=(K+1),NSUBS

             IF (TIME(J).GT.TIME(K)) THEN

!              Swap the values:
               RWORK=TIME(K)
               TIME(K)=TIME(J)
               TIME(J)=RWORK
               UWORK=ELAPSE(K)
               ELAPSE(K)=ELAPSE(J)
               ELAPSE(J)=UWORK
               IWORK=NENTER(K)
               NENTER(K)=NENTER(J)
               NENTER(J)=IWORK
               SWORK=SUBNAME(K)
               SUBNAME(K)=SUBNAME(J)
               SUBNAME(J)=SWORK

             ENDIF

           ENDDO

         ENDDO


!        Output timing information


      WRITE(*,'(''1'',//,20X,'' FLOW TRACE SUMMARY'',/)')
      WRITE(*,'(4X,''ROUTINE'',6X,''CPU'',6X,''%'',3X,
     &  ''CALLS'',2X,''AVERAGE'',4X,''ELAPSE'',4X,''%''
     &  ,6X,''AVERAGE'',1X,''CPU'')')
      WRITE(*,'(17X,''TIME'',4X,''CPU'',9X,''CPUTIME'',4X,
     &  ''TIME'',3X,''ELAPSE'',4X,''ELAPSE'',2X,''SPEEDUP'')')

         DO K=1,NSUBS

           SWORK=SUBNAME(K)
           TOTLAP=TIME(K)
           TOTELAP=ELAPSE(K)
           NCALLS=NENTER(K)
           AVELAP=TOTELAP/NCALLS
           P=100.0*TOTLAP/TOT
           PCENT=100.0*TOTELAP/TOTE
           AT=TIME(K)/NENTER(K)
           IF (AVELAP.EQ.0.0) THEN
             SPEEDUP=0.0
           ELSE
             SPEEDUP=AT/AVELAP
           ENDIF
           WRITE(*,'(/,T1,I3,T5,A8,T13,F10.4,T25,F5.2,T30,
     &       I5,T35,F10.4,T45,F10.4,T57,F5.2,T62,F10.4,T74,F5.2)')
     &       K,SWORK,TOTLAP,P,NENTER(K),AT,TOTELAP,PCENT,AVELAP, 
     &       SPEEDUP

         ENDDO

         SPEEDUP=TOT/TOTE
         WRITE(*,'(/,T3,''**TOTAL'',T12,F11.4,T44,F11.4,
     &     T74,F5.2)')TOT,TOTE,SPEEDUP

!      -----------------------------------------------------------------
       ELSEIF ((ICALL.EQ.3).AND.(ISTOP.EQ.0)) THEN

!      First call in subroutine

!        Switch off timer
         CPUEND=SECOND()
         CALL TIMEF(ELPEND)
      ELAPSE(K)=ELAPSE(K)+(ELPEND-ELPSTART)*1.E-3
         TIME(K)=TIME(K)+CPUEND-CPUSTART

!        Save name of calling subroutine
      J=J+1
      RETURNAM(J)=SUBNAME(K)

!        Check subroutine name
         DO K=1,NSUBS

           IF (SUBNAME(K).EQ.SUB) GOTO 10

         ENDDO

!        New subroutine entered
      NSUBS=NSUBS+1
         IF (NSUBS .LE. KLIMIT) THEN

      SUBNAME(NSUBS)=SUB
      K=NSUBS

         ELSE

           WRITE(*,'(''WARNING: More than''I4
     &       '' different subroutine calls to TIMER'')')KLIMIT
           ISTOP=1
           GOTO 9999

         ENDIF

 10      CONTINUE

!        Start timer for subroutine
      NENTER(K)=NENTER(K)+1
      CALL TIMEF(ELPSTART)
      CPUSTART=SECOND()
!      -----------------------------------------------------------------
       ELSEIF ((ICALL.EQ.4).AND.(ISTOP.EQ.0)) THEN

!      Return from subroutine

!        Stop timer
      CPUEND=SECOND()
      CALL TIMEF(ELPEND)
      ELAPSE(K)=ELAPSE(K)+(ELPEND-ELPSTART)*1.E-3
      TIME(K)=TIME(K)+CPUEND-CPUSTART

!        Find name of calling program
         DO K=1,NSUBS

           IF (SUBNAME(K).EQ.RETURNAM(J)) GOTO 11

         ENDDO

         WRITE(*,'(3X,''Calling prog:-'',1X,A8,1X,''not found
     &     ,now in'',1X,A8/3X,''TIMER being DISABLED for the rest
     &     of this run'')')RETURNAM(J),SUBNAME(J+1)
      ISTOP=1
         GOTO 9999
 11      CONTINUE

!        Start timer for calling program
      CALL TIMEF(ELPSTART)
      CPUSTART=SECOND()
      J=J-1

!      -----------------------------------------------------------------
       ELSEIF ((ICALL.LT.1).OR.(ICALL.GT.6)) THEN

!      If ICALL<1 or ICALL>6 then there is an error. 
!      If 4<ICALL<=6 then this call
!      to TIMER is ignored. These values of I are recognised by the
!      TIMER3A version.

         WRITE(*,'(3X,
     &     ''Illegal call to TIMER by subroutine'',1X,A8)')SUB

       ENDIF

 9999    CONTINUE

       RETURN
       END


