!
! Description:
!  This subroutine is part of the wavetrain diagnostic output code
!  developed by Anne Guillaume at MeteoFrance and ECMWF.
!  Introduced into the unified wave moel at UM4.1
!
! Method:
!
!
!
! Current Code Owner: Martin Holt
!
! History:
! Version   Date     Comment
! -------   ----     -------
! UM4.1    June 1996 Code introduced to UM.  M Holt
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!
!- End of header

      SUBROUTINE WTRAIN(PSPEC,KBLO,KJS,KJL,KANG,KFRE,PFWIND,
     %                  PDWIND,
     %                  PFREQ,PFBIN,PTHETA,PRES,KDANG,PDMAX,
     %                  PECUT,PEMINR,PEMAXR,PDTMIN,KWTMAX,
     %                  PMISS,PSWH,PERIO,PDIR,KWTOT
     %                  ,KFLAGWS,PMCOEF,KREOSP,KWTRA,df
     %                  )
C
C**** *WTRAIN* - FIND WAVE TRAINS AND COMPUTE INTEGRATED PARAMETERS.
C
C     A.GUILLAUME      ECMWF                02/07/92
C     A.GUILLAUME      ECMWF  save memory space 2/94
C     M.Holt        UKMO   included array DF  - use in calculation of
C                          spectral integrated parameters (ukmo freqs)
C
C*    PURPOSE.
C     --------
C
C       FIND WAVE TRAINS.
C
C**   INTERFACE.
C     ----------
C
C       *CALL* *WTRAIN(PSPEC,KBLO,KJS,KJL,KANG,KFRE,PFWIND,PDWIND,
C                      PFREQ,PFBIN,PTHETA,PRES,KDANG,PDMAX,PECUT,
C                      PEMINR,PEMAXR,PDTMIN,KWTMAX,PSWH,PERIO,PDIR,
C                      KWTOT,KFLAGWS,PMCOEF,KREOSP,KWTRA,df)
C
C     I/      *PSPEC*   - SPECTRUM.
C     I/      *KBLO*    - DIMENSION OF ONE BLOCK.
C     I/      *KJS*     - INDEX OF FIRST POINT OF BLOCK TO USE.
C     I/      *KJL*     - INDEX OF LAST POINT OF BLOCK TO USE.
C     I/      *KANG*    - NUMBER OF DIRECTIONS.
C     I/      *KFRE*    - NUMBER OF FREQUENCIES.
C     I/      *PFWIND*  - WIND SPEED
C     I/      *PDWIND*  - WIND DIRECTION (IN RADIAN)
C     I/      *PFREQ*   - FREQUENCY MATRIX.
C     I/      *PFBIN*   - PFREQ(IF+1)=PFREQ(IF)*(1+PFBIN)
C     I/      *PTHETA*  - DIRECTION MATRIX (RADIAN)
C     I/ *PRES*    - INTERPOLLATION PRECISION (TYPICALLY PRES=1000.)
C     I/   *KDANG*   - MAX NUMBER OF DIRECTIONS IN SPREADING OF THE WT
C                         (TYPICALLY TO ACHIEVE 60DEG, KDANG=2 WHEN
C                          KANG=12)
C     I/   *PDMAX*   - MAX ANGULAR DISTANCE BETWEEN WIND AND WINDSEA.
C                         (TYPICALLY PI/3)
C     I/   *PECUT*   - WAVE TRAINS WITH ENERGY LESS THAN PECUT*ETOT
C                         ARE DISCARDED AT THE END(TYPICALLY, 0.04)
C     I/   *PEMINR*  - FOR MERGING WAVE TRAINS WITH CLOSE PERIODS
C                         (TYPICALLY 1./(1.+3.*PFBIN) )
C     I/   *PEMAXR*  - FOR MERGING WAVE TRAINS WITH CLOSE PERIODS
C                         (TYPICALLY (1.+3.*PFBIN) )
C     I/   *PDTMIN*  - FOR MERGING WAVE TRAINS WITH CLOSE DIRECTIONS
C                         (TYPICALLY PI/4)
C     I/   *KWTMAX*  - MAX NB OF WAVE TRAINS (TYPICALLY 5)
C     I/    *PMISS*   - MISSING VALUE, SHOULD BE NEGATIVE.
C      /O   *PSWH*    - SWH OF WAVE TRAINS.
C      /O   *PERIO*   - MEAN PERIOD  OF WAVE TRAINS.
C      /O   *PDIR*    - MEAN DIRECTION OF WAVE TRAINS.
C      /O   *KWTOT*   - FINAL NB OF WAVE TRAINS
C     I/    *KFLAGWS* - FLAG VALUE TO ISOLATE WINDSEA
C (DONE IF KFLAGWS.EQ.1,MUST BE SET TO 0 OTHERWISE,TO SAVE MEMORY SPACE)
C     I/   *PMCOEF*  - TUNING FACTOR FOR FINDING WINDSEA (0.9, 0.8..)
C     I/   *KREOSP*  - FLAG VALUE TO REORGANIZE WAVE TRAIN INDEX MATRIX
C                         DONE IF KREOSP=1
C        NOTE there are calls to wtreorg with KREOSP=1
C                             hard wired  in the arg list
C /O *KWTRA*   - WAVE TRAIN INDEX MATRIX (ONLY USEFUL LATER IF KREOSP=1)
C    I/      *df*      - array of frequency intervals (ie as for UKMO)
C
C     METHOD.
C     -------
C
C       NONE.
C
C     EXTERNALS.
C     ----------
C
C       FINDPIC
C       TRHOU
C       VAGDIRT
C       VTOTT
C       WTRAIN1
C       WTRAIN2
C
C     REFERENCE.
C     ----------
C
C       NONE.
C
      DIMENSION PFREQ(KFRE),df(kfre),PTHETA(KANG)
      DIMENSION PSPEC(KBLO,KANG,KFRE)
      DIMENSION KWTRA(KJL-KJS+1,KANG,KFRE)
      DIMENSION PFWIND(1),PDWIND(1),KWTOT(KBLO)
      DIMENSION PSWH(KBLO,KWTMAX),
     %          PDIR(KBLO,KWTMAX),PERIO(KBLO,KWTMAX)
C WORKING ARRAYS :
C               *ZETOF*   - 1)TOTAL ENERGY 2)MIN ENERGY TO DISCARD WT
C               *ZWORK*   -
C
      DIMENSION ZETOF(KBLO)
      DIMENSION ZWORK(KJL-KJS+1,KANG,KFRE)
C
C*    *PARAMETER* OF GLOBAL CONSTANTS.
C
CCC      PARAMETER (G = 9.806, PI = 3.14159265358978, CIRC = 40000000.,
CCC     1           ZPI = 2.*PI, RAD = PI/180., DEG = 180./PI,
CCC     2           R = CIRC/ZPI)

C*L------------------COMDECK C_G----------------------------------------
C G IS MEAN ACCEL DUE TO GRAVITY AT EARTH'S SURFACE
      REAL G

      PARAMETER(G=9.80665)
C*----------------------------------------------------------------------
C*L------------------COMDECK C_PI---------------------------------------
CLL
CLL 4.0 19/09/95  New value for PI. Old value incorrect
CLL               from 12th decimal place. D. Robinson
CLL
      REAL PI,PI_OVER_180,RECIP_PI_OVER_180

      PARAMETER(
     & PI=3.14159265358979323846, ! Pi
     & PI_OVER_180 =PI/180.0,   ! Conversion factor degrees to radians
     & RECIP_PI_OVER_180 = 180.0/PI ! Conversion factor radians to
     &                              ! degrees
     & )
C*----------------------------------------------------------------------


C
C..FUNCION IN LINE
      IDELTA(I,J)=(ISIGN(1,I-J)+ISIGN(1,J-I))/2
      XPI=2.*PI/FLOAT(KANG)
C
C     ---------------------------------------------------------------
C
         ZPI=2.*PI
         RAD=PI_OVER_180
         DEG=RECIP_PI_OVER_180
C*    0. INITIALIZE KWTOT.
C        -----------------
C
      DO 10 J=KJS,KJL
      KWTOT(J)=KWTMAX
10    CONTINUE
C
C     ---------------------------------------------------------------
C
C*    1. COMPUTE TOTAL ENERGY.
C        ---------------------
C
c     array zetof holds the gridpoint energy scaled by pecut. this array
c     is passed into regroup and used as array pemin
c
100   CONTINUE

      CALL VTOTT(PSPEC,KBLO,KJS,KJL,KANG,KFRE,PFREQ,PFBIN,ZETOF,df)
      DO 101 J=KJS,KJL
      ZETOF(J)=PECUT*ZETOF(J)
101   CONTINUE
C
C     ---------------------------------------------------------------
C
C*    2. FIND WAVE TRAINS.
C        ----------------
C
200   CONTINUE
      CALL WTRAIN2(PSPEC,KBLO,KJS,KJL,KANG,KFRE,PRES,KDANG,
     %             KWTRA,KWTMAX-1)
CAG   PRINT*,'AFTER WTRAIN2'
CAG   PRINT 201,KWTRA
201   FORMAT('KWTRA AFTER WTRAIN2',/,(24I2))
C
C     ---------------------------------------------------------------
C
C*    3. COMPUTE INTEGRATED PARAMETERS AND
C        CLASSIFY WT.SPECTRA BY PERIOD.
C        ------------------------------
C
300   CONTINUE
      DO 301 IWT=1,KWTMAX
      DO 311 JFRE=1,KFRE
      DO 311 JANG=1,KANG
      DO 311 J=KJS,KJL
      ZWORK(J-KJS+1,JANG,JFRE)=
     %   PSPEC(J,JANG,JFRE)*IDELTA(KWTRA(J-KJS+1,JANG,JFRE),IWT)
311   CONTINUE
      CALL VINTPAR(ZWORK,PSWH(KJS,IWT),PERIO(KJS,IWT),
     %             PDIR(KJS,IWT),KJL-KJS+1,1,KJL-KJS+1,KANG,KFRE,
     %             PFREQ,PFBIN,PTHETA,PMISS,df)
301   CONTINUE
CAG   PRINT*,'BEFORE WTREORG'
CAG   PRINT 312,PERIO,PSWH,PDIR
312   FORMAT('PERIO =',5F9.2,/,'PSWH  =',5F9.2,/,'PDIR  =',5F9.2,/)
C
CCMH  note the hardwired arguments here
c
      CALL WTREORG(PERIO,PSWH,PDIR,KBLO,KJS,KJL,
     %             KWTMAX-1,KWTOT,PFWIND,PDWIND,0.,PMISS,
     %             0,1.,1,KANG,KFRE,KWTRA)
C
CAG   PRINT*,'AFTER WTREORG PERIOD'
CAG   PRINT 312,PERIO,PSWH,PDIR
CAG   PRINT 321,KWTRA
321   FORMAT('KWTRA AFTER WTREORG',/,(24I2))
C
C     ---------------------------------------------------------------
C
C*    4. REDUCE NB WAVE TRAINS BY MERGING CLOSE ONES.
C        -------------------------------------------
C
400   CONTINUE
ccc   print*,'before calling regroup aray pfreq   kfre'
ccc   print*,kfre,pfreq
      CALL REGROUP(PSPEC,KWTRA,PSWH,PERIO,PDIR,KBLO,KJS,KJL,
     %             KANG,KFRE,KWTMAX-1,KWTOT,
     %             PEMINR,PEMAXR,PDTMIN,PFREQ,PFBIN,PTHETA,
     %             ZETOF,PMISS,df)
CAG   PRINT*,'AFTER REGROUP'
CAG   PRINT 312,PERIO,PSWH,PDIR
CAG   PRINT 401,KWTRA
401   FORMAT('KWTRA AFTER REGROUP',/,(24I2))
C
C     ---------------------------------------------------------------
C
C*    5. COMPUTE SWH PSWH.
C        -----------------
C
500   CONTINUE
ccmh why kwtmax-1 ??
      DO 501 IWT=1,KWTMAX-1
      DO 501 J=KJS,KJL
ccmh  PSWH(J,IWT)=4.004*SQRT(AMAX1(0.,PSWH(J,IWT)*XPI*PFBIN/2.))
ccmh  use this when peto is filled using UKMO df() in vTOTT ?
CCmh  pswh is set with peto in the call to regroup / vintpar / vtott
      PSWH(J,IWT)=4.004*SQRT(MAX(0.,PSWH(J,IWT)*XPI))
501   CONTINUE
c
CAG   PRINT*,'AFTER swh'
CAG   PRINT 312,PERIO,PSWH,PDIR
C
C     ---------------------------------------------------------------
C
C*    6. FIND WIND SEA AND CLASSE LES WT PAR PSWH DECROISSANTS.
C        -----------------------------------------------------
C
600   CONTINUE
      CALL WTREORG(PSWH,PERIO,PDIR,KBLO,KJS,KJL,
     %             KWTMAX,KWTOT,PFWIND,PDWIND,PDMAX,PMISS,
     %                   0,PMCOEF,KREOSP,KANG,KFRE,KWTRA)
C
ccmh  but i have set kflagws to zero so could comment out these calls ?
      CALL WTREORG(PSWH,PERIO,PDIR,KBLO,KJS,KJL,
     %             KWTMAX,KWTOT,PFWIND,PDWIND,PDMAX,PMISS,
     %             KFLAGWS,PMCOEF,KREOSP,KANG,KFRE,KWTRA)
ccc   PRINT*,'AFTER WTREORG WINDSEA'
ccc   PRINT 312,PERIO,PSWH,PDIR
ccc   PRINT 601,KWTRA
601   FORMAT('KWTRA AFTER WTREORG',/,(24I2))
C
C     ---------------------------------------------------------------
C
C*    7. CONVERT DIRECTION IN DEGREE.
C        ---------------------------
C
700   CONTINUE
      DO 701 IWT=1,KWTMAX
      DO 701 J=KJS,KJL
C
c     first convert to degrees in 0 to 360
C     WAM / UM convention is zero=north / incr clockwise
c
      IF (PDIR(J,IWT).NE.PMISS)then
       PDIR(J,IWT)=MOD(PDIR(J,IWT)*180./PI,360.)
ccc    PDIR(J,IWT)=MOD(270. - PDIR(J,IWT),360.)
      endif
701   CONTINUE

      RETURN
      END
