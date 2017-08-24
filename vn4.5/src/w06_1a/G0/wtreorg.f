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

      SUBROUTINE WTREORG(PSWH,PERIO,PDIR,KBLO,KJS,KJL,
     %                   KWTMAX,KWTOT,PFWIND,PDWIND,PDMAX,PMISS,
     %                   KFLAGWS,PMCOEF,KREOSP,KANG,KFRE,KWTRA)
C
C**** *WTREORG* - ROUTINE TO FIND WINDSEA AND CLASSIFY WAVE TRAINS
C
C     A.GUILLAUME      ECMWF                01/07/92
C     A.GUILLAUME      ECMWF save memory space 02/94
C
C     PURPOSE.
C     --------
C
C          *WTREORG*
C
C**   INTERFACE.
C     ----------
C
C          *CALL* * WTREORG(PSWH,PERIO,PDIR,KBLO,KJS,KJL,
C                           KWTMAX,KWTOT,PFWIND,PDWIND,PDMAX,PMISS,
C                           KFLAGWS,PMCOEF,KREOSP,KANG,KFRE,KWTRA)
C
C       I/O     *PSWH*    - PSWH OF WAVE TRAINS.
C       I/O     *PERIO*   - MEAN PERIOD  OF WAVE TRAINS.
C       I/O     *PDIR*    - MEAN DIRECTION OF WAVE TRAINS.
C       I/      *KBLO*    - DIMENSION OF ONE BLOCK.
C       I/      *KJS*     - INDEX OF FIRST POINT OF BLOCK TO USE.
C       I/      *KJL*     - INDEX OF LAST POINT OF BLOCK TO USE.
C       I/      *KWTMAX*  - MAX NUMBER OF WAVE TRAINS.
C       I/O     *KWTOT*   - FINAL NB OF WAVE TRAINS
C       I/      *PFWIND*  - WIND SPEED TO COMPUTE P.M.FREQUENCY
C       I/      *PDWIND*  - WIND DIRECTION (RADIAN)
C    I/      *PDMAX*   - MAX ANGULAR DISTANCE BETWEEN WIND AND WINDSEA.
C                   NO WINDSEA FURTHER THAN ZTEMAX FROM WIND DIRECTION
C       I/      *PMISS*   - MISSING VALUE
C  I/  *KFLAGWS* - FLAG VALUE TO ISOLATE WINDSEA (DONE IF KFLAGWS.EQ.1,
C                    MUST BE SET TO 0 OTHERWISE,TO SAVE MEMORY SPACE)
C  I/      *PMCOEF*  - TUNING FACTOR FOR FINDING WINDSEA (0.9, 0.8..)
C  I/      *KREOSP*  - FLAG VALUE TO REORGANIZE WAVE TRAIN INDEX MATRIX
C                           DONE IF KREOSP=1
C
C       I/      *KANG*    - NUMBER OF DIRECTIONS.
C       I/      *KFRE*    - NUMBER OF FREQUENCIES.
C  I/O     *KWTRA*   - WAVE TRAIN INDEX MATRIX (ONLY USED IF KREOSP=1)
C
C     METHOD.
C     -------
C
C
C     EXTERNALS.
C     ----------
C
C          NONE.
C
C     REFERENCE.
C     ----------
C
C          NONE.
C
      DIMENSION PSWH(KBLO,KWTMAX),PERIO(KBLO,KWTMAX),
     %          PDIR(KBLO,KWTMAX),PFWIND(1),PDWIND(1),
     %          KWTOT(KBLO),KWTRA(KJL-KJS+1,KANG,KFRE)
C..WORKING ARRAYS:
      DIMENSION ZTET1(KBLO),ZTET2(KBLO),ZTEMAX(KBLO)
      DIMENSION ZPEMIN(KBLO*KFLAGWS),IWDSEA(KBLO*KFLAGWS)
C
C*    *PARAMETER* OF GLOBAL CONSTANTS.
C
CCC   PARAMETER (G = 9.806, PI = 3.14159265358978, CIRC = 40000000.,
CCc  1           ZPI = 2.*PI, RAD = PI/180., DEG = 180./PI,
CCC  2           R = CIRC/ZPI)

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
C..FUNCTION IN LINE
C
      IDELTA(I,J)=(ISIGN(1,I-J)+ISIGN(1,J-I))/2

         ZPI=2.*PI
         RAD=PI_OVER_180
         DEG=RECIP_PI_OVER_180


      IWTR1=1
CAG   WRITE(6,*)' KFLAGWS ',KFLAGWS
CCMH UKMO calls with kflagws=0 so comment out lines to help memory probs
      IF(KFLAGWS.NE.1) GO TO 500
      IWTR1=2
C
C          1. COMPUTE P.M.FREQUENCY.
C             ----------------------
C
100   CONTINUE
ccc   DO 101 J=KJS,KJL
ccc   ZPEMIN(J)=1./0.13/G*PFWIND(J)/PMCOEF
ccc   ZTEMAX(J)=PDMAX
101   CONTINUE
C
C          2. SHIFT WAVE TRAINS TO PUT WIND SEA IN 1ST.
C             ----------------------------------------
C
200   CONTINUE
ccc   DO 201 JWTR=KWTMAX-1,1,-1
cc    DO 201 J=KJS,KJL
cc    PSWH(J,JWTR+1)=PSWH(J,JWTR)
cc    PDIR(J,JWTR+1)=PDIR(J,JWTR)
cc    PERIO(J,JWTR+1)=PERIO(J,JWTR)
201   CONTINUE
CCMH  IF(KREOSP.EQ.1) THEN
cc       DO 202 JFRE=1,KFRE
cc       DO 202 JANG=1,KANG
cc       DO 202 JWTR=KWTMAX,1,-1
cc       DO 202 J=1,KJL-KJS+1
cc       KWTRA(J,JANG,JFRE)=(JWTR+1)*IDELTA(KWTRA(J,JANG,JFRE),JWTR)
cc   %      +KWTRA(J,JANG,JFRE)*(1-IDELTA(KWTRA(J,JANG,JFRE),JWTR))
ccc202      CONTINUE
CAG   WRITE(6,*)'IN WTREORG, AFTER WINDSEA SHIFT'
CAG   PRINT 1000,PERIO,PSWH,PDIR
CAG   PRINT 2000,KWTRA
1000  FORMAT('PERIO =',5F9.2,/,'PSWH  =',5F9.2,/,'PDIR  =',5F9.2,/)
2000  FORMAT(/,(24I2))
cc    ENDIF
cc    DO 203 J=KJS,KJL
cc    PSWH(J,1)=0.
cc    PDIR(J,1)=PMISS
cc    PERIO(J,1)=PMISS
cc    KWTOT(J)=KWTOT(J)+1
203   CONTINUE
C
C       3. FIND WINDSEA WT.
C         THE CLOSER TO WIND DIRECTION WITH PERIOD LESS THAN PM PERIOD
C         ------------------------------------------------------------
C
300   CONTINUE
cc    DO 301 J=KJS,KJL
cc    IWDSEA(J)=1
301   CONTINUE
cc    DO 302 JWTR=2,KWTMAX
cc    DO 303 J=KJS,KJL
cc    ZTET1(J)=MOD(PDIR(J,JWTR)-PDWIND(J)+ZPI,ZPI)
cc    ZTET2(J)=MOD(PDWIND(J)-PDIR(J,JWTR)+ZPI,ZPI)
cc    ZTET1(J)=AMIN1(ABS(ZTET1(J)),ABS(ZTET2(J)))
303   CONTINUE
cc    DO 304 J=KJS,KJL
cc    IF((PERIO(J,JWTR).LT.ZPEMIN(J))
cc   %          .AND.
cc   %    (PERIO(J,JWTR).GT.0.)
cc   %          .AND.
cc   %    (ZTET1(J).LT.ZTEMAX(J))) THEN
cc       IWDSEA(J)=JWTR
cc       ZTEMAX(J)=ZTET1(J)
cc    ENDIF
304   CONTINUE
302   CONTINUE
C
C          4. PUT FIND WINDSEA WT IN 1ST.
C             ---------------------------
C
400   CONTINUE
cc    DO 410 J=KJS,KJL
cc    PSWH(J,1)=PSWH(J,IWDSEA(J))
cc    PDIR(J,1)=PDIR(J,IWDSEA(J))
cc    PERIO(J,1)=PERIO(J,IWDSEA(J))
410   CONTINUE
cc    DO 420 J=KJS,KJL
cc    PSWH(J,IWDSEA(J))=0.
cc    PDIR(J,IWDSEA(J))=PMISS
cc    PERIO(J,IWDSEA(J))=PMISS
420   CONTINUE
cc    IF(KREOSP.EQ.1) THEN
cc       DO 430 JFRE=1,KFRE
cc       DO 430 JANG=1,KANG
cc       DO 430 J=1,KJL-KJS+1
cc       KWTRA(J,JANG,JFRE)=IDELTA(KWTRA(J,JANG,JFRE),IWDSEA(J))
cc   %     +KWTRA(J,JANG,JFRE)*(1-IDELTA(KWTRA(J,JANG,JFRE),IWDSEA(J)))
430      CONTINUE
CAG      WRITE(6,*)'IN WTREORG, AFTER WINDSEA IN 1'
CAG      PRINT 2000,KWTRA
CAG      WRITE(6,*)'IN WTREORG, AFTER SHIFT '
CAG      PRINT 1000,PERIO,PSWH,PDIR
CAG      PRINT 2000,KWTRA
cc    ENDIF
C
C          5. CLASSEMENT PAR ENERGIE DES AUTRES WT.
C             -------------------------------------
C
500   CONTINUE
c
      DO 501 ITR1=IWTR1,KWTMAX
      DO 501 ITR2=ITR1+1,KWTMAX
      DO 502 J=KJS,KJL
C     IF(PSWH(J,ITR1).LT.PSWH(J,ITR2)) ZTEMAX(J)=0
C     ELSE                             ZTEMAX(J)=1
      ZTEMAX(J)=AMAX1(0.,SIGN(1.,PSWH(J,ITR1)-PSWH(J,ITR2)))
502   CONTINUE
      DO 503 J=KJS,KJL
C     IF(PSWH(J,ITR1).LT.PSWH(ITR2)) ZTET1(J)=PSWH(J,ITR2)
C                                    ZTET2(J)=PSWH(J,ITR1)
C     ELSE                           ZTET1(J)=PSWH(J,ITR1)
C                                    ZTET2(J)=PSWH(J,ITR2)
C     ENDIF
      ZTET1(J)=PSWH(J,ITR1)*ZTEMAX(J)
     %         +PSWH(J,ITR2)*(1.-ZTEMAX(J))
      ZTET2(J)=PSWH(J,ITR2)*ZTEMAX(J)
     %         +PSWH(J,ITR1)*(1.-ZTEMAX(J))
503   CONTINUE
      DO 504 J=KJS,KJL
      PSWH(J,ITR1)=ZTET1(J)
      PSWH(J,ITR2)=ZTET2(J)
504   CONTINUE
      DO 505 J=KJS,KJL
         ZTET1(J)=PDIR(J,ITR1)*ZTEMAX(J)
     %            +PDIR(J,ITR2)*(1.-ZTEMAX(J))
         ZTET2(J)=PDIR(J,ITR2)*ZTEMAX(J)
     %            +PDIR(J,ITR1)*(1.-ZTEMAX(J))
505   CONTINUE
      DO 506 J=KJS,KJL
      PDIR(J,ITR1)=ZTET1(J)
      PDIR(J,ITR2)=ZTET2(J)
506   CONTINUE
      DO 507 J=KJS,KJL
      ZTET1(J)=PERIO(J,ITR1)*ZTEMAX(J)
     %         +PERIO(J,ITR2)*(1.-ZTEMAX(J))
      ZTET2(J)=PERIO(J,ITR2)*ZTEMAX(J)
     %         +PERIO(J,ITR1)*(1.-ZTEMAX(J))
507      CONTINUE
      DO 508 J=KJS,KJL
      PERIO(J,ITR1)=ZTET1(J)
      PERIO(J,ITR2)=ZTET2(J)
508   CONTINUE
      IF(KREOSP.EQ.1) THEN
         DO 509 J=1,KJL-KJS+1
         JJZZ=NINT(ZTEMAX(J+KJS-1))
         ITR11=ITR1*JJZZ+ITR2*(1-JJZZ)
         ITR22=ITR2*JJZZ+ITR1*(1-JJZZ)
         DO 509 JFRE=1,KFRE
         DO 509 JANG=1,KANG
         KWTRA(J,JANG,JFRE)=ITR11*IDELTA(KWTRA(J,JANG,JFRE),ITR1)
     %                     +ITR22*IDELTA(KWTRA(J,JANG,JFRE),ITR2)
     %   +KWTRA(J,JANG,JFRE)*(1-IDELTA(KWTRA(J,JANG,JFRE),ITR1))*
     %                       (1-IDELTA(KWTRA(J,JANG,JFRE),ITR2))
509      CONTINUE
CAG   WRITE(6,*)' AFTER classification, ITR1,ITR2 ',ITR1,ITR2
CAG   PRINT 599,KWTRA
599   FORMAT('KWTRA IN WTREORG, AFTER classification ',/,(24I2))
      ENDIF
501   CONTINUE
      RETURN
      END
