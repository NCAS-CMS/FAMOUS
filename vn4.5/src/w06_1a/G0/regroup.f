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

      SUBROUTINE REGROUP(PSPEC,KWTRA,PETO,PERIO,PDIR,KBLO,KJS,KJL,
     %                   KANG,KFRE,KWTMAX,KWTOT,
     %                   PEMINR,PEMAXR,PDTMIN,PFREQ,PFBIN,PTHETA,
     %                   PEMIN,PMISS,df)
C
C**** *REGROUP* - ROUTINE TO REDUCE NUMBER OF WAVE TRAINS
C
C     A.GUILLAUME      ECMWF                26/06/92
C     A.GUILLAUME      reduce memory space  09/02/94
C     M. HOLT          set pmiss for zero ht 20/2/96
C
C     PURPOSE.
C     --------
C
C          *REGROUP* REDUCES THE NB OF WT
C
C**   INTERFACE.
C     ----------
C
C          *CALL* *REGROUP(PSPEC,KWTRA,PETO,PERIO,PDIR,KBLO,KJS,KJL,
C                          KANG,KFRE,KWTMAX,KWTOT,
C                          PEMINR,PEMAXR,PDTMIN,PFREQ,PFBIN,PTHETA,
C                          PEMIN,PMISS)
C
C       I       *PSPEC*   - WAVE SPECTRUM.
C I/O   *KWTRA*   - WAVE TRAIN INDEX ASSOCIATED WITH EACH BIN OF PSPEC.
C       I/O     *PETO*    - TOTAL ENERGY OF WAVE TRAINS.
C       I/O     *PERIO*   - MEAN PERIOD  OF WAVE TRAINS.
C       I/O     *PDIR*    - MEAN DIRECTION OF WAVE TRAINS.
C       I/      *KBLO*    - DIMENSION OF ONE BLOCK.
C       I/      *KJS*     - INDEX OF FIRST POINT OF BLOCK TO USE.
C       I/      *KJL*     - INDEX OF LAST POINT OF BLOCK TO USE.
C       I/      *KANG*    - NUMBER OF DIRECTIONS.
C       I/      *KFRE*    - NUMBER OF FREQUENCIES.
C       I/      *KWTMAX*  - MAX NUMBER OF WAVE TRAINS.
C        /O     *KWTOT*   - NUMBER OF WAVE TRAINS.
C       I/      *PEMINR*  - FOR MERGING WAVE TRAINS WITH CLOSE PERIODS
C       I/      *PEMAXR*  - FOR MERGING WAVE TRAINS WITH CLOSE PERIODS
C     I/      *PDTMIN*  - FOR MERGING WAVE TRAINS WITH CLOSE DIRECTIONS
C       I/      *PFREQ*   - FREQUENCY MATRIX.
C       I/      *PFBIN*   - PFREQ(IF+1)=PFREQ(IF)*(1+PFBIN)
C       I/      *PTHETA*  - DIRECTION MATRIX.
C       I/      *PEMIN*   - MIN ENERGY CUT-OFF
C       I/      *PMISS*   - MISSING VALUE
C
C     METHOD.
C     -------
C
C
C     EXTERNALS.
C     ----------
C          VTOTT
C          VAGDIRT
C
C     REFERENCE.
C     ----------
C
C          NONE.
C
      DIMENSION PSPEC(KBLO,KANG,KFRE),KWTRA(KJL-KJS+1,KANG,KFRE)
      DIMENSION PETO(KBLO,KWTMAX),PERIO(KBLO,KWTMAX),
     %          PDIR(KBLO,KWTMAX),KWTOT(KBLO),PEMIN(KBLO)
      DIMENSION PFREQ(KFRE),PTHETA(KANG),df(kfre)
C..WORKING ARRAYS:
      DIMENSION ZWORK(KJL-KJS+1,KANG,KFRE),ZPER(KBLO)
      DIMENSION ZTET1(KBLO),ZTET2(KBLO)
C
C*    *PARAMETER* OF GLOBAL CONSTANTS.
C
CCC   PARAMETER (G = 9.806, PI = 3.14159265358978, CIRC = 40000000.,
CCC  1           ZPI = 2.*PI, RAD = PI/180., DEG = 180./PI,
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
      IDELTA(I,J)=(ISIGN(1,I-J)+ISIGN(1,J-I))/2

         ZPI=2.*PI
         RAD=PI_OVER_180
         DEG=RECIP_PI_OVER_180

C
C          1. MERGING LOOP.
C             -------------
C
ccc   WRITE(6,*)'in routine     regroup array pfreq   kfre'
ccc   WRITE(6,*)kfre,pfreq
100   CONTINUE
      DO 101 JWTR=KWTMAX,1,-1
      DO 101 JWTR2=JWTR-1,1,-1
      DO 102 J=KJS,KJL
      ZPER(J)=PERIO(J,JWTR2)/PERIO(J,JWTR)
      ZTET1(J)=MOD(PDIR(J,JWTR)-PDIR(J,JWTR2)+ZPI,ZPI)
      ZTET2(J)=MOD(PDIR(J,JWTR2)-PDIR(J,JWTR)+ZPI,ZPI)
      ZTET1(J)=AMIN1(ABS(ZTET1(J)),ABS(ZTET2(J)))
102   CONTINUE
C
C          1.1 REGROUPEMENT PAR PERIODE ET PAR DIRECTION.
C              ------------------------------------------
C
110   CONTINUE
      DO 111 J=KJS,KJL
CAG   WRITE(6,*)'J IN REGROUP ',J,' JWTR,JWTR2 ',JWTR,JWTR2
CAG   WRITE(6,*)' ZTET1(J),ZTET2(J),ZPER(J) ',
CAG  %       ZTET1(J),ZTET2(J),ZPER(J)
CAG   PRINT 115,KWTRA
115   FORMAT('KWTRA IN REGROUP',(24I2))
      IF((ZPER(J).LT.PEMAXR).AND.(ZPER(J).GT.PEMINR)) THEN
         IF(ZTET1(J).LT.PDTMIN) THEN
            DO 112 JFRE=1,KFRE
            DO 112 JANG=1,KANG
               KWTRA(J-KJS+1,JANG,JFRE)=
     %            JWTR2*IDELTA(KWTRA(J-KJS+1,JANG,JFRE),JWTR)
     %           +KWTRA(J-KJS+1,JANG,JFRE)
     %            *(1-IDELTA(KWTRA(J-KJS+1,JANG,JFRE),JWTR))
112         CONTINUE
            PETO(J,JWTR)=0.
            PERIO(J,JWTR)=PMISS
            PDIR(J,JWTR)=PMISS
         ENDIF
      ENDIF
111   CONTINUE
C
C          1.2. COMPUTE INTEGRATED PARAMETERS OF NEW WAVE TRAIN.
C               -----------------------------------------------
C
120   CONTINUE
      DO 123 JFRE=1,KFRE
      DO 123 JANG=1,KANG
      DO 123 J=1,KJL-KJS+1
      ZWORK(J,JANG,JFRE)=
     %PSPEC(J+KJS-1,JANG,JFRE)*IDELTA(KWTRA(J,JANG,JFRE),JWTR2)
123   CONTINUE
      CALL VINTPAR(ZWORK,PETO(KJS,JWTR2),PERIO(KJS,JWTR2),
     %             PDIR(KJS,JWTR2),KJL-KJS+1,1,KJL-KJS+1,KANG,KFRE,
     %             PFREQ,PFBIN,PTHETA,PMISS,df)
101   CONTINUE
CAG   WRITE(6,*)' NB TRAINS AVANT REGROUPEMENT ',KWTOT
C
C          2. COUNT AND REORGANIZE WAVE TRAINS.
C             ---------------------------------
C
200   CONTINUE
      DO 201 J=KJS,KJL
      KWTOT(J)=0
201   CONTINUE
      DO 202 J=KJS,KJL
      DO 203 JWTR=1,KWTMAX
CCMH   * note that pemin is proportion of total energy so
CCMH   * if ice point we have
CCMH   * total energy is 2pi*emin or zero ? so the pemin test is passed
      IF((PETO(J,JWTR).NE.PMISS)
     %   .AND.(PETO(J,JWTR).GT.1.e-2)
     %   .AND.(PETO(J,JWTR).GT.PEMIN(J))) THEN
         KWTOT(J)=KWTOT(J)+1
         PETO(J,KWTOT(J))=PETO(J,JWTR)
         PDIR(J,KWTOT(J))=PDIR(J,JWTR)
         PERIO(J,KWTOT(J))=PERIO(J,JWTR)
         DO 204 JFRE=1,KFRE
         DO 204 JANG=1,KANG
         KWTRA(J-KJS+1,JANG,JFRE)=
     %   KWTOT(J)*IDELTA(KWTRA(J-KJS+1,JANG,JFRE),JWTR)
     %   +KWTRA(J-KJS+1,JANG,JFRE)
     %   *(1-IDELTA(KWTRA(J-KJS+1,JANG,JFRE),JWTR))
204      CONTINUE
      ELSE
         PETO(J,JWTR)=0.
         PDIR(J,JWTR)=PMISS
         PERIO(J,JWTR)=PMISS
         DO 205 JFRE=1,KFRE
         DO 205 JANG=1,KANG
         KWTRA(J-KJS+1,JANG,JFRE)=
     %   (KWTMAX+1)*IDELTA(KWTRA(J-KJS+1,JANG,JFRE),JWTR)
     %   +KWTRA(J-KJS+1,JANG,JFRE)
     %   *(1-IDELTA(KWTRA(J-KJS+1,JANG,JFRE),JWTR))
205      CONTINUE
      ENDIF
203   CONTINUE
      DO 206 JWTR=KWTOT(J)+1,KWTMAX
      PETO(J,JWTR)=0.
      PDIR(J,JWTR)=PMISS
      PERIO(J,JWTR)=PMISS
206   CONTINUE
202   CONTINUE
CAG   WRITE(6,*)' NB TRAINS APRES REGROUPEMENT ',KWTOT
      RETURN
      END
