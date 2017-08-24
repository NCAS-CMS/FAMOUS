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

      SUBROUTINE VAGDIRT(PSPEC,KBLO,KJS,KJL,KANG,KFRE,PFREQ,PFBIN,
     %                   PTHETA,PMISS,PTHETM,df)
C
C**** *VAGDIRT* - ROUTINE TO COMPUTE MEAN WAVE DIRECTION.
C
C     A.GUILLAUME      ECMWF               13/3/92.
C
C
C     PURPOSE.
C     --------
C           *VAGDIRT* CACULATES THE MEAN DIRECTIONS OF WAVE FIELD.
C                     DIRECTIONS ARE GIVEN IN RADIAN.
C
C**   INTERFACE.
C     ----------
C          *CALL* *VAGDIRT(PSPEC,KBLO,KJS,KJL,KANG,KFRE,PFREQ,PFBIN,
C                        PTHETA,PMISS,PTHETM)*
C
C       I/      *PSPEC*   - SPECTRUM.
C       I/      *KBLO*    - DIMENSION OF ONE BLOCK.
C       I/      *KJS*     - INDEX OF FIRST POINT OF BLOCK TO USE.
C       I/      *KJL*     - INDEX OF LAST POINT OF BLOCK TO USE.
C       I/      *KANG*    - NUMBER OF DIRECTIONS.
C       I/      *KFRE*    - NUMBER OF FREQUENCIES.
C       I/      *PFREQ*   - FREQUENCY ARRAY.
C       I/      *PFBIN*   - PFREQ(IF+1)=PFREQ(IF)*(1+PFBIN)
C       I/      *PTHETA   - DIRECTIONS ARRAY.
C    I/      *PMISS    - MISSING VALUE WHEN PTHETM CANNOT BE COMPUTED.
C        /O     *PTHETM*  - MEAN WAVE DIRECTIONS.
C
C     METHOD.
C     -------
C
C     EXTERNALS.
C     ----------
C
C     REFERENCES.
C     -----------
C
      DIMENSION PSPEC(KBLO,KANG,KFRE),PTHETM(KBLO)
      DIMENSION PFREQ(KFRE), PTHETA(KANG),df(kfre)
C WORKING ARRAYS
      DIMENSION ZXX(KBLO),ZYY(KBLO),ZZZ(KBLO)
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


         ZPI=2.*PI
         RAD=PI_OVER_180
         DEG=RECIP_PI_OVER_180

      DO 1 J=KJS,KJL
      ZXX(J)=0.
      ZYY(J)=0.
      PTHETM(J)=PMISS
1     CONTINUE
      DO 2 JANG=1,KANG
      DO 3 J=KJS,KJL
      ZZZ(J)=0.
3     CONTINUE
      DO 4 JFRE=1,KFRE
      DO 4 J=KJS,KJL

CCC      ZZZ(J)=ZZZ(J)+PSPEC(J,JANG,JFRE)
CCC     %             *PFREQ(JFRE)*(1.+1./(1.+PFBIN))
       zzz(j)=zzz(j)+pspec(j,jang,jfre)*df(jfre)
4     CONTINUE
      DO 6 J=KJS,KJL
      ZXX(J)=ZXX(J)+ZZZ(J)*COS(PTHETA(JANG))
      ZYY(J)=ZYY(J)+ZZZ(J)*SIN(PTHETA(JANG))
6     CONTINUE
2     CONTINUE
      DO 7 J=KJS,KJL
      ZZZ(J)=SQRT(AMAX1(ZXX(J)*ZXX(J)+ZYY(J)*ZYY(J),0.))
7     CONTINUE
      DO 8 J=KJS,KJL
      IF(ZZZ(J).EQ.0.) GO TO 8
      PTHETM(J)=ACOS(AMIN1(1.,(AMAX1(-1.,ZXX(J)/ZZZ(J)))))
C      IN COMMENT, NON VECTORIALISED CODE OF THE NEXT TWO LINES.
C      IF(ZYY.LE.0.) PTHETM(IGR)=-PTHETM(IGR)+2*PI
      ZXX(J)=AMAX1(0.,SIGN(1.,ZYY(J)))
      PTHETM(J)=PTHETM(J)*ZXX(J)+(-PTHETM(J)+2*PI)*(1-ZXX(J))
8     CONTINUE
      RETURN
      END
