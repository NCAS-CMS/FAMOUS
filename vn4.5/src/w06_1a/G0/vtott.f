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

      SUBROUTINE VTOTT(PSPEC,KBLO,KJS,KJL,KANG,KFRE,PFREQ,
     %                 PFBIN,PETO,df)
C
C**** *VTOTT*   - ROUTINE TO CALCULATE TOTAL ENERGY OF SPECTRUM.
C
C     A.GUILLAUME      ECMWF              13/3/92.
C
C     PURPOSE.
C     --------
C          *VTOTT* CALCULATES THE TOTAL ENERGY OF WAVE SPECTRUM.
C
C**   INTERFACE.
C     ----------
C
C       *CALL* *VTOTT(PSPEC,KBLO,KJS,KJL,KANG,KFRE,PFREQ,PFBIN,PETO,df)
C
C       I/      *PSPEC*   - SPECTRUM.
C       I/      *KBLO*    - DIMENSION OF ONE BLOCK.
C       I/      *KJS*     - INDEX OF FIRST POINT OF BLOCK TO USE.
C       I/      *KJL*     - INDEX OF LAST POINT OF BLOCK TO USE.
C       I/      *KANG*    - NUMBER OF DIRECTIONS.
C       I/      *KFRE*    - NUMBER OF FREQUENCIES.
C       I/      *PFREQ*   - FREQUENCY ARRAY.
C       I/      *PFBIN*   - PFREQ(IF+1)=PFREQ(IF)*(1+PFBIN)
C        /O     *PETO*    - TOTAL ENERGY, EXACTLY TOTAL ENERGY IS:
C                           ETOT*2*3.1416*PFBIN/NTET/2.
C       I/      *df*      - array of frequency intervals
C     METHOD.
C     -------
C
C        PETO, USED TO DETERMINE THE TOTAL ENERGY ( TOTAL ENERGY=
C        ETOT*2*3.1416*PFBIN/NTET/2.), IS CALCULATED FROM THE SPECTRUM.
C        FREQUENCIES ARE ASSUMED TO BE IN GEOMETRIC PROGRESSION OF
C        FACTOR (1+PFBIN).
c
c           frequency interval df is calculated as
c
c           df = 0.5* [ Fi+1 - Fi-1]     which rearranges as
c
c           df = pfbin/2  * Fi * [1  +1/(1+pfbin)]
C
C     EXTERNALS.
C     ----------
C          NONE.
C
C     REFERENCE.
C     ----------
C
      DIMENSION PSPEC(KBLO,KANG,KFRE),PETO(KBLO)
      DIMENSION PFREQ(KFRE),df(kfre)

      DO 1 J=KJS,KJL
      PETO(J)=0.
1     CONTINUE
      DO 2 JFRE=1,KFRE
      DO 2 JANG=1,KANG
      DO 2 J=KJS,KJL

      if(pfbin.gt.0.0001) then
      PETO(J)=PETO(J)+PSPEC(J,JANG,JFRE)*PFREQ(JFRE)
     %        *(1.+1./(1.+PFBIN))
      else
c
CC      if using array DF this fills peto with SUM(e*df)
c
       peto(j)=peto(j)+pspec(j,jang,jfre)*df(jfre)
      endif
2     CONTINUE
      RETURN
      END
