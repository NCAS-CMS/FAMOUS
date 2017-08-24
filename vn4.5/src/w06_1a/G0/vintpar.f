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

      SUBROUTINE VINTPAR(PSPEC,PETO,PERIO,PDIR,KBLO,KJS,KJL,KANG,KFRE,
     %                   PFREQ,PFBIN,PTHETA,PMISS,df)
C
C**** *VINTPAR* - ROUTINE TO COMPUTE INTEGRATED PARAMETERS.
C
C     A.GUILLAUME      ECMWF                02/07/92
C
C     PURPOSE.
C     --------
C
C          *VINTPAR* COMPUTES ETOT, MEAN DIR AND MEAN PERIOD.
C
C**   INTERFACE.
C     ----------
C
C          *CALL* *VINTPAR(PSPEC,PETO,PERIO,PDIR,KBLO,KJS,KJL,KANG,KFRE,
C                          PFREQ,PFBIN,PTHETA,PMISS)
C
C       I/      *PSPEC    - WAVE SPECTRUM.
C        /O     *PETO*    - TOTAL ENERGY OF WAVE TRAINS.
C        /O     *PERIO*   - MEAN PERIOD  OF WAVE TRAINS.
C        /O     *PDIR*    - MEAN DIRECTION OF WAVE TRAINS.
C       I/      *KBLO*    - DIMENSION OF ONE BLOCK.
C       I/      *KJS*     - INDEX OF FIRST POINT OF BLOCK TO USE.
C       I/      *KJL*     - INDEX OF LAST POINT OF BLOCK TO USE.
C       I/      *KANG*    - NUMBER OF DIRECTIONS.
C       I/      *KFRE*    - NUMBER OF FREQUENCIES.
C       I/      *PFREQ*   - FREQUENCY MATRIX.
C       I/      *PFBIN*   - PFREQ(IF+1)=PFREQ(IF)*(1+PFBIN)
C       I/      *PTHETA*  - DIRECTION MATRIX.
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
      DIMENSION PSPEC(KBLO,KANG,KFRE)
      DIMENSION PETO(KBLO),PERIO(KBLO),
     %          PDIR(KBLO)
      DIMENSION PFREQ(KFRE),PTHETA(KANG),df(kfre)
C..WORKING ARRAYS:
      DIMENSION ZWORK(KJL-KJS+1,KANG,KFRE),ZETO(KJL-KJS+1)
C
C          1. COMPUTE INTEGRATED PARAMETERS OF NEW WAVE TRAIN.
C             -----------------------------------------------
C
100   CONTINUE
c
c     this call puts sum (E*df) into array PETO
c
      CALL VTOTT(PSPEC,KBLO,KJS,KJL,KANG,KFRE,PFREQ,PFBIN,PETO,df)
      CALL VAGDIRT(PSPEC,KBLO,KJS,KJL,KANG,KFRE,
     %             PFREQ,PFBIN,PTHETA,PMISS,PDIR,df)
      DO 101 JFRE=1,KFRE
      DO 101 JANG=1,KANG
      DO 101 J=KJS,KJL
      ZWORK(J-KJS+1,JANG,JFRE)=PSPEC(J,JANG,JFRE)
     %                   /PFREQ(JFRE)
101   CONTINUE
c
c     this call puts sum (E/f *df) into array ZETO
c
      CALL VTOTT(ZWORK,KJL-KJS+1,1,KJL-KJS+1,KANG,KFRE,PFREQ,PFBIN,
     %           ZETO,df)
      DO 102 J=KJS,KJL
      IF(PETO(J).EQ.0.) THEN
         PERIO(J)=PMISS
      ELSE
cc    nb here if NOT using pfbin then peto holds sum(E*df)
cc    nb but zeto holds sum(E/f * df) calc  using same method in VTOT
         PERIO(J)=ZETO(J-KJS+1)/PETO(J)
      ENDIF
102   CONTINUE
      RETURN
      END
