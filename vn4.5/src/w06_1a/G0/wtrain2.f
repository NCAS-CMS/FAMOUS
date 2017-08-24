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

      SUBROUTINE WTRAIN2(PSPEC,KBLO,KJS,KJL,KANG,KFRE,PRES,KDANG,
     %                   KWTRA,KWTMAX)
C
C**** *WTRAIN2* - WAVE TRAIN 2ND ROUTINE (FIND WAVE TRAINS)
C
C     A.GUILLAUME
C     A.GUILLAUME   ECMWF    02/94 modified to save memory space
C
C
C*    PURPOSE.
C     --------
C
C       FIND WAVE TRAINS.
C
C**   INTERFACE.
C     ----------
C
C   *CALL* *WTRAIN2(PSPEC,KBLO,KJS,KJL,KANG,KFRE,PFREQ,KNRJ,PRES,KDANG,
C                       PSPEWT,KWTMAX)*
C
C       I/      *PSPEC*   - SPECTRUM.
C       I/      *KBLO*    - DIMENSION OF ONE BLOCK.
C       I/      *KJS*     - INDEX OF FIRST POINT OF BLOCK TO USE.
C       I/      *KJL*     - INDEX OF LAST POINT OF BLOCK TO USE.
C       I/      *KANG*    - NUMBER OF DIRECTIONS.
C       I/      *KFRE*    - NUMBER OF FREQUENCIES.
C       I/      *PRES*    - INTERPOLLATION PRECISION
C I/      *KDANG*   - MAX NUMBER OF DIRECTIONS IN SPREADING OF THE WT
C                         (TYPICALLY TO ACHIEVE 60DEG, KDANG=2 WHEN
C                            KANG=12)
C        /O     *KWTRA*   - WAVE TRAIN INDEX OF EACH SPECTRUM BIN.
C       I/      *KWTMAX*  - MAX NB OF WAVE TRAINS.
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
C       WTRAIN1
C
C     REFERENCE.
C     ----------
C
C       NONE.
C
      DIMENSION PSPEC(KBLO,KANG,KFRE),KWTRA(KJL-KJS+1,KANG,KFRE)
C WORKING ARRAYS :
      DIMENSION KPICA(KJL-KJS+1),KPICF(KJL-KJS+1)
      DIMENSION KNRJ(KJL-KJS+1,KANG,KFRE),KNR(KJL-KJS+1,KANG,KFRE)
C
C*    0. INITIALIZE KWTRA.
C        ----------------
C
10    CONTINUE
      DO 11 JANG=1,KANG
      DO 11 JFRE=1,KFRE
      DO 11 J=1,KJL-KJS+1
      KWTRA(J,JANG,JFRE)=KWTMAX+1
11    CONTINUE
C
C*    1. RE-SCALE SPECTRA TO INTEGER VALUES.
C        -----------------------------------
C
100   CONTINUE

      CALL WTRAIN1(PSPEC,KBLO,KJS,KJL,KANG,KFRE,PRES,KNRJ)
CAG   WRITE(6,*)' IN WTRAIN2, AFTER WTRAIN1 '
C
C*    2. START LOOP TO FIND WAVE TRAINS.
C        -------------------------------
C
200   CONTINUE
      DO 201 IWT=1,KWTMAX
      CALL FINDPIC(KNRJ,KJL-KJS+1,1,KJL-KJS+1,KANG,KFRE,KPICA,KPICF)


      CALL TRHOU(KNRJ,KJL-KJS+1,1,KJL-KJS+1,KANG,KFRE,KPICA,KPICF,
     %           KDANG,KNR)
CAG   WRITE(6,*)' IN WTRAIN2, AFTER TRHOU, IWT ',IWT
      DO 202 JFRE=1,KFRE
      DO 202 JANG=1,KANG
      DO 202 J=1,KJL-KJS+1
      KWTRA(J,JANG,JFRE)=IWT*KNR(J,JANG,JFRE)
     %         +(1-KNR(J,JANG,JFRE))*KWTRA(J,JANG,JFRE)
202   CONTINUE
      DO 203 JFRE=1,KFRE
      DO 203 JANG=1,KANG
      DO 203 J=1,KJL-KJS+1
      KNRJ(J,JANG,JFRE)=KNRJ(J,JANG,JFRE)
     %                      *(1-KNR(J,JANG,JFRE))
203   CONTINUE
201   CONTINUE
      RETURN
      END
