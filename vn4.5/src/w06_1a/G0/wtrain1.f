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

      SUBROUTINE WTRAIN1(PSPEC,KBLO,KJS,KJL,KANG,KFRE,PRES,ISPEC)
C
C**** *WTRAIN1* - WAVE TRAIN 1ST ROUTINE.
C
C     A.GUILLAUME
C
C
C*    PURPOSE.
C     --------
C
C       COMPUTE MAX OF SPECTRA AND REDUCE TO INTEGER VALUES.
C
C**   INTERFACE.
C     ----------
C
C       *CALL* *WTRAIN1(PSPEC,KBLO,KJS,KJL,KANG,KFRE,PRES,ISPEC)
C
C       I/      *PSPEC*   - SPECTRUM.
C       I/      *KBLO*    - DIMENSION OF ONE BLOCK.
C       I/      *KJS*     - INDEX OF FIRST POINT OF BLOCK TO USE.
C       I/      *KJL*     - INDEX OF LAST POINT OF BLOCK TO USE.
C       I/      *KANG*    - NUMBER OF DIRECTIONS.
C       I/      *KFRE*    - NUMBER OF FREQUENCIES.
C       I/      *PRES*    - INTERPOLLATION PRECISION
C        /O     *ISPEC*   - INTEGER SPECTRA
C
C     METHOD.
C     -------
C
C       NONE.
C
C     EXTERNALS.
C     ----------
C
C       NONE.
C
C     REFERENCE.
C     ----------
C
C       NONE.
C
C
      DIMENSION PSPEC(KBLO,KANG,KFRE)
      DIMENSION ISPEC(KJL-KJS+1,KANG,KFRE)
C..WORKING ARRAYS:
      DIMENSION ZFMAX(KJL-KJS+1)
C
C*    1. INITIALISE ZFMAX.
C        ----------------
C
100   CONTINUE
      DO 101 J=1,KJL-KJS+1
      ZFMAX(J) = 0.1E-20
101   CONTINUE
C
C*    2. FIND MAX OF SPECTRA.
C        -------------------
C
200   CONTINUE
      DO 201 JFRE=1,KFRE
      DO 201 JANG=1,KANG
      DO 201 J=1,KJL-KJS+1
      ZFMAX(J)=AMAX1(ZFMAX(J),PSPEC(J+KJS-1,JANG,JFRE))
201   CONTINUE
C
C*    3. RE_SCALE SPECTRA TO INTEGER VALUES.
C        -----------------------------------
C
300   CONTINUE
      DO 301 JFRE=1,KFRE
      DO 301 JANG=1,KANG
      DO 301 J=1,KJL-KJS+1
      ISPEC(J,JANG,JFRE)=INT(PSPEC(J+KJS-1,JANG,JFRE)/ZFMAX(J)*PRES)
301   CONTINUE
      RETURN
      END
