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

      SUBROUTINE FINDPIC(KNRJ,KBLO,KJS,KJL,KANG,KFRE,KPICA,KPICF)
C
C**** *FINDPIC* - ROUTINE TO FIND PEAK OF 2D-SPECTRA
C
C     A.GUILLAUME      METEO-FRANCE         20/08/91
C     A.GUILLAUME      ECMWF                30/06/92
C
C     PURPOSE.
C     --------
C
C          *FINDPIC* FIND PEAK DIRECTION AND FREQUENCY
C
C**   INTERFACE.
C     ----------
C
C          *CALL* *FINDPIC(KNRJ,KBLO,KJS,KJL,KANG,KFRE,KPICA,KPICF)
C
C       I/      *KNRJ*    - 2D WAVE SPECTRA (*FACTOR TO BE INTEGER)
C       I/      *KBLO*    - DIMENSION OF ONE BLOCK.
C       I/      *KJS*     - INDEX OF FIRST POINT OF BLOCK TO USE.
C       I/      *KJL*     - INDEX OF LAST POINT OF BLOCK TO USE.
C       I/      *KANG*    - NUMBER OF DIRECTIONS.
C       I/      *KFRE*    - NUMBER OF FREQUENCIES.
C        /O     *KPICA*   - INDEX OF PEAK DIRECTION.
C        /O     *KPICF*   - INDEX OF PEAK FREQUENCY.
C
C     METHOD.
C     -------
C
C
C     EXTERNALS.
C     ----------
C          NONE.
C
C     REFERENCE.
C     ----------
C
C          NONE.
C
      DIMENSION KNRJ(KBLO,KANG,KFRE)
      DIMENSION KPICA(KBLO),KPICF(KBLO)
C..WORKING ARRAYS:
      DIMENSION INRJMX(KBLO)
C
C          1. CHERCHE LE MAX D'ENERGIE.
C             ------- -- --- ----------
C
100   CONTINUE
      DO 101 J=KJS,KJL
      INRJMX(J)=0
101   CONTINUE
      DO 102 JFRE=1,KFRE
      DO 102 JANG=1,KANG
      DO 102 J=KJS,KJL
      INRJMX(J)=MAX0(KNRJ(J,JANG,JFRE),INRJMX(J))
102   CONTINUE
C
C          2. ISOLE UN PIC.
C             ----- -- ---
C
200   CONTINUE
      DO 201 J=KJS,KJL
      KPICF(J)=1
      KPICA(J)=1
201   CONTINUE
      DO 202 JFRE=1,KFRE
      DO 202 JANG=1,KANG
      DO 202 J=KJS,KJL
      IF(KNRJ(J,JANG,JFRE).EQ.INRJMX(J)) THEN
         KPICA(J)=JANG
         KPICF(J)=JFRE
      ENDIF
202   CONTINUE
      RETURN
      END
