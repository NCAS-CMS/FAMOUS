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
C Modelling at the the above address.
C ******************************COPYRIGHT******************************
C
C $Header: /u/um1/vn4.1/mods/source/RCS/anf1f401,v 1.2 1996/06/21 10:13:
!+ Perform multiple fast fourier transforms by calling FTRANS
! Subroutine Interface:
      SUBROUTINE FOURIER(A,TRIGS,IFAX,INC,JUMP,N,LOT,ISIGN)

      IMPLICIT NONE
! Description:
!
!   SUBROUTINE 'FOURIER' - MULTIPLE FAST REAL PERIODIC TRANSFORM
!   UNIFIED MODEL RE-WRITE OF ECMWF ROUTINE FFT991
!
!   REAL TRANSFORM OF LENGTH N PERFORMED BY REMOVING REDUNDANT
!   OPERATIONS FROM COMPLEX TRANSFORM OF LENGTH N
!
!   INPUT INFORMATION:
!   A IS THE ARRAY CONTAINING INPUT & OUTPUT DATA
!   TRIGS IS A PREVIOUSLY PREPARED LIST OF TRIG FUNCTION VALUES
!   IFAX IS A PREVIOUSLY PREPARED LIST OF FACTORS OF N
!   INC IS THE INCREMENT WITHIN EACH DATA 'VECTOR'
!       (E.G. INC=1 FOR CONSECUTIVELY STORED DATA)
!   JUMP IS THE INCREMENT BETWEEN THE START OF EACH DATA VECTOR
!   N IS THE LENGTH OF THE DATA VECTORS
!   LOT IS THE NUMBER OF DATA VECTORS
!   ISIGN = +1 FOR TRANSFORM FROM SPECTRAL TO GRIDPOINT
!         = -1 FOR TRANSFORM FROM GRIDPOINT TO SPECTRAL
!
!   ORDERING OF COEFFICIENTS:
!       A(0),B(0),A(1),B(1),A(2),B(2),...,A(N/2),B(N/2)
!       WHERE B(0)=B(N/2)=0; (N+2) LOCATIONS REQUIRED
!
!   ORDERING OF DATA:
!       X(0),X(1),X(2),...,X(N-1), 0 , 0 ; (N+2) LOCATIONS REQUIRED
!
!   N MUST BE COMPOSED OF FACTORS 2,3 & 5 BUT DOES NOT HAVE TO BE EVEN
!
!   DEFINITION OF TRANSFORMS:
!
!   ISIGN=+1: X(J)=SUM(K=0,...,N-1)(C(K)*EXP(2*I*J*K*PI/N))
!       WHERE C(K)=A(K)+I*B(K) AND C(N-K)=A(K)-I*B(K)
!
!   ISIGN=-1: A(K)=(1/N)*SUM(J=0,...,N-1)(X(J)*COS(2*J*K*PI/N))
!             B(K)=-(1/N)*SUM(J=0,...,N-1)(X(J)*SIN(2*J*K*PI/N))
!---------------------------------------------------------------------
!
! Current code owner: M.H.Mawson
!
! History:
! Version   Date        Comment
! =======   ====        =======
! 4.1       June '96    Original code at 4.1. Based on modifications by
!                       Ken Hawick on public domain provided software.
!                       This is primarily for workstation usage.
!  4.5  07/05/98  Blocking size increased from 64 to 512 to give
!                 much better vector performance on Fujitsu VPP700
!                                               RBarnes@ecmwf.int
!
! Code description:
!   FORTRAN 77 + common Fortran 90 extensions.
!   Written to UM programming standards version 7.
!   DOCUMENTATION:        NIL.
!END------------------------------------------------------------------

! Subroutine arguments
!
!   Scaler arguments with intent(in):
      INTEGER
     &     INC,        ! IN Increment between elements of data vector
     &     JUMP,       ! IN Increment between start of each data vector
     &     N,          ! IN Length of data vector in grid-point space
     &                 !      without extra zeroes
     &     LOT,        ! IN Number of data vectors
     &     ISIGN,      ! IN Determines type of transform
     &     IFAX(10)    ! IN List of factors of n

!   Array arguments with intent(in):
      REAL TRIGS(N)    ! IN Trigonometrical functions

!   Array arguments with intent(out):
      REAL A(JUMP*LOT) ! INOUT Data


      REAL WORK((N+2)*64) ! General workspace

! local scalers:
      INTEGER NFAX,       ! NUMBER OF FACTORS
     &        NX,         ! N+1 EXCEPT WHERE N IS ODD THEN HOLDS N
     &        NBLOX,      ! NUMBER OF BLOCKS LOT IS SPLIT INTO
     &        NB,         ! DO LOOP COUNTER
     &        ISTART,     ! START ADDRESS FOR A BLOCK
     &        NVEX,       ! NUMBER OF ELEMENTS IN VECTOR
     &        IA,         ! USED TO PASS ISTART TO FTRANS
     &        IX,         ! VARIABLE USED FOR ADDRESSING
     &        LA,         ! VARIABLE USED FOR ADDRESSING
     &        IGO,        ! A CONTROL VARIABLE
     &        K,          ! DO LOOP COUNTER
     &        IFAC,       ! HOLDS CURRENT FACTOR
     &        IERR        ! HOLDS ERROR STATUS

      INTEGER I,J,II,IZ,JJ,IBASE,JBASE  ! loop/indexing variables.


! Function and subroutine calls:
      EXTERNAL FTRANS

!- End of Header --------------------------------------------------

C------------------------------------------------------
C Section 1. Set up information for sections 2 and 3:
C------------------------------------------------------

C Set number of factors and NX:
      NFAX=IFAX(1)
      NX=N+1
      IF (MOD(N,2).EQ.1) NX=N

C Calculate number of blocks of 64 data vectors are to be
C split into:
      NBLOX=1+(LOT-1)/64
      NVEX=LOT-(NBLOX-1)*64

C------------------------------------------------------
C Section 2. ISIGN=+1, spectral to gridpoint transformation
C------------------------------------------------------

      IF (ISIGN.EQ.1) THEN  ! spectral-to-gridpoint transform:
        ISTART=1
        DO NB=1,NBLOX
          IA=ISTART
          I=ISTART
          DO J=1,NVEX
            A(I+INC)=0.5*A(I)
            I=I+JUMP
          ENDDO
          IF (MOD(N,2).NE.1) THEN
            I=ISTART+N*INC
            DO J=1,NVEX
              A(I)=0.5*A(I)
              I=I+JUMP
            ENDDO
          END IF
          IA=ISTART+INC
          LA=1
          IGO=1

          DO K=1,NFAX
            IFAC=IFAX(K+1)
            IERR=-1
            IF (IGO.EQ.1) THEN   !  Invoke Fourier Synthesis pass
              CALL FTRANS(-1,A(IA),A(IA+LA*INC),WORK(1),WORK(IFAC*LA+1),
     &                     TRIGS,INC,1,JUMP,NX,NVEX,N,IFAC,LA,IERR)
            ELSE
              CALL FTRANS(-1,WORK(1),WORK(LA+1),A(IA),A(IA+IFAC*LA*INC),
     &                     TRIGS,1,INC,NX,JUMP,NVEX,N,IFAC,LA,IERR)
            END IF
C           IF (IERR.NE.0) GO TO 400
            LA=IFAC*LA
            IGO=-IGO
            IA=ISTART
          ENDDO

C If necessary, copy results back to A:
          IF (MOD(NFAX,2).NE.0) THEN
            IBASE=1
            JBASE=IA
            DO JJ=1,NVEX
              I=IBASE
              J=JBASE
              DO II=1,N
                A(J)=WORK(I)
                I=I+1
                J=J+INC
              ENDDO
              IBASE=IBASE+NX
              JBASE=JBASE+JUMP
            ENDDO
          END IF

C Fill in zeros at end:
          IX=ISTART+N*INC
          DO J=1,NVEX
            A(IX)=0.0
            A(IX+INC)=0.0
            IX=IX+JUMP
          ENDDO
          ISTART=ISTART+NVEX*JUMP
          NVEX=64
        ENDDO


      ELSE  ! isign=-1, gridpoint-to-spectral transform

C------------------------------------------------------
C Section 3: ISIGN=-1, gridpoint to spectral transform
C------------------------------------------------------

        ISTART=1
        DO NB=1,NBLOX
          IA=ISTART
          LA=N
          IGO=+1

          DO K=1,NFAX
            IFAC=IFAX(NFAX+2-K)
            LA=LA/IFAC
            IERR=-1
            IF (IGO.EQ.1) THEN ! Invoke Fourier analysis pass
              CALL FTRANS(1,A(IA),A(IA+IFAC*LA*INC),WORK(1),WORK(LA+1),
     &                     TRIGS,INC,1,JUMP,NX,NVEX,N,IFAC,LA,IERR)
            ELSE
              CALL FTRANS(1,WORK(1),WORK(IFAC*LA+1),A(IA),A(IA+LA*INC),
     &                      TRIGS,1,INC,NX,JUMP,NVEX,N,IFAC,LA,IERR)
            END IF
C           IF (IERR.NE.0) GO TO 500
            IGO=-IGO
            IA=ISTART+INC
          ENDDO

C If necessary, copy results back to A:
          IF (MOD(NFAX,2).NE.0) THEN
            IBASE=1
            JBASE=IA
            DO JJ=1,NVEX
              I=IBASE
              J=JBASE
              DO II=1,N
                A(J)=WORK(I)
                I=I+1
                J=J+INC
              ENDDO
              IBASE=IBASE+NX
              JBASE=JBASE+JUMP
            ENDDO
          END IF

C Shift A(0) and fill in zero imaginary parts:
          IX=ISTART
          DO J=1,NVEX
            A(IX)=A(IX+INC)
            A(IX+INC)=0.0
            IX=IX+JUMP
          ENDDO
          IF (MOD(N,2).NE.1) THEN
            IZ=ISTART+(N+1)*INC
            DO J=1,NVEX
              A(IZ)=0.0
              IZ=IZ+JUMP
            ENDDO
          END IF

          ISTART=ISTART+NVEX*JUMP
          NVEX=64
        ENDDO
      END IF

C Error messages:
C 400 CONTINUE
C     IF(IERR.NE.0) THEN
C       IF(IERR.EQ.1) THEN
C         WRITE(6,410) NVEX
C 410     FORMAT(16H1VECTOR LENGTH =,I4,17H, GREATER THAN 64)
C       ELSE IF(IERR.EQ.2) THEN
C         WRITE(6,420) IFAC
C 420     FORMAT( 9H1FACTOR =,I3,17H, NOT CATERED FOR)
C       ELSE IF(IERR.EQ.3) THEN
C         WRITE(6,430) IFAC
C 430     FORMAT(9H1FACTOR =,I3,31H, ONLY CATERED FOR IF LA*IFAC=N)
C       ELSE
C         WRITE(4,440) IFAC
C 440     FORMAT(' UNRECOGNISED ERROR MESSAGE, CODE ',I3)
C       END IF
C     END IF

C End of routine FOURIER

      RETURN
      END

!- End of subroutine code-----------------------------------------

C-----------------------------------------------------------------------
C Subroutine FTRANS
C
C $Header: /u/um1/vn4.1/mods/source/RCS/anf1f401,v 1.2 1996/06/21 10:13:
C-----------------------------------------------------------------------
C  Fourier transform:
C
!+ Public Domain provided Fourier transform routine.
! Subroutine Interface:
      SUBROUTINE FTRANS(ICTL,A,B,C,D,TRIGS,INC1,INC2,INC3,INC4,LOT,N,
     &                  IFAC, LA,IERR)
! Description:
! Calculates Fourier Transforms.
!
! Current code owner: Public Domain.
!
! History:
! Version   Date        Comment
! =======   ====        =======
! 4.1       June '96    Original code at 4.1.
!                       Public domain provided software.
!                       This is primarily for workstation usage.
!
! Code description:
!   FORTRAN 77 + common Fortran 90 extensions.
!   Written to UM programming standards version 7.
!   DOCUMENTATION:        NIL.
!END------------------------------------------------------------------

      IMPLICIT NONE
!
! Subroutine arguments

      INTEGER ICTL   ! Control:  1 = analysis; -1 = synthesis


      INTEGER INC1,  ! ADDRESSING INCREMENT FOR A
     &        INC2,  ! ADDRESSING INCREMENT FOR C
     &        INC3,  ! INCREMENT BETWEEN INPUT VECTORS A
     &        INC4,  ! INCREMENT BETWEEN INPUT VECTORS C
     &        LOT,   ! NUMBER OF VECTORS
     &        N,     ! LENGTH OF THE VECTORS
     &        IFAC,  ! CURRENT FACTOR OF N
     &        LA,    ! N/(PRODUCT OF FACTORS USED SO FAR)
     &        IERR   ! Error INDICATOR:
                     !   0 - PASS COMPLETED WITHOUT ERROR
                     !   1 - LOT GREATER THAN 64
                     !   2 - IFAC NOT CATERED FOR
                     !   3 - IFAC ONLY CATERED FOR IF LA=N/IFAC
      REAL A(N),     ! First real input vector
     &     B(N),
     &     C(N),     ! First real output vector
     &     D(N),
     &     TRIGS(N)  ! Precalculated list of sines & cosines

C  for Fourier analysis:
C     A IS FIRST REAL INPUT VECTOR
C         EQUIVALENCE B(1) WITH A(IFAC*LA*INC1+1)
C     C IS FIRST REAL OUTPUT VECTOR
C         EQUIVALENCE D(1) WITH C(LA*INC2+1)
C
C  or for synthesis:
C     A IS FIRST REAL INPUT VECTOR
C         EQUIVALENCE B(1) WITH A (LA*INC1+1)
C     C IS FIRST REAL OUTPUT VECTOR
C         EQUIVALENCE D(1) WITH C(IFAC*LA*INC2+1)

C-----------------------------------------------------------------------

      INTEGER IINK,
     &        JINK,
     &        IJUMP, JUMP,
     &        KSTOP,
     &        IBASE,
     &        JBASE,
     &        IBAD,
     &        IGO,
     &        I, J, K, L, M,
     &            KB, KC, KD, KE, KF,
     &        IA, IB, IC, ID, IE, IF, IG, IH,
     &        JA, JB, JC, JD, JE, JF, JG, JH,
     &        IJK

      REAL A0,  A1,  A2,  A3,  A4,  A5,  A6
      REAL A10, A11
      REAL A20, A21

      REAL B0,  B1,  B2,  B3,  B4,  B5,  B6
      REAL B10, B11
      REAL B20, B21

      REAL C1, C2, C3, C4, C5
      REAL S1, S2, S3, S4, S5

      REAL Z, ZQRT5, ZSIN36, ZSIN45, ZSIN60, ZSIN72
      REAL QQRT5, SSIN36, SSIN45, SSIN60, SSIN72

      REAL AA10(64),AA11(64),AA20(64),AA21(64),
     &     BB10(64),BB11(64),BB20(64),BB21(64)

!      DOUBLE PRECISION SIN36, SIN45, SIN72, SIN60, QRT5
      REAL SIN36, SIN45, SIN72, SIN60, QRT5

      DATA SIN36/0.587785252292473/,SIN72/0.951056516295154/,
     &     QRT5/0.559016994374947/,SIN60/0.866025403784437/

!- End of Header --------------------------------------------------

      IF( ICTL .EQ. 1 )THEN  !  Do Fourier Analysis:


      M=N/IFAC
      IINK=LA*INC1
      JINK=LA*INC2
      IJUMP=(IFAC-1)*IINK
      KSTOP=(N-IFAC)/(2*IFAC)

      IBAD=1
      IF (LOT.GT.64) GO TO 910
      IBASE=0
      JBASE=0
      IGO=IFAC-1
      IF (IGO.EQ.7) IGO=6
      IBAD=2
      IF (IGO.LT.1.OR.IGO.GT.6) GO TO 910
      GO TO (200,300,400,500,600,800),IGO


C     CODING FOR FACTOR 2
 200  CONTINUE
      IA=1
      IB=IA+IINK
      JA=1
      JB=JA+(2*M-LA)*INC2

      IF (LA.EQ.M) GO TO 290

      DO 220 L=1,LA
         I=IBASE
         J=JBASE
CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
         DO 210 IJK=1,LOT
            C(JA+J)=A(IA+I)+A(IB+I)
            C(JB+J)=A(IA+I)-A(IB+I)
            I=I+INC3
            J=J+INC4
 210     CONTINUE
         IBASE=IBASE+INC1
         JBASE=JBASE+INC2
 220  CONTINUE
      JA=JA+JINK
      JINK=2*JINK
      JB=JB-JINK
      IBASE=IBASE+IJUMP
      IJUMP=2*IJUMP+IINK
      IF (JA.EQ.JB) GO TO 260
      DO 250 K=LA,KSTOP,LA
         KB=K+K
         C1=TRIGS(KB+1)
         S1=TRIGS(KB+2)
         JBASE=0
         DO 240 L=1,LA
            I=IBASE
            J=JBASE
CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
            DO 230 IJK=1,LOT
               C(JA+J)=A(IA+I)+(C1*A(IB+I)+S1*B(IB+I))
               C(JB+J)=A(IA+I)-(C1*A(IB+I)+S1*B(IB+I))
               D(JA+J)=(C1*B(IB+I)-S1*A(IB+I))+B(IA+I)
               D(JB+J)=(C1*B(IB+I)-S1*A(IB+I))-B(IA+I)
               I=I+INC3
               J=J+INC4
 230        CONTINUE
            IBASE=IBASE+INC1
            JBASE=JBASE+INC2
 240     CONTINUE
         IBASE=IBASE+IJUMP
         JA=JA+JINK
         JB=JB-JINK
 250  CONTINUE
      IF (JA.GT.JB) GO TO 900
 260  CONTINUE
      JBASE=0
      DO 280 L=1,LA
         I=IBASE
         J=JBASE
CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
         DO 270 IJK=1,LOT
            C(JA+J)=A(IA+I)
            D(JA+J)=-A(IB+I)
            I=I+INC3
            J=J+INC4
 270     CONTINUE
         IBASE=IBASE+INC1
         JBASE=JBASE+INC2
 280  CONTINUE
      GO TO 900

 290  CONTINUE
      Z=1.0/FLOAT(N)
      DO 294 L=1,LA
         I=IBASE
         J=JBASE
CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
         DO 292 IJK=1,LOT
            C(JA+J)=Z*(A(IA+I)+A(IB+I))
            C(JB+J)=Z*(A(IA+I)-A(IB+I))
            I=I+INC3
            J=J+INC4
 292     CONTINUE
         IBASE=IBASE+INC1
         JBASE=JBASE+INC2
 294  CONTINUE
      GO TO 900

C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C Coding for factor 3:
 300  CONTINUE
      IA=1
      IB=IA+IINK
      IC=IB+IINK
      JA=1
      JB=JA+(2*M-LA)*INC2
      JC=JB

      IF (LA.EQ.M) GO TO 390

      DO 320 L=1,LA
         I=IBASE
         J=JBASE
CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
         DO 310 IJK=1,LOT
            C(JA+J)=A(IA+I)+(A(IB+I)+A(IC+I))
            C(JB+J)=A(IA+I)-0.5*(A(IB+I)+A(IC+I))
            D(JB+J)=SIN60*(A(IC+I)-A(IB+I))
            I=I+INC3
            J=J+INC4
 310     CONTINUE
         IBASE=IBASE+INC1
         JBASE=JBASE+INC2
 320  CONTINUE
      JA=JA+JINK
      JINK=2*JINK
      JB=JB+JINK
      JC=JC-JINK
      IBASE=IBASE+IJUMP
      IJUMP=2*IJUMP+IINK
      IF (JA.EQ.JC) GO TO 360
      DO 350 K=LA,KSTOP,LA
         KB=K+K
         KC=KB+KB
         C1=TRIGS(KB+1)
         S1=TRIGS(KB+2)
         C2=TRIGS(KC+1)
         S2=TRIGS(KC+2)
         JBASE=0
         DO 340 L=1,LA
            I=IBASE
            J=JBASE
CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
            DO 330 IJK=1,LOT
               A1=(C1*A(IB+I)+S1*B(IB+I))+(C2*A(IC+I)+S2*B(IC+I))
               B1=(C1*B(IB+I)-S1*A(IB+I))+(C2*B(IC+I)-S2*A(IC+I))
               A2=A(IA+I)-0.5*A1
               B2=B(IA+I)-0.5*B1
               A3=SIN60*
     $              ((C1*A(IB+I)+S1*B(IB+I))-(C2*A(IC+I)+S2*B(IC+I)))
               B3=SIN60*
     $              ((C1*B(IB+I)-S1*A(IB+I))-(C2*B(IC+I)-S2*A(IC+I)))
               C(JA+J)=A(IA+I)+A1
               D(JA+J)=B(IA+I)+B1
               C(JB+J)=A2+B3
               D(JB+J)=B2-A3
               C(JC+J)=A2-B3
               D(JC+J)=-(B2+A3)
               I=I+INC3
               J=J+INC4
 330        CONTINUE
            IBASE=IBASE+INC1
            JBASE=JBASE+INC2
 340     CONTINUE
         IBASE=IBASE+IJUMP
         JA=JA+JINK
         JB=JB+JINK
         JC=JC-JINK
 350  CONTINUE
      IF (JA.GT.JC) GO TO 900
 360  CONTINUE
      JBASE=0
      DO 380 L=1,LA
         I=IBASE
         J=JBASE
CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
         DO 370 IJK=1,LOT
            C(JA+J)=A(IA+I)+0.5*(A(IB+I)-A(IC+I))
            D(JA+J)=-SIN60*(A(IB+I)+A(IC+I))
            C(JB+J)=A(IA+I)-(A(IB+I)-A(IC+I))
            I=I+INC3
            J=J+INC4
 370     CONTINUE
         IBASE=IBASE+INC1
         JBASE=JBASE+INC2
 380  CONTINUE
      GO TO 900

 390  CONTINUE
      Z=1.0/FLOAT(N)
      ZSIN60=Z*SIN60
      DO 394 L=1,LA
         I=IBASE
         J=JBASE
CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
         DO 392 IJK=1,LOT
            C(JA+J)=Z*(A(IA+I)+(A(IB+I)+A(IC+I)))
            C(JB+J)=Z*(A(IA+I)-0.5*(A(IB+I)+A(IC+I)))
            D(JB+J)=ZSIN60*(A(IC+I)-A(IB+I))
            I=I+INC3
            J=J+INC4
 392     CONTINUE
         IBASE=IBASE+INC1
         JBASE=JBASE+INC2
 394  CONTINUE
      GO TO 900

C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C Coding for factor 4
 400  CONTINUE
      IA=1
      IB=IA+IINK
      IC=IB+IINK
      ID=IC+IINK
      JA=1
      JB=JA+(2*M-LA)*INC2
      JC=JB+2*M*INC2
      JD=JB

      IF (LA.EQ.M) GO TO 490

      DO 420 L=1,LA
         I=IBASE
         J=JBASE
CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
         DO 410 IJK=1,LOT
            C(JA+J)=(A(IA+I)+A(IC+I))+(A(IB+I)+A(ID+I))
            C(JC+J)=(A(IA+I)+A(IC+I))-(A(IB+I)+A(ID+I))
            C(JB+J)=A(IA+I)-A(IC+I)
            D(JB+J)=A(ID+I)-A(IB+I)
            I=I+INC3
            J=J+INC4
 410     CONTINUE
         IBASE=IBASE+INC1
         JBASE=JBASE+INC2
 420  CONTINUE
      JA=JA+JINK
      JINK=2*JINK
      JB=JB+JINK
      JC=JC-JINK
      JD=JD-JINK
      IBASE=IBASE+IJUMP
      IJUMP=2*IJUMP+IINK
      IF (JB.EQ.JC) GO TO 460
      DO 450 K=LA,KSTOP,LA
         KB=K+K
         KC=KB+KB
         KD=KC+KB
         C1=TRIGS(KB+1)
         S1=TRIGS(KB+2)
         C2=TRIGS(KC+1)
         S2=TRIGS(KC+2)
         C3=TRIGS(KD+1)
         S3=TRIGS(KD+2)
         JBASE=0
         DO 440 L=1,LA
            I=IBASE
            J=JBASE
CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
            DO 430 IJK=1,LOT
               A0=A(IA+I)+(C2*A(IC+I)+S2*B(IC+I))
               A2=A(IA+I)-(C2*A(IC+I)+S2*B(IC+I))
               A1=(C1*A(IB+I)+S1*B(IB+I))+(C3*A(ID+I)+S3*B(ID+I))
               A3=(C1*A(IB+I)+S1*B(IB+I))-(C3*A(ID+I)+S3*B(ID+I))
               B0=B(IA+I)+(C2*B(IC+I)-S2*A(IC+I))
               B2=B(IA+I)-(C2*B(IC+I)-S2*A(IC+I))
               B1=(C1*B(IB+I)-S1*A(IB+I))+(C3*B(ID+I)-S3*A(ID+I))
               B3=(C1*B(IB+I)-S1*A(IB+I))-(C3*B(ID+I)-S3*A(ID+I))
               C(JA+J)=A0+A1
               C(JC+J)=A0-A1
               D(JA+J)=B0+B1
               D(JC+J)=B1-B0
               C(JB+J)=A2+B3
               C(JD+J)=A2-B3
               D(JB+J)=B2-A3
               D(JD+J)=-(B2+A3)
               I=I+INC3
               J=J+INC4
 430        CONTINUE
            IBASE=IBASE+INC1
            JBASE=JBASE+INC2
 440     CONTINUE
         IBASE=IBASE+IJUMP
         JA=JA+JINK
         JB=JB+JINK
         JC=JC-JINK
         JD=JD-JINK
 450  CONTINUE
      IF (JB.GT.JC) GO TO 900
 460  CONTINUE
      SIN45=SQRT(0.5)
      JBASE=0
      DO 480 L=1,LA
         I=IBASE
         J=JBASE
CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
         DO 470 IJK=1,LOT
            C(JA+J)=A(IA+I)+SIN45*(A(IB+I)-A(ID+I))
            C(JB+J)=A(IA+I)-SIN45*(A(IB+I)-A(ID+I))
            D(JA+J)=-A(IC+I)-SIN45*(A(IB+I)+A(ID+I))
            D(JB+J)=A(IC+I)-SIN45*(A(IB+I)+A(ID+I))
            I=I+INC3
            J=J+INC4
 470     CONTINUE
         IBASE=IBASE+INC1
         JBASE=JBASE+INC2
 480  CONTINUE
      GO TO 900
C
 490  CONTINUE
      Z=1.0/FLOAT(N)
      DO 494 L=1,LA
         I=IBASE
         J=JBASE
CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
         DO 492 IJK=1,LOT
            C(JA+J)=Z*((A(IA+I)+A(IC+I))+(A(IB+I)+A(ID+I)))
            C(JC+J)=Z*((A(IA+I)+A(IC+I))-(A(IB+I)+A(ID+I)))
            C(JB+J)=Z*(A(IA+I)-A(IC+I))
            D(JB+J)=Z*(A(ID+I)-A(IB+I))
            I=I+INC3
            J=J+INC4
 492     CONTINUE
         IBASE=IBASE+INC1
         JBASE=JBASE+INC2
 494  CONTINUE
      GO TO 900

C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C Coding for factor 5
 500  CONTINUE
      IA=1
      IB=IA+IINK
      IC=IB+IINK
      ID=IC+IINK
      IE=ID+IINK
      JA=1
      JB=JA+(2*M-LA)*INC2
      JC=JB+2*M*INC2
      JD=JC
      JE=JB

      IF (LA.EQ.M) GO TO 590

      DO 520 L=1,LA
         I=IBASE
         J=JBASE
CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
         DO 510 IJK=1,LOT
            A1=A(IB+I)+A(IE+I)
            A3=A(IB+I)-A(IE+I)
            A2=A(IC+I)+A(ID+I)
            A4=A(IC+I)-A(ID+I)
            A5=A(IA+I)-0.25*(A1+A2)
            A6=QRT5*(A1-A2)
            C(JA+J)=A(IA+I)+(A1+A2)
            C(JB+J)=A5+A6
            C(JC+J)=A5-A6
            D(JB+J)=-SIN72*A3-SIN36*A4
            D(JC+J)=-SIN36*A3+SIN72*A4
            I=I+INC3
            J=J+INC4
 510     CONTINUE
         IBASE=IBASE+INC1
         JBASE=JBASE+INC2
 520  CONTINUE
      JA=JA+JINK
      JINK=2*JINK
      JB=JB+JINK
      JC=JC+JINK
      JD=JD-JINK
      JE=JE-JINK
      IBASE=IBASE+IJUMP
      IJUMP=2*IJUMP+IINK
      IF (JB.EQ.JD) GO TO 560
      DO 550 K=LA,KSTOP,LA
         KB=K+K
         KC=KB+KB
         KD=KC+KB
         KE=KD+KB
         C1=TRIGS(KB+1)
         S1=TRIGS(KB+2)
         C2=TRIGS(KC+1)
         S2=TRIGS(KC+2)
         C3=TRIGS(KD+1)
         S3=TRIGS(KD+2)
         C4=TRIGS(KE+1)
         S4=TRIGS(KE+2)
         JBASE=0
         DO 540 L=1,LA
            I=IBASE
            J=JBASE
CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
            DO 530 IJK=1,LOT
               A1=(C1*A(IB+I)+S1*B(IB+I))+(C4*A(IE+I)+S4*B(IE+I))
               A3=(C1*A(IB+I)+S1*B(IB+I))-(C4*A(IE+I)+S4*B(IE+I))
               A2=(C2*A(IC+I)+S2*B(IC+I))+(C3*A(ID+I)+S3*B(ID+I))
               A4=(C2*A(IC+I)+S2*B(IC+I))-(C3*A(ID+I)+S3*B(ID+I))
               B1=(C1*B(IB+I)-S1*A(IB+I))+(C4*B(IE+I)-S4*A(IE+I))
               B3=(C1*B(IB+I)-S1*A(IB+I))-(C4*B(IE+I)-S4*A(IE+I))
               B2=(C2*B(IC+I)-S2*A(IC+I))+(C3*B(ID+I)-S3*A(ID+I))
               B4=(C2*B(IC+I)-S2*A(IC+I))-(C3*B(ID+I)-S3*A(ID+I))
               A5=A(IA+I)-0.25*(A1+A2)
               A6=QRT5*(A1-A2)
               B5=B(IA+I)-0.25*(B1+B2)
               B6=QRT5*(B1-B2)
               A10=A5+A6
               A20=A5-A6
               B10=B5+B6
               B20=B5-B6
               A11=SIN72*B3+SIN36*B4
               A21=SIN36*B3-SIN72*B4
               B11=SIN72*A3+SIN36*A4
               B21=SIN36*A3-SIN72*A4
               C(JA+J)=A(IA+I)+(A1+A2)
               C(JB+J)=A10+A11
               C(JE+J)=A10-A11
               C(JC+J)=A20+A21
               C(JD+J)=A20-A21
               D(JA+J)=B(IA+I)+(B1+B2)
               D(JB+J)=B10-B11
               D(JE+J)=-(B10+B11)
               D(JC+J)=B20-B21
               D(JD+J)=-(B20+B21)
               I=I+INC3
               J=J+INC4
 530        CONTINUE
            IBASE=IBASE+INC1
            JBASE=JBASE+INC2
 540     CONTINUE
         IBASE=IBASE+IJUMP
         JA=JA+JINK
         JB=JB+JINK
         JC=JC+JINK
         JD=JD-JINK
         JE=JE-JINK
 550  CONTINUE
      IF (JB.GT.JD) GO TO 900
 560  CONTINUE
      JBASE=0
      DO 580 L=1,LA
         I=IBASE
         J=JBASE
CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
         DO 570 IJK=1,LOT
            A1=A(IB+I)+A(IE+I)
            A3=A(IB+I)-A(IE+I)
            A2=A(IC+I)+A(ID+I)
            A4=A(IC+I)-A(ID+I)
            A5=A(IA+I)+0.25*(A3-A4)
            A6=QRT5*(A3+A4)
            C(JA+J)=A5+A6
            C(JB+J)=A5-A6
            C(JC+J)=A(IA+I)-(A3-A4)
            D(JA+J)=-SIN36*A1-SIN72*A2
            D(JB+J)=-SIN72*A1+SIN36*A2
            I=I+INC3
            J=J+INC4
 570     CONTINUE
         IBASE=IBASE+INC1
         JBASE=JBASE+INC2
 580  CONTINUE
      GO TO 900
C
 590  CONTINUE
      Z=1.0/FLOAT(N)
      ZQRT5=Z*QRT5
      ZSIN36=Z*SIN36
      ZSIN72=Z*SIN72
      DO 594 L=1,LA
         I=IBASE
         J=JBASE
CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
         DO 592 IJK=1,LOT
            A1=A(IB+I)+A(IE+I)
            A3=A(IB+I)-A(IE+I)
            A2=A(IC+I)+A(ID+I)
            A4=A(IC+I)-A(ID+I)
            A5=Z*(A(IA+I)-0.25*(A1+A2))
            A6=ZQRT5*(A1-A2)
            C(JA+J)=Z*(A(IA+I)+(A1+A2))
            C(JB+J)=A5+A6
            C(JC+J)=A5-A6
            D(JB+J)=-ZSIN72*A3-ZSIN36*A4
            D(JC+J)=-ZSIN36*A3+ZSIN72*A4
            I=I+INC3
            J=J+INC4
 592     CONTINUE
         IBASE=IBASE+INC1
         JBASE=JBASE+INC2
 594  CONTINUE
      GO TO 900

C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C Coding for factor 6
 600  CONTINUE
      IA=1
      IB=IA+IINK
      IC=IB+IINK
      ID=IC+IINK
      IE=ID+IINK
      IF=IE+IINK
      JA=1
      JB=JA+(2*M-LA)*INC2
      JC=JB+2*M*INC2
      JD=JC+2*M*INC2
      JE=JC
      JF=JB
C
      IF (LA.EQ.M) GO TO 690
C
      DO 620 L=1,LA
         I=IBASE
         J=JBASE
CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
         DO 610 IJK=1,LOT
            A11=(A(IC+I)+A(IF+I))+(A(IB+I)+A(IE+I))
            C(JA+J)=(A(IA+I)+A(ID+I))+A11
            C(JC+J)=(A(IA+I)+A(ID+I)-0.5*A11)
            D(JC+J)=SIN60*((A(IC+I)+A(IF+I))-(A(IB+I)+A(IE+I)))
            A11=(A(IC+I)-A(IF+I))+(A(IE+I)-A(IB+I))
            C(JB+J)=(A(IA+I)-A(ID+I))-0.5*A11
            D(JB+J)=SIN60*((A(IE+I)-A(IB+I))-(A(IC+I)-A(IF+I)))
            C(JD+J)=(A(IA+I)-A(ID+I))+A11
            I=I+INC3
            J=J+INC4
 610     CONTINUE
         IBASE=IBASE+INC1
         JBASE=JBASE+INC2
 620  CONTINUE
      JA=JA+JINK
      JINK=2*JINK
      JB=JB+JINK
      JC=JC+JINK
      JD=JD-JINK
      JE=JE-JINK
      JF=JF-JINK
      IBASE=IBASE+IJUMP
      IJUMP=2*IJUMP+IINK
      IF (JC.EQ.JD) GO TO 660
      DO 650 K=LA,KSTOP,LA
         KB=K+K
         KC=KB+KB
         KD=KC+KB
         KE=KD+KB
         KF=KE+KB
         C1=TRIGS(KB+1)
         S1=TRIGS(KB+2)
         C2=TRIGS(KC+1)
         S2=TRIGS(KC+2)
         C3=TRIGS(KD+1)
         S3=TRIGS(KD+2)
         C4=TRIGS(KE+1)
         S4=TRIGS(KE+2)
         C5=TRIGS(KF+1)
         S5=TRIGS(KF+2)
         JBASE=0
         DO 640 L=1,LA
            I=IBASE
            J=JBASE
CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
            DO 630 IJK=1,LOT
               A1=C1*A(IB+I)+S1*B(IB+I)
               B1=C1*B(IB+I)-S1*A(IB+I)
               A2=C2*A(IC+I)+S2*B(IC+I)
               B2=C2*B(IC+I)-S2*A(IC+I)
               A3=C3*A(ID+I)+S3*B(ID+I)
               B3=C3*B(ID+I)-S3*A(ID+I)
               A4=C4*A(IE+I)+S4*B(IE+I)
               B4=C4*B(IE+I)-S4*A(IE+I)
               A5=C5*A(IF+I)+S5*B(IF+I)
               B5=C5*B(IF+I)-S5*A(IF+I)
               A11=(A2+A5)+(A1+A4)
               A20=(A(IA+I)+A3)-0.5*A11
               A21=SIN60*((A2+A5)-(A1+A4))
               B11=(B2+B5)+(B1+B4)
               B20=(B(IA+I)+B3)-0.5*B11
               B21=SIN60*((B2+B5)-(B1+B4))
               C(JA+J)=(A(IA+I)+A3)+A11
               D(JA+J)=(B(IA+I)+B3)+B11
               C(JC+J)=A20-B21
               D(JC+J)=A21+B20
               C(JE+J)=A20+B21
               D(JE+J)=A21-B20
               A11=(A2-A5)+(A4-A1)
               A20=(A(IA+I)-A3)-0.5*A11
               A21=SIN60*((A4-A1)-(A2-A5))
               B11=(B5-B2)-(B4-B1)
               B20=(B3-B(IA+I))-0.5*B11
               B21=SIN60*((B5-B2)+(B4-B1))
               C(JB+J)=A20-B21
               D(JB+J)=A21-B20
               C(JD+J)=A11+(A(IA+I)-A3)
               D(JD+J)=B11+(B3-B(IA+I))
               C(JF+J)=A20+B21
               D(JF+J)=A21+B20
               I=I+INC3
               J=J+INC4
 630        CONTINUE
            IBASE=IBASE+INC1
            JBASE=JBASE+INC2
 640     CONTINUE
         IBASE=IBASE+IJUMP
         JA=JA+JINK
         JB=JB+JINK
         JC=JC+JINK
         JD=JD-JINK
         JE=JE-JINK
         JF=JF-JINK
 650  CONTINUE
      IF (JC.GT.JD) GO TO 900
 660  CONTINUE
      JBASE=0
      DO 680 L=1,LA
         I=IBASE
         J=JBASE
CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
         DO 670 IJK=1,LOT
            C(JA+J)=(A(IA+I)+0.5*(A(IC+I)-A(IE+I))) +
     $           SIN60*(A(IB+I)-A(IF+I))
            D(JA+J)=-(A(ID+I)+0.5*(A(IB+I)+A(IF+I))) -
     $           SIN60*(A(IC+I)+A(IE+I))
            C(JB+J)=A(IA+I)-(A(IC+I)-A(IE+I))
            D(JB+J)=A(ID+I)-(A(IB+I)+A(IF+I))
            C(JC+J)=(A(IA+I)+0.5*(A(IC+I)-A(IE+I))) -
     $           SIN60*(A(IB+I)-A(IF+I))
            D(JC+J)=-(A(ID+I)+0.5*(A(IB+I)+A(IF+I))) +
     $           SIN60*(A(IC+I)+A(IE+I))
            I=I+INC3
            J=J+INC4
 670     CONTINUE
         IBASE=IBASE+INC1
         JBASE=JBASE+INC2
 680  CONTINUE
      GO TO 900
C
 690  CONTINUE
      Z=1.0/FLOAT(N)
      ZSIN60=Z*SIN60
      DO 694 L=1,LA
         I=IBASE
         J=JBASE
CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
         DO 692 IJK=1,LOT
            A11=(A(IC+I)+A(IF+I))+(A(IB+I)+A(IE+I))
            C(JA+J)=Z*((A(IA+I)+A(ID+I))+A11)
            C(JC+J)=Z*((A(IA+I)+A(ID+I))-0.5*A11)
            D(JC+J)=ZSIN60*((A(IC+I)+A(IF+I))-(A(IB+I)+A(IE+I)))
            A11=(A(IC+I)-A(IF+I))+(A(IE+I)-A(IB+I))
            C(JB+J)=Z*((A(IA+I)-A(ID+I))-0.5*A11)
            D(JB+J)=ZSIN60*((A(IE+I)-A(IB+I))-(A(IC+I)-A(IF+I)))
            C(JD+J)=Z*((A(IA+I)-A(ID+I))+A11)
            I=I+INC3
            J=J+INC4
 692     CONTINUE
         IBASE=IBASE+INC1
         JBASE=JBASE+INC2
 694  CONTINUE
      GO TO 900

C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C Coding for factor 8:
 800  CONTINUE
      IBAD=3
      IF (LA.NE.M) GO TO 910
      IA=1
      IB=IA+IINK
      IC=IB+IINK
      ID=IC+IINK
      IE=ID+IINK
      IF=IE+IINK
      IG=IF+IINK
      IH=IG+IINK
      JA=1
      JB=JA+LA*INC2
      JC=JB+2*M*INC2
      JD=JC+2*M*INC2
      JE=JD+2*M*INC2
      Z=1.0/FLOAT(N)
      ZSIN45=Z*SQRT(0.5)

      DO 820 L=1,LA
         I=IBASE
         J=JBASE
CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
         DO 810 IJK=1,LOT
            C(JA+J)=Z*(((A(IA+I)+A(IE+I))+(A(IC+I)+A(IG+I)))+
     *           ((A(ID+I)+A(IH+I))+(A(IB+I)+A(IF+I))))
            C(JE+J)=Z*(((A(IA+I)+A(IE+I))+(A(IC+I)+A(IG+I)))-
     *           ((A(ID+I)+A(IH+I))+(A(IB+I)+A(IF+I))))
            C(JC+J)=Z*((A(IA+I)+A(IE+I))-(A(IC+I)+A(IG+I)))
            D(JC+J)=Z*((A(ID+I)+A(IH+I))-(A(IB+I)+A(IF+I)))
            C(JB+J)=Z*(A(IA+I)-A(IE+I))
     *           +ZSIN45*((A(IH+I)-A(ID+I))-(A(IF+I)-A(IB+I)))
            C(JD+J)=Z*(A(IA+I)-A(IE+I))
     *           -ZSIN45*((A(IH+I)-A(ID+I))-(A(IF+I)-A(IB+I)))
            D(JB+J)=ZSIN45*((A(IH+I)-A(ID+I))+(A(IF+I)-A(IB+I)))
     *           +Z*(A(IG+I)-A(IC+I))
            D(JD+J)=ZSIN45*((A(IH+I)-A(ID+I))+(A(IF+I)-A(IB+I)))
     *           -Z*(A(IG+I)-A(IC+I))
            I=I+INC3
            J=J+INC4
 810     CONTINUE
         IBASE=IBASE+INC1
         JBASE=JBASE+INC2
 820  CONTINUE

 900  CONTINUE
      IBAD=0
 910  CONTINUE
      IERR=IBAD

C-----------------------------------------------------------------------
      ELSE  !  Do Fourier Synthesis:

      M=N/IFAC
      IINK=LA*INC1
      JINK=LA*INC2
      JUMP=(IFAC-1)*JINK
      KSTOP=(N-IFAC)/(2*IFAC)

      IBAD=1
      IF (LOT.GT.64) GO TO 1910
      IBASE=0
      JBASE=0
      IGO=IFAC-1
      IF (IGO.EQ.7) IGO=6
      IBAD=2
      IF (IGO.LT.1.OR.IGO.GT.6) GO TO 1910
      GO TO (1200,1300,1400,1500,1600,1800),IGO

C Coding for factor 2:

 1200 CONTINUE
      IA=1
      IB=IA+(2*M-LA)*INC1
      JA=1
      JB=JA+JINK

      IF (LA.EQ.M) GO TO 1290

      DO L=1,LA
         I=IBASE
         J=JBASE
CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
         DO IJK=1,LOT
            C(JA+J)=A(IA+I)+A(IB+I)
            C(JB+J)=A(IA+I)-A(IB+I)
            I=I+INC3
            J=J+INC4
         ENDDO
         IBASE=IBASE+INC1
         JBASE=JBASE+INC2
      ENDDO
      IA=IA+IINK
      IINK=2*IINK
      IB=IB-IINK
      IBASE=0
      JBASE=JBASE+JUMP
      JUMP=2*JUMP+JINK
      IF (IA.EQ.IB) GO TO 1260
      DO K=LA,KSTOP,LA
         KB=K+K
         C1=TRIGS(KB+1)
         S1=TRIGS(KB+2)
         IBASE=0
         DO L=1,LA
            I=IBASE
            J=JBASE
CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
            DO IJK=1,LOT
               C(JA+J)=A(IA+I)+A(IB+I)
               D(JA+J)=B(IA+I)-B(IB+I)
               C(JB+J)=C1*(A(IA+I)-A(IB+I))-S1*(B(IA+I)+B(IB+I))
               D(JB+J)=S1*(A(IA+I)-A(IB+I))+C1*(B(IA+I)+B(IB+I))
               I=I+INC3
               J=J+INC4
            ENDDO
            IBASE=IBASE+INC1
            JBASE=JBASE+INC2
         ENDDO
         IA=IA+IINK
         IB=IB-IINK
         JBASE=JBASE+JUMP
      ENDDO
      IF (IA.GT.IB) GO TO 1900
 1260 CONTINUE
      IBASE=0
      DO L=1,LA
         I=IBASE
         J=JBASE
CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
         DO IJK=1,LOT
            C(JA+J)=A(IA+I)
            C(JB+J)=-B(IA+I)
            I=I+INC3
            J=J+INC4
         ENDDO
         IBASE=IBASE+INC1
         JBASE=JBASE+INC2
      ENDDO
      GO TO 1900

 1290 CONTINUE
      DO L=1,LA
         I=IBASE
         J=JBASE
CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
         DO IJK=1,LOT
            C(JA+J)=2.0*(A(IA+I)+A(IB+I))
            C(JB+J)=2.0*(A(IA+I)-A(IB+I))
            I=I+INC3
            J=J+INC4
         ENDDO
         IBASE=IBASE+INC1
         JBASE=JBASE+INC2
      ENDDO
      GO TO 1900

C Coding for factor 3:
 1300 CONTINUE
      IA=1
      IB=IA+(2*M-LA)*INC1
      IC=IB
      JA=1
      JB=JA+JINK
      JC=JB+JINK

      IF (LA.EQ.M) GO TO 1390

      DO L=1,LA
         I=IBASE
         J=JBASE
CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
         DO IJK=1,LOT
            C(JA+J)=A(IA+I)+A(IB+I)
            C(JB+J)=(A(IA+I)-0.5*A(IB+I))-(SIN60*(B(IB+I)))
            C(JC+J)=(A(IA+I)-0.5*A(IB+I))+(SIN60*(B(IB+I)))
            I=I+INC3
            J=J+INC4
         ENDDO
         IBASE=IBASE+INC1
         JBASE=JBASE+INC2
      ENDDO
      IA=IA+IINK
      IINK=2*IINK
      IB=IB+IINK
      IC=IC-IINK
      JBASE=JBASE+JUMP
      JUMP=2*JUMP+JINK
      IF (IA.EQ.IC) GO TO 1360
      DO K=LA,KSTOP,LA
         KB=K+K
         KC=KB+KB
         C1=TRIGS(KB+1)
         S1=TRIGS(KB+2)
         C2=TRIGS(KC+1)
         S2=TRIGS(KC+2)
         IBASE=0
         DO L=1,LA
            I=IBASE
            J=JBASE
CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
            DO IJK=1,LOT
               C(JA+J)=A(IA+I)+(A(IB+I)+A(IC+I))
               D(JA+J)=B(IA+I)+(B(IB+I)-B(IC+I))
               C(JB+J)=
     &              C1*((A(IA+I)-0.5*(A(IB+I)+A(IC+I)))-
     &              (SIN60*(B(IB+I)+B(IC+I))))
     &              -S1*((B(IA+I)-0.5*(B(IB+I)-B(IC+I)))+
     &              (SIN60*(A(IB+I)-A(IC+I))))
               D(JB+J)=
     &              S1*((A(IA+I)-0.5*(A(IB+I)+A(IC+I)))-
     &              (SIN60*(B(IB+I)+B(IC+I))))
     &              +C1*((B(IA+I)-0.5*(B(IB+I)-B(IC+I)))+
     &              (SIN60*(A(IB+I)-A(IC+I))))
               C(JC+J)=
     &              C2*((A(IA+I)-0.5*(A(IB+I)+A(IC+I)))+
     &              (SIN60*(B(IB+I)+B(IC+I))))
     &              -S2*((B(IA+I)-0.5*(B(IB+I)-B(IC+I)))-
     &              (SIN60*(A(IB+I)-A(IC+I))))
               D(JC+J)=
     &              S2*((A(IA+I)-0.5*(A(IB+I)+A(IC+I)))+
     &              (SIN60*(B(IB+I)+B(IC+I))))
     &              +C2*((B(IA+I)-0.5*(B(IB+I)-B(IC+I)))-
     &              (SIN60*(A(IB+I)-A(IC+I))))
               I=I+INC3
               J=J+INC4
            ENDDO
            IBASE=IBASE+INC1
            JBASE=JBASE+INC2
         ENDDO
         IA=IA+IINK
         IB=IB+IINK
         IC=IC-IINK
         JBASE=JBASE+JUMP
      ENDDO
      IF (IA.GT.IC) GO TO 1900

 1360 CONTINUE
      IBASE=0
      DO L=1,LA
         I=IBASE
         J=JBASE
CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
         DO IJK=1,LOT
            C(JA+J)=A(IA+I)+A(IB+I)
            C(JB+J)=(0.5*A(IA+I)-A(IB+I))-(SIN60*B(IA+I))
            C(JC+J)=-(0.5*A(IA+I)-A(IB+I))-(SIN60*B(IA+I))
            I=I+INC3
            J=J+INC4
         ENDDO
         IBASE=IBASE+INC1
         JBASE=JBASE+INC2
      ENDDO
      GO TO 1900

 1390  CONTINUE
      SSIN60=2.0*SIN60
      DO L=1,LA
         I=IBASE
         J=JBASE
CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
         DO IJK=1,LOT
            C(JA+J)=2.0*(A(IA+I)+A(IB+I))
            C(JB+J)=(2.0*A(IA+I)-A(IB+I))-(SSIN60*B(IB+I))
            C(JC+J)=(2.0*A(IA+I)-A(IB+I))+(SSIN60*B(IB+I))
            I=I+INC3
            J=J+INC4
         ENDDO
         IBASE=IBASE+INC1
         JBASE=JBASE+INC2
      ENDDO
      GO TO 1900

C Coding for factor 4:
 1400 CONTINUE
      IA=1
      IB=IA+(2*M-LA)*INC1
      IC=IB+2*M*INC1
      ID=IB
      JA=1
      JB=JA+JINK
      JC=JB+JINK
      JD=JC+JINK

      IF (LA.EQ.M) GO TO 1490

      DO L=1,LA
         I=IBASE
         J=JBASE
CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
         DO IJK=1,LOT
            C(JA+J)=(A(IA+I)+A(IC+I))+A(IB+I)
            C(JB+J)=(A(IA+I)-A(IC+I))-B(IB+I)
            C(JC+J)=(A(IA+I)+A(IC+I))-A(IB+I)
            C(JD+J)=(A(IA+I)-A(IC+I))+B(IB+I)
            I=I+INC3
            J=J+INC4
         ENDDO
         IBASE=IBASE+INC1
         JBASE=JBASE+INC2
      ENDDO
      IA=IA+IINK
      IINK=2*IINK
      IB=IB+IINK
      IC=IC-IINK
      ID=ID-IINK
      JBASE=JBASE+JUMP
      JUMP=2*JUMP+JINK
      IF (IB.EQ.IC) GO TO 1460
      DO K=LA,KSTOP,LA
         KB=K+K
         KC=KB+KB
         KD=KC+KB
         C1=TRIGS(KB+1)
         S1=TRIGS(KB+2)
         C2=TRIGS(KC+1)
         S2=TRIGS(KC+2)
         C3=TRIGS(KD+1)
         S3=TRIGS(KD+2)
         IBASE=0
         DO L=1,LA
            I=IBASE
            J=JBASE
CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
            DO IJK=1,LOT
               C(JA+J)=(A(IA+I)+A(IC+I))+(A(IB+I)+A(ID+I))
               D(JA+J)=(B(IA+I)-B(IC+I))+(B(IB+I)-B(ID+I))
               C(JC+J)=
     &              C2*((A(IA+I)+A(IC+I))-(A(IB+I)+A(ID+I)))
     &              -S2*((B(IA+I)-B(IC+I))-(B(IB+I)-B(ID+I)))
               D(JC+J)=
     &              S2*((A(IA+I)+A(IC+I))-(A(IB+I)+A(ID+I)))
     &              +C2*((B(IA+I)-B(IC+I))-(B(IB+I)-B(ID+I)))
               C(JB+J)=
     &              C1*((A(IA+I)-A(IC+I))-(B(IB+I)+B(ID+I)))
     &              -S1*((B(IA+I)+B(IC+I))+(A(IB+I)-A(ID+I)))
               D(JB+J)=
     &              S1*((A(IA+I)-A(IC+I))-(B(IB+I)+B(ID+I)))
     &              +C1*((B(IA+I)+B(IC+I))+(A(IB+I)-A(ID+I)))
               C(JD+J)=
     &              C3*((A(IA+I)-A(IC+I))+(B(IB+I)+B(ID+I)))
     &              -S3*((B(IA+I)+B(IC+I))-(A(IB+I)-A(ID+I)))
               D(JD+J)=
     &              S3*((A(IA+I)-A(IC+I))+(B(IB+I)+B(ID+I)))
     &              +C3*((B(IA+I)+B(IC+I))-(A(IB+I)-A(ID+I)))
               I=I+INC3
               J=J+INC4
            ENDDO
            IBASE=IBASE+INC1
            JBASE=JBASE+INC2
         ENDDO
         IA=IA+IINK
         IB=IB+IINK
         IC=IC-IINK
         ID=ID-IINK
         JBASE=JBASE+JUMP
      ENDDO
      IF (IB.GT.IC) GO TO 1900
 1460 CONTINUE
      IBASE=0
      SIN45=SQRT(0.5)
      DO L=1,LA
         I=IBASE
         J=JBASE
CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
         DO IJK=1,LOT
            C(JA+J)=A(IA+I)+A(IB+I)
            C(JB+J)=SIN45*((A(IA+I)-A(IB+I))-(B(IA+I)+B(IB+I)))
            C(JC+J)=B(IB+I)-B(IA+I)
            C(JD+J)=-SIN45*((A(IA+I)-A(IB+I))+(B(IA+I)+B(IB+I)))
            I=I+INC3
            J=J+INC4
         ENDDO
         IBASE=IBASE+INC1
         JBASE=JBASE+INC2
      ENDDO
      GO TO 1900

 1490 CONTINUE
      DO L=1,LA
         I=IBASE
         J=JBASE
CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
         DO IJK=1,LOT
            C(JA+J)=2.0*((A(IA+I)+A(IC+I))+A(IB+I))
            C(JB+J)=2.0*((A(IA+I)-A(IC+I))-B(IB+I))
            C(JC+J)=2.0*((A(IA+I)+A(IC+I))-A(IB+I))
            C(JD+J)=2.0*((A(IA+I)-A(IC+I))+B(IB+I))
            I=I+INC3
            J=J+INC4
         ENDDO
         IBASE=IBASE+INC1
         JBASE=JBASE+INC2
      ENDDO
      GO TO 1900

C Coding for factor 5:

 1500 CONTINUE
      IA=1
      IB=IA+(2*M-LA)*INC1
      IC=IB+2*M*INC1
      ID=IC
      IE=IB
      JA=1
      JB=JA+JINK
      JC=JB+JINK
      JD=JC+JINK
      JE=JD+JINK

      IF (LA.EQ.M) GO TO 1590

      DO L=1,LA
         I=IBASE
         J=JBASE
CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
         DO IJK=1,LOT
            C(JA+J)=A(IA+I)+(A(IB+I)+A(IC+I))
            C(JB+J)=((A(IA+I)-0.25*(A(IB+I)+A(IC+I)))+QRT5*
     &           (A(IB+I)-A(IC+I))) - (SIN72*B(IB+I)+SIN36*B(IC+I))
            C(JC+J)=((A(IA+I)-0.25*(A(IB+I)+A(IC+I)))-QRT5*
     &           (A(IB+I)-A(IC+I))) - (SIN36*B(IB+I)-SIN72*B(IC+I))
            C(JD+J)=((A(IA+I)-0.25*(A(IB+I)+A(IC+I)))-QRT5*
     &           (A(IB+I)-A(IC+I))) + (SIN36*B(IB+I)-SIN72*B(IC+I))
            C(JE+J)=((A(IA+I)-0.25*(A(IB+I)+A(IC+I)))+QRT5*
     &           (A(IB+I)-A(IC+I))) + (SIN72*B(IB+I)+SIN36*B(IC+I))
            I=I+INC3
            J=J+INC4
         ENDDO
         IBASE=IBASE+INC1
         JBASE=JBASE+INC2
      ENDDO
      IA=IA+IINK
      IINK=2*IINK
      IB=IB+IINK
      IC=IC+IINK
      ID=ID-IINK
      IE=IE-IINK
      JBASE=JBASE+JUMP
      JUMP=2*JUMP+JINK
      IF (IB.EQ.ID) GO TO 1560
      DO K=LA,KSTOP,LA
         KB=K+K
         KC=KB+KB
         KD=KC+KB
         KE=KD+KB
         C1=TRIGS(KB+1)
         S1=TRIGS(KB+2)
         C2=TRIGS(KC+1)
         S2=TRIGS(KC+2)
         C3=TRIGS(KD+1)
         S3=TRIGS(KD+2)
         C4=TRIGS(KE+1)
         S4=TRIGS(KE+2)
         IBASE=0
         DO L=1,LA
            I=IBASE
            J=JBASE
CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
            DO IJK=1,LOT

               AA10(IJK)=(A(IA+I)-0.25*((A(IB+I)+A(IE+I))+
     &              (A(IC+I)+A(ID+I))))
     &              +QRT5*((A(IB+I)+A(IE+I))-(A(IC+I)+A(ID+I)))
               AA20(IJK)=(A(IA+I)-0.25*((A(IB+I)+A(IE+I))+
     &              (A(IC+I)+A(ID+I))))
     &              -QRT5*((A(IB+I)+A(IE+I))-(A(IC+I)+A(ID+I)))
               BB10(IJK)=(B(IA+I)-0.25*((B(IB+I)-B(IE+I))+
     &              (B(IC+I)-B(ID+I))))
     &              +QRT5*((B(IB+I)-B(IE+I))-(B(IC+I)-B(ID+I)))
               BB20(IJK)=(B(IA+I)-0.25*((B(IB+I)-B(IE+I))+
     &              (B(IC+I)-B(ID+I))))
     &              -QRT5*((B(IB+I)-B(IE+I))-(B(IC+I)-B(ID+I)))
               AA11(IJK)=SIN72*(B(IB+I)+B(IE+I))+
     &              SIN36*(B(IC+I)+B(ID+I))
               AA21(IJK)=SIN36*(B(IB+I)+B(IE+I))-
     &              SIN72*(B(IC+I)+B(ID+I))
               BB11(IJK)=SIN72*(A(IB+I)-A(IE+I))+
     &              SIN36*(A(IC+I)-A(ID+I))
               BB21(IJK)=SIN36*(A(IB+I)-A(IE+I))-
     &              SIN72*(A(IC+I)-A(ID+I))

               C(JA+J)=A(IA+I)+((A(IB+I)+A(IE+I))+(A(IC+I)+A(ID+I)))
               D(JA+J)=B(IA+I)+((B(IB+I)-B(IE+I))+(B(IC+I)-B(ID+I)))
               C(JB+J)=C1*(AA10(IJK)-AA11(IJK))-S1*(BB10(IJK)+BB11(IJK))
               D(JB+J)=S1*(AA10(IJK)-AA11(IJK))+C1*(BB10(IJK)+BB11(IJK))
               C(JE+J)=C4*(AA10(IJK)+AA11(IJK))-S4*(BB10(IJK)-BB11(IJK))
               D(JE+J)=S4*(AA10(IJK)+AA11(IJK))+C4*(BB10(IJK)-BB11(IJK))
               C(JC+J)=C2*(AA20(IJK)-AA21(IJK))-S2*(BB20(IJK)+BB21(IJK))
               D(JC+J)=S2*(AA20(IJK)-AA21(IJK))+C2*(BB20(IJK)+BB21(IJK))
               C(JD+J)=C3*(AA20(IJK)+AA21(IJK))-S3*(BB20(IJK)-BB21(IJK))
               D(JD+J)=S3*(AA20(IJK)+AA21(IJK))+C3*(BB20(IJK)-BB21(IJK))

               I=I+INC3
               J=J+INC4
            ENDDO
            IBASE=IBASE+INC1
            JBASE=JBASE+INC2
         ENDDO
         IA=IA+IINK
         IB=IB+IINK
         IC=IC+IINK
         ID=ID-IINK
         IE=IE-IINK
         JBASE=JBASE+JUMP
      ENDDO
      IF (IB.GT.ID) GO TO 1900
 1560  CONTINUE
      IBASE=0
      DO L=1,LA
         I=IBASE
         J=JBASE
CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
         DO IJK=1,LOT
            C(JA+J)=(A(IA+I)+A(IB+I))+A(IC+I)
            C(JB+J)=(QRT5*(A(IA+I)-A(IB+I))+
     &           (0.25*(A(IA+I)+A(IB+I))-A(IC+I)))
     &           -(SIN36*B(IA+I)+SIN72*B(IB+I))
            C(JE+J)=-(QRT5*(A(IA+I)-A(IB+I))+
     &           (0.25*(A(IA+I)+A(IB+I))-A(IC+I)))
     &           -(SIN36*B(IA+I)+SIN72*B(IB+I))
            C(JC+J)=(QRT5*(A(IA+I)-A(IB+I))-
     &           (0.25*(A(IA+I)+A(IB+I))-A(IC+I)))
     &           -(SIN72*B(IA+I)-SIN36*B(IB+I))
            C(JD+J)=-(QRT5*(A(IA+I)-A(IB+I))-
     &           (0.25*(A(IA+I)+A(IB+I))-A(IC+I)))
     &           -(SIN72*B(IA+I)-SIN36*B(IB+I))
            I=I+INC3
            J=J+INC4
         ENDDO
         IBASE=IBASE+INC1
         JBASE=JBASE+INC2
      ENDDO
      GO TO 1900

 1590 CONTINUE
      QQRT5=2.0*QRT5
      SSIN36=2.0*SIN36
      SSIN72=2.0*SIN72
      DO L=1,LA
         I=IBASE
         J=JBASE
CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
         DO IJK=1,LOT
            C(JA+J)=2.0*(A(IA+I)+(A(IB+I)+A(IC+I)))
            C(JB+J)=(2.0*(A(IA+I)-0.25*(A(IB+I)+A(IC+I)))
     &           +QQRT5*(A(IB+I)-A(IC+I)))-(SSIN72*B(IB+I)+
     &           SSIN36*B(IC+I))
            C(JC+J)=(2.0*(A(IA+I)-0.25*(A(IB+I)+A(IC+I)))
     &           -QQRT5*(A(IB+I)-A(IC+I)))-(SSIN36*B(IB+I)-
     &           SSIN72*B(IC+I))
            C(JD+J)=(2.0*(A(IA+I)-0.25*(A(IB+I)+A(IC+I)))
     &           -QQRT5*(A(IB+I)-A(IC+I)))+(SSIN36*B(IB+I)-
     &           SSIN72*B(IC+I))
            C(JE+J)=(2.0*(A(IA+I)-0.25*(A(IB+I)+A(IC+I)))
     &           +QQRT5*(A(IB+I)-A(IC+I)))+(SSIN72*B(IB+I)+
     &           SSIN36*B(IC+I))
            I=I+INC3
            J=J+INC4
         ENDDO
         IBASE=IBASE+INC1
         JBASE=JBASE+INC2
      ENDDO
      GO TO 1900

C Coding for factor 6:
 1600 CONTINUE
      IA=1
      IB=IA+(2*M-LA)*INC1
      IC=IB+2*M*INC1
      ID=IC+2*M*INC1
      IE=IC
      IF=IB
      JA=1
      JB=JA+JINK
      JC=JB+JINK
      JD=JC+JINK
      JE=JD+JINK
      JF=JE+JINK

      IF (LA.EQ.M) GO TO 1690

      DO L=1,LA
         I=IBASE
         J=JBASE
CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
         DO IJK=1,LOT
            C(JA+J)=(A(IA+I)+A(ID+I))+(A(IB+I)+A(IC+I))
            C(JD+J)=(A(IA+I)-A(ID+I))-(A(IB+I)-A(IC+I))
            C(JB+J)=((A(IA+I)-A(ID+I))+0.5*(A(IB+I)-A(IC+I)))
     &           -(SIN60*(B(IB+I)+B(IC+I)))
            C(JF+J)=((A(IA+I)-A(ID+I))+0.5*(A(IB+I)-A(IC+I)))
     &           +(SIN60*(B(IB+I)+B(IC+I)))
            C(JC+J)=((A(IA+I)+A(ID+I))-0.5*(A(IB+I)+A(IC+I)))
     &           -(SIN60*(B(IB+I)-B(IC+I)))
            C(JE+J)=((A(IA+I)+A(ID+I))-0.5*(A(IB+I)+A(IC+I)))
     &           +(SIN60*(B(IB+I)-B(IC+I)))
            I=I+INC3
            J=J+INC4
         ENDDO
         IBASE=IBASE+INC1
         JBASE=JBASE+INC2
      ENDDO
      IA=IA+IINK
      IINK=2*IINK
      IB=IB+IINK
      IC=IC+IINK
      ID=ID-IINK
      IE=IE-IINK
      IF=IF-IINK
      JBASE=JBASE+JUMP
      JUMP=2*JUMP+JINK
      IF (IC.EQ.ID) GO TO 1660
      DO K=LA,KSTOP,LA
         KB=K+K
         KC=KB+KB
         KD=KC+KB
         KE=KD+KB
         KF=KE+KB
         C1=TRIGS(KB+1)
         S1=TRIGS(KB+2)
         C2=TRIGS(KC+1)
         S2=TRIGS(KC+2)
         C3=TRIGS(KD+1)
         S3=TRIGS(KD+2)
         C4=TRIGS(KE+1)
         S4=TRIGS(KE+2)
         C5=TRIGS(KF+1)
         S5=TRIGS(KF+2)
         IBASE=0
         DO L=1,LA
            I=IBASE
            J=JBASE
CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
            DO IJK=1,LOT

               AA11(IJK)= (A(IE+I)+A(IB+I))+(A(IC+I)+A(IF+I))
               AA20(IJK)=(A(IA+I)+A(ID+I))-0.5*AA11(IJK)
               AA21(IJK)=SIN60*((A(IE+I)+A(IB+I))-(A(IC+I)+A(IF+I)))
               BB11(IJK)= (B(IB+I)-B(IE+I))+(B(IC+I)-B(IF+I))
               BB20(IJK)=(B(IA+I)-B(ID+I))-0.5*BB11(IJK)
               BB21(IJK)=SIN60*((B(IB+I)-B(IE+I))-(B(IC+I)-B(IF+I)))

               C(JA+J)=(A(IA+I)+A(ID+I))+AA11(IJK)
               D(JA+J)=(B(IA+I)-B(ID+I))+BB11(IJK)
               C(JC+J)=C2*(AA20(IJK)-BB21(IJK))-S2*(BB20(IJK)+AA21(IJK))
               D(JC+J)=S2*(AA20(IJK)-BB21(IJK))+C2*(BB20(IJK)+AA21(IJK))
               C(JE+J)=C4*(AA20(IJK)+BB21(IJK))-S4*(BB20(IJK)-AA21(IJK))
               D(JE+J)=S4*(AA20(IJK)+BB21(IJK))+C4*(BB20(IJK)-AA21(IJK))

               AA11(IJK)=(A(IE+I)-A(IB+I))+(A(IC+I)-A(IF+I))
               BB11(IJK)=(B(IE+I)+B(IB+I))-(B(IC+I)+B(IF+I))
               AA20(IJK)=(A(IA+I)-A(ID+I))-0.5*AA11(IJK)
               AA21(IJK)=SIN60*((A(IE+I)-A(IB+I))-(A(IC+I)-A(IF+I)))
               BB20(IJK)=(B(IA+I)+B(ID+I))+0.5*BB11(IJK)
               BB21(IJK)=SIN60*((B(IE+I)+B(IB+I))+(B(IC+I)+B(IF+I)))

               C(JD+J)=C3*((A(IA+I)-A(ID+I))+AA11(IJK))-
     &              S3*((B(IA+I)+B(ID+I))-BB11(IJK))
               D(JD+J)=S3*((A(IA+I)-A(ID+I))+AA11(IJK))+
     &              C3*((B(IA+I)+B(ID+I))-BB11(IJK))
               C(JB+J)=C1*(AA20(IJK)-BB21(IJK))-S1*(BB20(IJK)-AA21(IJK))
               D(JB+J)=S1*(AA20(IJK)-BB21(IJK))+C1*(BB20(IJK)-AA21(IJK))
               C(JF+J)=C5*(AA20(IJK)+BB21(IJK))-S5*(BB20(IJK)+AA21(IJK))
               D(JF+J)=S5*(AA20(IJK)+BB21(IJK))+C5*(BB20(IJK)+AA21(IJK))

               I=I+INC3
               J=J+INC4
            ENDDO
            IBASE=IBASE+INC1
            JBASE=JBASE+INC2
         ENDDO
         IA=IA+IINK
         IB=IB+IINK
         IC=IC+IINK
         ID=ID-IINK
         IE=IE-IINK
         IF=IF-IINK
         JBASE=JBASE+JUMP
      ENDDO
      IF (IC.GT.ID) GO TO 1900
 1660 CONTINUE
      IBASE=0
      DO L=1,LA
         I=IBASE
         J=JBASE
CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
         DO IJK=1,LOT
            C(JA+J)=A(IB+I)+(A(IA+I)+A(IC+I))
            C(JD+J)=B(IB+I)-(B(IA+I)+B(IC+I))
            C(JB+J)=(SIN60*(A(IA+I)-A(IC+I)))-
     &           (0.5*(B(IA+I)+B(IC+I))+B(IB+I))
            C(JF+J)=-(SIN60*(A(IA+I)-A(IC+I)))-
     &           (0.5*(B(IA+I)+B(IC+I))+B(IB+I))
            C(JC+J)=SIN60*(B(IC+I)-B(IA+I))+
     &           (0.5*(A(IA+I)+A(IC+I))-A(IB+I))
            C(JE+J)=SIN60*(B(IC+I)-B(IA+I))-
     &           (0.5*(A(IA+I)+A(IC+I))-A(IB+I))
            I=I+INC3
            J=J+INC4
         ENDDO
         IBASE=IBASE+INC1
         JBASE=JBASE+INC2
      ENDDO
      GO TO 1900

 1690 CONTINUE
      SSIN60=2.0*SIN60
      DO L=1,LA
         I=IBASE
         J=JBASE
CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
         DO IJK=1,LOT
            C(JA+J)=(2.0*(A(IA+I)+A(ID+I)))+(2.0*(A(IB+I)+A(IC+I)))
            C(JD+J)=(2.0*(A(IA+I)-A(ID+I)))-(2.0*(A(IB+I)-A(IC+I)))
            C(JB+J)=(2.0*(A(IA+I)-A(ID+I))+(A(IB+I)-A(IC+I)))
     &           -(SSIN60*(B(IB+I)+B(IC+I)))
            C(JF+J)=(2.0*(A(IA+I)-A(ID+I))+(A(IB+I)-A(IC+I)))
     &           +(SSIN60*(B(IB+I)+B(IC+I)))
            C(JC+J)=(2.0*(A(IA+I)+A(ID+I))-(A(IB+I)+A(IC+I)))
     &           -(SSIN60*(B(IB+I)-B(IC+I)))
            C(JE+J)=(2.0*(A(IA+I)+A(ID+I))-(A(IB+I)+A(IC+I)))
     &           +(SSIN60*(B(IB+I)-B(IC+I)))
            I=I+INC3
            J=J+INC4
         ENDDO
         IBASE=IBASE+INC1
         JBASE=JBASE+INC2
      ENDDO
      GO TO 1900

C Coding for factor 8:
 1800 CONTINUE
      IBAD=3
      IF (LA.NE.M) GO TO 1910
      IA=1
      IB=IA+LA*INC1
      IC=IB+2*LA*INC1
      ID=IC+2*LA*INC1
      IE=ID+2*LA*INC1
      JA=1
      JB=JA+JINK
      JC=JB+JINK
      JD=JC+JINK
      JE=JD+JINK
      JF=JE+JINK
      JG=JF+JINK
      JH=JG+JINK
      SSIN45=SQRT(2.0)

      DO L=1,LA
         I=IBASE
         J=JBASE
CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
         DO IJK=1,LOT
            C(JA+J)=2.0*(((A(IA+I)+A(IE+I))+A(IC+I))+
     &           (A(IB+I)+A(ID+I)))
            C(JE+J)=2.0*(((A(IA+I)+A(IE+I))+A(IC+I))-
     &           (A(IB+I)+A(ID+I)))
            C(JC+J)=2.0*(((A(IA+I)+A(IE+I))-A(IC+I))-
     &           (B(IB+I)-B(ID+I)))
            C(JG+J)=2.0*(((A(IA+I)+A(IE+I))-A(IC+I))+
     &           (B(IB+I)-B(ID+I)))
            C(JB+J)=2.0*((A(IA+I)-A(IE+I))-B(IC+I))
     &           +SSIN45*((A(IB+I)-A(ID+I))-(B(IB+I)+B(ID+I)))
            C(JF+J)=2.0*((A(IA+I)-A(IE+I))-B(IC+I))
     &           -SSIN45*((A(IB+I)-A(ID+I))-(B(IB+I)+B(ID+I)))
            C(JD+J)=2.0*((A(IA+I)-A(IE+I))+B(IC+I))
     &           -SSIN45*((A(IB+I)-A(ID+I))+(B(IB+I)+B(ID+I)))
            C(JH+J)=2.0*((A(IA+I)-A(IE+I))+B(IC+I))
     &           +SSIN45*((A(IB+I)-A(ID+I))+(B(IB+I)+B(ID+I)))
            I=I+INC3
            J=J+INC4
         ENDDO
         IBASE=IBASE+INC1
         JBASE=JBASE+INC2
      ENDDO

 1900 CONTINUE
      IBAD=0
 1910 CONTINUE
      IERR=IBAD

      ENDIF    !  end of Fourier Synthesis

      RETURN   ! end of FTRANS
      END

!- End of subroutine code-----------------------------------------
