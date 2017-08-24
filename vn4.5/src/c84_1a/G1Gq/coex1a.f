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
C Modelling at the above address.
C ******************************COPYRIGHT******************************
C
CLL   SUBROUTINE COEX,COEX2,CMPS,XPND,INSTIN,EXTRIN -----------------
CLL
CLL   PURPOSE:   PACK TO AND UNPACK FROM WGDOS FORMAT
CLL
CLL   (Note that for optimal performance the routines INSTIN and
CLL    EXTRIN are inline expanded (enable option 8 on fpp)
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  Date
CLL   3.2  19/04/93  Code for ne real missing data indicator (TCJ).
CLL   4.2  Nov. 96   T3E migration: WHENNE,WHENEQ replaced by 
CLL                   portable fortran code.   S.J.Swarbrick
CLL   4.3  Apr. 97 T3E migration: Calling the functions CRI2IBM/IBM2CRI
CLL                directly replaces the calls to the (now unsupported)
CLL                routines USICTI, USICTC, USSCTC, USSCTI. Ian Edmond
CLL
CLL   PROGRAMMING STANDARD: UNIFIED MODEL DOCUMENTATION PAPER NO. 4,
CLL   STANDARD B, VERSION 2, DATED 18/01/90
CLL
CLL  Logical component number: S72
CLL
CLL   SYSTEM TASK: P7
CLL   OWNED BY P J SMITH
CLL
CLL   DOCUMENTATION:  ??????
CLLEND-------------------------------------------------------------

      SUBROUTINE COEX(FIELD,M,ICOMP,N,IX,IY,NUM,ISC,OCO,RMDI,LWORD)
      INTEGER N,ICOMP(N),M,IX,IY,NUM,ISC,LWORD
      REAL FIELD(M),RMDI
      LOGICAL OCO
      INTEGER JJ,IFAIL
C     INTEGER IRC(IY)    ! this was cause of errors as IY not set
      EXTERNAL COEX2
C
C     CRAY VERSION C.1.1  16/11/90  P J SMITH
C
C
C                     OCO=.TRUE.                 OCO=.FALSE.
C
C      FIELD   =  FIELD FOR COMPRESSING     RETURNED EXPANDED DATA
C          M   =  SIZE OF FIELD             SIZE OF FIELD
C      ICOMP   =  RETURNED COMPRESSED DATA  FIELD FOR EXPANSION
C          N   =  SIZE OF COMP                 -------
C         IX   =  X DIMENSION OF FIELD      RETURNED X DIMENSION
C         IY   =  Y DIMENSION OF FIELD      RETURNED Y DIMENSION
C        NUM   =  TOTAL NO. OF COMPRESSED      -------
C                 (32 BIT) WORDS RETURNED
C        ISC   =  ACCURACY IN POWER OF 2
C        OCO   =  .TRUE. FOR COMPRESSION    .FALSE. FOR EXPANSION
C
C        USERS OF THIS ROUTINE SHOULD ENSURE THAT THE ARRAY 'COMP'
C      WILL BE BIG ENOUGH TO HOLD THE COMPRESSED DATA.
C        IF INPUT ARRAY IS ONE DIMENSIONAL PUT IY=1 AND IX=DIMENSION,
C      WHEN USING WITH PRINTFILE FIELDS USE IX=192,IY=121  NOT IX=23232,
C      IY=1 THIS WILL MAKE THE COMPRESSION MORE EFFICIENT.
C      FOR MOST PRINTFILE FIELDS USING AN ACCURACY IF AN EIGTH (ISC=-3)
C      A COMPRESSION FACTOR OF 4 CAN BE ACHIEVED. FOR DDDFFF FIELDS
C      USERS ARE ADVISED TO SPLIT THE FIELD INTO U AND V
C      COMPONENTS.
C
C     CRAY ROUTINE     64 BIT WORDS - CRAY FULLWORD
C
      IF (OCO) THEN
C     WRITE(6,*) 'CRAY COEX - PACKING'
C
C         COMPRESSION OF DATA
C
C
C         INSERT SCALE FACTOR, COLS AND ROWS INTO HEADER
C
          IER=CRI2IBM(2,1,ICOMP(1),32,ISC,1,64,32)
          IER=CRI2IBM(2,1,ICOMP(2),0,IX,1,64,16)
          IER=CRI2IBM(2,1,ICOMP(2),16,IY,1,64,16)
C
C         CALL AUTOTASKED/DYNAMIC ARRAY PART OF COEX
C
          CALL COEX2(FIELD,M,ICOMP,N,IX,IY,NUM,ISC,OCO,RMDI)
C         CALL COEX2(FIELD,M,ICOMP,N,IX,IY,NUM,ISC,OCO,RMDI,IRC)
C
C         CHECK RETURN CODES
C
C         IFAIL=0
C         DO JJ = 1,IY
C             IF (IRC(JJ).NE.0) IFAIL=IFAIL+1
C         ENDDO
C         IF(IFAIL.GT.0)THEN
C         DO JJ = 1,IY
C             IF (IRC(JJ).NE.0) THEN
C              WRITE(6,*)'RETURN CODE',IRC(JJ),'FROM COEX FOR ROW',JJ
C             ENDIF
C         ENDDO
C         ENDIF
C
C         END OF COMPRESSION SECTION
C
      ELSE
C     WRITE(6,*) 'CRAY COEX - UNPACKING'
C
C         EXPANSION SECTION
C
C
C         EXTRACT SCALE FACTOR, COLS AND ROWS FROM HEADER
C
          IER=IBM2CRI(2,1,ICOMP(1),32,ISC,1,64,32)
          IER=IBM2CRI(2,1,ICOMP(2),0,IX,1,64,16)
          IER=IBM2CRI(2,1,ICOMP(2),16,IY,1,64,16)
C
C         CALL AUTOTASKED/DYNAMIC ARRAY PART OF COEX
C
          CALL COEX2(FIELD,M,ICOMP,N,IX,IY,NUM,ISC,OCO,RMDI)
C         CALL COEX2(FIELD,M,ICOMP,N,IX,IY,NUM,ISC,OCO,RMDI,IRC)
C
C         CHECK RETURN CODES (NOT YET IMPLEMENTED)
C
C         DO JJ = 1,IY
C             IF (IRC(JJ).NE.0) THEN
C                 WRITE(6,*)' NON-ZERO RETURN CODE FOR ROW ',JJ,IRC(JJ)
C             ENDIF
C         ENDDO
      END IF
C
C
      RETURN

      END
      SUBROUTINE COEX2(FIELD,M,ICOMP,N,IX,IY,NUM,ISC,OCO,RMDI)
C     SUBROUTINE COEX2(FIELD,M,ICOMP,N,IX,IY,NUM,ISC,OCO,RMDI,IRC)
      INTEGER N,ICOMP(N),M,IX,IY,NUM,ISC
C     INTEGER IRC(IY)
      REAL FIELD(M),RMDI
      LOGICAL OCO
      INTEGER JJ,KK,ILEN,IST,ICX,JCX,JBIT,JLEN,NOB,NOC,IER2
      INTEGER IC(IY),ICB(IY),NOP(IY),JC(IY),IBIT(IY),ISTART(IY)
      INTEGER JCOMP(IX,IY),IERR1(IY),IERR2(IY),IERR3(IY)
      REAL ACC,APREC
      REAL BASE(IY)
      EXTERNAL CMPS,XPND,STRMOV
C
C     CRAY VERSION C.1.1  16/11/90  P J SMITH
C     CRAY ROUTINE     64 BIT WORDS - CRAY FULLWORD
C
C                     OCO=.TRUE.                 OCO=.FALSE.
C
C      FIELD   =  FIELD FOR COMPRESSING     RETURNED EXPANDED DATA
C          M   =  SIZE OF FIELD             SIZE OF FIELD
C      ICOMP   =  RETURNED COMPRESSED DATA  FIELD FOR EXPANSION
C          N   =  SIZE OF COMP                 -------
C         IX   =  X DIMENSION OF FIELD      X DIMENSION OF FIELD
C         IY   =  Y DIMENSION OF FIELD      Y DIMENSION OF FIELD
C        NUM   =  TOTAL NO. OF COMPRESSED      -------
C                 (32 BIT) WORDS RETURNED
C        ISC   =  ACCURACY IN POWER OF 2    ACCURACY IN POWER OF 2
C        OCO   =  .TRUE. FOR COMPRESSION    .FALSE. FOR EXPANSION
C        IRC   =  RETURN CODE FOR EACH ROW  RETURN CODE FOR EACH ROW
C
C     INITIALISE TEMP. VARIABLES/ARRAYS
C
      IF (.NOT. OCO) THEN
        DO JJ=1,IY
          NOP(JJ)=0
          IBIT(JJ)=0
        ENDDO
        DO JJ=1,IY
          BASE(JJ)=0.0
          IBIT(JJ)=0.0
        ENDDO
      ENDIF

      DO JJ=1,IY
          DO JCX=1,IX
             JCOMP(JCX,JJ) = 0
          END DO
      END DO
      IF (OCO) THEN
C
C         COMPRESSION OF DATA
C
          IC(1) = 1
          ACC   = ISC
          APREC = 2.**ACC
C
C         PUT PACKED DATA FOR EACH ROW INTO TEMP ARRAY JCOMP
C
CFPP$ CNCALL
          DO 10 JJ = 1,IY
C             IP      = POSITION IN INPUT ARRAY
              IP      = (JJ-1)*IX + 1
              CALL CMPS(IX,FIELD(IP),JCOMP(2,JJ),NOP(JJ),APREC,
     +                                          IBIT(JJ),BASE(JJ),RMDI)
   10     CONTINUE
C
C         ADD BASE VALUE, NUMBER OF BITS, AND LENGTH
C
CFPP$ CNCALL
          DO 11 JJ = 1,IY
              IERR1(JJ)=CRI2IBM(3,1, JCOMP(1,JJ),0, BASE(JJ),1,64,32)
              IERR2(JJ)=CRI2IBM(2,1,JCOMP(1,JJ),32,IBIT(JJ),1,64,16)
              IERR3(JJ)=CRI2IBM(2,1,JCOMP(1,JJ),48,NOP(JJ),1,64,16)
   11     CONTINUE
C
C         CHECK ROW HEADER AND SET RETURN CODES
C
CFPP$ NOCONCUR
C         DO 12 JJ = 1,IY
C             IF (IERR1(JJ).NE.0) IRC(JJ) = IRC(JJ) + 1
C             IF (IERR2(JJ).NE.0) IRC(JJ) = IRC(JJ) + 2
C             IF (IERR3(JJ).NE.0) IRC(JJ) = IRC(JJ) + 4
C             IF (JCOMP(1,JJ).EQ.0) THEN
C                 IF (IBIT(JJ).NE.0) IRC(JJ) = IRC(JJ) + 8
C                 IF (NOP(JJ) .NE.0) IRC(JJ) = IRC(JJ) + 16
C             ENDIF
C  12     CONTINUE
C
C         CALCULATE POSITIONS IN OUTPUT ARRAY FOR PACKED ROWS
C         (FIRST PACKED ROW STARTS AT WORD 1; BIT 31)
C
          IC(1)     = 2
          ICB(1)    = -1
          ISTART(1) = 5
CFPP$ NOCONCUR
          DO 20 JJ = 2,IY
              IF (MOD(NOP(JJ-1),2).EQ.1) THEN
                  IC(JJ ) = IC(JJ-1) + NOP(JJ-1)/2 + 1
                  ICB(JJ) = -ICB(JJ-1)
                  IF (ICB(JJ).GT.0) IC(JJ) = IC(JJ) + 1
              ELSE
                  IC(JJ)  = IC(JJ-1) + (NOP(JJ-1)+1)/2 + 1
                  ICB(JJ) = ICB(JJ-1)
              ENDIF
              ISTART(JJ)  = 5
              IF(ICB(JJ).EQ.1) ISTART(JJ) = 1
   20     CONTINUE
C
C         MOVE TEMP. ARRAY INTO OUTPUT ARRAY
C
CFPP$ NOCONCUR
          DO 30 JJ=1,IY
              NOB  = NOP(JJ)*4 + 8
C CHECK IF PACKED FIELD GREATER THAN UN PACKED FIELD
C             IF(NOB.GT.IX*8)IRC(JJ)=IRC(JJ)+32
              IST  = ISTART(JJ)
              ICX  = IC(JJ)
              CALL STRMOV(JCOMP(1,JJ),1,NOB,ICOMP(ICX),IST)
   30     CONTINUE
C
C         INSERT TOTAL LENGTH OF THIS FIELD
          NUM = IC(IY)*2 + NOP(IY)
          IF (ICB(IY).LT.0) NUM = NUM + 1
          IER2=CRI2IBM(2,1,ICOMP(1),0,NUM,1,64,32)
C
C         END OF COMPRESSION SECTION
C
      ELSE
C
C         EXPANSION SECTION
C
          ACC   = ISC
          APREC = 2.**ACC
          ICX   = 2
          JCX   = -1
CFPP$ NOCONCUR
          DO 40 JJ = 1,IY
C
C             MOVE PACKED ROWS INTO TEMP ARRAYS
C
              IF (JCX.LT.0) THEN
C
C                 EXTRACT BASE, NO. BITS, NO 32 BIT WORDS
C
                  IERR=IBM2CRI(3,1,ICOMP(ICX),32,BASE(JJ),1,64,32)
                  IERR=IBM2CRI(2,1,ICOMP(ICX+1),0,IBIT(JJ),1,64,16)
                  IERR=IBM2CRI(2,1,ICOMP(ICX+1),16,NOP(JJ),1,64,16)
C                 SAVE START POSITION OF ROW
                  IC(JJ)     = ICX
                  ISTART(JJ) = 5
              ELSE
                  IERR=IBM2CRI(3,1,ICOMP(ICX),0,BASE(JJ),1,64,32)
                  IERR=IBM2CRI(2,1,ICOMP(ICX),32,IBIT(JJ),1,64,16)
                  IERR=IBM2CRI(2,1,ICOMP(ICX),48,NOP(JJ),1,64,16)
C                 SAVE START POSITION OF ROW
                  IC(JJ)     = ICX
                  ISTART(JJ) = 1
              END IF
C
C             CALCULATE START OF NEXT ROW
C
              IF (MOD(NOP(JJ),2).EQ.1) THEN
                  ICX   = ICX + NOP(JJ)/2 + 1
                  JCX   = -JCX
                  IF (JCX.GT.0) ICX = ICX + 1
              ELSE
                  ICX   = ICX + (NOP(JJ)+1)/2 + 1
              END IF
   40     CONTINUE
C
C         MOVE EACH PACKED ROW INTO TEMP ARRAY JCOMP
C
CFPP$ NOCONCUR
          DO 50 JJ = 1,IY
              ICX  = IC(JJ)
              IST  = ISTART(JJ)
              NOB  = NOP(JJ)*4 + 8
              CALL STRMOV(ICOMP(ICX),IST,NOB,JCOMP(1,JJ),1)
   50     CONTINUE
C
C         CALCULATE START OF EACH ROW IN FIELD
C
          ICX = 1
          DO 60 JJ = 1,IY
              IC(JJ) = ICX
              ICX    = ICX + IX
   60     CONTINUE
C
C         UNPACK DATA INTO FIELD
C
CFPP$ CNCALL
          DO 70 JJ=1,IY
              ICX    = IC(JJ)
              CALL XPND(IX,JCOMP(2,JJ),FIELD(ICX),APREC,IBIT(JJ),
     +                                            BASE(JJ),NOP(JJ),RMDI)
   70     CONTINUE
      END IF
      RETURN

      END
      SUBROUTINE CMPS(IX,FIELD,ICOMP,NOP,APREC,IBIT,BASE,RMDI)
CFPP$ NOCONCUR R
CFPP$ EXPAND (INSTIN)
      INTEGER IX,NOP,ICOMP(IX-1),IBIT
      REAL FIELD(IX),APREC,BASE,RMDI
      INTEGER I,IBASE,II,IMAX,I2,JBIT,JJ,JWORD,NMIS,NPOINT,NZERO
      INTEGER IMAP(IX),IMIS(IX),IZERO(IX),ITEMP(IX)
      INTEGER IGATH1(IX),IGATH2(IX)
      REAL ATEMP(IX),BPREC
      LOGICAL OBTMIS,OBTZER
C
C     CRAY VERSION C.1.1  16/11/90  P J SMITH
C
C         IX   =  LENTH OF FIELD
C      FIELD   =  FIELD FOR COMPRESSING
C      ICOMP   =  RETURNED COMPRESSED DATA
C        NOP   =  NUMBER OF WORDS OF COMP FILLED BY THIS CALL
C      APREC   =  PRECISION
C       IBIT   =  NUMBER OF BITS INTO WHICH DATA PACKED
C       BASE   =  REFERENCE (MIN) VALUE FOR PACKED DATA
C       RMDI   =  MISSING DATA INDICATOR
C
C     INITIALISE VARIABLES/TEMP ARRAYS
C
      OBTMIS = .FALSE.
      OBTZER = .FALSE.
      BPREC  = 1./APREC
      DO I = 1,IX
          IMAP(I)  = 1
      END DO
      DO I = 1,IX-1
          ICOMP(I) = 0
      END DO
C
C     SCAN FIELD FOR MISSING DATA AND STORE RESULT IN IMIS,
C     SCAN FIELD FOR ZERO VALUES AND STORE RESULT IN IZERO,
C     SCAN FIELD FOR MINIMUM VALUE (IGNORE MISSING DATA INDICATOR)
C
      BASE  = 999999.
      NMIS  = 0
      NZERO = 0
      JJ    = 0

      DO I = 1,IX
       IF (FIELD(I).NE.RMDI) THEN
         IMIS(I)=0
       ELSE
         IMIS(I)=1
       ENDIF
      END DO

      DO I = 1,IX
       IF (FIELD(I).NE.0.0) THEN
         IZERO(I)=0
       ELSE
         IZERO(I)=1
       ENDIF
      END DO

C
C GET NO. OF NON-RMDI POINTS + COMPRESS INDEX TO REMOVE THEM
C
       JJ  =0
       DO I=1,IX
         IF (FIELD(I).NE.RMDI) THEN
           JJ        =JJ+1
           IGATH1(JJ)=I
         END IF
       END DO
C
      NMIS=IX-JJ
C
      IF(JJ.NE.0)THEN
C REMOVE MISSING DATA
CDIR$ IVDEP
       DO I =1,JJ
       ATEMP(I)=FIELD(IGATH1(I))
       END DO
C
C GET NO. OF NON-ZERO (NON-RMDI) POINTS + COMPRESS INDEX TO REMOVE THEM
C
       NPOINT=0
       DO I  =1,JJ
         IF (ATEMP(I).NE.0.0) THEN
           NPOINT        =NPOINT+1
           IGATH2(NPOINT)=I
         END IF
       END DO
C
       NZERO=JJ-NPOINT
C
C SET BASE VALUE
       DO I =1,JJ
       IF(ATEMP(I).LT.BASE) BASE=ATEMP(I)
       END DO
      ENDIF
C
C     CHECK IF A BITMAP FOR MISSING DATA IS NEEDED,
C
      IF (NMIS.GT.0) THEN
          OBTMIS = .TRUE.
      ELSE
          OBTMIS = .FALSE.
      END IF
C
C     ROUND BASE TO PRECISION REQUIRED
C
      IF (BASE.NE.0.0) THEN
          BASE  = BASE*BPREC
          IBASE = NINT(BASE)
          BASE  = IBASE*APREC
      ELSE
          IBASE=0
      END IF
! Fujitsu vectorization directive
!OCL NOVREC
C
C     FIND DIFFERENCE FROM BASE AND SCALE
C     FIND MAXIMUM DIFFERENCE
C
      IMAX = 0
      DO 20 I = 1,JJ
          ATEMP(I) = ATEMP(I)*BPREC
          ITEMP(I) = NINT(ATEMP(I)) - IBASE
          IF(ITEMP(I).LT.0) ITEMP(I) = 0
          IF (IMAX.LT.ITEMP(I)) IMAX = ITEMP(I)
   20 CONTINUE
C
C     FIND NUMBER OF BITS REQUIRED TO STORE MAX DIFFERENCE
C
      IBIT  = 0
      IF (IMAX.GT.0) THEN
          I2    = 1
          DO WHILE(IMAX.GE.I2)
              IBIT  = IBIT + 1
              I2    = I2*2
          ENDDO
      ENDIF
C
C     SET START POSITION IN OUTPUT ARRAY
C
      JWORD = 1
      JBIT  = 63
C
C     IF BIT-MAPPING FOR MISSING DATA THEN INSERT BITMAP
C
      IF (OBTMIS) THEN
          DO 30 I = 1,IX
              IF (IMIS(I).EQ.1) THEN
                  ICOMP(JWORD) = IBSET(ICOMP(JWORD),JBIT)
                  IMAP(I)      = 0
              END IF
              IF (JBIT.EQ.0) THEN
                  JBIT  = 63
                  JWORD = JWORD + 1
              ELSE
                  JBIT  = JBIT - 1
              END IF
   30     CONTINUE
      END IF
C
C     IF WORTHWHILE USE BIT MAP AND COMPRESS OUT ZEROS.
C
      IF (IBIT.GT.0) THEN
          IF (NZERO.GT.IX/IBIT) THEN
              OBTZER = .TRUE.
              DO 40 I = 1,IX
                  IF (IZERO(I).EQ.1) THEN
                      ICOMP(JWORD) = IBCLR(ICOMP(JWORD),JBIT)
                      IMAP(I)      = 0
                  ELSE
                      ICOMP(JWORD) = IBSET(ICOMP(JWORD),JBIT)
                  END IF
                  IF (JBIT.EQ.0) THEN
                      JBIT  = 63
                      JWORD = JWORD + 1
                  ELSE
                      JBIT  = JBIT - 1
                  END IF
   40         CONTINUE
          ELSE
              OBTZER = .FALSE.
          END IF
C
C         IF BIT MAP INCLUDED FILL TO END OF CURRENT 32 BIT WORD
C         AND SET POINTERS TO START OF NEXT 32 BIT BOUNDARY.
C
          IF (OBTZER .OR. OBTMIS) THEN
              DO 50 I = 0,JBIT
                  ICOMP(JWORD) = IBSET(ICOMP(JWORD),I)
   50         CONTINUE
              IF (JBIT.NE.63) THEN
                  IF (JBIT.GE.31) THEN
                      JBIT = 31
                  ELSE
                      JWORD = JWORD + 1
                      JBIT = 63
                  ENDIF
              ENDIF
          ELSE
              JBIT = 63
          END IF
C
C         IF BIT MAPPING ZEROS - COMPRESS OUT UNWANTED ZERO VALUES
C        (OTHERWISE PACK ALL NON RMDI DATA (JJ) )
C
          IF (OBTZER) THEN
CDIR$ IVDEP
           DO I= 1,NPOINT
           ITEMP(I)=ITEMP(IGATH2(I))
           END DO
          ELSE
           NPOINT=JJ
          END IF
C
C         MOVE INTO OUTPUT ARRAY USING MINIMUM NUMBER OF BITS REQUIRED
C
          DO 70 I = 1,NPOINT
              CALL INSTIN(ICOMP(JWORD),2,JBIT,IBIT,ITEMP(I),0)
              JBIT = JBIT - IBIT
              IF (JBIT.LT.0) THEN
                  JWORD = JWORD + 1
                  JBIT  = JBIT + 64
              END IF
   70     CONTINUE
      ELSEIF(IBIT.EQ.0) THEN
C
C         IF BIT MAP INCLUDED FILL TO END OF CURRENT 32 BIT WORD
C
          IF (OBTMIS) THEN
              DO 80 I = 0,JBIT
                  ICOMP(JWORD) = IBSET(ICOMP(JWORD),I)
   80         CONTINUE
          END IF
      END IF
C
C     CALCULATE LENGTH IN 32 BIT WORDS
C
      NOP = JWORD*64 - JBIT - 1
      NOP = (NOP+31)/32
C
C     SET FLAGS TO INDICATE BIT MAPS
C
      IF (OBTZER) IBIT = IBIT + 128
      IF (OBTMIS) IBIT = IBIT + 32

      RETURN
! Fujitsu vectorization directive
!OCL NOVREC

      END
      SUBROUTINE XPND(IX,ICOMP,FIELD,APREC,IBIT,BASE,NOP,RMDI)
CFPP$ NOCONCUR R
CFPP$ EXPAND (EXTRIN)
      INTEGER IX,NOP,ICOMP(NOP+1),IBIT
      REAL FIELD(IX),APREC,BASE,RMDI
      INTEGER I,JWORD,JBIT,JJ,NPOINT,NZERO
      INTEGER IMAP(IX),IMIS(IX),IZERO(IX),IMIN(IX),ITEMP(IX)
      INTEGER IGATH(IX)
      REAL ATEMP(IX)
      LOGICAL OBTMIN,OBTMIS,OBTZER,OBTMAP,BTEST
      EXTERNAL INSTIN,EXTRIN
C
C     CRAY VERSION C.1.1  16/11/90  P J SMITH
C
C         IX   =  LENTH OF FIELD
C      ICOMP   =  DATA FOR EXPANSION
C      FIELD   =  FIELD OF EXPANDED DATA
C      APREC   =  PRECISION
C       IBIT   =  NUMBER OF BITS INTO WHICH DATA PACKED
C       BASE   =  REFERENCE (MIN) VALUE FOR PACKED DATA
C        NOP   =  SIZE OF COMP
C       RMDI   =  MISSING DATA INDICATOR
C
C
C     INITIALISE VARIABLES/TEMP ARRAYS
C
      OBTMAP   = .FALSE.
      OBTMIS   = .FALSE.
      OBTMIN   = .FALSE.
      OBTZER   = .FALSE.
      DO I = 1,IX
          IMAP(I)  = 1
          IMIS(I)  = 0
          IMIN(I)  = 0
          IZERO(I) = 0
      END DO
C
C     CHECK IF BITMAP USED FOR ZERO VALUES
C
      IF (IBIT.GE.128) THEN
          OBTZER = .TRUE.
          OBTMAP = .TRUE.
          IBIT   = IBIT-128
      ELSE
          OBTMIN = .FALSE.
      ENDIF
C
C     CHECK IF BITMAP USED FOR MINIMUM VALUES (NOT CURRENTLY USED)
C
      IF (IBIT.GE.64) THEN
          OBTMIN = .TRUE.
          OBTMAP = .TRUE.
          IBIT   = IBIT - 64
      ELSE
          OBTMIN = .FALSE.
      ENDIF
C
C     CHECK IF BITMAP USED FOR MISSING DATA
C
      IF (IBIT.GE.32) THEN
          OBTMIS = .TRUE.
          OBTMAP = .TRUE.
          IBIT   = IBIT - 32
      ELSE
          OBTMIS = .FALSE.
      ENDIF
C
C     SET START POSITION IN ICOMP
C
      JWORD = 1
      JBIT  = 63
C
C     EXTRACT BITMAPS
C
      IF (OBTMIS) THEN
C
C         EXTRACT MISSING DATA BITMAP
C
          DO 10 I=1,IX
              IF (BTEST(ICOMP(JWORD),JBIT)) THEN
                  IMIS(I) = 1
                  IMAP(I) = 0
              ELSE
                  IMIS(I) = 0
              ENDIF
              IF(JBIT.GT.0) THEN
                  JBIT    = JBIT - 1
              ELSE
                  JBIT    = 63
                  JWORD   = JWORD + 1
              ENDIF
   10     CONTINUE
      ENDIF
      IF (OBTMIN) THEN
C
C         EXTRACT MINIMUM VALUE BITMAP (NOT USED AT PRESENT)
C
          DO 20 I=1,IX
              IF(BTEST(ICOMP(JWORD),JBIT)) THEN
                  IMIN(I) = 1
                  IMAP(I) = 0
              ELSE
                  IMIN(I) = 0
              ENDIF
              IF(JBIT.GT.0) THEN
                  JBIT    = JBIT - 1
              ELSE
                  JBIT    = 63
                  JWORD   = JWORD + 1
              ENDIF
   20     CONTINUE
      ENDIF
      IF (OBTZER) THEN
C
C         EXTRACT ZERO VALUE BITMAP
C
          DO 30 I=1,IX
              IF(BTEST(ICOMP(JWORD),JBIT)) THEN
                  IZERO(I)= 0
              ELSE
                  IZERO(I)= 1
                  IMAP(I) = 0
              ENDIF
              IF(JBIT.GT.0) THEN
                  JBIT    = JBIT - 1
              ELSE
                  JBIT    = 63
                  JWORD   = JWORD + 1
              ENDIF
   30     CONTINUE
      ENDIF
C
C     IF BIT MAP USED FIND NUMBER OF POINTS STORED
C     AND RESET POINTERS TO BEGINNING OF 32 BIT BOUNDARY
C
      IF (OBTMAP) THEN
          NPOINT = 0
          DO 40 I=1,IX
              IF (IMAP(I).EQ.1) NPOINT = NPOINT + 1
   40     CONTINUE
          IF (JBIT.NE.63) THEN
              IF (JBIT.GE.31) THEN
                  JBIT  = 31
              ELSE
                  JBIT  = 63
                  JWORD = JWORD + 1
              ENDIF
          ENDIF
      ELSE
          NPOINT = IX
      ENDIF
      IF (IBIT.GT.0) THEN
C
C         UNPACK SCALED VALUES TO TEMP ARRAY
C
          DO 50 I=1,NPOINT
              CALL EXTRIN(ICOMP(JWORD),2,JBIT,IBIT,ITEMP(I),0)
              JBIT = JBIT - IBIT
              IF (JBIT.LT.0) THEN
                  JWORD = JWORD + 1
                  JBIT  = JBIT + 64
              ENDIF
   50     CONTINUE
C
C         ADD DIFFERENCES TO MINIMUM VALUE AND UNSCALE
C
          DO 60 I=1,NPOINT
              ATEMP(I)=ITEMP(I)*APREC+BASE
   60     CONTINUE
C
C         MOVE INTO UNPACKED ARRAY FIELD
C
C FIRST GET GATHER INDEX
C
          JJ  =0
          DO I=1,IX
            IF (IMAP(I).EQ.1) THEN
              JJ       =JJ+1
              IGATH(JJ)=I
            END IF
          END DO
C                                                                       
          DO I=1,IX
          FIELD(I)=0.
          END DO

          DO I=1,JJ
          FIELD(IGATH(I)) = ATEMP(I)
          END DO
C
C         IF MINIMUMS BIT MAPPED FILL ZEROS IN FIELD WITH BASE
C
          IF (OBTMIN) THEN
              DO 80 I=1,IX
                  IF(IMIN(I).EQ.1) FIELD(I) = BASE
   80         CONTINUE
          ENDIF
C
C         IF MISSING DATA BIT MAPPED FILL ZEROS IN FIELD WITH MDI
C
          IF (OBTMIS) THEN
              DO 90 I=1,IX
                  IF(IMIS(I).EQ.1) FIELD(I) = RMDI
   90         CONTINUE
          ENDIF

      ELSEIF (IBIT.EQ.0) THEN

C
C         ALL POINTS IN ROW HAVE SAME VALUE E.G. POLE ROW
C
          DO 100 I=1,IX
              FIELD(I)=BASE
  100     CONTINUE
C
C         IF MISSING DATA BIT MAPPED FILL ZEROS IN FIELD WITH MDI
C
          IF (OBTMIS) THEN
              DO 110 I=1,IX
                  IF(IMIS(I).EQ.1) FIELD(I) = RMDI
  110         CONTINUE
          ENDIF
      ENDIF
C
      RETURN

      END
      SUBROUTINE INSTIN(ICOMP,IWORD,ISTART,NBIT,INUM,ISIGN)
CFPP$ NOCONCUR R
      INTEGER IWORD,ICOMP(IWORD),ISTART,NBIT,INUM,ISIGN
      INTEGER ISNUM,ISCOMP,NUM
      EXTERNAL MOVBIT
C
C     CRAY VERSION C.1.1  16/11/90  P J SMITH
C
C
C     SUBROUTINE TO INSERT AN INTEGER VALUE INTO NBITS OF ICOMP
C     =========================================================
C
C     ICOMP  = FULLWORD ARRAY TO WHICH AN NBIT INTEGER IS TO BE ADDED
C     IWORD  = DIMESION OF ICOMP (ENABLES INTEGER TO CROSS WORDS)
C     ISTART = POSITION IN ICOMP WHERE AN NBIT INTEGER TO BE ADDED
C     NBIT   = NUMBER OF BITS INTO WHICH INTEGER IS TO BE STORED
C     INUM   = FULLWORD INTEGER TO BE ADDED TO ICOMP
C     ISIGN  = INDICATOR FOR STORAGE IN ICOMP: 0= +VE INTEGER
C                                              1= SIGNED INTEGER
C                                              2= SIGN BIT & +VE INT.
C
C     START POSITIONS IN ICOMP
C     |               |               |               |
C     |       ^       |       ^       |       ^       |       ^       ^
C     6666555555555544444444443333333333222222222211111111110000000000
C     3210987654321098765432109876543210987654321098765432109876543210
C
C     0000000001111111111222222222233333333334444444444555555555566666
C     1234567890123456789012345678901234567890123456789012345678901234
C
C
      ISNUM  = 64 - NBIT + 1
      ISCOMP = 64 - ISTART
      NUM    = NBIT
C
C     MOVE SIGN BIT
C
      IF (ISIGN.EQ.2) THEN
          CALL MOVBIT(INUM,1,1,ICOMP(1),ISCOMP)
          NUM    = NBIT - 1
          ISNUM  = ISNUM + 1
          ISCOMP = ISCOMP + 1
          INUM   = -INUM
      ENDIF
C
C     MOVE INTEGER
C
      CALL MOVBIT(INUM,ISNUM,NUM,ICOMP(1),ISCOMP)
C
      RETURN

      END
      SUBROUTINE EXTRIN(ICOMP,IWORD,ISTART,NBIT,INUM,ISIGN)
CFPP$ NOCONCUR R
      INTEGER IWORD,ICOMP(IWORD),ISTART,NBIT,INUM,ISIGN
      INTEGER ISCOMP,ISNUM,NUM,IBIT
      LOGICAL BTEST
      INTEGER IBSET,IBCLR
      EXTERNAL MOVBIT
C
C     CRAY VERSION C.1.1  16/11/90  P J SMITH
C
C
C     SUBROUTINE TO EXTRACT AN INTEGER VALUE FROM NBITS OF ICOMP
C     ==========================================================
C
C     ICOMP  = ARRAY FROM WHICH AN NBIT INTEGER IS TO BE RETRIEVED
C     IWORD  = DIMESION OF ICOMP (ENABLES INTEGER TO CROSS WORDS)
C     ISTART = POSITION IN ICOMP WHERE AN NBIT INTEGER STARTS
C     NBIT   = NUMBER OF BITS IN WHICH INTEGER HAS BEEN STORED
C     INUM   = INTEGER EXTRACTED FROM ICOMP
C     ISIGN  = INDICATOR FOR STORAGE IN ICOMP: 0= +VE INTEGER
C                                              1= SIGNED INTEGER
C                                              2= SIGN BIT & +VE INT.
C
C     START POSITIONS IN ICOMP
C     |               |               |               |
C     |       ^       |       ^       |       ^       |       ^       ^
C     6666555555555544444444443333333333222222222211111111110000000000
C     3210987654321098765432109876543210987654321098765432109876543210
C
C     0000000001111111111222222222233333333334444444444555555555566666
C     1234567890123456789012345678901234567890123456789012345678901234
C
C
      INUM   = 0
      ISCOMP = 64 - ISTART
      ISNUM  = 64 - NBIT + 1
      NUM    = NBIT
C
C     MOVE INTEGER
C
      IF (ISIGN.EQ.0) THEN
C
C         POSITIVE INTEGER
C
          CALL MOVBIT(ICOMP,ISCOMP,NUM,INUM,ISNUM)

      ELSEIF(ISIGN.EQ.1) THEN
C
C         SIGNED INTEGER
C
C         SIGN
          CALL MOVBIT(ICOMP,ISCOMP,1,INUM,1)
          ISCOMP = ISCOMP + 1
          ISNUM  = ISNUM +1
          NUM    = NUM -1
C         INTEGER
          CALL MOVBIT(ICOMP,ISCOMP,NUM,INUM,ISNUM)
C         SET UNDIFINED IF INUM NEGATIVE
          IF (INUM.LT.0) THEN
              DO IBIT=NUM,63
                  INUM = IBSET(INUM,IBIT)
              END DO
          ENDIF

      ELSEIF(ISIGN.EQ.2) THEN
C
C         SIGN BIT PLUS POSITIVE INTEGER
C
          ISCOMP = ISCOMP + 1
          ISNUM  = ISNUM +1
          NUM    = NUM -1
C         INTEGER
          CALL MOVBIT(ICOMP,ISCOMP,NUM,INUM,ISNUM)
C         SIGN
          IF (BTEST(ICOMP(1),ISTART)) INUM = -INUM
      ENDIF

      RETURN

      END
