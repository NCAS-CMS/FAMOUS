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
CLL  SUBROUTINE PP_FILE -----------------------------------------
CLL
CLL  Purpose:- To output a field to a PP_FILE
CLL
CLL  Tested under compiler CFT77
CLL  Tested under OS version 5.1
CLL
CLL TJ, RR      <- programmer of some or all of previous code or changes
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  Date
CLL   3.2    13/07/93 Changed CHARACTER*(*) to CHARACTER*(80) for
CLL                   portability.  Author Tracey Smith.
CLL   3.2  19/04/93  Code for new real missing data indicator (TCJ).
CLL  3.4   04/08/94  No packing indicator change from -26 to -99  PJS
!LL  4.1   22/11/96  Modify I/O calls for MPP use  P.Burton
!LL  4.3   30/04/97  Added code to use UM_SECTOR_SIZE to make transfers
!LL                  well-formed.
!LL                  B. Carruthers  Cray Research.
!LL  4.4   16/06/97  Add processing after the write, so
!LL                  that all the processors know the answer
!LL                    Author: Bob Carruthers, Cray Rsearch.
!LL  4.5   28/05/98  Code for parallel processing in COEX Packing
!LL                    Author: Paul Burton & Bob Carruthers
CLL
CLL  Programming standard: U M DOC  Paper NO. 4,
CLL
CLL  Logical components covered C4
CLL
CLL  Project TASK: C4
CLL
CLL  External documentation  C4
CLL
CLLEND-------------------------------------------------------------

C
C*L  ARGUMENTS:---------------------------------------------------
      SUBROUTINE PP_FILE(PPFIELD,LENBUF,NUM_WORDS,RMDI,COMP_ACCRCY,
     1PPHORIZ_OUT,UNITPP,IWA,N_COLS_OUT,N_ROWS_OUT,PACKING,
     2PACKING_TYPE,ICODE,CMESSAGE)
      IMPLICIT NONE


      CHARACTER*(80) CMESSAGE !OUT OUT MESSAGE FROM ROUTINE
C
      LOGICAL
     *  PACKING            !IN OVERALL Packing switch (T if pckng reqd)

      INTEGER
     *  ICODE              !IN    RETURN CODE FROM ROUTINE
     *, LENBUF             !IN     LENGTH OFF PP BUFFER
     *, UNITPP             !IN     OUTPUT PP UNIT NUMBER
     *, LEN_IO             !NOT USED, BUT NEEDED FOR BUFFOUT CALL
C
      INTEGER
     *  N_ROWS_OUT    !IN   PPHORIZ_OUT=N_ROWS_OUT*N_COLS_OUT
     *, N_COLS_OUT    !IN    PPHORIZ_OUT=N_COLS_OUT*N_ROWS_OUT
     *, NUM_OUT       !IN    NUMBER OF COMPRESSED (32 BIT) WORDS
     *, COMP_ACCRCY   !IN    PACKING ACCURACY IN POWER OF 2
     *, U_ROWS        !IN    NO OF U,V, ROWS
     *, P_ROWS        !IN    PRESS/TEMP ROWS
     *, PPHORIZ_OUT   !IN    SIZE OF OUTPUT FIELD
     *, NUM_WORDS     !IN    NUMBER OF 64 BIT WORDS WORTH OF DATA
     *, PACKING_TYPE  ! OUT set to 1 if WGDOS packing else set to zero.
C
      REAL
     *  PPFIELD(PPHORIZ_OUT)   !INOUT ARRAY TO STORE PPDATA
     *, BUFOUT(LENBUF)         !OUTPUT PP BUFFER (ROUNDED UP)
     *, RMDI                   !IN     Missing data indicator
C
C

C
cdir$ cache_align bufout
!====================== COMDECK CNTL_IO ========================
! Description:
!
!     Defines the sector size for well-formed transfers on Cray
!     Research systems.  Disk addresses must start on a sector
!     boundary, and transfers must be a number of sectors.  Disk
!     word addresses start at 0.
!
!     On the T3E, well-formed transfers must also start on a
!     cache-line boundary in memory.
!
!   4.3    30/04/97  New deck       B. Carruthers, Cray Research
!   4.4    27/10/97  Remove DATA statement. C.P. Jones
!
C
      INTEGER UM_SECTOR_SIZE    ! Sector size on disk for I/O
C
      COMMON / CNTL_IO / UM_SECTOR_SIZE
C*---------------------------------------------------------------------

C*L  WORKSPACE USAGE:-------------------------------------------------
C   DEFINE LOCAL WORKSPACE ARRAYS: 1 REAL ARRAY
C   AT FULL FIELD LENGTH
C
C*---------------------------------------------------------------------
C
C*L EXTERNAL SUBROUTINES CALLED---------------------------------------
      EXTERNAL SETPOS,COEX,BUFFOUT
C*------------------------------------------------------------------
CL  MAXIMUM VECTOR LENGTH ASSUMED IS (ROWS-1) * ROWLENGTH
CL---------------------------------------------------------------------
C----------------------------------------------------------------------
C    DEFINE LOCAL VARIABLES
      INTEGER
     *  ML            !     LONGITUDE COUNTER
     *, JL            !     LATITUDE COUNTER
     *, IWA           !     RECORD NUMBER
     *, II            !     COUNTER
     *, LENGTH_FULLWRD!     LENGTH IN BITS OF FULLWORD VAR
     *, LEN_BUF_WORDS !     NUM_WORDS ROUNDED BY 512 AND ACTUALLY

      INTEGER
     *  JJ            !     Local counter

      REAL
     *  IX            !     RETURN VALUE FROM UNIT COMMAND

      LOGICAL
     *  UV                 !
C
C
C    REMEMBER THAT BUFFER OUT STARTS AT ADDRESS 0 THUS IPPLOOK GOES
C    FROM 0 to 262143 ie THE NEXT ADDRESS SHOULD BE IWA=262144 to
C    IWA=325119 then IWA=325120 to 388095 then 388096 etc
C
C======================================================================
      LENGTH_FULLWRD=64   !   LENGTH IN BITS OF FULLWORD VAR
CL   At this point packing,if required,will be done using the WGDOS
CL   method of packing.
      PACKING_TYPE=0
C Note the value of -26 corresponds to -15 (F) in ppxref.
C The packing acuracy is scaled to allow greater accuracy.
C Packing will only be attempted if there are at least 2 points per row
C in the PPfield.
C
      IF(PACKING.AND.COMP_ACCRCY.GT.-99.AND.N_COLS_OUT.GE.2)
     &   PACKING_TYPE=1
C
      IF(PACKING_TYPE.EQ.1)THEN
        CALL COEX(PPFIELD,PPHORIZ_OUT,BUFOUT,LENBUF,N_COLS_OUT,
     &  N_ROWS_OUT,NUM_OUT,COMP_ACCRCY,.TRUE.,RMDI,LENGTH_FULLWRD)

        NUM_WORDS=(NUM_OUT+1)/2 ! Round up to the nearest 64 Bit CRAY Wd
C  COEX returns the number of IBM words needed to hold the packed data
C                             ~~~
        LEN_BUF_WORDS=((NUM_WORDS+um_sector_size-1)/um_sector_size)*
     2   um_sector_size
      ELSE  ! No packing required.
        DO 1 JJ=1,PPHORIZ_OUT
        BUFOUT(JJ) = PPFIELD(JJ)
    1   CONTINUE
        NUM_WORDS=PPHORIZ_OUT
        LEN_BUF_WORDS=LENBUF
      ENDIF
      DO JJ=NUM_WORDS+1,LEN_BUF_WORDS
        BUFOUT(JJ)= 0.0
      ENDDO
      CALL SETPOS_single(UNITPP,IWA,ICODE)
      CALL BUFFOUT_single(UNITPP,BUFOUT(1),LEN_BUF_WORDS,LEN_IO,IX)
C     WRITE(6,102) IWA,LEN_BUF_WORDS
  100 FORMAT(//,32X,'   ARRAY BUFOUT AT END OF PPOUT ',//,32(10F8.0/))
  102 FORMAT(' FROM PP_FILE    IWA  LEN_BUF_WORDS ',2I12)
C
      IF (IX.NE.-1.0.OR.LEN_IO.NE.LEN_BUF_WORDS) THEN
        CALL IOERROR('Buffer out Data Field',IX,LEN_IO,
     &                LEN_BUF_WORDS)
        CMESSAGE='PPFILE  : I/O error - PP Data Field Output'
        ICODE=7
        RETURN
      ENDIF
  999 CONTINUE
      RETURN
      END
