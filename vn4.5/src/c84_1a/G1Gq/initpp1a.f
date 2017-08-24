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
CLL  Routine: INIT_PP  -------------------------------------------------
CLL
CLL  Purpose: Initialises direct access PP files at the start of
CLL           the run.  NB: Sequential PP files need no initialisation.
CLL
CLL  Tested under compiler:   cft77
CLL  Tested under OS version: UNICOS 5.1
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  Date
CLL   3.1  12/02/93  Modify args to allow correct setting of PP_FIXHD(5)
CLL   3.2    13/07/93 Changed CHARACTER*(*) to CHARACTER*(80) for
CLL                   portability.  Author Tracey Smith.
CLL   3.3  05/10/93  Flush buffer of each pp file after initialisation
CLL                  to ensure complete file ready for re-start in the
CLL                  event of a 'hard' failure. R. Rawlins
!LL  4.3   30/04/97  Added code to use UM_SECTOR_SIZE to make transfers
!LL                  well-formed.
!LL                  B. Carruthers  Cray Research.
CLL
CLL  Programming standard: UM Doc Paper 3, version 1 (15/1/90)
CLL
CLL  Logical components covered: D401
CLL
CLL  Project task:
CLL
CLL  External documentation: On-line UM document C61 - Zonal mean
CLL                          calculations.
CLL
CLL  -------------------------------------------------------------------
C*L  Interface and arguments: ------------------------------------------
C
      SUBROUTINE INIT_PP ( FTN_UNIT,FILE_TYPE_LETTER,
     &                     LEN1_LOOKUP,PP_LEN2_LOOKUP,FIXHD,
     1                     INTHD,REALHD,LEVDEPC,LEN_FIXHD,LEN_INTHD,
     2                     LEN_REALHD,LEN1_LEVDEPC,LEN2_LEVDEPC,
     3                     ICODE,CMESSAGE)
C
      IMPLICIT NONE
C
      CHARACTER*1
     &    FILE_TYPE_LETTER    ! IN  - File type (p-PP, b-bndry)
      INTEGER
     1    FTN_UNIT            ! IN  - Fortran unit number
     2,   LEN1_LOOKUP         ! IN  - Size of PP header
     3,   PP_LEN2_LOOKUP      ! IN  - Max allowable fields
     4,   LEN_FIXHD           ! IN    LENGTH OF FIXED CONSTANTS
     5,   LEN_INTHD           ! IN    LENGTH OF INTEGER CONSTANTS
     6,   LEN_REALHD          ! IN    LENGTH OF REAL CONSTANTS
     7,   LEN1_LEVDEPC        ! IN    LENGTH OF 1st Dim of lev depndt
     8,   LEN2_LEVDEPC        ! IN    LENGTH OF 2nd Dim of lev depndt
     9,   ICODE               ! OUT - Error exit code
     A,   PP_LEN_INTHD        ! OUT - Length of PP FILE integer header
     B,   PP_LEN_REALHD       ! OUT - Length of PP FILE real header
     C,   PP_LEN2_LEVDEPC     ! OUT - Length of 2nd dim of PP lev dep
C
      PARAMETER
     1   (PP_LEN_INTHD=15
     2,   PP_LEN_REALHD=6
     3,   PP_LEN2_LEVDEPC=4)
C
      INTEGER
     *    FIXHD(LEN_FIXHD)          ! IN    ARRAY OF FIXED CONSTANTS
     *,   INTHD(LEN_INTHD)          ! IN    ARRAY OF integer CONSTANTS
     *,   LEVDEPC(LEN1_LEVDEPC*LEN2_LEVDEPC)  ! IN LEV DEP CONSTANTS
     *,   PP_FIXHD(LEN_FIXHD)       ! OUT   ARRAY of fixed constants
     *,   PP_INTHD(PP_LEN_INTHD)    ! OUT   ARRAY of integer constants
     *,   PP_LEVDEPC(LEN1_LEVDEPC*PP_LEN2_LEVDEPC) ! OUT level dep cts
C
      REAL
     *    REALHD(LEN_REALHD)        ! IN    ARRAY OF REAL CONSTANTS
     *,   PP_REALHD(PP_LEN_REALHD)  ! OUT   ARRAY OF REAL CONSTANTS
C
      CHARACTER*80
     1    CMESSAGE            ! OUT - Error message
C
C*----------------------------------------------------------------------
C
C  External subroutines
C
      EXTERNAL SETPOS,IOERROR,POSERROR,BUFFOUT,FLUSH_BUFFER
C
C  Local variables
C
      INTEGER IPPLOOK(LEN1_LOOKUP,PP_LEN2_LOOKUP)
C
cdir$ cache_align pp_fixhd, pp_inthd, pp_realhd, pp_levdepc, ipplook
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
      INTEGER
     1       II,JJ,IWA,IX,LEN_IO,START_BLOCK  !
      REAL A_IO
! History:
! Version  Date  Comment
!  3.4   18/5/94 Add PP missing data indicator. J F Thomson
C*L------------------COMDECK C_MDI-------------------------------------
      REAL    RMDI_PP   ! PP missing data indicator (-1.0E+30)
      REAL    RMDI_OLD  ! Old real missing data indicator (-32768.0)
      REAL    RMDI      ! New real missing data indicator (-2**30)
      INTEGER IMDI      ! Integer missing data indicator

      PARAMETER(
     & RMDI_PP  =-1.0E+30,
     & RMDI_OLD =-32768.0,
     & RMDI =-32768.0*32768.0,
     & IMDI =-32768)
C*----------------------------------------------------------------------
CL----------------------------------------------------------------------
CL 1. Reserve space
CL
      DO 1 II=1,PP_LEN2_LOOKUP
      DO 2 JJ=1,LEN1_LOOKUP
      IPPLOOK(JJ,II)=-99
    2 CONTINUE
    1 CONTINUE
CL----------------------------------------------------------------------
CL 1.1 Set up FIXED header record for the PP FILE
CL
      DO 3 II=1,LEN_FIXHD
      PP_FIXHD(II)=FIXHD(II)
    3 CONTINUE
      IF (FILE_TYPE_LETTER.EQ.'p') THEN
        PP_FIXHD(5)=3
      ELSEIF (FILE_TYPE_LETTER.EQ.'b') THEN
        PP_FIXHD(5)=5
      ELSE
        ICODE=100
        CMESSAGE='INIT_PP  : Unknown output file type letter'
        RETURN
      ENDIF
      PP_FIXHD(101)=PP_LEN_INTHD
      PP_FIXHD(105)=PP_FIXHD(100)+PP_FIXHD(101)
      PP_FIXHD(106)=PP_LEN_REALHD
      PP_FIXHD(110)=PP_FIXHD(105)+PP_FIXHD(106)
      PP_FIXHD(111)=LEN1_LEVDEPC
      PP_FIXHD(112)=PP_LEN2_LEVDEPC
      PP_FIXHD(115)=0
      PP_FIXHD(120)=0
      PP_FIXHD(125)=0
      PP_FIXHD(130)=0
      PP_FIXHD(135)=0
      PP_FIXHD(140)=0
      PP_FIXHD(142)=0
      PP_FIXHD(144)=0
      PP_FIXHD(150)=PP_FIXHD(110)+ PP_FIXHD(111)*PP_FIXHD(112)
      PP_FIXHD(151)=LEN1_LOOKUP
      PP_FIXHD(152)=PP_LEN2_LOOKUP
      pp_fixhd(160)=     ! make sure the data starts on a sector bndry
     2 ((pp_fixhd(150)+pp_len2_lookup*len1_lookup+um_sector_size-1)/
     3 um_sector_size)*um_sector_size+1
CL----------------------------------------------------------------------
CL 1.2 Set up INTEGER constants record for the PP FILE
CL
      IF(PP_FIXHD(5).LE.2) THEN !  set all values initially to MDI
        DO II=1,PP_LEN_INTHD
          PP_INTHD(II)=INTHD(21)
        ENDDO
      ELSE
        DO II=1,PP_LEN_INTHD
          PP_INTHD(II)=IMDI
        ENDDO
      ENDIF

      PP_INTHD(6)=INTHD(6)
      PP_INTHD(7)=INTHD(7)
      PP_INTHD(8)=INTHD(8)
      PP_INTHD(9)=INTHD(9)
      PP_INTHD(10)=INTHD(10)
      PP_INTHD(13)=INTHD(13)
CL----------------------------------------------------------------------
CL 1.3 Set up REAL constants record for the PP FILE
CL
      PP_REALHD(1)=REALHD(1)
      PP_REALHD(2)=REALHD(2)
      PP_REALHD(3)=REALHD(3)
      PP_REALHD(4)=REALHD(4)
      PP_REALHD(5)=REALHD(5)
      PP_REALHD(6)=REALHD(6)
CL----------------------------------------------------------------------
CL 1.4 Set up LEVEL DEPENDANT constants record for the PP FILE
CL
      DO 5 II=1,LEN1_LEVDEPC*PP_LEN2_LEVDEPC
      PP_LEVDEPC(II)=LEVDEPC(II)
    5 CONTINUE
CL----------------------------------------------------------------------
CL 2.1 BUFFER OUT Header Records starting with the FIXED LENGTH
CL
      CALL BUFFOUT(FTN_UNIT,PP_FIXHD(1),LEN_FIXHD,LEN_IO,A_IO)
        IF(A_IO.NE.-1.0.OR.LEN_IO.NE.LEN_FIXHD) THEN
           CALL IOERROR('bufferout of fixed length header',A_IO,LEN_IO,
     &                    LEN_FIXHD)
           CMESSAGE='INIT_PP:I/O error'
           ICODE=1
           RETURN
        ENDIF
      START_BLOCK=LEN_FIXHD+1
CL----------------------------------------------------------------------
CL 2.2 BUFFER OUT Integer Constants
CL

      IF(FIXHD(100).GT.0) THEN  ! Any integer constants to output ?

C Check for error in file pointers

C        WRITE(6,*)  'START_BLOCK FIXHD(100)'
C        WRITE(6,*)   START_BLOCK
C        WRITE(6,*)   FIXHD(100)
C        WRITE(6,*)   FTN_UNIT
         IF(FIXHD(100).NE.START_BLOCK) THEN  ! Check start address
            CALL POSERROR('integer constants',START_BLOCK,100,
     &      PP_FIXHD(100))
            CMESSAGE='INIT_PP:  Addressing conflict'
            ICODE=2
            RETURN
         END IF

         CALL BUFFOUT (FTN_UNIT,PP_INTHD(1),PP_FIXHD(101),LEN_IO,A_IO)

C Check for I/O errors

         IF(A_IO.NE.-1.0.OR.LEN_IO.NE.PP_FIXHD(101)) THEN
            CALL IOERROR('buffer out of integer constants',A_IO,LEN_IO
     &     ,PP_FIXHD(101))
            CMESSAGE='INIT_PP: I/O Error'
            ICODE=3
            RETURN
         END IF

         START_BLOCK=START_BLOCK+PP_FIXHD(101)

      END IF

CL----------------------------------------------------------------------
CL 2.3 BUFFER OUT Real Constants
CL

      IF(PP_FIXHD(105).GT.0) THEN   ! Any real constants to output ?

C Check for error in file pointers

        IF(PP_FIXHD(105).NE.START_BLOCK) THEN
          CALL POSERROR('real constants',START_BLOCK,100,PP_FIXHD(105))
          CMESSAGE='INIT_PP: Addressing conflict'
          ICODE=4
          RETURN
        END IF

        CALL BUFFOUT(FTN_UNIT,PP_REALHD(1),PP_FIXHD(106),LEN_IO,A_IO)

C Check for I/O errors

        IF(A_IO.NE.-1.0.OR.LEN_IO.NE.PP_FIXHD(106)) THEN
          CALL IOERROR('buffer out of real constants',A_IO,LEN_IO
     &                 ,PP_FIXHD(106))
          CMESSAGE='INIT_PP: I/O Error'
          ICODE=5
          RETURN
        END IF

        START_BLOCK=START_BLOCK+PP_FIXHD(106)

      END IF

CL----------------------------------------------------------------------
CL 2.4 BUFFER OUT Level Dependant Constants.
CL

      IF(PP_FIXHD(112).GT.0) THEN ! Any level dependant constants ?

C Check for error in file pointers

         IF(PP_FIXHD(110).NE.START_BLOCK) THEN
            CALL POSERROR('real constants',START_BLOCK,100,
     &                     PP_FIXHD(110))
            CMESSAGE='INIT_PP: Addressing conflict'
            ICODE=6
            RETURN
         END IF

         CALL BUFFOUT (FTN_UNIT,PP_LEVDEPC(1)
     &              ,PP_FIXHD(111)*PP_FIXHD(112),LEN_IO,A_IO)

C Check for I/O errors

         IF(A_IO.NE.-1.0.OR.LEN_IO.NE.(PP_FIXHD(111)*PP_FIXHD(112)
     &        ))THEN
           CALL IOERROR('buffer out of lev dep constants',A_IO,LEN_IO
     &            ,PP_FIXHD(111))
           CMESSAGE='INIT_PP: I/O Error'
           ICODE=7
           RETURN
         END IF

         START_BLOCK=START_BLOCK+ PP_FIXHD(111)*PP_FIXHD(112)

      END IF
CL----------------------------------------------------------------------
CL 2.5 BUFFER OUT Lookup Table
CL
C     IWA= 0
C     CALL SETPOS(FTN_UNIT,3,IWA,ICODE)
           IF(PP_FIXHD(152).GT.0) THEN

C Check for error in file pointers

             IF(PP_FIXHD(150).NE.START_BLOCK) THEN
               CALL POSERROR('lookup table',START_BLOCK,100,
     &              PP_FIXHD(150))
               CMESSAGE='INIT_PP: Addressing conflict'
               ICODE=8
               RETURN
             END IF

      CALL BUFFOUT (FTN_UNIT,
     *           IPPLOOK(1,1),LEN1_LOOKUP*PP_LEN2_LOOKUP,LEN_IO,A_IO)
C
C Check for I/O errors

            IF(A_IO.NE.-1.0.OR.LEN_IO.NE.(PP_FIXHD(151)*PP_FIXHD(152)))
     &          THEN
              CALL IOERROR('buffer out of PP LOOKUP TABLE ',A_IO,LEN_IO
     &            ,PP_FIXHD(152))
              CMESSAGE='INIT_PP: I/O Error'
              ICODE=9
              RETURN
            END IF
C
C Clear file buffer : force last buffer to be written to file
C  to avoid problems with continuation runs following hard failures.
C
      CALL FLUSH_BUFFER(FTN_UNIT,ICODE)
      IF(ICODE.NE.0) THEN
         CMESSAGE='INIT_PP: Problem flushing buffer'
         ICODE=10
         RETURN
      ENDIF
C
            START_BLOCK=START_BLOCK+(PP_FIXHD(151)*PP_FIXHD(152))
C
          END IF
      RETURN
      END

