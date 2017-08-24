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
CLL  SUBROUTINE PR_FIXHD---------------------------------------
CLL
CLL  Purpose: Prints out fixed length header record and checks
CLL           validity of information.
CLL
CLL AD, DR, SI  <- programmer of some or all of previous code or changes
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  Date
CLL
CLL  3.1   22/12/92     Allow use by ancillary field headers
CLL                     Author A. Dickinson    Reviewer C. Wilson
CLL   3.2    13/07/93 Changed CHARACTER*(*) to CHARACTER*(80) for
CLL                   portability.  Author Tracey Smith.
CLL  3.2   19/04/93     Code for new real missing data indicator.
CLL                     Author T.Johns         Reviewer A.Dickinson
CLL  3.4   5/07/94  Allowed for new PF model vert coord (density*r*r)
CLL                 and for Varobs file.    Author Colin Parrett.
CLL  4.0  21/11/95   Allowed for Covariance files.   C.Parrett
CLL  4.1  09/04/96  Introduce wave sub-model.  RTHBarnes.
CLL  4.1   16/04/96  Allowed for OPS Obstore files.  Author Colin Parret
!LL  4.1  23/05/96   Removed check of data length for MPP code  P.Burton
CLL  4.1  22/05/96   Print out UM Version Number   D. Robinson
CLL
CLL  Programming standard:
CLL           Unified Model Documentation Paper No 3
CLL           Version No 2  09/09/91
CLL
CLL  System component: C25
CLL
CLL  System task: F3
CLL
CLL  Documentation:
CLL           Unified Model Documentation Paper No F3
CLL           Version No 5 9/2/90
CLL
CLL------------------------------------------------------------
C*L Arguments:-------------------------------------------------
      SUBROUTINE PR_FIXHD
     *(FIXHD,LEN_FIXHD,LEN_INTHD,LEN_REALHD,LEN1_LEVDEPC
     *,LEN2_LEVDEPC,LEN1_ROWDEPC,LEN2_ROWDEPC,LEN1_COLDEPC,LEN2_COLDEPC
     *,LEN1_FLDDEPC,LEN2_FLDDEPC,LEN_EXTCNST,LEN_DUMPHIST,LEN_CFI1
     *,LEN_CFI2,LEN_CFI3,LEN1_LOOKUP,LEN2_LOOKUP,LEN_DATA
     *,ICODE,CMESSAGE)

      IMPLICIT NONE

      INTEGER
     * LEN_FIXHD     !IN Length of fixed length header
     *,LEN_INTHD     !IN Length of integer header
     *,LEN_REALHD    !IN Length of real header
     *,LEN1_LEVDEPC  !IN 1st dim of level dep consts
     *,LEN2_LEVDEPC  !IN 2nd dim of level dep consts
     *,LEN1_ROWDEPC  !IN 1st dim of row dep consts
     *,LEN2_ROWDEPC  !IN 2nd dim of row dep consts
     &,LEN1_COLDEPC  !IN 1st dim of column dep consts
     &,LEN2_COLDEPC  !IN 2nd dim of column dep consts
     &,LEN1_FLDDEPC  !IN 1st dim of field dep consts
     &,LEN2_FLDDEPC  !IN 2nd dim of field dep consts
     &,LEN_EXTCNST   !IN Length of extra constants
     &,LEN_DUMPHIST  !IN Length of history block
     &,LEN_CFI1      !IN Length of comp field index 1
     &,LEN_CFI2      !IN Length of comp field index 2
     &,LEN_CFI3      !IN Length of comp field index 3
     &,LEN1_LOOKUP   !IN 1st dim of lookup
     &,LEN2_LOOKUP   !IN 2nd dim of lookup

      INTEGER
     * FIXHD(LEN_FIXHD) !IN Fixed length header
     &,LEN_DATA         !IN Length of real data
     *,ICODE          !OUT Return code; successful=0
     *                !                 error > 0

      CHARACTER*(80)
     * CMESSAGE       !OUT Error message if ICODE > 0

C -------------------------------------------------------------
C Workspace usage:---------------------------------------------
C None
C -------------------------------------------------------------
C*L External subroutines called:-------------------------------
C None
C*-------------------------------------------------------------
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
C Local variables:---------------------------------------------
      INTEGER I
C--------------------------------------------------------------

      ICODE=0
      CMESSAGE=' '

      WRITE(6,'('' '')')
      WRITE(6,'('' FIXED LENGTH HEADER'')')
      WRITE(6,'('' -------------------'')')

      WRITE(6,'('' Dump format version'',I6)')FIXHD(1)
      WRITE(6,'('' UM Version No      '',I6)')FIXHD(12)

      IF(FIXHD(2).EQ.1)THEN
        WRITE(6,'('' Atmospheric data'')')
      ELSE IF(FIXHD(2).EQ.2)THEN
        WRITE(6,'('' Oceanic data'')')
      ELSE IF (FIXHD(2).EQ.4) THEN  
        WRITE(6,'('' Wave sub-model data'')')  
      ELSE
        WRITE(6,'('' ***FATAL ERROR*** Invalid data type: FIXHD(2)='',
     *I9)')FIXHD(2)
       ICODE=4
       CMESSAGE='PR_FIXHD: Consistency check'
       RETURN
      ENDIF

      IF(FIXHD(3).EQ.1)THEN
        WRITE(6,'('' On hybrid levels'')')
      ELSEIF(FIXHD(3).EQ.2)THEN
        WRITE(6,'('' On sigma levels'')')
      ELSEIF(FIXHD(3).EQ.3)THEN
        WRITE(6,'('' On pressure levels'')')
      ELSEIF(FIXHD(3).EQ.4)THEN
        WRITE(6,'('' On depth levels'')')
      ELSEIF(FIXHD(3).EQ.5)THEN
        WRITE(6,'('' Charney-Phillips on radius levels'')')
      ELSEIF (FIXHD(3).EQ.6) THEN
        WRITE(6,'('' Wave model direction levels and frequency pseudo-le
     &vels'')')   
      ELSEIF(FIXHD(3).EQ.IMDI.AND.FIXHD(5).EQ.4)THEN
        WRITE(6,'('' Missing data indicator used for vert coord type
     *'')')
      ELSE
        WRITE(6,'('' ***FATAL ERROR*** Invalid level type: FIXHD(3)='',
     *I9)')FIXHD(3)
       ICODE=4
       CMESSAGE='PR_FIXHD: Consistency check'
       RETURN
      ENDIF

      IF(FIXHD(4).EQ.0)THEN
        WRITE(6,'('' Over global domain'')')
      ELSEIF(FIXHD(4).EQ.1)THEN
        WRITE(6,'('' Over N. Hemispheric domain'')')
      ELSEIF(FIXHD(4).EQ.2)THEN
        WRITE(6,'('' Over S. Hemispheric domain'')')
      ELSEIF(FIXHD(4).EQ.3)THEN
        WRITE(6,'('' Over LA domain with no wrap around'')')
      ELSEIF(FIXHD(4).EQ.4)THEN
        WRITE(6,'('' Over LA domain with wrap around'')')
      ELSEIF(FIXHD(4).EQ.103)THEN
        WRITE(6,'('' Over rotated LA domain'')')
      ELSE
        WRITE(6,'('' ***FATAL ERROR*** Invalid domain: FIXHD(4)='',
     *I9)')FIXHD(4)
       ICODE=4
       CMESSAGE='PR_FIXHD: Consistency check'
       RETURN
      ENDIF

      IF(FIXHD(5).EQ.1)THEN
        WRITE(6,'('' Instantaneous dump'')')
      ELSEIF(FIXHD(5).EQ.2)THEN
        WRITE(6,'('' Meaned dump'')')
      ELSEIF(FIXHD(5).EQ.3)THEN
        WRITE(6,'('' FIELDS file'')')
      ELSEIF(FIXHD(5).EQ.4)THEN
        WRITE(6,'('' Ancillary dataset'')')
      ELSEIF(FIXHD(5).EQ.5)THEN
        WRITE(6,'('' Boundary dataset'')')
      ELSEIF(FIXHD(5).EQ.6)THEN
        WRITE(6,'('' AC Observation File'')')
      ELSEIF(FIXHD(5).EQ.7)THEN
        WRITE(6,'('' Var Observation File'')')
      ELSEIF(FIXHD(5).EQ.8)THEN
        WRITE(6,'('' Cx file (model columns at ob locations)'')')
       ELSEIF(FIXHD(5).EQ.9)THEN
         WRITE(6,'('' Covariance File'')')
      ELSE IF (FIXHD(5).EQ.10) THEN
        WRITE (6, '('' OPS Obstore file '')')
      ELSE
        WRITE(6,'('' ***FATAL ERROR*** Invalid dump type: FIXHD(5)='',
     *I9)')FIXHD(5)
       ICODE=4
       CMESSAGE='PR_FIXHD: Consistency check'
       RETURN
      ENDIF

      WRITE(6,'('' Exp No ='',I6,'' Run Id ='',I6)') FIXHD(7),FIXHD(6)

      IF(FIXHD(8).EQ.1)THEN
        WRITE(6,'('' Gregorian calendar'')')
      ELSEIF(FIXHD(8).EQ.2)THEN
        WRITE(6,'('' 360-day calendar'')')
      ELSEIF(FIXHD(8).EQ.IMDI.AND.FIXHD(5).EQ.4)THEN
        WRITE(6,'('' Missing data indcator used as calendar'')')
      ELSE
        WRITE(6,'('' *** FATAL ERROR *** Invalid calendar type: FIXHD(8)
     *='',I9)')FIXHD(8)
       ICODE=4
       CMESSAGE='PR_FIXHD: Consistency check'
       RETURN
      ENDIF

      IF(FIXHD(9).EQ.1)THEN
        WRITE(6,'('' Arakawa A grid'')')
      ELSEIF(FIXHD(9).EQ.2)THEN
        WRITE(6,'('' Arakawa B grid'')')
      ELSEIF(FIXHD(9).EQ.3)THEN
        WRITE(6,'('' Arakawa C grid'')')
      ELSEIF(FIXHD(9).EQ.4)THEN
        WRITE(6,'('' Arakawa D grid'')')
      ELSEIF(FIXHD(9).EQ.5)THEN
        WRITE(6,'('' Arakawa E grid'')')
      ELSEIF(FIXHD(9).EQ.IMDI.AND.FIXHD(5).EQ.4)THEN
        WRITE(6,'('' Missing data indicator used for grid type'')')
      ELSEIF(FIXHD(9).EQ.IMDI.AND.FIXHD(5).EQ.5)THEN
        WRITE(6,'('' Missing data indicator used for grid type'')')
      ELSE
        WRITE(6,'('' *** FATAL ERROR *** Invalid grid type: FIXHD(9)=''
     *  ,I9)')FIXHD(9)
       ICODE=4
       CMESSAGE='PR_FIXHD: Consistency check'
       RETURN
      ENDIF

      WRITE(6,'(''                 Year  Month Day Hour Min  Sec DayNo
     *'')')
      WRITE(6,'('' Data time     ='',7I5)')(FIXHD(I),I=21,27)
      WRITE(6,'('' Validity time ='',7I5)')(FIXHD(I),I=28,34)
      WRITE(6,'('' Creation time ='',7I5)')(FIXHD(I),I=35,41)

      WRITE(6,'(''                      Start   1st dim  2nd dim  1st pa
     *rm  2nd parm'')')
      WRITE(6,'('' Integer Consts   '',2I9,9X,I9)')FIXHD(100),
     *FIXHD(101),LEN_INTHD
      WRITE(6,'('' Real Consts      '',2I9,9X,I9)')FIXHD(105),
     *FIXHD(106),LEN_REALHD
      WRITE(6,'('' Level Dep Consts '',5I9)')FIXHD(110),
     *FIXHD(111),FIXHD(112),LEN1_LEVDEPC,LEN2_LEVDEPC
      WRITE(6,'('' Row Dep Consts   '',5I9)')FIXHD(115),
     *FIXHD(116),FIXHD(117),LEN1_ROWDEPC,LEN2_ROWDEPC
      WRITE(6,'('' Column Dep Consts'',5I9)')FIXHD(120),
     *FIXHD(121),FIXHD(122),LEN1_COLDEPC,LEN2_COLDEPC
      WRITE(6,'('' Fields of Consts '',5I9)')FIXHD(125),
     *FIXHD(126),FIXHD(127),LEN1_FLDDEPC,LEN2_FLDDEPC
      WRITE(6,'('' Extra Consts     '',2I9,9X,I9)')FIXHD(130),
     *FIXHD(131),LEN_EXTCNST
      WRITE(6,'('' History Block    '',2I9,9X,I9)')FIXHD(135),
     *FIXHD(136),LEN_DUMPHIST
      WRITE(6,'('' CFI No 1         '',2I9,9X,I9)')FIXHD(140),
     *FIXHD(141),LEN_CFI1
      WRITE(6,'('' CFI No 2         '',2I9,9X,I9)')FIXHD(142),
     *FIXHD(143),LEN_CFI2
      WRITE(6,'('' CFI No 3         '',2I9,9X,I9)')FIXHD(144),
     *FIXHD(145),LEN_CFI3
      WRITE(6,'('' Lookup Tables    '',5I9)')FIXHD(150),
     *FIXHD(151),FIXHD(152),LEN1_LOOKUP,LEN2_LOOKUP
      WRITE(6,'('' Model Data       '',2I9,9X,I9)')FIXHD(160),
     *FIXHD(161),LEN_DATA

C Check model parameters against header record entries

      IF(FIXHD(101).GT.0)THEN
        IF(LEN_INTHD.NE.FIXHD(101))THEN
        WRITE(6,'('' *ERROR* Integer Consts'')')
        WRITE(6,'('' Parameter and header values dont match'')')
        ICODE=4
        CMESSAGE='PR_FIXHD: Consistency check'
        RETURN
        ENDIF
      ENDIF

      IF(FIXHD(105).GT.0)THEN
        IF(LEN_REALHD.NE.FIXHD(106))THEN
        WRITE(6,'('' *ERROR* Real Consts'')')
        WRITE(6,'('' Parameter and header values dont match'')')
        ICODE=4
        CMESSAGE='PR_FIXHD: Consistency check'
        RETURN
        ENDIF
      ENDIF

      IF(FIXHD(110).GT.0)THEN
       IF(LEN1_LEVDEPC.NE.0)THEN
        IF(LEN1_LEVDEPC.NE.FIXHD(111).OR.LEN2_LEVDEPC.NE.FIXHD(112))THEN
        WRITE(6,'('' *ERROR* Level Dep Consts'')')
        WRITE(6,'('' Parameter and header values dont match'')')
        ICODE=4
        CMESSAGE='PR_FIXHD: Consistency check'
        RETURN
        ENDIF
       ENDIF
      ENDIF

      IF(FIXHD(115).GT.0)THEN
       IF(LEN1_ROWDEPC.NE.0)THEN
        IF(LEN1_ROWDEPC.NE.FIXHD(116).OR.LEN2_ROWDEPC.NE.FIXHD(117))THEN
        WRITE(6,'('' *ERROR* Row Dep Consts'')')
        WRITE(6,'('' Parameter and header values dont match'')')
        ICODE=4
        CMESSAGE='PR_FIXHD: Consistency check'
        RETURN
        ENDIF
       ENDIF
      ENDIF

      IF(FIXHD(120).GT.0)THEN
       IF(LEN1_COLDEPC.NE.0)THEN
        IF(LEN1_COLDEPC.NE.FIXHD(121).OR.LEN2_COLDEPC.NE.FIXHD(122))THEN
        WRITE(6,'('' *ERROR* Column Dep Consts'')')
        WRITE(6,'('' Parameter and header values dont match'')')
        ICODE=4
        CMESSAGE='PR_FIXHD: Consistency check'
        RETURN
        ENDIF
       ENDIF
      ENDIF

      IF(FIXHD(125).GT.0)THEN
       IF(LEN1_FLDDEPC.NE.0)THEN
        IF(LEN1_FLDDEPC.NE.FIXHD(126).OR.LEN2_FLDDEPC.NE.FIXHD(127))THEN
        WRITE(6,'('' *ERROR* Fields of Consts'')')
        WRITE(6,'('' Parameter and header values dont match'')')
        ICODE=4
        CMESSAGE='PR_FIXHD: Consistency check'
        RETURN
        ENDIF
       ENDIF
      ENDIF

      IF(FIXHD(130).GT.0)THEN
       IF(LEN_EXTCNST.NE.0)THEN
        IF(LEN_EXTCNST.NE.FIXHD(131))THEN
        WRITE(6,'('' *ERROR* Extra Consts'')')
        WRITE(6,'('' Parameter and header values dont match'')')
        ICODE=4
        CMESSAGE='PR_FIXHD: Consistency check'
        RETURN
        ENDIF
       ENDIF
      ENDIF

      IF(FIXHD(135).GT.0)THEN
       IF(LEN_DUMPHIST.NE.0)THEN
        IF(LEN_DUMPHIST.NE.FIXHD(136))THEN
        WRITE(6,'('' *ERROR* History File'')')
        WRITE(6,'('' Parameter and header values dont match'')')
        ICODE=4
        CMESSAGE='PR_FIXHD: Consistency check'
        RETURN
        ENDIF
       ENDIF
      ENDIF

      IF(FIXHD(140).GT.0)THEN
       IF(LEN_CFI1.NE.0)THEN
        IF(LEN_CFI1.NE.FIXHD(141))THEN
        WRITE(6,'('' *ERROR* CFI No 1'')')
        WRITE(6,'('' Parameter and header values dont match'')')
        ICODE=4
        CMESSAGE='PR_FIXHD: Consistency check'
        RETURN
        ENDIF
       ENDIF
      ENDIF

      IF(FIXHD(142).GT.0)THEN
       IF(LEN_CFI2.NE.0)THEN
        IF(LEN_CFI2.NE.FIXHD(143))THEN
        WRITE(6,'('' *ERROR* CFI No 2'')')
        WRITE(6,'('' Parameter and header values dont match'')')
        ICODE=4
        CMESSAGE='PR_FIXHD: Consistency check'
        RETURN
        ENDIF
       ENDIF
      ENDIF

      IF(FIXHD(144).GT.0)THEN
       IF(LEN_CFI3.NE.0)THEN
        IF(LEN_CFI3.NE.FIXHD(145))THEN
        WRITE(6,'('' *ERROR* CFI No 3'')')
        WRITE(6,'('' Parameter and header values dont match'')')
        ICODE=4
        CMESSAGE='PR_FIXHD: Consistency check'
        RETURN
        ENDIF
       ENDIF
      ENDIF

      IF(FIXHD(150).GT.0)THEN
        IF(LEN2_LOOKUP.NE.IMDI)THEN
        IF(LEN1_LOOKUP.NE.FIXHD(151).OR.LEN2_LOOKUP.NE.FIXHD(152))THEN
        WRITE(6,'('' *ERROR* Lookup Table'')')
        WRITE(6,'('' Parameter and header values dont match'')')
        ICODE=4
        CMESSAGE='PR_FIXHD: Consistency check'
        RETURN
        ENDIF
        ENDIF
      ENDIF


      RETURN
      END
