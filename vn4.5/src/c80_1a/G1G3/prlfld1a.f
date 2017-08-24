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
CLL  SUBROUTINE PR_LFLD----------------------------------------
CLL
CLL  Purpose: Prints out selected values from logical data
CLL           using information from associated PP header.
CLL
CLL  Written by D. Robinson
CLL
CLL  Model            Modification history:
CLL version  date
CLL   3.3  22/11/93  New routine (adapted from deck PRIFLD1A)
CLL
CLL  System component: R30/W30
CLL
CLL  System task: F3
CLL
CLL  Programming standard:
CLL           Unified Model Documentation Paper No 3
CLL           Version No 1 15/1/90
CLL
CLL  Documentation:
CLL           Unified Model Documentation Paper No F3
CLL           Version No 5 9/2/90
CLL
CLL------------------------------------------------------------
C*L Arguments:-------------------------------------------------
      SUBROUTINE PR_LFLD(LOOKUP,RLOOKUP,LEN1_LOOKUP,LD1,K)

      IMPLICIT NONE

      INTEGER
     * K             !IN Field number ie position in 2nd dim
     *               !   of LOOKUP
     *,LEN1_LOOKUP   !IN First dimension of LOOKUP table
     *,LOOKUP(LEN1_LOOKUP,*)  !IN Integer equivalence of PP LOOKUP

      REAL
     * RLOOKUP(LEN1_LOOKUP,*) !IN Real equivalence of PP LOOKUP

      LOGICAL
     * LD1(*)        !IN Kth field in data array
C -------------------------------------------------------------
C*L External subroutines called:-------------------------------
C None
C--------------------------------------------------------------
C*L Local control constants:-----------------------------------
      INTEGER
     * NS_PTS        !PARAM No of points down to print
     *,EW_PTS        !PARAM No of points across to print
      PARAMETER(NS_PTS=6,EW_PTS=5)
C -------------------------------------------------------------
C Workspace usage:---------------------------------------------
      REAL LON(EW_PTS)     ! Longitudes printed out
      INTEGER I(EW_PTS)    ! Index of values printed out
      CHARACTER*12 DASH(EW_PTS)  !Stores dashed lines
C*-------------------------------------------------------------
C Local variables:---------------------------------------------
      INTEGER
     * N_ROWS      ! No of rows in field
     *,N_COLS      ! No of colums in field
     *,ROW         ! Row number
     *,R_INC,F_INC ! No of rows/points between printed lines
     *,J,L         ! Loop counts
     *,EW_PRINT    ! No of E-W values printed out
     *,POS_MIN     ! Position of Minimum value of field
     *,POS_MAX     ! Position of Maximum value of field
     *,F_MIN       ! Minimum value of field
     *,F_MAX       ! Maximum value of field

      REAL
     * LAT         ! Latitude
C--------------------------------------------------------------

C*L------------------ COMDECK LOOKADD ----------------------------------
CLL
CLL Purpose : Contains information about the format
CLL           of the PP header
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  Date
CLL   4.0  12/09/95   Change NPERIODS to LBUSER3, BRSVD1 to BULEV,
CLL                   BRSVD2 to BHULEV and definitions for BRLEV and
CLL                   BHRLEV. Corresponding changes made to STWORK1A
CLL                   and PPHEAD1A. (Andrew Brady)   
CLL  4.0  12/10/95  Change item 45 from lbuser7 to model_code. RTHBarnes
CLL
CLL Programming standard :
CLL
CLL Logical components covered : F092
CLL
CLL Project task :
CLL
CLL External documentation:
CLL
CLLEND -----------------------------------------------------------------
C
      INTEGER
C Validity time
     &       LBYR,       ! Year
     &       LBMON,      ! Month
     &       LBDAT,      ! Day of month
     &       LBHR,       ! Hour
     &       LBMIN,      ! Minute
     &       LBDAY       ! Day number

C Data time

      INTEGER
     &       LBYRD,      ! Year
     &       LBMOND,     ! Month
     &       LBDATD,     ! Day of month
     &       LBHRD,      ! Hour
     &       LBMIND,     ! Minute
     &       LBDAYD      ! Day number

      INTEGER
     &       LBTIM,      ! Time indicator
     &       LBFT,       ! Forcast period (hours)
     &       LBLREC,     ! Length of data record
     &       LBCODE,     ! Grid type code
     &       LBHEM,      ! Hemisphere indicator
     &       LBROW,      ! Number of rows in grid
     &       LBNPT,      ! Number of points per row
     &       LBEXT,      ! Length of extra data
     &       LBPACK,     ! Packing method indicator
     &       LBREL       ! Header release number

      INTEGER
     &       LBFC,       ! Field code
     &       LBCFC,      ! Second field code
     &       LBPROC,     ! Processing code
     &       LBVC,       ! Vertical coordinate type
     &       LBRVC,      ! Coordinate type for reference level
     &       LBEXP,      ! Experiment number
     &       LBEGIN,     ! Start record
     &       LBNREC,     ! No of records-Direct access only
     &       LBPROJ,     ! Met-O-8 projection number
     &       LBTYP,      ! Met-O-8 field type
     &       LBLEV,      ! Met-O-8 level code
     &       LBRSVD1,    ! Reserved for future PP-package use
     &       LBRSVD2,    ! Reserved for future PP-package use
     &       LBRSVD3,    ! Reserved for future PP-package use
     &       LBRSVD4,    ! Reserved for future PP-package use
     &       LBSRCE      ! =1111 to indicate following apply to UM
      INTEGER
     &       DATA_TYPE,  ! Indicator for real/int or timeseries
     &       NADDR,      ! Start address in DATA_REAL or DATA_INT
     &       LBUSER3,    ! Free for user-defined function   
     &       ITEM_CODE,  ! Stash code
     &       LBPLEV,     ! Pseudo-level indicator (if defined)
     &       LBUSER6,    ! Free for user-defined function
     &       MODEL_CODE ! internal model identifier
      INTEGER
     &       BULEV,      ! Upper level boundary (Bk for ATMOS)
     &       BHULEV,     ! Upper level boundary (Ak for ATMOS)   
     &       BRSVD3,     ! Reserved for future PP-package use
     &       BRSVD4,     ! Reserved for future PP-package use
     &       BDATUM,     ! Datum value
     &       BACC,       ! (Packed fields) Packing accuracy
     &       BLEV,       ! Level
     &       BRLEV,      ! Lower level boundary (Bk for ATMOS)   
     &       BHLEV,      ! (Hybrid levels) A-level of value
     &       BHRLEV,     ! Lower level boundary (Ak for ATMOS)   
     &       BPLAT,      ! Real latitude of 'pseudo' N Pole
     &       BPLON,      ! Real longitude of 'pseudo' N Pole
     &       BGOR,       ! Grid orientation
     &       BZY,        ! Zeroth latitude
     &       BDY,        ! Latitude interval
     &       BZX,        ! Zeroth longitude
     &       BDX,        ! Longitude interval
     &       BMDI,       ! Missing data indicator
     &       BMKS        ! M,K,S scaling factor

C Mapping of MPP_LOOKUP; analogous to mapping in PP header

      INTEGER
     &       P_NADDR,    ! Address on local PE
     &       P_LBLREC    ! Local length of record

      PARAMETER (
     &       P_NADDR=1,
     &       P_LBLREC=2)
C*----------------------------------------------------------------------
C NADDR IS LOCATION IN PP-HEADER (LOOKUP) FOR START POSN OF VARIABLE
C ITEM_CODE is the location in PP header for a code defined as
C           (section number)*1000+item number
C DATA_TYPE is the location in the PP header defining data as REAL or
C           INTEGER.
C LBNPT is the location defining the number of points per row
C
      PARAMETER(
C Validity time
     &       LBYR=1,
     &       LBMON=2,
     &       LBDAT=3,
     &       LBHR=4,
     &       LBMIN=5,
     &       LBDAY=6,
C Data time
     &       LBYRD=7,
     &       LBMOND=8,
     &       LBDATD=9,
     &       LBHRD=10,
     &       LBMIND=11,
     &       LBDAYD=12)

      PARAMETER (
     &       LBTIM=13,
     &       LBFT=14,
     &       LBLREC=15,
     &       LBCODE=16,
     &       LBHEM=17,
     &       LBROW=18,
     &       LBNPT=19,
     &       LBEXT=20,
     &       LBPACK=21,
     &       LBREL=22,
     &       LBFC=23,
     &       LBCFC=24,
     &       LBPROC=25,
     &       LBVC=26,
     &       LBRVC=27)

      PARAMETER (
     &       LBEXP=28,
     &       LBEGIN=29,
     &       LBNREC=30,
     &       LBPROJ=31,
     &       LBTYP=32,
     &       LBLEV=33,
     &       LBRSVD1=34,
     &       LBRSVD2=35,
     &       LBRSVD3=36,
     &       LBRSVD4=37,
     &       LBSRCE=38,
     &       DATA_TYPE=39,
     &       NADDR=40,
     &       LBUSER3=41,    
     &       ITEM_CODE=42,
     &       LBPLEV=43,
     &       LBUSER6=44,
     &       MODEL_CODE=45)

      PARAMETER (
     &       BULEV=46,
     &       BHULEV=47, 
     &       BRSVD3=48,
     &       BRSVD4=49,
     &       BDATUM=50,
     &       BACC=51,
     &       BLEV=52,
     &       BRLEV=53,
     &       BHLEV=54,
     &       BHRLEV=55,
     &       BPLAT=56,
     &       BPLON=57,
     &       BGOR=58,
     &       BZY=59,
     &       BDY=60,
     &       BZX=61,
     &       BDX=62,
     &       BMDI=63,
     &       BMKS=64)

C
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

CL Internal structure: None

C Initialise string used to create table boundaries
      DO 50 J=1,EW_PTS
50    DASH(J)='------------'

      IF(LOOKUP(LBCODE,K).EQ.IMDI) THEN
C       IF LBCODE IS MISSING DATA, ASSUME THAT THE FIELD IN DUMP
C       HAS NOT BEEN WRITTEN TO BY STASH.
C       THIS SHOULD ONLY OCCUR TO DIAGNOSTIC PARTS OF THE DUMP BEFORE
C       FIRST WRITE BY STASH TO THAT AREA/HEADER.
        WRITE(6,*) 'MESSAGE FROM PR_LFLD'
        WRITE(6,*) 'LBCODE NOT SET; ASSUME DATA NOT SET. NO PRINT'
        RETURN
      END IF

C No of rows and columns in field
      N_ROWS=LOOKUP(LBROW,K)
      N_COLS=LOOKUP(LBNPT,K)


      IF(N_COLS.NE.0.AND.N_COLS.NE.IMDI)THEN

C No of E-W values to be printed
      EW_PRINT=MIN(N_COLS,EW_PTS)

C Calculate longitudes and addresses of values to be printed from 1st ro
      I(1)=1
      LON(1)=RLOOKUP(BZX,K)+RLOOKUP(BDX,K)
      DO 100 J=1,EW_PTS-2
      I(J+1)=I(J)+N_COLS/(EW_PTS-1)
      LON(J+1)=LON(J)+RLOOKUP(BDX,K)*(N_COLS/(EW_PTS-1))
100   CONTINUE
      I(EW_PTS)=N_COLS
      LON(EW_PTS)=RLOOKUP(BZX,K)+RLOOKUP(BDX,K)*N_COLS

C Initialise row and field pointers
      ROW=1
      LAT=RLOOKUP(BZY,K)+RLOOKUP(BDY,K)
      R_INC=N_ROWS/(NS_PTS-1)
      F_INC=R_INC*N_COLS

C Print 1st row
      WRITE(6,'(14X,9A12)')(DASH(J),J=1,EW_PRINT)
      WRITE(6,'('' FIELD NO'',I4,'':''9(F10.3,2X))')
     *K,(LON(J),J=1,EW_PRINT)
      WRITE(6,'(14X,9A12)')(DASH(J),J=1,EW_PRINT)

C Print remaining rows except last
      DO 200 L=1,NS_PTS-1
      WRITE(6,'(1X,I3,'':'',F8.3,'':'',3X,9(L9,3X))')ROW,LAT,
     *(LD1(I(J)),J=1,EW_PRINT)
      DO 300 J=1,EW_PTS
      I(J)=I(J)+F_INC
300   CONTINUE
      ROW=ROW+R_INC
      LAT=LAT+R_INC*RLOOKUP(BDY,K)
200   CONTINUE

C Calculate addresses used to print values for last row
      I(1)=1+(N_ROWS-1)*N_COLS
      DO 400 J=1,EW_PTS-2
      I(J+1)=I(J)+N_COLS/(EW_PTS-1)
400   CONTINUE
      I(EW_PTS)=N_ROWS*N_COLS

C Set row pointers to last row
      LAT=RLOOKUP(BZY,K)+RLOOKUP(BDY,K)*N_ROWS
      ROW=N_ROWS

C Print last row
      WRITE(6,'(1X,I3,'':'',F8.3,'':'',3X,9(L9,3X))')ROW,LAT,
     *(LD1(I(J)),J=1,EW_PRINT)
      WRITE(6,'(14X,9A12)')(DASH(J),J=1,EW_PRINT)
      ELSE

C Print out summary of non standard fields

      EW_PRINT=MIN(EW_PTS,LOOKUP(LBLREC,K))
      WRITE(6,'(14X,9A12)')(DASH(J),J=1,EW_PRINT)
      WRITE(6,'('' FIELD NO'',I4,'':  DATA NOT ON MODEL GRID''
     *,'' SO FIRST FEW VALUES PRINTED'')')K
      WRITE(6,'(14X,9A12)')(DASH(J),J=1,EW_PRINT)
      WRITE(6,'(1X,3X,'':'',8X,'':'',3X,9(L9,3X))')
     *(LD1(J),J=1,EW_PRINT)
      WRITE(6,'(14X,9A12)')(DASH(J),J=1,EW_PRINT)

      ENDIF

      WRITE(6,'('' '')')

      RETURN
      END

