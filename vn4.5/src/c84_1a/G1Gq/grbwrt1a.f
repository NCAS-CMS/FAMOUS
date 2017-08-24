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
CLL  Routine: GRBWRT------------------------------------------------
CLL
CLL  Purpose: This routine acts as an interface between the model and
CLL  GRIB format output routines.
CLL
CLL  Author:   D.M.Goddard        Date:           23 December 1993
CLL  Reviewer:                    Date of review:
CLL
CLL  Tested under compiler:   cft77
CLL  Tested under OS version: UNICOS 5.1
CLL
CLL  Code version no: 1           Date: 15 October 1993
CLL
CLL  Modification History:
CLL  3.4   11/10/94 : Correct setting of reals in lookup table
CLL                   and add return code and message to PP2GRIB call
CLL                   R A Stratton.
!    4.0   10/03/95 : Allow alternative grib packing to be used and
!                     improve error traping. R A Stratton.
!LL  4.3   06/02/97  Modify I/O calls for MPP use  P.Burton
CLL
CLL  Programming standard: UM Doc Paper 3, version 1 (15/1/90)
CLL
CLL  Logical components covered: ...
CLL
CLL  Project task: ...
CLL
CLL  External documentation: On-line UM document ??? - ??????????
CLL
CLL  -------------------------------------------------------------------
C*L  Interface and arguments: ------------------------------------------
      SUBROUTINE GRIB_FILE(LEN1_LOOKUP,LEN2_LOOKUP,LOOKUP,RLOOKUP,IENT,
     &                     FIELD,PPHORIZ_OUT,LENBUF,NUM_CRAY_WORDS,
     &                     UNITPP,IWA,GRIB_PACKING,ICODE,CMESSAGE)

      INTEGER
     &     LEN1_LOOKUP !  IN   first dimension of LOOKUP
     &    ,LEN2_LOOKUP !  IN   second dimension of LOOKUP
     &    ,LENBUF      !  IN   No of points in output field
     &    ,IENT        !  IN   level indicator for processing LOOKUP.
     &    ,IWA         !  IN   Record number
     &    ,PPHORIZ_OUT !  IN
     &    ,UNITPP      !  IN   Output PP unit number
     &    ,GRIB_PACKING !  IN  Packing profile for grib
     &    ,LEN_FIELD
     &    ,ICODE          !  OUT  Return code
     &    ,NUM_CRAY_WORDS !  OUT  Number of cray words output in grib
     &    ,LOOKUP(LEN1_LOOKUP,LEN2_LOOKUP) ! Integer lookup headers
      REAL
     &     FIELD(PPHORIZ_OUT)  ! IN   Unpacked output array
     &    ,RLOOKUP(LEN1_LOOKUP,LEN2_LOOKUP) ! REAL lookup headers
      CHARACTER
     &     CMESSAGE*(*)     ! OUT  Will contain any error messages
C
C LOCAL VARIABLES
C
      INTEGER
     &     ILABEL(45)       ! Integer part of LOOKUP for level IENT
     &    ,LEN_IO
     &    ,IX
      REAL
     &     RLABEL(19)       ! Real part of LOOKUP for level IENT
     &    ,WORK_ARRAY(LENBUF) ! GRIB packed output array
     &    ,BUFOUT(LENBUF)   ! Output PP BUFFER
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

CL
CL 1. Fill arrays ILABEL and RLABEL
CL
      DO J=1,45
        ILABEL(J)=LOOKUP(J,IENT)
      ENDDO
      DO J=1,19
        RLABEL(J)=RLOOKUP(J+45,IENT)
      ENDDO
CL
CL 2. Convert data to GRIB code
CL
      CALL PP2GRIB(FIELD,WORK_ARRAY,LENBUF,NUM_CRAY_WORDS,GRIB_PACKING,
     &             ILABEL,RLABEL,ICODE,CMESSAGE)
      IF (ICODE.NE.0) THEN
        RETURN
      ENDIF
C     WRITE(6,*) NUM_CRAY_WORDS,LENBUF
C     write(6,*) (ilabel(j),j=1,45)
C     write(6,*) (rlabel(j),j=1,19)
CL
CL 3. Put coded data into BUFOUT for output
CL
      DO I=1,NUM_CRAY_WORDS
        BUFOUT(I)=WORK_ARRAY(I)
      ENDDO
      DO I=NUM_CRAY_WORDS+1,LENBUF
        BUFOUT(I)=0.0
      ENDDO
CL
CL 4. Update lookup for this field
CL
      DO J=1,45
        LOOKUP(J,IENT)=ILABEL(J)
      ENDDO
      DO J=1,19
        RLOOKUP(J+45,IENT)=RLABEL(J)
      ENDDO
      LOOKUP(LBLREC,IENT)=NUM_CRAY_WORDS
      LOOKUP(LBEGIN,IENT)=IWA
      LOOKUP(LBNREC,IENT)=NUM_CRAY_WORDS
      LOOKUP(DATA_TYPE,IENT)=1
      LOOKUP(NADDR,IENT)=IWA
CL
CL 5. Output BUFOUT
CL
      CALL SETPOS_single(UNITPP,IWA,ICODE)
      CALL BUFFOUT_single(UNITPP,BUFOUT(1),NUM_CRAY_WORDS,LEN_IO,IX)
      RETURN
      END
