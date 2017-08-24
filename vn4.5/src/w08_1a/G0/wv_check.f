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

      SUBROUTINE CHECK (IREFRA, ML, KL, IINPC,
C*    *PARAMETER*  FOR ARRAY DIMENSIONS.  parall
C
     & NANG, NFRE, NGX, NGY, NBLO, NIBLO, NOVER,
     & NIBLD, NBLD, NIBLC, NBLC,
C*    *COMMON* *GRIDPAR*  GENERAL GRID INFORMATION.
C
     & DELPHI, DELLAM, SINPH, COSPH, IGL, IJS, IJL2, IJLS, IJL, IJLT,
C
C*    *COMMON* *MAP*  LON/LAT INDEX OF EACH SEA POINT.

     & IXLG, KXLT, NX, NY, IPER, icase, AMOWEP, AMOSOP, AMOEAP,
     & AMONOP, XDELLA, XDELLO,

     & icode)

C*    *PARAMETER* OF GLOBAL CONSTANTS.
C
      PARAMETER (G = 9.806, PI = 3.14159265358978, CIRC = 40000000.,
     1           ZPI = 2.*PI, RAD = PI/180., DEG = 180./PI,
     2           R = CIRC/ZPI)
C
C*     VARIABLE.   TYPE.     PURPOSE.
C      ---------   -------   --------
C      *G*         REAL      ACCELLERATION OF GRAVITY.
C      *PI*        REAL      PI.
C      *CIRC*      REAL      EARTH CIRCUMFERENCE (METRES).
C      *RAD*       REAL      PI / 180.
C      *DEG*       REAL      180. / PI.
C      *ZPI*       REAL      2. * PI.
C      *R*         REAL      EARTH RADIUS        (METRES).
C
C include parwvsh as value of ndepth is printed out
C*    *COMMON* *SHALLOW*   SHALLOW WATER TABLES.
C     ndepth = length of shallow water tables
      integer ndepth
      PARAMETER (NDEPTH = 52)
C

C*    *COMMON* *GRIDPAR*  GENERAL GRID INFORMATION.
C
      real DELPHI  ! grid increment for latitude (in metres)
      real DELLAM  ! grid increment for long. at equator (in metres)
      real SINPH(NGY), COSPH(NGY) ! sin / cos of latitude

      integer IGL   ! number of blocks
      integer IJS(NBLO)  ! index of first point of second lat
      integer IJL2(NBLO) ! index of last  point of second lat
      integer IJLS(NBLO) ! index of first point of lat before last
      integer IJL(NBLO)  ! index of last  point of lat before last
      integer IJLT(NBLO) ! total number of gridpoints in a block
C
C*    *COMMON* *MAP*  LON/LAT INDEX OF EACH SEA POINT.
C
      integer IXLG(NIBLO,NBLO) ! longitude grid index
      integer KXLT(NIBLO,NBLO) ! latitude grid index
      integer NY, NX           ! number of lats / longs in grid
      integer IPER             ! set to 1 for a periodic grid
      integer ICASE            ! set to 1 for lat-long grid

      real AMOWEP    ! most western longitude (deg)
      real AMOSOP    ! most southern latitude (deg)
      real AMOEAP    ! most eastern longitude (deg)
      real AMONOP    ! most northern latitdue (deg)
      real XDELLA    ! grid increment for latitude (deg)
      real XDELLO    ! grid increment for longitude (deg)
C
C*    *PARAMETER*  FOR ARRAY DIMENSIONS.  parall
       INTEGER
     & NANG,       ! number of direction components
     & NFRE,       ! number of frequency components
     & NGX,        ! number of cols in LS mask grid
     & NGY,        ! number of rows in LS mask grid
     & NBLO,       ! max number of blocks
     & NIBLO,      ! max number datapoints per block
     & NOVER,      ! max number datapoints in overlap row
     & NIBLD, NBLD, NIBLC, NBLC
C

C need calls to ARGWVSH because a value for ndepth is printed
C ----------------------------------------------------------------------
C
C**** *CHECK* - ROUTINE TO CHECK CONSISTENCY BETWEEN COMPUTED BLOCKS.
C
C     H.GUNTHER            ECMWF       04/04/1990
C
C*    PURPOSE.
C     -------
C
C       *CHECK* CHECKS CONSISTENCY BETWEEN BLOCK INDICES.
C
C**   INTERFACE.
C     ----------
C
C       *CALL* *CHECK (IREFRA, ML, KL, IINPC)*
C          *IREFRA*  - REFRACTION OPTION.
C          *ML*      - NUMBER OF FREQUENCIES.
C          *KL*      - NUMBER OF DIRECTIONS.
C          *IINPC*   - NUMBER INPUT POINTS FROM A PREVIOUS COARSE GRID.
C
C     METHOD.
C     -------
C
C       NONE.
C
C     EXTERNALS.
C     ----------
C
C       *ABORT*     - TERMINATES PROCESSING.
C       *OUTPP*     - WRITE OUT A GRID.
C
C     REFERENCE.
C     ----------
C
C       NONE.
C
C ----------------------------------------------------------------------
C
      CHARACTER*1 LST(NGX,NGY)
      DIMENSION GRID(NGX,NGY)
C
C*     VARIABLE.   TYPE.     PURPOSE.
C      ---------   -------   --------
C      *LSTAB*     CHARACTER  LAND SEA TABLE  L = LAND
C                                             S = SEA
C                                             + = SEA AND OUTPUT POINT.
C      *GRID*      REAL       ARRAY FOR GRIDDED PRINT OUTPUT.
C
C ----------------------------------------------------------------------
C
C*    1. COMPARE LENGTH OF OVERLAPPING LAT.
C        -----------------------------------
C
      iu06=6
 1000 CONTINUE
      DO 1001 IG=1,IGL-1
         IU1 = IJS(IG+1)-1
         IO2 = IJL(IG)-IJLS(IG)+1
         IF (IU1.NE.IO2) THEN
            WRITE (IU06,*) ' *****************************************'
            WRITE (IU06,*) ' *                                       *'
            WRITE (IU06,*) ' *      FATAL ERROR IN SUB. CHECK        *'
            WRITE (IU06,*) ' *      =========================        *'
            WRITE (IU06,*) ' *                                       *'
            WRITE (IU06,*) ' * LENGTH OF FIRST LAT. IN BLOCK IG+1    *'
            WRITE (IU06,*) ' * IS NOT EQUAL TO SECOND TO LAST OF     *'
            WRITE (IU06,*) ' * BLOCK IG                              *'
            WRITE (IU06,*) ' * BLOCK  NUMBER IS IG = ', IG
            WRITE (IU06,*) ' * LENGTH IN BLOCK IG   IS IU1 = ', IU1
            WRITE (IU06,*) ' * LENGTH IN BLOCK IG+1 IS IO2 = ', IO2
            WRITE (IU06,*) ' *                                       *'
            WRITE (IU06,*) ' *****************************************'

            icode=-1
            call abort
         ENDIF
         IU2 = IJL2(IG+1)-IJS(IG+1)+1
         IO1 = IJLT(IG)-IJL(IG)
         IF (IU2.NE.IO1) THEN
            WRITE (IU06,*) ' *****************************************'
            WRITE (IU06,*) ' *                                       *'
            WRITE (IU06,*) ' *      FATAL ERROR IN SUB. CHECK        *'
            WRITE (IU06,*) ' *      =========================        *'
            WRITE (IU06,*) ' *                                       *'
            WRITE (IU06,*) ' * LENGTH OF SECOND LAT. IN BLOCK IG+1   *'
            WRITE (IU06,*) ' * IS NOT EQUAL TO LAST OF BLOCK IG      *'
            WRITE (IU06,*) ' * BLOCK  NUMBER IS IG = ', IG
            WRITE (IU06,*) ' * LENGTH IN BLOCK IG   IS IU2 = ', IU2
            WRITE (IU06,*) ' * LENGTH IN BLOCK IG+1 IS IO1 = ', IO1
            WRITE (IU06,*) ' *                                       *'
            WRITE (IU06,*) ' *****************************************'
            icode=-1
            call abort
         ENDIF
 1001 CONTINUE
C
C ----------------------------------------------------------------------
C
C*    2. GENERATE LAND SEA TABLE FROM INDEX ARRAYS.
C        ------------------------------------------
C
 2000 CONTINUE
      DO 2001 K=1,NY
      DO 2001 I=1,NX
         LST(I,K) = 'L'
 2001 CONTINUE

      IERR = 0
      DO 2002 IG=1,IGL
         DO 2003 IJ=IJS(IG),IJL(IG)
            IF (IXLG(IJ,IG).NE.0.OR.KXLT(IJ,IG).NE.0)
     1          LST(IXLG(IJ,IG),KXLT(IJ,IG)) = 'S'
 2003    CONTINUE
 2002 CONTINUE
C
C*    2.4 PRINT LAND SEA MAP.
C         -------------------
C
 2400 CONTINUE
      ILEN = 120
      IPAGE = (NX+ILEN-1)/ILEN
      IF (IPAGE.GT.1) THEN
         LAST = (NX-ILEN*(IPAGE-1)+IPAGE-2)/(IPAGE-1)
         IF (LAST.LE.10) THEN
            ILEN = ILEN + 10
            IPAGE = (NX+ILEN-1)/ILEN
         ENDIF
      ENDIF
      DO 2401 L=1,IPAGE
         IA = (L-1)*ILEN
         IE = MIN(IA+ILEN,NX)
         IA = IA + 1
         BMOWEP = AMOWEP +REAL(IA-1)*XDELLO
         BMOEAP = AMOWEP +REAL(IE-1)*XDELLO
         WRITE (IU06,'(1H1,'' LAND SEA MAP OF FULL GRID '',
     1               ''   L = LAND  S = SEA '',
     2               ''                PAGE: '',I2)') L
         WRITE (IU06,'(2X,''LONGITUDE IS FROM '',F7.2,'' TO '',F7.2)')
     1              BMOWEP, BMOEAP
         WRITE (IU06,'(2X,''LATITUDE  IS FROM '',F7.2,'' TO '',F7.2)')
     1              AMONOP, AMOSOP
         WRITE (IU06,'(2X,130I1)') (MOD(I,10),I=IA,IE)
         DO 2402 K=NY,1,-1
            WRITE (IU06,'(1X,I1,130A1)') MOD(K,10),(LST(I,K),I=IA,IE)
 2402    CONTINUE
         WRITE (IU06,'(2X,130I1)') (MOD(I,10),I=IA,IE)
 2401 CONTINUE
C
C ----------------------------------------------------------------------
C
C*    5. OUTPUT OF OVERALL GRID INFORMATION.
C        -----------------------------------
C
 5000 CONTINUE
ccmh  WRITE (IU06,'(1H1,'' GRID SUMMERY:'')')
      WRITE (IU06,'(1H1,'' GRID SUMMARY:'')')
      WRITE (IU06,*) ' NUMBER OF BLOCKS GENERATED IS IGL ....: ', IGL
      IJFLAT = 0
      IJLLAT = 0
      IJMAX = 0
      ISEA = 0
      IPOI = 0
      DO 5001 IG=1,IGL
         IJFLAT= MAX(IJFLAT,IJS(IG)-1)
         IJLLAT= MAX(IJLLAT,IJLT(IG)-IJL(IG))
         IJMAX = MAX(IJMAX,IJLT(IG))
         IPOI  = IPOI + IJLT(IG)
         ISEA  = ISEA + IJL(IG)-IJS(IG) + 1
 5001 CONTINUE
      IOV = IPOI-ISEA
      WRITE (IU06,*) ' MAXIMUM NUMBER OF POINTS IN A BLOCK ..: ', IJMAX
      WRITE (IU06,*) ' TOTAL NUMBER OF POINT IN ALL BLOCKS ..: ', IPOI
      WRITE (IU06,*) ' TOTAL NUMBER OF SEA POINTS ...........: ', ISEA
      WRITE (IU06,*) ' TOTAL NUMBER OF POINTS IN OVERLAP.....: ', IOV
      WRITE (IU06,*) ' MAXIMUM LENGTH OF FIRST LAT OF A BLOCK: ', IJFLAT
      WRITE (IU06,*) ' MAXIMUM LENGTH OF LAST  LAT OF A BLOCK: ', IJLLAT
C
C ----------------------------------------------------------------------
      NIBLC = 1
      NIBLD = NIBLO
      NBLC = 1
      NBLD = NBLO
C
C*    6. OUTPUT OF OPTIMAL DIMENSIONS.
C        -----------------------------
C
 6000 CONTINUE
      WRITE (IU06,'(//,'' DIMENSIONS OF ARRAYS, WHICH ARE USED'',
     1             '' IN PRESET AND CHIEF '',/)')
      WRITE (IU06,'(''                                     DEFINED'',
     1           ''      USED'',''  REQUIRED'')')
      WRITE (IU06,'('' NUMBER OF DIRECTIONS        NANG '', 3I10)')
     1           NANG, KL, KL
      WRITE (IU06,'('' NUMBER OF FREQUENCIES       NFRE '', 3I10)')
     1           NFRE, ML, ML
      WRITE (IU06,'('' NUMBER LONGITUDE GRID POINTS NGX '', 3I10)')
     1           NGX, NX, NX
      WRITE (IU06,'('' NUMBER LATITUDE GRID POINTS  NGY '', 3I10)')
     1           NGY, NY, NY
      WRITE (IU06,'('' NUMBER OF BLOCKS            NBLO '', 3I10)')
     1           NBLO, IGL, NBLO
      WRITE (IU06,'('' MAXIMUM BLOCK LENGTH       NIBLO '', 3I10)')
     1           NIBLO, IJMAX, NIBLO
      WRITE (IU06,'('' NUMBER POINTS IN FIRST LAT NOVER '', 3I10)')
     1           NOVER, MAX(1,IJFLAT), MAX(1,IJFLAT)
      IF (IREFRA.EQ.0) THEN
         WRITE (IU06,'('' DEPTH REFRAC. BLOCK LENGTH NIBLD '', 3I10)')
     1           NIBLD, 1, 1
         WRITE (IU06,'('' DEPTH REFRAC. BLOCK NUMBER  NBLD '', 3I10)')
     1           NBLD, 1, 1
         WRITE (IU06,'('' CURR. REFRAC. BLOCK LENGTH NIBLC '', 3I10)')
     1           NIBLC, 1, 1
         WRITE (IU06,'('' CURR. REFRAC. BLOCK NUMBER  NBLC '', 3I10)')
     1           NBLC, 1, 1
      ELSE IF (IREFRA.EQ.1) THEN
         WRITE (IU06,'('' DEPTH REFRAC. BLOCK LENGTH NIBLD '', 3I10)')
     1           NIBLD, IJMAX, NIBLO
         WRITE (IU06,'('' DEPTH REFRAC. BLOCK NUMBER  NBLD '', 3I10)')
     1           NBLD, IGL, NBLO
         WRITE (IU06,'('' CURR. REFRAC. BLOCK LENGTH NIBLC '', 3I10)')
     1           NIBLC, 1, 1
         WRITE (IU06,'('' CURR. REFRAC. BLOCK NUMBER  NBLC '', 3I10)')
     1           NBLC, 1, 1
      ELSE IF (IREFRA.EQ.2) THEN
         WRITE (IU06,'('' DEPTH REFRAC. BLOCK LENGTH NIBLD '', 3I10)')
     1           NIBLD, IJMAX, NIBLO
         WRITE (IU06,'('' DEPTH REFRAC. BLOCK NUMBER  NBLD '', 3I10)')
     1           NBLD, IGL, NBLO
         WRITE (IU06,'('' CURR. REFRAC. BLOCK LENGTH NIBLC '', 3I10)')
     1           NIBLC, IJMAX, NIBLO
         WRITE (IU06,'('' CURR. REFRAC. BLOCK NUMBER  NBLC '', 3I10)')
     1           NBLC, IGL, NBLO
      ENDIF

      WRITE (IU06,'('' SHALLOW WATER TABLE LEN.  NDEPTH '', 3I10)')
     1           NDEPTH, NDEPTH, NDEPTH

      WRITE (IU06,'(/,'' THE DIMENSIONS IN PRESET AND CHIEF HAVE TO '',
     1             '' BE THE VALUES IN COLUMN - REQUIRED - '')')
      WRITE (IU06,'(  '' IF YOU WANT TO USE THE OPTIMAL DIMENSION'',
     1             '' LENGTH IN THE WAMODEL, THEN  '',/,
     2             '' RERUN PREPROC WITH THE DIMENSION'',
     3             '' GIVEN AS -USED-'')')

      RETURN
      END
