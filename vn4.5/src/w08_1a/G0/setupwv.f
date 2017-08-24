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

      subroutine SETUPWV(ia1,ml,kl,irefra,ishallo,
C*    *PARAMETER*  FOR ARRAY DIMENSIONS.  parall
C
     & NANG, NFRE, NGX, NGY, NBLO, NIBLO, NOVER,
     & NIBLD, NBLD, NIBLC, NBLC,
C*    *COMMON* *MAP*  LON/LAT INDEX OF EACH SEA POINT.

     & IXLG, KXLT, NX, NY, IPER, icase, AMOWEP, AMOSOP, AMOEAP,
     & AMONOP, XDELLA, XDELLO,

C*    *COMMON* *GRIDPAR*  GENERAL GRID INFORMATION.
C
     & DELPHI, DELLAM, SINPH, COSPH, IGL, IJS, IJL2, IJLS, IJL, IJLT,
C
C*    *COMMON* *FREDIR* - FREQUENCY AND DIRECTION GRID.
C
     & FR, DFIM, GOM, C, DELTH, DELTR, TH, COSTH, SINTH,
C
C*    *COMMON* *INDNL* - INDICES AND WEIGHTS USED IN THE COMPUTATION
C                        OF THE NONLINEAR TRANSFER RATE.
C
     & IKP, IKP1, IKM, IKM1, K1W, K2W, K11W, K21W, AF11, FKLAP,
     & FKLAP1, FKLAM, FKLAM1, ACL1, ACL2,  CL11, CL21, DAL1, DAL2, FRH,

C*    *COMMON* *SHALLOW*   SHALLOW WATER TABLES.
C
     & DEPTH, DEPTHA, DEPTHD,TCGOND, TFAK, TSIHKD, INDEP,

C*    *COMMON* *COUPL* - PARAMETERS FOR COUPLING.
C
     & BETAMAX, ZALP, ALPHA, XKAPPA, XNLEV,
C
C*    *COMMON* *TABLE* - TABLE FOR TOTAL STRESS AND HIGH FREQ STRESS.
C
     & TAUT, DELTAUW, DELU, TAUHFT, DELUST, DELALP,

C*    *COMMON* *CURRENT* - CURRENT FIELD.
C
     & U, V,

C*    *COMMON* *REFDOT* - DEPTH AND CURRENT PART OF THETA DOT.
C
     & THDD, THDC, SIDC,

C*    *COMMON* *UBUF*  GRID POINT DEPENDENT CONSTANTS
C
     & KLAT, KLON,

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
C*    *COMMON* *SHALLOW*   SHALLOW WATER TABLES.
C     ndepth = length of shallow water tables
      integer ndepth
      PARAMETER (NDEPTH = 52)
C
C*    *COMMON* *TABLE* - TABLE FOR TOTAL STRESS AND HIGH FREQ STRESS.
C     ! table dimensions !
      INTEGER    ITAUMAX, JUMAX, IUSTAR, IALPHA
      PARAMETER (ITAUMAX=100, JUMAX=100, IUSTAR=100, IALPHA=100)
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
C*    *COMMON* *FREDIR* - FREQUENCY AND DIRECTION GRID.
C
      real FR(nfre)    ! frequencies (Hz)
      real DFIM(nfre)  ! frequency interval * direction interval
      real GOM(nfre)   ! deep water group velocity
      real C(nfre)     ! deep water phase velocity
      real DELTH       ! angular increment of spectrum (radians)
      real DELTR       ! delth times radius of earth (m)
      real TH(nang)    ! directions in radians
      real COSTH(nang), SINTH(nang)
C
C*    *COMMON* *INDNL* - INDICES AND WEIGHTS USED IN THE COMPUTATION
C                        OF THE NONLINEAR TRANSFER RATE.
C
      integer IKP(NFRE+4), IKP1(NFRE+4)
      integer IKM(NFRE+4), IKM1(NFRE+4)
! IKP: freq. index storing energy increments into bins. for wave 3
! IKM: freq. index storing energy increments into bins. for wave 4
      integer K1W(NANG,2), K2W(NANG,2)
      integer K11W(NANG,2),K21W(NANG,2)
! K1W angular index array for storing incrfements into bins wave3
! K2W angular index array for storing incrfements into bins wave4
! K?1W holds K?W(.,1)-1 and K?W(.,2)+1

      real AF11(NFRE+4) ! weight for DIA. is multiplied by freq **11
      real FKLAP(NFRE+4), FKLAP1(NFRE+4) ! weight for interpolation
      real FKLAM(NFRE+4), FKLAM1(NFRE+4) ! '+lambda' terms wave 3 / 4
      real ACL1, ACL2,  CL11, CL21 ! angular weight '1+lambda' terms
      real DAL1, DAL2              ! 1/acl1 1/acl2
      real FRH(30)                 ! tail frequency ratio **5
C
C*    *COMMON* *SHALLOW*   SHALLOW WATER TABLES.
C
      real DEPTH(NIBLO, NBLO)  ! water depth (metres)
      real DEPTHA, DEPTHD      ! min depth and increment for tables (m)
      real TCGOND(NDEPTH,NFRE) ! shallow water group velocity table
      real TFAK(NDEPTH,NFRE)   ! wave number table
      real TSIHKD(NDEPTH,NFRE) ! table for omega /sinh(2kd)

      integer INDEP(NIBLO)     ! depth index for gridpoint :one block

C*    *COMMON* *COUPL* - PARAMETERS FOR COUPLING.
C
      real BETAMAX      ! parameter for wind input
      real ZALP         ! shifts growth curve
      real ALPHA        ! charnock constant
      real XKAPPA       ! von karman constant
      real XNLEV        ! assumed height of input winds
C
C*    *COMMON* *TABLE* - TABLE FOR TOTAL STRESS AND HIGH FREQ STRESS.
C
      real TAUT(0:ITAUMAX,0:JUMAX)   ! stress table
      real DELTAUW                   ! wave stress increment
      real DELU                      ! wind increment
      real TAUHFT(0:IUSTAR,0:IALPHA) ! high freq. stress table
      real DELUST                    ! ustar increment
      real DELALP                    ! alpha increment

C*    *COMMON* *CURRENT* - CURRENT FIELD.
C
      real U(0:NIBLC,NBLC)  ! u-component of current
      real V(0:NIBLC,NBLC)  ! v component of current
C
C*    *COMMON* *REFDOT* - DEPTH AND CURRENT PART OF THETA DOT.
C
      real THDD(NIBLD,NANG,nblo) ! depth gradient part of theta dot
      real THDC(NIBLC,NANG,nblo) ! current gradient part of theta dot
c
CCNEW MH array SIDC to replace use of SL as workspace in WAM
cccc later make nangc / nfrec etc and  either set to nang or 1
      real SIDC(NIBLC,NANG,nfre,nblo)
ccc                        ! current gradient part of SIGMA dot
C
C*    *COMMON*  *SOURCE* - SOURCE FUNCTION AND TAIL FLAG.
C
      real SL(0:NIBLO,NANG,NFRE) ! total source function array
      real FCONST(NIBLO,NFRE)  ! tail flag=1/0 for prognostic/diagnostic
CCREFRA
c! for propagation with refraction
c! sl     = sigma dot term
c! fconst = sigma/sinh2kd
CREFRA
C
C*    *COMMON* *UBUF*  GRID POINT DEPENDENT CONSTANTS
C
      integer KLAT(NIBLO,2,nblo) ! index of gridpoint south and north
      integer KLON(NIBLO,2,nblo) ! index of gridpoint west and east
c                                ! land points marked by zero
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

        logical ia1(ngx,ngy)   ! in land-sea mask land=T

C ----------------------------------------------------------------------
C   based on program preproc of WAM cycle 4
C   version 1.0   M Holt 3/5/95
C
C**** *PROGRAM PREPROC* - PREPARE DATA (BUT NOT WINDS) FOR INPUT
C                         TO WAM WAVE MODELS.
C
C*    PURPOSE.
C     --------
C
C       TO ARRANGE A GRID FOR THE WAM WAVE MODEL AND COMPUTE
C       ALL FIXED MODEL PARAMETERS WHICH ARE STORED IN DIFFERENT
C       COMMON BLOCKS.
C
C     METHOD.
C     -------
C
C       [A REPRESENTATIVE TOPOGRAPHIC DATA SET ON LAT-LONG
C       COORDINATES CONTAINING THE MODEL SQUARE BOX REGION IS
C       READ IN.]
C
C       [THE MODEL REGION IS EXTRACTED AND INTERPOLATED
C       ONTO GIVEN LAT-LONG GRID INCREMENTS (SEE SUB TOPOAR).
C       THE PROGRAM CHECKS FOR A PERIODIC LATITUDE GRID. IF THE
C       GRID IS NOT PERIODIC A CLOSED BASIN IS ASSUMED.
C       THE PROGRAM DOES NOT DISTINGUISH BETWEEN DEEP AND SHALLOW
C       WATER.]
C
C   UM - grid details read in from model dump are formatted into blocks
C
C       -BLOCK STRUCTURE :
C        GRID POINTS ARE COLLECTED INTO A 1-DIMENSIONAL ARRAY,
C        BLOCKS OF MAXIMALLY NIBLO ELEMENTS,  GRID POINTS
C        (ONLY SEAPOINTS) ARE COUNTED ALONG LINES OF LATITUDES
C        FROM WEST TO EAST WORKING FROM SOUTH TO NORTH.
C        BLOCKS OVERLAP OVER TWO LATITUDE LINES,TO COMPUTE NORTH
C        -SOUTH ADVECTION TERMS.
C
C       -NESTED GRIDS: THE GRID GENERATED CAN BE A
C         - COARSE GRID WHICH MEANS OUTPUT OF SPECTRA
C                       FOR A FOLLOW UP FINE GRID RUN.
C         - FINE   GRID WHICH MEANS INPUT OF SPECTRA
C                       FROM  AN EARLIER COARSE GRID RUN.
C         - COARSE AND FINE GRID
C
C       - REFRACTION: CONTROLLED BY THE REFRACTION OPTION
C         A CURRENT FIELD IS READ, INTERPOLATED TO THE MODEL
C         GRID AND STORED IN THE GRID OUTPUT FILE.
C
C       - PARAMETERS FOR ARRAY DIMENSIONS: THE PRORAM CHECKS
C         ALL DIMENSIONS INTERNALLY. ONLY THE BLOCK LENGTH
C         (NIBLO) IS USED FOR THE SET UP OF THE GRID, ALL
C         THE OTHER PARAMETERS HAVE TO BE LARGE ENOUGH TO
C         GET A SUCCESFULL RUN OF THE JOB. AT THE END OF
C         THE OUTPUT PROTOCOLL A LIST IS PRINTED FOR THE
C         OPTIMAL SETTINGS OF THE DIMENSION.
C
C**   INTERFACE.
C     ----------
C
C       *PROGRAM* *SETUP*
C
C
C       ALL UNITS ARE DEFINE IN SECTION 1. OF THIS PROGRAM.
C
C  arrays from these WAM common blocks are filled and returned to MAIN
C
C       COMMON BLOCKS COUPLE, CURRENT, FREDIR, INDNL, GRIDPAR, MAP,
C       COUT, TABLE, AND SHALLOW ARE WRITTEN TO UNIT IU07 AND/OR IU17.
C       ALL FREQUENCY AND DIRECTION DEPENDENT ARRAYS ARE WRITTEN FROM
C       1 TO THE USED NUMBER OF FREQUENCIES AND THE USED NUMBER OF
C       DIRECTIONS.
C       OTHER ARRAYS ARE WRITTEN ACCORDING TO THEIR DIMENSIONS.
C
C     EXTERNALS.
C     ----------
C
C       *AKI*       - COMPUTES WAVE NUMBER.
C       *CHECK*     - CHECKS CONSISTENCY OF BLOCK OVERLAPS.
C       *JAFU*      - ANGULAR INDEX OF NON LINEAR INTERACTION
C       *MBLOCK*    - PREPARES ONE BLOCK
C       *MFREDIR*   - COMPUTES FREQUENCY/DIRECTION COMMON FREDIR
C       *MGRID*     - ARRANGES GRID FOR MODEL.
C       *MTABS*     - COMPUTES TABLES USED FOR SHALLOW WATER
C       *MUBUF*     - COMPUTES COMMON UBUF.
C       *NLWEIGT*   - COMPUTES NON LINEAR WEIGHTS IN COMMON INDNL
C       *STRESS*    - STRESS TABLE.
C       *TAUHF*     - HIGH FREQUENCY STRESS TABLE.
C
C     REFERENCE.
C     ----------
C
C       NONE.
C
      iu06=6
C
C ----------------------------------------------------------------------
Cxx
C*    3. INITIALISE TOTAL NUMBER OF BLOCKS,
C*       AND GRID INCREMENTS IN RADIENS AND METRES.
C        ------------------------------------------
C
 3000 CONTINUE
      IGL=0
      DELPHI =  XDELLA*CIRC/360.
      DELLAM =  XDELLO*CIRC/360.
      DO 3001 K=1,NY
         XLAT = (AMOSOP + REAL(K-1)*XDELLA)*RAD
         SINPH(K) = SIN(XLAT)
         COSPH(K) = COS(XLAT)
 3001 CONTINUE
C
C ----------------------------------------------------------------------
C
C*    4. COMPUTE GRID INDEPENDENT COMMON BLOCKS.
C        ---------------------------------------
C
 4000 CONTINUE
C
C*    4.1 COMMON FREDIR (FREQUENCY/DIRECTION CONST).
C         ------------------------------------------
C
 4100 CONTINUE

C       array FR is read from the dump
C       here call mfredir to calculate the other arrays
C
        WRITE(6,*)'calling mfredir from setupwv'
        CALL MFREDIR (fr(1),
C*    *PARAMETER*  FOR ARRAY DIMENSIONS.  parall
C
     & NANG, NFRE, NGX, NGY, NBLO, NIBLO, NOVER,
     & NIBLD, NBLD, NIBLC, NBLC,
C*    *COMMON* *FREDIR* - FREQUENCY AND DIRECTION GRID.
C
     & FR, DFIM, GOM, C, DELTH, DELTR, TH, COSTH, SINTH,
C
     & icode)
C
C
C*    4.2 COMMON INDNL (WEIGHT OF NON-LINEAR INTERACTION).
C         ------------------------------------------------
C
 4200 CONTINUE
C
      WRITE(6,*)'calling nlweigt'
      call nlweigt (ml,kl,
C*    *PARAMETER*  FOR ARRAY DIMENSIONS.  parall
C
     & NANG, NFRE, NGX, NGY, NBLO, NIBLO, NOVER,
     & NIBLD, NBLD, NIBLC, NBLC,
C*    *COMMON* *FREDIR* - FREQUENCY AND DIRECTION GRID.
C
     & FR, DFIM, GOM, C, DELTH, DELTR, TH, COSTH, SINTH,
C
C*    *COMMON* *INDNL* - INDICES AND WEIGHTS USED IN THE COMPUTATION
C                        OF THE NONLINEAR TRANSFER RATE.
C
     & IKP, IKP1, IKM, IKM1, K1W, K2W, K11W, K21W, AF11, FKLAP,
     & FKLAP1, FKLAM, FKLAM1, ACL1, ACL2,  CL11, CL21, DAL1, DAL2, FRH,

     & icode)
      WRITE(6,*)'after nlweigt'
C
C*    4.3 COMMON SHALLOW (SHALLOW WATER TABLES).
C         --------------------------------------
C
      WRITE(6,*)'calling mtabs'
      CALL MTABS (ml,kl,
C*    *PARAMETER*  FOR ARRAY DIMENSIONS.  parall
C
     & NANG, NFRE, NGX, NGY, NBLO, NIBLO, NOVER,
     & NIBLD, NBLD, NIBLC, NBLC,
C*    *COMMON* *SHALLOW*   SHALLOW WATER TABLES.
C
     & DEPTH, DEPTHA, DEPTHD,TCGOND, TFAK, TSIHKD, INDEP,

C*    *COMMON* *FREDIR* - FREQUENCY AND DIRECTION GRID.
C
     & FR, DFIM, GOM, C, DELTH, DELTR, TH, COSTH, SINTH,
C
     & icode)
C
C*    4.4 COMMON COUPLE.
C         --------------
C
      BETAMAX = 1.20
      ZALP    = 0.0110
      ALPHA   = 0.0100
      XKAPPA  = 0.41
      XNLEV   = 10.0
C
C*    4.4 COMMON TABLE (STRESS TABLES).
C         -----------------------------
C
      WRITE(6,*)'calling stress'
      CALL STRESS (
C*    *COMMON* *COUPL* - PARAMETERS FOR COUPLING.
C
     & BETAMAX, ZALP, ALPHA, XKAPPA, XNLEV,
C
C*    *COMMON* *TABLE* - TABLE FOR TOTAL STRESS AND HIGH FREQ STRESS.
C
     & TAUT, DELTAUW, DELU, TAUHFT, DELUST, DELALP,

     & icode)

      WRITE(6,*)'calling tauhf'
      CALL TAUHF (FR(ML),
C*    *PARAMETER*  FOR ARRAY DIMENSIONS.  parall
C
     & NANG, NFRE, NGX, NGY, NBLO, NIBLO, NOVER,
     & NIBLD, NBLD, NIBLC, NBLC,
C*    *COMMON* *COUPL* - PARAMETERS FOR COUPLING.
C
     & BETAMAX, ZALP, ALPHA, XKAPPA, XNLEV,
C
C*    *COMMON* *TABLE* - TABLE FOR TOTAL STRESS AND HIGH FREQ STRESS.
C
     & TAUT, DELTAUW, DELU, TAUHFT, DELUST, DELALP,

     & icode)
C
C ----------------------------------------------------------------------
C
C*    5. GENERATE OUTPUT GRID INFORMATION.
C        ---------------------------------
C
 5000 CONTINUE
C
C*    5.2 COMPUTATION OF BLOCKS.
C         ----------------------
C
 5200 CONTINUE

      WRITE(6,*)'calling mgrid'
      CALL MGRID (IA1,
C*    *PARAMETER*  FOR ARRAY DIMENSIONS.  parall
C
     & NANG, NFRE, NGX, NGY, NBLO, NIBLO, NOVER,
     & NIBLD, NBLD, NIBLC, NBLC,
C*    *COMMON* *MAP*  LON/LAT INDEX OF EACH SEA POINT.

     & IXLG, KXLT, NX, NY, IPER, icase, AMOWEP, AMOSOP, AMOEAP,
     & AMONOP, XDELLA, XDELLO,

C*    *COMMON* *GRIDPAR*  GENERAL GRID INFORMATION.
C
     & DELPHI, DELLAM, SINPH, COSPH, IGL, IJS, IJL2, IJLS, IJL, IJLT,
C
     & icode)
C
C ----------------------------------------------------------------------
C
C*    8. GENERATE AND WRITE COMMON UBUF.
C        -------------------------------
C
 8000 CONTINUE
C
         WRITE(6,*)'calling mubuf'
         CALL MUBUF (IA1,
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

C*    *COMMON* *UBUF*  GRID POINT DEPENDENT CONSTANTS
C
     & KLAT, KLON,

     & icode)


c TO ADD
c   propdot needs a current field if irefra=2
c   read current field here (steady state currents only for WAM !
c
      if(irefra.eq.2)then
       WRITE(6,*)'dummy call read currents still to add'
      endif
C
C     this subroutine call included from initmdl
C     CC call propdot moved to main after depths are set

C ----------------------------------------------------------------------
C
C*    10. CONSISTENCY CHECK OF COMPUTED BLOCK PARAMETERS AND
C*        OUTPUT OF NECESSARY DIMENSIONS.
C         --------------------------------------------------
C
 9100 CONTINUE

      WRITE(6,*)'calling subroutine check'
      CALL CHECK (IREFRA, ML, KL, IINPC,
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

      if(icode.ne.0) then
       WRITE(6,*)'calling abort in setupwv'
       WRITE(6,*)'icode ',icode,' returned from subroutine check'
       call abort
      endif

      return
      END
