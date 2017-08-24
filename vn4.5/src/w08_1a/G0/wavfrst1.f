! subroutine WAV_FOR_STEP
!
! Description:
!   called by WAV_STEP: interfaces the UM control with wave model
!   (WAM derived) subroutines
!
!
! Current Code Owner: Martin Holt
!
! History:
! Version   Date     Comment
! -------   ----     -------
! 4.1       June 1996 Original code. M Holt
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!- End of header

         SUBROUTINE WAV_FOR_STEP(ishallo, irefra,
     & energy, mdata, idelt, idelpro,
C*    *PARAMETER*  FOR ARRAY DIMENSIONS.  parall
C
     & NANG, NFRE, NGX, NGY, NBLO, NIBLO, NOVER,
     & NIBLD, NBLD, NIBLC, NBLC,
C*    *COMMON* *FREDIR* - FREQUENCY AND DIRECTION GRID.
C
     & FR, DFIM, GOM, C, DELTH, DELTR, TH, COSTH, SINTH,
C
C*    *COMMON*  *WIND* - VARIABLES USED FOR WIND COMPUTATIONS.
C
     & U10NEW, U10OLD, THWNEW, THWOLD, USNEW, USOLD, Z0NEW,
     & Z0OLD, TAUW,

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

C*    *COMMON* *INDNL* - INDICES AND WEIGHTS USED IN THE COMPUTATION
C                        OF THE NONLINEAR TRANSFER RATE.
C
     & IKP, IKP1, IKM, IKM1, K1W, K2W, K11W, K21W, AF11, FKLAP,
     & FKLAP1, FKLAM, FKLAM1, ACL1, ACL2,  CL11, CL21, DAL1, DAL2, FRH,

C*    *COMMON* *UBUF*  GRID POINT DEPENDENT CONSTANTS
C
     & KLAT, KLON,

C*    *COMMON* *GRIDPAR*  GENERAL GRID INFORMATION.
C
     & DELPHI, DELLAM, SINPH, COSPH, IGL, IJS, IJL2, IJLS, IJL, IJLT,
C
C*    *COMMON* *MAP*  LON/LAT INDEX OF EACH SEA POINT.

     & IXLG, KXLT, NX, NY, IPER, icase, AMOWEP, AMOSOP, AMOEAP,
     & AMONOP, XDELLA, XDELLO,

C*    *COMMON* *REFDOT* - DEPTH AND CURRENT PART OF THETA DOT.
C
     & THDD, THDC, SIDC,

C*    *COMMON* *CURRENT* - CURRENT FIELD.
C
     & U, V,


     & len_pd,len_sd,len_s2,len_p2,
     & icode,cmessage)

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
C*    *COMMON*  *MEANPA* - INTEGRATED PARAMETERS OF A BLOCK.
C
      real EMEAN(NIBLO)  ! total energy
      real FMEAN(NIBLO)  ! mean frequency
      real THQ(NIBLO)    ! mean wave direction (radians)
      real AKMEAN(NIBLO) ! mean wave number
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
C*    *COMMON*  *WIND* - VARIABLES USED FOR WIND COMPUTATIONS.
C
C   U10NEW etc changed to dimension (niblo,nblo) MH 9/6/95

      real U10NEW(NIBLO,nblo) ! new wind speed m/s
      real U10OLD(NIBLO,NBLO) ! intermediate storage windspeed
      real THWNEW(NIBLO,nblo) ! wind direction rads / oceanographic
      real THWOLD(NIBLO,NBLO) ! intermediate storage wind direction
      real USNEW (NIBLO,nblo) ! new friction velocity ustar
      real USOLD (NIBLO,NBLO) ! intermediate storage ustar
      real Z0NEW (NIBLO,nblo) ! new roughness length (m)
      real Z0OLD (NIBLO,NBLO) ! intermediate storage Z0
      real TAUW(NIBLO,NBLO)   ! wave stress in (m/s)**2
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
C*    *COMMON* *UBUF*  GRID POINT DEPENDENT CONSTANTS
C
      integer KLAT(NIBLO,2,nblo) ! index of gridpoint south and north
      integer KLON(NIBLO,2,nblo) ! index of gridpoint west and east
c                                ! land points marked by zero
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
C*    *COMMON* *CURRENT* - CURRENT FIELD.
C
      real U(0:NIBLC,NBLC)  ! u-component of current
      real V(0:NIBLC,NBLC)  ! v component of current
C
c       local diagnostic arrays passed through argument list
c THIS set for full grid - pass down to WAMODEL
c
       real sinp (len_sd)
       real snl (len_sd)
       real sds (len_sd)
       real sbf (len_sd)
       real stl (len_sd)
c
c       local diagnostic arrays passed through argument list
c THIS set for full grid - pass down to WAMODEL
c
       real sadv (len_pd)
c
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

      real energy(mdata,nang,nfre)

      INTEGER ICODE            ! OUT return code
      CHARACTER*80 CMESSAGE    ! OUT message accompanying return code

c     local arrays:

C*    *COMMON* *SPE3* - A BLOCK OF SPECTRA.
C
      real FL3(0:NIBLO,NANG,NFRE)
C
C*    *COMMON* *SPE1* - A BLOCK OF SPECTRA.
C
      real FL1(0:NIBLO,NANG,NFRE)
C
c       local diagnostic arrays passed through argument list
c THIS set for one block only - pass down to IMPLSCH
c
       real sin2 (len_s2)
       real snl2 (len_s2)
       real sds2 (len_s2)
       real sbf2 (len_s2)
       real stl2 (len_s2)
c
c       local diagnostic arrays passed through argument list
c THIS set for one block only - pass down to PROPAGS
c
       real sadv2 (len_p2)
c

      real temp(mdata)
      real over(nover,nang,nfre,nblo) ! array to hold overlapping
C                                     ! rows of energy at time t
C ----------------------------------------------------------------------
C
C**** *WAMODEL* - 3-G WAM MODEL - TIME INTEGRATION OF WAVE FIELDS.
C
C     S.D. HASSELMANN  MPI       1.12.85
C
C     G. KOMEN         KNMI         6.86  MODIFIED FOR SHALLOW WATER
C     P. JANSSEN                          ASPECTS.
C
C     S.D. HASSELMANN  MPI       15.2.87  MODIFIED FOR CYBER 205.
C
C     P. LIONELLO      ISDGM      6.3.87  MODIFIED TO OUTPUT SWELL.
C
C     S.D. HASSELMANN  MPI        1.6.87  ALL VERSIONS COMBINED INTO
C                                         ONE MODEL. DEEP AND SHALLOW
C                                         WATER , CRAY AND CYBER 205
C                                         VERSION.
C
C     CYCLE_2 MODICIFATIONS:
C     ----------------------
C
C     L. ZAMBRESKY     GKSS        10.87  OPTIMIZED FOR CRAY, CYBER 205
C     H. GUNTHER
C
C     A. SPEIDEL       MPI          4.88  VARIABLE DIMENSIONS, INTERNAL
C                                         CHECKS (CFL-CRITERION).
C
C     A. SPEIDEL       MPI         11.88  CHANGES FOR CRAY-2.
C
C     K. HUBBERT       POL          6.89  DEPTH AND CURRENT REFRACTION.
C                                         PRECALCULATION OF TERMS IN
C                                         *PROPDOT*.
C                                         SOLVE WAVE ACTION EQUATION
C                                         FOR CURRENT REFRACTION.
C
C     CYCLE_3 MODICIFATIONS:
C     ----------------------
C
C     R. PORTZ , S.D. HASSELMANN   MPI          1990
C
C      - RESTRUCTURE MODEL TO CALL THE ACTUAL INTEGRATION IN TIME
C        AS A SUBROUTINE: WAMODEL. A SHELL PROGRAM "WAMSHELL" READS
C        OUTPUT FROM PREPROC AND COMPUTES THE WIND ARRAYS FOR THE
C        INTEGRATION PERIOD FROM PREWIND, WHICH HAS BEEN INCORPORATED
C        AS A SUBROUTINE.
C      - ALL INTERMEDIATE AND RESTART I/O IS DONE IN THE SUBROUTINE
C        WAMODEL AND INPREST.
C      - THE COMMON BLOCK IN THE PREPROCESSOR AND MODEL ARE MADE
C        COMPATIBLE.
C      - THE COMPUTATION OF SEVERAL PARAMETERS HAS BEEN TRANSFERRED
C        FROM THE MODEL TO PREPROC.
C      - DEPTH AND CURRENT REFRACTION HAS BEEN INCORPORATED INTO THE
C        MODEL.
C      - OPEN BOUNDARIES ARE INCORPORATED IN THE MODEL.
C      - SEVERAL MINOR ERRORS HAVE BEEN REMOVED.
C      - THE BUFFERED I/O FOR THE CYBER 205 HAS BEEN CHANGED INTO A
C        BINARY READ AND WRITE.
C
C     CYCLE_4 MODICIFATIONS:
C     ----------------------
C
C     L. ZAMBRESKY   GKSS/ECMWF   6.89  ECMWF SUB VERSION
C                                       BASED ON CYCLE_2.
C
C     H. GUNTHER     GKSS/ECMWF 10.89  ECMWF SUB VERSION REORGANIZED.
C                                      - COMMON BLOCK STRUCTURE.
C                                      - BLOCKING STRUCTURE.
C                                      - TIME COUNTING.
C                                      - GRIDDED OUTPUT FIELDS.
C                                      - HEADERS ADDED TO OUTPUT FILES.
C                                      - ERRORS IN PROPAGATION CORRECTED
C
C     P.A.E.M. JANSSEN KNMI      1990  COUPLED MODEL.
C
C     H. GUNTHER     GKSS/ECMWF  8.91  LOGARITHMIC DEPTH TABLES.
C                                      MPI CYCLE_3 AND ECMWF VERSIONS
C                                      COMBINED INTO CYCLE_4.
C
CSHALLOW
C          DIFFERENCES FOR SHALLOW WATER RUNS TO DEEP WATER RUNS
C          ARE ENCLOSED IN COMMENT LINES : 'CSHALLOW'.
CSHALLOW
CNEST
C          DIFFERENCES FOR NESTED GRID RUNS TO NORMAL RUNS
C          ARE ENCLOSED IN COMMENT LINES : 'CNEST'.
CNEST
CREFRA
C          DIFFERENCES FOR REFRACTION RUNS TO NORMAL RUNS
C          ARE ENCLOSED IN COMMENT LINES : 'CREFRA'.
CREFRA
C
C*    PURPOSE.
C     --------
C
C       COMPUTATION OF THE 2-D FREQUENCY-DIRECTION WAVE SPECTRUM AT ALL
C       GRID POINTS FOR A GIVEN INITIAL SPECTRUM AND FORCING SURFACE
C       STRESS FIELD.
C
C**   INTERFACE.
C     ----------
C
C       *CALL* *WAMODEL (NADV)*
C         *NADV*     INTEGER   NUMBER OF ADVECTION ITERATIONS.
C
C       (in original WAM this is number of adsvection steps before
C         next wind input - not required by UM wave)
C
C     METHOD.
C     -------
C
C       GRID POINTS ARE LAT - LONG,VECTORIZATION IS ACHIEVED BY RUNNING
C       THROUGH THE GRID POINTS IN AN INNER LOOP ORGANIZED AS 1-D ARRAY
C       IN BLOCKS,-ALL COMPUTATIONS ARE CARRIED OUT FOR ONE BLOCK AT A
C       TIME (SEE "BLOCK STRUCTURE" BELOW)
C
C       ALL COMPONENTS OF THE SPECTRUM ARE COMPUTED PROGNOSTICALLY FROM
C       THE SPECTRAL TRANSPORT EQUATION UP TO A VARIABLE CUT-OFF
C       FREQUENCY = MAX(4*FPM,2.5*FMEAN),WHERE FPM IS THE
C       PIERSON MOSKOVITZ FREQUENCY AND FMEAN IS THE MEAN FREQUENCY,
C       BEYOND THE PROGNOSTIC CUTOFF A DIAGNOSTIC F**-5 TAIL IS ATTACHED
C       CONTINUOUSLY FOR EACH DIRECTION,
C
C       SOURCE FUNCTIONS ARE TAKEN FROM KOMEN ET AL(1984)
C
C       THE NONLINEAR TRANSFER IS PARAMETERIZED BY THE DISCRETE INTER-
C       ACTION APPROXIMATION OF HASSELMANN ET AL (1985B)
C
C       THE SOURCE FUNCTION AND THE ADVECTION TERM ARE INTEGRATED ON TWO
C       DIFFERENT TIME STEP LEVELS AND WITH DIFFERENT METHODS,-THE
C       ADVECTION TIME STEP IS A MULTIPLE OF THE SOURCE FUNCTION
C       TIME STEP.
C
C       THE SOURCE FUNCTIONS ARE INTEGRATED IMPLICITLY ACCORDING TO
C       HASSELMANN AND HASSELMANN (1985A),-THE RELEVANT FUNCTIONAL
C       DERIVATIVES OF THE INDIVIDUAL SOURCE FUNCTIONS REQUIRED FOR THE
C       SOLUTION OF THE IMPLICIT EQUATION ARE COMPUTED WITHIN THE SOURCE
C       FUNCTION SUBS,- THE TIME STEP IS TYPICALLY 20 MIN,
C
C       THE ADVECTION IS INTEGRATED BY A FIRST ORDER UPWIND SCHEME,ALSO
C       ACCORDING TO HASSELMANN AND HASSELMANN (1985A),-THE ADVECTIVE
C       TIMESTEP IS DEPENDENT ON THE FREQUENCY AND SPATIAL GRID IN
C       ACCORDANCE WITH CFL,
C
C       WINDS ARE READ IN EVERY WIND TIME STEP.IF THE WIND TIME STEP IS
C       GREATER THAN THE SOURCE TERM TIME STEP DELTWIND/DELTSOURCE STEPS
C       ARE INTEGRATED WITH CONSTANT WINDS,
C       WIND TIME STEP,PROPAGATION TIME STEP AND SOURCE TERM TIME STEP
C       SHOULD HAVE INTEGER RATIOS, THEY ARE GIVEN IN SECONDS AT
C       FULL MINUTES.
C
CNEST
C       ZERO ENERGY INFLUX IS ASSUMED AT COAST LINES. OPEN BOUNDARIES
C       ARE INCORPORATED IN THE MODEL, IF IT RUNS AS A NESTED GRID.
CNEST
C
C       BLOCK STRUCTURE (SEE PREPROC FOR DETAILS):
C       SEA POINTS ARE COLLECTED INTO A 1-DIMENSIONAL ARRAY.
C       BLOCKS OF MAXIMALLY NIBLO ELEMENTS.
C       SEA POINTS ARE COUNTED ALONG LINES OF LATITUDES FROM LEFT COAST
C       TO RIGHT COAST WORKING FROM SOUTH TO NORTH.
C       BLOCKS OVERLAP OVER TWO LATITUDE LINES,TO COMPUTE NORTH-SOUTH
C       ADVECTION TERMS, SEE ALSO COMMON GRIDPAR AND UBUF.
C
C       THE WIND FILES FOR THE BLOCKED WINDS CREATED BY PREWIND ARE
C       READ AND DELETED IN SUB IMPLSCH (IU17 AND IU18). THE FILE
c       NAMES ARE CREATED IN SUB CREWFN AND AN IMPLICIT OPEN IS USED.
C
C       THE FILE HANDLING SUBS OPENFIL, GSFILE AND CREWFN ARE COMPUTER
C       DEPENDENT AND MAY BE ADOPTED BY THE USER.
C       THE PROGRAM CLOSES AND DELETES ALL WORK FILES.
C
C       ALL PARAMETERS HAVE TO BE THE VALUES GIVEN AT THE END OF THE
C       PREPROC OUTPUT IN COLUMN 'REQUIRED'.
C
C     EXTERNALS.
C     ----------
C
C       *AIRSEA*    - SURFACE LAYER STRESS.
CREFRA
C       *DOTDC*     - READ COMMON REFDOT.
CREFRA
C       *FEMEAN*    - COMPUTATION OF MEAN FREQUENCY AT EACH GRID POINT.
C
C       *IMPLSCH*   - IMPLICIT SCHEME FOR INTEGRATION OF SOURCE
C                     FUNCTIONS IN TIME AND INPUT OF WINDS.
CREFRA
C       *INTPOL*    - MAP SPECTRUM FROM SIGMA TO OMEGA SPACE.
CREFRA
CSHALLOW
C       *SBOTTOM*   - COMPUTES BOTTOM DISSIPATION SOURCE TERM AND
C                     LINEAR CONTRIBUTION TO FUNCTIONAL MATRIX.
CSHALLOW
C       *SDISSIP*   - COMPUTATION OF DISSIPATION SOURCE FUNCTION
C                     AND LINEAR CONTRIBUTION OF DISSIPATION TO
C                     FUNCTIONAL MATRIX IN IMPLICIT SCHEME.
C       *SEMEAN*    - COMPUTATION OF TOTAL ENERGY AT EACH GRID POINT.
C
C       *SINPUT*    - COMPUTATION OF INPUT SOURCE FUNCTION, AND
C                     LINEAR CONTRIBUTION OF INPUT SOURCE FUNCTION
C                     TO FUNCTIONAL MATRIX IN IMPLICIT SCHEME.
C       *SNONLIN*   - COMPUTATION OF NONLINEAR TRANSFER RATE AND
C                     DIAGONAL LINEAR CONTRIBUTION OF NONLINEAR SOURCE
C                     FUNCTION TO FUNCTIONAL MATRIX.
C
C       *STRESSO*   - COMPUTATION OF WAVE STRESS.
C
C     REFERENCE.
C     ----------
C
C       SNYDER, R.L., F.W. DOBSON, J.A. ELLIOT, AND R.B. LONG:
C          ARRAY MEASUREMENTS OF ATMOSPHERIC PRESSURE FLUCTUATIONS
C          ABOVE SURFACE GRAVITY WAVES. J.FLUID MECH. 102, 1-59 ,1981.
C       G. KOMEN, S. HASSELMANN, K. HASSELMANN:
C          ON THE EXISTENCE OF A FULLY DEVELOPED WIND SEA SPECTRUM.
C          JPO,1984.
C       S. HASSELMANN, K. HASSELMANN, J.H. ALLENDER, T.P. BARNETT:
C          IMPROVED METHODS OF COMPUTING AND PARAMETERIZING THE
C          NONLINEAR ENERGY TRANSFER IN A GRAVITY WAVE SPECTRUM.
C          JPO, 1985.
C       S. HASSELMANN, K. HASSELMANN: A GLOBAL WAVE MODEL,
C          WAM REPORT,JUNE,30/1985.
C       P. JANSSEN, G. KOMEN: A SHALLOW WATER EXTENSION OF THE
C          3-G WAM-MODEL. WAM REPORT 1985.
C       THE WAMDI GROUP: THE WAM MODEL - A THIRD GENERATION OCEAN
C          WAVE PREDICTION MODEL. JPO, VOL. 18, NO. 12, 1988.
C       P.A.E.M JANSSEN: JPO, 1989 AND 1991.
C       K. HASSELMANN: TRANSPORT EQUATION OF FINITE DEPTH SURFACE
C          WAVE SPECTRUM IN TIME DPENDANT CURRENT AND DEPTH FIELD USING
C          NONCANONICAL SPACIAL (SPHERICAL) AND WAVE NUMBER (FRQUENCY-
C          DIRECTION) COORDINATES. WAM REPROT 1988.
C
C ----------------------------------------------------------------------
      iu06=6

c  note required to add the time controls here previously done within
c  implsch - need do over ratio idelpro to idelt.
c
c   eg idelpro=1200 idelt=1200 says n_srce-step = 1
c      idelpro=3600 idelt=1200 says n-srce-step = 3
c
       n_srce_step=int(idelpro/idelt)


C
C*    1.5  LOOP FOR BLOCKS OF LATITUDES.
C          -----------------------------
C
       if(igl.gt.1) then   ! fill array of overlap energies
        nstart_ov=1
        do ig=1,igl

         nend_blok=nstart_ov + ijlt(ig) -1
         nst_ov=nend_blok - (ijlt(ig)-ijls(ig))

         do m=1,nfre
          do k=1,nang
           ifill=1
c          fill from only the first row of each block.
           do ip=nst_ov,nst_ov+ijl(ig)-ijls(ig)
            over(ifill,k,m,ig)=energy(ip,k,m)
            ifill=ifill+1
           enddo
          enddo
         enddo
         nstart_ov=nst_ov
        enddo
       endif

c      initialise index for extracting blocked data from energy array
c
       nstart=1

       DO 1500 IG=1,IGL
C
C*    1.5.2 INPUT NEIGHBOURING GRID POINT INDICES (COMMON BLOCK UBUF).
C           ----------------------------------------------------------
CSHALLOW
C
C*    1.5.3 COMPUTE SHALLOW WATER TABLE INDICES.
C           ------------------------------------
C     calculate indep for the present block
c
            IF (ISHALLO.NE.1) THEN
               DO 1530 IJ=1,IJLT(IG)
                  XD = LOG(DEPTH(IJ,IG)/DEPTHA)/LOG(DEPTHD)+1.
                  ID = NINT(XD)
                  ID = MAX(ID,1)
                  INDEP(IJ) = MIN(ID,NDEPTH)
 1530          CONTINUE
            ENDIF
CSHALLOW
C
C
C*    1.5.4 COUPLING WITH NEIGHBOURING BLOCKS IG +- 1 AND START
C*          INPUT OF SPECTRA FOR BLOCK IG+1.
C           ----------------------------------------------------
c            here fill fl1 for this block; no need to use fl2 -
c            select data from appropriate part of array 'energy'
c            FL1 requires data from 1 to ijlt in each block.
c
c            the wam routines fillbl splitbl add on / take off
c            the overlapping rows - indexed by ijs ijl relative to
c            start of each block
c
CCC    note if not calling propags then fill fl3 here:
C
c      set index for end of present block on data grid
c
       nend=nstart + ijlt(ig)-1

       do l=1,nfre
        do k=1,nang

         do ip=1,niblo
          fl1(ip,k,l)=0.
         enddo

         ifill=1
         do ip=nstart,nend
          fl1(ifill,k,l)=energy(ip,k,l)
          ifill=ifill+1
         enddo
        enddo
       enddo

c            if block number greater than one then copy overlap of
c            first row values at time t.

       if(ig.gt.1) then
        do m=1,nfre
         do k=1,nang
          do ip=1,ijs(ig)-1
           FL1(ip,k,m)=over(ip,k,m,ig)
          enddo
         enddo
        enddo
       endif

C*    1.5.5 COMPUTATION OF PROPAGATION.
C           ---------------------------

      CALL PROPAGS(FL1, FL3, IG, irefra, ishallo, idelpro,
C*    *PARAMETER*  FOR ARRAY DIMENSIONS.  parall
C
     & NANG, NFRE, NGX, NGY, NBLO, NIBLO, NOVER,
     & NIBLD, NBLD, NIBLC, NBLC,
C*    *COMMON* *FREDIR* - FREQUENCY AND DIRECTION GRID.
C
     & FR, DFIM, GOM, C, DELTH, DELTR, TH, COSTH, SINTH,
C
C*    *COMMON* *GRIDPAR*  GENERAL GRID INFORMATION.
C
     & DELPHI, DELLAM, SINPH, COSPH, IGL, IJS, IJL2, IJLS, IJL, IJLT,
C
C*    *COMMON* *MAP*  LON/LAT INDEX OF EACH SEA POINT.

     & IXLG, KXLT, NX, NY, IPER, icase, AMOWEP, AMOSOP, AMOEAP,
     & AMONOP, XDELLA, XDELLO,

C*    *COMMON* *REFDOT* - DEPTH AND CURRENT PART OF THETA DOT.
C
     & THDD, THDC, SIDC,

C*    *COMMON* *SHALLOW*   SHALLOW WATER TABLES.
C
     & DEPTH, DEPTHA, DEPTHD,TCGOND, TFAK, TSIHKD, INDEP,

C*    *COMMON* *CURRENT* - CURRENT FIELD.
C
     & U, V,

C*    *COMMON* *UBUF*  GRID POINT DEPENDENT CONSTANTS
C
     & KLAT, KLON,

c* argument list for source term diagnostic arrays for PROPAGS
c len_p2 (source diagnostics) =1 or nang*nfre*niblo as required
c
c THIS set for one block only - pass from WAMODEL to PROPAGS
c
     & sadv2, len_p2,
c
     & icode)



C*    1.5.6 INTEGRATION OF SOURCE TERMS OVER SUB TIME STEPS BETWEEN
C*          PROPAGATION TIME STEPS.
C           -------------------------------------------------------
C
       do istep=1,n_srce_step

       CALL IMPLSCH (FL3, FL1, IJS(IG), IJL(IG),
     &    IG, IGL, ishallo,idelt,
C*    *PARAMETER*  FOR ARRAY DIMENSIONS.  parall
C
     & NANG, NFRE, NGX, NGY, NBLO, NIBLO, NOVER,
     & NIBLD, NBLD, NIBLC, NBLC,
C*    *COMMON* *FREDIR* - FREQUENCY AND DIRECTION GRID.
C
     & FR, DFIM, GOM, C, DELTH, DELTR, TH, COSTH, SINTH,
C
C*    *COMMON*  *MEANPA* - INTEGRATED PARAMETERS OF A BLOCK.
C
     & EMEAN, FMEAN, THQ, AKMEAN,

C*    *COMMON*  *SOURCE* - SOURCE FUNCTION AND TAIL FLAG.
C
     & SL, FCONST,

C*    *COMMON*  *WIND* - VARIABLES USED FOR WIND COMPUTATIONS.
C
     & U10NEW, U10OLD, THWNEW, THWOLD, USNEW, USOLD, Z0NEW,
     & Z0OLD, TAUW,

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

C*    *COMMON* *INDNL* - INDICES AND WEIGHTS USED IN THE COMPUTATION
C                        OF THE NONLINEAR TRANSFER RATE.
C
     & IKP, IKP1, IKM, IKM1, K1W, K2W, K11W, K21W, AF11, FKLAP,
     & FKLAP1, FKLAM, FKLAM1, ACL1, ACL2,  CL11, CL21, DAL1, DAL2, FRH,

c* argument list for source term diagnostic arrays for implsch
c stl to hold increments from spectral tail calculations
c other arrays source terms as standard notation
c len_s2 (source diagnostics) =1 or nang*nfre*niblo as required
c
c THIS set for one block only - pass from WAMODEL to implsch
c
     & sin2, snl2, sds2, sbf2, stl2, len_s2,
c
     & icode)


       enddo

CNEST
C
C   original code in implsch here extracted coarse mesh boundary outputs
c   and inserted fine mesh boundary inputs. This is best done at the top
c   level. code has been removed from here. ALSO will expect to use UM
c   routines that are available in preference to WAM supplied routines
C
CNEST
C
C   copy the energy for time t+dt back to main array.
c
c   Both propags and implsch only work on points ijs to ijl
c   within each block. but propags accesses rows up to IJLT
c   using indices in klat / klon.
c
c   ALSO  note that array FL3 is the output from IMPLSCH with
c   values at t+dt
c
       n11=nstart -1 +ijs(ig)
       n22=nstart -1 +ijl(ig)

       do l=1,nfre
        do k=1,nang
c
cc     pick out from point ijs in the block not from point 1
c
         ifill=0
         do ip=n11,n22
          energy(ip,k,l)=FL3(ifill+ijs(ig),k,l)
          ifill=ifill+1
         enddo
        enddo
       enddo


c
c copy blocks of diagnostics into full array
c

       if(len_s2.eq.niblo*nang*nfre.and.
     &   len_sd.eq.mdata*nang*nfre) then
c
c       set istart to account for blocks already copied
c
        istart=0
        if(ig.gt.1) then
         do ii=1,ig-1
          istart=istart + ijl(ii) - ijs(ii)+1
         enddo
        endif

        do l=1,nfre
         do m=1,nang

          nstar1=((l-1)*nang + m-1)*mdata + istart - ijs(ig) +1
          nstar2=((l-1)*nang + m-1)*niblo

          do ip=ijs(ig),ijl(ig)
           if(nstar2+ip.gt.len_s2)then
             WRITE(6,*)'error in nstar2 values:',len_s2,nstar2+ip
             icode=99
             cmessage='WAV_FOR_STEP error in nstar2'
             return
           endif
           if(nstar1+ip.gt.len_sd)then
             WRITE(6,*)'error in nstar1 values:',len_sd,nstar1+ip
             icode=99
             cmessage='WAV_FOR_STEP error in nstar1'
             return
           endif
           sinp(nstar1+ip) = sin2(nstar2+ip)
           snl(nstar1+ip) = snl2(nstar2+ip)
           sds(nstar1+ip) = sds2(nstar2+ip)
           sbf(nstar1+ip) = sbf2(nstar2+ip)
           stl(nstar1+ip) = stl2(nstar2+ip)
           sadv(nstar1+ip)= sadv2(nstar2+ip)
          enddo

         enddo
        enddo
       endif

C*    BRANCHING BACK TO 1.5 FOR NEXT BLOCK OF LATITUDES
C
C      update nstart for next block:
c
       nstart = nend - (ijlt(ig)-ijls(ig))
c

 1500 CONTINUE
C
      RETURN
      END
