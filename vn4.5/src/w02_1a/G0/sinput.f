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

      SUBROUTINE SINPUT (F, FL, IJS, IJL, IG, ishallo,
C*    *PARAMETER*  FOR ARRAY DIMENSIONS.  parall
C
     & NANG, NFRE, NGX, NGY, NBLO, NIBLO, NOVER,
     & NIBLD, NBLD, NIBLC, NBLC,
C*    *COMMON* *COUPL* - PARAMETERS FOR COUPLING.
C
     & BETAMAX, ZALP, ALPHA, XKAPPA, XNLEV,
C
C*    *COMMON* *FREDIR* - FREQUENCY AND DIRECTION GRID.
C
     & FR, DFIM, GOM, C, DELTH, DELTR, TH, COSTH, SINTH,
C
C*    *COMMON*  *MEANPA* - INTEGRATED PARAMETERS OF A BLOCK.
C
     & EMEAN, FMEAN, THQ, AKMEAN,

C*    *COMMON* *SHALLOW*   SHALLOW WATER TABLES.
C
     & DEPTH, DEPTHA, DEPTHD,TCGOND, TFAK, TSIHKD, INDEP,

C*    *COMMON*  *SOURCE* - SOURCE FUNCTION AND TAIL FLAG.
C
     & SL, FCONST,

C*    *COMMON*  *WIND* - VARIABLES USED FOR WIND COMPUTATIONS.
C
     & U10NEW, U10OLD, THWNEW, THWOLD, USNEW, USOLD, Z0NEW,
     & Z0OLD, TAUW,

     & icode)

C*    *COMMON* *SHALLOW*   SHALLOW WATER TABLES.
C     ndepth = length of shallow water tables
      integer ndepth
      PARAMETER (NDEPTH = 52)
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
c       values used in airsea / wave stress routines
      PARAMETER (ROAIR = 1.225, ROWATER = 1000.)
      PARAMETER (XEPS = ROAIR/ROWATER, XINVEPS = 1./XEPS)
c

C*    *COMMON* *COUPL* - PARAMETERS FOR COUPLING.
C
      real BETAMAX      ! parameter for wind input
      real ZALP         ! shifts growth curve
      real ALPHA        ! charnock constant
      real XKAPPA       ! von karman constant
      real XNLEV        ! assumed height of input winds
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
C*    *COMMON* *SHALLOW*   SHALLOW WATER TABLES.
C
      real DEPTH(NIBLO, NBLO)  ! water depth (metres)
      real DEPTHA, DEPTHD      ! min depth and increment for tables (m)
      real TCGOND(NDEPTH,NFRE) ! shallow water group velocity table
      real TFAK(NDEPTH,NFRE)   ! wave number table
      real TSIHKD(NDEPTH,NFRE) ! table for omega /sinh(2kd)

      integer INDEP(NIBLO)     ! depth index for gridpoint :one block

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

C ----------------------------------------------------------------------
C
C**** *SINPUT* - COMPUTATION OF INPUT SOURCE FUNCTION.
C
C     P.A.E.M. JANSSEN    KNMI      AUGUST    1990
C
C     OPTIMIZED BY : H. GUENTHER
C
C*    PURPOSE.
C     ---------
C
C       COMPUTE INPUT SOURCE FUNCTION AND STORE ADDITIVELY INTO NET
C       SOURCE FUNCTION ARRAY, ALSO COMPUTE FUNCTIONAL DERIVATIVE OF
C       INPUT SOURCE FUNCTION.
C
C**   INTERFACE.
C     ----------
C
C       *CALL* *SINPUT (F, FL, IJS, IJL, IG)*
C          *F*   - SPECTRUM.
C          *FL*  - DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE.
C          *IJS* - INDEX OF FIRST GRIDPOINT.
C          *IJL* - INDEX OF LAST GRIDPOINT.
C          *IG*  - BLOCK NUMBER.
C
C     METHOD.
C     -------
C
C       SEE REFERENCE.
C
C     EXTERNALS.
C     ----------
C
C       NONE.
C
C     REFERENCE.
C     ----------
C
C       P. JANSSEN, J.P.O., 1989.
C       P. JANSSEN, J.P.O., 1991
C
C ----------------------------------------------------------------------
Ccc what to do about this ??
CDIR$ VFUNCTION EXPHF
CDIR$ VFUNCTION ALOGHF
C
C ----------------------------------------------------------------------
C
      DIMENSION F(0:NIBLO,NANG,NFRE), FL(0:NIBLO,NANG,NFRE)
C
C ----------------------------------------------------------------------
C
      DIMENSION TEMP(NIBLO,NANG), TEMP1(NIBLO,NANG),
     1          UCO(NIBLO), ZCO(NIBLO), UCN(NIBLO), ZCN(NIBLO),
     2          UFAC1(NIBLO), UFAC2(NIBLO), CM(NIBLO)
C
C ----------------------------------------------------------------------
C
C*    1. PRECALCULATED ANGULAR DEPENDENCE.
C        ---------------------------------
C
 1000 CONTINUE
      DO 1001 K=1,NANG
         TKD=TH(K)
         DO 1002 IJ=IJS,IJL
            TEMP (IJ,K) = COS(TKD-THWOLD(IJ,IG))
            TEMP1(IJ,K) = COS(TKD-THWNEW(IJ,ig))
 1002    CONTINUE
 1001 CONTINUE
C
C ----------------------------------------------------------------------
C
C*    2. LOOP OVER FREQUENCIES.
C        ----------------------
C
      CONST1  = XEPS*BETAMAX/XKAPPA**2

      DO 2001 M=1,NFRE
         FAC = ZPI*FR(M)
         CONST=FAC*CONST1
C
C*      INVERSE OF PHASE VELOCITIES.
C       ----------------------------
C
         IF (ISHALLO.EQ.1) THEN
            DO 2002 IJ=IJS,IJL
               CM(IJ) = FAC/G
 2002       CONTINUE
         ELSE
            DO 2003 IJ=IJS,IJL
               CM(IJ) = TFAK(INDEP(IJ),M)/FAC
 2003       CONTINUE
         END IF
C
C*      PRECALCULATE FREQUENCY DEPENDENCE.
C       ----------------------------------
C
         DO 2004 IJ=IJS,IJL
            UCO(IJ) = USOLD(IJ,IG)*CM(IJ) + ZALP
            ZCO(IJ) = ALOG(G*Z0OLD(IJ,IG)*CM(IJ)**2)
            UCN(IJ) = USNEW(IJ,ig)*CM(IJ) + ZALP
            ZCN(IJ) = ALOG(G*Z0NEW(IJ,ig)*CM(IJ)**2)
 2004    CONTINUE
C
C*    2.1 LOOP OVER DIRECTIONS.
C         ---------------------
C
         DO 2101 K=1,NANG
            DO 2102 IJ=IJS,IJL
               UFAC1(IJ) = 0.
               UFAC2(IJ) = 0.
 2102       CONTINUE
            DO 2103 IJ=IJS,IJL
               IF (TEMP(IJ,K).GT.0.01) THEN
                  X    = TEMP(IJ,K)*UCO(IJ)
                  ZARG = XKAPPA/X
                  ZLOG = ZCO(IJ) + ZARG
                  IF (ZLOG.LT.0.) THEN
                     UFAC1(IJ) = CONST*EXP(ZLOG)*ZLOG**4*X**2
                  ENDIF
               ENDIF
 2103       CONTINUE
C
            DO 2104 IJ=IJS,IJL
               IF (TEMP1(IJ,K).GT.0.01) THEN
                  X    = TEMP1(IJ,K)*UCN(IJ)
                  ZARG = XKAPPA/X
                  ZLOG = ZCN(IJ) + ZARG
                  IF (ZLOG.LT.0.) THEN
                     UFAC2(IJ) = CONST*EXP(ZLOG)*ZLOG**4*X**2
                  ENDIF
               ENDIF
 2104       CONTINUE
C
C*    2.2 ADDING INPUT SOURCE TERM TO NET SOURCE FUNCTION.
C         ------------------------------------------------
C
            DO 2201 IJ=IJS,IJL
               SL(IJ,K,M) = 0.5*(UFAC1(IJ)+UFAC2(IJ))*F(IJ,K,M)
               FL(IJ,K,M) = UFAC2(IJ)
 2201       CONTINUE
 2101    CONTINUE
 2001 CONTINUE

      RETURN
      END
