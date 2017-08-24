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

      SUBROUTINE STRESSO (F, IJS, IJL, IG, IGL,
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
C*    *COMMON*  *SOURCE* - SOURCE FUNCTION AND TAIL FLAG.
C
     & SL, FCONST,

C*    *COMMON* *TABLE* - TABLE FOR TOTAL STRESS AND HIGH FREQ STRESS.
C
     & TAUT, DELTAUW, DELU, TAUHFT, DELUST, DELALP,

C*    *COMMON*  *WIND* - VARIABLES USED FOR WIND COMPUTATIONS.
C
     & U10NEW, U10OLD, THWNEW, THWOLD, USNEW, USOLD, Z0NEW,
     & Z0OLD, TAUW,

     & icode)

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
C*    *COMMON* *TABLE* - TABLE FOR TOTAL STRESS AND HIGH FREQ STRESS.
C
      real TAUT(0:ITAUMAX,0:JUMAX)   ! stress table
      real DELTAUW                   ! wave stress increment
      real DELU                      ! wind increment
      real TAUHFT(0:IUSTAR,0:IALPHA) ! high freq. stress table
      real DELUST                    ! ustar increment
      real DELALP                    ! alpha increment

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
C**** *STRESSO* - COMPUTATION OF WAVE STRESS.
C
C     H. GUNTHER      GKSS/ECMWF  NOVEMBER  1989 CODE MOVED FROM SINPUT.
C     P.A.E.M. JANSSEN      KNMI  AUGUST    1990
C
C*    PURPOSE.
C     --------
C
C       COMPUTE NORMALIZED WAVE STRESS FROM INPUT SOURCE FUNCTION
C
C**   INTERFACE.
C     ----------
C
C       *CALL* *STRESSO (F, IJS, IJL, IG, IGL)*
C          *F*   - WAVE SPECTRUM.
C          *IJS* - INDEX OF FIRST GRIDPOINT.
C          *IJL* - INDEX OF LAST GRIDPOINT.
C          *IG*  - ACTUAL BLOCK NUMBER.
C          *IGL* - NUMBER OF BLOCKS.
C
C     METHOD.
C     -------
C
C       THE INPUT SOURCE FUNCTION IS INTEGRATED OVER FREQUENCY
C       AND DIRECTIONS.
C       BECAUSE ARRAY *SL* IS USED, ONLY THE INPUT SOURCE
C       HAS TO BE STORED IN *SL* (CALL FIRST SINPUT, THEN
C       STRESSO, AND THEN THE REST OF THE SOURCE FUNCTIONS)
C
C     EXTERNALS.
C     -----------
C
C       NONE.
C
C     REFERENCE.
C     ----------
C
C       R SNYDER ET AL,1981.
C       G. KOMEN, S. HASSELMANN AND K. HASSELMANN, JPO, 1984.
C       P. JANSSEN, JPO, 1985
C
C ----------------------------------------------------------------------
C
      DIMENSION F(0:NIBLO,NANG,NFRE)
C
C ----------------------------------------------------------------------
C
      DIMENSION CONSTF(NFRE), TAUHF(NIBLO), TEMP(NIBLO),
     1          XSTRESS(NIBLO), YSTRESS(NIBLO)
C ----------------------------------------------------------------------
C
C*    1. PRECOMPUTE FREQUENCY SCALING.
C        -----------------------------
C
      CONST0  = DELTH*(ZPI)**4*FR(NFRE)**5/G**2
C
      DO 1000 M=1,NFRE
         CONSTF(M) =ZPI*XINVEPS*FR(M)*DFIM(M)
 1000 CONTINUE
C
C*    2. COMPUTE WAVE STRESS OF ACTUEL BLOCK.
C        ------------------------------------
C
C*    2.1 PRESET STRESS ARRAYS.
C         ---------------------
C
      DO 2100 IJ=IJS,IJL
         XSTRESS(IJ)=0.
         YSTRESS(IJ)=0.
 2100 CONTINUE
C
C*    2.2 INTEGRATE INPUT SOURCE FUNCTION OVER FREQUENCY AND DIRECTIONS.
C         --------------------------------------------------------------
C
      DO 2200 K=1,NANG
         DO 2210 M=1,NFRE
            CONST1=CONSTF(M)*SINTH(K)
            CONST2=CONSTF(M)*COSTH(K)
            DO 2220 IJ=IJS,IJL
               XSTRESS(IJ)=XSTRESS(IJ)+SL(IJ,K,M)*CONST1
               YSTRESS(IJ)=YSTRESS(IJ)+SL(IJ,K,M)*CONST2
 2220       CONTINUE
 2210    CONTINUE
 2200 CONTINUE
C
C*    2.3 CALCULATE HIGH-FREQUENCY CONTRIBUTION TO STRESS.
C     ----------------------------------------------------
C
      DO 2300 IJ=IJS,IJL
         TEMP(IJ) = 0.
 2300 CONTINUE
C
      DO 2310 K=1,NANG
         TKD = TH(K)
         DO 2320 IJ=IJS,IJL
            COSW     = MAX(COS(TKD-THWNEW(IJ,ig)),0.)
            TEMP(IJ) = TEMP(IJ)+F(IJ,K,NFRE)*COSW**3
 2320    CONTINUE
 2310 CONTINUE
C
      DO 2330 IJ=IJS,IJL
         UST   = USNEW(IJ,ig)
         XI    = UST / DELUST
         I     = MIN (IUSTAR-1, INT(XI))
         DELI1 = MIN (1. ,XI-FLOAT(I))
         DELI2   = 1. - DELI1
C
         XJ    = (G*Z0NEW(IJ,ig)/UST**2-ALPHA) / DELALP
         J     = MIN (IALPHA-1, INT(XJ))
         DELJ1 = MIN (1. ,XJ-FLOAT(J))
         DELJ2   = 1. - DELJ1
C
         TAU1 = (TAUHFT(I,J  )*DELI2 + TAUHFT(I+1,J  )*DELI1)*DELJ2
     1        + (TAUHFT(I,J+1)*DELI2 + TAUHFT(I+1,J+1)*DELI1)*DELJ1
C
         TAUHF(IJ) = CONST0*TEMP(IJ)*UST**2*TAU1
 2330 CONTINUE
C
      DO 2340 IJ=IJS,IJL
         XSTRESS(IJ) = XSTRESS(IJ)+TAUHF(IJ)*SIN(THWNEW(IJ,ig))
         YSTRESS(IJ) = YSTRESS(IJ)+TAUHF(IJ)*COS(THWNEW(IJ,ig))
         TAUW(IJ,IG) = SQRT(XSTRESS(IJ)**2+YSTRESS(IJ)**2)
 2340 CONTINUE

      RETURN
      END
