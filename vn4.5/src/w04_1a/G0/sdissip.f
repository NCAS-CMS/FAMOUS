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

      SUBROUTINE SDISSIP (F, FL, IJS, IJL, ishallo,
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

C*    *COMMON* *SHALLOW*   SHALLOW WATER TABLES.
C
     & DEPTH, DEPTHA, DEPTHD,TCGOND, TFAK, TSIHKD, INDEP,

C*    *COMMON*  *SOURCE* - SOURCE FUNCTION AND TAIL FLAG.
C
     & SL, FCONST,

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
      PARAMETER (CDIS = 4.5)
      PARAMETER (CONSD = -0.5*CDIS*ZPI**9/G**4)
      PARAMETER (CONSS = -0.5*CDIS*ZPI)

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
C**** *SDISSIP* - COMPUTATION OF DISSIPATION SOURCE FUNCTION.
C
C     S.D.HASSELMANN.
C     MODIFIED TO SHALLOW WATER : G. KOMEN , P. JANSSEN
C     OPTIMIZATION : L. ZAMBRESKY
C
C*    PURPOSE.
C     --------
C       COMPUTE DISSIPATION SOURCE FUNCTION AND STORE ADDITIVELY INTO
C       NET SOURCE FUNCTION ARRAY. ALSO COMPUTE FUNCTIONAL DERIVATIVE
C       OF DISSIPATION SOURCE FUNCTION.
C
C**   INTERFACE.
C     ----------
C
C       *CALL* *SDISSIP (F, FL, IJS, IJL)*
C          *F*   - SPECTRUM.
C          *FL*  - DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE
C          *IJS* - INDEX OF FIRST GRIDPOINT
C          *IJL* - INDEX OF LAST GRIDPOINT
C
C     METHOD.
C     -------
C
C       SEE REFERENCES.
C
C     EXTERNALS.
C     ----------
C
C       NONE.
C
C     REFERENCE.
C     ----------
C
C       G.KOMEN, S. HASSELMANN AND K. HASSELMANN, ON THE EXISTENCE
C          OF A FULLY DEVELOPED WINDSEA SPECTRUM, JGR, 1984.
C
C ----------------------------------------------------------------------
C
      DIMENSION TEMP1(NIBLO), SDS(NIBLO)
C
C ----------------------------------------------------------------------
C
      DIMENSION F(0:NIBLO,NANG,NFRE), FL(0:NIBLO,NANG,NFRE)
C ----------------------------------------------------------------------
C
C*    1. ADDING DISSIPATION AND ITS FUNCTIONAL DERIVATIVE TO NET SOURCE
C*       FUNCTION AND NET SOURCE FUNCTION DERIVATIVE.
C        --------------------------------------------------------------
C
 2000 CONTINUE
      IF (ISHALLO.EQ.1) THEN
         DO 2001 IJ=IJS,IJL
            SDS(IJ) = CONSD*EMEAN(IJ)**2*FMEAN(IJ)**9
 2001    CONTINUE

         DO 2002 M=1,NFRE
            DO 2003 IJ=IJS,IJL
               X         = (FR(M)/FMEAN(IJ))**2
               TEMP1(IJ) = SDS(IJ)*( X + X**2)
 2003       CONTINUE
            DO 2004 K=1,NANG
               DO 2005 IJ=IJS,IJL
                  SL(IJ,K,M) = SL(IJ,K,M)+TEMP1(IJ)*F(IJ,K,M)
                  FL(IJ,K,M) = FL(IJ,K,M)+TEMP1(IJ)
 2005          CONTINUE
 2004       CONTINUE
 2002    CONTINUE
      ELSE
CSHALLOW
         DO 2101 IJ=IJS,IJL
            SDS(IJ) = CONSS*FMEAN(IJ)*EMEAN(IJ)**2*AKMEAN(IJ)**4
 2101    CONTINUE

         DO 2102 M=1,NFRE
            DO 2103 IJ=IJS,IJL
               X         = TFAK(INDEP(IJ),M)/AKMEAN(IJ)
               TEMP1(IJ) = SDS(IJ)*( X + X**2)
 2103       CONTINUE
            DO 2104 K=1,NANG
               DO 2105 IJ=IJS,IJL
                  SL(IJ,K,M) = SL(IJ,K,M)+TEMP1(IJ)*F(IJ,K,M)
                  FL(IJ,K,M) = FL(IJ,K,M)+TEMP1(IJ)
 2105          CONTINUE
 2104       CONTINUE
 2102    CONTINUE
CSHALLOW
      ENDIF
C
      RETURN
      END
