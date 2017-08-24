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

      SUBROUTINE SBOTTOM (F, FL, IJS, IJL, IG,
C*    *PARAMETER*  FOR ARRAY DIMENSIONS.  parall
C
     & NANG, NFRE, NGX, NGY, NBLO, NIBLO, NOVER,
     & NIBLD, NBLD, NIBLC, NBLC,
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
      PARAMETER (CONST = -2.0*0.038/G)

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

CSHALLOW
C ----------------------------------------------------------------------
C
C**** *SBOTTOM* - COMPUTATION OF BOTTOM FRICTION.
C
C     G.J.KOMEN AND Q.D.GAO
C     OPTIMIZED BY L.F. ZAMBRESKY
C
C*    PURPOSE.
C     --------
C
C       COMPUTATION OF BOTTOM FRICTION DISSIPATION
C
C**   INTERFACE.
C     ----------
C
C       *CALL* *SBOTTOM (F, FL, IJS, IJL, IG)*
C          *F*   - SPECTRUM.
C          *FL*  - DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE
C          *IJS* - INDEX OF FIRST GRIDPOINT
C          *IJL* - INDEX OF LAST GRIDPOINT
C          *IG*  - BLOCK NUMBER
C
C     METHOD.
C     -------
C
C       SEE REFERENCES.
C
C     REFERENCES.
C     -----------
C
C       HASSELMANN ET AL, D. HYDR. Z SUPPL A12(1973) (JONSWAP)
C       BOUWS AND KOMEN, JPO 13(1983)1653-1658
C
C ----------------------------------------------------------------------
C
      DIMENSION F(0:NIBLO,NANG,NFRE), FL(0:NIBLO,NANG,NFRE)
      DIMENSION SBO(NIBLO)
C
C ----------------------------------------------------------------------
C
      DO 1050 M=1,NFRE
         DO 1010 IJ=IJS,IJL
            ARG = 2.* DEPTH(IJ,IG)*TFAK(INDEP(IJ),M)
            ARG = MIN(ARG,50.)
            SBO(IJ) = CONST*TFAK(INDEP(IJ),M)/SINH(ARG)
 1010    CONTINUE

         DO 1040 K=1,NANG
            DO 1030 IJ=IJS,IJL
               SL(IJ,K,M) = SL(IJ,K,M)+SBO(IJ)*F(IJ,K,M)
               FL(IJ,K,M) = FL(IJ,K,M)+SBO(IJ)
 1030      CONTINUE
 1040    CONTINUE
 1050 CONTINUE

      RETURN
      END
