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

      SUBROUTINE FEMEAN (F, IJS, IJL, ishallo,
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
C**** *FEMEAN* - COMPUTATION OF MEAN FREQUENCY AT EACH GRID POINT.
C
C     S.D. HASSELMANN
C     MODIFIED : P.JANSSEN (INTEGRATION OF F**-4 TAIL)
C     OPTIMIZED BY : L. ZAMBRESKY AND H. GUENTHER
C
C
C*    PURPOSE.
C     --------
C
C       COMPUTE MEAN FREQUENCY AT EACH GRID POINT.
C
C**   INTERFACE.
C     ----------
C
C       *CALL* *FEMEAN (F, IJS, IJL)*
C              *F*   - SPECTRUM.
C              *IJS* - INDEX OF FIRST GRIDPOINT
C              *IJL* - INDEX OF LAST GRIDPOINT
C
C     METHOD.
C     -------
C
C       NONE.
C
C     EXTERNALS.
C     ----------
C
C       NONE.
C
C     REFERENCE.
C     ----------
C
C       NONE.
C
C ----------------------------------------------------------------------
C
      DIMENSION F(0:NIBLO,NANG,NFRE)
C
C ----------------------------------------------------------------------
C
      DIMENSION TEMP1(NIBLO), TEMP2(NIBLO)
C
C ----------------------------------------------------------------------
C
C*    1. INITIALISE MEAN FREQUENCY ARRAY AND TAIL FACTOR.
C        ------------------------------------------------
C
 1000 CONTINUE
      DO 1001 IJ=IJS,IJL
         FMEAN(IJ) = 0.
 1001 CONTINUE
      DELT25 = 0.20*DELTH
C
C ----------------------------------------------------------------------
C
C*    2. INTEGRATE OVER FREQUENCIES AND DIRECTIONS.
C        ------------------------------------------
C
 2000 CONTINUE

      IF (ISHALLO.EQ.1) THEN
C
C*    2.1 DEEP WATER INTEGRATION.
C         -----------------------
C
 2100 CONTINUE
         DO 2101 M=1,NFRE
            FD = DFIM(M)/FR(M)
            DO 2102 IJ=IJS,IJL
               TEMP2(IJ) = 0.
 2102       CONTINUE
            DO 2103 K=1,NANG
               DO 2104 IJ=IJS,IJL
                  TEMP2(IJ) = TEMP2(IJ)+F(IJ,K,M)
 2104          CONTINUE
 2103       CONTINUE
            DO 2105 IJ=IJS,IJL
               FMEAN(IJ) = FMEAN(IJ)+FD*TEMP2(IJ)
 2105       CONTINUE
 2101    CONTINUE
CSHALLOW
      ELSE
C
C*    2.2 SHALLOW WATER INTEGRATION.
C         --------------------------
C
 2200 CONTINUE
         DO 2201 IJ=IJS,IJL
            AKMEAN(IJ) = 0.
 2201    CONTINUE
         DO 2202 M=1,NFRE
            FD=DFIM(M)/FR(M)
            DO 2203 IJ=IJS,IJL
               TEMP1(IJ) = DFIM(M)/SQRT(TFAK(INDEP(IJ),M))
               TEMP2(IJ) = 0.
 2203       CONTINUE
            DO 2204 K=1,NANG
               DO 2205 IJ=IJS,IJL
                  TEMP2(IJ) = TEMP2(IJ)+F(IJ,K,M)
 2205          CONTINUE
 2204       CONTINUE
            DO 2206 IJ=IJS,IJL
               FMEAN(IJ) = FMEAN(IJ)+FD*TEMP2(IJ)
               AKMEAN(IJ) = AKMEAN(IJ)+TEMP1(IJ)*TEMP2(IJ)
 2206       CONTINUE
 2202    CONTINUE
C
C        ADD TAIL TO MEAN WAVENUMBER AND NORMALIZE WITH TOTAL ENERGY.
C
         DEL2 = DELT25*SQRT(G)/ZPI
         DO 2207 IJ=IJS,IJL
            AKMEAN(IJ) = AKMEAN(IJ)+DEL2*TEMP2(IJ)
            AKMEAN(IJ) = (EMEAN(IJ)/AKMEAN(IJ))**2
 2207    CONTINUE
      ENDIF
CSHALLOW
C
C*    3. ADD TAIL CORRECTION TO MEAN FREQUENCY AND
C*       NORMALIZE WITH TOTAL ENERGY.
C        ------------------------------------------
C
 3000 CONTINUE


      DO 3001 IJ=IJS,IJL
         FMEAN(IJ) = FMEAN(IJ)+DELT25*TEMP2(IJ)
         FMEAN(IJ) = EMEAN(IJ)/FMEAN(IJ)
 3001 CONTINUE

      RETURN
      END
