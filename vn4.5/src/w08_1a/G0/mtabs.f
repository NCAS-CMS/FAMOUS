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

      SUBROUTINE MTABS (ml,kl,
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

C*    *COMMON* *SHALLOW*   SHALLOW WATER TABLES.
C
      real DEPTH(NIBLO, NBLO)  ! water depth (metres)
      real DEPTHA, DEPTHD      ! min depth and increment for tables (m)
      real TCGOND(NDEPTH,NFRE) ! shallow water group velocity table
      real TFAK(NDEPTH,NFRE)   ! wave number table
      real TSIHKD(NDEPTH,NFRE) ! table for omega /sinh(2kd)

      integer INDEP(NIBLO)     ! depth index for gridpoint :one block

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

      iu06=6
C ----------------------------------------------------------------------
C
C**** *MTABS* - ROUTINE TO COMPUTE TABLES USED FOR SHALLOW WATER.
C
C     H.GUNTHER            ECMWF       04/04/1990
C
C*    PURPOSE.
C     -------
C
C       TO COMPUTE TABLES USED FOR SHALLOW WATER.
C
C**   INTERFACE.
C     ----------
C
C       *CALL* *MTABS (ML, KL)*
C          *ML*      - NUMBER OF FREQUENCIES.
C          *KL*      - NUMBER OF DIRECTIONS.
C
C     METHOD.
C     -------
C
C       TABLES FOR GROUP VELOCITY, WAVE NUMBER AND OMEGA/SINH(2KD)
C       ARE COMPUTED AT ALL FREQUENCIES AND FOR A DEPTH TABLE
C       OF LENGTH NDEPTH, STARTING AT DEPTHA METERS AND INCREMENTED
C        BY DEPTHD METRES.
C
C     EXTERNALS.
C     ----------
C
C       *AKI*       - FUNCTION TO COMPUTE WAVE NUMBER.
C
C     REFERENCE.
C     ----------
C
C       NONE.
C
C ----------------------------------------------------------------------
C
      DEPTHA = 5.
      DEPTHD = 1.1
C
C ----------------------------------------------------------------------
C
C*    1. GROUP VELOCITY AND WAVE NUMBER.
C        -------------------------------
C
 1000 CONTINUE
C
C*    1.1 LOOP OVER FREQUENCIES.
C         ----------------------
C
 1100 CONTINUE
      GH = G/(4.*PI)
      DO 1101 M=1,ML
         OM=ZPI*FR(M)
C
C*    1.1.1 LOOP OVER DEPTH.
C           -----------------
C
 1110 CONTINUE
         DO 1111 JD=1,NDEPTH
            AD=DEPTHA*DEPTHD**(JD-1)
            AK=AKI(OM,AD)
            TFAK(JD,M)=AK
            AKD=AK*AD
            IF(AKD.LE.10.0) THEN
               TCGOND(JD,M) = 0.5*SQRT(G*TANH(AKD)/AK)*
     1                       (1.0+2.0*AKD/SINH(2.*AKD))
               TSIHKD(JD,M) = OM/SINH(2.*AKD)
            ELSE
               TCGOND(JD,M) = GH/FR(M)
               TSIHKD(JD,M) = 0.
            ENDIF
 1111    CONTINUE
 1101 CONTINUE
C
C ----------------------------------------------------------------------
C
C*    2. PRINT TABLES.
C        -------------
C
 2000 CONTINUE
      NAN  = 10
      NSTP = NDEPTH/NAN
      NSTP = MAX(NSTP,1)
      DEPTHE = DEPTHA*DEPTHD**(NDEPTH-1)
      WRITE (IU06,'(1H1, '' SHALLOW WATER TABLES:'',/)')
      WRITE (IU06,'(''  LOGARITHMIC DEPTH FROM: DEPTHA = '',F5.1,
     1  '' TO DEPTHE  = '',F5.1, ''IN STEPS OF DEPTHD = '',F5.1)')
     2    DEPTHA, DEPTHE, DEPTHD
      WRITE (IU06,'(''  PRINTED IN STEPS OF '',I3,'' ENTRIES'',/)') NSTP
      DO 2001 JD=1,NDEPTH,NSTP
         AD=DEPTHA*DEPTHD**(JD-1)
         WRITE (IU06,'(1X,''DEPTH = '',F7.1,'' METRES '')') AD
         WRITE (IU06,'(1X,''GROUP VELOCITY IN METRES/SECOND'')')
         WRITE (IU06,'(1x,13F10.5)') (TCGOND(JD,M),M=1,ML)
         WRITE (IU06,'(1X,''WAVE NUMBER IN 1./METRES'')')
         WRITE (IU06,'(1x,13F10.5)') (TFAK(JD,M),M=1,ML)
         WRITE (IU06,'(1X,''OMEGA/SINH(2KD) IN 1./SECOND'')')
         WRITE (IU06,'(1x,13F10.5)') (TSIHKD(JD,M),M=1,ML)
 2001 CONTINUE

      RETURN
      END
C
      REAL FUNCTION AKI (OM, BETA)

C ----------------------------------------------------------------------
C
C**** *AKI* - FUNCTION TO COMPUTE WAVE NUMBER.
C
C     G. KOMEN, P. JANSSEN   KNMI        01/06/1986
C
C*    PURPOSE.
C     -------
C
C       *AKI* COMPUTES THE WAVE NUMBER AS FUNCTION OF
C             CIRCULAR FREQUENCY AND WATER DEPTH.
C
C**   INTERFACE.
C     ----------
C
C       *FUNCTION* *AKI (OM, BETA)*
C          *OM*      - CIRCULAR FREQUENCY.
C          *BETA*    - WATER DEPTH.
C
C     METHOD.
C     -------
C
C       NEWTONS METHOD TO SOLVE THE DISPERSION RELATION IN SHALLOW
C       WATER.
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
C
C*    *PARAMETER*  RELATIVE ERROR LIMIT OF NEWTON'S METHOD.
C
      PARAMETER (EBS = 0.0001)
C
C ----------------------------------------------------------------------
C
C*    1. START VALUE:  MAXIMUM FROM DEEP  AND EXTREM SHALLOW WATER
C                      WAVE NUMBER.
C        ---------------------------------------------------------
C
 1000 CONTINUE
      AKM1=OM**2/(4.*G)
      AKM2=OM/(2.*SQRT(G*BETA))
      AO=AMAX1(AKM1,AKM2)
C
C ----------------------------------------------------------------------
C
C*    2. ITERATION LOOP.
C        ---------------
C
 2000 CONTINUE
      AKP = AO
      BO = BETA*AO
      IF (BO.GT.60.) THEN
        AKI = OM**2/G
      ELSE
        TH = G*AO*TANH(BO)
        STH = SQRT(TH)
        AO = AO+(OM-STH)*STH*2./(TH/AO+G*BO/COSH(BO)**2)
        IF (ABS(AKP-AO).GT.EBS*AO) GO TO 2000
        AKI = AO
      ENDIF

      RETURN
      END
