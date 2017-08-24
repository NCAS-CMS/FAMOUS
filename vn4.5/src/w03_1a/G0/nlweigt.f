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

      SUBROUTINE NLWEIGT (ml,kl,
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
C*    *PARAMETER*  FOR DISCRETE APPROXIMATION OF NONLINEAR TRANSFER
C
      PARAMETER (ALAMD=0.25, CON=3000., DELPHI1=-11.48, DELPHI2=33.56)
      PARAMETER (CO = 1.1)

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
C**** *NLWEIGT* - COMPUTATION OF INDEX ARRAYS AND WEIGHTS FOR THE
C                 COMPUTATION OF THE NONLINEAR TRANSFER RATE.
C
C     SUSANNE HASSELMANN JUNE 86.
C
C     H. GUNTHER   ECMWF/GKSS  DECEMBER 90 - CYCLE_4 MODIFICATIONS.
C                                            4 FREQUENCIES ADDED.
C
C*    PURPOSE.
C     --------
C
C       COMPUTATION OF PARAMETERS USED IN DISCRETE INTERACTION
C       PARAMETERIZATION OF NONLINEAR TRANSFER.
C
C**   INTERFACE.
C     ----------
C
C       *CALL* *NLWEIGT (ML, KL)*
C          *ML*     INTEGER   NUMBER OF FREQUENCIES.
C          *KL*     INTEGER   NUMBER OF DIRECTIONS.
C
C     METHOD.
C     -------
C
C       NONE.
C
C     EXTERNALS.
C     ----------
C
C       *JAFU*      - FUNCTION FOR COMPUTATION OF ANGULAR INDICES
C                     OF K(F,THET).
C
C     REFERENCE.
C     ----------
C       S. HASSELMANN AND K. HASSELMANN, JPO, 1985 B.
C
C
C ----------------------------------------------------------------------
C
C*    *PARAMETER*  FOR DISCRETE APPROXIMATION OF NONLINEAR TRANSFER
C
C*     VARIABLE.   TYPE.     PURPOSE.
C      ---------   -------   --------
C      *ALAMD*     REAL      LAMBDA
C      *CON*       REAL      WEIGHT FOR DISCRETE APPROXIMATION OF
C                            NONLINEAR TRANSFER
C      *DELPHI1*   REAL
C      *DELPHI2*   REAL
C
C ----------------------------------------------------------------------
C
C*     VARIABLE.   TYPE.     PURPOSE.
C      ---------   -------   --------
C      *CO*        REAL      FREQUENCY RATIO.
C
C ----------------------------------------------------------------------
C
      DIMENSION JA1(NANG,2), JA2(NANG,2), FRLON(2*NFRE+2)
C
      iu06=6
C ----------------------------------------------------------------------
C
C*    1. COMPUTATION FOR ANGULAR GRID.
C        -----------------------------
C
 1000 CONTINUE
C
      DELTHA = DELTH*DEG
      CL1 = DELPHI1/DELTHA
      CL2 = DELPHI2/DELTHA
C
C*    1.1 COMPUTATION OF INDICES OF ANGULAR CELL.
C         ---------------------------------------
C
      KLP1 = KL+1
      IC = 1
      DO 1001 KH=1,2
         KLH = KL
         IF (KH.EQ.2) KLH=KLP1
         DO 1002 K=1,KLH
            KS = K
            IF (KH.GT.1) KS=KLP1-K+1
            IF (KS.GT.KL) GO TO 1002
            CH = IC*CL1
            JA1(KS,KH) = JAFU(CH,K,KLP1)
            CH = IC*CL2
            JA2(KS,KH) = JAFU(CH,K,KLP1)
 1002    CONTINUE
         IC = -1
 1001 CONTINUE
C
C*    1.2 COMPUTATION OF ANGULAR WEIGHTS.
C         -------------------------------
C
      ICL1 = CL1
      CL1  = CL1-ICL1
      ICL2 = CL2
      CL2  = CL2-ICL2
      ACL1 = ABS(CL1)
      ACL2 = ABS(CL2)
      CL11 = 1.-ACL1
      CL21 = 1.-ACL2
      AL11 = (1.+ALAMD)**4
      AL12 = (1.-ALAMD)**4
      DAL1 = 1./AL11
      DAL2 = 1./AL12
C
C*    1.3 COMPUTATION OF ANGULAR INDICES.
C         -------------------------------
C
      ISG = 1
      DO 1301 KH=1,2
         CL1H = ISG*CL1
         CL2H = ISG*CL2
         DO 1302 K=1,KL
            KS = K
            IF (KH.EQ.2) KS = KL-K+2
            IF(K.EQ.1) KS = 1
            K1 = JA1(K,KH)
            K1W(KS,KH) = K1
            IF (CL1H.LT.0.) THEN
               K11 = K1-1
               IF (K11.LT.1) K11 = KL
            ELSE
              K11 = K1+1
              IF (K11.GT.KL) K11 = 1
            ENDIF
            K11W(KS,KH) = K11
            K2 = JA2(K,KH)
            K2W(KS,KH) = K2
            IF (CL2H.LT.0) THEN
               K21 = K2-1
               IF(K21.LT.1) K21 = KL
            ELSE
               K21 = K2+1
               IF (K21.GT.KL) K21 = 1
            ENDIF
            K21W(KS,KH) = K21
 1302    CONTINUE
         ISG = -1
 1301 CONTINUE
C
C*    2. COMPUTATION FOR FREQUENCY GRID.
C        -------------------------------
C
 2000 CONTINUE
C
      DO 2001 M=1,ML
         FRLON(M) = FR(M)
 2001 CONTINUE
      DO 2002 M=ML+1,2*ML+2
         FRLON(M) = CO*FRLON(M-1)
 2002 CONTINUE
      F1P1 = ALOG10(CO)
      DO 2003 M=1,ML+4
         FRG = FRLON(M)
         AF11(M) = CON * FRG**11
         FLP = FRG*(1.+ALAMD)
         FLM = FRG*(1.-ALAMD)
         IKN = IFIX(ALOG10(1.+ALAMD)/F1P1+.000001)
         IKN = M+IKN
         IKP(M) = IKN
         FKP = FRLON(IKP(M))
         IKP1(M) = IKP(M)+1
         FKLAP(M) = (FLP-FKP)/(FRLON(IKP1(M))-FKP)
         FKLAP1(M) = 1.-FKLAP(M)
         IF (FRLON(1).GE.FLM) THEN
            IKM(M) = 1
            IKM1(M) = 1
            FKLAM(M) = 0.
            FKLAM1(M) = 0.
         ELSE
            IKN = IFIX(ALOG10(1.-ALAMD)/F1P1+.0000001)
            IKN = M+IKN-1
            IF (IKN.LT.1) IKN = 1
            IKM(M) = IKN
            FKM = FRLON(IKM(M))
            IKM1(M) = IKM(M)+1
            FKLAM(M) = (FLM-FKM)/(FRLON(IKM1(M))-FKM)
            FKLAM1(M) = 1.-FKLAM(M)
         ENDIF
 2003 CONTINUE
C
C*    3. COMPUTE TAIL FREQUENCY RATIOS.
C        ------------------------------
C
 3000 CONTINUE
      IE = MIN(30,ML+3)
      DO 3001 I=1,IE
         M = ML+I-1
         FRH(I) = (FRLON(ML)/FRLON(M))**5
 3001 CONTINUE
C
C*    4. PRINTER PROTOCOL.
C        -----------------
C
 4000 CONTINUE
      WRITE(IU06,'(1H1,'' NON LINEAR INTERACTION PARAMETERS:'')')
      WRITE(IU06,'(1H0,'' COMMON INDNL: CONSTANTS'')')
      WRITE(IU06,'(1X,''    ACL1       ACL2   '',
     1             ''    CL11       CL21   '',
     2             ''    DAL1       DAL2'')')
      WRITE(IU06,'(1X,6F11.8)') ACL1, ACL2, CL11, CL21, DAL1, DAL2

      WRITE(IU06,'(1H0,'' COMMON INDNL: FREQUENCY ARRAYS'')')
      WRITE(IU06,'(1X,'' M   IKP IKP1  IKM IKM1'',
     1          ''   FKLAP       FKLAP1 '',
     2          ''   FKLAM       FKLAM1     AF11'')')
      DO 4001 M=1,ML+4
         WRITE(IU06,'(1X,I2,4I5,4F11.8,E11.3)')
     1      M, IKP(M), IKP1(M), IKM(M), IKM1(M),
     2      FKLAP(M), FKLAP1(M), FKLAM(M), FKLAM1(M), AF11(M)
 4001 CONTINUE

      WRITE(IU06,'(1H0,'' COMMON INDNL: ANGULAR ARRAYS'')')
      WRITE(IU06,'(1X,''  |--------KH = 1----------|'',
     1              ''|--------KH = 2----------|'')')
      WRITE(IU06,'(1X,'' K   K1W   K2W  K11W  K21W'',
     1              ''   K1W   K2W  K11W  K21W'')')
      DO 4002 K=1,KL
      WRITE(IU06,'(1X,I2,8I6)') K,(K1W(K,KH), K2W(K,KH), K11W(K,KH),
     2              K21W(K,KH),KH=1,2)
 4002 CONTINUE
      WRITE(IU06,'(1H0,'' COMMON INDNL: TAIL ARRAY FRH'')')
      WRITE(IU06,'(1X,8F10.7)') (FRH(M),M=1,30)
      RETURN
      END
C
      INTEGER FUNCTION JAFU (CL, J, IAN)

C ----------------------------------------------------------------------
C
C**** *JAFU* - FUNCTION TO COMPUTE THE INDEX ARRAY FOR THE
C              ANGLES OF THE INTERACTING WAVENUMBERS.
C
C     S. HASSELMANN        MPIFM        01/12/1985.
C
C*    PURPOSE.
C     --------
C
C       INDICES DEFINING BINS IN FREQUENCY AND DIRECTION PLANE INTO
C       WHICH NONLINEAR ENERGY TRANSFER INCREMENTS ARE STORED. NEEDED
C       FOR COMPUTATION OF THE NONLINEAR ENERGY TRANSFER.
C
C**   INTERFACE.
C     ----------
C
C       *FUNCTION* *JAFU (CL, J, IAN)*
C          *CL*  - WEIGHTS.
C          *J*   - INDEX IN ANGULAR ARRAY.
C          *IAN* - NUMBER OF ANGLES IN ARRAY.
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
C        S. HASSELMANN AND K. HASSELMANN,JPO, 1985 B.
C
C ----------------------------------------------------------------------
C
      IDPH = CL
      JA = J+IDPH
      IF (JA.LE.0)   JA = IAN+JA-1
      IF (JA.GE.IAN) JA = JA-IAN+1
      JAFU = JA

      RETURN
      END
