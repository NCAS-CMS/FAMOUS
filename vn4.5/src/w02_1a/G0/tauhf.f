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

      SUBROUTINE TAUHF (FRMAX,
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
C*    *COMMON* *TABLE* - TABLE FOR TOTAL STRESS AND HIGH FREQ STRESS.
C     ! table dimensions !
      INTEGER    ITAUMAX, JUMAX, IUSTAR, IALPHA
      PARAMETER (ITAUMAX=100, JUMAX=100, IUSTAR=100, IALPHA=100)
C

      PARAMETER (JTOT = 50)

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
C**** *TAUHF* - COMPUTATION OF HIGH-FREQUENCY STRESS.
C
C     PETER A.E.M. JANSSEN    KNMI      OCTOBER 90
C
C*    PURPOSE.
C     ---------
C
C       COMPUTE HIGH-FREQUENCY WAVE STRESS
C
C**   INTERFACE.
C     ----------
C
C       *CALL* *TAUHF (FRMAX)*
C          *FRMAX - LAST MODEL FREQUENCY FR(ML).
C
C     METHOD.
C     -------
C
C       SEE REFERENCE FOR WAVE STRESS CALCULATION.
C
C     EXTERNALS.
C     ----------
C
C       NONE.
C
C     REFERENCE.
C     ----------
C
C       FOR QUASILINEAR EFFECT SEE PETER A.E.M. JANSSEN,1990.
C
C ----------------------------------------------------------------------
C*    1. PRELIMINARY CALCULATIONS.
C        -------------------------
C
       USTARM = 5.
       ALPHAM = 10.*ALPHA
       DELUST = USTARM/FLOAT(IUSTAR)
       DELALP = ALPHAM/FLOAT(IALPHA)
C
       CONST1 = BETAMAX/XKAPPA**2
       OMEGAC = ZPI*FRMAX
C
       DO 1100 L=0,IALPHA
          DO 1200 K=0,IUSTAR
             TAUHFT(K,L) = 0.
 1200     CONTINUE
 1100 CONTINUE
C
C*    2. CALCULATE HIGH-FREQUENCY CONTRIBUTION TO STRESS.
C        ------------------------------------------------
C
      X0 = 0.05
      DO 2100 L=0,IALPHA
         DO 2200 K=0,IUSTAR
            UST      = MAX(FLOAT(K)*DELUST,0.000001)
            Z0       = UST**2*(ALPHA+FLOAT(L)*DELALP)/G
            OMEGACC  = MAX(OMEGAC,X0*G/UST)
            YC       = OMEGACC*SQRT(Z0/G)
            DELY     = MAX((1.-YC)/FLOAT(JTOT),0.)
            DO 2300 J=1,JTOT
               Y        = YC+FLOAT(J-1)*DELY
               OMEGA    = Y*SQRT(G/Z0)
               CM       = G/OMEGA
               ZX       = UST/CM +ZALP
               ZARG     = MIN(XKAPPA/ZX,20.)
               ZMU      = MIN(G*Z0/CM**2*EXP(ZARG),1.)
C
               ZLOG         = MIN(ALOG(ZMU),0.)
               ZBETA        = CONST1*ZMU*ZLOG**4
               TAUHFT(K,L)  = TAUHFT(K,L)+ZBETA/Y*DELY
 2300       CONTINUE
 2200    CONTINUE
 2100 CONTINUE

      RETURN
      END
