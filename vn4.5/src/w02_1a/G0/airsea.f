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
C
C
              SUBROUTINE AIRSEA (U10, TAUW, US, Z0, IJS, IJL,
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
C**** *AIRSEA* - DETERMINE TOTAL STRESS IN SURFACE LAYER.
C
C     P.A.E.M. JANSSEN    KNMI      AUGUST    1990
C
C*    PURPOSE.
C     --------
C
C       COMPUTE TOTAL STRESS.
C
C**   INTERFACE.
C     ----------
C
C       *CALL* *AIRSEA (U10, TAUW, US, Z0, IJS, IJL)*
C          *U10*  - INPUT BLOCK OF WINDSPEEDS U10.
C          *TAUW* - INPUT BLOCK OF WAVE STRESSES.
C          *US*   - OUTPUT BLOCK OF SURFACE STRESSES.
C          *ZO*   - OUTPUT BLOCK OF ROUGHNESS LENGTH.
C          *IJS*  - INDEX OF FIRST GRIDPOINT.
C          *IJL*  - INDEX OF LAST GRIDPOINT.
C
C     METHOD.
C     -------
C
C       USE TABLE TAUT(TAUW,U) AND LINEAR INTERPOLATION.
C
C     EXTERNALS.
C     ----------
C
C       NONE.
C
C     REFERENCE.
C     ---------
C
C       NONE.
C
C ----------------------------------------------------------------------
C
      DIMENSION U10(NIBLO), TAUW(NIBLO), US(NIBLO), Z0(NIBLO)
C
C ----------------------------------------------------------------------
C
C*    1. DETERMINE TOTAL STRESS FROM TABLE.
C        ----------------------------------
C
      if(deltauw.le.0.or.delu.le.0) then
      WRITE(6,*)'error in airsea deltauw or delu is zero: setting icode'
       icode=1
       goto 999
      endif

 1000 CONTINUE
      DO 1001 IJ=IJS,IJL
         XI      = TAUW(IJ)/DELTAUW
         I       = MIN ( ITAUMAX-1, INT(XI) )
         DELI1   = XI - FLOAT(I)
         DELI2   = 1. - DELI1
         XJ      = U10(IJ)/DELU
         J       = MIN ( JUMAX-1, INT(XJ) )
         DELJ1   = XJ - FLOAT(J)
         DELJ2   = 1. - DELJ1
         DELTOLD = (TAUT(I,J  )*DELI2 + TAUT(I+1,J  )*DELI1)*DELJ2
     1           + (TAUT(I,J+1)*DELI2 + TAUT(I+1,J+1)*DELI1)*DELJ1
         US(IJ) = SQRT(DELTOLD)
 1001 CONTINUE
C
C*    2. DETERMINE ROUGHNESS LENGTH.
C        ---------------------------
C
 2000 CONTINUE
      DO 2001 IJ=IJS,IJL
         UST     = MAX (US(IJ), 0.000001)
         X       = MIN (TAUW(IJ)/UST**2, 0.999)
         Z0(IJ)  = ALPHA*UST**2 / (G*SQRT(1.-X))
 2001 CONTINUE
C
  999 continue
      RETURN
      END
