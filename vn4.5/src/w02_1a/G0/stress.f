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

      SUBROUTINE STRESS(
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


C ----------------------------------------------------------------------
C
C**** *STRESS* - COMPUTATION OF TOTAL STRESS.
C
C     P.A.E.M. JANSSEN    KNMI      AUGUST    1990
C
C*    PURPOSE.
C     ---------
C
C       TO GENERATE STRESS TABLE TAU(TAUW,U10).
C
C**   INTERFACE.
C     ----------
C
C       *CALL* *STRESS*
C
C     METHOD.
C     -------
C
C       A STEADY STATE WIND PROFILE IS ASSUMED.
C       THE WIND STRESS IS COMPUTED USING THE ROUGHNESSLENGTH
C
C                  Z1=Z0/SQRT(1-TAUW/TAU)
C
C       WHERE Z0 IS THE CHARNOCK RELATION , TAUW IS THE WAVE-
C       INDUCED STRESS AND TAU IS THE TOTAL STRESS.
C       WE SEARCH FOR STEADY-STATE SOLUTIONS FOR WHICH TAUW/TAU < 1.
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
C ----------------------------------------------------------------------
C
      PARAMETER (XM=0.50, XNU=0.00001, G=9.806, NITER=10, EPS1=0.00001)
C
C*     VARIABLE.   TYPE.     PURPOSE.
C      ---------   -------   --------
C      *XM*        REAL      POWER OF TAUW/TAU IN ROUGHNESS LENGTH.
C      *XNU*       REAL      KINEMATIC VISCOSITY OF AIR.
C      *G*         REAL      ACCELERATION OF GRAVITY.
C      *NITER*     INTEGER   NUMBER OF ITERATIONS TO OBTAIN TOTAL STRESS
C      *EPS1*      REAL      SMALL NUMBER TO MAKE SURE THAT A SOLUTION
C                            IS OBTAINED IN ITERATION WITH TAU>TAUW.
C
C ----------------------------------------------------------------------
C
C*    1.DETERMINE TOTAL STRESS.
C       -----------------------
C
C*    1.1 INITIALISE CONSTANTS.
C         ---------------------
C
      UMAX    = 50.
      TAUWMAX = 5.
      DELU    = UMAX/FLOAT(JUMAX)
      DELTAUW = TAUWMAX/FLOAT(ITAUMAX)
C
C*    1.2 DETERMINE STRESS.
C         -----------------
C
      DO 1000 I=0,ITAUMAX
         DO 1100 J=0,JUMAX
            ZTAUW   = FLOAT(I)*DELTAUW
            UTOP    = FLOAT(J)*DELU
            CDRAG   = 0.0012875
            WCD     = SQRT(CDRAG)
            USTOLD  = UTOP*WCD
            TAUOLD  = MAX(USTOLD**2, ZTAUW+EPS1)
C
            DO 1200 ITER=1,NITER
               X      = ZTAUW/TAUOLD
               UST    = SQRT(TAUOLD)
               Z0     = ALPHA*UST**2/(G)/(1.-X)**XM
               ZNU    = 0.1*XNU/UST
               Z0     = MAX(ZNU,Z0)
               F      = UST-XKAPPA*UTOP/(ALOG(XNLEV/Z0))
               DELF   = 1.-XKAPPA*UTOP/(ALOG(XNLEV/Z0))**2*2./UST*
     *                  (1.-(XM+1)*X)/(1.-X)
               UST    = UST-F/DELF
               TAUOLD =  MAX(UST**2., ZTAUW+EPS1)
 1200       CONTINUE
            TAUT(I,J)  = TAUOLD
 1100    CONTINUE
C
C*    END DO LOOP OVER INDICES OF TAU-TABLE
C
 1000 CONTINUE

      RETURN
      END
