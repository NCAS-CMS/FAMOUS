C ******************************COPYRIGHT******************************
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
C*LL  SUBROUTINE EX_COEF------------------------------------------------
CLL
CLL  Purpose: To calculate exchange coefficients for boundary layer
CLL           subroutine KMKH.
CLL
CLL  Suitable for single column use (via *IF definition IBM).
CLL
CLL  Version 1 written by Fiona Hewer, May 1992.
CLL
CLL  Model            Modification history:
CLL version  Date
CLL
CLL   3.4  18/10/94   *DECK inserted into UM version 3.4. S Jackson
CLL
CLL   4.3  04/02/97   Reduction in mixing length above b.l. top
CLL                   under control of logical switch L_MIXLEN &
CLL                   removal of enhanced stable momentum mixing
CLL                   under control of logical switch L_MOM.
CLL                                                      R.N.B.Smith
!LL   4.5  23/10/98   Prevent TIMER from performing barrier  P.Burton
CLL
CLL  Programming standard: Unified Model Documentation Paper No 4,
CLL                        Version 2, dated 18/1/90.
CLL
CLL  System component covered: Part of P243.
CLL
CLL  Project task:
CLL
CLL  Documentation: UMDP No.24
CLL
CLL---------------------------------------------------------------------
C*
C*L  Arguments :-
      SUBROUTINE EX_COEF (
     + P_FIELD,U_FIELD,P1,P_POINTS,BL_LEVELS,
     + CCB,CCT,NTML,L_MOM,L_MIXLEN,
     + AKH,BKH,CCA,DZL,PSTAR,RDZ,RI,TV,U_P,V_P,ZH,ZLB,Z0M,H_BLEND,
     + RHOKM,RHOKH,LTIMER
     +)
      IMPLICIT NONE
      LOGICAL LTIMER
C
      LOGICAL
     & L_MOM       ! IN Switch for convective momentum transport.
     &,L_MIXLEN    ! IN Switch for reducing the turbulent mixing
C                  !    length above the top of the boundary layer.
C
      INTEGER
     + P_FIELD     ! IN No. of P-grid points in whole field.
     +,U_FIELD     ! IN No. of U-grid points in whole field.
     +,P1          ! IN First P-grid point to be processed.
     +,P_POINTS    ! IN No. of P-grid points to be processed.
     +,BL_LEVELS   ! IN maximum number of boundary layer levels
      INTEGER
     + CCB(P_FIELD)              ! IN  Convective Cloud Base.
     +,CCT(P_FIELD)              ! IN  Convective Cloud Top.
     &,NTML(P_FIELD)             ! IN  Number of turbulently mixed
C                                !     layers.
      REAL
     + AKH(BL_LEVELS)            ! IN Hybrid "A" for layer interfaces.
C                                !    AKH(K) is value for lower boundary
C                                !    of layer K.
     +,BKH(BL_LEVELS)            ! IN Hybrid "B" for layer interfaces.
C                                !    BKH(K) is value for lower boundary
C                                !    of layer K.
     +,CCA(P_FIELD)              ! IN Convective Cloud Amount.
     +,DZL(P_FIELD,BL_LEVELS)    ! IN Layer depths (m).  DZL(,K) is the
C                                !    distance from layer boundary K-1/2
C                                !    to layer boundary K+1/2.  For K=1
C                                !    the lower boundary is the surface.
     +,PSTAR(P_FIELD)            ! IN Surface pressure (Pa).
     +,RDZ(P_FIELD,BL_LEVELS)    ! IN Reciprocal of distance between
C                                !    hybrid levels (m-1).  1/RDZ(,K) is
C                                !    the vertical distance from level
C                                !    K-1 to level K, except that for
C                                !    K=1 it is just the height of the
C                                !    lowest atmospheric level.
     +,RI(P_FIELD,2:BL_LEVELS)   ! IN Richardson number for lower
C                                !    interface of layer.
     +,TV(P_FIELD,BL_LEVELS)     ! IN Virtual temperature (K) - from
C                                !    SUBROUTINE Z.
     +,U_P(P_FIELD,BL_LEVELS)    ! IN Westerly wind component on P-grid
C                                !    (m per s).
     +,V_P(P_FIELD,BL_LEVELS)    ! IN Southerly wind component on P-grid
C                                !    (m per s).
     +,ZH(P_FIELD)               ! IN Boundary layer height (m).
     +,ZLB(P_FIELD,BL_LEVELS)    ! IN ZLB(,K) is height above surface of
C                                !    lower boundary of layer K (m).
     +,Z0M(P_FIELD)              ! IN Roughness length for momentum (m).
     +,H_BLEND(P_FIELD)          ! IN Blending height for effective
C                                !    roughness length scheme
      REAL
     + RHOKM(U_FIELD,BL_LEVELS)  ! INOUT Layer K-1 - to - layer K
C                                !       exchange coefficient for
C                                !       momentum, on UV-grid with first
C                                !       and last rows set to "missing
C                                !       data".Output as GAMMA*RHOKM*
C                                !       RDZUV for P244 (IMPL_CAL).
     +,RHOKH(P_FIELD,BL_LEVELS)  ! INOUT Layer K-1 - to - layer K
C                                !       exchange coefficient for FTL,
C                                !       on P-grid.Output as GAMMA*
C                                !       *RHOKH*RDZ for P244(IMPL_CAL)
C*
C*L---------------------------------------------------------------------
      EXTERNAL TIMER
C*
C*L---------------------------------------------------------------------
C    Local and other symbolic constants :-
C*L------------------COMDECK C_LHEAT------------------------------------
C LC IS LATENT HEAT OF CONDENSATION OF WATER AT 0DEGC
C LF IS LATENT HEAT OF FUSION AT 0DEGC
      REAL LC,LF

      PARAMETER(LC=2.501E6,
     &          LF=0.334E6)
C*----------------------------------------------------------------------
C*L------------------COMDECK C_R_CP-------------------------------------
C R IS GAS CONSTANT FOR DRY AIR
C CP IS SPECIFIC HEAT OF DRY AIR AT CONSTANT PRESSURE
C PREF IS REFERENCE SURFACE PRESSURE
      REAL R,CP,KAPPA,PREF  ! PREFTOKI

      PARAMETER(R=287.05,
     &          CP=1005.,
     &          KAPPA=R/CP,
     &          PREF=100000.)
C*----------------------------------------------------------------------

C*L------------------COMDECK C_EPSLON-----------------------------------
C EPSILON IS RATIO OF MOLECULAR WEIGHTS OF WATER AND DRY AIR
      REAL EPSILON,C_VIRTUAL

      PARAMETER(EPSILON=0.62198,
     &          C_VIRTUAL=1./EPSILON-1.)
C*----------------------------------------------------------------------

C*L-----------------COMDECK C_VKMAN-------------------------------------
C VKMAN IS VON KARMAN'S CONSTANT
      REAL VKMAN

      PARAMETER(VKMAN=0.4)
C*----------------------------------------------------------------------

      REAL EH,EM,G0,DH,DM,LAMBDA_MIN,A_LAMBDA
      PARAMETER (
     + EH=25.0                  ! Used in calc of stability function FH.
     +,EM=4.0                   ! Used in calc of stability function FM.
     +,G0=10.0                  ! Used in stability function calcs.
     +,DH=G0/EH                 ! Used in calc of stability function FH.
     +,DM=G0/EM                 ! Used in calc of stability function FM.
     +,LAMBDA_MIN=40.0          ! Minimum value of length scale LAMBDA.
     +,A_LAMBDA=2.0             ! used in calc of LAMBDA_EFF
     +)
C*
C
C  Define local storage.
C
C  Scalars.
C
      REAL
     + DVDZM   ! Modulus of wind shear gradient across lower level bdy.
     +,DZU     ! Westerly wind shear between adjacent levels.
     +,DZV     ! Southerly wind shear between adjacent levels.
     +,DVMOD2  ! Square of modulus of wind shear between adjacent levels
     +,ELH     ! Mixing length for heat and moisture at lower layer bdy.
     +,ELM     ! Mixing length for momentum at lower layer boundary.
     +,FH      ! (Value of) stability function for heat & moisture.
     +,FM      ! (Value of) stability function for momentum transport.
     +,RHO     ! Air density at lower layer boundary.
     +,RTMRI   ! Temporary in stability function calculation.
     +,VKZ     ! Temporary in calculation of ELH.
     +,WK      ! Temporary in calculation of RHO.
     +,WKM1    ! Temporary in calculation of RHO.
     +,LAMBDAM ! Asymptotic mixing length for turbulent transport
C              ! of momentum.
     +,LAMBDAH ! Asymptotic mixing length for turbulent transport
C              ! of heat/moisture.
     +,LAMBDA_EFF ! Effective mixing length used with effective
C                 ! roughness length scheme.
      INTEGER
     + I       ! Loop counter (horizontal field index).
     +,K       ! Loop counter (vertical level index).
     +,KM1     ! K-1.
C
C Layer interface K_LOG_LAYR-1/2 is the highest which requires log
C profile correction factors to the vertical finite differences.
C The value should be reassessed if the vertical resolution is changed.
C We could set K_LOG_LAYR = BL_LEVELS and thus apply the correction
C factors for all the interfaces treated by the boundary layer scheme;
C this would be desirable theoretically but expensive computationally
C because of the use of the log function.
C
      INTEGER    K_LOG_LAYR
      PARAMETER (K_LOG_LAYR=2)
C*
      IF (LTIMER) THEN
        CALL TIMER('EX_COEF   ',103)
      ENDIF
C
C-----------------------------------------------------------------------
CL 1.  Loop round "boundary" levels; calculate the stability-
CL     dependent turbulent mixing coefficients.
C-----------------------------------------------------------------------
C
      DO 2 K=2,BL_LEVELS
        KM1 = K-1
CMIC$ DO ALL VECTOR SHARED(BL_LEVELS, P_POINTS, P1, DM, DH,
CMIC$1   P_FIELD, PSTAR, TV, DZL, RDZ, ZLB, Z0M, LAMBDA_MIN, ZH,
CMIC$2   U_P, V_P, RI, K, KM1, U_FIELD, RHOKM, RHOKH, AKH, BKH, CCA,
CMIC$3   CCB, CCT, H_BLEND, NTML)
CMIC$4   PRIVATE(RHO, VKZ, ELM, ELH, DZU, DZV, DVMOD2, DVDZM,
CMIC$5   RTMRI, FM, FH, I, WKM1, WK, LAMBDAM, LAMBDAH, LAMBDA_EFF,
CMIC$6   A_LAMBDA)
CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
        DO 21 I=P1,P1+P_POINTS-1
C-----------------------------------------------------------------------
CL 2.1 Calculate asymptotic mixing lengths LAMBDAM and LAMBDAH
CL     (currently assumed equal).
C-----------------------------------------------------------------------
C
        LAMBDAM = MAX ( LAMBDA_MIN , 0.15*ZH(I) )
        LAMBDAH = LAMBDAM
        LAMBDA_EFF = MAX (LAMBDAM, A_LAMBDA*H_BLEND(I) )
        IF ( L_MIXLEN ) THEN
          IF ( K .GE. NTML(I) + 2 ) THEN
            LAMBDAM = LAMBDA_MIN
            LAMBDAH = LAMBDA_MIN
            IF (ZLB(I,K) .GT. A_LAMBDA*H_BLEND(I))
     &                                  LAMBDA_EFF = LAMBDA_MIN
          ENDIF
        ENDIF
C
C-----------------------------------------------------------------------
CL 2.2 Calculate mixing lengths ELH, ELM at layer interface K-1/2.
C-----------------------------------------------------------------------
C
C  Incorporate log profile corrections to the vertical finite
C  differences into the definitions of ELM and ELH.
C  To save computing logarithms for all K, the values of ELM and ELH
C  are unchanged for K > K_LOG_LAYR.
C
          IF (K .LE. K_LOG_LAYR) THEN
            VKZ = VKMAN / RDZ(I,K)
            ELM = VKZ / ( LOG( ( ZLB(I,K) + Z0M(I) + 0.5*DZL(I,K)   ) /
     &                          ( ZLB(I,K) + Z0M(I) - 0.5*DZL(I,KM1) ) )
     &                     + VKZ / LAMBDA_EFF )
            ELH = VKZ / ( LOG( ( ZLB(I,K) + Z0M(I) + 0.5*DZL(I,K)   ) /
     &                          ( ZLB(I,K) + Z0M(I) - 0.5*DZL(I,KM1) ) )
     &                     + VKZ / LAMBDAH )
          ELSE
            VKZ = VKMAN * ( ZLB(I,K) + Z0M(I) )
            ELM = VKZ / (1.0 + VKZ/LAMBDA_EFF )
            ELH = VKZ / (1.0 + VKZ/LAMBDAH )
          ENDIF
C
C-----------------------------------------------------------------------
CL 2.2 Calculate air density RHO = P/(R*TV) for interface between
CL     "current" and previous (lower) layers (i.e. at level K-1/2).
C-----------------------------------------------------------------------
C Repeat of KMKH calculation, could be passed in from KMKH.
          WKM1 = 0.5 * DZL(I,KM1) * RDZ(I,K)
          WK = 0.5 * DZL(I,K) * RDZ(I,K)
          RHO =               ! Calculate rho at K-1/2, from P243.111 :-
     +     ( AKH(K) + BKH(K)*PSTAR(I) )    ! Pressure at K-1/2, P243.112
     +     /                               ! divided by ...
     +     ( R *                           ! R times ...
     +     ( TV(I,KM1)*WK + TV(I,K)*WKM1 )  ! TV at K-1/2, from P243.113
     +     )
C
C-----------------------------------------------------------------------
CL 2.3 Calculate wind shear and magnitude of gradient thereof across
CL     interface K-1/2.
C-----------------------------------------------------------------------
C Repeat of KMKH calculation, could be passed in from KMKH.
          DZU = U_P(I,K) - U_P(I,KM1)
          DZV = V_P(I,K) - V_P(I,KM1)
          DVMOD2 = MAX ( 1.0E-12 , DZU*DZU + DZV*DZV )
          DVDZM = SQRT (DVMOD2) * RDZ(I,K)
C
C-----------------------------------------------------------------------
CL 2.4 Calculate (values of) stability functions FH, FM.
C-----------------------------------------------------------------------
C
          IF (RI(I,K) .GE. 0.0) THEN
            RTMRI = 0.0
            FM = 1.0 / ( 1.0 + G0*RI(I,K) )
            FH = FM
C           !-----------------------------------------------------------
C           ! If convective cloud exists in layer K allow neutral mixing
C           ! of momentum between layers K-1 and K. This is to ensure
C           ! that a reasonable amount of momentum is mixed in the
C           ! presence of convection; it is not be required when
C           ! momentum transport is included in the convection scheme.
C           !-----------------------------------------------------------
            IF ( .NOT.L_MOM .AND. (CCA(I) .GT. 0.0) .AND.
     &           (K .GE. CCB(I)) .AND. (K .LT. CCT(I)) )
     &         FM = 1.0
          ELSE
            RTMRI = (ELM/ELH) * SQRT ( -RI(I,K) )
            FM = 1.0 - ( G0*RI(I,K) / ( 1.0 + DM*RTMRI ) )
            FH = 1.0 - ( G0*RI(I,K) / ( 1.0 + DH*RTMRI ) )
          ENDIF
C
C-----------------------------------------------------------------------
CL 2.5 Calculate exchange coefficients RHO*KM, RHO*KH for interface
CL     K-1/2.
C-----------------------------------------------------------------------
C
          RHOKM(I,K) = RHO * ELM * ELM * DVDZM * FM
          RHOKH(I,K) = RHO * ELH * ELM * DVDZM * FH
   21   CONTINUE
    2 CONTINUE
      IF (LTIMER) THEN
        CALL TIMER('EX_COEF   ',104)
      ENDIF
      RETURN
      END
