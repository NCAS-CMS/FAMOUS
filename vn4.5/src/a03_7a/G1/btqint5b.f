C *****************************COPYRIGHT******************************
C (c) CROWN COPYRIGHT 1997, METEOROLOGICAL OFFICE, All Rights Reserved.
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
! SUBROUTINE BTQ_INT
!!!  Purpose: To interpolate buoyancy parameters BT and BQ from full
!!!  levels to half levels
!!!
!!! ORIGINAL PROGRAMMER: J James
!!! CURRENT CODE OWNER:  RNB Smith
!!!
!!! HISTORY:
!!! DATE   VERSION   COMMENT
!!! ----   -------   -------
!!!
!!! 30/1/97  4.3     New Deck.  S.Jackson
!!! 8/9/97  4.4      L_BL_LSPICE specifies mixed phase precipitation
!!!                  scheme.                       D.Wilson
!!!
!!! CODE DESCRIPTION:
!!!   LANGUAGE: FORTRAN 77 + CRAY EXTENSIONS
!!!   THIS CODE IS WRITTEN TO UMDP3 PROGRAMMING STANDARDS.
!!!
!!! SYSTEM COMPONENT COVERED: P24
!!! SYSTEM TASK:
!!!---------------------------------------------------------------------
      SUBROUTINE BTQ_INT (
     & P_FIELD,P1
     &,P_POINTS,BL_LEVELS
     &,BQ,BT,BF,DZL
     &,RDZ
     &,QW,QCF,TL
     &,L_BL_LSPICE,LTIMER
     &  )

      IMPLICIT NONE

! ARGUMENTS WITH INTENT IN. IE: INPUT VARIABLES.

      LOGICAL LTIMER          ! IN Flag for TIMER diagnostics
     &,       L_BL_LSPICE             !IN
!                              TRUE  Use scientific treatment of mixed
!                                    phase precip scheme.
!                              FALSE Do not use mixed phase precip
!                                    considerations

      INTEGER
     & P_FIELD                ! IN No. of P-grid points in whole field.
     &,P1                     ! IN First P-grid point to be processed.
     &,P_POINTS               ! IN No. of P-grid points to be processed.
     &,BL_LEVELS              ! IN No. of atmospheric levels for which
!                                boundary layer fluxes are calculated.
!                                Assumed ! <=30 for dimensioning GAMMA()
!                                in common deck C_GAMMA
      REAL
     & DZL(P_FIELD,BL_LEVELS) ! IN Layer depths (m).  DZL(,K) is the
!                                  distance from layer boundary K-1/2
!                                  to layer boundary K+1/2.  For K=1
!                                  the lower boundary is the surface.
     &,RDZ(P_FIELD,BL_LEVELS) ! IN Reciprocal of distance between
!                                  hybrid levels (m-1).  1/RDZ(,K) is
!                                  the vertical distance from level
!                                  K-1 to level K, except that for
!                                  K=1 it is just the height of the
!                                  lowest atmospheric level.

      REAL  ! INOUT arrays,
     & QW(P_FIELD,BL_LEVELS)  ! INOUT Total water content (kg/kg air).
     &,QCF(P_FIELD,BL_LEVELS)  ! INOUT Ice water content (kg/kg air).
     &,TL(P_FIELD,BL_LEVELS)  ! INOUT Liquid/frozen water temperature
!                                     (K).
     &,BQ(P_FIELD,BL_LEVELS)  ! INOUT A buoyancy parameter
!                                (beta q tilde).
     &,BT(P_FIELD,BL_LEVELS)  ! INOUT A buoyancy parameter
!                                (beta T tilde).
     &,BF(P_FIELD,BL_LEVELS)  ! INOUT A buoyancy parameter
!                                (beta F tilde).

!-----------------------------------------------------------------------
!    External references :-
      EXTERNAL TIMER

!-----------------------------------------------------------------------
!    Local and other symbolic constants :-
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

C*L------------------COMDECK C_G----------------------------------------
C G IS MEAN ACCEL DUE TO GRAVITY AT EARTH'S SURFACE
      REAL G

      PARAMETER(G=9.80665)
C*----------------------------------------------------------------------
C*L------------------COMDECK C_O_DG_C-----------------------------------
C ZERODEGC IS CONVERSION BETWEEN DEGREES CELSIUS AND KELVIN
C TFS IS TEMPERATURE AT WHICH SEA WATER FREEZES
C TM IS TEMPERATURE AT WHICH FRESH WATER FREEZES AND ICE MELTS
      REAL ZERODEGC,TFS,TM

      PARAMETER(ZERODEGC=273.15,
     &          TFS=271.35,
     &          TM=273.15)
C*----------------------------------------------------------------------

C*L------------------COMDECK C_EPSLON-----------------------------------
C EPSILON IS RATIO OF MOLECULAR WEIGHTS OF WATER AND DRY AIR
      REAL EPSILON,C_VIRTUAL

      PARAMETER(EPSILON=0.62198,
     &          C_VIRTUAL=1./EPSILON-1.)
C*----------------------------------------------------------------------


      REAL ETAR,GRCP,LCRCP,LFRCP,LS,LSRCP
      PARAMETER (
     & ETAR=1.0/(1.0-EPSILON)   ! Used in buoyancy parameter BETAC.
     &,GRCP=G/CP                ! Used in DZTL, FTL calculations.
     &,LCRCP=LC/CP              ! Latent heat of condensation / CP.
     &,LFRCP=LF/CP              ! Latent heat of fusion / CP.
     &,LS=LC+LF                 ! Latent heat of sublimation.
     &,LSRCP=LS/CP              ! Latent heat of sublimation / CP.
     &)


!  Define local storage.

!  (b) Scalars.

      REAL
     & BQM     ! Temporary in calculation of Richardson number RI(,K)
     &,BTM     ! Temporary in calculation of Richardson number RI(,K)
     &,BFM     ! Temporary in calculation of Richardson number RI(,K)
     &,DZB     ! Temporary in calculation of Richardson number RI(,K).
     &,DZQW    ! Difference of QW between "current" and "lower" levels.
     &,DZTL    ! Liquid/ice static energy difference between adjacent
!                levels.
     &,DZQI    ! Difference of QCF between "current" and "lower" levels.
     &,WK      ! Temporary in calculation of RHO.
     &,WKM1    ! Temporary in calculation of RHO.

      INTEGER
     & I       ! Loop counter (horizontal field index).
     &,K       ! Loop counter (vertical level index).
!              ! mixing layer; set to BL_LEVELS-1.

      IF (LTIMER) THEN
        CALL TIMER('BTQ_INT ',3)
      ENDIF
!-----------------------------------------------------------------------
!! 1.  Loop round "boundary" levels.
!-----------------------------------------------------------------------

      DO K=2,BL_LEVELS
        DO I=P1,P1+P_POINTS-1

!-----------------------------------------------------------------------
!! 2 Calculate wind shear and other "jumps" between "current" and
!!   previous (lower) level.
!-----------------------------------------------------------------------

          DZTL = TL(I,K) - TL(I,K-1) + GRCP/RDZ(I,K)   ! Used in P243.C2
          DZQW = QW(I,K) - QW(I,K-1)                   ! Used in P243.C2
          DZQI = QCF(I,K) - QCF(I,K-1)

!-----------------------------------------------------------------------
!! 3. Calculate buoyancy parameters BT, BQ, DZB at interface between
!!     "current" and previous levels (i.e. at level K-1/2, if current
!!     level is level K).
!-----------------------------------------------------------------------

          WKM1 = 0.5 * DZL(I,K-1) * RDZ(I,K)         ! P243.C5 (2nd eqn)
          WK = 0.5 * DZL(I,K) * RDZ(I,K)             ! P243.C5 (1st eqn)

          BTM = WK*BT(I,K) + WKM1*BT(I,K-1)          ! P243.C3 (1st eqn)
          BQM = WK*BQ(I,K) + WKM1*BQ(I,K-1)          ! P243.C3 (2nd eqn)
          BFM = WK*BF(I,K) + WKM1*BF(I,K-1)

          IF (L_BL_LSPICE) THEN
            DZB = BTM*DZTL + BQM*DZQW - BFM*DZQI
          ELSE
            DZB = BTM*DZTL + BQM*DZQW
          ENDIF

!  For rationale of next IF block, see the discussion in the last
!  paragraph of Appendix C of the P243 documentation.

          IF (DZB.GT.0.0) THEN
            BTM = 0.5 * ( BT(I,K) + BT(I,K-1) )
            BQM = 0.5 * ( BQ(I,K) + BQ(I,K-1) )
          BFM = 0.5 * ( BF(I,K) + BF(I,K-1) )
          IF (L_BL_LSPICE) THEN
            DZB = BTM*DZTL + BQM*DZQW - BFM*DZQI
          ELSE
            DZB = BTM*DZTL + BQM*DZQW
          ENDIF
          ENDIF

          BT(I,K-1) = BTM
          BQ(I,K-1) = BQM
          BF(I,K-1) = BFM

        ENDDO ! p_points
      ENDDO ! bl_levels

      IF (LTIMER) THEN
        CALL TIMER('BTQ_INT ',4)
      ENDIF

      RETURN
      END
