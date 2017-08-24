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
!+ Convert thetal,qt to theta,q,cl,cf,cloud amount.
!
! Subroutine Interface:
      SUBROUTINE THLQT2THQ(P_FIELD,Q_LEVELS,
     &                     PSTAR,P_EXNER,
     &                     AKH,BKH,AK,BK,RHCRIT,
     &                     THETA,Q,QCF,QCL,RHCPT,
     &                     ICODE)

      IMPLICIT NONE
!
! Description:
!   Convert thetal,qt to theta,q,cl,cf,cloud amount.
!
! Method:
!   Called in ATM_STEP to split conserved variables for 2nd and
!   subsequent loops over dynamics if using long physics timestep
!   option. Uses potl.temp to & from temp conversion and routine LS_CLD
!   to recover cloud water, ice & amount from conserved variables.
!
! Current Code Owner: R.T.H.Barnes (FR)
!
! History:
! Version  Date         Comment
! -------  ----         -------
!  4.0  31/08/95  New routine. R.T.H.Barnes.
!  4.3  3/6/97    Interface and variable names changed to allow calling
!                 from atmdyn.                            P.Burton
!  4.5  13/05/98  Change to subroutine statement: new variable passed
!                 in, and altered call to GLUE_CLD.  S. Cusack
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 v6 programming standards.
!
! Subroutine arguments
      INTEGER
     &  P_FIELD   ! IN : size of horizonal field
     &, Q_LEVELS  ! IN : number of moist levels
     &, ICODE     ! OUT: return code

      REAL
     &  PSTAR(P_FIELD)              ! IN
     &, P_EXNER(P_FIELD,Q_LEVELS+1) ! IN
     &, AKH(Q_LEVELS+1)             ! IN
     &, BKH(Q_LEVELS+1)             ! IN
     &, AK(Q_LEVELS)                ! IN
     &, BK(Q_LEVELS)                ! IN
     &, RHCRIT(Q_LEVELS+1)          ! IN
     &, RHCPT(P_FIELD,Q_LEVELS)     ! IN
     &, THETA(P_FIELD,Q_LEVELS)     ! IN/OUT
     &, Q(P_FIELD,Q_LEVELS)         ! IN/OUT
     &, QCF(P_FIELD,Q_LEVELS)       ! IN/OUT
     &, QCL(P_FIELD,Q_LEVELS)       ! IN/OUT



! Constants
C*L -----------------COMDECK PHYSCONS----------------------------------
C
C   Purpose : contains physical constants required by the whole of the
C             model. It is made up of individual COMDECKS for sets of
C             of related constants, each routine can access one or
C             several of these COMDECKS seperately
C   System Component : F07
C   System task : Z
C  END
C*----------------------------------------------------------------------
C*L------------------COMDECK C_A----------------------------------------
C A IS MEAN RADIUS OF EARTH
      REAL A

      PARAMETER(A=6371229.)
C*----------------------------------------------------------------------

C*L------------------COMDECK C_G----------------------------------------
C G IS MEAN ACCEL DUE TO GRAVITY AT EARTH'S SURFACE
      REAL G

      PARAMETER(G=9.80665)
C*----------------------------------------------------------------------
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

C*L-----------------COMDECK C_VKMAN-------------------------------------
C VKMAN IS VON KARMAN'S CONSTANT
      REAL VKMAN

      PARAMETER(VKMAN=0.4)
C*----------------------------------------------------------------------

C*L------------------COMDECK C_EPSLON-----------------------------------
C EPSILON IS RATIO OF MOLECULAR WEIGHTS OF WATER AND DRY AIR
      REAL EPSILON,C_VIRTUAL

      PARAMETER(EPSILON=0.62198,
     &          C_VIRTUAL=1./EPSILON-1.)
C*----------------------------------------------------------------------

C*L------------------COMDECK C_OMEGA------------------------------------
C OMEGA IS MAGNITUDE OF EARTH'S ANGULAR VELOCITY
      REAL OMEGA

      PARAMETER(OMEGA=7.292116E-5)
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

C*L------------------COMDECK C_PI---------------------------------------
CLL
CLL 4.0 19/09/95  New value for PI. Old value incorrect
CLL               from 12th decimal place. D. Robinson
CLL
      REAL PI,PI_OVER_180,RECIP_PI_OVER_180

      PARAMETER(
     & PI=3.14159265358979323846, ! Pi
     & PI_OVER_180 =PI/180.0,   ! Conversion factor degrees to radians
     & RECIP_PI_OVER_180 = 180.0/PI ! Conversion factor radians to
     &                              ! degrees
     & )
C*----------------------------------------------------------------------

C*L------------------COMDECK C_KT_FT-----------------------------------
      REAL KT2MS,    ! Knots to m/s conversion
     &     FT2M      ! Feet to meters conversion

      PARAMETER(
     & KT2MS=1852.0/3600.0,
     & FT2M =0.3048)
C*----------------------------------------------------------------------
! History:
! Version  Date  Comment
!  3.4   18/5/94 Add PP missing data indicator. J F Thomson
C*L------------------COMDECK C_MDI-------------------------------------
      REAL    RMDI_PP   ! PP missing data indicator (-1.0E+30)
      REAL    RMDI_OLD  ! Old real missing data indicator (-32768.0)
      REAL    RMDI      ! New real missing data indicator (-2**30)
      INTEGER IMDI      ! Integer missing data indicator

      PARAMETER(
     & RMDI_PP  =-1.0E+30,
     & RMDI_OLD =-32768.0,
     & RMDI =-32768.0*32768.0,
     & IMDI =-32768)
C*----------------------------------------------------------------------
C*L------------------COMDECK C_LAPSE ----------------------------------
      REAL LAPSE,LAPSE_TROP
      PARAMETER(LAPSE=0.0065)     !  NEAR SURFACE LAPSE RATE
      PARAMETER(LAPSE_TROP=0.002) !  TROPOPAUSE LAPSE RATE
C*----------------------------------------------------------------------



! Local parameters:
! Local scalars:
      INTEGER   I,K   ! Loop counters over P_FIELD,Q_LEVELS
      REAL      PL,PU ! Lower and upper pressure values

! Local dynamic arrays:
      REAL   CF(P_FIELD,Q_LEVELS) ! cloud fraction from LS_CLD
      REAL   P_X_C(P_FIELD,Q_LEVELS) ! save P_EXNER_C for speed
      REAL   LS_GRID_QC(P_FIELD,Q_LEVELS) ! Qc from LS_CLD
      REAL   LS_BS(P_FIELD,Q_LEVELS) ! bs from LS_CLD

! Function & Subroutine calls:
      External   GLUE_CLD
C*L------------------ COMDECK P_EXNERC ---------------------------------
C statement function to define exner pressure at model levels
C from exner pressures at model half-levels
      REAL P_EXNER_C
      REAL R_P_EXNER_C
      REAL P_EXU_DUM,P_EXL_DUM,PU_DUM,PL_DUM,KAPPA_DUM
C arithmetic mean
      P_EXNER_C(P_EXU_DUM,P_EXL_DUM,PU_DUM,PL_DUM,KAPPA_DUM) =
     & 0.5*(P_EXU_DUM + P_EXL_DUM)
      R_P_EXNER_C(P_EXU_DUM,P_EXL_DUM,PU_DUM,PL_DUM,KAPPA_DUM) =
     & 2.0/(P_EXU_DUM + P_EXL_DUM)

C*------------------- --------------------------------------------------


!- End of header

! 1. Convert thetal and qt to theta and q

! 1.1 Convert thetal to temperaturel
      DO  K = 1,Q_LEVELS
        DO  I = 1,P_FIELD
          PU=PSTAR(I)*BKH(K+1)+AKH(K+1)
          PL=PSTAR(I)*BKH(K)+AKH(K)
          P_X_C(I,K) =
     &      P_EXNER_C(P_EXNER(I,K+1),P_EXNER(I,K),PU,PL,KAPPA)

          THETA(I,K)=THETA(I,K)*P_X_C(I,K)
        END DO ! I
      END DO ! K

! 1.2 Call LS_CLD to convert to temperature, q and cloud variables
      CALL GLUE_CLD(AK,BK,PSTAR,RHCRIT,Q_LEVELS,RHCPT,P_FIELD,P_FIELD,
     &              THETA,CF,Q,QCF,QCL,LS_GRID_QC,LS_BS,ICODE)

! 1.3 Convert temperature back to theta
      DO  K = 1,Q_LEVELS
        DO  I = 1,P_FIELD
          THETA(I,K)=THETA(I,K)/P_X_C(I,K)
        END DO ! I
      END DO ! K

      RETURN
      END

