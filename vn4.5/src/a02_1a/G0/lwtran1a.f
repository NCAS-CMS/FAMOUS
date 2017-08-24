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
CLL Subroutine LWTRAN   ----------------------------------------------
CLL
CLL            Purpose :
CLL  It calculates clear-sky transmissivities in each of the NBANDS
CLL  longwave spectral bands (and, optionally, additional diagnostic
CLL  ones) from the pathlengths for each effective absorbing gas.
CLL  (Where the absorption by a gas includes terms with different
CLL  pathlength scaling, like water vapour line & continuum, they are
CLL  different gases as far as LWTRAN is concerned.)  It uses look-up
CLL  tables derived from line data as described by Slingo and Wilderspin
CLL  (April 1986, Quart.J.R.Met.Soc., 112, 472, 371-386), or UMDP 23,
CLL  which incorporate a full angular integration.  Interpolation is
CLL  logarithmic in the pathlength, with values at half-decade intervals
CLL  from 10**-9 to 10**3 kg/m2.
CLL    The version of routines LWTRAN and LWLKIN
CLL  used in version 1A (gaseous effects treated as in Slingo &
CLL  Wilderspin, 1986) of the UM LW code.
CLL  LWLKIN must be CALLed to initialize TRTAB before LWTRAN is CALLed
CLL  (LWTRAN would normally be CALLed via LWMAST and LWRAD).
CLL
CLL        Author: William Ingram
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  Date
CLL   4.2    Sept.96  T3E migration: *DEF CRAY removed;
CLL                   *DEF T3E used for T3E library functions;
CLL                   dynamic allocation no longer *DEF controlled;
CLL                   cray HF functions replaced by T3E lib functions.
CLL                       S.J.Swarbrick
CLL
CLL Programming standard :
CLL  It conforms with standard A of version 3 (07/9/90) of UMDP 4, and
CLL  contains no 8X-deprecated features.
CLL  If UPDATE *DEF CRAY is off, the code is standard FORTRAN 77 except
CLL  for having ! comments (it then sets the "vector length" to 1) but
CLL  otherwise it includes an automatic array also.
CLL
CLL Logical components covered : P23
CLL  Component 232 (longwave radiation),
CLL  It is also intended to be easily extended to perform
CLL  some of the functions of D23 (radiation diagnostics), by diagnosing
CLL  additional transmissivities.
CLL
CLL Project task : P23 (radiation)
CLL
CLL External documentation:      UMDP 23.
CLL
CLLEND -----------------------------------------------------------------
C*L
      SUBROUTINE LWTRAN (PATH, TRTAB, DTRTAB, FLSTBD, KEXP,
     &     L,                                                           
     &     TRANS)
C*
      INTEGER NBANDS          ! Number of spectral bands in the longwave
      PARAMETER (NBANDS=6)    ! This run uses the standard set of
! longwave bands as described by Slingo and Wilderspin
! (April 1986, Quart.J.R.Met.Soc., 112, 472, 371-386) or UMDP 23.
      INTEGER NGASES
C Effective number of absorbing gases treated in the longwave
      PARAMETER (NGASES=4)    !  Standard set is water vapour line and
C continuum (which have to be treated separately because the pathlength
C scaling is different), ozone and carbon dioxide.
      INTEGER NTRANS, NDIATR
      PARAMETER (NDIATR=0)
      PARAMETER (NTRANS=NBANDS+NDIATR)
C Number of transmissivities to be calculated - one for each band is
C needed to construct the actual fluxes, but we also allow for NDIATR
C for diagnostic uses, such as possible "narrow-band" flux diagnostics.
      INTEGER IT              ! Dimension of look-up tables
!  Dimension of Lookup Tables
      PARAMETER (IT=25)
      REAL ZSTART,               ! Pathlength where tables start
     &     OFFSET                ! Effective offset for use of tables
      INTEGER LOG10S,            ! Log10 (ZSTART)
     &     NDEC                  ! number of entries per decade
      PARAMETER (LOG10S=-9, NDEC=2)
      PARAMETER (OFFSET=1-LOG10S*NDEC, ZSTART=10.**LOG10S)
C*L
      INTEGER!, INTENT(IN) ::
     &     L,                    ! Number of points                     
     &     FLSTBD(NGASES,2)      ! First & last band in which each
C                                ! effective gas is active
      REAL!, INTENT(IN) ::
     &     PATH (L,NGASES),      ! Scaled pathlengths for each gas
     &     TRTAB(IT,NTRANS,NGASES),  ! Transmissivity look-up table
     &     DTRTAB(IT,NTRANS,NGASES), ! and table of its differences
     &     KEXP(NTRANS)          !  k1 in Eq 2.3.6, used instead of
C                                !  (D)TRTAB for continuum absorption.
      REAL!, INTENT(OUT) ::
     &     TRANS(L,NTRANS)       ! Transmissivities
C*
CL    !  No EXTERNAL routines called
CL    !  Two workspace arrays of size L, one real (Y) & one integer (I) 
C
      REAL RLNR10,               ! NDEC/ln(10)
     &     TGAS,                 ! Transmissivity due to a single "gas"
     &     Y(L)                  ! Used in the interpolation
      INTEGER JTRANS, GAS, J,    ! Loop over transmissivity, gas & point
     &     I(L)                  ! INT(Y)
!
! Local workspace
      REAL LOGPATH(L,3)   
      REAL EXPPATH(FLSTBD(4,2)-FLSTBD(4,1)+1,L) 
      INTEGER fldiff
!
! No. of rows for exp_v32 function
      fldiff = FLSTBD(4,2) - FLSTBD(4,1) + 1
!
      RLNR10 = REAL(NDEC) / LOG(10.)     ! Cannot put this in a
C     !  PARAMETER statement in FORTRAN77, but the CRAY compiler's
C     !  optimizer will make it have the same effect as if it were.
C
C     !  First, initialize the transmissivities to 1 - we will assume
C     !  random overlap of lines of different gases, so that the total
C     !  transmissivity in each band is the product of the
C     !  transmissivities for the individual gases.
C
      DO 1 JTRANS=1, NTRANS
Cfpp$  Select(CONCUR)
       DO 1 J=1, L
        TRANS(J,JTRANS) = 1.
    1 CONTINUE
C
C     ! Then loop through those effective gases which use look-up tables
C
      DO J=1,L
        DO GAS=1,3
          logpath(j,gas)=log(path(j,gas))
        END DO
      END DO
!
      DO 2 GAS=1, 3
Cfpp$  Select(CONCUR)
       DO 20 J=1, L
        Y(J) = LOGPATH(J,GAS) * RLNR10 + OFFSET               
C
C OFFSET ALLLOWS FOR THE START-POINT OF THE TABLES (SAME FOR ALL GASES).
        I(J) = INT(Y(J))
        Y(J) = Y(J) - REAL(I(J))
        IF (I(J).GT.IT) I(J) = IT
   20  CONTINUE
       DO 22 JTRANS=FLSTBD(GAS,1), FLSTBD(GAS,2)
Cfpp$   Select(CONCUR)
        DO 22 J=1, L
         IF (I(J).GE.1) THEN
            TGAS = TRTAB(I(J),JTRANS,GAS) + Y(J)*DTRTAB(I(J),JTRANS,GAS)
          ELSE
            TGAS =
     &   1.  -  ( 1. - TRTAB(1,JTRANS,GAS) )  *  PATH(J,GAS) / ZSTART
         ENDIF
C        !   Assume random overlap of lines of different gases, so that
C        !   the total transmissivity is the product of the
C        !   transmissivities for the individual gases:
         TRANS(J,JTRANS) = TRANS(J,JTRANS) * TGAS
   22  CONTINUE
    2 CONTINUE
C
C     ! Currently H2O continuum is just exponential (2.3.4), & CFCs will
C     ! be too.  Again, transmissivities just multiply:
C
! Use exppath to store products....
      do jtrans = FLSTBD(4,1), FLSTBD(4,2)
        do j    = 1,L
          exppath(jtrans-FLSTBD(4,1)+1,j) = KEXP(JTRANS)*PATH(J,4) 
        end do
      end do

! ....then compute exp of products in exppath
      DO JTRANS=1,FLSTBD(4,2)-FLSTBD(4,1)+1
        DO J=1,L
          exppath(jtrans,j)=exp(exppath(jtrans,j))
        end do
      end do

      DO 3 JTRANS=FLSTBD(4,1), FLSTBD(4,2)
Cfpp$  Select(CONCUR)
       DO 30 J=1, L
        TRANS(J,JTRANS) =
     &          TRANS(J,JTRANS) * EXPPATH(JTRANS-FLSTBD(4,1)+1,J)   
   30  CONTINUE
    3 CONTINUE
C
      RETURN
      END
      SUBROUTINE LWLKIN (LWLUT)
      INTEGER NBANDS          ! Number of spectral bands in the longwave
      PARAMETER (NBANDS=6)    ! This run uses the standard set of
! longwave bands as described by Slingo and Wilderspin
! (April 1986, Quart.J.R.Met.Soc., 112, 472, 371-386) or UMDP 23.
      INTEGER NGASES
C Effective number of absorbing gases treated in the longwave
      PARAMETER (NGASES=4)    !  Standard set is water vapour line and
C continuum (which have to be treated separately because the pathlength
C scaling is different), ozone and carbon dioxide.
      INTEGER NTRANS, NDIATR
      PARAMETER (NDIATR=0)
      PARAMETER (NTRANS=NBANDS+NDIATR)
C Number of transmissivities to be calculated - one for each band is
C needed to construct the actual fluxes, but we also allow for NDIATR
C for diagnostic uses, such as possible "narrow-band" flux diagnostics.
      INTEGER IT              ! Dimension of look-up tables
!  Dimension of Lookup Tables
      PARAMETER (IT=25)
      REAL!, INTENT(OUT)
     &     LWLUT(IT,NTRANS,NGASES,2)
      REAL TRTAB(IT,NTRANS,NGASES),
     &     CO2(IT), WL(IT,6), WC(IT,4), O3(IT)
C     ! Equivalence arrays named after the various gases to the relevant
C     ! areas of TRTAB, to make the DATA statements easier to understand
C     ! and to change the order in which the gases are treated.
      EQUIVALENCE
     &     (CO2(1),TRTAB(1,3,1)), (WL(1,1),TRTAB(1,1,2)),
     &     (O3(1),TRTAB(1,4,3)),  (WC(1,1),TRTAB(1,2,4))
      INTEGER JTRANS, GAS, J     ! Loop over transmissivity, gas & point
      DATA CO2 /
     &   0.999999E+00, 0.999999E+00, 0.999998E+00, 0.999996E+00,
     &   0.999989E+00, 0.999967E+00, 0.999897E+00, 0.999677E+00,
     &   0.998989E+00, 0.996873E+00, 0.990651E+00, 0.973934E+00,
     &   0.935907E+00, 0.868782E+00, 0.777512E+00, 0.671052E+00,
     &   0.553106E+00, 0.430641E+00, 0.316395E+00, 0.218707E+00,
     &   0.137333E+00, 0.704820E-01, 0.259600E-01, 0.564402E-02,
     &   0.467002E-03 /
      DATA O3 /
     + 0.100000E+01, 0.100000E+01, 0.100000E+01, 0.999991E+00,
     + 0.999983E+00, 0.999949E+00, 0.999846E+00, 0.999505E+00,
     + 0.998456E+00, 0.995171E+00, 0.985324E+00, 0.957910E+00,
     + 0.892432E+00, 0.778131E+00, 0.657125E+00, 0.569070E+00,
     + 0.496698E+00, 0.433166E+00, 0.380384E+00, 0.346929E+00,
     + 0.325094E+00, 0.281032E+00, 0.198720E+00, 0.989509E-01,
     + 0.291727E-01 /
      DATA ((WL(J,JTRANS),J=1,IT),JTRANS=1,2) /
     &   0.999998E+00, 0.999997E+00, 0.999992E+00, 0.999977E+00,
     &   0.999928E+00, 0.999775E+00, 0.999296E+00, 0.997833E+00,
     &   0.993562E+00, 0.982165E+00, 0.955877E+00, 0.905200E+00,
     &   0.821497E+00, 0.698231E+00, 0.536604E+00, 0.355461E+00,
     &   0.192011E+00, 0.799150E-01, 0.236460E-01, 0.399601E-02,
     &   0.226021E-03, 0.202656E-05, 0.000000E+00, 0.000000E+00,
     &   0.000000E+00,
     &   0.999999E+00, 0.999999E+00, 0.999999E+00, 0.999999E+00,
     &   0.999998E+00, 0.999994E+00, 0.999985E+00, 0.999954E+00,
     &   0.999857E+00, 0.999552E+00, 0.998612E+00, 0.995824E+00,
     &   0.988135E+00, 0.969417E+00, 0.930882E+00, 0.863802E+00,
     &   0.760731E+00, 0.616712E+00, 0.439173E+00, 0.259404E+00,
     &   0.120350E+00, 0.407310E-01, 0.818002E-02, 0.621974E-03,
     &   0.798702E-05 /
      DATA ((WL(J,JTRANS),J=1,IT),JTRANS=3,4) /
     &   0.100000E+01, 0.100000E+01, 0.100000E+01, 0.999999E+00,
     &   0.999999E+00, 0.999999E+00, 0.999998E+00, 0.999996E+00,
     &   0.999990E+00, 0.999972E+00, 0.999912E+00, 0.999726E+00,
     &   0.999153E+00, 0.997459E+00, 0.992810E+00, 0.981556E+00,
     &   0.958210E+00, 0.915852E+00, 0.846374E+00, 0.742139E+00,
     &   0.600660E+00, 0.432756E+00, 0.265738E+00, 0.130967E+00,
     &   0.466160E-01,
     &   0.100000E+01, 0.100000E+01, 0.100000E+01, 0.100000E+01,
     &   0.100000E+01, 0.100000E+01, 0.100000E+01, 0.999999E+00,
     &   0.999999E+00, 0.999999E+00, 0.999999E+00, 0.999997E+00,
     &   0.999994E+00, 0.999983E+00, 0.999948E+00, 0.999838E+00,
     &   0.999495E+00, 0.998438E+00, 0.995319E+00, 0.986802E+00,
     &   0.966416E+00, 0.925201E+00, 0.853906E+00, 0.743495E+00,
     &   0.586898E+00 /
      DATA ((WL(J,JTRANS),J=1,IT),JTRANS=5,6) /
     &   0.100000E+01, 0.100000E+01, 0.100000E+01, 0.100000E+01,
     &   0.100000E+01, 0.100000E+01, 0.999999E+00, 0.999999E+00,
     &   0.999999E+00, 0.999998E+00, 0.999996E+00, 0.999990E+00,
     &   0.999971E+00, 0.999911E+00, 0.999721E+00, 0.999136E+00,
     &   0.997389E+00, 0.992502E+00, 0.980239E+00, 0.953619E+00,
     &   0.904250E+00, 0.823875E+00, 0.703543E+00, 0.538438E+00,
     &   0.342889E+00,
     &   0.999999E+00, 0.999999E+00, 0.999999E+00, 0.999997E+00,
     &   0.999993E+00, 0.999979E+00, 0.999935E+00, 0.999796E+00,
     &   0.999361E+00, 0.998036E+00, 0.994181E+00, 0.983981E+00,
     &   0.960933E+00, 0.918253E+00, 0.851627E+00, 0.759414E+00,
     &   0.646099E+00, 0.525366E+00, 0.413100E+00, 0.316496E+00,
     &   0.233320E+00, 0.160128E+00, 0.979180E-01, 0.511810E-01,
     &   0.232610E-01 /
C     !  Initialize unused parts to zero to prevent INDEF problems
      DATA ((TRTAB(J,JTRANS,1),J=1,IT),JTRANS=1,2),
     &     ((TRTAB(J,JTRANS,1),J=1,IT),JTRANS=4,6),
     &     ((TRTAB(J,JTRANS,3),J=1,IT),JTRANS=1,3),
     &     ((TRTAB(J,JTRANS,3),J=1,IT),JTRANS=5,6),
     &     ((TRTAB(J,JTRANS,4),J=1,IT),JTRANS=1,6)
     &     / IT*0., IT*0., IT*0., IT*0., IT*0., IT*0., IT*0., IT*0.,
     &     IT*0., IT*0., IT*0., IT*0., IT*0., IT*0., IT*0., IT*0. /
C
      DO 1 GAS=1, NGASES
       DO 1 JTRANS=1, NTRANS
        DO 1 J=1, IT
         LWLUT(J,JTRANS,GAS,1) = TRTAB(J,JTRANS,GAS)
    1 CONTINUE
C
      DO 2 GAS=1, NGASES
       DO 2 JTRANS=1, NTRANS
        DO 2 J=1, IT-1
         LWLUT(J,JTRANS,GAS,2) =
     &    LWLUT(J+1,JTRANS,GAS,1) - LWLUT(J,JTRANS,GAS,1)
    2 CONTINUE
C
C     ! Set the last element for each gas and band to zero, so that the
C     ! extrapolation done for any pathlength greater than the maximum
C     ! catered for just gives the greatest value in TRTAB.
C
      DO 3 GAS=1, NGASES
       DO 3 JTRANS=1, NTRANS
        LWLUT(IT,JTRANS,GAS,2) = 0.
    3 CONTINUE
C
      RETURN
      END
