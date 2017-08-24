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

      SUBROUTINE IMPLSCH (FL3, FL, IJS, IJL, IG, IGL, ishallo,
     & idelt,
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

C*    *COMMON*  *SOURCE* - SOURCE FUNCTION AND TAIL FLAG.
C
     & SL, FCONST,

C*    *COMMON*  *WIND* - VARIABLES USED FOR WIND COMPUTATIONS.
C
     & U10NEW, U10OLD, THWNEW, THWOLD, USNEW, USOLD, Z0NEW,
     & Z0OLD, TAUW,

C*    *COMMON* *SHALLOW*   SHALLOW WATER TABLES.
C
     & DEPTH, DEPTHA, DEPTHD,TCGOND, TFAK, TSIHKD, INDEP,

C*    *COMMON* *COUPL* - PARAMETERS FOR COUPLING.
C
     & BETAMAX, ZALP, ALPHA, XKAPPA, XNLEV,
C
C*    *COMMON* *TABLE* - TABLE FOR TOTAL STRESS AND HIGH FREQ STRESS.
C
     & TAUT, DELTAUW, DELU, TAUHFT, DELUST, DELALP,

C*    *COMMON* *INDNL* - INDICES AND WEIGHTS USED IN THE COMPUTATION
C                        OF THE NONLINEAR TRANSFER RATE.
C
     & IKP, IKP1, IKM, IKM1, K1W, K2W, K11W, K21W, AF11, FKLAP,
     & FKLAP1, FKLAM, FKLAM1, ACL1, ACL2,  CL11, CL21, DAL1, DAL2, FRH,

c* argument list for source term diagnostic arrays for implsch
c stl to hold increments from spectral tail calculations
c other arrays source terms as standard notation
c len_s2 (source diagnostics) =1 or nang*nfre*niblo as required
c
c THIS set for one block only - pass from WAMODEL to implsch
c
     & sin2, snl2, sds2, sbf2, stl2, len_s2,
c
     & icode)

C*    *COMMON* *SHALLOW*   SHALLOW WATER TABLES.
C     ndepth = length of shallow water tables
      integer ndepth
      PARAMETER (NDEPTH = 52)
C
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
      PARAMETER (GZPI28 = G/28./ZPI)

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
C*    *COMMON*  *SOURCE* - SOURCE FUNCTION AND TAIL FLAG.
C
      real SL(0:NIBLO,NANG,NFRE) ! total source function array
      real FCONST(NIBLO,NFRE)  ! tail flag=1/0 for prognostic/diagnostic
CCREFRA
c! for propagation with refraction
c! sl     = sigma dot term
c! fconst = sigma/sinh2kd
CREFRA
C
C*    *COMMON*  *WIND* - VARIABLES USED FOR WIND COMPUTATIONS.
C
C   U10NEW etc changed to dimension (niblo,nblo) MH 9/6/95

      real U10NEW(NIBLO,nblo) ! new wind speed m/s
      real U10OLD(NIBLO,NBLO) ! intermediate storage windspeed
      real THWNEW(NIBLO,nblo) ! wind direction rads / oceanographic
      real THWOLD(NIBLO,NBLO) ! intermediate storage wind direction
      real USNEW (NIBLO,nblo) ! new friction velocity ustar
      real USOLD (NIBLO,NBLO) ! intermediate storage ustar
      real Z0NEW (NIBLO,nblo) ! new roughness length (m)
      real Z0OLD (NIBLO,NBLO) ! intermediate storage Z0
      real TAUW(NIBLO,NBLO)   ! wave stress in (m/s)**2
C
C*    *COMMON* *SHALLOW*   SHALLOW WATER TABLES.
C
      real DEPTH(NIBLO, NBLO)  ! water depth (metres)
      real DEPTHA, DEPTHD      ! min depth and increment for tables (m)
      real TCGOND(NDEPTH,NFRE) ! shallow water group velocity table
      real TFAK(NDEPTH,NFRE)   ! wave number table
      real TSIHKD(NDEPTH,NFRE) ! table for omega /sinh(2kd)

      integer INDEP(NIBLO)     ! depth index for gridpoint :one block

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
c       local diagnostic arrays passed through argument list
c THIS set for one block only - pass down to IMPLSCH
c
       real sin2 (len_s2)
       real snl2 (len_s2)
       real sds2 (len_s2)
       real sbf2 (len_s2)
       real stl2 (len_s2)
c
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
C**** *IMPLSCH* - IMPLICIT SCHEME FOR TIME INTEGRATION OF SOURCE
C****             FUNCTIONS.
C
C     S.D.HASSELMANN.  MPI
C     H. GUENTHER AND L. ZAMBRESKY  OPTIMIZATION PERFORMED.
C     H. GUENTHER      GKSS/ECMWF   OCTOBER 1989  NEW WIND FIELD
C                                                 INTERFACE AND
C                                                 TIME COUNTING
C     P.A.E.M. JANSSEN KNMI         AUGUST  1990  COUPLED MODEL
C     H. GUENTHER      GKSS/ECMWF   JUNE    1991  NEW SEPARATION OF
C                                                  DIAG- AND PROGNOSTIC
C                                                  PART OF SPECTRUM.
C
C*    PURPOSE.
C     --------
C
C       THE IMPLICIT SCHEME ENABLES THE USE OF A TIMESTEP WHICH IS
C       LARGE COMPARED WITH THE CHARACTERISTIC DYNAMIC TIME SCALE.
C       THE SCHEME IS REQUIRED FOR THE HIGH FREQUENCIES WHICH
C       RAPIDLY ADJUST TO A QUASI-EQUILIBRIUM.
C
C**   INTERFACE.
C     ----------
C
C       *CALL* *IMPLSCH (FL3, FL, IJS, IJL, IG, IGL)*
C          *FL3*    - FREQUENCY SPECTRUM(INPUT AND OUTPUT).
C          *FL*     - DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE
C          *IJS*    - INDEX OF FIRST GRIDPOINT
C          *IJL*    - INDEX OF LAST GRIDPOINT
C          *IG*     - BLOCK NUMBER
C          *IGL*    - NUMBER OF BLOCKS
C
C     METHOD.
C     -------
C
C       THE SPECTRUM AT TIME (TN+1) IS COMPUTED AS
C       FN+1=FN+DELT*(SN+SN+1)/2., WHERE SN IS THE TOTAL SOURCE
C       FUNCTION AT TIME TN, SN+1=SN+(DS/DF)*DF - ONLY THE DIAGONAL
C       TERMS OF THE FUNCTIONAL MATRIX DS/DF ARE COMPUTED, THE
C       NONDIAGONAL TERMS ARE NEGLIGIBLE.
C       THE ROUTINE IS CALLED AFTER PROPAGATION FOR TIME PERIOD
C       BETWEEN TWO PROPAGATION CALLS - ARRAY FL3 CONTAINS THE
C       SPECTRUM AND FL IS USED AS AN INTERMEDIATE STORAGE FOR THE
C       DIAGONAL TERM OF THE FUNCTIONAL MATRIX.
C
C     EXTERNALS.
C     ---------
C
C       *FEMEAN*    - COMPUTATION OF MEAN FREQUENCY AT EACH GRID POINT.
CSHALLOW
C       *SBOTTOM*   - COMPUTES BOTTOM DISSIPATION SOURCE TERM AND
C                     LINEAR CONTRIBUTION TO FUNCTIONAL MATRIX.
CSHALLOW
C       *SDISSIP*   - COMPUTATION OF DISSIPATION SOURCE FUNCTION
C                     AND LINEAR CONTRIBUTION OF DISSIPATION TO
C                     FUNCTIONAL MATRIX IN IMPLICIT SCHEME.
C       *SEMEAN*    - COMPUTATION OF TOTAL ENERGY AT EACH GRID POINT.
C       *SINPUT*    - COMPUTATION OF INPUT SOURCE FUNCTION, AND
C                     LINEAR CONTRIBUTION OF INPUT SOURCE FUNCTION
C                     TO FUNCTIONAL MATRIX IN IMPLICIT SCHEME.
C       *SNONLIN*   - COMPUTATION OF NONLINEAR TRANSFER RATE AND
C                     DIAGONAL LINEAR CONTRIBUTION OF NONLINEAR SOURCE
C                     FUNCTION TO  FUNCTIONAL MATRIX.
C       *STRESSO*   - COMPUTATION NORMALISED WAVE STRESS.
C           !!!!!!! MAKE SURE THAT SINPUT IS CALLED FIRST, STRESSO
C           !!!!!!! NEXT, AND THEN THE REST OF THE SOURCE FUNCTIONS.
C
C     REFERENCE.
C     ----------
C
C       S. HASSELMANN AND K. HASSELMANN, "A GLOBAL WAVE MODEL",
C       30/6/85 (UNPUBLISHED NOTE)
C
C ----------------------------------------------------------------------
C
      DIMENSION FL(0:NIBLO,NANG,NFRE), FL3(0:NIBLO,NANG,NFRE)

cc    local array used when extracting source term diagnostics
      real temp2(0:niblo,nang,nfre)
C
C ----------------------------------------------------------------------
C
      DIMENSION MIJ(NIBLO), MFMF(NIBLO), GADIAG(NIBLO),
     1          TEMP(NIBLO,NFRE), DELFL(NFRE)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   all these are local arrays - since SL is initialised to zero in
c   this subroutine
cc
cc  It seems that sl is been used simply as a temp
cc  work space in WAM
cc
cc    EQUIVALENCE (SL(1,3,1), MIJ(1))
cc    EQUIVALENCE (SL(1,5,1), MFMF(1))
cc    EQUIVALENCE (SL(1,7,1), GADIAG(1))
cc    EQUIVALENCE (SL(1,9,1), TEMP(1,1))
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc
         DELT = IDELT
         DELT5 = 0.5*DELT
C ----------------------------------------------------------------------
C
C*    1. INITIALISATION.
C        ---------------
C
 1000 CONTINUE
C
C ----------------------------------------------------------------------
c initialisation of local diagnostic arrays

       do i=1,len_s2
        sin2 (i)=0.
        snl2 (i)=0.
        sds2 (i)=0.
        sbf2 (i)=0.
        stl2 (i)=0.
       enddo

C ----------------------------------------------------------------------
C
C*    2. COMPUTATION OF IMPLICIT INTEGRATION.
C        ------------------------------------
C
C         INTEGRATION IS DONE FOR A BLOCK
C         OF LATITUDES BETWEEN PROPAGATION CALLS.
C
 2000 CONTINUE
C ----------------------------------------------------------------------
C
C*    2.2 COMPUTE MEAN PARAMETERS.
C         ------------------------
C
 2200 CONTINUE

         CALL SEMEAN(FL3, IJS, IJL,
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

     & icode)

         CALL FEMEAN(FL3, IJS, IJL, ishallo,
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
C
C ----------------------------------------------------------------------
C
C*    2.3 COMPUTATION OF SOURCE FUNCTIONS.
C         --------------------------------
C
 2300 CONTINUE
C
C*    2.3.1 INITIALISE SOURCE FUNCTION AND DERIVATIVE ARRAY.
C           ------------------------------------------------
C
         DO 2311 M=1,NFRE
         DO 2311 K=1,NANG
         DO 2311 IJ=0,NIBLO
            SL(IJ,K,M) = 0.
            FL(IJ,K,M) = 0.
 2311    CONTINUE
C
C*    2.3.2 ADD SOURCE FUNCTIONS AND WAVE STRESS.
C           -------------------------------------
C
         CALL SINPUT (FL3, FL, IJS, IJL, IG, ishallo,
C*    *PARAMETER*  FOR ARRAY DIMENSIONS.  parall
C
     & NANG, NFRE, NGX, NGY, NBLO, NIBLO, NOVER,
     & NIBLD, NBLD, NIBLC, NBLC,
C*    *COMMON* *COUPL* - PARAMETERS FOR COUPLING.
C
     & BETAMAX, ZALP, ALPHA, XKAPPA, XNLEV,
C
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

C*    *COMMON*  *SOURCE* - SOURCE FUNCTION AND TAIL FLAG.
C
     & SL, FCONST,

C*    *COMMON*  *WIND* - VARIABLES USED FOR WIND COMPUTATIONS.
C
     & U10NEW, U10OLD, THWNEW, THWOLD, USNEW, USOLD, Z0NEW,
     & Z0OLD, TAUW,

     & icode)

c extract diagnostics if required
         if(len_s2.eq.nang*nfre*niblo) then
           WRITE(6,*)'extracting diagnostics Sinput'
           do l=1,nfre
            do m=1,nang
             nstart=((l-1)*nang + m-1)*niblo
             do ip=ijs,ijl
              sin2(nstart+ip)=sl(ip,m,l)*delt
              temp2(ip,m,l)=sl(ip,m,l)
             enddo
            enddo
           enddo
         endif

         CALL STRESSO (FL3, IJS, IJL, IG, igl,
C*    *PARAMETER*  FOR ARRAY DIMENSIONS.  parall
C
     & NANG, NFRE, NGX, NGY, NBLO, NIBLO, NOVER,
     & NIBLD, NBLD, NIBLC, NBLC,
C*    *COMMON* *COUPL* - PARAMETERS FOR COUPLING.
C
     & BETAMAX, ZALP, ALPHA, XKAPPA, XNLEV,
C
C*    *COMMON* *FREDIR* - FREQUENCY AND DIRECTION GRID.
C
     & FR, DFIM, GOM, C, DELTH, DELTR, TH, COSTH, SINTH,
C
C*    *COMMON*  *SOURCE* - SOURCE FUNCTION AND TAIL FLAG.
C
     & SL, FCONST,

C*    *COMMON* *TABLE* - TABLE FOR TOTAL STRESS AND HIGH FREQ STRESS.
C
     & TAUT, DELTAUW, DELU, TAUHFT, DELUST, DELALP,

C*    *COMMON*  *WIND* - VARIABLES USED FOR WIND COMPUTATIONS.
C
     & U10NEW, U10OLD, THWNEW, THWOLD, USNEW, USOLD, Z0NEW,
     & Z0OLD, TAUW,

     & icode)

         CALL SNONLIN (FL3, FL, IJS, IJL, IG, ishallo,
C*    *PARAMETER*  FOR ARRAY DIMENSIONS.  parall
C
     & NANG, NFRE, NGX, NGY, NBLO, NIBLO, NOVER,
     & NIBLD, NBLD, NIBLC, NBLC,
C*    *COMMON* *INDNL* - INDICES AND WEIGHTS USED IN THE COMPUTATION
C                        OF THE NONLINEAR TRANSFER RATE.
C
     & IKP, IKP1, IKM, IKM1, K1W, K2W, K11W, K21W, AF11, FKLAP,
     & FKLAP1, FKLAM, FKLAM1, ACL1, ACL2,  CL11, CL21, DAL1, DAL2, FRH,

C*    *COMMON*  *MEANPA* - INTEGRATED PARAMETERS OF A BLOCK.
C
     & EMEAN, FMEAN, THQ, AKMEAN,

C*    *COMMON* *SHALLOW*   SHALLOW WATER TABLES.
C
     & DEPTH, DEPTHA, DEPTHD,TCGOND, TFAK, TSIHKD, INDEP,

C*    *COMMON*  *SOURCE* - SOURCE FUNCTION AND TAIL FLAG.
C
     & SL, FCONST,

     & icode)

c extract diagnostics if required
         if(len_s2.eq.nang*nfre*niblo) then
           WRITE(6,*)'extracting diagnostics Snl'
           do l=1,nfre
            do m=1,nang
             nstart=((l-1)*nang + m-1)*niblo
             do ip=ijs,ijl
              snl2(nstart+ip)=(sl(ip,m,l) - temp2(ip,m,l))*delt
              temp2(ip,m,l)=sl(ip,m,l)
             enddo
            enddo
           enddo
         endif

         CALL SDISSIP (FL3 ,FL, IJS, IJL, ishallo,
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

C*    *COMMON*  *SOURCE* - SOURCE FUNCTION AND TAIL FLAG.
C
     & SL, FCONST,

     & icode)

c extract diagnostics if required
         if(len_s2.eq.nang*nfre*niblo) then
           WRITE(6,*)'extracting diagnostics Sds'
           do l=1,nfre
            do m=1,nang
             nstart=((l-1)*nang + m-1)*niblo
             do ip=ijs,ijl
              sds2(nstart+ip)=(sl(ip,m,l) - temp2(ip,m,l))*delt
              temp2(ip,m,l)=sl(ip,m,l)
             enddo
            enddo
           enddo
         endif

CSHALLOW
         IF(ISHALLO.NE.1) then
          CALL SBOTTOM (FL3, FL, IJS, IJL, IG,
C*    *PARAMETER*  FOR ARRAY DIMENSIONS.  parall
C
     & NANG, NFRE, NGX, NGY, NBLO, NIBLO, NOVER,
     & NIBLD, NBLD, NIBLC, NBLC,
C*    *COMMON* *SHALLOW*   SHALLOW WATER TABLES.
C
     & DEPTH, DEPTHA, DEPTHD,TCGOND, TFAK, TSIHKD, INDEP,

C*    *COMMON*  *SOURCE* - SOURCE FUNCTION AND TAIL FLAG.
C
     & SL, FCONST,

     & icode)

c extract diagnostics if required
         if(len_s2.eq.nang*nfre*niblo) then
           WRITE(6,*)'extracting diagnostics Sbf'
           do l=1,nfre
            do m=1,nang
             nstart=((l-1)*nang + m-1)*niblo
             do ip=ijs,ijl
              sbf2(nstart+ip)=(sl(ip,m,l) - temp2(ip,m,l))*delt
              temp2(ip,m,l)=sl(ip,m,l)
             enddo
            enddo
           enddo
         endif
        ENDIF
CSHALLOW
C ----------------------------------------------------------------------
C
C*    2.4 COMPUTATION OF NEW SPECTRA.
C         ---------------------------
C
C       INCREASE OF SPECTRUM IN A TIME STEP IS LIMITED TO A FINITE
C       FRACTION OF A TYPICAL F**(-5) EQUILIBRIUM SPECTRUM.
C
 2400 CONTINUE

         DO 2401 M=1,NFRE
cc
CCMH  note this term 1200 limits delt to be a multiple of 20 mins
CCMH  or else some constant in here is hardwired to 20 minutes
cc
            DELFL(M) = 0.62E-04*FR(M)**(-5.)*DELT/1200.
            DO 2402 K=1,NANG
               DO 2403 IJ=IJS,IJL
                  GTEMP1 = MAX((1.-DELT5*FL(IJ,K,M)),1.)
                  GTEMP2 = DELT*SL(IJ,K,M)/GTEMP1
                  FLHAB = ABS(GTEMP2)
                  FLHAB = MIN(FLHAB,DELFL(M))
                  FL3(IJ,K,M) = FL3(IJ,K,M) + SIGN(FLHAB,GTEMP2)
                  FL3(IJ,K,M) = MAX(FL3(IJ,K,M),0.)
 2403          CONTINUE
 2402       CONTINUE
 2401    CONTINUE
C
C ----------------------------------------------------------------------
C
C*    2.5 REPLACE DIAGNOSTIC PART OF SPECTRA BY A F**(-5) TAIL.
C         -----------------------------------------------------
C
 2500 CONTINUE
C
C*    2.5.1 COMPUTE MEAN PARAMETERS.
C           ------------------------
C
         CALL SEMEAN(FL3, IJS, IJL,
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

     & icode)

         CALL FEMEAN(FL3, IJS, IJL, ishallo,
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
C
C*    2.5.2 COMPUTE LAST FREQUENCY INDEX OF PROGNOSTIC PART OF SPECTRUM.
C*          FREQUENCIES LE MAX(4*F(PM) , 2.5*FMEAN).
C           ------------------------------------------------------------
C
         FPMH = 2.5/FR(1)
         FPM = 4.*GZPI28/FR(1)
cc
ccmh note from elsewhere that 24.1598 is 1./(log10(1.1))
cc
         DO 2521 IJ=IJS,IJL
            FPM4 = FPM/(USNEW(IJ,ig)+0.1E-9)
            MIJ(IJ) = ALOG10(FPM4)*24.1589+2.
            FPM4 = FMEAN(IJ)*FPMH
            MFMF(IJ) = ALOG10(FPM4)*24.1589+1.
 2521    CONTINUE

         DO 2522 IJ=IJS,IJL
            MIJ(IJ) = MAX(MFMF(IJ),MIJ(IJ))
            MIJ(IJ) = MIN(MIJ(IJ),NFRE)
 2522    CONTINUE
C
C*    2.5.3 COMPUTE TAIL ENERGY RATIOS.
C           ---------------------------
C
         DO 2531 M=1,NFRE
            DELFL(M) = (1./FR(M))**5.
 2531    CONTINUE
         DO 2532 IJ=IJS,IJL
            GADIAG(IJ) = FR(MIJ(IJ))**5.
 2532    CONTINUE
C
C*    2.5.4 MERGE TAIL INTO SPECTRA.
C           ------------------------
C
         DO 2541 M=1,NFRE
            DO 2542 IJ=IJS,IJL
               FCONST(IJ,M) = 0.
               TEMP(IJ,M) = GADIAG(IJ)*DELFL(M)
 2542       CONTINUE
 2541    CONTINUE
         DO 2543 IJ=IJS,IJL
            J = MIJ(IJ)
            DO 2544 M=1,J
               FCONST(IJ,M) = 1.
               TEMP(IJ,M) = 0.
 2544       CONTINUE
 2543    CONTINUE
C
         DO 2545 K=1,NANG
            DO 2546 IJ=IJS,IJL
               GADIAG(IJ) = FL3(IJ,K,MIJ(IJ))
 2546       CONTINUE
            DO 2547 M=1,NFRE
               DO 2548 IJ=IJS,IJL
                   FL3(IJ,K,M) = GADIAG(IJ)*TEMP(IJ,M)
     1                         + FL3(IJ,K,M)*FCONST(IJ,M)
 2548          CONTINUE
 2547       CONTINUE
 2545    CONTINUE

c extract diagnostics if required
         if(len_s2.eq.nang*nfre*niblo) then
           WRITE(6,*)'extracting diagnostics Stail'
           do l=1,nfre
            do m=1,nang
             nstart=((l-1)*nang + m-1)*niblo
             do ip=ijs,ijl
              stl2(nstart+ip)=(sl(ip,m,l) - temp2(ip,m,l))*delt
             enddo
            enddo
           enddo
         endif
C
      RETURN
      END
