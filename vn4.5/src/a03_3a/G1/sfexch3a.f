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
C*LL  SUBROUTINE SF_EXCH------------------------------------------------
CLL
CLL  Purpose: Calculate coefficients of turbulent exchange between
CLL           the surface and the lowest atmospheric layer, and
CLL           "explicit" fluxes between the surface and this layer.
CLL
CLL  Suitable for Single Column use if *IF definition IBM is selected.
CLL
CLL          Canopy evaporation made implicit
CLL     with respect to canopy water content (requiring TIMESTEP to be
CLL     passed in).
CLL
CLL CW, FH      <- programmer of some or all of previous code or changes
CLL
CLL  Model            Modification history:
CLL version  Date
CLL   3.4  18/10/94   *DECK inserted into UM version 3.4. S Jackson
CLL
CLL   4.0  05/01/95   rhostar*cD*vshr before interpolation output
CLL                   via argument list.        R.N.B.Smith
CLL
CLL   4.0  30/12/94   Changes to formulation of cH and 10m wind
CLL                   diagnostic;  z0h set to z0(veg).
CLL                                                 R.N.B.Smith
CLL   4.1   08/05/96  decks A03_2C and A03_3B removed
CLL                                     S D Jackson
CLL
CLL   4.1  08/05/96   Logical switch for rapidly mixing boundary layer
CLL                                               S D Jackson
CLL
!LL   4.1  12/7/96    Added MPP code  P.Burton
CLL   4.2   Oct. 96   T3E migration - *DEF CRAY removed, WHENIMD removed
CLL                                     S J Swarbrick
!LL   4.3  14/01/97   MPP code : Corrected setting of polar rows
!LL                                                     P.Burton
CLL  4.3  09/06/97  Add swapbounds for RHOKM_1 & CDR10M.  
CLL                                   D.Sexton/RTHBarnes
CLL  4.4  08/09/97  L_BL_LSPICE specifies mixed phase precipitation
CLL                 scheme.                     D. Wilson
CLL  4.5    Jul. 98  Kill the IBM specific lines (JCThil)
CLL  Programming standard: Unified Model Documentation Paper No 4,
CLL                        Version 2, dated 18/1/90.
CLL
CLL  System component covered: Part of P243.
CLL
CLL  Project task:
CLL
CLL  Documentation: UM Documentation Paper No 24, section P243.
CLL                 See especially sub-section (ix).
CLL
CLLEND------------------------------------------------------------------
C*
C*L  Arguments ---------------------------------------------------------
      SUBROUTINE SF_EXCH (
     + P_POINTS,LAND_PTS,U_POINTS,ROW_LENGTH,P_ROWS,U_ROWS
     +,LAND_INDEX,P1,GATHER
     +,AK_1,BK_1,TIMESTEP
     +,CANOPY,CATCH,CF_1,ICE_FRACT,LAND_MASK,PSTAR,Q_1
     +,QCF_1,QCL_1,RESIST,ROOTD,SMC,SMVCCL,SMVCWT,LYING_SNOW,T_1,TSTAR
     +,U_1,V_1,U_1_P,V_1_P,U_0,V_0,Z0V,SIL_OROG,Z1,Z0MSEA,HO2R2_OROG
     &,BQ_1,BT_1,BF_1,CD,CH,EA,ES,ESL,FQW_1,FQW_LEAD,FTL_1,FTL_LEAD
     +,TAUX_1,TAUY_1,QW_1,RHOKE,RHOKEA,RHOKES,RHOKESL,RHOKH_1,RHOKM_1
     +,RIB,TL_1,TSTAR_NL,VSHR,Z0H,Z0M,Z0M_EFF,H_BLEND
     &,T1_SD,Q1_SD
     &,RHO_CD_MODV1
     +,CDR10M,CHR1P5M,CER1P5M,FME
     +,SU10,SV10,SQ1P5,ST1P5,SFME
     +,NRML
     &,L_Z0_OROG,L_RMBL,L_BL_LSPICE,ERROR,LTIMER
     +)
      IMPLICIT NONE
C
C  Input variables.  All fields are on P grid except where noted.
C  Fxxx in a comment indicates the file from which the data are taken.
C
C
C       GENERAL NOTES ABOUT GRID-DEFINITION INPUT VARIABLES.
C       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C  For global data :-
C
C  An Arakawa B-grid is assumed in which each pole is represented by a
C  row of P-grid points.  These polar rows are omitted in the input and
C  output of the present subroutine, so that the argument P_ROWS is two
C  less than the total number of P-rows in the grid. Land specific
C  variables that are required as INput by the higher level routine,
C  BDY_LAYR, are stored on P_grid land points only and land pts on polar
C  rows are not input or output by this routine ; diagnostic variables
C  must be defined on land and sea points for post processing.
C  If defined variable IBM is selected then land point calculations are
C  performed using the array LAND_INDEX to select land points. But note
C  that elements of LAND_INDEX define land points on the full field
C  (ie including polar rows).
C
C  Entire fields of UV-grid values are taken as input, but the two
C  polemost rows are (a) not updated, in the case of INOUT fields, or
C  (b) set to zero, in the case of OUT fields.
C
C  For limited-area data :-
C
C  The above applies, but for "polar rows", etc., read "rows at the
C  north and south boundaries of the area", etc.  E.g. if you want to
C  do calculations in UV-rows n to m inclusive, the input data will be
C  on P-rows n to m+1, and UV-rows n-1 to m+1.  P-rows n to m will
C  then be updated.  Land specific variables are processed as for global
C  data.
C
C  For both cases, the following equalities apply amongst the input
C  grid-definition variables :-
C
C            P_POINTS = P_ROWS * ROW_LENGTH
C            U_POINTS = U_ROWS * ROW_LENGTH
C              U_ROWS = P_ROWS + 1
C            LAND_PTS <= P_POINTS
C
C  An error condition is returned if the input variables don't satisfy
C  these equalities.  (There is of course redundancy here; a compromise
C  between economy, clarity and easy dimensioning is intended.)
C
C  NB: All this has severe implications for batching/macrotasking;
C      effectively it can't be done on a shared-memory machine without
C      either rewriting this routine or using expensive synchronizations
C      (or other messy and/or undesirable subterfuges).
C
C
      LOGICAL LTIMER
C
      INTEGER              !    Variables defining grid.
     + P_POINTS            ! IN Number of P-grid points to be processed.
     +,LAND_PTS            ! IN Number of land points to be processed.
     +,U_POINTS            ! IN Number of UV-grid points.
     +,ROW_LENGTH          ! IN No. of points in latitude row (inclusive
C                          !    of endpoints for ltd. area model).
     +,P_ROWS              ! IN Number of rows of data on P-grid.
     +,U_ROWS              ! IN Number of rows of data on UV-grid.
     +,LAND_INDEX(LAND_PTS)! IN Index for compressed land point array;
C                          !    ith element holds position in the FULL
C                          !    field of the ith land pt to be processed
     +,P1                  ! IN First P-point to be processed.
      LOGICAL
     + GATHER              ! IN If true then leads variables are comp-
C                          !    ressed for sea-ice calculations. This
C                          !    saves duplicating calculations if there
C                          !    are a relatively few of sea-ice points.
C                          !    Set to false for a limited area run
C                          !    with a high proportion of sea-ice.
      REAL
     + AK_1                ! IN Hybrid "A" for lowest model layer.
     +,BK_1                ! IN Hybrid "B" for lowest model layer.
     +,TIMESTEP            ! IN Timestep in seconds for EPDT calc.
     +,CANOPY(LAND_PTS)    ! IN Surface water (kg per sq metre).  F642.
     +,CATCH(LAND_PTS)     ! IN Surface capacity (max. surface water)
C                          !    (kg per sq metre).  F6416.
     +,CF_1(P_POINTS)      ! IN Cloud fraction for lowest atmospheric
C                          !    layer (decimal fraction).
     +,ICE_FRACT(P_POINTS) ! IN Fraction of gridbox which is sea-ice.
     +,LYING_SNOW(P_POINTS)! IN Lying snow amount (kg per sq metre).
     +,PSTAR(P_POINTS)     ! IN Surface pressure (Pascals).
     +,Q_1(P_POINTS)       ! IN Specific humidity for lowest atmospheric
C                          !    layer (kg water per kg air).
     +,QCF_1(P_POINTS)     ! IN Cloud ice for lowest atmospheric layer
C                          !    (kg water per kg air).
     +,QCL_1(P_POINTS)     ! IN Cloud liquid water for lowest atm layer
C                          !    (kg water per kg air).
     +,RESIST(LAND_PTS)    ! IN "Stomatal" resistance to evaporation
C                          !    (seconds per metre).  F6415.
     +,ROOTD(LAND_PTS)     ! IN "Root depth" (metres).  F6412.
     +,SMC(LAND_PTS)       ! IN Soil moisture content (kg per sq m).
C                          !    F621.
     +,SMVCCL(LAND_PTS)    ! IN Critical volumetric SMC (cubic metres
C                          !    per cubic metre of soil).  F6232.
     +,SMVCWT(LAND_PTS)    ! IN Volumetric wilting point (cubic m of
C                          !    water per cubic m of soil).  F6231.
C
C    Note: (SMVCCL - SMVCWT) is the critical volumetric available soil
C          moisture content.                            ~~~~~~~~~
C
      REAL                !    (Split to avoid > 19 continuations.)
     + T_1(P_POINTS)      ! IN Temperature for lowest atmospheric layer
C                         !    (Kelvin).
     +,TSTAR(P_POINTS)    ! IN Mean gridsquare surface temperature (K).
     +,U_1(U_POINTS)      ! IN West-to-east wind component for lowest
C                         !    atmospheric layer (m/s).  On UV grid.
     +,V_1(U_POINTS)      ! IN South-to-north wind component for lowest
C                         !    atmospheric layer (m/s).  On UV grid.
     +,U_1_P(P_POINTS)    ! IN West-to-east wind component for lowest
C                         !    atmospheric layer (m/s).  On P grid.
     +,V_1_P(P_POINTS)    ! IN South-to-north wind component for lowest
C                         !    atmospheric layer (m/s).  On P grid.
     +,U_0(U_POINTS)      ! IN West-to-east component of ocean surface
C                         !    current (m/s; ASSUMED zero over land).
C                         !    UV grid.  F615.
     +,V_0(U_POINTS)      ! IN South-to-north component of ocean surface
C                         !    current (m/s; ASSUMED zero over land).
C                         !    UV grid.  F616.
     +,Z0V(P_POINTS)      ! IN Vegetative roughness length (m).  F6418.
     +,SIL_OROG(LAND_PTS) ! IN Silhouette area of unresolved orography
C                         !    per unit horizontal area
     +,Z1(P_POINTS)       ! IN Height of lowest atmospheric level (m).
     &,HO2R2_OROG(LAND_PTS) ! IN Peak to trough height of unresolved
C                         !    orography devided by 2SQRT(2) (m).
      LOGICAL
     + LAND_MASK(P_POINTS) ! IN .TRUE. for land; .FALSE. elsewhere. F60.
     +,SU10                ! IN STASH flag for 10-metre W wind.
     +,SV10                ! IN STASH flag for 10-metre S wind.
     +,SQ1P5               ! IN STASH flag for 1.5-metre sp humidity.
     +,ST1P5               ! IN STASH flag for 1.5-metre temperature.
     +,SFME                ! IN STASH flag for wind mixing energy flux.
     +,L_Z0_OROG           ! IN .TRUE. to use orographic roughness.
     +,L_RMBL              ! IN .TRUE. to use rapidly mixing boundary
C                          !    scheme
     &,L_BL_LSPICE               ! IN
!                              TRUE  Use scientific treatment of mixed
!                                    phase precip scheme.
!                              FALSE Do not use mixed phase precip
!                                    considerations
C
C  Modified (INOUT) variables.
C
      REAL
     + Z0MSEA(P_POINTS)   ! INOUT Sea-surface roughness length for
C                         !       momentum (m).  F617.
C
C  Output variables.
C
      REAL
     + BQ_1(P_POINTS)   ! OUT A buoyancy parameter for lowest atm level
C                       !     ("beta-q twiddle").
     +,BT_1(P_POINTS)   ! OUT A buoyancy parameter for lowest atm level.
C                       !     ("beta-T twiddle").
     &,BF_1(P_POINTS)
!        OUT A buoyancy parameter for lowest atm level.
!            ("beta-F twiddle").
     +,CD(P_POINTS)     ! OUT Bulk transfer coefficient for momentum.
     +,CH(P_POINTS)     ! OUT Bulk transfer coefficient for heat and/or
C                       !     moisture.
     +,CDR10M(U_POINTS)  ! OUT Reqd for calculation of 10m wind (u & v).
C                        !     NBB: This is output on the UV-grid, but
C                        !     with the first and last rows set to a
C                        !     "missing data indicator".
C                        !     Sea-ice leads ignored. See 3.D.7 below.
     +,CHR1P5M(P_POINTS) ! OUT Reqd for calculation of 1.5m temperature.
C                        !     Sea-ice leads ignored. See 3.D.7 below.
     +,CER1P5M(P_POINTS) ! OUT Reqd for calculation of 1.5m sp humidity.
C                        !     Sea-ice leads ignored. See 3.D.7 below.
     &,RHO_CD_MODV1(P_POINTS)
C                       ! OUT rhostar*cD*vshr before horizontal
C                       !     interpolation output as a diagnostic.
     +,EA(P_POINTS)     ! OUT "Explicit" evaporation with only aero-
C                       !     dynamic resistance (+ve), or condensation
C                       !     (-ve), averaged over grid box (kg/m2/s).
     +,ES(P_POINTS)     ! OUT "Explicit" evapotranspiration through a
C                       !     non-aerodynamic resistance.  Always non-
C                       !     negative (kg/m2/s).
     +,ESL(P_POINTS)    ! OUT ES without fractional weighting factor
C                       !     "L" is for "local".
      REAL                !     (Split to avoid > 19 continuations.)
     + FQW_1(P_POINTS)    ! OUT "Explicit" surface flux of QW (i.e.
C                         !     evaporation), on P-grid (kg/m2/s).
     +,FQW_LEAD(P_POINTS) ! OUT FQW for leads fraction of gridsquare.
C                         !     Missing data indicator at non sea-ice.
     +,FTL_1(P_POINTS)    ! OUT "Explicit" surface flux of TL = H/CP.
C                         !     (sensible heat / CP).
     +,FTL_LEAD(P_POINTS) ! OUT FTL for leads fraction of gridsquare,
C                         !     Missing data indicator at non sea-ice.
     +,TAUX_1(U_POINTS)   ! OUT "Explicit" x-component of surface
C                         !     turbulent stress; on UV-grid; first and
C                         !     last rows set to a "missing data
C                         !     indicator". (Newtons per square metre)
     +,TAUY_1(U_POINTS)   ! OUT "Explicit" y-component of surface
C                         !     turbulent stress; on UV-grid; first and
C                         !     last rows set to a "missing data
C                         !     indicator". (Newtons per square metre)
     +,QW_1(P_POINTS)     ! OUT Total water content of lowest
C                         !     atmospheric layer (kg per kg air).
C
      REAL ! Surface exchange coefficients;passed to subroutine IMPL_CAL
     + RHOKE(P_POINTS)  ! OUT For FQW, then *GAMMA(1) for implicit calcs
     +,RHOKEA(P_POINTS) ! OUT For EA, then *GAMMA(1) for implicit calcs.
     +,RHOKES(P_POINTS)  ! OUT For ES, then *GAMMA(1) for implicit calcs
     +,RHOKESL(P_POINTS) ! OUT For ESL,then *GAMMA(1) for implicit calcs
     +,RHOKH_1(P_POINTS) ! OUT For FTL,then *GAMMA(1) for implicit calcs
     +,RHOKM_1(U_POINTS) ! OUT For momentum, then *GAMMA(1) for implicit
C                        !     calculations. NBB: This is output on the
C                        !     UV-grid, but with the first and last
C                        !     rows set to a "missing data indicator".
      REAL ! End of surface exchange coefficients.
     + RIB(P_POINTS)      ! OUT Bulk Richardson number for lowest layer.
     +,TL_1(P_POINTS)     ! OUT Liquid/frozen water temperature for
C                         !     lowest atmospheric layer (K).
     +,TSTAR_NL(P_POINTS) ! OUT TSTAR No Leads: surface temperature
C                         !     over sea-ice fraction of gridsquare.
C                         !     =TSTAR over non sea-ice points.
     +,VSHR(P_POINTS)     ! OUT Magnitude of surface-to-lowest-lev. wind
     +,Z0H(P_POINTS)      ! OUT Roughness length for heat and moisture m
     +,Z0M(P_POINTS)      ! OUT Roughness length for momentum (m).
     +,Z0M_EFF(P_POINTS)  ! OUT Effective roughness length for momentum
     +,H_BLEND(P_POINTS)  ! OUT Blending height
     &,T1_SD(P_POINTS)    ! OUT Standard deviation of turbulent
C                         !     fluctuations of surface layer
C                         !     temperature (K).
     &,Q1_SD(P_POINTS)    ! OUT Standard deviation of turbulent
C                         !     fluctuations of surface layer
C                         !     specific humidity (kg/kg).
     +,FME(P_POINTS)      ! OUT Wind mixing energy flux (Watts/sq m).
      INTEGER
     + NRML(P_POINTS)     ! OUT 1 if surface layer unstable, else 0.
     +,ERROR              ! OUT 1 if grid definition faulty; else 0.
C*
C*L  Symbolic constants ------------------------------------------------
C
C   (1) UM-wide common parameters.
C
C*L------------------COMDECK C_O_DG_C-----------------------------------
C ZERODEGC IS CONVERSION BETWEEN DEGREES CELSIUS AND KELVIN
C TFS IS TEMPERATURE AT WHICH SEA WATER FREEZES
C TM IS TEMPERATURE AT WHICH FRESH WATER FREEZES AND ICE MELTS
      REAL ZERODEGC,TFS,TM

      PARAMETER(ZERODEGC=273.15,
     &          TFS=271.35,
     &          TM=273.15)
C*----------------------------------------------------------------------

C*L------------------COMDECK C_LHEAT------------------------------------
C LC IS LATENT HEAT OF CONDENSATION OF WATER AT 0DEGC
C LF IS LATENT HEAT OF FUSION AT 0DEGC
      REAL LC,LF

      PARAMETER(LC=2.501E6,
     &          LF=0.334E6)
C*----------------------------------------------------------------------
C*L------------------COMDECK C_G----------------------------------------
C G IS MEAN ACCEL DUE TO GRAVITY AT EARTH'S SURFACE
      REAL G

      PARAMETER(G=9.80665)
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

C
C   (2) Boundary Layer local parameters.
C
C*L-----------COMDECK C_CHARNK FOR SUBROUTINE SF_EXCH----------
C CHARNOCK is a constant in the Charnock formula for sea-surface
C          roughness length for momentum (Z0MSEA).
      REAL CHARNOCK

      PARAMETER(CHARNOCK = 0.012)
C*----------------------------------------------------------------------
C*L-----------COMDECK C_DENSTY FOR SUBROUTINE SF_EXCH----------
C RHOSEA    = density of sea water (kg/m3)
C RHO_WATER = density of pure water (kg/m3)
      REAL RHOSEA,RHO_WATER

      PARAMETER(RHOSEA = 1000.0,
     &          RHO_WATER = 1000.0)
C*----------------------------------------------------------------------
C*L------------------COMDECK C_GAMMA------------------------------------
C GAMMA holds the implicit weighting coeff. for up to 30 B.L. levels.
C It is only required for the the number of B.L. levels actually used,
C so it does not need to be set up to 30 when less BL levels are used.
      REAL GAMMA(30)       ! Max of 30 Boundary Layer levels assumed.
C
      DATA GAMMA / 2 * 2.0 , 1.5 , 27 * 1.0 /
C*----------------------------------------------------------------------
C*L-----------COMDECK C_HT_M FOR SUBROUTINE SF_EXCH----------
C Z10M  = height of 10m level for diagnostic calculations (m).
C Z1P5M = height of 1.5m level for diagnostic calculations (m).
      REAL Z10M,Z1P5M

      PARAMETER(Z10M  = 10.0,
     &          Z1P5M = 1.5)
C*----------------------------------------------------------------------
C*L-----------COMDECK C_ROUGH FOR SUBROUTINE SF_EXCH----------
C Z0FSEA = roughness length for free convective heat and moisture
C          transport over the sea (m).
C Z0HSEA = roughness length for free heat and moisture transport
C          over the sea (m).
C Z0MIZ  = roughness length for heat, moisture and momentum over
C          the Marginal Ice Zone (m).
C Z0SICE = roughness length for heat, moisture and momentum over
C          sea-ice (m).
      REAL Z0FSEA,Z0HSEA,Z0MIZ,Z0SICE

      PARAMETER(Z0FSEA = 1.3E-3,
     &          Z0HSEA = 1.0E-4,
     &          Z0MIZ  = 1.0E-1,
     &          Z0SICE = 3.0E-3)
C*----------------------------------------------------------------------
CLL  Model           Modification history :
CLL version  Date
CLL   3.4  18/10/94   *COMDECK inserted into UM version 3.4. S Jackson
CLL
C*L------------------COMDECK C_SURF------------------------------------
      REAL    RI_CRIT   ! Critical Richardson number, where Z0M_EFF=Z0M.
C                       ! Linear interpolation between RIB=0 and RI_CRIT

      REAL    OROG_DRAG_PARAM    ! Tunable parameter in calculation of
C                                ! Effective roughness length for
C                                ! momentum
      PARAMETER(
     & RI_CRIT=1.0,
     & OROG_DRAG_PARAM=0.3)
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
C
C   (3) Derived local parameters.
C
      REAL ETAR,GRCP,LCRCP,LFRCP,LS,LSRCP,H_BLEND_MIN,H_BLEND_MAX
      PARAMETER (
     + ETAR=1./(1.-EPSILON)  ! Used in calc of buoyancy parameter BETAC.
     +,GRCP=G/CP             ! Used in calc of dT across surface layer.
     +,LCRCP=LC/CP           ! Evaporation-to-dT conversion factor.
     +,LFRCP=LF/CP           ! Freezing-to-dT conversion factor.
     +,LS=LF+LC              ! Latent heat of sublimation.
     +,LSRCP=LS/CP           ! Sublimation-to-dT conversion factor.
     +,H_BLEND_MIN=0.0       ! Minimum blending height.
     &,H_BLEND_MAX=1000.0    ! Maximum blending height (m).
     +)
C*
! MPP Common block
! ------------------------ Comdeck PARVARS -------------------------
! Parameters and common blocks required by the MPP-UM
!========================== COMDECK PARPARM ====================
!   Description:
!
!   This COMDECK contains PARAMETERs for the MPP-UM
!
!   Two sets of parameters are set up -
!     i)  for the MPP-UM itself.
!     ii) for the interface to the Message Passing Software.
!
!   History:
!
!   Model    Date     Modification history
!  version
!   4.1      27/1/96  New comdeck based on first section of
!                     old PARVARS.   P.Burton
!   4.2      21/11/96 Add new field type parameter and
!                     magic number used in addressing to indicate
!                     if a calculation is for local data, or data
!                     on the dump on disk (ie. global data)  P.Burton
!   4.2      18/11/96 Moved MaxFieldSize to comdeck AMAXSIZE and
!                     removed Maxbuf.  P.Burton
!   4.2      18/7/96  Removed some unused variables      P.Burton
!   4.4      11/07/97 Reduced MAXPROC to 256 to save memory  P.Burton
!
! ---------------------- PARAMETERS ---------------------
!
! =======================================================
! Parameters needed for the MPP-UM
! =======================================================

      INTEGER   Ndim_max        ! maximum number of spatial dimensions
      PARAMETER (Ndim_max = 3 ) ! 3d data


      INTEGER
     &   fld_type_p           ! indicates a grid on P points
     &,  fld_type_u           ! indicates a grid on U points
     &,  fld_type_unknown     ! indicates a non-standard grid.
      PARAMETER (
     &   fld_type_p=1
     &,  fld_type_u=2
     &,  fld_type_unknown=-1)

      INTEGER
     &   local_data
     &,  global_dump_data
      PARAMETER (
     &   local_data=1        ! Used in addressing to indicate if
     &,  global_dump_data=2) ! calculation is for a local or
!                            ! global (ie. disk dump) size

! =======================================================
! Parameters needed for the Message Passing Software
! =======================================================


      INTEGER
     &   Maxproc              ! Max number of processors
      PARAMETER (
     &   MAXPROC = 256)

      INTEGER
     &   PNorth       ! North processor address in the neighbour array
     &,  PEast        ! East  processor address in the neighbour array
     &,  PSouth       ! South processor address in the neighbour array
     &,  PWest        ! West  processor address in the neighbour array
     &,  NoDomain     ! Value in neighbour array if the domain has
     &                !  no neighbor in this direction. Otherwise
     &                !  the value will be the tid of the neighbor
      PARAMETER (
     &   PNorth   = 1
     &,  PEast    = 2
     &,  PSouth   = 3
     &,  PWest    = 4
     &,  NoDomain = -1)

      INTEGER
     &   BC_STATIC            ! Static boundary conditions
     &,  BC_CYCLIC            ! Cyclic boundary conditions
      PARAMETER (
     &   BC_STATIC = 1
     &,  BC_CYCLIC = 2)

! ---------------------- End of comdeck PARPARM ---------------------
!========================== COMDECK PARCOMM ====================
!
! *** NOTE : This comdeck requires comdeck PARPARM to be *CALLed
!            first.
!
!   Description:
!
!   This COMDECK contains COMMON blocks for the MPP-UM
!
!
!   Two COMMON blocks are defined:
!     i)  UM_PARVAR holds information required by the
!         Parallel Unified Model itself
!     ii) MP_PARVAR holds information required by the interface to
!         the Message Passing Software used by the PUM
!
!   Key concepts used in the inline documentation are:
!     o GLOBAL data - the entire data domain processed by the UM
!     o LOCAL data - the fragment of the GLOBAL data which is
!       stored by this particular process
!     o PERSONAL data - the fragment of the LOCAL data which is
!       updated by this particular process
!     o HALO data - a halo around the PERSONAL data which forms
!       the LOCAL data
!
!     Acronyms used:
!     LPG - Logical Process Grid, this is the grid of logical
!           processors; each logical processor handles one of the
!           decomposed parts of the global data. It does not
!           necessarily represent a physical grid of processors.
!
!   History:
!
!   4.1      27/1/96  New comdeck based on second section of
!                     old PARVARS.   P.Burton
!   4.2     19/08/96  Removed some unused variables, and added
!                     current_decomp_type variable to allow use
!                     of flexible decompositions.
!                     Added nproc_max to indicate the max. number
!                     of processors used for MPP-UM
!                                                      P.Burton
!
! -------------------- COMMON BLOCKS --------------------
!
! =======================================================
! Common block for the Parallel Unified Model
! =======================================================

      INTEGER
     &   first_comp_pe       ! top left pe in LPG
     &,  last_comp_pe        ! bottom right pe in LPG
     &,  current_decomp_type ! current decomposition type
     &,  Offx                ! halo size in EW direction
     &,  Offy                ! halo size in NS direction
     &,  glsize(Ndim_max)    ! global data size
     &,  lasize(Ndim_max)    ! local data size
     &,  blsizep(Ndim_max)   ! personal p data area
     &,  blsizeu(Ndim_max)   ! personal u data area
     &,  datastart(Ndim_max) ! position of personal data in global data
     &                       !   (in terms of standard Fortran array
     &                       !    notation)
     &,  gridsize(Ndim_max)  ! size of the LPG in each dimension
     &,  gridpos(Ndim_max)   ! position of this process in the LPG
!                            ! 0,1,2,...,nproc_x-1 etc.

      LOGICAL
     &    atbase             ! process at the bottom of the LPG
     &,   attop              ! process at the top of the LPG
     &,   atleft             ! process at the left of the LPG
     &,   atright            ! process at the right of the LPG
! NB: None of the above logicals are mutually exclusive

      COMMON /UM_PARVAR/
     &                  first_comp_pe,last_comp_pe
     &,                 current_decomp_type,Offx, Offy
     &,                 glsize,lasize,blsizep,blsizeu
     &,                 datastart,gridsize,gridpos
     &,                 atbase,attop,atleft,atright

! =======================================================
! Common block for the Message Passing Software
! =======================================================

      INTEGER
     &  bound(Ndim_max)           ! type of boundary (cyclic or static)
     &                            !  in each direction
     &, g_lasize(Ndim_max,0:maxproc)
!                                 ! global copy of local data size
     &, g_blsizep(Ndim_max,0:maxproc)
!                                 ! global copy of personal p data area
     &, g_blsizeu(Ndim_max,0:maxproc)
!                                 ! global copy of personal u data area
     &, g_datastart(Ndim_max,0:maxproc)
!                                 ! global copy of datastart
     &, g_gridpos(Ndim_max,0:maxproc)
!                                 ! global copy of gridpos
     &, nproc                     ! number of processors in current
!                                 ! decomposition
     &, nproc_max                 ! maximum number of processors
     &, nproc_x                   ! number of processors in x-direction
     &, nproc_y                   ! number of processors in y-direction
     &, mype                      ! number of this processor
     &                            !  (starting from 0)
     &, neighbour(4)              ! array with the tids of the four
     &                            ! neighbours in the horizontal plane
     &, gc_proc_row_group         ! GID for procs along a proc row
     &, gc_proc_col_group         ! GID for procs along a proc col
     &, gc_all_proc_group         ! GID for all procs

      COMMON /MP_PARVAR/
     &                  bound
     &,                 g_lasize,g_blsizep,g_blsizeu
     &,                 g_datastart,g_gridpos
     &,                 nproc,nproc_max,nproc_x,nproc_y,mype
     &,                 neighbour,gc_proc_row_group
     &,                 gc_proc_col_group, gc_all_proc_group



! ---------------------- End of comdeck PARCOMM -----------------------
! --------------------- End of comdeck PARVARS ---------------------
C*L
C   External subprograms called.
C
      EXTERNAL FCDCH,QSAT,QSAT_WAT,SFL_INT
      EXTERNAL P_TO_UV,UV_TO_P
      EXTERNAL TIMER
C*
C
C   Define local storage.
C
C   (a) Workspace.
C
C*L  Workspace ---------------------------------------------------------
C  25 blocks of real workspace are required, as follows.
      REAL
     + CD_LEAD(P_POINTS)  ! Bulk transfer coefficient for momentum
C                         !  over sea-ice leads.Missing data over non
C                         !  sea-ice points.(Temporary store for Z0MIZ)
     +,CD_MIZ(P_POINTS)   ! Bulk transfer coefficient for momentum
C                         !  over the sea-ice Marginal Ice Zone.
C                         !  Missing data indicator over non sea-ice.
     +,CH_LEAD(P_POINTS)  ! Bulk transfer coefficient for heat and
C                         !  or moisture over sea ice leads.
C                         !  Missing data indicator over non sea-ice.
     +,CH_MIZ(P_POINTS)   ! Bulk transfer coefficient for heat and
C                         !  or moisture over the Marginal Ice Zone.
C                         !  Missing data indicator over non sea-ice.
     +,CD_STD(P_POINTS)   ! Local drag coefficient for
C                         !  calculation of interpolation coefficients
     +,DQ(P_POINTS)       ! Sp humidity difference between surface
C                         !  and lowest atmospheric level (Q1 - Q*).
C                         !  Holds value over sea-ice where ICE_FRACT
C                         !  >0 i.e. Leads contribution not included.
     &,DQI(P_POINTS)
!        Ice water difference between surface
!        and lowest atmospheric level (Q1 - Q*).
!        Holds value over sea-ice where ICE_FRACT
!        >0 i.e. Leads contribution not included.
     +,DQ_LEAD(P_POINTS)  ! DQ for leads fraction of gridsquare.
C                         !  Missing data indicator over non sea-ice.
     &,DQI_LEAD(P_POINTS)
!        DQI for leads fraction of gridsquare.
!        Missing data indicator over non sea-ice.
     +,DTEMP(P_POINTS)    ! Liquid/ice static energy difference
C                         !  between surface and lowest atmospheric
C                         !  level, divided by CP (a modified
C                         !  temperature difference).
C                         !  Holds value over sea-ice where ICE_FRACT
C                         !  >0 i.e. Leads contribution not included.
     +,DTEMP_LEAD(P_POINTS) ! DTEMP for leads fraction of gridsquare.
C                           !  Missing data indicator over non sea-ice.
     +,FRACA(P_POINTS)    ! Fraction of surface moisture flux with
C                         !  only aerodynamic resistance.
     +,PSIS(P_POINTS)     ! Soil moisture availability factor.
     +,QSL(P_POINTS)      ! Saturated sp humidity at liquid/ice
C                         !  temperature and pressure of lowest
C                         !  atmospheric level.
     +,QSTAR(P_POINTS)    ! Surface saturated sp humidity. Holds
C                         !  value over sea-ice where ICE_FRACT > 0.
C                         !  i.e. Leads contribution not included.
     +,QSTAR_LEAD(P_POINTS) ! QSTAR for sea-ice leads.
C                           !  Missing data indicator over non sea-ice.
     +,RESFS(P_POINTS)    ! Combined soil, stomatal and aerodynamic
C                         !  resistance factor = PSIS/(1+RS/RA) for
C                         !  fraction 1 - FRACA of the gridbox
     +,RESFT(P_POINTS)    ! Total resistance factor for moisture
C                         !  transfer to/from the surface
C                         !   = FRACA + (1-FRACA)*RESFS so that
C                         !  1/RT = RESFT/RA
     +,RIB_LEAD(P_POINTS) ! Bulk Richardson no. for sea-ice leads at
C                         !  lowest layer. At non sea-ice points holds
C                         !  RIB for FCDCH calculation, then set to
C                         !  to missing data indicator.
     +,RS(P_POINTS)       ! Stomatal resistance to evaporation.
     +,U_0_P(P_POINTS)    ! West-to-east component of ocean surface
C                         !  current (m/s; zero over land if U_0 OK).
C                         !  P grid.  F615.
     +,V_0_P(P_POINTS)    ! South-to-north component of ocean surface
C                         !  current (m/s; zero over land if V_0 OK).
C                         !  P grid.  F616.
     +,Z0F(P_POINTS)      ! Roughness length for free-convective heat
C                         !  and moisture transport.
     +,Z0FS(P_POINTS)     ! Roughness length for free-convective heat
C                         !  and moisture transport over sea.
     +,Z0HS(P_POINTS)     ! Roughness length for heat and moisture
C                         !  transport over sea.
     &,WIND_PROFILE_FACTOR(P_POINTS)
C                         ! For transforming effective surface transfer
C                         ! coefficients to those excluding form drag.
C
C  Workspace (reqd for compression).
      INTEGER
     + SICE_INDEX(P_POINTS) ! Index vector for gather to sea-ice points
      LOGICAL ITEST(P_POINTS)  ! Used as 'logical' for compression.
C*
C
C   (b) Scalars.
C
      INTEGER
     + I           ! Loop counter (horizontal field index).
     +,J           ! Offset counter within I-loop.
     +,L           ! Loop counter (land point field index).
     +,NSICE       ! Number of sea-ice points.
     +,SI          ! Loop counter (sea-ice field index).
      REAL
     + AL          ! Temporary in calculation of buoyancy parameters.
     +,ALPHAL      ! Temporary in calculation of buoyancy parameters.
     +,BETAC       ! Temporary in calculation of buoyancy parameters.
     +,BETAQ       ! Temporary in calculation of buoyancy parameters.
     +,BETAT       ! Temporary in calculation of buoyancy parameters.
     +,CHN         ! Neutral-stability value of CH, used as a first
C                  ! approximation to the "true" CH.
      REAL         ! (Split to avoid having > 19 continuation lines.)
     + RHOSTAR     ! Surface air density in kg per cubic metre.
     +,TAU         ! Magnitude of surface wind stress over sea.
     +,SMCRIT      ! "Critical" available SMC in kg per sq m.
     +,USHEAR      ! U-component of surface-to-lowest-level wind shear.
     +,VSHEAR      ! V-component of surface-to-lowest-level wind shear.
     +,VSHR2       ! Square of magnitude of surface-to-lowest-level
C                  ! wind shear.
     +,ZETAM       ! Temporary in calculation of CHN.
     +,ZETAH       ! Temporary in calculation of CHN.
     +,ZETA1       ! Temporary in calculation of land roughness
     +,ZETA2       ! Temporary in calculation of land roughness
     +,ZETA3       ! Temporary in calculation of land roughness
     +,Z0          ! Roughness length over land.
     +,EPDT        ! "Potential" Evaporation * Timestep
     &,VS          ! Surface layer friction velocity
     &,VSF1_CUBED  ! Cube of surface layer free convective scaling
C                  ! velocity
     &,WS1         ! Turbulent velocity scale for surface layer
C
C-----------------------------------------------------------------------
CL  0.  Check that the scalars input to define the grid are consistent.
C-----------------------------------------------------------------------
C
      IF (LTIMER) THEN
        CALL TIMER('SFEXCH  ',3)
      ENDIF

      ERROR=0
      IF (                                                          
     +     U_POINTS .NE. (U_ROWS*ROW_LENGTH) .OR.
     +     P_POINTS .NE. (P_ROWS*ROW_LENGTH) .OR.
     +     LAND_PTS .GT.  P_POINTS )  THEN
        ERROR=1
        GOTO6
      ENDIF
C
C-----------------------------------------------------------------------
CL  1.  Construct SICE_INDEX for compression onto sea points in
CL      sea-ice leads calculations.
C-----------------------------------------------------------------------
C
        DO I = 1,P_POINTS
          ITEST(I) = .FALSE.
          IF (ICE_FRACT(I).GT.0.0 .AND. .NOT.LAND_MASK(I))
     &      ITEST(I) = .TRUE.
        ENDDO
C
C  Routine whenimd is functionally equivalent to WHENILE, so ITEST is
C  1 for "False", 0 for "True".
C
C
        NSICE = 0
        DO I=1,P_POINTS                             
          IF(ITEST(I))THEN                        
            NSICE = NSICE + 1                      
            SICE_INDEX(NSICE) = I                      
          END IF                                  
        END DO            
C
C-----------------------------------------------------------------------
CL  2.  Calculate QSAT values required later and components of ocean
CL      current.
C        Done here to avoid loop splitting.
C        QSTAR 'borrowed' to store P at level 1 (just this once).
C        PSIS 'borrowed' to store leads and non sea-ice surface temp.
C-----------------------------------------------------------------------
C
C-----------------------------------------------------------------------
CL  2.1 IF (GATHER) THEN
CL       Calculate temperatures and pressures for QSAT calculations.
CL       Calculate QSAT values. For sea-ice points, separate values
CL       are required for the leads (QSTAR_LEAD) and sea-ice (QSTAR)
CL       fractions respectively. QSTAR_LEAD = missing data, elsewhere.
CL       Use RS to store compressed PSTAR for this section only.
CL       NB Unlike QSTAR, TSTAR values at sea-ice points are gridsq.
CL       means and so include the leads contribution.
CL      ELSE
CL       As above with QSTAR_LEAD done on full field.
CL      ENDIF
C-----------------------------------------------------------------------
      IF (GATHER) THEN
        DO I = 1,P_POINTS
          IF (L_BL_LSPICE) THEN
            TL_1(I) = T_1(I) - LCRCP*QCL_1(I)                 ! P243.9
          ELSE
            TL_1(I) = T_1(I) - LCRCP*QCL_1(I) - LSRCP*QCF_1(I) !P243.9
          ENDIF
          TSTAR_NL(I) = TSTAR(I)
          QSTAR_LEAD(I) = 1.0E30               ! Missing data indicator
          QSTAR(I) = AK_1 + BK_1*PSTAR(I)
        ENDDO
        IF (NSICE.GT.0) THEN
CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
          DO SI = 1,NSICE
            I = SICE_INDEX(SI)
            TSTAR_NL(I) = (TSTAR(I)-(1.0-ICE_FRACT(I)) *TFS)
     +                / ICE_FRACT(I)                          ! P2430.1
            PSIS(SI) = TFS
            RS(SI) = PSTAR(I)
          ENDDO
        ENDIF
        IF (L_BL_LSPICE) THEN
          CALL QSAT_WAT(QSL,TL_1,QSTAR,P_POINTS)
        ELSE
          CALL QSAT(QSL,TL_1,QSTAR,P_POINTS)
        ENDIF

        CALL QSAT(QSTAR,TSTAR_NL,PSTAR,P_POINTS)
C            ...values at sea-ice points contain ice contribution only
        IF (NSICE.GT.0) CALL QSAT(QSTAR_LEAD,PSIS,RS,NSICE)
C            ...values at sea-ice points only
      ELSE
C-----------------------------------------------------------------------
CL  2.1  Single Column Model selected.
CL       Calculate temperatures and pressures for QSAT calculations.
CL       If there is sea-ice, separate values of surface saturated
CL       specific humidity are required for the leads (QSTAR_LEAD)
CL       and sea-ice (QSTAR) fractions respectively.
CL       NB Unlike QSTAR, TSTAR values at sea-ice points are gridsq.
CL       means and so include the leads contribution.
C-----------------------------------------------------------------------
        DO I = 1,P_POINTS
          IF (L_BL_LSPICE) THEN
            TL_1(I) = T_1(I) - LCRCP*QCL_1(I)                 ! P243.9
          ELSE
            TL_1(I) = T_1(I) - LCRCP*QCL_1(I) - LSRCP*QCF_1(I) !P243.9
          ENDIF
          TSTAR_NL(I) = TSTAR(I)
C Set to missing data at non sea-ice points after QSAT.
          PSIS(I) = TSTAR(I)
          QSTAR(I) = AK_1 + BK_1*PSTAR(I)
          IF ( ICE_FRACT(I).GT.0.0 .AND. .NOT.LAND_MASK(I) ) THEN
            TSTAR_NL(I) = (TSTAR(I)-(1.0-ICE_FRACT(I)) *TFS)
     +                / ICE_FRACT(I)                          ! P2430.1
            PSIS(I) = TFS
          ENDIF
        ENDDO
        IF (L_BL_LSPICE) THEN
          CALL QSAT_WAT(QSL,TL_1,QSTAR,P_POINTS)
        ELSE
          CALL QSAT(QSL,TL_1,QSTAR,P_POINTS)
        ENDIF

        CALL QSAT(QSTAR,TSTAR_NL,PSTAR,P_POINTS)
C          ...values at sea-ice points contain ice contribution only
        IF (NSICE.GT.0) CALL QSAT(QSTAR_LEAD,PSIS,PSTAR,P_POINTS)
C          ...values at sea-ice points contain leads contribution only
        DO I=1,P_POINTS
          IF ( .NOT.(ICE_FRACT(I).GT.0.0 .AND. .NOT.LAND_MASK(I)) )
     +      QSTAR_LEAD(I) = 1.0E30
        ENDDO
      ENDIF                 ! End of IF (GATHER) THEN... ELSE.          
C-----------------------------------------------------------------------
CL  2.2  Set components of ocean surface current.
C-----------------------------------------------------------------------
      CALL UV_TO_P(U_0,U_0_P,U_POINTS,P_POINTS,ROW_LENGTH,U_ROWS)
      CALL UV_TO_P(V_0,V_0_P,U_POINTS,P_POINTS,ROW_LENGTH,U_ROWS)
C
C-----------------------------------------------------------------------
CL  3.  Loop round gridpoints to be processed, performing calculations
CL      BEFORE call to FCDCH which necessitates splitting of loop.
C
C       Note that there are a number of redundant operations in this
C       loop, mainly setting "default" values before entering an IF-
C       block.  These are necessary to prevent the flow control
C       becoming too complicated for the compiler to vectorize, and
C       should therefore be left in.
C
C       Note too that this section has been split into several small
C       loops instead of being kept as one larger loop.  This is also
C       done to enable the compiler to vectorize.
C-----------------------------------------------------------------------
C
C-----------------------------------------------------------------------
CL  3.1 Fix roughness lengths for the various surface types.
C-----------------------------------------------------------------------
      DO I = 1,P_POINTS
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
CL  3.1.1 Liquid sea. Overwrite sea-ice and land in 3.1.2, 3.1.3.
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        Z0M(I) = Z0MSEA(I)                                     ! P243.B5
        Z0H(I) = Z0HSEA                                        !    "
        Z0M_EFF(I) = Z0MSEA(I)
        Z0F(I) = Z0FSEA                                        !    "
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
CL  3.1.2 Sea ice: Z0MIZ set on all points for input to FCDCH routine
CL        in CD_MIZ,CH_MIZ calculations. Similarily Z0HSEA,Z0FSEA for
CL        CD_LEAD,CH_LEAD calculations. Z0SICE for CD,CH over sea-ice.
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        CD_LEAD(I) = Z0MIZ
        Z0HS(I) = Z0HSEA
        Z0FS(I) = Z0FSEA
        IF (ICE_FRACT(I).GT.0.0 .AND. .NOT.LAND_MASK(I)) THEN
          Z0M(I) = Z0SICE                                      ! P243.B4
          Z0H(I) = Z0SICE                                      !    "
          Z0M_EFF(I) = Z0SICE
          Z0F(I) = Z0SICE                                      !    "
        ENDIF
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
CL  3.1.2a Specify blending height for all points. Set to minimum value
CL         so that LAMBDA_EFF = LAMBDA over the sea in KMKH.
CL         This avoids having to pass land-sea mask into KMKH.
CL         Also set the wind profile factor to the default value of 1.0
CL         for non-land points.
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        H_BLEND(I) = H_BLEND_MIN
        WIND_PROFILE_FACTOR(I) = 1.0
C
      ENDDO
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
CL  3.1.3 Land.  Reduce roughness if there is any snow lying.
CL        Eqns P243.B1, B2.
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
      DO L = 1,LAND_PTS
        I = LAND_INDEX(L) - (P1-1)
        IF (LYING_SNOW(I) .LT. 5.0E3) THEN ! Only reduce non-orographic
C                                          ! roughness for land points
C                                          ! without permanent snow.
C
          Z0 = Z0V(I) - 4.0E-4 * LYING_SNOW(I)
          ZETA1 = MIN( 5.0E-4 , Z0V(I) )
          Z0M(I) = MAX( ZETA1 , Z0 )
          Z0H(I) = MIN( Z0V(I) , Z0M(I) )
          Z0F(I) = Z0H(I)
        ELSE                 ! for permanent land-ice Z0V is appropriate
          Z0M(I) = Z0V(I)         ! P243.B1, based on P243.B2 (2nd case)
          Z0H(I) = Z0V(I)         !    "   ,   "    "    "    ( "    " )
          Z0F(I) = Z0V(I)         !    "   ,   "    "    "    ( "    " )
        ENDIF
C
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
CL  3.1.4 Orographic roughness. Calculate Z0M_EFF in neutral case.
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        IF (L_Z0_OROG) THEN
C
C         ! Set blending height, effective roughness length and
C         ! wind profile factor at land points.
C

          IF (SIL_OROG(L) .EQ. RMDI) THEN
             ZETA1 = 0.0
          ELSE
             ZETA1 = 0.5 * OROG_DRAG_PARAM * SIL_OROG(L)
          ENDIF
          ZETA2 = LOG ( 1.0 + Z1(I)/Z0M(I) )
          ZETA3 = 1.0 / SQRT ( ZETA1/(VKMAN*VKMAN) + 1.0/(ZETA2*ZETA2) )
          ZETA2 = 1.0 / EXP(ZETA3)
          H_BLEND(I) = MAX ( Z1(I) / (1.0 - ZETA2) ,
     &                       HO2R2_OROG(L) * SQRT(2.0) )
          H_BLEND(I) = MIN ( H_BLEND_MAX, H_BLEND(I) )

          Z0M_EFF(I) = H_BLEND(I) / EXP ( VKMAN / SQRT ( ZETA1 +
     &                 (VKMAN / LOG ( H_BLEND(I) / Z0M(I) ) ) *
     &                 (VKMAN / LOG ( H_BLEND(I) / Z0M(I) ) ) ) )

          WIND_PROFILE_FACTOR(I) = LOG ( H_BLEND(I) / Z0M_EFF(I) ) /
     &                             LOG ( H_BLEND(I) / Z0M(I) )

        ELSE ! Orographic roughness not represented so
C            ! leave blending height and wind profile factor at their
C            ! default values and set effective roughness length to its
C            ! value based on land type.
C
          Z0M_EFF(I) = Z0M(I)
        ENDIF

      ENDDO
C-----------------------------------------------------------------------
CL  3.2 Calculate buoyancy parameters for the lowest model level.
C-----------------------------------------------------------------------
      DO I=1,P_POINTS
        IF (L_BL_LSPICE) THEN
          QW_1(I) = Q_1(I) + QCL_1(I)                         ! P243.10
        ELSE
          QW_1(I) = Q_1(I) + QCL_1(I) + QCF_1(I)              ! P243.10
        ENDIF
        BETAT = 1.0 / T_1(I)                         ! P243.19 (1st eqn)
        BETAQ = C_VIRTUAL /
     +     ( 1.0 + C_VIRTUAL*Q_1(I) - QCL_1(I) - QCF_1(I) )
C                                                  ... P243.19 (2nd eqn)
C
       IF (TL_1(I).GT.TM.OR.L_BL_LSPICE) THEN
          ALPHAL = (EPSILON * LC * QSL(I)) / (R * TL_1(I) * TL_1(I))
C                                       ... P243.20 (Clausius-Clapeyron)
C
          AL = 1.0 / (1.0 + LCRCP*ALPHAL)                      ! P243.21
          BETAC = CF_1(I) * AL * (LCRCP*BETAT - ETAR*BETAQ)
C                                                  ... P243.19 (3rd eqn)
C
        ELSE
          ALPHAL = (EPSILON * LS * QSL(I)) / (R * TL_1(I) * TL_1(I))
C                                       ... P243.20 (Clausius-Clapeyron)
C
          AL = 1.0 / (1.0 + LSRCP*ALPHAL)                      ! P243.21
          BETAC = CF_1(I) * AL * (LSRCP*BETAT - ETAR*BETAQ)
C                                                  ... P243.19 (3rd eqn)
C
        ENDIF
        BT_1(I) = BETAT - ALPHAL*BETAC               ! P243.18 (1st eqn)
        BQ_1(I) = BETAQ + BETAC                      ! P243.18 (2nd eqn)
        BF_1(I) = BETAQ *EPSILON*ETAR               ! P243.18 (2nd eqn)
      ENDDO
C-----------------------------------------------------------------------
CL  3.3 Calculate temperature (strictly, liquid/ice static energy) and
CL      humidity jumps, and wind shear, across the surface layer.
CL      Separate values are required for the leads and ice fractions
CL      of sea-ice grid-squares.
C-----------------------------------------------------------------------
      IF (GATHER) THEN
        DO I=1,P_POINTS
          DTEMP(I) = TL_1(I) - TSTAR_NL(I)
     &                 + GRCP * ( Z1(I) + Z0M_EFF(I) - Z0H(I) )
C                                                             ! P243.118
          DQ(I) = QW_1(I) - QSTAR(I)                          ! P243.119
          DTEMP_LEAD(I) = 1.0E30                ! Missing data indicator
          DQ_LEAD(I) = 1.0E30                   ! Missing data indicator
          USHEAR = U_1_P(I) - U_0_P(I)
          VSHEAR = V_1_P(I) - V_0_P(I)
          VSHR2 = MAX (1.0E-6 , USHEAR*USHEAR + VSHEAR*VSHEAR)
          VSHR(I) = SQRT(VSHR2)
C                                ... P243.120 (previous 4 lines of code)
        ENDDO
C-----------------------------------------------------------------------
CL 3.3.1 Calculate leads values by looping round sea-ice points only.
C        Avoids an if test in the above loop, so code can run faster.
C-----------------------------------------------------------------------
CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
        DO SI=1,NSICE
          I = SICE_INDEX(SI)
          DTEMP_LEAD(I) = TL_1(I)-TFS + GRCP*(Z1(I)+Z0MSEA(I)-Z0HS(I))
          DQ_LEAD(I) = QW_1(I) - QSTAR_LEAD(SI)
        ENDDO
      ELSE
      DO I=1,P_POINTS
        USHEAR = U_1_P(I) - U_0_P(I)
        VSHEAR = V_1_P(I) - V_0_P(I)
        VSHR2 = MAX (1.0E-6 , USHEAR*USHEAR + VSHEAR*VSHEAR)
        VSHR(I) = SQRT(VSHR2)
C                                ... P243.120 (previous 4 lines of code)
        DTEMP(I) = TL_1(I) - TSTAR_NL(I)
     &                 + GRCP * ( Z1(I) + Z0M_EFF(I) - Z0H(I) )
C                                                             ! P243.118
        DQ(I) = QW_1(I) - QSTAR(I)                            ! P243.119
        IF ( ICE_FRACT(I).GT.0.0 .AND. .NOT.LAND_MASK(I) ) THEN
          DTEMP_LEAD(I) = TL_1(I)-TFS + GRCP*(Z1(I)+Z0MSEA(I)-Z0HS(I))
          DQ_LEAD(I) = QW_1(I) - QSTAR_LEAD(I)
        ELSE
          DTEMP_LEAD(I) = 1.0E30          ! Missing data indicator
          DQ_LEAD(I) = 1.0E30             ! Missing data indicator
        ENDIF
C
      ENDDO
      ENDIF                   ! End of IF (GATHER) THEN... ELSE...
C-----------------------------------------------------------------------
CL  3.4 Evaporation over land surfaces without snow is limited by
CL      soil moisture availability and stomatal resistance.
C       Set FRACA (= fA in the documentation) according to P243.68,
C       PSIS according to P243.65, and RESFS (= fS) according to P243.75
C       and P243.61, using neutral-stability value of CH (as explained
C       in section (v) of the P243 documentation).
C-----------------------------------------------------------------------
      DO I=1,P_POINTS
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
CL  3.4.1 Set parameters (workspace) to values relevant for non-land
CL        points.  Only aerodynamic resistance applies.
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        FRACA(I) = 1.0
        PSIS(I) = 0.0
        RESFT(I) = 1.0
        RS(I) = 0.0
        RESFS(I) = 0.0
      ENDDO
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
CL  3.4.2 Over-write workspace for other points.
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
      DO L=1,LAND_PTS
        I = LAND_INDEX(L) - (P1-1)
C
C  Calculate the soil moisture availability factor, PSIS.
C
        SMCRIT = RHO_WATER * ROOTD(L) * (SMVCCL(L)-SMVCWT(L))
C                                                            ... P243.66
C
        PSIS(I) = 0.0
        IF (SMC(L).GE.SMCRIT .AND. SMCRIT.GT.0.0)
     +    PSIS(I) = 1.0
        IF (SMC(L).LT.SMCRIT .AND. SMC(L).GT.0.0)
     +    PSIS(I) = SMC(L)/SMCRIT
C
C  For snow-covered land or land points with negative moisture flux
C  set the fraction of the flux with only aerodynamic resistance to 1
C  (no surface/stomatal resistance to evap from snow or condensation).
C
        FRACA(I) = 1.0
C
C  When there is positive moisture flux from snow-free land, calculate
C  the fraction of the flux from the "canopy".
C
        IF (DQ(I).LT.0.0 .AND. LYING_SNOW(I).LE.0.0) FRACA(I) = 0.0
        IF (DQ(I).LT.0.0.AND.LYING_SNOW(I).LE.0.0.AND.CATCH(L).GT.0.0)
     +    FRACA(I) = CANOPY(L)/CATCH(L)
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
CL  3.4.2.1 Calculate neutral stability value of CH (CHN), as an
CL          approximation to CH.
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        ZETAM = LOG ( (Z1(I) + Z0M_EFF(I))/Z0M_EFF(I) )
        ZETAH = LOG ( (Z1(I) + Z0M_EFF(I))/Z0H(I) )
        CHN = (VKMAN/ZETAH) * (VKMAN/ZETAM) * WIND_PROFILE_FACTOR(I)
C                                 ... P243.41 (previous 3 lines of code)
C
        RS(I) = RESIST(L)
        RESFS(I) = PSIS(I) / ( 1.0 + CHN*VSHR(I)*RS(I) )
        RESFT(I) = FRACA(I) + (1.0 - FRACA(I)) * RESFS(I)
      ENDDO         ! Evaporation over land points only, section 3.4.2
C-----------------------------------------------------------------------
CL  3.5 Calculate bulk Richardson numbers for the surface layer.
CL      At sea-ice points RIB contains value for ice only (not leads).
CL      Initialise RIB_LEAD to RIB so that it contains sensible
CL      values at non sea ice points for the FCDCH calculation below.
C-----------------------------------------------------------------------
      IF (GATHER) THEN
        DO I=1,P_POINTS
          RIB(I) = G * Z1(I) *
     +                 ( BT_1(I)*DTEMP(I) + BQ_1(I)*RESFT(I)*DQ(I) ) /
     +                 ( VSHR(I)*VSHR(I) )
          RIB_LEAD(I) = RIB(I)
        ENDDO
C-----------------------------------------------------------------------
CL  3.5.1  Calculate bulk Richardson no. for leads at sea-ice points
CL         only.
C-----------------------------------------------------------------------
CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
        DO SI = 1,NSICE
          I = SICE_INDEX(SI)
          RIB_LEAD(I) = G * Z1(I) *
     +                      ( BT_1(I) * DTEMP_LEAD(I) +
     +                        BQ_1(I) * RESFT(I) * DQ_LEAD(I) ) /
     +                      ( VSHR(I) * VSHR(I) )
C                            ... P2430.2, for sea-ice leads.
        ENDDO
      ELSE
      DO I=1,P_POINTS
        RIB(I) = G * Z1(I) *
     +               ( BT_1(I)*DTEMP(I) + BQ_1(I)*RESFT(I)*DQ(I) ) /
     +               ( VSHR(I)*VSHR(I) )
C                            ... P243.43 (G times middle line is surface
C                                layer buoyancy difference, P243.25)
C
        RIB_LEAD(I) = RIB(I)
        IF ( ICE_FRACT(I).GT.0.0 .AND. .NOT.LAND_MASK(I) ) THEN
          RIB_LEAD(I) = G * Z1(I) *
     +                      ( BT_1(I) * DTEMP_LEAD(I) +
     +                        BQ_1(I) * RESFT(I) * DQ_LEAD(I) ) /
     +                      ( VSHR(I) * VSHR(I) )
C                            ... P2430.2, for sea-ice leads.
        ENDIF
      ENDDO
      ENDIF           ! End of IF (GATHER) THEN... ELSE...
C
C-----------------------------------------------------------------------
CL  3.6 Calculate stability corrected effective roughness length.
CL  Simple linear interpolation when RIB between 0 and RIB_CRIT (>0) for
CL  form drag term.
C-----------------------------------------------------------------------

      IF (L_Z0_OROG) THEN
CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
        DO L = 1,LAND_PTS
          I = LAND_INDEX(L) - (P1-1)

          IF (RIB(I) .GE. RI_CRIT) THEN
            Z0M_EFF(I) = Z0M(I)
          ELSEIF (RIB(I) .GT. 0.0 ) THEN
            IF (SIL_OROG(L) .EQ. RMDI) THEN
               ZETA1 = 0.0
            ELSE
               ZETA1 = 0.5 * OROG_DRAG_PARAM * SIL_OROG(L) *
     &                       ( 1.0 - RIB(I) / RI_CRIT )
            ENDIF
            Z0M_EFF(I) = H_BLEND(I) / EXP ( VKMAN / SQRT ( ZETA1 +
     &                   (VKMAN / LOG ( H_BLEND(I) / Z0M(I) ) ) *
     &                   (VKMAN / LOG ( H_BLEND(I) / Z0M(I) ) ) ) )

          ENDIF

          WIND_PROFILE_FACTOR(I) = LOG ( H_BLEND(I) / Z0M_EFF(I) ) /
     &                             LOG ( H_BLEND(I) / Z0M(I) )


        ENDDO
      ENDIF

C
C-----------------------------------------------------------------------
CL  3.7 Calculate CD, CH via routine FCDCH.
CL  Calculate CD_MIZ,CH_MIZ,CD_LEAD,CH_LEAD on full field then set
CL  non sea-ice points to missing data (contain nonsense after FCDCH)
C   Unlike the QSAT calculations above, arrays are not compressed to
C   sea-ice points for FCDCH. This is because it would require extra
C   work space and initial tests showed that with with the extra
C   compression calculations required no time was saved.
C   NB CD_LEAD stores Z0MIZ for calculation of CD_MIZ,CH_MIZ.
C-----------------------------------------------------------------------
C
      CALL FCDCH(RIB,CD_LEAD,CD_LEAD,CD_LEAD,Z1,WIND_PROFILE_FACTOR,
     &           P_POINTS,CD_MIZ,CH_MIZ,CD_STD,LTIMER)
C                                           ! Marginal Ice Zone.P2430.9
      CALL FCDCH(RIB_LEAD,Z0MSEA,Z0HS,Z0FS,Z1,WIND_PROFILE_FACTOR,
     &           P_POINTS,CD_LEAD,CH_LEAD,CD_STD,LTIMER)
C                                           ! Sea-ice leads.P2430.8
      CALL FCDCH(RIB,Z0M_EFF,Z0H,Z0F,Z1,WIND_PROFILE_FACTOR,
     &           P_POINTS,CD,CH,CD_STD,LTIMER)
      DO I=1,P_POINTS
C       IF ( an ordinary sea points (no sea-ice) or a land point)
        IF (.NOT.(ICE_FRACT(I).GT.0.0 .AND. .NOT.LAND_MASK(I)) ) THEN
          CD_MIZ(I) = 1.0E30
          CH_MIZ(I) = 1.0E30
          CD_LEAD(I) = 1.0E30
          CH_LEAD(I) = 1.0E30
          RIB_LEAD(I) = 1.0E30
        ENDIF
      ENDDO
C
C-----------------------------------------------------------------------
CL  4.  Loop round gridpoints to be processed, performing calculations
CL      AFTER call to FCDCH which necessitates splitting of loop.
C-----------------------------------------------------------------------
CL  4.1 Recalculate RESFS using "true" CH.
C-----------------------------------------------------------------------
CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
      DO L = 1,LAND_PTS
        I = LAND_INDEX(L) - (P1-1)
          RESFS(I) = PSIS(I) / ( 1.0 + CH(I)*VSHR(I)*RS(I) ) ! P243.75
          EPDT = -PSTAR(I)/(R*TSTAR(I))*CH(I)*VSHR(I)*DQ(I)*TIMESTEP
          IF (EPDT .GT. 0.0 .AND. LYING_SNOW(I) .LE. 0.0
     &        .AND. CATCH(L) .GT. 0.0) THEN
            FRACA(I) = CANOPY(L)/(CATCH(L)+EPDT)
            RESFT(I) = FRACA(I) + (1.0 - FRACA(I)) * RESFS(I)
          ENDIF ! Positive evaporation from snow-free land
      ENDDO ! Loop over land-points
C-----------------------------------------------------------------------
CL  4.D Call SFL_INT to calculate CDR10M, CHR1P5M and CER1P5M -
CL      interpolation coefficients used in SF_EVAP and IMPL_CAL to
CL      calculate screen temperature, specific humidity and 10m winds.
C-----------------------------------------------------------------------
C
      IF (SU10 .OR. SV10 .OR. SQ1P5 .OR. ST1P5) THEN

        CALL SFL_INT (
     &  P_POINTS,U_POINTS,RIB,Z1,Z0M,Z0M_EFF,Z0H,Z0F,CD_STD,CD,CH,
     &  RESFT,WIND_PROFILE_FACTOR,
     &  CDR10M,CHR1P5M,CER1P5M,
     +  SU10,SV10,ST1P5,SQ1P5,LTIMER
     + )
      ENDIF
C-----------------------------------------------------------------------
CL  4.2 Now that diagnostic calculations are over, update CD and CH
CL      to their correct values (i.e. gridsquare means).
C-----------------------------------------------------------------------
      DO I = 1,P_POINTS
        IF ( ICE_FRACT(I).GT.0.0 .AND. .NOT.LAND_MASK(I) ) THEN
          IF ( ICE_FRACT(I).LT. 0.7 ) THEN
            CD(I) = ( ICE_FRACT(I)*CD_MIZ(I) +
     +                (0.7-ICE_FRACT(I))*CD_LEAD(I) ) / 0.7  ! P2430.5
            CH(I) = ( ICE_FRACT(I)*CH_MIZ(I) +
     +                (0.7-ICE_FRACT(I))*CH_LEAD(I) ) / 0.7  ! P2430.4
          ELSE
            CD(I) = ( (1.0-ICE_FRACT(I))*CD_MIZ(I) +
     +                (ICE_FRACT(I)-0.7)*CD(I) ) / 0.3       ! P2430.7
            CH(I) = ( (1.0-ICE_FRACT(I))*CH_MIZ(I) +
     +                (ICE_FRACT(I)-0.7)*CH(I) ) / 0.3       ! P2430.7
          ENDIF
        ENDIF
C-----------------------------------------------------------------------
CL  4.3 Calculate the surface exchange coefficients RHOK(*).
C-----------------------------------------------------------------------
        RHOSTAR = PSTAR(I) / ( R*TSTAR(I) )
C                        ... surface air density from ideal gas equation
C
        RHOKM_1(I) = RHOSTAR * CD(I) * VSHR(I)                ! P243.124
        RHOKH_1(I) = RHOSTAR * CH(I) * VSHR(I)                ! P243.125
        RHOKEA(I) = RHOKH_1(I) * FRACA(I)                     ! P243.126
        RHOKESL(I) = RHOKH_1(I) * RESFS(I)                    ! P243.127
        RHOKES(I) = RHOKESL(I) * (1.0 - FRACA(I))             ! P243.128
        RHOKE(I) = RHOKEA(I) + RHOKES(I)                      ! P243.129
C
C     RHOSTAR * CD * VSHR stored for diagnostic output before
C     horizontal interpolation.
C
        RHO_CD_MODV1(I) = RHOKM_1(I)
C
C-----------------------------------------------------------------------
CL  4.4 Calculate the "explicit" fluxes of sensible heat and moisture.
C-----------------------------------------------------------------------
C
C  Firstly for land points, sea points with no sea-ice (ICE_FRACT=0.0)
C  and the sea-ice parts of the sea-ice points.
        FTL_1(I) = -RHOKH_1(I) * DTEMP(I)                     ! P243.134
        FQW_1(I) = -RHOKE(I) * DQ(I)                          ! P243.135
        EA(I) = -RHOKEA(I) * DQ(I)                            ! P243.136
        ES(I) = -RHOKES(I) * DQ(I)                            ! P243.137
        ESL(I) = -RHOKESL(I) * DQ(I)                          ! P243.138
        FTL_LEAD(I) = 1.0E30                   ! Missing data indicator
        FQW_LEAD(I) = 1.0E30                   ! Missing data indicator
C
C  Secondly for sea ice points.
        IF ( (ICE_FRACT(I).GT.0.0) .AND. .NOT.LAND_MASK(I) ) THEN
          FTL_LEAD(I) = -RHOKH_1(I)*DTEMP_LEAD(I) * (1.0-ICE_FRACT(I))
C                                                           ...P2430.11
          FTL_1(I) = FTL_LEAD(I) + FTL_1(I)*ICE_FRACT(I)   ! P2430.13/15
          FQW_LEAD(I) = -RHOKE(I)*DQ_LEAD(I) * (1.0-ICE_FRACT(I))
C                                                           ...P2430.12
          FQW_1(I) = FQW_LEAD(I) + FQW_1(I)*ICE_FRACT(I)   ! P2430.14/16
        ENDIF
C-----------------------------------------------------------------------
CL  4.4.1 Set indicator for unstable suface layer (buoyancy flux +ve.)
CL        if required by logical L_RMBL.
C-----------------------------------------------------------------------
        IF ( L_RMBL.AND.BT_1(I)*FTL_1(I)+BQ_1(I)*FQW_1(I).GT.0.0) THEN
          NRML(I) = 1
        ELSE
          NRML(I) = 0
        ENDIF
C-----------------------------------------------------------------------
CL  4.5 Multiply surface exchange coefficients that are on the P-grid
CL      by GAMMA(1).Needed for implicit calculations in P244(IMPL_CAL).
CL      RHOKM_1 dealt with in section 4.1 below.
C-----------------------------------------------------------------------
        RHOKE(I) = RHOKE(I) * GAMMA(1)
        RHOKEA(I) = RHOKEA(I) * GAMMA(1)
        RHOKES(I) = RHOKES(I) * GAMMA(1)
        RHOKESL(I) = RHOKESL(I) * GAMMA(1)
        RHOKH_1(I) = RHOKH_1(I) * GAMMA(1)
C-----------------------------------------------------------------------
CL  4.5.1 Calculate the standard deviations of layer 1 turbulent
CL        fluctuations of temperature and humidity using approximate
CL        formulae from first order closure.
C-----------------------------------------------------------------------
        VS = SQRT ( RHOKM_1(I)/RHOSTAR * VSHR(I) )
        VSF1_CUBED = 1.25 * ( Z1(I) + Z0M(I) ) * G *
     &                ( BT_1(I)*FTL_1(I) + BQ_1(I)*FQW_1(I) ) / RHOSTAR
C       !---------------------------------------------------------------
C       ! Only calculate standard deviations for unstable surface layers
C       !---------------------------------------------------------------
        IF (VSF1_CUBED .GT. 0.0) THEN
          WS1 = ( VSF1_CUBED + VS * VS * VS ) ** (1.0/3.0)
          T1_SD(I) = MAX ( 0.0 , 1.93 * FTL_1(I) / (RHOSTAR * WS1) )
          Q1_SD(I) = MAX ( 0.0 , 1.93 * FQW_1(I) / (RHOSTAR * WS1) )
        ELSE
          T1_SD(I) = 0.0
          Q1_SD(I) = 0.0
        ENDIF
C-----------------------------------------------------------------------
CL  4.6 For sea points, calculate the wind mixing energy flux and the
CL      sea-surface roughness length on the P-grid, using time-level n
CL      quantities.  Use CD_LEAD at sea-ice points.
C-----------------------------------------------------------------------
        IF (.NOT.LAND_MASK(I)) THEN
          TAU = RHOKM_1(I) * VSHR(I)                          ! P243.130
          IF (ICE_FRACT(I) .GT. 0.0)
     &      TAU = RHOSTAR * CD_LEAD(I) * VSHR(I) * VSHR(I)
          IF (SFME) FME(I) = (1.0-ICE_FRACT(I)) * TAU * SQRT(TAU/RHOSEA)
C                                                             ! P243.96
          Z0MSEA(I) = MAX ( Z0HSEA ,
     &                      (CHARNOCK/G) * (TAU / RHOSTAR) )
C                                         ... P243.B6 (Charnock formula)
C                      TAU/RHOSTAR is "mod VS squared", see eqn P243.131
C
        ENDIF ! of IF (.NOT. LAND_MASK), land-points done in next loop.
      ENDDO ! Loop over points for sections 4.2 - 4.6
      DO L=1,LAND_PTS
      I = LAND_INDEX(L) - (P1-1)
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
CL  4.7 Set Z0MSEA to Z0V, FME to zero for land points.
C   (Former because UM uses same storage for Z0V
C   and Z0MSEA.)
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Z0MSEA(I) = Z0V(I)
      IF (SFME) FME(I) = 0.0
      ENDDO ! Loop over points for section 4.7
C
C-----------------------------------------------------------------------
CL  5.  Calculate "explicit" surface fluxes of momentum (on UV-grid).
C-----------------------------------------------------------------------
CL  5.1 Interpolate exchange coefficient to UV-grid, then mutiply
CL      by GAMMA(1) to be passed to subroutine IMPL_CAL (P244) which
CL      only uses RHOKM_1 when mulitplied by GAMMA(1).
C-----------------------------------------------------------------------
C
C  PSIS used purely as spare workspace here.
C
! RHOKM_1 contains incorrect data in halos. The P_TO_UV can interpolate 
! this into the real data, so first we must update east/west halos.
      CALL SWAPBOUNDS(RHOKM_1,ROW_LENGTH,U_POINTS/ROW_LENGTH,1,0,1)  

      CALL P_TO_UV(RHOKM_1,PSIS,P_POINTS,U_POINTS,ROW_LENGTH,P_ROWS)
      DO I=1,U_POINTS-2*ROW_LENGTH
        J = I+ROW_LENGTH
        RHOKM_1(J) = PSIS(I)
        TAUX_1(J) = RHOKM_1(J) * ( U_1(J) - U_0(J) )         ! P243.132
        TAUY_1(J) = RHOKM_1(J) * ( V_1(J) - V_0(J) )         ! P243.133
        RHOKM_1(J) = GAMMA(1) * RHOKM_1(J)
      ENDDO
C-----------------------------------------------------------------------
CL  5.2 Set first and last rows to "missing data indicator".
C-----------------------------------------------------------------------
      IF (attop) THEN
        DO I=1,ROW_LENGTH
          RHOKM_1(I) = 1.0E30
          TAUX_1(I) = 1.0E30
          TAUY_1(I) = 1.0E30
        ENDDO
      ENDIF

      IF (atbase) THEN
        DO I= (U_ROWS-1)*ROW_LENGTH + 1 , U_ROWS*ROW_LENGTH
          RHOKM_1(I) = 1.0E30
          TAUX_1(I) = 1.0E30
          TAUY_1(I) = 1.0E30
        ENDDO
      ENDIF
C-----------------------------------------------------------------------
CL  5.D Interpolate CDR10M to UV-grid.
C-----------------------------------------------------------------------
! CDR10M contains incorrect data in halos. The P_TO_UV can interpolate  
! this into the real data, so first we must update east/west halos.
      CALL SWAPBOUNDS(CDR10M,ROW_LENGTH,U_POINTS/ROW_LENGTH,1,0,1)  

      IF (SU10 .OR. SV10) THEN
        CALL P_TO_UV(CDR10M,PSIS,P_POINTS,U_POINTS,ROW_LENGTH,P_ROWS)
        DO I=1,U_POINTS-2*ROW_LENGTH
          J = I + ROW_LENGTH
          CDR10M(J) = PSIS(I)
        ENDDO
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
CL  5.D.1 Set first and last rows to "missing data indicator".
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        IF (attop) THEN
          DO I=1,ROW_LENGTH
            CDR10M(I) = 1.0E30
          ENDDO
        ENDIF

        IF (atbase) THEN
          DO I= (U_ROWS-1)*ROW_LENGTH + 1 , U_ROWS*ROW_LENGTH
            CDR10M(I) = 1.0E30
          ENDDO
        ENDIF
      ENDIF
    6 CONTINUE   ! Branch for error exit.
      IF (LTIMER) THEN
        CALL TIMER('SFEXCH  ',4)
      ENDIF
      RETURN
      END
