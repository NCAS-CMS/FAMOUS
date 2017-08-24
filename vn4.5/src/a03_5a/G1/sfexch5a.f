C *****************************COPYRIGHT******************************
C (c) CROWN COPYRIGHT 1996, METEOROLOGICAL OFFICE, All Rights Reserved.
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
C*LL  SUBROUTINE SF_EXCH------------------------------------------------
CLL
CLL  Purpose: Calculate coefficients of turbulent exchange between
CLL           the surface and the lowest atmospheric layer, and
CLL           "explicit" fluxes between the surface and this layer.
CLL
CLL  Suitable for Single Column use.
CLL
CLL          Canopy evaporation made implicit
CLL     with respect to canopy water content (requiring TIMESTEP to be
CLL     passed in).
CLL
CLL  Model            Modification history:
CLL version  Date
CLL   4.1  07/05/96   New deck. M.J.Woodage
CLL   4.2   Oct. 96   T3E migration - *DEF CRAY removed
CLL                                     S J Swarbrick
!LL   4.3  14/01/97   MPP code : Corrected setting of polar rows
!LL                                                     P.Burton
CLL   4.3  15/05/97   By-pass call to SF_STOM when land points=0 to
CLL                   prevent occasional failures with MPP. R.Rawlins
CLL  4.3  09/06/97  Add swapbounds for CDR10M.  D.Sexton/RTHBarnes
CLL  4.4  08/09/97  L_BL_LSPICE specifies mixed phase precipitation
CLL                 scheme                     D.Wilson
CLL  4.5  20/08/98    Option to include a thermal plant canopy 
CLL                                            M.Best
CLL  4.5    Jul. 98  Kill the IBM specific lines (JCThil)
CLL
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
     & P_POINTS,LAND_PTS,U_POINTS,ROW_LENGTH,P_ROWS,U_ROWS
     &,LAND_INDEX,P1,GATHER
     &,AK_1,BK_1
     &,CANOPY,CATCH,CO2,CF_1,SM_LEVELS,DZSOIL,HCONS,F_TYPE
     &,HT,LAI,PAR,GPP,NPP,RESP_P
     &,ICE_FRACT,LAND_MASK,LYING_SNOW
     &,PSTAR,Q_1,QCF_1,QCL_1,RADNET_C,GC,RESIST,ROOTD,SMC
     &,SMVCCL,SMVCWT
     &,T_1,TIMESTEP,TI,TS1,TSTAR
     &,U_1,V_1,U_1_P,V_1_P,U_0,V_0,V_ROOT,V_SOIL
     &,VFRAC,Z0V,SIL_OROG,Z1,CANCAP,Z0MSEA,HO2R2_OROG     
     &, ALPHA1,ASHTF,BQ_1,BT_1,BF_1,CD,CH
     &,EPOT,FQW_1,FSMC,FTL_1,E_SEA,H_SEA,TAUX_1,TAUY_1,QW_1
     &,FRACA,RESFS,F_SE,RESFT,RHOKE,RHOKH_1,RHOKM_1
     &,RHOKPM,RHOKPM_POT
     &,RIB,TL_1,VSHR,Z0H,Z0M,Z0M_EFF,H_BLEND
     &,T1_SD,Q1_SD
     &,RHO_CD_MODV1
     &,CDR10M,CHR1P5M,CER1P5M,FME
     &,SU10,SV10,SQ1P5,ST1P5,SFME
     &,RHO_ARESIST,ARESIST,RESIST_B
     &,NRML
     &,L_Z0_OROG,L_RMBL,L_BL_LSPICE,ERROR,LTIMER
     &)
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
     & P_POINTS            ! IN Number of P-grid points to be processed.
     &,LAND_PTS            ! IN Number of land points to be processed.
     &,U_POINTS            ! IN Number of UV-grid points.
     &,ROW_LENGTH          ! IN No. of points in latitude row (inclusive
C                          !    of endpoints for ltd. area model).
     &,P_ROWS              ! IN Number of rows of data on P-grid.
     &,U_ROWS              ! IN Number of rows of data on UV-grid.
     &,LAND_INDEX(LAND_PTS)! IN Index for compressed land point array;
C                          !    ith element holds position in the FULL
C                          !    field of the ith land pt to be processed
     &,P1                  ! IN First P-point to be processed.

      LOGICAL
     & GATHER              ! IN If true then leads variables are comp-
C                          !    ressed for sea-ice calculations. This
C                          !    saves duplicating calculations if there
C                          !    are a relatively few of sea-ice points.
C                          !    Set to false for a limited area run
C                          !    with a high proportion of sea-ice.
!---------------------------------------------------------------------
! Extra variables for the interactive stomatal resistance model
!---------------------------------------------------------------------
      INTEGER
     & SM_LEVELS           ! IN Number of soil moisture levels
     &,F_TYPE(LAND_PTS)    ! IN Plant functional type:
C                          !     1 - Broadleaf Tree
C                          !     2 - Needleleaf Tree
C                          !     3 - C3 Grass
C                          !     4 - C4 Grass

      REAL
     & HT(LAND_PTS)        ! IN Canopy height (m).
     &,LAI(LAND_PTS)       ! IN Leaf area index.
     &,PAR(P_POINTS)       ! IN Photosynthetically active radiation
C                          !    (W/m2).
     &,GPP(LAND_PTS)       ! OUT Gross Primary Productivity
C                          !    (kg C/m2/s).
     &,NPP(LAND_PTS)       ! OUT Net Primary Productivity
C                          !    (kg C/m2/s).
     &,RESP_P(LAND_PTS)    ! OUT Plant respiration rate (kg C/m2/s).

      REAL
     & AK_1                ! IN Hybrid "A" for lowest model layer.
     &,BK_1                ! IN Hybrid "B" for lowest model layer.
     &,CANOPY(LAND_PTS)    ! IN Surface water (kg per sq metre).  F642.
     &,CATCH(LAND_PTS)     ! IN Surface capacity (max. surface water)
C                          !    (kg per sq metre).  F6416.
     &,CF_1(P_POINTS)      ! IN Cloud fraction for lowest atmospheric
C                          !    layer (decimal fraction).
     &,CO2                 ! IN CO2 mixing ratio (kg CO2/kg air).
     &,DZSOIL(SM_LEVELS)   ! IN Thicknesses of the soil layers (m).
     &,HCONS(LAND_PTS)     ! IN Soil thermal conductivity including
C                          !    the effects of water and ice (W/m/K).
     &,ICE_FRACT(P_POINTS) ! IN Fraction of gridbox which is sea-ice.
     &,LYING_SNOW(P_POINTS)! IN Lying snow amount (kg per sq metre).
     &,PSTAR(P_POINTS)     ! IN Surface pressure (Pascals).
     &,Q_1(P_POINTS)       ! IN Specific humidity for lowest atmospheric
C                          !    layer (kg water per kg air).
     &,QCF_1(P_POINTS)     ! IN Cloud ice for lowest atmospheric layer
C                          !    (kg water per kg air).
     &,QCL_1(P_POINTS)     ! IN Cloud liquid water for lowest atm layer
C                          !    (kg water per kg air).
     &,GC(LAND_PTS)        ! IN Interactive canopy conductance
C                          !    to evaporation (m/s)
     &,RESIST(LAND_PTS)    ! IN Fixed "stomatal" resistance
C                          !    to evaporation (s/m)
     &,ROOTD(LAND_PTS)     ! IN "Root depth" (metres).  F6412.
     &,SMC(LAND_PTS)       ! IN Soil moisture content (kg per sq m).
C                          !    F621.
     &,SMVCCL(LAND_PTS)    ! IN Critical volumetric SMC (cubic metres
C                          !    per cubic metre of soil).  F6232.
     &,SMVCWT(LAND_PTS)    ! IN Volumetric wilting point (cubic m of
C                          !    water per cubic m of soil).  F6231.
C
C    Note: (SMVCCL - SMVCWT) is the critical volumetric available soil
C          moisture content.                            ~~~~~~~~~
C
     &,STHU(LAND_PTS,SM_LEVELS)! IN Unfrozen soil moisture content of
C                         !    each layer as a fraction of
C                         !    saturation.
C
      REAL                !    (Split to avoid > 19 continuations.)
     & T_1(P_POINTS)      ! IN Temperature for lowest atmospheric layer
C                         !    (Kelvin).
     &,TIMESTEP           ! IN Timestep in seconds for EPDT calc.
     &,TI(P_POINTS)       ! IN Temperature of sea-ice surface layer (K).
     &,TS1(LAND_PTS)      ! IN Temperature of top soil layer (K)
     &,TSTAR(P_POINTS)    ! IN Mean gridsquare surface temperature (K).
     &,U_1(U_POINTS)      ! IN West-to-east wind component for lowest
C                         !    atmospheric layer (m/s).  On UV grid.
     &,V_1(U_POINTS)      ! IN South-to-north wind component for lowest
C                         !    atmospheric layer (m/s).  On UV grid.
     &,U_1_P(P_POINTS)    ! IN West-to-east wind component for lowest
C                         !    atmospheric layer (m/s).  On P grid.
C                         !    (Same as U_1 for Single Column Model.)
     &,V_1_P(P_POINTS)    ! IN South-to-north wind component for lowest
C                         !    atmospheric layer (m/s).  On P grid.
C                         !    (Same as V_1 for Single Column Model.)
     &,U_0(U_POINTS)      ! IN West-to-east component of ocean surface
C                         !    current (m/s; ASSUMED zero over land).
C                         !    UV grid.  F615.
     &,V_0(U_POINTS)      ! IN South-to-north component of ocean surface
C                         !    current (m/s; ASSUMED zero over land).
C                         !    UV grid.  F616.
     &,V_ROOT(LAND_PTS)   ! IN Volumetric soil moisture concentration
C                         !    in the rootzone (m3 H2O/m3 soil).
     &,V_SOIL(LAND_PTS)   ! IN Volumetric soil moisture concentration
C                         !    in the top soil layer (m3 H2O/m3 soil).
     &,VFRAC(LAND_PTS)    ! IN Vegetated fraction.
     &,Z0V(P_POINTS)      ! IN Vegetative roughness length (m).  F6418.
     &,SIL_OROG(LAND_PTS) ! IN Silhouette area of unresolved orography
C                         !    per unit horizontal area
     &,Z1(P_POINTS)       ! IN Height of lowest atmospheric level (m).
     &,HO2R2_OROG(LAND_PTS) ! IN Peak to trough height of unresolved
C                         !    orography devided by 2SQRT(2) (m).
      LOGICAL
     & LAND_MASK(P_POINTS) ! IN .TRUE. for land; .FALSE. elsewhere. F60.
     &,SU10                ! IN STASH flag for 10-metre W wind.
     &,SV10                ! IN STASH flag for 10-metre S wind.
     &,SQ1P5               ! IN STASH flag for 1.5-metre sp humidity.
     &,ST1P5               ! IN STASH flag for 1.5-metre temperature.
     &,SFME                ! IN STASH flag for wind mixing energy flux.
     +,L_RMBL                    ! IN T to use rapidly mixing boundary
C                                !    scheme in IMPL_CAL
     &,L_BL_LSPICE               ! IN
!                              TRUE  Use scientific treatment of mixed
!                                    phase precip scheme.
!                              FALSE Do not use mixed phase precip
!                                    considerations
     &,L_Z0_OROG           ! IN .TRUE. to use orographic roughness.
C
C  Modified (INOUT) variables.
C
      REAL
     & CANCAP(P_POINTS)   ! INOUT Volumetric heat capacity of
C                         !       vegetation canopy (J/Kg/m3).
     &,RADNET_C(P_POINTS) ! INOUT Adjusted net radiation for vegetation 
C                         !       over land (W/m2).
     &,Z0MSEA(P_POINTS)   ! INOUT Sea-surface roughness length for      
C                         !       momentum (m).  F617.
C
C  Output variables.
C
      REAL
     & ALPHA1(P_POINTS) ! OUT Gradient of saturated specific humidity
C                       !     with respect to temperature between the
C                       !     bottom model layer and the surface
     &,ASHTF(P_POINTS)  ! OUT Coefficient to calculate surface
C                       !     heat flux into soil or sea-ice (W/m2/K).
     &,BQ_1(P_POINTS)   ! OUT A buoyancy parameter for lowest atm level
C                       !     ("beta-q twiddle").
     &,BT_1(P_POINTS)   ! OUT A buoyancy parameter for lowest atm level.
C                       !     ("beta-T twiddle").
     &,BF_1(P_POINTS)   
!        OUT A buoyancy parameter for lowest atm level. 
!            ("beta-F twiddle").
     &,CD(P_POINTS)     ! OUT Bulk transfer coefficient for momentum.
     &,CH(P_POINTS)     ! OUT Bulk transfer coefficient for heat and/or
C                       !     moisture.
     &,CDR10M(U_POINTS) ! OUT Reqd for calculation of 10m wind (u & v).
C                       !     NBB: This is output on the UV-grid, but
C                       !     with the first and last rows set to a
C                       !     "missing data indicator".
C                       !     Sea-ice leads ignored. See 3.D.7 below.
     &,CHR1P5M(P_POINTS)! OUT Reqd for calculation of 1.5m temperature.
C                       !     Sea-ice leads ignored. See 3.D.7 below.
     &,CER1P5M(P_POINTS)! OUT Reqd for calculation of 1.5m sp humidity.
C                       !     Sea-ice leads ignored. See 3.D.7 below.
     &,RHO_CD_MODV1(P_POINTS)
C                       ! OUT rhostar*cD*vshr before horizontal
C                       !     interpolation output as a diagnostic.
      REAL              !     (Split to avoid > 19 continuations.)
     & EPOT(P_POINTS)   ! OUT potential evaporation on P-grid
C                       !      (kg/m2/s).
     &,FQW_1(P_POINTS)  ! OUT "Explicit" surface flux of QW (i.e.
C                       !      evaporation), on P-grid (kg/m2/s).
     &,FTL_1(P_POINTS)  ! OUT "Explicit" surface flux of TL = H/CP.
C                       !     (sensible heat / CP).
     &,FSMC(LAND_PTS)   ! OUT soil moisture availability.
     &,FRACA(P_POINTS)  ! OUT Fraction of surface moisture flux with
C                       !     only aerodynamic resistance.
     &,E_SEA(P_POINTS)  ! OUT Evaporation from sea times leads
C                       !     fraction (kg/m2/s). Zero over land.
     &,H_SEA(P_POINTS)  ! OUT Surface sensible heat flux over sea
C                       !     times leads fraction (W/m2).
C                       !     Zero over land.
     &,TAUX_1(U_POINTS) ! OUT "Explicit" x-component of surface
C                       !     turbulent stress; on UV-grid; first and
C                       !     last rows set to a "missing data
C                       !     indicator". (Newtons per square metre)
     &,TAUY_1(U_POINTS) ! OUT "Explicit" y-component of surface
C                       !     turbulent stress; on UV-grid; first and
C                       !     last rows set to a "missing data
C                       !     indicator". (Newtons per square metre)
     &,QW_1(P_POINTS)   ! OUT Total water content of lowest
C                       !     atmospheric layer (kg per kg air).
     &,RESFS(P_POINTS)  ! OUT Combined soil, stomatal and aerodynamic
C                       !     resistance factor = PSIS/(1+RS/RA) for
C                       !     fraction (1-FRACA)
     &,F_SE(P_POINTS)   ! OUT Fraction of the evapotranspiration which
C                       !     is bare soil evaporation.
     &,RESFT(P_POINTS)  ! OUT Total resistance factor
C                       !     FRACA+(1-FRACA)*RESFS.
C
      REAL ! Surface exchange coefficients;passed to subroutine IMPL_CAL
     & RHOKE(P_POINTS)   ! OUT For FQW, then *GAMMA(1) for implicit calc
     &,RHOKH_1(P_POINTS) ! OUT For FTL,then *GAMMA(1) for implicit calcs
     &,RHOKM_1(U_POINTS) ! OUT For momentum, then *GAMMA(1) for implicit
C                        !     calculations. NBB: This is output on the
C                        !     UV-grid, but with the first and last
C                        !     rows set to a "missing data indicator".
     &,RHOKPM(P_POINTS)  ! OUT NB NOT * GAMMA for implicit calcs.
     &,RHOKPM_POT(P_POINTS)
C                         ! OUT Surface exchange coeff. for
C                               potential evaporation.
     &,Z0M_EFF(P_POINTS)  ! OUT Effective roughness length for momentum
     &,H_BLEND(P_POINTS)  ! OUT Blending height
     &,T1_SD(P_POINTS)    ! OUT Standard deviation of turbulent
C                         !     fluctuations of surface layer
C                         !     temperature (K).
     &,Q1_SD(P_POINTS)    ! OUT Standard deviation of turbulent
C                         !     fluctuations of surface layer
C                         !     specific humidity (kg/kg).
     &,RIB(P_POINTS)     ! OUT Bulk Richardson number for lowest layer.
     &,TL_1(P_POINTS)    ! OUT Liquid/frozen water temperature for
C                        !     lowest atmospheric layer (K).
     &,VSHR(P_POINTS)    ! OUT Magnitude of surface-to-lowest-lev. wind
     &,Z0H(P_POINTS)     ! OUT Roughness length for heat and moisture m
     &,Z0M(P_POINTS)     ! OUT Roughness length for momentum (m).
     &,FME(P_POINTS)     ! OUT Wind mixing energy flux (Watts/sq m).
     &,RHO_ARESIST(P_POINTS)  ! OUT, RHOSTAR*CD_STD*VSHR  for SCYCLE
     &,ARESIST(P_POINTS)      ! OUT, 1/(CD_STD*VSHR)      for SCYCLE
     &,RESIST_B(P_POINTS)     ! OUT, (1/CH-1/CD_STD)/VSHR for SCYCLE
C
      INTEGER
     & NRML(P_POINTS)    ! OUT 1 if surface layer unstable, else 0.
     &,ERROR             ! OUT 1 if grid definition faulty; else 0.
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
C RHO_WATER removed to avoid clash with declaration in C_DENSTY
C J.Smith 28/06/95
      REAL OMEGA1,RHO_SNOW,DEFF_SNOW,SNOW_HCON,SNOW_HCAP
      INTEGER PSOIL
      PARAMETER (
     + PSOIL=4                  ! No. of soil layers (must = NSOIL).
     +,OMEGA1=3.55088E-4        ! Tunable characteristic freq (rad/s).
     +,RHO_SNOW=250.0           ! Density of lying snow (kg per m**3).
     +,DEFF_SNOW=0.1            ! Depth of `effective' snow surface
C                               ! layer (m).
     +,SNOW_HCON=0.265          ! Thermal conductivity of lying snow
C                               ! (Watts per m per K).
     +,SNOW_HCAP=0.63E6         ! Thermal capacity of lying snow
C                               ! (J/K/m3)
     +)
      REAL KAPPAI
      PARAMETER (
     + KAPPAI=2.09          ! Thermal conductivity of sea-ice (W per
C                           ! m per K).
     +)
      REAL DE
      PARAMETER (
     + DE = 0.1             ! Effective thickness of sea-ice surface
C                           ! layer (m).
     +)
C-----------------------------------------------------------------------
C*L-----------COMDECK C_SICEHC FOR SUBROUTINE IMPL_CAL----------
C AI  = reciprocal effective areal heat capacity of sea-ice,
C          ( 1 / (J per sq m per K)).
      REAL AI

      PARAMETER(AI  = 4.8E-6)
C*----------------------------------------------------------------------

C
C   (3) Derived local parameters.
C
      REAL ETAR,GRCP,LCRCP,LFRCP,LS,LSRCP,H_BLEND_MIN,H_BLEND_MAX

      PARAMETER (
     & ETAR=1./(1.-EPSILON)  ! Used in calc of buoyancy parameter BETAC.
     &,GRCP=G/CP             ! Used in calc of dT across surface layer.
     &,LCRCP=LC/CP           ! Evaporation-to-dT conversion factor.
     &,LFRCP=LF/CP           ! Freezing-to-dT conversion factor.
     &,LS=LF+LC              ! Latent heat of sublimation.
     &,LSRCP=LS/CP           ! Sublimation-to-dT conversion factor.
     &,H_BLEND_MIN=0.0       ! Minimum blending height.
     &,H_BLEND_MAX=1000.0    ! Maximum blending height (m).
     &)
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
      EXTERNAL SF_ROUGH,SF_RIB,FCDCH,QSAT,SFL_INT,SF_FLUX,SF_STOM
     &,QSAT_WAT
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
     & CD_LEAD(P_POINTS)  ! Bulk transfer coefficient for momentum
C                         !  over sea-ice leads.Missing data over non
C                         !  sea-ice points.(Temporary store for Z0MIZ)
     &,CD_MIZ(P_POINTS)   ! Bulk transfer coefficient for momentum
C                         !  over the sea-ice Marginal Ice Zone.
C                         !  Missing data indicator over non sea-ice.
     &,CH_LEAD(P_POINTS)  ! Bulk transfer coefficient for heat and
C                         !  or moisture over sea ice leads.
C                         !  Missing data indicator over non sea-ice.
     &,CH_MIZ(P_POINTS)   ! Bulk transfer coefficient for heat and
C                         !  or moisture over the Marginal Ice Zone.
C                         !  Missing data indicator over non sea-ice.
     &,CD_STD(P_POINTS)   ! Local drag coefficient for
C                         !  calculation of interpolation coefficients
     &,DQ(P_POINTS)       ! Sp humidity difference between surface
C                         !  and lowest atmospheric level (Q1 - Q*).
C                         !  Holds value over sea-ice where ICE_FRACT
C                         !  >0 i.e. Leads contribution not included.
     &,DQI(P_POINTS)        
!        Ice water difference between surface
!        and lowest atmospheric level (Q1 - Q*).
!        Holds value over sea-ice where ICE_FRACT
!        >0 i.e. Leads contribution not included.
     &,DQ_LEAD(P_POINTS)  ! DQ for leads fraction of gridsquare.
C                         !  Missing data indicator over non sea-ice.
     &,DQI_LEAD(P_POINTS)   
!        DQI for leads fraction of gridsquare.
!        Missing data indicator over non sea-ice.
     &,DTEMP(P_POINTS)    ! Liquid/ice static energy difference
C                         !  between surface and lowest atmospheric
C                         !  level, divided by CP (a modified
C                         !  temperature difference).
C                         !  Holds value over sea-ice where ICE_FRACT
C                         !  >0 i.e. Leads contribution not included.
     &,DTEMP_LEAD(P_POINTS) ! DTEMP for leads fraction of gridsquare.
C                           !  Missing data indicator over non sea-ice.
     &,EPDT(P_POINTS)     ! "Potential" Evaporation * Timestep
     &,NL0(LAND_PTS)      ! Nitrogen concentration of the top leaf
C                         ! (kg N/kg C).
     &,PSIS(P_POINTS)     ! Soil moisture availability factor.
     &,PSTAR_ICE(P_POINTS)! Surface pressure over sea ice (Pa).
     &,QS1(P_POINTS)        ! Sat. specific humidity qsat(TL_1,PSTAR)
     &,QSL(P_POINTS)      ! Saturated sp humidity at liquid/ice
C                         !  temperature and pressure of lowest
C                         !  atmospheric level.
     &,QSTAR(P_POINTS)    ! Surface saturated sp humidity. Holds
C                         !  value over sea-ice where ICE_FRACT > 0.
C                         !  i.e. Leads contribution not included.
     &,QSTAR_LEAD(P_POINTS) ! QSTAR for sea-ice leads.
C                         ! Missing data indicator over non sea-ice.
     &,RHOSTAR(P_POINTS)  ! Surface air density in kg per cubic metre.
     &,RIB_LEAD(P_POINTS) ! Bulk Richardson no. for sea-ice leads at
C                         ! lowest layer. At non sea-ice points holds
C                         ! RIB for FCDCH calculation, then set to
C                         ! to missing data indicator.
     &,RA(P_POINTS)       ! Aerodynamic resistance.
     &,ROOT(LAND_PTS)     ! Root biomass (kg C/m2).
     &,TSTAR_NL(P_POINTS) ! TSTAR No Leads: surface temperature
C                         ! over sea-ice fraction of gridsquare.
C                         ! =TSTAR over non sea-ice points.
     &,U_0_P(P_POINTS)    ! West-to-east component of ocean surface
C                         ! current (m/s; zero over land if U_0 OK).
C                         ! P grid.  F615.
     &,V_0_P(P_POINTS)    ! South-to-north component of ocean surface
C                         ! current (m/s; zero over land if V_0 OK).
C                         ! P grid.  F616.
     &,WIND_PROFILE_FACTOR(P_POINTS)
C                         ! For transforming effective surface transfer
C                         ! coefficients to those excluding form drag.

     &,Z0F(P_POINTS)      ! Roughness length for free-convective heat
C                         ! and moisture transport.
     &,Z0FS(P_POINTS)     ! Roughness length for free-convective heat
C                         ! and moisture transport over sea.
     &,Z0HS(P_POINTS)     ! Roughness length for heat and moisture
C                         ! transport over sea.
C
C  Workspace (reqd for compression).
      INTEGER
     & SICE_INDEX(P_POINTS) ! Index vector for gather to sea-ice points
      LOGICAL ITEST(P_POINTS)  ! Used as 'logical' for compression.
C*
C
C   (b) Scalars.
C
      INTEGER
     & I           ! Loop counter (horizontal field index).
     &,J           ! Offset counter within I-loop.
     &,K           ! Offset counter within I-loop.
     &,L           ! Loop counter (land point field index).
     &,N           ! Loop counter (land point field index).
     &,NSICE       ! Number of sea-ice points.
     &,SI          ! Loop counter (sea-ice field index).
      REAL
     & TAU         ! Magnitude of surface wind stress over sea.
     &,VS          ! Surface layer friction velocity
     &,VSF1_CUBED  ! Cube of surface layer free convective scaling
C                  ! velocity
     &,WS1         ! Turbulent velocity scale for surface layer

!-------------------------------------------------------------------
! Extra work variables for the canopy (stomatal) conductance model.
!-------------------------------------------------------------------
      LOGICAL
     & INT_STOM              ! T for interactive stomatal resistance.
      PARAMETER (INT_STOM=.TRUE.)

C
C-----------------------------------------------------------------------
CL  0.  Check that the scalars input to define the grid are consistent.
C-----------------------------------------------------------------------
C
      IF (LTIMER) THEN
        CALL TIMER('SFEXCH  ',3)
      ENDIF

      ERROR=0
      IF ( U_ROWS .NE. (P_ROWS+1) .OR.
     &     U_POINTS .NE. (U_ROWS*ROW_LENGTH) .OR.
     &     P_POINTS .NE. (P_ROWS*ROW_LENGTH) .OR.
     &     LAND_PTS .GT.  P_POINTS )  THEN
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
          QSTAR_LEAD(I) = 1.0E30                 ! Missing data indicato
          QSTAR(I) = AK_1 + BK_1*PSTAR(I)
        ENDDO
        IF (NSICE.GT.0) THEN
CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
          DO SI = 1,NSICE
            I = SICE_INDEX(SI)
            TSTAR_NL(I) = (TSTAR(I)-(1.0-ICE_FRACT(I)) *TFS)
     &                / ICE_FRACT(I)                          ! P2430.1
            PSIS(SI) = TFS
            PSTAR_ICE(SI) = PSTAR(I)
          ENDDO
        ENDIF
        IF (L_BL_LSPICE) THEN
          CALL QSAT_WAT(QSL,TL_1,QSTAR,P_POINTS)
        ELSE
          CALL QSAT(QSL,TL_1,QSTAR,P_POINTS)
        ENDIF

        CALL QSAT(QSTAR,TSTAR_NL,PSTAR,P_POINTS)

C            ...values at sea-ice points contain ice contribution only
        IF (NSICE.GT.0) CALL QSAT(QSTAR_LEAD,PSIS,PSTAR_ICE,NSICE)
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
CL       Also initialise RIB to 0
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
     &                / ICE_FRACT(I)                          ! P2430.1
            PSIS(I) = TFS
          ENDIF
          RIB(I) = 0.0
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
     &      QSTAR_LEAD(I) = 1.0E30
        ENDDO
      ENDIF                ! End of IF (GATHER) THEN... ELSE.
C-----------------------------------------------------------------------
CL  2.2  Set components of ocean surface current.
C-----------------------------------------------------------------------
      CALL UV_TO_P(U_0,U_0_P,U_POINTS,P_POINTS,ROW_LENGTH,U_ROWS)
      CALL UV_TO_P(V_0,V_0_P,U_POINTS,P_POINTS,ROW_LENGTH,U_ROWS)
C
C-----------------------------------------------------------------------
CL  3. Calculation of transfer coefficients and surface layer stability
C-----------------------------------------------------------------------
C
C-----------------------------------------------------------------------
CL  3.1 Calculate neutral roughness lengths and blending height for
CL      surface
C-----------------------------------------------------------------------

      CALL SF_ROUGH (
     & P_POINTS,LAND_PTS,LAND_MASK,
     & P1,LAND_INDEX,
     & L_Z0_OROG,Z1,Z0MSEA,ICE_FRACT,
     & LYING_SNOW,Z0V,SIL_OROG,HO2R2_OROG,RIB,Z0M_EFF,Z0M,Z0H,
     & WIND_PROFILE_FACTOR,H_BLEND,CD_LEAD,Z0HS,Z0F,Z0FS,
     & LTIMER)


C-----------------------------------------------------------------------
CL  3.2 Calculate buoyancy parameters and bulk Richardson number for
CL      the lowest model level.
C-----------------------------------------------------------------------
C Calculate QSAT(TL1,P*)
C
      CALL QSAT(QS1,TL_1,PSTAR,P_POINTS)

      CALL SF_RIB (
     & P_POINTS,LAND_PTS,LAND_MASK,INT_STOM,
     & GATHER,P1,LAND_INDEX,
     & NSICE,SICE_INDEX,ICE_FRACT,
     & PSTAR,AK_1,BK_1,Q_1,QW_1,QCL_1,QCF_1,
     & CF_1,T_1,TL_1,QSL,QSTAR,QSTAR_LEAD,
     & QS1,TSTAR_NL,Z1,Z0M_EFF,Z0M,Z0H,Z0HS,Z0MSEA,
     & WIND_PROFILE_FACTOR,U_1_P,U_0_P,V_1_P,V_0_P,
     & ROOTD,SMVCCL,SMVCWT,SMC,V_SOIL,VFRAC,CANOPY,CATCH,
     & LYING_SNOW,GC,RESIST,RIB,RIB_LEAD,PSIS,VSHR,ALPHA1,
     & BT_1,BQ_1,BF_1,FRACA,RESFS,DQ,DQ_LEAD,DTEMP,
     & DTEMP_LEAD,L_BL_LSPICE,
     & LTIMER)

C-----------------------------------------------------------------------
CL  3.3 Calculate stability corrected effective roughness length.
CL  Simple linear interpolation when RIB between 0 and RIB_CRIT (>0) for
CL  form drag term.
C-----------------------------------------------------------------------


      CALL SF_ROUGH (
     & P_POINTS,LAND_PTS,LAND_MASK,
     & P1,LAND_INDEX,
     & L_Z0_OROG,Z1,Z0MSEA,ICE_FRACT,
     & LYING_SNOW,Z0V,SIL_OROG,HO2R2_OROG,RIB,Z0M_EFF,Z0M,Z0H,
     & WIND_PROFILE_FACTOR,H_BLEND,CD_LEAD,Z0HS,Z0F,Z0FS,
     & LTIMER)

C
C-----------------------------------------------------------------------
CL  3.4 Calculate CD, CH via routine FCDCH.
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
          CD_MIZ(I) = 1.E30
          CH_MIZ(I) = 1.E30
          CD_LEAD(I) = 1.E30
          CH_LEAD(I) = 1.E30
          RIB_LEAD(I) = 1.E30
        ENDIF
      ENDDO
C

C-----------------------------------------------------------------------
CL  4.  Loop round gridpoints to be processed, performing calculations
CL      AFTER call to FCDCH which necessitates splitting of loop.
C-----------------------------------------------------------------------
CL  4.1 Recalculate RESFS using "true" CH and EPDT
C-----------------------------------------------------------------------

CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
      DO L = 1,LAND_PTS
        I = LAND_INDEX(L) - (P1-1)
          EPDT(I) = -PSTAR(I)/(R*TSTAR(I))*CH(I)*VSHR(I)*DQ(I)*TIMESTEP
      ENDDO ! Loop over land-points

!-----------------------------------------------------------------------
! If the interactive surface resistance is requested call SF_STOM
!-----------------------------------------------------------------------
      IF (INT_STOM) THEN

!-----------------------------------------------------------------------
! Calculate the aerodynamic resistance
!-----------------------------------------------------------------------
        DO I=1,P_POINTS
          RA(I) = 1.0 / (CH(I) * VSHR(I))
        ENDDO

CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
        DO L = 1,LAND_PTS
          I = LAND_INDEX(L) - (P1-1)
!-----------------------------------------------------------------------
! For mesoscale model release assume uniform functional types and top
! leaf nitrogen concentrations. Assume that (fine) root biomass is
! equal to leaf biomass.
!-----------------------------------------------------------------------
          NL0(L) = 50.0E-3
          ROOT(L) = 0.04 * LAI(L)

        ENDDO ! Loop over land-points

        IF(LAND_PTS.GT.0) THEN    ! Omit if no land points

        CALL SF_STOM  (LAND_PTS,LAND_INDEX,P1,P_POINTS
     &,                F_TYPE,CO2,HT,PAR,LAI,NL0,PSTAR
     &,                Q_1,RA,ROOT,TSTAR,SMVCCL,V_ROOT,SMVCWT
     &,                VFRAC,GPP,NPP,RESP_P,GC,LTIMER,FSMC)

        ENDIF                     ! End test on land points

!-----------------------------------------------------------------------
! Convert carbon fluxes to gridbox mean values
!-----------------------------------------------------------------------
CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
        DO L = 1,LAND_PTS
          I = LAND_INDEX(L) - (P1-1)

            GPP(L) = VFRAC(L) * GPP(L)
            NPP(L) = VFRAC(L) * NPP(L)
            RESP_P(L) = VFRAC(L) * RESP_P(L)

        ENDDO ! Loop over land-points

      ENDIF  ! INT_STOM


      CALL SF_RESIST (
     & P_POINTS,LAND_PTS,LAND_MASK,INT_STOM,
     & P1,LAND_INDEX,
     & ROOTD,SMVCCL,SMVCWT,SMC,V_SOIL,VFRAC,CANOPY,CATCH,DQ,
     & EPDT,LYING_SNOW,GC,RESIST,VSHR,CH,PSIS,FRACA,RESFS,
     & F_SE,RESFT,LTIMER)


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
     &  SU10,SV10,ST1P5,SQ1P5,LTIMER
     & )
      ENDIF
C-----------------------------------------------------------------------
CL  4.2 Now that diagnostic calculations are over, update CD and CH
CL      to their correct values (i.e. gridsquare means).
C-----------------------------------------------------------------------
      DO I = 1,P_POINTS
        IF ( ICE_FRACT(I).GT.0.0 .AND. .NOT.LAND_MASK(I) ) THEN
          IF ( ICE_FRACT(I).LT. 0.7 ) THEN
            CD(I) = ( ICE_FRACT(I)*CD_MIZ(I) +
     &                (0.7-ICE_FRACT(I))*CD_LEAD(I) ) / 0.7  ! P2430.5
        CD_STD(I) = CD(I)         ! for SCYCLE: no orog. over sea+ice
            CH(I) = ( ICE_FRACT(I)*CH_MIZ(I) +
     &                (0.7-ICE_FRACT(I))*CH_LEAD(I) ) / 0.7  ! P2430.4
          ELSE
            CD(I) = ( (1.0-ICE_FRACT(I))*CD_MIZ(I) +
     &                (ICE_FRACT(I)-0.7)*CD(I) ) / 0.3       ! P2430.7
        CD_STD(I) = CD(I)         ! for SCYCLE: no orog. over sea+ice
            CH(I) = ( (1.0-ICE_FRACT(I))*CH_MIZ(I) +
     &                (ICE_FRACT(I)-0.7)*CH(I) ) / 0.3       ! P2430.7
          ENDIF
        ENDIF
C-----------------------------------------------------------------------
CL  4.3 Calculate the surface exchange coefficients RHOK(*).
C-----------------------------------------------------------------------
        RHOSTAR(I) = PSTAR(I) / ( R*TSTAR(I) )
C                        ... surface air density from ideal gas equation
C  Calculate resistances for use in Sulphur Cycle
C  (Note that CD_STD, CH and VSHR should never = 0)
         RHO_ARESIST(I) = RHOSTAR(I) * CD_STD(I) * VSHR(I)
             ARESIST(I) = RHOSTAR(I)/RHO_ARESIST(I)
             RESIST_B(I)= (CD_STD(I)/CH(I) - 1.0) * ARESIST(I)
!
        RHOKM_1(I) = RHOSTAR(I) * CD(I) * VSHR(I)             ! P243.124
        RHOKH_1(I) = RHOSTAR(I) * CH(I) * VSHR(I)             ! P243.125
        RHOKE(I) = RESFT(I) * RHOKH_1(I)
C
C     RHOSTAR * CD * VSHR stored for diagnostic output before
C     horizontal interpolation.
C
        RHO_CD_MODV1(I) = RHOKM_1(I)


      ENDDO


      CALL SF_FLUX (
     & P_POINTS,LAND_PTS,LAND_MASK,
     & P1,LAND_INDEX,
     & ALPHA1,DQ,DQ_LEAD,DTEMP,DTEMP_LEAD,DZSOIL,HCONS,ICE_FRACT,
     & LYING_SNOW,QS1,QW_1,RADNET_C,RESFT,RHOKE,RHOKH_1,TI,TL_1,TS1,
     & Z0H,Z0M_EFF,Z1,
     & ASHTF,E_SEA,EPOT,FQW_1,FTL_1,H_SEA,RHOKPM,RHOKPM_POT,
     & TSTAR,VFRAC,TIMESTEP,CANCAP,
     & LTIMER)

C-----------------------------------------------------------------------
CL  4.4.1 Set indicator for unstable suface layer (buoyancy flux +ve.).
CL        if required by logical L_RMBL
C-----------------------------------------------------------------------

      DO I=1,P_POINTS

        IF (L_RMBL.AND.BT_1(I)*FTL_1(I)+BQ_1(I)*FQW_1(I).GT.0.0 )THEN
          NRML(I) = 1
        ELSE
          NRML(I) = 0
        ENDIF
C-----------------------------------------------------------------------
CL  4.5 Multiply surface exchange coefficients that are on the P-grid
CL      by GAMMA(1).Needed for implicit calculations in P244(IMPL_CAL).
CL      RHOKM_1 dealt with in section 4.1 below.
C-----------------------------------------------------------------------
        RHOKH_1(I) = RHOKH_1(I) * GAMMA(1)
        RHOKE(I) = RHOKE(I) * GAMMA(1)
C-----------------------------------------------------------------------
CL  4.5.1 Calculate the standard deviations of layer 1 turbulent
CL        fluctuations of temperature and humidity using approximate
CL        formulae from first order closure.
C-----------------------------------------------------------------------
        VS = SQRT ( RHOKM_1(I)/RHOSTAR(I) * VSHR(I) )
        VSF1_CUBED = 1.25 * ( Z1(I) + Z0M(I) ) * G *
     &              ( BT_1(I)*FTL_1(I) + BQ_1(I)*FQW_1(I) ) / RHOSTAR(I)
C       !---------------------------------------------------------------
C       ! Only calculate standard deviations for unstable surface layers
C       !---------------------------------------------------------------
        IF (VSF1_CUBED .GT. 0.0) THEN
          WS1 = ( VSF1_CUBED + VS * VS * VS ) ** (1.0/3.0)
          T1_SD(I) = MAX ( 0.0 , 1.93 * FTL_1(I) / (RHOSTAR(I) * WS1) )
          Q1_SD(I) = MAX ( 0.0 , 1.93 * FQW_1(I) / (RHOSTAR(I) * WS1) )
        ELSE
          T1_SD(I) = 0.0
          Q1_SD(I) = 0.0
        ENDIF
C-----------------------------------------------------------------------
CL  4.6 For sea points, calculate the wind mixing energy flux and the
CL      sea-surface roughness length on the P-grid, using time-level n
CL      quantities.
C-----------------------------------------------------------------------
        IF (.NOT.LAND_MASK(I)) THEN
          TAU = RHOKM_1(I) * VSHR(I)                          ! P243.130
          IF (ICE_FRACT(I) .GT. 0.0)
     &      TAU = RHOSTAR(I) * CD_LEAD(I) * VSHR(I) * VSHR(I)
          IF (SFME) FME(I) = (1.0-ICE_FRACT(I)) * TAU * SQRT(TAU/RHOSEA)
C                                                             ! P243.96
          Z0MSEA(I) = MAX ( Z0HSEA ,
     &                      (CHARNOCK/G) * (TAU / RHOSTAR(I)) )
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
! RHOKM_1 contains duff data in halos. The P_TO_UV can interpolate this
! into the real data, so first we must update east/west halos
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
