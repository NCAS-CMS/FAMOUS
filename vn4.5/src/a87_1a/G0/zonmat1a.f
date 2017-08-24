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
CLL  Routine: ZONM_ATM -------------------------------------------------
CLL
CLL  Purpose: Calculates zonal mean, quarter global mean and global
CLL           mean values from atmospheric prognostic fields, and
CLL           prints formatted output summary on UNIT 6.
CLL           NB: zonal means are performed on p-grid or u-grid as
CLL           appropriate.  Quarter global means on p-grid include
CLL           overlap row(s).  Formatted output assumes up to 20
CLL           levels initially (general no. of levels to be allowed
CLL           in future).
CLL
CLL  Tested under compiler:   cft77
CLL  Tested under OS version: UNICOS 5.0
CLL
CLL  Author  N.Farnon
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  Date
CLL   3.2    13/07/93 Changed CHARACTER*(*) to CHARACTER*(80) for
CLL                   portability.  Author Tracey Smith.
CLL   3.4    14/09/94 Reduce workspace required by altered definitions
CLL          of IHYDRO, IRAD etc.
CLL   4.0    24/10/95 Simplify do loop structure to remove Cray
CLL          compiler optimisation mis-translation and improve
CLL          vectorisation. R.Rawlins.
!LL  4.3   5/3/97 Correct error in zonal mean print. R A Stratton.
CLL
CLL  Programming standard: UM Doc Paper 3, version 1 (15/1/90)
CLL
CLL  Logical components covered: D61,D64
CLL
CLL  Project task: C61
CLL
CLL  External documentation: UM document C61 - Zonal mean calculations.
CLL
CLLEND------------------------------------------------------------------
C
C*L  Interface and arguments: ------------------------------------------
C
      SUBROUTINE ZONM_ATM ( PMSL,
     1                 PSTAR, T, Q,
     2                 U, V, TSTAR, SOILT, SOILM, SNOWD, CANOPYW,
     3                 SH, EVAP, PPTN, LSRN, CVRN, SNOW, AICE, HICE,
     &                 TAUX, TAUY,
     *                 SUBL, SOEV, SOHF, THRF,SNML, SFRU,SBRU,
     *                 SDTR, SDSR, TDTR, TOLR, TOSW, TISW,
     *                  WFCA,
     *                 CLLIQ,CLICE,
     *                 SIHF,SIMH,SISH,SIST,SISS,SIEF,
     4                 LAND,
     5                 DELTA_AK, DELTA_BK,
     6                 COS_P_LATITUDE, COS_U_LATITUDE,
     7                 ROW_LENGTH, P_ROWS, U_ROWS, P_LEVELS, Q_LEVELS,
     8                 ST_LEVELS, SM_LEVELS, P_FIELD, U_FIELD, BANNER,
     &                 IHYDRO,IRAD,ICLOUD,IFLUXL,ISICE,
     *                IPRTWIND,IPRTKE,IPRTQ,IPRTTEMP,IPRTVAR,IPRTEXTRA,
     9                 ICODE, CMESSAGE)
C
      IMPLICIT NONE
C
C*L------------------COMDECK C_LHEAT------------------------------------
C LC IS LATENT HEAT OF CONDENSATION OF WATER AT 0DEGC
C LF IS LATENT HEAT OF FUSION AT 0DEGC
      REAL LC,LF

      PARAMETER(LC=2.501E6,
     &          LF=0.334E6)
C*----------------------------------------------------------------------
C
C Parameters set by user interface.
      INTEGER
     1       ROW_LENGTH,           ! IN  - Number of points per row
     2       P_ROWS,               ! IN  - Number of rows (p-grid)
     3       P_LEVELS,             ! IN  - Number of pressure levels
     4       Q_LEVELS,             ! IN  - Number of wet levels
!                               levels in multilevel hydrology scheme:
     5       ST_LEVELS,            ! IN  - Number of soil layers
     6       SM_LEVELS             ! IN  - Number of soil moisture
C Parameters derived from those above.
      INTEGER
     1       P_FIELD,              ! IN  - Size of p-grid 2D field
     2       U_ROWS,               ! IN  - Number of rows (u-grid)
     3       U_FIELD               ! IN  - Size of u-grid 2D field
      INTEGER
     &    IHYDRO,                  ! IN  - Hydrology printout
     &    IRAD,                    ! IN  - Radiation printout
     &    ICLOUD,                  ! IN  - Cloud water printout
     &    IFLUXL,                  ! IN  - Fluxes over land
     &    ISICE                    ! IN  - sea-ice
      LOGICAL
     1    LAND(P_FIELD)            ! IN  - Land-sea mask (p-grid)
      REAL
     1    PMSL(P_FIELD),           ! IN  - Pressure at mean sea-level
     2    PSTAR(P_FIELD),          ! IN  - Surface pressure
     3    T(P_FIELD,P_LEVELS),     ! IN  - Atmospheric temperature
     4    Q(P_FIELD,Q_LEVELS),     ! IN  - Atmospheric humidity
     5    U(U_FIELD,P_LEVELS),     ! IN  - Atmospheric u-velocity
     6    V(U_FIELD,P_LEVELS)      ! IN  - Atmospheric v-velocity
      REAL
     1    TSTAR(P_FIELD),          ! IN  - Surface temperature
     2    SOILT(P_FIELD,ST_LEVELS),! IN  - Soil temperature
     3    SOILM(P_FIELD),          ! IN  - Soil moisture
     4    SNOWD(P_FIELD),          ! IN  - Snow depth
     5    CANOPYW(P_FIELD),        ! IN  - Canopy water
     6    SH(P_FIELD),             ! IN  - Surface sensible heat
     7    EVAP(P_FIELD),           ! IN  - Surface evaporation
     8    PPTN(P_FIELD),           ! IN  - Surface total precipitation
     9    LSRN(P_FIELD),           ! IN  - Surface dynamic rain
     A    CVRN(P_FIELD),           ! IN  - Surface convective rain
     B    SNOW(P_FIELD),           ! IN  - Surface snowfall
     &    AICE(P_FIELD),           ! IN  - Seaice fraction (sea points)
     &    HICE(P_FIELD),           ! IN  - Seaice thickness (sea points)
     C    TAUX(U_FIELD),           ! IN  - Surface windstress (x)
     D    TAUY(U_FIELD),           ! IN  - Surface windstress (y)
     &    SUBL(P_FIELD),           ! IN  - Sublimation
     &    SOEV(P_FIELD),           ! IN  - Soil evaporation
     &    THRF(IHYDRO),            ! IN  - Throughfall
     &    SNML(IHYDRO),            ! IN  - Snow melt
     &    SFRU(P_FIELD),           ! IN  - surface runoff
     &    SBRU(P_FIELD),           ! IN  - Sub surface runoff
     &    SOHF(IFLUXL),            ! IN  - Soil heat flux
     &    SDTR(IRAD),           ! IN  - net surface downward rad
     &    SDSR(IRAD),           ! IN  - surface solar downward rad
     &    TDTR(IRAD),           ! IN  - TOA total downward rad
     &    TOLR(IRAD),           ! IN  - TOA outgoing lw rad
     &    TOSW(IRAD),           ! IN  - outgoing solar radiation
     &    TISW(IRAD),           ! IN  - incoming solar radiation
     &    WFCA(IFLUXL),           ! IN  - water flux into canopy
     &    CLLIQ(ICLOUD,Q_LEVELS),  ! IN  - cloud liquid water
     &    CLICE(ICLOUD,Q_LEVELS),  ! IN  - cloud Ice water
     &    SIHF(ISICE),           ! IN  - sea-ice heat flux
     &    SIMH(ISICE),           ! IN  - sea-ice melting heat
     &    SISH(ISICE),           ! IN  - sea-ice sensible heat
     &    SIST(ISICE),           ! IN  - sea-ice surface total rad
     &    SISS(ISICE),           ! IN  - sea-ice surface solar
     &    SIEF(ISICE)            ! IN  - sea-ice energy flux
      REAL
     1    DELTA_AK(P_LEVELS),      ! IN  - Hybrid layer thickness A
     2    DELTA_BK(P_LEVELS),      ! IN  - Hybrid layer thickness B
     3    COS_P_LATITUDE(P_FIELD), ! IN  - Cosine latitude (p-grid)
     4    COS_U_LATITUDE(U_FIELD)  ! IN  - Cosine latitude (u-grid)
      INTEGER
     &    IPRTWIND,                ! IN  - wind printout
     &    IPRTKE,                  ! IN  - Kinetic energy printout
     &    IPRTQ,                   ! IN  - Q printout
     &    IPRTTEMP,                ! IN  - temp printout
     &    IPRTVAR,                 ! IN  - variance printout
     &    IPRTEXTRA                ! IN  - extra printouts
      INTEGER
     1    ICODE                    ! OUT - Error exit code
      CHARACTER*80
     1    BANNER                   ! IN  - Description of data fields
      CHARACTER*(80)
     1    CMESSAGE                 ! OUT - Error message
C
C*----------------------------------------------------------------------
C
C  Subroutines called
C
      EXTERNAL ZONM,COLM,GLBM,P_TO_UV
C

C  Local variables
C
      REAL
     &    AMASK(P_FIELD),          ! LOC - Mask of all ones
     &    LMASK(P_FIELD),          ! LOC - Land mask (p-grid)
     &    SMASK(P_FIELD)           ! LOC - Sea mask (p-grid)
      LOGICAL
     1    LAPTS(P_ROWS,2),         ! LOC - Marks all rows
     2    LLPTS(P_ROWS,2),         ! LOC - Marks rows with land pts
     3    LSPTS(P_ROWS,2),         ! LOC - Marks rows with sea pts
     4    LQGAPTS(4,2),            ! LOC - Marks all pts for each 1/4
     5    LQGLPTS(4,2),            ! LOC - Marks land pts for each 1/4
     6    LQGSPTS(4,2),            ! LOC - Marks sea pts for each 1/4
     7    LGAPTS(2),               ! LOC - Marks all pts of globe
     8    LGLPTS(2),               ! LOC - Marks land pts of globe
     9    LGSPTS(2)                ! LOC - Marks sea pts of globe
      INTEGER
     1    I,                       ! LOC - Field index
     2    ROW,                     ! LOC - Row number
     3    J,                       ! LOC - Row length
     4    LEVEL,                   ! LOC - Vertical level
     8    LPTS(P_ROWS,2),          ! LOC - No of land points/row
     9    SPTS(P_ROWS,2),          ! LOC - No of sea points/row
     A    QGLPTS(4,2),             ! LOC - No of land pts/ 1/4 globe
     B    QGSPTS(4,2),             ! LOC - No of sea pts/ 1/4 globe
     C    GLPTS(2),                ! LOC - No of land pts/globe
     D    GSPTS(2),                ! LOC - No of sea pts/globe
     E    QUART,                   ! LOC - Index for each 1/4 globe
     F    IGRID                    ! LOC - Index for p or u-grid
      REAL
     1    P_MASS(P_FIELD,P_LEVELS),! IN  - Mass weighting (p-grid)
     2    U_MASS(U_FIELD,P_LEVELS),! LOC - Mass weighting (u-grid)
     3    S_PMASS(P_FIELD),        ! LOC - Dummy wgt for surf(p-grid)
     4    S_UMASS(U_FIELD),        ! LOC - Dummy wgt for surf(u-grid)
     5    PSTAR_U(U_FIELD)         ! LOC - wgt for atmos (u-grid)

C Additional local variables used in the calculation of Energy
C diagnostics

      REAL
     &    WORK1(P_FIELD,P_LEVELS),       ! 3D work field
     &    WORK2(P_FIELD)                 ! 2D work field
     &    ,WORK3(U_FIELD,P_LEVELS)       ! 3D work field

C Values used in the calculation of Zonal means

      REAL
     &    Z_TKE(U_ROWS),               ! Total Kinetic Energy
     &    Z_ZKE(U_ROWS),               ! Zonal Kinetic Energy
     &    Z_EKE(U_ROWS),               ! Eddy Kinetic Energy
     &    Z_VAR_T(P_ROWS),             ! Variance of T
     &    Z_VAR_Q(P_ROWS)              ! Variance of Q
      REAL
     &  Z_VAR_T_LEV(P_ROWS,P_LEVELS), ! Variance of T at each level
     &  Z_VAR_Q_LEV(P_ROWS,Q_LEVELS), ! Variance of Q at each level
     &  Z_TKE_LEV(U_ROWS,P_LEVELS),   ! Total Kinetic Energy per level
     &  Z_ZKE_LEV(U_ROWS,P_LEVELS),   ! Zonal Kinetic Energy per level
     &  Z_EKE_LEV(U_ROWS,P_LEVELS)    ! Eddy Kinetic Energy per level

C Values used in the calculation of Global means

      REAL
     &    G_TKE,                       ! Total Kinetic Energy
     &    G_ZKE,                       ! Zonal Kinetic Energy
     &    G_EKE,                       ! Eddy Kinetic Energy
     &    G_VAR_T,                     ! Variance of T
     &    G_VAR_Q                      ! Variance of Q
      REAL
     &  G_VAR_T_LEV(P_LEVELS),        ! Variance of T at each level
     &  G_VAR_Q_LEV(Q_LEVELS),        ! Variance of Q at each level
     &  G_TKE_LEV(P_LEVELS),          ! Total Kinetic Energy per level
     &  G_ZKE_LEV(P_LEVELS),          ! Zonal Kinetic Energy per level
     &  G_EKE_LEV(P_LEVELS)           ! Eddy Kinetic Energy per level
      REAL
     1    TATMOS(P_FIELD),         ! LOC - Col mean atmos temperature
     2    QATMOS(P_FIELD),         ! LOC - Col mean atmos humidity
     3    UATMOS(U_FIELD),         ! LOC - Col mean atmos temperature
     4    VATMOS(U_FIELD),         ! LOC - Col mean atmos humidity
     5    CLLIQ_F(ICLOUD),        ! LOC - Cloud liquid water
     6    CLICE_F(ICLOUD)         ! LOC - Cloud ice water
      INTEGER                      ! NOTE -Start & end pts for each 1/4
     1    START_ROW(2,4),          ! LOC  - Start of 1/4 global section
     &    END_ROW(2,4),            ! LOC  - End of 1/4 global section
     &    LAST_ROW                 ! LOC  - p_rows or u_rows
      REAL
     1    Z_TATMOS(P_ROWS),        ! LOC - Zonm of col mn atmos temperat
     2    Z_QATMOS(P_ROWS),        ! LOC - Zonm of col mean atmos humidi
     3    Z_UATMOS(U_ROWS),        ! LOC - Zonal mn of col mean atmos u-
     4    Z_VATMOS(U_ROWS)         ! LOC - Zonal mn of col mean atmos v-
      REAL
     1    Z_T(P_ROWS,P_LEVELS),    ! LOC - Atmospheric temperature
     2    Z_Q(P_ROWS,Q_LEVELS),    ! LOC - Atmospheric humidity
     3    Z_U(U_ROWS,P_LEVELS),    ! LOC - Atmospheric u-velocity
     4    Z_V(U_ROWS,P_LEVELS)     ! LOC - Atmospheric v-velocity
      REAL
     1    Z_PMSL(P_ROWS),          ! LOC - Pressure at mean sea-level
     2    Z_PSTAR(P_ROWS),         ! LOC - Surface pressure
     3    Z_TSTAR(P_ROWS),         ! LOC - Surface temperature
     4    Z_SOILT(P_ROWS,ST_LEVELS), ! LOC - Soil temperature
     5    Z_SOILM(P_ROWS),         ! LOC - Soil moisture
     6    Z_SNOWD(P_ROWS),         ! LOC - Snow depth
     7    Z_CANOPYW(P_ROWS),       ! LOC - Canopy water
     8    Z_SH(P_ROWS),            ! LOC - Surface sensible heat
     9    Z_EVAP(P_ROWS),          ! LOC - Surface evaporation
     A    Z_PPTN(P_ROWS),          ! LOC - Surface total precipitation
     B    Z_LSRN(P_ROWS),          ! LOC - Surface dynamic rain
     C    Z_CVRN(P_ROWS),          ! LOC - Surface convective rain
     D    Z_SNOW(P_ROWS),          ! LOC - Surface snowfall
     E    Z_TAUX(U_ROWS),          ! LOC - Surface windstress (x)
     F    Z_TAUY(U_ROWS),          ! LOC - Surface windstress (y)
     G    Z_CLLIQ(P_ROWS),         ! LOC - Cloud liquid water
     H    Z_CLICE(P_ROWS),         ! LOC - Cloud ice water
     I    Z_TDTR(P_ROWS),          ! LOC - TOA net total down rad
     J    Z_TOLR(P_ROWS),          ! LOC - TOA outgoing LW rad
     K    Z_TOSW(P_ROWS),          ! LOC - outgoing solar rad
     L    Z_TISW(P_ROWS),          ! LOC - incoming solar rad
     M    Z_ALBEDO(P_ROWS)         ! LOC - planetary albedo
C      land only
      REAL
     1    Z_L_TSTAR(P_ROWS),       ! LOC - Surface temperature
     2    Z_L_SH(P_ROWS),          ! LOC - Surface sensible heat
     3    Z_L_EVAP(P_ROWS),        ! LOC - Surface evaporation
     4    Z_L_PPTN(P_ROWS),        ! LOC - Surface total precipitation
     5    Z_L_SNOW(P_ROWS),        ! LOC - Surface snowfall
     6    Z_L_TAUX(U_ROWS),        ! LOC - Surface windstress (x)
     7    Z_L_TAUY(U_ROWS),        ! LOC - Surface windstress (y)
     8    Z_L_SUBL(P_ROWS),        ! LOC - sublimation
     9    Z_L_SFRU(P_ROWS),        ! LOC - surface runoff
     A    Z_L_SBRU(P_ROWS),        ! LOC - Subsurface runoff
     B    Z_L_SDTR(P_ROWS),        ! LOC - surface downward total rad
     C    Z_L_SDSR(P_ROWS)         ! LOC - surface solar downward rad
C       Sea only
      REAL
     1    Z_S_TSTAR(P_ROWS),       ! LOC - Surface temperature
     2    Z_S_SH(P_ROWS),          ! LOC - Surface sensible heat
     3    Z_S_EVAP(P_ROWS),        ! LOC - Surface evaporation
     4    Z_S_PPTN(P_ROWS),        ! LOC - Surface total precipitation
     5    Z_S_SNOW(P_ROWS),        ! LOC - Surface snowfall
     &    Z_S_AICE(P_ROWS),        ! LOC - Seaice fraction
     &    Z_S_HICE(P_ROWS),        ! LOC - Seaice thickness
     6    Z_S_TAUX(U_ROWS),        ! LOC - Surface windstress (x)
     7    Z_S_TAUY(U_ROWS),        ! LOC - Surface windstress (y)
     8    Z_S_SUBL(P_ROWS),        ! LOC - Sublimation
     9    Z_S_SDTR(P_ROWS),        ! LOC - surface downward total rad
     A    Z_S_SDSR(P_ROWS)         ! LOC - surface solar downward rad
      REAL
     1    G_Q_PMSL(4),             ! LOC - 1/4 Global mean PMSL
     2    G_Q_PSTAR(4),            ! LOC - 1/4 Global mean surface press
     3    G_Q_TSTAR(4),            ! LOC - Surface temperature
     4    G_Q_L_TSTAR(4),          ! LOC - Surface temperature
     5    G_Q_S_TSTAR(4),          ! LOC - Surface temperature
     6    G_Q_SOILT(4,ST_LEVELS),  ! LOC - Soil temperature
     7    G_Q_SOILM(4),            ! LOC - Soil moisture
     8    G_Q_SNOWD(4),            ! LOC - Snow depth
     9    G_Q_CANOPYW(4),          ! LOC - Canopy water
     A    G_Q_SH(4),               ! LOC - Surface sensible heat
     B    G_Q_EVAP(4),             ! LOC - Surface evaporation
     C    G_Q_PPTN(4),             ! LOC - Surface total precipitation
     D    G_Q_LSRN(4),             ! LOC - Surface dynamic rain
     E    G_Q_CVRN(4),             ! LOC - Surface convective rain
     F    G_Q_SNOW(4),             ! LOC - Surface snowfall
     &    G_Q_SUBL(4),           ! LOC -Sublimation
     &    G_Q_SDTR(4),           ! LOC -Net surface downward total rad
     &    G_Q_SDSR(4),           ! LOC -Net surface downward total rad
     &    G_Q_TDTR(4),           ! LOC -TOA  downward total rad
     &    G_Q_TOLR(4),           ! LOC -TOA outgoing lw rad
     &    G_Q_TOSW(4),           ! LOC -outgoing radiation
     &    G_Q_TISW(4),           ! LOC -incoming radiation
     &    G_Q_ALBEDO(4),         ! LOC -Planetary albedo
     &    G_Q_CLLIQL(4,Q_LEVELS),  ! LOC - Cloud liquid water
     &    G_Q_CLICEL(4,Q_LEVELS),  ! LOC - Cloud Ice water
     &    G_Q_CLLIQ(4),           ! LOC - Cloud liquid water
     &    G_Q_CLICE(4),           ! LOC - Cloud Ice water
     &    G_Q_SIHF(4),           ! LOC - Sea-ice heat flux
     &    G_Q_SIMH(4),           ! LOC - Sea-ice melting heat
     &    G_Q_SISH(4),           ! LOC - Sea-ice sensible heat
     &    G_Q_SIST(4),           ! LOC - Sea-ice surface total rad
     &    G_Q_SISS(4),           ! LOC - Sea-ice surface solar rad
     &    G_Q_SIEF(4)            ! LOC - Sea-ice energy flux
C     Land only quarter globe means
      REAL
     &    G_Q_L_SH(4),             ! LOC - Surface sensible heat
     &    G_Q_L_EVAP(4),           ! LOC - Surface evaporation
     &    G_Q_L_PPTN(4),           ! LOC - Surface total precipitation
     &    G_Q_L_SNOW(4),           ! LOC - Surface snowfall
     &    G_Q_L_SUBL(4),           ! LOC - Sublimation
     &    G_Q_L_SOEV(4),           ! LOC - Soil evaporation
     &    G_Q_L_THRF(4),           ! LOC - throughfall
     &    G_Q_L_SNML(4),           ! LOC - Snowmelt
     &    G_Q_L_SFRU(4),           ! LOC - Surface runoff
     &    G_Q_L_SBRU(4),           ! LOC - Subsurface runoff
     &    G_Q_L_SOHF(4),           ! LOC -Soil heat flux (top 2 layers)
     &    G_Q_L_SDTR(4),           ! LOC -Net surface downward total rad
     &    G_Q_L_SDSR(4)            ! LOC -Net surface downward total rad
C     Land only quarter globe means derived from above
      REAL
     &    G_Q_L_WFSS(4),           ! LOC - Net water flux surface snow
     &    G_Q_L_WFCA(4),           ! LOC - Net water flux canopy
     &    G_Q_L_WFSO(4),           ! LOC - Net water flux into soil
     &    G_Q_WAFL(4),             ! LOC - Net water flux atmosphere
     &    G_Q_L_ENFS(4),           ! LOC - Net energy flux into soil
     &    G_Q_ENFL(4)              ! LOC - Net energy flux atmosphere
C     Sea only quarter globe means
      REAL
     &    G_Q_S_SH(4),             ! LOC - Surface sensible heat
     &    G_Q_S_EVAP(4),           ! LOC - Surface evaporation
     &    G_Q_S_PPTN(4),           ! LOC - Surface total precipitation
     &    G_Q_S_SNOW(4),           ! LOC - Surface snowfall
     &    G_Q_S_SUBL(4),           ! LOC - Sublimation
     &    G_Q_S_SDTR(4),           ! LOC -Net surface downward total rad
     &    G_Q_S_SDSR(4)            ! LOC -Net surface downward total rad
C     Sea only quarter globe means derived from above
C     REAL
C    &    G_Q_S_SH(4)              ! LOC - Surface sensible heat
      REAL
     1    G_Q_T(4,P_LEVELS),       ! LOC - 1/4Global mean atmos temperat
     2    G_Q_Q(4,Q_LEVELS),       ! LOC - 1/4 Global mean atmos humidit
     3    G_Q_U(4,P_LEVELS),       ! LOC - Atmospheric u-velocity
     4    G_Q_V(4,P_LEVELS),       ! LOC - Atmospheric v-velocity
     5    G_Q_TAUX(4),             ! LOC - Surface windstress (x)
     6    G_Q_TAUY(4),             ! LOC - Surface windstress (y)
     7    G_Q_L_TAUX(4),           ! LOC - Surface windstress (x)
     8    G_Q_L_TAUY(4),           ! LOC - Surface windstress (y)
     9    G_Q_S_TAUX(4),           ! LOC - Surface windstress (x)
     A    G_Q_S_TAUY(4),           ! LOC - Surface windstress (y)
     B    G_Q_TKE(4),              ! LOC - Total Kinetic Energy
     C    G_Q_ZKE(4),              ! LOC - Zonal Kinetic Energy
     D    G_Q_EKE(4)               ! LOC - Eddy Kinetic Energy
      REAL
     1    G_Q_TATMOS(4),           ! LOC - 1/4 global mn atmos temperatu
     2    G_Q_QATMOS(4),           ! LOC - 1/4 global mean atmos humidit
     3    G_Q_UATMOS(4),           ! LOC - 1/4 global mean atmos u-winds
     4    G_Q_VATMOS(4)            ! LOC - 1/4 global mean atmos v-winds
      REAL
     1    G_PMSL,                  ! LOC - 1/4 Global mean PMSL
     2    G_PSTAR,                 ! LOC - 1/4 Global mean surface press
     3    G_TSTAR,                 ! LOC - Surface temperature
     4    G_L_TSTAR,               ! LOC - Surface temperature
     5    G_S_TSTAR,               ! LOC - Surface temperature
     6    G_SOILT(ST_LEVELS),      ! LOC - Soil temperature
     7    G_SOILM,                 ! LOC - Soil moisture
     8    G_SNOWD,                 ! LOC - Snow depth
     9    G_CANOPYW,               ! LOC - Canopy water
     A    G_SH,                    ! LOC - Surface sensible heat
     B    G_EVAP,                  ! LOC - Surface evaporation
     C    G_PPTN,                  ! LOC - Surface total precipitation
     D    G_LSRN,                  ! LOC - Surface dynamic rain
     E    G_CVRN,                  ! LOC - Surface convective rain
     F    G_SNOW                   ! LOC - Surface snowfall
      REAL
     1    G_T(P_LEVELS),           ! LOC - 1/4Global mean atmos temperat
     2    G_Q(Q_LEVELS),           ! LOC - 1/4 Global mean atmos humidit
     3    G_U(P_LEVELS),           ! LOC - Atmospheric u-velocity
     4    G_V(P_LEVELS),           ! LOC - Atmospheric v-velocity
     5    G_TAUX,                  ! LOC - Surface windstress (x)
     6    G_TAUY,                  ! LOC - Surface windstress (y)
     7    G_L_TAUX,                ! LOC - Surface windstress (x)
     8    G_L_TAUY,                ! LOC - Surface windstress (y)
     9    G_S_TAUX,                ! LOC - Surface windstress (x)
     A    G_S_TAUY,                ! LOC - Surface windstress (y)
     &    G_SUBL,                  ! LOC - Sublimation
     &    G_SDTR,                  ! LOC - Surface net total down rad
     &    G_SDSR,                  ! LOC - Surface net solar down rad
     &    G_TDTR,                  ! LOC - TOA net total down rad
     &    G_TOLR,                  ! LOC - TOA out going LW rad
     &    G_TOSW,                  ! LOC - outgoing solar radiation
     &    G_TISW,                  ! LOC - incoming solar radiation
     &    G_ALBEDO,                ! LOC - Planetary albedo
     &    G_CLLIQL(Q_LEVELS),      ! LOC - Cloud liquid water
     &    G_CLICEL(Q_LEVELS),      ! LOC - Cloud Ice water
     &    G_CLLIQ,                 ! LOC - Cloud liquid water
     &    G_CLICE,                 ! LOC - Cloud Ice water
     &    G_SIHF,           ! LOC - Sea-ice heat flux
     &    G_SIMH,           ! LOC - Sea-ice melting heat
     &    G_SISH,           ! LOC - Sea-ice sensible heat
     &    G_SIST,           ! LOC - Sea-ice surface total rad
     &    G_SISS,           ! LOC - Sea-ice surface solar rad
     &    G_SIEF            ! LOC - Sea-ice energy flux
      REAL
     1    G_TATMOS,                ! LOC - Global mn atmos temperature
     2    G_QATMOS,                ! LOC - Global mean atmos humidity
     3    G_UATMOS,                ! LOC - Global mean atmos u-winds
     4    G_VATMOS                 ! LOC - Global mean atmos v-winds
C     Land only global means
      REAL
     &    G_L_SH,                ! LOC - Surface sensible heat
     &    G_L_EVAP,              ! LOC - Surface evaporation
     &    G_L_PPTN,              ! LOC - Surface total precipitation
     &    G_L_SNOW,              ! LOC - Surface snowfall
     &    G_L_SUBL,              ! LOC - Sublimation
     &    G_L_SOEV,              ! LOC - Soil evaporation
     &    G_L_THRF,              ! LOC - throughfall
     &    G_L_SNML,              ! LOC - Snowmelt
     &    G_L_SFRU,              ! LOC - Surface runoff
     &    G_L_SBRU,              ! LOC - Subsurface runoff
     &    G_L_SOHF,              ! LOC -Soil heat flux (top 2 layers)
     &    G_L_SDTR,              ! LOC -Net surface downward total rad
     &    G_L_SDSR               ! LOC -Net surface downward total rad
C     Land only global means derived from above
      REAL
     &    G_L_WFSS,           ! LOC - Net water flux surface snow
     &    G_L_WFCA,           ! LOC - Net water flux canopy
     &    G_L_WFSO,           ! LOC - Net water flux into soil
     &    G_WAFL,             ! LOC - Net water flux atmosphere
     &    G_ENFL,             ! LOC - Net energy flux atmosphere
     &    G_L_ENFS            ! LOC - Net energy flux into soil
C     Sea only  global means
      REAL
     &    G_S_SH,             ! LOC - Surface sensible heat
     &    G_S_EVAP,           ! LOC - Surface evaporation
     &    G_S_PPTN,           ! LOC - Surface total precipitation
     &    G_S_SNOW,           ! LOC - Surface snowfall
     &    G_S_SUBL,           ! LOC - Sublimation
     &    G_S_SDTR,           ! LOC -Net surface downward total rad
     &    G_S_SDSR            ! LOC -Net surface downward total rad
C
      CHARACTER*130
     &    CTITLE          ! title for print
     &    ,CHEAD1         ! header for table
     &    ,CHEAD2         ! second header for table
CL----------------------------------------------------------------------

CL 1. Initializing the start and end pt for each 1/4 of globe
CL
CL 1.1 For p-grid
CL
      START_ROW(1,1)=1                    ! 1st quarter start
        END_ROW(1,1)=(P_ROWS-1)/3 + 1     ! 1st quarter end
      START_ROW(1,2)=(P_ROWS-1)/3 + 1     ! 2nd quarter start
        END_ROW(1,2)=(P_ROWS-1)/2 + 1     ! 2nd quarter end
      START_ROW(1,3)=(P_ROWS-1)/2 + 1     ! 3rd quarter start
        END_ROW(1,3)=2*(P_ROWS-1)/3 + 1   ! 3rd quarter end
      START_ROW(1,4)=2*(P_ROWS-1)/3 + 1   ! 4th quarter start
        END_ROW(1,4)=P_ROWS               ! 4th quarter end
CL
CL 1.2 For u-grid
CL
      START_ROW(2,1)=1                    ! 1st quarter start
        END_ROW(2,1)=(U_ROWS)/3           ! 1st quarter end
      START_ROW(2,2)=(U_ROWS)/3 + 1       ! 2nd quarter start
        END_ROW(2,2)=(U_ROWS)/2           ! 2nd quarter end
      START_ROW(2,3)=(U_ROWS)/2 + 1       ! 3rd quarter start
        END_ROW(2,3)=2*(U_ROWS)/3         ! 3rd quarter end
      START_ROW(2,4)=2*(U_ROWS)/3 + 1     ! 4th quarter start
        END_ROW(2,4)=U_ROWS               ! 4th quarter end
CL----------------------------------------------------------------------
CL 2. Calculate mass weights and set up masks for weighted sums
CL
      DO I=1,P_FIELD
        AMASK(I) = 1.0               ! Set mask of all land & sea pts =1
        IF (LAND(I)) THEN
          LMASK(I) = 1.0             ! Set mask of all land pts =1
        ELSE
          LMASK(I) = 0.0             ! Set mask of all other pts =0
        ENDIF
          SMASK(I) = 1.0-LMASK(I)    ! Set mask of all sea pts =1 or 0
      END DO
CL----------------------------------------------------------------------
CL 3. Calculate no of land/sea points on row-by-row basis
CL    and set logical arrays to denote active land/sea rows

      DO IGRID=1,2                  ! p-grid, then u-grid

        IF(IGRID.EQ.1) THEN
           LAST_ROW=P_ROWS
        ELSE
           LAST_ROW=U_ROWS
        ENDIF

        DO ROW=1,LAST_ROW           ! iterate over all rows in grid

          LPTS(ROW,IGRID) = 0       ! Initialize no. of land pts =0
          DO J=1,ROW_LENGTH          ! Loop through each pt/row
            LPTS(ROW,IGRID) = LPTS(ROW,IGRID)+
     &              LMASK((ROW-1)*ROW_LENGTH+J) ! Sum land pts for grid
          ENDDO                      ! J 1,ROW_LENGTH
          SPTS(ROW,IGRID) = ROW_LENGTH-LPTS(ROW,IGRID)  ! No. of sea pts
          LAPTS(ROW,IGRID) = .TRUE.  ! land or sea pts/row?= always true
          LLPTS(ROW,IGRID) = LPTS(ROW,IGRID).GT.0   ! Any land pts/row=t
          LSPTS(ROW,IGRID) = SPTS(ROW,IGRID).GT.0   ! Any sea  pts/row=t

        ENDDO                       ! ROW   1, last row
      ENDDO                         ! IGRID 1,2

      DO 21 IGRID=1,2
      GLPTS(IGRID) = 0              ! Initialize no. of global pts =0
      DO 20 I=1,4                   ! Loop through each 1/4 of globe
        QGLPTS(I,IGRID) = 0         ! Init. no. of 1/4 global pts =0
        DO 30 ROW=START_ROW(IGRID,I),END_ROW(IGRID,I)
     &                              ! Loop from start to end pt for each
          QGLPTS(I,IGRID) = QGLPTS(I,IGRID)+LPTS(ROW,IGRID)! Land pts/qu
 30     CONTINUE
        QGSPTS(I,IGRID) = ROW_LENGTH * (END_ROW(IGRID,I)
     &          - START_ROW(IGRID,I)) - QGLPTS(I,IGRID)    ! No sea pts
        LQGAPTS(I,IGRID) = .TRUE.                 ! All land&seapts/quar
        LQGLPTS(I,IGRID) = QGLPTS(I,IGRID).GT.0   ! All land pts/quart=t
        LQGSPTS(I,IGRID) = QGSPTS(I,IGRID).GT.0   ! All sea pts/quart=tr
        GLPTS(IGRID) = GLPTS(IGRID) + QGLPTS(I,IGRID) ! No. of global la
 20   CONTINUE

      GSPTS(IGRID) = ROW_LENGTH * P_ROWS - GLPTS(IGRID) ! No. of global
      LGAPTS(IGRID) = .TRUE.                  ! All land&sea pts=true fo
      LGLPTS(IGRID) = GLPTS(IGRID) .GT. 0     ! Logical land pts=true fo
      LGSPTS(IGRID) = GSPTS(IGRID) .GT. 0     ! Logical sea pts=true for
 21   CONTINUE

CL----------------------------------------------------------------------
CL 4. Compute mass weights from P* and hybrid coordinates
CL    Allow for negative sign of delta_ak and delta_bk
CL
      DO LEVEL=1,P_LEVELS
        DO I=1,P_FIELD
          P_MASS(I,LEVEL) = -DELTA_AK(LEVEL)-DELTA_BK(LEVEL)*PSTAR(I)
        END DO
CL
CL 4.1 Interpolate onto U-Grid
CL
        CALL P_TO_UV(P_MASS(1,LEVEL),U_MASS(1,LEVEL),P_FIELD,U_FIELD,
     &               ROW_LENGTH,P_ROWS)
      END DO

CL
CL 4.2 Interpolate PSTAR onto U-Grid
CL
      CALL P_TO_UV(PSTAR,PSTAR_U,P_FIELD,U_FIELD,ROW_LENGTH,P_ROWS)
CL
CL 4.3 Set dummy weighting for surface variables to one
CL     Following two DO loops labelled due to fpp translation problem
CL
      DO 430 I=1,P_FIELD
        S_PMASS(I)=1.0                ! Dummy weight for surface vars(p-
 430  CONTINUE

      DO 431 I=1,U_FIELD
        S_UMASS(I)=1.0                ! Dummy weight for surface vars(u-
 431  CONTINUE

CL----------------------------------------------------------------------
CL 5. Compute column means of primary atmospheric variables
CL
      CALL COLM(T,TATMOS,P_MASS,P_ROWS,ROW_LENGTH,P_LEVELS)
      CALL COLM(Q,QATMOS,P_MASS,P_ROWS,ROW_LENGTH,Q_LEVELS)
      CALL COLM(U,UATMOS,U_MASS,U_ROWS,ROW_LENGTH,P_LEVELS)
      CALL COLM(V,VATMOS,U_MASS,U_ROWS,ROW_LENGTH,P_LEVELS)
CL----------------------------------------------------------------------
CL 6. Compute zonal means
CL
CL 6.1 Mass weighted on p-grid for surface variables (all points)
CL
        CALL ZONM(PMSL,Z_PMSL,AMASK,S_PMASS,LAPTS(1,1),ROW_LENGTH,
     &                                                        P_ROWS)
        CALL ZONM(PSTAR,Z_PSTAR,AMASK,S_PMASS,LAPTS(1,1),ROW_LENGTH,
     &                                                        P_ROWS)
        CALL ZONM(TSTAR,Z_TSTAR,AMASK,S_PMASS,LAPTS(1,1),ROW_LENGTH,
     &                                                        P_ROWS)
        CALL ZONM(SNOWD,Z_SNOWD,AMASK,S_PMASS,LAPTS(1,1),ROW_LENGTH,
     &                                                        P_ROWS)
        CALL ZONM(SH,Z_SH,AMASK,S_PMASS,LAPTS(1,1),ROW_LENGTH,P_ROWS)

        CALL ZONM(EVAP,Z_EVAP,AMASK,S_PMASS,LAPTS(1,1),ROW_LENGTH,
     &                                                        P_ROWS)
        CALL ZONM(PPTN,Z_PPTN,AMASK,S_PMASS,LAPTS(1,1),ROW_LENGTH,
     &                                                        P_ROWS)
        CALL ZONM(LSRN,Z_LSRN,AMASK,S_PMASS,LAPTS(1,1),ROW_LENGTH,
     &                                                        P_ROWS)
        CALL ZONM(CVRN,Z_CVRN,AMASK,S_PMASS,LAPTS(1,1),ROW_LENGTH,
     &                                                        P_ROWS)
        CALL ZONM(SNOW,Z_SNOW,AMASK,S_PMASS,LAPTS(1,1),ROW_LENGTH,
     &                                                        P_ROWS)
CL
CL 6.2 Mass weighted on p-grid for column meaned atmos variables (all po
CL
        CALL ZONM(TATMOS,Z_TATMOS,AMASK,PSTAR,LAPTS(1,1),ROW_LENGTH,
     &                                                        P_ROWS)
        CALL ZONM(QATMOS,Z_QATMOS,AMASK,PSTAR,LAPTS(1,1),ROW_LENGTH,
     &                                                        P_ROWS)
CL
CL 6.3 Mass weighted on u-grid for column meaned atmos variables (all po
CL
        CALL ZONM(UATMOS,Z_UATMOS,AMASK,PSTAR_U,LAPTS(1,1),ROW_LENGTH,
     &                                                        U_ROWS)
        CALL ZONM(VATMOS,Z_VATMOS,AMASK,PSTAR_U,LAPTS(1,1),ROW_LENGTH,
     &                                                        U_ROWS)
CL
CL 6.4 Mass weighted on u-grid for surface variables (all points)
CL
        CALL ZONM(TAUX,Z_TAUX,AMASK,S_UMASS,LAPTS(1,2),ROW_LENGTH,
     &                                                        U_ROWS)
        CALL ZONM(TAUY,Z_TAUY,AMASK,S_UMASS,LAPTS(1,2),ROW_LENGTH,
     &                                                        U_ROWS)
CL
CL 6.5 Mass weighted on p-levels for atmospheric variables (all points)
CL
      DO LEVEL=1,P_LEVELS
        CALL ZONM(U(1,LEVEL),Z_U(1,LEVEL),AMASK,U_MASS(1,LEVEL),
     &            LAPTS(1,2),ROW_LENGTH,U_ROWS)
        CALL ZONM(V(1,LEVEL),Z_V(1,LEVEL),AMASK,U_MASS(1,LEVEL),
     &            LAPTS(1,2),ROW_LENGTH,U_ROWS)
        CALL ZONM(T(1,LEVEL),Z_T(1,LEVEL),AMASK,P_MASS(1,LEVEL),
     &            LAPTS(1,1),ROW_LENGTH,P_ROWS)


CL
CL 6.6 Mass weighted on p_levels for kinetic energy
CL
        DO I=1,U_FIELD
          WORK3(I,LEVEL)=
     &        U(I,LEVEL)*U(I,LEVEL)+V(I,LEVEL)*V(I,LEVEL)
        END DO

        CALL ZONM(WORK3(1,LEVEL),Z_TKE_LEV(1,LEVEL),AMASK,
     &                 U_MASS(1,LEVEL),LAPTS(1,2),ROW_LENGTH,U_ROWS)

        DO I=1,U_ROWS
          Z_ZKE_LEV(I,LEVEL)=Z_U(I,LEVEL)*Z_U(I,LEVEL)+
     &                       Z_V(I,LEVEL)*Z_V(I,LEVEL)
          Z_EKE_LEV(I,LEVEL)=Z_TKE_LEV(I,LEVEL)-Z_ZKE_LEV(I,LEVEL)
        END DO

      END DO

C      Total KE in vertical mean sense, mass weighted

      CALL COLM(WORK3,WORK2,U_MASS,U_ROWS,ROW_LENGTH,P_LEVELS)
      CALL ZONM(WORK2,Z_TKE,
     &          AMASK,PSTAR_U,LAPTS(1,2),ROW_LENGTH,U_ROWS)

      DO LEVEL=1,P_LEVELS
       DO I=1,U_FIELD
        WORK3(I,LEVEL) =   Z_ZKE_LEV((I-1)/ROW_LENGTH+1,LEVEL)
       ENDDO
      ENDDO

      CALL COLM(WORK3,WORK2,U_MASS,U_ROWS,ROW_LENGTH,P_LEVELS)
      CALL ZONM(WORK2,Z_ZKE,AMASK,PSTAR_U,LAPTS(1,2),ROW_LENGTH,
     &          U_ROWS)
      DO I=1,U_ROWS
        Z_EKE(I)=Z_TKE(I)-Z_ZKE(I)
      END DO
CL
CL 6.7 Mass weighted on p-levels for temperature variance
CL
      DO LEVEL=1,P_LEVELS

        DO I=1,P_FIELD
          WORK1(I,LEVEL)=T(I,LEVEL)*T(I,LEVEL)
        END DO

        CALL ZONM(WORK1(1,LEVEL),WORK2,AMASK,P_MASS(1,LEVEL),
     &                 LAPTS(1,1),ROW_LENGTH,P_ROWS)

        DO I=1,P_ROWS
          Z_VAR_T_LEV(I,LEVEL)=WORK2(I)-Z_T(I,LEVEL)*Z_T(I,LEVEL)
        END DO

      END DO

C Vertical mean temperature variance, mass weighted

      CALL COLM(WORK1,WORK2,P_MASS,P_ROWS,ROW_LENGTH,P_LEVELS)
      CALL ZONM(WORK2,Z_VAR_T,AMASK,PSTAR,LAPTS(1,1),ROW_LENGTH,P_ROWS)

      DO I=1,P_ROWS
        Z_VAR_T(I)=Z_VAR_T(I)-Z_TATMOS(I)*Z_TATMOS(I)
      END DO
CL
CL 6.8 Mass weighted on q-levels for moisture and moisture variance
CL
      DO LEVEL=1,Q_LEVELS
        CALL ZONM(Q(1,LEVEL),Z_Q(1,LEVEL),AMASK,P_MASS(1,LEVEL),
     &            LAPTS(1,1),ROW_LENGTH,P_ROWS)

        DO I=1,P_FIELD
          WORK1(I,LEVEL)=Q(I,LEVEL)*Q(I,LEVEL)
        END DO

        CALL ZONM(WORK1(1,LEVEL),WORK2,
     &           AMASK,P_MASS(1,LEVEL),LAPTS(1,1),ROW_LENGTH,P_ROWS)

        DO I=1,P_ROWS
          Z_VAR_Q_LEV(I,LEVEL)=WORK2(I)-Z_Q(I,LEVEL)*Z_Q(I,LEVEL)
        END DO

      END DO

C Vertical mean moisture variance, mass weighted

      CALL COLM(WORK1,WORK2,P_MASS,P_ROWS,ROW_LENGTH,Q_LEVELS)
      CALL ZONM(WORK2,Z_VAR_Q,AMASK,PSTAR,LAPTS(1,1),ROW_LENGTH,P_ROWS)

      DO I=1,P_ROWS
        Z_VAR_Q(I)=Z_VAR_Q(I)-Z_QATMOS(I)*Z_QATMOS(I)
      END DO
CL
CL 6.9 Mass weighted on p/u-grid for surface variables (land points)
CL
      DO LEVEL=1,ST_LEVELS
        CALL ZONM(SOILT(1,LEVEL),Z_SOILT(1,LEVEL),LMASK,S_PMASS,
     &            LLPTS(1,1),ROW_LENGTH,P_ROWS)
      END DO

      CALL ZONM(SOILM,Z_SOILM,LMASK,S_PMASS,LLPTS(1,1),ROW_LENGTH,
     &                                                        P_ROWS)
      CALL ZONM(TSTAR,Z_L_TSTAR,LMASK,S_PMASS,LLPTS(1,1),ROW_LENGTH,
     &                                                        P_ROWS)
      CALL ZONM(CANOPYW,Z_CANOPYW,LMASK,S_PMASS,LLPTS(1,1),ROW_LENGTH,
     &                                                        P_ROWS)
      CALL ZONM(SH,Z_L_SH,LMASK,S_PMASS,LLPTS(1,1),ROW_LENGTH,P_ROWS)
      CALL ZONM(EVAP,Z_L_EVAP,LMASK,S_PMASS,LLPTS(1,1),ROW_LENGTH,
     &                                                        P_ROWS)
      CALL ZONM(PPTN,Z_L_PPTN,LMASK,S_PMASS,LLPTS(1,1),ROW_LENGTH,
     &                                                        P_ROWS)
      CALL ZONM(SNOW,Z_L_SNOW,LMASK,S_PMASS,LLPTS(1,1),ROW_LENGTH,
     &                                                        P_ROWS)
      CALL ZONM(TAUX,Z_L_TAUX,LMASK,S_UMASS,LLPTS(1,2),ROW_LENGTH,
     &                                                        U_ROWS)
      CALL ZONM(TAUY,Z_L_TAUY,LMASK,S_UMASS,LLPTS(1,2),ROW_LENGTH,
     &                                                        U_ROWS)
CL
CL 6.10 Mass weighted on p/u-grid for surface variables (sea points)
CL
      CALL ZONM(TSTAR,Z_S_TSTAR,SMASK,S_PMASS,LSPTS(1,1),ROW_LENGTH,
     &                                                        P_ROWS)
      CALL ZONM(SH,Z_S_SH,SMASK,S_PMASS,LSPTS(1,1),ROW_LENGTH,P_ROWS)
      CALL ZONM(EVAP,Z_S_EVAP,SMASK,S_PMASS,LSPTS(1,1),ROW_LENGTH,
     &                                                        P_ROWS)
      CALL ZONM(PPTN,Z_S_PPTN,SMASK,S_PMASS,LSPTS(1,1),ROW_LENGTH,
     &                                                        P_ROWS)
      CALL ZONM(SNOW,Z_S_SNOW,SMASK,S_PMASS,LSPTS(1,1),ROW_LENGTH,
     &                                                        P_ROWS)
      CALL ZONM(AICE,Z_S_AICE,SMASK,S_PMASS,LSPTS(1,1),ROW_LENGTH,
     &                                                        P_ROWS)
      CALL ZONM(HICE,Z_S_HICE,SMASK,S_PMASS,LSPTS(1,1),ROW_LENGTH,
     &                                                        P_ROWS)
      CALL ZONM(TAUX,Z_S_TAUX,SMASK,S_UMASS,LSPTS(1,2),ROW_LENGTH,
     &                                                        U_ROWS)
      CALL ZONM(TAUY,Z_S_TAUY,SMASK,S_UMASS,LSPTS(1,2),ROW_LENGTH,
     &                                                        U_ROWS)
CL----------------------------------------------------------------------
CL 7. Compute quarter global and global means
CL
      DO 90,I=1,4
        CALL GLBM(PMSL,G_Q_PMSL(I),START_ROW(1,I),END_ROW(1,I),
     &     COS_P_LATITUDE,AMASK,LQGAPTS(I,1),P_ROWS,ROW_LENGTH,.TRUE.)
        CALL GLBM(PSTAR,G_Q_PSTAR(I),START_ROW(1,I),END_ROW(1,I),
     &     COS_P_LATITUDE,AMASK,LQGAPTS(I,1),P_ROWS,ROW_LENGTH,.TRUE.)
        CALL GLBM(TSTAR,G_Q_TSTAR(I),START_ROW(1,I),END_ROW(1,I),
     &     COS_P_LATITUDE,AMASK,LQGAPTS(I,1),P_ROWS,ROW_LENGTH,.TRUE.)
        CALL GLBM(TSTAR,G_Q_L_TSTAR(I),START_ROW(1,I),END_ROW(1,I),
     &     COS_P_LATITUDE,LMASK,LQGLPTS(I,1),P_ROWS,ROW_LENGTH,.TRUE.)
        CALL GLBM(TSTAR,G_Q_S_TSTAR(I),START_ROW(1,I),END_ROW(1,I),
     &     COS_P_LATITUDE,SMASK,LQGSPTS(I,1),P_ROWS,ROW_LENGTH,.TRUE.)
        CALL GLBM(SOILM,G_Q_SOILM(I),START_ROW(1,I),END_ROW(1,I),
     &     COS_P_LATITUDE,LMASK,LQGLPTS(I,1),P_ROWS,ROW_LENGTH,.TRUE.)
        CALL GLBM(SNOWD,G_Q_SNOWD(I),START_ROW(1,I),END_ROW(1,I),
     &     COS_P_LATITUDE,AMASK,LQGAPTS(I,1),P_ROWS,ROW_LENGTH,.TRUE.)
        CALL GLBM(CANOPYW,G_Q_CANOPYW(I),START_ROW(1,I),END_ROW(1,I),
     &     COS_P_LATITUDE,LMASK,LQGLPTS(I,1),P_ROWS,ROW_LENGTH,.TRUE.)
        CALL GLBM(SH,G_Q_SH(I),START_ROW(1,I),END_ROW(1,I),
     &     COS_P_LATITUDE,AMASK,LQGAPTS(I,1),P_ROWS,ROW_LENGTH,.TRUE.)
        CALL GLBM(SH,G_Q_L_SH(I),START_ROW(1,I),END_ROW(1,I),
     &     COS_P_LATITUDE,LMASK,LQGAPTS(I,1),P_ROWS,ROW_LENGTH,.TRUE.)
        CALL GLBM(SH,G_Q_S_SH(I),START_ROW(1,I),END_ROW(1,I),
     &     COS_P_LATITUDE,SMASK,LQGAPTS(I,1),P_ROWS,ROW_LENGTH,.TRUE.)
        CALL GLBM(EVAP,G_Q_EVAP(I),START_ROW(1,I),END_ROW(1,I),
     &     COS_P_LATITUDE,AMASK,LQGAPTS(I,1),P_ROWS,ROW_LENGTH,.TRUE.)
        CALL GLBM(EVAP,G_Q_L_EVAP(I),START_ROW(1,I),END_ROW(1,I),
     &     COS_P_LATITUDE,LMASK,LQGAPTS(I,1),P_ROWS,ROW_LENGTH,.TRUE.)
        CALL GLBM(EVAP,G_Q_S_EVAP(I),START_ROW(1,I),END_ROW(1,I),
     &     COS_P_LATITUDE,SMASK,LQGAPTS(I,1),P_ROWS,ROW_LENGTH,.TRUE.)
        CALL GLBM(PPTN,G_Q_PPTN(I),START_ROW(1,I),END_ROW(1,I),
     &     COS_P_LATITUDE,AMASK,LQGAPTS(I,1),P_ROWS,ROW_LENGTH,.TRUE.)
        CALL GLBM(PPTN,G_Q_L_PPTN(I),START_ROW(1,I),END_ROW(1,I),
     &     COS_P_LATITUDE,LMASK,LQGAPTS(I,1),P_ROWS,ROW_LENGTH,.TRUE.)
        CALL GLBM(PPTN,G_Q_S_PPTN(I),START_ROW(1,I),END_ROW(1,I),
     &     COS_P_LATITUDE,SMASK,LQGAPTS(I,1),P_ROWS,ROW_LENGTH,.TRUE.)
        CALL GLBM(LSRN,G_Q_LSRN(I),START_ROW(1,I),END_ROW(1,I),
     &     COS_P_LATITUDE,AMASK,LQGAPTS(I,1),P_ROWS,ROW_LENGTH,.TRUE.)
        CALL GLBM(CVRN,G_Q_CVRN(I),START_ROW(1,I),END_ROW(1,I),
     &     COS_P_LATITUDE,AMASK,LQGAPTS(I,1),P_ROWS,ROW_LENGTH,.TRUE.)
        CALL GLBM(SNOW,G_Q_SNOW(I),START_ROW(1,I),END_ROW(1,I),
     &     COS_P_LATITUDE,AMASK,LQGAPTS(I,1),P_ROWS,ROW_LENGTH,.TRUE.)
        CALL GLBM(SNOW,G_Q_L_SNOW(I),START_ROW(1,I),END_ROW(1,I),
     &     COS_P_LATITUDE,LMASK,LQGAPTS(I,1),P_ROWS,ROW_LENGTH,.TRUE.)
        CALL GLBM(SNOW,G_Q_S_SNOW(I),START_ROW(1,I),END_ROW(1,I),
     &     COS_P_LATITUDE,SMASK,LQGAPTS(I,1),P_ROWS,ROW_LENGTH,.TRUE.)
        CALL GLBM(TATMOS,G_Q_TATMOS(I),START_ROW(1,I),END_ROW(1,I),
     &     COS_P_LATITUDE,AMASK,LQGAPTS(I,1),P_ROWS,ROW_LENGTH,.TRUE.)
        CALL GLBM(QATMOS,G_Q_QATMOS(I),START_ROW(1,I),END_ROW(1,I),
     &     COS_P_LATITUDE,AMASK,LQGAPTS(I,1),P_ROWS,ROW_LENGTH,.TRUE.)
        CALL GLBM(UATMOS,G_Q_UATMOS(I),START_ROW(2,I),END_ROW(2,I),
     &     COS_U_LATITUDE,AMASK,LQGAPTS(I,1),U_ROWS,ROW_LENGTH,.FALSE.)
        CALL GLBM(VATMOS,G_Q_VATMOS(I),START_ROW(2,I),END_ROW(2,I),
     &     COS_U_LATITUDE,AMASK,LQGAPTS(I,1),U_ROWS,ROW_LENGTH,.FALSE.)
        CALL GLBM(TAUX,G_Q_TAUX(I),START_ROW(2,I),END_ROW(2,I),
     &     COS_U_LATITUDE,AMASK,LQGAPTS(I,2),U_ROWS,ROW_LENGTH,.FALSE.)
        CALL GLBM(TAUY,G_Q_TAUY(I),START_ROW(2,I),END_ROW(2,I),
     &     COS_U_LATITUDE,AMASK,LQGAPTS(I,2),U_ROWS,ROW_LENGTH,.FALSE.)
        CALL GLBM(TAUX,G_Q_L_TAUX(I),START_ROW(2,I),END_ROW(2,I),
     &     COS_U_LATITUDE,LMASK,LQGLPTS(I,2),U_ROWS,ROW_LENGTH,.FALSE.)
        CALL GLBM(TAUY,G_Q_L_TAUY(I),START_ROW(2,I),END_ROW(2,I),
     &     COS_U_LATITUDE,LMASK,LQGLPTS(I,2),U_ROWS,ROW_LENGTH,.FALSE.)
        CALL GLBM(TAUX,G_Q_S_TAUX(I),START_ROW(2,I),END_ROW(2,I),
     &     COS_U_LATITUDE,SMASK,LQGSPTS(I,2),U_ROWS,ROW_LENGTH,.FALSE.)
        CALL GLBM(TAUY,G_Q_S_TAUY(I),START_ROW(2,I),END_ROW(2,I),
     &     COS_U_LATITUDE,SMASK,LQGSPTS(I,2),U_ROWS,ROW_LENGTH,.FALSE.)
90      CONTINUE

      DO LEVEL=1,ST_LEVELS
        DO I=1,4
          CALL GLBM(SOILT(1,LEVEL),G_Q_SOILT(I,LEVEL), START_ROW(1,I),
     &    END_ROW(1,I),COS_P_LATITUDE,LMASK, LQGLPTS(I,1),P_ROWS,
     &    ROW_LENGTH,.TRUE.)
        END DO
      END DO

CL 7.1 Winds and Kinetic energy at each level, mass weighted

      DO 120 LEVEL=1,P_LEVELS

        DO I=1,U_FIELD
          WORK3(I,LEVEL)=
     &        U(I,LEVEL)*U(I,LEVEL)+V(I,LEVEL)*V(I,LEVEL)
          WORK2(I)=Z_ZKE_LEV((I-1)/ROW_LENGTH+1,LEVEL)
        END DO

        CALL GLBM(WORK3(1,LEVEL),G_TKE_LEV(LEVEL),
     &       1,U_ROWS,COS_U_LATITUDE,
     &       U_MASS(1,LEVEL),LGAPTS(2),U_ROWS,ROW_LENGTH,.FALSE.)
        CALL GLBM(WORK2,G_ZKE_LEV(LEVEL),1,U_ROWS,COS_U_LATITUDE,
     &       U_MASS(1,LEVEL),LGAPTS(2),U_ROWS,ROW_LENGTH,.FALSE.)
        G_EKE_LEV(LEVEL)=G_TKE_LEV(LEVEL)-G_ZKE_LEV(LEVEL)

        DO 130,I=1,4
          CALL GLBM(U(1,LEVEL),G_Q_U(I,LEVEL),START_ROW(2,I),
     &        END_ROW(2,I),COS_U_LATITUDE,U_MASS(1,LEVEL),LQGAPTS(I,2),
     &                                    U_ROWS,ROW_LENGTH,.FALSE.)
          CALL GLBM(V(1,LEVEL),G_Q_V(I,LEVEL),START_ROW(2,I),
     &        END_ROW(2,I),COS_U_LATITUDE,U_MASS(1,LEVEL),LQGAPTS(I,2),
     &                                    U_ROWS,ROW_LENGTH,.FALSE.)
130     CONTINUE
120   CONTINUE

CL 7.2 Column mean energies, mass weighted

      CALL COLM(WORK3,WORK2,U_MASS,U_ROWS,ROW_LENGTH,P_LEVELS)
      DO I=1,4
        CALL GLBM(WORK2,G_Q_TKE(I),START_ROW(2,I),END_ROW(2,I),
     &    COS_U_LATITUDE,PSTAR_U,LQGAPTS(I,2),U_ROWS,ROW_LENGTH,.FALSE.)
      END DO
      CALL GLBM(WORK2,G_TKE,1,U_ROWS,
     &    COS_U_LATITUDE,PSTAR_U,LGAPTS(2),U_ROWS,ROW_LENGTH,.FALSE.)
      DO I=1,U_FIELD
        WORK2(I)=Z_ZKE((I-1)/ROW_LENGTH+1)
      END DO
      DO I=1,4
        CALL GLBM(WORK2,G_Q_ZKE(I),START_ROW(2,I),END_ROW(2,I),
     &   COS_U_LATITUDE,PSTAR_U,LQGAPTS(I,2),U_ROWS,ROW_LENGTH,.FALSE.)
        G_Q_EKE(I)=G_Q_TKE(I)-G_Q_ZKE(I)
      END DO
      CALL GLBM(WORK2,G_ZKE,1,U_ROWS,
     &  COS_U_LATITUDE,PSTAR_U,LGAPTS(2),U_ROWS,ROW_LENGTH,.FALSE.)
      G_EKE=G_TKE-G_ZKE

CL 7.3 Temperature and temperature variance at each level, mass weighted

      DO 140 LEVEL=1,P_LEVELS

        DO I=1,P_FIELD
          WORK1(I,LEVEL)=T(I,LEVEL)*T(I,LEVEL)
        END DO

        DO I=1,4
          CALL GLBM(T(1,LEVEL),G_Q_T(I,LEVEL),
     &          START_ROW(1,I),END_ROW(1,I),COS_P_LATITUDE,
     &          P_MASS(1,LEVEL),LQGAPTS(I,1),P_ROWS,ROW_LENGTH,.TRUE.)
        ENDDO

        CALL GLBM(T(1,LEVEL),G_T(LEVEL),1,P_ROWS,COS_P_LATITUDE,
     &            P_MASS(1,LEVEL),LGAPTS(1),P_ROWS,ROW_LENGTH,.TRUE.)
        CALL GLBM(WORK1(1,LEVEL),G_VAR_T_LEV(LEVEL),1,P_ROWS,
     &       COS_P_LATITUDE,P_MASS(1,LEVEL),LGAPTS(1),P_ROWS,ROW_LENGTH,
     &       .TRUE.)
        G_VAR_T_LEV(LEVEL)=G_VAR_T_LEV(LEVEL)-G_T(LEVEL)*G_T(LEVEL)

140   CONTINUE

CL 7.4 Column mean temperature and temperature variance, mass weighted

      CALL COLM(WORK1,WORK2,P_MASS,P_ROWS,ROW_LENGTH,P_LEVELS)
      CALL GLBM(TATMOS,G_TATMOS,1,P_ROWS,COS_P_LATITUDE,
     &            PSTAR,LGAPTS(1),P_ROWS,ROW_LENGTH,.TRUE.)
      CALL GLBM(WORK2,G_VAR_T,1,P_ROWS,
     &     COS_P_LATITUDE,PSTAR,LGAPTS(1),P_ROWS,ROW_LENGTH,.TRUE.)
      G_VAR_T=G_VAR_T-G_TATMOS*G_TATMOS

CL 7.5 Moisture and moisture variance at each level, mass weighted

      DO LEVEL=1,Q_LEVELS

        DO I=1,P_FIELD
          WORK1(I,LEVEL)=Q(I,LEVEL)*Q(I,LEVEL)
        END DO

        DO I=1,4
          CALL GLBM(Q(1,LEVEL),G_Q_Q(I,LEVEL),
     &         START_ROW(1,I),END_ROW(1,I),COS_P_LATITUDE,
     &         P_MASS(1,LEVEL),LQGAPTS(I,1),P_ROWS,ROW_LENGTH,.TRUE.)
        ENDDO

        CALL GLBM(Q(1,LEVEL),G_Q(LEVEL),1,P_ROWS,COS_P_LATITUDE,
     &            P_MASS(1,LEVEL),LGAPTS(1),P_ROWS,ROW_LENGTH,.TRUE.)
        CALL GLBM(WORK1(1,LEVEL),G_VAR_Q_LEV(LEVEL),1,P_ROWS,
     &       COS_P_LATITUDE,P_MASS(1,LEVEL),LGAPTS(1),P_ROWS,ROW_LENGTH,
     &       .TRUE.)
        G_VAR_Q_LEV(LEVEL)=G_VAR_Q_LEV(LEVEL)-G_Q(LEVEL)*G_Q(LEVEL)

      END DO

CL 7.6 Column mean moisture and moisture variance, mass weighted

      CALL COLM(WORK1,WORK2,P_MASS,P_ROWS,ROW_LENGTH,Q_LEVELS)
      CALL GLBM(QATMOS,G_QATMOS,1,P_ROWS,COS_P_LATITUDE,
     &            PSTAR,LGAPTS(1),P_ROWS,ROW_LENGTH,.TRUE.)
      CALL GLBM(WORK2,G_VAR_Q,1,P_ROWS,
     &         COS_P_LATITUDE,PSTAR,LGAPTS(1),P_ROWS,ROW_LENGTH,.TRUE.)
      G_VAR_Q=G_VAR_Q-G_QATMOS*G_QATMOS
CL
CL 7.7 Cloud water liquid and ice over whole atmosphere, put in rows
CL
      IF ((IPRTEXTRA.EQ.1) .AND. (ICLOUD.GT.1)) THEN
CL
      CALL COLM(CLLIQ,CLLIQ_F,P_MASS,P_ROWS,ROW_LENGTH,Q_LEVELS)
      CALL ZONM(CLLIQ_F,Z_CLLIQ,AMASK,PSTAR,LAPTS(1,1),ROW_LENGTH,
     &          P_ROWS)
      CALL COLM(CLICE,CLICE_F,P_MASS,P_ROWS,ROW_LENGTH,Q_LEVELS)
      CALL ZONM(CLICE_F,Z_CLICE,AMASK,PSTAR,LAPTS(1,1),ROW_LENGTH,
     &          P_ROWS)

      END IF
CL
CL----------------------------------------------------------------------
CL 8. More global means
CL
      CALL GLBM(PMSL,G_PMSL,1,P_ROWS,COS_P_LATITUDE,
     &            AMASK,LGAPTS(1),P_ROWS,ROW_LENGTH,.TRUE.)
      CALL GLBM(PSTAR,G_PSTAR,1,P_ROWS,COS_P_LATITUDE,
     &            AMASK,LGAPTS(1),P_ROWS,ROW_LENGTH,.TRUE.)
      CALL GLBM(TSTAR,G_TSTAR,1,P_ROWS,COS_P_LATITUDE,
     &            AMASK,LGAPTS(1),P_ROWS,ROW_LENGTH,.TRUE.)
      CALL GLBM(TSTAR,G_L_TSTAR,1,P_ROWS,COS_P_LATITUDE,
     &            LMASK,LGLPTS(1),P_ROWS,ROW_LENGTH,.TRUE.)
      CALL GLBM(TSTAR,G_S_TSTAR,1,P_ROWS,COS_P_LATITUDE,
     &            SMASK,LGSPTS(1),P_ROWS,ROW_LENGTH,.TRUE.)
      CALL GLBM(SOILM,G_SOILM,1,P_ROWS,COS_P_LATITUDE,
     &            LMASK,LGLPTS(1),P_ROWS,ROW_LENGTH,.TRUE.)
      CALL GLBM(SNOWD,G_SNOWD,1,P_ROWS,COS_P_LATITUDE,
     &            AMASK,LGAPTS(1),P_ROWS,ROW_LENGTH,.TRUE.)
      CALL GLBM(CANOPYW,G_CANOPYW,1,P_ROWS,COS_P_LATITUDE,
     &            LMASK,LGLPTS(1),P_ROWS,ROW_LENGTH,.TRUE.)
      CALL GLBM(SH,G_SH,1,P_ROWS,COS_P_LATITUDE,
     &            AMASK,LGAPTS(1),P_ROWS,ROW_LENGTH,.TRUE.)
      CALL GLBM(SH,G_L_SH,1,P_ROWS,COS_P_LATITUDE,
     &            LMASK,LGAPTS(1),P_ROWS,ROW_LENGTH,.TRUE.)
      CALL GLBM(SH,G_S_SH,1,P_ROWS,COS_P_LATITUDE,
     &            SMASK,LGAPTS(1),P_ROWS,ROW_LENGTH,.TRUE.)
      CALL GLBM(EVAP,G_EVAP,1,P_ROWS,COS_P_LATITUDE,
     &            AMASK,LGAPTS(1),P_ROWS,ROW_LENGTH,.TRUE.)
      CALL GLBM(EVAP,G_L_EVAP,1,P_ROWS,COS_P_LATITUDE,
     &            LMASK,LGAPTS(1),P_ROWS,ROW_LENGTH,.TRUE.)
      CALL GLBM(EVAP,G_S_EVAP,1,P_ROWS,COS_P_LATITUDE,
     &            SMASK,LGAPTS(1),P_ROWS,ROW_LENGTH,.TRUE.)
      CALL GLBM(PPTN,G_PPTN,1,P_ROWS,COS_P_LATITUDE,
     &            AMASK,LGAPTS(1),P_ROWS,ROW_LENGTH,.TRUE.)
      CALL GLBM(PPTN,G_L_PPTN,1,P_ROWS,COS_P_LATITUDE,
     &            LMASK,LGAPTS(1),P_ROWS,ROW_LENGTH,.TRUE.)
      CALL GLBM(PPTN,G_S_PPTN,1,P_ROWS,COS_P_LATITUDE,
     &            SMASK,LGAPTS(1),P_ROWS,ROW_LENGTH,.TRUE.)
      CALL GLBM(LSRN,G_LSRN,1,P_ROWS,COS_P_LATITUDE,
     &            AMASK,LGAPTS(1),P_ROWS,ROW_LENGTH,.TRUE.)
      CALL GLBM(CVRN,G_CVRN,1,P_ROWS,COS_P_LATITUDE,
     &            AMASK,LGAPTS(1),P_ROWS,ROW_LENGTH,.TRUE.)
      CALL GLBM(SNOW,G_SNOW,1,P_ROWS,COS_P_LATITUDE,
     &            AMASK,LGAPTS(1),P_ROWS,ROW_LENGTH,.TRUE.)
      CALL GLBM(SNOW,G_L_SNOW,1,P_ROWS,COS_P_LATITUDE,
     &            LMASK,LGAPTS(1),P_ROWS,ROW_LENGTH,.TRUE.)
      CALL GLBM(SNOW,G_S_SNOW,1,P_ROWS,COS_P_LATITUDE,
     &            SMASK,LGAPTS(1),P_ROWS,ROW_LENGTH,.TRUE.)
      CALL GLBM(UATMOS,G_UATMOS,1,U_ROWS,COS_U_LATITUDE,
     &            PSTAR_U,LGAPTS(2),U_ROWS,ROW_LENGTH,.FALSE.)
      CALL GLBM(VATMOS,G_VATMOS,1,U_ROWS,COS_U_LATITUDE,
     &            PSTAR_U,LGAPTS(2),U_ROWS,ROW_LENGTH,.FALSE.)
      CALL GLBM(TAUX,G_TAUX,1,U_ROWS,COS_U_LATITUDE,
     &            AMASK,LGAPTS(2),U_ROWS,ROW_LENGTH,.FALSE.)
      CALL GLBM(TAUY,G_TAUY,1,U_ROWS,COS_U_LATITUDE,
     &            AMASK,LGAPTS(2),U_ROWS,ROW_LENGTH,.FALSE.)
      CALL GLBM(TAUX,G_L_TAUX,1,U_ROWS,COS_U_LATITUDE,
     &            LMASK,LGLPTS(2),U_ROWS,ROW_LENGTH,.FALSE.)
      CALL GLBM(TAUY,G_L_TAUY,1,U_ROWS,COS_U_LATITUDE,
     &            LMASK,LGLPTS(2),U_ROWS,ROW_LENGTH,.FALSE.)
      CALL GLBM(TAUX,G_S_TAUX,1,U_ROWS,COS_U_LATITUDE,
     &            SMASK,LGSPTS(2),U_ROWS,ROW_LENGTH,.FALSE.)
      CALL GLBM(TAUY,G_S_TAUY,1,U_ROWS,COS_U_LATITUDE,
     &            SMASK,LGSPTS(2),U_ROWS,ROW_LENGTH,.FALSE.)

      DO LEVEL=1,P_LEVELS
        CALL GLBM(U(1,LEVEL),G_U(LEVEL),1,U_ROWS,COS_U_LATITUDE,
     &            U_MASS(1,LEVEL),LGAPTS(2),U_ROWS,ROW_LENGTH,.FALSE.)
        CALL GLBM(V(1,LEVEL),G_V(LEVEL),1,U_ROWS,COS_U_LATITUDE,
     &            U_MASS(1,LEVEL),LGAPTS(2),U_ROWS,ROW_LENGTH,.FALSE.)
      END DO
C
      DO LEVEL=1,ST_LEVELS
          CALL GLBM(SOILT(1,LEVEL),G_SOILT(LEVEL),1,P_ROWS,
     &        COS_P_LATITUDE,LMASK,LGLPTS(1),P_ROWS,ROW_LENGTH,.TRUE.)
      END DO
C
      CALL ZONM(SUBL,Z_L_SUBL,LMASK,S_PMASS,LLPTS(1,1),ROW_LENGTH,
     &                                                        P_ROWS)
      CALL ZONM(SUBL,Z_S_SUBL,SMASK,S_PMASS,LSPTS(1,1),ROW_LENGTH,
     &                                                        P_ROWS)
      DO I=1,4
        CALL GLBM(SUBL,G_Q_L_SUBL(I),START_ROW(1,I),END_ROW(1,I),
     &     COS_P_LATITUDE,LMASK,LQGAPTS(I,1),P_ROWS,ROW_LENGTH,.TRUE.)
        CALL GLBM(SUBL,G_Q_SUBL(I),START_ROW(1,I),END_ROW(1,I),
     &     COS_P_LATITUDE,AMASK,LQGAPTS(I,1),P_ROWS,ROW_LENGTH,.TRUE.)
        CALL GLBM(SUBL,G_Q_S_SUBL(I),START_ROW(1,I),END_ROW(1,I),
     &     COS_P_LATITUDE,SMASK,LQGAPTS(I,1),P_ROWS,ROW_LENGTH,.TRUE.)
      ENDDO
      CALL GLBM(SUBL,G_SUBL,1,P_ROWS,COS_P_LATITUDE,
     &     AMASK,LGAPTS(1),P_ROWS,ROW_LENGTH,.TRUE.)
      CALL GLBM(SUBL,G_L_SUBL,1,P_ROWS,COS_P_LATITUDE,
     &            LMASK,LGAPTS(1),P_ROWS,ROW_LENGTH,.TRUE.)
      CALL GLBM(SUBL,G_S_SUBL,1,P_ROWS,COS_P_LATITUDE,
     &            SMASK,LGAPTS(1),P_ROWS,ROW_LENGTH,.TRUE.)
      DO I=1,4
        CALL GLBM(SOEV,G_Q_L_SOEV(I),START_ROW(1,I),END_ROW(1,I),
     &     COS_P_LATITUDE,LMASK,LQGAPTS(I,1),P_ROWS,ROW_LENGTH,.TRUE.)
      ENDDO
      CALL GLBM(SOEV,G_L_SOEV,1,P_ROWS,COS_P_LATITUDE,
     &            LMASK,LGAPTS(1),P_ROWS,ROW_LENGTH,.TRUE.)
C
CL----------------------------------------------------------------------
CL
CL 8.a Additional zonal means ONLY calculated if field available
CL
      IF (IPRTEXTRA.EQ.1) THEN
       IF (IHYDRO.GT.1) THEN
        CALL ZONM(SFRU,Z_L_SFRU,LMASK,S_PMASS,LLPTS(1,1),ROW_LENGTH,
     &                                                        P_ROWS)
        CALL ZONM(SBRU,Z_L_SBRU,LMASK,S_PMASS,LLPTS(1,1),ROW_LENGTH,
     &                                                        P_ROWS)
        DO I=1,4
          CALL GLBM(THRF,G_Q_L_THRF(I),START_ROW(1,I),END_ROW(1,I),
     &      COS_P_LATITUDE,LMASK,LQGAPTS(I,1),P_ROWS,ROW_LENGTH,.TRUE.)
          CALL GLBM(SNML,G_Q_L_SNML(I),START_ROW(1,I),END_ROW(1,I),
     &      COS_P_LATITUDE,LMASK,LQGAPTS(I,1),P_ROWS,ROW_LENGTH,.TRUE.)
          CALL GLBM(SFRU,G_Q_L_SFRU(I),START_ROW(1,I),END_ROW(1,I),
     &      COS_P_LATITUDE,LMASK,LQGAPTS(I,1),P_ROWS,ROW_LENGTH,.TRUE.)
          CALL GLBM(SBRU,G_Q_L_SBRU(I),START_ROW(1,I),END_ROW(1,I),
     &      COS_P_LATITUDE,LMASK,LQGAPTS(I,1),P_ROWS,ROW_LENGTH,.TRUE.)
        ENDDO
        CALL GLBM(SFRU,G_L_SFRU,1,P_ROWS,COS_P_LATITUDE,
     &            LMASK,LGAPTS(1),P_ROWS,ROW_LENGTH,.TRUE.)
        CALL GLBM(SBRU,G_L_SBRU,1,P_ROWS,COS_P_LATITUDE,
     &            LMASK,LGAPTS(1),P_ROWS,ROW_LENGTH,.TRUE.)
        CALL GLBM(THRF,G_L_THRF,1,P_ROWS,COS_P_LATITUDE,
     &            LMASK,LGAPTS(1),P_ROWS,ROW_LENGTH,.TRUE.)
        CALL GLBM(SNML,G_L_SNML,1,P_ROWS,COS_P_LATITUDE,
     &            LMASK,LGAPTS(1),P_ROWS,ROW_LENGTH,.TRUE.)
       ENDIF
       IF (IRAD.GT.1) THEN   ! RADIATION
        CALL ZONM(SDTR,Z_L_SDTR,LMASK,S_PMASS,LLPTS(1,1),ROW_LENGTH,
     &                                                        P_ROWS)
        CALL ZONM(SDTR,Z_S_SDTR,SMASK,S_PMASS,LSPTS(1,1),ROW_LENGTH,
     &                                                        P_ROWS)
        CALL ZONM(SDSR,Z_L_SDSR,LMASK,S_PMASS,LLPTS(1,1),ROW_LENGTH,
     &                                                        P_ROWS)
        CALL ZONM(SDSR,Z_S_SDSR,SMASK,S_PMASS,LSPTS(1,1),ROW_LENGTH,
     &                                                        P_ROWS)
        CALL ZONM(TDTR,Z_TDTR,AMASK,S_PMASS,LAPTS(1,1),ROW_LENGTH,
     &                                                        P_ROWS)
        CALL ZONM(TOLR,Z_TOLR,AMASK,S_PMASS,LAPTS(1,1),ROW_LENGTH,
     &                                                        P_ROWS)
        CALL ZONM(TOSW,Z_TOSW,AMASK,S_PMASS,LAPTS(1,1),ROW_LENGTH,
     &                                                        P_ROWS)
        CALL ZONM(TISW,Z_TISW,AMASK,S_PMASS,LAPTS(1,1),ROW_LENGTH,
     &                                                        P_ROWS)
        DO I=1,P_ROWS
           IF (Z_TISW(I) .GT. 10.0E-10) THEN
              Z_ALBEDO(I)=Z_TOSW(I)/Z_TISW(I)
           ELSE
              Z_ALBEDO(I) = 0.0
           END IF
        END DO
        DO I=1,4
          CALL GLBM(SDTR,G_Q_L_SDTR(I),START_ROW(1,I),END_ROW(1,I),
     &     COS_P_LATITUDE,LMASK,LQGAPTS(I,1),P_ROWS,ROW_LENGTH,.TRUE.)
          CALL GLBM(SDTR,G_Q_SDTR(I),START_ROW(1,I),END_ROW(1,I),
     &     COS_P_LATITUDE,AMASK,LQGAPTS(I,1),P_ROWS,ROW_LENGTH,.TRUE.)
          CALL GLBM(SDTR,G_Q_L_SDTR(I),START_ROW(1,I),END_ROW(1,I),
     &     COS_P_LATITUDE,LMASK,LQGAPTS(I,1),P_ROWS,ROW_LENGTH,.TRUE.)
          CALL GLBM(SDTR,G_Q_S_SDTR(I),START_ROW(1,I),END_ROW(1,I),
     &     COS_P_LATITUDE,SMASK,LQGAPTS(I,1),P_ROWS,ROW_LENGTH,.TRUE.)
          CALL GLBM(SDSR,G_Q_SDSR(I),START_ROW(1,I),END_ROW(1,I),
     &     COS_P_LATITUDE,AMASK,LQGAPTS(I,1),P_ROWS,ROW_LENGTH,.TRUE.)
          CALL GLBM(SDSR,G_Q_L_SDSR(I),START_ROW(1,I),END_ROW(1,I),
     &     COS_P_LATITUDE,LMASK,LQGAPTS(I,1),P_ROWS,ROW_LENGTH,.TRUE.)
          CALL GLBM(SDSR,G_Q_S_SDSR(I),START_ROW(1,I),END_ROW(1,I),
     &     COS_P_LATITUDE,SMASK,LQGAPTS(I,1),P_ROWS,ROW_LENGTH,.TRUE.)
          CALL GLBM(TDTR,G_Q_TDTR(I),START_ROW(1,I),END_ROW(1,I),
     &     COS_P_LATITUDE,AMASK,LQGAPTS(I,1),P_ROWS,ROW_LENGTH,.TRUE.)
          CALL GLBM(TOLR,G_Q_TOLR(I),START_ROW(1,I),END_ROW(1,I),
     &     COS_P_LATITUDE,AMASK,LQGAPTS(I,1),P_ROWS,ROW_LENGTH,.TRUE.)
          CALL GLBM(TOSW,G_Q_TOSW(I),START_ROW(1,I),END_ROW(1,I),
     &     COS_P_LATITUDE,AMASK,LQGAPTS(I,1),P_ROWS,ROW_LENGTH,.TRUE.)
          CALL GLBM(TISW,G_Q_TISW(I),START_ROW(1,I),END_ROW(1,I),
     &     COS_P_LATITUDE,AMASK,LQGAPTS(I,1),P_ROWS,ROW_LENGTH,.TRUE.)
        G_Q_ALBEDO(I)=G_Q_TOSW(I)/G_Q_TISW(I)
        ENDDO
        CALL GLBM(SDTR,G_L_SDTR,1,P_ROWS,COS_P_LATITUDE,
     &            LMASK,LGAPTS(1),P_ROWS,ROW_LENGTH,.TRUE.)
        CALL GLBM(SDTR,G_SDTR,1,P_ROWS,COS_P_LATITUDE,
     &            AMASK,LGAPTS(1),P_ROWS,ROW_LENGTH,.TRUE.)
        CALL GLBM(SDTR,G_L_SDTR,1,P_ROWS,COS_P_LATITUDE,
     &            LMASK,LGAPTS(1),P_ROWS,ROW_LENGTH,.TRUE.)
        CALL GLBM(SDTR,G_S_SDTR,1,P_ROWS,COS_P_LATITUDE,
     &            SMASK,LGAPTS(1),P_ROWS,ROW_LENGTH,.TRUE.)
        CALL GLBM(SDSR,G_SDSR,1,P_ROWS,COS_P_LATITUDE,
     &            AMASK,LGAPTS(1),P_ROWS,ROW_LENGTH,.TRUE.)
        CALL GLBM(SDSR,G_L_SDSR,1,P_ROWS,COS_P_LATITUDE,
     &            LMASK,LGAPTS(1),P_ROWS,ROW_LENGTH,.TRUE.)
        CALL GLBM(SDSR,G_S_SDSR,1,P_ROWS,COS_P_LATITUDE,
     &            SMASK,LGAPTS(1),P_ROWS,ROW_LENGTH,.TRUE.)
        CALL GLBM(TDTR,G_TDTR,1,P_ROWS,COS_P_LATITUDE,
     &            AMASK,LGAPTS(1),P_ROWS,ROW_LENGTH,.TRUE.)
        CALL GLBM(TOLR,G_TOLR,1,P_ROWS,COS_P_LATITUDE,
     &            AMASK,LGAPTS(1),P_ROWS,ROW_LENGTH,.TRUE.)
        CALL GLBM(TOSW,G_TOSW,1,P_ROWS,COS_P_LATITUDE,
     &            AMASK,LGAPTS(1),P_ROWS,ROW_LENGTH,.TRUE.)
        CALL GLBM(TISW,G_TISW,1,P_ROWS,COS_P_LATITUDE,
     &            AMASK,LGAPTS(1),P_ROWS,ROW_LENGTH,.TRUE.)
        G_ALBEDO=G_TOSW/G_TISW
       ENDIF
CL
CL    Fluxes over land - except for energy flux ENFL
CL
       IF (IFLUXL.GT.1) THEN
        DO I=1,4
          CALL GLBM(SOHF,G_Q_L_SOHF(I),START_ROW(1,I),END_ROW(1,I),
     &     COS_P_LATITUDE,LMASK,LQGAPTS(I,1),P_ROWS,ROW_LENGTH,.TRUE.)
          G_Q_L_WFSS(I)=G_Q_L_SNOW(I)-G_Q_L_SUBL(I)-G_Q_L_SNML(I)
          CALL GLBM(WFCA,G_Q_L_WFCA(I),START_ROW(1,I),END_ROW(1,I),
     &     COS_P_LATITUDE,LMASK,LQGAPTS(I,1),P_ROWS,ROW_LENGTH,.TRUE.)
          G_Q_L_WFSO(I)=G_Q_L_THRF(I)+G_Q_L_SNML(I)-G_Q_L_SOEV(I)
     &            -G_Q_L_SFRU(I)-G_Q_L_SBRU(I)
          G_Q_WAFL(I)=G_Q_EVAP(I)*86400.-G_Q_PPTN(I)
          G_Q_L_ENFS(I)=G_Q_L_SDTR(I)-G_Q_L_SH(I) - LC*G_Q_L_EVAP(I)
     &     -LF*(G_Q_L_SUBL(I)+G_Q_L_SNML(I))/86400.
          G_Q_ENFL(I)=G_Q_TDTR(I)-G_Q_SDTR(I)+G_Q_SH(I)
     &   +(LC*G_Q_PPTN(I)+LF*G_Q_SNOW(I))/86400.
        ENDDO
        CALL GLBM(SOHF,G_L_SOHF,1,P_ROWS,COS_P_LATITUDE,
     &            LMASK,LGAPTS(1),P_ROWS,ROW_LENGTH,.TRUE.)
        G_L_WFSS=G_L_SNOW-G_L_SUBL-G_L_SNML
        CALL GLBM(WFCA,G_L_WFCA,1,P_ROWS,COS_P_LATITUDE,
     &            LMASK,LGAPTS(1),P_ROWS,ROW_LENGTH,.TRUE.)
        G_L_WFSO=G_L_THRF+G_L_SNML-G_L_SOEV-G_L_SFRU-G_L_SBRU
        G_WAFL=G_EVAP*86400.-G_PPTN
        G_L_ENFS=G_L_SDTR-G_L_SH - LC*G_L_EVAP
     &     -LF*(G_L_SUBL+G_L_SNML)/86400.
        G_ENFL=G_TDTR-G_SDTR+G_SH+(LC*G_PPTN+LF*G_SNOW)/86400.
       ENDIF
CL
CL   Cloud water quarter globe means
CL
       IF (ICLOUD.GT.1) THEN
        DO LEVEL=1,Q_LEVELS
         DO I=1,4
          CALL GLBM(CLLIQ(1,LEVEL),G_Q_CLLIQL(I,LEVEL),
     &         START_ROW(1,I),END_ROW(1,I),COS_P_LATITUDE,
     &         P_MASS(1,LEVEL),LQGAPTS(I,1),P_ROWS,ROW_LENGTH,.TRUE.)
          CALL GLBM(CLICE(1,LEVEL),G_Q_CLICEL(I,LEVEL),
     &         START_ROW(1,I),END_ROW(1,I),COS_P_LATITUDE,
     &         P_MASS(1,LEVEL),LQGAPTS(I,1),P_ROWS,ROW_LENGTH,.TRUE.)
         ENDDO

        CALL GLBM(CLLIQ(1,LEVEL),G_CLLIQL(LEVEL),1,P_ROWS,
     &      COS_P_LATITUDE,P_MASS(1,LEVEL),LGAPTS(1),P_ROWS,ROW_LENGTH,
     &      .TRUE.)
        CALL GLBM(CLICE(1,LEVEL),G_CLICEL(LEVEL),1,P_ROWS,
     &      COS_P_LATITUDE,P_MASS(1,LEVEL),LGAPTS(1),P_ROWS,ROW_LENGTH,
     &      .TRUE.)

        END DO
        CALL COLM(CLLIQ,WORK2,P_MASS,P_ROWS,ROW_LENGTH,Q_LEVELS)
        DO I=1,4
          CALL GLBM(WORK2,G_Q_CLLIQ(I),
     &              START_ROW(1,I),END_ROW(1,I),COS_P_LATITUDE,
     &              PSTAR,LQGAPTS(I,1),P_ROWS,ROW_LENGTH,.TRUE.)
        ENDDO
        CALL GLBM(WORK2,G_CLLIQ,1,P_ROWS,
     &        COS_P_LATITUDE,PSTAR,LGAPTS(1),P_ROWS,ROW_LENGTH,.TRUE.)
        CALL COLM(CLICE,WORK2,P_MASS,P_ROWS,ROW_LENGTH,Q_LEVELS)
        DO I=1,4
          CALL GLBM(WORK2,G_Q_CLICE(I),
     &              START_ROW(1,I),END_ROW(1,I),COS_P_LATITUDE,
     &              PSTAR,LQGAPTS(I,1),P_ROWS,ROW_LENGTH,.TRUE.)
        ENDDO
        CALL GLBM(WORK2,G_CLICE,1,P_ROWS,
     &        COS_P_LATITUDE,PSTAR,LGAPTS(1),P_ROWS,ROW_LENGTH,.TRUE.)
       ENDIF
CL
CL    Fluxes over Sea-ice
CL
       IF (ISICE.GT.1) THEN
        DO I=1,4
          CALL GLBM(SIHF,G_Q_SIHF(I),START_ROW(1,I),END_ROW(1,I),
     &     COS_P_LATITUDE,SMASK,LQGAPTS(I,1),P_ROWS,ROW_LENGTH,.TRUE.)
          CALL GLBM(SIMH,G_Q_SIMH(I),START_ROW(1,I),END_ROW(1,I),
     &     COS_P_LATITUDE,SMASK,LQGAPTS(I,1),P_ROWS,ROW_LENGTH,.TRUE.)
          CALL GLBM(SISS,G_Q_SISS(I),START_ROW(1,I),END_ROW(1,I),
     &     COS_P_LATITUDE,SMASK,LQGAPTS(I,1),P_ROWS,ROW_LENGTH,.TRUE.)
          CALL GLBM(SIST,G_Q_SIST(I),START_ROW(1,I),END_ROW(1,I),
     &     COS_P_LATITUDE,SMASK,LQGAPTS(I,1),P_ROWS,ROW_LENGTH,.TRUE.)
          CALL GLBM(SIEF,G_Q_SIEF(I),START_ROW(1,I),END_ROW(1,I),
     &     COS_P_LATITUDE,SMASK,LQGAPTS(I,1),P_ROWS,ROW_LENGTH,.TRUE.)
          CALL GLBM(SISH,G_Q_SISH(I),START_ROW(1,I),END_ROW(1,I),
     &     COS_P_LATITUDE,SMASK,LQGAPTS(I,1),P_ROWS,ROW_LENGTH,.TRUE.)
        ENDDO
        CALL GLBM(SIHF,G_SIHF,1,P_ROWS,COS_P_LATITUDE,
     &            SMASK,LGAPTS(1),P_ROWS,ROW_LENGTH,.TRUE.)
        CALL GLBM(SIMH,G_SIMH,1,P_ROWS,COS_P_LATITUDE,
     &            SMASK,LGAPTS(1),P_ROWS,ROW_LENGTH,.TRUE.)
        CALL GLBM(SISS,G_SISS,1,P_ROWS,COS_P_LATITUDE,
     &            SMASK,LGAPTS(1),P_ROWS,ROW_LENGTH,.TRUE.)
        CALL GLBM(SIST,G_SIST,1,P_ROWS,COS_P_LATITUDE,
     &            SMASK,LGAPTS(1),P_ROWS,ROW_LENGTH,.TRUE.)
        CALL GLBM(SIEF,G_SIEF,1,P_ROWS,COS_P_LATITUDE,
     &            SMASK,LGAPTS(1),P_ROWS,ROW_LENGTH,.TRUE.)
        CALL GLBM(SISH,G_SISH,1,P_ROWS,COS_P_LATITUDE,
     &            SMASK,LGAPTS(1),P_ROWS,ROW_LENGTH,.TRUE.)
       ENDIF
      ENDIF

CL----------------------------------------------------------------------
CL 9. Print results neatly formatted
CL
      WRITE(6,'(1X)')
      WRITE(6,'(A80/)') BANNER
      WRITE(6,115)
115   FORMAT(T2,'ROW SEA LAND  PMSL    PSTAR     SNOWD  CANOPYW LSRN',
     &T54,'CVRN TATMOS QATMOS SOILM  T*SEA  T*LAND SOILT1 SOILT2 ',
     &T108,'SOILT3 SOILT4',T122,'SRUN SUBRUN')
      WRITE(6,125)
125   FORMAT(T10,'ONLY',T16,'MB',T24,'MB',T34,'MM',T42,'MM',T50,'MM/DAY'
     &,T59,'DEG.K',T66,'GM/KG',T74,'MM',T80,'DEG.K',T87,'DEG.K',T94,
     &'DEG.K',T101,'DEG.K',T108,'DEG.K',T115,'DEG.K',T123,'MM/DAY'/)

      DO ROW=1,P_ROWS
       IF (IHYDRO.GT.1) THEN
        WRITE(6,135)ROW,SPTS(ROW,1),LPTS(ROW,1),(Z_PMSL(ROW)*0.01),
     &(Z_PSTAR(ROW)*0.01),Z_SNOWD(ROW),Z_CANOPYW(ROW),Z_LSRN(ROW),
     &Z_CVRN(ROW),Z_TATMOS(ROW),(Z_QATMOS(ROW)*1000.),Z_SOILM(ROW),
     &Z_S_TSTAR(ROW),Z_L_TSTAR(ROW),(Z_SOILT(ROW,LEVEL),
     & LEVEL=1,ST_LEVELS),Z_L_SFRU(ROW),Z_L_SBRU(ROW)
       ELSE
        WRITE(6,135)ROW,SPTS(ROW,1),LPTS(ROW,1),(Z_PMSL(ROW)*0.01),
     &(Z_PSTAR(ROW)*0.01),Z_SNOWD(ROW),Z_CANOPYW(ROW),Z_LSRN(ROW),
     &Z_CVRN(ROW),Z_TATMOS(ROW),(Z_QATMOS(ROW)*1000.),Z_SOILM(ROW),
     &Z_S_TSTAR(ROW),Z_L_TSTAR(ROW),(Z_SOILT(ROW,LEVEL),
     &                                       LEVEL=1,ST_LEVELS)
       END IF
        IF (MOD(ROW,5).EQ.0) WRITE(6,'(1X)')
      END DO

135   FORMAT(T2,I3,T6,I3,T10,I3,T15,F7.2,T23,F7.2,T31,F9.2,T41,F5.2,T47,
     &F5.3,T53,F5.3,T59,F6.2,T66,F5.3,T73,F6.2,T80,F6.2,T87,F6.2,T94,
     &F6.2,T101,F6.2,T108,F6.2,T115,F6.2,T122,2F5.2)
      WRITE(6,'(1X)')
      WRITE(6,'(A80/)') BANNER
      WRITE(6,215)
215   FORMAT(T2,'ROW SEA LAND',T16,'S.HEAT  (W/M2)',T34,'EVAP  (MM/DAY)'
     &,T52,'PPTN  (MM/DAY)',T70,'SNOW  (MM/DAY)',T87,'AICE'
     &,T93,'HICE(M)',T101,'SUBLM(mm/day)',T115,'Cloud water whole')
      WRITE(6,220)
220   FORMAT(T115,'atmos (g/kg*1000)')
      WRITE(6,225)
225   FORMAT(T16,'ALL',T22,'SEA',T27,'LAND',T34,'ALL',T40,'SEA',T45,
     &'LAND',T52,'ALL',T57,'SEA',T63,'LAND',T70,'ALL',T76,'SEA',T81,
     &'LAND',T88,'SEA',T95,'SEA',T102,'SEA',T108,'LAND',T115,'LIQUID',
     & T123,'ICE'/)

      DO ROW=1,P_ROWS
        IF (ICLOUD.GT.1) THEN
         WRITE(6,235)ROW,SPTS(ROW,1),LPTS(ROW,1),Z_SH(ROW),Z_S_SH(ROW),
     &Z_L_SH(ROW),(Z_EVAP(ROW)*86400.),(Z_S_EVAP(ROW)*86400.),
     &(Z_L_EVAP(ROW)*86400.),Z_PPTN(ROW),
     &Z_S_PPTN(ROW),Z_L_PPTN(ROW),Z_SNOW(ROW),Z_S_SNOW(ROW),
     &Z_L_SNOW(ROW),Z_S_AICE(ROW),Z_S_HICE(ROW),Z_S_SUBL(ROW),
     &Z_L_SUBL(ROW),1000000.*Z_CLLIQ(ROW),1000000.*Z_CLICE(ROW)
        ELSE
         WRITE(6,235)ROW,SPTS(ROW,1),LPTS(ROW,1),Z_SH(ROW),Z_S_SH(ROW),
     &Z_L_SH(ROW),(Z_EVAP(ROW)*86400.),(Z_S_EVAP(ROW)*86400.),
     &(Z_L_EVAP(ROW)*86400.),Z_PPTN(ROW),
     &Z_S_PPTN(ROW),Z_L_PPTN(ROW),Z_SNOW(ROW),Z_S_SNOW(ROW),
     &Z_L_SNOW(ROW),Z_S_AICE(ROW),Z_S_HICE(ROW),Z_S_SUBL(ROW),
     &Z_L_SUBL(ROW)
        END IF
        IF (MOD(ROW,5).EQ.0) WRITE(6,'(1X)')
      END DO
235   FORMAT(T2,I3,T6,I3,T10,I3,T15,F5.1,T21,F5.1,T27,F5.1,T33,F5.3,T39,
     &F5.3,T45,F5.3,T51,F5.3,T57,F5.3,T63,F5.3,T69,F5.3,T75,F5.3,T81,
     &F5.3,T87,F5.3,T94,F5.3,T100,2(1X,F6.3),T115,2(F6.2,1X))
C
      IF ((IPRTWIND.EQ.1) .AND. (IRAD.EQ.1)) THEN
        WRITE(6,'(1X)')
        WRITE(6,'(A80/)') BANNER
        WRITE(6,315)
315   FORMAT(T2,'ROW SEA LAND UATMOS VATMOS',T35,'TAUX  (N/M2)',T56,
     &'TAUY  (N/M2)')
        WRITE(6,325)
325   FORMAT(T18,'M/S',T25,'M/S',T31,'ALL',T38,'SEA',T45,'LAND',T52,
     &'ALL',T59,'SEA',T66,'LAND'/)

        DO ROW=1,U_ROWS
          WRITE(6,330)ROW,SPTS(ROW,2),LPTS(ROW,2),Z_UATMOS(ROW),
     &Z_VATMOS(ROW),Z_TAUX(ROW),Z_S_TAUX(ROW),Z_L_TAUX(ROW),Z_TAUY(ROW),
     &Z_S_TAUY(ROW),Z_L_TAUY(ROW)
        END DO
330     FORMAT(T2,I3,T6,I3,T10,I3,T15,F6.2,T22,F6.2,T29,6F7.3)
      END IF
C
      IF ((IPRTWIND.NE.1) .AND. (IRAD.GT.1)) THEN
        WRITE(6,'(1X)')
        WRITE(6,'(A80/)') BANNER
        WRITE(6,331)
331   FORMAT(T2,'ROW TOTAL SURFACE  SURFACE SOLAR  TOTAL TOA  LW TOA',
     &'    ALBEDO')
        WRITE(6,332)
332   FORMAT(T7,2('Land',4X,'Sea',4X),2X,2('rad',7X))

        DO ROW=1,U_ROWS
          WRITE(6,333)ROW,Z_L_SDTR(ROW),Z_S_SDTR(ROW),Z_L_SDSR(ROW),
     &Z_S_SDSR(ROW),Z_TDTR(ROW),Z_TOLR(ROW),Z_ALBEDO(ROW)
        END DO
333   FORMAT(T2,I3,1X,2(2(F6.1,1X),1X),1X,2(F7.2,2X),1X,F7.4)
      END IF

      IF ((IPRTWIND.EQ.1) .AND. (IRAD.GT.1)) THEN
        WRITE(6,'(1X)')
        WRITE(6,'(A80/)') BANNER
        WRITE(6,335)
335   FORMAT(T2,'ROW SEA LAND UATMOS VATMOS',T35,'TAUX  (N/M2)',T56,
     &'TAUY  (N/M2)',T72,'TOTAL SURFACE  SURFACE SOLAR  TOTAL TOA',
     &'  LW TOA    ALBEDO')
        WRITE(6,340)
340   FORMAT(T18,'M/S',T25,'M/S',T31,'ALL',T38,'SEA',T45,'LAND',T52,
     &'ALL',T59,'SEA',T66,'LAND',T73,2('Land',4X,'Sea',4X),
     & 2X,2('rad',7X))

        DO ROW=1,U_ROWS
          WRITE(6,345)ROW,SPTS(ROW,2),LPTS(ROW,2),Z_UATMOS(ROW),
     &Z_VATMOS(ROW),Z_TAUX(ROW),Z_S_TAUX(ROW),Z_L_TAUX(ROW),Z_TAUY(ROW),
     &Z_S_TAUY(ROW),Z_L_TAUY(ROW),
     &Z_L_SDTR(ROW),Z_S_SDTR(ROW),Z_L_SDSR(ROW),
     &Z_S_SDSR(ROW),Z_TDTR(ROW),Z_TOLR(ROW),Z_ALBEDO(ROW)


        END DO
345     FORMAT(T2,I3,T6,I3,T10,I3,T15,F6.2,T22,F6.2,T29,6F7.3,
     &  T72,2(2(F6.1,1X),1X),1X,2(F7.2,2X),1X,F7.4)
      END IF
C
      IF (IPRTWIND.EQ.1) THEN
        WRITE(6,'(1X)')
        WRITE(6,'(A80/)') BANNER
        WRITE(6,415)
415   FORMAT(T2,'ROW   U1',T14,'U2',T20,'U3',T26,'U4',T32,'U5',T38,'U6',
     &T44,'U7',T50,'U8',T56,'U9',T62,'U10',T68,'U11',T74,'U12',T80,
     &'U13',T86,'U14',T92,'U15',T98,'U16',T104,'U17',T110,'U18',T116,
     &'U19',T122,'U20')
        WRITE(6,425)
425   FORMAT(T7,'M/S',T13,'M/S',T19,'M/S',T25,'M/S',T31,'M/S',
     &T37,'M/S',T43,'M/S',T49,'M/S',T55,'M/S',T61,'M/S',T67,
     &'M/S',T73,'M/S',T79,'M/S',T85,'M/S',T91,'M/S',T97,
     &'M/S',T103,'M/S',T109,'M/S',T115,'M/S',T121,'M/S'/)
        DO ROW=1,U_ROWS
          WRITE(6,435)ROW,(Z_U(ROW,LEVEL), LEVEL=1,P_LEVELS)
          IF (MOD(ROW,5).EQ.0) WRITE(6,'(1X)')
        END DO
435     FORMAT(T1,I3,T5,F5.2,19F6.2)
c
        WRITE(6,'(1X)')
        WRITE(6,'(A80/)') BANNER
        WRITE(6,515)
515   FORMAT(T2,'ROW   V1',T14,'V2',T20,'V3',T26,'V4',T32,'V5',T38,'V6',
     &T44,'V7',T50,'V8',T56,'V9',T62,'V10',T68,'V11',T74,'V12',T80,
     &'V13',T86,'V14',T92,'V15',T98,'V16',T104,'V17',T110,'V18',T116,
     &'V19',T122,'V20')
        WRITE(6,525)
525   FORMAT(T7,'M/S',T13,'M/S',T19,'M/S',T25,'M/S',T31,'M/S',
     &T37,'M/S',T43,'M/S',T49,'M/S',T55,'M/S',T61,'M/S',T67,
     &'M/S',T73,'M/S',T79,'M/S',T85,'M/S',T91,'M/S',T97,
     &'M/S',T103,'M/S',T109,'M/S',T115,'M/S',T121,'M/S'/)

        DO ROW=1,U_ROWS
          WRITE(6,535)ROW,(Z_V(ROW,LEVEL), LEVEL=1,P_LEVELS)
          IF (MOD(ROW,5).EQ.0) WRITE(6,'(1X)')
        END DO
535     FORMAT(T1,I3,T5,F5.2,19F6.2)
      ENDIF

      IF (IPRTKE.EQ.1) THEN
      WRITE(6,'(1X)')
      WRITE(6,'(A80/)') BANNER

      WRITE(6,540)
540   FORMAT(T2,'ROW   TKE1',T14,'TKE2',T20, 'TKE3',T26,'TKE4',T32,
     &'TKE5',T38,'TKE6', T44,'TKE7',T50,'TKE8',T56,'TKE9',T62,
     &'TKE10',T68,'TKE11',T74,'TKE12',T80, 'TKE13',T86,'TKE14',T92,
     &'TKE15',T98, 'TKE16',T104,'TKE17',T110,'TKE18',T116,
     &'TKE19',T122,'TKE20')

      WRITE(6,543)
543   FORMAT(T7,'M2/S2',T13,'M2/S2',T19,'M2/S2',T25,'M2/S2',T31,'M2/S2',
     & T37,'M2/S2',T43,'M2/S2',T49,'M2/S2',T55,'M2/S2',T61,'M2/S2',T67,
     & 'M2/S2',T73,'M2/S2',T79,'M2/S2',T85,'M2/S2',T91,'M2/S2',T97,
     & 'M2/S2',T103,'M2/S2',T109,'M2/S2',T115,'M2/S2',T121,'M2/S2'/)

      DO ROW=1,U_ROWS
        WRITE(6,545)  ROW,(Z_TKE_LEV(ROW,LEVEL)*0.5,LEVEL=1,P_LEVELS)
        IF(MOD(ROW,5).EQ.0) WRITE(6,'(1X)')
      END DO
545   FORMAT(T1,I3,T7,F5.1,19F6.1)

      WRITE(6,'(1X)')
      WRITE(6,'(A80/)') BANNER

      WRITE(6,550)
550   FORMAT(T2,'ROW   ZKE1',T14,'ZKE2',T20, 'ZKE3',T26,'ZKE4',T32,
     &'ZKE5',T38,'ZKE6', T44,'ZKE7',T50,'ZKE8',T56,'ZKE9',T62,
     &'ZKE10',T68,'ZKE11',T74,'ZKE12',T80, 'ZKE13',T86,'ZKE14',T92,
     &'ZKE15',T98, 'ZKE16',T104,'ZKE17',T110,'ZKE18',T116,
     &'ZKE19',T122,'ZKE20')

      WRITE(6,543)

      DO ROW=1,U_ROWS
        WRITE(6,545) ROW,(Z_ZKE_LEV(ROW,LEVEL)*0.5,LEVEL=1,P_LEVELS)
        IF(MOD(ROW,5).EQ.0) WRITE(6,'(1X)')
      END DO

      WRITE(6,'(1X)')
      WRITE(6,'(A80/)') BANNER

      WRITE(6,560)
560   FORMAT(T2,'ROW   EKE1',T14,'EKE2',T20, 'EKE3',T26,'EKE4',T32,
     &'EKE5',T38,'EKE6', T44,'EKE7',T50,'EKE8',T56,'EKE9',T62,
     &'EKE10',T68,'EKE11',T74,'EKE12',T80, 'EKE13',T86,'EKE14',T92,
     &'EKE15',T98, 'EKE16',T104,'EKE17',T110,'EKE18',T116,
     &'EKE19',T122,'EKE20')

      WRITE(6,543)

      DO ROW=1,U_ROWS
        WRITE(6,545) ROW,(Z_EKE_LEV(ROW,LEVEL)*0.5,LEVEL=1,P_LEVELS)
        IF(MOD(ROW,5).EQ.0) WRITE(6,'(1X)')
      END DO
      ENDIF

      IF(IPRTTEMP.EQ.1) THEN
        WRITE(6,'(1X)')
        WRITE(6,'(A80/)') BANNER
        WRITE(6,615)
615   FORMAT(T2,'ROW  T1',T13,'T2',T19,'T3',T25,'T4',T31,'T5',T37,'T6',
     &T43,'T7',T49,'T8',T55,'T9',T61,'T10',T67,'T11',T73,'T12',T79,
     &'T13',T85,'T14',T91,'T15',T97,'T16',T103,'T17',T109,'T18',T115,
     &'T19',T121,'T20')
        WRITE(6,625)
625   FORMAT(T6,'DEG.K',T12,'DEG.K',T18,'DEG.K',T24,'DEG.K',T30,'DEG.K',
     &T36,'DEG.K',T42,'DEG.K',T48,'DEG.K',T54,'DEG.K',T60,'DEG.K',T66,
     &'DEG.K',T72,'DEG.K',T78,'DEG.K',T84,'DEG.K',T90,'DEG.K',T96,
     &'DEG.K',T102,'DEG.K',T108,'DEG.K',T114,'DEG.K',T120,'DEG.K'/)
       DO ROW=1,P_ROWS
         WRITE(6,635)ROW,(Z_T(ROW,LEVEL), LEVEL=1,P_LEVELS)
         IF (MOD(ROW,5).EQ.0) WRITE(6,'(1X)')
       END DO
635    FORMAT(T1,I3,T5,F5.1,19F6.1)

       IF(IPRTVAR.EQ.1) THEN
         WRITE(6,'(1X)')
         WRITE(6,'(A80/)') BANNER
         WRITE(6,640)
640   FORMAT(T2,'ROW  TVR1',T13,'TVR2',T19, 'TVR3',T25,'TVR4',T31,
     &'TVR5',T37,'TVR6', T43,'TVR7',T49,'TVR8',T55,'TVR9',T61, 'TVR10',
     &T67,'TVR11',T73,'TVR12',T79,'TVR13',T85,'TVR14',T91, 'TVR15',T97,
     &  'TVR16',T103,'TVR17',T109,'TVR18',T115, 'TVR19',T121,'TVR20')

         WRITE(6,643)
643   FORMAT(T8,'K2',T14,'K2',T20,'K2',T26,'K2',T32,'K2', T38,'K2',T44,
     &'K2',T50,'K2',T56,'K2',T62,'K2',T68, 'K2',T74,'K2',T80,'K2', T86,
     &'K2',T92,'K2',T98, 'K2',T104,'K2',T110,'K2',T116,'K2',T122,'K2'/)

         DO ROW=1,P_ROWS
          WRITE(6,646) ROW,(Z_VAR_T_LEV(ROW,LEVEL),LEVEL=1,P_LEVELS)
          IF(MOD(ROW,5).EQ.0) WRITE(6,'(1X)')
         END DO
646      FORMAT(T1,I3,T7,F5.0,7F6.0,12F6.1)
       ENDIF
      ENDIF

      IF(IPRTQ.EQ.1) THEN
        WRITE(6,'(1X)')
        WRITE(6,'(A80/)') BANNER
        WRITE(6,715)
715   FORMAT(T2,'ROW Q1',T13,'Q2',T20,'Q3',T27,'Q4',T34,'Q5',T41,'Q6',
     &T47,'Q7',T53,'Q8',T59,'Q9',T65,'Q10',T71,'Q11',T77,'Q12',T83,
     &'Q13',T89,'Q14',T95,'Q15',T101,'Q16',T107,'Q17',T113,'Q18',T119,
     &'Q19',T125,'Q20')
        WRITE(6,725)
725   FORMAT(T5,'GM/KG',T12,'GM/KG',T19,'GM/KG',T26,'GM/KG',T33,'GM/KG',
     &T40,'GM/KG',T46,'GM/KG',T52,'GM/KG',T58,'GM/KG',T64,'GM/KG',T70,
     &'GM/KG',T76,'GM/KG',T82,'GM/KG',T88,'GM/KG',T94,'GM/KG',T100,
     &'GM/KG',T106,'GM/KG',T112,'GM/KG',T118,'GM/KG',T124,'GM/KG'/)
        DO ROW=1,P_ROWS
          WRITE(6,735)ROW,(Z_Q(ROW,LEVEL)*1000., LEVEL=1,Q_LEVELS)
          IF (MOD(ROW,5).EQ.0) WRITE(6,'(1X)')
        END DO
735     FORMAT(T1,I3,T5,F6.3,4F7.3,15F6.3)
        IF (IPRTVAR.EQ.1) THEN
         WRITE(6,'(1X)')
         WRITE(6,'(A80/)') BANNER

         WRITE(6,647)
647   FORMAT(T2,'ROW  QVR1',T13,'QVR2',T19, 'QVR3',T25,'QVR4',T31,
     &'QVR5',T37,'QVR6', T43,'QVR7',T49,'QVR8',T55,'QVR9',T60, 'QVR10',
     &T66,'QVR11',T72,'QVR12',T78,'QVR13',T84,'QVR14',T90, 'QVR15',T96,
     &  'QVR16',T102,'QVR17',T108,'QVR18',T114, 'QVR19',T120,'QVR20')

         WRITE(6,648)
 648  FORMAT(T6,'G2KG2',T12,'G2KG2',T18,'G2KG2',T24, 'G2KG2',T30,
     &'G2KG2', T36,'G2KG2',T42, 'G2KG2',T48,'G2KG2',T54,'G2KG2',T60,
     &'G2KG2', T66, 'G2KG2',T72,'G2KG2',T78,'G2KG2', T84, 'G2KG2',T90,
     &'G2KG2',T96, 'G2KG2',T102,'G2KG2',
     & T108,'G2KG2',T114,'G2KG2',T120,'G2KG2'/)

        DO ROW=1,P_ROWS
        WRITE(6,649)ROW,(Z_VAR_Q_LEV(ROW,LEVEL)*1.E6 ,LEVEL=1,Q_LEVELS)
         IF(MOD(ROW,5).EQ.0) WRITE(6,'(1X)')
        END DO
649     FORMAT(T1,I3,T6,F5.1,5F6.1,4F6.2,10F6.3)

       ENDIF
      ENDIF
CL
CL 9.2 ** Quarter Global Means Output **
CL
      WRITE(6,'(1X)')
      WRITE(6,'(A80/)') BANNER
      WRITE(6,15)
15    FORMAT(T40,'QUARTER GLOBAL MEAN STATISTICS'/)
      WRITE(6,116)
116   FORMAT(T3,'ROWS',T12,'PMSL',T20,'PSTAR',T28,'SNOWD',T37,'CANOPYW',
     &T47,'SH (W/m2)',T62,'LSRN',T68,'CVRN',
     &T76,'SOILM',T82,'T*ALL',T89,'T*SEA',T96,'T*LAND',T103,'SOILT1',
     &T110,'SOILT2',T117,'SOILT3',T124,'SOILT4')
      WRITE(6,126)
126   FORMAT(T4,'TO',T13,'MB',T21,'MB',T28,'MM',T38,'MM',T45,'All   Land
     &  Sea ',T62,'MM/',T68,'MM/',T76,'MM',T82,'DEG.K'
     &,T89,'DEG.K',T96,'DEG.K',T103,'DEG.K',T110,'DEG.K',T117,'DEG.K',
     & T124,'DEGK')
      WRITE(6,136)
136   FORMAT(T62,'DAY',T68,'DAY')
      DO 310, QUART=1,4
      WRITE(6,146)START_ROW(1,QUART),END_ROW(1,QUART),
     &                              (G_Q_PMSL(QUART)*0.01),
     &(G_Q_PSTAR(QUART)*0.01),G_Q_SNOWD(QUART),G_Q_CANOPYW(QUART),
     &G_Q_SH(QUART),G_Q_L_SH(QUART),G_Q_S_SH(QUART),G_Q_LSRN(QUART),
     &G_Q_CVRN(QUART),G_Q_SOILM(QUART),G_Q_TSTAR(QUART),
     &G_Q_S_TSTAR(QUART),G_Q_L_TSTAR(QUART),
     &(G_Q_SOILT(QUART,LEVEL),LEVEL=1,ST_LEVELS)
310   CONTINUE
146   FORMAT(T2,I3,':',T6,I3,T10,F7.2,T18,F7.2,T26,F9.3,T36,F7.3,T44,
     &F5.2,T50,F5.2,T56,F5.2,T62,F5.3,T68,F5.3,T76,F5.1,T82,
     &F6.2,T89,F6.2,T96,F6.2,T103,F6.2,T110,F6.2,T117,F6.2,T124,F6.2,
     &T131,F6.2)
      WRITE(6,117)(G_PMSL*0.01),(G_PSTAR*0.01),G_SNOWD,G_CANOPYW,G_SH,
     &G_L_SH,G_S_SH,G_LSRN,G_CVRN,G_SOILM,G_TSTAR,
     &G_S_TSTAR,G_L_TSTAR,(G_SOILT(LEVEL), LEVEL=1,ST_LEVELS)
117   FORMAT(T3,'GLOBAL ',T10,F7.2,T18,F7.2,T26,F9.3,T36,F7.3,T44,
     &F5.2,T50,F5.2,T56,F5.2,T62,F5.3,T68,F5.3,T76,F5.1,T82,
     &F6.2,T89,F6.2,T96,F6.2,T103,F6.2,T110,F6.2,T117,F6.2,T124,F6.2)

      WRITE(6,'(1X)')
      CTITLE   =' Rows     Land   Sea   Evaporation(mm/day)  Rainfall(mm
     &/day)  Snowfall(mm/day) Sublimation(mm/day) Soil Evap'
      CHEAD1   ='  to                   All   Land  Sea     All   Land
     &Sea   All   Land  Sea   All    Land  Sea    (mm/day)'
       WRITE(6,'(1X)')
       WRITE(6,'(A130)') CTITLE
       WRITE(6,'(A130)') CHEAD1
      DO 1310, QUART=1,4
      WRITE(6,1146)START_ROW(1,QUART),END_ROW(1,QUART),QGLPTS(QUART,1),
     &QGSPTS(QUART,1),(G_Q_EVAP(QUART)*86400.),
     &(G_Q_L_EVAP(QUART)*86400.),
     &(G_Q_S_EVAP(QUART)*86400.),G_Q_PPTN(QUART),G_Q_L_PPTN(QUART),
     &G_Q_S_PPTN(QUART),G_Q_SNOW(QUART),G_Q_L_SNOW(QUART),
     &G_Q_S_SNOW(QUART),G_Q_SUBL(QUART),G_Q_L_SUBL(QUART),
     &G_Q_S_SUBL(QUART),G_Q_L_SOEV(QUART)
1310   CONTINUE
1146   FORMAT(T2,I3,':',T6,I3,2I6,T22,4(1X,3(1X,F5.3)),2X,F5.3)
      WRITE(6,1117) GLPTS(1),GSPTS(1),(G_EVAP*86400.),(G_L_EVAP*86400.),
     &(G_S_EVAP*86400.),G_PPTN,G_L_PPTN,G_S_PPTN,G_SNOW,G_L_SNOW,
     &G_S_SNOW,G_SUBL,G_L_SUBL,G_S_SUBL,G_L_SOEV
1117   FORMAT(T2,'GLOBAL ',2I6,T22,4(1X,3(1X,F5.3)),2X,F5.3)

      IF(IPRTWIND.EQ.1) THEN
        WRITE(6,'(1X)')
        WRITE(6,216)
216   FORMAT(T3,'ROWS',T10,'TATMOS QATMOS',T38,'UATMOS VATMOS',T57,
     &'TAUX  (N/M2)',T78,'TAUY  (N/M2)',T95,'TKE',T102,'ZKE',T109,
     &'EKE')
        WRITE(6,226)
226   FORMAT(T10,'DEG.K',T17,'GM/KG',T40,'M/S',T47,'M/S',T53,'ALL',T60,
     &'SEA',T67,'LAND',T74,'ALL',T81,'SEA',T88,'LAND',T94,'M2/S2',
     &T101,'M2/S2',T108,'M2/S2'/)

        DO 320, QUART=1,4
      WRITE(6,236)START_ROW(1,QUART),END_ROW(1,QUART),G_Q_TATMOS(QUART),
     &(G_Q_QATMOS(QUART)*1000.),START_ROW(2,QUART),END_ROW(2,QUART),
     &G_Q_UATMOS(QUART),G_Q_VATMOS(QUART),G_Q_TAUX(QUART),
     &G_Q_S_TAUX(QUART),G_Q_L_TAUX(QUART),G_Q_TAUY(QUART),
     &G_Q_S_TAUY(QUART),G_Q_L_TAUY(QUART),
     &G_Q_TKE(QUART)*0.5,G_Q_ZKE(QUART)*0.5,G_Q_EKE(QUART)*0.5
320     CONTINUE

236   FORMAT(T2,I3,':',T6,I3,T10,F6.2,F6.2,T30,I3,':',T34,I3,T38,F6.2,
     &F6.2,6F7.3,3F7.2)
      WRITE(6,127)G_TATMOS,(G_QATMOS*1000.),G_UATMOS,G_VATMOS,G_TAUX,
     &          G_S_TAUX,G_L_TAUX,G_TAUY,G_S_TAUY,G_L_TAUY,
     &          G_TKE*0.5,G_ZKE*0.5,G_EKE*0.5
127   FORMAT(T3,'GLOBAL ',T10,F6.2,F6.2,T31,'GLOBAL',T38,F6.2,F6.2,
     &6F7.3,3F7.2)
        WRITE(6,'(1X)')
        WRITE (6,316)
316   FORMAT(T3,'ROW',T11,'U1',T17,'U2',T23,'U3',T29,'U4',T35,'U5',
     &T41,'U6',T47,'U7',T53,'U8',T59,'U9',T65,'U10',T71,'U11',T77,'U12',
     &T83,'U13',T89,'U14',T95,'U15',T101,'U16',T107,'U17',T113,'U18',
     &T119,'U19',T125,'U20')
        WRITE(6,326)
326   FORMAT(T11,'M/S',T17,'M/S',T23,'M/S',T29,'M/S',T35,'M/S',
     &T41,'M/S',T47,'M/S',T53,'M/S',T59,'M/S',T65,'M/S',T71,
     &'M/S',T77,'M/S',T83,'M/S',T89,'M/S',T95,'M/S',T101,
     &'M/S',T107,'M/S',T113,'M/S',T119,'M/S',T125,'M/S'/)
        DO QUART=1,4
          WRITE(6,336)START_ROW(2,QUART),END_ROW(2,QUART),
     &                           (G_Q_U(QUART,LEVEL), LEVEL=1,P_LEVELS)
        END DO
336     FORMAT(T2,I3,':',T6,I3,T10,F5.2,19F6.2)
        WRITE(6,137)(G_U(LEVEL), LEVEL=1,P_LEVELS)
137     FORMAT(T3,'GLOBAL ',T10,F5.2,19F6.2)
        WRITE(6,'(1X)')
        WRITE (6,416)
416   FORMAT(T3,'ROW',T11,'V1',T17,'V2',T23,'V3',T29,'V4',T35,'V5',
     &T41,'V6',T47,'V7',T53,'V8',T59,'V9',T65,'V10',T71,'V11',T77,'V12',
     &T83,'V13',T89,'V14',T95,'V15',T101,'V16',T107,'V17',T113,'V18',
     &T119,'V19',T125,'V20')
        WRITE(6,426)
426   FORMAT(T11,'M/S',T17,'M/S',T23,'M/S',T29,'M/S',T35,'M/S',
     &T41,'M/S',T47,'M/S',T53,'M/S',T59,'M/S',T65,'M/S',T71,
     &'M/S',T77,'M/S',T83,'M/S',T89,'M/S',T95,'M/S',T101,
     &'M/S',T107,'M/S',T113,'M/S',T119,'M/S',T125,'M/S'/)
        DO QUART=1,4
         WRITE(6,436)START_ROW(2,QUART),END_ROW(2,QUART),
     &                           (G_Q_V(QUART,LEVEL), LEVEL=1,P_LEVELS)
        END DO
436     FORMAT(T2,I3,':',T6,I3,T10,F5.2,19F6.2)
        WRITE(6,147)(G_V(LEVEL), LEVEL=1,P_LEVELS)
147     FORMAT(T3,'GLOBAL ',T10,F5.2,19F6.2)
        WRITE(6,1140)
      ENDIF
      IF(IPRTKE.EQ.1) THEN
 1140   FORMAT(T3,'KINETIC ENERGY IN M2/S2')
        WRITE(6,1141)(G_TKE_LEV(LEVEL)*0.5,LEVEL=1,P_LEVELS)
 1141   FORMAT(T3,'TKE',T10,F5.1,19F6.1)
        WRITE(6,1142)(G_ZKE_LEV(LEVEL)*0.5,LEVEL=1,P_LEVELS)
 1142   FORMAT(T3,'ZKE',T10,F5.1,19F6.1)
        WRITE(6,1143)(G_EKE_LEV(LEVEL)*0.5,LEVEL=1,P_LEVELS)
 1143   FORMAT(T3,'EKE',T10,F5.1,19F6.1)
      END IF

      IF(IPRTTEMP.EQ.1) THEN
        WRITE(6,'(1X)')
        WRITE (6,516)
516     FORMAT(T3,'ROWS',T10,' T1     T2     T3     T4     T5     T6'
     &,T53,'T7     T8     T9     T10')
        WRITE(6,526)
526   FORMAT(T10,'DEG.K  DEG.K  DEG.K  DEG.K  DEG.K  DEG.K  DEG.K',
     &T59,'DEG.K  DEG.K  DEG.K'/)
        DO QUART=1,4
          WRITE(6,536)START_ROW(1,QUART),END_ROW(1,QUART),
     &                                 (G_Q_T(QUART,LEVEL), LEVEL=1,10)
        END DO
536     FORMAT(T2,I3,':',T6,I3,T9,10F7.2)
        WRITE(6,157)(G_T(LEVEL), LEVEL=1,10)
157     FORMAT(T3,'GLOBAL ',T9,10F7.2)
        WRITE(6,1150)
 1150   FORMAT(T2,'TEMPERATURE VARIANCE IN DEG. K2')
        WRITE(6,1151) (G_VAR_T_LEV(LEVEL),LEVEL=1,10)
 1151   FORMAT(T9,10F7.2)
        WRITE(6,'(1X)')
        WRITE (6,546)
546   FORMAT(T3,'ROWS',T10,' T11    T12    T13    T14    T15    T16'
     &,T53,'T17    T18    T19    T20')
        WRITE(6,556)
556   FORMAT(T10,'DEG.K  DEG.K  DEG.K  DEG.K  DEG.K  DEG.K  DEG.K',
     &T59,'DEG.K  DEG.K  DEG.K'/)
        DO QUART=1,4
          WRITE(6,566)START_ROW(1,QUART),END_ROW(1,QUART),
     &                          (G_Q_T(QUART,LEVEL), LEVEL=11,P_LEVELS)
        END DO

566     FORMAT(T2,I3,':',T6,I3,T9,10F7.2)
        WRITE(6,167)(G_T(LEVEL), LEVEL=11,P_LEVELS)
167     FORMAT(T3,'GLOBAL ',T9,10F7.2)
        WRITE(6,1160)
 1160   FORMAT(T3,'TEMPERATURE VARIANCE IN DEG. K2')
        WRITE(6,1161) (G_VAR_T_LEV(LEVEL),LEVEL=11,P_LEVELS)
 1161   FORMAT(T9,10F7.2)
      ENDIF

      IF(IPRTQ.EQ.1) THEN
        WRITE(6,'(1X)')
        WRITE (6,616)
616   FORMAT(T3,'ROW',T11,'Q1',T18,'Q2',T25,'Q3',T32,'Q4',T39,'Q5',
     &T45,'Q6',T51,'Q7',T57,'Q8',T63,'Q9',T69,'Q10',T75,'Q11',T81,'Q12',
     &T87,'Q13',T93,'Q14',T99,'Q15',T105,'Q16',T111,'Q17',T117,'Q18',
     &T123,'Q19',T129,'Q20')
        WRITE(6,626)
626   FORMAT(T10,'GM/KG',T17,'GM/KG',T24,'GM/KG',T31,'GM/KG',T38,'GM/KG'
     &,T44,'GM/KG',T50,'GM/KG',T56,'GM/KG',T62,'GM/KG',T68,'GM/KG',T74,
     &'GM/KG',T80,'GM/KG',T86,'GM/KG',T92,'GM/KG',T98,'GM/KG',T104,
     &'GM/KG',T110,'GM/KG',T116,'GM/KG',T122,'GM/KG',T128,'GM/KG'/)

        DO QUART=1,4
          WRITE(6,636)START_ROW(1,QUART),END_ROW(1,QUART),
     &                   (G_Q_Q(QUART,LEVEL)*1000., LEVEL=1,Q_LEVELS)
        END DO
636     FORMAT(T1,I3,':',T5,I3,T8,5F7.3,15F6.3)
        WRITE(6,177)(G_Q(LEVEL)*1000., LEVEL=1,Q_LEVELS)
177     FORMAT(T2,'GLOBAL ',T8,5F7.3,15F6.3)
        WRITE(6,1170)
 1170   FORMAT(T3,'MOISTURE VARIANCE IN G2/KG2')
        WRITE(6,1171) (G_VAR_Q_LEV(LEVEL)*1.E6,LEVEL=1,Q_LEVELS)
 1171   FORMAT(T8,5F7.3,15F6.3)
      ENDIF

      WRITE(6,'(1X)')
      WRITE(6,1180)
 1180 FORMAT(T3,'GLOBAL MEANS')
      WRITE(6,1181) G_TKE*0.5, G_ZKE*0.5, G_EKE*0.5
 1181 FORMAT(T5,'TKE',T10,F8.1,T19,'ZKE',T24,F8.1,T33,'EKE',F8.1,
     &      '(M2 S-2)')
      WRITE(6,1182) G_VAR_T,G_VAR_Q*1.E6
 1182 FORMAT(T5,'T_VARIANCE',T15,F8.1,T23,' K2  Q_VARIANCE',T38,F8.1,
     &         T47,' G2/KG2')
CL----------------------------------------------------------------------
C     Extra variables added and only printed if available
C
      IF (IPRTEXTRA.EQ.1) THEN
       IF(IHYDRO.GT.1) THEN
      CTITLE   ='Quarter-globe and Global means over land  from hydrolog
     &y section'
      CHEAD1   ='          Thro Fall  Snowmelt  Sur runoff  Subsur Runof
     &f'
      CHEAD2   ='           mm/day     mm/day     mm/day      mm/day'
        WRITE(6,'(1X)')
        WRITE(6,'(A130)') CTITLE
        WRITE(6,'(A130)') CHEAD1
        WRITE(6,'(A130)') CHEAD2

        DO QUART=1,4
         WRITE(6,2004)START_ROW(1,QUART),END_ROW(1,QUART),
     &      G_Q_L_THRF(QUART),G_Q_L_SNML(QUART),G_Q_L_SFRU(QUART),
     &      G_Q_L_SBRU(QUART)
        END DO
 2004   FORMAT(T1,I3,':',T5,I3,T8,4(3x,F8.3))
        WRITE(6,2006) G_L_THRF,G_L_SNML,G_L_SFRU,G_L_SBRU
 2006   FORMAT(T2,'GLOBAL ',T8,4(3x,F8.3))
       END IF

       IF(IRAD.GT.1) THEN
      CTITLE   ='   Quarter-globe and Global means - Radiation in W/m2'
      CHEAD1   ='          Total Surface        Surface Solar      Total
     & TOA  LW TOA   SW TOA   Albedo'
      CHEAD2   ='         All    Land   Sea    All    land   Sea     rad
     &       rad      rad'
        WRITE(6,'(1X)')
        WRITE(6,'(A130)') CTITLE
        WRITE(6,'(A130)') CHEAD1
        WRITE(6,'(A130)') CHEAD2
        DO I=1,4
          WRITE(6,2014)START_ROW(1,I),END_ROW(1,I),
     &      G_Q_SDTR(I),G_Q_L_SDTR(I),G_Q_S_SDTR(I),
     &      G_Q_SDSR(I),G_Q_L_SDSR(I),G_Q_S_SDSR(I),
     &      G_Q_TDTR(I),G_Q_TOLR(I),G_Q_TOSW(I),G_Q_ALBEDO(I)
        END DO
 2014   FORMAT(T1,I3,':',T5,I3,T8,6(1X,F6.1),T52,3(F7.2,2X),F7.4)
        WRITE(6,2016) G_SDTR,G_L_SDTR,G_S_SDTR,G_SDSR,G_L_SDSR,
     &      G_S_SDSR,G_TDTR,G_TOLR,G_TOSW,G_ALBEDO
 2016   FORMAT(T2,'GLOBAL ',T8,6(1X,F6.1),T52,3(F7.2,2X),F7.4)
       END IF

       IF (IFLUXL.GT.1) THEN
      CTITLE   ='     Quarter-globe and Global means - Fluxes over land
     &and over whole atmosphere'
      CHEAD1   ='        Soil Flux  Snow Bud  Can Bud  Soil Bud  Water(l
     &+s) Soil En  Energy Flux (l+s)'
      CHEAD2   ='         W/m2      mm/day    mm/day    mm/day    mm/day
     &     W/m2     W/m2'
        WRITE(6,'(1X)')
        WRITE(6,'(A130)') CTITLE
        WRITE(6,'(A130)') CHEAD1
        WRITE(6,'(A130)') CHEAD2
        DO I=1,4
         WRITE(6,2024)START_ROW(1,I),END_ROW(1,I),
     &      G_Q_L_SOHF(I),G_Q_L_WFSS(I),G_Q_L_WFCA(I),G_Q_L_WFSO(I),
     &      G_Q_WAFL(I),G_Q_L_ENFS(I),G_Q_ENFL(I)
        END DO
 2024   FORMAT(T1,I3,':',T5,I3,T8,7(F8.3,2x))
        WRITE(6,2026) G_L_SOHF,G_L_WFSS,G_L_WFCA,G_L_WFSO,G_WAFL,
     &      G_L_ENFS,G_ENFL
 2026   FORMAT(T2,'GLOBAL ',T8,7(F8.3,2x))
       END IF
C
C  Cloud water printouts
C
       IF(ICLOUD.GT.1) THEN
       CHEAD1  ='ROW       1     2     3     4     5     6     7     8
     &   9    10    11    12    13    14    15    16    17    18    19
     &  20'
       CHEAD2  =' to '
       CTITLE  ='         Cloud liquid water - all levels (g/kg x1000)'
         WRITE(6,'(1X)')
         WRITE(6,'(A130)') CTITLE
         WRITE(6,'(A130)') CHEAD1
         WRITE(6,'(A130)') CHEAD2
         DO I=1,4
        WRITE(6,3024)START_ROW(1,I),END_ROW(1,I),(1000000.*G_Q_CLLIQL(I,
     &LEVEL),LEVEL=1,Q_LEVELS)
         END DO
         WRITE(6,3026) (1000000.*G_CLLIQL(LEVEL),LEVEL=1,Q_LEVELS)
      CTITLE   ='         Cloud ice water    - all levels (g/kg x1000)'
         WRITE(6,'(1X)')
         WRITE(6,'(A130)') CTITLE
         WRITE(6,'(A130)') CHEAD1
         WRITE(6,'(A130)') CHEAD2
         DO I=1,4
        WRITE(6,3024)START_ROW(1,I),END_ROW(1,I),(1000000.*G_Q_CLICEL(I,
     &LEVEL),LEVEL=1,Q_LEVELS)
         END DO
         WRITE(6,3026) (1000000.*G_CLICEL(LEVEL),LEVEL=1,Q_LEVELS)
         CTITLE  ='  Cloud water whole atmosphere (g/kg x1000)'
         CHEAD1  =' ROW      Liquid    Ice'
         CHEAD2  ='  to     '
         WRITE(6,'(1X)')
         WRITE(6,'(A130)') CTITLE
         WRITE(6,'(A130)') CHEAD1
         WRITE(6,'(A130)') CHEAD2
         DO I=1,4
          WRITE(6,3030)START_ROW(1,I),END_ROW(1,I),
     &    (1000000.*G_Q_CLLIQ(I)),(1000000.*G_Q_CLICE(I))
         END DO
         WRITE(6,3032) (1000000.*G_CLLIQ),(1000000.*G_CLICE)
 3024    FORMAT(T1,I3,':',T5,I3,T8,20(F6.2))
 3026    FORMAT(T2,'GLOBAL ',T8,20(F6.2))
 3030    FORMAT(T1,I3,':',T5,I3,T8,2(3X,F6.2))
 3032    FORMAT(T2,'GLOBAL ',T8,2(3X,F6.2))
       ENDIF
       IF(ISICE.GT.1) THEN
       CTITLE  ='    Quarter-globe means  Sea-ice meaned over all sea po
     &ints'
      CHEAD1   ='        Heat Flux  Melt ht   S Heat  Sur Rad  Sur SW
     &Energy Fl'
      CHEAD2   ='         W/m2      W/m2      W/m2     W/m2     W/m2
     & W/m2'
        WRITE(6,'(1X)')
        WRITE(6,'(A130)') CTITLE
        WRITE(6,'(A130)') CHEAD1
        WRITE(6,'(A130)') CHEAD2
        DO I=1,4
         WRITE(6,3124)START_ROW(1,I),END_ROW(1,I),G_Q_SIHF(I),
     &   G_Q_SIMH(I),G_Q_SISH(I),G_Q_SIST(I),G_Q_SISS(I),G_Q_SIEF(I)
        END DO
         WRITE(6,3126) G_SIHF,G_SIMH,G_SISH,G_SIST,
     &     G_SISS,G_SIEF
 3124    FORMAT(T1,I3,':',T5,I3,T8,6(1X,F8.3))
 3126    FORMAT(T2,'GLOBAL ',T8,6(1X,F8.3))
       ENDIF
      ENDIF
       WRITE(6,'(1X)')
CL----------------------------------------------------------------------
 999  CONTINUE
      RETURN
CL----------------------------------------------------------------------
      END
