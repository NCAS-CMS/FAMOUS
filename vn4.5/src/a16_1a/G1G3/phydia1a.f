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
CLL  SUBROUTINE PHY_DIAG------------------------------------------------
CLL
CLL  PURPOSE:   Calculation of physical diagnostics. These diagnostics
CLL             are those based on temperatures humidities and pressure.
CLL             Those diagnostics involving the winds are in the other
CLL             diagnostic routine DYN_DIAG
CLL             Date: 08/06/92
CLL     Exner at model level made consistent with the geopotential eqn.
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  Date
CLL  3.1     10/01/93 Get model half Hts for all dependant routines.
CLL  3.1     09/02/93
CLL  Modified by C.Wilson
CLL  Add alternate call to V_INT_T for temperatures on pressure levels
CLL  to call V_INT_TP if QV_INT_TP is true
CLL  so as to reduce bias in stratospheric temperature output
CLL   3.2    13/07/93 Changed CHARACTER*(*) to CHARACTER*(80) for
CLL                   portability.  Author Tracey Smith.
CLL  3.2  07/06/93  Correct indexing of WORK1 for use in DIA_THADV
CLL   3.4    21/9/94  Calculate model level geopotential heights
CLL                   (in m).       Author S.A.Woltering.
CLL  3.4  15/07/94  Redimension WORK1 from P_TO_UV for input to
CLL                 DIA_THADV using new utility routine.  Rick Rawlins.
CLL  4.0  22/11/95  Correction to Model level geopotential heights.
CLL                                 Author Steve Woltering
CLL  4.1  15/05/96  Put tracer fields on pressure levels.   D.Podd
CLL  4.4  16/10/97  Extra swapbounds for P_EXNER_HALF. D Robinson.
!LL  4.5  15/04/98  Added start-end arguments to V_INT, V_INT_T and
!LL                 V_INT_Z routines. S.D.Mullerworth
CLL  4.5   3/09/98  Corrected relative humidity diagnostic. 
CLL                                                    D Wilson.
CLL  4.5  05/06/98 New argument L_VINT_TP. Use to call V_INT_TP or
CLL                V_INT_T. Also L_LSPICE for correcting relative
CLL                humidity diagnostic. D. Robinson.
CLL
CLL
CLL  Programming standard: U M DOC  Paper NO. 4,
CLL
CLL  SYSTEM COMPONENTS D471 D421 D481 D482 D432 D423 D422 D431 D472
CLL                    D435 D484 D441
CLL
CLL  SYSTEM TASK: D4
CLL
CLL  External documentation
CLL
CLLEND------------------------------------------------------------------
C
C*L  ARGUMENTS:---------------------------------------------------------
      SUBROUTINE PHY_DIAG(
     &  FIRST_ROW , TOP_ROW_START , P_LAST_ROW , U_LAST_ROW,
     &  P_BOT_ROW_START , U_BOT_ROW_START , upd_P_ROWS , upd_U_ROWS,
     &  FIRST_FLD_PT , LAST_P_FLD_PT , LAST_U_FLD_PT,
     &  FIRST_VALID_PT , LAST_P_VALID_PT , LAST_U_VALID_PT,
     &  VALID_P_ROWS, VALID_U_ROWS,
     &  START_POINT_NO_HALO, START_POINT_INC_HALO,
     &  END_P_POINT_NO_HALO, END_P_POINT_INC_HALO,
     &  END_U_POINT_NO_HALO, END_U_POINT_INC_HALO,
     &  FIRST_ROW_PT ,  LAST_ROW_PT , tot_P_ROWS , tot_U_ROWS,
     &  GLOBAL_ROW_LENGTH, GLOBAL_P_FIELD, GLOBAL_U_FIELD,
     &  MY_PROC_ID , NP_PROC_ID , SP_PROC_ID ,
     &  GC_ALL_GROUP, GC_ROW_GROUP, GC_COL_GROUP, N_PROCS,
     &  EW_Halo , NS_Halo, halo_4th, extra_EW_Halo, extra_NS_Halo,
     &  LOCAL_ROW_LENGTH, FIRST_GLOBAL_ROW_NUMBER,
     &  at_top_of_LPG,at_right_of_LPG, at_base_of_LPG,at_left_of_LPG,
C   primary data in
     &  PSTAR,U,V,Q,THETA,OROG,P_EXNER_HALF,LAND,TSTAR,TRACERS,
C   primary data constants
     &  U_ROWS,P_ROWS,ROW_LENGTH,P_LEVELS,Q_LEVELS,P_FIELD,
     &  U_FIELD,AK,BK,AKH,BKH,EW_SPACE,NS_SPACE,SEC_U_LATITUDE,
     &  TR_LEVELS,TR_VARS,TR_P_FIELD_DA,
C   Input pressure for required interpolation
     &  T_P_PRESS,HTS_PRESS,REL_HUMID_PRESS,WBPT_PRESS,TH_ADV_PRESS,
     &  TR_PRESS,
C   Input pressure for indices for product field
     &  H2_IND,
C   DIAGNOSTICS OUT
     &  HEIGHTS,  T_P, REL_HUMID_P, WBPT,     SNPROB, MINUS20_ICAO,
     &  MINUS20_P,  FREEZE_ICAO, FREEZE_Z, FREEZE_P, CONTRAIL_LOWERHT,
     &  CONTRAIL_UPPERHT, TROP_P,   TROP_T,   TROP_Z, TROP_ICAO,
     &  MINUS20_Z,  TH_ADV_SINGLE,     DUCT_HEIGHT,     DUCT_INT,
     &  P_MSL,      TH_ADV_MEAN,  HEIGHTS2,  MODHT,
     &  TRACERWORK,INT16,PT_TRACER,
C   diagnostic lengths
     &  T_P_LEVS,HTS_LEVS,REL_HUMID_P_LEVS,WBPT_LEVS,TH_ADV_P_LEVS,
     &  H2_P_LEVS,
     &  TR_PRESS_LEVS,NUM_STASH_LEVELSDA,
C   diagnostic logical indicators
     &  QMODEL_HALF_HEIGHT,QHEIGHTS,    QT_P,          QREL_HUMID_P,
     &  QWBPT,             QSNPROB,     QMINUS20_ICAO, QMINUS20_P,
     &  QFREEZE_ICAO,      QFREEZE_Z,   QFREEZE_P,  QCONTRA_LOWERHT,
     &  QCONTRA_UPPERHT,   QTROP_P,     QTROP_T,       QTROP_Z,
     &  QTROP_ICAO,        QMINUS20_Z,  QTH_ADV_SINGLE,   QDUCT_HEIGHT,
     &  QDUCT_INT,         QPMSL,       QTH_ADV_MEAN, QHTS2,  QMODHT,
     &  SF_TRACER,
     &  Z_REF, L_VINT_TP, L_LSPICE,
     &  ICODE,CMESSAGE)
C
      IMPLICIT NONE
C*L------------------COMDECK C_EPSLON-----------------------------------
C EPSILON IS RATIO OF MOLECULAR WEIGHTS OF WATER AND DRY AIR
      REAL EPSILON,C_VIRTUAL

      PARAMETER(EPSILON=0.62198,
     &          C_VIRTUAL=1./EPSILON-1.)
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
! Comdeck TYPFLDPT
! Variables which point to useful positions in a horizontal field

      INTEGER
     &  FIRST_ROW        ! First updatable row on field
     &, TOP_ROW_START    ! First point of north-pole (global) or
!                        ! Northern (LAM) row
!                        ! for processors not at top of LPG, this
!                        ! is the first point of valid data
!                        ! (ie. Northern halo).
     &, P_LAST_ROW       ! Last updatable row on pressure point field
     &, U_LAST_ROW       ! Last updatable row on wind point field
     &, P_BOT_ROW_START  ! First point of south-pole (global) or
!                        ! Southern (LAM) row on press-point field
     &, U_BOT_ROW_START  ! First point of south-pole (global) or
!                        ! Southern (LAM) row on wind-point field
!                        ! for processors not at base of LPG, this
!                        ! is the start of the last row of valid data
!                        ! (ie. Southern halo).
     &, upd_P_ROWS       ! number of P_ROWS to be updated
     &, upd_U_ROWS       ! number of U_ROWS to be updated
     &, FIRST_FLD_PT     ! First point on field
     &, LAST_P_FLD_PT    ! Last point on pressure point field
     &, LAST_U_FLD_PT    ! Last point on wind point field
! For the last three variables, these indexes are the start points
! and end points of "local" data - ie. missing the top and bottom
! halo regions.
     &, FIRST_VALID_PT   ! first valid point of data on field
     &, LAST_P_VALID_PT  ! last valid point of data on field
     &, LAST_U_VALID_PT  ! last valid point of data on field
     &, VALID_P_ROWS     ! number of valid rows of P data
     &, VALID_U_ROWS     ! number of valid rows of U data
     &, START_POINT_NO_HALO
!                        ! first non-polar point of field (misses
!                        ! halo for MPP code)
     &, START_POINT_INC_HALO
!                        ! first non-polar point of field (includes
!                        ! halo for MPP code)
     &, END_P_POINT_NO_HALO
!                        ! last non-polar point of P field (misses
!                        ! halo for MPP code)
     &, END_P_POINT_INC_HALO
!                        ! last non-polar point of P field (includes
!                        ! halo for MPP code)
     &, END_U_POINT_NO_HALO
!                        ! last non-polar point of U field (misses
!                        ! halo for MPP code)
     &, END_U_POINT_INC_HALO
!                        ! last non-polar point of U field (includes
!                        ! halo for MPP code)
     &, FIRST_ROW_PT     ! first data point along a row
     &, LAST_ROW_PT      ! last data point along a row
! For the last two variables, these indexes are the start and
! end points along a row of the "local" data - ie. missing out
! the east and west halos
     &, tot_P_ROWS         ! total number of P_ROWS on grid
     &, tot_U_ROWS         ! total number of U_ROWS on grid
     &, GLOBAL_ROW_LENGTH  ! length of a global row
     &, GLOBAL_P_FIELD     ! size of a global P field
     &, GLOBAL_U_FIELD     ! size of a global U field
!

     &, MY_PROC_ID         ! my processor id
     &, NP_PROC_ID         ! processor number of North Pole Processor
     &, SP_PROC_ID         ! processor number of South Pole Processor
     &, GC_ALL_GROUP       ! group id of group of all processors
     &, GC_ROW_GROUP       ! group id of group of all processors on this
!                          ! processor row
     &, GC_COL_GROUP       ! group id of group of all processors on this
!                          ! processor column
     &, N_PROCS            ! total number of processors

     &, EW_Halo            ! Halo size in the EW direction
     &, NS_Halo            ! Halo size in the NS direction

     &, halo_4th           ! halo size for 4th order calculations
     &, extra_EW_Halo      ! extra halo size required for 4th order
     &, extra_NS_Halo      ! extra halo size required for 4th order
     &, LOCAL_ROW_LENGTH   ! size of local row
     &, FIRST_GLOBAL_ROW_NUMBER
!                          ! First row number on Global Grid    

! Variables which indicate if special operations are required at the
! edges.
      LOGICAL
     &  at_top_of_LPG    ! Logical variables indicating if this
     &, at_right_of_LPG  ! processor is at the edge of the Logical
     &, at_base_of_LPG   ! Processor Grid and should process its edge
     &, at_left_of_LPG   ! data differently.

! End of comdeck TYPFLDPT

      INTEGER
     *  P_FIELD            !IN    1ST DIMENSION OF FIELD OF PSTAR
     *, U_FIELD            !IN    1ST DIMENSION OF FIELD OF U,V
     *, U_ROWS             !IN    NUMBER OF ROWS FOR U,V FIELD
     *, P_ROWS             !IN    NUMBER OF ROWS FOR P,T FIELD
     *, ROW_LENGTH         !IN    NUMBER OF POINTS PER ROW
     *, LEVELS             !IN    NUMBER OF MODEL LEVELS
     *, P_LEVELS           !IN    NUMBER OF PRESSURE LEVELS
     *, Q_LEVELS           !IN    NUMBER OF WET LEVELS
     *, ICODE              !IN    RETURN CODE   :  IRET=0   NORMAL EXIT
     *, NUM_STASH_LEVELSDA !IN    NUMBER OF STASH LEVELS
     *, TR_LEVELS          !IN    NUMBER OF TRACER LEVELS
     *, TR_VARS            !IN    NUMBER OF TRACERS
     *, TR_P_FIELD_DA      !IN    P_FIELD for DA of tracer arrays
     *, INT16              !IN    Dummy variable for STASH_MAXLEN(16)
     *, PT_TRACER(TR_VARS+1)
     *                     !IN    POINTERS TO TRACERS IN STASHWORK
     *, TR_PRESS_LEVS(TR_VARS+1)
     *                     !IN    NO OF LEVS ON WHICH TO INTERP TRACERS
     *, T_P_LEVS           !IN    NO OF LEVS ON WHICH TO INTERP T_P
     *, REL_HUMID_P_LEVS   !IN    NO OF LEVS ON WHICH TO INTERP REL HUM
     *, HTS_LEVS           !IN    NO OF LEVS ON WHICH TO INTERP HEIGHTS
     *, H2_P_LEVS          !IN    NO OF LEVS ON WHICH TO INTERP H**2
     *, WBPT_LEVS          !IN    NO OF LEVS ON WHICH TO CALCULATE WBPT
     *, TH_ADV_P_LEVS      !IN    NO OF LEVS ON WHICH TO CALCULATE THADV
     *, Z_REF              !IN    LEVEL OF MODEL USED TO CALCULATE PMSL
     *, BOTTOM_QC_LEVEL    ! Tropopause quality control level
     *, TOP_QC_LEVEL       !      "        "      "       "
      INTEGER
     &   H2_IND(H2_P_LEVS) !IN  indices for pressure levels for h**2
C
C
      LOGICAL
     * QT_P           !IN  LOGICAL FLAG FOR   PRESS INTER TEMPERTATURE
     *,QHEIGHTS       !IN     "     "    "    HEIGHTS ON ANY P SURFACE
     *,QREL_HUMID_P   !IN     "     "    "    HUMIDITIES ANY P SURFACE
     *,QWBPT          !IN     "     "    "    WET BULB POTENTIAL TEMP
     *,QSNPROB        !IN     "     "    "    SNOW PROBABILITY
     *,QMINUS20_ICAO  !IN     "     "    "    -20 deg C LEVEL ICAO HT
     *,QMINUS20_P     !IN     "     "    "    -20 deg C LEVEL PRESSURE
     *,QFREEZE_ICAO   !IN     "     "    "    FREEZING LEVEL ICAO HT
     *,QFREEZE_Z      !IN     "     "    "    FREEZING LEVEL HT
     *,QFREEZE_P      !IN     "     "    "    FREEZING LEVEL PRESSURE
     *,QCONTRA_LOWERHT!IN     "     "    "    HEIGHT OF LOWER CONTRAIL
     *,QCONTRA_UPPERHT!IN     "     "    "    HEIGHT OF UPPER CONTRAIL
     *,QMODEL_HALF_HEIGHT!IN  "     "    "    HEIGHTS ON MODEL LAYER
     *,LAND(P_FIELD)  !IN    LAND SEA MASK SET TRUE FOR LAND
      LOGICAL
     * QTROP_P        !IN     "     "    "    TROPOPAUSE PRESSURE
     *,QTROP_T        !IN     "     "    "    TROPOPAUSE TEMPERATURE
     *,QTROP_Z        !IN     "     "    "    TROPOPAUSE HEIGHT
     *,QTROP_ICAO     !IN     "     "    "    TROPOPAUSE ICAO HEIGHT
     *,QMINUS20_Z     !IN     "     "    "    -20 deg C LEVEL HEIGHT
     *,QTH_ADV_SINGLE !IN     "     "    "    THERMAL ADV SINGLE LEV
     *,QDUCT_HEIGHT   !IN     "     "    "    RADIO DUCT HEIGHT
     *,QDUCT_INT      !IN     "     "    "    RADIO DUCT INTENSITY
     *,QPMSL          !IN     "     "    "    PRESSURE AT MEAN SEA LVL
     *,QTH_ADV_MEAN   !IN     "     "    "    THERMAL ADVECTION MEAN
     *,QHTS2          !IN     "     "    "    heights**2 on pressure lev
     *,QMODHT         !IN     "     "    "    Model level geopot heights
     *,SF_TRACER(TR_VARS+1) !IN   LOGICAL FLAGS FOR TRACERS
     &, L_VINT_TP     !IN  Switch to control Vert Int for Output of
                      !    Temperature on model levels. D Robinson.
     &, L_LSPICE      !IN  Switch for New cloud/precip microphysics.
C
      CHARACTER*80 CMESSAGE
C
      REAL
     * PSTAR(P_FIELD)         !IN    PRIMARY MODEL ARRAY FOR PSTAR FIELD
     *,TSTAR(P_FIELD)         !IN    PRIMARY MODEL ARRAY FOR TSTAR FIELD
     *,OROG(P_FIELD)         !IN    PRIMARY MODEL OROGRAPHY
     *,P_EXNER_HALF(P_FIELD,P_LEVELS+1) !IN  EXNER PRESS ON 1/2 LVLS
     *,THETA(P_FIELD,P_LEVELS) !IN PRIMARY MODEL ARRAY FOR THETA FIELD
     *,U(U_FIELD,P_LEVELS)    !INT PRIMARY MODEL ARRAY FOR U FIELD
     *,V(U_FIELD,P_LEVELS)    !IN PRIMARY MODEL ARRAY FOR V FIELD
     *,Q(P_FIELD,Q_LEVELS)    !IN PRIMARY MODEL ARRAY FOR HUMIDITY
     *,TRACERS(TR_P_FIELD_DA,TR_LEVELS,TR_VARS+1)
     *                        !IN PRIMARY MODEL ARRAY FOR TRACERS
     *,HEIGHTS(P_FIELD,HTS_LEVS) ! OUTPUT HEIGHTS ON ANY P SURFACE
     *,T_P(P_FIELD,T_P_LEVS) ! OUTPUT TEMPS ON ANY P SURFACE
     *,P_MSL(P_FIELD)  ! OUTPUT PMSL
     *,TROP_T(P_FIELD)  ! OUTPUT TEMPS OF TROPOPAUSE
     *,TROP_P(P_FIELD)  ! OUTPUT PRESSURE OF TROPOPAUSE
     *,TROP_Z(P_FIELD)  ! OUTPUT HEIGHT OF TROPOPAUSE PRESSURE SURFACE
     *,TROP_ICAO(P_FIELD)  ! OUTPUT ICAO HT OF TROPOPAUSE PRESSURE
     *,HEIGHTS2(P_FIELD,H2_P_LEVS) ! OUTPUT height**2 on any P surface
     *,MODHT(P_FIELD,P_LEVELS)  ! OUTPUT GEOPOT. MODEL LEV HEIGHTS
C
      REAL
     * REL_HUMID_P(P_FIELD,REL_HUMID_P_LEVS)
     *                             ! OUTPUT HUMIDITIES ON ANY P SURFACE
     *,WBPT(P_FIELD,WBPT_LEVS)     ! WET BULB POTTEMP ON ANY P SURFACE
     *,SNPROB(P_FIELD)             ! SNOW PROBABILITY
     *,MINUS20_ICAO(P_FIELD)       ! MINUS 20 LEVEL ICAO HEIGHT
     *,MINUS20_P(P_FIELD)          ! MINUS 20 LEVEL PRESSURE
     *,FREEZE_Z(P_FIELD)           ! HEIGHT OF THE FREEZING LEVEL
     *,FREEZE_ICAO(P_FIELD)        ! ICAO HEIGHT OF THE FREEZING LEVEL
     *,FREEZE_P(P_FIELD)           ! PRESSURE OF THE FREEZING LEVEL
     *,CONTRAIL_UPPERHT(P_FIELD)   ! UPPER VALUE OF THE CONTRAIL
     *,CONTRAIL_LOWERHT(P_FIELD)   ! LOWER VALUE OF THE CONTRAIL
     *,MINUS20_Z(P_FIELD)          ! MINUS 20 LEVEL TRUE HEIGHT
     *,TH_ADV_SINGLE(P_FIELD,TH_ADV_P_LEVS)! THERMAL ADV SINGLE LEVELS
     *,DUCT_HEIGHT(P_FIELD)        ! RADIO DUCT HEIGHT
     *,DUCT_INT(P_FIELD)           ! RADIO DUCT INTENSITY
     *,TH_ADV_MEAN(P_FIELD)        ! THERMAL ADVECTION MEAN 850/700/500
     *,TRACER_P(TR_P_FIELD_DA)     ! TRACER FIELD ON ANY P SURFACE
     *,TRACERWORK(INT16)       ! STASHWORK FOR SUBSTITUTION OF TRACERS
C----------------------------------------------------------------------
C            AK,BK  DEFINE HYBRID VERTICAL COORDINATES P=A+BP*,
C       DELTA_AK,DELTA_BK  DEFINE LAYER PRESSURE THICKNESS PD=AD+BDP*,
C----------------------------------------------------------------------
      REAL
     * AKH(P_LEVELS+1)             !IN    value at layer boundary
     *,BKH(P_LEVELS+1)             !IN    value at layer boundary
     *,AK (P_LEVELS)               !IN    VALUE AT LAYER CENTRE
     *,BK (P_LEVELS)               !IN    VALUE AT LAYER CENTRE
     *,EW_SPACE                    !IN    LONGITUDE GRID SPACING
     *,NS_SPACE                    !IN    LATITUDE GRID SPACING
     *,SEC_U_LATITUDE(U_FIELD)     !IN    1/COS(LAT) AT U POINTS
     *,TIMESTEP                    !IN    TIMESTEP
     *,T_P_PRESS(T_P_LEVS)         !IN Pressures reqd for Temp inter
     *,HTS_PRESS(HTS_LEVS)         !IN     "       "   "   HTS
     *,REL_HUMID_PRESS(REL_HUMID_P_LEVS)!IN "   "   "   REL Humidity
     *,WBPT_PRESS(WBPT_LEVS)       !IN     "       "   "   W.B.P.T
     *,TH_ADV_PRESS(TH_ADV_P_LEVS) !IN     "       "   "   Thermal adv
     *,TR_PRESS(TR_VARS+1,NUM_STASH_LEVELSDA)
     *                             !IN Pressures reqd for Tracer interp
C*----------------------------------------------------------------------

C*L  WORKSPACE USAGE:---------------------------------------------------
C   REAL ARRAYS REQUIRED AT FULL FIELD LENGTH  4*(P_LEVELS)+7+Q_LEVELS
C-----------------------------------------------------------------------
      REAL
     * MODEL_HALF_HEIGHT(P_FIELD,P_LEVELS+1) !OUT HEIGHTS OF MODEL HALF
     *, WORK1(P_FIELD,P_LEVELS+1)  ! Will hold Rel Humid and P_HALF
     *                             ! and temperatures for FREEZE
     *                             ! and pressures for THADV
     *,WORK2(P_FIELD)         ! Workspace for Mean thermal advection
     *, PHI_STAR(P_FIELD)     ! Geopotential
     *, P(P_FIELD,P_LEVELS)   ! Pressure ARRAY
     *, PP                    ! Intermediate PRESSURE
     *, P_EXNER_FULL_L        ! Full Exner pressure for a point/level
     *, BOTTOM_QC_PRESS       ! pressure of the bottom QC level
     *, TOP_QC_PRESS          ! pressure of the top QC level
     *, PZ(P_FIELD)           ! Pressure surface in which results reQD
     *, T(P_FIELD)            ! Temperature array
     *, Q_SAT(P_FIELD) !
     *, DUMMY1
     *, DUMMY2
C
C*---------------------------------------------------------------------
C
C*L EXTERNAL SUBROUTINES CALLED------------%--------------------------
      EXTERNAL V_INT_Z,V_INT_ZH,PMSL,V_INT_T,P_TO_UV,TROP,V_INT,QSAT
      EXTERNAL QCTROP,THETAW,ICAO_HT,SNOWPR,FREEZE,DUCT,CONTRAIL
      EXTERNAL DIA_THADV
      EXTERNAL V_INT_TP
C*------------------------------------------------------------------
CL  MAXIMUM VECTOR LENGTH ASSUMED IS (ROWS+1) * ROWLENGTH
C----------------------------------------------------------------------
C*L  DEFINE LOCAL VARIABLES
      INTEGER
     *  P_POINTS      !     NUMBER OF P POINTS NEEDED
     *, ROWS_P        !     NUMBER OF P ROWS   NEEDED
     *, U_POINTS      !     NUMBER OF U POINTS UPDATED
     *, START_U       !     START POSITION FOR U POINTS UPDATED
     *, VERT_POINTS   !     NUMBER OF POINTS NON-ZERO DIFFUSION COEFFS
     *, ISL           !     LOCAL COUNTER
     *, BL            !     BOTTOM LEVEL REQD
     *, TL            !     TOP LEVEL REQD
     *, T_REF         !     REF level for below sfce Temp extrapola
     *, IPOINT        !     Reference point for checking output
     &, P_FLD_VALID   !     No of P points excluding halos & unused rows

      REAL
     *  PRESS_REQD    !     The pressure required for THETAW and THADV
     *, T0            !     Temperature required for FREEZE
     *, CP_OVER_G     !     Used in model level heights calculation
C-----------------------------------------------------------------------
C R IS GAS CONSTANT FOR DRY AIR
C CP IS SPECIFIC HEAT OF DRY AIR AT CONSTANT PRESSURE
C PREF IS REFERENCE SURFACE PRESSURE
C-----------------------------------------------------------------------
      INTEGER    K,I,II,IK,ITR ! LOOP COUNTERS IN ROUTINE

      REAL
     &    PU,PL
        PARAMETER (CP_OVER_G = CP / G)
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
C*L------------------ COMDECK P_EXNERC ---------------------------------
C statement function to define exner pressure at model levels
C from exner pressures at model half-levels
      REAL P_EXNER_C
      REAL R_P_EXNER_C
      REAL P_EXU_DUM,P_EXL_DUM,PU_DUM,PL_DUM,KAPPA_DUM
C consistent with geopotential see eqn 26 DOC PAPER 10
      P_EXNER_C(P_EXU_DUM,P_EXL_DUM,PU_DUM,PL_DUM,KAPPA_DUM) =
     & (P_EXU_DUM*PU_DUM - P_EXL_DUM*PL_DUM)/
     & ( (PU_DUM-PL_DUM)*(KAPPA_DUM + 1) )
      R_P_EXNER_C(P_EXU_DUM,P_EXL_DUM,PU_DUM,PL_DUM,KAPPA_DUM) =
     & ( (PU_DUM-PL_DUM)*(KAPPA_DUM + 1) )/
     & (P_EXU_DUM*PU_DUM - P_EXL_DUM*PL_DUM)

C*------------------- --------------------------------------------------

C*----------------------------------------------------------------------

      T_REF=2  ! This value is used in the vertical interpolation and is
C         set to the 2nd model level. This is not dependant on vert res
      ICODE=0

! Size of p-fields excluding halos and unused rows
      P_FLD_VALID=LAST_P_FLD_PT-FIRST_FLD_PT+1

          CALL SWAPBOUNDS(P_EXNER_HALF,ROW_LENGTH,P_ROWS,
     &    Offx,Offy,P_LEVELS+1)
C-----------------------------------------------------------------------
C     Calculation of pressure at ETA levels
C-----------------------------------------------------------------------
      DO K=1,P_LEVELS
        DO I=1,P_FIELD
          P(I,K)=AK(K)+BK(K)*PSTAR(I)
        ENDDO
      ENDDO
C-----------------------------------------------------------------------
      DO I=1,P_FIELD
        PHI_STAR(I)=OROG(I)*G
      ENDDO
C     set logical if model_half_heights needed for later
      IF(QHEIGHTS   .OR.
     *   QMODHT     .OR.
     *   QTROP_Z    .OR.
     *   QSNPROB    .OR.
     *   QFREEZE_Z  .OR.
     *   QMINUS20_Z )    QMODEL_HALF_HEIGHT=.TRUE.
C
      IF(QMODEL_HALF_HEIGHT) THEN
        CALL V_INT_ZH(P_EXNER_HALF,THETA,Q,PHI_STAR,
     *  MODEL_HALF_HEIGHT,P_FIELD,P_LEVELS,Q_LEVELS)
      ELSEIF(QHEIGHTS) THEN
        WRITE(6,100)
        WRITE(6,100)
  100   FORMAT('  WARNING MODEL_HALF_HEIGHTS NOT SWITCHED ON')
        ICODE=1
        CMESSAGE='PHY_DIAG:to calculate HTs,model half HTs needed first'
        GOTO 1000
      ENDIF
C-----------------------------------------------------------------------
CL--------------------------Section 0-----------------------------------
CL------------------GEOPOTENTIAL MODEL LEVEL HEIGHTS ITEM 281-----------
CL
CL   section 0 calculation of geopotential height of a model level
CL   surface. Heights in M
C-----------------------------------------------------------------------
      IF (QMODHT) THEN
        DO K=1,Q_LEVELS
          DO I=FIRST_FLD_PT,LAST_P_FLD_PT
            PU = AKH(K+1)+BKH(K+1)*PSTAR(I)
            PL = AKH(K)+BKH(K)*PSTAR(I)
            P_EXNER_FULL_L = P_EXNER_C( P_EXNER_HALF(I,K+1),
     &                 P_EXNER_HALF(I,K),PU,PL,KAPPA )
            MODHT(I,K) = MODEL_HALF_HEIGHT(I,K) + CP_OVER_G*
     &     (1.0+C_VIRTUAL*Q(I,K))*THETA(I,K)
     &    *(P_EXNER_HALF(I,K) - P_EXNER_FULL_L)
          ENDDO
        ENDDO
        IF(P_LEVELS .GT. Q_LEVELS) THEN
          DO K=Q_LEVELS+1,P_LEVELS
            DO I=FIRST_FLD_PT,LAST_P_FLD_PT
              PU = AKH(K+1)+BKH(K+1)*PSTAR(I)
              PL = AKH(K)+BKH(K)*PSTAR(I)
              P_EXNER_FULL_L = P_EXNER_C( P_EXNER_HALF(I,K+1),
     &                   P_EXNER_HALF(I,K),PU,PL,KAPPA )
              MODHT(I,K) = MODEL_HALF_HEIGHT(I,K) + CP_OVER_G*
     &        THETA(I,K)
     &      *(P_EXNER_HALF(I,K) - P_EXNER_FULL_L)
            ENDDO
          ENDDO
        END IF
      ENDIF
C-----------------------------------------------------------------------
CL--------------------------Section 1-----------------------------------
CL------------------GEOPOTENTIAL HEIGHTS ITEM 202-----------------------
CL
CL   section 2 calculation of geopotential height of an arbitary
CL   pressure surface. Pressure input is in Pascals and heights in M
C-----------------------------------------------------------------------
C   First check the height codes are valid and point to STASH_PRESSURE
C-----------------------------------------------------------------------
      IF(QHEIGHTS)  THEN
        DO K=1,HTS_LEVS
          DO I=FIRST_FLD_PT,LAST_P_FLD_PT
            PZ(I)=HTS_PRESS(K)*100.0    ! convert to pascals
          ENDDO
          CALL V_INT_Z(PZ,P(1,Z_REF),PSTAR,P_EXNER_HALF,THETA,Q,
     &    MODEL_HALF_HEIGHT,HEIGHTS(1,K),P_FIELD,P_LEVELS,Q_LEVELS,
     &    Z_REF,AKH,BKH,FIRST_FLD_PT,LAST_P_FLD_PT)
        ENDDO
      ENDIF
C-----------------------------------------------------------------------
CL--------------------------Section 2-----------------------------------
CL------------------MEAN SEA LEVEL PRESSURE ITEM 222--------------------
C-----------------------------------------------------------------------
      IF(QPMSL) THEN
        CALL PMSL(P_MSL,P(1,Z_REF),PSTAR,P_EXNER_HALF,THETA,Q,PHI_STAR
     &    ,P_FIELD,P_LEVELS,Q_LEVELS,Z_REF,AKH,BKH,FIRST_FLD_PT
     &    ,LAST_P_FLD_PT)
      ENDIF
C-----------------------------------------------------------------------
CL--------------------------Section 3-----------------------------------
CL------------------TEMPERATURES ON PRESSURE ITEM 203-------------------
CL
CL   Pressure input is in Pascals and potential temperature is input
CL                                                          in Kelvin
C-----------------------------------------------------------------------
      IF(QT_P) THEN
        DO K=1,T_P_LEVS
          DO I=1,P_FIELD
            PZ(I)=T_P_PRESS(K)*100.0    ! convert to pascals
          ENDDO
!         From Vn 4.5, L_VINT_TP is set through UMUI
          IF (L_VINT_TP) THEN
           CALL V_INT_TP(T_P(1,K),PZ,P(1,T_REF),PSTAR,P_EXNER_HALF,
     &     THETA,P_FIELD,P_LEVELS,T_REF,AKH,BKH)
          ELSE
           CALL V_INT_T(T_P(1,K),PZ,P(1,T_REF),PSTAR,P_EXNER_HALF,
     &        THETA,P_FIELD,P_LEVELS,T_REF,AKH,BKH,FIRST_FLD_PT
     &        ,LAST_P_FLD_PT)
          ENDIF
        ENDDO
      ENDIF
C-----------------------------------------------------------------------
CL--------------------------Section 4-----------------------------------
CL------------------RELATIVE HUM ON PRESSURE ITEM 204-------------------
C-----------------------------------------------------------------------
      IF (QREL_HUMID_P) THEN
C-----------------------------------------------------------------------
C Calculation of temperature at theta points
C Calculation of relative humidity at pressure points
C-----------------------------------------------------------------------
        DO K=1,Q_LEVELS
          DO I=FIRST_FLD_PT,LAST_P_FLD_PT
            PU=PSTAR(I)*BKH(K+1) + AKH(K+1)
            PL=PSTAR(I)*BKH(K) + AKH(K)
            P_EXNER_FULL_L=
     &        P_EXNER_C(P_EXNER_HALF(I,K+1),P_EXNER_HALF(I,K),
     &                                                     PU,PL,KAPPA)
            T(I)=THETA(I,K)*P_EXNER_FULL_L
          ENDDO
          CALL QSAT(Q_SAT(FIRST_FLD_PT),T(FIRST_FLD_PT)
     &      ,P(FIRST_FLD_PT,K),P_FLD_VALID)
          DO I=FIRST_FLD_PT,LAST_P_FLD_PT
            WORK1(I,K)=Q(I,K)/Q_SAT(I)*100.0   ! Rel humidity in WORK1
! Relative humidity should NOT be limited to 100% but do so in old
! versions of the cloud scheme to maintain bit comparison.
            IF(WORK1(I,K).GT.100.0 .AND. .NOT. L_LSPICE) THEN
               WORK1(I,K)=100.0
            ELSEIF(WORK1(I,K).LT.0.) THEN
               WORK1(I,K)=0.0
            ENDIF
          ENDDO
        ENDDO
C-----------------------------------------------------------------------
CL    Interpolate relative humidity on a pressure surface
C-----------------------------------------------------------------------
        DO K=1,REL_HUMID_P_LEVS
          DO I=FIRST_FLD_PT,LAST_P_FLD_PT
            PZ(I)=REL_HUMID_PRESS(K)*100.0   ! convert to Pascals
          ENDDO
          CALL V_INT(P,PZ,WORK1,REL_HUMID_P(1,K),
     &    P_FIELD,Q_LEVELS,DUMMY1,DUMMY2,.FALSE.
     &    ,FIRST_FLD_PT,LAST_P_FLD_PT)
        ENDDO
      ENDIF
C-----------------------------------------------------------------------
C------------------TROPOPAUSE  ITEMS 214 215 216 217--------------------
CL----------Section 5 calculation of tropopause temp/ht and pressure----
C-----------------------------------------------------------------------
      IF (QTROP_T.AND.QTROP_P.AND.QTROP_Z) THEN
        IF(QMODEL_HALF_HEIGHT) THEN
C-----------------------------------------------------------------------
C      Call tropopause program with quality control
C-----------------------------------------------------------------------
          BOTTOM_QC_PRESS=60000.0  ! pressure of the bottom QC level
          TOP_QC_PRESS=10000.0     ! pressure of the top QC level
          BOTTOM_QC_LEVEL=1
          TOP_QC_LEVEL=P_LEVELS
          DO K=1,P_LEVELS
            PP= AK(K) + BK(K) * 100000 ! Pstar assumed to be 1000MB
            IF(BOTTOM_QC_PRESS.LT.PP)  THEN
              BOTTOM_QC_LEVEL=K
            ENDIF
            IF(TOP_QC_PRESS.LT.PP)  THEN
              TOP_QC_LEVEL=K
            ENDIF
          ENDDO
          DO K=1,P_LEVELS+1
            DO I=1,P_FIELD
              WORK1(I,K)=AKH(K)+BKH(K)*PSTAR(I)  ! Calculate P_HALF
            ENDDO
          ENDDO
          CALL QCTROP(THETA,WORK1,P_EXNER_HALF
     &   ,MODEL_HALF_HEIGHT,TROP_T,TROP_P,TROP_Z,P_FIELD,P_LEVELS
     &   ,P(1,Z_REF),P(1,T_REF),PSTAR,Q,Q_LEVELS,Z_REF,T_REF,Z_REF
     &   ,BOTTOM_QC_LEVEL,TOP_QC_LEVEL,AKH,BKH
     &   ,FIRST_FLD_PT,LAST_P_FLD_PT)
C         Last Z_REF used as Min TROP L

C-----------------------------------------------------------------------
C This piece of code is to ensure sensible values of TROP ht and TEMP.
C Although,the TROP height is also produced it can be rubbish so for
C the moment it is best to call V_INTZ
C-----------------------------------------------------------------------
          DO I=FIRST_FLD_PT,LAST_P_FLD_PT
            IF(TROP_T(I).LT.0) THEN
              TROP_T(I)=199.0
            ENDIF
            IF(TROP_P(I).LT.0) THEN
              TROP_P(I)=111.0
            ENDIF
          ENDDO
C         CALL V_INT_Z(TROP_P,P(1,Z_REF),PSTAR,P_EXNER_HALF,THETA,Q,
C    &    MODEL_HALF_HEIGHT,TROP_Z,P_FIELD,P_LEVELS,Q_LEVELS,Z_REF,
C    &    AKH,BKH,FIRST_FLD_PT,LAST_P_FLD_PT)
C-----------------------------------------------------------------------
C         Calculates ICAO heights from Trop pressure field
C-----------------------------------------------------------------------
          IF(QTROP_ICAO) THEN
            CALL ICAO_HT(TROP_P(FIRST_FLD_PT),P_FLD_VALID
     &        ,TROP_ICAO(FIRST_FLD_PT))
          ENDIF
C
        ELSE
          WRITE(6,444)
          WRITE(6,444)
 444      FORMAT(' Subroutine TROP not called No MODEL_HALF_HEIGHTS')
        ENDIF
      ELSEIF(QTROP_T.NEQV.QTROP_P.OR.QTROP_T.NEQV.QTROP_Z.OR.QTROP_P.
     &  NEQV.QTROP_Z) THEN
        WRITE(6,*)'Subroutine TROP not called-tropopause temperature,'
        WRITE(6,*)' pressure and height all need to be selected'
      ENDIF
C-----------------------------------------------------------------------
CL-------------------SECTION 6 ITEM 205---------------------------------
CL          Calculation of Wet Bulb Potential Temperature
C-----------------------------------------------------------------------
      IF (QWBPT) THEN
        DO K=1,WBPT_LEVS
          PRESS_REQD=WBPT_PRESS(K)*100.0   ! convert to Pascals
          ICODE=0
          CALL THETAW
     1    (PRESS_REQD,THETA,Q,P,PSTAR,P_EXNER_HALF,WBPT(1,K),
     2    AK,BK,AKH,BKH,P_FIELD,P_LEVELS,Q_LEVELS,T_REF,
     &    FIRST_FLD_PT,LAST_P_FLD_PT,
     3    ICODE,CMESSAGE)
          IF(ICODE.NE.0) THEN
            ICODE=1
            CMESSAGE='PHY_DIAG error in calculation of WBPT'
            GOTO 1000
          ENDIF
        ENDDO
      ENDIF
C-----------------------------------------------------------------------
CL---------------------SECTION 7 ITEM 206-------------------------------
CL             Calculation of Snow probability
C-----------------------------------------------------------------------
      IF (QSNPROB)THEN
        CALL SNOWPR
     1  (P,PSTAR,P_EXNER_HALF,THETA,Q,MODEL_HALF_HEIGHT,
     2  P_FIELD,P_LEVELS,Q_LEVELS,Z_REF,AKH,BKH,
     &  SNPROB,FIRST_FLD_PT,LAST_P_FLD_PT)
      ENDIF
C-----------------------------------------------------------------------
CL--------------------SECTION 8 ITEMS 207/8/18--------------------------
CL     Calculation of heights and pressure of -20 isotherm
C-----------------------------------------------------------------------
      IF (QMINUS20_P.AND.QMINUS20_Z) THEN
        IF(QMODEL_HALF_HEIGHT)THEN
C-----------------------------------------------------------------------
C         Set T0 to 253.16K (-20 deg C)
C-----------------------------------------------------------------------
          T0=253.16
C-----------------------------------------------------------------------
C         Calculate temperatures at full levels (stored in WORK1)
C-----------------------------------------------------------------------
          DO K=1,P_LEVELS
            DO I=FIRST_FLD_PT,LAST_P_FLD_PT
              PU=PSTAR(I)*BKH(K+1) + AKH(K+1)
              PL=PSTAR(I)*BKH(K) + AKH(K)
              P_EXNER_FULL_L=
     &          P_EXNER_C(P_EXNER_HALF(I,K+1),P_EXNER_HALF(I,K),
     &                                                     PU,PL,KAPPA)
             WORK1(I,K)=THETA(I,K)*P_EXNER_FULL_L
            ENDDO
          ENDDO
C-----------------------------------------------------------------------
          CALL FREEZE
     1    (T0,P,THETA,WORK1,P_EXNER_HALF,PSTAR,Q,MODEL_HALF_HEIGHT,
     2    OROG,MINUS20_Z,MINUS20_P,
     &    P_FIELD,P_LEVELS,Q_LEVELS,Z_REF,AKH,BKH
     &    ,FIRST_FLD_PT,LAST_P_FLD_PT)
C-----------------------------------------------------------------------
          IF(QMINUS20_ICAO)THEN
            CALL ICAO_HT(MINUS20_P(FIRST_FLD_PT),P_FLD_VALID
     &        ,MINUS20_ICAO(FIRST_FLD_PT))
          ENDIF
C-----------------------------------------------------------------------
        ELSE
          WRITE(6,556)
          WRITE(6,556)
  556     FORMAT(' Subroutine FREEZE not called-no MODEL_HALF_HEIGHTs')
        ENDIF
      ELSEIF(QMINUS20_P.NEQV.QMINUS20_Z) THEN
        WRITE(6,*)'Subroutine FREEZE not called-both -20 level pressure'
        WRITE(6,*)'and height need to be selected'
      ENDIF
C-----------------------------------------------------------------------
CL-------------------SECTION 9 ITEMS 209-211----------------------------
CL       Calculation of Freezing level heights and pressure
C-----------------------------------------------------------------------
      IF (QFREEZE_Z.AND.QFREEZE_P)THEN
        IF(QMODEL_HALF_HEIGHT)THEN
C-----------------------------------------------------------------------
C         Set T0 to 273.16K (0 deg C)
C-----------------------------------------------------------------------
          T0=273.16
C-----------------------------------------------------------------------
C         Calculate temperatures at full levels (stored in WORK1)
C-----------------------------------------------------------------------
          DO K=1,P_LEVELS
            DO I=FIRST_FLD_PT,LAST_P_FLD_PT
              PU=PSTAR(I)*BKH(K+1) + AKH(K+1)
              PL=PSTAR(I)*BKH(K) + AKH(K)
              P_EXNER_FULL_L=
     &          P_EXNER_C(P_EXNER_HALF(I,K+1),P_EXNER_HALF(I,K),
     &                                                     PU,PL,KAPPA)
             WORK1(I,K)=THETA(I,K)*P_EXNER_FULL_L
            ENDDO
          ENDDO
C-----------------------------------------------------------------------
          CALL FREEZE
     1    (T0,P,THETA,WORK1,P_EXNER_HALF,PSTAR,Q,MODEL_HALF_HEIGHT,
     2    OROG,FREEZE_Z,FREEZE_P,
     &    P_FIELD,P_LEVELS,Q_LEVELS,Z_REF,AKH,BKH
     &    ,FIRST_FLD_PT,LAST_P_FLD_PT)
          IF(QFREEZE_ICAO)THEN
            CALL ICAO_HT(FREEZE_P(FIRST_FLD_PT),P_FLD_VALID
     &        ,FREEZE_ICAO(FIRST_FLD_PT))
          ENDIF
C-----------------------------------------------------------------------
        ELSE
          WRITE(6,555)
          WRITE(6,555)
  555     FORMAT(' Subroutine FREEZE not called-no MODEL_HALF_HEIGHTs')
        ENDIF
      ELSEIF(QFREEZE_P.NEQV.QFREEZE_Z) THEN
        WRITE(6,*)'Subroutine FREEZE not called - both freezing level'
        WRITE(6,*)'pressure and height need to be selected'
      ENDIF
C-----------------------------------------------------------------------
CL-----------------SECTION 10 ITEMS 220/1-------------------------------
CL         Calculation of duct height and intensity
C-----------------------------------------------------------------------
      IF(QDUCT_INT.AND.QDUCT_HEIGHT)THEN
        CALL DUCT
     1  (PSTAR,TSTAR,THETA,Q,U,V,
     2   AK,BK,AKH,BKH,LAND,DUCT_HEIGHT,DUCT_INT,P_EXNER_HALF,
     3  P_LEVELS,Q_LEVELS,ROW_LENGTH,P_ROWS,U_ROWS,P_FIELD,U_FIELD)
      ELSEIF(QDUCT_INT.NEQV.QDUCT_HEIGHT)THEN
        WRITE(6,*)'Subroutine DUCT not called - both duct height and'
        WRITE(6,*)'intensity need to be selected'
      ENDIF
C-----------------------------------------------------------------------
CL-----------------SECTION 11 ITEMS 212/3-------------------------------
CL         Calculation of contrail lower and upper height
C-----------------------------------------------------------------------
      IF(QCONTRA_LOWERHT.AND.QCONTRA_UPPERHT)THEN
C-----------------------------------------------------------------------
C         Calculate temperatures at full levels (stored in WORK1)
C-----------------------------------------------------------------------
        DO K=1,P_LEVELS
          DO I=1,P_FIELD
            PU=PSTAR(I)*BKH(K+1) + AKH(K+1)
            PL=PSTAR(I)*BKH(K)   + AKH(K)
            P_EXNER_FULL_L=
     &        P_EXNER_C(P_EXNER_HALF(I,K+1),P_EXNER_HALF(I,K),
     &                                                     PU,PL,KAPPA)
            WORK1(I,K)=THETA(I,K)*P_EXNER_FULL_L
          ENDDO
        ENDDO
C-----------------------------------------------------------------------
        CALL CONTRAIL
     1  (P,WORK1,PSTAR,P_EXNER_HALF,THETA,CONTRAIL_UPPERHT,
     2  CONTRAIL_LOWERHT,P_FIELD,P_LEVELS,FIRST_FLD_PT,LAST_P_FLD_PT)
      ELSEIF(QCONTRA_UPPERHT.NEQV.QCONTRA_LOWERHT)THEN
      WRITE(6,*)'Subroutine CONTRAIL not called - both lower height and'
        WRITE(6,*)'upper height need to be selected'
      ENDIF
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
CL-----------------SECTION 12 ITEMS 219/223-----------------------------
CL         Calculation of thermal advection
C-----------------------------------------------------------------------
      IF(QTH_ADV_SINGLE.OR.QTH_ADV_MEAN)THEN
C-----------------------------------------------------------------------
CL    Interpolate P field from P/T points to U/V points for all levels.
C-----------------------------------------------------------------------
        DO K=1,P_LEVELS
          CALL P_TO_UV
     1    (P(1,K),WORK1(1,K),P_FIELD,U_FIELD,ROW_LENGTH,
     2     P_ROWS)
        ENDDO
      ENDIF
C       Change effective dimension of WORK1 for input to DIA_THADV
C       (only first U_FIELD points defined for each level from P_TO_UV)
       CALL CHANGE_DIMENS(WORK1,P_FIELD,U_FIELD,P_LEVELS,ICODE)
C-----------------------------------------------------------------------
CL    Item 219 single level thermal advection---------------------------
C-----------------------------------------------------------------------
      IF(QTH_ADV_SINGLE)THEN
C-----------------------------------------------------------------------
CL    Loop over required levels
C-----------------------------------------------------------------------
        DO K=1,TH_ADV_P_LEVS
          PRESS_REQD=TH_ADV_PRESS(K)*100.0     ! Convert to pascals
          CALL DIA_THADV
     1    (U,V,THETA,P,WORK1,PSTAR,PRESS_REQD,P_EXNER_HALF,
     2    TH_ADV_SINGLE(1,K),P_FIELD,U_FIELD,P_LEVELS,ROW_LENGTH,
     3    P_ROWS,SEC_U_LATITUDE,EW_SPACE,NS_SPACE,AKH,BKH)

        ENDDO
      ENDIF
C-----------------------------------------------------------------------
CL    Item 223 mean thermal advection over levels 850/700/500mb-
C-----------------------------------------------------------------------
      IF(QTH_ADV_MEAN)THEN
C-----------------------------------------------------------------------
C     Initialise array
C-----------------------------------------------------------------------
        DO I=FIRST_FLD_PT,LAST_P_FLD_PT
          TH_ADV_MEAN(I)=0.0
        ENDDO
C-----------------------------------------------------------------------
CL    Loop over levels 850/700/500mb
C-----------------------------------------------------------------------
        DO K=1,3
          IF(K.EQ.1)PRESS_REQD=85000.0
          IF(K.EQ.2)PRESS_REQD=70000.0
          IF(K.EQ.3)PRESS_REQD=50000.0
          CALL DIA_THADV
     1    (U,V,THETA,P,WORK1,PSTAR,PRESS_REQD,P_EXNER_HALF,WORK2,
     2    P_FIELD,U_FIELD,P_LEVELS,ROW_LENGTH,P_ROWS,
     3    SEC_U_LATITUDE,EW_SPACE,NS_SPACE,AKH,BKH)
          DO I=FIRST_FLD_PT,LAST_P_FLD_PT
            TH_ADV_MEAN(I)=TH_ADV_MEAN(I)+WORK2(I)
          ENDDO
        ENDDO
        DO I=FIRST_FLD_PT,LAST_P_FLD_PT
          TH_ADV_MEAN(I)=TH_ADV_MEAN(I)/3.0
        ENDDO
      ENDIF
C-----------------------------------------------------------------------
CL    Section 13 ITEMS 224 products of other fields
C-----------------------------------------------------------------------
CL   Item 224 - height**2 on pressure levels
CL   Must be requested on identical levels to height on pressure
CL   levels.
C
      IF (QHTS2) THEN
          DO K=1,H2_P_LEVS
            DO I=FIRST_FLD_PT,LAST_P_FLD_PT
              HEIGHTS2(I,K) = HEIGHTS(I,H2_IND(K))*HEIGHTS(I,H2_IND(K))
            ENDDO
          ENDDO
      ENDIF
C-----------------------------------------------------------------------
CL-----------------SECTION 14 ITEMS 226-254-----------------------------
CL         Interpolate tracers onto pressure levels
C-----------------------------------------------------------------------
      IF (TR_VARS.GT.0) THEN
        DO ITR=1,TR_VARS
          IF (SF_TRACER(ITR)) THEN
            DO K=1,TR_PRESS_LEVS(ITR)
              DO I=FIRST_FLD_PT,LAST_P_FLD_PT
                PZ(I)=TR_PRESS(ITR,K)*100.0    ! convert to Pascals
C-----------------------------------------------------------------------
C   Extract current tracer field from STASHWORK/TRACERWORK
C-----------------------------------------------------------------------
                TRACER_P(I)=TRACERWORK( PT_TRACER(ITR) - 1 + I +
     &                                ( P_FIELD * (K-1) ) )
              ENDDO
              CALL V_INT(P,PZ,TRACERS(1,1,ITR),TRACER_P,
     &                   P_FIELD,TR_LEVELS,DUMMY1,DUMMY2,.FALSE.
     &                   ,FIRST_FLD_PT,LAST_P_FLD_PT)
C-----------------------------------------------------------------------
C   Copy interpolated tracer field back to STASHWORK/TRACERWORK
C-----------------------------------------------------------------------
              DO I=FIRST_FLD_PT,LAST_P_FLD_PT
                TRACERWORK( PT_TRACER(ITR) - 1 + I +
     &                    ( P_FIELD * (K-1) ) ) = TRACER_P(I)
              ENDDO
            ENDDO
          ENDIF
        ENDDO
      ENDIF
C-----------------------------------------------------------------------
1000  CONTINUE
      RETURN
      END
