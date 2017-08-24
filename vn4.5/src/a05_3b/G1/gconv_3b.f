C ******************************COPYRIGHT******************************
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
C
!   Subroutine GLUE_CONV--------------------------------------
!
!   Level 3 control routine
!
!   Purpose: Calls CONVECT to calculate and add convection increments.
!            It is an extra level of Control routine to avoid the
!            use of *IF DEF around calls to different CONVECT routines
!            as per proposal of S. Foreman 22/8/94 for plug compatible
!            routines.
!       Test version with tracers, momentum transports, CAPE closure.
!
!   Code Description:
!    Language FORTRAN 77 + extensions.
!    This code is written to UMDP3 v6 programming standards.
!
!         Author: P.Inness        Reviewer: C. Wilson
!
!   Modification History:
!
!   Version      Date
!
!     4.3        03/02/97   New deck   sdjackson
!     4.4        26/09/97   Provision for 3D conv cloud amount and rain
!                           out before calculation of water path. J.M.G
!LL  4.4  Oct 97    Add halo mask to stop redundant calculations
!LL                                               Alan Dickinson
!
!   System components covered : P27
!
!   System task : P0
!
!   Documentation:
!
!  END -----------------------------------------------------------------

      SUBROUTINE GLUE_CONV(
     &     NP_FIELD,NPNTS,NLEV,NBL,
     &     TH,Q,PSTAR,BLAND,U,V,TRACER,
     &     DTHBYDT,DQBYDT,DUBYDT,DVBYDT,
     &     RAIN,SNOW,
     &     CCA,ICCB,ICCT,CCLWP,CCW,
     &     ICCBPxCCA,ICCTPxCCA,GBMCCWP,GBMCCW,
     &     LCBASE,LCTOP,LCCA,LCCLWP,
     &     CAPE,
     &     EXNER,AK,BK,AKM12,BKM12,DELAK,DELBK,
     &     TIMESTEP,
     &     T1_SD,Q1_SD,
     &     L_MOM,L_TRACER,L_CAPE,NTRA,TRLEV,L_XSCOMP,
     &     L_SDXS
     &     ,l_halo
     &    ,N_CCA_LEV, L_3D_CCA, L_CCW, MPARWTR
     &    , ANVIL_FACTOR, TOWER_FACTOR, UD_FACTOR, L_CLOUD_DEEP
     &    , L_PHASE_LIM, UP_FLUX,FLG_UP_FLX,DWN_FLUX,FLG_DWN_FLX
     &    , ENTRAIN_UP, FLG_ENTR_UP,DETRAIN_UP,FLG_DETR_UP,ENTRAIN_DWN
     &    , FLG_ENTR_DWN,DETRAIN_DWN,FLG_DETR_DWN
     & )
      IMPLICIT NONE

!----------------------------------------------------------------------
!Some of the following variables are 'dummy' as they are used in other
! CONVECT versions
!-----------------------------------------------------------------------
!----------------------------------------------------------------------
! IN variables
!---------------------------------------------------------------------
      INTEGER NP_FIELD            ! LENGTH OF DATA (ALSO USED TO
                                  ! SPECIFY STARTING POINT OF
                                  ! DATA PASSED IN)
C
      INTEGER NPNTS               ! IN FULL VECTOR LENGTH
C
      INTEGER NLEV                ! IN NUMBER OF MODEL LAYERS
C
      INTEGER NBL                 ! IN NUMBER OF BOUNDARY LAYER LEVELS
C
      INTEGER NTRA                ! NUMBER OF TRACER FIELDS
C
      INTEGER TRLEV               ! NUMBER OF MODEL LEVELS ON WHICH
                                  ! TRACERS ARE INCLUDED
      INTEGER N_CCA_LEV           ! IN Number of convective cloud levels
C
      LOGICAL BLAND(NP_FIELD)     ! IN LAND/SEA MASK
C
      LOGICAL L_TRACER            ! IN SWITCH FOR INCLUSION OF TRACERS
C
      LOGICAL L_MOM               ! IN SWITCH FOR INCLUSION OF
                                  !    MOMENTUM TRANSPORTS
C
      LOGICAL L_XSCOMP            ! IN SWITCH FOR ALLOWING COMPENSATING
                                  !    COOLING AND DRYING OF THE
                                  !    ENVIRONMENT IN INITIATING LAYER
C
      LOGICAL L_SDXS              ! IN SWITCH FOR ALLOWING PARCEL EXCESS
                                  !    TO BE SET TO S.D. OF TURBULENT
                                  !    FLUCTUATIONS IN LOWEST MODEL
                                  !    LAYER
      LOGICAL L_CAPE              ! IN SWITCH FOR USE OF CAPE CLOSURE
C
      LOGICAL L_3D_CCA            ! IN SWITCH FOR USE OF 3D CLOUD AMOUNT
C
      LOGICAL L_CCW               ! IN IF .TRUE. THEN PRECIP NOT INC. IN
!                                 !    CONV. CLOUD WATER PATH.
!
      LOGICAL L_CLOUD_DEEP        ! IN If TRUE then use depth criterion
!                                 !    to determine existence of anvils
!
      LOGICAL L_PHASE_LIM         ! IN DUMMY variable added to match
!                                 !    argument lists.
C
      REAL PSTAR(NP_FIELD)        ! IN SURFACE PRESSURE (PA)
C
      REAL EXNER(NP_FIELD,NLEV+1) ! IN EXNER RATIO
C
      REAL AK(NLEV),              ! IN HYBRID CO-ORDINATE COEFFICIENTS
     *     BK(NLEV)               !    DEFINE PRESSURE AT MID-POINT
                                  !    OF LAYER K
C
      REAL AKM12(NLEV+1),         ! IN HYBRID CO-ORDINATE COEFFICIENTS
     *     BKM12(NLEV+1)          !    TO DEFINE PRESSURE AT
                                  !    LEVEL K-1/2
C
      REAL DELAK(NLEV),           ! IN DIFFERENCE IN HYBRID CO-ORDINATE
     *     DELBK(NLEV)            !    COEFFICIENTS ACROSS LAYER K
C
      REAL TIMESTEP               ! IN MODEL TIMESTEP (SECS)
C
      REAL U(NP_FIELD,NLEV)       ! IN MODEL U FIELD (M/S)
C
      REAL V(NP_FIELD,NLEV)       ! IN MODEL V FIELD (M/S)
C
      REAL T1_SD(NP_FIELD)        ! IN Standard deviation of turbulent
C                                 !    fluctuations of layer 1
C                                 !    temperature (K).
      REAL Q1_SD(NP_FIELD)        ! IN Standard deviation of turbulent
C                                 !    fluctuations of layer 1
C                                 !    humidity (kg/kg).
      REAL MPARWTR              ! IN Reservoir of convective cloud water
!                               !    left in a layer after conv. precip.
      REAL ANVIL_FACTOR         ! IN used in calculation of cloud amount
     &    ,TOWER_FACTOR         !    on model levels if L_3D_CCA = .T.
!
      REAL UD_FACTOR            ! IN Factor used in calculation of CCWP
!                               !    for radiation, if L_CCW is true.
!                                                                       
      LOGICAL l_halo(NP_FIELD)    ! IN Mask for halos
      LOGICAL FLG_UP_FLX          ! STASH FLAG FOR UPDRAUGHT MASS FLUX
C
      LOGICAL FLG_DWN_FLX         ! STASH FLAG FOR DOWNDRAGHT MASS FLUX
!
      LOGICAL FLG_ENTR_UP         ! STASH FLAG FOR UPDRAUGHT ENTRAINMENT
!
      LOGICAL FLG_ENTR_DWN        ! STASH FLAG FOR DOWNDRAUGHT ENTRAINMN
!
      LOGICAL FLG_DETR_UP         ! STASH FLAG FOR UPDRAUGHT DETRAINMENT
!
      LOGICAL FLG_DETR_DWN        ! STASH FLAG FOR DOWNDRAUGHT DETRAINMN
!                                                                       
!----------------------------------------------------------------------
! INOUT variables
!----------------------------------------------------------------------
      REAL TH(NP_FIELD,NLEV)      ! INOUT
                                  ! IN MODEL POTENTIAL TEMPERATURE (K)
                                  ! OUT MODEL POTENTIAL TEMPERATURE
                                  !     AFTER CONVECTION (K)
C
      REAL Q(NP_FIELD,NLEV)       ! INOUT
                                  ! IN MODEL MIXING RATIO (KG/KG)
                                  ! OUT MODEL MIXING RATIO AFTER
                                  !     AFTER CONVECTION (KG/KG)
C
      REAL TRACER(NP_FIELD,TRLEV,NTRA)! INOUT
                                      ! IN  MODEL TRACER FIELDS (KG/KG)
                                      ! OUT MODEL TRACER FIELDS AFTER
                                      !     CONVECTION (KG/KG)
C

!----------------------------------------------------------------------
! OUT variables
!----------------------------------------------------------------------
      REAL DTHBYDT(NP_FIELD,NLEV) ! OUT INCREMENTS TO POTENTIAL
                                  !     TEMPERATURE DUE TO CONVECTION
                                  !     (K/S)
C
      REAL DQBYDT(NP_FIELD,NLEV)  ! OUT INCREMENTS TO MIXING RATIO
                                  !     DUE TO CONVECTION (KG/KG/S)
C
      REAL DUBYDT(NP_FIELD,NLEV)  ! OUT INCREMENTS TO U DUE TO
                                  !     CONVECTIVE MOMENTUM TRANSPORT
                                  !     (M/S**2)
C
      REAL DVBYDT(NP_FIELD,NLEV)  ! OUT INCREMENTS TO V DUE TO
                                  !     CONVECTIVE MOMENTUM TRANSPORT
                                  !     (M/S**2)
C
      REAL RAIN(NP_FIELD)         ! OUT SURFACE CONVECTIVE RAINFALL
                                  !     (KG/M**2/S)
C
      REAL SNOW(NP_FIELD)         ! OUT SURFACE CONVECTIVE SNOWFALL
                                  !     (KG/M**2/S)
C
      REAL CCA(NP_FIELD,N_CCA_LEV) ! OUT CONVECTIVE CLOUD AMOUNT (%)
!                                  !     2D or 3D depending on L_3D_CCA
C
      INTEGER ICCB(NP_FIELD)      ! OUT CONVECTIVE CLOUD BASE LEVEL
C
      INTEGER ICCT(NP_FIELD)      ! OUT CONVECTIVE CLOUD TOP LEVEL
C
      REAL CCLWP(NP_FIELD)        ! OUT CONDENSED WATER PATH (KG/M**2)
C
      REAL CCW(NP_FIELD,NLEV)     ! OUT CONVECTIVE CLOUD LIQUID WATER
                                  ! (G/KG) ON MODEL LEVELS
C
      REAL ICCBPxCCA(NP_FIELD)    ! OUT CONV. CLD BASE PRESSURE x CCA
C
      REAL ICCTPxCCA(NP_FIELD)    ! OUT CONV. CLD TOP PRESSURE x CCA
C
      REAL GBMCCWP(NP_FIELD)      ! OUT GRIDBOX MEAN CCWP
C
      REAL GBMCCW(NP_FIELD,NLEV)  ! OUT GRIDBOX MEAN CCW
C
      REAL LCCA(NP_FIELD)         ! OUT LOWEST CONV.CLOUD AMOUNT (%)
C
      INTEGER LCBASE(NP_FIELD)    ! OUT LOWEST CONV.CLOUD BASE LEVEL
C
      INTEGER LCTOP(NP_FIELD)     ! OUT LOWEST CONV.CLOUD TOP LEVEL
C
      REAL LCCLWP(NP_FIELD)       ! OUT CONDENSED WATER PATH (KG/M**2)
                                  !     FOR LOWEST CONV.CLOUD
C
      REAL CAPE(NPNTS)            ! OUT MODEL VALUES OF CONVECTIVE
                                  !     AVAILABLE POTENTIAL ENERGY
!
      REAL UP_FLUX(NP_FIELD,NLEV)     ! OUT UPDRAUGHT MASS FLUX
!
      REAL DWN_FLUX(NP_FIELD,NLEV)    ! OUT DOWNDRAUGHT MASS FLUX
!
      REAL ENTRAIN_UP(NP_FIELD,NLEV)  ! OUT FRACTIOAL ENTRAINMENT RATE F
                                      !     UPDRAUGHTS
!
      REAL DETRAIN_UP(NP_FIELD,NLEV)   ! OUT FRACTIONAL DETRAINMEN RATE 
                                       ! UPDRAUGHTS
!
      REAL ENTRAIN_DWN(NP_FIELD,NLEV)  ! OUT FRACTIONAL DETRAINMENT RATE
                                       !     DOWNDRAUGHTS
!
      REAL DETRAIN_DWN(NP_FIELD,NLEV)  ! OUT FRACTIONAL DETRAINMENT RATE
 
                                  !     FOR DIAGNOSTIC OUTPUT

! External subroutines called

      EXTERNAL
     &      CONVECT

! Local variables


!--------------- SECTION 5 CONVECTION ------------------------

!  5.2 Call CONVECT to calculate and add convection increments

      CALL CONVECT(
     &     NP_FIELD,NPNTS,NLEV,NBL,
     &     TH,Q,PSTAR,BLAND,U,V,TRACER,
     &     DTHBYDT,DQBYDT,DUBYDT,DVBYDT,
     &     RAIN,SNOW,
     &     CCA,ICCB,ICCT,CCLWP,CCW,
     &     ICCBPxCCA,ICCTPxCCA,GBMCCWP,GBMCCW,
     &     LCBASE,LCTOP,LCCA,LCCLWP,
     &     CAPE,
     &     EXNER,AK,BK,AKM12,BKM12,DELAK,DELBK,
     &     TIMESTEP,
     &     T1_SD,Q1_SD,
     &     L_MOM,L_TRACER,L_CAPE,NTRA,TRLEV,L_XSCOMP,
     &     L_SDXS
     &    ,N_CCA_LEV,L_3D_CCA,L_CCW,MPARWTR
     &    ,ANVIL_FACTOR ,TOWER_FACTOR
     &     ,l_halo
     &     ,UD_FACTOR, L_CLOUD_DEEP
     &     ,UP_FLUX,FLG_UP_FLX,DWN_FLUX,FLG_DWN_FLX,ENTRAIN_UP
     &     ,FLG_ENTR_UP,DETRAIN_UP,FLG_DETR_UP,ENTRAIN_DWN
     &     ,FLG_ENTR_DWN,DETRAIN_DWN,FLG_DETR_DWN
     & )

      RETURN
      END
