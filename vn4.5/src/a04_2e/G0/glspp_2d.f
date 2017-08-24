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
!   SUBROUTINE GLUE_LSPP------------------------------------------------
!
!   Level 3 control routine
!
!   Purpose: Calculate large-scale (dynamical) precipitation. LS_PPNC is
!            the gather/scatter routine which then calls LSP_EVAP,
!            LSP_FRMT,LSP_FORM.
!            GLUE is an extra level of control routine to avoid using
!            *IF DEF around calls to different LS_PPN versions, as per
!            S. Foreman's 22/8/94 proposal for plug compatibility.
!
!   A04_2D : Moisture distribution about mean gridbox supersaturation.
!            Long timestep and explicit short timestep ice fallouts.
!   A04_2E : MPP version of 2D scheme (loss of exact bit comparison).
!
!   Called by : LSPP_CTL1
!
!   Code description: Language FORTRAN 77 + extensions.
!
!   Programming standard: Unified Model Documentation Paper No 3,
!                         Version 6.
!
!   Author: Andrew Bushell     Reviewer: C. Wilson
!
!   Modification History from UM Version 4.0:
!    Version      Date
!
!   4.1  06/06/96  Add Sulphur Cycle tracers for wet scavenging in
!                  Version 2D of LS_PPN                    M. Woodage
!   4.4  01/07/97  Changes to use of arguments if 2A cloud scheme is
!                  chosen. MUST USE 1A cloud scheme with 2D, 2E ppn.
!   4.4  08/09/97  Add RHCRIT argument   D.Wilson
!    4.5  02/04/98   Add NH3 to argument list and pass to LS_PPN
!                     (For S Cycle)               M Woodage
!    4.5  12/03/98   Add aged soot to argument list and pass to
!                    LS_PPN                       Luke Robinson.
!   4.5  03/09/97  Added dummy 3D rain and snow variables. D.Wilson
!   4.5  01/05/98  Restrict murk aerosol calculations to aerosol
!                  levels=boundary levels. P.Clark
!   4.5  13/05/98  Altered call to glue routine (dummy variables).
!                                  S. Cusack
!
!   System components covered :
!
!   System task :
!
!   Documentation: UMDP No.
!
!  END -----------------------------------------------------------------
!
      SUBROUTINE GLUE_LSPP(
     & AK,BK,CF,DELTA_AK,DELTA_BK,PSTAR,TIMESTEP
     &,BLAND,CW_SEA,CW_LAND
     &,LS_GRID_QC,LS_BS
     &,RHCRIT,DUMMY, L_DUMMY
     &,Q_LEVELS,PFIELD
     &,POINTS,K1STPT,LSPICE_DIM1,LSPICE_DIM2,A_LEVELS,Q,QCF,QCL,T
     &,SO2,L_SULPC_SO2
     &,NH3,L_SULPC_NH3
     &,SO4_AIT,SO4_ACC,SO4_DIS
     &  ,AGED_SOOT                     !INOUT
     &  ,L_SOOT
     &,AEROSOL,L_MURK,LSRAIN,LSSNOW,LSRAIN3D,LSSNOW3D
     &,LSCAV_SO2,LSCAV_SO4AIT,LSCAV_SO4ACC,LSCAV_SO4DIS
     &,LSCAV_NH3
     &  ,LSCAV_AGEDSOOT                !INOUT
     &,ERROR
     & )
      IMPLICIT NONE
!-----------------------------------------------------------------------
! Some of the following variables are dummy, for use in other LS_PPN
! versions.
!-----------------------------------------------------------------------
! IN variables
!-----------------------------------------------------------------------
      INTEGER Q_LEVELS         ! No. of "wet" levels in the model.
!
      INTEGER POINTS           ! No. of gridpoints being processed.
!
      INTEGER PFIELD           ! No. of points in global field (at one
!                                vertical level).
      INTEGER K1STPT           ! First gridpoint processed within
!                                within complete field.
      INTEGER A_LEVELS         ! No.of aerosol levels used
      INTEGER LSPICE_DIM1       ! Dimension of dummy LSRAIN3D 
!                                 and LSSNOW3D.
      INTEGER LSPICE_DIM2       ! Dimension of dummy LSRAIN3D
!                                 and LSSNOW3D.
      REAL CF(PFIELD,Q_LEVELS) ! Cloud fraction.
!
      REAL PSTAR(PFIELD)       ! Surface pressure (Pa).
!
      REAL AK(Q_LEVELS)        ! Hybrid co-ordinate for centre of layer.
!
      REAL BK(Q_LEVELS)        ! Hybrid co-ordinate for centre of layer.
!
      REAL RHCRIT(Q_LEVELS)    ! Critical humidity for cloud formation.
!
      REAL DUMMY(PFIELD,Q_LEVELS)    ! Dummy variable in this version
!
      REAL DELTA_AK(Q_LEVELS)  ! Change of hybrid co-ord across layer.
!                                (Upper minus lower).
      REAL DELTA_BK(Q_LEVELS)  ! Change of hybrid co-ord across layer.
!                                (Upper minus lower).
      REAL LS_GRID_QC(PFIELD,Q_LEVELS)
! WARNING: Input contents of this argument are dependent upon ls cloud
! scheme chosen,
! 1A: Grid-box mean cloud condensate at processed levels (kg/ kg air).
! 2A: Liquid cloud fraction on model levels.
! This glue routine is only compatible with the 1A choice.
      REAL LS_BS(PFIELD,Q_LEVELS)
! WARNING: Input contents of this argument are dependent upon ls cloud
! scheme chosen,
! 1A: Maximum moisture fluctuation /6*sigma on levels (kg per kg air).
! 2A: Frozen cloud fraction on model levels.
! This glue routine is only compatible with the 1A choice.
      REAL TIMESTEP            ! Timestep (sec).
!
      REAL CW_SEA              ! Threshold cloud liquid water content
!                                over sea for conversion to ppn
!                                (kg water per m**3)
      REAL CW_LAND             ! Threshold cloud liquid water content
!                                over land for conversion to ppn
!                                (kg water per m**3)
      LOGICAL BLAND(PFIELD)    ! Land/sea mask
!
      LOGICAL L_MURK           ! Aerosol needs scavenging.
!
      LOGICAL L_SULPC_SO2   ! Sulphur Cycle on, tracers to be scavenged
     &       ,L_SULPC_NH3         ! indicates if NH3 present
!
      LOGICAL L_SOOT        ! soot on, tracers to be scavenged
!
      LOGICAL L_DUMMY ! Dummy variable
!
!-----------------------------------------------------------------------
! INOUT variables
!-----------------------------------------------------------------------
      REAL Q(PFIELD,Q_LEVELS)        ! Specific humidity
!                                      (kg water/kg air).
      REAL QCF(PFIELD,Q_LEVELS)      ! Cloud ice (kg per kg air).
!
      REAL QCL(PFIELD,Q_LEVELS)      ! Cloud liquid water (kg/ kg air).
!
      REAL T(PFIELD,Q_LEVELS)        ! Temperature (K).
!
      REAL AEROSOL(PFIELD,Q_LEVELS)  ! Aerosol (K).
      REAL                  ! Sulphur Cycle tracers for wet scavenging
     &    SO2(PFIELD,Q_LEVELS)
     &    ,NH3(PFIELD,Q_LEVELS)
     &   ,SO4_AIT(PFIELD,Q_LEVELS)
     &   ,SO4_ACC(PFIELD,Q_LEVELS)
     &   ,SO4_DIS(PFIELD,Q_LEVELS)
     &   ,AGED_SOOT(PFIELD,Q_LEVELS)   !INOUT
!
!
!-----------------------------------------------------------------------
! OUT variables
!-----------------------------------------------------------------------
      REAL LSRAIN(PFIELD)   ! Surface rainfall rate (kg per sq m per s).
!
      REAL LSSNOW(PFIELD)   ! Surface snowfall rate (kg per sq m per s).
      REAL LSRAIN3D(LSPICE_DIM1,LSPICE_DIM2) ! Dummy
      REAL LSSNOW3D(LSPICE_DIM1,LSPICE_DIM2) ! Dummy 
      REAL                  ! column totals of scavenged S Cycle tracers
     &    LSCAV_SO2(PFIELD)
     &    ,LSCAV_NH3(PFIELD)
!
     &   ,LSCAV_SO4AIT(PFIELD)
     &   ,LSCAV_SO4ACC(PFIELD)
     &   ,LSCAV_SO4DIS(PFIELD)
     &   ,LSCAV_AGEDSOOT(PFIELD)       !INOUT
!
!
      INTEGER ERROR         ! Return code - 0 if OK,
!                                           1 if bad arguments.
!
!    External subroutine called ----------------------------------------
      EXTERNAL  LS_PPN
!----------------------------------------------------------------------
!
      CALL LS_PPN(
! Input data not changed on output
     & AK,BK,CF
     &,DELTA_AK,DELTA_BK,PSTAR,TIMESTEP
     &,BLAND,CW_SEA,CW_LAND
     &,LS_GRID_QC,LS_BS
! Size and control data
     &,Q_LEVELS,PFIELD,POINTS,K1STPT,A_LEVELS
! Input data changed on output
     &,Q,QCF,QCL,T
     &,SO2,L_SULPC_SO2
     &,NH3,L_SULPC_NH3
     &,SO4_AIT,SO4_ACC,SO4_DIS
     &   ,AGED_SOOT                    !INOUT
     &   ,L_SOOT
     &,AEROSOL,L_MURK
! Output data
     &,LSRAIN,LSSNOW
     &,LSCAV_SO2,LSCAV_SO4AIT,LSCAV_SO4ACC,LSCAV_SO4DIS
     &,LSCAV_NH3
     &   ,LSCAV_AGEDSOOT               !INOUT
     &,ERROR)
!
      RETURN
      END
