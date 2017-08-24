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
!
!!!  SUBROUTINE SF_LBEST-----------------------------------------------
!!!
!!!  Purpose: Calculate combination height for surface tiles, and
!!!           heat moisture and wind variables at that height
!!!
!!! Simon Jackson <- programmer of some or all of previous code changes
!!!
!!!  Model            Modification history:
!!! version  Date
!!!
!!!   4.3  7/2/97     New deck. S Jackson
!!!
!!!  Programming standard: Unified Model Documentation Paper No 4,
!!!                        Version ?, dated ?.
!!!
!!!  System component covered: P24.
!!!
!!!  Project task:
!!!
!!!  Documentation: UMDP 24.
!!!---------------------------------------------------------------------

!!  Arguaments --------------------------------------------------------

      SUBROUTINE SF_LBEST (
     & P_POINTS,P_FIELD,P1,H_BLEND_OROG,
     & QCL_1,QCF_1,QSTAR_GB,Q_1,TSTAR_GB,T_1,U_1,V_1,
     & Z0M_EFF,Z0H,Z0M,Z1_UV,Z1_TQ,H_BLEND,HEAT_BLEND_FACTOR,
     & Q_BLEND,QW_BLEND,T_BLEND,TL_BLEND,U_BLEND,V_BLEND,
     & WIND_BLEND_FACTOR,LTIMER
     & )

      IMPLICIT NONE

      INTEGER              !    Variables defining grid.
     & P_POINTS            ! IN Number of P-grid points to be processed.
     &,P_FIELD             ! IN Total number of p points.
     &,P1                  ! IN First p point to be processed.
     &,LAND_MASK(P_FIELD)  ! IN Land/sea mask

      LOGICAL
     & LTIMER

      REAL
     & H_BLEND_OROG(P_FIELD)! IN Blending height for effective
!                               roughness lengths
     &,QCL_1(P_FIELD)      ! IN Level 1 cloud water liquid water content
     &,QCF_1(P_FIELD)      ! IN Lev 1 cloud fraction
     &,QSTAR_GB(P_FIELD)   ! IN Mean surface QSTAR
     &,Q_1(P_FIELD)        ! IN Level 1 Q
     &,TSTAR_GB(P_FIELD)   ! IN Mean surface temperature
     &,T_1(P_FIELD)        ! IN Level 1 temperature
     &,U_1(P_FIELD)        ! IN U wind component for lowest
!                               atmospheric layer (m/s).  On P grid.
     &,V_1(P_FIELD)        ! IN V wind component for lowest
!                               atmospheric layer (m/s).  On P grid.
     &,Z0M_EFF(P_FIELD)    ! IN Effective roughness length for momentum
     &,Z0H(P_FIELD)        ! IN Roughness length for heat and moisture m
     &,Z0M(P_FIELD)        ! IN Roughness length for momentum (m).
     &,Z1_UV(P_FIELD)      ! IN Height of level 1 on UV levels
     &,Z1_TQ(P_FIELD)      ! IN Height of level 1 on TQ levels

! Output
      REAL
     & H_BLEND(P_FIELD)    ! OUT Blending height for times
     &,HEAT_BLEND_FACTOR(P_FIELD)
!                            OUT Heat Blending factor
     &,Q_BLEND(P_FIELD)    ! OUT Est of Q at blending height
!                                using neutral profile
     &,QW_BLEND(P_FIELD)   ! OUT Est of QW at blending height
!                                using neutral profile
     &,T_BLEND(P_FIELD)    ! OUT Est of temperature at blending height
!                                using neutral profile
     &,TL_BLEND(P_FIELD)   ! OUT Est of TL at blending height
!                                using neutral profile
     &,U_BLEND(P_FIELD)    ! OUT U component of wind at blending height
     &,V_BLEND(P_FIELD)    ! OUT U component of wind at blending height
     &,WIND_BLEND_FACTOR(P_FIELD)
!                            OUT Wind Blending factor
!  Work Variables

      INTEGER
     & I            ! Loop counter

      REAL
     & LAPSE        ! Atmospheric lapse rate in surface layer

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
C*L-----------------COMDECK C_VKMAN-------------------------------------
C VKMAN IS VON KARMAN'S CONSTANT
      REAL VKMAN

      PARAMETER(VKMAN=0.4)
C*----------------------------------------------------------------------

!!----------------------------------------------------------------------
!!!-----------COMDECK C_ROUGH FOR SUBROUTINE SF_EXCH----------
! Z0HSEA = roughness length for heat and moisture transport
!          over the sea (m).
! Z0MIZ  = roughness length for heat, moisture and momentum over
!          the Marginal Ice Zone (m).
! Z0SICE = roughness length for heat, moisture and momentum over
!          sea-ice (m).
      REAL Z0HSEA,Z0MIZ,Z0SICE

      PARAMETER(Z0HSEA = 4.0E-5,
     &          Z0MIZ  = 1.0E-1,
     &          Z0SICE = 3.0E-3)
!!----------------------------------------------------------------------
      REAL    RI_CRIT   ! Critical Richardson number, where Z0M_EFF=Z0M.
!                       ! Linear interpolation between RIB=0 and RI_CRIT
                                                                       
      REAL    OROG_DRAG_PARAM    ! Tunable parameter in calculation of
!                                ! Effective roughness length for 
!                                ! momentum                     
      PARAMETER(                                                
     & RI_CRIT=0.5,                                              
     & OROG_DRAG_PARAM=0.3)                                         
!*----------------------------------------------------------------------
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

      REAL LCRCP,LS,LSRCP
      PARAMETER (
     & LCRCP=LC/CP           ! Evaporation-to-dT conversion factor.
     &,LS=LF+LC              ! Latent heat of sublimation.
     &,LSRCP=LS/CP           ! Sublimation-to-dT conversion factor.
     & )

      EXTERNAL TIMER

      IF (LTIMER) THEN
        CALL TIMER('SF_LBEST',3)
      ENDIF

!#################################################################
! Start of code
!#################################################################

      DO I=P1,P1+P_POINTS-1

        H_BLEND(I) = Z1_UV(I) + Z0M_EFF(I)

        IF (H_BLEND(I) .NE. Z1_UV(I) + Z0M_EFF(I)) THEN   
          WIND_BLEND_FACTOR(I) = LOG ( H_BLEND(I) / Z0M_EFF(I) ) /
     &                LOG ( (Z1_UV(I) + Z0M_EFF(I)) / Z0M_EFF(I) )

        U_BLEND(I) = U_1(I) * WIND_BLEND_FACTOR(I)
        V_BLEND(I) = V_1(I) * WIND_BLEND_FACTOR(I)

        ELSE

          WIND_BLEND_FACTOR(I) = 1.0

          U_BLEND(I) = U_1(I)
          V_BLEND(I) = V_1(I)

        ENDIF

        H_BLEND(I) = Z1_TQ(I) + Z0M_EFF(I)

        IF (H_BLEND(I) .NE. Z1_TQ(I) + Z0M_EFF(I)) THEN
          HEAT_BLEND_FACTOR(I) = LOG ( H_BLEND(I) / Z0H(I) ) /
     &                LOG ( (Z1_TQ(I) + Z0M_EFF(I)) / Z0H(I) )

          T_BLEND(I) = TSTAR_GB(I) - G/CP * (H_BLEND(I) - Z0H(I))
     &               + (T_1(I) + G/CP * (Z1_TQ(I) + Z0M_EFF(I)- Z0H(I))
     &                    - TSTAR_GB(I) ) * HEAT_BLEND_FACTOR(I)

         Q_BLEND(I) = QSTAR_GB(I) + ( Q_1(I) - QSTAR_GB(I) ) *
     &                 HEAT_BLEND_FACTOR(I)
        ELSE

          HEAT_BLEND_FACTOR(I) = 1.0

          T_BLEND(I) = T_1(I)
          Q_BLEND(I) = Q_1(I)

        ENDIF

        TL_BLEND(I) = T_BLEND(I) - LCRCP*QCL_1(I) - LSRCP*QCF_1(I)
                                                           ! P243.9

        QW_BLEND(I) = Q_BLEND(I) + QCL_1(I) + QCF_1(I)     ! P243.10

      ENDDO


      IF (LTIMER) THEN
        CALL TIMER('SF_LBEST',4)
      ENDIF

      RETURN
      END


