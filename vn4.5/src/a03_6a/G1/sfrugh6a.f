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
!!!  SUBROUTINE SF_ROUGH-----------------------------------------------
!!!
!!!  Purpose: Calculate roughness lengths, blending height and wind
!!!           profile factor
!!!
!!! SJ         <- programmerof some or all of previous code changes
C Modification History:
C Version Date     Change
C  4.5    Jul. 98  Kill the IBM specific lines (JCThil)
!!!
!!!--------------------------------------------------------------------

!  Arguaments --------------------------------------------------------

      SUBROUTINE SF_ROUGH (
     & P_FIELD,P_POINTS,LAND_FIELD,LAND_PTS,LAND_MASK,L_LAND,P1,LAND1,
     & LAND_INDEX,
     & L_Z0_OROG,Z1_UV,Z0MSEA,ICE_FRACT,
     & LYING_SNOW,Z0V,SIL_OROG,HO2R2_OROG,RIB,Z0M_EFF,Z0M,Z0H,
     & WIND_PROFILE_FACTOR,H_BLEND_OROG,MIZ_RUF,Z0HS,LTIMER
     & )

      IMPLICIT NONE

      INTEGER               !    Variables defining grid.
     & P_POINTS             ! IN Number of P-grid points to be processed
     &,P_FIELD              ! IN Number of points on P-grid.
     &,P1                   ! IN First P-point to be processed.
     &,LAND1                ! IN First land point to be processed.
     &,LAND_PTS             ! IN Number of land points to be processed.
     &,LAND_FIELD           ! IN Number of land points.

     &,LAND_INDEX(LAND_FIELD)!IN Index for compressed land point array;
!                               i'th element holds position in the FULL
!                               field of the ith land pt to be
!                               processed
      LOGICAL
     & LAND_MASK(P_FIELD)   ! IN .TRUE. for land; .FALSE. elsewhere. F60
     &,L_LAND               ! IN .TRUE. to calculate land points only
!                                 This saves time when tiling
     &,L_Z0_OROG            ! IN .TRUE. to use orographic roughness.
     &,LTIMER               ! IN .TRUE. for timer diagnostics

      REAL
     & HO2R2_OROG(LAND_FIELD)!IN Peak to trough height of unresolved
!                                orography devided by 2SQRT(2) (m).
     &,ICE_FRACT(P_FIELD)   ! IN Fraction of gridbox which is sea-ice.
     &,LYING_SNOW(P_FIELD)  ! IN Lying snow amount (kg per sq metre).
     &,RIB(P_FIELD)         ! IN Bulk Richardson number for lowest layer
     &,SIL_OROG(LAND_FIELD) ! IN Silhouette area of unresolved orography
!                                per unit horizontal area
     &,Z0V(P_FIELD)         ! IN Vegetative roughness length (m).  F6418
     &,Z1_UV(P_FIELD)       ! IN Height of lowest atmospheric level (m).

!  Modified (INOUT) variables.

      REAL
     & Z0MSEA(P_FIELD)      ! INOUT Sea-surface roughness length for
!                                  momentum (m).  F617.

!  Output variables.

      REAL
     & MIZ_RUF(P_FIELD)     ! OUT Surface roughness length for the
!                                 marginal ice zone at sea-ice points.
     &,H_BLEND_OROG(P_FIELD)!OUT Blending height
     &,WIND_PROFILE_FACTOR(P_FIELD)
!                           ! OUT For transforming effective surface
!                              transfer coefficients to those excluding
!                              form drag.
     &,Z0M_EFF(P_FIELD)     ! OUT Effective roughness length for
!                              momentum (m)
     &,Z0H(P_FIELD)         ! OUT Roughness length for heat and moisture
     &,Z0M(P_FIELD)         ! OUT Roughness length for momentum (m).
     &,Z0HS(P_FIELD)        ! OUT Roughness length for heat and moisture
!                              transport over sea.

!  Work Variables

      INTEGER
     & I            ! Loop counter
     &,L            ! Another loop counter - this time for land points

      REAL
     & RIB_FN       ! Interpolation coefficient for 0 < RIB < RI_CRIT
     &,ZETA1        ! Work space
     &,Z0           ! yet more workspace

!   Common parameters

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

!   Local parameters

      REAL H_BLEND_MIN,H_BLEND_MAX
      PARAMETER (
     & H_BLEND_MIN=0.0        ! Minimun value of blending height
     &,H_BLEND_MAX=1000.0     ! Maximum value of blending height
     & )


      EXTERNAL TIMER

!-----------------------------------------------------------------------
!!  1 Fix roughness lengths for the various surface types.
!-----------------------------------------------------------------------
      IF (LTIMER) THEN
        CALL TIMER('SF_ROUGH',3)
      ENDIF

      IF(.NOT.L_LAND) THEN ! sea points as well as land points
        DO I = P1,P1+P_POINTS-1
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!  1.1 Liquid sea. Overwrite sea-ice and land in 3.1.2, 3.1.3.
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          Z0M(I) = Z0MSEA(I)                                  ! P243.B5
          Z0H(I) = Z0HSEA                                     !    "
          Z0M_EFF(I) = Z0MSEA(I)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!  1.2 Sea ice: Z0MIZ set on all points for input to FCDCH routine
!!        in CD_MIZ,CH_MIZ calculations. Similarily Z0HSEA
!!        CD_LEAD,CH_LEAD calculations. Z0SICE for CD,CH over sea-ice.
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          MIZ_RUF(I) = Z0MIZ
          Z0HS(I) = Z0HSEA
          IF (ICE_FRACT(I).GT.0.0 .AND. .NOT.LAND_MASK(I)) THEN
            Z0M(I) = Z0SICE                                   ! P243.B4
            Z0H(I) = Z0SICE                                   !    "
            Z0M_EFF(I) = Z0SICE
          ENDIF
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!  1.2a Specify blending height for all points. Set to minimum value
!!         so that LAMBDA_EFF = LAMBDA over the sea in KMKH.
!!         This avoids having to pass land-sea mask into KMKH.
!!         Also set the wind profile factor to the default value of 1.0
!!         for non-land points.
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

          H_BLEND_OROG(I) = H_BLEND_MIN
          WIND_PROFILE_FACTOR(I) = 1.0

        ENDDO
      ENDIF ! End of L_LAND block

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!  1.3 Land.  Reduce roughness if there is any snow lying.
!!        Eqns P243.B1, B2.
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
      DO L = LAND1,LAND1+LAND_PTS-1
        I = LAND_INDEX(L)

! Only reduce non-orographic roughness for land points without permanent
! snow.
        IF (LYING_SNOW(I) .LT. 5.0E3) THEN

          Z0 = Z0V(I) - 4.0E-4 * LYING_SNOW(I)
          ZETA1 = MIN( 5.0E-4 , Z0V(I) )
          Z0M(I) = MAX( ZETA1 , Z0 )
          Z0H(I) = MIN( 0.1*Z0V(I) , Z0M(I) )
        ELSE                 ! for permanent land-ice Z0V is appropriate
          Z0M(I) = Z0V(I)         ! P243.B1, based on P243.B2 (2nd case)
          Z0H(I) = Z0V(I)         !    "   ,   "    "    "    ( "    " )
        ENDIF

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!  1.4 Orographic roughness. Calculate Z0M_EFF in neutral case.
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        IF (L_Z0_OROG) THEN

! Set blending height, effective roughness length and
! wind profile factor at land points.

          H_BLEND_OROG(I) = MAX ( Z1_UV(I) + Z0M(I) ,
     &                            HO2R2_OROG(L) * SQRT(2.0) )
          H_BLEND_OROG(I) = MIN ( H_BLEND_MAX, H_BLEND_OROG(I) )


! Apply simple stability correction to form drag if RIB is stable

          IF (SIL_OROG(L) .EQ. RMDI) THEN
             ZETA1 = 0.0
          ELSE
             RIB_FN =  ( 1.0 - RIB(I) / RI_CRIT )
             IF (RIB_FN.GT.1.0) RIB_FN = 1.0
             IF (RIB_FN.LT.0.0) RIB_FN = 0.0
             ZETA1 = 0.5 * OROG_DRAG_PARAM * SIL_OROG(L) * RIB_FN
          ENDIF

          Z0M_EFF(I) = H_BLEND_OROG(I) / EXP ( VKMAN / SQRT ( ZETA1 +
     &                 (VKMAN / LOG ( H_BLEND_OROG(I) / Z0M(I) ) ) *
     &                 (VKMAN / LOG ( H_BLEND_OROG(I) / Z0M(I) ) ) ) )


          IF (RIB(I).GT.RI_CRIT) Z0M_EFF(I)=Z0M(I)

          WIND_PROFILE_FACTOR(I) = LOG( H_BLEND_OROG(I) / Z0M_EFF(I) ) /
     &                             LOG( H_BLEND_OROG(I) / Z0M(I) )

        ELSE
! Orographic roughness not represented so leave blending height and
! wind profile factor at their default values and set effective
! roughness length to its value based on land type.

          Z0M_EFF(I) = Z0M(I)
        ENDIF

      ENDDO

      IF (LTIMER) THEN
        CALL TIMER('SF_ROUGH',4)
      ENDIF

      RETURN
      END

