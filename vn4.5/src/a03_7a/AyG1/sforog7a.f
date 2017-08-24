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
!!!  SUBROUTINE SF_OROG------------------------------------------------
!!!
!!!  Purpose: Calculate roughness lengths, blending height and wind
!!!           profile factor
!!!
!!! SJ, RE        <- programmer of some or all of previous code changes
C Modification History:
C Version Date     Change
C  4.5    Jul. 98  Kill the IBM specific lines (JCThil)
!!!
!!!--------------------------------------------------------------------
      SUBROUTINE SF_OROG(
     & P_FIELD,LAND_FIELD,TILE_PTS,LAND_INDEX,TILE_INDEX,
     & L_Z0_OROG,LTIMER,
     & HO2R2_OROG,RIB,SIL_OROG,Z0M,Z1,
     & WIND_PROFILE_FACTOR,Z0M_EFF
     & )

      IMPLICIT NONE

      INTEGER
     & P_FIELD              ! IN Number of points on P-grid.
     &,LAND_FIELD           ! IN Number of land points.
     &,TILE_PTS             ! IN Number of tile points.
     &,LAND_INDEX(P_FIELD)  ! IN Index of land points.
     &,TILE_INDEX(LAND_FIELD)! IN Index of tile points.

      LOGICAL
     & LTIMER               ! IN .TRUE. for timer diagnostics
     &,L_Z0_OROG            ! IN .TRUE. to use orographic roughness.

      REAL
     & HO2R2_OROG(LAND_FIELD)!IN Peak to trough height of unresolved
!                            !   orography divided by 2SQRT(2) (m).
     &,RIB(LAND_FIELD)      ! IN Bulk Richardson number for lowest layer
     &,SIL_OROG(LAND_FIELD) ! IN Silhouette area of unresolved orography
!                           !    per unit horizontal area
     &,Z0M(LAND_FIELD)      ! IN Roughness length for momentum (m).
     &,Z1(P_FIELD)          ! IN Height of lowest atmospheric level (m).

!  Output variables.

      REAL
     & WIND_PROFILE_FACTOR(LAND_FIELD)
!                           ! OUT For transforming effective surface
!                           !     transfer coefficients to those
!                           !     excluding form drag.
     &,Z0M_EFF(LAND_FIELD)  ! OUT Effective roughness length for
!                                 momentum (m)

!  Work Variables

      INTEGER
     & I            ! Horizontal field index
     &,J            ! Tile field index
     &,L            ! Land field index

      REAL
     & H_BLEND_OROG ! Blending height
     &,RIB_FN       ! Interpolation coefficient for 0 < RIB < RI_CRIT
     &,ZETA1        ! Work space
     &,ZETA2        ! More work space
     &,ZETA3        ! Even more work space

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

!   Local parameters

      REAL H_BLEND_MIN,H_BLEND_MAX
      PARAMETER (
     & H_BLEND_MIN=0.0       ! Minimun value of blending height
     &,H_BLEND_MAX=1000.0    ! Maximum value of blending height
     & )

      EXTERNAL TIMER

      IF (LTIMER) THEN
        CALL TIMER('SF_OROG ',3)
      ENDIF

! Set blending height, effective roughness length and
! wind profile factor at land points.
      DO J=1,TILE_PTS
        L = TILE_INDEX(J)
        I = LAND_INDEX(L)

        WIND_PROFILE_FACTOR(L) = 1.0
        Z0M_EFF(L) = Z0M(L)

        IF (L_Z0_OROG) THEN

          ZETA1 = 0.5 * OROG_DRAG_PARAM * SIL_OROG(L)
          ZETA2 = LOG ( 1.0 + Z1(I)/Z0M(L) )
          ZETA3 = 1.0 / SQRT ( ZETA1/(VKMAN*VKMAN) + 1.0/(ZETA2*ZETA2) )
          ZETA2 = 1.0 / EXP(ZETA3)
          H_BLEND_OROG = MAX ( Z1(I) / (1.0 - ZETA2) ,
     &                       HO2R2_OROG(L) * SQRT(2.0) )
          H_BLEND_OROG = MIN ( H_BLEND_MAX, H_BLEND_OROG )

! Apply simple stability correction to form drag if RIB is stable

          IF (SIL_OROG(L) .EQ. RMDI) THEN
             ZETA1 = 0.0
          ELSE
             RIB_FN =  ( 1.0 - RIB(L) / RI_CRIT )
             IF (RIB_FN.GT.1.0) RIB_FN = 1.0
             IF (RIB_FN.LT.0.0) RIB_FN = 0.0
             ZETA1 = 0.5 * OROG_DRAG_PARAM * SIL_OROG(L) * RIB_FN
          ENDIF

          Z0M_EFF(L) = H_BLEND_OROG / EXP ( VKMAN / SQRT ( ZETA1 +
     &                 (VKMAN / LOG ( H_BLEND_OROG / Z0M(L) ) ) *
     &                 (VKMAN / LOG ( H_BLEND_OROG / Z0M(L) ) ) ) )
          IF (RIB(L).GT.RI_CRIT) Z0M_EFF(L)=Z0M(L)

          WIND_PROFILE_FACTOR(L) = LOG( H_BLEND_OROG / Z0M_EFF(L) ) /
     &                             LOG( H_BLEND_OROG / Z0M(L) )

        ENDIF

      ENDDO

      IF (LTIMER) THEN
        CALL TIMER('SF_OROG ',4)
      ENDIF

      RETURN
      END

!!!  SUBROUTINE SF_OROG_GB --------------------------------------------
!!!
!!!  Purpose: Calculate effective roughness length and blending height
!!!
!!! SJ, RE        <- programmer of some or all of previous code changes
!!!
!!!--------------------------------------------------------------------
      SUBROUTINE SF_OROG_GB(
     & P_FIELD,P1,P_POINTS,LAND_FIELD,LAND1,LAND_PTS,LAND_INDEX,
     & LAND_MASK,L_Z0_OROG,HO2R2_OROG,RIB,SIL_OROG,Z0M,Z1,
     & H_BLEND_OROG,Z0M_EFF,LTIMER
     & )

      IMPLICIT NONE

      INTEGER
     & P_FIELD              ! IN Number of points on P-grid.
     &,P1                   ! IN First P-point to be processed.
     &,P_POINTS             ! IN Number of P-grid points to be processed
     &,LAND_FIELD           ! IN Number of land points.
     &,LAND1                ! IN First land point to be processed.
     &,LAND_PTS             ! IN Number of land points to be processed.
     &,LAND_INDEX(P_FIELD)  ! IN Index of land points.

      LOGICAL
     & LAND_MASK(P_FIELD)    ! IN .TRUE. for land; .FALSE. elsewhere.
     &,LTIMER               ! IN .TRUE. for timer diagnostics
     &,L_Z0_OROG            ! IN .TRUE. to use orographic roughness.


      REAL
     & HO2R2_OROG(LAND_FIELD)!IN Peak to trough height of unresolved
!                            !   orography divided by 2SQRT(2) (m).
     &,RIB(P_FIELD)         ! IN GBM Bulk Richardson number for lowest
!                           !    layer
     &,SIL_OROG(LAND_FIELD) ! IN Silhouette area of unresolved orography
!                           !    per unit horizontal area
     &,Z0M(LAND_FIELD)      ! IN GBM Roughness length for momentum (m).
     &,Z1(P_FIELD)          ! IN Height of lowest atmospheric level (m).

!  Output variables.

      REAL
     & H_BLEND_OROG(P_FIELD)! OUT Blending height
     &,Z0M_EFF(P_FIELD)     ! OUT Effective roughness length for
!                                 momentum (m)

!  Work Variables

      INTEGER
     & I            ! Horizontal field index
     &,L            ! Land field index

      REAL
     & RIB_FN       ! Interpolation coefficient for 0 < RIB < RI_CRIT
     &,ZETA1        ! Work space
     &,ZETA2        ! More work space
     &,ZETA3        ! Even more work space

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

!   Local parameters

      REAL H_BLEND_MIN,H_BLEND_MAX
      PARAMETER (
     & H_BLEND_MIN=0.0       ! Minimun value of blending height
     &,H_BLEND_MAX=1000.0    ! Maximum value of blending height
     & )

      EXTERNAL TIMER

      IF (LTIMER) THEN
        CALL TIMER('SF_OROG ',3)
      ENDIF

! Set blending height, effective roughness length and
! wind profile factor at land points.

      DO L = LAND1,LAND1+LAND_PTS-1
        I = LAND_INDEX(L)

        H_BLEND_OROG(I) = H_BLEND_MIN
        Z0M_EFF(I) = Z0M(L)

        IF (L_Z0_OROG) THEN

          ZETA1 = 0.5 * OROG_DRAG_PARAM * SIL_OROG(L)
          ZETA2 = LOG ( 1.0 + Z1(I)/Z0M(L) )
          ZETA3 = 1.0 / SQRT ( ZETA1/(VKMAN*VKMAN) + 1.0/(ZETA2*ZETA2) )
          ZETA2 = 1.0 / EXP(ZETA3)
          H_BLEND_OROG(I) = MAX ( Z1(I) / (1.0 - ZETA2) ,
     &                       HO2R2_OROG(L) * SQRT(2.0) )
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
     &                 (VKMAN / LOG ( H_BLEND_OROG(I) / Z0M(L) ) ) *
     &                 (VKMAN / LOG ( H_BLEND_OROG(I) / Z0M(L) ) ) ) )
          IF (RIB(I).GT.RI_CRIT) Z0M_EFF(I) = Z0M(L)

        ENDIF

      ENDDO

      IF (LTIMER) THEN
        CALL TIMER('SF_OROG ',4)
      ENDIF

      RETURN
      END
