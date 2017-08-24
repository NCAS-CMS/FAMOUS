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
C*LL  SUBROUTINE SF_ROUGH-----------------------------------------------
CLL
CLL  Purpose: Calculate roughness lengths, blending height and wind
CLL           profile factor
CLL
C Modification History:
C Version Date     Change
C  4.5    Jul. 98  Kill the IBM specific lines (JCThil)
CLL
CLLEND------------------------------------------------------------------
C*
C*L  Arguaments --------------------------------------------------------
      SUBROUTINE SF_ROUGH (
     & P_POINTS,LAND_PTS,LAND_MASK,
     & P1,LAND_INDEX,
     & L_Z0_OROG,Z1,Z0MSEA,ICE_FRACT,
     & LYING_SNOW,Z0V,SIL_OROG,HO2R2_OROG,RIB,Z0M_EFF,Z0M,Z0H,
     & WIND_PROFILE_FACTOR,H_BLEND,CD_LEAD,Z0HS,Z0F,Z0FS,
     & LTIMER)

      IMPLICIT NONE
C
      INTEGER              !    Variables defining grid.
     & P_POINTS            ! IN Number of P-grid points to be processed.
     &,LAND_PTS            ! IN Number of land points to be processed.
C
     &,LAND_INDEX(LAND_PTS)! IN Index for compressed land point array;
C                          !    i'th element holds position in the FULL
C                          !    field of the ith land pt to be
C                          !    processed
     &,P1                  ! IN First P-point to be processed.
      LOGICAL
     & LAND_MASK(P_POINTS) ! IN .TRUE. for land; .FALSE. elsewhere. F60.
     &,L_Z0_OROG           ! IN .TRUE. to use orographic roughness.
C
      REAL
     & HO2R2_OROG(LAND_PTS)! IN Peak to trough height of unresolved
C                          !    orography devided by 2SQRT(2) (m).
     &,ICE_FRACT(P_POINTS) ! IN Fraction of gridbox which is sea-ice.
     &,LYING_SNOW(P_POINTS)! IN Lying snow amount (kg per sq metre).
     &,RIB(P_POINTS)       ! IN Bulk Richardson number for lowest layer.
     &,SIL_OROG(LAND_PTS)  ! IN Silhouette area of unresolved orography
C                          !    per unit horizontal area
     &,Z0V(P_POINTS)       ! IN Vegetative roughness length (m).  F6418.
     &,Z1(P_POINTS)        ! IN Height of lowest atmospheric level (m).
C
C  Modified (INOUT) variables.
C
      REAL
     & Z0MSEA(P_POINTS)    ! INOUT Sea-surface roughness length for
C                          !       momentum (m).  F617.
C
C  Output variables.
C
      REAL
     & CD_LEAD(P_POINTS)  ! Bulk transfer coefficient for momentum
C                         !  over sea-ice leads.Missing data over non
C                         !  sea-ice points.(Temporary store for Z0MIZ)
     &,H_BLEND(P_POINTS)   ! OUT Blending height
     &,WIND_PROFILE_FACTOR(P_POINTS)
C                          ! For transforming effective surface transfer
C                          ! coefficients to those excluding form drag.
     &,Z0M_EFF(P_POINTS)   ! OUT Effective roughness length for momentum
     &,Z0F(P_POINTS)       ! Roughness length for free-convective heat
C                          ! and moisture transport.
     &,Z0H(P_POINTS)       ! OUT Roughness length for heat and moisture
     &,Z0M(P_POINTS)       ! OUT Roughness length for momentum (m).
     &,Z0FS(P_POINTS)      ! Roughness length for free-convective heat
C                          ! and moisture transport over sea.
     &,Z0HS(P_POINTS)      ! Roughness length for heat and moisture
C                          ! transport over sea.

      LOGICAL LTIMER       ! Logical switch for TIMER diags
C
C  Work Variables
C
      INTEGER
     & I            ! Loop counter
     &,L            ! Another loop counter - this time for land points

      REAL
     & RIB_FN       ! Interpolation coefficient for 0 < RIB < RI_CRIT
     &,ZETA1        ! Work space
     &,ZETA2        ! More work space
     &,ZETA3        ! Even more work space
     &,Z0           ! yet more workspace
C
C   Common parameters

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

C   Local parameters

      REAL H_BLEND_MIN,H_BLEND_MAX
      PARAMETER (
     & H_BLEND_MIN=0.0       ! Minimun value of blending height
     &,H_BLEND_MAX=1000.0     ! Maximum value of blending height
     & )

      IF (LTIMER) THEN
        CALL TIMER('SFROUGH ',3)
      ENDIF
C-----------------------------------------------------------------------
CL  1 Fix roughness lengths for the various surface types.
C-----------------------------------------------------------------------
      DO I = 1,P_POINTS
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
CL  1.1 Liquid sea. Overwrite sea-ice and land in 3.1.2, 3.1.3.
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        Z0M(I) = Z0MSEA(I)                                     ! P243.B5
        Z0H(I) = Z0HSEA                                        !    "
        Z0M_EFF(I) = Z0MSEA(I)
        Z0F(I) = Z0FSEA                                        !    "
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
CL  1.2 Sea ice: Z0MIZ set on all points for input to FCDCH routine
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
CL  1.2a Specify blending height for all points. Set to minimum value
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
CL  1.3 Land.  Reduce roughness if there is any snow lying.
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
CL  1.4 Orographic roughness. Calculate Z0M_EFF in neutral case.
C- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        IF (L_Z0_OROG) THEN
C
C         ! Set blending height, effective roughness length and
C         ! wind profile factor at land points.
C

          ZETA1 = 0.5 * OROG_DRAG_PARAM * SIL_OROG(L)
          ZETA2 = LOG ( 1.0 + Z1(I)/Z0M(I) )
          ZETA3 = 1.0 / SQRT ( ZETA1/(VKMAN*VKMAN) + 1.0/(ZETA2*ZETA2) )
          ZETA2 = 1.0 / EXP(ZETA3)
          H_BLEND(I) = MAX ( Z1(I) / (1.0 - ZETA2) ,
     &                       HO2R2_OROG(L) * SQRT(2.0) )
          H_BLEND(I) = MIN ( H_BLEND_MAX, H_BLEND(I) )


! Apply simple stability correction to form drag if RIB is stable

          IF (SIL_OROG(L) .EQ. RMDI) THEN
             ZETA1 = 0.0
          ELSE
             RIB_FN =  ( 1.0 - RIB(I) / RI_CRIT )
             IF (RIB_FN.GT.1.0) RIB_FN = 1.0
             IF (RIB_FN.LT.0.0) RIB_FN = 0.0
             ZETA1 = 0.5 * OROG_DRAG_PARAM * SIL_OROG(L) * RIB_FN
          ENDIF


          Z0M_EFF(I) = H_BLEND(I) / EXP ( VKMAN / SQRT ( ZETA1 +
     &                 (VKMAN / LOG ( H_BLEND(I) / Z0M(I) ) ) *
     &                 (VKMAN / LOG ( H_BLEND(I) / Z0M(I) ) ) ) )


          IF (RIB(I).GT.RI_CRIT) Z0M_EFF(I)=Z0M(I)

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

      IF (LTIMER) THEN
        CALL TIMER('SFROUGH ',4)
      ENDIF
      RETURN
      END

