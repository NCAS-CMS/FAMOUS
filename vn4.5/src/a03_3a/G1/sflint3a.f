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
CLL  SUBROUTINE SFL_INT------------------------------------------------
CLL
CLL  Purpose: To calculate interpolation coefficients for 10m winds
CLL           and 1.5m temperature/specific humidity diagnostics
CLL           using a generalisation of the method of Geleyn (1988).
CLL
CLL  Suitable for single column use (via *IF definition IBM).
CLL
CLL  Model            Modification history:
CLL version  Date
CLL
CLL   3.4  18/10/94   *DECK inserted into UM version 3.4. S Jackson
CLL
CLL   4.0  30/12/94   Modified calculation of 10m wind interpolation
CLL                   factor when effective roughness length used;
CLL                   10m wind assumed to lie on "local" profile at
CLL                   height z0m+10 metres above the surface.
CLL                                                    R.N.B.Smith
CLL   4.2   Oct. 96   T3E migration - *DEF CRAY removed
CLL                                     S J Swarbrick
CLL
CLL  Programming standard: Unified Model Documentation Paper No 4,
CLL                        Version 2, dated 18/1/90.
CLL
CLL  Logical component covered: Part of P243.
CLL
CLL  System Task:
CLL
CLL  External Documentation: UMDP No.24
CLL
CLL---------------------------------------------------------------------
C*L  Arguments :-
      SUBROUTINE SFL_INT (
     & P_POINTS,U_POINTS,RIB,Z1,Z0M,Z0M_EFF,Z0H,Z0F,CD_STD,CD,CH
     &,RESFT,WIND_PROFILE_FACTOR
     &,CDR10M,CHR1P5M,CER1P5M
     +,SU10,SV10,ST1P5,SQ1P5,LTIMER
     +)
      IMPLICIT NONE
      INTEGER
     + P_POINTS    ! IN No. of P-grid points to be processed.
     +,U_POINTS    ! IN No. of UV-grid points to be processed.
      REAL
     + RIB(P_POINTS)     ! IN Bulk Richardson number for
C                        !    lowest layer.
     +,Z1(P_POINTS)      ! IN Height of lowest model level (m).
     +,Z0M(P_POINTS)     ! IN Roughness length for momentum (m).
     +,Z0M_EFF(P_POINTS) ! IN Effective roughness length for
C                        !    momentum (m).
     +,Z0H(P_POINTS)     ! IN Roughness length for heat and
C                        !    moisture (m).
     +,Z0F(P_POINTS)     ! IN Roughness length in the free
C                        !    convective limit (m).
     &,CD_STD(P_POINTS)  ! IN Surface drag coefficient for shear stress
C                        !    only, i.e. without orographic part of drag
     &,CD(P_POINTS)      ! IN Effective surface drag coefficient,
C                        !    including orographic part of drag
     &,CH(P_POINTS)      ! IN Surface transfer coefficient for heat and
C                        !    moisture.
     +,RESFT(P_POINTS)   ! IN Total resistance factor for moisture
C                        !    transfer from the surface.
     &,WIND_PROFILE_FACTOR(P_POINTS)
C                        ! IN Factor for converting effective friction
C                        !    velocity to local one.
      LOGICAL
     + SU10                      ! IN 10m U-wind diagnostic flag
     +,SV10                      ! IN 10m V-wind diagnostic flag
     +,ST1P5                     ! IN screen temp diagnostic flag
     +,SQ1P5                     ! IN screen specific humidity
C                                !    diagnostic flag
     +,LTIMER                    ! IN TIMER diagnostics flag
C Output variables
C
      REAL
     + CDR10M(U_POINTS)  ! OUT interpolation coefficicent for 10m wind
     +,CHR1P5M(P_POINTS) ! OUT Interpolation coefficient for 1.5m
C                        !     temperature
     +,CER1P5M(P_POINTS) ! OUT Interpolation coefficient for 1.5m
C                        !     specific humidity
C*
C*L---------------------------------------------------------------------
      EXTERNAL TIMER
C*
C*L---------------------------------------------------------------------
C    Local and other symbolic constants :-
C*L-----------------COMDECK C_VKMAN-------------------------------------
C VKMAN IS VON KARMAN'S CONSTANT
      REAL VKMAN

      PARAMETER(VKMAN=0.4)
C*----------------------------------------------------------------------

      REAL Z1P5M,Z10M
      PARAMETER (
     + Z1P5M = 1.5  ! for diagnosis of screen values of temperature
C                   ! and humidity (assumed to be at 1.5m).
     +,Z10M = 10.0  ! for diagnosis of 10m winds.
     +)
C
C  Define local storage.
C
C  (a) Local work arrays.
C
      REAL
     + Z1E(P_POINTS)    ! Level 1 height + Z0M_EFF
     +,SQRTCD(P_POINTS) ! Square root of CD
C
C  (b) Scalars.
C
      REAL
     + Z1S              ! Level 1 height + Z0M_EFF - Z0H
     +,Z1P5ME           ! Z1P5M + Z0H
     +,Z10ME            ! Z10M + Z0M
     +,SQRTCD_K         ! Temporary storage in calc of 1.5 amd 10m diags
     +,Z_OVER_Z1        ! Temporary storage in calc of 1.5 amd 10m diags
     +,CDNZ             ! Neutral drag coef. for momentum @ 10m
     +,CDNZ1            ! Neutral drag coef. for momentum @ level1
     +,CHNZ             ! Neutral drag coef. for heat/moisture @ 1.5m
     +,CHNZ1            ! Neutral drag coef. for heat/moisture @ level1
     &,CDTEMP1          ! Workspace in calc of interpolation coeffs.
     &,CDTEMP2          ! Workspace in calc of interpolation coeffs.
     &,CDTEMP3          ! Workspace in calc of interpolation coeffs.

      INTEGER
     + I       ! Loop counter (horizontal field index).
C*
      IF (LTIMER) THEN
        CALL TIMER('SFL_INT   ',3)
      ENDIF
C
C-----------------------------------------------------------------------
CL 1. This routine uses a generalised formulation of Geleyn (1988)
CL   to interpolate to screen and 10m height using surface and bottom
CL   level values as well as roughness lengths for heat and momentum.
CL    Start of main loop. Set up variables needed later (eg height of
CL   envelope orography).
C-----------------------------------------------------------------------

      DO I=1,P_POINTS
        SQRTCD(I) = SQRT(CD(I))
        Z1E(I) = Z1(I) + Z0M_EFF(I)
      ENDDO

C-----------------------------------------------------------------------
CL 2. If selected calculate interpolation coefficient for 10m winds.
C-----------------------------------------------------------------------
C
      IF(SU10.OR.SV10) THEN
        DO I=1,P_POINTS
          Z10ME = Z10M + Z0M(I)
          CDNZ = LOG( Z10ME / Z0M(I) )
          CDNZ1 = LOG( Z1E(I) / Z0M_EFF(I) )
          SQRTCD_K = SQRTCD(I) / VKMAN
          Z_OVER_Z1 = Z10M  / Z1(I)

          IF (RIB(I).GE.0.0) THEN
C
C Stable case
C
            CDR10M(I) = Z_OVER_Z1 + SQRTCD_K *
     +                  (CDNZ - Z_OVER_Z1*CDNZ1)
          ELSE
C
C Unstable Case
C
            CDTEMP1 = EXP( 1.0 / SQRTCD_K )
            CDTEMP2 = Z1E(I) / Z0M_EFF(I)
            CDTEMP3 = LOG( ( ( Z1E(I) - Z0M(I) ) * CDTEMP1 -
     &                       ( Z0M_EFF(I) - Z0M(I) ) * CDTEMP2 ) /
     &                     ( ( Z1E(I) - Z10ME ) * CDTEMP1 -
     &                       ( Z0M_EFF(I) - Z10ME ) * CDTEMP2 ) )
            CDR10M(I) = SQRTCD_K * ( CDNZ + CDTEMP3 )
          ENDIF
          CDR10M(I) = CDR10M(I) * WIND_PROFILE_FACTOR(I)
        ENDDO
      ENDIF
C
C-----------------------------------------------------------------------
CL 3. If selected calculate interpolation coefficient for 1.5m screen
CL      temperature and specific humidity.
C-----------------------------------------------------------------------
C
      IF (ST1P5.OR.SQ1P5) THEN
        DO I=1,P_POINTS
C
C variables to be used later
          Z1S = Z1E(I) - Z0H(I)
          Z1P5ME = Z1P5M + Z0H(I)
          CHNZ = LOG( Z1P5ME / Z0H(I) )
          CHNZ1 = LOG( Z1E(I) / Z0H(I) )
          SQRTCD(I) = SQRT(CD_STD(I))
          SQRTCD_K =0.0
          IF (SQRTCD(I) .GT. 0.0) SQRTCD_K = CH(I) / (SQRTCD(I) * VKMAN)
          Z_OVER_Z1 =  Z1P5M  / Z1S

C
C Stable case
C
          IF (RIB(I).GE.0.0) THEN
            CHR1P5M(I) = Z_OVER_Z1 + SQRTCD_K *
     +                  (CHNZ - Z_OVER_Z1 * CHNZ1)
          ELSE
C
C Unstable Case
C
            CDTEMP1 = EXP(1.0 / SQRTCD_K)
            CDTEMP2 = ( Z1P5M * Z1E(I) +
     &                  Z0H(I) * ( Z1S - Z1P5M ) * CDTEMP1 ) /
     &                ( Z1S * Z1P5ME )
            CHR1P5M(I) = 1.0 - SQRTCD_K * LOG ( CDTEMP2 )
          ENDIF

C
            CER1P5M(I) = ( CHR1P5M(I) - 1.0 ) * RESFT(I)      ! P243.123
        ENDDO
      ENDIF
C
      IF (LTIMER) THEN
        CALL TIMER('SFL_INT ',4)
      ENDIF
      RETURN
      END
