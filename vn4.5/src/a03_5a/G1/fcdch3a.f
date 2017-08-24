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
CLL   SUBROUTINE FCDCH--------------------------------------------------
CLL
CLL  Purpose: Calculate bulk transfer coefficients at one or more
CLL           gridpoints, according to formulae derived by R N B Smith,
CLL           October 1989.
CLL
CLL  Model            Modification history:
CLL version  Date
CLL
CLL   3.4  18/10/94   *DECK inserted into UM version 3.4. S Jackson
CLL
CLL   4.0  30/12/94   Revised stability functions and calcs. for
CLL                   removing form drag effects from cH and cD(std);
CLL                   cD(std) used in definition of the Prandtl number.
CLL                                                      R.N.B.Smith
CLL  4.5  12/05/98  Optimize use of sqrt.  RBarnes@ecmwf.int
CLL
CLL  Programming standard: Unified Model Documentation Paper No 4,
CLL                        Version 2, dated 18/1/90.
CLL
CLL  System component covered: Part of P243.
CLL
CLL  Documentation: UM Documentation Paper No 24, section P243.
CLL                 See especially sub-section (iv).
CLL
C*L  Arguments:---------------------------------------------------------
      SUBROUTINE FCDCH(
     & RIB,Z0M,Z0H,Z0F,Z1,WIND_PROFILE_FACTOR,POINTS,CD,CH,CD_STD,LTIMER
     &)
      IMPLICIT NONE

      LOGICAL LTIMER

      INTEGER POINTS ! IN Number of gridpoints treated.

      REAL
     + RIB(POINTS)   ! IN Bulk Richardson number.
     +,Z0M(POINTS)   ! IN Roughness length for momentum transport (m).
     +,Z0H(POINTS)   ! IN Roughness length for heat and moisture (m).
     +,Z0F(POINTS)   ! IN Roughness length for free-convective heat and
     +               !    moisture transport (m).
     +,Z1(POINTS)    ! IN Height of centre of lowest model layer(m).
     &,WIND_PROFILE_FACTOR(POINTS)
C                    ! IN for adjusting the surface transfer
C                    !    coefficients to remove form drag effects.

      REAL
     & CD(POINTS)    ! OUT Surface drag coefficient including form drag.
     +,CH(POINTS)    ! OUT Bulk transfer coefficient for heat/moisture.
     &,CD_STD(POINTS)! OUT Surface drag coefficient excluding form drag.

C*L  Workspace usage----------------------------------------------------
C    No work areas are required.
C
C*----------------------------------------------------------------------
C*L  No external subprograms are called.

      EXTERNAL TIMER

C*----------------------------------------------------------------------
C  Common and local physical constants.
C*L-----------------COMDECK C_VKMAN-------------------------------------
C VKMAN IS VON KARMAN'S CONSTANT
      REAL VKMAN

      PARAMETER(VKMAN=0.4)
C*----------------------------------------------------------------------

      REAL ALPHAR,HETGEN,CZ,DM,THIRD
      PARAMETER (
     & ALPHAR=5.0,  ! Tunable parameter in FM and FH calculation.
     & HETGEN=0.0,  ! Tunable parameter to represent 'the degree of
C                   ! heterogeneity' of the surface; must be > or = 0.0
C                   ! and < or = 1.0
     + CZ=4.0,      ! Tunable parameter in unstable Fh, Fm calculations,
C                   ! equal to (3h)**-1.5 in the documentation.
     + DM=2.0,      ! Tunable parameter in unstable Fm calculation.
     + THIRD=1./3.  ! One third.
     +)
C
C  Define local variables (more or less in order of first appearance).
C
      INTEGER I       ! Loop counter; horizontal field index.
      REAL
     + KARMAN2        ! Square of von Karman's constant.
     +,ZETAM          ! See documentation for definition.
     +,ZETAH          ! See documentation for definition.
     &,CDN            ! CD for neutral conditions.
     &,CHN            ! CH for neutral conditions.
     &,CDN_STD        ! CD_STD for neutral conditions.
     &,PRANDTL        ! Prandtl number at neutrality.
     +,RFZ            ! Temporary in calculation of FM and FH.
     &,RIF            ! Flux Richardson number.
     &,AM             ! Temporary in calculation of FM and FH.
     &,AH             ! Temporary in calculation of FM and FH.
     +,BM             ! Temporary in calculation of FM and FH.
     +,BH             ! Temporary in calculation of FM and FH.
     &,BM_STD         ! Temporary in calculation of FM_STD.
     &,FM             ! Stability factor for CD.
     &,FH             ! Stability factor for CH.
     &,FM_STD         ! Stability factor for CD_STD.
     &,temp_sqrt      ! sqrt temporary value

      IF (LTIMER) THEN
        CALL TIMER('FCDCH   ',3)
      ENDIF

      KARMAN2=VKMAN*VKMAN
      DO 1 I=1,POINTS
C
C-----------------------------------------------------------------------
CL 1. Calculate neutral CD, CH.
C-----------------------------------------------------------------------
C
C  (A) Store ZETAM, ZETAH.
C
        ZETAM = LOG( (Z1(I) + Z0M(I)) / Z0M(I) )
        ZETAH = LOG( (Z1(I) + Z0M(I)) / Z0H(I) )
C
C  (B) Calculate neutral CD, CH.  Eqns P243.40, P243.41.
C
        CDN = KARMAN2 / ( ZETAM * ZETAM )
        CHN = KARMAN2 / ( ZETAH * ZETAM ) * WIND_PROFILE_FACTOR(I)
        CDN_STD = CDN * WIND_PROFILE_FACTOR(I) * WIND_PROFILE_FACTOR(I)

        PRANDTL = CDN_STD / CHN
C
C  (C) Calculate temporary quantities.
C
        AM = 2.0 * ALPHAR / PRANDTL
        AH = AM
C
C-----------------------------------------------------------------------
CL 2. Calculate functions Fm, Fh.
C-----------------------------------------------------------------------
        RFZ=0.0
        BM=0.0
        BH=0.0
        BM_STD=0.0

        RIF = RIB(I) / PRANDTL
C
C  Case 1: stable boundary layer (RIB > 0).
C
        IF (RIB(I) .GT. 0.0) THEN
          IF ( 1.0/RIF .GT. HETGEN*ALPHAR ) THEN
            FM = 1.0 - HETGEN * ALPHAR * RIF
            FM = ( FM * FM ) /
     &            ( 1.0 + 2.0 * (1.0-HETGEN) * ALPHAR * RIF )
            FH = FM
            FM_STD = FM
          ELSE
            FM = 0.0
            FH = 0.0
            FM_STD = 0.0
          ENDIF
        ELSE
C
C  Case 2: unstable boundary layer (RIB < or = 0).
C
C  (A) Store 1/Fz in RFZ.  Eqn P243.51, as approximated by P243.52.
C
          RFZ = CZ * SQRT ( Z1(I) / Z0F(I) )
C
C  (B) Store BM, BH and BM_STD.
C
          BM = DM * AM * CDN * RFZ
          BH = AH * CHN * RFZ
          BM_STD = DM * AM * CDN_STD * RFZ
C
C  (C) Finally calculate FM, FH and FM_STD.
C
          temp_sqrt = SQRT(-RIB(I))
          FM = 1.0 - AM * RIB(I) / ( 1.0 + BM * temp_sqrt )
          FH = 1.0 - AH * RIB(I) / ( 1.0 + BH * temp_sqrt )
          FM_STD = 1.0 - AM * RIB(I) / ( 1.0 + BM_STD * temp_sqrt )
        ENDIF
C
C-----------------------------------------------------------------------
CL 3. Calculate output coefficients.  Eqns P243.53, P243.54.
C-----------------------------------------------------------------------
C
        CD(I) = CDN * FM
        CH(I) = CHN * FH
        CD_STD(I) = CDN_STD * FM_STD

    1 CONTINUE

      IF (LTIMER) THEN
        CALL TIMER('FCDCH   ',4)
      ENDIF

      RETURN
      END
