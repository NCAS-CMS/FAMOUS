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
!+ Subroutine to calculate differences in source functions.
!
! Method:
!       Using the polynomial fit to the Planck function, values
!       of this function at the boundaries of layers are found
!       and differences across layers are determined. If the
!       Planckian is being taken to vary quadratically across
!       the layer second differences are found.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.0             27-07-95                Original Code
!                                               (J. M. Edwards)
!       4.1             14-03-95                Explicit formulation
!                                               for sea ice introduced.
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE DIFF_PLANCK_SOURCE(N_PROFILE, N_LAYER
     &   , N_DEG_FIT, THERMAL_COEFFICIENT
     &   , T_REF_PLANCK, T_LEVEL, T_GROUND
     &   , PLANCK_SOURCE, DIFF_PLANCK, PLANCK_GROUND
     &   , L_IR_SOURCE_QUAD, T, DIFF_PLANCK_2
     &   , N_FRAC_ICE_POINT, I_FRAC_ICE_POINT, ICE_FRACTION
     &   , PLANCK_FREEZE_SEA
     &   , NPD_PROFILE, NPD_LAYER, NPD_THERMAL_COEFF
     &   )
!
!
      IMPLICIT NONE
!
!
!     DUMMY ARRAY SIZES
      INTEGER   !, INTENT(IN)
     &     NPD_PROFILE
!             MAXIMUM NUMBER OF PROFILES
     &   , NPD_LAYER
!             MAXIMUM NUMBER OF LAYERS
     &   , NPD_THERMAL_COEFF
!             NUMBER OF THERMAL COEFFICIENTS
!
!     COMDECKS INCLUDED.
C*L------------------COMDECK C_O_DG_C-----------------------------------
C ZERODEGC IS CONVERSION BETWEEN DEGREES CELSIUS AND KELVIN
C TFS IS TEMPERATURE AT WHICH SEA WATER FREEZES
C TM IS TEMPERATURE AT WHICH FRESH WATER FREEZES AND ICE MELTS
      REAL ZERODEGC,TFS,TM

      PARAMETER(ZERODEGC=273.15,
     &          TFS=271.35,
     &          TM=273.15)
C*----------------------------------------------------------------------

!
!     DUMMY ARGUMENTS.
      INTEGER   !, INTENT(IN)
     &     N_PROFILE
!             NUMBER OF PROFILES
     &   , N_LAYER
!             NUMBER OF LAYERS
     &   , N_DEG_FIT
!             DEGREE OF FITTING FUNCTION
      LOGICAL   !, INTENT(IN)
     &     L_IR_SOURCE_QUAD
!             IR-SOURCE QUADRATIC
      REAL      !, INTENT(IN)
     &     THERMAL_COEFFICIENT(0: NPD_THERMAL_COEFF-1)
!             COEFFICIENTS OF FIT TO PLANCK FNC
     &   , T_REF_PLANCK
!             PLANCKIAN REFERENCE TEMPERATURE
     &   , T_LEVEL(NPD_PROFILE, 0: NPD_LAYER)
!             TEMPERATURES ON LEVELS
     &   , T_GROUND(NPD_PROFILE)
!             TEMPERATURES AT GROUND
      REAL      !, INTENT(OUT)
     &     PLANCK_SOURCE(NPD_PROFILE, 0: NPD_LAYER)
!             PLANCK FUNCTION ON LEVELS
     &   , DIFF_PLANCK(NPD_PROFILE, NPD_LAYER)
!             DIFFERENCES IN PLANCKIAN FNC
     &   , DIFF_PLANCK_2(NPD_PROFILE, NPD_LAYER)
!             TWICE 2ND DIFFERENCES IN PLANCKIAN FUNCTION
     &   , T(NPD_PROFILE, 0: NPD_LAYER)
!             TEMPERATURES AT CENTRES OF LAYERS
     &   , PLANCK_GROUND(NPD_PROFILE)
!             PLANCKIAN FUNCTION AT GROUND
!
!     ARGUMENTS RELATING TO SEA ICE.
      INTEGER   !, INTENT(IN)
     &     N_FRAC_ICE_POINT
!             NUMBER OF POINTS WITH FRACTIONAL ICE COVER
     &   , I_FRAC_ICE_POINT(NPD_PROFILE)
!             INDICES OF POINTS WITH FRACTIONAL ICE COVER
      REAL      !, INTENT(IN)
     &     ICE_FRACTION(NPD_PROFILE)
!             ICE FRACTION
!
      REAL      !, INTENT(OUT)
     &     PLANCK_FREEZE_SEA
!             PLANCK FUNCTION OVER FREEZING SEA
!
!
!     LOCAL VARIABLES.
      INTEGER
     &     I
!             LOOP VARIABLE
     &   , J
!             LOOP VARIABLE
     &   , L
!             LOOP VARIABLE
      REAL
     &     T_RATIO(NPD_PROFILE)
!             TEMPERATURE RATIO
!
!     VARIABLES FOR FRACTIONAL ICE COVER.
      INTEGER
     &     LG
!             GATHERED LOOP VARIABLE
      REAL
     &     T_ICE_G(NPD_PROFILE)
!             TEMPERATURE OF ICE SURFACE (GATHERED OVER SEA-ICE POINTS)
     &   , PLANCK_GROUND_G(NPD_PROFILE)
!             PLANCKIAN OF ICE SURFACE (GATHERED OVER SEA-ICE POINTS)


!
!
!
!     CALCULATE THE CHANGE IN THE THERMAL SOURCE FUNCTION
!     ACROSS EACH LAYER FOR THE INFRA-RED PART OF THE SPECTRUM.
      DO L=1, N_PROFILE
         T_RATIO(L)=T_LEVEL(L, 0)/T_REF_PLANCK
         PLANCK_SOURCE(L, 0)
     &      =THERMAL_COEFFICIENT(N_DEG_FIT)
      ENDDO
      DO J=N_DEG_FIT-1, 0, -1
         DO L=1, N_PROFILE
            PLANCK_SOURCE(L, 0)
     &         =PLANCK_SOURCE(L, 0)
     &         *T_RATIO(L)+THERMAL_COEFFICIENT(J)
         ENDDO
      ENDDO
      DO I=1, N_LAYER
         DO L=1, N_PROFILE
            T_RATIO(L)=T_LEVEL(L, I)/T_REF_PLANCK
            PLANCK_SOURCE(L, I)
     &         =THERMAL_COEFFICIENT(N_DEG_FIT)
         ENDDO
         DO J=N_DEG_FIT-1, 0, -1
            DO L=1, N_PROFILE
               PLANCK_SOURCE(L, I)
     &            =PLANCK_SOURCE(L, I)
     &            *T_RATIO(L)+THERMAL_COEFFICIENT(J)
            ENDDO
         ENDDO
         DO L=1, N_PROFILE
            DIFF_PLANCK(L, I)=PLANCK_SOURCE(L, I)
     &         -PLANCK_SOURCE(L, I-1)
         ENDDO
      ENDDO
!     CALCULATE THE SECOND DIFFERENCE IF REQUIRED.
      IF (L_IR_SOURCE_QUAD) THEN
         DO I=1, N_LAYER
!           USE THE SECOND DIFFERENCE FOR TEMPORARY STORAGE.
!           OF THE PLANCKIAN AT THE MIDDLE OF THE LAYER.
            DO L=1, N_PROFILE
               T_RATIO(L)=T(L, I)/T_REF_PLANCK
               DIFF_PLANCK_2(L, I)
     &            =THERMAL_COEFFICIENT(N_DEG_FIT)
            ENDDO
            DO J=N_DEG_FIT-1, 0, -1
               DO L=1, N_PROFILE
                  DIFF_PLANCK_2(L, I)
     &               =DIFF_PLANCK_2(L, I)
     &               *T_RATIO(L)+THERMAL_COEFFICIENT(J)
               ENDDO
            ENDDO
            DO L=1, N_PROFILE
               DIFF_PLANCK_2(L, I)=2.0E+00*(PLANCK_SOURCE(L, I)
     &            +PLANCK_SOURCE(L, I-1)-2.0E+00*DIFF_PLANCK_2(L, I))
            ENDDO
         ENDDO
      ENDIF
!
!     SOURCE AT THE SURFACE.
      DO L=1, N_PROFILE
         T_RATIO(L)=T_GROUND(L)/T_REF_PLANCK
         PLANCK_GROUND(L)=THERMAL_COEFFICIENT(N_DEG_FIT)
      ENDDO
      DO J=N_DEG_FIT-1, 0, -1
         DO L=1, N_PROFILE
            PLANCK_GROUND(L)=PLANCK_GROUND(L)*T_RATIO(L)
     &         +THERMAL_COEFFICIENT(J)
         ENDDO
      ENDDO
!
!     WHERE THERE IS FRACTIONAL SEA-ICE THE FORMULATION MUST BE
!     EXTENDED, BUT IT IS CONVENIENT TO CARRY OUT THE ABOVE OPERATION
!     AT ALL POINTS TO AVOID THE USE OF INDIRECT ADDRESSING.
!
!     CALCULATE THE SOURCE FUNCTION OVER OPEN SEA, ADOPTING THE MODEL'S
!     CONVENTION THAT THE TEMPAERTURE THERE IS FIXED.
      T_RATIO(1)=TFS/T_REF_PLANCK
      PLANCK_FREEZE_SEA=THERMAL_COEFFICIENT(N_DEG_FIT)
      DO J=N_DEG_FIT-1, 0, -1
         PLANCK_FREEZE_SEA=PLANCK_FREEZE_SEA*T_RATIO(1)
     &      +THERMAL_COEFFICIENT(J)
      ENDDO
!
!     DETERMINE THE TEMPERATURE OF THE ICE.
      DO L=1, N_FRAC_ICE_POINT
         LG=I_FRAC_ICE_POINT(L)
         T_ICE_G(L)=(T_GROUND(LG)
     &      -TFS*(1.0E+00-ICE_FRACTION(LG)))/ICE_FRACTION(LG)
      ENDDO
!
!     CALCULATE THE SOURCE FUNCTION AT POINTS WITH FRACTIONAL ICE.
      DO L=1, N_FRAC_ICE_POINT
         T_RATIO(L)=T_ICE_G(L)/T_REF_PLANCK
         PLANCK_GROUND_G(L)=THERMAL_COEFFICIENT(N_DEG_FIT)
      ENDDO
      DO J=N_DEG_FIT-1, 0, -1
         DO L=1, N_FRAC_ICE_POINT
            PLANCK_GROUND_G(L)=PLANCK_GROUND_G(L)*T_RATIO(L)
     &         +THERMAL_COEFFICIENT(J)
         ENDDO
      ENDDO
!
!     DETERMINE THE OVERALL PLANCKIAN FUNCTION OF THE SURFACE.
      DO L=1, N_FRAC_ICE_POINT
         LG=I_FRAC_ICE_POINT(L)
         PLANCK_GROUND(LG)=ICE_FRACTION(LG)*PLANCK_GROUND_G(L)
     &      +PLANCK_FREEZE_SEA*(1.0E+00-ICE_FRACTION(LG))
      ENDDO
!
!
!
      RETURN
      END
