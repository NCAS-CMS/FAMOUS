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
!+ Subroutine to find the optical depth and single scattering albedo.
!
! Method:
!       Depending on the treatment of scattering, the optical and
!       and single scattering albedo are determined from the
!       extinctions supplied.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.0             27-07-95                Original Code
!                                               (J. M. Edwards)
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE SINGLE_SCATTERING(I_SCATTER_METHOD_BAND
     &   , N_PROFILE, N_LAYER, I_TOP
     &   , D_MASS
     &   , K_GREY_TOT, K_EXT_SCAT, K_GAS_ABS
     &   , TAU, OMEGA
     &   , NPD_PROFILE, NPD_LAYER
     &   )
!
!
!
      IMPLICIT NONE
!
!
!     SIZES OF ARRAYS.
      INTEGER   !, INTENT(IN)
     &     NPD_PROFILE
!             MAXIMUM NUMBER OF PROFILES
     &   , NPD_LAYER
!             MAXIMUM NUMBER OF LAYERS
!
!     INCLUDE COMDECKS.
!     ------------------------------------------------------------------
!     COMDECK FOR TWO-STREAM RADIATION CODE.
!
!     MODULE FOR SETTING MACHINE PRECISION.
!
      REAL
     &     TOL_MACHINE
!             MACHINE TOLERANCE
     &   , SQRT_TOL_MACHINE
!             SQRT OF MACHINE TOLERANCE
!
!
!     THE PRECISION SHOULD BE ABOUT 2/2^(SIZE OF SIGNIFICAND)
!
!     THE IEEE-FORMAT USES 53 BITS FOR THE SIGNIFICAND
!     IN DOUBLE PRECISION
!
!     THE CRAY FORMAT USES 47 BITS IN SINGLE PRECISION.
!
      PARAMETER(
     &     TOL_MACHINE=1.42E-14
     &   , SQRT_TOL_MACHINE=1.19E-7
     &   )
!
!     ------------------------------------------------------------------
!     ------------------------------------------------------------------
!     COMDECK FOR TWO-STREAM RADIATION CODE.
!
!     MODULE TO SET THE METHODS OF TREATING SCATTERING.
!
      INTEGER
     &     IP_SCATTER_FULL
!             FULL TREATMENT OF SCATTERING
     &   , IP_NO_SCATTER_ABS
!             SCATTERING IGNORED COMPLETELY.
     &   , IP_NO_SCATTER_EXT
!             SCATTERING TREATED AS ABSORPTION
!
      PARAMETER(
     &     IP_SCATTER_FULL=1
     &   , IP_NO_SCATTER_ABS=2
     &   , IP_NO_SCATTER_EXT=3
     &   )
!
!     ------------------------------------------------------------------
!
!     DUMMY VARIABLES.
      INTEGER   !, INTENT(IN)
     &     I_SCATTER_METHOD_BAND
!             TREATMENT OF SCATTERING IN THIS BAND
!
!                       Atmospheric Properties
      INTEGER   !, INTENT(IN)
     &     N_PROFILE
!             NUMBER OF PROFILES
     &   , N_LAYER
!             NUMBER OF LAYERS
     &   , I_TOP
!             TOP LAYER TO CONSIDER
      REAL      !, INTENT(IN)
     &     D_MASS(NPD_PROFILE, NPD_LAYER)
!             MASS THICKNESS OF EACH LAYER
!
!                       Optical Propeties
      REAL      !, INTENT(IN)
     &     K_GREY_TOT(NPD_PROFILE, NPD_LAYER)
!             ABSORPTIVE EXTINCTION
     &   , K_EXT_SCAT(NPD_PROFILE, NPD_LAYER)
!             SCATTERING EXTINCTION
     &   , K_GAS_ABS(NPD_PROFILE, NPD_LAYER)
!             GASEOUS EXTINCTION
!
!                       Single Scattering Propeties
      REAL      !, INTENT(OUT)
     &     TAU(NPD_PROFILE, NPD_LAYER)
!             OPTICAL DEPTH
     &   , OMEGA(NPD_PROFILE, NPD_LAYER)
!             SINGLE SCATTERING ALBEDO
!
!
!
!     LOCAL VARIABLES.
      INTEGER
     &     L
!             LOOP VARIABLE
     &   , I
!             LOOP VARIABLE
!
!
!
!     THE MACHINE TOLERANCE IS ADDED TO THE DENOMINATOR IN THE
!     EXPRESSION FOR OMEGA TO PREVENT DIVISION BY ZERO: THIS IS
!     SIGNIFICANT ONLY IF THE TOTAL EXTINCTION IS SMALL, AND THUS
!     WILL NOT SENSIBLY AFFECT ANY PHYSICAL RESULTS.
!
      IF (I_SCATTER_METHOD_BAND.EQ.IP_SCATTER_FULL) THEN
!
         DO I=I_TOP, N_LAYER
            DO L=1, N_PROFILE
               TAU(L, I)=(K_GREY_TOT(L, I)+K_GAS_ABS(L, I))
     &            *D_MASS(L, I)
               OMEGA(L, I)=K_EXT_SCAT(L, I)
     &            /(K_GREY_TOT(L, I)+K_GAS_ABS(L, I)+TOL_MACHINE)
            ENDDO
         ENDDO
!
      ELSE IF (I_SCATTER_METHOD_BAND.EQ.IP_NO_SCATTER_ABS) THEN
!
!        THE SCATTERING EXTINCTION IS IGNORED COMPLETELY, SO
!        ONLY THE ABSORPTIVE CONTRIBUTIONS TO THE SINGLE
!        SCATTERING PROPETIES ARE INCLUDED.
!
         DO I=I_TOP, N_LAYER
            DO L=1, N_PROFILE
               TAU(L, I)=(K_GREY_TOT(L, I)+K_GAS_ABS(L, I)
     &            -K_EXT_SCAT(L, I))*D_MASS(L, I)
               OMEGA(L, I)=0.0
            ENDDO
         ENDDO
!
      ELSE IF (I_SCATTER_METHOD_BAND.EQ.IP_NO_SCATTER_EXT) THEN
!
!        THE SCATTERING EXTINCTION IS ADDED ON TO THE ABSORPTION.
!
         DO I=I_TOP, N_LAYER
            DO L=1, N_PROFILE
               TAU(L, I)=(K_GREY_TOT(L, I)+K_GAS_ABS(L, I))
     &            *D_MASS(L, I)
               OMEGA(L, I)=0.0E+00
            ENDDO
         ENDDO
!
      ENDIF
!
!
!
      RETURN
      END
