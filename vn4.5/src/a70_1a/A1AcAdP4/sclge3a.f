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
!+ Subroutine to set geometry of clouds.
!
! Method:
!       For use in multi-column mode arrays are set for each layer
!       pointing to profiles which have non-negligible clear or
!       cloudy fractions. The topmost cloudy layers are also
!       detected.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.0             27-07-95                Original Code
!                                               (J. M. Edwards)
!       4.2             Nov. 96   T3E migration: CALL WHENFGT,WHENFLE 
!                                  replaced by portable fortran code.
!                                                S.J.Swarbrick
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE SET_CLOUD_GEOMETRY(N_PROFILE, N_LAYER
     &   , L_GLOBAL_CLOUD_TOP, N_CLOUD_TOP_GLOBAL
     &   , W_CLOUD
     &   , N_CLOUD_PROFILE, I_CLOUD_PROFILE
     &   , N_CLOUD_TOP
     &   , N_FREE_PROFILE, I_FREE_PROFILE
     &   , NPD_PROFILE, NPD_LAYER
     &   )
!
!
!
      IMPLICIT NONE
!
!
!     SIZES OF DUMMY ARRAYS.
      INTEGER   !, INTENT(IN)
     &     NPD_PROFILE
!             MAXIMUM NUMBER OF PROFILES
     &   , NPD_LAYER
!             MAXIMUM NUMBER OF LAYERS
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
!     MODULE FOR SETTING MACHINE-DEPENDENT TOLERANCES.
!     (THE COMDECK PRMCH3A MUST ALWAYS BE INCLUDED BEFORE THIS COMDECK.)
!
      REAL
     &     TOL_DIV
!             TOLERANCE FOR DIVISION
     &   , TOL_TEST
!             TOLERANCE FOR TESTING EQUALITY
!
      PARAMETER(
     &     TOL_DIV=3.2E+01*TOL_MACHINE
     &   , TOL_TEST=1.6E+01*TOL_MACHINE
     &   )
!
!     ------------------------------------------------------------------
!
!     DUMMY ARGUMENTS.
      INTEGER   !, INTENT(IN)
     &     N_PROFILE
!             NUMBER OF PROFILES
     &   , N_LAYER
!             NUMBER OF LAYERS
      LOGICAL   !, INTENT(IN)
     &     L_GLOBAL_CLOUD_TOP
!             FLAG TO USE A GLOBAL VALUE FOR THE TOPS OF CLOUDS
      INTEGER   !, INTENT(IN)
     &     N_CLOUD_TOP_GLOBAL
!             GLOBAL TOPMOST CLOUDY LAYER
!
      INTEGER   !, INTENT(OUT)
     &     N_CLOUD_TOP
!             TOPMOST CLOUDY LAYER
     &   , N_FREE_PROFILE(NPD_LAYER)
!             NUMBER OF FREE PROFILES
     &   , I_FREE_PROFILE(NPD_PROFILE, NPD_LAYER)
!             CLOUD-FREE PROFILES
     &   , N_CLOUD_PROFILE(NPD_LAYER)
!             NUMBER OF CLOUDY PROFILES
     &   , I_CLOUD_PROFILE(NPD_PROFILE, NPD_LAYER)
!             PROFILES CONTAINING CLOUDS
      REAL      !, INTENT(IN)
     &     W_CLOUD(NPD_PROFILE, NPD_LAYER)
!             AMOUNTS OF CLOUD
!
!
!     LOCAL VARIABLES.
      INTEGER
     &     I,J,L                                                        
!             LOOP VARIABLE
!
!
!     SUBROUTINES CALLED:
!
!
!
      DO I=1, N_LAYER
!
         J                 =1
         N_CLOUD_PROFILE(I)=0
         DO L   =1,N_PROFILE
           IF (W_CLOUD(L,I).GT.TOL_TEST) THEN
             I_CLOUD_PROFILE(J,I)=L
             J                   =J+1
             N_CLOUD_PROFILE(  I)=N_CLOUD_PROFILE(I)+1
           END IF
         END DO
!
!
         J                =1
         N_FREE_PROFILE(I)=0
         DO L   =1,N_PROFILE
           IF (W_CLOUD(L,I).LE.TOL_TEST) THEN
             I_FREE_PROFILE(J,I)=L
             J                  =J+1
             N_FREE_PROFILE(  I)=N_FREE_PROFILE(I)+1
           END IF
         END DO
!
      ENDDO
!
      IF (L_GLOBAL_CLOUD_TOP) THEN
         N_CLOUD_TOP=N_CLOUD_TOP_GLOBAL
      ELSE
         N_CLOUD_TOP=1
         DO WHILE ( (N_CLOUD_TOP.LT.N_LAYER).AND.
     &              (N_CLOUD_PROFILE(N_CLOUD_TOP).EQ.0) )
            N_CLOUD_TOP=N_CLOUD_TOP+1
         ENDDO
      ENDIF
!
!
!
      RETURN
      END
