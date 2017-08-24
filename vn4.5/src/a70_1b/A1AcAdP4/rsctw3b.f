C ******************************COPYRIGHT******************************
C (c) CROWN COPYRIGHT 1998, METEOROLOGICAL OFFICE, All Rights Reserved.
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
!+ Subroutine to rescale optical depth and albedo.
!
! Method:
!       The standard rescaling formulae are applied.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.5             11-06-98                Optimised version
!                                               (P. Burton)
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE RESCALE_TAU_OMEGA(N_PROFILE
     &   , I_LAYER_FIRST, I_LAYER_LAST
     &   , TAU, OMEGA, FORWARD_SCATTER
     &   , NPD_PROFILE, NPD_LAYER
     &   )
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
!
!     DUMMY ARGUMENTS.
      INTEGER   !, INTENT(IN)
     &     N_PROFILE
!             NUMBER OF PROFILES
     &   , I_LAYER_FIRST
!             FIRST LAYER TO RESCALE
     &   , I_LAYER_LAST
!             FIRST LAYER TO RESCALE
      REAL      !, INTENT(IN)
     &     FORWARD_SCATTER(NPD_PROFILE, NPD_LAYER)
!             FORWARD SCATTERING
      REAL      !, INTENT(INOUT)
     &     TAU(NPD_PROFILE, NPD_LAYER)
!             OPTICAL DEPTH
     &   , OMEGA(NPD_PROFILE, NPD_LAYER)
!             ALBEDO OF SINGLE SCATTERING
!
!     LOCAL VARIABLES.
      INTEGER
     &     I
!             LOOP VARIABLE
     &   , L
!             LOOP VARIABLE

!     Temporary scalars for the rescheduled divide
      REAL
     &  TEMP1,TEMP2,TEMP3,TEMP4,TEMP5



      DO I=I_LAYER_FIRST, I_LAYER_LAST
         IF(N_PROFILE.GT.0) THEN
            TEMP1=1.0E+00 - OMEGA(1,I)*FORWARD_SCATTER(1,I)
            TEMP2=1.0E+00 - FORWARD_SCATTER(1,I)
            DO L=1,N_PROFILE-1
               TEMP5=TEMP2/TEMP1
               TAU(L,I)=TAU(L,I)*TEMP1
               TEMP3=1.0-OMEGA(L+1,I)*FORWARD_SCATTER(L+1,I)
               TEMP4=1.0-FORWARD_SCATTER(L+1,I)
               TEMP2=TEMP4
               TEMP1=TEMP3
               OMEGA(L,I) = OMEGA(L,I)*TEMP5
            ENDDO
            TEMP5=TEMP2/TEMP1
            TAU(N_PROFILE,I) = TAU(N_PROFILE,I)*TEMP1
            OMEGA(N_PROFILE,I)= OMEGA(N_PROFILE,I)*TEMP5
         END IF
      END DO

      RETURN
      END
