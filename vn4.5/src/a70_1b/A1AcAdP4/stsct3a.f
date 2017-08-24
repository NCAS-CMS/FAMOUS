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
!+ Subroutine to set the scattering flag in this a band.
!
! Method:
!       Straightforward.
!
! Current Owner of Code: J. M. Edwards
!
! History:
!       Version         Date                    Comment
!       4.0             27-07-95                Original Code
!                                               (J. M. Edwards)
!                                               (J. M. Edwards)
!
! Description of Code:
!   FORTRAN 77 with extensions listed in documentation.
!
!- ---------------------------------------------------------------------
      SUBROUTINE SET_SCATTERING(I_SCATTER_METHOD
     &   , L_SWITCH_SCATTER
     &   , I_SCATTER_METHOD_BAND
     &   )
!
!
!
      IMPLICIT NONE
!
!
!
!     INCLUDE COMDECKS
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
!
!     DUMMY ARGUMENTS.
      INTEGER   !, INTENT(IN)
     &     I_SCATTER_METHOD
!             METHOD OF TREATING SCATTERING
      LOGICAL   !, INTENT(IN)
     &     L_SWITCH_SCATTER
!             SCATTERING SWITCH FOR THE BAND
      INTEGER   !, INTENT(OUT)
     &     I_SCATTER_METHOD_BAND
!             SCATTERING FLAG IN THIS BAND
!
!
!
      IF (I_SCATTER_METHOD.EQ.IP_SCATTER_FULL) THEN
!
!        PERFORM A FULL SCATTERING CALCULATION IN THIS BAND
         I_SCATTER_METHOD_BAND=IP_SCATTER_FULL
!
      ELSE IF (I_SCATTER_METHOD.EQ.IP_NO_SCATTER_ABS) THEN
!
!        SCATTERING EXTINCTION IS TO BE NEGLECTED IF SPECIFIED IN THIS
!        BAND. OTHERWISE A FULL SCATTERING CALCULATION IS REQUIRED.
         IF (L_SWITCH_SCATTER) THEN
            I_SCATTER_METHOD_BAND=IP_NO_SCATTER_ABS
         ELSE
            I_SCATTER_METHOD_BAND=IP_SCATTER_FULL
         ENDIF
!
      ELSE IF (I_SCATTER_METHOD.EQ.IP_NO_SCATTER_EXT) THEN
!
!        SCATTERING EXTINCTION IS TO BE TREATED AS ABSORPTION IF
!        SPECIFIED IN THIS BAND. OTHERWISE A FULL SCATTERING
!        CALCULATION IS REQUIRED.
         IF (L_SWITCH_SCATTER) THEN
            I_SCATTER_METHOD_BAND=IP_NO_SCATTER_ABS
         ELSE
            I_SCATTER_METHOD_BAND=IP_SCATTER_FULL
         ENDIF
!
      ENDIF
!
!
!
      RETURN
      END
