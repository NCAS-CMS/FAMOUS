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
!+ Subroutine to increment a sum of fluxes.
!
! Method:
!       The arrays holding the summed fluxes are incremented
!       by a weighted sum of the variables suffixed with _INCR.
!       Arguments specify which arrays are to be incremented.
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
      SUBROUTINE AUGMENT_FLUX(N_PROFILE, N_LAYER, N_AUGMENT
     &   , ISOLIR, L_CLEAR
     &   , WEIGHT_INCR
     &   , FLUX_DIRECT, FLUX_TOTAL
     &   , FLUX_DIRECT_INCR, FLUX_TOTAL_INCR
     &   , FLUX_DIRECT_CLEAR, FLUX_TOTAL_CLEAR
     &   , FLUX_DIRECT_INCR_CLEAR, FLUX_TOTAL_INCR_CLEAR
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
!     INCLUDE COMDECKS
!     ------------------------------------------------------------------
!     COMDECK FOR TWO-STREAM RADIATION CODE.
!
!     MODULE TO SET FLAGS FOR DIFFERENT PORTIONS
!     OF THE SPECTRUM.
!
      INTEGER
     &     IP_SOLAR
!             SOLAR REGION
     &   , IP_INFRA_RED
!             INFRA-RED REGION
!
      PARAMETER(
     &     IP_SOLAR=1
     &   , IP_INFRA_RED=2
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
     &   , N_AUGMENT
!             LENGTH OF VECTOR TO AUGMENT
     &   , ISOLIR
!             SPECTRAL REGION
      LOGICAL   !, INTENT(IN)
     &     L_CLEAR
!             CLEAR FLUXES CALCULATED
      REAL  !, INTENT(IN)
     &     WEIGHT_INCR
!             WEIGHT TO APPLY TO FLUXES
     &   , FLUX_DIRECT_INCR(NPD_PROFILE, 0: NPD_LAYER)
!             DIRECT FLUX IN BAND
     &   , FLUX_TOTAL_INCR(NPD_PROFILE, 2*NPD_LAYER+2)
!             TOTAL FLUX IN BAND
     &   , FLUX_DIRECT_INCR_CLEAR(NPD_PROFILE, 0: NPD_LAYER)
!             CLEAR DIRECT FLUX IN BAND
     &   , FLUX_TOTAL_INCR_CLEAR(NPD_PROFILE, 2*NPD_LAYER+2)
!             CLEAR TOTAL FLUX IN BAND
      REAL  !, INTENT(INOUT)
     &     FLUX_DIRECT(NPD_PROFILE, 0: NPD_LAYER)
!             DIRECT FLUX
     &   , FLUX_TOTAL(NPD_PROFILE, 2*NPD_LAYER+2)
!             TOTAL FLUX
     &   , FLUX_DIRECT_CLEAR(NPD_PROFILE, 0: NPD_LAYER)
!             CLEAR DIRECT FLUX
     &   , FLUX_TOTAL_CLEAR(NPD_PROFILE, 2*NPD_LAYER+2)
!             CLEAR TOTAL FLUX
!
!     LOCAL ARGUMENTS.
      INTEGER
     &     I
!             LOOP VARIABLE
     &   , L
!             LOOP VARIABLE
!
!
!     INCREMENT THE ACTUAL FLUXES.

! There are four possible cases

      IF(ISOLIR.EQ.IP_SOLAR.AND.L_CLEAR) THEN

        DO I=0,N_LAYER
          DO L=1,N_PROFILE
            FLUX_DIRECT(L,I) = FLUX_DIRECT(L,I)+
     &        WEIGHT_INCR*FLUX_DIRECT_INCR(L,I)
            FLUX_DIRECT_CLEAR(L, I)=FLUX_DIRECT_CLEAR(L, I)+
     &        WEIGHT_INCR*FLUX_DIRECT_INCR_CLEAR(L, I)
          END DO
        END DO

        DO I=1, N_AUGMENT
          DO L=1, N_PROFILE
            FLUX_TOTAL(L, I)=FLUX_TOTAL(L, I)+
     &        WEIGHT_INCR*FLUX_TOTAL_INCR(L, I)
            FLUX_TOTAL_CLEAR(L, I)=FLUX_TOTAL_CLEAR(L, I)+
     &        WEIGHT_INCR*FLUX_TOTAL_INCR_CLEAR(L, I)
          ENDDO
        ENDDO

      ELSE IF(ISOLIR.EQ.IP_SOLAR.AND..NOT.L_CLEAR) THEN

        DO I=0, N_LAYER
          DO L=1, N_PROFILE
            FLUX_DIRECT(L, I)=FLUX_DIRECT(L, I)+
     &        WEIGHT_INCR*FLUX_DIRECT_INCR(L, I)
          ENDDO
        ENDDO

        DO I=1, N_AUGMENT
          DO L=1, N_PROFILE
            FLUX_TOTAL(L, I)=FLUX_TOTAL(L, I)+
     &        WEIGHT_INCR*FLUX_TOTAL_INCR(L, I)
          ENDDO
        ENDDO

      ELSE IF(ISOLIR.NE.IP_SOLAR.AND.L_CLEAR) THEN

        DO I=1, N_AUGMENT
          DO L=1, N_PROFILE
            FLUX_TOTAL(L, I)=FLUX_TOTAL(L, I)+
     &        WEIGHT_INCR*FLUX_TOTAL_INCR(L, I)
            FLUX_TOTAL_CLEAR(L, I)=FLUX_TOTAL_CLEAR(L, I)+
     &        WEIGHT_INCR*FLUX_TOTAL_INCR_CLEAR(L, I)
          ENDDO
        ENDDO

      ELSE IF(ISOLIR.NE.IP_SOLAR.AND..NOT.L_CLEAR) THEN

        DO I=1, N_AUGMENT
          DO L=1, N_PROFILE
            FLUX_TOTAL(L, I)=FLUX_TOTAL(L, I)+
     &        WEIGHT_INCR*FLUX_TOTAL_INCR(L, I)
          ENDDO
        ENDDO

      END IF

!
      RETURN
      END
