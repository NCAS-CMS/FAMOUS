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
!  SUBROUTINE KMKH---------------------------------------------------
!!!
!!!  Purpose: To set the turbulent mixing coefficients KM and KH
!!!           (Note: should be used after any vertical interpolation
!!!                  but before any horizontal interpolation.)
!!!
!!!
!!!  Model            Modification history
!!! version  date
!!!
!!!   4.4  09/09/97   New deck R.N.B.Smith
!!!  4.5    Jul. 98  Kill the IBM specific lines (JCThil)
!!!
!!!  Programming standard:
!!!
!!!  System component covered: Part of P243.
!!!
!!!  Documentation: UMDP No.24
!!!
!!!---------------------------------------------------------------------

! Arguments :-

      SUBROUTINE KMKH (
     & P_FIELD,P1,P_POINTS,BL_LEVELS
     &,RHOKM,RHO_KM,RHOKH
     &,RHOKMZ,RHOKHZ
     &,NTML,CUMULUS,RHOKM_TOP,RHOKH_TOP
     &,UNSTABLE,NTDSC,DSC
     &,LTIMER
     & )

      IMPLICIT NONE

      LOGICAL LTIMER             ! IN Flag for TIMER diagnostics

      INTEGER
     & P_FIELD                ! IN No. of P-grid points in whole field
     &,P1                     ! IN First P-grid point to be processed.
     &,P_POINTS               ! IN No. of P-grid points to be
!                                  processed.
     &,BL_LEVELS              ! IN No. of atmospheric levels for
!                                  which boundary layer fluxes are
!                                  calculated.
      LOGICAL
     & CUMULUS(P_FIELD)       ! IN Flag for cumulus the the b.l.
     &,UNSTABLE(P_FIELD)      ! IN Flag for unstable surface layer.

     &,DSC(P_FIELD)           ! IN Flag set if decoupled stratocumulus
!                             !    layer found.
      INTEGER
     & NTML(P_FIELD)          ! IN Number of model levels in the
!                                   turbulently mixed layer.
     &,NTDSC(P_FIELD)         ! IN Top level for turbulent mixing in
!                             !    cloud layer.                         
                                                                        
      REAL
     & RHOKMZ(P_FIELD,2:BL_LEVELS)
!                             ! IN Non-local turbulent mixing
!                                  coefficient for momentum.
     &,RHOKHZ(P_FIELD,2:BL_LEVELS)
!                             ! IN Non-local turbulent mixing
!                                  coefficient for heat and moisture.   
     &,RHOKM_TOP(P_FIELD,2:BL_LEVELS)
!                             ! IN Non-local top-down turbulent
!                             !    mixing coefficient for momentum.
     &,RHOKH_TOP(P_FIELD,2:BL_LEVELS)
!                             ! IN Non-local top-down turbulent
!                             !    mixing coefficient for heat
!                             !    and moisture. 
      REAL
     & RHOKM(P_FIELD,BL_LEVELS)
!                             ! INOUT Layer K-1 - to - layer K
!                                     turbulent mixing coefficient
!                                     for momentum.
     &,RHOKH(P_FIELD,BL_LEVELS)
!                             ! INOUT Layer K-1 - to - layer K
!                                     turbulent mixing coefficient
!                                     for heat and moisture.

      REAL
     & RHO_KM(P_FIELD,2:BL_LEVELS)
!                             ! OUT RHO * KM before horizontal
!                                   interpolation to UV-grid.


!!----------------------------------------------------------------------
!    External references :-
      EXTERNAL TIMER

!!----------------------------------------------------------------------

!  Define local storage.

      INTEGER
     & I             ! Loop counter (horizontal field index).
     &,K             ! Loop counter (vertical level index).


      IF (LTIMER) THEN
        CALL TIMER('KMKH    ',3)
      ENDIF

!-----------------------------------------------------------------------
!! Set KM and KH to be the maximum of the local and non-local values and
!! store RHO_KM on P-grid for output.
!-----------------------------------------------------------------------
      DO K=2,BL_LEVELS
        DO I=P1,P1+P_POINTS-1
          IF(CUMULUS(I) .AND. K.GE.NTML(I)) THEN
            RHOKH(I,K) = 0.0
            RHOKM(I,K) = 0.0
          ENDIF
          IF(UNSTABLE(I) .AND. K.GT.NTML(I)) THEN
            RHOKH(I,K) = 0.0
            RHOKM(I,K) = 0.0
          ENDIF
          IF(CUMULUS(I) .AND. K.GE.NTML(I) .AND. K.LT.NTML(I)+2) THEN
            RHOKHZ(I,K)=0.0
            RHOKMZ(I,K)=0.0
            RHOKH_TOP(I,K)=0.0
            RHOKM_TOP(I,K)=0.0
          ENDIF
          IF(DSC(I) .AND. K.GT.NTDSC(I)) THEN
            RHOKH(I,K) = 0.0
            RHOKM(I,K) = 0.0
          ENDIF
!
          RHOKH(I,K) = MAX( RHOKH(I,K) , RHOKHZ(I,K)+RHOKH_TOP(I,K) )
          RHOKM(I,K) = MAX( RHOKM(I,K) , RHOKMZ(I,K)+RHOKM_TOP(I,K) )
!
          RHO_KM(I,K) = RHOKM(I,K)
        ENDDO ! P_POINTS
      ENDDO ! BL_LEVELS

      IF (LTIMER) THEN
        CALL TIMER('KMKH    ',4)
      ENDIF

      RETURN
      END
