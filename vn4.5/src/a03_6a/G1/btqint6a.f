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
! SUBROUTINE BTQ_INT
!!!  Purpose: To interpolate buoyancy parameters BT and BQ from full
!!!  levels to half levels
!!!
!!! ORIGINAL PROGRAMMER: J. James
!!! CURRENT CODE OWNER:  R.N.B. Smith
!!!
!!! HISTORY:
!!! DATE   VERSION   COMMENT
!!! ----   -------   -------
!!!
!!! 10/9/97  4.4     New Deck.  R.N.B.Smith
!!!
!!! CODE DESCRIPTION:
!!!   LANGUAGE: FORTRAN 77 + CRAY EXTENSIONS
!!!   THIS CODE IS WRITTEN TO UMDP3 PROGRAMMING STANDARDS.
!!!
!!! SYSTEM COMPONENT COVERED: P24
!!!---------------------------------------------------------------------
      SUBROUTINE BTQ_INT (
     & P_FIELD,P1,P_POINTS,BL_LEVELS
     &,DZL,RDZ,BQ,BT,BQ_CLD,BT_CLD,A_QS,A_DQSDT
     &,BQM,BTM,BQM_CLD,BTM_CLD,A_QSM,A_DQSDTM
     &,LTIMER
     &  )

      IMPLICIT NONE

! ARGUMENTS WITH INTENT IN. IE: INPUT VARIABLES.

      LOGICAL LTIMER          ! IN Flag for TIMER diagnostics

      INTEGER
     & P_FIELD                ! IN No. of P-grid points in whole field.
     &,P1                     ! IN First P-grid point to be processed.
     &,P_POINTS               ! IN No. of P-grid points to be processed.
     &,BL_LEVELS              ! IN No. of atmospheric levels for which
!                                boundary layer fluxes are calculated.
!                                Assumed ! <=30 for dimensioning GAMMA()
!                                in common deck C_GAMMA
      REAL
     & DZL(P_FIELD,BL_LEVELS) ! IN Layer depths (m).  DZL(,K) is the
!                                  distance from layer boundary K-1/2
!                                  to layer boundary K+1/2.  For K=1
!                                  the lower boundary is the surface.
     &,RDZ(P_FIELD,BL_LEVELS) ! IN Reciprocal of distance between
!                                  hybrid levels (m-1).  1/RDZ(,K) is
!                                  the vertical distance from level
!                                  K-1 to level K, except that for
!                                  K=1 it is just the height of the
!                                  lowest atmospheric level.
     &,BQ(P_FIELD,BL_LEVELS)  ! IN A buoyancy parameter for clear air
!                             !    on p,T,q-levels (full levels).
     &,BT(P_FIELD,BL_LEVELS)  ! IN A buoyancy parameter for clear air
!                             !    on p,T,q-levels (full levels).
     &,BQ_CLD(P_FIELD,BL_LEVELS)
!                             ! IN A buoyancy parameter for cloudy air
!                             !    on p,T,q-levels (full levels).
     &,BT_CLD(P_FIELD,BL_LEVELS)
!                             ! IN A buoyancy parameter for cloudy air
!                             !    on p,T,q-levels (full levels).
     &,A_QS(P_FIELD,BL_LEVELS)
!                             ! IN Saturated lapse rate factor
!                             !    on p,T,q-levels (full levels).
     &,A_DQSDT(P_FIELD,BL_LEVELS)
!                             ! IN Saturated lapse rate factor
!                             !    on p,T,q-levels (full levels).

      REAL  ! OUT arrays,
     & BQM(P_FIELD,BL_LEVELS) ! OUT A buoyancy parameter for clear air
!                             !    on intermediate levels (half levels):
!                             !    (*,K) elements are k+1/2 values.
     &,BTM(P_FIELD,BL_LEVELS) ! OUT A buoyancy parameter for clear air
!                             !    on intermediate levels (half levels):
!                             !    (*,K) elements are k+1/2 values.
     &,BQM_CLD(P_FIELD,BL_LEVELS)
!                             ! OUT A buoyancy parameter for cloudy air
!                             !    on intermediate levels (half levels):
!                             !    (*,K) elements are k+1/2 values.
     &,BTM_CLD(P_FIELD,BL_LEVELS)
!                             ! OUT A buoyancy parameter for cloudy air
!                             !    on intermediate levels (half levels):
!                             !    (*,K) elements are k+1/2 values.
     &,A_QSM(P_FIELD,BL_LEVELS)
!                             ! OUT Saturated lapse rate factor
!                             !    on intermediate levels (half levels):
!                             !    (*,K) elements are k+1/2 values.
     &,A_DQSDTM(P_FIELD,BL_LEVELS)
!                             ! OUT Saturated lapse rate factor
!                             !    on intermediate levels (half levels):
!                             !    (*,K) elements are k+1/2 values.


!-----------------------------------------------------------------------
!    External references :-
      EXTERNAL TIMER

!-----------------------------------------------------------------------
!    Local and other symbolic constants :-

!  Define local storage.

!  (b) Scalars.

      REAL
     & WK      ! Temporary in weighting factor.
     &,WKM1    ! Temporary in weighting factor.

      INTEGER
     & I       ! Loop counter (horizontal field index).
     &,K       ! Loop counter (vertical level index).

      IF (LTIMER) THEN
        CALL TIMER('BTQ_INT ',3)
      ENDIF
!-----------------------------------------------------------------------
!! 1.  Loop round levels.
!-----------------------------------------------------------------------

      DO K=2,BL_LEVELS
        DO I=P1,P1+P_POINTS-1

!-----------------------------------------------------------------------
!! 1.1 Calculate buoyancy parameters at half levels,
!!     i.e. at level K-1/2, if current level is level K.
!-----------------------------------------------------------------------

          WKM1 = 0.5 * DZL(I,K-1) * RDZ(I,K)         ! P243.C5 (2nd eqn)
          WK = 0.5 * DZL(I,K) * RDZ(I,K)             ! P243.C5 (1st eqn)

          BTM(I,K-1) = WKM1*BT(I,K) + WK*BT(I,K-1)
          BQM(I,K-1) = WKM1*BQ(I,K) + WK*BQ(I,K-1)
          BTM_CLD(I,K-1) = WKM1*BT_CLD(I,K) + WK*BT_CLD(I,K-1)
          BQM_CLD(I,K-1) = WKM1*BQ_CLD(I,K) + WK*BQ_CLD(I,K-1)
          A_QSM(I,K-1) = WKM1*A_QS(I,K) + WK*A_QS(I,K-1)
          A_DQSDTM(I,K-1) = WKM1*A_DQSDT(I,K) + WK*A_DQSDT(I,K-1)

        ENDDO ! p_points
      ENDDO ! bl_levels

      IF (LTIMER) THEN
        CALL TIMER('BTQ_INT ',4)
      ENDIF

      RETURN
      END
