!+ Calculates total flux across a horizontal field
!
! Subroutine Interface
      SUBROUTINE FLUX_DIAG(FLUX,AREA,FIELD_SIZE,START_POINT,NPNTS,
     &                     CONV_FAC,TSTEP,FLUXFLD)
      IMPLICIT NONE
!
! Description:
! Part of the energy correction suite of routines:
! Sums scaled version of FLUX from START_POINT for NPNTS points
! NB: This routine has been rewritten and the interface changed
! slightly to make it suitable for use on MPP machines.
! Instead of just passing in the part of FLUX to be summed, the
! whole field is passed in, and the part to be summed is defined
! via START_POINT and NPNTS
!
! Method:
! The field is scaled and then summed via DO_SUMS which provides
! a suitable sum for the platform
!
! Current code owner : Paul Burton
!
! History
!  Model    Date      Modification history from model version 4.1
!  version
!    4.1    13/11/95  Modified interface to make more suitable for
!                     MPP use. P.Burton
!LL  4.4  05/09/97   Net flux now accumulated in prognostic field
!LL                  S.D. Mullerworth
!
! Subroutine Arguments:

      INTEGER FIELD_SIZE,        ! IN size of FLUX
     &        START_POINT,       ! IN local point to start sum at
     &        NPNTS              ! IN number of points to sum

      REAL    FLUX(FIELD_SIZE),  ! IN flux to be summed
     &        AREA(FIELD_SIZE),  ! IN area of grid box
     &        CONV_FAC,          ! IN conversion factor to translate
     &                           !    flux into energy units
     &        TSTEP              ! IN timestep
      REAL   
     &        FLUXFLD(FIELD_SIZE) ! INOUT sum of fluxes

! Parameters and COMMON
C*L------------------COMDECK C_A----------------------------------------
C A IS MEAN RADIUS OF EARTH
      REAL A

      PARAMETER(A=6371229.)
C*----------------------------------------------------------------------


! Local variabels
      INTEGER END_POINT,   ! local point to end sum at
     &        I            ! loop counter

      REAL
     &  WORK(FIELD_SIZE)        ! work space
     &  ,FACTOR                 ! Multiplication factor

      END_POINT=START_POINT+NPNTS-1
      FACTOR=(A*A)*(CONV_FAC*TSTEP)

      DO I=START_POINT,END_POINT
        FLUXFLD(I)=FLUXFLD(I)+(AREA(I)*FLUX(I))*FACTOR
      ENDDO

      RETURN
      END
