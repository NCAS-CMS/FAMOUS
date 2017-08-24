C *****************************COPYRIGHT******************************
C (c) CROWN COPYRIGHT 1996, METEOROLOGICAL OFFICE, All Rights Reserved.
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
!+ Contains subroutines CALC_ENERGY_SUM_ARRAY and CALC_MASS_SUM_ARRAY
!  called from deck EMDIAG1B
!  These routines calculate the numbers that must eventually be
!  globally summed.
!
!
! Subroutine Interface:
      SUBROUTINE CALC_ENERGY_SUM_ARRAY(VAR,AREA,MASS,RS,CONV_FACT,
     &                                 FIELD_SIZE,START_POINT,END_POINT,
     &                                 SUM_ARRAY)
      IMPLICIT NONE
!
! Description:
! Part of the energy correction suite of routines:
! Calculates the energy of each grid box based on field VAR, between the
! points START_POINT and END_POINT and adds this in to the corresponding
! grid box on SUM_ARRAY
!
! Method:
! Look at the code!
!
! Current code owner : Paul Burton
!
! History
!  Model    Date      Modification history from model version 4.1
!  version
!    4.1    7/11/95   New DECK created to make EMDIAG suitable for
!                     MPP use. P.Burton
!
! Subroutine Arguments:

      INTEGER FIELD_SIZE,        ! IN size of filed VAR and SUM_ARRAY
     &        START_POINT,       ! IN point to start at
     &        END_POINT          ! IN point to end at

      REAL VAR(FIELD_SIZE),      ! IN energy variable
     &     AREA(FIELD_SIZE),     ! IN area of each grid box
     &     MASS(FIELD_SIZE),     ! IN mass of each grid box
     &     RS(FIELD_SIZE),       ! IN radius at this level
     &     SUM_ARRAY(FIELD_SIZE) ! INOUT contains array which will be
!                                !       globally summed

      REAL CONV_FACT             ! IN conversion factor to translate
!                                !    flux into energy units

! Parameters
C*L------------------COMDECK C_G----------------------------------------
C G IS MEAN ACCEL DUE TO GRAVITY AT EARTH'S SURFACE
      REAL G

      PARAMETER(G=9.80665)
C*----------------------------------------------------------------------

! Local variables

      INTEGER I ! loop variable

! loop from START_POINT to END_POINT and add energy at grid box
! to the SUM_ARRAY

       DO I=START_POINT,END_POINT
         SUM_ARRAY(I)=SUM_ARRAY(I)+
     &                RS(I)*RS(I)*AREA(I)*MASS(I)*VAR(I)*CONV_FACT/G
       ENDDO

       RETURN
       END


! Subroutine Interface:
      SUBROUTINE CALC_MASS_SUM_ARRAY(VAR,AREA,RS,
     &                               FIELD_SIZE,START_POINT,END_POINT,
     &                               SUM_ARRAY)
      IMPLICIT NONE
!
! Description:
! Part of the energy correction suite of routines:
! Calculates the mass of each grid box based on field VAR, between the
! points START_POINT and END_POINT and adds this in to the corresponding
! grid box on SUM_ARRAY
!
! Method:
! Look at the code!
!
! Current code owner : Paul Burton
!
! History
!  Model    Date      Modification history from model version 4.1
!  version
!    4.1    7/11/95   New DECK created to make EMDIAG suitable for
!                     MPP use. P.Burton
!
! Subroutine Arguments:

      INTEGER FIELD_SIZE,        ! IN size of filed VAR and SUM_ARRAY
     &        START_POINT,       ! IN point to start at
     &        END_POINT          ! IN point to end at

      REAL VAR(FIELD_SIZE),      ! IN energy variable
     &     AREA(FIELD_SIZE),     ! IN area of each grid box
     &     RS(FIELD_SIZE),       ! IN radius at this level
     &     SUM_ARRAY(FIELD_SIZE) ! INOUT contains array which will be
!                                !       globally summed

! Parameters
C*L------------------COMDECK C_G----------------------------------------
C G IS MEAN ACCEL DUE TO GRAVITY AT EARTH'S SURFACE
      REAL G

      PARAMETER(G=9.80665)
C*----------------------------------------------------------------------

! Local variables

      INTEGER I ! loop variable

! loop from START_POINT to END_POINT and add energy at grid box
! to the SUM_ARRAY

       DO I=START_POINT,END_POINT
         SUM_ARRAY(I)=SUM_ARRAY(I)+
     &                RS(I)*RS(I)*AREA(I)*VAR(I)/G
       ENDDO

       RETURN
       END
