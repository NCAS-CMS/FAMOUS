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
!    SUBROUTINE INFILT------------------------------------------------
!
! Subroutine Interface:
      SUBROUTINE INFILT (NPNTS,SOIL_PTS,SOIL_INDEX,B,KS
     &,                  INFIL_FAC,STHF,INFIL
C LOGICAL LTIMER
     +,LTIMER
     +)

      IMPLICIT NONE
!
! Description:
!     Calculates the maximum infiltration rate at the soil
!     surface                                            (Cox, 6/95)
!
! Documentation : UM Documentation Paper 25
!
! Current Code Owner : David Gregory
!
! History:
! Version   Date     Comment
! -------   ----     -------
!  4.1                New deck.   Peter Cox
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!
! System component covered: P25
! System Task: P25
!


! Global variables:

! Subroutine arguments
!   Scalar arguments with intent(IN) :
      INTEGER
     & NPNTS                ! IN Number of gridpoints.
     &,SOIL_PTS             ! IN Number of soil points.

!   Array arguments with intent(IN) :
      INTEGER
     & SOIL_INDEX(NPNTS)    ! IN Array of soil points.

      REAL
     & B(NPNTS)             ! IN Clapp-Hornberger exponent.
     &,KS(NPNTS)            ! IN Saturated hydraulic conductivity
!                           !    (kg/m2/s).
     &,INFIL_FAC(NPNTS)     ! IN Infiltration enhancement factor.
     &,STHF(NPNTS)          ! IN Frozen soil moisture content of
!                           !    the layer as a fraction of
!                           !    saturation.
     &,SATHH(NPNTS)         ! IN Saturated soil water pressure (m).

C
      LOGICAL LTIMER        ! Logical switch for TIMER diags
!   Array arguments with intent(OUT) :
      REAL
     & INFIL(NPNTS)         ! OUT Maximum infiltration rate at the
!                           !     soil surface (kg/m2/s).

! Local scalars:
      INTEGER
     & I,J                  ! WORK Loop counter.

! Local arrays:
      REAL
     & STHUMAX(NPNTS)       ! WORK Maximum unfrozen soil moisture of
!                           !      the top layer as a fraction of
!                           !      saturation.

      IF (LTIMER) THEN
        CALL TIMER('INFILT  ',103)
      ENDIF

!-----------------------------------------------------------------------
! Initialise the infiltration rate (for land ice)
!-----------------------------------------------------------------------
      DO I=1,NPNTS
        INFIL(I)=KS(I)
      ENDDO

!-----------------------------------------------------------------------
! Diagnose the maximum unfrozen water
!-----------------------------------------------------------------------
!     DO J=1,SOIL_PTS
!       I=SOIL_INDEX(J)
!-----------------------------------------------------------------------
! Neglect the effect of frozen water on the infiltration rate - to avoid
! excessive runoff in the Boreal regions               (P.M.Cox 26/2/96)
!       STHUMAX(I)=1.0-STHF(I)       ! Old code
!-----------------------------------------------------------------------
!       STHUMAX(I)=1.0
!     ENDDO

!----------------------------------------------------------------------
! Calculate the hydraulic conductivity of the top layer.
!----------------------------------------------------------------------
!     CALL HYD_CON (NPNTS,SOIL_PTS,SOIL_INDEX,B,KS,STHUMAX,INFIL,
!    &              LTIMER)
!
!----------------------------------------------------------------------
! Include infiltration enhancement by vegetation
!----------------------------------------------------------------------
CDIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
      DO J=1,SOIL_PTS
        I=SOIL_INDEX(J)
        INFIL(I)=INFIL(I)*INFIL_FAC(I)
      ENDDO

      IF (LTIMER) THEN
        CALL TIMER('INFILT  ',104)
      ENDIF

      RETURN
      END
