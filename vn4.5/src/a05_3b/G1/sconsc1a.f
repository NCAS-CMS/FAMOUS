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
!     SUBROUTINE SCONSCV -----------------------------------------------
!
!    Purpose: Scavenge Sulphur Cycle tracers by convective precipitation
!             assuming loss rate = cost * conv ppn rate at surface
!             for all levels below cloud top.
!
!             Called from CONV_CT1 if Sulphur Cycle is on
!
!  Current code owner:  M. Woodage
!
! History:
! Version    Date     Comment
! -------    ----     -------
!  4.1      06/06/96  Original Code          M. Woodage
!  4.3  17/04/97    Tidy DEFS and code so that blank source is not
CLL                 produced (A. Brady)
!  4.4   30/09/97    Add logical to control below cloud scavenging. 
!                    Use conv cloud amount to adjust amount scavenged
!                    assuming CCA=0.05.                     (M Woodage)
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!
!  System component covered:
!
!  System task:
!
! Documentation:  Not yet available
!
!-------------------------------------------------------------------
!
      SUBROUTINE SCONSCV(TRACER,
     &                   TIMESTEP,
     &                   TR_LEVS,
     &                   NPFLD,
     &                   FIRST_POINT,LAST_POINT,
     &                   CCLDBASE,CCLDTOP,
     &               L_SCAV_BELOW_CLOUD,CCA,
     &                   RAINRATE,SNOWRATE,
     &                   K_RAIN,K_SNOW,
     &                   ACCU_SCAV_TR,
     &                   AKDIFF,BKDIFF,
     &                   PLEVS,P_STAR)
!
      IMPLICIT NONE
!
      INTEGER NPFLD                 ! IN, no. of pts in a 2_D field
      INTEGER CCLDBASE(NPFLD),      ! IN, convective cloud base
     &        CCLDTOP(NPFLD),       ! IN, convective cloud top
     &        PLEVS,                ! IN, no. of p_levels
     &        TR_LEVS,              ! IN, no. of tracer levels
     &        FIRST_POINT,          ! IN, first point for calcs to be do
     &        LAST_POINT            ! IN, last point for calcns to be do
!
!
      REAL TIMESTEP,                ! IN, timestep in secs
     &     K_RAIN,                  ! IN, scavenging rate coeff for rain
     &     K_SNOW                   ! IN, scavenging rate coeff for snow
!
      LOGICAL L_SCAV_BELOW_CLOUD    !IN, control for scavenging levels
!
      REAL  CCA(NPFLD)       ! IN, convective cloud amount (fraction)
      REAL RAINRATE(NPFLD),         ! IN conv rain rate at surface kg/m2
     &     SNOWRATE(NPFLD),         ! IN conv snow rate at surface kg/m2
     &     TRACER(NPFLD,TR_LEVS),   ! IN and OUT,   tracer
     &     ACCU_SCAV_TR(NPFLD),     ! OUT, column total of scvnged trcr
     &     AKDIFF(PLEVS),           ! IN, for layer thickness calcn:
     &     BKDIFF(PLEVS),           ! IN,       "
     &     P_STAR(NPFLD)            ! IN,       "
!
!  Local variables
!
      INTEGER I,K                    ! loop counters
      INTEGER START_LEVEL    ! lowest level for scavenging
!
      REAL TERMR,                    ! to assist calcn of scav rate
     &     TERMS,                    !
     &     RDZ,                      ! mass p.u.area of air in layer
     &     DELTA_TR                  ! tracer increment due to scvnging
!
      REAL
     &     TOTRATE,                  ! total scav rate
     &     INVTRAT                   ! 1/(1+TOTRATE)
!
C*L------------------COMDECK C_G----------------------------------------
C G IS MEAN ACCEL DUE TO GRAVITY AT EARTH'S SURFACE
      REAL G

      PARAMETER(G=9.80665)
C*----------------------------------------------------------------------
!
!
! Initialise ACCU_SCAV_TR array to zero before adding accumulations
      DO I=1,NPFLD
      ACCU_SCAV_TR(I)=0.0
      END DO
!
! Calculate total scavenging rate
!
      DO I=FIRST_POINT,LAST_POINT         ! leave out polar rows
!
      IF (CCLDTOP(I).GT.0 .AND. CCA(I).GT.0.0) THEN
!
! Set up START_LEVEL for scavenging
      IF (L_SCAV_BELOW_CLOUD) THEN
        START_LEVEL = 1
      ELSE
        START_LEVEL = CCLDBASE(I)
      END IF
!
        IF (RAINRATE(I).LE.0.0) THEN    ! check for negative ppn
          TERMR=0.0
         ELSE
          TERMR=K_RAIN*RAINRATE(I)
        ENDIF
!
        IF (SNOWRATE(I).LE.0.0) THEN
          TERMS=0.0
         ELSE
          TERMS=K_SNOW*SNOWRATE(I)
        ENDIF
!
! Calculate TOTRATE, *3600.0 because K_RAIN and K_SNOW are derived for
!  ppn rates in mm/hr, but model values are kg/m2/s (cf CON_SCAV)
        TOTRATE=(TERMR+TERMS)*3600.0*TIMESTEP
! Increase TOTRATE to obtain rate in cloudy part of grid box
! Assume CCA=0.05
       TOTRATE=TOTRATE / 0.05
        INVTRAT=1.0/(1.0+TOTRATE)
!
! Do scavenging, leaving out N and S polar rows
! Calculate amount of tracer scavenged and add to column total
!
      DO K = START_LEVEL,CCLDTOP(I)
!
! Calculate proportion of tracer mixing ratio scavenged out
       DELTA_TR=TRACER(I,K)*(1.0-INVTRAT)
! Reduce DELTA_TR to allow for non_cloudy part of grid box
      DELTA_TR = DELTA_TR * 0.05
!
! Calculate mass of air per unit area in layer for conversion of tracer
!  mixing ratio increment to mass p.u.a. for STASH
       RDZ=(-AKDIFF(K)-BKDIFF(K)*P_STAR(I))/G
!
! Increment column total mass p.u.a. of scavenged tracer
       ACCU_SCAV_TR(I)=ACCU_SCAV_TR(I)+DELTA_TR*RDZ
!
! Decrement tracer mixing ratio
          TRACER(I,K)=TRACER(I,K)-DELTA_TR
!
          END DO                      ! end K loop
        END IF
!
      END DO                          ! end of I loop
!
!
      RETURN
      END
!
