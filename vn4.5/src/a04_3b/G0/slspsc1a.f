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
!
!   SUBROUTINE SLSPSCV ----------------------------------------------
!
! Purpose: This subroutine carries out wet scavenging of S Cycle tracers
!          assuming   loss rate = cost * ls ppn rate in layer above.
!
!          Called by LS_PPNC Version 2D if Sulphur Cycle is on.
!
! Current code owner: M. Woodage
!
! History:
! Version   Date     Comment
! -------   ----     ------
!   4.1   06/06/96   Original code            M. Woodage
!
! Code description:
!   Language: FORTRAN 77  + common extensions
!   This code is written to UMDP3 v6 programming standards.
!
! System Components covered:
!
! System task:
!
! Documentation:  Not yet available
!
!-------------------------------------------------------------------
!
      SUBROUTINE SLSPSCV(TRACER,        ! INOUT, S Cycle tracer
     &                   LSCAV_TR,      ! INOUT, accumulated scavngd tr
     &                   K_RAIN,        ! IN, scavenging coeff for rain
     &                   K_SNOW,        ! IN, scavenging coeff for snow
     &                   RDZ,           ! IN, AIR MASS P.U.AREA OF LAYER
     &                   TSTEP,         ! IN
     &                   NGPTS,         ! IN, NO. OF GATHERED POINTS
     &                   RAINRATE,      ! IN, LS RAIN IN LAYER ABOVE
     &                   SNOWRATE)      ! IN, LS SNOW IN LAYER ABOVE
!
      IMPLICIT NONE
!
      INTEGER NGPTS                ! IN, no points to be processed
!
      REAL TSTEP                   ! IN, timestep in secs
      REAL K_RAIN                  ! IN, scav rate const for rain
      REAL K_SNOW                  ! IN, scav rate const for snow
      REAL RDZ(NGPTS)              ! IN, -(DELTA_AK+DELTA_BK*PSTAR)/G
      REAL RAINRATE(NGPTS)         ! IN, rain rate mm/hr=kg/m2/s*3600
      REAL SNOWRATE(NGPTS)         ! IN, snow rate  "
      REAL TRACER(NGPTS)           ! INOUT, tracer to be scavenged
      REAL LSCAV_TR(NGPTS)         ! INOUT, accumulated scavenged tr
!
!  Local variables
      INTEGER I,J                   ! LOOP COUNTERS
!
      REAL TERMR                    ! local vars to assist calcn
      REAL TERMS                    !
      REAL DELTA_TR                 ! tracer increment due to scavnging
      REAL TOTRATE                  ! total LSPPN scavenging rate
      REAL INVTRAT                  ! 1/(1+TOTRATE)
!
! Calculate total scavenging rate array
!
      DO I=1,NGPTS
!
        IF (RAINRATE(I).LE.0.0) THEN      ! Check for negative ppn
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
! Calculate TOTRATE, *3600.0 because K_RAIN and K_SNOW values are
!  suitable for ppn rate in mm/hr, but model rates are in kg/m2/s
        TOTRATE=(TERMR+TERMS)*3600.0*TSTEP
        INVTRAT=1.0/(1.0+TOTRATE)
!
! Do scavenging
! Calculate proportion of tracer mixing ratio scavenged out
        DELTA_TR=TRACER(I)*(1.0-INVTRAT)
!
! Increment accumulated scavenged tracer in column, multiplying by  RDZ
! to convert mmr to mass per unit area.
        LSCAV_TR(I)=LSCAV_TR(I)+DELTA_TR*RDZ(I)
!
! Decrement tracer mixing ratio
        TRACER(I)=TRACER(I)-DELTA_TR
!
        END DO                             ! END OF I LOOP
!
      RETURN
      END
!
