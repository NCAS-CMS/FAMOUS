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
!    SUBROUTINE RAINOUT ------------------------------------
!
! Purpose: This subroutine removes dissolved tracer aerosol assuming
!          the amount in grid box is reduced in the same proportion as
!          the reduction in the total condensed water (liquid + ice)
!          due to ppn. (It is assumed that the concn of tracer
!          is the same in every droplet and ice particle)
!
!           Called by LSPP_CTL
!
! Current Code Owner: D L Roberts
!
! History:
! Version    Date    Comments
! -------    ----    --------
!   4.1    23/05/96  Original code               D L Roberts
!   4.3    17/03/96  Include layer thickness calcn for diagnostics.
!                    Disallow rainout if ppn does not reach surface.
!                                                       M Woodage
!   4.5   07/07/98  Generalise routine for any tracer (originally 
!                   written for dissolved sulphate aerosol and 
!                   called RAINOUT_SULPHATE)           M Woodage
!
! Code Description:
!   Language: FORTRAN77 + common extensions
!   This code is written to UMDP3 v6 programming standards
!
! System components covered:
!
! System task:
!
! Documentation: Not yet available
!
!-----------------------------------------------------------------
!
      SUBROUTINE RAINOUT(
     &       QCF             ! IN
     &      ,QCL             ! IN
     &      ,QPREVIOUS       ! IN
     &     ,LS_RAIN      ! IN
     &     ,LS_SNOW      ! IN
     &      ,TRACER          ! IN/OUT
     &      ,FIRST_POINT     ! IN
     &      ,LAST_POINT      ! IN
     &      ,P_FIELD         ! IN
     &      ,Q_LEVELS        ! IN
     &      ,RNOUT_TRACER    ! OUT
     &      ,AKDIFF,BKDIFF,PSTAR   !IN
     &      )
!
      IMPLICIT NONE
!
      INTEGER
     &        Q_LEVELS,         !IN, no. of wet levels
     &        P_FIELD,          !IN, no. of pts in full 2_D field
     &        FIRST_POINT,      !IN, first point in 2D domain
     &        LAST_POINT        !IN, last point in 2D domain
!
      REAL TRACER(P_FIELD,Q_LEVELS),  !INOUT mmr of dissolved tracer
     &     QCL(P_FIELD,Q_LEVELS),      !IN cloud liquid water (mmr)
     &     QCF(P_FIELD,Q_LEVELS),      !IN cloud frozen water (mmr)
     &    LS_RAIN(P_FIELD), ! IN Large-scale rain at the surface
     &    LS_SNOW(P_FIELD), ! IN Large-scale snow at the surface
     &     QPREVIOUS(P_FIELD,Q_LEVELS) !IN, total condensed water
!                                            before precipitation
      REAL AKDIFF(Q_LEVELS),     !IN, for layer thickness calcn
     &     BKDIFF(Q_LEVELS),     !IN,     "
     &     PSTAR(P_FIELD)        !IN,     "
      REAL RNOUT_TRACER(P_FIELD)      ! OUT tracer removed kg/m2/ts
!
!  Local variables
!
      INTEGER I,LEVEL         ! loop variables
!
      REAL
     &     QREMAIN,           ! total condensed water after precipn.
     &     FRACTION           ! fraction of water remaining
      REAL DELTA_TRACER       ! amount tracer removed from grid box
      REAL SURF_PRECIP(P_FIELD) ! total precipn at surface
!
      REAL RDZ               ! mass p.u.area of air in layer
!
C*L------------------COMDECK C_G----------------------------------------
C G IS MEAN ACCEL DUE TO GRAVITY AT EARTH'S SURFACE
      REAL G

      PARAMETER(G=9.80665)
C*----------------------------------------------------------------------
!
!     WANT TO RESTRICT CALCULATIONS TO POINTS WHERE THERE IS SOME
!     CONDENSED WATER, AND
!     TO POINTS WHERE SOME PRECIPITATION ACTUALLY REACHES
!     THE SURFACE.
!     THERE ARE THEN THREE CASES TO CONSIDER.
!     (A) CONDENSED WATER CONTENT HAS ACTUALLY INCREASED. IN THIS
!          CASE WE LEAVE THE DISSOLVED TRACER UNCHANGED.
!     (B) CONDENSED WATER CONTENT HAS DECREASED BUT IS NON-NEGATIVE.
!          IN THIS CASE WE REDUCE THE DISSOLVED TRACER IN THE
!         SAME RATIO.
!     (C) CONDENSED WATER CONTENT HAS DECREASED TO LESS THAN ZERO.
!          IN THIS CASE WE REMOVE ALL THE DISSOLVED TRACER.
!         (MAYBE THIS CASE SHOULD NOT OCCUR. HOWEVER IT COSTS
!         ALMOST NOTHING TO TRAP IT.)
!
!   Initialise RNOUT_TRACER to zero before doing rainout
        DO I = FIRST_POINT,LAST_POINT
          RNOUT_TRACER(I) = 0.0
        SURF_PRECIP(I) = LS_RAIN(I) + LS_SNOW(I)
        END DO
!
      DO LEVEL = 1,Q_LEVELS           !  loops over wet levels
        DO I = FIRST_POINT,LAST_POINT ! loop over domain on a level.
!
          IF( ( QPREVIOUS(I,LEVEL) .GT. 0.0 ) .AND.
     &        ( SURF_PRECIP(I) .GT. 0.0     ) )  THEN
!
            QREMAIN = QCF(I,LEVEL) + QCL(I,LEVEL)
            IF( QREMAIN .GT. QPREVIOUS(I,LEVEL) ) THEN
              FRACTION = 1.0
            ELSE IF ( QREMAIN .GE. 0.0 ) THEN
              FRACTION = QREMAIN/QPREVIOUS(I,LEVEL)
            ELSE
              FRACTION = 0.0
            ENDIF
!
! Calculate amount TRACER removed per grid box
            DELTA_TRACER = TRACER(I,LEVEL) * (1.0-FRACTION)
!
! Calculate mass of air per unit area in layer for conversion of
! tracer mixing ratio increment to mass p.u.area for STASH
      RDZ=(-AKDIFF(LEVEL)-BKDIFF(LEVEL)*PSTAR(I))/G
!
!  Accumulate removed TRACER for each level
            RNOUT_TRACER(I) = RNOUT_TRACER(I) + DELTA_TRACER*RDZ
!
! Decrement TRACER
            TRACER(I,LEVEL) = TRACER(I,LEVEL) - DELTA_TRACER
!
          ENDIF            ! END QPREVIOUS condition
!
        END DO             ! END OF I LOOP
      END DO               ! END OF LEVEL LOOP
!
      RETURN
      END
!
