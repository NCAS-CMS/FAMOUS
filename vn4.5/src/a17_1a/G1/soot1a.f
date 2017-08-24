C *****************************COPYRIGHT******************************
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
C*LL  SUBROUTINES NEW2OLD and SOOTSCAV ---------------------------------
!LL  Purpose:
!LL        NEW2OLD converts a proportion of the fresh soot to an aged
!LL        variety. This conversion takes place as an exponential decay
!LL        with an e-folding time of 1.6 days.
!LL
!LL        SOOTSCAV causes a fraction of the aged soot to become
!LL        scavenged by cloud droplets, creating the third mode of
!LL        soot, soot in cloud water.
!LL
!LL  Modification History from Version 4.4
!LL     Version    Date
!LL       4.5      Jun 1998          New Deck       Luke Robinson
!LL
!LL  Programming standard: Unified Model Documentation Paper No 3
!LL  Code Description:
!LL  Language: FORTRAN77 + common extensions
!LL
!LL  Logical component covered:
!LL
!LL  Project task:
!LL
!LL  Documentation: Not yet available.
!LL
C*L  Arguments:---------------------------------------------------------
       SUBROUTINE NEW2OLD(NPTS,
     &                    FIRST_POINT,
     &                    LAST_POINT,
     &                    P_FIELD,
     &                    SootBefore,
     &                    SootAfter,
     &                    TimeStep)

! Converts fresh soot to aged using an exponential decay of former.

       INTEGER                     !,INTENT(IN)
     &         NPTS,               ! No. of points in 3D array.
     &         FIRST_POINT,
     &         LAST_POINT,
     &         P_FIELD             ! No. of points in each level.

       real TimeStep               !,INTENT(IN)

       real                        !,INTENT(INOUT)
     &         SootAfter(NPTS),    ! Aged soot.
     &         SootBefore(NPTS)    ! Fresh soot.

!     LOCAL VARIABLES.
!
       INTEGER i,j,k
       real Delta                  ! Amount of soot converted to aged.
       real rate                   ! Decay rate.

       parameter(rate=7.1E-6)

! This loop cycles through all points on all levels,
! but avoids the polar points.

       do k = 1, NPTS-P_FIELD + 1, P_FIELD
         do j = FIRST_POINT, LAST_POINT
             i = k + j - 1
             Delta         = (rate * TimeStep) * SootBefore(i)
             SootBefore(i) = SootBefore(i) - Delta
             SootAfter(i)  = SootAfter(i) + Delta
         enddo
       enddo

       return
       end
!===================================================================

      SUBROUTINE SOOTSCAV(Soot,
     &             SOOTINCLOUD,
     &             CLOUDF,
     &             NPNTS,QPNTS,NPFLD,TSTEP,
     &             QCL,QCF,
     &             FIRST_POINT,LAST_POINT
     &             )

! Performs nucleation scavenging, creating soot in cloud water
! from a proportion of the aged soot.

!-------------------------------------------------------------------

      INTEGER
     &        NPNTS,            !IN no. of pts in 3_D array on P_LEVS
     &        QPNTS,            !IN no. of pts in 3_D array on Q_LEVS
     &        NPFLD,            !IN no. of pts in 2_D field
     &        FIRST_POINT,      !IN first point for calcns to be done
     &        LAST_POINT        !IN last  point for calcns to be done
!
      REAL
     &     CLEARF(QPNTS),         ! IN clear air fraction (1-CLOUDF)
     &     CLOUDF(QPNTS),         ! IN cloud fraction (range 0 TO 1)
     &     DELTASOOTEVAP(NPNTS),  ! amount of soot released
!                                   when cloud water evaporates.
     &     DELTAST_NUCL(NPNTS),   ! amount of soot converted from aged
!                                   to SOOTINCLOUD by nucleation.
     &     EVAPTIME,       ! timescale for cloud droplets to evaporate
     &     QCTOTAL(QPNTS), ! total condensed water amount.(QCL+QCF)
     &     QCL(QPNTS),         !IN  cloud liquid water (mmr)
     &     QCF(QPNTS),         !IN  cloud frozen water (mmr)
!
     &     SOOT(NPNTS),        !INOUT mass mix rat of SOOT
     &     SOOTINCLOUD(NPNTS), !OUT mass mix rat of soot
!                               suspended in cloudwater.
     &     TSTEP                  !IN physics timestep

! Contains constants required for soot conversion and nucleation
! scavenging by cloud droplets. No diffusional scavenging is performed.
!
      REAL
     &     CLOUDTAU,       ! air parcel lifetime in cloud
     &     EVAPTAU,        ! timescale for suspended soot to evaporate
     &     NUCTAU,         ! timescale for accumulation mode particles
     &     THOLD           ! Cloud liquid water threshold
!                           for nucleation scavenging to occur.
!
      PARAMETER (
     &        CLOUDTAU = 1.08E4,           ! secs (=3 hours)
     &        EVAPTAU = 300.0,             ! secs  (=5 mins)
     &        NUCTAU = 30.0,               ! secs
     &        THOLD = 1.0E-8               ! kg/kg
     &          )
!-------------------------------------------------------------------
!
!--------------------------------------------------------------------
! Initialise DELTA increments to 0.0
!--------------------------------------------------------------------
!
      DO J=1,NPNTS-NPFLD+1,NPFLD      ! J loops over all model points
        DO I=FIRST_POINT,LAST_POINT   ! I loop omits N & S polar rows
          K=I+J-1
          DELTASOOTEVAP(K) = 0.0
          DELTAST_NUCL(K) = 0.0
        END DO
      END DO


!-------------------------------------------------------------------
! Calculate the total water content and the clear air fraction.
!-------------------------------------------------------------------

      DO J=1,QPNTS-NPFLD+1,NPFLD     ! J loops over all model points
        DO I=FIRST_POINT,LAST_POINT  ! I loop omits N & S polar rows
          K=I+J-1
          QCTOTAL(K) = QCL(K) + QCF(K)
          CLEARF(K)=1.0-CLOUDF(K)
        END DO
      END DO
!
!-------------------------------------------------------------------
!    Release of aerosol from evaporating cloud droplets:
!    if no condensed water (liquid + ice) in grid box, release
!    soot as aged soot.
!--------------------------------------------------------------------
!
      DO J=1,QPNTS-NPFLD+1,NPFLD      ! J loops over all model points
        DO I=FIRST_POINT,LAST_POINT   ! I loop omits N & S polar rows
          K=I+J-1
!     If cloud fraction less than 0.95, release some in clear  air.
          IF ( QCTOTAL(K) .LT. THOLD ) THEN
            DELTASOOTEVAP(K) = SOOTINCLOUD(K)
          ELSE IF  ( CLOUDF(K).LT.0.95 ) THEN
            EVAPTIME = EVAPTAU + 0.5*CLOUDTAU
            DELTASOOTEVAP(K) = ( 1.0 - EXP(-TSTEP/EVAPTIME) )
     &          *SOOTINCLOUD(K)
          ELSE
            DELTASOOTEVAP(K) = 0.0
          ENDIF
        END DO
      END DO
!
!     Also release any dissolved aerosol in a non-wet level,
!     where it should not be.
!
      IF (QPNTS.LT.NPNTS) THEN         ! ie. if dry points exist.
      DO J=QPNTS+1,NPNTS-NPFLD+1,NPFLD ! J loop omits wet points
        DO I=FIRST_POINT,LAST_POINT    ! I loop omits N & S polar rows
          K=I+J-1
          DELTASOOTEVAP(K) = SOOTINCLOUD(K)
        END DO
      END DO
      ENDIF
!
!
!-------------------------------------------------------------------
! Nucleation of aerosol forming SOOTINCLOUD (i.e. soot acting as CCN)
!-------------------------------------------------------------------
!
!    THIS CODE ASSUMES THAT THE PARAMETER NUCTAU, WHICH IS THE
!    TIMESCALE FOR NUCLEATION ONCE A PARTICLE ENTERS A CLOUD, IS
!    VERY SHORT COMPARED WITH CLOUDTAU.
!
      DO J=1,QPNTS-NPFLD+1,NPFLD     ! J loops over all model points
        DO I=FIRST_POINT,LAST_POINT  ! I loop omits N AND S polar rows
          K=I+J-1
          IF ((QCTOTAL(K) .GE. THOLD) .AND. (CLOUDF(K).GT.0.0)) THEN
            NUCTIME=NUCTAU + ( (CLEARF(K)*CLOUDTAU)/(2.0*CLOUDF(K)) )
            DELTAST_NUCL(K) = ( 1.0 - EXP(-TSTEP/NUCTIME) )*SOOT(K)
          ENDIF
        END DO
      END DO
!
!-------------------------------------------------------------------
!   UPDATE soot.
!--------------------------------------------------------------------
!
      DO J=1,QPNTS-NPFLD+1,NPFLD      ! J loops over all model points
        DO I=FIRST_POINT,LAST_POINT   ! I loop omits N & S polar rows
          K=I+J-1
            SOOT(K) = SOOT(K)
     &              + DELTASOOTEVAP(K)
     &              - DELTAST_NUCL(K)
            SOOTINCLOUD(K) = SOOTINCLOUD(K)
     &                     - DELTASOOTEVAP(K)
     &                     + DELTAST_NUCL(K)
        END DO
      END DO
!
!--------------------------------------------------------------------
      RETURN
      END
