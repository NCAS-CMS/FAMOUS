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
CLL  SUBROUTINE TROPIN------------------------------------------------
CLL
CLL  Purpose:  Finds the tropopause & returns index of where it is.
CLL      Not suitable for single-column model use, as it does
CLL             horizontal filling-in, & includes dynamical allocation.
CLL
CLL    Based on routine TROP, but
CLL      1) taking temperature rather than theta as input
CLL      2) rather than returning the temperature, pressure and height
CLL      of a continuously-varying tropopause found by extrapolating
CLL      lapse rates, it returns the index of the adjacent layer
CLL      boundary (N meaning the bottom of the Nth model layer)
CLL      3) "filling in" where no tropopause is found.
CLL
CLL        Author:  William Ingram
CLL
CLL  Model                Modification history:
CLL version  Date
CLL   4.2   23/9/96       New code, based on TROP
CLL
CLL    Note that the definition of "tropopause" matches TROP, which
CLL      is not quite the WMO definition, though the same critical
CLL      lapse rate is used (unless comdeck C_LAPSE is altered).
CLL      For details of the interpolation assumptions, see UMDP S1
CLL      section 3.2.2 or Swinbank & Wilson (1990: SRFRTN 48).
CLL      Any physical changes to one routine should be considered for
CLL                                         mirroring in the other.
CLLEND----------------------------------------------------------------
C
C*L  Arguments:-------------------------------------------------------
      SUBROUTINE TROPIN(PSTAR, T, P_EXNER_HALF, IT, L1, L2, ROW_LENGTH,
     &     P_LEVELS, MIN_TROP_LEVEL, MAX_TROP_LEVEL, AKH, BKH, WRAP)

      IMPLICIT NONE

      INTEGER!, INTENT (IN)
     &     L1,                          !   Number of points in arrays
     &     L2,                          !   Number of points to process
     &     ROW_LENGTH,                  !   Number of points per row
     &     P_LEVELS,                    !   Number of model levels
     &     MIN_TROP_LEVEL,
     &     MAX_TROP_LEVEL
C     ! Limits on where the tropopause may be deemed to be - between
C     !  the MIN_TROP_LEVELth and MAX_TROP_LEVELth layers (with the
C     !  convention used here for layer boundaries, the actual index
C     !  IT returned has MIN_TROP_LEVEL < IT =< MAX_TROP_LEVEL.)

      REAL!, INTENT(IN)
     &     PSTAR(L1),                   !  Surface pressure
     &     T(L1,P_LEVELS),              !  Temperature at layer centres
     &     P_EXNER_HALF(L1,P_LEVELS+1), !  Pexner at layer boundaries
     &     AKH(P_LEVELS+1),             !  Hybrid coordinate A & B
     &     BKH(P_LEVELS+1)              !   values for layer boundaries

      LOGICAL!, INTENT (IN)
     &     WRAP
C     ! Do the rows wrap round (so that the last point of a row is
C     !  geographically beside the first point of the same row) ?

      INTEGER!, INTENT (OUT)
     &     IT(L1)
C     ! Integer indexing the tropopause, taken to be @ a layer boundary
C     !   with the convention that N means the bottom of layer N.

C Workspace usage:-----------------------------------------------------

      REAL LAPSE_RATE(L2,
     & MIN_TROP_LEVEL+1:MAX_TROP_LEVEL+1-1/(1+P_LEVELS-MAX_TROP_LEVEL))
C     ! Lapse rate between layer centres (with top level MAX_TROP_LEVEL
C     !    if MAX_TROP_LEVEL = P_LEVELS, MAX_TROP_LEVEL+1 otherwise,
C     !    assuming MAX_TROP_LEVEL =< P_LEVELS as it should be).

      LOGICAL LTROP(L2)
C     !  Logical array which indicates whether we are still seeking a
C     !    tropopause (if not, it has already been found lower down)

C*---------------------------------------------------------------------
C Define local variables:----------------------------------------------

      INTEGER I, J,         ! Loopers over level & point
     &     POINT,           ! Point counter at ends of rows
     &     JP1,
C     !  J+1, except where this would cause out-of-bounds reference
     &     NNEIGH,          ! Number of well-defined tropopauses among
C     ! the 8 nearest neighbours of a point without one of its own
     &     FILLIN,
C     !  Used to fill in points without a clearly-defined tropopause
     &     DTI              ! Default tropopause index for points where
C     ! not even one nearest neighbour has a well-defined tropopause

      REAL PJP1, PJ, PJM1,  ! Pressures at half levels J+1/J/J-1
     &     P_EXNER_FULL_J,  ! Exner pressures at centres of current
     &     P_EXNER_FULL_JM1,!                      layer & one below
     &     DEL_EXNER_J,     ! Differences of Exner pressure across
     &     DEL_EXNER_JM1,   !                           half layers
     &     DENOM            ! Denominator in lapse rate expression
C*L------------------COMDECK C_G----------------------------------------
C G IS MEAN ACCEL DUE TO GRAVITY AT EARTH'S SURFACE
      REAL G

      PARAMETER(G=9.80665)
C*----------------------------------------------------------------------
C*L------------------COMDECK C_R_CP-------------------------------------
C R IS GAS CONSTANT FOR DRY AIR
C CP IS SPECIFIC HEAT OF DRY AIR AT CONSTANT PRESSURE
C PREF IS REFERENCE SURFACE PRESSURE
      REAL R,CP,KAPPA,PREF  ! PREFTOKI

      PARAMETER(R=287.05,
     &          CP=1005.,
     &          KAPPA=R/CP,
     &          PREF=100000.)
C*----------------------------------------------------------------------

C*L------------------COMDECK C_LAPSE ----------------------------------
      REAL LAPSE,LAPSE_TROP
      PARAMETER(LAPSE=0.0065)     !  NEAR SURFACE LAPSE RATE
      PARAMETER(LAPSE_TROP=0.002) !  TROPOPAUSE LAPSE RATE
C*----------------------------------------------------------------------

      REAL CP_OVER_G, P_EXNER_500, P_EXNER_50
      PARAMETER ( CP_OVER_G = CP / G )

C----------------------------------------------------------------------
C*L------------------ COMDECK P_EXNERC ---------------------------------
C statement function to define exner pressure at model levels
C from exner pressures at model half-levels
      REAL P_EXNER_C
      REAL R_P_EXNER_C
      REAL P_EXU_DUM,P_EXL_DUM,PU_DUM,PL_DUM,KAPPA_DUM
C consistent with geopotential see eqn 26 DOC PAPER 10
      P_EXNER_C(P_EXU_DUM,P_EXL_DUM,PU_DUM,PL_DUM,KAPPA_DUM) =
     & (P_EXU_DUM*PU_DUM - P_EXL_DUM*PL_DUM)/
     & ( (PU_DUM-PL_DUM)*(KAPPA_DUM + 1) )
      R_P_EXNER_C(P_EXU_DUM,P_EXL_DUM,PU_DUM,PL_DUM,KAPPA_DUM) =
     & ( (PU_DUM-PL_DUM)*(KAPPA_DUM + 1) )/
     & (P_EXU_DUM*PU_DUM - P_EXL_DUM*PL_DUM)

C*------------------- --------------------------------------------------


C     ! 1.  Set up local constants and initialise arrays

      P_EXNER_500 = (500./1000.)**KAPPA
      P_EXNER_50  =  (50./1000.)**KAPPA

      DTI = ( MIN_TROP_LEVEL + MAX_TROP_LEVEL ) / 2

      DO I=1, L2
        LTROP(I) = .TRUE.
      ENDDO

CL    ! Compute lapse rate between full levels: equation 3.16, UMDP S1

      DO J=MIN_TROP_LEVEL+1, MIN(MAX_TROP_LEVEL+1,P_LEVELS)
        DO I=1, L2
C         ! Exner pressure at full levels
          PJP1 = AKH(J+1) + BKH(J+1) * PSTAR(I)
          PJ   = AKH(J)   + BKH(J)   * PSTAR(I)
          PJM1 = AKH(J-1) + BKH(J-1) * PSTAR(I)
          P_EXNER_FULL_J   = P_EXNER_C
     &      (P_EXNER_HALF(I,J+1), P_EXNER_HALF(I,J),  PJP1, PJ, KAPPA)
          P_EXNER_FULL_JM1 = P_EXNER_C
     &      (P_EXNER_HALF(I,J),  P_EXNER_HALF(I,J-1), PJ, PJM1, KAPPA)
C         ! Exner pressure difference across half layers
          DEL_EXNER_J   = P_EXNER_HALF(I,J) - P_EXNER_FULL_J
          DEL_EXNER_JM1 = P_EXNER_FULL_JM1  - P_EXNER_HALF(I,J)
C         ! Denominator
          DENOM = ( T(I,J-1) * DEL_EXNER_JM1 / P_EXNER_FULL_JM1
     &                + T(I,J) * DEL_EXNER_J / P_EXNER_FULL_J )
C         !  Lapse rate between level j-1 and j
          LAPSE_RATE(I,J) = ( T(I,J-1) - T(I,J) )/( CP_OVER_G * DENOM )
        ENDDO
      ENDDO

CL    ! 2.  Find level of tropopause, where it is well defined

      DO J=MIN_TROP_LEVEL+1, MAX_TROP_LEVEL

C 'J+1' level for lapse rate test; allows J iteration up to P_LEVELS
        JP1=MIN(J+1,P_LEVELS)

        DO I=1, L2
C         ! Exner pressure at full levels
          PJP1 = AKH(J+1) + BKH(J+1) * PSTAR(I)
          PJ   = AKH(J)   + BKH(J)   * PSTAR(I)
          PJM1 = AKH(J-1) + BKH(J-1) * PSTAR(I)
          P_EXNER_FULL_J = P_EXNER_C
     &       (P_EXNER_HALF(I,J+1), P_EXNER_HALF(I,J), PJP1, PJ, KAPPA)
          P_EXNER_FULL_JM1 = P_EXNER_C
     &       (P_EXNER_HALF(I,J), P_EXNER_HALF(I,J-1), PJ, PJM1, KAPPA)

C         ! Not-quite-WMO criteria for interval containing tropopause
C         ! (where 'interval' stretches between layer centres j and j-1)

          IF (  P_EXNER_FULL_JM1 .GT. P_EXNER_50  .AND.
     &            P_EXNER_FULL_J .LT. P_EXNER_500 .AND.
     &           LAPSE_RATE(I,J) .LT. LAPSE_TROP  .AND.
     &         LAPSE_RATE(I,JP1) .LT. LAPSE_TROP  .AND. LTROP(I) )
     &    THEN
            LTROP(I)=.FALSE.
            IT(I) = J
          ENDIF
        ENDDO
      ENDDO


CL    ! 4.  Fill in where the above criteria did not find a tropopause

C     !  Run through all internal points and where no tropopause was
C     !   found, set the level to the average of those found at the 8
C     !   surrounding points.  If none of these did find one, then DTI,
C     !   the middle of the permitted range, is set.
C     !   (Using integers means rounding down - probably the best thing
C     !   overall as the lack of vertical resolution to pick out a
C     !   tropopause precisely is likely to mean they'll be diagnosed
C     !   too high, if anything, at neighbouring points.  Similarly
C     !   it's not worth worrying about level number being non-linear
C     !   in height or pressure at those points where ipso facto the
C     !   result is so arbitrary.)

      DO J=1, L2/ROW_LENGTH-2
        DO I=2, ROW_LENGTH-1
          POINT = I+J*ROW_LENGTH
          IF ( LTROP(POINT) ) THEN
            FILLIN = 0
            NNEIGH = 0
            IF ( .NOT. LTROP(POINT-1-ROW_LENGTH) ) THEN
              FILLIN = FILLIN + IT(POINT-1-ROW_LENGTH)
              NNEIGH = NNEIGH + 1
            ENDIF
            IF ( .NOT. LTROP(POINT+1-ROW_LENGTH) ) THEN
              FILLIN = FILLIN + IT(POINT+1-ROW_LENGTH)
              NNEIGH = NNEIGH + 1
            ENDIF
            IF ( .NOT. LTROP(POINT-1+ROW_LENGTH) ) THEN
              FILLIN = FILLIN + IT(POINT-1+ROW_LENGTH)
              NNEIGH = NNEIGH + 1
            ENDIF
            IF ( .NOT. LTROP(POINT+1+ROW_LENGTH) ) THEN
              FILLIN = FILLIN + IT(POINT+1+ROW_LENGTH)
              NNEIGH = NNEIGH + 1
            ENDIF
            IF ( .NOT. LTROP(POINT-1) ) THEN
              FILLIN = FILLIN + IT(POINT-1)
              NNEIGH = NNEIGH + 1
            ENDIF
            IF ( .NOT. LTROP(POINT+1) ) THEN
              FILLIN = FILLIN + IT(POINT+1)
              NNEIGH = NNEIGH + 1
            ENDIF
            IF ( .NOT. LTROP(POINT-ROW_LENGTH) ) THEN
              FILLIN = FILLIN + IT(POINT-ROW_LENGTH)
              NNEIGH = NNEIGH + 1
            ENDIF
            IF ( .NOT. LTROP(POINT+ROW_LENGTH) ) THEN
              FILLIN = FILLIN + IT(POINT+ROW_LENGTH)
              NNEIGH = NNEIGH + 1
            ENDIF
            IF ( NNEIGH .EQ. 0 ) THEN
               IT(POINT) = DTI
             ELSE
               IT(POINT) = FILLIN / NNEIGH
            ENDIF
          ENDIF
        ENDDO
      ENDDO

C     !  If the grid wraps round, the first & last points on each row
C     !    are in fact internal points & can be filled in similarly if
C     !    a few indices are altered.  (Indeed, they must be if the
C     !    model is to give the same answers regardless of
C     !    decomposition - there may be no physical justification for
C     !    being so scrupulous.)
C     !  If not, again such points return DTI if no tropopause was
C     !    found on the WMO criteria.

      IF ( WRAP ) THEN
         DO J=1, L2/ROW_LENGTH-2
           POINT = 1+J*ROW_LENGTH
           IF ( LTROP(POINT) ) THEN
             FILLIN = 0
             NNEIGH = 0
             IF ( .NOT. LTROP(POINT-1) ) THEN
               FILLIN = FILLIN + IT(POINT-1)
               NNEIGH = NNEIGH + 1
             ENDIF
             IF ( .NOT. LTROP(POINT+1-ROW_LENGTH) ) THEN
               FILLIN = FILLIN + IT(POINT+1-ROW_LENGTH)
               NNEIGH = NNEIGH + 1
             ENDIF
             IF ( .NOT. LTROP(POINT-1+2*ROW_LENGTH) ) THEN
               FILLIN = FILLIN + IT(POINT-1+2*ROW_LENGTH)
               NNEIGH = NNEIGH + 1
             ENDIF
             IF ( .NOT. LTROP(POINT+1+ROW_LENGTH) ) THEN
               FILLIN = FILLIN + IT(POINT+1+ROW_LENGTH)
               NNEIGH = NNEIGH + 1
             ENDIF
             IF ( .NOT. LTROP(POINT-1+ROW_LENGTH) ) THEN
               FILLIN = FILLIN + IT(POINT-1+ROW_LENGTH)
               NNEIGH = NNEIGH + 1
             ENDIF
             IF ( .NOT. LTROP(POINT+1) ) THEN
               FILLIN = FILLIN + IT(POINT+1)
               NNEIGH = NNEIGH + 1
             ENDIF
             IF ( .NOT. LTROP(POINT-ROW_LENGTH) ) THEN
               FILLIN = FILLIN + IT(POINT-ROW_LENGTH)
               NNEIGH = NNEIGH + 1
             ENDIF
             IF ( .NOT. LTROP(POINT+ROW_LENGTH) ) THEN
               FILLIN = FILLIN + IT(POINT+ROW_LENGTH)
               NNEIGH = NNEIGH + 1
             ENDIF
             IF ( NNEIGH .EQ. 0 ) THEN
                IT(POINT) = DTI
              ELSE
                IT(POINT) = FILLIN / NNEIGH
             ENDIF
           ENDIF
           POINT = (1+J)*ROW_LENGTH
           IF ( LTROP(POINT) ) THEN
             FILLIN = 0
             NNEIGH = 0
             IF ( .NOT. LTROP(POINT-1-ROW_LENGTH) ) THEN
               FILLIN = FILLIN + IT(POINT-1-ROW_LENGTH)
              NNEIGH = NNEIGH + 1
             ENDIF
             IF ( .NOT. LTROP(POINT+1-2*ROW_LENGTH) ) THEN
               FILLIN = FILLIN + IT(POINT+1-2*ROW_LENGTH)
               NNEIGH = NNEIGH + 1
             ENDIF
             IF ( .NOT. LTROP(POINT-1+ROW_LENGTH) ) THEN
               FILLIN = FILLIN + IT(POINT-1+ROW_LENGTH)
               NNEIGH = NNEIGH + 1
             ENDIF
             IF ( .NOT. LTROP(POINT+1) ) THEN
               FILLIN = FILLIN + IT(POINT+1)
               NNEIGH = NNEIGH + 1
             ENDIF
             IF ( .NOT. LTROP(POINT-1) ) THEN
               FILLIN = FILLIN + IT(POINT-1)
               NNEIGH = NNEIGH + 1
             ENDIF
             IF ( .NOT. LTROP(POINT+1-ROW_LENGTH) ) THEN
               FILLIN = FILLIN + IT(POINT+1-ROW_LENGTH)
               NNEIGH = NNEIGH + 1
             ENDIF
             IF ( .NOT. LTROP(POINT-ROW_LENGTH) ) THEN
               FILLIN = FILLIN + IT(POINT-ROW_LENGTH)
               NNEIGH = NNEIGH + 1
             ENDIF
             IF ( .NOT. LTROP(POINT+ROW_LENGTH) ) THEN
               FILLIN = FILLIN + IT(POINT+ROW_LENGTH)
               NNEIGH = NNEIGH + 1
             ENDIF
             IF ( NNEIGH .EQ. 0 ) THEN
                IT(POINT) = DTI
              ELSE
                IT(POINT) = FILLIN / NNEIGH
             ENDIF
           ENDIF
         ENDDO
       ELSE ! if not WRAP
         DO J=1, L2/ROW_LENGTH-2
           IF ( LTROP(1+J*ROW_LENGTH) )   IT(1+J*ROW_LENGTH)   = DTI
           IF ( LTROP((1+J)*ROW_LENGTH) ) IT((1+J)*ROW_LENGTH) = DTI
         ENDDO
      ENDIF
C
      DO J=1, ROW_LENGTH
        IF ( LTROP(J) )      IT(J)      = DTI
        IF ( LTROP(L2+1-J) ) IT(L2+1-J) = DTI
      ENDDO
C
      RETURN
      END
