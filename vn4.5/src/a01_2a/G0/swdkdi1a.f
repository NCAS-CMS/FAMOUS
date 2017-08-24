C ******************************COPYRIGHT******************************
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
C
CLL  Subroutine SWDKDI -------------------------------------------------
CLL
CLL  Its function is to calculate SW diagnostics (those that are not
CLL  naturally zero at night points) in the case when the entire domain
CLL  is in darkness.  It simply reproduces the relevant code from SWRAD,
CLL  which cannot be CALLed itself as it dynamically allocates arrays
CLL  dimensioned by the number of sunlit points.  For more information,
CLL  see that routine.
CLL
CLL  Model           Modification history:
CLL version  Date
CLL   4.3  16/10/96  Written by William Ingram & reviewed by Paul Burton
CLL
CLLEND --------------------------------------------------------------
C*L
      SUBROUTINE SWDKDI (ABIN, BBIN, LCAIN, CCAIN,
     &     LCA3L, LCA3ON, TCASW, TCASWO, LCLD3,
     &     NDO, NLEVS, NCLDS, L1)
      IMPLICIT NONE
      EXTERNAL SWDTCA
C*
C
C     !   Dimensions:
      INTEGER!, INTENT(IN) ::
     &     L1,                       ! Number of points in input arrays
     &     NDO,                      ! Number of points to be treated
     &     NLEVS,                    ! Number of levels
     &     NCLDS                     ! Number of possibly cloudy levels
C     !  Physical inputs:
      REAL!,
     &     ABIN(NLEVS+1), BBIN(NLEVS+1), ! As and Bs at layer boundaries
     &     LCAIN(L1,1/(NCLDS+1)+NCLDS),  ! Layer cloud fractional cover
     &     CCAIN(L1)                     ! Convective Cloud Amount
C     !  Control quantities:
      LOGICAL!, INTENT(IN) ::
     &     LCLD3,                  ! Is the 3-cloud trick on (2A SW) ?
     &     TCASWO,                 !      Is TCASW wanted ?
     &     LCA3ON                  !      And LCA3L ?
C     ! Note that if LCLD3, LCA3L is needed to calculate TCASW & so
C     !  will be calculated whenever TCASWO or LCA3ON - so space must
C     !  then be available (via "implied diagnostics" in the std UM).
C     !  And outputs:
      REAL!, INTENT(OUT) ::
     &     TCASW(L1),              !   Total cloud amount in SW
C     ! (i.e. fraction of the grid-box with cloud at some level)
     &     LCA3L(L1,NCLDS)         ! Diagnostic of layer cloud amount
C     ! restricted to 3 layers, calculated at all points on SW timesteps
C*
C     !  Constants:
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

C*L
CL    !  Dynamically allocated workspace:
      INTEGER INDEX(NDO)
C     !  Index for maximum(input)/only(used) cloud cover for a "type"
      REAL MAXCLD(NDO)                   !  Maximum cloud cover
     &                                   !                  for a "type"
C*
      INTEGER LEVEL, J,                  ! Loopers over level and point
     &     TYPE,                         !       & cloud "type" (H/M/L)
     &     RANGE(3,2),                   ! The range of level numbers
C     !  (counting down from the highest potentially cloudy level) for
C     !  the 3 cloud "types" - i.e. the RANGE(n,1)th to RANGE(n,2)th
C     !  potentially cloudy levels are assigned to the nth cloud type.
C     !  The values are set by comparing model eta values with BOUNDS.
     &     FSTLEV,                       ! The equivalent of RANGE for
     &     LSTLEV,                       !  a particular cloud type, but
C                                        !  counting up from the surface
     &     NCLEAR                        ! NLEVS-NCLDS
      REAL BOUNDS(2),                    ! Eta values that define where
C     ! cloud changes from "high" to "medium", & from "medium" to "low"
     &     ETA,                          ! Eta at the layer boundary
C     !                                  !    currently being checked
     &     ETALST                        !       & the previous one
      LOGICAL SET                        ! Has RANGE been set yet ?
      DATA BOUNDS / .37, .79 /
      DATA SET / .FALSE. /
      SAVE RANGE, SET                    ! SET must be specified too as
C     !   FORTRAN requires a variable initialized by a DATA statement to
C     !   have the SAVE attribute only if its value has not changed.

CL    ! If LCLD3 is on, the first time into the routine, find where
CL    ! cloud type boundaries will lie in terms of the numbering of this
CL    !  run's eta levels:
C
      IF ( LCLD3 .AND. .NOT. SET ) THEN
        NCLEAR = NLEVS - NCLDS
        RANGE(1,1) = 1
        LEVEL = NCLEAR + 1
        DO J=1, 2
  101     ETA = BBIN(NLEVS+2-LEVEL) + ABIN(NLEVS+2-LEVEL) / PREF
          IF ( ETA .LT. BOUNDS(J) ) THEN
             LEVEL  = LEVEL + 1
             ETALST = ETA
C            ! This assumes the vertical resolution is not too crude in
C            !    the troposphere - but it would have to be rather worse
C            !    even than the old 11-layer Cyber climate model.
             GO TO 101
           ELSE
C            ! This has found the first layer boundary below BOUNDS -
C            !   is this or the previous one closer ?
             IF ( BOUNDS(J)-ETALST .LT. ETA-BOUNDS(J) ) LEVEL = LEVEL-1
             RANGE(J+1,1) = LEVEL - NCLEAR
             RANGE(J,2)   = RANGE(J+1,1) - 1
          ENDIF
        ENDDO
        RANGE(3,2) = NCLDS
        SET = .TRUE.
      ENDIF
C
C     ! Next IF block calculates the diagnostic LCA3L.
C     ! This must be done if this diagnostic is wanted in its own right
C     !   or, when the 3-cloud trick is on, if TCASW is, as then the
C     !   latter is calculated from LCA3L.
C
      IF ( LCA3ON .OR. TCASWO .AND. LCLD3 ) THEN
        DO TYPE=1, 3
          FSTLEV = NCLDS + 1 - RANGE(TYPE,2)
          LSTLEV = NCLDS + 1 - RANGE(TYPE,1)
Cfpp$     Select(CONCUR)
          DO J=1, NDO
            MAXCLD(J) = LCAIN(J,FSTLEV)
            INDEX(J)  = FSTLEV
          ENDDO
          DO LEVEL=FSTLEV+1, LSTLEV
Cfpp$       Select(CONCUR)
            DO 163 J=1, NDO
              IF ( MAXCLD(J) .LT. LCAIN(J,LEVEL) ) THEN
                MAXCLD(J) = LCAIN(J,LEVEL)
                INDEX(J) = LEVEL
              ENDIF
  163       CONTINUE                             ! Next J
          ENDDO                                  ! Next LEVEL
          DO LEVEL=FSTLEV, LSTLEV
Cfpp$       Select(CONCUR)
            DO 164 J=1, NDO
              IF ( LEVEL .EQ. INDEX(J) ) THEN
                 LCA3L(J,LEVEL) = MAXCLD(J)
               ELSE
                 LCA3L(J,LEVEL) = 0.
              ENDIF
  164       CONTINUE                            ! Next J
          ENDDO                                 ! Next LEVEL
        ENDDO                                   ! Next TYPE
      ENDIF   !  LCA3ON .OR. TCASWO .AND. LCLD3
C
C
CL    !  If wanted, diagnose total cloud amount as seen by the SW:
C
      IF ( TCASWO ) THEN
        IF ( LCLD3 ) THEN
           CALL SWDTCA (LCA3L, CCAIN, NCLDS, L1, NDO, TCASW)
         ELSE
           CALL SWDTCA (LCAIN, CCAIN, NCLDS, L1, NDO, TCASW)
        ENDIF
      ENDIF
C
      RETURN
      END
