C ******************************COPYRIGHT******************************
C (c) CROWN COPYRIGHT 1995, METEOROLOGICAL OFFICE, All Rights Reserved.
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
CLL   SUBROUTINE TRAC_VERT_ADV
CLL
CLL   PURPOSE: CALCULATES VERTICAL ADVECTION INCREMENTS TO A FIELD AT A
CLL            SINGLE MODEL LEVEL USING A POSITIVE DEFINITE SCHEME.
CLL            IN CALCULATING THE INCREMENTS THE TEST FOR THE
CLL            DIRECTION OF THE WIND HAS BEEN REVERSED TO TAKE INTO
CLL            ACCOUNT THE CHANGE IN SIGN INTRODUCED BY MASS
CLL            WEIGHTING.
CLL   SUITABLE FOR SINGLE COLUMN USE.
CLL
CLL M.MAWSON    <- PROGRAMMER OF SOME OR ALL OF PREVIOUS CODE OR CHANGES
CLL
CLL  MODEL            MODIFICATION HISTORY FROM MODEL VERSION 3.0:
CLL VERSION  DATE
CLL 4.0   1/8/95 TOP LEVEL FLUXES AND LIMITERS MODIFIED
CLL              FOR TRACER VERTICAL ADVECTION.  T DAVIES
!LL 4.5  23/6/98 Optimisation changes
!LL              D.Salmond, B.Carruthers and P.Burton
CLL
CLL   PROGRAMMING STANDARD: UNIFIED MODEL DOCUMENTATION PAPER NO. 4,
CLL                         STANDARD B.
CLL
CLL   SYSTEM COMPONENTS COVERED: P123
CLL
CLL   SYSTEM TASK: P1
CLL
CLL   DOCUMENTATION: U. M. Doc Paper 11 by M.J.P. Cullen
CLL
CLLEND-----------------------------------------------------------------

C
C*L   ARGUMENTS:-------------------------------------------------------
      SUBROUTINE TRAC_VERT_ADV
     &                       (FIELD,ETADOT,PSTAR, P_FIELD,
     &                        ADVECTION_TIMESTEP,START_LEVEL,END_LEVEL,
     &                        FIRST_POINT,POINTS,P_LEVELS,
     &                        TR_FIRST_LEVEL,TR_LAST_LEVEL,
     &                        RS,AK,BK,DELTA_AK,DELTA_BK,
     &                        FIELD_LOWER_BOUNDARY,
     &                        L_TRACER_THETAL_QT,L_SUPERBEE)

      IMPLICIT NONE

      LOGICAL
     & L_SUPERBEE          !IN True then use SUPERBEE limiter,
     &                     !   False then use VAN LEER limiter.
     & , L_TRACER_THETAL_QT   ! TRUE IF USING TRACER ADVECTION FOR
     &                        ! THETAL AND QT

      INTEGER
     & P_FIELD             !IN  DIMENSION OF FIELDS ON PRESSURE GRID.
     &,FIRST_POINT         !IN  FIRST POINT TO BE UPDATED
     &,POINTS              !IN  NO. OF POINTS TO BE UPDATED
     &,P_LEVELS            !IN  NUMBER OF P LEVELS
     &,TR_FIRST_LEVEL      !IN  BOTTOM LEVEL FOR WHICH TRACER DEFINED
     &,TR_LAST_LEVEL       !IN  TOP LEVEL FOR WHICH TRACER DEFINED
     &,START_LEVEL         !IN  BOTTOM LEVEL FOR WHICH TRACER ADVECTED
     &,END_LEVEL           !IN  TOP LEVEL FOR WHICH TRACER ADVECTED

      REAL
     & PSTAR(P_FIELD)             !IN  SURFACE PRESSURE
     &,ETADOT(P_FIELD,P_LEVELS)   !IN  ADVECTING ETADOT FIELD,
     &                            !    MASS-WEIGHTED
     &,FIELD(P_FIELD,TR_LAST_LEVEL+1-TR_FIRST_LEVEL)   !IN  FIELD TO BE
     &                                                 !    ADVECTED.
     &,ADVECTION_TIMESTEP         !IN

      REAL
     & AK(P_LEVELS)                    !IN
     &,BK(P_LEVELS)                    !IN
     &,DELTA_AK(P_LEVELS)              !IN
     &,DELTA_BK(P_LEVELS)              !IN
     &,RS(P_FIELD,P_LEVELS)            !IN
     &,FIELD_LOWER_BOUNDARY(P_FIELD)   !IN

C*---------------------------------------------------------------------

C*L   DEFINE ARRAYS AND VARIABLES USED IN THIS ROUTINE-----------------
C DEFINE LOCAL ARRAYS: 25 ARE REQUIRED

      REAL
     & FLUX_DELTA_T(P_FIELD,-1:+1)       !  FLUX * ADVECTION TIMESTEP
     &,B1(P_FIELD,-1:+1)                 !  ARGUMENT OF B_TERM
     &,B2(P_FIELD)                 !  ARGUMENT OF B_TERM
     &,B_TERM(P_FIELD)             !
     &,COURANT(P_FIELD,-1:+1)            !  COURANT NUMBER
     &,ABS_COURANT(P_FIELD,-1:+1)    !  ABSOLUTE VALUE OF COURANT NUMBER
     &,COURANT_MW(P_FIELD,-1:+1)         !  MASS WEIGHTED COURANT NUMBER
     &,RS_SQUARED_DELTAP(P_FIELD)  !
     &,MW(P_FIELD,-1:+1)                 !  MASS WEIGHT
     &,FIELD_INC(P_FIELD,-1:+1)          !
     &,RS_SQUARED_DELTAP_RECIP(P_FIELD,-1:+1)
     &,B_SWITCH(P_FIELD)           !  ENTROPY CONDITION SWITCH.

      LOGICAL
     &  COURANT_MW_GE_0(P_FIELD)
C*---------------------------------------------------------------------
C DEFINE LOCAL VARIABLES

      INTEGER
     &  START_P_UPDATE       ! FIRST P POINT TO BE UPDATED
     &,END_P_UPDATE         ! LAST P POINT TO BE UPDATED
     &,TR_LEVEL             ! LEVEL INDEX FOR TRACER FIELD
     &,LEVEL                ! LEVEL INDEX

C COUNT VARIABLES FOR DO LOOPS ETC.

      INTEGER
     & I,IM,I0,IP,ITM,II
c
      integer saved_levels(3), saved_levels_loc(2)
c
      equivalence (saved_levels(1), im)
      equivalence (saved_levels(2), i0)
      equivalence (saved_levels(3), ip)

C*L   NO EXTERNAL SUBROUTINE CALLS:------------------------------------
C*---------------------------------------------------------------------

C CALL COMDECK FOR RADIUS OF EARTH
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


CL---------------------------------------------------------------------
CL    INTERNAL STRUCTURE.
CL---------------------------------------------------------------------
CL
CL---------------------------------------------------------------------
CL    SECTION 0.     INITIALISATION
CL---------------------------------------------------------------------

      START_P_UPDATE  = FIRST_POINT
      END_P_UPDATE    = START_P_UPDATE + POINTS - 1

      IM=-1
      I0=0
      IP=+1

      saved_levels_loc(1)=1
      saved_levels_loc(2)=3

CL
CL---------------------------------------------------------------------
CL    SECTION 1.     CALCULATE FIELD INCREMENTS FOR VERTICAL ADVECTION
CL---------------------------------------------------------------------

      DO I=START_P_UPDATE,END_P_UPDATE
        RS_SQUARED_DELTAP(I) = RS(I,1)*RS(I,1)*
     &                         (DELTA_AK(1)+DELTA_BK(1)*PSTAR(I))
        RS_SQUARED_DELTAP_RECIP(I,IP)=1./RS_SQUARED_DELTAP(I)
        MW(I,IP) = RS_SQUARED_DELTAP(I)
      ENDDO
      DO I=START_P_UPDATE,END_P_UPDATE
        COURANT_MW(I,IP) = 0.0
        COURANT(I,IP) = 0.0
      END DO

C  LOOP OVER LAYER BOUNDARIES, EXCEPT TOP AND BOTTOM ZONAL ADVECTIVE
C  FLUX ASSUMED AT TOP AND BOTTOM TRACER LEVEL

      DO LEVEL=START_LEVEL-1,END_LEVEL

        IF(LEVEL.NE.START_LEVEL-1) THEN
           ITM=IM
           IM=I0
           I0=IP
           IP=ITM
        ENDIF

        TR_LEVEL=LEVEL+1-TR_FIRST_LEVEL

CL Copy values of B1,FLUX_DELTA_T,COURANT_MW into position for
CL next level

C----------------------------------------------------------------------
CL    SECTION 1.1
C----------------------------------------------------------------------

C----------------------------------------------------------------------
CL    SECTION 1.2    CALCULATE COURANT NUMBER
C----------------------------------------------------------------------

        IF(LEVEL.GT.0) THEN
          IF(LEVEL.EQ.START_LEVEL-1) THEN
            DO I=START_P_UPDATE,END_P_UPDATE
              RS_SQUARED_DELTAP(I)=RS(I,LEVEL+1)*RS(I,LEVEL+1)*
     &                  (DELTA_AK(LEVEL+1)+DELTA_BK(LEVEL+1)*PSTAR(I))
              RS_SQUARED_DELTAP_RECIP(I,IP)=1./RS_SQUARED_DELTAP(I)
            end do
            do i=start_p_update,end_p_update
              COURANT_MW(I,IP)=ETADOT(I,LEVEL+1) *ADVECTION_TIMESTEP
              COURANT(I,IP) = COURANT_MW(I,IP)*
     &                     RS_SQUARED_DELTAP_RECIP(I,IP)
            END DO
          ELSE IF(LEVEL.LT.P_LEVELS) THEN
            DO I=START_P_UPDATE,END_P_UPDATE
              MW(I,IP) = RS_SQUARED_DELTAP(I)
              RS_SQUARED_DELTAP(I)=RS(I,LEVEL+1)*RS(I,LEVEL+1)*
     &                  (DELTA_AK(LEVEL+1)+DELTA_BK(LEVEL+1)*PSTAR(I))
              MW(I,IP) = 0.5*(MW(I,IP)+RS_SQUARED_DELTAP(I))
              RS_SQUARED_DELTAP_RECIP(I,IP)=1./RS_SQUARED_DELTAP(I)
            end do
            do i=start_p_update,end_p_update
              COURANT_MW(I,IP)=ETADOT(I,LEVEL+1) *ADVECTION_TIMESTEP
              COURANT(I,IP) = COURANT_MW(I,IP) / MW(I,IP)
            END DO
          ELSE
            DO I=START_P_UPDATE,END_P_UPDATE
              COURANT(I,IP) = 0.
            END DO
          END IF
        END IF

C----------------------------------------------------------------------
CL    SECTION 1.3    CALCULATE FLUX_DELTA_T AND B1
C----------------------------------------------------------------------

        IF(LEVEL.EQ.START_LEVEL-1) THEN
          DO I=START_P_UPDATE,END_P_UPDATE
            FLUX_DELTA_T(I,IP)=(FIELD(I,TR_LEVEL+1)
!    &                          - FIELD_LOWER_BOUNDARY(I))*
     &                                                   )*
     &                            COURANT_MW(I,IP)
          END DO
        ELSE IF(LEVEL.EQ.P_LEVELS.OR.LEVEL.EQ.TR_LAST_LEVEL) THEN
          DO I=START_P_UPDATE,END_P_UPDATE
            FLUX_DELTA_T(I,IP)=0
          END DO
        ELSE
          DO I=START_P_UPDATE,END_P_UPDATE
            FLUX_DELTA_T(I,IP) = (FIELD(I,TR_LEVEL+1)-
     &                      FIELD(I,TR_LEVEL)) * COURANT_MW(I,IP)
          END DO
        END IF

        DO I= START_P_UPDATE,END_P_UPDATE
          ABS_COURANT(I,IP) = ABS(COURANT(I,IP))
          B1(I,IP)=FLUX_DELTA_T(I,IP)*0.5*
     &                (1.0-ABS_COURANT(I,IP))
        END DO

CL End of calculation for first pass through loop.
CL Remaining calculations can only be carried out if B1 is available
CL for 3 consecutive levels.

        IF(LEVEL.GT.START_LEVEL) THEN

          DO I=START_P_UPDATE,END_P_UPDATE
            COURANT_MW_GE_0(I)=(COURANT_MW(I,I0).GE.0.0)
          ENDDO

C----------------------------------------------------------------------
CL    SECTION 1.4    CALCULATE B2
C----------------------------------------------------------------------

          IF(LEVEL.EQ.END_LEVEL) THEN
            DO I=START_P_UPDATE,END_P_UPDATE
              B2(I) = FLUX_DELTA_T(I,IM)*0.5*(MW(I,I0)/MW(I,IM)-
     &                                           ABS_COURANT(I,IM))
              B_SWITCH(I) = SIGN(1.,COURANT(I,I0)*COURANT(I,IM))
              IF (COURANT_MW_GE_0(I)) THEN
                B2(I) = 0.0
                B_SWITCH(I) = 1.0
              END IF
            END DO
          ELSE
            DO I=START_P_UPDATE,END_P_UPDATE
              ii=max(0, nint(sign(1.0, courant_mw(i, i0))))
              ii=saved_levels(saved_levels_loc(ii+1))
              B2(I)=FLUX_DELTA_T(I,ii)*0.5*(MW(I,I0)/MW(I,ii)-
     &              ABS_COURANT(I,ii))
            end do
            do i=start_p_update,end_p_update
              ii=max(0, nint(sign(1.0, courant_mw(i, i0))))
              ii=saved_levels(saved_levels_loc(ii+1))
              B_SWITCH(I) = SIGN(1.,COURANT(I,I0)*COURANT(I,ii))
            END DO
          END IF

C----------------------------------------------------------------------
CL    SECTION 1.5    CALCULATE B_TERM
C----------------------------------------------------------------------

          IF(L_SUPERBEE) THEN
CL    SUPERBEE LIMITER.

            DO I=START_P_UPDATE,END_P_UPDATE
              IF(ABS(B2(I)).GT.1.0E-8) THEN
                B_SWITCH(I) = B_SWITCH(I) * B1(I,I0)/B2(I)
                IF(B_SWITCH(I).GT.0.5.AND.B_SWITCH(I).LT.2.0) THEN
                  B_TERM(I) = B2(I) * MAX(B_SWITCH(I),1.0)
                ELSE IF (B_SWITCH(I).LE.0.0) THEN
                  B_TERM(I) = 0.0
                ELSE
                  B_TERM(I) = 2.0 * B2(I) * MIN(B_SWITCH(I),1.0)
                END IF
              ELSE
                B_SWITCH(I) = 0.0
                B_TERM(I) = 0.0
              END IF
            END DO

          ELSE

CL    VAN LEER LIMITER.

C LOOP OVER ALL POINTS
            DO I=START_P_UPDATE,END_P_UPDATE
              B_TERM(I) = 0.0
              IF (B1(I,I0)*B2(I)*B_SWITCH(I).GT.0.0)
     &           B_TERM(I) = 2.0*B1(I,I0)*B2(I)*B_SWITCH(I)/
     &                       (B1(I,I0)+B2(I)*B_SWITCH(I))
            END DO

          END IF

C----------------------------------------------------------------------
CL    SECTION 1.6    CALCULATE INCREMENTS TO FIELD
C----------------------------------------------------------------------

CDIR$ IVDEP
      IF( L_TRACER_THETAL_QT) THEN
CL  FOR TRACER ADVECTION OF THETAL & QT MODIFY TOP LEVEL LIMITER TO
CL  LOOK LIKE CENTRED SCHEME INCREMENT.
       IF(LEVEL.GE.P_LEVELS)THEN
        DO I=START_P_UPDATE,END_P_UPDATE
         B_TERM(I)=0.5*FLUX_DELTA_T(I,I0)
        END DO
       ENDIF
      ENDIF

          DO I=START_P_UPDATE,END_P_UPDATE
            IF (COURANT_MW_GE_0(I)) THEN
              FIELD_INC(I,I0) = B_TERM(I)-FLUX_DELTA_T(I,I0)
              FIELD_INC(I,IP) = -B_TERM(I)
            ELSE
              FIELD_INC(I,I0) = - B_TERM(I)
              FIELD_INC(I,IP) = B_TERM(I)-FLUX_DELTA_T(I,I0)
            ENDIF
          END DO


C----------------------------------------------------------------------
CL    SECTION 1.7    UPDATE FIELD
C----------------------------------------------------------------------

          DO I=START_P_UPDATE,END_P_UPDATE
            FIELD(I,TR_LEVEL-1) = FIELD(I,TR_LEVEL-1) + FIELD_INC(I,I0)*
     &                            RS_SQUARED_DELTAP_RECIP(I,IM)
            FIELD(I,TR_LEVEL) = FIELD(I,TR_LEVEL) +
     &                          FIELD_INC(I,IP) *
     &                          RS_SQUARED_DELTAP_RECIP(I,I0)
          END DO

        END IF

CL End loop over levels

      END DO

CL    END OF ROUTINE TRAC_VERT_ADV

      RETURN
      END
