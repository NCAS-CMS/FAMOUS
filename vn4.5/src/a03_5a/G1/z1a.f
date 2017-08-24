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
C*LL  SUBROUTINE Z -----------------------------------------------------
CLL
CLL   Purpose: Calculate virtual temperature at one model level, and
CLL            depth of layer containing this level, and height of
CLL            the top of this layer.
CLL
CLL        Suitable for Single Column use.
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  date
CLL   3.4   20/04/94  DEF TIMER replaced by LOGICAL LTIMER
CLL                   Argument LTIMER added
CLL                                                 S.J.Swarbrick
CLL
CLL
CLL   Programming standard: Unified Model Documentation Paper No 4,
CLL                         Version 2, dates 18/1/90.
CLL
CLL   System component covered: ancillary to P24.
CLL
CLL   Documentation: UM Documentation Paper 24, section P243.
CLL                  See especially Appendix A.
CLL
C*----------------------------------------------------------------------
C*L
C-----------------------------------------------------------------------
CL  Arguments :-
C-----------------------------------------------------------------------
      SUBROUTINE Z (
     + POINTS
     +,EXNER_LOWER,EXNER_UPPER,PSTAR,AKH,BKH,Q,QCF,QCL,T,Z_LOWER
     +,TV,Z_UPPER,DELTA_Z,DELTA_Z_LOWER,LTIMER
     +)
      IMPLICIT NONE
      LOGICAL LTIMER
      INTEGER
     + POINTS              ! IN No of gridpoints being processed.
      REAL
     + EXNER_LOWER(POINTS) ! IN Exner function for lower boundary of
C                          !    this layer.
     +,EXNER_UPPER(POINTS) ! IN Exner function for upper boundary of
C                          !    this layer.
     +,PSTAR(POINTS)       ! IN surface pressure (Pa)
     +,AKH(2)              ! IN AK value at bottom and top of this layer
     +,BKH(2)              ! IN BK value at bottom and top of this layer
     +,Q(POINTS)           ! IN Sp humidity at this level (kg water
C                          !    per kg of air).
     +,QCF(POINTS)         ! IN Cloud ice at this level (kg per
C                          !    kg of air).
     +,QCL(POINTS)         ! IN Cloud liquid water at this level (kg
C                          !    per kg of air).
     +,T(POINTS)           ! IN Temperature at this level (K).
     +,Z_LOWER(POINTS)     ! IN Height above surface of lower boundary
C                          !    of this layer (metres).
      REAL
     + TV(POINTS)          ! OUT Virtual temperature for this level
C                          !     (K).
     +,Z_UPPER(POINTS)     ! OUT Height above surface of upper boundary
C                          !     of this layer (metres).
     +,DELTA_Z(POINTS)     ! OUT Depth of this layer (metres).
     +,DELTA_Z_LOWER(POINTS) ! OUT Depth of lower half layer (metres).
C*
C*L
      EXTERNAL TIMER
C*
C-----------------------------------------------------------------------
C*L Local and other parameters.
C-----------------------------------------------------------------------
C*L------------------COMDECK C_EPSLON-----------------------------------
C EPSILON IS RATIO OF MOLECULAR WEIGHTS OF WATER AND DRY AIR
      REAL EPSILON,C_VIRTUAL

      PARAMETER(EPSILON=0.62198,
     &          C_VIRTUAL=1./EPSILON-1.)
C*----------------------------------------------------------------------

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

      REAL CPRG
      PARAMETER (
     + CPRG=CP/G           ! CP upon G.
     +)
C*
C-----------------------------------------------------------------------
C  Declare local variable.
C-----------------------------------------------------------------------
      INTEGER
     + I                   ! Loop counter; horizontal field index.
C-----------------------------------------------------------------------
CL  No significant structure.
C-----------------------------------------------------------------------

      REAL
     &    PU,PL
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


      IF (LTIMER) THEN
      CALL TIMER('Z       ',3)
      ENDIF
      DO 1 I=1,POINTS
C
C  Calculate virtual temperature.  Cf eqn P243.A2 (which calculates
C  virtual potential temperature).
C          ~~~~~~~~~
        TV(I) = T(I) * ( 1.0 + C_VIRTUAL*Q(I) - QCF(I) - QCL(I) )
C
C  Calculate layer depth, eqn P243.A3.
C
        PU=PSTAR(I)*BKH(2) + AKH(2)
        PL=PSTAR(I)*BKH(1) + AKH(1)
        DELTA_Z(I) = CPRG * ( TV(I) /
     +   P_EXNER_C( EXNER_UPPER(I),EXNER_LOWER(I),PU,PL,KAPPA) !Exner k
     +                      ) *
     +    ( EXNER_LOWER(I) - EXNER_UPPER(I) )  ! -(Exner k+1/2 - k-1/2)
C
C  Calculate lower half layer depth, eqn P243.A6.
C
        DELTA_Z_LOWER(I) = CPRG *  TV(I) *
     +    ( EXNER_LOWER(I) /
     +   P_EXNER_C( EXNER_UPPER(I),EXNER_LOWER(I),PU,PL,KAPPA) -1.)
C
C  Calculate height of top of layer, eqn P243.A4.
C
        Z_UPPER(I) = Z_LOWER(I) + DELTA_Z(I)
    1 CONTINUE
      IF (LTIMER) THEN
      CALL TIMER('Z       ',4)
      ENDIF
      RETURN
      END
