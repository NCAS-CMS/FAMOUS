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
CLL  SUBROUTINE V_INT_TP----------------------------------------------
CLL
CLL  Purpose:  To calculate the temperature along an
CLL            arbitrary pressure level. Assumes input data is
CLL            on model levels.
CLL
CLL A.Dickinson <- programmer of some or all of previous code or changes
CLL D.Robinson  <- programmer of some or all of previous code or changes
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  Date
CLL 3.1     09/02/93
CLL              : NEW ALternative routine for OUTPUT temperatures
CLL              : on pressure levels. Defines EXNER at full levels
CLL              : so as to minimise errors for isothermal atmosphere
CLL              : Reduces large biases above 100hPA in stratosphere
CLL              : C Wilson 08/02/93
CLL 4.2     01/07/96 
CLL              : Revised for CRAY T3E. Faster version of P_EXNER_C
CLL              : introduced. Unneccessary calculation of pressure
CLL              : removed. 
CLL              : New arguments START and END introduced to
CLL              : facilitate the removal of duplicate calculations
CLL              : when using domain decomposition in MPP mode.     
CLL              : Author: A. Dickinson    Reviewer: F. Rawlins     
CLL  4.5    09/01/98  CRAY T3E optimisation: replace rtor_v by powr_v
CLL                                                    Deborah Salmond
CLL                                                                     
CLL
CLL  Documentation: The interpolation formulae are described in
CLL                 unified model on-line documentation paper S1.
CLL
CLL  -----------------------------------------------------------------
C
C*L  Arguments:-------------------------------------------------------
      SUBROUTINE V_INT_TP
     *(T,P,PL,PSTAR,P_EXNER_HALF,THETA,POINTS,P_LEVELS,L,AKH,BKH)
!     *,START,END)

      IMPLICIT NONE

      INTEGER
     * POINTS    !IN Number of points per level
     *,P_LEVELS  !IN Number of model levels
     *,L         !IN Reference level for below-surface T extrapolation
                 ! Use L=2
     *,START     !IN First point to be processed in POINTS dimension
     *,END       !IN Last point to be processed in POINTS dimension

      REAL
     * T(POINTS)     !OUT Temperature along input pressure surface
     *,P(POINTS)     !IN Pressure surface on which results required
     *,PL(POINTS)    !IN Pressure at reference level L
     *,PSTAR(POINTS) !IN Surface pressure
     *,P_EXNER_HALF(POINTS,P_LEVELS+1) !IN Exner pressure at model
     *                                 !   half levels
     *,THETA(POINTS,P_LEVELS) !IN Potential temperature at full levels
     *,AKH(P_LEVELS+1) !IN Hybrid coords. A values at half levels.
     *,BKH(P_LEVELS+1) !IN Hybrid coords. B values at half levels.

C Workspace usage:-----------------------------------------------------
C
      REAL
     * P_EXNER(POINTS) 
C External subroutines called:-----------------------------------------
C None
C*---------------------------------------------------------------------
C Define local variables:----------------------------------------------
      REAL PTOP,PBOT,PKP2,PKP1,PK,PKM1,PPP1,PP,PPM1,P1,P2,P3
      INTEGER I,K,IC
      REAL P_EXNER_FULL_1,P_EXNER_FULL_2,TK,TERM1,TERM2,
     * P_EXNER_FULL_K,P_EXNER_FULL_KP1,P_EXNER_FULL_KM1
     *,P_EXNER_FULL_L,P_EXNER_FULL_LM1
C----------------------------------------------------------------------
C Constants from comdecks:---------------------------------------------
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

C*L------------------COMDECK C_LAPSE ----------------------------------
      REAL LAPSE,LAPSE_TROP
      PARAMETER(LAPSE=0.0065)     !  NEAR SURFACE LAPSE RATE
      PARAMETER(LAPSE_TROP=0.002) !  TROPOPAUSE LAPSE RATE
C*----------------------------------------------------------------------
C----------------------------------------------------------------------

CL 1. Define local constants

      REAL CP_OVER_G,LAPSE_R_OVER_G
      PARAMETER(CP_OVER_G=CP/G)
      PARAMETER(LAPSE_R_OVER_G=LAPSE*R/G)

C*L------------------ COMDECK P_EXNRC2 ---------------------------------
C statement function to define exner pressure at model levels
C from exner pressures at model half-levels
C this version for output and LASS/GLOSS processing
C Revised formulation chosen to have minimum error with isothermal
C profile
      REAL P_EXNER_C                                          
      REAL P_EXU_DUM,P_EXL_DUM,dummy1,dummy2,dummy3             
! 3 dummy args to prevent compile errors (S.J.Swarbrick)
      P_EXNER_C(P_EXU_DUM,P_EXL_DUM,dummy1,dummy2,dummy3) =     
     & (P_EXU_DUM - P_EXL_DUM)/                                    
     & (  ALOG(P_EXU_DUM/P_EXL_DUM) )                 
C                                         

C*------------------- --------------------------------------------------



      IC=0
      START=1
      END  =POINTS
CL 2. Special cases   (i) Below ground
CL                   (ii) Bottom of bottom layer
CL                  (iii) Top of top layer and above

      DO 200 I=START,END

C Convert target pressure to Exner pressure

      P_EXNER(I)=(P(I)/PREF)**KAPPA

CL (i) Below ground: equation (3.15)
CL A lapse rate of 6.5 deg/km is assumed. L is a
CL reference level - usually the first level above the model's
CL boundary layer.

      IF(P(I).GT.PSTAR(I))THEN

        P_EXNER_FULL_L = P_EXNER_C
     +  (P_EXNER_HALF(I,L+1),P_EXNER_HALF(I,L),
     +   dummy1,dummy2,dummy3)
        T(I)=THETA(I,L)*P_EXNER_FULL_L
     *  *(P(I)/PL(I))**LAPSE_R_OVER_G
        IC=IC+1

CL (ii) Bottom layer: equation (3.14)

      ELSE IF(P_EXNER(I).GT.P_EXNER_HALF(I,2)) THEN

        P_EXNER_FULL_1 = P_EXNER_C
     +  (P_EXNER_HALF(I,2),P_EXNER_HALF(I,1),
     +   dummy1,dummy2,dummy3)
        P_EXNER_FULL_2 = P_EXNER_C
     +  (P_EXNER_HALF(I,3),P_EXNER_HALF(I,2),
     +   dummy1,dummy2,dummy3)

        TK=THETA(I,1)*P_EXNER_FULL_1

        TERM1=(TK-THETA(I,2)*P_EXNER_FULL_2)
     *       *THETA(I,1)*(P_EXNER(I)-P_EXNER_FULL_1)

        TERM2=THETA(I,2)*(P_EXNER_HALF(I,2)-P_EXNER_FULL_2)
     *   +THETA(I,1)*(P_EXNER_FULL_1-P_EXNER_HALF(I,2))

       T(I)=TK+TERM1/TERM2
       IC=IC+1

      ENDIF


CL (iii) Top layer and above: equation (3.13)

      IF(P_EXNER(I).LE.P_EXNER_HALF(I,P_LEVELS))THEN

       P_EXNER_FULL_L   = P_EXNER_C (P_EXNER_HALF(I,P_LEVELS+1),
     +                    P_EXNER_HALF(I,P_LEVELS),
     +                    dummy1,dummy2,dummy3)
       P_EXNER_FULL_LM1 = P_EXNER_C (P_EXNER_HALF(I,P_LEVELS),
     +                    P_EXNER_HALF(I,P_LEVELS-1),
     +                    dummy1,dummy2,dummy3)

       TK=THETA(I,P_LEVELS)*P_EXNER_FULL_L

       TERM1=(TK-THETA(I,P_LEVELS-1)*P_EXNER_FULL_LM1)
     *       *THETA(I,P_LEVELS)*(P_EXNER_FULL_L-P_EXNER(I))

       TERM2=THETA(I,P_LEVELS)*(P_EXNER_HALF(I,P_LEVELS)-P_EXNER_FULL_L)
     * +THETA(I,P_LEVELS-1)*(P_EXNER_FULL_LM1-P_EXNER_HALF(I,P_LEVELS))

       T(I)=TK+TERM1/TERM2
       IC=IC+1

      ENDIF

200   CONTINUE

C Loop over levels

      DO 300 K=2,P_LEVELS-1

      IF(IC.EQ.END-START+1)GOTO 400


CL 3. Middle levels: equation (3.12)
CL Two alternatives are used depending on whether P_EXNER(I) falls in
CL the top or bottom half of layer k.

      DO 310 I=1,POINTS


      P_EXNER_FULL_K   = P_EXNER_C (P_EXNER_HALF(I,K+1),
     +                   P_EXNER_HALF(I,K),dummy1,dummy2,dummy3)

C Top half of layer k.

      IF(P_EXNER(I).GT.P_EXNER_HALF(I,K+1))THEN
        IF(P_EXNER(I).LE.P_EXNER_FULL_K)THEN

       P_EXNER_FULL_KP1 = P_EXNER_C (P_EXNER_HALF(I,K+2),
     +                    P_EXNER_HALF(I,K+1),dummy1,dummy2,dummy3)

       TK=THETA(I,K)*P_EXNER_FULL_K

       TERM1=(THETA(I,K+1)*P_EXNER_FULL_KP1-TK)
     *       *THETA(I,K)*(P_EXNER_FULL_K-P_EXNER(I))

       TERM2=THETA(I,K)*(P_EXNER_FULL_K-P_EXNER_HALF(I,K+1))
     *   +THETA(I,K+1)*(P_EXNER_HALF(I,K+1)-P_EXNER_FULL_KP1)

       T(I)=TK+TERM1/TERM2
       IC=IC+1

       ENDIF
      ENDIF

C Bottom half of layer k.


      IF(P_EXNER(I).GT.P_EXNER_FULL_K)THEN
       IF(P_EXNER(I).LE.P_EXNER_HALF(I,K))THEN

       P_EXNER_FULL_KM1 = P_EXNER_C (P_EXNER_HALF(I,K),
     +                    P_EXNER_HALF(I,K-1),dummy1,dummy2,dummy3)

       TK=THETA(I,K)*P_EXNER_FULL_K

       TERM1=(THETA(I,K-1)*P_EXNER_FULL_KM1-TK)
     *       *THETA(I,K)*(P_EXNER(I)-P_EXNER_FULL_K)

       TERM2=THETA(I,K)*(P_EXNER_HALF(I,K)-P_EXNER_FULL_K)
     *   +THETA(I,K-1)*(P_EXNER_FULL_KM1-P_EXNER_HALF(I,K))

       T(I)=TK+TERM1/TERM2
       IC=IC+1

       ENDIF
      ENDIF

310   CONTINUE

300   CONTINUE

400   CONTINUE

      RETURN
      END

