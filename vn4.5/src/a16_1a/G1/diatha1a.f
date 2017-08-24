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
CLL  SUBROUTINE THADV---------------------------------------------------
CLL
CLL  PURPOSE:   Calculates thermal advection at pressure levels
CLL  Tested under compiler CFT77
CLL  Tested under OS version 5.1
CLL
CLL D.Robinson  <- programmer of some or all of previous code or changes
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  Date
CLL   4.4  23/06/97   Initialise TH_TP. D. Robinson.
!LL   4.5  15/04/98   Added start-end arguments to V_INT and
!LL                   V_INT_T routines. S.D.Mullerworth
CLL
CLL  Logical components covered: D424
CLL
CLL  Project TASK:
CLL
CLL  Programming standard: U M DOC  Paper NO. 4,
CLL
CLL  External documentation:
CLL
CLLEND------------------------------------------------------------------
C
C*L  ARGUMENTS:---------------------------------------------------------
      SUBROUTINE DIA_THADV(
C data in
     & U,V,THETA,PPT,PUV,PSTAR,PRESS_REQD,P_EXNER_HALF,
C data out
     & TH_TP,
C constants in
     & P_FIELD,U_FIELD,P_LEVELS,ROW_LENGTH,P_ROWS,SEC_U_LATITUDE,
     & EW_SPACE,NS_SPACE,AKH,BKH)
C*
C*L
C-----------------------------------------------------------------------
      IMPLICIT NONE
C-----------------------------------------------------------------------
      EXTERNAL  V_INT_T,V_INT,UV_TO_P_FULL
C-----------------------------------------------------------------------
      INTEGER
     * P_FIELD                       ! IN  NO OF POINTS ON P/T GRID
     *,U_FIELD                       ! IN  NO OF POINTS ON U/V GRID
     *,P_LEVELS                      ! IN  NO OF MODEL LEVELS
     *,P_ROWS                        ! IN  NO OF ROWS ON P/T GRID
     *,ROW_LENGTH                    ! IN  NO OF COLUMNS
C-----------------------------------------------------------------------
      REAL
     * U(U_FIELD,P_LEVELS)              ! IN  U FIELD   AT FULL LEVELS
     *,V(U_FIELD,P_LEVELS)              ! IN  V FIELD   AT FULL LEVELS
     *,THETA(P_FIELD,P_LEVELS)          ! IN  POTENTIAL TEMPERATURE " "
     *,PRESS_REQD                       ! IN  PRESSURE AT WHICH TO
     *                                  !    CALCULATE THERMAL ADVECTION
     *,TH_TP(P_FIELD)                   ! OUT  THERMAL ADVECTION AT
     *                                  !      REQUIRED LEVEL
     *,PPT(P_FIELD,P_LEVELS)            ! IN  PRESS FIELD AT P/T POINTS
     *,PUV(U_FIELD,P_LEVELS)            ! IN  PRESS FIELD AT U/V POINTS
     *,PSTAR(P_FIELD)                   ! IN  SURFACE PRESSURE FIELD
     *,P_EXNER_HALF(P_FIELD,P_LEVELS+1) ! IN  EXNER PRESSURE AT MODEL
     *                                  !     HALF LEVELS
     *,SEC_U_LATITUDE(U_FIELD)          ! IN  1/COS(LAT) AT U/V POINTS
     *,EW_SPACE                         ! IN  LONGITUDE GRID SPACING
     *,NS_SPACE                         ! IN  LATITUDE GRID SPACING
     *                                  !     BOTH THE ABOVE IN DEGREES
     *,AKH,BKH                          ! IN  Hybrid Coords. A and B
     *                                  !     values on half levels.
C*
C*L
C-----------------------------------------------------------------------
C Local Variables
C-----------------------------------------------------------------------
      INTEGER
     * I,J,K               !  LOOP COUNTERS
     *,T_REF             !  REFERENCE LEVEL FOR VERTICAL INTERPOLATION
     *                   !  OF TEMPERATURES
C-----------------------------------------------------------------------
      REAL
     * LONG_SI_OVER_TWO_A         ! LONGITUDE STEP INVERSE (IN RADIANS)
     *                            ! DIVIDED BY (2*EARTH'S RADIUS(A))
     *,LAT_SI_OVER_TWO_A          ! LATITUDE STEP INVERSE (IN RADIANS)
     *                            ! DIVIDED BY (2*EARTH'S RADIUS(A))
     *,T_P(P_FIELD)               ! TEMPERATURE FIELD ON REQUIRED LEVEL
     *,U_P(U_FIELD)               ! U FIELD ON REQUIRED LEVEL
     *,V_P(U_FIELD)               ! V FIELD ON REQUIRED LEVEL
     *,TH_UV(U_FIELD)             ! THERMAL ADVECTION AT U/V POINTS
     *,PPT_REQ(P_FIELD)           ! PRESSURE SURFACE REQUIRED ARRAY
     *,DUMMY1                     ! DUMMY ARGUMENT1
     *,DUMMY2                     ! DUMMY ARGUMENT2
     *,TERM1                      ! 1st TERM USED IN CALCULATION
     *,TERM2                      ! 2nd TERM USED IN CALCULATION
C-----------------------------------------------------------------------
C Constants required : A=radius of Earth,
C                      RECIP_PI_OVER_180=180/Pi
C-----------------------------------------------------------------------
C*L------------------COMDECK C_A----------------------------------------
C A IS MEAN RADIUS OF EARTH
      REAL A

      PARAMETER(A=6371229.)
C*----------------------------------------------------------------------

C*L------------------COMDECK C_PI---------------------------------------
CLL
CLL 4.0 19/09/95  New value for PI. Old value incorrect
CLL               from 12th decimal place. D. Robinson
CLL
      REAL PI,PI_OVER_180,RECIP_PI_OVER_180

      PARAMETER(
     & PI=3.14159265358979323846, ! Pi
     & PI_OVER_180 =PI/180.0,   ! Conversion factor degrees to radians
     & RECIP_PI_OVER_180 = 180.0/PI ! Conversion factor radians to
     &                              ! degrees
     & )
C*----------------------------------------------------------------------

      T_REF=2                  ! SET TO SECOND MODEL LEVEL
C*----------------------------------------------------------------------
C-----------------------------------------------------------------------
CL    Initialise TH_TP
C-----------------------------------------------------------------------
      DO I=1,P_FIELD
        TH_TP(I)=0.0
      ENDDO
C-----------------------------------------------------------------------
CL    Set up Pressure required array
C-----------------------------------------------------------------------
      DO I=1,P_FIELD
        PPT_REQ(I)=PRESS_REQD
      ENDDO
C-----------------------------------------------------------------------
CL    Interpolate the temperature field onto the required level
C-----------------------------------------------------------------------
      CALL V_INT_T(T_P,PPT_REQ,PPT(1,T_REF),PSTAR,P_EXNER_HALF,
     &THETA,P_FIELD,P_LEVELS,T_REF,AKH,BKH,1,P_FIELD)
C-----------------------------------------------------------------------
CL    Interpolate the U field onto the required level
C-----------------------------------------------------------------------
      CALL V_INT(PUV,PPT_REQ,U,U_P,
     &U_FIELD,P_LEVELS,DUMMY1,DUMMY2,.FALSE.,1,U_FIELD)
C-----------------------------------------------------------------------
CL    Interpolate the V field onto the required level
C-----------------------------------------------------------------------
       CALL V_INT(PUV,PPT_REQ,V,V_P,
     & U_FIELD,P_LEVELS,DUMMY1,DUMMY2,.FALSE.,1,U_FIELD)
C-----------------------------------------------------------------------
CL    Calculation of thermal  advection
C-----------------------------------------------------------------------
C
C     *T4     _      *T1     Y = Latitude
C             ^              X = Longitude
C            dY              dY= Latitude interval in radians
C             ^              dX= Longitude interval in radians
C     <---dX--*U,V--->       A = Radius of Earth in metres
C             ^              T1-4= Temperatures in Kelvin  (T_P)
C             ^              U,V= Wind components in m/s   (U_P,V_P)
C             ^
C     *T3     _      *T2
C
C                                      Term1               Term2
C                      d(Temp)    ( U*(T1+T2-T3-T4)   V*(T4+T1-T2-T3))
C  Thermal advection = ------- = -( --------------- + ---------------)
C   Degrees/second     d(Time)    (  Cos(Y)*dX*2*A       2*dY*A      )
C
C-----------------------------------------------------------------------
CL    Calculate 1/(dX*2*A) and 1/(dY*2*A)
C-----------------------------------------------------------------------
      LONG_SI_OVER_TWO_A = 0.5*RECIP_PI_OVER_180/EW_SPACE/A
      LAT_SI_OVER_TWO_A  = 0.5*RECIP_PI_OVER_180/NS_SPACE/A
C-----------------------------------------------------------------------
      DO I=1,U_FIELD-1
C-----------------------------------------------------------------------
C         Calculate terms for all cases
C-----------------------------------------------------------------------
        TERM1=U_P(I)*(T_P(I+1)+T_P(I+1+ROW_LENGTH)-T_P(I+ROW_LENGTH)
     *  -T_P(I))*SEC_U_LATITUDE(I)*LONG_SI_OVER_TWO_A
C
        TERM2=V_P(I)*(T_P(I)+T_P(I+1)-T_P(I+1+ROW_LENGTH)
     *  -T_P(I+ROW_LENGTH))*LAT_SI_OVER_TWO_A
        TH_UV(I)=-(TERM1+TERM2)
      ENDDO
C-----------------------------------------------------------------------
CL     For regional model first and last values in each row recalculated
CL                  by extrapolation of adjacent two points
C-----------------------------------------------------------------------
      DO I=ROW_LENGTH,U_FIELD,ROW_LENGTH
        J=I-ROW_LENGTH+1
        TH_UV(I)=2.*TH_UV(I-1)-TH_UV(I-2)
        TH_UV(J)=2.*TH_UV(J+1)-TH_UV(J+2)
      ENDDO

C-----------------------------------------------------------------------
CL    Interpolate thermal advection from U/V to P/T points
C-----------------------------------------------------------------------
      CALL UV_TO_P_FULL
     1(TH_UV,TH_TP,U_FIELD,P_FIELD,ROW_LENGTH,P_ROWS)
C-----------------------------------------------------------------------

      RETURN
      END
