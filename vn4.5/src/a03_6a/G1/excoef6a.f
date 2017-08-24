C *****************************COPYRIGHT******************************
C (c) CROWN COPYRIGHT 1997, METEOROLOGICAL OFFICE, All Rights Reserved.
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
!!!  SUBROUTINE EX_COEF------------------------------------------------
!!!
!!!  Purpose: To calculate exchange coefficients for boundary layer
!!!           subroutine KMKH.
!!!
!!!  Suitable for single column use (via *IF definition IBM).
!!!
!!!  Model            Modification history:
!!! version  Date
!!!
!!!  4.4     10/9/97  New deck  R.N.B.Smith
!LL   4.5  23/10/98   Prevent TIMER from performing barrier  P.Burton
!!!
!!!  Programming standard:
!!!
!!!  System component covered: Part of P243.
!!!
!!!  Project task:
!!!
!!!  Documentation: UMDP No.24
!!!
!!!---------------------------------------------------------------------

!!  Arguments :-
      SUBROUTINE EX_COEF (
     & P_FIELD,P1,P_POINTS,BL_LEVELS
     &,CCB,CCT,NTML,L_MOM
     &,CCA,DZL,RDZ,DB,U_P,V_P
     &,RHO,ZH,ZLB,Z0M,H_BLEND
     &,CUMULUS,Z_LCL
     &,RHOKM,RHOKH,LTIMER
     & )

      IMPLICIT NONE

      LOGICAL LTIMER
      LOGICAL
     & L_MOM       ! IN Switch for convective momentum transport.

      INTEGER
     & P_FIELD     ! IN No. of P-grid points in whole field.
     &,P1          ! IN First P-grid point to be processed.
     &,P_POINTS    ! IN No. of P-grid points to be processed.
     &,BL_LEVELS   ! IN maximum number of boundary layer levels

      INTEGER
     & CCB(P_FIELD)              ! IN  Convective Cloud Base.
     &,CCT(P_FIELD)              ! IN  Convective Cloud Top.
     &,NTML(P_FIELD)             ! IN  Number of turbulently mixed
                                 !     layers.

      REAL
     & CCA(P_FIELD)              ! IN Convective Cloud Amount.
     &,DZL(P_FIELD,BL_LEVELS)    ! IN Layer depths (m).  DZL(,K) is the
!                                     distance from layer boundary K-1/2
!                                     to layer boundary K+1/2.  For K=1
!                                     the lower boundary is the surface.
     &,RDZ(P_FIELD,BL_LEVELS)    ! IN Reciprocal of distance between
!                                     hybrid levels (m-1).  1/RDZ(,K) is
!                                     the vertical distance from level
!                                     K-1 to level K, except that for
!                                     K=1 it is just the height of the
!                                     lowest atmospheric level.
     &,RHO(P_FIELD,BL_LEVELS)    ! IN Density at levels where Richardson
!                                     number is calculated.
     &,DB(P_FIELD,2:BL_LEVELS)   ! IN Buoyancy jump across lower
!                                     interface of layer.
     &,U_P(P_FIELD,BL_LEVELS)    ! IN Westerly wind component
!                                     horizontally interpolated to
!                                     P-grid. (m/s)
     &,V_P(P_FIELD,BL_LEVELS)    ! IN Southerly wind component
!                                     horizontally interpolated to
!                                     P-grid. (m/s)
     &,ZLB(P_FIELD,BL_LEVELS)    ! IN ZLB(,K) is height above surface of
!                                     lower boundary of layer K (m).
     &,Z0M(P_FIELD)              ! IN Roughness length for momentum (m).
     &,H_BLEND(P_FIELD)          ! IN Blending height for effective
!                                     roughness length scheme
     &,Z_LCL(P_FIELD)            ! IN Height of lifting condensation
!                                     level (m).
      LOGICAL
     & CUMULUS(P_FIELD)          ! IN Flag for boundary layer cumulus.

      REAL
     & RHOKM(P_FIELD,BL_LEVELS)  ! INOUT Layer K-1 - to - layer K
!                                        exchange coefficient for
!                                        momentum, on UV-grid with first
!                                        and last rows set to "missing
!                                        data".
     &,RHOKH(P_FIELD,BL_LEVELS)  ! INOUT Layer K-1 - to - layer K
!                                        exchange coefficient for FTL,
!                                        on P-grid.
     &,ZH(P_FIELD)               ! INOUT Mixing layer height (m). 

!-----------------------------------------------------------------------
      EXTERNAL TIMER

!-----------------------------------------------------------------------

!    Local and other symbolic constants :-
C*L------------------COMDECK C_LHEAT------------------------------------
C LC IS LATENT HEAT OF CONDENSATION OF WATER AT 0DEGC
C LF IS LATENT HEAT OF FUSION AT 0DEGC
      REAL LC,LF

      PARAMETER(LC=2.501E6,
     &          LF=0.334E6)
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

C*L------------------COMDECK C_EPSLON-----------------------------------
C EPSILON IS RATIO OF MOLECULAR WEIGHTS OF WATER AND DRY AIR
      REAL EPSILON,C_VIRTUAL

      PARAMETER(EPSILON=0.62198,
     &          C_VIRTUAL=1./EPSILON-1.)
C*----------------------------------------------------------------------

C*L-----------------COMDECK C_VKMAN-------------------------------------
C VKMAN IS VON KARMAN'S CONSTANT
      REAL VKMAN

      PARAMETER(VKMAN=0.4)
C*----------------------------------------------------------------------


      REAL EH,EM,G0,DH,DM,LAMBDA_MIN,A_LAMBDA
      PARAMETER (
     & EH=25.0                  ! Used in calc of stability function FH.
     &,EM=4.0                   ! Used in calc of stability function FM.
     &,G0=10.0                  ! Used in stability function calcs.
     &,DH=G0/EH                 ! Used in calc of stability function FH.
     &,DM=G0/EM                 ! Used in calc of stability function FM.
     &,LAMBDA_MIN=40.0          ! Minimum value of length scale LAMBDA.
     &,A_LAMBDA=2.0             ! used in calc of LAMBDA_EFF
     &)


!  Define local storage.

!  Arrays.

      REAL
     & DVDZM(P_FIELD,2:BL_LEVELS)  ! Modulus of wind shear.
     &,RI(P_FIELD,2:BL_LEVELS)     ! Local Richardson number.
     &,ZH_LOCAL(P_FIELD)           ! Mixed layer depth determined
!                                    from local Richardson number
!                                    profile.

      INTEGER
     & NTML_LOCAL(P_FIELD)         ! Number of model layers in the 
!                                    turbulently mixed layer as
!                                    determined from the local
!                                    Richardson number profile.

      LOGICAL
     & TOPBL(P_FIELD)              ! Flag for having reached the top
!                                    of the turbulently mixed layer. 

                                                                        
!  Scalars.

      REAL
     & DZU       ! Westerly wind shear between adjacent levels.         
     &,DZV       ! Southerly wind shear between adjacent levels.
     &,DVMOD2    ! Square of modulus of wind shear between adjacent
!                  levels
     &,ELH       ! Mixing length for heat & moisture at lower layer bdy.
     &,ELM       ! Mixing length for momentum at lower layer boundary.
     &,FH        ! (Value of) stability function for heat & moisture.
     &,FM        ! (Value of) stability function for momentum transport.
     &,RTMRI     ! Temporary in stability function calculation.
     &,VKZ       ! Temporary in calculation of ELH.
     &,LAMBDAM   ! Asymptotic mixing length for turbulent transport
!                  of momentum.
     &,LAMBDAH   ! Asymptotic mixing length for turbulent transport
!                  of heat/moisture.
     &,LAMBDA_EFF! Effective mixing length used with effective
!                  roughness length scheme.
      INTEGER
     & I         ! Loop counter (horizontal field index).
     &,K         ! Loop counter (vertical level index).
     &,KM1       ! K-1.
     &,MBL       ! Maximum number of model layers below mixed layer top.

! Layer interface K_LOG_LAYR-1/2 is the highest which requires log
! profile correction factors to the vertical finite differences.
! The value should be reassessed if the vertical resolution is changed.
! We could set K_LOG_LAYR = BL_LEVELS and thus apply the correction
! factors for all the interfaces treated by the boundary layer scheme;
! this would be desirable theoretically but expensive computationally
! because of the use of the log function.

      INTEGER    K_LOG_LAYR
      PARAMETER (K_LOG_LAYR=2)

      IF (LTIMER) THEN
        CALL TIMER('EX_COEF ',103)
      ENDIF

!  Set MBL, "maximum number of boundary layer levels" for the purposes  
!  of mixed layer height calculation.        
                                                                        
      MBL = BL_LEVELS - 1                                               

!-----------------------------------------------------------------------
!! 0. Initialise flag for having reached top of turbulently mixed layer 
!!    and also the number of turbulently mixed layers.   
!-----------------------------------------------------------------------
                                                                        
      DO I=P1,P1+P_POINTS-1                                             
        TOPBL(I) = .FALSE.                                              
        NTML_LOCAL(I) = 1  
      ENDDO                                                             
                                                                        
!
!-----------------------------------------------------------------------
!  1.1 Loop over levels calculating local wind shears and Richardson
!      numbers.
!-----------------------------------------------------------------------
!
      DO K=2,BL_LEVELS
        KM1 = K-1                                                       
        DO I=P1,P1+P_POINTS-1
                                                                        
          DZU = U_P(I,K) - U_P(I,KM1)                                   
          DZV = V_P(I,K) - V_P(I,KM1)                                   
          DVMOD2 = MAX ( 1.0E-12 , DZU*DZU + DZV*DZV )                  
          DVDZM(I,K) = SQRT (DVMOD2) * RDZ(I,K)       
                                                                        
          RI(I,K) = DB(I,K) / ( RDZ(I,K) * DVMOD2 )      

!-----------------------------------------------------------------------
!! 1.2 If either a stable layer (Ri > 1) or the maximum boundary layer  
!!     height has been reached, set boundary layer height (ZH) to       
!!     the height of the lower boundary of the current layer   
!-----------------------------------------------------------------------
                                                                        
          IF ( .NOT.TOPBL(I) .AND. ( RI(I,K).GT.1.0 .OR. K.GT.MBL
     &         .OR. (CUMULUS(I) .AND. ZLB(I,K).GE.Z_LCL(I)) ) ) THEN
            TOPBL(I) = .TRUE.                                           
            ZH_LOCAL(I) = ZLB(I,K)                        
            NTML_LOCAL(I) = K-1    
          ENDIF                                                         
        ENDDO  ! Loop over points
      ENDDO  ! Loop over levels
!   
!-----------------------------------------------------------------------
!! 2.  Loop round "boundary" levels; calculate the stability   
!!     dependent turbulent mixing coefficients.
!-----------------------------------------------------------------------

      DO K=2,BL_LEVELS
        KM1 = K-1
        DO I=P1,P1+P_POINTS-1
!-----------------------------------------------------------------------
!! 2.1 Calculate asymptotic mixing lengths LAMBDAM and LAMBDAH
!!     (currently assumed equal).
!-----------------------------------------------------------------------

        LAMBDAM = MAX ( LAMBDA_MIN , 0.15*ZH_LOCAL(I) )     
        LAMBDAH = LAMBDAM
        LAMBDA_EFF = MAX (LAMBDAM, A_LAMBDA*H_BLEND(I) )
        IF ( K .GE. NTML_LOCAL(I) + 1) THEN       
          LAMBDAM = LAMBDA_MIN
          LAMBDAH = LAMBDA_MIN
          IF (ZLB(I,K) .GT. A_LAMBDA*H_BLEND(I)) LAMBDA_EFF = LAMBDA_MIN
        ENDIF

!-----------------------------------------------------------------------
!! 2.2 Calculate mixing lengths ELH, ELM at layer interface K-1/2.
!-----------------------------------------------------------------------

!  Incorporate log profile corrections to the vertical finite
!  differences into the definitions of ELM and ELH.
!  To save computing logarithms for all K, the values of ELM and ELH
!  are unchanged for K > K_LOG_LAYR.

          IF (K .LE. K_LOG_LAYR) THEN
            VKZ = VKMAN / RDZ(I,K)
            ELM = VKZ / ( LOG( ( ZLB(I,K) + Z0M(I) + 0.5*DZL(I,K)   ) /
     &                          ( ZLB(I,K) + Z0M(I) - 0.5*DZL(I,KM1) ) )
     &                     + VKZ / LAMBDA_EFF )
            ELH = VKZ / ( LOG( ( ZLB(I,K) + Z0M(I) + 0.5*DZL(I,K)   ) /
     &                          ( ZLB(I,K) + Z0M(I) - 0.5*DZL(I,KM1) ) )
     &                     + VKZ / LAMBDAH )
          ELSE
            VKZ = VKMAN * ( ZLB(I,K) + Z0M(I) )
            ELM = VKZ / (1.0 + VKZ/LAMBDA_EFF )
            ELH = VKZ / (1.0 + VKZ/LAMBDAH )
          ENDIF

!-----------------------------------------------------------------------
!! 2.4 Calculate stability functions FH, FM.     
!-----------------------------------------------------------------------

          IF (RI(I,K) .GE. 0.0) THEN      
            RTMRI = 0.0
            FM = 1.0 / ( 1.0 + G0 * RI(I,K) )   
            FH = 1.0 / ( 1.0 + G0 * RI(I,K) )   
!
!           !-----------------------------------------------------------
!           ! If convective cloud exists in layer K allow neutral mixing
!           ! of momentum between layers K-1 and K. This is to ensure
!           ! that a reasonable amount of momentum is mixed in the
!           ! presence of convection; it is not be required when
!           ! momentum transport is included in the convection scheme.
!           !-----------------------------------------------------------
!
            IF ( .NOT.L_MOM .AND. (CCA(I) .GT. 0.0) .AND.
     &           (K .GE. CCB(I)) .AND. (K .LT. CCT(I)) )
     &         FM = 1.0
          ELSE
            RTMRI = (ELM/ELH) * SQRT ( -RI(I,K) )     
            FM = 1.0 - ( G0*RI(I,K) / ( 1.0 + DM*RTMRI ) )  
            FH = 1.0 - ( G0*RI(I,K) / ( 1.0 + DH*RTMRI ) )      
          ENDIF

!-----------------------------------------------------------------------
!! 2.5 Calculate exchange coefficients RHO*KM, RHO*KH for interface
!!     K-1/2.
!-----------------------------------------------------------------------

          RHOKM(I,K) = RHO(I,K) * ELM * ELM * DVDZM(I,K) * FM    
          RHOKH(I,K) = RHO(I,K) * ELH * ELM * DVDZM(I,K) * FH    

        ENDDO ! p_points
      ENDDO ! bl_levels

      DO I = P1,P1+P_POINTS-1   
        IF ( NTML_LOCAL(I) .GT. NTML(I)+1 .OR. 
     &     ( NTML(I) .EQ. 1 .AND. NTML_LOCAL(I) .GT. NTML(I) ) ) THEN
             ZH(I) = ZH_LOCAL(I)   
             NTML(I) = NTML_LOCAL(I)    
        ENDIF
      ENDDO      

      IF (LTIMER) THEN
        CALL TIMER('EX_COEF ',104)
      ENDIF

      RETURN
      END
