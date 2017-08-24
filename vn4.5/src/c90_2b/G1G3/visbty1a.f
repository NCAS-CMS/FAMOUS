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
CLL  SUBROUTINE VISBTY -----------------------------------------------
CLL
CLL     PURPOSE:
CLL Process fields of temperature, specific humidity, cloud liquid 
CLL water or ThetaL and qt to give visibility in metres.
CLL Calculated at model level (eg bottom eta level 25m)
CLL or level within surface layer eg screen ht ( 1.5M )
CLL
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  Date
CLL  3.1  23/10/92  New deck author P.Smith
CLL  3.1  20/01/93  New deck - used as mods at 2.7 & 2.8/3.0.
CLL                 Interfacing done by R.T.H.Barnes.
CLL  3.2  29/04/93  CCN and derived constants moved to MODECK C_VISBTY
CLL                 Programmer Pete Clark.
CLL  3.4  07/06/94  Aerosol field introduced. Programmer Pete Clark.
CLL  4.0 05/09/95  Lower limit to aerosol introduced.  Pete Clark.
!    4.4  09/01/97  Only liquid water and not cloud ice is now used
!                   to calculate visibility. Damian Wilson.
!LL  4.5  29/04/98  Scheme replaced with NIMROD visibility diagnostic
CLL
CLL  Programming standard: U M Doc. Paper No. 4
CLL
CLL  Logical components covered :
CLL
CLL  Project task:
CLL
CLL  External documentation
CLL    Forecasting Research Scientific Paper NO.4
CLL    Diagnosis of visibility in the UK Met Office Mesoscale Model
CLL    and the use of a visibility analysis to constrain initial
CLL    conditions.  SP Ballard, BJ Wright, BW Golding    1992
!      NIMROD diagnostic: 
!      Wright, B. J., 1997: Improvements to the Nimrod Visibility
!         Analysis/Forecast System. FR-Div. Tech. Rep., No. 217.  
!      Wright, B. J., 1997: A New Visibility Analysis/Forecast System 
!         for Nimrod. Met. Office FR Tech Rep., No. 222. 
CLL
CLLEND----------------------------------------------------------------
C
C*L  Arguments:-------------------------------------------------------
      SUBROUTINE VISBTY
     +           (AK,BK,PSTAR,T,Q,QCL,QCF             !INPUT    
     &           ,AEROSOL, PROB, RHCRIT, L_MURK       !INPUT            
     +           ,P_FIELD                             !INPUT
     +           ,VISIBILITY)                         !OUTPUT
      IMPLICIT NONE
C---------------------------------------------------------------------
C Workspace usage:----------------------------------------------------
C 3 real arrays of size P_FIELD                                         
C*--------------------------------------------------------------------
C*L-------------------------------------------------------------------
C input variables-----------------------------------------------------
C---------------------------------------------------------------------
      INTEGER
     *        P_FIELD                   ! IN NO. points in field.
      REAL
     &        AK,BK                     ! IN Ak and Bk of level
     &       ,PSTAR(P_FIELD)            ! IN Surface pressure
     &       ,T(P_FIELD)                ! IN Temperature        
     &       ,Q(P_FIELD)                ! IN Qt           
     &       ,QCL(P_FIELD)              ! IN cloud water array.         
     &       ,QCF(P_FIELD)              ! IN cloud ice array.           
     &       ,AEROSOL(P_FIELD)          ! IN Aerosol mixing ratio(ug/kg)
     &       ,PROB                    ! IN Probability level ( e.g 0.5
                                      !    corresponds to median).
     &       ,RHCRIT                  ! IN Citical RH (determines
                                      !    width of distribiution)
      LOGICAL
     &   L_MURK                        ! IN : Aerosol present
C---------------------------------------------------------------------
C output variables----------------------------------------------------
C---------------------------------------------------------------------
      REAL
     &      VISIBILITY(P_FIELD)         ! OUT visibility array.
C*--------------------------------------------------------------------
C*L-------------------------------------------------------------------
C Local varables:-----------------------------------------------------
C---------------------------------------------------------------------
      REAL
     &       QT(P_FIELD)              ! total of cloud water and vapour
     &      ,P(P_FIELD)               ! pressure of level    
     &      ,Qs(P_FIELD)              ! saturation vapour pressure     
C*L  External subroutine called ----------------------------------------
      EXTERNAL QSAT_WAT                                    
C*--------------------------------------------------------------------
C constants for visibility calculation used to be set here but now
C set in MODECK.
C PI needed to set new constants.
C---------------------------------------------------------------------
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

C*L-----------COMDECK C_DENSTY FOR SUBROUTINE SF_EXCH----------
C RHOSEA    = density of sea water (kg/m3)
C RHO_WATER = density of pure water (kg/m3)
      REAL RHOSEA,RHO_WATER

      PARAMETER(RHOSEA = 1000.0,
     &          RHO_WATER = 1000.0)
C*----------------------------------------------------------------------
CLL================ COMDECK C_VISBTY ===========================
CLL Description:
CLL   This COMDECK contains declarations for constants used to diagnose
CLL visibility. Constants are set as PARAMTERs.
CLL
CLL
CLL  Model            Modification history:
CLL version  Date
CLL  3.2    29/04/93  CCN Parameters moved here from VISBTY so that
CLL                   they can also be used to compute fog fraction.
CLL                   Programmer: Pete Clark.
CLL  4.0 05/09/95  Variable AEROMAX used as upper limit to aerosol in 
CLL                assimilation introduced. Programmer Pete Clark.
CLL  4.5 01/05/98  Completely re-written for NIMROD style diagnostic.
CLL
CLLEND----------------------------------------------------------------

C Define Parameters:
      REAL
     &  N0              ! Standard number density of the aerosol (/m3)
     &, B0              ! Activation parameter
     &, radius0         ! Radius of standard aerosol particle (m)
     &, rho             ! Density of the the aerosol (Kg/m3)
     &, rho_a           ! Density of air (Kg/m3)
     &, m0              ! Standard aerosol mass mixing ratio (Kg/Kg)
     &, power           ! Aerosol particle radius/mass loading power
     &, Beta0           ! Scattering coefficient normalisation
     &, LiminalContrast ! Liminal contrast
     &, LnLiminalContrast ! Natural log of Liminal contrast
     &, VisFactor       ! Constant incorporating the scattering 
     &                  !  coefficient, normalisation and 
     &                  !  transformation to visibility
     &                  !   ( = ln(liminal contrast) / Beta0 )
     &, RecipVisAir     ! Recipirical of the clean air visibility
     &, FourThirds      ! 4/3
     &, A0              ! Constant involving surface energy of water
     &, VISFOG          ! Visibility defining fog               
     &, VISMIST         ! Visibility defining mist            
     &, AERO0           ! Minimum allowed aerosol            
     &, AEROMAX         ! maximum allowed aerosol            
      PARAMETER (
     &                    N0 = 500.0E6
     &,                   B0 = 0.5
     &,              radius0 = 0.16E-6
     &,                  rho = 1700.0
     &,                rho_a = 1.0
     &,           FourThirds = 4.0/3.0
     &,                   m0 = FourThirds * Pi 
     &                         * radius0 * radius0 * radius0
     &                         * (rho/rho_a) * N0
     &,                power = 1.0/6.0
     &,                Beta0 = 1.5 * Pi
     &,      LiminalContrast = 0.02
     &,    lnLiminalContrast = -3.912023005
     &,            VisFactor = -LnLiminalContrast / Beta0 
     &,          RecipVisAir = 1.0E-5
     &,                   A0 = 1.2E-9 
     &)
        PARAMETER (VISFOG=1000.0, VISMIST=5000.0)
        PARAMETER (AERO0=0.1)
        PARAMETER (AEROMAX=1000.0)
!
! Local parameter variables
!
      REAL
     &  OneThird        ! 1/3
     &, RHmax           ! Maximum value of relative humidity 
     &                  !  which is allowed to feed into the
     &                  !  calculation of the 'fog' droplet radius
     &, RHmin           ! Minimum value of relative humidity which 
     &                  !  is allowed to feed into the calculation
     &                  !   of the 'fog' droplet radius
     &, Weight          ! Weighting on new value for iterative 
     &                  !  solution of droplet radius
     &, Delta_radius_star! Convergence required for iterative 
     &                  !  solution of droplet radius
     &, N               ! Local number density
     &, qt_limit        ! Smallest Qt value allowed
     &, radius_star_min !
     &, radius_star_max !
     &, radius_star_factor!  

      INTEGER
     &  Niterations     !  Maximum number of iteration used to  
     &                  !   estimate the water droplet radius
      PARAMETER ( OneThird = 1.0/3.0
     &,              RHmin = 0.001
     &,              RHmax = 0.99
     &,             Weight = 0.75
     &,  Delta_radius_star = 0.001
     &,        Niterations = 20 
     &,           qt_limit = 0.0001
     &,    radius_star_min = 1.0
     &,    radius_star_max = 1000.0
     &, radius_star_factor = 4.0 )
!
! Local workspace variables    
!
       INTEGER
     &  Point            !  Loop variable for points    
     &, Iteration        !  Loop variable iterations used to estimate
     &                   !   the water droplet radius

      REAL
     &  m_over_m0        !  Ratio of  aerosol mass mixing ratio and
     &                   !   the standard aerosol mass mixing ratio
     &, RecipVis         !  Recipirical of the visibility
     &, radius_dry       !  Radius of dry aerosol particle (m)
     &, radius           !  Radius of fog droplets (m)
     &, radius_star1     !  Previous estimate of water droplet radius
     &                   !   divided by the dry radius
     &, radius_star2     !  Current best estimate of water droplet 
     &                   !   radius divided by the dry radius
     &, radius_act       !  Activation droplet radius
     &, radius_star_act  !  Activation droplet rad divided by dry rad
     &, A                !  A0 divided by the dry radius
     &, RH_lim           !  Limited RH value (fractional)
     &, Fn               !  Value of droplet radius function
     &, Deriv            !  Derivative of droplet radius function
     &, radius_star_diff !  Absolute value of radius_star1 minus 
     &                   !    radius_star2
     &, RHterm           !  Relative humidity term in function to be 
     &                   !   minimised to find the droplet radius
     &, qLterm           !  Liquid water term in function to be 
     &                   !   minimised to find the droplet radius
     &, RHderiv          !  Derivative of relative humidity term 
     &, qLderiv          !  Derivative of liquid water term
     &, bs               !  Width of distribution in total water 
     &                   !   mixing ratio space (kg/kg)
     &, qt_mod           !  Modified total water value based on the 
     &                   !   probability of the value occurring 
     &                   !   assuming a triangular distriubtion
     &                   !   of width bs.
     &, qt_mod_factor    !  Factor to multiply bs to modify qt
!     Check Prob is legal
      IF ( Prob .LT. 0.0 .OR.  Prob .GT. 1.0 ) THEN
        Write(6,*)"INVALID PROBABILITY VALUE in VISBTY",Prob
        Prob=MIN(MAX(Prob,0.0),1.0)
      ENDIF
!     Create factor to multiply bs by to modify qt
      IF ( Prob .EQ. 0.5 ) THEN
        qt_mod_factor = 0.0
      ELSE IF ( Prob .GE. 0.0 .AND. Prob .LT. 0.5 ) THEN 
        qt_mod_factor = ( 1.0 - SQRT( 2.0 * Prob ) )
      ELSE IF ( Prob .GE. 0.5 .AND. Prob .LE. 1.0 ) THEN 
        qt_mod_factor = - ( 1.0 - SQRT( 2.0 * (1.0-Prob) ) )
      END IF
! ----------------------------------------------------------------------
! For the new cloud and precipitation scheme only use the liquid content
! 1. Calculate total of water vapour and liquid water contents, P and
! limit aerosol       
! ----------------------------------------------------------------------
      DO Point=1,P_FIELD                                                
        QT(Point) = Q(Point)+QCL(Point)                              
      END DO
   
!     Make sure aerosol greater than lower bound
      DO Point=1,P_FIELD                                                
        AEROSOL(Point)=MAX(AEROSOL(Point),AERO0)                        
      ENDDO

!     Calculate pressure
      DO Point=1,P_FIELD
        P(Point)=AK+BK*PSTAR(Point)
      ENDDO

      CALL QSAT_WAT (QS,T,P,P_FIELD)

      DO Point = 1 , P_FIELD
         
!-------------------------------------------------------------------
!* 2. Calculate the ratio of the aerosol mass mixing ratio to the 
!*    standard mass mixing ratio, m_over_m0, and the aerosol number 
!*    density, N, the dry radius, radius_dry:
!*                      p 
!*                  (m )  
!*           r = r0 (--) 
!*            d     (m0)
!*
!*
!*    And the activation radius:
!*
!*                             1/2
!*                  (       3 )  
!*                  ( 3 B0 r  )  
!*           r    = ( ------d-) 
!*            act   (   A0    )
!*
!*    and A (A0 divided by the dry radius).
! N.B. AEROSOL is in ug/kg, m in kg/kg
! If not available, use 10 ug/kg
!-------------------------------------------------------------------

        if (L_MURK) then  
          m_over_m0 = max(Aerosol(Point)/m0*1.0E-9, 0.0001)
        else
          m_over_m0 = max(10.0/m0*1.0E-9, 0.0001)
        endif

        N = N0 * m_over_m0**(1.0-3*power)
      
        radius_dry = radius0 * (m_over_m0)**power
        A = A0 / radius_dry

        radius_act = SQRT( (3 * B0 * radius_dry**3) / A0 )
        radius_star_act =  radius_act/radius_dry 

!-------------------------------------------------------------
!* 3. Calculate the width of the total water
!*    distribution and a modified value of total water, based
!*    on a probability.
!-------------------------------------------------------------

        bs = (1.0-RHcrit) * qs(Point)

        qt_mod = MAX( qt_limit, qt(Point)+ qt_mod_factor* bs)


!====================================================================
!* 4.  Use Newton-Raphson to iteratively improve on a first-guess
!*     droplet radius, using the droplet growth equation and the
!*     geometric relation between liquid water and droplet radius.
!====================================================================
!* 4.1 Calculate a first guess relative humidity, qt/qs, but limit it
!*     to be in the range 0.001 -> 0.999.
!*     From this calculate a first-guess normalised radius using a
!*     simplified version of the droplet growth equation:
!*
!*                              1/3
!*                (       B0   )
!*           r  = ( 1 - ------ )
!*            *   (     ln(RH) )
!*
!----------------------------------------------------------------------

        RH_lim = MIN( MAX( qt_mod/qs(Point), RHmin ) , RHmax )
        radius_star2 = (1.0-B0/LOG(RH_lim))**OneThird 

!----------------------------------------------------------------------
!* 4.2 Initialise the iteration counter, the normalised radius 
!*     difference, and the updated normalised radius value.
!----------------------------------------------------------------------

        Iteration = 0
        radius_star_diff = 1.0
        radius_star1 = radius_star2

        Do While ( Iteration .LT. Niterations .AND. 
     &               radius_star_diff .GT. Delta_radius_star )

!----------------------------------------------------------------------
!* 4.3 Update the iteration counter and the normalised radius value.
!----------------------------------------------------------------------

          Iteration = Iteration + 1 
          radius_star1 = Weight * radius_star2 
     &                 + ( 1.0 - Weight ) * radius_star1

!----------------------------------------------------------------------
!* 4.4 Calculate the relative humidity term:
!*
!*                      ( A        B0   )
!*          RHterm = exp( --  -  ------ )
!*                      ( r       3     )
!*                      (  *     r  - 1 )
!*                      (         *     )
!*
!*      and its derivative with respect to the normalised radius:
!*
!*                    (                 2    )
!*                    (   A       3 B0 r     )
!*          RHderiv = ( - --  +  -------*- 2 ) * RHterm
!*                    (    2     (  3     )  )
!*                    (   r      ( r  - 1 )  )
!*                    (    *     (  *     )  )
!*
!----------------------------------------------------------------------

          If ( radius_star1 .LT. radius_star_act ) then
            RHterm  = EXP( A/radius_star1 
     &                     - B0/(radius_star1**3-1.0) )* qs(Point)
            RHderiv = - RHterm * ( -A/(radius_star1**2)
     &                + (3.0*B0*radius_star1**2)
     &                /(radius_star1**3-1.0)**2 )
          Else
            RHterm  = EXP( A/radius_star_act
     &                     - B0/(radius_star_act**3-1.0) ) * qs(Point)
            RHderiv = 0.0
          Endif


!----------------------------------------------------------------------
!* 4.5 Calculate the liquid water mixing ratio term:
!*
!*                                          
!*                   4             3 (  3     )
!*          qLterm = - Pi rho_w N r  ( r  - 1 )
!*                   3             d (  *     )       
!*
!*      and its derivative with respect to the normalised radius:
!*                                          
!*                                  3  2
!*          qLderiv = 4 Pi rho_w N r  r 
!*                                  d  * 
!*      
!----------------------------------------------------------------------

          qLterm  = N * FourThirds * Pi * RHO_WATER * radius_dry**3
     &                * ( radius_star1**3 - 1.0 )
          qLderiv  = - N * 4.0 * Pi * RHO_WATER 
     &                * radius_dry**3 * radius_star1**2 

!----------------------------------------------------------------------
!* 4.6 Calculate the function, Fn, and its derivative, Deriv, and 
!*     an improved estimate of the normalised radius, 
!*     using Newton Raphson:
!*
!*          Fn = qt - RHterm - qLterm
!*
!*          Deriv = RHderiv + qLderiv
!*      
!*                          Fn
!*          r      = r  -  -----
!*           * new    *    Deriv 
!* 
!*     The new estimate of the normalised radius is limited lie between
!*     prescribed maximum and minimum values and within a factor of the
!*     previous value to ensure that the soultion does not diverge.     
!----------------------------------------------------------------------

          Fn    = qt_mod - RHterm - qLterm
          Deriv = RHderiv + qLderiv

          radius_star2 = radius_star1 - Fn/Deriv

          IF ( radius_star2 .LT. radius_star_min )  
     &        radius_star2 = radius_star_min
          IF ( radius_star2 .GT. radius_star_max )  
     &        radius_star2 = radius_star_max
          IF ( radius_star2 .GT. radius_star_factor * radius_star1 )  
     &        radius_star2 = radius_star_factor * radius_star1
          IF ( radius_star2 .LT. radius_star1 / radius_star_factor )  
     &        radius_star2 = radius_star1 / radius_star_factor
 
!---------------------------------------------------------------------
!* 4.7 Calculate difference between the old and the new values of the 
!*     normalised radius.
!---------------------------------------------------------------------

          radius_star_diff = ABS( radius_star1 - radius_star2 )

        END DO

!---------------------------------------------------------------------
!* 5.  Calculate the radius from the final normalised radius.
!---------------------------------------------------------------------
 
        radius = radius_star2 * radius_dry

!---------------------------------------------------------------------
!* 6. Calculate the visibility, Vis, using the equation:
!*
!*                 ln(liminal contrast)
!*           Vis = -------------2------
!*                     Beta0 N r
!*
!*    (An extra term RecipVisAir is included in the recipical of 
!*     visibility to limit visibilities to 100km in clean air).
!---------------------------------------------------------------------

        RecipVis = (N * radius**2) / VisFactor + RecipVisAir
        Visibility(Point) = 1/RecipVis

      END DO
  
      RETURN
      END
C
CLL  SUBROUTINE VISTOQT ----------------------------------------------- 
CLL
CLL     PURPOSE:
CLL Invert relationship between aerosol, visibility and water
CLL content   
CLL This is needed for fog probability calculation.
!             Since 4.5, adopted NIMROD based code:
CLL
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  Date
CLL  3.4  07/06/94  First written. Programmer Pete Clark.
CLL  4.0 05/09/95  Diagnosed equivalent RH constrained to be >1%. 
CLL                Lower limit to aerosol introduced.  Pete Clark.
!    4.5  30/04/98 NIMROD code adopted. 
!                   Pete Clark responsible for UM implementation of
!                   Bruce Wright's NIMROD code.
CLL
CLL  Programming standard: Unified Model Documentation Paper No 3,
CLL                        Version 7, dated 11/3/93.
CLL
CLL  Logical components covered :
CLL
CLL  Project task:
CLL
CLL  External documentation
CLL    Forecasting Research Scientific Paper NO.4
CLL    Diagnosis of visibility in the UK Met Office Mesoscale Model
CLL    and the use of a visibility analysis to constrain initial
CLL    conditions.  SP Ballard, BJ Wright, BW Golding    1992
!      Wright, B. J., 1997: Improvements to the Nimrod Visibility
!         Analysis/Forecast System. FR-Div. Tech. Rep., No. 217.  
!      Wright, B. J., 1997: A New Visibility Analysis/Forecast System 
!         for Nimrod. Met. Office FR Tech Rep., No. 222. 
CLL
CLLEND----------------------------------------------------------------
C
C*L  Arguments:-------------------------------------------------------
      SUBROUTINE VISTOQT                                                
     &           (VISIBILITY                          !INPUT
     &           ,Qs                                  !INPUT
     &           ,AEROSOL                             !INPUT
     &           ,L_MURK                              !INPUT
     &           ,Npoints                             !INPUT
     &           ,qt )                                !OUTPUT  
      IMPLICIT NONE
C---------------------------------------------------------------------
C Workspace usage:----------------------------------------------------
C None
C*--------------------------------------------------------------------
C*L-------------------------------------------------------------------
C input variables-----------------------------------------------------
C---------------------------------------------------------------------
      INTEGER
     &        Npoints             ! IN NO. points in field.            
      REAL
     &        VISIBILITY          ! IN visibility 
                                  !   NB Original code had Visibility
                                  !      as an array - this feature was
                                  !      not used so has been removed   
     &       ,Qs(Npoints)         !  Saturated humidity mixing ratio
     &       ,AEROSOL(Npoints)    ! IN Aerosol mixing ratio(ug/kg)     
      LOGICAL
     &        L_MURK                    ! IN : Aerosol present
C---------------------------------------------------------------------
C output variables----------------------------------------------------
C---------------------------------------------------------------------
      REAL
     &        qt(Npoints)         ! OUT Total water mixing ratio (kg/kg)
C*--------------------------------------------------------------------
C*L-------------------------------------------------------------------
C---------------------------------------------------------------------
C*--------------------------------------------------------------------
C constants for visibility calculation used to be set here but now
C set in COMDECK.
C PI needed to set new constants.
C---------------------------------------------------------------------
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

C*L-----------COMDECK C_DENSTY FOR SUBROUTINE SF_EXCH----------
C RHOSEA    = density of sea water (kg/m3)
C RHO_WATER = density of pure water (kg/m3)
      REAL RHOSEA,RHO_WATER

      PARAMETER(RHOSEA = 1000.0,
     &          RHO_WATER = 1000.0)
C*----------------------------------------------------------------------
CLL================ COMDECK C_VISBTY ===========================
CLL Description:
CLL   This COMDECK contains declarations for constants used to diagnose
CLL visibility. Constants are set as PARAMTERs.
CLL
CLL
CLL  Model            Modification history:
CLL version  Date
CLL  3.2    29/04/93  CCN Parameters moved here from VISBTY so that
CLL                   they can also be used to compute fog fraction.
CLL                   Programmer: Pete Clark.
CLL  4.0 05/09/95  Variable AEROMAX used as upper limit to aerosol in 
CLL                assimilation introduced. Programmer Pete Clark.
CLL  4.5 01/05/98  Completely re-written for NIMROD style diagnostic.
CLL
CLLEND----------------------------------------------------------------

C Define Parameters:
      REAL
     &  N0              ! Standard number density of the aerosol (/m3)
     &, B0              ! Activation parameter
     &, radius0         ! Radius of standard aerosol particle (m)
     &, rho             ! Density of the the aerosol (Kg/m3)
     &, rho_a           ! Density of air (Kg/m3)
     &, m0              ! Standard aerosol mass mixing ratio (Kg/Kg)
     &, power           ! Aerosol particle radius/mass loading power
     &, Beta0           ! Scattering coefficient normalisation
     &, LiminalContrast ! Liminal contrast
     &, LnLiminalContrast ! Natural log of Liminal contrast
     &, VisFactor       ! Constant incorporating the scattering 
     &                  !  coefficient, normalisation and 
     &                  !  transformation to visibility
     &                  !   ( = ln(liminal contrast) / Beta0 )
     &, RecipVisAir     ! Recipirical of the clean air visibility
     &, FourThirds      ! 4/3
     &, A0              ! Constant involving surface energy of water
     &, VISFOG          ! Visibility defining fog               
     &, VISMIST         ! Visibility defining mist            
     &, AERO0           ! Minimum allowed aerosol            
     &, AEROMAX         ! maximum allowed aerosol            
      PARAMETER (
     &                    N0 = 500.0E6
     &,                   B0 = 0.5
     &,              radius0 = 0.16E-6
     &,                  rho = 1700.0
     &,                rho_a = 1.0
     &,           FourThirds = 4.0/3.0
     &,                   m0 = FourThirds * Pi 
     &                         * radius0 * radius0 * radius0
     &                         * (rho/rho_a) * N0
     &,                power = 1.0/6.0
     &,                Beta0 = 1.5 * Pi
     &,      LiminalContrast = 0.02
     &,    lnLiminalContrast = -3.912023005
     &,            VisFactor = -LnLiminalContrast / Beta0 
     &,          RecipVisAir = 1.0E-5
     &,                   A0 = 1.2E-9 
     &)
        PARAMETER (VISFOG=1000.0, VISMIST=5000.0)
        PARAMETER (AERO0=0.1)
        PARAMETER (AEROMAX=1000.0)
!
! Local parameter variables
!
      REAL
     &  qt_limit       ! Smallest total water mixing ratio value allowed
      PARAMETER (   qt_limit = 0.0001 )
!
! Local workspace variables    
!
       INTEGER
     &  Point            !  Loop variable for points    

      REAL
     &  qL               !  Liquid water mixing ratio (Kg/Kg).
     &, radius_dry       !  Dry particle radius for aerosol (m)
     &, radius           !  Radius of fog droplets (m)
     &, radius_star      !  Water droplet radius divided by dry radius
     &, radius_act       !  Activation droplet radius
     &, radius_star_act  !  Activation droplet radius divided by the dry
     &, radius_star_used !  Water droplet radius divided by the dry 
     &                   !   radius actually used for the relative 
     &                   !   humidity calculation
     &, RH               !  Relative humidity derived from visibility
     &, A                !  A0 divided by the dry radius
     &, m_over_m0        !  Ratio of the aerosol mass mixing ratio and
     &                   !   the standard aerosol mass mixing ratio
     &, N                !  Number density of aerosol particles (/m3)
  
  
!**

      Do Point = 1 , Npoints

!---------------------------------------------------------------------
!* 1.  Calculate the ratio of the aerosol mass mixing ratio to the 
!*     standard mass mixing ratio, m_over_m0, and the aerosol number 
!*     density, N:
!*
!*                      (1-3p)
!*                  (m ) 
!*           N = N0 (--) 
!*                  (m0)
!* 
!*     And the dry radius, radius_dry:
!*                      p 
!*                  (m )  
!*           r = r0 (--) 
!*            d     (m0)
!*
!*     And A (A0 divided by the dry radius).
!*
!*     And the activation radius:
!*
!*                             1/2
!*                  (       3 )  
!*                  ( 3 B0 r  )  
!*           r    = ( ------d-) 
!*            act   (   A0    )
!*
! N.B. AEROSOL is in ug/kg, m in kg/kg
! If not available, use 10 ug/kg
!-------------------------------------------------------------------

        if (L_MURK) then  
          m_over_m0 = max(Aerosol(Point)/m0*1.0E-9, 0.0001)
        else
          m_over_m0 = max(10.0/m0*1.0E-9, 0.0001)
        endif

        N = N0 * m_over_m0**(1.0-3*power)
        
        radius_dry = radius0 * (m_over_m0)**power
        A = A0 / radius_dry

        radius_act = SQRT( (3 * B0 * radius_dry**3) / A0 )
        radius_star_act =  radius_act/radius_dry 
      
!----------------------------------------------------------------------
!* 2.  Calculate a water droplet radius, from the visibility:
!*
!*                                1/2
!*               ( ln( epsilon ) ) 
!*           r = (---------------)
!*               (  Vis N Beta0  )  
!*
!*    (An extra term RecipVisAir is included in the recipical of 
!*     visibility to limit visibilities to 100km in clean air).
!----------------------------------------------------------------------

        radius = SQRT( (VisFactor/N) * 
     &           ((1.0/Visibility) - RecipVisAir) )

!----------------------------------------------------------------------
!* 3.  Provided the diagnosed radius is greater than the dry radius,
!*     calculate the normalised droplet radius, and the saturated 
!*     humidity mixing ratio.
!----------------------------------------------------------------------

        If ( radius .GT. radius_dry ) then

          radius_star = radius / radius_dry

!----------------------------------------------------------------------
!* 5.  Calculate the corresponding liquid water mixing ratio:
!*
!*                                         
!*               4            (  3     3 )
!*          qL = - Pi rho_w N ( r  - r   )
!*               3            (       d  )       
!*
!----------------------------------------------------------------------

          qL = FourThirds * Pi * rho_water * N  * 
     &         ( radius**3 - radius_dry**3 )

!----------------------------------------------------------------------
!* 6.  Calculate the relative humidity:
!*
!*                  ( A        B0   )
!*          RH = exp( --  -  ------ )
!*                  ( r       3     )
!*                  (  *     r  - 1 )
!*                  (         *     )
!*
!----------------------------------------------------------------------

          If ( radius_star .LT. radius_star_act ) then
            RH = EXP( A/radius_star
     &                - B0 /( radius_star **3 - 1.0 ) )
          Else
            RH = EXP( A/radius_star_act  
     &                  - B0 /( radius_star_act **3 - 1.0 ) )
          Endif

!----------------------------------------------------------------------
!* 7.  Calculate the total water mixing ratio:  qt = RH * qs(T) + qL
!----------------------------------------------------------------------

          qt(Point) = MAX( RH * qs(Point) + qL, qt_limit )
        
!----------------------------------------------------------------------
!* 8. If the droplet radius is less than the dry radius, then set the 
!*    total water mixing ratio to the minimum value.
!----------------------------------------------------------------------

        Else

          qt(Point) = qt_limit

        End if


      End Do

      RETURN
      END
C
CLL  SUBROUTINE FOG_FR------------------------------------------------
CLL
CLL  Purpose: Calculates fog fraction, using the large scale cloud
CLL           scheme. The fog fraction is similar to the cloud fraction
CLL           except it records the fraction of a grid box with RH
CLL           greater than that required for the critical visibility
CLL           (e.g 1 km).
!LL           Since 4.5, adopted NIMROD based code:
!LL           Calculates the fraction of a gridsquare with visibility 
!LL           less than threshold, Vis_thresh, given the total water
!LL           mixing ratio, qt, temperature, T, pressure p, and the 
!LL           aerosol mass mixing ratio, m, assuming a triangular
!LL           distribution of states about the median, characterised by 
!LL           a critical relative humdity value, RHcrit.
CLL           NB:  Throughout, levels are counted from the bottom up,
CLL           i.e. the lowest level under consideration is level 1, the
CLL           next lowest level 2, and so on.
CLL
CLL           Suitable for single-column use.
CLL
CLL   Model         Modification history:
CLL  version  Date
CLL
CLL     3.2 04/05/93 Created by Pete Clark
!LL     3.4 04/08/95 LS_CLD replaced by GLUE_CLD. Andrew Bushell.
CLL   4.2    Oct. 96  T3E migration: *DEF CRAY removed (dynamic
CLL                    allocation now unconditional)
CLL                                   S.J.Swarbrick
!       4.4 01/07/97 Calculation is now based on liquid water
!                    and not on ice. Calculates liquid fog fraction.
!                    This scheme is severely tied to the ideas behind
!                    the 1A cloud scheme. Damian Wilson.
!LL     4.5 30/04/98 NIMROD code adopted. Calculation is still based
!LL                  on liquid water and not ice, but arguments changed
!LL                  to pass T, q, qcl and qcf separately. This is to
!LL                  make code independent of cloud scheme and capable
!LL                  of future development to include ice properly.
!LL                  Pete Clark responsible for UM implementation of
!LL                  Bruce Wright's NIMROD code.
CLL
CLL Programming standard:  Unified Model Documentation Paper No 3,
CLL                        Version 5, dated 08/12/92.
CLL
CLL Documentation:  
!LL    Wright, B. J., 1997: Improvements to the Nimrod Visibility
!LL       Analysis/Forecast System. FR-Div. Tech. Rep., No. 217.  
!LL    Wright, B. J., 1997: A New Visibility Analysis/Forecast System 
!LL       for Nimrod. Met. Office FR Tech Rep., No. 222. 
CLL
CLLEND----------------------------------------------------------------
C
C*L
C*LArguments:---------------------------------------------------------
      SUBROUTINE FOG_FR(
     + AK,BK,PSTAR,RHCRIT,LEVELS,POINTS,PFIELD,
     & T,AEROSOL,L_MURK,Q,QCL,QCF,VIS,FF,NVIS,                          
     & ERROR
     +)
      IMPLICIT NONE
      INTEGER
     + LEVELS              ! IN No. of levels being processed.
     +,POINTS              ! IN No. of gridpoints being processed.
     +,PFIELD              ! IN No. of points in global field (at one
C                          !    vertical level).
     &,NVIS                ! IN No. of visibility thresholds
      REAL
     + PSTAR(PFIELD)       ! IN Surface pressure (Pa).
     +,RHCRIT(LEVELS)      ! IN Critical relative humidity.  See the
C                          !    the paragraph incorporating eqs P292.11
C                          !    to P292.14; the values need to be tuned
C                          !    for the given set of levels.
     +,AK(LEVELS)          ! IN Hybrid "A" co-ordinate.
     +,BK(LEVELS)          ! IN Hybrid "B" co-ordinate.
     +,Q(PFIELD,LEVELS)    ! IN Specific Humidity  
C                          !    (kg per kg air).
     &,QCL(PFIELD,LEVELS)  ! Cloud liquid water content at              
C                          !     processed levels (kg per kg air).      
     &,QCF(PFIELD,LEVELS)  ! Cloud ice content at processed levels      
C                          !    (kg per kg air).                        
     &,T(PFIELD,LEVELS)    ! IN Temperature (K).                     
     &,AEROSOL(PFIELD,LEVELS) ! IN Aerosol mixing ratio(ug/kg)
     &,VIS(NVIS)              ! Visibility thresholds                   
      LOGICAL
     &   L_MURK               ! IN : Aerosol present

      REAL
     + FF(PFIELD,LEVELS,NVIS)   ! OUT Vis prob at processed levels      
C                          !     (decimal fraction).
      INTEGER ERROR        ! OUT 0 if OK; 1 if bad arguments.
C
C*--------------------------------------------------------------------
C*L  Workspace usage----------------------------------------------------
      REAL                 ! "Automatic" arrays on Cray.
     & P(POINTS)                                                        
     &,QT(POINTS)          ! total of cloud water and vapour       
     &,QS(POINTS)          ! Saturated spec humidity for temp T
     &,qt_thresh(POINTS)   ! modified qt 
     &,bs      
C*L  External subroutine called ----------------------------------------
      EXTERNAL QSAT_WAT,VISTOQT                                   
C* Local, including SAVE'd, storage------------------------------------
C
      INTEGER K,I,J     ! Loop counters: K - vertical level index.      
C                       !                I - horizontal field index.
                        !                J - Vis threshold index.  
C*--------------------------------------------------------------------
C*  Local and other physical constants----------------------------------
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

CLL================ COMDECK C_VISBTY ===========================
CLL Description:
CLL   This COMDECK contains declarations for constants used to diagnose
CLL visibility. Constants are set as PARAMTERs.
CLL
CLL
CLL  Model            Modification history:
CLL version  Date
CLL  3.2    29/04/93  CCN Parameters moved here from VISBTY so that
CLL                   they can also be used to compute fog fraction.
CLL                   Programmer: Pete Clark.
CLL  4.0 05/09/95  Variable AEROMAX used as upper limit to aerosol in 
CLL                assimilation introduced. Programmer Pete Clark.
CLL  4.5 01/05/98  Completely re-written for NIMROD style diagnostic.
CLL
CLLEND----------------------------------------------------------------

C Define Parameters:
      REAL
     &  N0              ! Standard number density of the aerosol (/m3)
     &, B0              ! Activation parameter
     &, radius0         ! Radius of standard aerosol particle (m)
     &, rho             ! Density of the the aerosol (Kg/m3)
     &, rho_a           ! Density of air (Kg/m3)
     &, m0              ! Standard aerosol mass mixing ratio (Kg/Kg)
     &, power           ! Aerosol particle radius/mass loading power
     &, Beta0           ! Scattering coefficient normalisation
     &, LiminalContrast ! Liminal contrast
     &, LnLiminalContrast ! Natural log of Liminal contrast
     &, VisFactor       ! Constant incorporating the scattering 
     &                  !  coefficient, normalisation and 
     &                  !  transformation to visibility
     &                  !   ( = ln(liminal contrast) / Beta0 )
     &, RecipVisAir     ! Recipirical of the clean air visibility
     &, FourThirds      ! 4/3
     &, A0              ! Constant involving surface energy of water
     &, VISFOG          ! Visibility defining fog               
     &, VISMIST         ! Visibility defining mist            
     &, AERO0           ! Minimum allowed aerosol            
     &, AEROMAX         ! maximum allowed aerosol            
      PARAMETER (
     &                    N0 = 500.0E6
     &,                   B0 = 0.5
     &,              radius0 = 0.16E-6
     &,                  rho = 1700.0
     &,                rho_a = 1.0
     &,           FourThirds = 4.0/3.0
     &,                   m0 = FourThirds * Pi 
     &                         * radius0 * radius0 * radius0
     &                         * (rho/rho_a) * N0
     &,                power = 1.0/6.0
     &,                Beta0 = 1.5 * Pi
     &,      LiminalContrast = 0.02
     &,    lnLiminalContrast = -3.912023005
     &,            VisFactor = -LnLiminalContrast / Beta0 
     &,          RecipVisAir = 1.0E-5
     &,                   A0 = 1.2E-9 
     &)
        PARAMETER (VISFOG=1000.0, VISMIST=5000.0)
        PARAMETER (AERO0=0.1)
        PARAMETER (AEROMAX=1000.0)
!
!  (a) Scalars effectively expanded to workspace by the Cray (using
!      vector registers).
C-----------------------------------------------------------------------
C  Check input arguments for potential over-writing problems.
C-----------------------------------------------------------------------
      ERROR=0
      IF(POINTS.GT.PFIELD)THEN
        ERROR=1
        GOTO1000
      ENDIF
C
C
C-----------------------------------------------------------------------
CL Subroutine structure :
CL Loop round levels to be processed.
C-----------------------------------------------------------------------
C
      DO K=1,LEVELS
C
C-----------------------------------------------------------------------
CL 1. Calculate Pressure and initialise temporary arrays
C-----------------------------------------------------------------------
C
        DO I=1,POINTS
          P(I)=AK(K)+PSTAR(I)*BK(K)
          QT(I)=Q(I,K)+QCL(I,K) 
        ENDDO ! Loop over points

!-----------------------------------------------------------------------
!* 2.  Calculate total water threshold corresponding to visibility 
!      Since Qs is needed more than once, pre-calculate and pass it
!-----------------------------------------------------------------------

        CALL QSAT_WAT (Qs,T(1,K),P,POINTS)

        DO J=1,NVIS

          Call VISTOQT( VIS(J), Qs, AEROSOL(1,K), L_MURK, 
     &                points, qt_thresh )


!-----------------------------------------------------------------------
!* 3.  Calculate the width of the distribution in total water space, bs:
!*
!*           bs = ( 1 - RHcrit ) * qs(T)
!*
!-----------------------------------------------------------------------

          Do I = 1 , points

            bs = (1.0-RHcrit(K)) * qs(I)

!=======================================================================
!* 4.  Calculate the fraction of states in a triangular
!*     distribution which exceed the total water threshold.
!=======================================================================

!-----------------------------------------------------------------------
!* 4.1 If total water threshold value is less than the total water value
!*     minus the width of the distribution, then all of the states have 
!*     a total water value exceeding the threshold, so set the 
!*     visibility fraction to 1.0
!-----------------------------------------------------------------------

            if ( qt_thresh(I) .LE. qt(I)-bs ) then

              FF(I,K,J) = 1.0

!-----------------------------------------------------------------------
!* 4.2 If total water threshold value is greater than the total water 
!*     value minus the width of the distribution, but less than the 
!*     total water value then the visibility fraction, VF, is given by:
!*    
!*                                                    2
!*                             ( qt       - qt + bs  )
!*            VF = 1.0 - 0.5 * (    thresh           )
!*                             ( ------------------- )
!*                             (          bs         )
!*
!-----------------------------------------------------------------------

             Else if ( qt_thresh(I) .GT. qt(I)-bs .AND.
     &                 qt_thresh(I) .LE. qt(I) ) then

               FF(I,K,J) = 1.0 - 0.5 * 
     &              (( qt_thresh(I) - qt(I) + bs )/ bs)**2

!-----------------------------------------------------------------------
!* 4.3 If total water threshold value is greater than the total water 
!*     value, but less than the total water value plus the width of the 
!*     distribution, then the visibility fraction, VF, is given by:
!*
!*                                              2
!*                       ( qt + bs - qt        )
!*            VF = 0.5 * (             thresh  )
!*                       ( ------------------- )
!*                       (          bs         )
!*
!-----------------------------------------------------------------------

             Else if ( qt_thresh(I) .GT. qt(I) .AND.   
     &                 qt_thresh(I) .LE. qt(I)+bs    ) then

                FF(I,K,J)= 0.5 * (( qt(I) + bs - qt_thresh(I))/bs)**2

!-----------------------------------------------------------------------
!* 4.4 If total water threshold value is greater than the total water 
!*     value plus the width of the distribution, then non of the states 
!*     have a total water value exceeding the threshold, so set the 
!*     visibility fraction to 0.0
!-----------------------------------------------------------------------

             Else

               FF(I,K,J) = 0.0

            End if

          End Do ! Loop over Points I

        End Do ! Loop over VIS J

      ENDDO ! Loop over levels
C
 1000 CONTINUE ! Error exit
      RETURN
      END
