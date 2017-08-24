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
!  SUBROUTINE LSP_ICE------------------------------------------------
!
!  Purpose: Form or augment ice at the expense of cloud water or
!           vapour in one model layer.
!           Also perform flux divergence of falling ice and rain,
!           Evaporation and melting of snow,
!           Evaporation of rain, formation of rain.
!           This is the principal subroutine of the 3B large scale
!           precipitation scheme.
!
! S Ballard   <- programmer
! D Wilson    <- programmer
!
!  Model            Modification history from model version 4.5:
! version  Date
!  4.5     Feb 98  New deck                      Damian Wilson
!
!  Programming standard: Unified Model Documentation Paper No 4,
!                        Version 1dated  12/9/89.
!
!  Logical component covered: Part of P26.
!
!  System task:
!
!  Documentation: Unified Model Documentation Paper No 26.
!
!  Arguments:-----------------------------------------------------------
      SUBROUTINE LSP_ICE(
     &  P,RHODZ,TIMESTEPFIXED,POINTS,
     &  RHCPT,
     &  SO4_ACC,SO4_DIS,
     &  QCF,QCL,Q,RAIN,SNOW,VF,
     &  T,CFLIQ,CFICE,BLAND,CX,CONSTP)
      IMPLICIT NONE
!
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

C*L------------------COMDECK C_O_DG_C-----------------------------------
C ZERODEGC IS CONVERSION BETWEEN DEGREES CELSIUS AND KELVIN
C TFS IS TEMPERATURE AT WHICH SEA WATER FREEZES
C TM IS TEMPERATURE AT WHICH FRESH WATER FREEZES AND ICE MELTS
      REAL ZERODEGC,TFS,TM

      PARAMETER(ZERODEGC=273.15,
     &          TFS=271.35,
     &          TM=273.15)
C*----------------------------------------------------------------------

C*L------------------COMDECK C_LHEAT------------------------------------
C LC IS LATENT HEAT OF CONDENSATION OF WATER AT 0DEGC
C LF IS LATENT HEAT OF FUSION AT 0DEGC
      REAL LC,LF

      PARAMETER(LC=2.501E6,
     &          LF=0.334E6)
C*----------------------------------------------------------------------
! --------------------------COMDECK C_LSPMIC----------------------------
! SPECIFIES MICROPHYSICAL PARAMETERS FOR AUTOCONVERSION, HALLETT MOSSOP
! PROCESS, ICE NUCLEATION. ALSO SPECIFIES NUMBER OF ITERATIONS OF
! THE MICROPHYSICS.
! ----------------------------------------------------------------------
!
! input variables
!
        INTEGER               !, INTENT(IN)      
     &    LSITER
!                 Number of iterations in microphysics.
     &,   ADV_TYPE
!                 Vertical advection method.
!
        REAL                  !, INTENT(IN)
     &    AUTOLIM_SEA          
!           Liquid Water limit for autoconversion over sea
     &,   AUTOLIM_LAND          
!           Liquid Water limit for autoconversion over land
     &,   AUTORATE_SEA
!           Rate constant for autoconversion over sea
     &,   AUTORATE_LAND
!           Rate constant for autoconversion over land
     &,   M0               
!           Nucleation mass
     &,   QCFMIN
!           Minimum allowed QCF after microphysics
     &,   TNUC             
!           Maximum Temp for ice nuclei nucleation (deg C)
     &,   THOMO
!           Maximum Temp for homogenous nucleation (deg C)
     &,   HM_T_MIN
!           Min temp for production of Hallett Mossop splinters (deg C)
     &,   HM_T_MAX
!           Max temp for production of Hallett Mossop splinters (deg C)
     &,   HM_DECAY
!           Residence distance for Hallett Mossop splinters (1/deg C)
     &,   HM_RQCL
!           Reciprocal of scaling liquid water content for HM process
!
! ----------------------------------------------------------------------
!      AUTOCONVERSION TERMS
! ----------------------------------------------------------------------
!
       REAL
     &    INHOMOG_RATE  ! Inhomogeneity factor for autoconversion rate
     &,   INHOMOG_LIM   ! Inhomogeneity factor for autoconversion limit
     &,   EC_AUTO       ! Collision collection coefficient
     &,   N_DROP_LAND    ! Droplet concentration over land
     &,   N_DROP_LAND_CR ! (N_DROP_LAND)^(-1/3)
     &,   N_DROP_SEA     ! Droplet concentration over sea
     &,   N_DROP_SEA_CR  ! (N_DROP_SEA)^(-1/3)
! WARNING. BE AWARE THAT DROPLET CONCENTRATION IS ALSO DEFINED IN THE
! RADIATION SCHEME. ARE YOU HAPPY THAT THE VALUES ARE CONSISTENT?
     &,   R_THRESH      ! Threshold droplet radius for autoconversion
!
        PARAMETER(INHOMOG_RATE=1.0
     &,           INHOMOG_LIM=1.0
     &,           EC_AUTO=0.55
! The compilation is being picky about using non integer powers in
! parameter statements. The best I can do at the moment is to
! directly define 1/cubed roots of droplet concentrations. 
! Just be careful you remember to change all of these.
     &,           N_DROP_LAND   =6.0E8
     &,           N_DROP_LAND_CR=1.18563E-3
     &,           N_DROP_SEA    =1.5E8
     &,           N_DROP_SEA_CR =1.88207E-3
!
     &,           R_THRESH=7.0E-6)
!
! The numbers 5907.24 and 4188.79 represent combinations of
! physical constants. Do NOT change them.
! PLEASE REFER TO UMDP 26 EQUATIONS P26.129 TO P26.136 FOR
! AN EXPLANATION OF THE AUTOCONVERSION PARAMETERIZATION.
! REMEMBER THAT THE PARAMETERIZATION IS VERY ROUGH.
!
        PARAMETER(AUTORATE_LAND=5907.24*EC_AUTO*INHOMOG_RATE
     &                          *N_DROP_LAND_CR
     &,           AUTORATE_SEA =5907.24*EC_AUTO*INHOMOG_RATE
     &                          *N_DROP_SEA_CR           
     &,           AUTOLIM_LAND =4188.79*R_THRESH**3
     &                          *N_DROP_LAND*INHOMOG_LIM
     &,           AUTOLIM_SEA  =4188.79*R_THRESH**3
     &                          *N_DROP_SEA *INHOMOG_LIM)
!        PARAMETER(AUTOLIM_SEA=2.155E-4
!     &,           AUTOLIM_LAND=8.621E-4
!     &,           AUTORATE_SEA=6.11
!     &,           AUTORATE_LAND=3.85)
!
! ----------------------------------------------------------------------
!     ITERATIONS OF MICROPHYSICS
! ----------------------------------------------------------------------
!
        PARAMETER(LSITER=1  
!         Advise 1 iteration for every 10 minutes
!         or less of timestep.
     &,           ADV_TYPE=2)
!         ADV_TYPE=1: Original formulation
!         ADV_TYPE=2: Revised formulation has better fall through layers
!
! ----------------------------------------------------------------------
!     NUCLEATION OF ICE
! ----------------------------------------------------------------------
!
! Note that the assimilation scheme uses temperature thresholds
! in its calculation of qsat.
        PARAMETER(M0=1.0E-12
     &,           QCFMIN=1.0E-8
     &,           TNUC=-10.0
     &,           THOMO=-40.0)
!
! ----------------------------------------------------------------------
!     HALLETT MOSSOP PROCESS
! ----------------------------------------------------------------------
!
        PARAMETER(HM_T_MIN=-8.0
! Switch off Hallett Mossop in this version but allow functionality
     &,           HM_T_MAX=-273.0
!    &,           HM_T_MAX=-3.0
     &,           HM_DECAY=1.0/7.0
     &,           HM_RQCL=1.0/0.1E-3)
!
! --------------------------COMDECK C_LSPDIF---------------------------
! input variables
!
        REAL                  !, INTENT(IN)
     &    APB1,APB2,APB3
!           Terms in deposition and sublimation
     &,   APB4,APB5,APB6
!           Terms in evap of melting snow and rain
     &,   TW1,TW2,TW3
!           Numerical fit to wet bulb temperature
     &,   TW4,TW5
!           Numerical fit to wet bulb temperature
!
! The APB and TW parameters represent diffusional growth constants
! and wet bulb temperature parameters. Do not change them.
        PARAMETER(APB1=7.14E11
     &,           APB2=1.16E8
     &,           APB3=2.416E2
     &,           APB4=5.57E11
     &,           APB5=1.03E8
     &,           APB6=2.04E2
     &,           TW1=1329.31
     &,           TW2=0.0074615
     &,           TW3=0.85E5
     &,           TW4=40.637
     &,           TW5=275.0)
!
!
      INTEGER         !, INTENT(IN)
     & POINTS
!        Number of points to be processed.
!
      REAL            !, INTENT(IN)
     & TIMESTEPFIXED
!        Timestep of physics in model (s).
     &, RHCPT(POINTS)
!        Critical humidity of all points for cloud formation.
!
      REAL            !, INTENT(IN)
     &  CFLIQ(POINTS)
!         Liquid cloud fraction in this layer.
     &, CFICE(POINTS)
!         Frozen cloud fraction in this layer.
     &, P(POINTS)
!         Air pressure at this level (Pa).
     &, RHODZ(POINTS)
!         Air mass p.u.a. in this layer (kg per sq m).
     &, SO4_ACC(POINTS)
!         Sulphur cycle variable
     &, SO4_DIS(POINTS)
!         Sulphur cycle variable
!
      REAL            !, INTENT(INOUT)
     & Q(POINTS)
!        Specific humidity at this level (kg wat per kg air).
     &,QCF(POINTS)
!        Cloud ice (kg water per kg air).
     &,QCL(POINTS)
!        Cloud liquid water (kg water per kg air).
     &,T(POINTS)
!        Temperature at this level (K).
     &,RAIN(POINTS)
!        On input: Rate of rainfall entering this layer from above.
!        On output: Rate of rainfall leaving this layer.
!                   (kg per sq m per s).
     &,SNOW(POINTS)
!        On input: Rate of snowfall entering this layer from above.
!        On Output: Rate of snowfall leaving this layer.
!                    (kg per sq m per s).
     &,VF(POINTS)
!        On input: Fall speed of ice into layer from above.
!        On output: Fall speed of ice into layer below.
!                   (m per s).
!
      LOGICAL         !, INTENT(IN)
     &  BLAND(POINTS)
!         Land/sea mask
!
!  Workspace usage: 3 real arrays---------------------------------------
      REAL
     &  QS(POINTS)
!         Saturated sp humidity for (T,p) in layer
     &, QSL(POINTS)
!         Saturated sp humidity for (T,p) in layer
!         wrt water at all temps
     &, SNOWT(POINTS)
!         Cumulative fall out of snow within iterations.
!
! external subprograms are called --------------------------------------
      EXTERNAL QSAT,QSAT_WAT
!
!  Local (derived) physical constants ----------------------------------
      REAL LCRCP,LFRCP,LSRCP,CONW,RHO1
      PARAMETER(
     &  LCRCP=LC/CP
!         Latent heat of condensation / Cp (K).
     &, LFRCP=LF/CP
!         Latent heat of fusion / Cp (K).
     &, LSRCP=LCRCP+LFRCP
!         Sum of the above (S for Sublimation).
     &, CONW=R/(EPSILON*LC)
!         Constant in wet bulb temperature calculation.
     &, RHO1=1.0
!         Reference density of air (kg/m3)
     &  )
!
! ----------------------------------------------------------------------
!  1   Define local scalars.
! ----------------------------------------------------------------------
      INTEGER
     &  I
! Loop counter (horizontal field index).
     &, J
! Counter for the iterations
!
!       Reals effectively expanded to workspace by the Cray (using
!       vector registers).
!
      REAL
!       Real workspace.  At end of DO loop, contains :-
     &  RHO(POINTS)
!         Density of air in the layer.
     &, RHOR(POINTS)
!         1.0/RHO to speed up calculations.
     &, VTEMP
!         Virtual temperature as at start of loop.
     &, TEMPC
!         temperature degree C as at start of loop.
     &, ESI(POINTS)
!         saturation vapour pressure (wrt ice below zero)
     &, ESW(POINTS)
!         saturation vapour pressure (wrt water at all T)
     &, DQI
!         increment to/from ice/snow
     &, DQIL
!         increment to/from cloud water
     &, DPR
!         increment to/from rain
     &, CFICETEMP
!         fraction of ice inferred for fall speed calculations.
     &, FQI
!         fallspeed for ice
     &, DHI(POINTS)
!         CFL limit
     &, DHIR(POINTS)
!         1.0/DHI
     &, DHILSITERR(POINTS)
!         1.0/(DHI*LSITER)
     &, FQIRQI
!         saved flux of ice out of layer
     &, FQIRQI2
!         saved flux of ice out of layer from layer above
     &, QUP
!         updated ice for long timestep
     &, QCLNEW
!         updated liquid cloud in implicit calculations
     &, TEMP7
!         term in melting
     &, PR02
!         term in evaporation of rain
     &, PR04
!         square of pr02
     &, QC
!         term in autoconversion of cloud to rain
     &, APLUSB
!         denominator in deposition or evaporation of ice
     &, CORR(POINTS)
!         density correction for fall speed
     &, ROCOR(POINTS)
!         density correction for fall speed
     &, VR1
!         Mean fall speed of rain
     &, VS1
!         Mean fall speed of snow
     &, LAMR1
!         Inverse lambda in rain exponential distribution
     &, LAMR2
!         Inverse lambda in rain exponential distribution
     &, LAMFAC1
!         Expression containing calculations with lambda
     &, LAMS1
!         Inverse lambda in snow exponential distribution
     &, FV1
!         Mean velocity difference between rain and snow
     &, TIMESTEP
!         Timestep of each iteration
     &, CORR2(POINTS)
!         Temperature correction of viscosity etc.
     &, RHNUC
!         Relative humidity required for nucleation
     &, TCG(POINTS)
!         Temperature Factor for X1I in Cox Golding calculation
     &, TCGI
!         Inverse of TCG
     &, RATEQ
!         Constant effecting rate of deposition/sublimation of ice
     &, RATEQCF
!         Constant representing effect of sub grid distribution of ice
     &, RATEQS(POINTS)
!         Critical humidity for ice deposition
     &, RATEQSL(POINTS)
!         Critical humidity for rain evaporation
     &, HM_NORMALIZE
!         Normalization for Hallett Mossop process
     &, HM_RATE
!         Increase in deposition due to Hallett Mossop process
! Obtain the size for CX and CONSTP
! Sets up the size of arrays for CX and CONSTP
      REAL CX(16),CONSTP(16)
!
!
! ----------------------------------------------------------------------
!  2.1 Start the microphysics calculations
! ----------------------------------------------------------------------
! Set up the iterations
       TIMESTEP=TIMESTEPFIXED/LSITER
! Set up sub grid scale constants. For no sub grid scale variability
! use the set up RATEQ=1.0, RATEQCF=0.0, RATEQS=1.0.
! Ideally represent these in a comdeck but require RHCRIT
       RATEQ=1.0
       RATEQCF=0.0
! Set up SNOWT(I) to be zero for all I.
! Points_do1:
       DO I=1,POINTS
         SNOWT(I)=0.0
       END DO ! Points_do1
! Set up Hallett Mossop calculation
       HM_NORMALIZE=1.0/(1.0-EXP((HM_T_MIN-HM_T_MAX)*HM_DECAY))
! ----------------------------------------------------------------------
!  2.2 Start iterating.
! ----------------------------------------------------------------------
! Iters_do1:
       DO J=1,LSITER
! ----------------------------------------------------------------------
!  2.3  Calculate sat humidity mixing ratios
! ----------------------------------------------------------------------
! Qsat with respect to ice
       CALL QSAT(QS,T,P,POINTS)
!
! Qsat with respect to liquid water
       CALL QSAT_WAT(QSL,T,P,POINTS)
! ----------------------------------------------------------------------
!  2.4 Start loop over points.
! ----------------------------------------------------------------------
! Points_do2:
       DO I=1,POINTS
!
! ----------------------------------------------------------------------
!  3.1 Calculate density of air, RHO, via virtual temperature.
! ----------------------------------------------------------------------
         VTEMP=T(I)*(1.+C_VIRTUAL*Q(I)-QCL(I)-QCF(I)) ! Virtual Temp
         RHO(I)=P(I)/(R*VTEMP)
         RHOR(I)=1.0/RHO(I)
! Correction factor of fall speeds etc. due to density.
         CORR(I)=(RHO1*RHOR(I))**0.4
! Correction factor in viscosity etc. due to temperature.
         CORR2(I)=(T(I)/273.0)**1.5 * (393.0/(T(I)+120.0))
       ENDDO
       DO I=1,POINTS
         TEMPC=T(I)-ZERODEGC
! Combined correction factor
         ROCOR(I)=SQRT(RHO(I)*CORR(I)*CORR2(I))
! Calculate a temperature factor for N0snow. CX(13)=1.0 if there is a
! temperature dependence, and 0.0 if there is not.
         TCG(I)=EXP(-CX(13)*TEMPC/8.18)
       ENDDO
       DO I=1,POINTS
! ----------------------------------------------------------------------
!  3.2 Set T in deg C and saturated vapour pressures in N/m2
! ----------------------------------------------------------------------
         ESI(I)=QS(I)*P(I)/EPSILON
         ESW(I)=QSL(I)*P(I)/EPSILON
         TCGI=1.0/TCG(I)
         TEMPC=T(I)-ZERODEGC
! ----------------------------------------------------------------------
! 3.3 Calculate RATEQS and RATEQSL as a func of RHCRIT and cloud fracs.
! ----------------------------------------------------------------------
         RATEQS(I)=RHCPT(I)+CFICE(I)*(1.0-RHCPT(I))
         RATEQSL(I)=RHCPT(I)+CFLIQ(I)*(1.0-RHCPT(I))
!
! ----------------------------------------------------------------------
!  4   Check that ice cloud fraction is sensible.
! ----------------------------------------------------------------------
         CFICETEMP=CFICE(I)
! The possibility exists in multiple iterations that ice cloud fraction
! is equal to zero but nucleation and deposition from the last iteration
! has produced a finite ice content. Hence this section produces a fix
! which will stop the scheme crashing. Only need to use for more than 1
! iteration
         IF (LSITER.GT.1) THEN
           IF (QCF(I).GT.0.0 .AND. CFICE(I).LE.0.1) THEN
             CFICETEMP=MAX(CFLIQ(I),0.1)
           END IF
         END IF
!
! ----------------------------------------------------------------------
!  5   Falling ice is advected downwards
! ----------------------------------------------------------------------
! Estimate fall speed out of this layer. We want to avoid advecting
! very small amounts of snow between layers, as this can cause numerical
! problems in other routines, so if QCF is smaller than a single
! nucleation mass per metre cubed don't advect it.
         IF (QCF(I).GT.M0) THEN
!
! Estimate the mean fall speed across the entire gridbox.
! Use a top hat distrubution within the gridbox
           FQI=CONSTP(3)*CORR(I)*
     &         (RHO(I)*QCF(I)*CONSTP(1)*TCGI/CFICETEMP)**CX(3)
         ELSE
! QCF is smaller than zero so set fall speed to zero
           FQI=0.0
! Endif for calculation of fall speed
         END IF
! Calculate CFL quantity of timestep over level separation.
         DHI(I)=TIMESTEP*RHO(I)/RHODZ(I)
! Define DHIR and DHILSITERR(I) to speed up calculations.
         DHIR(I)=1.0/DHI(I)
         DHILSITERR(I)=1.0/(DHI(I)*LSITER)
! ----------------------------------------------------------------------
! Choice of advection methods to use
! ----------------------------------------------------------------------
         IF (ADV_TYPE .EQ. 1) THEN
! ----------------------------------------------------------------------
!  5.1 Original advection scheme
! ----------------------------------------------------------------------
! See if fall speed is small enough that not all the ice falls out of
! the layer.
           IF(VF(I).LE.DHIR(I))THEN
! short timestep solution
             VF(I)=FQI
             IF(VF(I).LE.DHIR(I))THEN
! flux out is just represented by the fall speed estimated above
               FQIRQI=FQI*RHO(I)*QCF(I)
             ELSE
! cannot allow more ice to leave than was already there
               FQIRQI=RHO(I)*QCF(I)*DHIR(I)
             ENDIF
! calculate new ice content in this layer by flux divergence
             QCF(I)=QCF(I)+(SNOW(I)-FQIRQI)*DHI(I)*RHOR(I)
           ELSE
! long timestep case
             QUP = SNOW(I)*RHOR(I)/VF(I)
             FQIRQI = SNOW(I) - (RHODZ(I)*(QUP-QCF(I))/TIMESTEP)
             QCF(I) = QUP
! VF must be as least as great as the fall velocity of the current layer
             IF(VF(I).LT.FQI) VF(I)=FQI
!
! END IF for VF(I).LE.DHIR(I)
!
           END IF
! ----------------------------------------------------------------------
!  5.2 Modified advection scheme treats better ice falling across layers
! ----------------------------------------------------------------------
           ELSEIF (ADV_TYPE .EQ. 2) THEN
! Fall of ice OUT of the layer
! FQIRQI is the flux out
           FQIRQI=RHO(I)*QCF(I)*MIN(FQI,DHIR(I))
! QCF(I) is what remains in the layer
           QCF(I)=QCF(I)-FQIRQI*DHI(I)*RHOR(I)
! Fall of ice INTO the layer
! QUP is ice content from flux in which remains in layer
           IF (VF(I).GT.DHIR(I)) THEN
             QUP=SNOW(I)*RHOR(I)/VF(I)
! FQIRQI2 is flux straight through the layer
             FQIRQI2=SNOW(I)-RHO(I)*QUP*DHIR(I)
           ELSE
             QUP=SNOW(I)*RHOR(I)*DHI(I)
             FQIRQI2=0.0
           ENDIF
! QCF is updated ice content in the layer
           QCF(I)=QCF(I)+QUP
! FQIRQI is updated flux out of the layer
           FQIRQI=FQIRQI+FQIRQI2
! Now update fall speed out of layer. This is a weighted average
! of fall speed from layer and excess fall speed from fall
! through the layer.
!          VF(I)=MAX(FQI,VF(I)-DHIR(I))
           IF (FQIRQI2.GT.0.0) THEN
             VF(I)=FQI + FQIRQI2/FQIRQI * (VF(I)-DHIR(I)-FQI)
           ELSE
             VF(I)=FQI
           ENDIF
! End of advection calculations for type 2
! ----------------------------------------------------------------------
!  5.3 Other advection methods aren't written yet!
! ----------------------------------------------------------------------
         ELSE
! Error: Advection type does not exist
!
! ENDIF for advection method
         ENDIF
! ----------------------------------------------------------------------
!  5.4 Snow is used to save fall out of layer
!      for calculation of fall into next layer
! ----------------------------------------------------------------------
         SNOWT(I)=SNOWT(I)+FQIRQI/LSITER
!
! ----------------------------------------------------------------------
!      Transfer processes only active at T less than 0 deg C
! ----------------------------------------------------------------------
         IF(T(I).LT.ZERODEGC) THEN
!
! ----------------------------------------------------------------------
!  6.1 Homogenous nucleation takes place at temperatures less than THOMO
! ----------------------------------------------------------------------
            IF (T(I).LT.(ZERODEGC+THOMO)) THEN
! Turn all liquid to ice
              QCF(I)=QCF(I)+QCL(I)
              T(I)=T(I)+LFRCP*QCL(I)
              QCL(I)=0.0
            END IF
! ----------------------------------------------------------------------
!  6.2 Heteorgenous nucleation occurs for temps less than TNUC deg C
! ----------------------------------------------------------------------
           IF (T(I).LT.(ZERODEGC+TNUC)) THEN
! Calculate number of active ice nucleii
             DQI=MIN(0.01*EXP(-0.6*TEMPC),1.0E5)
! Each nucleus can grow to arbitary mass of M0 kg
             DQI=M0*DQI*RHOR(I)
! RHNUC represents how much moisture is available for ice formation.
             RHNUC=(188.92+2.81*(T(I)-ZERODEGC)
     &       +0.013336*(T(I)-ZERODEGC)**2-10.0)*0.01
             RHNUC=MIN(RHNUC,1.0)
! Predict transfer of mass to ice.
             DQI=MAX(MIN(DQI,Q(I)+QCL(I)
     &           -RATEQS(I)*MAX(QSL(I)*RHNUC,QS(I))),0.0)
             QCF(I)=QCF(I)+DQI
! This comes initially from liquid water
             DQIL=MIN(DQI,QCL(I))
             QCL(I)=QCL(I)-DQIL
             T(I)=T(I)+LFRCP*DQIL
! If more moisture is required then nucleation removes from vapour.
             DQI=DQI-DQIL
             T(I)=T(I)+LSRCP*DQI
             Q(I)=Q(I)-DQI
! END IF for nucleation
           END IF
!
! ----------------------------------------------------------------------
!  7   Deposition/Sublimation of snow - explicit.
!      Hallett Mossop process enhances growth.
! ----------------------------------------------------------------------
           IF(QCF(I).GT.M0) THEN
! Calculate transfer rate as a function of QCF and T
             PR02=RHO(I)*QCF(I)*CONSTP(1)*TCGI
             APLUSB=(APB1-APB2*T(I))*ESI(I)
             APLUSB=APLUSB+(T(I)**3)*P(I)*APB3
             DQI=TCG(I)*CONSTP(5)*T(I)**2*ESI(I)*RATEQ*
     &      (MIN((Q(I)+QCL(I)),QSL(I))-RATEQCF*QCF(I)-RATEQS(I)*QS(I))*
     &       (0.65*CONSTP(13)*CORR2(I)*PR02**CX(1)+CONSTP(6)*ROCOR(I)*
     &       PR02**CX(2))/(QS(I)*APLUSB*RHO(I))
! Limits depend on whether deposition or sublimation occurs
             IF (DQI.GT.0.0) THEN
! Deposition is occuring.
! Hallett Mossop Enhancement
               IF ( (T(I)-ZERODEGC) .GE. HM_T_MAX) THEN
! Temperature is greater than maximum threshold for HM.
                 HM_RATE=0.0
               ELSEIF ((T(I)-ZERODEGC) .LT. HM_T_MAX
! Temperature is between HM thresholds
     &           .AND. (T(I)-ZERODEGC) .GT. HM_T_MIN) THEN
                 HM_RATE=(1.0-EXP( (T(I)-ZERODEGC-HM_T_MAX)*HM_DECAY) )
     &             *HM_NORMALIZE
               ELSE
! Temperature is less than minimum threshold for HM.
                 HM_RATE=EXP( (T(I)-ZERODEGC-HM_T_MIN)*HM_DECAY)
               ENDIF
! Calculate enhancement factor for HM process.
               HM_RATE=1.0+HM_RATE*QCL(I)*HM_RQCL
! Calculate Transfer. Limit is available moisture.
               DQI=MIN(DQI*TIMESTEP*HM_RATE,
     &         (Q(I)+QCL(I)-RATEQCF*QCF(I)-RATEQS(I)*QS(I))
     &          /(1.0+RATEQCF))
             ELSE
! Sublimation is occuring. Limits are spare moisture capacity and QCF
               DQI=MAX(MAX(DQI*TIMESTEP,
     &           (Q(I)+QCL(I)-RATEQCF*QCF(I)-RATEQS(I)*QS(I))
     &                            /(1.0+RATEQCF)),-QCF(I))
             END IF
! Adjust ice content
             QCF(I)=QCF(I)+DQI
             DQIL=MAX(MIN(DQI,QCL(I)),0.0)
! Adjust liquid content (deposits before vapour by Bergeron Findeison
!  process).
             QCL(I)=QCL(I)-DQIL
             T(I)=T(I)+LFRCP*DQIL
             DQI=DQI-DQIL
! Adjust vapour content
             Q(I)=Q(I)-DQI
             T(I)=T(I)+LSRCP*DQI
! END IF for QCF.GT.M0.
           END IF
!
! ----------------------------------------------------------------------
!  8   Riming of snow by cloud water -implicit in QCL
! ----------------------------------------------------------------------
           IF (QCF(I).GT.M0.AND.QCL(I).GT.0.0) THEN
               QCLNEW=QCL(I)/(1.0+CONSTP(4)*TCG(I)*CORR(I)*TIMESTEP*
     &         (RHO(I)*QCF(I)*CONSTP(1)*TCGI)**CX(4))
! Recalculate water contents
               QCF(I)=QCF(I)+(QCL(I)-QCLNEW)
               T(I)=T(I)+LFRCP*(QCL(I)-QCLNEW)
               QCL(I)=QCLNEW
! END IF for QCF.GT.M0.AND.QCL(I).GT.0.0
           END IF
!
! ----------------------------------------------------------------------
!  9   Capture of rain by snow - implicit in rain
! ----------------------------------------------------------------------
           IF (RAIN(I).GT.0.0.AND.QCF(I).GT.M0) THEN
! Calculate velocities
             VR1=CORR(I)*CONSTP(11)/6.0*
     &              (RAIN(I)/(CONSTP(8)*CORR(I)))**CX(5)
             VS1=CONSTP(3)*CORR(I)*(RHO(I)*QCF(I)*CONSTP(1)*TCGI)**CX(3)
! Estimate the mean absolute differences in velocities.
             FV1=MAX(ABS(VR1-VS1),(VR1+VS1)/8.0)
! Calculate functions of slope parameter lambda
             LAMR1=(RAIN(I)/(CONSTP(8)*CORR(I)))**(CX(10))
             LAMS1=(RHO(I)*QCF(I)*CONSTP(1)*TCGI)**(-CX(6))
             LAMFAC1=CONSTP(16)*(LAMR1**6.0*LAMS1**CX(16)) +
     &               CONSTP(15)*(LAMR1**5.0*LAMS1**CX(15)) +
     &               CONSTP(14)*(LAMR1**4.0*LAMS1**CX(14))
! Calculate transfer
             DPR=TCG(I)*CONSTP(9)*LAMS1**(-CX(8))*LAMR1**(-CX(9))*FV1*
     &       LAMFAC1*TIMESTEP*RHOR(I)
             DPR=MIN(DPR,RAIN(I)*(DHI(I)*LSITER)*RHOR(I))
! Adjust ice and rain contents
             QCF(I)=QCF(I)+DPR
             RAIN(I)=RAIN(I)-DPR*RHO(I)*DHILSITERR(I)
             T(I)=T(I)+LFRCP*DPR
!      Endif for RAIN.GT.0.0 in capture term
           END IF
! ----------------------------------------------------------------------
!      End of transfer processes only active at T less than 0 deg C
! ----------------------------------------------------------------------
         END IF
! ----------------------------------------------------------------------
!  10  Evaporate melting snow - implicit in subsaturation
! ----------------------------------------------------------------------
         IF(QCF(I).GT.M0.AND.T(I).GT.ZERODEGC)THEN
! Calculate transfer as a function of QCF, T and specific humidity
           PR02=RHO(I)*QCF(I)*CONSTP(1)*TCGI
           PR04=((APB4-APB5*T(I))*ESW(I)+APB6*P(I)*T(I)**3)
           DPR=TCG(I)*RATEQ*CONSTP(5)*T(I)**2*ESW(I)*TIMESTEP*
     &     (0.65*CONSTP(13)*CORR2(I)*PR02**CX(1)
     &      +CONSTP(6)*ROCOR(I)*PR02**CX(2))/(QSL(I)*RHO(I)*PR04)
           DPR=DPR*MAX(
     &         (RATEQS(I)*QSL(I)+RATEQCF*QCF(I)-Q(I)-QCL(I)),0.0)
     &                /(1.0+DPR*(1.0+RATEQCF))
! Extra check to see we don't get a negative QCF
           DPR=MIN(DPR,QCF(I))
! Update values of ice and vapour
           QCF(I)=QCF(I)-DPR
           Q(I)=Q(I)+DPR
           T(I)=T(I)-DPR*LSRCP
         END IF
!
! ----------------------------------------------------------------------
!  11  Melting of snow - explicit
!      USE WET BULB TEMP (DEG.C) IN SNOW MELT CALC.
!      Use a numerical approximation.
! ----------------------------------------------------------------------
         IF(QCF(I).GT.M0.AND.T(I).GT.ZERODEGC)THEN
           TEMPC=T(I)-ZERODEGC
! An approximate calculation of wet bulb temperature
           TEMP7=TEMPC-RATEQ*
     &            (RATEQS(I)*QSL(I)+RATEQCF*QCF(I)-Q(I)-QCL(I))
     &           *(TW1+TW2*(P(I)-TW3) - TW4*(T(I)-TW5) )
           TEMP7=MAX(TEMP7,0.0)
! End of wet bulb temp formulations.
           PR02=RHO(I)*QCF(I)*CONSTP(1)*TCGI
           DPR=TCG(I)*CONSTP(7)*TIMESTEP*
     &            (0.65*CONSTP(13)*CORR2(I)*PR02**CX(1)
     &         + CONSTP(6)*ROCOR(I)*PR02**CX(2))*RHOR(I)
! Solve implicitly in terms of temperature
           DPR=TEMP7*(1.0-1.0/(1.0+DPR*LFRCP))/LFRCP
           DPR=MIN(DPR,QCF(I))
! Update values of ice and Rain
           QCF(I)=QCF(I)-DPR
           RAIN(I)=RAIN(I)+DPR*RHO(I)*DHILSITERR(I)
           T(I)=T(I)-LFRCP*DPR
! ENDIF for melting snow
         END IF
       ENDDO
!
! ----------------------------------------------------------------------
!  12  Evaporation of rain - implicit in subsaturation
! ----------------------------------------------------------------------
       DO I=1,POINTS
         IF(RAIN(I).GT.0.0)THEN
           PR04=((APB4-APB5*T(I))*ESW(I)+APB6*P(I)*T(I)**3)
! Define LAMR1 and LAMR2
           LAMR1=RAIN(I)/(CONSTP(8)*CORR(I))
           LAMR2=LAMR1**(CX(12)*CX(10))
           LAMR1=LAMR1**(CX(11)*CX(10))
! New, consistent evaporation method, with rain fall speed relationship.
           DPR=CONSTP(2)*T(I)**2*ESW(I)*TIMESTEP
           DPR=DPR*( (0.78*CORR2(I)*LAMR2)
     &               + (CONSTP(12)*ROCOR(I)*LAMR1) )
           DPR=DPR*RATEQ/(QSL(I)*RHO(I)*PR04)
! Calculate transfers.
           DPR=DPR*MAX((RATEQSL(I)*QSL(I)-Q(I)-QCL(I)),0.0)/(1.0+DPR)
           DPR=DPR*RHO(I)*DHILSITERR(I)
           DPR=MIN(DPR,RAIN(I))
! Update values of rain and vapour
           RAIN(I)=RAIN(I)-DPR
           Q(I)=Q(I)+DPR*DHI(I)*LSITER*RHOR(I)
           T(I)=T(I)-DPR*LCRCP*DHI(I)*LSITER*RHOR(I)
! END IF for evaporation of rain.
         END IF
!
! ----------------------------------------------------------------------
!  13  Accretion of cloud on rain - implicit in liquid water content
! ----------------------------------------------------------------------
         IF(RAIN(I).GT.0.0.AND.QCL(I).GT.0.0)THEN
! New accretion formulation.
           PR02=RAIN(I)/(CONSTP(8)*CORR(I))
           QCLNEW=QCL(I)/
     &            ((1.0+CONSTP(10)*CORR(I)*TIMESTEP*PR02**CX(7)))
! Now calculate increments to rain.
           RAIN(I)=RAIN(I)+(QCL(I)-QCLNEW)*RHO(I)*DHILSITERR(I)
           QCL(I)=QCL(I)-(QCL(I)-QCLNEW)
! END IF for accretion of cloud on rain.
         END IF
       ENDDO
!
! ----------------------------------------------------------------------
!  14  Autoconversion of cloud to rain - explicit
! ----------------------------------------------------------------------
       DO I=1,POINTS
         IF (QCL(I).GT.0.0.AND.CFLIQ(I).GT.0.0) THEN
! Use a liquid cloud fraction here as this term is very non-linear
! The section below is a simple way of proceeding.
           IF (BLAND(I)) THEN
! Land point
             QC=MIN(AUTOLIM_LAND*CFLIQ(I)*RHOR(I),QCL(I))
           ELSE
! Sea point
             QC=MIN(AUTOLIM_SEA*CFLIQ(I)*RHOR(I),QCL(I))
           END IF
           IF (BLAND(I)) THEN
! Land point
             DPR=MIN(AUTORATE_LAND*(RHO(I)*QCL(I)/CFLIQ(I))**1.333
     &                 *TIMESTEP*QCL(I)/CORR2(I),QCL(I)-QC)
           ELSE
! Sea point
             DPR=MIN(AUTORATE_SEA*(RHO(I)*QCL(I)/CFLIQ(I))**1.333
     &                 *TIMESTEP*QCL(I)/CORR2(I),QCL(I)-QC)
           END IF
! End of calculation of autoconversion amount DPR
           QCL(I)=QCL(I)-DPR
           RAIN(I)=RAIN(I)+DPR*RHO(I)*DHILSITERR(I)
! ENDIF for autoconversion.
         END IF
! ----------------------------------------------------------------------
!  15  Now continue the loops over points and iterations.
! ----------------------------------------------------------------------
! Continue DO loop over points
       END DO ! Points_do2
! Continue DO loop over iterations
       END DO ! Iters_do1
!
! Copy contents of SNOWT to SNOW, to fall into next layer down
! Points_do3
      DO I=1,POINTS
        SNOW(I)=SNOWT(I)
! ----------------------------------------------------------------------
!  16 Melt any SNOW which has reached here, as long as T is large enough
! ----------------------------------------------------------------------
! Only use if long timestep case. In which case melt the excess snow
! which falls straight through a layer. Use DQI variable to save space.
! DQI APPROXIMATELY represents the excess SNOW.
           DQI=MIN(VF(I)*DHI(I)-1.0,1.0)
           IF (DQI.GT.0.0) THEN
! Long timestep case
             TEMPC=T(I)-ZERODEGC
             IF (SNOW(I).GT.0.0.AND.T(I).GT.ZERODEGC) THEN
! Numerical approximation of wet bulb temperature.
               TEMP7=TEMPC-RATEQ*(RATEQS(I)*QSL(I)
     &           +RATEQCF*QCF(I)-Q(I)-QCL(I))*(TW1+TW2*
     &           (P(I)-TW3) - TW4*(T(I)-TW5) )
               TEMP7=MAX(TEMP7,0.0)
! End of wet bulb calculation
               DPR=TEMP7/(LFRCP*LSITER)
               DPR=MIN(DPR,SNOW(I)*DHI(I)*RHOR(I)*DQI)
! Update values of snow and rain
               SNOW(I)=SNOW(I)-DPR*RHO(I)*DHIR(I)
               RAIN(I)=RAIN(I)+DPR*RHO(I)*DHIR(I)
               T(I)=T(I)-LFRCP*DPR*LSITER
! END IF for long timestep
             END IF
! END IF for melting of excess snow.
           END IF
! ----------------------------------------------------------------------
!  17  Remove any small amount of QCF which is left over to be tidy.
!      If QCF is less than QCFMIN and isn't growing
!      by deposition (assumed to be given by RHCPT) then remove it.
! ----------------------------------------------------------------------
!           DQI=M0*RHOR(I)*MIN( 0.01*EXP(-0.6*TEMPC),1.0E5 )
!           IF (QCF(I).LT.MIN( MAX(M0*RHOR(I),DQI),1.0E-5*QS(I) ) ) THEN
           IF (QCF(I).LT.QCFMIN.AND.
     &     (T(I).GT.ZERODEGC .OR. (Q(I)+QCL(I) .LE. RHCPT(I)*QS(I))
     &     .OR. QCF(I).LT.0.0)  )  THEN
             Q(I)=Q(I)+QCF(I)
             T(I)=T(I)-LSRCP*QCF(I)
             QCF(I)=0.0
           END IF
! END DO for melting of excess snow loop over points.
      END DO ! Points_do3
! ----------------------------------------------------------------------
!  18  End of the LSP_ICE subroutine
! ----------------------------------------------------------------------
      RETURN
      END
