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
CLL
CLL   SUBROUTINE SLABCNTL
CLL   -------------------
CLL
CLL   THIS ROUTINE IS FOR USE WITH THE 'SLAB' OCEAN MODEL ONLY.
CLL
CLL   THE ROUTINE CALLS TRANSA2S TO CONVERT THE ATMOSPHERE
CLL   VARIABLES TO THE FORMAT REQUIRED FOR THE SLAB MODEL,
CLL   SLAB TO UPDATE THE SEA SURFACE TEMPERATURE, SLABICE TO
CLL   UPDATE THE ICE CONCENTRATION AND DEPTH, AND TRANSS2A TO
CLL   CONVERT THE SLAB VARIABLES BACK TO THE ATMOSPHERE FORMAT.
CLL
CLL   THIS ROUTINE FORMS PART OF SYSTEM COMPONENT P40.
CLL   IT CAN BE COMPILED BY CFT77, BUT DOES NOT CONFORM TO ANSI
CLL   FORTRAN77 STANDARDS, BECAUSE OF THE INLINE COMMENTS.
CLL   DYNAMIC ALLOCATION AND THE USE OF ENDDO.
CLL   IT ADHERES TO THE STANDARDS OF DOCUMENTATION PAPER 3, VERSION 5.
CLL
CLL   ALL QUANTITIES IN THIS ROUTINE ARE IN S.I. UNITS UNLESS
CLL   OTHERWISE STATED.
CLL
CLL   CALLED BY: SLABSTEP
CLL   WRITTEN BY C.A.SENIOR (15/1/91)
CLL   MODIFIED BY A.B.KEEN (02/02/93)
CLL   MODIFIED BY C.A.SENIOR (22/03/93)
CLL   MODIFIED BY C.A.SENIOR (30/03/93)
CLL   MODIFIED BY A.B.KEEN (22/04/93)
CLL   MODIFIED BY A.B.KEEN (27/04/93)
CLL   MODIFIED BY C.A.SENIOR (08/07/93)
CLL   MODIFIED BY C.A.SENIOR (24/02/94)
CLL   Modified at version 3.4 by J.Thomson, C.Senior, R.E.Carnell
CLL   Changes to satisfy slab code review.
CLL   Addition of ice dynamics and alternative version of thermodyn.
CLL   Slabtemp advection added.
CLL   Version    Description of change
CLL     4.0      Vertical SST advection added (R.Carnell)
CLL              Corrections made to ice model and new diagnostics
CLL              added. J.F.Crossley
!LL  4.4   04/08/97  Add missing ARGOINDX to various argument lists.
!LL                  D. Robinson.
!LL
CLL   VERSION NUMBER 1.1
CLL   REVIEWER: W.INGRAM (02/03/93)
CLL   DOCUMENTATION: UM DOCUMENTATION PAPER 58; THE SLAB OCEAN MODEL
CLL
CLLEND---------------------------------------------------------------
C*L
      SUBROUTINE SLABCNTL(
!========================== COMDECK ARGOINDX ==========================
     &  J_1, J_2, J_3, J_JMT, J_JMTM1, J_JMTM2, J_JMTP1, JST, JFIN 
     &, J_FROM_LOC, J_TO_LOC, JMT_GLOBAL, JMTM1_GLOBAL, JMTM2_GLOBAL  
     &, JMTP1_GLOBAL, J_OFFSET, O_MYPE, O_EW_HALO, O_NS_HALO        
     &, J_PE_JSTM1, J_PE_JSTM2, J_PE_JFINP1, J_PE_JFINP2, O_NPROC,
     &  imout,jmout,J_PE_IND_MED,NMEDLEV,lev_med,lev_hud,imout_hud,
     &  jmout_hud,J_PE_IND_HUD,med_topflow,
     & L1,U_FIELD,ICOLS,JROWS,LAND,DT,DZ1,
     + SOLARIN,BLUEIN,EVAP,LONGWAVE,SENSIBLE,HEATCONV,
     + SNOWLS,SNOWCONV,TSTARATM,SLABTEMP,HICEATM,HSNOWATM,
     + AICEATM,SUBLIMA,TOPMELTZ,BOTMELTZ,
     + UICE,VICE,
     + UCURRENT,VCURRENT,WSX,WSY,
     + H0,AMXSOUTH,AMXNORTH,AICEMIN,HICEMIN,
     + TCLIM,HCLIM,CALIB,HICESLB,
     + AINC_DYN,HINC_DYN,HSINC_DYN,HINC_DIFF,
     + HINC_ADV,HSINC_ADV,AREAS,
     + AINC_THERM,HINC_THERM,HSINC_THERM,OIFLUX,
     + PRESSURE,PMAX,LEADFLUX,ATMSFLUX,DTADV,DTDIFF,CARYHEAT,
     + SNOWSLAB,SNOWLEAD,DTICE,
     + EDDYDIFF,epsilon,Ah,HCLIMIT,Ah_ice,
     + Pstar_ice_strength,kappa_ice_strength,cdw,tol_icav,tol_ifree,
     + weight_ifree,nmax_icav,nmax_ifree,
     + L_THERM,L_IDYN,L_IDRIF,LGLOBAL,L_SLBADV,
     + COS_P_LATITUDE,COS_U_LATITUDE,SEC_P_LATITUDE,
     + SIN_U_LATITUDE,CORIOLIS,
     + ADJHCONV,wtsfc,wtbase,
     + DELTA_LONG,DELTA_LAT,BASE_LAT)
C
!========================== COMDECK TYPOINDX ===========================
!
! Description:
!
!       This comdeck contains all the indices and row-wise loop
!       control variables required by the ocean MPP code.
!
! History:
!
!=======================================================================

      INTEGER J_1     ! Local value of loop control for J = 1, n
     &,       J_2     !   "     "    "   "     "     "  J = 2, n
     &,       J_3     !   "     "    "   "     "     "  J = 3, n
     &,       J_JMT   !   "     "    "   "     "     "  J = n, JMT
     &,       J_JMTM1 !   "     "    "   "     "     "  J = n, JMTM1
     &,       J_JMTM2 !   "     "    "   "     "     "  J = n, JMTM2
     &,       J_JMTP1 !   "     "    "   "     "     "  J = n, JMTP1
     &,       JST     ! First row this process considers (no halo)
     &,       JFIN    ! Last   "    "     "     "        "
     &,       J_FROM_LOC
     &,       J_TO_LOC
     &,       JMT_GLOBAL       ! Global value of JMT
     &,       JMTM1_GLOBAL     ! Global value of JMT - 1
     &,       JMTM2_GLOBAL     ! Global value of JMT - 2
     &,       JMTP1_GLOBAL     ! Global value of JMT + 1
     &,       J_OFFSET         ! Global value of JST - 1
     &,       O_MYPE           ! MYPE for ocean arg lists
     &,       O_EW_HALO        ! EW HALO for ocean arg lists
     &,       O_NS_HALO        ! NS HALO for ocean arg lists
     &,       J_PE_JSTM1
     &,       J_PE_JSTM2
     &,       J_PE_JFINP1
     &,       J_PE_JFINP2
     &,       O_NPROC
     &,       imout(4),jmout(4)! i,j indices for pts in Med outflow
     &,       J_PE_IND_MED(4)  ! no for each PE in Med outflow
     &,       NMEDLEV          ! no of levels for Med outflow
     &,      lev_med   ! level at which deep advective Med outflow
c                        exits the Mediterranean
     &,      lev_hud   ! level at which deep advective flow
c                        enters the Hudson Bay
     &,      imout_hud(4)  ! zonal index for Hudson Bay outflow
     &,      jmout_hud(4)  ! merid index for Hudson Bay outflow
     &,      J_PE_IND_HUD(4)  ! PE's involved in HB outflow
     &,      med_topflow   ! last level for which there is inflow to 
C                          ! Mediterranean




C
      INTEGER L1,           ! IN SIZE OF DATA VECTORS (P GRID)
     + U_FIELD,             ! IN SIZE OF U FIELDS
     + ICOLS,               ! IN NUMBER OF POINTS IN A ROW
     + JROWS                ! IN NUMBER OF POINTS IN A COLUMN
C
      LOGICAL LAND(L1)      ! IN ATMOSPHERIC MODEL LAND-SEA MASK
     +                      !    FALSE AT OCEAN POINTS
     +,CALIB                ! IN TRUE IF SLAB CALIBRATION EXPT
     +,LGLOBAL              ! IN TRUE if global model (hence cyclic)
     +,L_THERM              ! IN TRUE if coupled model ice thermo
     +,L_IDYN               ! IN TRUE if cav fluid ice dynamics
     +,L_IDRIF              ! IN TRUE if ice depth advection
     +,L_SLBADV             ! IN TRUE if slabtemp advection
C
      REAL DT            ! IN TIMESTEP FOR UPDATING THE SLAB OCEAN IN S
     +,DZ1               ! IN THICKNESS OF THE SLAB OCEAN IN METRES.
     +,DELTA_LONG        ! IN EW grid spacing (degrees)
     +,DELTA_LAT         ! IN NS grid spacing (degrees)
     +,BASE_LAT          ! IN latitude of first row (degrees)
     +,AH                ! IN Diffusion Coefficent for slab temperature
      REAL AH_ICE        ! IN Diffusion coeff. for ice depth.
     &,epsilon           ! IN Minimum depth of gbm ice (metres)
     &,cdw               ! IN quadratic water stress coefficient
     &,Pstar_ice_strength! IN Parameter in ice strength calculation
     &,kappa_ice_strength! IN Parameter in ice srength calculation
     &,Weight_ifree      ! IN Weighting for free drift relaxation.
C
      REAL
     + SOLARIN(L1)   ! IN NET DOWNWARD SHORTWAVE FLUX FROM THE
     +               !    ATMOSPHERE (ALL FREQUENCIES).
     +,BLUEIN(L1)    ! IN NET DOWNWARD SHORTWAVE FLUX FROM THE
     +               !    ATMOSPHERE (BAND 1, SEA POINTS)
     +,EVAP(L1)      ! IN SURFACE EVAPORATION FROM THE WATER
     +               !    FRACTION OF ALL OCEAN POINTS. AT SEA-ICE
     +               !    POINTS, THIS IS WEIGHTED BY THE
     +               !    FRACTIONAL LEAD AREA. (KG M-2 S-1)
     +,LONGWAVE(L1)  ! IN NET DOWNWARD LONGWAVE HEAT FLUX.
     +,SENSIBLE(L1)  ! IN SENSIBLE HEAT FLUX (+VE UPWARD) FOR
     +               !    THE WATER FRACTION OF ALL OCEAN POINTS.
     +               !    AREA-WEIGHTED AT SEA-ICE POINTS.
     +,HEATCONV(L1)  ! IN HEAT CONVERGENCE RATE, IN W M-2.
     +,SNOWLS(L1)    ! IN LARGE-SCALE SNOWFALL RATE (KG M-2 S-1)
     +,SNOWCONV(L1)  ! IN CONVECTIVE SNOWFALL RATE (KG M-2 S-1)
     +,TSTARATM(L1)  ! INOUT SEA SURFACE TEMP FROM ATMOS MODEL (K)
     +,SLABTEMP(L1)  ! INOUT TEMPERATURE OF THE SLAB OCEAN (C)
     +,HICEATM(L1)   ! INOUT EQUIVALENT ICE DEPTH FROM ATMOS MODEL (M)
     +,HSNOWATM(L1)  ! INOUT SNOW DEPTH FROM ATMOS MODEL(KG M-2)
     +,AICEATM(L1)   ! INOUT ICE CONCENTRATION FROM ATMOS MODEL
     +,UICE(U_FIELD) ! INOUT X COMPONENT OF ICE VELOCITY (m/s)
     +,VICE(U_FIELD) ! INOUT Y COMPONENT OF ICE VELOCITY (m/s)
     +,SUBLIMA(L1)   ! IN ACCUMULATED SUBLIMATION, IN KG M-2
     +,TOPMELTZ(L1)  ! IN RATE OF MELTING OF SNOW IN W M-2.
     +               !    (THIS CAN BE TRANSFERRED TO ICE.)
     +,BOTMELTZ(L1)  ! IN DIFFUSIVE HEAT FLUX THROUGH ICE. IN W M-2
     +               !    IF THIS IS +VE, ICE MELTS AT THE BASE.
     +               !    IF IT IS -VE, ICE ACCRETES THERE.
     +,HICESLB(L1)   ! OUT MEAN ICE DEPTH OVER WHOLE GRID BOX
     +               !     IN M.
     +,TCLIM(L1)     ! IN CLIMATOLOGICAL SEA SURFACE TEMPS K
     +,HCLIM(L1)     ! IN CLIMATOLOGICAL SEA-ICE DEPTHS M
     +,ADJHCONV(L1)  ! OUT REDISTRIBUTED HEAT CONVERGENCES
     +,wtsfc(l1)     ! OUT w x slab temp at surface
     +,wtbase(l1)    ! OUT w x slab temp at base
     +,COS_P_LATITUDE(L1)     ! IN COS LATITUDE ON P GRID
     +,COS_U_LATITUDE(U_FIELD)! IN COS LATITUDE ON UV GRID
     +,SEC_P_LATITUDE(L1)     ! IN 1.0/COS LATITUDE ON P GRID
     +,SIN_U_LATITUDE(U_FIELD)! IN SIN LATITUDE ON UV GRID
     +,CORIOLIS(L1)           ! IN 2 * OMEGA * SIN(LAT) ON P GRID

C
       REAL
     + UCURRENT(U_FIELD)    ! IN X COMPONENT OF SURFACE CURRENT (M/S)
     +,VCURRENT(U_FIELD)    ! IN Y COMPONENT OF SURFACE CURRENT (M/S)
     +,WSX(U_FIELD)         ! IN X COMPONENT OF SURFACE STRESS (N/M2)
     +,WSY(U_FIELD)         ! IN Y COMPONENT OF SURFACE STRESS (N/M2)
C
       REAL
     + AINC_THERM(L1)       ! OUT ice fraction inc (thermodynamics)
     +,HINC_THERM(L1)       ! OUT ice depth inc (thermodynamics)
     +,HSINC_THERM(L1)      ! OUT snow depth inc *ice fract (therm)
     +,AINC_DYN(L1)         ! OUT ice fraction inc (dynamics)
     +,HINC_DYN(L1)         ! OUT ice depth inc (dynamics)
     +,HSINC_DYN(L1)        ! OUT snow depth inc *ice fract (dynamics)
     +,HINC_DIFF(L1)        ! OUT ice depth inc (diffusion)
     +,HINC_ADV(L1)         ! OUT ice depth inc (advection)
     +,HSINC_ADV(L1)        ! OUT snow depth inc *ice fract (advection)
     +,AREAS(L1)            ! OUT grid box areas
     +,OIFLUX(L1)           ! OUT ocean to ice heat flux
     +,PRESSURE(L1)         ! OUT internal ice pressure
     +,PMAX(L1)             ! OUT ice strength
     +,ATMSFLUX(L1)         ! OUT net heat into slab through leads
     +,LEADFLUX(L1)         ! OUT net heat into ice through leads
     +,DTADV(L1)            ! OUT slab heating rate due to advection
     +,DTDIFF(L1)           ! OUT slab heating rate due to diffusion
     +,CARYHEAT(L1)         ! OUT negative heat flux due to slab
     +                      ! temperatures falling below freezing.
     +                      ! W M-2
     +,SNOWSLAB(L1)         ! OUT snowfall rate melting in slab.
     +,SNOWLEAD(L1)         ! OUT snowfall rate melting in leads.
     +,DTICE(L1)            ! OUT slab heating rate due to ice melt etc.
C
       REAL
     + H0            ! IN MINIMUM LOCAL ICE DEPTH (M)
     +,AMXSOUTH      ! IN MAX ICE CONC - SOUTHERN HEMISPHERE
     +,AMXNORTH      ! IN MAX ICE CONC - NORTHEN HEMISPHERE
     +,AICEMIN       ! IN MIN ICE CONC - GLOBAL
     +,HICEMIN       ! IN MIN GRIDBOX MEAN ICE DEPTH (M)
     +,HCLIMIT       ! IN LIMIT FOR REDISTRIBUTING HEAT CONVERGENCES
C
C
C     Include COMDECKS
C
C*L------------------COMDECK C_LHEAT------------------------------------
C LC IS LATENT HEAT OF CONDENSATION OF WATER AT 0DEGC
C LF IS LATENT HEAT OF FUSION AT 0DEGC
      REAL LC,LF

      PARAMETER(LC=2.501E6,
     &          LF=0.334E6)
C*----------------------------------------------------------------------
! History:
! Version  Date  Comment
!  3.4   18/5/94 Add PP missing data indicator. J F Thomson
C*L------------------COMDECK C_MDI-------------------------------------
      REAL    RMDI_PP   ! PP missing data indicator (-1.0E+30)
      REAL    RMDI_OLD  ! Old real missing data indicator (-32768.0)
      REAL    RMDI      ! New real missing data indicator (-2**30)
      INTEGER IMDI      ! Integer missing data indicator

      PARAMETER(
     & RMDI_PP  =-1.0E+30,
     & RMDI_OLD =-32768.0,
     & RMDI =-32768.0*32768.0,
     & IMDI =-32768)
C*----------------------------------------------------------------------
C
C     VARIABLES LOCAL TO THIS ROUTINE
C
      INTEGER J           ! LOOP COUNTER
     +,L2                 ! NUMBER OF DATA POINTS TO BE PROCESSED
     +,JROWSM1            ! NUMBER OF P ROWS MINUS ONE
     +,HALFROWS           ! HALF THE NUMBER OF ROWS
     +,NPOINTS            ! NUMBER OF POINTS IN EACH HEMISPHERE
     +,SPOINTS            ! TOTAL NO OF POINTS - NPOINTS
     +,SPTS1              ! FIRST POINT IN SOUTHERN HEMISPHERE
     +,NPTSP1             ! NPOINTS PLUS 1
C
      REAL
     + SNOWRATE(L1)       ! RATE OF SNOWFALL, IN KG M-2 S-1.
     +,SUBLIMZ(L1)        ! SUBLIMATION RATE (KG M_2 S_1)
     +,HSNOWSLB(L1)       ! SNOW DEPTH, NOT AVERAGED OVER GRID
     +                    !    BOX, JUST OVER THE ICE PORTION IN M.
     +,TCLIMC(L1)         ! CLIMATOLOGICAL SEA SURFACE TEMPS C
     +,AICESLB(L1)        ! ICE CONCENTRATION.
      REAL ONEEM6         ! SMALL VALUE FOR ICE FRACTION INEQUALITY.
      PARAMETER (ONEEM6 = 1.0E-06)
C
      LOGICAL ICY(L1)     ! TRUE IF BOX CONTAINS ICE.
     +                    ! THIS IS RESET IN SLABICE
     +,OPENSEA(L1)        ! TRUE IF BOX IS ICE_FREE OCEAN POINT.
     +,NEWICE(L1)         ! TRUE IF ICE IS FORMING.
C
C
C    1. INITIALISE VARIABLES
C
C       INITIALISE THE AMOUNT OF DATA TO PROCESS
C
      L2 = L1
      JROWSM1 = JROWS-1
C
C       CONVERT SUBLIMATION TO A RATE (USE LOCAL ARRAY FOR SAFETY)
C       AND INITIALISE ADJUSTED HEAT CONVERGENCE ARRAY
C
      DO J=1,L2
          SUBLIMZ(J)  = SUBLIMA(J)/DT
          ADJHCONV(J) = HEATCONV(J)
      END DO
C
C       INITIALISE THE ARRAY ICY AND THE ARRAY OPENSEA
C       AND THE ARRAYS NEWICE AND OIFLUX
C
      DO J=1,L2
        IF ( .NOT. LAND(J)) THEN
          IF ( (AICEATM(J) .GE. (AICEMIN-oneem6))) THEN
            ICY(J)     = .TRUE.
            OPENSEA(J) = .FALSE.
          ELSE
            ICY(J)     = .FALSE.
            OPENSEA(J) = .TRUE.
          ENDIF
          newice(j)     = .false.
          oiflux(j)     = 0.0
          ainc_therm(j) = 0.0
          hinc_therm(j) = 0.0
          hsinc_therm(j)= 0.0
          ainc_dyn(j)   = 0.0
          hinc_dyn(j)   = 0.0
          hsinc_dyn(j)  = 0.0
          hinc_diff(j)  = 0.0
          hinc_adv(j)   = 0.0
          hsinc_adv(j)  = 0.0
          areas(j)      = 0.0
          leadflux(j)   = 0.0
          atmsflux(j)   = 0.0
          dtadv(j)      = 0.0
          dtdiff(j)     = 0.0
          caryheat(j)   = 0.0
          snowslab(j)   = 0.0
          snowlead(j)   = 0.0
          dtice(j)      = 0.0
        ELSE
          ICY(J)     = .FALSE.
          OPENSEA(J) = .FALSE.
          newice(j)  = .false.
          oiflux(j)     = rmdi
          ainc_therm(j) = rmdi
          hinc_therm(j) = rmdi
          hsinc_therm(j)= rmdi
          ainc_dyn(j)   = rmdi
          hinc_dyn(j)   = rmdi
          hsinc_dyn(j)  = rmdi
          hinc_diff(j)  = rmdi
          hinc_adv(j)   = rmdi
          hsinc_adv(j)  = rmdi
          areas(j)      = rmdi
          leadflux(j)   = rmdi
          atmsflux(j)   = rmdi
          dtadv(j)      = rmdi
          dtdiff(j)     = rmdi
          caryheat(j)   = rmdi
          snowslab(j)   = rmdi
          snowlead(j)   = rmdi
          dtice(j)      = rmdi
        ENDIF
      END DO
      DO J=1,L2
        IF ( opensea(j) ) THEN
          wtsfc(j)      = 0.0
          wtbase(j)     = 0.0
        ELSE
          wtsfc(j)      = rmdi
          wtbase(j)     = rmdi
        ENDIF
      END DO
C
C-----------------------------------------------------------------------
C
C     2. CALL SLBHCADJ TO ADJUST HEAT CONVERGENCES
C        IN NON-CALIBRATION EXPERIMENT
C
C
      HALFROWS = JROWS/2
      NPOINTS  = ICOLS*HALFROWS
      SPOINTS  = L1 - NPOINTS
      NPTSP1   = NPOINTS+1
      SPTS1    = NPTSP1+ICOLS
C
      IF (.NOT.CALIB) THEN
C
C     ** NORTHERN HEMISPHERE
        CALL SLBHCADJ(L1,NPOINTS,
     +                ADJHCONV(1),
     +                COS_P_LATITUDE(1),
     +                ICY(1),
     +                HCLIMIT,
     +                OPENSEA(1))
C
C     ** SOUTHERN HEMISPHERE
        CALL SLBHCADJ(L1,NPOINTS,
     +                ADJHCONV(SPTS1),
     +                COS_P_LATITUDE(SPTS1),
     +                ICY(SPTS1),
     +                HCLIMIT,
     +                OPENSEA(SPTS1))
C
      ENDIF
C
C-----------------------------------------------------------------------
C
C
C     3. COMPUTE ATMSFLUX AND SNOWOUT FROM ATMOSPHERE FLUXES AND
C        SNOWRATES.
C     IT SHOULD BE POINTED OUT THAT AT SEA-ICE POINTS, THE
C     RADIATIVE FLUXES, THE SENSIBLE
C     HEAT FLUX AND THE EVAPORATION WERE ALREADY WEIGHTED BY THE
C     FRACTIONAL AREA OF LEADS WHEN THEY WERE DIAGNOSED, SO NO
C     SPECIAL CODE IS NECESSARY HERE.
C
      IF (l_therm.or.l_idyn.or.l_idrif) then
        DO J=1,L2
          IF (.NOT.LAND(J)) THEN
            ATMSFLUX(J) = (SOLARIN(J)-BLUEIN(J)) + LONGWAVE(J)
     +                  - ( SENSIBLE(J) + LC*EVAP(J) )
            LEADFLUX(J) = ATMSFLUX(J) * AICEATM(J)
            ATMSFLUX(J) = ATMSFLUX(J) * ( 1.0 - AICEATM(J) )
     +                  + BLUEIN(J)
          ENDIF
        END DO
      ELSE
        DO J=1,L2
          ATMSFLUX(J) = SOLARIN(J) + LONGWAVE(J)
     +                 - ( SENSIBLE(J) + LC*EVAP(J) )
          LEADFLUX(J) = ATMSFLUX(J)
        END DO
      ENDIF
C
      DO J=1,L2
          SNOWRATE(J) = SNOWLS(J) + SNOWCONV(J)
          SNOWSLAB(J) = SNOWRATE(J)
          IF (ICY(J)) THEN
            SNOWSLAB(J) = 0.0
            SNOWLEAD(J) = SNOWRATE(J)
          ENDIF
      END DO
C
C
C
C-----------------------------------------------------------------------
C
C    4. CALL TRANSA2S TO CONVERT ATMOSPHERE VARIABLES TO FORMAT
C       REQUIRED FOR SLAB MODEL
C
        CALL TRANSA2S(L1,
     +                L1,
     +                L_THERM,
     +                LAND,
     +                TSTARATM,
     +                SLABTEMP,
     +                HICEATM,
     +                HICESLB,
     +                HICEMIN,
     +                HSNOWATM,
     +                HSNOWSLB,
     +                AICEATM,
     +                AICESLB,
     +                AICEMIN,
     +                TCLIM,
     +                TCLIMC,
     +                HCLIM)
C
C-----------------------------------------------------------------------
C
C     5. CALL UMSLAB TO UPDATE SLAB OCEAN TEMPERATURE
C
C
        CALL UMSLAB(
!========================== COMDECK ARGOINDX ==========================
     &  J_1, J_2, J_3, J_JMT, J_JMTM1, J_JMTM2, J_JMTP1, JST, JFIN 
     &, J_FROM_LOC, J_TO_LOC, JMT_GLOBAL, JMTM1_GLOBAL, JMTM2_GLOBAL  
     &, JMTP1_GLOBAL, J_OFFSET, O_MYPE, O_EW_HALO, O_NS_HALO        
     &, J_PE_JSTM1, J_PE_JSTM2, J_PE_JFINP1, J_PE_JFINP2, O_NPROC,
     &  imout,jmout,J_PE_IND_MED,NMEDLEV,lev_med,lev_hud,imout_hud,
     &  jmout_hud,J_PE_IND_HUD,med_topflow,
     &              ICY,
     +              ATMSFLUX,
     +              ADJHCONV,
     +              SNOWSLAB,
     +              CARYHEAT,
     +              SLABTEMP,
     +              wtsfc,
     +              wtbase,
     +              LAND,
     +              AICESLB,
     +              NEWICE,
     +              OIFLUX,
     +              DTADV,DTDIFF,
     +              L1,L1,DT,DZ1,u_field,
     +              JROWS,
     +              JROWSM1,
     +              ICOLS,
     +              TCLIMC,
     +              HCLIM,
     +              CALIB,
     +              L_THERM,L_IDYN,L_IDRIF,
     +              L_SLBADV,LGLOBAL,
     +              EDDYDIFF,
     +              DELTA_LONG,DELTA_LAT,BASE_LAT,
     +              AH,
     +              OPENSEA,
     +              ucurrent,
     +              vcurrent,
     +              COS_P_LATITUDE,COS_U_LATITUDE,
     +              SEC_P_LATITUDE,SIN_U_LATITUDE)
C
C-----------------------------------------------------------------------
C
C     6. CALL SLAB_ICEDRIFT OR SLAB_ICEDYN TO UPDATE SEA_ICE VARIABLES
C        (DYNAMIC INCREMENTS)
C
      IF (L_IDRIF.OR.L_IDYN) THEN
C Initialise diagnostic arrays to hold increments
        DO J=1,L2
          IF (.NOT. LAND(J)) THEN
            AINC_DYN(J)  = AICESLB(J)
            HINC_DYN(J)  = HICESLB(J)
            HSINC_DYN(J) = HSNOWSLB(J)*AICESLB(J)
          ENDIF
        END DO
      ENDIF
      IF (L_IDRIF) THEN
C
        call slab_icedrift(
!========================== COMDECK ARGOINDX ==========================
     &  J_1, J_2, J_3, J_JMT, J_JMTM1, J_JMTM2, J_JMTP1, JST, JFIN 
     &, J_FROM_LOC, J_TO_LOC, JMT_GLOBAL, JMTM1_GLOBAL, JMTM2_GLOBAL  
     &, JMTP1_GLOBAL, J_OFFSET, O_MYPE, O_EW_HALO, O_NS_HALO        
     &, J_PE_JSTM1, J_PE_JSTM2, J_PE_JFINP1, J_PE_JFINP2, O_NPROC,
     &  imout,jmout,J_PE_IND_MED,NMEDLEV,lev_med,lev_hud,imout_hud,
     &  jmout_hud,J_PE_IND_HUD,med_topflow,
     &   l1,l2,icols,jrows,jrowsm1,land
     &  ,LGLOBAL,aicemin,amxnorth,amxsouth,Ah_ice
     &  ,delta_lat,delta_long,base_lat,dt,cos_u_latitude
     &  ,cos_p_latitude,sec_p_latitude,sin_u_latitude
     &  ,wsx,wsy,ucurrent,vcurrent,aiceslb,hiceslb,hsnowslb,icy,newice
     &  ,hinc_diff,hinc_adv,hsinc_adv,areas
     &  )
      ELSEIF (L_IDYN) THEN
        call slab_icedyn(
!========================== COMDECK ARGOINDX ==========================
     &  J_1, J_2, J_3, J_JMT, J_JMTM1, J_JMTM2, J_JMTP1, JST, JFIN 
     &, J_FROM_LOC, J_TO_LOC, JMT_GLOBAL, JMTM1_GLOBAL, JMTM2_GLOBAL  
     &, JMTP1_GLOBAL, J_OFFSET, O_MYPE, O_EW_HALO, O_NS_HALO        
     &, J_PE_JSTM1, J_PE_JSTM2, J_PE_JFINP1, J_PE_JFINP2, O_NPROC,
     &  imout,jmout,J_PE_IND_MED,NMEDLEV,lev_med,lev_hud,imout_hud,
     &  jmout_hud,J_PE_IND_HUD,med_topflow,
     &   icols,jrows,jrowsm1,LGLOBAL,delta_lat,delta_long,dt
     &  ,amxsouth,amxnorth,aicemin,Pstar_ice_strength
     &  ,kappa_ice_strength,cdw,tol_ifree
     &  ,nmax_ifree,weight_ifree,tol_icav,nmax_icav
     &  ,land,cos_u_latitude,cos_p_latitude,sec_p_latitude
     &  ,sin_u_latitude,coriolis,wsx,wsy,ucurrent,vcurrent
     &  ,aiceslb,hiceslb,hsnowslb,icy,newice,opensea,uice,vice
     &  ,pmax,pressure
     &  )
      ENDIF
C
C Copy increments into diagnostic arrays
C
      IF (L_IDRIF.OR.L_IDYN) THEN
        area  = 0.0
        DO J=1,L2
          IF (.NOT. LAND(J)) THEN
            AINC_DYN(J)  = AICESLB(J)  - AINC_DYN(J)
            HINC_DYN(J)  = HICESLB(J)  - HINC_DYN(J)
            HSINC_DYN(J) = HSNOWSLB(J)*AICESLB(J) - HSINC_DYN(J)
            area=area+areas(j)
          ENDIF
        END DO
      ENDIF
C
C-----------------------------------------------------------------------
C
C     7. CALL SLABICE TO UPDATE SEA_ICE VARIABLES
C        NEED SEPARATE CALLS FOR NORTHERN AND SOUTHERN HEMISPHERES
C
C Initialise diagnostic arrays to hold increments
C
      DO J=1,L2
        IF (.NOT. LAND(J)) THEN
          AINC_THERM(J)  = AICESLB(J)
          HINC_THERM(J)  = HICESLB(J)
          HSINC_THERM(J) = HSNOWSLB(J)*AICESLB(J)
          DTICE(J)       = SLABTEMP(J)
        ENDIF
      END DO
C     ** NORTHERN HEMISPHERE
        CALL SLABICE(ICY(1),
     +               LEADFLUX(1),
     +               NEWICE(1),
     +               OIFLUX(1),
     +               ATMSFLUX(1),
     +               ADJHCONV(1),
     +               HICESLB(1),
     +               HSNOWSLB(1),
     +               SNOWLEAD(1),
     +               SUBLIMZ(1),
     +               AICESLB(1),
     +               CARYHEAT(1),
     +               TOPMELTZ(1),
     +               BOTMELTZ(1),
     +               SLABTEMP(1),
     +               L1,NPOINTS,DT,DZ1,H0,
     +               L_THERM,L_IDYN,L_IDRIF,
     +               AMXNORTH,AICEMIN,
     +               TCLIMC(1),HCLIM(1),CALIB)
C
C     ** SOUTHERN HEMISPHERE
        CALL SLABICE(ICY(NPTSP1),
     +               LEADFLUX(NPTSP1),
     +               NEWICE(NPTSP1),
     +               OIFLUX(NPTSP1),
     +               ATMSFLUX(NPTSP1),
     +               ADJHCONV(NPTSP1),
     +               HICESLB(NPTSP1),
     +               HSNOWSLB(NPTSP1),
     +               SNOWLEAD(NPTSP1),
     +               SUBLIMZ(NPTSP1),
     +               AICESLB(NPTSP1),
     +               CARYHEAT(NPTSP1),
     +               TOPMELTZ(NPTSP1),
     +               BOTMELTZ(NPTSP1),
     +               SLABTEMP(NPTSP1),
     +               L1,SPOINTS,DT,DZ1,H0,
     +               L_THERM,L_IDYN,L_IDRIF,
     +               AMXSOUTH,AICEMIN,
     +               TCLIMC(NPTSP1),HCLIM(NPTSP1),CALIB)
C
C Copy increments into diagnostic arrays
C
      DO J=1,L2
        IF (.NOT. LAND(J)) THEN
          AINC_THERM(J)  = AICESLB(J)  - AINC_THERM(J)
          HINC_THERM(J)  = HICESLB(J)  - HINC_THERM(J)
          HSINC_THERM(J) = HSNOWSLB(J)*AICESLB(J) - HSINC_THERM(J)
          DTICE(J)       = ( SLABTEMP(J) - DTICE(J) ) / DT
        ENDIF
      END DO
C
C
C-----------------------------------------------------------------------
C
C    8. CALL TRANSS2A TO CONVERT SLAB VARIABLES TO FORMAT
C       REQUIRED FOR ATMOSPHERE
C
        CALL TRANSS2A(L1,L1,
     +                LAND,
     +                SLABTEMP,
     +                TSTARATM,
     +                AICESLB,
     +                AICEATM,
     +                HICESLB,
     +                HICEATM,
     +                HICEMIN,
     +                HSNOWSLB,
     +                HSNOWATM,
     +                AICEMIN)
C
C-----------------------------------------------------------------------
C
C    9. FOR CALIBRATION EXPERIMENT INITIALISE HEATCONV
C
      IF (CALIB) THEN
      DO J=1,L2
          HEATCONV(J) = ADJHCONV(J)
      END DO
      ENDIF
C
      RETURN
      END
