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
!    Subroutine GRAVSETT ----------------------------------------------
!
! Purpose: To perform gravitational settlement of tracer particles
!          down to the lowest layer of the model.
!          This version allows tracers to fall through 1 or 2 layers.
!
! Current owners of code:                 S Woodward, M Woodage
!
! History:
! Version    Date     Comment
! -------    ----     -------
!   4.4    03/10/97   Original code        S Woodward, M Woodage
!
! Code description:
!  Language: FORTRAN77 + extensions
!  Programming standard: UMDP 3 Vn 6
!
! System components covered:
!
! System task:
!
!Documentation: Not yet available
!
!-----------------------------------------------------------------------
!
      SUBROUTINE GRAVSETT(
     & PFIELD,NLEVS,TRACFLD,DIAM,RHOP,PSTAR,AK,BK,DELTA_AK,DELTA_BK,T,
     & FIRST_POINT,LAST_POINT,DT,DRYDEP)
!
!
      IMPLICIT NONE
!
      INTEGER NLEVS              !IN number of model levels
      INTEGER PFIELD             !IN number of grid points
      INTEGER FIRST_POINT        !IN first point for calcns to be done
      INTEGER LAST_POINT         !IN last point for calcns to be done
!
      REAL DIAM                  !IN tracer particle diameter
      REAL RHOP                  !IN tracer particle density
      REAL PSTAR(PFIELD)         !IN Pstar
      REAL AK(NLEVS)             !IN A vals on levs
      REAL BK(NLEVS)             !IN B vals on levs
      REAL DELTA_AK(NLEVS)       !IN A(lev+1/2) - A(lev-1/2)
      REAL DELTA_BK(NLEVS)       !IN B(lev+1/2) - B(lev-1/2)
      REAL T(PFIELD,NLEVS)       !IN temperature
      REAL DT                    !IN timestep s
!
      REAL TRACFLD(PFIELD,NLEVS) !IN/OUT tracer field
!
      REAL DRYDEP(PFIELD) !OUT deposition flux from layer2 (kg m-2 s-1)
!

! Include COMDECKS
!
! Parameters and Common blocks
! ------------------------ Comdeck PARVARS -------------------------
! Parameters and common blocks required by the MPP-UM
!========================== COMDECK PARPARM ====================
!   Description:
!
!   This COMDECK contains PARAMETERs for the MPP-UM
!
!   Two sets of parameters are set up -
!     i)  for the MPP-UM itself.
!     ii) for the interface to the Message Passing Software.
!
!   History:
!
!   Model    Date     Modification history
!  version
!   4.1      27/1/96  New comdeck based on first section of
!                     old PARVARS.   P.Burton
!   4.2      21/11/96 Add new field type parameter and
!                     magic number used in addressing to indicate
!                     if a calculation is for local data, or data
!                     on the dump on disk (ie. global data)  P.Burton
!   4.2      18/11/96 Moved MaxFieldSize to comdeck AMAXSIZE and
!                     removed Maxbuf.  P.Burton
!   4.2      18/7/96  Removed some unused variables      P.Burton
!   4.4      11/07/97 Reduced MAXPROC to 256 to save memory  P.Burton
!
! ---------------------- PARAMETERS ---------------------
!
! =======================================================
! Parameters needed for the MPP-UM
! =======================================================

      INTEGER   Ndim_max        ! maximum number of spatial dimensions
      PARAMETER (Ndim_max = 3 ) ! 3d data


      INTEGER
     &   fld_type_p           ! indicates a grid on P points
     &,  fld_type_u           ! indicates a grid on U points
     &,  fld_type_unknown     ! indicates a non-standard grid.
      PARAMETER (
     &   fld_type_p=1
     &,  fld_type_u=2
     &,  fld_type_unknown=-1)

      INTEGER
     &   local_data
     &,  global_dump_data
      PARAMETER (
     &   local_data=1        ! Used in addressing to indicate if
     &,  global_dump_data=2) ! calculation is for a local or
!                            ! global (ie. disk dump) size

! =======================================================
! Parameters needed for the Message Passing Software
! =======================================================


      INTEGER
     &   Maxproc              ! Max number of processors
      PARAMETER (
     &   MAXPROC = 256)

      INTEGER
     &   PNorth       ! North processor address in the neighbour array
     &,  PEast        ! East  processor address in the neighbour array
     &,  PSouth       ! South processor address in the neighbour array
     &,  PWest        ! West  processor address in the neighbour array
     &,  NoDomain     ! Value in neighbour array if the domain has
     &                !  no neighbor in this direction. Otherwise
     &                !  the value will be the tid of the neighbor
      PARAMETER (
     &   PNorth   = 1
     &,  PEast    = 2
     &,  PSouth   = 3
     &,  PWest    = 4
     &,  NoDomain = -1)

      INTEGER
     &   BC_STATIC            ! Static boundary conditions
     &,  BC_CYCLIC            ! Cyclic boundary conditions
      PARAMETER (
     &   BC_STATIC = 1
     &,  BC_CYCLIC = 2)

! ---------------------- End of comdeck PARPARM ---------------------
!========================== COMDECK PARCOMM ====================
!
! *** NOTE : This comdeck requires comdeck PARPARM to be *CALLed
!            first.
!
!   Description:
!
!   This COMDECK contains COMMON blocks for the MPP-UM
!
!
!   Two COMMON blocks are defined:
!     i)  UM_PARVAR holds information required by the
!         Parallel Unified Model itself
!     ii) MP_PARVAR holds information required by the interface to
!         the Message Passing Software used by the PUM
!
!   Key concepts used in the inline documentation are:
!     o GLOBAL data - the entire data domain processed by the UM
!     o LOCAL data - the fragment of the GLOBAL data which is
!       stored by this particular process
!     o PERSONAL data - the fragment of the LOCAL data which is
!       updated by this particular process
!     o HALO data - a halo around the PERSONAL data which forms
!       the LOCAL data
!
!     Acronyms used:
!     LPG - Logical Process Grid, this is the grid of logical
!           processors; each logical processor handles one of the
!           decomposed parts of the global data. It does not
!           necessarily represent a physical grid of processors.
!
!   History:
!
!   4.1      27/1/96  New comdeck based on second section of
!                     old PARVARS.   P.Burton
!   4.2     19/08/96  Removed some unused variables, and added
!                     current_decomp_type variable to allow use
!                     of flexible decompositions.
!                     Added nproc_max to indicate the max. number
!                     of processors used for MPP-UM
!                                                      P.Burton
!
! -------------------- COMMON BLOCKS --------------------
!
! =======================================================
! Common block for the Parallel Unified Model
! =======================================================

      INTEGER
     &   first_comp_pe       ! top left pe in LPG
     &,  last_comp_pe        ! bottom right pe in LPG
     &,  current_decomp_type ! current decomposition type
     &,  Offx                ! halo size in EW direction
     &,  Offy                ! halo size in NS direction
     &,  glsize(Ndim_max)    ! global data size
     &,  lasize(Ndim_max)    ! local data size
     &,  blsizep(Ndim_max)   ! personal p data area
     &,  blsizeu(Ndim_max)   ! personal u data area
     &,  datastart(Ndim_max) ! position of personal data in global data
     &                       !   (in terms of standard Fortran array
     &                       !    notation)
     &,  gridsize(Ndim_max)  ! size of the LPG in each dimension
     &,  gridpos(Ndim_max)   ! position of this process in the LPG
!                            ! 0,1,2,...,nproc_x-1 etc.

      LOGICAL
     &    atbase             ! process at the bottom of the LPG
     &,   attop              ! process at the top of the LPG
     &,   atleft             ! process at the left of the LPG
     &,   atright            ! process at the right of the LPG
! NB: None of the above logicals are mutually exclusive

      COMMON /UM_PARVAR/
     &                  first_comp_pe,last_comp_pe
     &,                 current_decomp_type,Offx, Offy
     &,                 glsize,lasize,blsizep,blsizeu
     &,                 datastart,gridsize,gridpos
     &,                 atbase,attop,atleft,atright

! =======================================================
! Common block for the Message Passing Software
! =======================================================

      INTEGER
     &  bound(Ndim_max)           ! type of boundary (cyclic or static)
     &                            !  in each direction
     &, g_lasize(Ndim_max,0:maxproc)
!                                 ! global copy of local data size
     &, g_blsizep(Ndim_max,0:maxproc)
!                                 ! global copy of personal p data area
     &, g_blsizeu(Ndim_max,0:maxproc)
!                                 ! global copy of personal u data area
     &, g_datastart(Ndim_max,0:maxproc)
!                                 ! global copy of datastart
     &, g_gridpos(Ndim_max,0:maxproc)
!                                 ! global copy of gridpos
     &, nproc                     ! number of processors in current
!                                 ! decomposition
     &, nproc_max                 ! maximum number of processors
     &, nproc_x                   ! number of processors in x-direction
     &, nproc_y                   ! number of processors in y-direction
     &, mype                      ! number of this processor
     &                            !  (starting from 0)
     &, neighbour(4)              ! array with the tids of the four
     &                            ! neighbours in the horizontal plane
     &, gc_proc_row_group         ! GID for procs along a proc row
     &, gc_proc_col_group         ! GID for procs along a proc col
     &, gc_all_proc_group         ! GID for all procs

      COMMON /MP_PARVAR/
     &                  bound
     &,                 g_lasize,g_blsizep,g_blsizeu
     &,                 g_datastart,g_gridpos
     &,                 nproc,nproc_max,nproc_x,nproc_y,mype
     &,                 neighbour,gc_proc_row_group
     &,                 gc_proc_col_group, gc_all_proc_group



! ---------------------- End of comdeck PARCOMM -----------------------
! --------------------- End of comdeck PARVARS ---------------------
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

C*L------------------COMDECK C_G----------------------------------------
C G IS MEAN ACCEL DUE TO GRAVITY AT EARTH'S SURFACE
      REAL G

      PARAMETER(G=9.80665)
C*----------------------------------------------------------------------
!-------------------COMDECK C_SULCHM--------------------------------
! Parameters for Sulphur Cycle Chemistry
      REAL
     &     EVAPTAU,       ! timescale for dissolved SO4 to evaporate
     &     NUCTAU,        ! timescale for accumulation mode particles
!                           to nucleate once they enter a cloud.
     &     DIFFUSE_AIT,   ! diffusion coefficient of Aitken particles
     &     K_SO2OH_HI,            ! high pressure reaction rate limit
     &     K_DMS_OH,              ! reaction rate for DMS+OH  cc/mcl/s
     &     BRAT_SO2,              ! branching ratio for SO2 in DMS oxidn
     &     BRAT_MSA,              ! branching ratio for MSA in DMS oxidn
     &     AVOGADRO,             ! no. of molecules in 1 mole
     &     RMM_H2O2,             ! relative molecular mass H2O2 kg/mole
     &     RMM_AIR,              ! relative molecular mass dry air
     &     RMM_W,                ! relative molecular mass water
     &     RELM_S_H2O2,          ! rel atomic mass sulphur/RMM_H2O2
     &     RELM_S_2N,         ! rel atomic mass Sulphur/2*Nitrogen
     &     PARH,                ! power of temp dependence of K_SO2OH_LO
     &     K1,                  ! parameters for calcn of K_SO2OH_LO
     &     T1,                  !
     &     FC,                  ! parameters for interpolation between
     &     FAC1,                !   LO and HI reaction rate limits
     &     K2,K3,K4,            ! parameters for calcn of K_HO2_HO2
     &     T2,T3,T4,            !
     &     CLOUDTAU,              ! air parcel lifetime in cloud
     &     CHEMTAU,               ! chem lifetime in cloud before oxidn
     &     O3_MIN,            ! min mmr of O3 required for oxidn
     &     THOLD                  ! threshold for cloud liquid water
!
!
      PARAMETER (
     &           EVAPTAU = 300.0,             ! secs  (=5 mins) 
     &             NUCTAU = 30.0,         ! secs
     &       DIFFUSE_AIT = 1.7134E-9,        ! sq m/s
     &        K_SO2OH_HI = 1.5E-12,    ! cc/mcl/s from STOCHEM model
     &           K_DMS_OH = 9.1E-12,      ! cc/mcl/s
     &          BRAT_SO2 = 0.9,   
     &           BRAT_MSA = 1.0-BRAT_SO2,
     &           AVOGADRO = 6.022E23,     ! per mole
     &           RMM_H2O2 = 3.40E-2,      ! kg/mole
     &            RMM_AIR = 2.896E-2,     ! kg/mole
     &              RMM_W = 1.8E-2,       ! kg/mole
     &        RELM_S_H2O2 = 3.206/3.40,
     &           RELM_S_2N = 3.206/2.80,
     &               PARH = 3.3,
     &                K1 = 3.0E-31,    ! (cc/mcl)2/s from STOCHEM
     &                 T1 = 300.0,        ! K
     &                FC = 0.6,        ! from STOCHEM model
     &              FAC1 = 1.0317,  ! 0.75-1.27*LOG10(FC) from STOCHEM
     &                 K2 = 2.3E-13,      ! cc/mcl/s
     &                 K3 = 1.9E-33,      ! (cc/mcl)2/s
     &                 K4 = 1.4E-21,      ! cc/mcl
     &                 T2 = 600.0,        ! K
     &                 T3 = 890.0,        ! K
     &                 T4 = 2200.0,       ! K
     &           CLOUDTAU = 1.08E4,       ! secs (=3 hours)
     &            CHEMTAU = 9.0E2,        ! secs (=15 mins)
     &              O3_MIN = 1.6E-8,    !(kg/kg, equiv. 10ppbv)
     &              THOLD = 1.0E-8        ! kg/kg
     &          )
!
      REAL RAD_AIT,         ! median radius of Aitken mode particles
     &     DIAM_AIT,        !   "    diameter    " 
     &     RAD_ACC,         ! median radius of acccumulation mode
     &     DIAM_ACC,        !   "    diameter    "
     &     CHI,             ! mole fraction of S in particle
     &     RHO_SO4,         ! density of  SO4 particle
     &     SIGMA,           ! standard devn of particle size distn
     &     E_PARM,          ! param relating size distns of Ait & Acc
     &     NUM_STAR         ! threshold concn of accu mode particles
                            !  below which PSI=1
!
      PARAMETER (
     &           RAD_AIT = 24.0E-9,          ! m
     &          DIAM_AIT = 2.0*RAD_AIT,    
     &           RAD_ACC = 95.0E-9,          ! m
     &          DIAM_ACC = 2.0*RAD_ACC,    
     &               CHI = 32.0/132.0,
     &           RHO_SO4 = 1769.0,            ! kg/m3
     &             SIGMA = 1.4,
     &            E_PARM = 0.9398, 
     &          NUM_STAR = 1.0E6             ! m-3
     &          )
!
!*---------------------------------------------------------------------
!
! External subroutines called
      EXTERNAL VGRAV
!
! Local variables
!
      INTEGER K                  !LOC loop counter for levels
      INTEGER J                  !LOC loop counter for points
!
      REAL P(PFIELD)             !LOC  Pressure
      REAL VRHOCDT(PFIELD)       !LOC  v*rho*tracer*deltat @lev
      REAL RHOK2(PFIELD)         !LOC  rho(lev+2)
      REAL RHOK1(PFIELD)         !LOC  rho(lev+1)
      REAL RHOK(PFIELD)          !LOC  rho(lev)
      REAL DZK(PFIELD)           !LOC thickness of layer lev
      REAL DZK1(PFIELD)          !LOC thickness of layer lev+1
      REAL DZK2(PFIELD)          !LOC thickness of layer lev+2
      REAL V(PFIELD,NLEVS) !LOC deposition velocity (vstokes corrected)
      REAL MASSOUT2K2(PFIELD)    !LOC flux falling 2 levs from lev k+2
      REAL MASSOUT1K2(PFIELD)    !LOC flux falling 1 levs from lev k+2
      REAL MASSOUT2K1(PFIELD)    !LOC flux falling 2 levs from lev k+1
      REAL MASSOUT1K1(PFIELD)    !LOC flux falling 1 levs from lev k+1
      REAL MASSOUT2K(PFIELD)     !LOC flux falling 2 levs from lev k
      REAL MASSOUT1K(PFIELD)     !LOC flux falling 1 levs from lev k
      REAL DUMMY1(PFIELD,NLEVS)  !LOC
      REAL DUMMY2(PFIELD,NLEVS)  !LOC
!
!
! Calculate settlement velocity
!
      CALL VGRAV(PFIELD,NLEVS,DIAM,RHOP,PSTAR,AK,BK,T,V,DUMMY1,DUMMY2,
     &           FIRST_POINT,LAST_POINT)
!
! Calculate new tracer mixing ratios
!
! Initialise deposition flux to zero
      DO J = 1,PFIELD
        DRYDEP(J)=0.
      ENDDO
!
! Level 1 (K at start of loop)
!
      DO J = FIRST_POINT,LAST_POINT
        P(J)=AK(1)+BK(1)*PSTAR(J)
        RHOK(J)=P(J)/(R*T(J,1))
        DZK(J)=-(DELTA_AK(1)+DELTA_BK(1)*PSTAR(J))/(RHOK(J)*G)
        MASSOUT2K(J)=0.
        MASSOUT1K(J)=0.
      ENDDO                !J
!
! Level 2 (K+1 at start of loop)
!   NB  deposit tracer direct to ground from lev 2 if V high enough
!
      DO J = FIRST_POINT,LAST_POINT
!
        P(J)=AK(2)+BK(2)*PSTAR(J)
        RHOK1(J)=P(J)/(R*T(J,2))
        DZK1(J)=-(DELTA_AK(2)+DELTA_BK(2)*PSTAR(J))/(RHOK1(J)*G)
!
!   check for deposition :
        IF (V(J,2)*DT .GT. DZK(J)) THEN
!       some tracer deposited onto ground
!
          IF ( V(J,2)*DT .GT. DZK1(J)+DZK(J) ) THEN
!         all deposited to ground
            MASSOUT2K1(J)=RHOK1(J)*TRACFLD(J,2)*DZK1(J)
            MASSOUT1K1(J)=0.
          ELSE IF ( V(J,2)*DT .GT. DZK1(J) ) THEN
!         some deposited to ground, some to layer 1
            MASSOUT2K1(J)=RHOK1(J)*TRACFLD(J,2)*(V(J,2)*DT-DZK(J))
            MASSOUT1K1(J)=RHOK1(J)*TRACFLD(J,2)*DZK1(J)-MASSOUT2K1(J)
          ELSE
!         some deposited to ground, some to layer1, some left in layer2
            MASSOUT2K1(J)=RHOK1(J)*TRACFLD(J,2)*(V(J,2)*DT-DZK(J))
            MASSOUT1K1(J)=RHOK1(J)*TRACFLD(J,2)*DZK(J)
          ENDIF
!
          DRYDEP(J)=MASSOUT2K1(J)/DT
!
        ELSE
!       only falls into layer 1
          MASSOUT2K1(J)=0.
          IF ( V(J,2)*DT .GT. DZK1(J)) THEN
!         all falls into layer 1
            MASSOUT1K1(J)=RHOK1(J)*TRACFLD(J,2)*DZK1(J)
          ELSE
!         some to layer 1 , some left in layer2
            MASSOUT1K1(J)=RHOK1(J)*TRACFLD(J,2)*V(J,2)*DT
          ENDIF
!
        ENDIF
!
      ENDDO                 !END J LOOP
!
! Main loop through levels, from bottom up
!
      DO K = 1,NLEVS-2
!
        DO J = FIRST_POINT,LAST_POINT
!
          P(J)=AK(K+2)+BK(K+2)*PSTAR(J)
          RHOK2(J)=P(J)/(R*T(J,K+2))
          DZK2(J)=-(DELTA_AK(K+2)+DELTA_BK(K+2)*PSTAR(J))/
     &             (RHOK2(J)*G)
!
!       Calculate mass of tracer falling between levels
!
!        limit fall to 2 levs
         IF (V(J,K+2)*DT.GT.(DZK1(J)+DZK(J)))
     &       V(J,K+2)=(DZK1(J)+DZK(J))/DT
!
!         check how far tracer falls:
          IF ( V(J,K+2)*DT .GT. DZK1(J) ) THEN
!         it falls through more than 1 layer
             IF ( V(J,K+2)*DT .GT. (DZK2(J)+DZK1(J)) ) THEN
!            all into layer k
               MASSOUT2K2(J)=RHOK2(J)*TRACFLD(J,K+2)*DZK2(J)
               MASSOUT1K2(J)=0.
             ELSE IF ( V(J,K+2)*DT .GT. DZK2(J) ) THEN
!            some into k+1, some into k
               MASSOUT2K2(J)=
     &           RHOK2(J)*TRACFLD(J,K+2)*(V(J,K+2)*DT-DZK1(J))
               MASSOUT1K2(J)=
     &           RHOK2(J)*TRACFLD(J,K+2)*DZK2(J)-MASSOUT2K2(J)
             ELSE
!            some left in k+2, some into k+1, some into k
               MASSOUT2K2(J)=
     &           RHOK2(J)*TRACFLD(J,K+2)*(V(J,K+2)*DT-DZK1(J))
               MASSOUT1K2(J)=RHOK2(J)*TRACFLD(J,K+2)*DZK1(J)
             ENDIF
!
           ELSE
!          falls no more than 1 layer
             MASSOUT2K2(J)=0.
             IF (V(J,K+2)*DT .GT. DZK2(J)) THEN
!            all falls into layer k+1
               MASSOUT1K2(J)=RHOK2(J)*TRACFLD(J,K+2)*DZK2(J)
             ELSE
!            some falls into k+1, some left in k+2
               MASSOUT1K2(J)=RHOK2(J)*TRACFLD(J,K+2)*V(J,K+2)*DT
             ENDIF
!
          ENDIF
!
! Update tracer field
!
          TRACFLD(J,K)=TRACFLD(J,K)+(MASSOUT2K2(J)+MASSOUT1K1(J)-
     &                 MASSOUT2K(J)-MASSOUT1K(J))/(RHOK(J)*DZK(J))
!
! Put k+2 vals in k+1's & k+1's in k's
          MASSOUT1K(J)=MASSOUT1K1(J)
          MASSOUT1K1(J)=MASSOUT1K2(J)
          MASSOUT2K(J)=MASSOUT2K1(J)
          MASSOUT2K1(J)=MASSOUT2K2(J)
          DZK(J)=DZK1(J)
          DZK1(J)=DZK2(J)
          RHOK(J)=RHOK1(J)
          RHOK1(J)=RHOK2(J)
!
        ENDDO            !END J LOOP
!
      ENDDO              !END K LOOP
!
! Top 2 levels
!
      DO J=FIRST_POINT,LAST_POINT
!
         TRACFLD(J,NLEVS-1)=TRACFLD(J,NLEVS-1)+
     &    (MASSOUT1K1(J)-MASSOUT2K(J)-MASSOUT1K(J))/(RHOK(J)*DZK(J))
         TRACFLD(J,NLEVS)=TRACFLD(J,NLEVS)-
     &   (MASSOUT2K1(J)+MASSOUT1K1(J))/
     &                 (RHOK1(J)*DZK1(J))
!
      ENDDO        !J
!
      RETURN
      END
