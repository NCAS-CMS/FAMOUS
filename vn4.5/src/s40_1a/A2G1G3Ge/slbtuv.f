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
C*LL
CLL   SUBROUTINE SLAB_T_UV
CLL   -------------------
CLL
CLL   INTERPOLATES U AND V CURRENTS ONTO ARAKAWA C GRID
CLL   BEFORE CALLING SLAB TEMPERTURE ADVECTION
CLL   CALLED FROM UMSLAB
CLL
CLL   IT CAN BE COMPILED BY CFT77, BUT DOES NOT CONFORM TO ANSI
CLL   FORTRAN77 STANDARDS, BECAUSE OF THE INLINE COMMENTS.
CLL
CLL   ALL QUANTITIES IN THIS ROUTINE ARE IN S.I. UNITS UNLESS
CLL   OTHERWISE STATED.
CLL
CLL   WRITTEN BY R.E.CARNELL (05/07/94)
CLL
CLL
CLL  MODEL            MODIFICATION HISTORY SINCE INSERTION IN UM 3.4:
CLL VERSION  DATE
CLL
CLL   4.0           Added vertical SST advection. R.CArnell(J.Crossley)
!LL   4.4  04/08/97 Add missing ARGOINDX to various argument lists.
!LL                 D. Robinson.
CLL
CLL    ADHERES TO THE STANDARDS OF MET. DYNAMICS SUBROUTINES.
CLLEND---------------------------------------------------------------
C*L
      subroutine slab_t_uv(
!========================== COMDECK ARGOINDX ==========================
     &  J_1, J_2, J_3, J_JMT, J_JMTM1, J_JMTM2, J_JMTP1, JST, JFIN 
     &, J_FROM_LOC, J_TO_LOC, JMT_GLOBAL, JMTM1_GLOBAL, JMTM2_GLOBAL  
     &, JMTP1_GLOBAL, J_OFFSET, O_MYPE, O_EW_HALO, O_NS_HALO        
     &, J_PE_JSTM1, J_PE_JSTM2, J_PE_JFINP1, J_PE_JFINP2, O_NPROC,
     &  imout,jmout,J_PE_IND_MED,NMEDLEV,lev_med,lev_hud,imout_hud,
     &  jmout_hud,J_PE_IND_HUD,med_topflow,
C in : model data and atmos/ancillary fields.
     & l1,l2,icols,jrows,jrowsm1,landmask
     &,lglobal,u_field,dz1
     &,delta_lat,delta_long,base_lat,timestep
     &,cos_u_latitude,cos_p_latitude
     &,sec_p_latitude
     &,sin_u_latitude
     &,ucurrent,vcurrent
C
C inout : primary variables
     &,slabtemp
     &,opensea
C
C out   : diagnostics
     &,wtsfc
     &,wtbase
     & )
C
      implicit none
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




      integer
     & l1                           ! in size of data vectors
     &,l2                           ! in amount of data to be processed
     &,icols                        ! in number of columns EW
     &,jrows                        ! in number of rows NS
     &,jrowsm1                      ! in number of rows NS - 1
     &,dz1                          ! in depth of slab ocean
     &,u_field                      ! in points in u_field
      logical
     & lglobal                      ! in true if model is global
     &,landmask(icols,jrows)        ! in mask true at land points
     &,opensea(icols,jrows)         ! in true if open sea (no ice)
      real
     & delta_lat                    ! in meridional grid spacing deg
     &,delta_long                   ! in zonal grid spacing in degrees
     &,base_lat                     ! in base latitude in degrees
     &,timestep                     ! in slab timestep in seconds
     &,cos_p_latitude(icols,jrows)  ! in cosine of latitude on p grid
     &,cos_u_latitude(icols,jrowsm1)! in cosine of latitude on uv grid
     &,sec_p_latitude(icols,jrows)  ! in secont of latitude on p grid
     &,sin_u_latitude(icols,jrowsm1)! in sine of latitude on uv grid
     &,ucurrent(icols,jrowsm1)      ! in zonal surface current (M/S)
     &,vcurrent(icols,jrowsm1)      ! in meridional sfc current (M/S)
     &,slabtemp(icols,jrows)        ! inout slab ocean temperature C
     &,wtsfc(icols,jrows)           ! out w x slab temp at surface
     &,wtbase(icols,jrows)          ! out w x slab temp at base
C
C Global UM parameters
C*L------------------COMDECK C_A----------------------------------------
C A IS MEAN RADIUS OF EARTH
      REAL A

      PARAMETER(A=6371229.)
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

C
C variables local to this subroutine are now defined
C
      integer
     & i,j                          ! loop counters
     &,icolsm1                      ! number of tracer columns - 1
      logical
     & ocean(icols,jrows)           ! true for ocean points on p grid
      real
     & uv                           ! workspace scalar
     &,dlat_rad                     ! grid spacing in radians
     &,dlon_rad                     ! grid spacing in radians
     &,tdiff                        ! slabt_old-slabtemp
C
      real
     & hmask(icols,jrows)           ! 0.0 for land 1.0 for sea at p pts
     &,tmask(icols,jrows)           ! 1.0 for opensea 0.0 ice/land p pts
     &,umask(icols,jrowsm1)         ! 0.0 for uv land 1.0 for sea uv pts
C
      real
     & ucurrent_c(icols,jrowsm1)    ! u current on C grid
     &,vcurrent_c(icols,jrowsm1)    ! v current on C grid
     &,slabt_old(icols,jrows)       ! initial slabtemp
     &,slabt_work(icols,jrows)      ! slabtemp  (no mdi)
C*
C start executable code
C
C initialise various constants.
      icolsm1  = icols-1
      dlat_rad = delta_lat * pi_over_180
      dlon_rad = delta_long * pi_over_180
C
C First set up land sea and ice-free sea masks
C
      do j = 1,jrows
        do i = 1,icols
          ocean(i,j) = .not.landmask(i,j)
          hmask(i,j) = 0.0
          if (ocean(i,j))       hmask(i,j) = 1.0
          tmask(i,j) = 0.0
          if (opensea(i,j))     tmask(i,j) = 1.0
          slabt_work(i,j) = slabtemp(i,j)
          if (slabt_work(i,j).eq.rmdi) slabt_work(i,j) = 0.0
        end do
      end do
C
C Calculate Arakawa B grid velocity mask.
      do j = 1,jrowsm1
        do i = 1,icolsm1
          umask(i,j) = 1.0
          uv = hmask(i,j)+hmask(i+1,j)+hmask(i,j+1)+hmask(i+1,j+1)
          if (uv.lt.3.5) umask(i,j) = 0.0
        end do
        if (lglobal) then
          umask(icols,j) = 1.0
          uv = hmask(icols,j)+hmask(icols,j+1)+hmask(1,j)+hmask(1,j+1)
          if (uv.lt.3.5) umask(icols,j) = 0.0
         else
          umask(icols,j) = 0.0     ! what should i do here ?
        endif
      end do
      do j=1,jrowsm1
        do i=1,icols
          ucurrent(i,j) = ucurrent(i,j)*umask(i,j)
          vcurrent(i,j) = vcurrent(i,j)*umask(i,j)
        end do
      end do
c
C Interpolate currents to C grid.
c
      call uv_to_cu(
!========================== COMDECK ARGOINDX ==========================
     &  J_1, J_2, J_3, J_JMT, J_JMTM1, J_JMTM2, J_JMTP1, JST, JFIN 
     &, J_FROM_LOC, J_TO_LOC, JMT_GLOBAL, JMTM1_GLOBAL, JMTM2_GLOBAL  
     &, JMTP1_GLOBAL, J_OFFSET, O_MYPE, O_EW_HALO, O_NS_HALO        
     &, J_PE_JSTM1, J_PE_JSTM2, J_PE_JFINP1, J_PE_JFINP2, O_NPROC,
     &  imout,jmout,J_PE_IND_MED,NMEDLEV,lev_med,lev_hud,imout_hud,
     &  jmout_hud,J_PE_IND_HUD,med_topflow,
     &     ucurrent,ucurrent_c,jrowsm1,icols)
      call uv_to_cv(
!========================== COMDECK ARGOINDX ==========================
     &  J_1, J_2, J_3, J_JMT, J_JMTM1, J_JMTM2, J_JMTP1, JST, JFIN 
     &, J_FROM_LOC, J_TO_LOC, JMT_GLOBAL, JMTM1_GLOBAL, JMTM2_GLOBAL  
     &, JMTP1_GLOBAL, J_OFFSET, O_MYPE, O_EW_HALO, O_NS_HALO        
     &, J_PE_JSTM1, J_PE_JSTM2, J_PE_JFINP1, J_PE_JFINP2, O_NPROC,
     &  imout,jmout,J_PE_IND_MED,NMEDLEV,lev_med,lev_hud,imout_hud,
     &  jmout_hud,J_PE_IND_HUD,med_topflow,
     &     vcurrent,vcurrent_c,jrowsm1,icols)
c
C Copy initial slabt temp to workspace
c
      do j=1,jrows
        do i=1,icols
          slabt_old(i,j) = slabtemp(i,j)
        end do
      end do
C
C Call slab_temp_advect to advect slab temperature.
C
      call slab_temp_advect(
     & L1,u_field,landmask
     &,icols,jrows,jrowsm1,lglobal,dlat_rad,dlon_rad,timestep,a,dz1
     &,ucurrent_c,vcurrent_c,tmask,opensea,cos_p_latitude
     &,sec_p_latitude,cos_u_latitude
     &,sin_u_latitude,slabt_work,wtsfc,wtbase
     &          )
C
C copy work variables to primary space
C
      do j=1,jrows
        do i=1,icols
          if (tmask(i,j) .eq. 1.0) slabtemp(i,j)=slabt_work(i,j)
        end do
      end do
      return
      end
