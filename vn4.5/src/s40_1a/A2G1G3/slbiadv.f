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
!+ C grid upwind ice advection scheme for use in slab models.
!
! Subroutine Interface:
      subroutine slab_ice_advect(
     & icols,jrows,jrowsm1,lglobal
     &,delta_lat,delta_long,timestep,a
     &,uice,vice,hmask,umask,icy
     &,cos_p_latitude,cos_u_latitude,sin_u_latitude
     &,aice,hice,hsnow,areas
     & )

      IMPLICIT NONE
!
! Description:
!   Advects ice thickness, compactness and snow depth on Arakawa
!   C grid for slab ocean model using upwind advection.
!
! Method:
! First copy primary variables to wrokspace to avoid data dependency.
! Each point in the grid is treated in turn with extra loops for the
! first and last rows. The area of each grid square and the length of
! each face is calculated using spherical geometry. The sign of the
! velocities on each face is used to determine the upwind value of
! the primary fields being advected (ice depth, ice fraction and snow
! depth) and increments due to advection accross each face are
! determined (and masked so that land is not advected). These facial
! increments are combined to give a change in each property for the
! grid box. The same process is then carried out for the polar rows.
!
! Current Code Owner: J F Thomson
!
! History:
! Version   Date     Comment
! -------   ----     -------
!  3.4     26/06/94  Original code. J.F.Thomson
!  4.0               Diagnose grid box areas and correct masking.
!                     J.F.Crossley
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!
! System component covered:
! System Task:              P40
!
! Declarations:
!   These are of the form:-
!
! Global variables (*CALLed COMDECKs etc...):

! Subroutine arguments
!   Scalar arguments with intent(in):
      integer icols          ! number of columns.
      integer jrows          ! number of rows.
      integer jrowsm1        ! number of rows minus 1.
      logical lglobal        ! true for a global model.
      real delta_lat         ! zonal grid spacing in degrees.
      real delta_long        ! meridional grid spacing in degrees.
      real timestep          ! timestep in seconds.
      real a                 ! radius of earth in metres.

!   Array  arguments with intent(in):
      real uice(icols,jrowsm1)           ! zonal sea ice velocity.
      real vice(icols,jrowsm1)           ! meridional sea ice velocity.
      real hmask(icols,jrows)            ! 0.0 for ha land 1.0 for sea.
      real umask(icols,jrowsm1)          ! 0.0 for uv land 1.0 for sea.
      logical icy(icols,jrows)           ! true if ice is present.
      real cos_p_latitude(icols,jrows)   ! cos(latitude) at p points
      real cos_u_latitude(icols,jrowsm1) ! cos(latitude) at u points
      real sin_u_latitude(icols,jrowsm1) ! sin(latitude) at u points

!   Scalar arguments with intent(InOut):

!   Array  arguments with intent(InOut):
      real aice(icols,jrows)             ! ice compactness.
      real hice(icols,jrows)             ! ice thickness.
      real hsnow(icols,jrows)            ! snow depth.

!   Scalar arguments with intent(out):

!   Array  arguments with intent(out):
      real areas(icols,jrows)            ! grid box areas

! Local parameters:
      real fract_pole           ! fraction of polar grid area/next row
      parameter ( fract_pole = 0.125 )

! Local scalars:
      real um,vm,ump,vmp        ! velocity components on each box face.
      real qinx,qinxp,qiny,qinyp! advective ice depth increments.
      real qinyp_pole, qiny_pole! advective ice depth incs (poles)
      real qax,qaxp,qay,qayp    ! advective ice fraction increments.
      real qayp_pole, qay_pole  ! advective ice fraction incs (poles)
      real qsx,qsxp,qsy,qsyp    ! advective snow depth increments
      real qsyp_pole, qsy_pole  ! advective snow depth incs (poles)
      real qijh,qija,qijs       ! grid box mean advective incs.
      real area                 ! grid box area in square metres.
      real facey,faceyp,facex,facexp ! lengths of faces in metres.
      logical lcheck            ! true if a surrounding point is icy.
      integer icolsm1           ! number of columns minus 1
      integer icolsm2           ! number of columns minus 2
      integer jrowsm2           ! number of rows minus 2
      integer i,i1,i2           ! loop indices
      integer j                 ! loop index

! Local dynamic arrays:
      real hice_init(icols,jrows) ! initial ice thickness.
      real aice_init(icols,jrows) ! initial ice compactness.
      real hsno_init(icols,jrows) ! initial snow depth.

! Function & Subroutine calls:

!- End of header
! start executable code
!
      icolsm1=icols-1
      icolsm2=icols-2
      jrowsm2=jrows-2
      qinyp_pole=0.0
      qiny_pole=0.0
      qayp_pole=0.0
      qay_pole=0.0
      qsyp_pole=0.0
      qsy_pole=0.0
!
! Copy initial ice thickness, compactness and snow depth to workspace
      do j=1,jrows
        do i=1,icols
          hice_init(i,j) = hice(i,j)
          aice_init(i,j) = aice(i,j)
          hsno_init(i,j) = hsnow(i,j)
        end do
      end do
!
! Loop over grid advecting primary variables.
      do j=1,jrowsm2
        do i=1,icols
          i1=i+1
          i2=i+2
          if (i1.gt.icols) i1=i1-icols
          if (i2.gt.icols) i2=i2-icols
          if ( hmask(i1,j+1).gt.0.1 ) then
            area   = ABS( a**2 * delta_long * (sin_u_latitude(i,j+1)
     &               -sin_u_latitude(i,j)) )
            areas(i1,j+1) = area
          endif
          lcheck = (icy(i,j+1).or.icy(i1,j+1).or.icy(i1,j)
     &              .or.icy(i1,j+2).or.icy(i2,j+1))
          if ( hmask(i1,j+1).gt.0.1 .and. lcheck ) then
!
            facex  = ABS ( a * delta_lat )
            facexp = facex
            facey  = ABS ( a*cos_u_latitude(i,j)*delta_long )
            faceyp = ABS ( a*cos_u_latitude(i,j+1)*delta_long )
!
            um  = uice(i,j)
            vm  = vice(i,j)
            ump = uice(i1,j)
            vmp = vice(i,j+1)
!
            if (um.ge.0.)  qinx  = um*hice_init(i,j+1)*hmask(i,j+1)
            if (um.lt.0.)  qinx  = um*hice_init(i1,j+1)*hmask(i1,j+1)
            if (vm.ge.0.)  qiny  =-vm*hice_init(i1,j+1)*hmask(i1,j+1)
            if (vm.lt.0.)  qiny  =-vm*hice_init(i1,j)*hmask(i1,j)
            if (ump.ge.0.) qinxp =-ump*hice_init(i1,j+1)*hmask(i1,j+1)
            if (ump.lt.0.) qinxp =-ump*hice_init(i2,j+1)*hmask(i2,j+1)
            if (vmp.ge.0.) qinyp = vmp*hice_init(i1,j+2)*hmask(i1,j+2)
            if (vmp.lt.0.) qinyp = vmp*hice_init(i1,j+1)*hmask(i1,j+1)
            if (j.eq.1)    qiny_pole=qiny_pole-qiny*facey
     &                               *timestep/(area*fract_pole)
            if (j.eq.1) areas(i1,j) = area*fract_pole
            if (j.eq.jrowsm2)qinyp_pole=qinyp_pole-qinyp*faceyp
     &                               *timestep/(area*fract_pole)
            if (j.eq.jrowsm2) areas(i1,j+2) = area*fract_pole
!
            if (um.ge.0.)  qax  = um*aice_init(i,j+1)*hmask(i,j+1)
            if (um.lt.0.)  qax  = um*aice_init(i1,j+1)*hmask(i1,j+1)
            if (vm.ge.0.)  qay  =-vm*aice_init(i1,j+1)*hmask(i1,j+1)
            if (vm.lt.0.)  qay  =-vm*aice_init(i1,j)*hmask(i1,j)
            if (ump.ge.0.) qaxp =-ump*aice_init(i1,j+1)*hmask(i1,j+1)
            if (ump.lt.0.) qaxp =-ump*aice_init(i2,j+1)*hmask(i2,j+1)
            if (vmp.ge.0.) qayp = vmp*aice_init(i1,j+2)*hmask(i1,j+2)
            if (vmp.lt.0.) qayp = vmp*aice_init(i1,j+1)*hmask(i1,j+1)
            if (j.eq.1)    qay_pole=qay_pole-qay*facey
     &                               *timestep/(area*fract_pole)
            if (j.eq.jrowsm2)qayp_pole=qayp_pole-qayp*faceyp
     &                               *timestep/(area*fract_pole)
!
            if (um.ge.0.)  qsx  = qax*hsno_init(i,j+1)
            if (um.lt.0.)  qsx  = qax*hsno_init(i1,j+1)
            if (vm.ge.0.)  qsy  = qay*hsno_init(i1,j+1)
            if (vm.lt.0.)  qsy  = qay*hsno_init(i1,j)
            if (ump.ge.0.) qsxp = qaxp*hsno_init(i1,j+1)
            if (ump.lt.0.) qsxp = qaxp*hsno_init(i2,j+1)
            if (vmp.ge.0.) qsyp = qayp*hsno_init(i1,j+2)
            if (vmp.lt.0.) qsyp = qayp*hsno_init(i1,j+1)
!
            qijh = hice(i1,j+1)
            qija = aice(i1,j+1)
            qijs = hsnow(i1,j+1)*qija
!
       hice(i1,j+1) = qijh + ( qinx*facex + qinxp*facexp + qiny*facey
     &      + qinyp*faceyp ) * timestep/area
       aice(i1,j+1) = qija + ( qax*facex + qaxp*facexp + qay*facey
     &      + qayp*faceyp ) * timestep/area
            if (aice(i1,j+1).gt.0.0) then
       hsnow(i1,j+1)= ( qijs + ( qsx*facex + qsxp*facexp + qsy*facey
     &        + qsyp*faceyp ) * timestep/area )/aice(i1,j+1)
              if (j.eq.1)    qsy_pole=qsy_pole-qsy*facey
     &                                *timestep/(area*0.125)
              if (j.eq.jrowsm2)qsyp_pole=qsyp_pole-qsyp*faceyp
     &                                 *timestep/(area*0.125)
            endif
!
          endif
        end do
      end do
!
! Special code for polar rows.
      qiny_pole=qiny_pole/icols
      qinyp_pole=qinyp_pole/icols
      qay_pole=qay_pole/icols
      qayp_pole=qayp_pole/icols
      qsy_pole=qsy_pole/icols
      qsyp_pole=qsyp_pole/icols
! First row.
      do i=1,icols
        i1=i+1
        i2=i+2
        if (i1.gt.icols) i1=i1-icols
        if (i2.gt.icols) i2=i2-icols
        lcheck = (icy(i,1).or.icy(i1,1).or.icy(i1,2)
     &            .or.icy(i2,1))
        if (hmask(i1,1).gt.0.1 .and. lcheck) then
!
          qijh = hice(i1,1)
          qija = aice(i1,1)
          qijs = hsnow(i1,1)*qija
!
          hice(i1,1) = qijh + qiny_pole
          aice(i1,1) = qija + qay_pole
          if (aice(i1,1).gt.0.0)
     &     hsnow(i1,1)= ( qijs + qsy_pole ) / aice(i1,1)
        endif
      end do
! Last row.
      do i=1,icols
        i1=i+1
        i2=i+2
        if (i1.gt.icols) i1=i1-icols
        if (i2.gt.icols) i2=i2-icols
        lcheck = (icy(i,jrows).or.icy(i1,jrows).or.icy(i1,jrowsm1)
     &            .or.icy(i2,jrows))
        if (hmask(i1,jrows).gt.0.1 .and. lcheck) then
!
          qijh = hice(i1,jrows)
          qija = aice(i1,jrows)
          qijs = hsnow(i1,jrows)*qija
!
          hice(i1,jrows) = qijh + qinyp_pole
          aice(i1,jrows) = qija + qayp_pole
          if (aice(i1,jrows).gt.0.0)
     &      hsnow(i1,jrows)= ( qijs + qsyp_pole ) / aice(i1,jrows)
        endif
      end do
!
      return
      end
