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
!+ Relaxation routine to correct ice velocities for slab model.
!
! Subroutine Interface:
      SUBROUTINE slab_icecavrx(
     &           imt,jmt,jmtm1,delta_lat,delta_long
     &          ,tol,nmax
     &          ,cos_p_latitude,cos_u_latitude,sec_p_latitude
     &          ,hmask,umask,umc,vmc,ccalc,cavrow
     &          ,ax,ay,bx,by,pmax
     &          ,uice,vice,pressure
     &          )

      IMPLICIT NONE
!
! Description:
!   DYNAMIC SEA ICE MODEL SUBROUTINE TO CORRECT C GRID FREE DRIFT
!   VELOCITIES USING THE CAVITATING FLUID ICE RHEOLOGY.
!
! Method:
! Calculate various constants to improve vectorisation of main loop.
! Loop over maximum number of iterations.
!   Copy initial velocities in this loop to workspace.
!   Loop over each point for which calculations are required.
!     Calculate divergence of velocity field.
!     Calculate pressure fiedl increment required prevent convergence,
! or reduce pressure within a diverging grid square to zero.
!     Compare resulting pressure with ice strength and reduce if
!     necessary.
!     Calculate velocity increments for each face implied by pressure
!     correction, masking land and land boundaries.
!   End loop over grid points.
!   Loop over grid points and find maximum velocity component
!   correction for this iteration. Compare with tolerance.
!   If tolerance exceeds largest velocity correction jump out of loop.
! End loop over maximum iterations and print warning if max itertations
! was required.
!
! Current Code Owner: J.F.Thomson
!
! History:
! Version   Date     Comment
! -------   ----     -------
!  3.4      6/94     Original code. J.Thomson
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
!     INTEGER      ExampleVariable      !Description of variable
!
! Global variables (*CALLed COMDECKs etc...):
C*L------------------COMDECK C_A----------------------------------------
C A IS MEAN RADIUS OF EARTH
      REAL A

      PARAMETER(A=6371229.)
C*----------------------------------------------------------------------


! Subroutine arguments
!   Scalar arguments with intent(in):
      INTEGER imt                      ! number of columns
      INTEGER jmt                      ! number of rows
      INTEGER jmtm1                    ! number of rows minus 1
      REAL    delta_lat                ! zonal grid spacing in radians
      REAL    delta_long               ! merid. grid spacing in radians
      REAL    tol                      ! tolerance for correction scheme
      INTEGER nmax                     ! maximum iterations

!   Array  arguments with intent(in):
      REAL    cos_p_latitude(imt,jmt)  ! cos latitude on p pts
      REAL    cos_u_latitude(imt,jmtm1)! cos latitude on p pts
      REAL    sec_p_latitude(imt,jmt)  ! 1/(cos latitude) on p pts
      REAL    hmask(imt,jmt)           ! 1. for sea 0. for land p pts
      REAL    umask(imt,jmtm1)         ! 1. for sea 0. for land uv pts
      REAL    umc(imt,jmtm1)           ! 1. for sea 0. for land cu pts
      REAL    vmc(imt,jmtm1)         ! 1. for sea 0. for land cv pts
      REAL    ax(imt,jmtm1)            ! cdw * sin(psi) for C x pts
      REAL    ay(imt,jmtm1)            ! cdw * sin(psi) for C y pts
      REAL    bx(imt,jmtm1)            ! mf+cdw*cos(psi) for C x pts
      REAL    by(imt,jmtm1)            ! mf+cdw*cos(psi) for C y pts
     &                                 ! m is gbm ice depth * rhoice
     &                                 ! f is coriolis parameter
      REAL    pmax(imt,jmt)            ! ice strength (max pressure)
      LOGICAL ccalc(imt,jmt)           ! true if dynamics calcs needed
      LOGICAL cavrow(jmt)              ! true if dynamics calcs needed

!   Scalar arguments with intent(InOut):

!   Array  arguments with intent(InOut):
      REAL    uice(imt,jmtm1)          ! x ice vel. component (C grid)
      REAL    vice(imt,jmtm1)          ! y ice vel. component (C grid)
      REAL    pressure(imt,jmt)        ! internal ice pressure

!   Scalar arguments with intent(out):

!   Array  arguments with intent(out):

! Local parameters:
      REAL    zero                     ! 0.0
      PARAMETER ( zero = 0.0 )
      REAL    cos_80_deg               ! cos ( 80 degrees )
      PARAMETER ( cos_80_deg = 0.1736 )

! Local scalars:
      INTEGER imtm1                    ! number of columns minus 1
      INTEGER imtm2                    ! number of columns minus 2
      INTEGER i,j,n,istart             ! loop counters
      INTEGER i1,j1                    ! pointers
      INTEGER iter                     ! iteration count
      INTEGER niter                    ! iteration count for polar rows
      REAL    errm,erru,errv           ! error terms in iteration
      REAL    pcorr,delp               ! pressure corrections
      REAL    rdlat                    ! recip of merid grid spacing
      REAL    rdlon                    ! recip of zonal grid spacing
      REAL    t1                       ! intermediate value in calc
      REAL    pdiff                    ! pressure gradient
      REAL    duij,duip1j,dvij,dvijp1  ! velocity increments
      REAL    div                      ! divergence at a point

! Local dynamic arrays:
      REAL    uwork(imt,jmtm1)         ! work array for u ice velocity.
      REAL    vwork(imt,jmtm1)         ! work array for v ice velocity.
      REAL con1(imt,jmt),con2(imt,jmt),con3(imt,jmt),dmlt(imt,jmt)
      REAL con5(imt,jmt),con6(imt,jmt),con7(imt,jmt),con8(imt,jmt)
      REAL con9(imt,jmt)

! Function & Subroutine calls:
      External timer

!- End of header
      call timer('icecavrx',3)
      imtm1=imt-1
      imtm2=imt-2
      rdlat=1.0/delta_lat
      rdlon=1.0/delta_long

! Initialise arrays for use in cavitating fluid calculations
      do j=1,jmt              ! start j loop
        do i=1,imt            ! start i loop
          con1(i,j)=1.0/a*cos_p_latitude(i,j)
          con2(i,j)=zero
          con3(i,j)=zero
          dmlt(i,j)=zero
          con5(i,j)=zero
          con6(i,j)=zero
          con7(i,j)=zero
          con8(i,j)=zero
          con9(i,j)=zero
        end do                ! end i loop
      end do                  ! end j loop
      do j=1,jmtm1            ! start j loop
        j1=j+1
        do i=1,imt            ! start i loop
          i1=i+1
          if (i1.gt.imt) i1=i1-imt
          if (abs(ax(i,j)).gt.1.0e-4) then
            con2(i,j) = (1.0/ax(i,j))*sec_p_latitude(i1,j1)
            con5(i,j) = (con1(i1,j1)*rdlon)/ax(i,j)
          endif
          if (abs(ay(i,j)).gt.1.0e-4) then
            con3(i,j) = cos_u_latitude(i,j)/ay(i,j)
            con7(i,j) = rdlat/(a*ay(i,j))
          endif
          if (abs(ax(i1,j)).gt.1.0e-4)
     &      con6(i,j) = (con1(i1,j1)*rdlon)/ax(i1,j)
          if (abs(ay(i,j1)).gt.1.0e-4)
     &      con8(i,j) = rdlat/(a*ay(i,j1))
        end do                ! end i loop
      end do                  ! end j loop
      do j=1,jmtm1            ! start j loop
        j1=j+1
        do i=1,imt            ! start i loop
          i1=i+1
          if (i1.gt.imt) i1=i1-imt
          t1 =
     &    ( ( con2(i1,j)+con2(i,j) )
     &    *rdlon*rdlon
     &    + (con3(i,j1)+con3(i,j))*rdlat*rdlat )
          if ( abs(t1) .gt. 1.0e-4 )
     &      con9(i,j) =  a*a*cos_p_latitude(i1,j1)/t1
        end do                ! end i loop
      end do                  ! end j loop

! Cavitating Fluid Correction Scheme
! Repeat correction up to nmax times
      do iter=1,nmax          ! start iter loop
        errm = zero
! Copy initial velocities to workspace
        do j=1,jmtm1          ! start j loop
          do i=1,imt        ! start i loop
            uwork(i,j) = uice(i,j) * umc(i,j)
            vwork(i,j) = vice(i,j) * vmc(i,j)
          end do            ! end i loop
        end do              ! end j loop
! Loop over computational grid
        do j=1,jmt-2        ! start j loop
          j1=j+1
          if (cavrow(j1)) then
            niter=1
            if ( cos_u_latitude(1,j1) .lt. cos_80_deg ) niter = 20
            do n=1,niter      ! start n loop
              do istart=2,3              ! start istart loop
                do i=istart,imtm1,2      ! start i loop
                  i1 = i+1
                  delp = zero
                  pcorr = zero
                  div = zero
                  div = ( (uice(i1,j)-uice(i,j)) * rdlon
     &            + ( vice(i,j1)*cos_u_latitude(i,j1)
     &            - vice(i,j)*cos_u_latitude(i,j) ) * rdlat )
     &            * con1(i1,j1)
                  pcorr = -div * con9(i,j) * hmask(i1,j1)
                  if ( div .le. zero ) then
                    delp = pcorr
                  endif
! Test ice pressure against ice strength, pmax.
                  pdiff = pmax(i1,j1) - ( pressure(i1,j1) + delp )
                  if (pdiff.lt.zero) delp=pmax(i1,j1)-pressure(i1,j1)
                  if (div.gt.zero .and. pressure(i1,j1).le.zero) then
                    delp = zero
                  elseif (div.gt.zero.and.pressure(i1,j1).gt.zero) then
                    delp = max(pcorr,-pressure(i1,j1))
                  endif
! Calculate velocity corrections.
                  duij   = -delp*con5(i,j)
                  duip1j =  delp*con6(i,j)
                  dvij   = -delp*con7(i,j)
                  dvijp1 =  delp*con8(i,j)
! Add corrections to velocities.
                  uice(i,j)  = (uice(i,j)+duij) * umc(i,j)
                  uice(i1,j) = (uice(i1,j)+duip1j) * umc(i1,j)
                  vice(i,j)  = (vice(i,j)+dvij) * vmc(i,j)
                  vice(i,j1) = (vice(i,j1)+dvijp1) * vmc(i,j1)
! Calculate new ice pressure.
                  pressure(i1,j1) = pressure(i1,j1) + delp
! End loop
                end do                 ! end i loop
              end do                 ! end istart loop
            end do                   ! end n loop
          endif
        end do                         ! end j loop
! Check maximum error and jump out of loop if this is within the
! tolerance set above.
        errm = zero
        do j = 1,jmt-2                 ! start j loop
          if (cavrow(j+1)) then
            do i = 2,imtm1               ! start i loop
              erru = abs(uice(i,j)-uwork(i,j))
              errv = abs(vice(i,j)-vwork(i,j))
              errm = max(errm,erru,errv)
            end do                       ! end i loop
          endif
        end do                         ! end j loop
        if ( errm .lt. tol ) go to 999
! End loop over maximum iterations
      end do                           ! end iter loop
      if (iter.ge.nmax) write(6,*)
     &   'Maximum number of iterations exceeded in correction scheme'
 999  continue
      write(6,*) 'Number of iterations in correction scheme = ',iter

      call timer('icecavrx',4)
      return
      end
