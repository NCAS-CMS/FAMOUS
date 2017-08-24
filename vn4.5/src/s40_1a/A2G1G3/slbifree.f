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
!+ Subroutine calculates free drift ice velocity for SLAB ocean.
!
! Subroutine Interface:
      SUBROUTINE slab_icefreec(
     & imt,jmt,jmtm1,lglobal
     &,tol,nmax,weight
     &,delta_lat,delta_long
     &,cos_p_latitude,umask,ccalc
     &,umc,vmc,pressure
     &,ax,ay,bx,by
     &,x_stress,y_stress
     &,uice,vice
     & )

      IMPLICIT NONE
!
! Description:
!     DYNAMIC SEA ICE MODEL SUBROUTINE TO CALCULATE MODIFIED FREE
!     DRIFT ICE VELOCITIES IN C GRID WITH GIVEN PRESSURE.
!
! Method:
! Loop over maximum number of iterations.
!   Copy initial velocities for this sweep to workspace.
!   Loop over each point in grid requiring free drift calculation
!   using chequerboard scheme to improve vectorisation.
!     C grid free drift underrelaxation scheme.
!     See documentation (to be released - contact J. Thomson)
!     for details of calculation.
!   End loop over grid.
!   Calculate maximum velocity change in this iteration and compare
!   with tolerance. If tolerance exceeds max change, jump out of
!   loop.
! End loop over max iterations, printing warning if max iterations
! reached
!
! Current Code Owner: J.F.Thomson
!
! History:
! Version   Date     Comment
! -------   ----     -------
! 3.4       6/94     Original code. J.Thomson
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
      INTEGER     imt            ! number of tracer columns.
      INTEGER     jmt            ! number of tracer rows.
      INTEGER     jmtm1          ! number of velocity rows.
      LOGICAL     lglobal        ! true for global models.
      REAL        delta_lat      ! meridional grid spacing in radians.
      REAL        delta_long     ! zonal grid spacing in radians.
      REAL        weight         ! weighting term for C grid calc.
      REAL        tol            ! tolerance for C grid calc.
      INTEGER     nmax           ! maximum iterations

!   Array  arguments with intent(in):
      REAL    cos_p_latitude(imt,jmt) ! cosine of p grid points.
      REAL    umask(imt,jmtm1)        ! in 1.0 for uv land 0.0 for sea.
      REAL    umc(imt,jmtm1)          ! in 1.0 for cu land 0.0 for sea.
      REAL    vmc(imt,jmtm1)          ! in 1.0 for cv land 0.0 for sea.
      LOGICAL ccalc(imt,jmt)          ! true if C grid calcs required.
      REAL    pressure(imt,jmtm1)     ! internal ice pressure.
      REAL    x_stress(imt,jmtm1)     ! zonal stress on ice (C grid)
      REAL    y_stress(imt,jmtm1)     ! merid. stress on ice (C grid)
      REAL    ax(imt,jmtm1)           ! cdw * sin(psi) for C x pts
      REAL    ay(imt,jmtm1)           ! cdw * sin(psi) for C y pts
      REAL    bx(imt,jmtm1)           ! mf + cdw * cos(psi) for C x pts
      REAL    by(imt,jmtm1)           ! mf + cdw * cos(psi) for C y pts
     &                                ! m = gbm ice depth * rhoice
     &                                ! f = coriolis parameter

!   Scalar arguments with intent(InOut):

!   Array  arguments with intent(InOut):
      REAL    uice(imt,jmtm1)         ! zonal sea ice velocity.
      REAL    vice(imt,jmtm1)         ! meridional sea ice velocity.

!   Scalar arguments with intent(out):

!   Array  arguments with intent(out):

! Local parameters:
      REAL    zero                    ! 0.0
      PARAMETER ( zero = 0.0 )
      REAL    fr1                     ! 1/4 as used in relaxation
      PARAMETER ( fr1  = 0.25 )
      REAL    fr2                     ! 1/16 as used in relaxation
      PARAMETER ( fr2  = 0.0625 )

! Local scalars:
      REAL    r1                      ! term in C grid calc
      REAL    r2                      ! term in C grid calc
      REAL    rdet                    ! term in C grid calc
      REAL    det                     ! term in C grid calc
      REAL    errm                    ! error term in C grid calc
      REAL    erru                    ! error term in C grid calc
      REAL    errv                    ! error term in C grid calc
      REAL    dpdx                    ! pressure gradient
      REAL    dpdy                    ! pressure gradient
      INTEGER imtm1                   ! number of columns minus 1
      INTEGER imtm2                   ! number of columns minus 2
      INTEGER i,j,istart,iter         ! loop counters

! Local dynamic arrays:
      REAL    uwork(imt,jmtm1)        ! work array for u ice velocity.
      REAL    vwork(imt,jmtm1)        ! work array for v ice velocity.

! Function & Subroutine calls:
!     External

!- End of header
      imtm1=imt-1
      imtm2=imt-2

! Iterative C grid 'free drift' velocity calculations
! Begin loop over maximum number of iterations
      do iter=1,nmax

! Copy initial velocities to work space
        do j=1,jmtm1
          do i=1,imt
            uwork(i,j) = uice(i,j)
            vwork(i,j) = vice(i,j)
          end do
        end do
!
! Calculate velocity components using w as a weight.
!
        do j=2,jmt-2
          do istart=2,3
          do i=istart,imtm1,2
          if(ccalc(i+1,j+1)) then
            dpdx = ( pressure(i+1,j+1) - pressure(i,j+1) )
     &             / ( a * cos_p_latitude(i+1,j+1) * delta_long )
            dpdy = ( pressure(i+1,j+1) - pressure(i+1,j) )
     &             / ( a * delta_lat )
            r1 = dpdx - x_stress(i,j) - fr1 * bx(i,j)
     &  *(vice(i-1,j) + vice(i-1,j+1) + vice(i,j+1) )
            r2 = dpdy - y_stress(i,j) + fr1 * by(i,j)
     &  *(uice(i+1,j) + uice(i,j-1) + uice(i+1,j-1) )
            det = ( ax(i,j)*ay(i,j) + bx(i,j)*by(i,j)*fr2 )
            if ( abs(det) .lt. 1.0e-6 ) then
              rdet = zero
            else
              rdet = 1. / det
            endif
            uice(i,j) = (weight*rdet*( -ay(i,j)*r1 - bx(i,j)*r2*fr1 )
     &                  + (1.0-weight) * uwork(i,j) ) * umc(i,j)
            vice(i,j) = (weight*rdet*( -ax(i,j)*r2 + by(i,j)*r1*fr1 )
     &                  + (1.0-weight) * vwork(i,j) ) * vmc(i,j)
          else
            uice(i,j)=0.0
            vice(i,j)=0.0
          endif
          end do
          end do
        end do
        do j=2,jmt-2
          if (lglobal) then
            if(ccalc(2,j+1)) then
              dpdx = ( pressure(2,j+1) - pressure(1,j+1) )
     &               / ( a * cos_p_latitude(2,j+1) * delta_long )
              dpdy = ( pressure(2,j+1) - pressure(2,j) )
     &               / ( a * delta_lat )
              r1 = dpdx - x_stress(1,j) - fr1 * bx(1,j)
     &    *(vice(imt,j) + vice(imt,j+1) + vice(1,j+1) )
              r2 = dpdy - y_stress(1,j) + fr1 * by(1,j)
     &    *(uice(2,j) + uice(1,j-1) + uice(2,j-1) )
              det = ( ax(1,j)*ay(1,j) + bx(1,j)*by(1,j)*fr2 )
              if ( abs(det) .lt. 1.0e-6 ) then
                rdet = zero
              else
                rdet = 1. / det
              endif
              uice(1,j) = (weight*rdet*(-ay(1,j)*r1 - bx(1,j)*r2*fr1)
     &                  + (1.0-weight) * uwork(1,j) ) * umc(1,j)
              vice(1,j) = (weight*rdet*(-ax(1,j)*r2 + by(1,j)*r1*fr1)
     &                  + (1.0-weight) * vwork(1,j) ) * vmc(1,j)
            else
              uice(1,j)=0.0
              vice(1,j)=0.0
            endif
            if(ccalc(1,j+1)) then
              dpdx = ( pressure(1,j+1) - pressure(imt,j+1) )
     &               / ( a * cos_p_latitude(1,j+1) * delta_long )
              dpdy = ( pressure(1,j+1) - pressure(1,j) )
     &               / ( a * delta_lat )
              r1 = dpdx - x_stress(imt,j) - fr1 * bx(imt,j)
     &    *(vice(imtm1,j) + vice(imtm1,j+1) + vice(imt,j+1) )
              r2 = dpdy - y_stress(imt,j) + fr1 * by(imt,j)
     &    *(uice(1,j) + uice(imt,j-1) + uice(1,j-1) )
              det = (ax(imt,j)*ay(imt,j) + bx(imt,j)*by(imt,j)*fr2)
              if ( abs(det) .lt. 1.0e-6 ) then
                rdet = zero
              else
                rdet = 1. / det
              endif
          uice(imt,j) = (weight*rdet*(-ay(imt,j)*r1 - bx(imt,j)*r2*fr1)
     &                  + (1.0-weight) * uwork(imt,j) ) * umc(imt,j)
          vice(imt,j) = (weight*rdet*(-ax(imt,j)*r2 + by(imt,j)*r1*fr1)
     &                  + (1.0-weight) * vwork(imt,j) ) * vmc(imt,j)
            else
              uice(imt,j)=0.0
              vice(imt,j)=0.0
            endif
          else
            uice(1,j)=0.0
            vice(1,j)=0.0
            uice(imt,j)=0.0
            vice(imt,j)=0.0
          endif
        end do
! Check maximum error and jump out of loop if this is within the
! tolerance set above.
        errm = zero
        do j = 2,jmt-2
          do i = 1,imt
            erru = abs(uice(i,j)-uwork(i,j))
            errv = abs(vice(i,j)-vwork(i,j))
            if (erru .gt. errm) errm=erru
            if (errv .gt. errm) errm=errv
          end do
        end do
        if (errm .le. tol) go to 888
! End loop over maximum iterations
      end do
      if (iter.ge.nmax) write(6,*)
     &   'Maximum number of iterations exceeded in free drift calc.'
 888  continue
      write(6,*) 'Number of iterations in free drift calc = ',iter
! Copy velocities from adjacent rows into rows 1 and jmtm1
        do i=1,imt
          uice(i,1)     = uice(i,2)*umc(i,1)
          uice(i,jmtm1) = uice(i,jmt-2)*umc(i,jmtm1)
          vice(i,1)     = vice(i,2)*vmc(i,1)
          vice(i,jmtm1) = vice(i,jmt-2)*vmc(i,jmtm1)
        end do
      return
      end
