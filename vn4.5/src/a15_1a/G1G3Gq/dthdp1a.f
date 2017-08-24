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
CLL Subroutine dthe_dp -----------------------------------------------
CLL
CLL Purpose: To compute dtheta_dp on model half-levels
CLL          for use in potential vorticity calculation.
CLL          These are the natural definition of static stability
CLL          given that theta is held on model levels.
CLL          The subroutine calculates both dtheta/dp and the
CLL          half-levels (called e_levels) between the top and
CLL          bottom full levels. Therefore there are p_levels-1
CLL          e_levels and dtheta/dp's.
CLL
CLL Not suitable for single column use.
CLL
CLL  13/9/93 Written by Terry Davies
CLL
CLL  Model            Modification history:
CLL version  Date
CLL   3.3   14/12/93  Original version
CLL
CLL Programming Standard: UM DOC Paper3, Version 4 (05/02/92)
CLL
CLL Logical Component Covered: D415
CLL
CLL System Task: D4
CLL
CLL Documentation: U.M.D.P No 13. Derivation and Calculation of
CLL                Unified Model Potential Vorticity.
CLL                by Simon Anderson and Ian Roulstone.
CLL
CLLEND

C*L ARGUMENTS: ---------------------------------------------------------
      subroutine dthe_dp
     1                  (pstar,theta,p_field,p_levels,
     2                   ak,bk,akh,bkh,n_levels,
     3                   e_levels,dthe_dph)

      implicit none

C Input variables ------------------------------------------------------

      integer
     & p_field                 !IN    Size of field on pressure points.
     &,p_levels                !IN    Number of pressure levels.
     &,n_levels                !IN    Number of half-levels (p_levels-1)
C                                     for dtheta/dp and e_levels.

      real
     & pstar(p_field)          !IN    Surface pressure field.
     &,theta(p_field,p_levels) !IN    Theta field on p_levels

      real
     & ak(p_levels)            !IN    A coefficient of hybrid
     &                         !      coordinates at full levels.
     &,bk(p_levels)            !IN    B coefficient of hybrid
     &                         !      coordinates at full levels.
     &,akh(p_levels+1)          !IN    A coefficient of hybrid
     &                         !      coordinates at half levels.
     &,bkh(p_levels+1)         !IN    B coefficient of hybrid
     &                         !      coordinates at half levels.

C Output variables -----------------------------------------------------

      real
     & e_levels(n_levels)    !OUT   Model half-levels over range.
     &,dthe_dph(p_field,n_levels)   !OUT dtheta/dp on half levels


C*----------------------------------------------------------------------
C*L Workspace Usage: 2  arrays are required.
      real
     & pressure(p_field,2)   ! Pressure on model levels
     &,thetae(p_field,2)    ! thetate on model levels
C                            ! 2 levels used, IL, IU used to switch

C*----------------------------------------------------------------------
C*L External subroutine calls:

C*----------------------------------------------------------------------
C*L Call comdecks to get required variables:

C*----------------------------------------------------------------------
C*L Define local variables.
      integer i,j,ii       ! Loop variables.
     & ,il,iu,is           ! Pointers for pressure arrays

C ----------------------------------------------------------------------
CL Section 1 Compute model half-levels over required range
CL ~~~~~~~~~
C ----------------------------------------------------------------------

CL Section 1.1
CL ~~~~~~~~~~~ Half-levels in e_levels

      do i=1,p_levels-1
      ii=i+1
      e_levels(i)=0.00001*akh(ii)+bkh(ii)
      end do

CL Section 1.2 Calculate pressure at bottom of range
CL             Store in pressure(p_field,1)

      do j=1,p_field
      pressure(j,1)=ak(1)+bk(1)*pstar(j)
      end do

C ----------------------------------------------------------------------
CL Section 2 Calculate dtheta/dp on model half-levels
CL ~~~~~~~~~
C ----------------------------------------------------------------------

      is=1
      il=2
      do i=1,p_levels-1
      iu=il
      il=is
      ii=i+1
       do j=1,p_field
       pressure(j,iu)=ak(ii)+bk(ii)*pstar(j)
       dthe_dph(j,i)=(theta(j,ii-1)-theta(j,ii))/
     2           (pressure(j,il)-pressure(j,iu))
       end do
      is=iu
      end do

      return
      end

