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
CLL SUBROUTINE th_pvint -----------------------------------------------
CLL
CLL Purpose: Interpolates fields onto a desired pv surface.
CLL          Note that routine assumes that pv is monotonic.
CLL
CLL Suitable for single column use.
CLL
CLL  Model            Modification history:
CLL Version   Date
CLL   3.1   21/01/93  Written By Simon Anderson.
CLL   3.1   18/01/93  New deck at the release of Version 3.1.
CLL   3.2   28/07/93  Change subroutine name to uppercase for
CLL                   portability.    Tracey Smith
CLL    3.3  14/12/93  Change to use linear interpolation only
CLL                   in vertical.   Terry Davies
CLL
CLL Programming standard UM DOC Paper 3, Version 4(05/02/92) A
CLL
CLL Logical Component Covered: D415
CLL
CLL Project Task: D4
CLL
CLL Documentation: U.M.D.P No.13. Derivation and Calculation of
CLL                Unified Model Potential Vorticity.
CLL                By Simon Anderson and Ian Roulstone.
CLL
CLLEND------------------------------------------------------------------

C*L ARGUMENTS:----------------------------------------------------------
      SUBROUTINE TH_PVINT
     1                   (theta_on_press,pvort_p,
     2                    p_field,theta_pv_p_levs,rmdi,des_pv,
     3                    theta_on_pv)

      implicit none

C Input variables ------------------------------------------------------

      integer
     & p_field                 !IN  Points in horizontal p field.
     &,theta_pv_p_levs         !IN  Number of desired pressure levels.

      real
     & theta_on_press(p_field,theta_pv_p_levs)
     &                         !IN  Calculated field of
     &                         !    theta on pressure levels.
     &,pvort_p(p_field,theta_pv_p_levs)
     &                         !IN  Calculated field of potential
     &                         !    vorticity on pressure levels.

      real
     & rmdi                    !IN  Real missing data indicator.
     &,des_pv                  !IN  Desired pv surface in pv units we
     &                         !    want variables interpolated onto.


C Output variables -----------------------------------------------------

      real
     & theta_on_pv(p_field)    !OUT Interpolated theta field on
     &                         !    pv=des_pv surface.

C*----------------------------------------------------------------------
C*L Workspace usage:
      integer
     & base_level_p(p_field)   ! The pressure level below the desired pv
     &                         ! surface for each atmospheric column
     &                         ! (set to theta_pv_p_levs if Des_pv
     &                         !  is above the top level, or if
     &                         !  level is not found).
     &                         ! Calculated at p-points.
      real
     & des_pv_mks              ! Desired pv surface in mks units we
     &                         ! want variables interpolated onto.

C*----------------------------------------------------------------------
C*L External subroutine calls:   None.

C*----------------------------------------------------------------------
C*L Local variables:
      integer i,j,iii,level    ! Loop counts.

      real zthp1,zthp2,zpv1,zpv2

      logical l_down           ! Logical used to control Do While.

C ----------------------------------------------------------------------
CL Section 1 : Find the value of the base level. This is the largest
CL ~~~~~~~~~   value of the pressure level such that pv(level) < Des_pv,
CL             calculated on p-points.
C ----------------------------------------------------------------------

C Firstly, we must divide the desired pv field by 1,000,000 in
C order to get actual values of potential vorticity.

      des_pv_mks = des_pv * 1.0E-06

      do 110 i = 1,p_field
        base_level_p(i) = theta_pv_p_levs
 110  continue

C Find first level, searching down from above, above which Des_pv
C value can be found. Set to top level if none found.

      do 120 i=1,p_field
        l_down = .true.
        level = theta_pv_p_levs
        do while (l_down)
          level = level - 1
          if(pvort_p(i,level).ne.rmdi .and.
     &       pvort_p(i,level+1).ne.rmdi ) then
            if(pvort_p(i,level).lt.des_pv_mks.and.
     &         pvort_p(i,level+1).gt.des_pv_mks ) then
              base_level_p(i) = level
              l_down = .false.
            end if
          end if
        if (level.le.1) l_down = .false.
        end do
 120  continue

C When this loop is done, base_level_p will be the pressure level
C with the highest value of pv LESS than the desired pv value.
C If for any reason, the desired pv value does not lie within any two
C desired pressure levels then base_level_p is set to theta_pv_p_levs.


C ----------------------------------------------------------------------
CL Section 2 : This section will interpolate variables held on
CL ~~~~~~~~~   the p-grid onto the desired potential vorticity surface.
CL             Used for THETA_ON_PV at each point.
C ----------------------------------------------------------------------

C----  Given Theta as a function of pv and Des_pv lying between
C---- Zpv1 and Zpv2, calculate theta at DEs_pv by fitting
C---- linear Lagrange polynomial and evaluating
C----  Potential Vorticity = Des_pv.
C----
C----  If surface is set to theta_pv_p_levs, then we set the value
C----  to missing data.
C----
      do 210 i=1,p_field

        if(base_level_p(i).eq.theta_pv_p_levs) then
          theta_on_pv(i)=rmdi
        else

          iii=base_level_p(i)
      zpv1=pvort_p(i,iii)
      zthp1=theta_on_press(i,iii)
      zpv2=pvort_p(i,iii+1)
      zthp2=theta_on_press(i,iii+1)

      theta_on_pv(i)=(des_pv_mks-zpv2)*zthp1/(zpv1-zpv2)+
     &               (des_pv_mks-zpv1)*zthp2/(zpv2-zpv1)

        end if

 210  continue

      return
      end

