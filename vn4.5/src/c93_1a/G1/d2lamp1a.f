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
!+ 2D Limited area p.d.e. solver, using zero Dirichlet b.c.'s.
!
! Subroutine Interface:
      SUBROUTINE DEL_SQUARED_LAM_P
     & (P_field, row_len, LonOffset, LonPts,
     &  LatOffset, LatPts, EarthRadiusInv, DLonInv, DLatInv,
     &  cos_p_lat, sec_p_lat, cos_v_lat, source, solution)

      IMPLICIT NONE
!
! Description:
!  This routine solves a 2D Poisson equation on a section of a sphere,
!  defined by coordinate boundaries, using a multigrid method. Both the
!  source and solution arrays are on pressure grids.
!  Zero Dirichlet boundary conditions are imposed.
!
! Method:
!  Uses the multigrid solver MG_CNTL developed by Mark Mawson and
!  documented by him under the title "Numerical solution of Elliptic
!  Equations using Multigrid Methods". It is used here with version = 4.
!
!  The first guess for the solution is taken to be an array of zeros.
!
!  Provision is made for solving on a subset of the full limited area
!  grid thus allowing the problem to be solved on a grid for which the
!  multigrid method is efficient.  The subset region is specified by 4
!  variables: LonOffset; LonPts; LatOffset; LatPts.  Ideally, both
!  LonPts & LatPts (the length of the subset grid in each dimension)
!  SHOULD be of the form 2**(l-1) +1:  where l is an integer >= 1. They
!  MUST be of the form m*2**(l-1) +1:  where m is an integer >= 2,
!  and l is an integer >= 1.
!
!  IMPORTANT: This routine expects the source and solution arrays to be
!  of the same extent as the full limited area pressure grid (values of
!  the solution outside the subset region are set to zero).
!
! Current Code Owner: Phil Andrews
! History:      Model
! Date:        Version:  Comment:
!   3/11/93      3.3     Original version. (Phil Andrews).
!   7/12/93      3.3     Altered to call version 2.0 of the multigrid
!                        code. (Phil Andrews).
!   7/12/93      3.3     Split subroutine declaration line to
!                        avoid lexcon making the line too
!                        long.  Tracey Smith
CLL vn4.2  07/11/96:Allow this routine to be used in C93_2A (Farnon)
!
! Code Description:
!   Language: FORTRAN 77 + extensions
!
! System component covered:
! System Task:
!
! Declarations: these are of the form:-
!     INTEGER      ExampleVariable      !Description of variable
!
! Subroutine arguments:
!   Scalar arguments with Intent (In):
      INTEGER        P_field      ! Number of pressure points in full,
                                  ! unsubsetted, LAM grid.
      INTEGER        row_len      ! Row length of full LAM grid.
      INTEGER        LonOffset    ! Longitude offset of start of subset
                                  ! of LAM grid from start of full LAM
                                  ! grid.
      INTEGER        LonPts       ! Longitude extent of subset region.
      INTEGER        LatOffset    ! Latitude offset of start of subset
                                  ! of LAM grid from start of full LAM
                                  ! grid.
      INTEGER        LatPts       ! Latitude extent of subset region.

      REAL           EarthRadiusInv ! 1.0 / radius of the Earth.
      REAL           DLonInv      ! 1.0/ longitude step.
      REAL           DLatInv      ! 1.0/ latitude step.

! Array arguments with Intent (In):
      REAL           cos_p_lat(P_field) ! cos of latitude on full LAM
      REAL           sec_p_lat(P_field) ! sec of latitude on full LAM
      REAL           cos_V_lat(P_field - row_len) ! cos of latitude
                                  ! on V wind component grid.
      REAL           source(P_field)! Array holding the right
                                  ! hand side of Poisson's equn.

! Array arguments with Intent (Out):
      REAL           solution(P_field)! Is on the full LAM grid, but
                                  ! values outside the subset region are
                                  ! set to zero.

! Local parameters:
! N.B. These are all control parameters for the multigrid solver.
      INTEGER        version      ! Selects mode of operation of the
        PARAMETER   (version = 4) ! multigrid solver.
      INTEGER        k_bc         ! Select vertical bc type(s).
        PARAMETER   (k_bc = 1)
      INTEGER        MxnGrds      ! Max number of grids.
        PARAMETER   (MxnGrds = 100)
      INTEGER        levels       ! extent in k direction
        PARAMETER   (levels = 1)  ! i.e. only one level
      INTEGER        MaxIts       ! Max number of multigrid iteratons.
        PARAMETER   (MaxIts = 100)
      INTEGER        Iprint       ! Selects diagnostics from multigrid.
        PARAMETER   (Iprint = 0)
      INTEGER        Ksmooth      ! Contols choice of smoother
        PARAMETER   (Ksmooth = 16)! i-line zebra
      INTEGER        NPre         ! No of smooths on each grid on way
        PARAMETER   (NPre = 3)    ! down.
      INTEGER        NPost        ! No of smooths on each grid on way
        PARAMETER   (NPost = 3)   ! back up.
      INTEGER        NCoarse      ! No of smooths on coarsest grid.
        PARAMETER   (NCoarse = 2)
      INTEGER        KRestrict    ! Selects restriction algorithm
        PARAMETER   (KRestrict = 3) ! (Full weighting).
      INTEGER        NcGc         ! Selects cycle to use.
        PARAMETER   (NcGc = 1)    ! V cycle.

      REAL           Relax        ! Relaxation parameter for Jacobi
        PARAMETER   (Relax = 1)   ! smoother.
      REAL           Tol_Res      ! Sets factor to reduce initial
        PARAMETER   (Tol_Res = 1.0E-8)! residual by.
      REAL           WorstSmoothingRate ! min. acceptable improvement
        PARAMETER   (WorstSmoothingRate = 0.95) ! in PSR in one cycle.

! Local scalars:
      INTEGER        x            ! Loop counter in Longitude direction.
      INTEGER        y            ! Loop counter in Latitude  direction.
      INTEGER        w            ! Index for full LAM arrays

! Local dynamic arrays:
      REAL           Z_Q(1)       ! Z values at k levels.
      REAL           Z_Mid(2)     ! Z values at k + 1/2 levels.
      REAL           rhs(LonPts, LatPts, levels) ! Source terms on
                                  ! subset of LAM grid.
      REAL           Q(LonPts, LatPts, levels) ! On input to mg_cntl is
                                  ! first guess at the solution plus the
                                  ! boundary conditions at the edges
                                  ! (set to zero in this routine).
                                  ! On output is the solution.
      REAL           coeff_A  (LonPts, LatPts, levels)  ! Coeffs of the
      REAL           coeff_B  (LonPts, LatPts, levels)  ! pde being
      REAL           coeff_C1 (LonPts, LatPts, levels)  ! solved. See
      REAL           coeff_C2 (LonPts, LatPts, levels)  ! documentation
      REAL           coeff_DEF(LonPts, LatPts, levels)  ! for the
      REAL           coeff_D  (LonPts, LatPts, levels)  ! multigrid
      REAL           coeff_E  (LonPts, LatPts, levels)  ! solver.
      REAL           coeff_F  (LonPts, LatPts, levels)  !
      REAL           coeff_G  (LonPts, LatPts, levels)  !
      REAL           cosPlat(LonPts, LatPts) ! cos latitude at pressure
                                  ! points on subset of LAM grid.
      REAL           secPlat(LonPts, LatPts) ! sec latitude at pressure.
                                  ! points on subset of LAM grid.
      REAL           cosVlat(LonPts, LatPts) ! cos latitude at
                                  ! V wind grid points on subset of LAM.

! Function & Subroutine calls:
      External mg_cntl
!-
C*
! 1.0 Initialize
      Z_Q(1)   = 0.0
      Z_Mid(1) = 0.0
      Z_Mid(2) = 0.0

! Subset source array. Initalize the coeff_, trig & rhs arrays.
      Do y = 1, LatPts
        Do x = 1, LonPts
          w = (x + LonOffset) + ((y + LatOffset -1) * row_len)

          cosPlat(x, y)      = cos_p_lat(w)
          secPlat(x, y)      = sec_p_lat(w)
          rhs(x, y, 1)       = source(w)
          Q(x, y, 1)         = 0.0
          coeff_A(x, y, 1)   = 1.0
          coeff_B(x, y, 1)   = 1.0
          coeff_C1(x, y, 1)  = 0.0
          coeff_C2(x, y, 1)  = 0.0
          coeff_DEF(x, y, 1) = 0.0
          coeff_D(x, y, 1)   = 0.0
          coeff_E(x, y, 1)   = 0.0
          coeff_F(x, y, 1)   = 0.0
          coeff_G(x, y, 1)   = 0.0

        End do
      End do

      Do y = 1, LatPts
        Do x = 1, LonPts - 1
          w = (x + LonOffset) + ((y + LatOffset -1) * row_len)

          cosVlat(x, y) = cos_V_lat(w)

        End do

        cosVlat(LonPts, y) = 0.0 ! This extra row shouldnt be used!

      End do

! 2.0 Call the multigrid routine
      Call mg_cntl (MxnGrds, LonPts, LatPts, levels,     ! Intent In
     &  MaxIts, Tol_Res, Iprint, Ksmooth,                !   "     "
     &  NPre, NPost, NCoarse, Relax, Krestrict,          !   "     "
     &  NcGc, coeff_A, coeff_B, coeff_C1, coeff_C2,      !   "     "
     &  coeff_DEF, coeff_D, coeff_E, coeff_F, coeff_G,   !   "     "
     &  Q,                                               ! Intent InOut
     &  rhs, cosPlat, secPlat, cosVlat, EarthRadiusInv,  ! Intent In
     &  DLonInv, DLatInv, WorstSmoothingRate,            !   "     "
     &  version, k_bc, Z_Q, Z_Mid)                       !   "     "

! 3.0 Put the answer into the full LAM array.
! Initialize solution to zero:
      Do w = 1, P_field
        solution(w) = 0.0

      End do

! Now copy Q into solution
      Do y = 1, LatPts
        Do x = 1, LonPts
          w = x + LonOffset + ((y + LatOffset -1) * row_len)

          solution(w) = Q(x, y, 1)

        End do
      End do

! Thats it
      Return
      End
