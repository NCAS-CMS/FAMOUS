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
C*LL  SUBROUTINE LSP_FOCWWIL--------------------------------------------
!LL
!LL  Purpose: Calculate from temperature the Fraction Of Cloud Water
!LL           Which Is Liquid.
!LL     NOTE: Operates within range 0 to -9 deg.C based upon MRF
!LL           observational analysis. Not robust to changes in TM or T0C
!LL
!LL A.Bushell   <- programmer of some or all of previous code or changes
!LL
!LL  Model
!LL version  Date     Modification history from model version 4.0:
!LL
!LL   4.0    27/09/95 Subroutine created from in-line COMDECK.
!LL
!LL
!LL  Programming standard: Unified Model Documentation Paper No 4,
!LL                        Version 1, dated 12/9/89.
!LL
!LL  Logical component covered: Part of P26.
!LL
!LL  System task:
!LL
!LL  Documentation: Unified Model Documentation Paper No 26: Eq 26.50.
!LL
!LL  Called by components P26, P23.
C*
C*L  Arguments:---------------------------------------------------------
      SUBROUTINE LSP_FOCWWIL(
     & T,POINTS,ROCWWIL
     &)
      IMPLICIT NONE
      INTEGER          ! Input integer scalar :-
     & POINTS          ! IN Number of points to be processed.
      REAL             ! Input real arrays :-
     & T(POINTS)       ! IN Temperature at this level (K).
      REAL             ! Updated real arrays :-
     & ROCWWIL(POINTS) ! OUT Ratio Of Cloud Water Which Is Liquid.
C*L   External subprogram called :-
!     EXTERNAL None.
C*
!-----------------------------------------------------------------------
!  Common, then local, physical constants.
!-----------------------------------------------------------------------
C*L------------------COMDECK C_O_DG_C-----------------------------------
C ZERODEGC IS CONVERSION BETWEEN DEGREES CELSIUS AND KELVIN
C TFS IS TEMPERATURE AT WHICH SEA WATER FREEZES
C TM IS TEMPERATURE AT WHICH FRESH WATER FREEZES AND ICE MELTS
      REAL ZERODEGC,TFS,TM

      PARAMETER(ZERODEGC=273.15,
     &          TFS=271.35,
     &          TM=273.15)
C*----------------------------------------------------------------------

      REAL
     & TSTART          ! Temperature at which ROCWWIL reaches 1.
     &,TRANGE          ! Temperature range over which 0 < ROCWWIL < 1.
      PARAMETER(TSTART=TM,
     &          TRANGE=9.0)
!-----------------------------------------------------------------------
!  Define local scalars.
!-----------------------------------------------------------------------
!  (a) Reals effectively expanded to workspace by the Cray (using
!      vector registers).
      REAL             ! Real workspace. At end of DO loop, contains :-
     & TFOC            ! T(I) within DO loop. Allows routines to call
!                        LSP_FOCWWIL(WORK1, POINTS, WORK1) to save space
!  (b) Others.
      INTEGER I       ! Loop counter (horizontal field index).
!
      DO  I = 1, POINTS
!
        TFOC = T(I)
!-----------------------------------------------------------------------
!L 0. Calculate fraction of cloud water which is liquid (FL),
!L    according to equation P26.50.
!-----------------------------------------------------------------------
        IF (TFOC .LE. (TSTART - TRANGE)) THEN
!       Low temperatures, cloud water all frozen------------------------
          ROCWWIL(I) = 0.0
!
        ELSE IF (TFOC .LT. TSTART) THEN
!       Intermediate temperatures---------------------------------------
          ROCWWIL(I) = (TFOC - TSTART + TRANGE) / TRANGE
!
        ELSE
!       High temperatures, cloud water all liquid-----------------------
          ROCWWIL(I) = 1.0
!
        END IF
!
      END DO ! Loop over points
!
      RETURN
      END
