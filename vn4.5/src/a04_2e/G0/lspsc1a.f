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
C*LL  SUBROUTINE LSP_SCAV-----------------------------------------------
CLL
CLL  Purpose: Scavenge aerosol by large scale precipitation.
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  Date
CLL  3.4  15/08/94  New routine. Pete Clark.
CLL
CLL  Programming standard: Unified Model Documentation Paper No 3,
CLL                        Version 7, dated 11/3/93.
CLL
CLL  Logical component covered: Part of P26.
CLL
CLL  System task:
CLL
CLL  Documentation: Unified Model Documentation Paper No 26.
C*
C*L  Arguments:---------------------------------------------------------
      SUBROUTINE LSP_SCAV(
     & TIMESTEP,POINTS,RAIN,SNOW,AEROSOL
     &)
      IMPLICIT NONE
      INTEGER         ! Input integer scalar :-
     & POINTS         ! IN Number of points to be processed.
      REAL            ! Input real scalar :-
     & TIMESTEP       ! IN Timestep (s).
      REAL            ! Input real arrays :-
     & RAIN(POINTS)   ! IN Rate of rainfall in this layer from
C                     !       above
C*                    !       (kg per sq m per s).
     &,SNOW(POINTS)   ! IN Rate of snowfall in this layer from
C                     !       above
C*                    !       (kg per sq m per s).
      REAL            ! Updated real arrays :-
     & AEROSOL(POINTS) ! INOUT Aerosol mixing ratio
C*L   External subprogram called :-
C     EXTERNAL None
C-----------------------------------------------------------------------
C  Define local scalars.
C-----------------------------------------------------------------------
C  (a) Reals effectively expanded to workspace by the Cray (using
C      vector registers).
      REAL            ! Real workspace.
     & KRAIN,KSNOW
      PARAMETER(KRAIN=1.0E-4,KSNOW=1.0E-4)
      REAL            ! Real workspace.
     & RRAIN,RSNOW
C  (b) Others.
      INTEGER I       ! Loop counter (horizontal field index).
C
C Overall rate = KRAIN*(R) where R is in mm/hr=kg/m2/s*3600.0
      RRAIN=KRAIN*TIMESTEP*3600.0
      RSNOW=KSNOW*TIMESTEP*3600.0
      DO I=1,POINTS
        AEROSOL(I)=AEROSOL(I)/(1.0+RRAIN*RAIN(I)+RSNOW*SNOW(I))
      END DO
      RETURN
      END
