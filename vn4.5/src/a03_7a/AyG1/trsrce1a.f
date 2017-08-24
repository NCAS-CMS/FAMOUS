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
C*LL  SUBROUTINE TRSRCE ----------------------------------------------
CLL
CLL  Purpose: Adds source increments to a single level of the aerosol
CLL  field.
CLL
CLL  Suitable for single-column use.
CLL
CLL Pete Clark  <- programmer of some or all of previous code or changes
CLL
CLL  Model            Modification history from model version 3.4:
CLL version  Date
CLL   4.0  05/09/95 Model hour and minute added to call to TRSRCE.
CLL                 Diurnal variation of emissions, hard wired to
CLL                 peak at midday, with 10% amplitude, assumed.
CLL                 Programmer Pete Clark.
CLL   4.1  02/05/96 Amplitude of diurnal variation added to input
CLL                 argument list.             M.Woodage
CLL
CLL  Programming standard: Unified Model Documentation Paper No 3,
CLL                        Version 7, dated 11/3/93.
CLL
CLL
C*L  Arguments:---------------------------------------------------------
      SUBROUTINE TRSRCE(
     & DAKH,DBKH,
     & POINTS,PFIELD,
     & PSTAR,TR,TRSCE,TIMESTEP,I_HOUR,I_MINUTE,AMP,ERROR
     &)
C
      IMPLICIT NONE
      INTEGER
     & POINTS              ! IN No. of gridpoints being processed.
     &,PFIELD              ! IN No. of points in global field (at one
C                          !    vertical level).
     &,I_HOUR            ! IN Local time hour
     &,I_MINUTE          ! IN Local time minute
C
      REAL
     & DAKH,DBKH         ! IN Layer thickness in terms of ak and bk
     &,PSTAR(PFIELD)
     &,TR(PFIELD)        ! INOUT Tracer field
C                               !       (kg per kg air).
     &,TRSCE(PFIELD)     ! IN Tracer source field
C                               !       (kg per kg air per s).
     &,TIMESTEP          ! IN Timestep in seconds
     &,AMP               ! IN Amplitude of diurnal variation of emission
C
      INTEGER ERROR      ! OUT Error return code.
C*L------------------COMDECK C_G----------------------------------------
C G IS MEAN ACCEL DUE TO GRAVITY AT EARTH'S SURFACE
      REAL G

      PARAMETER(G=9.80665)
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
C*L  Workspace usage----------------------------------------------------
C*L  External subroutine called ----------------------------------------
C     None
C     EXTERNAL
C* Local, including SAVE'd, storage------------------------------------
C
C  (a) Scalars effectively expanded to workspace by the Cray (using
C      vector registers).
      REAL
     & DM,FACTOR,TS
      PARAMETER (FACTOR=1.0) ! Factor to multiply source term.
      REAL TIME, TZERO
      PARAMETER (TZERO=12.0) ! Time of maximum emissions
C
C  (b) Others.
      INTEGER I       ! Loop counters: K - vertical level index.
C-----------------------------------------------------------------------
C  Check input arguments for potential over-writing problems.
C-----------------------------------------------------------------------
      ERROR=0
      IF(POINTS.GT.PFIELD)THEN
        ERROR=1
        GOTO 9999
      ENDIF
C
C-----------------------------------------------------------------------
CL Subroutine structure :
CL Loop over field adding source term*timestep/level mass/unit area.
C-----------------------------------------------------------------------
C
C
C
C     Allow for source varying with time of day
C
      TS=1.0
      IF (AMP .GT. 0.0) THEN
        TS = 1.0 +
     &  AMP*COS( (REAL(I_HOUR) + REAL(I_MINUTE)/60.0-TZERO) * PI/12.0)
      ENDIF
      TS = TS * TIMESTEP * FACTOR
      DO I=1,POINTS
        DM = ABS((DAKH + PSTAR(I) * DBKH ) / G) ! DAKH,DBKH -ve
        TR(I) = TR(I) + TRSCE(I)*TS/DM
      ENDDO ! Loop over points
C
 9999 CONTINUE ! Error exit
      RETURN
      END
