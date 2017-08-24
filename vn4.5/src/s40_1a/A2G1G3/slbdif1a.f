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
CLL
C    SUBROUTINE SLABDIFF
C    -------------------
C
CLL   THIS ROUTINE IS FOR USE WITH THE 'SLAB' OCEAN MODEL ONLY.
CLL
CLL   DEL2 DIFFUSION OF SLAB TEMPERATURE, WITH CHECK FOR STABILITY.
CLL   DIFFUSION COEFFICIENT USED IS DEPENDENT ON RESOLUTION
CLL   SUGGESTED VALUES ARE
CLL
CLL   4.0E4  FOR 5 X 7.5 DEGREE LAT/LONG
CLL   2.0E4  FOR 2.5 X 3.75 DEGREE LAT/LONG
CLL
CLL   THIS ROUTINE FORMS PART OF SYSTEM COMPONENT P40.
CLL   IT CAN BE COMPILED BY CFT77, BUT DOES NOT CONFORM TO ANSI
CLL   FORTRAN77 STANDARDS, BECAUSE OF THE INLINE COMMENTS, THE
CLL   USE OF ENDDO AND DYNAMIC ALLOCATION.
CLL   IT ADHERES TO THE STANDARDS OF DOCUMENTATION PAPER 3, VERSION 5.
CLL
CLL   ALL QUANTITIES IN THIS ROUTINE ARE IN S.I. UNITS UNLESS
CLL   OTHERWISE STATED.
CLL
CLL   CALLED BY: UMSLAB
CLL
CLL   WRITTEN BY C.A.SENIOR (09/9/93)
CLL   MODIFIED BY C.A.SENIOR (27/10/93) to include stability test
CLL   MODIFIED BY C.A.SENIOR (14/12/93) to update to version 3.2
CLL   MODIFIED BY C.A.SENIOR (17/12/93) after review
CLL   MODIFIED BY C.A.SENIOR (25/02/94) after further review
CLL   VERSION NUMBER 1.1
CLL   REVIEWER: W.J.INGRAM
CLL
CLLEND---------------------------------------------------------------
C
      SUBROUTINE SLABDIFF(SLABTEMP,
     +                    OPENSEA,
     +                    L1,L2,
     +                    JROWS,
     +                    ICOLS,
     +                    AHDT,
     +                    DELTA_LONG,DELTA_LAT,BASE_LAT,
     +                    COS_P_LATITUDE,COS_U_LATITUDE,SEC_P_LATITUDE)

      INTEGER L1              ! IN SIZE OF DATA VECTORS
     +,L2                     ! IN AMOUNT OF DATA TO BE PROCESSED
     +,JROWS                  ! IN NO OF ROWS N-S
     +,ICOLS                  ! IN NO OF COLUMNS E-W
      REAL
     + SLABTEMP(ICOLS,JROWS)  ! INOUT HEAT CONTENT OF THE SLAB
     +,DT                     ! IN TIMESTEP FOR UPDATING SLAB MODEL
     +,DELTA_LONG             ! IN EW GRID SPACING (DEGREES)
     +,DELTA_LAT              ! IN NS GRID SPACING (DEGREES)
     +,BASE_LAT               ! IN LATITUDE OF FIRST ROW (DEGREES)
     +,AHDT                   ! DIFFUSION COEFFICENT * TIMESTEP
C
      REAL
     + COS_U_LATITUDE(ICOLS,JROWS)  ! COSINE OF LATITUDE ON U GRID
     +,COS_P_LATITUDE(ICOLS,JROWS)  ! COSINE OF LATITUDE ON P GRID
     +,SEC_P_LATITUDE(ICOLS,JROWS)  ! 1/COS_P_LATITUDE
C
      LOGICAL
     + OPENSEA(ICOLS,JROWS)   ! IN TRUE IF BOX CONTAINS ICE FREE SEA
     +                        ! POINTS, FALSE AT LAND/SEA_ICE POINTS
C
C Include Comdecks
C*L------------------COMDECK C_A----------------------------------------
C A IS MEAN RADIUS OF EARTH
      REAL A

      PARAMETER(A=6371229.)
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
C Local variables
C
      INTEGER
     + ICOLSM1                ! ICOLS MINUS 1
     +,JROWSM1                ! JROWS MINUS 1.
     +,J                      ! LOOPCOUNTER
     +,I                      ! LOOPCOUNTER
      REAL
     + ZMASK(ICOLS,JROWS)     ! MASK,1=OPEN SEA,0=LAND+SEA-ICE
     +,DYT                    ! GRID SPACING N-S
     +,DXT                    ! GRID SPACING E-W
     +,DYTR                   ! 1/DYT
     +,DXT2R                  ! 1/(2*DXT)
     +,DIFFUS(ICOLS,JROWS)    ! DIFFUSION INCREMENT
     +,ROWCOEF                ! COS WEIGHTED COEFF. E-W
     +,COLCOEF_J              ! COS WEIGHTED COEFF. N-S AT ROW J
     +,COLCOEF_JM1            ! COS WEIGHTED COEFF. N-S AT ROW J-1
     +,TEMPA(ICOLS)           ! E-W DIFUSION INCREMENTS
     +,STABLT                 ! STABILITY COEFFICENT
C
      PARAMETER ( STABLT= 0.1) ! STABILITY COEFFICIENT
C
C INITIALISE CONSTANTS.
C
      JROWSM1 = JROWS-1
      ICOLSM1 = ICOLS-1
C
C SET UP REAL OPEN SEA MASK AND EXCLUDE POLAR ROWS
C
      DO J = 2,JROWSM1
       DO I = 1,ICOLS
        IF ( OPENSEA(I,J) ) THEN
         ZMASK(I,J) = 1.0
        ELSE
         ZMASK(I,J) = 0.0
        ENDIF
       ENDDO
      ENDDO
      DO I = 1,ICOLS
       ZMASK(I,1)     = 0
       ZMASK(I,JROWS) = 0
      ENDDO
C
C CALCULATE GRID SPACINGS IN RADIANS
C
C
      DYT   = DELTA_LAT * A / RECIP_PI_OVER_180
      DXT   = DELTA_LONG * A / RECIP_PI_OVER_180
      DYTR  = 1. / DYT
      DXT2R = .5/DXT
C
C CALCULATE COEFFICIENTS
C
      DO J = 2,JROWSM1
       ROWCOEF = 4.0 * AHDT * SEC_P_LATITUDE(1,J)
     &                 * SEC_P_LATITUDE(1,J) * DXT2R
C
C CHECK STABILITY AND RESET ROWCOEF IF UNSTABLE
C
       IF (ROWCOEF * DXT2R .GT. STABLT) ROWCOEF = STABLT / DXT2R
C
       COLCOEF_J   = AHDT * COS_U_LATITUDE(1,J)
     &               * DYTR * DYTR * SEC_P_LATITUDE(1,J)
       COLCOEF_JM1 = AHDT * COS_U_LATITUDE(1,J-1)
     &               * DYTR * DYTR * SEC_P_LATITUDE(1,J)
C
C
C CALCULATE DIFFUSION INCREMENTS USING DEL2
C
C      1. E-W INCREMENTS
C
       DO I = 2,ICOLS
        TEMPA(I) = DXT2R * ( SLABTEMP(I,J) - SLABTEMP(I-1,J) )
       END DO
       TEMPA (1) = DXT2R * (SLABTEMP(1,J) - SLABTEMP(ICOLS,J) )
C
C      2. ADD IN N-S INCREMENTS
C
       DO I = 2,ICOLSM1
        DIFFUS(I,J) = ROWCOEF
     &  * ( ZMASK(I+1,J) * TEMPA(I+1) - ZMASK(I-1,J) * TEMPA(I))
     &  + COLCOEF_J
     &  * ZMASK(I,J+1) * ( SLABTEMP(I,J+1) - SLABTEMP(I,J) )
     &  + COLCOEF_JM1
     &  * ZMASK(I,J-1) * ( SLABTEMP(I,J-1) - SLABTEMP(I,J) )
       END DO
C
C  CALCULATE DIFFUSION INCREMENTS AT 1ST AND LAST COLUMNS
C
       DIFFUS(1,J) = ROWCOEF
     &  * ( ZMASK(2,J) * TEMPA(2) - ZMASK(ICOLS,J) * TEMPA(1) )
     &  + COLCOEF_J
     &  * ZMASK(1,J+1) * ( SLABTEMP(1,J+1) - SLABTEMP(1,J) )
     &  + COLCOEF_JM1
     &  * ZMASK(1,J-1) * ( SLABTEMP(1,J-1) - SLABTEMP(1,J) )
C
       DIFFUS(ICOLS,J) = ROWCOEF
     & * ( ZMASK(1,J) * TEMPA(1) - ZMASK(ICOLSM1,J) * TEMPA(ICOLS))
     & + COLCOEF_J
     & * ZMASK(ICOLS,J+1) * (SLABTEMP(ICOLS,J+1)
     &                                     -SLABTEMP(ICOLS,J))
     & + COLCOEF_JM1
     & * ZMASK(ICOLS,J-1) * (SLABTEMP (ICOLS,J-1)
     &                                     - SLABTEMP(ICOLS,J) )
      END DO
C
C ADD IN DIFFUSION INCREMENTS.
C MULTIPLY BY ZMASK SO DIFFUSION ONLY ADDED AT OPEN SEA POINTS
C DO NOT RESET POLAR ROWS TO MAINTAIN CONSERVANCY
C
      DO J = 2,JROWSM1
       DO I = 1,ICOLS
        SLABTEMP(I,J) = SLABTEMP(I,J) + DIFFUS(I,J) * ZMASK(I,J)
       END DO
      END DO
      RETURN
      END
