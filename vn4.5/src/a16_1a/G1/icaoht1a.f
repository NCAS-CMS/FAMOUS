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
CLL  SUBROUTINE ICAO_HT-----------------------------------------
CLL
CLL  PURPOSE:   Performs an interpolation from pressure levels to the
CLL             I.C.A.O standard atmosphere heights.
CLL             e.g. tropopause heights or max wind heights in
CLL             thousands of feet
CLL  Tested under compiler CFT77
CLL  Tested under OS version 5.1
CLL
CLL  Author J.T.Heming
CLL  Code version 1.0         Date 04/12/90
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  date
CLL
CLL  Logical components covered D4
CLL
CLL  Project TASK: D4
CLL
CLL  Programming standard: U M DOC  Paper NO. 4,
CLL
CLL  External documentation
CLL
CLLEND------------------------------------------------------------------
C
C*L  ARGUMENTS:---------------------------------------------------------
      SUBROUTINE ICAO_HT(
     & PRESS_IN,POINTS,ICAO_HEIGHT)
C-----------------------------------------------------------------------
      IMPLICIT NONE
C-----------------------------------------------------------------------
      INTEGER
     * POINTS           ! IN No. of points to be processed
C-----------------------------------------------------------------------
      REAL
     * PRESS_IN(POINTS) ! IN Pressures to be converted to ICAO heights
C-----------------------------------------------------------------------
      REAL
     * ICAO_HEIGHT(POINTS)  ! OUT ICAO height in thousands of feet
C*----------------------------------------------------------------------
! History:
! Version  Date  Comment
!  3.4   18/5/94 Add PP missing data indicator. J F Thomson
C*L------------------COMDECK C_MDI-------------------------------------
      REAL    RMDI_PP   ! PP missing data indicator (-1.0E+30)
      REAL    RMDI_OLD  ! Old real missing data indicator (-32768.0)
      REAL    RMDI      ! New real missing data indicator (-2**30)
      INTEGER IMDI      ! Integer missing data indicator

      PARAMETER(
     & RMDI_PP  =-1.0E+30,
     & RMDI_OLD =-32768.0,
     & RMDI =-32768.0*32768.0,
     & IMDI =-32768)
C*----------------------------------------------------------------------
C
C*L WORKSPACE USAGE-----------------------------------------------------
C   DEFINE LOCAL WORKSPACE ARRAYS
      REAL
     * PHST(15)
     *,HST(15)
     *,PRESSURE(POINTS) ! Local pressure array
C*----------------------------------------------------------------------
C
C*L EXTERNAL SUBROUTINES CALLED-----------------------------------------
C   NONE
C*----------------------------------------------------------------------
C
C-----------------------------------------------------------------------
C*L   DEFINE LOCAL VARIABLES
C-----------------------------------------------------------------------
      INTEGER
     * I,K                ! LOOP COUNTERS
C-----------------------------------------------------------------------
      REAL
     * TERM1
     *,TERM2
C-----------------------------------------------------------------------
      LOGICAL
     * FOUND(POINTS)
C*
C-----------------------------------------------------------------------
C*L   DATA FOR LOCAL ARRAYS
C-----------------------------------------------------------------------
      DATA PHST/101325.,95000.,85000.,70000.,50000.,40000.,30000.,
     + 25000.,20000.,15000.,10000.,7000.,5000.,3000.,999./
      DATA HST/0.,2000.,5000.,10000.,18000.,24000.,30000.,34000.,
     + 39000.,45000.,53000.,61000.,68000.,78000.,99000./
C*----------------------------------------------------------------------
C     LOOP OVER ALL POINTS
C-----------------------------------------------------------------------
      DO I=1,POINTS
        PRESSURE(I)=PRESS_IN(I)
C-----------------------------------------------------------------------
CL 1. If the pressure is less than 10mb it is assumed to be 10mb
CL    If pressure is greater than 1013.25mb it's assumed to be 1013.25mb
C-----------------------------------------------------------------------
        IF(PRESSURE(I).LE.1000.AND.PRESSURE(I).GE.0.) PRESSURE(I)=1000.
        IF(PRESSURE(I).GT.101325.) PRESSURE(I)=101325.
        FOUND(I)=.FALSE.
      ENDDO
C-----------------------------------------------------------------------
CL 2. FIND THE FIRST VALUE OF PHST THAT IS LESS THAN THE PRESSURE
C-----------------------------------------------------------------------
      DO 1 K=1,14
        DO I=1,POINTS
C-----------------------------------------------------------------------
CL 3. If the pressure is set to missing data indicator, set ICAO height
CL    to missing data indicator
C-----------------------------------------------------------------------
          IF(PRESSURE(I).EQ.RMDI.AND.(.NOT.FOUND(I)))THEN
            ICAO_HEIGHT(I)=RMDI
            FOUND(I)=.TRUE.
          ELSEIF((PHST(K+1).LT.PRESSURE(I)).AND.(.NOT.FOUND(I))) THEN
C-----------------------------------------------------------------------
CL 4. Calculate the ICAO height
C-----------------------------------------------------------------------
            TERM1=(3*PHST(K)-PRESSURE(I))*(PRESSURE(I)-PHST(K))
            TERM2=(3*PHST(K)-PHST(K+1))*(PHST(K+1)-PHST(K))
            ICAO_HEIGHT(I)=((HST(K+1)-HST(K))*TERM1/TERM2+HST(K))*0.001
            FOUND(I)=.TRUE.
C-----------------------------------------------------------------------
          ENDIF
        ENDDO
 1    CONTINUE
C-----------------------------------------------------------------------
      RETURN
      END
