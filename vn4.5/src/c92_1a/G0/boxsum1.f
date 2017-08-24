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
CLL  SUBROUTINE BOX_SUM
CLL
CLL  NOT SUITABLE FOR SINGLE COLUMN USE
CLL
CLL
CLL  SUITABLE FOR ROTATED GRIDS
CLL
CLL  ORIGINAL VERSION FOR CRAY Y-MP/IBM
CLL  WRITTEN 12/07/91 BY C. WILSON
CLL
CLL  CODE REVIEWED BY R.SMITH ??/??/??
CLL
CLL  VERSION NO. 2 DATED 16/09/91
CLL         COSMOS DSN MS15.CWUM.JOBS(BOXSUM2)
CLL  PROGRAMMING STANDARD: UNIFIED MODEL DOCUMENTATION PAPER NO. 4,
CLL  VERSION 1, DATED 12/09/89
! History:
! Version   Date     Comment
! -------   ----     ------- 
! 4.0      12/04/95  Imported into Unified model. D.M. Goddard          
! 4.1      12/06/96  Corrections for zonal means. D.M. Goddard
! 4.4      30/09/97  Corrections for portable model. D.M. Goddard
CLL
CLL  SYSTEM TASK:  S1 (part,extension for area mean interpolation)
CLL
CLL  PURPOSE:
CLL  Routine sums contributions from gridboxes for source data on a
CLL  regular lat-long grid to form means for gridboxes of a regular
CLL  lat-long grid specified as target.
CLL  Both grids are defined with the same pole and orientation;
CLL  the original data must be interpolated onto a rotated
CLL  grid ,if the target grid is a rotated grid, BEFORE calling this
CLL  routine.
CLL  The algorithms are general and will cope with either a finer
CLL  or coarser resolution source grid.
CLL
CLL  DOCUMENTATION:  UNIFIED MODEL DOCUMENTATION S1
CLL                  BY A.DICKINSON/C WILSON VERSION ??DATED ??/??/91
CLL
CLL
CLLEND-------------------------------------------------------------
C
C*L  ARGUMENTS:---------------------------------------------------
      SUBROUTINE BOX_SUM
     1  (SOURCE_ROW_LENGTH,SOURCE_ROWS,ROW_LENGTH,ROWS,
     2   LONG_L,COLAT_T,I_L,J_T,GLOBAL,BOXSUM,SOURCE)

      IMPLICIT NONE

      INTEGER
     *  SOURCE_ROW_LENGTH  !IN    Number of points per row (source data)
     *                     !      on rotated grid if necessary
     *, SOURCE_ROWS        !IN    Number of rows of source data
     *                     !      on rotated grid if necessary
     *, ROW_LENGTH         !IN    Number of points per row target area
     *, ROWS               !IN    Number of rows of target area
C
      INTEGER
     1 I_L(ROW_LENGTH+1)!IN Index of first source gridbox to overlap
     *                  !   with left hand side of target gridbox
     2,J_T(ROWS+1)      !IN Index of first source gridbox to overlap
     *                  !   top of target gridbox
C
CL N.B.I_L(I) is the first source gridbox to overlap LH side of target
CL            box  I of a row
CL     I_L(I+1) is the last source gridbox to overlap RH side of target
CL            box  I of a row
CL     J_T(J) is the first source gridbox to overlap top of target
CL            box on row J
CL     J_T(J+1) is the last source gridbox to overlap bottom of target
CL            box on row J
CL
CL REAL value of:-
CL     I_L(I) is also used to measure the 'longitude' of the RHS of the
CL            source gridbox
CL     J_T(J) is also used to measure the 'colatitude' of the bottom
CL            of the source gridbox
C

      REAL
     1 SOURCE(SOURCE_ROW_LENGTH,SOURCE_ROWS)!IN  source data
      REAL BOXSUM(ROW_LENGTH,ROWS) !OUT Sum of data on target grid
     2,    LONG_L(ROW_LENGTH +1)    !IN Left longitude of gridbox (in
     +                              ! units of souce gridbox EW length)
     3,    COLAT_T(ROWS +1)         !IN Colatitude of top of gridbox (in
     +                              ! units of source gridbox NS length)

      LOGICAL GLOBAL       !IN    true if global area required
C*---------------------------------------------------------------------

C*L  WORKSPACE USAGE:-------------------------------------------------
C   Define local workspace arrays:
C    1 real of length row_length
      REAL
     1     EW_SUM(ROW_LENGTH)     ! summed WE source data
     2,    EW_WEIGHT(ROW_LENGTH)  ! summed WE weights for source data
     3,    BOX_WEIGHT(ROW_LENGTH)  ! summed weights for target boxes
C
C*---------------------------------------------------------------------
C
C*L EXTERNAL SUBROUTINES CALLED---------------------------------------
C None
C*------------------------------------------------------------------

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

C    DEFINE LOCAL VARIABLES

      INTEGER I,J,I1,I2,IT,J1,J2,JT,K ! loop counters

      REAL RH_BOX

CL *********************************************************************
CL 1.0 Sum source boxes (whole and partial) contributions to target box
CL *********************************************************************

      DO 150 J=1,ROWS
        J1 = J_T(J)
        J2 = J_T(J+1)

        DO I=1,ROW_LENGTH
          BOX_WEIGHT(I)=0.0
        ENDDO

        DO 140 JT=J1,J2

CL *********************************************************************
CL 1.1 Sum  EW (whole and partial) contributions to target grid boxes
CL *********************************************************************

          DO I=1,ROW_LENGTH
            EW_SUM(I)=0.0
            EW_WEIGHT(I)=0.0
          ENDDO

          DO 120 I=1,ROW_LENGTH
            I1 = I_L(I)
            I2 = I_L(I+1)
            IF(I1.GT.I2.AND.GLOBAL) THEN
C   If grid box spans zero longitude need to split summation

              DO 101 IT=I1,SOURCE_ROW_LENGTH
                IF(IT.EQ.I1) THEN
C   Left side partial contribution
                  RH_BOX = LONG_L(I+1)
                  IF(RH_BOX.LT.LONG_L(I)) RH_BOX=RH_BOX
     &                                        +SOURCE_ROW_LENGTH
CXX               IF(NINT(SOURCE(I1,JT)).NE.IMDI) THEN
                  IF(NINT(SOURCE(I1,JT)).NE.NINT(RMDI)) THEN
CXX               IF(SOURCE(I1,JT).GE.0.0) THEN
                    EW_WEIGHT(I) =
     &                EW_WEIGHT(I) + (MIN(REAL(I1),RH_BOX) - LONG_L(I))
                    EW_SUM(I) = (MIN(REAL(I1),RH_BOX) - LONG_L(I))
     &                      *SOURCE(I1,JT) + EW_SUM(I)
                  ELSEIF(NINT(SOURCE(I1,JT)).EQ.NINT(RMDI)) THEN
C                   missing data
                  ELSE
                    write (6,*) '-ve source 1 ? I1 JT data',I1,JT,
     +              source(i1,jt)
                  ENDIF

                ELSE

C   Whole contributions
CXX               IF(NINT(SOURCE(IT,JT)).NE.IMDI) THEN
                  IF(NINT(SOURCE(IT,JT)).NE.NINT(RMDI)) THEN
CXX               IF(SOURCE(IT,JT).GE.0.0) THEN
                    EW_WEIGHT(I) = EW_WEIGHT(I)+ 1.0
                    EW_SUM(I) = EW_SUM(I) + SOURCE(IT,JT)
                  ELSEIF(NINT(SOURCE(IT,JT)).EQ.NINT(RMDI)) THEN
C                   missing data
                  ELSE
                    write (6,*) '-ve source 2 ? IT JT data',IT,JT,
     +              source(it,jt)
                  ENDIF

                ENDIF

101           CONTINUE

              DO 102 IT=1,I2
                IF(IT.EQ.I2) THEN
C   Right side partial contribution
CXX               IF(NINT(SOURCE(I2,JT)).NE.IMDI) THEN
                  IF(NINT(SOURCE(I2,JT)).NE.NINT(RMDI)) THEN
CXX               IF(SOURCE(I2,JT).GE.0.0) THEN
                    EW_WEIGHT(I) = EW_WEIGHT(I)+(LONG_L(I+1)-(I2-1))
                    EW_SUM(I)=EW_SUM(I) +
     &                               (LONG_L(I+1)-(I2-1))*SOURCE(I2,JT)
                  ELSEIF(NINT(SOURCE(I2,JT)).EQ.NINT(RMDI)) THEN
C                   missing data
                  ELSE
                    write (6,*) '-ve source 3 ? I2 JT data',I2,JT,
     +              source(i2,jt)
                  ENDIF
                ELSE

C   Whole contributions
CXX               IF(NINT(SOURCE(IT,JT)).NE.IMDI) THEN
                  IF(NINT(SOURCE(IT,JT)).NE.NINT(RMDI)) THEN
CXX               IF(SOURCE(IT,JT).GE.0.0) THEN
                    EW_WEIGHT(I) = EW_WEIGHT(I)+ 1.0
                    EW_SUM(I) = EW_SUM(I) + SOURCE(IT,JT)
                  ELSEIF(NINT(SOURCE(IT,JT)).EQ.NINT(RMDI)) THEN
C                   missing data
                  ELSE
                    write (6,*) '-ve source 4 ? IT JT data',IT,JT,
     +              source(it,jt)
                  ENDIF

                ENDIF

102           CONTINUE

          ELSE IF(I1.LT.I2)THEN ! no zero meridian crossing
            DO 110 IT=I1,I2
              IF(IT.EQ.I1) THEN
C   Left side partial contribution
CXX             IF(NINT(SOURCE(I1,JT)).NE.IMDI) THEN
                IF(NINT(SOURCE(I1,JT)).NE.NINT(RMDI)) THEN
CXX             IF(SOURCE(I1,JT).GE.0.0) THEN
                  EW_WEIGHT(I) = EW_WEIGHT(I) +
     &                      (MIN(REAL(I1),LONG_L(I+1)) - LONG_L(I))
                  EW_SUM(I) = (MIN(REAL(I1),LONG_L(I+1)) - LONG_L(I))
     &                      *SOURCE(I1,JT) + EW_SUM(I)
                  ELSEIF(NINT(SOURCE(I1,JT)).EQ.NINT(RMDI)) THEN
C                   missing data
                  ELSE
                    write (6,*) '-ve source 5 ? I1 JT data',I1,JT,
     +              source(i1,jt)
                ENDIF

              ELSE IF(IT.EQ.I2) THEN
C   Right side partial contribution
CXX             IF(NINT(SOURCE(I2,JT)).NE.IMDI) THEN
                IF(NINT(SOURCE(I2,JT)).NE.NINT(RMDI)) THEN
CXX             IF(SOURCE(I2,JT).GE.0.0) THEN
                  EW_WEIGHT(I) = EW_WEIGHT(I)+ (LONG_L(I+1)-(I2-1))
                  EW_SUM(I)=EW_SUM(I)+(LONG_L(I+1)-(I2-1))*SOURCE(I2,JT)
                  ELSEIF(NINT(SOURCE(I2,JT)).EQ.NINT(RMDI)) THEN
C                   missing data
                  ELSE
                    write (6,*) '-ve source 6 ? I2 JT data',I2,JT,
     +              source(i2,jt)
                ENDIF

              ELSE

C   Whole contributions

CXX             IF(NINT(SOURCE(IT,JT)).NE.IMDI) THEN
                IF(NINT(SOURCE(IT,JT)).NE.NINT(RMDI)) THEN
CXX             IF(SOURCE(IT,JT).GE.0.0) THEN
                  EW_WEIGHT(I) = EW_WEIGHT(I)+ 1.0
                  EW_SUM(I) = EW_SUM(I) + SOURCE(IT,JT)
                  ELSEIF(NINT(SOURCE(IT,JT)).EQ.NINT(RMDI)) THEN
C                   missing data
                  ELSE
                    write (6,*) '-ve source 7 ? IT JT data',IT,JT,
     +              source(it,jt)
                ENDIF

              ENDIF
110           CONTINUE

          ELSE 
!Zonal mean no need to average in EW direction
            DO K=1,ROW_LENGTH
              EW_WEIGHT(K)=1.0
              EW_SUM(K)=SOURCE(1,JT)
            END DO

          

          ENDIF

120       CONTINUE

CL *********************************************************************
CL 1.3 Add summed EW  box contributions into rows J  target grid boxes
CL *********************************************************************

          IF(JT.EQ.J1) THEN
C   Top row
            DO 130 I=1,ROW_LENGTH
              BOXSUM(I,J) = (MIN(REAL(J1),COLAT_T(J+1)) - COLAT_T(J))
     &                       *EW_SUM(I)
              BOX_WEIGHT(I) = (MIN(REAL(J1),COLAT_T(J+1)) - COLAT_T(J))
     &                       *EW_WEIGHT(I) + BOX_WEIGHT(I)
130         CONTINUE

          ELSE IF(JT.EQ.J2) THEN
C   Bottom of row J
            DO 131 I=1,ROW_LENGTH
              BOXSUM(I,J) = BOXSUM(I,J) +
     1         (1-(J2 - COLAT_T(J+1)))*EW_SUM(I)
              BOX_WEIGHT(I) = BOX_WEIGHT(I) +
     &         (1-(J2 - COLAT_T(J+1)))*EW_WEIGHT(I)
131         CONTINUE

          ELSE
C   Whole contributions to row J
            DO 132 I=1,ROW_LENGTH
              BOXSUM(I,J) = BOXSUM(I,J) + EW_SUM(I)
              BOX_WEIGHT(I) = BOX_WEIGHT(I) + EW_WEIGHT(I)
132         CONTINUE

          ENDIF

140     CONTINUE

          DO 142 I=1,ROW_LENGTH
            IF(BOX_WEIGHT(I).NE.0.0) THEN
              BOXSUM(I,J) = BOXSUM(I,J) / BOX_WEIGHT(I)
            ELSE
              BOXSUM(I,J) = RMDI
            ENDIF
142       CONTINUE

150   CONTINUE

      RETURN
      END
