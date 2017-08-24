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
CLL  SUBROUTINE V_INT-------------------------------------------------
CLL
CLL  Purpose:  Performs vertical interpolation from one arbitrary set
CLL            of pressure levels to another. The technique used is
CLL            linear interpolation in log(p). When interpolating
CLL            wind components there is an option (controlled by
CLL            MAX_WIND) for including data from max wind modelling.
CLL
CLL  Written by A. Dickinson
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  date
CLL   3.1  23/02/93    DO ALL directive inserted before loop 240
CLL                    Author: A. Dickinson    Reviewer: F. Rawlins
CLL
CLL   4.2  01/07/96    Revised for CRAY T3E. Vector gather replaced
CLL                    by algorithm in which each term in the 
CLL                    interpolation formula is first collected into
CLL                    a separate array and the interpolation
CLL                    calculation carried out after the loop over level
CLL                    New arguments START and END introduced to
CLL                    facilitate the removal of duplicate calculations
CLL                    when using domain decomposition in MPP mode.     
CLL                    Author: A. Dickinson    Reviewer: F. Rawlins     
!LL   4.5  14/04/98    Use assumption that neighbouring points are
!LL                    likely to be on or near same level. Jump out
!LL                    of loop-over-levels once level found. Results
!LL                    in a 40 percent speedup on 19 levels for
!LL                    non-vector machines. S.D.Mullerworth
CLL
CLL                                                                     
CLL Programming standard :
CLL
CLL Logical components covered : S111
CLL
CLL Project task :
CLL
CLL  Documentation: The interpolation formulae are described in
CLL                 unified model on-line documentation paper S1.
CLL
CLLEND -----------------------------------------------------------------
C
C*L  ARGUMENTS:-------------------------------------------------------
      SUBROUTINE V_INT(P_IN,P_OUT,DATA_IN,DATA_OUT,POINTS,LEVELS
     *               ,DATA_MAXW,P_MAXW,MAX_WIND
     &               ,START,END)

      IMPLICIT NONE

      INTEGER
     * POINTS ! Number of points to be processed.
     *,LEVELS ! Number of levels in source data.
     *,START  ! Start position at each level
     *,END    ! Last point to be processed at each level            

      REAL
     * P_IN(POINTS,LEVELS)   !IN 3-D field of pressures at which
     *                       ! source data is stored.
     *,P_OUT(POINTS)         !IN Array of pressure values to be
     *                       ! interpolated to.
     *,DATA_IN(POINTS,LEVELS)!IN Source data as 3-D field.
     *,DATA_OUT(POINTS)      !OUT Result of interpolation.
     *,DATA_MAXW(POINTS)     !IN Max wind data.
     *,P_MAXW(POINTS)        !IN Pressure of max wind data.

      LOGICAL
     * MAX_WIND !IN Switch to include max winds if required.

C Workspace usage:-----------------------------------------------------
      REAL                                             
     * P1(POINTS)            ! Upper input pressure \
     *,P2(POINTS)            ! Lower input pressure  \ Used in interp-
     *,D1(POINTS)            ! Upper input data      / olation formula
     *,D2(POINTS)            ! Lower input data     /              
C External subroutines called:-----------------------------------------
C None                                                  
C*---------------------------------------------------------------------
C Define local variables:----------------------------------------------
      INTEGER I,J
     &  ,LAST                   ! Stores level of preceding point
      REAL ALPHA
C----------------------------------------------------------------------

! Initialise LAST to any value between 1 and LEVELS
      LAST=2
      DO I=START,END

! Start from same level as last point. First check whether this point
! is above or below, then continue search in appropriate direction
        IF(P_OUT(I).GE.P_IN(I,LAST))THEN

! These next two loops exit immediately once level found.
! GOTO cuts out needless looping once level is found, reducing the
! cost of the routine by about 40 percent for 19 level runs.
          DO J=LAST,2,-1
            IF(P_OUT(I).LT.P_IN(I,J-1))THEN
              GOTO 240
            ENDIF
          ENDDO
        ELSE
          DO J=LAST+1,LEVELS
            IF(P_OUT(I).GE.P_IN(I,J))THEN
              GOTO 240
            ENDIF
          ENDDO
        ENDIF
 240    CONTINUE

! At this point, J is:
!    1         for below bottom level.
!    LEVELS+1  for above top level
!    Otherwise J is the level just above the point

        IF (J.GT.1.AND.J.LE.LEVELS)THEN
! Between top and bottom level
          P1(I)=P_IN(I,J)
          P2(I)=P_IN(I,J-1)
          D1(I)=DATA_IN(I,J)
          D2(I)=DATA_IN(I,J-1)
          LAST=J
        ELSE
! Special case; above top or below bottom.
! Set output field to top/bottom-most input field
          IF(J.EQ.LEVELS+1)J=LEVELS
          P1(I)=P_OUT(I)
          P2(I)=1.0
          D1(I)=DATA_IN(I,J)
          D2(I)=0.0
          LAST=J
        ENDIF
      ENDDO ! DO I=START,END

! If there is an extra level of winds from max wind modelling, include
! these in the interpolation. Repeat the level-finding logic because
! there are no calls with MAX_WIND=.TRUE. in UM so do not want to slow
! down the above loop by including the MAX_WIND test in the above.

      IF (MAX_WIND)THEN
        DO I=START,END

! If max wind level between current levels, redo interpolation
! incorporating max wind info.

! Start from same level as last point. First check whether this point
! is above or below, then check all levels above/below in turn
          IF(P_OUT(I).GE.P_IN(I,LAST))THEN
! Below LAST level.
! These loops exit immediately once level found.
! GOTO cuts out needless looping once level is found, reducing the
! cost of the routine by about 40 percent for 19 level runs.
            DO J=LAST,2,-1
              IF(P_OUT(I).LT.P_IN(I,J-1))THEN
                GOTO 340
              ENDIF
            ENDDO
          ELSE
            DO J=LAST+1,LEVELS
              IF(P_OUT(I).GE.P_IN(I,J))THEN
                GOTO 340
              ENDIF
            ENDDO
          ENDIF
 340      CONTINUE

          IF(J.GT.1.AND.J.LE.LEVELS)THEN
            IF(P_MAXW(I).LT.P_IN(I,J-1).AND.P_MAXW(I).GE.P_IN(I,J))THEN

              IF(P_OUT(I).LT.P_MAXW(I))THEN

! (i)  p(maxwind) > p(out) >= p(j)

                P2(I)=P_MAXW(I)
                D2(I)=DATA_MAXW(I)

              ELSE

! (ii) p(j-1) > p(out) >= p(maxwind)

                P1(I)=P_MAXW(I)
                D1(I)=DATA_MAXW(I)

              ENDIF
            ENDIF
          ENDIF

        ENDDO                   ! DO I=START,END

      ENDIF
 
CL 3. Compute equation (3.3)
   
C Compute alpha, the interpolation weight given by equation (3.4)
      DO I=START,END
          ALPHA=ALOG(P_OUT(I)/P2(I))                              
     *         /ALOG(P1(I)/P2(I))    
C Then apply equation (3.3)                                             
          DATA_OUT(I)=ALPHA*D1(I)+(1.-ALPHA)*D2(I) 
      ENDDO     

      RETURN
      END
