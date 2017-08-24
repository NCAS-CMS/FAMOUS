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
CLL  SUBROUTINE H_INT_CO----------------------------------------------
CLL
CLL  Purpose:  Calculates bi-linear horizontal interpolation
CLL            coefficients and gather indices for interpolating
CLL            between generalised latitude-longitude grids (eg
CLL            global, regional or rotated lat-lon grid) in which the
CLL            gridlength may vary with latitude and/or longitude. The
CLL            interpolation is carried out by subroutine
CLL            H_INT. Gather indices point to bottom left hand
CLL            corner and bottom right hand corner of each grid box on
CLL            source grid enclosing a target point. Two indices are
CLL            needed to cater for east-west (lambda direction) cyclic
CLL            boundaries when the source data is global. If a target po
CLL            falls outside the domain of the source data, one sided
CLL            differencing is used. The source latitude coordinates
CLL            must be supplied in decreasing order. The source long-
CLL            itude coordinates must be supplied in increasing order,
CLL            starting at any value, but not wrapping round. The
CLL            target points may be specified in any order.
CLL
CLL A.Dickinson <- programmer of some or all of previous code or changes
CLL J.Gregory   <- programmer of some or all of previous code or changes
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  Date
CLL   3.4   24/06/94  New checks to keep within array bounds for
CLL                   non_cyclic cases. D. Robinson
!     4.5   29/07/98  Optimisation changes for T3E
!                     Author D.M. Goddard
CLL
CLL  Programming standard:
CLL           Unified Model Documentation Paper No 3
CLL           Version No 1 15/1/90
CLL
CLL  System component: S121,S122
CLL
CLL  System task: S1
CLL
CLL  Documentation:
CLL            The interpolation formulae are described in
CLL            unified model on-line documentation paper S1.
CLL
CLL  -----------------------------------------------------------------
C
C*L  Arguments:-------------------------------------------------------
      SUBROUTINE H_INT_CO
     *(INDEX_B_L,INDEX_B_R,WEIGHT_T_R,WEIGHT_B_R,WEIGHT_T_L,WEIGHT_B_L
     *,LAMBDA_SRCE,PHI_SRCE,LAMBDA_TARG,PHI_TARG
     *,POINTS_LAMBDA_SRCE,POINTS_PHI_SRCE,POINTS,CYCLIC)

      IMPLICIT NONE

      INTEGER
     * POINTS_LAMBDA_SRCE !IN Number of lambda points on source grid
     *,POINTS_PHI_SRCE    !IN Number of phi points on source grid
     *,POINTS             !IN Total number of points on target grid
     *,INDEX_B_L(POINTS)  !OUT Index of bottom lefthand corner
     *                    !    of source gridbox
     *,INDEX_B_R(POINTS)  !OUT Index of bottom righthand corner
     *                    !    of source gridbox

      REAL
     * LAMBDA_TARG(POINTS) !IN Lambda coords of target grid in degrees
     *                     !   using same rotation as source grid
     *,PHI_TARG(POINTS)    !IN Phi coords of target grid in degrees
     *                     !   using same rotation as source grid
     *,WEIGHT_T_R(POINTS)  !OUT Weight applied to value at top right
     *                     !    hand corner of source gridbox
     *,WEIGHT_B_L(POINTS)  !OUT Weight applied to value at bottom left
     *                     !    hand corner of source gridbox
     *,WEIGHT_B_R(POINTS)  !OUT Weight applied to value at bottom right
     *                     !    hand corner of source gridbox
     *,WEIGHT_T_L(POINTS)  !OUT Weight applied to value at top left
     *                     !    hand corner of source gridbox
     *,LAMBDA_SRCE(POINTS_LAMBDA_SRCE) !IN Lambda coords of source grid
     *                     !    in degrees
     *,PHI_SRCE(POINTS_PHI_SRCE) !IN Phi coords of target grid in degree

      LOGICAL
     * CYCLIC              !IN =T, then source data is cyclic
     *                     !   =F, then source data is non-cyclic

C Local arrays:---------------------------------------------------------
      REAL
     * T_LAMBDA(POINTS) !Local value of target longitude

      INTEGER
     * IXP1(POINTS)    !Longitudinal index plus 1
     *,IX(POINTS)      !Longitudinal index
     *,IY(POINTS)      !Latitudinal index
C External subroutines called:------------------------------------------
C None
C*----------------------------------------------------------------------
C*L  Local variables:---------------------------------------------------
      REAL
     * A              !Longitudinal weight
     *,B              !Latitudinal weight

      INTEGER
     * I,J            !Loop indices
C ----------------------------------------------------------------------

CL    1. Initialise arrays

      IF(CYCLIC)THEN

C     1.1 Cyclic case

        DO I=1,POINTS

C  Scale target longitude so that it falls between LAMBDA_SRCE(1)
C  and LAMBDA_SRCE(1)+360
          T_LAMBDA(I)=MOD(LAMBDA_TARG(I)-LAMBDA_SRCE(1)+720.,360.)
     &    +LAMBDA_SRCE(1)

C Initialise longitudinal & latitudinal indices
          IX(I)=0
          IY(I)=1
        ENDDO

      ELSE

C     1.2 Non cyclic case

        DO I=1,POINTS

C  Assign target longitude to local array for use below
          T_LAMBDA(I)=LAMBDA_TARG(I)

C Initialise longitudinal & latitudinal indices
          IX(I)=1
          IY(I)=1

        ENDDO

      ENDIF

CL 2. Calculate lat and lon index of bottom left hand corner of
CL    source grid box enclosing each target point.

C Longitude
       DO I=1,POINTS
         DO J=1,POINTS_LAMBDA_SRCE
           IF(LAMBDA_SRCE(J).LE.T_LAMBDA(I))THEN
             IX(I)=J
           ELSE
             GOTO 120
           ENDIF
         ENDDO
 120   CONTINUE
       ENDDO

C Latitude
       DO I=1,POINTS
         DO J=1,POINTS_PHI_SRCE
           IF(PHI_SRCE(J).GE.PHI_TARG(I))THEN
             IY(I)=J
           ELSE
             GOTO 140
           ENDIF
         ENDDO
 140   CONTINUE
       ENDDO



CL  3. Correct 1-D indices for wrap around etc and then calculate
CL     2-D indices of bottom left and bottom right hand corner
CL     of each grid box.

      IF(CYCLIC)THEN
C     3.1 Cyclic case

      DO I=1,POINTS

C Set index for cyclic wrap around
        IF(IX(I).LT.1)THEN
          IX(I)=POINTS_LAMBDA_SRCE
          T_LAMBDA(I)=T_LAMBDA(I)+360.
        ENDIF

C Set index for one sided difference if target point to north or
C south of source area.
        IY(I)=MAX(IY(I),1)
        IY(I)=MIN(IY(I),POINTS_PHI_SRCE-1)

C 2-D indices
        INDEX_B_L(I)=IX(I)+IY(I)*POINTS_LAMBDA_SRCE
        INDEX_B_R(I)=INDEX_B_L(I)+1

C Correct for cyclic boundaries if target point outside source grid.

        IXP1(I)=IX(I)+1
        IF(IX(I).EQ.POINTS_LAMBDA_SRCE)THEN
          INDEX_B_R(I)=INDEX_B_R(I)-POINTS_LAMBDA_SRCE
          IXP1(I)=IXP1(I)-POINTS_LAMBDA_SRCE
        ENDIF

      ENDDO

      ELSE

C     3.2 Non cyclic case

      DO I=1,POINTS

C Set index for one sided difference if outside source area
        IX(I)=MAX(IX(I),1)
        IX(I)=MIN(IX(I),POINTS_LAMBDA_SRCE-1)
        IF (IX(I).LT.1) THEN ! IX(I) < 1 if POINTS_LAMBDA_SRCE = 1
          IX(I)=1
        ENDIF

        IXP1(I)=IX(I)+1
        IXP1(I)=MIN(IXP1(I),POINTS_LAMBDA_SRCE)

C Set index for one sided difference if outside source area
        IY(I)=MAX(IY(I),1)
        IY(I)=MIN(IY(I),POINTS_PHI_SRCE-1)
        IF (IY(I).LT.1) THEN ! IY(I) < 1 if POINTS_PHI_SRCE = 1
          IY(I)=1
        ENDIF


C 2-D indices
        INDEX_B_L(I)=IX(I)+IY(I)*POINTS_LAMBDA_SRCE
        INDEX_B_R(I)=INDEX_B_L(I)+1

      ENDDO

      ENDIF

CL 4. Compute interpolation weights

      DO I=1,POINTS

C Calculate basic weights (equation 2.2)
        A=(AMOD(360.+LAMBDA_SRCE(IXP1(I))-LAMBDA_SRCE(IX(I)),360.))
        IF(A.NE.0.)THEN
          A=(T_LAMBDA(I)-LAMBDA_SRCE(IX(I)))/A
        ELSE
          A=0.
        ENDIF

        B=(PHI_TARG(I)-PHI_SRCE(IY(I)+1))/
     *  (PHI_SRCE(IY(I))-PHI_SRCE(IY(I)+1))

C Calculate bi-linear interpolation weights as per equation (2.1)

        WEIGHT_T_R(I)=A*B
        WEIGHT_B_L(I)=(1.-A)*(1.-B)
        WEIGHT_T_L(I)=(1.-A)*B
        WEIGHT_B_R(I)=A*(1.-B)

      ENDDO

      RETURN
      END
