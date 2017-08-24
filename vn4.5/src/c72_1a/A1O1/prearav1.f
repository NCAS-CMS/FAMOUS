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
      SUBROUTINE PRE_AREAVER(GAPS_LAMBDA_SRCE,LAMBDA_SRCE
     &,GAPS_PHI_SRCE,PHI_SRCE,CYCLIC_SRCE,LROW_SRCE,WANT,MASK_SRCE
     &,GAPS_LAMBDA_TARG,LAMBDA_TARG,GAPS_PHI_TARG,PHI_TARG
     &,CYCLIC_TARG,SPHERICAL
     &,MAXL,COUNT_TARG,BASE_TARG,INDEX_SRCE,WEIGHT,ICODE,CMESSAGE)
CLL
CLL   Subroutine PRE_AREAVER
CLL
CLL   Calculate weights for area-averaging data on the source grid to
CLL   data on the target grid.
CLL
CLL  J.Gregory  <- programmer of some or all of previous code or changes
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  Date
CLL   3.2    13/07/93 Changed CHARACTER*(*) to CHARACTER*(80) for
CLL                   portability.  Author Tracey Smith.
CLL   3.3   3.9.93    Correct out-of-bounds errors for IXL and IYT
CLL   4.1   22.5.96   J.M.Gregory  Add argument SPHERICAL to produce
CLL                   correct weights for a spherical surface.
CLL
CLL   Standard, paper 4 version 4 (14.12.90)
CLL
CL
CL    The source grid and target grid are each specified by a lambda set
CL    and a phi set of coordinates delimiting the boxes. These sets are
CL    supplied in 1D arrays aa_bb for coordinate aa=LAMBDA,PHI on grid
CL    bb=SRCE,TARG.  The number of gaps is specified by GAPS_aa_bb,
CL    which is equal to the number of lines IF (CYCLIC_bb) and aa is
CL    LAMBDA, otherwise one less. (By "gap" we mean the interval
CL    spanning a box in one coordinate only. The total number of boxes,
CL    or grid points, is the product of the numbers of gaps in the two
CL    coordinates.) Whether the axes are cyclic is not known until
CL    run-time, so the dimensions of the arrays LAMBDA_bb are not known
CL    at compile-time, and they are dimensioned assumed-size. There are
CL    no restrictions on the base meridian of lambda, and it does not
CL    have to be the same for source and target grids. The lambda
CL    coordinates should be in increasing order (running from left to
CL    right), the phi decreasing (top to bottom). The coordinates must
CL    be given in degrees for cyclic axes, because a range of 360 is
CL    assumed, or IF (SPHERICAL), when trigonometric functions are
CL    used. IF (SPHERICAL), the weights computed are for a spherical
CL    surface, assuming that LAMBDA is longitude and PHI latitude.
CL    Otherwise, LAMBDA and PHI are treated as Cartesian axes on a
CL    plane.
CL
CL    The array MASK_SRCE is the land/sea mask for the source grid. The
CL    logical value indicating points to be used should be supplied in
CL    WANT. The first dimension of MASK_SRCE should be supplied in
CL    LROW_SRCE. Specifying this separately allows for the possibility
CL    of extra rows and columns in MASK_SRCE which are to be ignored.
CL
CL    The arrays COUNT_TARG and BASE_TARG should be dimensioned in the
CL    calling program to the number of boxes in the target array.
CL
CL    The arrays INDEX_SRCE and WEIGHT are returned in the form which
CL    the area-averaging routine DO_AREAVER expects. They are continuous
CL    lists comprising consecutive groups of entries. There is a group
CL    for each target point, for which the number of entries is spec-
CL    ified by COUNT_TARG, and the groups appear in the normal order of
CL    grid points. The size required for INDEX_SRCE and WEIGHT depends
CL    on how many source boxes go into each target box, on average, and
CL    is not known at compile-time. The maximum that could be needed is
CL    (GAPS_LAMBDA_SRCE+GAPS_LAMBDA_TARG)*(GAPS_PHI_SRCE+GAPS_PHI_TARG)
CL    and the size to which the arrays are actually dimensioned should
CL    be supplied in MAXL. The size used is returned in MAXL. It is the
CL    responsibility of the calling routine to provide enough space.
CL
      IMPLICIT NONE
C*L
      INTEGER
     & GAPS_LAMBDA_SRCE        !IN number of lambda gaps in source grid
     &,GAPS_PHI_SRCE           !IN number of phi gaps in source grid
     &,GAPS_LAMBDA_TARG        !IN number of lambda gaps in target grid
     &,GAPS_PHI_TARG           !IN number of phi gaps in target grid
     &,LROW_SRCE               !IN first dimension of MASK_SRCE
     &,MAXL                    !INOUT maximum entries in output lists
     &,COUNT_TARG(GAPS_LAMBDA_TARG,GAPS_PHI_TARG)
C                              !OUT no. of source boxes in target box
     &,BASE_TARG(GAPS_LAMBDA_TARG,GAPS_PHI_TARG)
C                              !OUT first index in list for target box
     &,INDEX_SRCE(MAXL)        !OUT list of source box indices
     &,ICODE                   !OUT return code
      LOGICAL
     & CYCLIC_SRCE             !IN source grid is cyclic
     &,CYCLIC_TARG             !IN target grid is cyclic
     &,WANT                    !IN indicator of wanted points in mask
     &,MASK_SRCE(LROW_SRCE,*)  !IN land/sea mask for source grid
     &,SPHERICAL               !IN calculate weights for a sphere
      REAL
     & LAMBDA_SRCE(*)          !IN source lambda line coordinates
     &,PHI_SRCE(*)             !IN source phi line coordinates
     &,LAMBDA_TARG(*)          !IN target lambda line coordinates
     &,PHI_TARG(*)             !IN target phi line coordinates
     &,WEIGHT(MAXL)            !OUT list of weights for source boxes
      CHARACTER*(80) CMESSAGE  !OUT error message
C*
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
      INTEGER
     & LINES_LAMBDA_SRCE       ! number of source lambda lines
     &,LINES_PHI_SRCE          ! number of source phi lines
     &,LINES_LAMBDA_TARG       ! number of target lambda lines
     &,LINES_PHI_TARG          ! number of target phi lines
     &,COUNT_LAMBDA(GAPS_LAMBDA_TARG)
C                              ! number of source lambda gaps per target
     &,COUNT_PHI(GAPS_PHI_TARG)
C                              ! number of source phi gaps per target
     &,INDEX_LAMBDA(GAPS_LAMBDA_SRCE+GAPS_LAMBDA_TARG)
C                              ! source lambda gap indices
     &,INDEX_PHI(GAPS_PHI_SRCE+GAPS_PHI_TARG)
C                              ! source phi gap indices
     &,IX1,IY1,IX2,IY2         ! working SRCE/TARG LAMBDA/PHI indices
     &,IX1N,IX1W               ! working indices
     &,IXL(GAPS_LAMBDA_TARG+1) ! source line on the left of target line
     &,IX2N                    ! target gap to the right of IX2
     &,IYT(GAPS_PHI_TARG+1)    ! source line above target line
     &,IXP,IYP,IP              ! pointers into lists
     &,IX,IY,I                 ! loop indices
      REAL
     & BASLAM                  ! minimum lambda for TEMP coordinates
     &,BTARG                   ! width of target gap
     &,BLO,BHI                 ! limits of gap
     &,TEMP_SRCE(GAPS_LAMBDA_SRCE+1)
C                              ! adjusted version of LAMBDA_SRCE
     &,TEMP_TARG(GAPS_LAMBDA_TARG+1)
C                              ! adjusted version of LAMBDA_TARG
     &,FRAC_LAMBDA(GAPS_LAMBDA_SRCE+GAPS_LAMBDA_TARG)
C                              ! fractions of target lambda gaps
     &,FRAC_PHI(GAPS_PHI_SRCE+GAPS_PHI_TARG)
C                              ! fractions of target phi gaps
     &,SUM                     ! sum of weights
C
CL    1  Set up the lambda coordinates to make them easier to handle.
C
CL    1.1  Produce in TEMP_SRCE a monotonically increasing set of angles
CL    equivalent to LAMBDA_SRCE i.e. equal under modulo 360.
C
      IF (CYCLIC_SRCE) THEN
        LINES_LAMBDA_SRCE=GAPS_LAMBDA_SRCE
      ELSE
        LINES_LAMBDA_SRCE=GAPS_LAMBDA_SRCE+1
      ENDIF
      BASLAM=LAMBDA_SRCE(1)
      DO 10 IX1=1,LINES_LAMBDA_SRCE
        IF (LAMBDA_SRCE(IX1).LT.BASLAM) THEN
          TEMP_SRCE(IX1)=LAMBDA_SRCE(IX1)+360.
        ELSE
          TEMP_SRCE(IX1)=LAMBDA_SRCE(IX1)
        ENDIF
   10 CONTINUE ! Labelled DO for the sake of fpp
C
CL    1.2  Produce in TEMP_TARG a set of angles equivalent to
CL    LAMBDA_TARG i.e. equal under modulo 360, but all in the range
CL    BASLAM to BASLAM+360, where BASLAM=min(TEMP_LAMBDA).
C
      IF (CYCLIC_TARG) THEN
        LINES_LAMBDA_TARG=GAPS_LAMBDA_TARG
      ELSE
        LINES_LAMBDA_TARG=GAPS_LAMBDA_TARG+1
      ENDIF
      DO 20 IX2=1,LINES_LAMBDA_TARG
        TEMP_TARG(IX2)=MOD(LAMBDA_TARG(IX2)-BASLAM,360.)
        IF (TEMP_TARG(IX2).LT.0.) TEMP_TARG(IX2)=TEMP_TARG(IX2)+360.
        TEMP_TARG(IX2)=TEMP_TARG(IX2)+BASLAM
   20 CONTINUE ! Labelled DO for the sake of fpp
C
CL    2  For each target lambda line, find the index of the next source
CL    lambda line to the left.
C
      DO IX2=1,LINES_LAMBDA_TARG
        DO IX1=1,LINES_LAMBDA_SRCE
          IF (TEMP_TARG(IX2).GE.TEMP_SRCE(IX1)) IXL(IX2)=IX1
        ENDDO
      ENDDO
C
CL    3  Find which source lambda gaps cover each target gap and the
CL    fractions they contribute.
C
C     At this point IXL(target_line) gives the index of the next source
C     lambda line to the left of the target lambda line, wrapping round
C     if the source grid is cyclic. This is also the index of the source
C     gap in which the target line falls. Similarly, the index of the
C     target line is also that of the target gap of which it is the
C     left-hand limit. Therefore also IXL(target_gap+1, wrapping round
C     if reqd.), is the index of the source gap which contains the
C     right-hand limit of the target gap. For each target gap, we loop
C     over all source gaps and find the fraction covered by each. Record
C     the fraction and the source index in cumulative lists. If the
C     source grid is not cyclic, parts of the target gap lying outside
C     the source grid are neglected.
C
      IXP=0
      DO IX2=1,GAPS_LAMBDA_TARG
        IX=0
        IX2N=MOD(IX2,LINES_LAMBDA_TARG)+1
        BTARG=TEMP_TARG(IX2N)-TEMP_TARG(IX2)
        IF (BTARG.LT.0.) THEN
          BTARG=BTARG+360.
          IX1N=IXL(IX2N)+LINES_LAMBDA_SRCE
        ELSE
          IX1N=IXL(IX2N)
        ENDIF
        DO IX1W=IXL(IX2),IX1N
          IX1=MOD(IX1W-1,LINES_LAMBDA_SRCE)+1
          IF (CYCLIC_SRCE.OR.IX1.NE.LINES_LAMBDA_SRCE) THEN
            IF (IX1W.EQ.IXL(IX2)) THEN
              BLO=TEMP_TARG(IX2)
            ELSE
              BLO=TEMP_SRCE(IX1)
            ENDIF
            IF (IX1W.EQ.IX1N) THEN
              BHI=TEMP_TARG(IX2N)
            ELSE
              BHI=TEMP_SRCE(MOD(IX1,LINES_LAMBDA_SRCE)+1)
            ENDIF
            IF (BHI.LT.BLO) BHI=BHI+360.
            IF (ABS(BHI-BLO).GT.1E-7*ABS(BTARG)) THEN
              IX=IX+1
              INDEX_LAMBDA(IXP+IX)=IX1
              FRAC_LAMBDA(IXP+IX)=(BHI-BLO)/BTARG
            ENDIF
          ENDIF
        ENDDO
        COUNT_LAMBDA(IX2)=IX
        IXP=IXP+COUNT_LAMBDA(IX2)
      ENDDO
C
CL    4  For each target phi line, find the index of the next source phi
CL    line above. Comments as for section 2, without wrap-round.
C
      LINES_PHI_SRCE=GAPS_PHI_SRCE+1
      LINES_PHI_TARG=GAPS_PHI_TARG+1
      DO IY2=1,LINES_PHI_TARG
        IYT(IY2)=0
        DO IY1=1,LINES_PHI_SRCE
          IF (PHI_TARG(IY2).LE.PHI_SRCE(IY1)) IYT(IY2)=IY1
        ENDDO
      ENDDO
C
CL    5  Find which source phi gaps cover each target gap and the
CL    fractions they contribute. Comments as for section 3, without
CL    wrap-round.
C
      IYP=0
      DO IY2=1,GAPS_PHI_TARG
        IY=0
        IF (SPHERICAL) THEN
C     Contain angle between +-90. There is no real area outside
C     these limits on a sphere.
          BTARG=SIN(MAX(MIN(PHI_TARG(IY2+1),90.),-90.)*PI_OVER_180)
     &    -SIN(MAX(MIN(PHI_TARG(IY2),90.),-90.)*PI_OVER_180)
        ELSE
           BTARG=PHI_TARG(IY2+1)-PHI_TARG(IY2)
        ENDIF
        DO IY1=IYT(IY2),IYT(IY2+1)
          IF (.NOT.(IY1.EQ.0.OR.IY1.EQ.LINES_PHI_SRCE)) THEN
            IF (IY1.EQ.IYT(IY2)) THEN
              BLO=PHI_TARG(IY2)
            ELSE
              BLO=PHI_SRCE(IY1)
            ENDIF
            IF (IY1.EQ.IYT(IY2+1)) THEN
              BHI=PHI_TARG(IY2+1)
            ELSE
              BHI=PHI_SRCE(IY1+1)
            ENDIF
            IF (SPHERICAL) THEN
              BLO=MAX(MIN(BLO,90.),-90.)
              BHI=MAX(MIN(BHI,90.),-90.)
            ENDIF
            IF (ABS(BHI-BLO).GT.1E-7*ABS(BTARG)) THEN
              IY=IY+1
              INDEX_PHI(IYP+IY)=IY1
C             Both numerator and denominator in the following are -ve.
              IF (SPHERICAL) THEN
                FRAC_PHI(IYP+IY)
     &          =(SIN(BHI*PI_OVER_180)-SIN(BLO*PI_OVER_180))/BTARG
              ELSE
                FRAC_PHI(IYP+IY)=(BHI-BLO)/BTARG
              ENDIF
            ENDIF
          ENDIF
        ENDDO
        COUNT_PHI(IY2)=IY
        IYP=IYP+COUNT_PHI(IY2)
      ENDDO
C
CL    6  For each target box, loop over contributing source boxes and
CL    calculate the weights for each one, ignoring land boxes.
C
C     After the first pass for each target box, go back and normalise
C     the weights to compensate for land source boxes and any outside
C     the source grid. Record the source box index and the weight in
C     cumulative lists.
C
      IP=0
      IYP=0
      DO IY2=1,GAPS_PHI_TARG
        IXP=0
        DO IX2=1,GAPS_LAMBDA_TARG
          I=0
          SUM=0
          DO IY=IYP+1,IYP+COUNT_PHI(IY2)
            DO IX=IXP+1,IXP+COUNT_LAMBDA(IX2)
              IF (MASK_SRCE(INDEX_LAMBDA(IX),INDEX_PHI(IY))
     &        .EQV.WANT) THEN
                I=I+1
                INDEX_SRCE(IP+I)=INDEX_LAMBDA(IX)
     &          +(INDEX_PHI(IY)-1)*GAPS_LAMBDA_SRCE
                WEIGHT(IP+I)=FRAC_LAMBDA(IX)*FRAC_PHI(IY)
                SUM=SUM+WEIGHT(IP+I)
              ENDIF
            ENDDO
          ENDDO
          COUNT_TARG(IX2,IY2)=I
          BASE_TARG(IX2,IY2)=IP
          DO I=1,COUNT_TARG(IX2,IY2)
            WEIGHT(IP+I)=WEIGHT(IP+I)/SUM
          ENDDO
          IP=IP+COUNT_TARG(IX2,IY2)
          IXP=IXP+COUNT_LAMBDA(IX2)
        ENDDO
        IYP=IYP+COUNT_PHI(IY2)
      ENDDO
      MAXL=IP
C
      ICODE=0
      CMESSAGE=' '
      RETURN
      END
