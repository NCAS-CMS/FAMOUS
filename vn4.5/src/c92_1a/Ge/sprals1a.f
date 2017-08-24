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
CLL SUBROUTINE SPIRAL_S-----------------------------------------------
CLL
CLL Written by C.P. Jones
CLL
CLL Reviewed by ??
CLL
CLL Programming standard:
CLL    Unified Model Documentation Paper No 3
CLL
CLL System component: S121
CLL
CLL System task: S1
CLL
CLL Purpose:
CLL   Attempts to set a value at points which are unresolved when
CLL   interpolating between one grid and another.  A value is set
CLL   by finding the mean of surrounding points which do have data
CLL   set within a search radius determined by NSEARCH.
CLL
CLL  Modification History:
CLL
CLL  Model
CLL  Version  Date
CLL  4.0      02/11/95  The Fortran for squaring an array needed to be
CLL     redefined for the DecAlpha because of lexcon. (N.Farnon)
CLL  4.4      23/09/97 amendment to calculation of row number and
CLL                    to check for missing data (M. J. Bell)       
CLL Documentation:
CLL   UMDP S1
CLL
CLL -------------------------------------------------------------
      SUBROUTINE SPIRAL_S(LAND_SEA_MASK,INDEX_UNRES,NO_POINT_UNRES,
     &           POINTS_PHI,POINTS_LAMBDA,DATA_FIELD,NSEARCH,SEA_LAND,
     &           CYCLIC)

      IMPLICIT NONE

C*L ARGUMENTS:---------------------------------------------------

      INTEGER
     & POINTS_PHI       !IN number of rows in grid
     &,POINTS_LAMBDA    !IN number of columns in grid
     &,NSEARCH          !IN number of points in each direction to search
     &,NO_POINT_UNRES   !INOUT number of unresolved points
     &,LAND_SEA_MASK(POINTS_LAMBDA*POINTS_PHI)
     &                  !IN land sea mask
     &,INDEX_UNRES(POINTS_LAMBDA*POINTS_PHI)
     &                  !INOUT index to unresolved pts
     &,SEA_LAND         !IN =0 for sea field  =1/-1 for land field

      REAL
     & DATA_FIELD(POINTS_LAMBDA*POINTS_PHI) !INOUT field

      LOGICAL
     & CYCLIC           ! IN =T if data covers complete latitude circle

C*L Parameters
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

C*L LOCAL VARIABLES

      INTEGER
     & I,J,JJ,JJJ,K,JK  ! indices
     &,IPT,IROW,ICOL  ! coordinate of unresolved point
     &,IPOINT,IUNRES  ! do loop variable
     &,NPOINTS        ! number of points in serach box
     &,IR((1+2*NSEARCH)*(1+2*NSEARCH))
     &                   ! row numbers of points to serach
     &,IC((1+2*NSEARCH)*(1+2*NSEARCH))
     &                   ! col numbers of points to search
     &,IND_SEARCH((1+2*NSEARCH)*(1+2*NSEARCH))
     &                   ! index to points to search
CLL
     &,NOT_YET_SET                  ! number of points still to set
     &,IND_YET_SET(POINTS_LAMBDA*POINTS_PHI) ! index of points
     &               ! still unresolved after calling this subroutine
     &,ISUM_MASK     ! number of surrounding points which have data

      REAL
     & SUM_DATA      ! sum of data surrounding unresolved points
     &, RMDI_TOL ! values within this tolerance counted as missing  

C*L   EXTERNAL ROUTINES
C     None
C---------------------------------------------------------------------
      RMDI_TOL  = ABS (RMDI) * 0.0001
C

C Calculate number of points in search box
      NPOINTS=(1+2*NSEARCH)**2 ! number of grid points in search box

C Loop around unresolved points
      NOT_YET_SET=0
      DO 10 IUNRES=1,NO_POINT_UNRES

C find unresolved point coordinate in terms of rows and cols
        IPT=INDEX_UNRES(IUNRES)
        IROW= (IPT - 1)/POINTS_LAMBDA   +   1   
        ICOL=IPT-(IROW-1)*POINTS_LAMBDA

C calculate surrounding points' coords in terms of rows and cols
        JJJ=1
        DO 20 J=-NSEARCH,NSEARCH
          DO 30 JJ=JJJ,JJJ+2*NSEARCH
            IR(JJ)=IROW+J
 30       CONTINUE
        JJJ=JJJ+1+2*NSEARCH
 20     CONTINUE

        JJJ=1+2*NSEARCH
        JK=1
        DO 40 J=-NSEARCH,NSEARCH
          DO 50 JJ=0,2*NSEARCH
            IC(JK+JJ*JJJ)=ICOL+J
 50       CONTINUE
        JK=JK+1
 40     CONTINUE

C Check that col and rows are in range of grid
        DO 70 IPOINT=1,NPOINTS
          IF(IC(IPOINT).GT.POINTS_LAMBDA) THEN
            IF(CYCLIC) THEN
              IC(IPOINT)=IC(IPOINT)-POINTS_LAMBDA
            ELSE
              IC(IPOINT)=POINTS_LAMBDA
            ENDIF
          ENDIF
          IF(IC(IPOINT).LT.1) THEN
            IF(CYCLIC) THEN
              IC(IPOINT)=IC(IPOINT)+POINTS_LAMBDA
            ELSE
              IC(IPOINT)=1
            ENDIF
          ENDIF
          IF(IR(IPOINT).LT.1) IR(IPOINT)=1
          IF(IR(IPOINT).GT.POINTS_PHI) IR(IPOINT)=POINTS_PHI
 70     CONTINUE

c Form index search array
        DO 80 IPOINT=1,NPOINTS
          IND_SEARCH(IPOINT)=(IR(IPOINT)-1)*POINTS_LAMBDA+IC(IPOINT)
 80     CONTINUE

C search for data around this point. If no data is found the point
C remains unresolved

        ISUM_MASK=0   ! number of points with data found
        SUM_DATA=0.0  ! sum of data of surrounding grid points

        DO 90 IPOINT=1,NPOINTS
          IF(IABS(LAND_SEA_MASK(IND_SEARCH(IPOINT))).EQ.IABS(SEA_LAND)
     &    .AND.DATA_FIELD(IND_SEARCH(IPOINT)) .GT. RMDI+RMDI_TOL) THEN
            SUM_DATA=SUM_DATA+DATA_FIELD(IND_SEARCH(IPOINT))
            ISUM_MASK=ISUM_MASK+1
          ENDIF
 90     CONTINUE

       IF(ISUM_MASK .GT. 0) THEN
C data found - take mean
          DATA_FIELD(IPT)=SUM_DATA/REAL(ISUM_MASK)
        ELSE
C data not found - point remains unresolved
          NOT_YET_SET=NOT_YET_SET+1
          IND_YET_SET(NOT_YET_SET)=IPT
        ENDIF

 10   CONTINUE


C amend output array with points remaining unresolved
      IF(NOT_YET_SET.GT.0) THEN
        DO 100 IPOINT=1,NOT_YET_SET
          INDEX_UNRES(IPOINT)=IND_YET_SET(IPOINT)
 100    CONTINUE
        NO_POINT_UNRES=NOT_YET_SET
      ELSE
        NO_POINT_UNRES=0
      ENDIF
      RETURN
      END
