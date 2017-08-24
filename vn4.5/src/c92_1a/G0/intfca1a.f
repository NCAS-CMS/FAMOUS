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
CLL  Routine: INTF_COAST_AJ------------------------------------------
CLL
CLL  Purpose: To calculate a suitable value of NSEARCH for call to
CLL           SPIRAL_S and then call SPIRAL_S
CL
CLL  Author: D.M.Goddard      Date: 17 November 1993
CLL  Reviewer:  Date of review:
CLL
CLL  Tested under compiler: cft77
CLL  Tested under OS version: UNICOS 7
CLL
CLL  Code version no: 1       Date: 17 November 1993
CLL
CLL  Modification History:
CLL
CLL  Model
CLL  Version    Date
CLL  4.0        26/07/95  Make local copy of land sea mask to work
CLL                       to correct error of mask not being reset
CLL                       after call to SPIRAL_S.
CLL  4.0        02/11/95  The Fortran for squaring an array needed to be
CLL       redefined for the DecAlpha because of lexcon. (N.Farnon)
CLL  Author: C.P.Jones           Reviewer  D.M.Goddard
!    4.5        30/07/97  Skip put of loop when nearest neighbour found.
!                         Author D.M. Goddard
CLL
CLL Programming standard: UM Doc Paper 3, version
CLL
CLL Logucal component number:
CLL
CLL Project task: S1
CLL
CLL
CLL Documentation:
CLL   UMDP S1
CLL
CLL -------------------------------------------------------------
      SUBROUTINE INTF_COAST_AJ
     &           (LAND_SEA_MASK,INDEX_UNRES,NO_POINT_UNRES,
     &           POINTS_PHI,POINTS_LAMBDA,DATA_FIELD,SEA_LAND,
     &           CYCLIC,MAXDIM)

      IMPLICIT NONE

C*L ARGUMENTS:---------------------------------------------------

      INTEGER
     & POINTS_PHI       !IN number of rows in grid
     &,POINTS_LAMBDA    !IN number of columns in grid
     &,NO_POINT_UNRES   !INOUT number of unresolved points
     &,LAND_SEA_MASK(POINTS_LAMBDA*POINTS_PHI)
     &                  !IN land sea mask
     &,INDEX_UNRES(POINTS_LAMBDA*POINTS_PHI)
     &                  !INOUT index to unresolved pts
     &,SEA_LAND         !IN =0 for sea field  =1/-1 for land field

      REAL
     & DATA_FIELD(POINTS_LAMBDA*POINTS_PHI) !IN field

      LOGICAL
     & CYCLIC           ! IN =T if data covers complete latitude circle

C*L LOCAL VARIABLES

      INTEGER
     & I,J,JJ,JJJ,K,JK,NSCH ! indices
     &,IPT,IROW,ICOL  ! coordinate of unresolved point
     &,IPOINT,IUNRES  ! do loop variable
     &,NPOINTS        ! number of points in serach box
     &,MAXDIM        ! largest dimension of field
     &,IR((1+2*MAXDIM)*(1+2*MAXDIM))
     &                   ! row numbers of points to serach
     &,IC((1+2*MAXDIM)*(1+2*MAXDIM))
     &                   ! col numbers of points to search
     &,IND_SEARCH((1+2*MAXDIM)*(1+2*MAXDIM))
     &                   ! index to points to search
     &               ! still unresolved after calling this subroutine
     &,ISUM_MASK     ! number of surrounding points which have data
     &,ISEARCH       ! largest dimension of field
     &,NSEARCH       ! minimum search radius required.
     &,LAND_SEA_TEMP(POINTS_LAMBDA*POINTS_PHI) ! local copy of mask


C*L   EXTERNAL ROUTINES
      EXTERNAL SPIRAL_S

      NSEARCH=0

C Take local copy of land sea mask to work with
      DO IPOINT=1,POINTS_LAMBDA*POINTS_PHI
        LAND_SEA_TEMP(IPOINT)=LAND_SEA_MASK(IPOINT)
      ENDDO

C toggle land sea mask to exclude unresolved points from meaning process
      DO IUNRES=1,NO_POINT_UNRES
        IF(SEA_LAND.EQ.0)THEN
          LAND_SEA_TEMP(INDEX_UNRES(IUNRES))=1
        ELSE
          LAND_SEA_TEMP(INDEX_UNRES(IUNRES))=0
        ENDIF
      ENDDO


C Loop around unresolved points
      DO IUNRES=1,NO_POINT_UNRES

C find unresolved point coordinate in terms of rows and cols
        IPT=INDEX_UNRES(IUNRES)
        IROW=INT(REAL(IPT)/REAL(POINTS_LAMBDA)+1)
        ICOL=IPT-(IROW-1)*POINTS_LAMBDA

        ISEARCH=MAXDIM
        DO I=1,MAXDIM

C Calculate number of points in search box
          NPOINTS=(1+2*I)**2 ! number of grid points in search box

C calculate surrounding points' coords in terms of rows and cols
          JJJ=1
          DO J=-I,I
            DO JJ=JJJ,JJJ+2*I
              IR(JJ)=IROW+J
            ENDDO
            JJJ=JJJ+1+2*I
          ENDDO

          JJJ=1+2*I
          JK=1
          DO J=-I,I
            DO JJ=0,2*I
              IC(JK+JJ*JJJ)=ICOL+J
            ENDDO
            JK=JK+1
          ENDDO

C Check that col and rows are in range of grid
          DO IPOINT=1,NPOINTS
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
          ENDDO

c Form index search array
          DO IPOINT=1,NPOINTS
            IND_SEARCH(IPOINT)=(IR(IPOINT)-1)*POINTS_LAMBDA+IC(IPOINT)
          ENDDO

C search for data around this point. If no data is found the point
C remains unresolved

          ISUM_MASK=0   ! number of points with data found

          DO IPOINT=1,NPOINTS
           IF(IABS(LAND_SEA_TEMP(IND_SEARCH(IPOINT))).EQ.IABS(SEA_LAND)
     *       .AND.DATA_FIELD(IND_SEARCH(IPOINT)).GE.0.0)THEN
             ISUM_MASK=ISUM_MASK+1
           ENDIF
         ENDDO

          IF(ISUM_MASK.GT.0)THEN
            ISEARCH=MIN0(ISEARCH,I)
            GOTO 100
          END IF

        ENDDO

 100    CONTINUE
        NSEARCH=MAX0(ISEARCH,NSEARCH)

      ENDDO


      CALL SPIRAL_S(LAND_SEA_TEMP,INDEX_UNRES,NO_POINT_UNRES,
     *                 POINTS_PHI,POINTS_LAMBDA,DATA_FIELD,NSEARCH,
     *                 SEA_LAND,CYCLIC)


      RETURN
      END
