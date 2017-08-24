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
!+ Controls horizontal interpolation
!
! Subroutine Interface:
      SUBROUTINE H_INT_CTL(IDIM,LEN_FIELD_OUT
     &,                    ROW_LENGTH_IN,ROW_LENGTH_OUT
     &,                    ROWS_IN,ROWS_OUT,AW_AREA_BOX
     &,                    GLOBAL,H_INT_TYPE
     &,                    AW_INDEX_TARG_LHS,AW_INDEX_TARG_TOP
     &,                    BL_INDEX_B_L,BL_INDEX_B_R
     &,                    AW_COLAT_T,AW_LONG_L,DATA_IN
     &,                    WEIGHT_T_R,WEIGHT_B_R,WEIGHT_T_L,WEIGHT_B_L
     &,                    DATA_OUT)

      IMPLICIT NONE
!
! Description:
!   <Say what this routine does>
!
! Method:
!   <Say how it does it: include references to external documentation>
!   <If this routine is very complex, then include a "pseudo code"
!    description of it to make its structure and method clear>
!
! Current Code Owner: D.M. Goddard
!
! History:
! Version   Date     Comment
! -------   ----     -------
! 4.0      160395   Original code. D.M. Goddard
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!
! System component covered: S121
! System Task:              S1
!
! Declarations:
!   These are of the form:-
!     INTEGER      ExampleVariable      !Description of variable
!
! Global variables (*CALLed COMDECKs etc...):

! Subroutine arguments
!   Scalar arguments with intent(in):
      INTEGER      IDIM             !Second dimension of index arrays
      INTEGER      LEN_FIELD_OUT    !No of points on target grid
      INTEGER      ROWS_IN          !No of rows on source grid
      INTEGER      ROWS_OUT         !No of rows on target grid
      INTEGER      ROW_LENGTH_IN    !No of pts per row on source grid
      INTEGER      ROW_LENGTH_OUT   !No of pts per row on target grid
      LOGICAL      GLOBAL           !True if global area required
      LOGICAL      H_INT_TYPE       !=T Area weighted interpolation;
                                    !=F Bi-linear interpolation

!   Array  arguments with intent(in):
      REAL         AW_AREA_BOX      !area of grid box in sq units of
                                    !  source grid
      INTEGER      AW_INDEX_TARG_LHS(ROW_LENGTH_OUT+1)
                                    !Index of source box overlapping
                                    !lhs of target grid-box
      INTEGER      AW_INDEX_TARG_TOP(ROWS_OUT+1)
                                    !Index of source box overlapping
                                    !top of target grid-box
      INTEGER      BL_INDEX_B_L(LEN_FIELD_OUT)
                                    !Gather index for bottom l.h.c of
                                    !source grid box. 1=P-pts; 2=UV-pts
      INTEGER      BL_INDEX_B_R(LEN_FIELD_OUT)
                                    !Gather index for bottom r.h.c of
                                    !source grid box. 1=P-pts; 2=UV-pts
      REAL         AW_COLAT_T(ROWS_OUT+1)
                                    !Colatitude of top of target grd-box
                                    ! (in units of DELTA_LAT_SRCE)
      REAL         AW_LONG_L(ROW_LENGTH_OUT+1)
                                    !Left longitude of target grid-box
                                    ! (in units of DELTA_LONG_SRCE)
      REAL         DATA_IN(ROW_LENGTH_IN*ROWS_IN)
                                    !Data before interpolation
      REAL         WEIGHT_T_R(LEN_FIELD_OUT) !\ Weights used in
      REAL         WEIGHT_B_R(LEN_FIELD_OUT) ! \bilinear horizontal
      REAL         WEIGHT_T_L(LEN_FIELD_OUT) ! /interpolation
      REAL         WEIGHT_B_L(LEN_FIELD_OUT) !/ 1=P-pts; 2=U-pts
                                             !  3=V-pts; 4=zonal
                                             !             means
!   Array  arguments with intent(out):
      REAL         DATA_OUT(ROW_LENGTH_OUT*ROWS_OUT)
                                    !Data after interpolation

!   ErrorStatus <Delete if ErrorStatus not used>
!     INTEGER      ErrorStatus          ! Error flag (0 = OK)


! Function & Subroutine calls:
      External H_INT_BL,H_INT_AW

!- End of header


      IF(.NOT.H_INT_TYPE)THEN

  ! 1: Bi-linear interpolation requested
        CALL TIMER('H_INT_BL',3)

        CALL H_INT_BL(ROWS_IN,ROW_LENGTH_IN,LEN_FIELD_OUT
     &,               BL_INDEX_B_L,BL_INDEX_B_R,DATA_IN
     &,               WEIGHT_B_L,WEIGHT_B_R
     &,               WEIGHT_T_L,WEIGHT_T_R
     &,               DATA_OUT)

        CALL TIMER('H_INT_BL',4)
      ELSE

  ! 2: Area weighted interpolation
        CALL TIMER('H_INT_AW',3)

        CALL H_INT_AW(ROWS_IN,ROWS_OUT
     &,               ROW_LENGTH_IN,ROW_LENGTH_OUT,GLOBAL
     &,               AW_INDEX_TARG_LHS,AW_INDEX_TARG_TOP
     &,               AW_COLAT_T,AW_LONG_L,DATA_IN,DATA_OUT)

        CALL TIMER('H_INT_AW',4)

      ENDIF   !

      RETURN
      END
