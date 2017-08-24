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
      SUBROUTINE DO_AREAVER(GAPS_LAMBDA_SRCE,GAPS_PHI_SRCE,LROW_SRCE
     &,INVERT_SRCE,DATA_SRCE,GAPS_LAMBDA_TARG,GAPS_PHI_TARG,COUNT_TARG
     &,BASE_TARG,LROW_TARG,WANT,MASK_TARG,INDEX_SRCE,WEIGHT,ADJUST
     &,DATA_TARG,ICODE,CMESSAGE)
CLL   Subroutine DO_AREAVER -------------------------------------------
CLL
CLL Purpose:
CLL
CLL   Perform area-averaging to transform data from the source grid to
CLL   the target grid, or adjust the values on the source grid to have
CLL   the area-averages supplied on the target grid. The latter mode
CLL   is intended for adjusting values obtained by interpolating from
CLL   "target" to "source" in order to conserve the area-averages.
CLL   This mode should be used ONLY if each source box belongs in
CLL   exactly one target box. ADJUST=0 selects normal area-averaging,
CLL   ADJUST=1 selects adjustment by addition (use this mode for fields
CLL   which may have either sign), ADJUST=2 selects adjustment by
CLL   multiplication (for fields which are positive-definite or
CLL   negative-definite).
CLL
CLL   The shape of the source and target grids are specified by their
CLL   dimensions GAPS_aa_bb, which give the number of gaps in the
CLL   aa=LAMBDA,PHI coordinate in the bb=SRCE,TARG grid. (The product
CLL   of GAPS_LAMBDA_bb and GAPS_PHI_bb is the number of boxes in the
CLL   bb grid.)
CLL
CLL   The input and output data are supplied as 2D arrays DATA_SRCE and
CLL   DATA_TARG, whose first dimensions should also be supplied. Speci-
CLL   fying these sizes separately from the actual dimensions of the
CLL   grids allows for columns and rows in the arrays to be ignored.
CLL   A target land/sea mask should be supplied in MASK_TARG, with the
CLL   value indicating wanted points specified in WANT. Points which
CLL   are unwanted or which lie outside the source grid are not altered
CLL   in DATA_TARG. DATA_SRCE can optionally be supplied with its rows
CLL   in reverse order (i.e. with the first row corresponding to
CLL   minimum LAMBDA).
CLL
CLL   The arrays COUNT_TARG, BASE_TARG, INDEX_SRCE and WEIGHT should be
CLL   supplied as returned by PRE_AREAVER q.v.
CLL
CLL   Programming Standard, paper 4 version 4 (14.12.90)
CLL
CLL Modification history:
CLL
CLL Logical components covered :
CLL
CLL Project task :
CLL
CLL External documentation: Unified Model documentation paper No:
CLL                         Version:
CLL
CLLEND -----------------------------------------------------------------
C
      IMPLICIT NONE
C*L
      INTEGER
     & GAPS_LAMBDA_SRCE        !IN number lambda gaps in source grid
     &,GAPS_PHI_SRCE           !IN number phi gaps in source grid
     &,LROW_SRCE               !IN first dimension of source arrays
     &,GAPS_LAMBDA_TARG        !IN number lambda gaps in target grid
     &,GAPS_PHI_TARG           !IN number phi gaps in target grid
     &,LROW_TARG               !IN first dimension of target arrays
     &,COUNT_TARG(GAPS_LAMBDA_TARG,GAPS_PHI_TARG)
C                              !IN no. of source boxes in target box
     &,BASE_TARG(GAPS_LAMBDA_TARG,GAPS_PHI_TARG)
C                              !IN first index in list for target box
     &,INDEX_SRCE(*)           !IN list of source box indices
     &,ADJUST                  !IN selects normal or adjust mode
     &,ICODE                   !OUT return code
      LOGICAL
     & INVERT_SRCE             !IN DATA_SRCE rows in reverse order
     &,WANT                    !IN indicator of wanted points in mask
     &,MASK_TARG(LROW_TARG,*)  !IN land/sea mask for target grid
C     NB alternative intents below apply for normal/adjust mode
      REAL
     & DATA_SRCE(LROW_SRCE,*)  !IN/INOUT data on source grid
     &,WEIGHT(*)               !IN list of weights for source boxes
     &,DATA_TARG(LROW_TARG,*)  !INOUT/IN data on target grid
      CHARACTER
     & CMESSAGE*(*)            !OUT error message
C*
      INTEGER
     & IP                      ! pointer into lists
     &,I                       ! loop index
     &,IX1(GAPS_LAMBDA_SRCE*GAPS_PHI_SRCE)
C                              ! working SRCE LAMBDA indices
     &,IY1(GAPS_LAMBDA_SRCE*GAPS_PHI_SRCE)
C                              ! working SRCE PHI indices
     &,IX2,IY2                 ! working TARG LAMBDA/PHI indices
      REAL
     & TEMP_TARG               ! workspace for area-average
     &,DELTA                   ! additive adjustment
     &,RATIO                   ! multiplicative adjustment
C
CL    Loop over all target boxes and calculate values as required.
C
C     The weights and source box indices are recorded in continuous
C     lists. COUNT_TARG indicates how many consecutive entries in these
C     lists apply to each target box.
C
      DO IY2=1,GAPS_PHI_TARG
        DO IX2=1,GAPS_LAMBDA_TARG
          IF ((MASK_TARG(IX2,IY2).EQV.WANT)
     &    .AND.COUNT_TARG(IX2,IY2).NE.0) THEN
            TEMP_TARG=0.
            DO I=1,COUNT_TARG(IX2,IY2)
              IP=BASE_TARG(IX2,IY2)+I
              IX1(I)=MOD(INDEX_SRCE(IP)-1,GAPS_LAMBDA_SRCE)+1
              IY1(I)=(INDEX_SRCE(IP)-1)/GAPS_LAMBDA_SRCE+1
              IF (INVERT_SRCE) IY1(I)=GAPS_PHI_SRCE-IY1(I)+1
              TEMP_TARG=TEMP_TARG+WEIGHT(IP)*DATA_SRCE(IX1(I),IY1(I))
            ENDDO
            IF (ADJUST.EQ.0) THEN
              DATA_TARG(IX2,IY2)=TEMP_TARG
            ELSEIF (ADJUST.EQ.1) THEN
              DELTA=DATA_TARG(IX2,IY2)-TEMP_TARG
              DO I=1,COUNT_TARG(IX2,IY2)
                DATA_SRCE(IX1(I),IY1(I))=DATA_SRCE(IX1(I),IY1(I))+DELTA
              ENDDO
            ELSEIF (ADJUST.EQ.2.AND.TEMP_TARG.NE.0.) THEN
              RATIO=DATA_TARG(IX2,IY2)/TEMP_TARG
              DO I=1,COUNT_TARG(IX2,IY2)
                DATA_SRCE(IX1(I),IY1(I))=DATA_SRCE(IX1(I),IY1(I))*RATIO
              ENDDO
            ENDIF
          ENDIF
        ENDDO
      ENDDO
C
      ICODE=0
      CMESSAGE=' '
      RETURN
      END
