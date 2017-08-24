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
CLL  SUBROUTINE CLOUD_COVER--------------------------------------------
CLL
CLL     PURPOSE:
CLL This routine combines model level cloud into three categories;
CLL High,Medium and Low. The criteria for the upper and lower
CLL boundaries were provided by C.F.O and are input to this routine.
CLL It also provides total cloud cover as it would be with either
CLL random or maximum/random cloud overlap.
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  Date
CLL   3.1   8/2/92    Total cloud cover diagnostics added.
CLL   3.2    13/07/93 Changed CHARACTER*(*) to CHARACTER*(80) for
CLL                   portability.  Author Tracey Smith.
CLL   3.3    10/12/93 Change *DEF to include A18_1A  Bruce M
CLL   3.3  15/12/93         Calculation of TCA_MAXRAN corrected to allow
CLL                                    for zero convective cloud cover.
!     4.4  31/07/97   Add cloud frequency diagnostic. Andrew Bushell
CLL
CLL  Programming standard: U M DOC  Paper NO. 4,
CLL
CLL  Logical components covered : D433
CLL
CLL  Project task: D4
CLL
CLL  External documentation  UMDP
CLL
CLLEND-------------------------------------------------------------
C
C*L ARGUMENTS:-----------------------------------------------------
      SUBROUTINE CLOUD_COVER
     * (MODEL_CLOUD, CONVECTIVE_CLOUD, CCB, CCT,
     * LOW_CLOUD, MED_CLOUD, HIGH_CLOUD, TCA_RANDOM, TCA_MAXRAN,
     & LAYER_CLOUD_PRESENT,
     * LOW_BOT_LEVEL,LOW_TOP_LEVEL,
     * MED_BOT_LEVEL,MED_TOP_LEVEL,HIGH_BOT_LEVEL,HIGH_TOP_LEVEL,
     * FLAG_LOW, FLAG_MED, FLAG_HIGH, FLAG_TRA, FLAG_TMR,
     & FLAG_LCP,CLOUD_LEVELS,POINTS,Q_LEVELS,ICODE,CMESSAGE)
C
      IMPLICIT NONE
      CHARACTER*(80) CMESSAGE  ! OUT Return message
      INTEGER
     * ICODE          ! OUT Return code  0 normal exit GT 0 error
     *,POINTS         ! IN  NO of Points to be processed.
     &,Q_LEVELS       ! IN  NO of HUMIDITY LEVELS
     &,CLOUD_LEVELS   ! IN  NO of RADIATION CLOUD LEVELS
     *,LOW_BOT_LEVEL  ! IN  The lowest level of the LOW category
     *,LOW_TOP_LEVEL  ! IN  The top level of the LOW category
     *,MED_BOT_LEVEL  ! IN  The lowest level of the MED  category
     *,MED_TOP_LEVEL  ! IN  The top level of the MED category
     *,HIGH_BOT_LEVEL ! IN  The lowest level of the Highest category
     *,HIGH_TOP_LEVEL ! IN  The top level of the Highest category
     *,CCB(POINTS)    ! IN  Convective cloud base
     *,CCT(POINTS)    ! IN  Convective cloud top
C
      REAL
     * MODEL_CLOUD(POINTS,CLOUD_LEVELS)   ! IN The model's layer cloud
     *,CONVECTIVE_CLOUD(POINTS)           ! IN Its convective cloud
     *,LOW_CLOUD(POINTS)   ! OUT Cloudy fraction in the lowest lyr
     *,MED_CLOUD(POINTS)   ! OUT Cloudy fraction for the medium lyr
     *,HIGH_CLOUD(POINTS)  ! OUT Cloudy fraction for the highest lyr
     *,TCA_RANDOM(POINTS)  ! OUT Total cloud amount with random overlap
     *,TCA_MAXRAN(POINTS)  ! OUT         and with maximum/random overlap
     &,LAYER_CLOUD_PRESENT(POINTS,CLOUD_LEVELS)   ! OUT Lcld frequency
C
      LOGICAL
     * FLAG_LOW         ! True if required FALSE if not
     *,FLAG_MED         !                "
     *,FLAG_HIGH        !                "
     *,FLAG_TRA         !                "
     *,FLAG_TMR         !                "
     &,FLAG_LCP         !                "
C*---------------------------------------------------------------
C
C*L WORKSPACE USAGE----------------------------------------------
C    One dynamically allocated array MAX_CONTIG:
      REAL MAX_CONTIG(POINTS)        ! Maximum total cloud cover in the
C     ! layer currently being considered and those below it through
C     ! which cloud extends contiguously.
C*---------------------------------------------------------------
C
C*L EXTERNAL SUBROUTINES CALLED----------------------------------
C   None.
C*---------------------------------------------------------------
C
C
C----------------------------------------------------------------
C   DEFINE LOCAL VARIABLES
C----------------------------------------------------------------
      INTEGER
     * I,J            !  LOOP COUNTERS
      REAL TOCLE      ! Total cloud in this layer
C

C------------------------------------------------------------------
CL  1. Calculate the cloudy fraction of each of the layers
C------------------------------------------------------------------

CL  1. Calculate the cloudy fraction for LOW

      IF(FLAG_LOW) THEN

        IF(LOW_BOT_LEVEL.GT.LOW_TOP_LEVEL)THEN
          ICODE=1
          CMESSAGE='CLDCVR :CLOUD_TYPE LOW Bottom lvl above top level'
          GOTO 9999
        ENDIF

        IF(LOW_TOP_LEVEL.GT.CLOUD_LEVELS) THEN
          ICODE=1
          CMESSAGE='CLDCVR :CLOUD_TYPE LOW TOP level above CLOUD_LEVELS'
          GOTO 9999
        ENDIF

C  Intialise the LOW_CLOUD ARRAY to the bottom layer amount.
        DO  I=1,POINTS
          LOW_CLOUD(I)=MODEL_CLOUD(I,LOW_BOT_LEVEL)
        ENDDO

C  Find the max fraction over the required layers
CL *** Following loop labelled due to fmp mistranslation
C
        DO 100 J=LOW_BOT_LEVEL+1,LOW_TOP_LEVEL
        DO I=1,POINTS
        IF(LOW_CLOUD(I).LT.MODEL_CLOUD(I,J)) THEN
          LOW_CLOUD(I)=MODEL_CLOUD(I,J)
        ENDIF
        ENDDO
 100    CONTINUE
      ENDIF
CL  2. Calculate the cloudy fraction for MEDIUM

      IF(FLAG_MED) THEN

        IF(MED_BOT_LEVEL.GT.MED_TOP_LEVEL)THEN
          ICODE=1
          CMESSAGE='CLDCVR:CLOUD_TYPE MED Bottom level above top level'
          GOTO 9999
        ENDIF

        IF(MED_TOP_LEVEL.GT.CLOUD_LEVELS) THEN
          ICODE=1
          CMESSAGE='CLDCVR:CLOUD_TYPE MED TOP level above CLOUD_LEVELS'
          GOTO 9999
        ENDIF

C  Intialise the MED_CLOUD ARRAY to the bottom layer amount.
        DO  I=1,POINTS
          MED_CLOUD(I)=MODEL_CLOUD(I,MED_BOT_LEVEL)
        ENDDO

C  Find the max fraction over the required layers
CL *** Following loop labelled due to fmp mistranslation
C
        DO 200 J=MED_BOT_LEVEL+1,MED_TOP_LEVEL
        DO I=1,POINTS
        IF(MED_CLOUD(I).LT.MODEL_CLOUD(I,J)) THEN
          MED_CLOUD(I)=MODEL_CLOUD(I,J)
        ENDIF
        ENDDO
 200    CONTINUE

      ENDIF
CL  4. Calculate the cloudy fraction for high

      IF(FLAG_HIGH) THEN

        IF(HIGH_BOT_LEVEL.GT.HIGH_TOP_LEVEL)THEN
          ICODE=1
          CMESSAGE='CLDCVR: CLOUD_TYPE HIGH Bottom lvl above top level'
          GOTO 9999
        ENDIF

        IF(HIGH_TOP_LEVEL.GT.CLOUD_LEVELS) THEN
          ICODE=1
          CMESSAGE='CLDCVR:CLOUD_TYPE HIGH TOP lvl above CLOUD_LEVELS'
          GOTO 9999
        ENDIF

C  Intialise the HIGH_CLOUD ARRAY to the bottom layer amount.
        DO  I=1,POINTS
          HIGH_CLOUD(I)=MODEL_CLOUD(I,HIGH_BOT_LEVEL)
        ENDDO

C  Find the max fraction over the required layers
CL *** Following loop labelled due to fmp mistranslation
C
        DO 300 J=HIGH_BOT_LEVEL+1,HIGH_TOP_LEVEL
        DO I=1,POINTS
        IF(HIGH_CLOUD(I).LT.MODEL_CLOUD(I,J)) THEN
          HIGH_CLOUD(I)=MODEL_CLOUD(I,J)
        ENDIF
        ENDDO
 300    CONTINUE


      ENDIF

C-----------------------------------------------------------------------
CL  2. Calculate the total cloud amounts
C-----------------------------------------------------------------------

C     These diagnostics are calculated just as in LWDCSF - see that
C       routine for comments on how it is done.

      IF ( FLAG_TRA ) THEN

        DO I=1, POINTS
          TCA_RANDOM(I) = 1. - CONVECTIVE_CLOUD(I)
        ENDDO

        DO 500 J=1, CLOUD_LEVELS
          DO I=1, POINTS
            TCA_RANDOM(I) = TCA_RANDOM(I) * ( 1. - MODEL_CLOUD(I,J) )
          ENDDO
  500   CONTINUE

        DO I=1, POINTS
          TCA_RANDOM(I) = 1. - TCA_RANDOM(I)
        ENDDO
C
      ENDIF

      IF ( FLAG_TMR ) THEN

        DO I=1, POINTS
          TCA_MAXRAN(I) = 1.
          MAX_CONTIG(I) = 0.
        ENDDO

        DO 600 J=1, CLOUD_LEVELS
          DO I=1, POINTS
          IF ( MODEL_CLOUD(I,J) .EQ. 0. .AND.
     &  ( CONVECTIVE_CLOUD(I) .EQ. 0. .OR.
     &                   J .LT. CCB(I) .OR. J .GE. CCT(I) )  ) THEN
               TCA_MAXRAN(I) = TCA_MAXRAN(I) * ( 1. - MAX_CONTIG(I) )
               MAX_CONTIG(I) = 0.
             ELSE
               TOCLE = MODEL_CLOUD(I,J)
               IF ( J .GE. CCB(I) .AND. J .LT. CCT(I) )
     &            TOCLE = TOCLE + CONVECTIVE_CLOUD(I) * ( 1. - TOCLE )
               MAX_CONTIG(I) = MAX(MAX_CONTIG(I),TOCLE)
            ENDIF
          ENDDO
  600   CONTINUE

        DO I=1, POINTS
          TCA_MAXRAN(I) = 1. - TCA_MAXRAN(I) * ( 1. - MAX_CONTIG(I) )
        ENDDO

      ENDIF
!-----------------------------------------------------------------------
!   3. Calculate the layer cloud frequency
!-----------------------------------------------------------------------
      IF ( FLAG_LCP ) THEN
        DO J=1, Q_LEVELS
          DO I=1, POINTS
            IF (MODEL_CLOUD(I,J) .LE. 0.)  THEN
              LAYER_CLOUD_PRESENT(I,J)=0.
            ELSE
              LAYER_CLOUD_PRESENT(I,J)=1.
            ENDIF
          END DO  ! I_do_lcp
        END DO  ! J_do_lcp
      ENDIF  ! Flag_lcp_if
!
 9999 CONTINUE
      RETURN
      END
