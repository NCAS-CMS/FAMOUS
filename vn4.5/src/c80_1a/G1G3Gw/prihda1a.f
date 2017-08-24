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
CLL  SUBROUTINE PR_INHDA---------------------------------------
CLL
CLL  Purpose: Prints out integer constants record and checks
CLL           validity of information.
CLL
CLL  Written by A. Dickinson
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  date
CLL   3.2    13/07/93 Changed CHARACTER*(*) to CHARACTER*(80) for
CLL                   portability.  Author Tracey Smith.
CLL   4.5    19/01/98 Remove SOIL_VARS and VEG_VARS. D. Robinson
CLL
CLL  Programming standard:  Unified Model Documentation Paper No 3
CLL                         Version No 1 15/1/90
CLL
CLL  System component: C25
CLL
CLL  System task: F3
CLL
CLL  Documentation: Unified Model Documentation Paper No F3
CLL                 Version No 5 9/2/90
CLL------------------------------------------------------------
C*L Arguments:-------------------------------------------------
      SUBROUTINE PR_INHDA
     *(INTHD,LEN_INTHD,ROW_LENGTH,P_ROWS
     &,P_LEVELS,Q_LEVELS,TR_LEVELS
     &,ST_LEVELS,SM_LEVELS,BL_LEVELS
     & , TR_VARS
     *,ICODE,CMESSAGE)

      IMPLICIT NONE

      INTEGER
     * LEN_INTHD      !IN Length of integer header
     *,ROW_LENGTH     !IN No of points along a model row
     *,P_ROWS         !IN No of model rows
     *,P_LEVELS       !IN No of model levels
     *,Q_LEVELS       !IN No of moist levels
     *,TR_LEVELS      !IN No of tracer levels
     &,ST_LEVELS      !IN No of soil temperature levels
     &,SM_LEVELS      !IN No of soil moisture levels
     *,BL_LEVELS      !IN No of b. layer levels
     *,TR_VARS        !IN No of tracer variables

      INTEGER
     & INTHD(LEN_INTHD) !IN Integer header
     *,ICODE          !OUT Return code; successful=0
     *                !                 error > 0

      CHARACTER*(80)
     * CMESSAGE       !OUT Error message if ICODE > 0

C -------------------------------------------------------------
C Workspace usage:---------------------------------------------
C None
C -------------------------------------------------------------
C*L External subroutines called:-------------------------------
C None
C*-------------------------------------------------------------

CL Internal structure: None
      ICODE=0
      CMESSAGE=' '

      WRITE(6,'('' '')')
      WRITE(6,'('' INTEGER CONSTANTS'')')
      WRITE(6,'('' -----------------'')')
      WRITE(6,'('' Number of timesteps since start of run -'',I7)')
     *INTHD(1)
      WRITE(6,'('' Meaning interval for mean fields (hrs) -'',I7)')
     *INTHD(2)
      WRITE(6,'('' Number of dumps used to generate mean  -'',I7)')
     *INTHD(3)
      WRITE(6,'('' No of hrs between neighbouring contiguous sections of
     * means -'',I7)')INTHD(4)
      WRITE(6,'('' No of hrs between end of one contiguous section and s
     *tart of next -'',I7)')INTHD(5)

      IF(ROW_LENGTH.NE.INTHD(6).OR.P_ROWS.NE.INTHD(7))THEN
      WRITE(6,'('' **FATAL ERROR** specifying model dimensions'')')
      WRITE(6,'('' Program dimensions ='',I7,'' x'',I7)')
     *ROW_LENGTH,P_ROWS
      WRITE(6,'('' File dimensions ='',I7,'' x'',I7)')
     *INTHD(6),INTHD(7)
      WRITE(6,'('' ***************'')')
        ICODE=4
        CMESSAGE='PR_INHDA: Consistency check'
        RETURN
      ELSE
      WRITE(6,'('' Number of E-W x N-S points -'',I7,'' x'',I7)')
     *INTHD(6),INTHD(7)
      ENDIF

      IF(P_LEVELS.NE.INTHD(8).OR.Q_LEVELS.NE.INTHD(9))THEN
      WRITE(6,'('' **FATAL ERROR** specifying no of model levels'')')
      WRITE(6,'('' Programmed levels (wet) ='',I2,''('',I2,'')'')')
     *P_LEVELS,Q_LEVELS
      WRITE(6,'('' File levels (wet) ='',I2,''('',I2,'')'')')
     *INTHD(8),INTHD(9)
      WRITE(6,'('' ***************'')')
        ICODE=4
        CMESSAGE='PR_INHDA: Consistency check'
        RETURN
      ELSE
      WRITE(6,'('' Number of levels -'',I7)')
     *INTHD(8)
      WRITE(6,'('' Number of wet levels -'',I7)')
     *INTHD(9)
      ENDIF
      IF((TR_LEVELS.NE.INTHD(12).AND.INTHD(14).NE.0)
     &    .OR.ST_LEVELS.NE.INTHD(10).OR.SM_LEVELS.NE.INTHD(28)
     &    .OR.BL_LEVELS.NE.INTHD(13))THEN
      WRITE(6,'('' **FATAL ERROR** specifying model level info'')')
      WRITE(6,'('' Programmed tracer, soil temperature,'',
     &            '' soil moisture,  b.l. levels ='',4I3)')
     &          TR_LEVELS,ST_LEVELS,SM_LEVELS,BL_LEVELS           
      WRITE(6,'('' File tracer, soil temperature, soil moisture,'',
     &          ''  b.l. levels ='',4I3)')
     &          INTHD(12),INTHD(10),INTHD(28),INTHD(13)
      WRITE(6,'('' ***************'')')
        ICODE=4
        CMESSAGE='PR_INHDA: Consistency check'
        RETURN
      ELSE
      WRITE(6,'('' Number of soil temperature levels -'',I7)') 
     &INTHD(10)
      WRITE(6,'('' Number of soil moisture levels -'',I7)')
     &INTHD(28)
      WRITE(6,'('' Number of tracer levels -'',I7)')
     *INTHD(12)
      WRITE(6,'('' Number of boundary layer levels -'',I7)')
     *INTHD(13)
      ENDIF

      IF (TR_VARS.NE.INTHD(14)) THEN
      WRITE(6,'('' **FATAL ERROR** specifying number of variables'')')
      WRITE (6,'('' Programmed TR_VARS = '',I3)') TR_VARS
      WRITE (6,'('' File       TR_VARS = '',I3)') INTHD(14)
      WRITE(6,'('' ***************'')')
        ICODE=4
        CMESSAGE='PR_INHDA: Consistency check'
        RETURN
      ELSE
      WRITE(6,'('' Number of passive tracers advected -'',I7)')
     *INTHD(14)
      ENDIF

      RETURN
      END
