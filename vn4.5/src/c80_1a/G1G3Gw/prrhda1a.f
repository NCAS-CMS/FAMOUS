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
CLL  SUBROUTINE PR_REHDA---------------------------------------
CLL
CLL  Purpose: Prints out real constants record and checks
CLL           validity of information.
CLL
CLL  Written by A. Dickinson 28/12/89
CLL
CLL  Model            Modification history from model version 3.0:
CLL version  date
CLL
CLL  Programming standard: Unified Model Documentation Paper No 3
CLL                        Version No 1 15/1/90
CLL
CLL  System component: C25
CLL
CLL  System task: F3
CLL
CLL  Documentation: Unified Model Documentation Paper No F3
CLL           Version No 5 9/2/90
CLL
CLL------------------------------------------------------------
C*L Arguments:-------------------------------------------------
      SUBROUTINE PR_REHDA(REALHD,LEN_REALHD)

      IMPLICIT NONE

      INTEGER
     * LEN_REALHD !IN Length of real header

      REAL
     * REALHD(LEN_REALHD) !IN Real header

C -------------------------------------------------------------
C Workspace usage:---------------------------------------------
C None
C -------------------------------------------------------------
C*L External subroutines called:-------------------------------
C None
C*-------------------------------------------------------------
C Local variables:---------------------------------------------
      INTEGER I
C--------------------------------------------------------------

CL Internal structure: None

      WRITE(6,'('' '')')
      WRITE(6,'('' REAL CONSTANTS'')')
      WRITE(6,'('' --------------'')')
      WRITE(6,'('' E-W grid spacing in degrees -'',e12.4)')
     *REALHD(1)
      WRITE(6,'('' N-S grid spacing in degress -'',e12.4)')
     *REALHD(2)
      WRITE(6,'('' Latitude of first row in degrees -'',e12.4)')
     *REALHD(3)
      WRITE(6,'('' Longitude of first point in a row in degrees -'',
     *e12.4)')REALHD(4)
      WRITE(6,'('' Real latitude of pseudo North Pole in degrees - '',
     *e12.4)')REALHD(5)
      WRITE(6,'('' Real longitude of pseudo North Pole in degrees - '',
     *e12.4)')REALHD(6)
      WRITE(6,'('' Grid orientation in degrees - '',
     *e12.4)')REALHD(7)
      WRITE(6,'(8X,''                   Year        Day       Hour
     *Minute     Second'')')
      WRITE(6,'('' Atmosphere time = '',
     *5F12.4)')(REALHD(I),I=8,12)

      WRITE(6,'('' Mass, energy, energy drift = '',3e12.4)')
     *REALHD(19),REALHD(20),REALHD(21)

      RETURN
      END
