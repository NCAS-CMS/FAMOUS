C ******************************COPYRIGHT******************************
C (c) CROWN COPYRIGHT 1998, METEOROLOGICAL OFFICE, All Rights Reserved.
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
CLL Subroutine INTF_NEW_FILES -----------------------------------------
C
CLL  Purpose: To test whether a new boundary data file needs to be
CLL           opened
CLL
CLL  Model            Modification history from model version 4.5
CLL version  Date
CLL  4.5   3/09/98    New deck added M.J.Bell
CLL
CLLEND ---------------------------------------------------------------
C ---------------------------------------------------------------------
      subroutine intf_new_files(first_unit, last_unit, max_n_intf, im,
     #    TYPE_LETTER_1, FT_OUTPUT, INTF_FREQ_HR, FT_STEPS, STEP,
     #    FT_FIRSTSTEP, INTERFACE_STEPS,
     #    LNEWBND )

CL  Purpose: determines new output interface files for a submodel.
CL
      implicit none
      integer first_unit ! IN first unit to test
      integer last_unit  ! IN last unit to test
      integer max_n_intf ! IN number of interface files
      integer im         ! IN  sub-model identifier
      character*1 TYPE_LETTER_1(20:last_unit) ! IN
      character*1 FT_OUTPUT(20:last_unit)     ! IN
      integer INTF_FREQ_HR(max_n_intf)     ! IN
      integer FT_STEPS(20:last_unit)          ! IN
      integer STEP                         ! IN model step no.
      integer FT_FIRSTSTEP(20:last_unit)      ! IN
      integer INTERFACE_STEPS(max_n_intf)  ! IN
      logical LNEWBND(max_n_intf)          ! OUT
C-----------------------------------------------------------------------
CL Declaration of local variables
      integer iunit
!     logical ll_intf_type   ! OUT T => file is an output interface file
      integer jintf          ! boundary file area number
C-----------------------------------------------------------------------
      do iunit = first_unit, last_unit

       if (type_letter_1(iunit).eq.'b') then  !  Boundary file

         IF (FT_OUTPUT(IUNIT).EQ.'Y') THEN ! Intf. data output?
          call intf_area ( im, iunit, JINTF)

          IF ( INTF_FREQ_HR(JINTF) .GT. 0) THEN

            IF (STEP.EQ.0) THEN

              LNEWBND(JINTF) = .TRUE. ! New intf data file required at
                                      ! first entry to ININTF1
            ELSE

              IF (FT_STEPS(IUNIT).EQ.0) LNEWBND(JINTF) =. FALSE. !False
                                         ! if incomplete single file

              IF (FT_STEPS(IUNIT).GT.0) LNEWBND(JINTF) = .NOT.(
C               step = first timestep to get boundary data
     *          (STEP-FT_FIRSTSTEP(IUNIT).EQ.0 .OR.
C               step = timestep to start new file
     *          MOD(STEP-FT_FIRSTSTEP(IUNIT),FT_STEPS(IUNIT)).NE.0)
     *          .AND.
     &     STEP.GT.FT_FIRSTSTEP(IUNIT)-INTERFACE_STEPS(JINTF))
C                 ! False if incomplete file in sequence
            END IF  ! STEP

           ENDIF  !  INTF_FREQ_HR
C
          ELSE    !  FT_OUTPUT(IUNIT)
C
C  Possible place for setting switches for reading in interface data
C
          ENDIF   !  FT_OUTPUT
        ENDIF  !  TYPE_LETTER_1

      ENDDO   ! IUNIT

      return
      end
