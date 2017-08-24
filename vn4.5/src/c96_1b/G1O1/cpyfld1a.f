C *****************************COPYRIGHT******************************
C (c) CROWN COPYRIGHT 1996, METEOROLOGICAL OFFICE, All Rights Reserved.
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
!+ Parallel UM: Copies one field to another, possibly with different
!               halo sizes
!
! Subroutine interface:
      SUBROUTINE COPY_FIELD(ORIG_FIELD,DEST_FIELD,
     &                      ORIG_FIELD_SIZE,DEST_FIELD_SIZE,
     &                      ORIG_ROW_LENGTH,ORIG_N_ROWS,LEVELS,
     &                      ORIG_Offx,ORIG_Offy,
     &                      DEST_Offx,DEST_Offy,
     &                      L_SWAP)

      IMPLICIT NONE
!
! Description:
! This routine copies one field into another, allowing for a
! different halo size in the two fields. If L_SWAP is true
! it will update the halos on the destination field
!
! Method
! Data is copied from ORIG_FIELD to DEST_FIELD with the Offx and
! Offy values of each field taken into account to work out the offsets.
! Data will be copied from the halo areas of ORIG_FIELD into
! corresponding halo areas of DEST_FIELD.
! A call to SWAPBOUNDS will update the halos if L_SWAP is .TRUE.
!
! Current Code Owner : Paul Burton
!
! History:
!  Model    Date      Modification history from model version 4.1
!  version
!    4.1    23/11/95   New DECK created for the Parallel Unified
!                      Model. P.Burton
!    4.3    25/07/97   Remove initialisation of Dest_Field
!                      if L_SWAP is true   D.Salmond
!
! Subroutine Arguments:

      INTEGER ORIG_FIELD_SIZE, ! IN horizontal size of ORIG_FIELD
     &        DEST_FIELD_SIZE, ! IN horizontal size of DEST_FIELD
     &        ORIG_ROW_LENGTH, ! IN row length of ORIG_FIELD
     &        ORIG_N_ROWS,     ! IN number of rows in ORIG_FIELD
     &        LEVELS,          ! IN number of levels in both fields
     &        ORIG_Offx,       ! IN halo size of ORIG in X direction
     &        ORIG_Offy,       ! IN halo size of ORIG in Y direction
     &        DEST_Offx,       ! IN halo size of DEST in X direction
     &        DEST_Offy        ! IN halo size of DEST in Y direction

      LOGICAL L_SWAP           ! IN do a halo swap of DEST_FIELD?

      REAL    ORIG_FIELD(ORIG_FIELD_SIZE,LEVELS),
     &                                       ! IN Field to copy from
     &        DEST_FIELD(DEST_FIELD_SIZE,LEVELS)
     &                                       ! OUT Field to copy to

! Local variables

      INTEGER  DEST_ROW_LENGTH,  ! row length of DEST_FIELD
     &         DEST_N_ROWS,      ! number of rows in DEST_FIELD
     &         MIN_ROW_LENGTH,   ! smallest row length in ORIG_FIELD
     &                           ! or DEST_FIELD
     &         MIN_N_ROWS,       ! smallest number of rows in ORIG_FIELD
     &                           ! or DEST_FIELD
     &         ORIG_Off_X,       ! X offset in ORIG_FIELD for copy loop
     &         ORIG_Off_Y,       ! Y offset in ORIG_FIELD for copy loop
     &         DEST_Off_X,       ! X offset in DEST_FIELD for copy loop
     &         DEST_Off_Y,       ! Y offset in DEST_FIELD for copy loop
     &         ORIG_INDEX,       ! point in horizontal ORIG_FIELD
     &         DEST_INDEX        ! point in horizontal DEST_FIELD

      INTEGER I,J,K  ! loop counters (column,row,level)

! ------------------------------------------------------------------

! Calculate the shape of DEST_FIELD
      DEST_ROW_LENGTH = ORIG_ROW_LENGTH - ORIG_Offx*2 + DEST_Offx*2
      DEST_N_ROWS     = ORIG_N_ROWS     - ORIG_Offy*2 + DEST_Offy*2

! Set DEST_FIELD to some "safe" value for all locations that aren't
! set in the copy
      IF(.NOT.L_SWAP)THEN
      DO K=1,LEVELS
        DO I=1,DEST_FIELD_SIZE
          DEST_FIELD(I,K)=0.0
        ENDDO
      ENDDO
      ENDIF

! Calculate the smallest size in each horizontal dimension
      MIN_ROW_LENGTH  = MIN(ORIG_ROW_LENGTH,DEST_ROW_LENGTH)
      MIN_N_ROWS      = MIN(ORIG_N_ROWS,DEST_N_ROWS)

! Calculate the offsets for the copy loop
      ORIG_Off_X=(ORIG_ROW_LENGTH-MIN_ROW_LENGTH)/2
      ORIG_Off_Y=(ORIG_N_ROWS-MIN_N_ROWS)/2
      DEST_Off_X=(DEST_ROW_LENGTH-MIN_ROW_LENGTH)/2
      DEST_Off_Y=(DEST_N_ROWS-MIN_N_ROWS)/2

! Copy from ORIG_FIELD to DEST_FIELD
      DO K=1,LEVELS
        DO J=1,MIN_N_ROWS
          DO I=1,MIN_ROW_LENGTH
            DEST_INDEX=I+DEST_Off_X+(J+DEST_Off_Y-1)*DEST_ROW_LENGTH
            ORIG_INDEX=I+ORIG_Off_X+(J+ORIG_Off_Y-1)*ORIG_ROW_LENGTH
            DEST_FIELD(DEST_INDEX,K) = ORIG_FIELD(ORIG_INDEX,K)
          ENDDO
        ENDDO
      ENDDO

      IF (L_SWAP) THEN
! Do a swap to update halos
        CALL SWAPBOUNDS(DEST_FIELD,DEST_ROW_LENGTH,DEST_N_ROWS,
     &                  DEST_Offx,DEST_Offy,LEVELS)
      ENDIF

      RETURN
      END
