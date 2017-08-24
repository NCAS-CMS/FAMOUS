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
!+ MPP version of routine to merge LBCs into full fields
!
! Subroutine Interface
      SUBROUTINE MPP_MERGEFLD(ROW_LENGTH,ROWS,
     &                        LENRIM,RIMWIDTH,RIMWEIGHTS,
     &                        FLD_TYPE,RIM,FIELD)
      IMPLICIT NONE
!
! Description:
! Merges Lateral Boundary Conditions (LBCs) with full field data
! (This is an MPP version of the MERGEFLD routine, upon which
!  this code is based)
!
! Method:
! Each boundary is treated seperately. If this processor
! lies on a particular boundary, the relevant part of that
! LBC is merged into the data held on the processor, using
! the weighting factor RIMWEIGHTS.
! The appropriate element of the RIMWEIGHTS array is shown
! in the diagram below (which assumes RIMWIDTH=4). Take
! note of what happens at the corners! It also shows the
! structure of the LBC arrays: Northern, Eastern, Southern
! and Western Boundaries.
!
!               1 1 1 1 1 1 1 1 1 1 1 1
! Northern -->  1 2 2 2 2 2 2 2 2 2 2 1
! Boundary      1 2 3 3 3 3 3 3 3 3 2 1
!               1 2 3 4 4 4 4 4 4 3 2 1
!              -------------------------
!               1 2 3 4 . . . . 4 3 2 1
! Western --->  1 2 3 4 . . . . 4 3 2 1  <-- Eastern
! Boundary      1 2 3 4 . . . . 4 3 2 1      Boundary
!               1 2 3 4 . . . . 4 3 2 1
!              -------------------------
!               1 2 3 4 4 4 4 4 4 3 2 1
! Southern -->  1 2 3 3 3 3 3 3 3 3 2 1
! Boundary      1 2 2 2 2 2 2 2 2 2 2 1
!               1 1 1 1 1 1 1 1 1 1 1 1
!
! IMPORTANT NOTE:
! This routine assumes that the rim fits entirely onto the
! processors at the edge of the LPG. If the rim overlaps
! into processors inside the LPG then this code will not
! work.
!
! Current code owner : Paul Burton
!
! History
!  Model    Date       Modification history from model version 4.1
!  version
!    4.1    16/1/96    New Deck for MPP code   P.Burton
!
! Subroutine Arguments:

      INTEGER
     &   ROW_LENGTH    ! IN length of local rows (inc. MPP halos)
     &,  ROWS          ! IN number of local rows (inc. MPP halos)
     &,  LENRIM        ! IN length of local LBC data
     &,  RIMWIDTH      ! IN width of boundary zone
     &,  FLD_TYPE      ! IN indicates type of field (P or U)

      REAL
     &   RIMWEIGHTS(RIMWIDTH)   ! IN weights to be given to boundary
!                               !    zone values
     &,  RIM(LENRIM)            ! IN boundary data to merge
     &,  FIELD(ROW_LENGTH*ROWS) ! INOUT field to merge boundary data
!                               !       into

! Parameters and COMMON
! ------------------------ Comdeck PARVARS -------------------------
! Parameters and common blocks required by the MPP-UM
!========================== COMDECK PARPARM ====================
!   Description:
!
!   This COMDECK contains PARAMETERs for the MPP-UM
!
!   Two sets of parameters are set up -
!     i)  for the MPP-UM itself.
!     ii) for the interface to the Message Passing Software.
!
!   History:
!
!   Model    Date     Modification history
!  version
!   4.1      27/1/96  New comdeck based on first section of
!                     old PARVARS.   P.Burton
!   4.2      21/11/96 Add new field type parameter and
!                     magic number used in addressing to indicate
!                     if a calculation is for local data, or data
!                     on the dump on disk (ie. global data)  P.Burton
!   4.2      18/11/96 Moved MaxFieldSize to comdeck AMAXSIZE and
!                     removed Maxbuf.  P.Burton
!   4.2      18/7/96  Removed some unused variables      P.Burton
!   4.4      11/07/97 Reduced MAXPROC to 256 to save memory  P.Burton
!
! ---------------------- PARAMETERS ---------------------
!
! =======================================================
! Parameters needed for the MPP-UM
! =======================================================

      INTEGER   Ndim_max        ! maximum number of spatial dimensions
      PARAMETER (Ndim_max = 3 ) ! 3d data


      INTEGER
     &   fld_type_p           ! indicates a grid on P points
     &,  fld_type_u           ! indicates a grid on U points
     &,  fld_type_unknown     ! indicates a non-standard grid.
      PARAMETER (
     &   fld_type_p=1
     &,  fld_type_u=2
     &,  fld_type_unknown=-1)

      INTEGER
     &   local_data
     &,  global_dump_data
      PARAMETER (
     &   local_data=1        ! Used in addressing to indicate if
     &,  global_dump_data=2) ! calculation is for a local or
!                            ! global (ie. disk dump) size

! =======================================================
! Parameters needed for the Message Passing Software
! =======================================================


      INTEGER
     &   Maxproc              ! Max number of processors
      PARAMETER (
     &   MAXPROC = 256)

      INTEGER
     &   PNorth       ! North processor address in the neighbour array
     &,  PEast        ! East  processor address in the neighbour array
     &,  PSouth       ! South processor address in the neighbour array
     &,  PWest        ! West  processor address in the neighbour array
     &,  NoDomain     ! Value in neighbour array if the domain has
     &                !  no neighbor in this direction. Otherwise
     &                !  the value will be the tid of the neighbor
      PARAMETER (
     &   PNorth   = 1
     &,  PEast    = 2
     &,  PSouth   = 3
     &,  PWest    = 4
     &,  NoDomain = -1)

      INTEGER
     &   BC_STATIC            ! Static boundary conditions
     &,  BC_CYCLIC            ! Cyclic boundary conditions
      PARAMETER (
     &   BC_STATIC = 1
     &,  BC_CYCLIC = 2)

! ---------------------- End of comdeck PARPARM ---------------------
!========================== COMDECK PARCOMM ====================
!
! *** NOTE : This comdeck requires comdeck PARPARM to be *CALLed
!            first.
!
!   Description:
!
!   This COMDECK contains COMMON blocks for the MPP-UM
!
!
!   Two COMMON blocks are defined:
!     i)  UM_PARVAR holds information required by the
!         Parallel Unified Model itself
!     ii) MP_PARVAR holds information required by the interface to
!         the Message Passing Software used by the PUM
!
!   Key concepts used in the inline documentation are:
!     o GLOBAL data - the entire data domain processed by the UM
!     o LOCAL data - the fragment of the GLOBAL data which is
!       stored by this particular process
!     o PERSONAL data - the fragment of the LOCAL data which is
!       updated by this particular process
!     o HALO data - a halo around the PERSONAL data which forms
!       the LOCAL data
!
!     Acronyms used:
!     LPG - Logical Process Grid, this is the grid of logical
!           processors; each logical processor handles one of the
!           decomposed parts of the global data. It does not
!           necessarily represent a physical grid of processors.
!
!   History:
!
!   4.1      27/1/96  New comdeck based on second section of
!                     old PARVARS.   P.Burton
!   4.2     19/08/96  Removed some unused variables, and added
!                     current_decomp_type variable to allow use
!                     of flexible decompositions.
!                     Added nproc_max to indicate the max. number
!                     of processors used for MPP-UM
!                                                      P.Burton
!
! -------------------- COMMON BLOCKS --------------------
!
! =======================================================
! Common block for the Parallel Unified Model
! =======================================================

      INTEGER
     &   first_comp_pe       ! top left pe in LPG
     &,  last_comp_pe        ! bottom right pe in LPG
     &,  current_decomp_type ! current decomposition type
     &,  Offx                ! halo size in EW direction
     &,  Offy                ! halo size in NS direction
     &,  glsize(Ndim_max)    ! global data size
     &,  lasize(Ndim_max)    ! local data size
     &,  blsizep(Ndim_max)   ! personal p data area
     &,  blsizeu(Ndim_max)   ! personal u data area
     &,  datastart(Ndim_max) ! position of personal data in global data
     &                       !   (in terms of standard Fortran array
     &                       !    notation)
     &,  gridsize(Ndim_max)  ! size of the LPG in each dimension
     &,  gridpos(Ndim_max)   ! position of this process in the LPG
!                            ! 0,1,2,...,nproc_x-1 etc.

      LOGICAL
     &    atbase             ! process at the bottom of the LPG
     &,   attop              ! process at the top of the LPG
     &,   atleft             ! process at the left of the LPG
     &,   atright            ! process at the right of the LPG
! NB: None of the above logicals are mutually exclusive

      COMMON /UM_PARVAR/
     &                  first_comp_pe,last_comp_pe
     &,                 current_decomp_type,Offx, Offy
     &,                 glsize,lasize,blsizep,blsizeu
     &,                 datastart,gridsize,gridpos
     &,                 atbase,attop,atleft,atright

! =======================================================
! Common block for the Message Passing Software
! =======================================================

      INTEGER
     &  bound(Ndim_max)           ! type of boundary (cyclic or static)
     &                            !  in each direction
     &, g_lasize(Ndim_max,0:maxproc)
!                                 ! global copy of local data size
     &, g_blsizep(Ndim_max,0:maxproc)
!                                 ! global copy of personal p data area
     &, g_blsizeu(Ndim_max,0:maxproc)
!                                 ! global copy of personal u data area
     &, g_datastart(Ndim_max,0:maxproc)
!                                 ! global copy of datastart
     &, g_gridpos(Ndim_max,0:maxproc)
!                                 ! global copy of gridpos
     &, nproc                     ! number of processors in current
!                                 ! decomposition
     &, nproc_max                 ! maximum number of processors
     &, nproc_x                   ! number of processors in x-direction
     &, nproc_y                   ! number of processors in y-direction
     &, mype                      ! number of this processor
     &                            !  (starting from 0)
     &, neighbour(4)              ! array with the tids of the four
     &                            ! neighbours in the horizontal plane
     &, gc_proc_row_group         ! GID for procs along a proc row
     &, gc_proc_col_group         ! GID for procs along a proc col
     &, gc_all_proc_group         ! GID for all procs

      COMMON /MP_PARVAR/
     &                  bound
     &,                 g_lasize,g_blsizep,g_blsizeu
     &,                 g_datastart,g_gridpos
     &,                 nproc,nproc_max,nproc_x,nproc_y,mype
     &,                 neighbour,gc_proc_row_group
     &,                 gc_proc_col_group, gc_all_proc_group



! ---------------------- End of comdeck PARCOMM -----------------------
! --------------------- End of comdeck PARVARS ---------------------

! Local variables

      INTEGER
     &  I    ! loop index East-West   : NB I and J take no account
     &, J    ! loop index North-South :    of the MPP halos
     &, IRIM ! position in RIM data array
     &, IFLD ! position in FIELD data array
     &, LBC_ROW_LEN  ! length of rows in LBC array
     &, EW_START_ROW  ! row to start merging in LBCs at : east/west
     &, EW_END_ROW    ! row to end merging in LBCs at   : boundaries
     &, S_START_ROW   ! row to start memrging in LBCs at : southern
     &, S_END_ROW    ! row to end merging in LBCs at    : boundary

      REAL
     &  RWT  ! rimweight of point being processed

! 0.0 Calculate the length of the LBC rows, and where to start and
!     end merging in the side LBC rows, and the southern LBC rows

        LBC_ROW_LEN=ROW_LENGTH-2*Offx
!       For most cases the length of the LBC row is just the
!       processor's row length minus the two MPP halos.
!       However, for U fields, the last point on each global row
!       is dropped in the LBC field:
        IF ((FLD_TYPE .EQ. fld_type_u) .AND. (atright))
     &    LBC_ROW_LEN=LBC_ROW_LEN-1

        EW_START_ROW=Offy+1
        EW_END_ROW=ROWS-Offy
!       This is a first-guess. If the processor is at the top or
!       bottom of the grid, it starts getting more complicated...

        IF (attop) THEN
!         At the top of the grid, the first RIMWIDTH rows are
!         dealt with at the Northern boundary
          EW_START_ROW=EW_START_ROW+RIMWIDTH
        ENDIF

        IF (atbase) THEN
!         At the bottom of the grid, the last RIMWIDTH rows are
!         dealt with at the Southern boundary
          EW_END_ROW=EW_END_ROW-RIMWIDTH

!         If we're at the base then we need to know where to
!         apply the Southern boundary rows
          S_START_ROW=ROWS-Offy-RIMWIDTH+1

          IF (FLD_TYPE .EQ. fld_type_u) THEN
!           For U fields we loose a row at the bottom of the
!           global grid
            EW_END_ROW=EW_END_ROW-1
            S_START_ROW=S_START_ROW-1
          ENDIF
          S_END_ROW=S_START_ROW+RIMWIDTH-1
        ENDIF


! 1.0  Merge Northern Boundary

      IRIM=1
      IF (attop) THEN  ! if this processor at Northern boundary

        DO J=1,RIMWIDTH
          DO I=1,LBC_ROW_LEN

            IFLD=I+Offx+(J+Offy-1)*ROW_LENGTH  ! point in FIELD

            RWT=RIMWEIGHTS(J)  ! simplistic first guess at RIMWEIGHT
!                              ! If we are at a corner it starts to
!                              ! get a bit more involved.....
            IF ((atleft) .AND. (I .LT. RIMWIDTH)) THEN
!             This is a left hand corner
              RWT=RIMWEIGHTS(MIN(I,J))
            ELSEIF ((atright) .AND.
     &              (I .GT. LBC_ROW_LEN-RIMWIDTH)) THEN
!             This is a right hand corner
              RWT=RIMWEIGHTS(MIN(LBC_ROW_LEN-I+1,J))
            ENDIF

            FIELD(IFLD)=RIM(IRIM)*RWT + FIELD(IFLD)*(1.0-RWT)

            IRIM=IRIM+1

          ENDDO ! I loop over points along row
        ENDDO ! J loop over rows in boundary area

      ENDIF ! If I'm a processor at the Northern boundary area

      IRIM=LBC_ROW_LEN*RIMWIDTH+1

! 2.0  Merge Eastern Boundary

      IF (atright) THEN  ! if this processor at Eastern boundary

        DO J=EW_START_ROW,EW_END_ROW
          DO I=1,RIMWIDTH

            IFLD=I+Offx+LBC_ROW_LEN-RIMWIDTH+(J-1)*ROW_LENGTH

            RWT=RIMWEIGHTS(RIMWIDTH+1-I)

            FIELD(IFLD)=RIM(IRIM)*RWT + FIELD(IFLD)*(1.0-RWT)

            IRIM=IRIM+1

          ENDDO ! I loop over points along rim
        ENDDO ! J loop over rows in boundary area

      ENDIF ! If I'm a processor at the Eastern boundary area

      IRIM=LBC_ROW_LEN*RIMWIDTH + (ROWS-2*Offy)*RIMWIDTH + 1

! 3.0  Merge Southern Boundary

      IF (atbase) THEN  ! if this processor at Southern boundary

        DO J=S_START_ROW,S_END_ROW
          DO I=1,LBC_ROW_LEN

            IFLD=I+Offx+(J-1)*ROW_LENGTH  ! point in FIELD

!            RWT=RIMWEIGHTS(J+1-S_START_ROW)
            RWT=RIMWEIGHTS(RIMWIDTH-(J-S_START_ROW))

!           And now the special case of the corners...
            IF ((atleft) .AND. (I .LT. RIMWIDTH)) THEN
!             This is a left hand corner
              RWT=RIMWEIGHTS(MIN(I,RIMWIDTH-(J-S_START_ROW)))
            ELSEIF ((atright) .AND.
     &              (I .GT. LBC_ROW_LEN-RIMWIDTH)) THEN
!             This is a right hand corner
              RWT=RIMWEIGHTS(MIN(LBC_ROW_LEN-I+1,
     &                           RIMWIDTH-(J-S_START_ROW)))
            ENDIF

            FIELD(IFLD)=RIM(IRIM)*RWT + FIELD(IFLD)*(1.0-RWT)

            IRIM=IRIM+1
          ENDDO ! I loop over points along row
        ENDDO ! J loop over rows in boundary area

      ENDIF ! If I'm a processor at the Southern boundary area

      IRIM=2*LBC_ROW_LEN*RIMWIDTH + (ROWS-2*Offy)*RIMWIDTH + 1

! 4.0 Merge Western Boundary

      IF (atleft) THEN  ! if this processor at Western boundary

        DO J=EW_START_ROW,EW_END_ROW
          DO I=1,RIMWIDTH

            IFLD=I+Offx+(J-1)*ROW_LENGTH

            RWT=RIMWEIGHTS(I)

            FIELD(IFLD)=RIM(IRIM)*RWT + FIELD(IFLD)*(1.0-RWT)

            IRIM=IRIM+1

          ENDDO ! I loop over points along rim
        ENDDO ! J loop over rows in boundary area

      ENDIF ! If I'm a processor at the Western boundary area


      RETURN

      END
