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
C
!+ Initialise model for submodel and internal model coupling
!
! Subroutine Interface:
      SUBROUTINE UM_Submodel_Init(ErrorStatus)

      IMPLICIT NONE
!
! Description:
!   UM_Submodel_Init initialises the model with information specifying
!   internal model and submodel partitions for the run, which is
!   required for control of coupling when more than one internal model
!   is present.
!
! Method:
!   The routine reads information from the user interface, providing
!   lists of internal models and their associated submodel data
!   partitions. This is required in both the reconfiguration and the
!   model as a prior step to calculating addressing in STASH_PROC.
!
! Current Code Owner: R. Rawlins
! History:
! Version   Date     Comment
! -------   ----     -------
! 3.5    07/04/95   Original code. R. Rawlins.
!LL 4.3-4.4   16/09/97 D1 addressing change and subsequent correction
!LL                    S.D.Mullerworth
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!
! System component covered:
! System Task:
! Declarations:
!
!
! Global variables (*CALLed COMDECKs etc...):
!
! Description:
!    Describes the number and identity of submodels available
!    within the system, and those included in the current
!    experiment.  Parameters set by the User Interface give
!    the relevant array sizes; other submodel configuration
!    information is either read from NAMELIST input, or
!    derived from dump header information.
!
! Current Code Owner: R. Rawlins
!
! History:
! Version  Date     Comment
! -------  ----     -------
! pre 3.0           Original code. T. Johns
! 3.5    07/04/95   Expansion for stage 1 of submodel project, allowing
!                   flexible specification of internal models within
!                   submodel partitions. R. Rawlins
!
! Declarations:
!
!  1. Internal model and submodel dump partition identifiers - fixed
!     for all experiments.
!
! Description:
!    Hold parameters defining internal model identifiers and submodel
!    data partition (ie main D1 data array and consequent dump), both
!    short and long form.
!
! Current Code Owner: R. Rawlins
!
! History:
! Version  Date     Comment
! -------  ----     -------
! pre 3.0           Original code. T. Johns
! 3.3    26/10/93   M. Carter. Part of an extensive mod that:
!                    1.Removes the limit on primary STASH item numbers.
!                    2.Removes the assumption that (section,item)
!                      defines the sub-model.
!                    3.Thus allows for user-prognostics.
!                    Add index to submodel home dump.
! 3.5    13/03/95   Expansion for stage 1 of submodel project, allowing
!                   flexible specification of internal models within
!                   submodel partitions. R. Rawlins
!
! Declarations:
!
!   Hold parameters defining internal model identifiers and submodel
!   data partition (ie main D1 data array and consequent dump), both
!   short and long form
      INTEGER
     *   A_IM,ATMOS_IM        ! Atmosphere internal model
     *  ,O_IM,OCEAN_IM        ! Ocean      internal model
     *  ,S_IM, SLAB_IM        ! Slab       internal model
     *  ,W_IM, WAVE_IM        ! Wave       internal model
     *  ,I_IM,SEAICE_IM       ! Sea-ice    internal model
     *  ,N_IM,NATMOS_IM       ! New dynamics (Charney-Phillips grid)
!                               atmosphere internal model
!
      PARAMETER(
     *   A_IM=1,ATMOS_IM=1       ! Atmosphere internal model
     *  ,O_IM=2,OCEAN_IM=2       ! Ocean      internal model
     *  ,S_IM=3, SLAB_IM=3       ! Slab       internal model
     *  ,W_IM=4, WAVE_IM=4       ! Wave       internal model
     *  ,I_IM=5,SEAICE_IM=5      ! Sea-ice    internal model
     *  ,N_IM=6,NATMOS_IM=6      ! New dynamics (Charney-Phillips grid)
!                                  atmosphere internal model
     *)
!
      INTEGER
     *   A_SM,ATMOS_SM        ! Atmosphere submodel partition
     *  ,O_SM,OCEAN_SM        ! Ocean      submodel partition
     *  ,W_SM, WAVE_SM        ! Wave       submodel partition
     *  ,N_SM,NATMOS_SM       ! New dynamics (Charney-Phillips grid)
!                                  atmosphere internal model
!
      PARAMETER(
     *   A_SM=1,ATMOS_SM=1    ! Atmosphere submodel partition
     *  ,O_SM=2,OCEAN_SM=2    ! Ocean      submodel partition
     *  ,W_SM=4, WAVE_SM=4    ! Wave       submodel partition
     *  ,N_SM=6,NATMOS_SM=6   ! New dynamics (Charney-Phillips grid)
!                                  atmosphere internal model
     *)
!
C

!
!  2. Maximum internal model/submodel array sizes for this version.
!
!
! Description:
!    Describes the number and identity of submodels available
!    within the system, and those included in the current
!    experiment.  Parameters set by the User Interface give
!    the relevant array sizes; other submodel configuration
!    information is either read from NAMELIST input, or
!    derived from dump header information.
!
! Current Code Owner: R. Rawlins
!
! History:
! Version  Date     Comment
! -------  ----     -------
! 3.5    13/07/95   Original code. D.M. Goddard
! 4.0     3/11/95   Reduce max internal model, submodel from 10 to 4
!                   to save space in model. At 4.0 the max no of 
!                   supported models is 3, 1 slot is reserved for
!                   expansion. Rick Rawlins.         
!  4.1  21/02/96  Wave model introduced as 4th sub-model.  RTHBarnes
!
! Declarations:
!
!
!  1. Maximum internal model/submodel array sizes for this version.
!
      INTEGER
     * N_INTERNAL_MODEL_MAX      ! Max no. of internal models
     *,N_SUBMODEL_PARTITION_MAX  ! Max no. of submodel dump partitions
     *,INTERNAL_ID_MAX           ! Max value of internal model id
     *,SUBMODEL_ID_MAX           ! Max value of submodel dump id

      PARAMETER(
     * N_INTERNAL_MODEL_MAX=4,                                          
     * N_SUBMODEL_PARTITION_MAX=4,
     * INTERNAL_ID_MAX=N_INTERNAL_MODEL_MAX,
     * SUBMODEL_ID_MAX=N_SUBMODEL_PARTITION_MAX)
!
!  3. Lists of internal models and their submodel dump partitions -
!     initialised by the user interface - experiment specific.
      INTEGER
     * N_INTERNAL_MODEL          ! No. of internal models
     *,N_SUBMODEL_PARTITION      ! No. of submodel partitions
     *,INTERNAL_MODEL_LIST(N_INTERNAL_MODEL_MAX) ! Internal models
     *,SUBMODEL_FOR_IM    (N_INTERNAL_MODEL_MAX) ! Submodel identifier
     *                           ! for each internal model in list
     &,SUBMODEL_FOR_SM(N_INTERNAL_MODEL_MAX) ! Submodel number for
!                                  each submodel id
!
! Namelist for information in 3.
      NAMELIST/NSUBMODL/N_INTERNAL_MODEL,N_SUBMODEL_PARTITION
     *,INTERNAL_MODEL_LIST,SUBMODEL_FOR_IM
!
!  4. Lists calculated in model from user interface supplied arrays -
!     - experiment specific.
      INTEGER
     * N_INTERNAL_FOR_SM(SUBMODEL_ID_MAX)  ! No of internal models in
!              each submodel partition indexed by sm identifier
     *,SUBMODEL_PARTITION_LIST(N_SUBMODEL_PARTITION_MAX)    ! List of
!              submodel partition identifiers
     *,SUBMODEL_PARTITION_INDEX(INTERNAL_ID_MAX)  ! Submodel partition
!              identifier indexed by internal model identifier
     *,INTERNAL_MODEL_INDEX(INTERNAL_ID_MAX)      ! Sequence number of
!              internal model indexed by internal model identifier:
!              required to map from id to STASH internal model sequence
      LOGICAL
     * LAST_IM_IN_SM(INTERNAL_ID_MAX)      ! Last internal model within
!                                a submodel partition if .TRUE.,
!                                indexed by internal model id.
! Common block for information in 3. and 4.
      COMMON/SUBMODL/N_INTERNAL_MODEL,N_SUBMODEL_PARTITION,
     *     INTERNAL_MODEL_LIST,SUBMODEL_FOR_IM,SUBMODEL_FOR_SM,
     *     N_INTERNAL_FOR_SM,SUBMODEL_PARTITION_LIST,
     *     SUBMODEL_PARTITION_INDEX,
     *     INTERNAL_MODEL_INDEX,
     *     LAST_IM_IN_SM

!
!  5. Time information specifying coupling frequencies between internal
!     models and submodels, and multipliers, indexed by sequence of
!     internal models and submodels (ie left to right along node tree).
!     {Not required at this release}.
!
! Namelists for information in 5. {Not required at this release}
!
!
!  6. Lists of coupling nodes defining coupling frequencies between
!     internal models and between submodel partitions. (Not defined
!     yet at this release).
!CALL CNODE
!
!  7. Variables dealing with general coupling switches at the control
!     level. {These will require revision at the next release when
!     coupling between internal models is dealt with more generally.
!     Logicals below are set in routine SETGRCTL.}

      LOGICAL
     * new_im   ! new internal model next group of timesteps if .true.
     *,new_sm   ! new submodel dump  next group of timesteps if .true.

      COMMON/CSUBMGRP/new_im,new_sm

      INTEGER SUBMODEL_IDENT
      COMMON/SUBMODID/SUBMODEL_IDENT                                    
C
C
C An alternative common block required by TYPD1
C     COMDECK CALTSUBM:
C     COMDECK TYPD1 needs access to N_SUBMODEL_PARTITION/_MAX
C     in CSUBMODL. However, they are not always called in the same
C     decks and in the right order. Therefore, copy the values to
C     another comdeck and *CALL it from TYPD1

      INTEGER ALT_N_SUBMODEL_PARTITION
      INTEGER ALT_N_SUBMODEL_PARTITION_MAX

      PARAMETER(ALT_N_SUBMODEL_PARTITION_MAX=4)

      COMMON/CALTSUBM/ALT_N_SUBMODEL_PARTITION
! Subroutine arguments
!   Scalar arguments with intent(in):

!   Array  arguments with intent(in):

!   Scalar arguments with intent(InOut):

!   Array  arguments with intent(InOut):

!   Scalar arguments with intent(out):

!   Array  arguments with intent(out):

!   ErrorStatus
      INTEGER      ErrorStatus          ! Error flag (0 = OK)

! Local parameters:

! Local scalars:
      INTEGER
     * s               ! submodel loop
     *,i               ! internal model loop
     *,sm              ! submodel identifier
     *,im              ! internal model identifier
     *,sm_prev         ! previous submodel identifier
     *,im_prev         ! previous internal model identifier

! Local dynamic arrays:

! Function & Subroutine calls: None

!- End of header
!
! 1. Initialise lists before obtaining values for this experiment.
!
      do i=1,N_INTERNAL_MODEL_MAX
         INTERNAL_MODEL_LIST(i)      = 0
         SUBMODEL_FOR_IM(i)          = 0
      enddo   ! i over internal model list

      do im=1,INTERNAL_ID_MAX
         SUBMODEL_PARTITION_INDEX(im) = 0
         INTERNAL_MODEL_INDEX(im) = 0
         LAST_IM_IN_SM(im)=.false.
      enddo   ! im over internal model ids

      do s=1,N_SUBMODEL_PARTITION_MAX
         SUBMODEL_PARTITION_LIST(s)= 0
         SUBMODEL_FOR_SM(s)=0
      enddo  ! s over submodel list

      do sm=1,SUBMODEL_ID_MAX
         N_INTERNAL_FOR_SM(sm)      = 0
      enddo  ! sm over submodel ids

!
! 2. Obtain internal model and submodel identifiers from umui
!    generated namelist.
!
      read(5,NSUBMODL)
!
!
! 3. Check umui supplied values.
!
!
! 3.1 Check for umui supplied dimensions against parameter maxima.

      if(N_INTERNAL_MODEL.gt.N_INTERNAL_MODEL_MAX) then
         write(6,*) 'UM_Submodel_In: FATAL ERROR. Too many internal ',
     *   'models =',N_INTERNAL_MODEL,
     *   ' :You need to increase N_INTERNAL_MODEL_MAX'
         ErrorStatus=1       ! Set error flag
      endif
!
! 3.2 Check umui suppiled values are valid
!
      do i=1,N_INTERNAL_MODEL ! loop over internal models

        im = INTERNAL_MODEL_LIST(i) ! internal model identifier
        if(im.le.0.or.im.gt.INTERNAL_ID_MAX) then
         write(6,*) 'UM_Submodel_In: FATAL ERROR. Illegal internal ',
     *   'model identifier=',im,
     *   ' :Check values in namelist NSUBMODL supplied by umui'
         ErrorStatus=1       ! Set error flag
        endif

        sm = SUBMODEL_FOR_IM(i)     ! submodel for this internal model
        if(sm.le.0.or.sm.gt.SUBMODEL_ID_MAX) then
         write(6,*) 'UM_Submodel_In: FATAL ERROR. Illegal submodel ',
     *   'dump identifier=',sm,
     *   ' :Check values in namelist NSUBMODL supplied by umui'
         ErrorStatus=1       ! Set error flag
        endif

      enddo ! i=1,N_INTERNAL_MODEL
!
! 4. Form internal model and submodel description arrays.
!
      sm_prev = 0             ! Null value of submodel identifier
      N_SUBMODEL_PARTITION=0  ! Count no. of submodel partitions

      do i=1,N_INTERNAL_MODEL ! loop over internal models

        im = INTERNAL_MODEL_LIST(i) ! internal model identifier
        sm = SUBMODEL_FOR_IM(i)     ! submodel for this internal model
        INTERNAL_MODEL_INDEX(im)=i  ! sequence no. for STASH arrays

        if(sm.ne.sm_prev) then  ! new submodel

           N_SUBMODEL_PARTITION = N_SUBMODEL_PARTITION+1
           SUBMODEL_PARTITION_LIST(N_SUBMODEL_PARTITION) = sm

!   Since this is a new submodel, the previous internal model must be
!   the last internal model in its submodel partition.
           IF(N_SUBMODEL_PARTITION.GT.1) THEN ! Not first dump
              LAST_IM_IN_SM(im_prev) = .true.
           ENDIF

        endif                   ! test on new submodel
        SUBMODEL_FOR_SM(IM) = N_SUBMODEL_PARTITION

        SUBMODEL_PARTITION_INDEX(im)=sm
        N_INTERNAL_FOR_SM(sm)=N_INTERNAL_FOR_SM(sm)+1

        im_prev=im
        sm_prev=sm

      enddo ! i=1,N_INTERNAL_MODEL

      LAST_IM_IN_SM(im) = .true.  ! last im in list is last im in sm

!
! 5. Check calculated dimensions against parameter maxima.

      if(N_SUBMODEL_PARTITION.gt.N_SUBMODEL_PARTITION_MAX) then
         write(6,*) 'UM_Submodel_In: FATAL ERROR. Too many submodels =',
     *   N_SUBMODEL_PARTITION,
     *   ' You need to increase N_SUBMODEL_PARTITION_MAX'
         ErrorStatus=1       ! Set error flag
      endif
!
C     Need a copy of No of submodels for use by TYPD1.
      ALT_N_SUBMODEL_PARTITION=N_SUBMODEL_PARTITION

      if (ALT_N_SUBMODEL_PARTITION_MAX.NE.N_SUBMODEL_PARTITION_MAX)THEN
        write(6,*)'UM_Submodel_In: Mismatch in parameters '
        WRITE(6,*)'N_SUBMODEL_PARTITION_MAX and '
        WRITE(6,*)'ALT_N_SUBMODEL_PARTITION_MAX. '
        WRITE(6,*)'They should be identical '
        ErrorStatus=1
        endif
      return
      end
