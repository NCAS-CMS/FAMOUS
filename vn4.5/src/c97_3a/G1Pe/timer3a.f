C ******************************COPYRIGHT******************************
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
C
CLL SUBROUTINE TIMER ------------------------------------------------
CLL
CLL                    Purpose:
CLL Allows the recording of time spent in any section of the program
CLL Two types of timings are supported:
CLL non-inclusive : if a timed section of code (1) contains another
CLL                 timed section of code (2), then the timing for
CLL                 section (1) will not include the time spent in
CLL                 section (2). This is the normal use for the timer
CLL                 routine in the UM up to vn3.4
CLL inclusive     : allows the user to measure the time taken between
CLL                 any two points in the code, irrespective of any
CLL                 other calls to the timer routine within the timed
CLL                 section
CLL
CLL NB: Non-inclusive timers DO INCLUDE any inclusive timer sections
CLL     contained within them. If this section of code should not be
CLL     included, then also time it with a non-inclusive timer
CLL
CLL Timer now also records the time spent in itself
CLL Parameters:
CLL section_name - 20 byte character string containing name of
CLL                timer reference
CLL
CLL action:
CLL  1 -> first call to timer (timer initialisation)
CLL  2 -> last call to timer (prints out the collected data)
CLL  3 -> non-inclusive start timer
CLL  4 -> non-inclusive end timer
CLL  5 -> inclusive start timer
CLL  6 -> inclusive end timer
CLL
CLL Timer should be called with action=1 before the first executable
CLL statement, and with action=2 after the last executable statement.
CLL
CLL   Model               Modification History
CLL  version    Date
CLL   4.1       16/08/94  Based on old UM timer routine
!LL   4.2       08/10/96  Corrected intermediate timer error in
!LL                       elapsed times.
!LL                       Corrected size of message arg. in
!LL                       TIMER_OUTPUT.
!LL                       P.Burton
!LL   4.3       23/01/97  Added better overview for MPP runs
!LL                       Corrected T3E wallclock time calculation
!LL                       P.Burton
!LL   4.5       17/04/98  Added barrier to allow imbalance to be
!LL                       included in the correct timer section
!LL                                                     P.Burton
!LL   4.5       09/07/98  Replaced missing array index for
!LL                       last_ni_wallclock_time_elapsed.
!LL                                                     P. Burton
CLL
CLL  Author : Paul Burton
CLL

      SUBROUTINE TIMER(section_name,action)

      IMPLICIT NONE

! Arguments:
      CHARACTER*(*) section_name  ! IN reference name for timed section
      INTEGER action            ! IN what action to take

! Local variables:
! ni prefix = non-inclusive timings
! in prefix = inclusive timings


      INTEGER max_timers      ! maximum number of timings to be handled
         PARAMETER ( max_timers=300 )

      CHARACTER*20 ni_timer_name(max_timers),
     &                        ! names of timer references
     &             in_timer_name(max_timers)
     &                        ! names of timer references

      INTEGER ni_number_of_times_timed(max_timers),
     &                        ! number of times that a section of code
     &                        ! has been timed
     &        in_number_of_times_timed(max_timers)
     &                        ! number of times that a section of code
     &                        ! has been timed

      INTEGER old_timer(max_timers)
     &                        ! the reference of the timer stopped when
     &                        ! a new one is started

      REAL ni_cpu_time_elapsed(max_timers),
     &                        ! total amount of cpu time spent in
     &                        ! a section of code
     &     ni_wallclock_time_elapsed(max_timers),
     &                        ! total amount of wallclock time
     &     in_cpu_time_elapsed(max_timers),
     &                        ! total amount of time spent in a section
     &                        ! of code
     &     in_wallclock_time_elapsed(max_timers)
     &                        ! total amount of wallclock time

      REAL ni_cpu_time_started,
     &                        ! for non-inclusive timer - cpu time
     &                        ! of starting
     &     ni_wallclock_time_started,
     &                        ! wallclock time of starting
     &     in_cpu_time_started(max_timers),
     &                        ! for inclusive timer - cpu time
     &                        !of starting
     &     in_wallclock_time_started(max_timers)
     &                        ! wallclock time of starting

      INTEGER current_timer   ! for non-inclusive timer - current
     &                        ! section of code being timed
      INTEGER ni_number_of_timers,
     &                        ! number of timers currently known about
     &        in_number_of_timers
     &                        ! number of timers currently known about

      LOGICAL in_timer_running(max_timers)
     &                        ! is a particular timer running?

      INTEGER section_ref     ! reference of the current section



      REAL cpu_time_into_timer,
     &                        ! cpu time at which timer routine entered
     &     wallclock_time_into_timer ! wallclock "

      LOGICAL timer_on        ! set to FALSE if an error occurs

! Saved variables:
      SAVE ni_timer_name,in_timer_name,
     &     ni_number_of_times_timed,in_number_of_times_timed,
     &     ni_cpu_time_elapsed,in_cpu_time_elapsed,
     &     ni_wallclock_time_elapsed,in_wallclock_time_elapsed,
     &     ni_cpu_time_started,in_cpu_time_started,
     &     ni_wallclock_time_started,in_wallclock_time_started,
     &     current_timer ,
     &     ni_number_of_timers, in_number_of_timers,
     &     in_timer_running,
     &     old_timer, timer_on

! Magic numbers (action types):
      INTEGER first_call_to_timer,
     &        last_call_to_timer,
     &        non_incl_start_timer,
     &        non_incl_end_timer,
     &        incl_start_timer,
     &        incl_end_timer,
     &        intermediate_output

      PARAMETER ( first_call_to_timer = 1,
     &            last_call_to_timer = 2,
     &            non_incl_start_timer = 3,
     &            non_incl_end_timer = 4,
     &            incl_start_timer = 5,
     &            incl_end_timer = 6,
     &            intermediate_output = 7 )

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
      INTEGER info
! Loop counters etc.
      INTEGER I



      EXTERNAL get_cpu_time,get_wallclock_time
      REAL get_cpu_time,get_wallclock_time
! ----------------------------------------------------------------------

      IF (( action .EQ. first_call_to_timer) .OR.
     &    ( action .EQ. non_incl_start_timer) .OR.
     &    ( action .EQ. incl_start_timer)) THEN
        CALL GC_GSYNC(nproc,info)
      ENDIF
      IF (action .GT. 100) action=action-100
! The following line is useful for general debugging purposes
! It prints out the name of every timed routine on entry and exit
!         WRITE(6,*) section_name,' action= ',action

! start up the timer timer

      cpu_time_into_timer = get_cpu_time()
      wallclock_time_into_timer = get_wallclock_time()
      in_number_of_times_timed(1) = in_number_of_times_timed(1) + 1

! check the length of the section_name

      IF (LEN(section_name) .GT. 20) THEN
        WRITE(6,*) 'TIMER has detected a non-fatal ERROR'
        WRITE(6,*) 'Section name ',section_name,' is too long.'
        WRITE(6,*) 'Maximum of 20 characters is allowed'
        WRITE(6,*) section_name,' will be truncated to 20 chars.'
      ENDIF


! diagnose what action to take:

      IF (action .EQ. first_call_to_timer) THEN

! First call to timer - do initialisation

        DO I=1,max_timers
          ni_timer_name(I)            = '                    '
          in_timer_name(I)            = '                    '

          ni_number_of_times_timed(I) = 0
          in_number_of_times_timed(I) = 0

          ni_cpu_time_elapsed(I)      = 0.
          in_cpu_time_elapsed(I)      = 0.
          ni_wallclock_time_elapsed(I)= 0.
          in_wallclock_time_elapsed(I)= 0.

          in_timer_running(I)=.FALSE.
        ENDDO

        timer_on = .TRUE.

        current_timer = 1
        ni_number_of_timers = 1
        in_number_of_timers = 1
        in_timer_name(1) = 'TIMER'

! and start the timer running

        ni_cpu_time_started = get_cpu_time()
        ni_wallclock_time_started = get_wallclock_time()
        ni_number_of_times_timed(current_timer) = 1
        old_timer(current_timer) = 0
        ni_timer_name(current_timer) = section_name

! ----------------------------------------------------------------------

      ELSEIF (timer_on .AND.
     &        ( (action .EQ. last_call_to_timer) .OR.
     &          (action .EQ. intermediate_output) ) )THEN

! Last call to timer - or intermediate output required, so
! print out table of results

      IF (action .EQ. last_call_to_timer) THEN
! the only active timer should be no.1

        IF (current_timer .NE. 1) THEN
          WRITE(6,*) 'TIMER has detected an ERROR'
          WRITE(6,*) 'Attempted to print results without switching ',
     &             'off all running non-inclusive timers.'
          WRITE(6,*) '** TIMER SWITCHED OFF FOR REST OF RUN **'
          timer_on = .FALSE.
          GOTO 9999
        ENDIF

! Make sure there are no inclusive timers still running

        section_ref = 0
        DO I=1,in_number_of_timers
          IF (in_timer_running(I)) section_ref = I
        ENDDO

        IF (section_ref .NE.0) THEN
          WRITE(6,*) 'TIMER has detected an ERROR'
          WRITE(6,*) 'Attempted to print results without switching ',
     &             'off all running inclusive timers.'
          WRITE(6,*) '** TIMER SWITCHED OFF FOR REST OF RUN **'
          timer_on = .FALSE.
          GOTO 9999
        ENDIF

! Just to make sure that timer isn't called again
        timer_on = .FALSE.

! and switch off the top level non-inclusive timer

        ni_cpu_time_elapsed(current_timer) =
     &      ni_cpu_time_elapsed(current_timer) +
     &      get_cpu_time() - ni_cpu_time_started

        ni_wallclock_time_elapsed(current_timer) =
     &      ni_wallclock_time_elapsed(current_timer) +
     &      get_wallclock_time() - ni_wallclock_time_started

      ENDIF ! If this is the final call to timer

      CALL TIMER_OUTPUT(
     &  in_number_of_timers, ni_number_of_timers,
     &  in_cpu_time_elapsed, ni_cpu_time_elapsed,
     &  in_wallclock_time_elapsed, ni_wallclock_time_elapsed,
     &  in_number_of_times_timed, ni_number_of_times_timed,
     &  in_timer_name, ni_timer_name,
     &  action,section_name)



! that's it

! ----------------------------------------------------------------------

      ELSEIF (timer_on .AND.
     &        (action .EQ. non_incl_start_timer) ) THEN

! Start a non-inclusive timer running

! Switch off the current timer
      ni_cpu_time_elapsed(current_timer) =
     &    ni_cpu_time_elapsed(current_timer) +
     &    get_cpu_time() - ni_cpu_time_started

      ni_wallclock_time_elapsed(current_timer) =
     &    ni_wallclock_time_elapsed(current_timer) +
     &    get_wallclock_time() - ni_wallclock_time_started

! See if we're already keeping records for this section

      section_ref = 0
      DO I=1,ni_number_of_timers
        IF (ni_timer_name(I) .EQ. section_name) section_ref = I
      ENDDO

! Check to make sure that there is no timer already running for
! this section
! (NB an inclusive timer with the same reference name is allowed
!  to run simultaneously with this one)

      IF (section_ref .EQ. current_timer) THEN
        WRITE(6,*) 'TIMER has detected an ERROR'
        WRITE(6,*) 'Simultaneous non-inclusive timers attempted ',
     &             'for section ',section_name
        WRITE(6,*) '** TIMER SWITCHED OFF FOR REST OF RUN **'
        timer_on = .FALSE.
        GOTO 9999
      ENDIF

! calculate the section reference for the new timer

      IF (section_ref .EQ. 0) THEN
!       this is a new section
        section_ref = ni_number_of_timers+1
        ni_timer_name(section_ref) = section_name
        ni_number_of_timers = section_ref
      ENDIF

! check that max_timers isn't exceeded:
      IF (ni_number_of_timers .GT. max_timers) THEN
        WRITE(6,*) 'TIMER has detected an ERROR'
        WRITE(6,*) 'More than ',max_timers,' non-inclusive ',
     &             'timers is not allowed.'
        WRITE(6,*) '** TIMER SWITCHED OFF FOR REST OF RUN **'
        timer_on = .FALSE.
        GOTO 9999
      ENDIF

! set up old_timer so that when this new timer is stopped, the
! current timer (that we've just stopped) can be restarted

      old_timer(section_ref)=current_timer

! now start up the new timer

      current_timer = section_ref
      ni_number_of_times_timed(current_timer) =
     &  ni_number_of_times_timed(current_timer) + 1
      ni_cpu_time_started = get_cpu_time()
      ni_wallclock_time_started = get_wallclock_time()

! that's it

! ----------------------------------------------------------------------

      ELSEIF (timer_on .AND.
     &        (action .EQ. non_incl_end_timer) ) THEN

! Stop a non-inclusive timer

! Make sure that we're being asked to end a timer that's actually
! running.

      IF (ni_timer_name(current_timer) .NE. section_name) THEN
        WRITE(6,*) 'TIMER has detected an ERROR'
        WRITE(6,*) 'Attempted to stop a non-active ',
     &             'non-inclusive timer ',section_name
        WRITE(6,*) '** TIMER SWITCHED OFF FOR REST OF RUN **'
        timer_on = .FALSE.
        GOTO 9999
      ENDIF

! OK - so stop this timer:

      ni_cpu_time_elapsed(current_timer) =
     &    ni_cpu_time_elapsed(current_timer) +
     &    get_cpu_time() - ni_cpu_time_started

      ni_wallclock_time_elapsed(current_timer) =
     &    ni_wallclock_time_elapsed(current_timer) +
     &    get_wallclock_time() - ni_wallclock_time_started

! and now restart the old timer (ie. the one that was in
! operation at the time this one was started)

      IF (old_timer(current_timer) .EQ. 0) THEN
! this means I have just stopped the top level timer - there
! are no more to stop. This is an error - I should do this
! by calling the timer with action=2
         WRITE(6,*) 'TIMER has detected an ERROR'
         WRITE(6,*) 'The top-level timer has been stopped'
         WRITE(6,*) '** TIMER SWITCHED OFF FOR REST OF RUN **'
        timer_on = .FALSE.
        GOTO 9999
      ENDIF

      current_timer=old_timer(current_timer)
      ni_cpu_time_started=get_cpu_time()
      ni_wallclock_time_started=get_wallclock_time()

! ----------------------------------------------------------------------

      ELSEIF (timer_on .AND.
     &        (action .EQ. incl_start_timer) ) THEN

! Start an inclusive timer running

! See if we're already keeping records for this section

      section_ref = 0
      DO I=1,in_number_of_timers
        IF (in_timer_name(I) .EQ. section_name) section_ref = I
      ENDDO

! and calculate the section reference

      IF (section_ref .EQ. 0) THEN
!       this is a new one
        section_ref = in_number_of_timers + 1
        in_timer_name(section_ref) = section_name
        in_number_of_timers = section_ref
      ENDIF

! check that max_timers isn't exceeded:

      IF (in_number_of_timers .GT. max_timers) THEN
        WRITE(6,*) 'TIMER has detected an ERROR'
        WRITE(6,*) 'More than ',max_timers,' inclusive ',
     &             'timers is not allowed.'
        WRITE(6,*) '** TIMER SWITCHED OFF FOR REST OF RUN **'
        timer_on = .FALSE.
        GOTO 9999
      ENDIF

! Check to make sure that there is no timer already running for
! this section
! (NB a non-inclusive timer with the same reference name is allowed
!  to run simultaneously with this one)

      IF (in_timer_running(section_ref)) THEN
        WRITE(6,*) 'TIMER has detected an ERROR'
        WRITE(6,*) 'Inclusive timer already running for ',
     &             section_name
        WRITE(6,*) '** TIMER SWITCHED OFF FOR REST OF RUN **'
        timer_on = .FALSE.
        GOTO 9999
      ENDIF

! so now we can start the timer for this section
      in_number_of_times_timed(section_ref) =
     &  in_number_of_times_timed(section_ref) + 1
      in_timer_running(section_ref) = .TRUE.
      in_cpu_time_started(section_ref) = get_cpu_time()
      in_wallclock_time_started(section_ref) = get_wallclock_time()

! that's it


! ----------------------------------------------------------------------

      ELSEIF (timer_on .AND.
     &        (action .EQ. incl_end_timer) ) THEN

! Stop an inclusive timer

! Find out what the reference number of this timer is

      section_ref = 0
      DO I=1,in_number_of_timers
        IF (in_timer_name(I) .EQ. section_name) section_ref = I
      ENDDO

      IF (section_ref .EQ. 0) THEN
        WRITE(6,*) 'TIMER has detected an ERROR'
        WRITE(6,*) 'Attempting to stop a non-existent ',
     &             'inclusive timer ',section_name
        WRITE(6,*) '** TIMER SWITCHED OFF FOR REST OF RUN **'
        timer_on = .FALSE.
        GOTO 9999
      ENDIF

! Make sure this timer is actually running at the moment

      IF (.NOT. in_timer_running(section_ref)) THEN
        WRITE(6,*) 'TIMER has detected an ERROR'
        WRITE(6,*) 'Attempting to stop a non-running ',
     &             'inclusive timer ',section_name
        WRITE(6,*) '** TIMER SWITCHED OFF FOR REST OF RUN **'
        timer_on = .FALSE.
        GOTO 9999
      ENDIF

! now we can stop it
      in_cpu_time_elapsed(section_ref) =
     &   in_cpu_time_elapsed(section_ref) +
     &   get_cpu_time() - in_cpu_time_started(section_ref)

      in_wallclock_time_elapsed(section_ref) =
     &   in_wallclock_time_elapsed(section_ref) +
     &   get_wallclock_time() - in_wallclock_time_started(section_ref)

      in_timer_running(section_ref) = .FALSE.

! that's it

! ----------------------------------------------------------------------

      ELSEIF (timer_on) THEN

        WRITE(6,*) 'TIMER has detected an ERROR'
        WRITE(6,*) 'Incorrect action= ',action,' supplied by ',
     &             'section ',section_name
        WRITE(6,*) 'Non-fatal error - TIMER will continue'

      ENDIF

 9999 CONTINUE

! stop the timer timer
      in_cpu_time_elapsed(1) = in_cpu_time_elapsed(1) +
     &     get_cpu_time() - cpu_time_into_timer
      in_wallclock_time_elapsed(1) = in_wallclock_time_elapsed(1) +
     &     get_wallclock_time() - wallclock_time_into_timer

      RETURN
      END

!*********************************************************************


      SUBROUTINE TIMER_OUTPUT(
     &  in_number_of_timers, ni_number_of_timers,
     &  in_cpu_time_elapsed, ni_cpu_time_elapsed,
     &  in_wallclock_time_elapsed, ni_wallclock_time_elapsed,
     &  in_number_of_times_timed, ni_number_of_times_timed,
     &  in_timer_name, ni_timer_name,
     &  action,message)

      IMPLICIT NONE

! Arguments

      INTEGER
     &  in_number_of_timers  ! IN number of inclusive timers
     &, ni_number_of_timers  ! IN number of non-inclusive timers
     &, in_number_of_times_timed(in_number_of_timers)
!                            ! IN number of times timed - inclusive
     &, ni_number_of_times_timed(ni_number_of_timers)
!                            ! IN number of times timed - non-incl.
     &, action  ! final output or intermediate

      REAL
     &  in_cpu_time_elapsed(in_number_of_timers)
!                            ! IN elapsed inclusive CPU time
     &, ni_cpu_time_elapsed(ni_number_of_timers)
!                            ! IN elapsed non-inclusive CPU time
     &, in_wallclock_time_elapsed(in_number_of_timers)
!                            ! IN elapsed inclusive wallclock time
     &, ni_wallclock_time_elapsed(ni_number_of_timers)

      CHARACTER*20
     &  in_timer_name(in_number_of_timers)
!                            ! IN name of timed section - inclusive
     &, ni_timer_name(ni_number_of_timers)
!                            ! IN name of timed section - non-incl.

      CHARACTER*(*)
     &   message              ! IN message to print


! Local variables

      INTEGER max_timers
      PARAMETER(max_timers=300)

      INTEGER last_call_to_timer,intermediate_output
      PARAMETER(last_call_to_timer=2,
     &          intermediate_output=7)

      INTEGER
     &  number_of_timers
     &, local_number_of_times_timed(max_timers)

      REAL
     &  local_cpu_time_elapsed(max_timers)
     &, local_wallclock_time_elapsed(max_timers)

      CHARACTER*20
     &  local_timer_name(max_timers)

! Variables required for using intermediate timers
! They record the values on the last call to this routine
      INTEGER
     &  last_in_number_of_times_timed(max_timers)
     &, last_ni_number_of_times_timed(max_timers)

      REAL
     &  last_in_cpu_time_elapsed(max_timers)
     &, last_ni_cpu_time_elapsed(max_timers)
     &, last_in_wallclock_time_elapsed(max_timers)
     &, last_ni_wallclock_time_elapsed(max_timers)

      LOGICAL
     &  first_intermediate_timer_call

      DATA first_intermediate_timer_call /.TRUE./
      DATA last_in_number_of_times_timed /max_timers*0/
      DATA last_ni_number_of_times_timed /max_timers*0/
      DATA last_in_cpu_time_elapsed /max_timers*0.0/
      DATA last_ni_cpu_time_elapsed /max_timers*0.0/
      DATA last_in_wallclock_time_elapsed /max_timers*0.0/
      DATA last_ni_wallclock_time_elapsed /max_timers*0.0/

      SAVE
     &  last_in_number_of_times_timed, last_ni_number_of_times_timed,
     &  last_in_cpu_time_elapsed, last_ni_cpu_time_elapsed,
     &  last_in_wallclock_time_elapsed, last_ni_wallclock_time_elapsed,
     &  first_intermediate_timer_call



      INTEGER sortwork_int    ! work variable for sort
      REAL    sortwork_real   ! work variable for sort
      CHARACTER*20 sortwork_char ! work variable for sort

      REAL total_cpu_time,         ! total cpu time spent in program
     &     total_wallclock_time,   ! total wallclock time spent in
!                                  ! program
     &     average_cpu_elapsed,    ! average cpu elapsed time
     &     average_wallclock_elapsed, ! average wallclock elapsed time
     &     percent_of_cpu_total,   ! % of cpu time spent in a section
     &     percent_of_wallclock_total,
     &                            ! % of wallclock time spent in a
!                                 ! section
     &     speed_up               ! speed_up=cpu/wallclock

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

! These are the declarations for MPP timer

      INTEGER info,
     &  wallclock_max_pe(max_timers),wallclock_min_pe(max_timers),
     &  cpu_max_pe(max_timers),cpu_min_pe(max_timers)

      REAL wallclock_mean(max_timers),cpu_mean(max_timers),
     &     wallclock_median(max_timers),cpu_median(max_timers),
     &     wallclock_sd(max_timers),cpu_sd(max_timers),
     &     wallclock_max(max_timers),wallclock_min(max_timers),
     &     cpu_max(max_timers),cpu_min(max_timers),
     &     cpu_total(max_timers),speedup(max_timers),
     &     efficiency(max_timers)

      INTEGER
     &  summ_n_timers    ! number of routines for ni summary
     &, routine_id          ! routine id on this processor

      REAL
     &  wallclock_times(0:MAXPROC)  ! wallclock time from each proc
     &, cpu_times(0:MAXPROC)        ! cpu time from each proc
     &, total_cpu,max_wall       ! total cpu, maxumum wallclock times

      CHARACTER*20 summ_section(max_timers)  ! names of sections
      COMMON /MPP_TIMER/ summ_n_timers,total_cpu,max_wall,
     &                   wallclock_times,cpu_times,
     &                   summ_section

! Variables for loops etc.
      INTEGER I,J,K,timer_kind

! Check to see if this is an intermediate output, and the first
! time it has been called
      IF ((action .EQ. intermediate_output) .AND.
     &    (first_intermediate_timer_call  ) ) THEN
! Copy the arguments into the last_* arrays
        first_intermediate_timer_call=.FALSE.

        DO I=1,in_number_of_timers
          last_in_number_of_times_timed(I)=in_number_of_times_timed(I)
          last_in_cpu_time_elapsed(I)=in_cpu_time_elapsed(I)
          last_in_wallclock_time_elapsed(I)=
     &      in_wallclock_time_elapsed(I)
        ENDDO

        DO I=1,ni_number_of_timers
          last_ni_number_of_times_timed(I)=ni_number_of_times_timed(I)
          last_ni_cpu_time_elapsed(I)=ni_cpu_time_elapsed(I)
          last_ni_wallclock_time_elapsed(I)=
     &      ni_wallclock_time_elapsed(I)
        ENDDO

        GOTO 9999  ! jump to end - no output on first call
      ENDIF

      WRITE(6,*)
      WRITE(6,*) '******************************************'
      WRITE(6,*)

      DO timer_kind=1,2  ! 1 is non-inclusive and 2 is inclusive
! Copy arguments into local arrays
        IF (action .EQ. last_call_to_timer) THEN
          WRITE(6,*) 'END OF RUN - TIMER OUTPUT'
          WRITE(6,*) 'Timer information is for whole run'
          IF (timer_kind .EQ. 1) THEN  ! non-inclusive timer
            number_of_timers=ni_number_of_timers
            DO I=1,number_of_timers
              local_timer_name(I)=ni_timer_name(I)
              local_cpu_time_elapsed(I)=ni_cpu_time_elapsed(I)
              local_wallclock_time_elapsed(I)=
     &          ni_wallclock_time_elapsed(I)
              local_number_of_times_timed(I)=
     &          ni_number_of_times_timed(I)
            ENDDO
          ELSE ! timer_kind .EQ. 2 - inclusive timer
            number_of_timers=in_number_of_timers
            DO I=1,number_of_timers
              local_timer_name(I)=in_timer_name(I)
              local_cpu_time_elapsed(I)=in_cpu_time_elapsed(I)
              local_wallclock_time_elapsed(I)=
     &          in_wallclock_time_elapsed(I)
              local_number_of_times_timed(I)=
     &          in_number_of_times_timed(I)
            ENDDO
          ENDIF ! which timer kind this was
        ELSE  ! this is an intermediate output call
          WRITE(6,*) 'INTERMEDIATE TIMER OUTPUT :',message
          WRITE(6,*) 'Timer information is only for code executed ',
     &               'since last intermediate timer output.'
          IF (timer_kind .EQ. 1) THEN  ! non-inclusive timer
            number_of_timers=ni_number_of_timers
            DO I=1,number_of_timers
              local_timer_name(I)=ni_timer_name(I)
              local_cpu_time_elapsed(I)=ni_cpu_time_elapsed(I)-
     &                                  last_ni_cpu_time_elapsed(I)
              local_wallclock_time_elapsed(I)=
     &          ni_wallclock_time_elapsed(I)-
     &          last_ni_wallclock_time_elapsed(I)
              local_number_of_times_timed(I)=
     &          ni_number_of_times_timed(I)-
     &          last_ni_number_of_times_timed(I)
            ENDDO
          ELSE ! timer kind .EQ. 2 - inclusive timer
            number_of_timers=in_number_of_timers
            DO I=1,number_of_timers
              local_timer_name(I)=in_timer_name(I)
              local_cpu_time_elapsed(I)=in_cpu_time_elapsed(I)-
     &                                  last_in_cpu_time_elapsed(I)
              local_wallclock_time_elapsed(I)=
     &          in_wallclock_time_elapsed(I)-
     &          last_in_wallclock_time_elapsed(I)
              local_number_of_times_timed(I)=
     &          in_number_of_times_timed(I)-
     &          last_in_number_of_times_timed(I)
            ENDDO
          ENDIF  ! what timer type
        ENDIF  ! what action to perform

! Do work for non-inclusive timers

! Calculate the total time in the program (based on non-inclusive
! timers)
        IF (timer_kind .EQ. 1) THEN
          total_cpu_time = 0.0
          total_wallclock_time = 0.0
          DO I=1,number_of_timers
          total_cpu_time = total_cpu_time + local_cpu_time_elapsed(I)
          total_wallclock_time =
     &      total_wallclock_time + local_wallclock_time_elapsed(I)
          ENDDO

          WRITE(6,*) 'PE ',mype,' Elapsed CPU Time: ',
     &               total_cpu_time
          WRITE(6,*) 'PE ',mype,'  Elapsed Wallclock Time: ',
     &                total_wallclock_time

! Calculate the total cpu time over all processors and the
! maximum elapsed time - so allowing a speedup to be caclulated

          total_cpu=total_cpu_time
          max_wall=total_wallclock_time

          CALL GC_RSUM(1,nproc,info,total_cpu)
          CALL GC_RMAX(1,nproc,info,max_wall)

          max_wall=MAX(max_wall,0.000001)
          WRITE(6,*)
          WRITE(6,*) 'Total Elapsed CPU Time: ',
     &               total_cpu
          WRITE(6,*) 'Maximum Elapsed Wallclock Time: ',
     &               max_wall
          WRITE(6,*) 'Speedup: ',total_cpu/max_wall
          WRITE(6,*) '--------------------------------------------'

        ENDIF

! Sort subroutines into time order (based on wallclock time)

        DO I=1,number_of_timers-1
          DO J=(I+1),number_of_timers
            IF (local_wallclock_time_elapsed(J) .GT.
     &          local_wallclock_time_elapsed(I)) THEN

!             Swap the two entries
              sortwork_real = local_cpu_time_elapsed(I)
              local_cpu_time_elapsed(I) = local_cpu_time_elapsed(J)
              local_cpu_time_elapsed(J) = sortwork_real

              sortwork_real = local_wallclock_time_elapsed(I)
              local_wallclock_time_elapsed(I) =
     &          local_wallclock_time_elapsed(J)
              local_wallclock_time_elapsed(J) = sortwork_real

              sortwork_int = local_number_of_times_timed(I)
              local_number_of_times_timed(I) =
     &          local_number_of_times_timed(J)
              local_number_of_times_timed(J) = sortwork_int

              sortwork_char = local_timer_name(I)
              local_timer_name(I) = local_timer_name(J)
              local_timer_name(J) = sortwork_char

            ENDIF
          ENDDO
        ENDDO

 20     FORMAT(20X, A45,I4)
 21     FORMAT(3X,'ROUTINE',6X,'CALLS',2X,'TOT CPU',4X,
     &       'AVERAGE',3X,'TOT WALL',2X,'AVERAGE',2X,
     &       '% CPU',4X,'% WALL',4X,'SPEED-UP')
 22     FORMAT(3X,'ROUTINE',6X,'CALLS',2X,'TOT CPU',4X,
     &         'AVERAGE',3X,'TOT WALL',2X,'AVERAGE',2X,'SPEED-UP')
 23     FORMAT(/,I3,1X,A12,1X,I4,4(2X,F8.2),2(2X,F6.2),4X,F6.2)
        IF (timer_kind .EQ. 1) THEN
          WRITE(6,20) 'Non-Inclusive Timer Summary for PE ',mype
          WRITE(6,21)
        ELSE
          WRITE(6,20) 'Inclusive Timer Summary for PE ',mype
          WRITE(6,22)
        ENDIF

        DO I=1,number_of_timers
          IF (local_number_of_times_timed(I) .NE. 0) THEN
            average_cpu_elapsed =  local_cpu_time_elapsed(I)/
     &                             local_number_of_times_timed(I)
            average_wallclock_elapsed = local_wallclock_time_elapsed(I)/
     &                                  local_number_of_times_timed(I)
          ELSE
             average_cpu_elapsed = 0.0
             average_wallclock_elapsed = 0.0
          ENDIF

          IF (local_wallclock_time_elapsed(I) .GT. 0) THEN
            speed_up=local_cpu_time_elapsed(I)/
     &               local_wallclock_time_elapsed(I)
          ELSE
            speed_up=1.0
          ENDIF

          IF (timer_kind .EQ. 1) THEN  ! non-inclusive timer has some
!                                      ! extra output

            percent_of_cpu_total = 100.0*local_cpu_time_elapsed(I)/
     &                             total_cpu_time
            percent_of_wallclock_total =
     &        100.0*local_wallclock_time_elapsed(I)/
     &        total_wallclock_time


            WRITE(6,23) I,local_timer_name(I),
     &                  local_number_of_times_timed(I),
     &                  local_cpu_time_elapsed(I),average_cpu_elapsed,
     &                  local_wallclock_time_elapsed(I),
     &                  average_wallclock_elapsed,
     &                  percent_of_cpu_total,
     &                  percent_of_wallclock_total,speed_up

          ELSE ! inclusive timer has slightly less to output

            WRITE(6,23) I,local_timer_name(I),
     &                  local_number_of_times_timed(I),
     &                  local_cpu_time_elapsed(I),average_cpu_elapsed,
     &                  local_wallclock_time_elapsed(I),
     &                  average_wallclock_elapsed,speed_up

          ENDIF

        ENDDO



! And now to assemble an overall timing assesment on PE0
! Each PE sends it total wallclock and cpu time spent in each routine
! to PE0, which calculates the average, s.d., max and min, and
! sorts on the basis of the average wallclock time
!
!
! We'll use the list of routines that PE0 already has as the master
! list.

        IF (mype .EQ. 0) THEN
          WRITE(6,*)
          WRITE(6,*) 'MPP Timing information : '
          WRITE(6,*)  nproc,' processors in configuration ',nproc_x,
     &                ' x ',nproc_y

          summ_n_timers=number_of_timers
          DO I=1,summ_n_timers
            summ_section(I)=local_timer_name(I)
          ENDDO
        ENDIF

! tell everyone else how many routines to do summary on - and which
! routines they are
        CALL GC_IBCAST(3213,1,0,nproc,info,summ_n_timers)
        CALL GC_CBCAST(3214,20*summ_n_timers,0,nproc,info,
     &                 summ_section)


        DO I=1,summ_n_timers

! which section_ref is this for me?

          routine_id=0
          DO J=1,number_of_timers
            IF (local_timer_name(J) .EQ. summ_section(I))
     &        routine_id=J
          ENDDO

          IF (routine_id .GT. 0) THEN
            wallclock_times(mype)=
     &        local_wallclock_time_elapsed(routine_id)
            cpu_times(mype)=local_cpu_time_elapsed(routine_id)
          ELSE
            wallclock_times(mype)=0.0
            cpu_times(mype)=0.0
          ENDIF

! send my information to PE 0.
          CALL GC_GSYNC (nproc,info)

          CALL GC_RSEND(1000+mype,1,0,info,wallclock_times(mype),
     &                wallclock_times(mype))
          CALL GC_GSYNC(nproc,info)

          IF (mype .EQ. 0) THEN
            DO J=0,nproc-1
              CALL GC_RRECV(1000+J,1,J,info,wallclock_times(J),
     &                      wallclock_times(J))
            ENDDO
          ENDIF
          CALL GC_GSYNC(nproc,info)

          CALL GC_RSEND(10000+mype,1,0,info,cpu_times(mype),
     &                cpu_times(mype))
          CALL GC_GSYNC(nproc,info)

          IF (mype .EQ. 0) THEN
            DO J=0,nproc-1
              CALL GC_RRECV(10000+J,1,J,info,cpu_times(J),
     &                      cpu_times(J))
            ENDDO
          ENDIF

          IF (mype .EQ. 0) THEN
! collect all the information - and start calculating the statistics
            wallclock_mean(I)=0.0
            cpu_total(I)=0.0
            wallclock_max(I)=-1.0E30
            wallclock_min(I)=1.0E30
            cpu_max(I)=-1.0E30
            cpu_min(I)=1.0E30

            DO J=0,nproc-1

              wallclock_mean(I)=wallclock_mean(I)+wallclock_times(J)
              cpu_total(I)=cpu_total(I)+cpu_times(J)

              IF (wallclock_times(J).GT.wallclock_max(I)) THEN
                wallclock_max(I)=wallclock_times(J)
                wallclock_max_pe(I)=J
              ENDIF
              IF (wallclock_times(J).LT.wallclock_min(I)) THEN
                wallclock_min(I)=wallclock_times(J)
                wallclock_min_pe(I)=J
              ENDIF
              IF (cpu_times(J).GT.cpu_max(I)) THEN
                cpu_max(I)=cpu_times(J)
                cpu_max_pe(I)=J
              ENDIF
              IF (cpu_times(J).LT.cpu_min(I)) THEN
                cpu_min(I)=cpu_times(J)
                cpu_min_pe(I)=J
              ENDIF

            ENDDO ! loop over processors

            IF (wallclock_max(I) .GT. 0.0) THEN
              speedup(I)=cpu_total(I)/wallclock_max(I)
            ELSE
              speedup(I)=1.0
            ENDIF
            efficiency(I)=speedup(I)/nproc

! and calculate the statistics
! first calculate the means
            wallclock_mean(I)=wallclock_mean(I)/nproc
            cpu_mean(I)=cpu_total(I)/nproc
! To stop a divide by zero later:
            IF (wallclock_mean(I) .EQ. 0.0) wallclock_mean(I)=1.0E-20
            IF (cpu_mean(I) .EQ. 0.0) cpu_mean(I)=1.0E-20
! and now the standard deviation
            wallclock_sd(I)=0.0
            cpu_sd(I)=0.0
            DO J=0,nproc-1
              wallclock_sd(I)=wallclock_sd(I)+
     &          (wallclock_times(J)-wallclock_mean(I))*
     &          (wallclock_times(J)-wallclock_mean(I))
              cpu_sd(I)=cpu_sd(I)+(cpu_times(J)-cpu_mean(I))*
     &                      (cpu_times(J)-cpu_mean(I))
            ENDDO
            wallclock_sd(I)=SQRT(wallclock_sd(I)/nproc)
            cpu_sd(I)=SQRT(cpu_sd(I)/nproc)

! Calculate the median
            DO J=0,nproc-2
              DO K=J+1,nproc-1
                IF (wallclock_times(K) .GT. wallclock_times(J)) THEN
                  sortwork_real=wallclock_times(J)
                  wallclock_times(J)=wallclock_times(K)
                  wallclock_times(K)=sortwork_real
                ENDIF
                IF (cpu_times(K) .GT. cpu_times(J)) THEN
                  sortwork_real=cpu_times(J)
                  cpu_times(J)=cpu_times(K)
                  cpu_times(K)=sortwork_real
                ENDIF
              ENDDO
            ENDDO

            IF (MOD(nproc,2) .EQ. 0) THEN
              wallclock_median(I)=(wallclock_times((nproc/2)-1)+
     &                             wallclock_times(nproc/2))*0.5
              cpu_median(I)=(cpu_times((nproc/2)-1)+
     &                       cpu_times(nproc/2))*0.5
            ELSE
              wallclock_median(I)=wallclock_times(nproc/2)
              cpu_median(I)=cpu_times(nproc/2)
            ENDIF

          ENDIF ! am I PE 0?

        ENDDO ! loop over sections

! Sort and output the information on PE 0

        IF (mype .EQ. 0) THEN

          DO I=1,summ_n_timers-1
            DO J=(I+1),summ_n_timers
              IF (wallclock_max(J) .GT. wallclock_max(I)) THEN

! Swap the entries I and J

              sortwork_char=summ_section(I)
              summ_section(I)=summ_section(J)
              summ_section(J)=sortwork_char

              sortwork_real=wallclock_mean(I)
              wallclock_mean(I)=wallclock_mean(J)
              wallclock_mean(J)=sortwork_real

              sortwork_real=wallclock_median(I)
              wallclock_median(I)=wallclock_median(J)
              wallclock_median(J)=sortwork_real

              sortwork_real=wallclock_sd(I)
              wallclock_sd(I)=wallclock_sd(J)
              wallclock_sd(J)=sortwork_real

              sortwork_real=wallclock_max(I)
              wallclock_max(I)=wallclock_max(J)
              wallclock_max(J)=sortwork_real

              sortwork_real=wallclock_min(I)
              wallclock_min(I)=wallclock_min(J)
              wallclock_min(J)=sortwork_real

              sortwork_int=wallclock_min_pe(I)
              wallclock_min_pe(I)=wallclock_min_pe(J)
              wallclock_min_pe(J)=sortwork_int

              sortwork_int=wallclock_max_pe(I)
              wallclock_max_pe(I)=wallclock_max_pe(J)
              wallclock_max_pe(J)=sortwork_int

              sortwork_real=cpu_mean(I)
              cpu_mean(I)=cpu_mean(J)
              cpu_mean(J)=sortwork_real

              sortwork_real=cpu_median(I)
              cpu_median(I)=cpu_median(J)
              cpu_median(J)=sortwork_real

              sortwork_real=cpu_sd(I)
              cpu_sd(I)=cpu_sd(J)
              cpu_sd(J)=sortwork_real

              sortwork_real=cpu_max(I)
              cpu_max(I)=cpu_max(J)
              cpu_max(J)=sortwork_real

              sortwork_real=cpu_min(I)
              cpu_min(I)=cpu_min(J)
              cpu_min(J)=sortwork_real

              sortwork_real=cpu_total(I)
              cpu_total(I)=cpu_total(J)
              cpu_total(J)=sortwork_real

              sortwork_real=speedup(I)
              speedup(I)=speedup(J)
              speedup(J)=sortwork_real

              sortwork_real=efficiency(I)
              efficiency(I)=efficiency(J)
              efficiency(J)=sortwork_real

              sortwork_int=cpu_min_pe(I)
              cpu_min_pe(I)=cpu_min_pe(J)
              cpu_min_pe(J)=sortwork_int

              sortwork_int=cpu_max_pe(I)
              cpu_max_pe(I)=cpu_max_pe(J)
              cpu_max_pe(J)=sortwork_int

              ENDIF
            ENDDO
          ENDDO

! and write out the information
          WRITE(6,*)
          IF (timer_kind .EQ. 1) THEN
            WRITE(6,*) 'MPP : None Inclusive timer summary'
          ELSE
            WRITE(6,*) 'MPP : Inclusive timer summary'
          ENDIF

          WRITE(6,*)
          WRITE(6,*)  'WALLCLOCK  TIMES'
          WRITE(6,40)
          DO I=1,summ_n_timers

            WRITE(6,41) I,summ_section(I),
     &                 wallclock_mean(I),wallclock_median(I),
     &                 wallclock_sd(I),
     &                 (wallclock_sd(I)/wallclock_mean(I))*100.0,
     &                 wallclock_max(I),wallclock_max_pe(I),
     &                 wallclock_min(I),wallclock_min_pe(I)
          ENDDO

          WRITE(6,*)
          WRITE(6,*)  'CPU TIMES (sorted by wallclock times)'
          WRITE(6,40)
          DO I=1,summ_n_timers
            WRITE(6,41) I,summ_section(I),
     &                 cpu_mean(I),cpu_median(I),
     &                 cpu_sd(I),
     &                 (cpu_sd(I)/cpu_mean(I))*100.0,
     &                 cpu_max(I),cpu_max_pe(I),
     &                 cpu_min(I),cpu_min_pe(I)
          ENDDO

          WRITE(6,*)
          WRITE(6,*) 'PARALLEL SPEEDUP SUMMARY ',
     &             '(sorted by wallclock times)'
          WRITE(6,50)
          DO I=1,summ_n_timers
            WRITE(6,51) I,summ_section(I),cpu_total(I),
     &                  wallclock_max(I),speedup(I),
     &                  efficiency(I)
          ENDDO


 40       FORMAT(4X,'ROUTINE',11X,'MEAN',3X,'MEDIAN',7X,
     &           'SD',3X,'% of mean',6X,'MAX',3X,
     &           '(PE)',6X,'MIN',3X,'(PE)')

 41       FORMAT(/,I3,1X,A12,1X,
     &           3(1X,F8.2),5X,F6.2,'%'
     &           2(1X,F8.2,1X,'(',I4,')'))

 50       FORMAT(4X,'ROUTINE',11X,'CPU TOTAL',3X,
     &           'WALLCLOCK MAX',3X,'SPEEDUP',3X,
     &           'PARALLEL EFFICIENCY')

 51       FORMAT(/,I3,1X,A17,1X,
     &           1X,F8.2,8X,F8.2,2X,F8.2,14X,F8.2)

        ENDIF
        WRITE(6,*)


      ENDDO ! loop over timer kind

! Finally copy the timer info into the last_* arrays so that the
! intermediate timer can calculate the timings since this point

      DO I=1,in_number_of_timers
        last_in_number_of_times_timed(I)=in_number_of_times_timed(I)
        last_in_cpu_time_elapsed(I)=in_cpu_time_elapsed(I)
        last_in_wallclock_time_elapsed(I)=
     &    in_wallclock_time_elapsed(I)
      ENDDO

      DO I=1,ni_number_of_timers
        last_ni_number_of_times_timed(I)=ni_number_of_times_timed(I)
        last_ni_cpu_time_elapsed(I)=ni_cpu_time_elapsed(I)
        last_ni_wallclock_time_elapsed(I)=
     &    ni_wallclock_time_elapsed(I)
      ENDDO


 9999 CONTINUE

      RETURN

      END





!*******************************************************************
      REAL FUNCTION get_cpu_time()

! Gets the cpu time from the system.

      IMPLICIT NONE
      REAL SECOND

      get_cpu_time=SECOND()

      RETURN
      END

      REAL FUNCTION get_wallclock_time()

! Gets the wallclock time from the wallclock.

      IMPLICIT NONE
      REAL temp

      CALL TIMEF(temp)

      get_wallclock_time=temp
      RETURN
      END
