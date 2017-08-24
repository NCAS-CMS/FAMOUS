!
! Description:
!  This subroutine is part of the wavetrain diagnostic output code
!  developed by Anne Guillaume at MeteoFrance and ECMWF.
!  Introduced into the unified wave moel at UM4.1
!
! Method:
!
!
!
! Current Code Owner: Martin Holt
!
! History:
! Version   Date     Comment
! -------   ----     -------
! UM4.1    June 1996 Code introduced to UM.  M Holt
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!
!- End of header

      SUBROUTINE TRHOU ( KNRJ, KBLO, KJS, KJL, KANG, KFRE,
     +                   KPICA, KPICF, KDANG, KNR )

C
C Subroutine TRHOU
C
C Purpose: to identify a single wave train given the position of
C a maximum in energy as the starting point.
C
C The algorithm scans away from the starting point marking values
C within the wave train, ending the wave train when the values start
C to increase of the wave train becomes joined only by a single
C filament.
C
C Modification History
C
C Versions 1-4: A Guillame, ECMWF. 1991-1994. Original algorithm
C Version  5  : S Foreman & M Holt, UK Met Office, February 1996
C                                   Improve efficiency on C90
C
C Comments on coding
C As it stands, the code does not autotask effectively because of
C dependencies in the outer loop. Stripmining the inner loop is the
C only effective way of multitasking. In the wave model the stripmining
C is effectively done before entry to the routine, making further
C efficiencies difficult.
C
C

      IMPLICIT NONE

C
C Arguments
C
      INTEGER
     +       KBLO                      ! IN   Number of horizontal point
     +      ,KJS                       ! IN   Start point within array
     +      ,KJL                       ! IN   End point in array
     +      ,KANG                      ! IN   Number of angles
     +      ,KFRE                      ! IN   Number of frequencies
     +      ,KPICA (KBLO)              ! IN   Directions of the peaks
     +      ,KPICF (KBLO)              ! IN   Frequencies of the peaks
     +      ,KDANG    ! IN   Number of directions to search on each side
     +       ,KNRJ (KBLO, KANG, KFRE)  ! IN   Normalised energy array
     +      ,KNR (KBLO, KANG, KFRE)    ! OUT  Mask: 1=in train, 0=out


C
C Local variables
C

      INTEGER
     +       ang        ! Direction - loop counter
     +      ,ang_sweep  ! Sense of change in direction: loop counter
     +      ,ang_lim    ! Limit of loop on direction
     +      ,ang_start  ! Starting point for direction loop
     +      ,angle      ! Current direction
     +      ,angle_close ! Direction closer to peak
     +      ,angle_off  ! Offset to be used to look closer to peak
     +      ,angle_this ! 1 if ang_off = 0
     +      ,freq        ! Frequency loop counter
     +      ,freq_sweep  ! Sense of change of frequency: loop counter
     +      ,freq_lim    ! Limit of loop oover frequency
     +      ,freq_start  ! Starting frequency displacement in sweep
     +      ,frqcy       ! Current frequency
     +      ,frqcy_close ! Frequency closer to peak
     +      ,frqcy_off   ! Offset used to look closer to peak
     +      ,frqcy_this  ! 1 if frqcy_off = 0
     +      ,max_around  ! Maximum energy around current
c                        ! direction and frequency,
c                        ! looking towards the peak
c
C
C Local arrays
C

      INTEGER
     +       angle_is (-KANG:2*KANG) ! Mapping to array position
     +      ,frqcy_is (-KFRE:2*KFRE) ! Mapping to array position



      INTEGER  ! Debugging variables
     +       j           ! Loop counter over horizontal points
     +      ,k           ! Loop counter over freeq
     +      ,m           ! Loop counter over dir

c       WRITE(6,*)' '
c       WRITE(6,*)'array knrj on entry to routine trhou'
c       WRITE(6,*)'            spectrum'
c       WRITE(6,*) ' '
c       do k=1,kfre
c        write(6,100) (knrj(1,m,k),m=1,kang)
c       enddo
c 100   format(16I10)
c       WRITE(6,*)' '




C
C Fill the index arrays that will be used to convert from calculated
C indexes to the actual positions in the array
C
      DO frqcy = -KFRE, 2*KFRE

       frqcy_is(frqcy) = MAX (frqcy, 1)   ! Map points below range
       frqcy_is(frqcy) = MIN (frqcy_is(frqcy), KFRE) ! Points above

      END DO     ! frqcy - setting frequency mapping

      DO angle = -KANG, 2*KANG

       angle_is(angle) = MOD ( angle +(KANG-1), KANG ) +1

      END DO


C
C Initialise the output array to zero (no point in wave train)
C


      DO freq = 1, KFRE
         DO ang = 1, KANG
            DO j = KJS, KJL
               KNR(j,ang,freq) = 0
            END DO
         END DO
      END DO

C
C Initialise output array around the peaks -
C these values must be in the wave train unless there is no energy there
C

      DO freq = -1, 1
       DO ang = -1, 1
        DO j = KJS, KJL

         angle = angle_is( ang + KPICA(j) )
         frqcy = frqcy_is( freq + KPICF(j) )

         KNR(j,angle,frqcy) = MIN( KNRJ(j, angle, frqcy), 1 )

        END DO
       END DO
      END DO



C
C Sweep around the peaks that have already been identified
C Look each direction followed by each frequency.
C

      DO freq_sweep = -1, 1, 2

       freq_lim = (KFRE - 1) * freq_sweep

       IF (freq_sweep .EQ. -1) THEN
        freq_start = 0
       ELSE
        freq_start = 1
       END IF


       DO freq = freq_start, freq_lim, freq_sweep

        DO ang_sweep = -1, 1, 2

         ang_lim = KDANG * ang_sweep ! Sweep KDANG directions eacy side

         IF (freq .EQ. 0) THEN

          ang_start = ang_sweep ! Do not need to include the peak point

          frqcy_off = 0         ! Do not look other side of peak
          frqcy_this = 1        ! and don't mask by present value

         ELSE

          ang_start = 0  ! But do need to include other points

          frqcy_off = freq_sweep ! Look towards peak
          frqcy_this = 0          ! Use appropriate spectral value

         END IF


         DO ang = ang_start, ang_lim, ang_sweep

          IF (ang .EQ. 0) THEN

           angle_off = 0    ! Don't look other side of spectral peak
           angle_this = 1   ! Don't reject because point not set yet

          ELSE

           angle_off = ang_sweep ! Look towards peak
           angle_this = 0        ! Include points closer to peak

          END IF


          DO j = KJS, KJL

C
C Define the positions relative to the current direction and frequency
C

           angle = angle_is( KPICA (j) + ang )    ! *_is converts to
           frqcy = frqcy_is( KPICF (j) + freq )   ! position in array

C
C And the points closer to the peak
C

           angle_close = angle_is( angle - angle_off )
           frqcy_close = frqcy_is( frqcy - frqcy_off )


C
C Seek the maximum of adjacent values closer to peak
C

           max_around = MAX (

     *           KNR(j,angle_close,frqcy) * KNRJ(j,angle_close,frqcy)
     *         , KNR(j,angle,frqcy_close) * KNRJ(j,angle,frqcy_close)
     *         , KNR(j,angle_close,frqcy_close)
     *                            * KNRJ(j,angle_close,frqcy_close)
     *         )



C
C Accept the point only if the energy not more than the surroundings
C and if there is a continuous area to it
C Take care near frequency and direction of peak not to reject
C the point because it has not been examined yet - angle_this
C and frqcy_this allow for this.
C

           KNR(j,angle,frqcy) =
     +           MIN ( KNRJ(j,angle,frqcy), 1)
     *         * MIN ( MAX(max_around-KNRJ(j,angle,frqcy)+1, 0), 1)
     *         * MAX( KNR(j,angle,frqcy_close), frqcy_this)
     *         * MAX( KNR(j,angle_close,frqcy), angle_this)




          END DO ! j          - over points
         END DO  ! ang        - over directions
        END DO   ! ang_sweep  - over direction sign
       END DO    ! freq       - over frequency
      END DO     ! freq_sweep - over frequency sign


c       WRITE(6,*)' '
c       WRITE(6,*)'array knr on exit from routine trhou'
c       WRITE(6,*)'             spectrum'
c       WRITE(6,*) ' '
c       do k=1,kfre
c        write(6,100) (knr(1,m,k),m=1,kang)
c       enddo
c
c       WRITE(6,*)' '

      RETURN
      END
