      SUBROUTINE PIPE_Init_Model(cdmodnam,jm)
C****
C               *****************************
C               * adapted from oasis PIPE_Init
C               * -------------     ------- *
C               *****************************
C
C**** *PIPE_Init_Model*  - Initialize coupled mode communication for
C      one model.
C
C     Purpose:
C     -------
C     Initialize coupler - model communication using named
c     pipes(fifo).
C
C**   Interface:
C     ---------
C       *CALL*  *PIPE_Init_Model(cdmodnam,jm)
C
C     Input:
C     -----
C            cdmodnam : remote models names (char string)
c                       after this routine PIPE_Init_Model has run,
c                       each of the models will be able to :
c                       o write to the pipe Preadm?? when writing to
c                       'WT'//cdmodnam
c                       o read from the pipe Pwritm?? when reading
c                       from 'RD'//cdmodnam
c                       using unix named pipes.
c                       NB : take care of NOT using the
c                       same cdmodnam into the coupler and the models
c                       as they are simple unix files whose name must
c                       be unique...
C                jm   : number of the model according to the namelist
c                       of the Oasis coupler (intger).
C
C     Output:
C     ------
C                kinfo    : error status
C
C     Workspace:
C     ---------
C     None
C
C     Externals:
C     ---------
C     assign, mknod, fsigctl, ,empty, ferror, getfpe
C
C     Reference:
C     ---------
C     See OASIS manual (1995)
C
C     History:
C     -------
C       Version   Programmer     Date      Description
C       -------   ----------     ----      -----------
C       2.0       JC. Thil       96/05/03  created
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
C* ---------------------------- Argument declarations ---------------

C
      CHARACTER*6 cdmodnam
      INTEGER     jm
C
C* ---------------------------- Local declarations ------------------
C
      CHARACTER*90  clcmd
      INTEGER       mknod
      CHARACTER*8   cprnam, cpwnam
      integer       nulou
C
C* ---------------------------- Poema verses ------------------------
C
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
C*    1. Initializations
C        ---------------
C
      data nulou/6/
C
C* Zeroes error codes
C
      ierror = 0
      icumul = 0

C*    3. Create pipes to communicate between the ocean model and the
C        -------------------------------------------------------------
c        coupler.
C        --------
C
C
C* Assign business for reading
C
C ---> Options u and unblocked valid for sequential unformatted files
C      Option u results in an immediate system call (see man assign)
C
c excerpt of the Cray documentation :
c
c                The undefined file structure
c                (specified by assign -s u) should be specified
c                for a pipe by the sending process. The unblocked
c                structure (specified by assign -s unblocked)
c                should be specified for a pipe by the receiving
c                process.cc
c
c                The file structure for the pipe of the sending
c                (write) process should be set to undefined
c                (assign -s u), which issues a system call for
c                each write.  You can also select a file
c                specification of system (assign -F system) for
c                the sending process.

c     from that end  Pwritm?? is the pipe we read from
c     hence the -s unblocked option.
      cprnam = 'WT'//cdmodnam
      WRITE (UNIT = clcmd,FMT ='
     $  (''assign -s unblocked -a Pwritm'',i2.2,
     $  '' f:'', A8 )') jm, cprnam
      print *, clcmd
      CALL assign(clcmd, ierror)
      CALL empty(clcmd,90)
C
C*    Find error if any
C
      IF (ierror .EQ. 0) THEN
        write(nulou,*)
     $    'Assign command done for pipe ', 'Pwritm', jm,
     $    ' file ', cprnam
      ELSE
        write(nulou,*)
     $    'Problem with assign command for pipe file',
     $    cprnam
      ENDIF
C
C*    Assign business for writing
C
c     from that end  Preadm?? is the pipe we write to,
c     hence the -s u option.
      cpwnam = 'RD'//cdmodnam
      WRITE (UNIT = clcmd,FMT ='
     $  (''assign -s u -a Preadm'',i2.2,
     $  '' f:'',A8 )') jm, cpwnam
      print *, clcmd
      CALL assign(clcmd, ierror)
      CALL empty(clcmd,90)
C
C*    Find error if any
C
      IF (ierror .EQ. 0) THEN
        write(nulou,*)
     $    'Assign command done for pipe', 'Preadm', jm,
     $    ' file ', cpwnam
      ELSE
        write(nulou,*)
     $    'Problem with assign command for pipe file',
     $    cpwnam
      ENDIF


C
C*    Mknod business
C
C ---> mode 4480: named pipe with rw permission for user
C      decimal value of octal 0010600 (see man mknod)
C
C
C* Reading pipes
C
      WRITE (UNIT = clcmd,FMT ='(''Preadm'',i2.2)')jm
      ierror = mknod (clcmd, 4480, 0)
      CALL empty(clcmd,90)
C
C*    Find error if any
C
      IF (ierror .EQ. 0) THEN
        write(nulou,*)
     $    'Mknod command done for pipe Preadm', jm
      ELSE
        write(nulou,*)
     $    'Mknod command useless for pipe Preadm', jm
        write(nulou,*)
     $    'File exists. The return code number is', ierror
      ENDIF
C
C*    Writing pipes
C
      WRITE (UNIT = clcmd,FMT ='(''Pwritm'',i2.2)')jm
      ierror = mknod (clcmd, 4480, 0)
      CALL empty(clcmd,90)
C
C*    Find error if any
C
      IF (ierror .EQ. 0) THEN
        write(nulou,*)
     $    'Mknod command done for pipe Pwritm', jm
      ELSE
        write(nulou,*)
     $    'Mknod command useless for pipe Pwritm', jm
        write(nulou,*)
     $    'File exists. The return code number is', ierror
      ENDIF
C
C
      RETURN
      END

