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

! subroutine to convert fields from blocks to full grid of datapoints
!
! Subroutine Interface:

          subroutine unblok(field,blkfld,mdata,
C*    *PARAMETER*  FOR ARRAY DIMENSIONS.  parall
C
     & NANG, NFRE, NGX, NGY, NBLO, NIBLO, NOVER,
     & NIBLD, NBLD, NIBLC, NBLC,
C*    *COMMON* *GRIDPAR*  GENERAL GRID INFORMATION.
C
     & DELPHI, DELLAM, SINPH, COSPH, IGL, IJS, IJL2, IJLS, IJL, IJLT,
C
     & icode)

      IMPLICIT NONE
!
! Description:
! This subroutine converts the input array of blocked data used by WAM:
! fldout(niblo,nblo) to data on the full grid of datapoints:
! array(mdata)
!
! Method:
!   The WAM model grid and block indexing arrays are used
!
! Current Code Owner: Martin Holt
!
! History:
! Version   Date     Comment
! -------   ----     -------
! <1.0>    July 1995 Original code.    M Holt
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!
! System component covered: <UM wave model>
! System Task:              <UM wave model>
!
!- End of header

C*    *PARAMETER*  FOR ARRAY DIMENSIONS.  parall
       INTEGER
     & NANG,       ! number of direction components
     & NFRE,       ! number of frequency components
     & NGX,        ! number of cols in LS mask grid
     & NGY,        ! number of rows in LS mask grid
     & NBLO,       ! max number of blocks
     & NIBLO,      ! max number datapoints per block
     & NOVER,      ! max number datapoints in overlap row
     & NIBLD, NBLD, NIBLC, NBLC
C
C*    *COMMON* *GRIDPAR*  GENERAL GRID INFORMATION.
C
      real DELPHI  ! grid increment for latitude (in metres)
      real DELLAM  ! grid increment for long. at equator (in metres)
      real SINPH(NGY), COSPH(NGY) ! sin / cos of latitude

      integer IGL   ! number of blocks
      integer IJS(NBLO)  ! index of first point of second lat
      integer IJL2(NBLO) ! index of last  point of second lat
      integer IJLS(NBLO) ! index of first point of lat before last
      integer IJL(NBLO)  ! index of last  point of lat before last
      integer IJLT(NBLO) ! total number of gridpoints in a block
C

       integer mdata           ! array size - max number of datapoints
       integer icode           ! out return code
       real field(mdata)       ! out   field of values at data points
       real blkfld(niblo,nblo) ! in    field arranged by blocks

! local variables

       integer nstart  ! index on full grid of first point in block
       integer nend    ! index on full grid of last  point in block
       integer ip      ! loop counter over all data points
       integer ig      ! loop counter over all blocks
       integer ifill,i         ! index

       do ip=1,mdata
         field(ip)=0.
       enddo
c
c   copy only from points ijs to ijl in each block, to ensure that time
c   dependent arrays such as U10 , z0 etc are correctly handled at
c   overlap rows
c
c   nstart should = 1 as ijs(1) = 1

       nstart=ijs(1)

       do ig=1,igl

        nend  = nstart + ijl(ig) - ijs(ig)
        ifill=ijs(ig)

        do i=nstart,nend
         field(i)=blkfld(ifill,ig)
         ifill=ifill+1
        enddo

        nstart=nend + 1

       enddo

       return
      end                                                              
