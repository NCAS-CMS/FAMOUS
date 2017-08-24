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
!
! Description:
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
! UM4.1    June 1996 Original code.  M Holt
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!
!- End of header

      subroutine wavetr(energ,pswh,perio,pdir,kwtot,pfreq,dff,ptheta,
     +len1,ndata,kwtmax,nfreq,ntheta,rmdi,icode)
c
c   top level subroutine to interface wavetrain programs with UKMO wave
c   model.   This subroutine is called by waveh from fldout.
c
c     arguments
c     #########
c     energ     in  energy array - one-dimensional - length len1
c     nfreq     in  number of frequencies
c     ntheta    in  number of directions
c     ndata     in  number of data points
c     len1      in  length of energy array = nfreq*ntheta*ndata
c     pfreq     in  frequency array
c     dff       in  frequency intervals array
c     ptheta    in  direction array
c     rmdi      in  real - missing data indicator
c     kwtmax    in  max number of wavetrains searched for
c
c               OUT
c     pswh      out sig wave height (ndata by kwtmax)
c     perio     out wave period     (   "           )
c     pdir      out wave direction  (   "           )
C                                   radians TO/ zero=east
c     kwtot     out number of wave trains at each gridpoint
C     icode     return code from this subroutine
c
C   * local but passed to main wavetrain processing* results not used
c     pfwind    wind speeds (ndata)     )  removed from input list
c     pdwind    wind directions (ndata) )
C     not used at present. resized to (1)

      integer nfreq,ntheta,ndata,kblo,kjs,kjl,kdang,nblok
      integer kwtra(ndata,ntheta,nfreq),kwtot(ndata)
C
      real pfwind(1),pdwind(1),pfreq(nfreq),energ(len1)
      real pdmax,pecut,peminr,pemaxr,pdtmin,pfbin
      real pmiss,pres,ptheta(ntheta),dff(nfreq)
      real pswh(ndata,kwtmax),perio(ndata,kwtmax),pdir(ndata,kwtmax)
c
c

C*L------------------COMDECK C_PI---------------------------------------
CLL
CLL 4.0 19/09/95  New value for PI. Old value incorrect
CLL               from 12th decimal place. D. Robinson
CLL
      REAL PI,PI_OVER_180,RECIP_PI_OVER_180

      PARAMETER(
     & PI=3.14159265358979323846, ! Pi
     & PI_OVER_180 =PI/180.0,   ! Conversion factor degrees to radians
     & RECIP_PI_OVER_180 = 180.0/PI ! Conversion factor radians to
     &                              ! degrees
     & )
C*----------------------------------------------------------------------


C
c     # set blocking information #
c     # note - process ndata/3 points per call
C     # because of memory restriction
c     # if larger grid used may need ndata/4 or whatever.
c     # at present, for oper global this runs in 16Mw

      icode=-1

      kblo=ndata

      nblok=8
      kjs=1
      kjl=int(ndata/nblok)
c
      do ip=1,ndata
       kwtot(ip)=0
      enddo

c     # set arguments # - see wavetrain.prog for details #
c      kdang is max spread of direction bins per wave train #
      kdang=4
      kflagws=0
      pdmax=0.33333*pi
      pecut=0.001
      peminr=1/1.3
      pemaxr=1.3
      pdtmin=0.25*pi
      pfbin=0.0
      pres=1000.

c     # note the wavetrain routine requires pmiss negative
      pmiss=-32768.

      pmcoef=0.8
      kreosp=0

      do j=1,kwtmax
       do ip=1,ndata
        pswh(ip,j)=pmiss
        perio(ip,j)=pmiss
        pdir(ip,j)=pmiss
       enddo
      enddo
c
      do ii=1,nblok

c      do the first block of points *
      WRITE(6,*)' processing kjs  to kjl ',kjs,kjl
      WRITE(6,*)' kwtmax is ',kwtmax
      call wtrain(energ,kblo,kjs,kjl,ntheta,nfreq,pfwind,pdwind,
     +            pfreq,pfbin,ptheta,pres,kdang,pdmax,
     +            pecut,peminr,pemaxr,pdtmin,kwtmax,
     +            pmiss,pswh,perio,pdir,kwtot,
     +            kflagws,pmcoef,kreosp,kwtra,dff)
c
c
c      do the second block of points *
      kjs=kjl+1
      kjl=int((ii+1)*ndata/nblok)
      if(kjl.ge.ndata) kjl=ndata

      enddo


cc    here replace pmiss with mdi

      do j=1,kwtmax
       do ip=1,ndata
        if(pswh(ip,j).eq.pmiss) pswh(ip,j)=rmdi
        if(perio(ip,j).eq.pmiss) perio(ip,j)=rmdi
        if(pdir(ip,j).eq.pmiss) pdir(ip,j)=0.
       enddo
      enddo

         WRITE(6,*)'setting mdis for absent trains : routine wavetr'
        do ip=1,ndata
         jstart=min(kwtot(ip)+1,kwtmax)
         do j=jstart,kwtmax
          pdir(ip,j)=rmdi
          perio(ip,j)=rmdi
          pswh(ip,j)=0.
         enddo
        enddo
c
      icode=0
      return
      end
