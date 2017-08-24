! *****************************COPYRIGHT******************************
! (c) CROWN COPYRIGHT 1997, METEOROLOGICAL OFFICE, All Rights Reserved.
!
! Use, duplication or disclosure of this code is subject to the
! restrictions as set forth in the contract.
!
!                Meteorological Office
!                London Road
!                BRACKNELL
!                Berkshire UK
!                RG12 2SZ
!
! If no contract has been raised with this copy of the code, the use,
! duplication or disclosure of it is strictly prohibited.  Permission
! to do so must first be obtained in writing from the Head of Numerical
! Modelling at the above address.
! ******************************COPYRIGHT******************************
!+
!
! Subroutine Interface:
      SUBROUTINE FOURIER(A,TRIGS,IFAX,INC,JUMP,N,LOT,ISIGN)
      implicit none
!
! Description:
!
!     SUBROUTINE 'FOURIER' - MULTIPLE FAST REAL PERIODIC TRANSFORM
!     UNIFIED MODEL RE-WRITE OF ECMWF ROUTINE FFT991
!
!     REAL TRANSFORM OF LENGTH N PERFORMED BY REMOVING REDUNDANT
!     OPERATIONS FROM COMPLEX TRANSFORM OF LENGTH N
!
!     INPUT INFORMATION ...
!     A IS THE ARRAY CONTAINING INPUT & OUTPUT DATA
!     TRIGS IS A PREVIOUSLY PREPARED LIST OF TRIG FUNCTION VALUES
!     IFAX IS A PREVIOUSLY PREPARED LIST OF FACTORS OF N
!     INC IS THE INCREMENT WITHIN EACH DATA 'VECTOR'
!         (E.G. INC=1 FOR CONSECUTIVELY STORED DATA)
!     JUMP IS THE INCREMENT BETWEEN THE START OF EACH DATA VECTOR
!     N IS THE LENGTH OF THE DATA VECTORS
!     LOT IS THE NUMBER OF DATA VECTORS
!     ISIGN = +1 FOR TRANSFORM FROM SPECTRAL TO GRIDPOINT
!           = -1 FOR TRANSFORM FROM GRIDPOINT TO SPECTRAL
!
!     ORDERING OF COEFFICIENTS:
!         A(0),B(0),A(1),B(1),A(2),B(2),...,A(N/2),B(N/2)
!         WHERE B(0)=B(N/2)=0; (N+2) LOCATIONS REQUIRED
!
!     ORDERING OF DATA:
!         X(0),X(1),X(2),...,X(N-1), 0 , 0 ; (N+2) LOCATIONS REQUIRED
!
!     VECTORIZATION IS ACHIEVED ON CRAY BY DOING THE TRANSFORMS
!     IN PARALLEL
!
!     N MUST BE COMPOSED OF FACTORS 2,3 & 5 BUT DOES NOT HAVE TO BE
!     EVEN
!
!     DEFINITION OF TRANSFORMS:
!     -------------------------
!
!     ISIGN=+1: X(J)=SUM(K=0,...,N-1)(C(K)*EXP(2*I*J*K*PI/N))
!         WHERE C(K)=A(K)+I*B(K) AND C(N-K)=A(K)-I*B(K)
!
!     ISIGN=-1: A(K)=(1/N)*SUM(J=0,...,N-1)(X(J)*COS(2*J*K*PI/N))
!               B(K)=-(1/N)*SUM(J=0,...,N-1)(X(J)*SIN(2*J*K*PI/N))
!     ..................................................................
!     MULTITASKING NOTE: ERROR TRAPPING NOT WORKING IN PRESENT FORM IF
!                        CODE IS MULTITASKED.
!
!     NOT SUITABLE FOR SINGLE COLUMN USE.
!
!     VERSION FOR CRAY Y-MP
!     REWRITTEN TO UNIFIED MODEL PROGRAMMING STANDARDS FROM ECMWF
!     ORIGINAL CODE BY M.H.MAWSON
!     ORIGINAL CODE WRITER C. TEMPERTON
!
! History:
! Version   Date     Comment
! -------   ----     -------
!
!    4.2*    06/12/96 * Not included in the UM until 4.3.
!                     This version gives exactly the same results
!                     on the C90 and T3E as those returned in
!                     libmet.a on the C90.  It deals with inc being
!                     values other than 1.  It does not have the
!                     following error:
!                     On the grid to spectral transform it returns
!                     the the B terms as being positive, whereas
!                     FFT991 returns the negative for the B terms.
!                     On the inverse transform a change of sign in
!                     required again if a change is effected on the
!                     forward transform.
!                     Supplied by A Mills, Cray.  K Rogers
!    4.3    11/4/97   Corrections made to vn4.2
!                     Additional routine InitFourier is not to be
!                     used by this version of Fourier. It has been
!                     included for use with a later version
!                     of Fourier which takes advantage of storing the
!                     factor and trig tables.
!                     Supplied by A Mills, Cray.  I Edmond
!    4.5    09/7/98   Added direct call to 'INC_1_FOURIER'
!                     if 'INC' is unity.  INC_1_FOURIER is
!                     essentially FOURIE2A with some additional
!                     checks.  This call avoids any data copies.
!                       Author: Bob Carruthers, Cray Research
!
! This routine is compatable with routine in FOURIER
!
! It calls the Cray standard routines scfft and csfft.
!
! It must rearrange the data before and after the call to scfft and csff
! The rearrangement of the data is dependent on whether the transform
! is grid to spectral or spectral to grid.
!
! Be careful to understand all of this before changing anything!
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!
! System component covered: <appropriate code>
! System Task:              <appropriate code>
!
! Declarations:

! 1.0 Subroutine arguments
!   1.1 Scalar arguments with intent(in):
      integer
     & inc
     &,n
     &,lot
     &,jump
     &,isign

!   1.2 Array arguments with intent(in):
      integer ifax(10)

      real    trigs(3*n)

!   1.3 Array arguments with intent(in/out):
      real a(1+jump*(lot-1)+inc*(n+1))

!   2.1 Local scalars:
      integer
     &        k
     &       ,i
     &       ,l_sign
     &       ,isy

      real
     &     scale
     &    ,local

!   2.2 Local arrays:
      real
     &     l_work  (  4 + 4 * n) ! this will work for both MPP and PVP
     &    ,l_table (100 + 4 * n) ! this will work for both MPP and PVP
     &    ,l_data  (n+2)

!   3.0 Local parameters:
      integer forward
        parameter(forward=1)         ! grid to spectral in scfft
      integer backward
        parameter(backward=-forward) ! spectral to grid in csff
!- End of header

      isy  = 0

      call scfft(     0, n, scale, l_data, l_data, l_table, l_work, isy)
      if (isign.eq.+1) then
        l_sign = backward
        scale  = 1.0
        do k   = 1, lot
          do i = 1, n+1, 2
            l_data(i)    =  a( 1+(k-1)*jump +(i-1)*inc)
            l_data(i+1)  = -a( 1+(k-1)*jump + i   *inc)
          enddo

      call csfft(l_sign, n, scale, l_data, l_data, l_table, l_work, isy)
          l_data(n+1) = 0.0
          l_data(n+2) = 0.0
          do i = 1, n+2
            a( 1+(k-1)*jump + (i-1) *inc) = l_data(i)
          enddo
        enddo
      else
        l_sign = forward
        scale  = 1.0/n
        do k   = 1, lot
          do i = 1, n
            l_data(i) = a( 1+(k-1)*jump + (i-1) *inc)
          enddo
          l_data(n+1) = 0.0
          l_data(n+2) = 0.0

      call scfft(l_sign, n, scale, l_data, l_data, l_table, l_work, isy)
          do i = 1, n+1, 2
            a( 1+(k-1)*jump +(i-1)*inc) =  l_data(i)
            a( 1+(k-1)*jump + i   *inc) = -l_data(i+1)
          enddo
        enddo
      endif

      end


!+
!
! Subroutine Interface:
      SUBROUTINE InitFourier(N,IFAX,TRIGS)
      implicit none
!
! Description:
!
!C     SUBROUTINE 'SET66' - COMPUTES FACTORS OF N & TRIGONOMETRIC
!C     FUNCTIONS REQUIRED BY STAGGERED SINE TRANSFORM Var_FFT66
!C     AND STAGGERED COSINE TRANSFORM Var_FFT55
!C
!C     AUTHOR: CLIVE TEMPERTON JANUARY 1987
!C
! modified Andrew Lorenc Oct 1995:
! the "standard" TRIGS for the FFT routines are put in TRIGS(1:N).
! extra TRIGS for the shifted (co)sine transforms are put in
! TRIGS(N+1:N+N).
! extra TRIGS for the sine transforms are put in TRIGS(2N+1:2N+N).
!
! Parent module: VarMod_Fourier
!
! History:
! Version   Date     Comment
! -------   ----     -------
! 4.3       11/4/97   Original code. Alistair Mills
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!
! System component covered: <appropriate code>
! System Task:              <appropriate code>
!
! Declarations:
!
! 1.0 Subroutine arguments
!   1.1 Scalar arguments with intent(in):
      Integer
     & N

!   1.2 Array arguments with intent(in):
      Integer
     & IFAX(10)

!   1.3 Array arguments with intent(out):
      Real
     & TRIGS(3*N)

!   2.1 Local scalars:
      Integer
     & I
     &,K
     &,L
     &,NU
     &,NFAX
     &,IFAC

      Real
     & DEL
     &,NIL
     &,NHL
     &,ANGLE

!   2.2 Local arrays:
      Integer
     & JFAX(10)
     &,LFAX(7)

      Character*80 ErrorMessage

      Data LFAX/6,8,5,4,3,2,1/

!- End of header

!     TRIGS FOR REAL PERIODIC TRANSFORM
!     ---------------------------------
      DEL=4.0*ASIN(1.0)/FLOAT(N)
      NIL=0
      NHL=(N/2)-1
!DIR$ IVDEP
! Fujitsu vectorization directive
!OCL NOVREC
      DO 10 K=NIL,NHL
      ANGLE=FLOAT(K)*DEL
      TRIGS(2*K+1)=COS(ANGLE)
      TRIGS(2*K+2)=SIN(ANGLE)
   10 CONTINUE
!
!     EXTRA TRIGS FOR (CO)SINE TRANSFORM
!     ----------------------------------

      DEL=0.5*DEL
      DO K=1,NHL
       ANGLE=FLOAT(K)*DEL
       TRIGS(2*N+K)=SIN(ANGLE)
      END DO

!     EXTRA TRIGS FOR SHIFTED (CO)SINE TRANSFORM
!     --------------------------------------
      DEL=0.5*DEL
      DO 15 K=1,N
      ANGLE=FLOAT(K)*DEL
      TRIGS(N+K)=SIN(ANGLE)
   15 CONTINUE

!
!     NOW FIND FACTORS OF N
!     ---------------------
!
! AL 9/10/95.  Removed code following Clive Temperton's advice:
! the "fix for small N" section in SET66 was needed for obscure
! historical reasons in the version which called CAL routines.
! It's not necessary in the all-Fortran version, and removing it
! should make things go a bit faster.
!
!     FIND FACTORS OF N (8,6,5,4,3,2; ONLY ONE 8 ALLOWED)
!     LOOK FOR SIXES FIRST, STORE FACTORS IN DESCENDING ORDER
      NU=N
      IFAC=6
      K=0
      L=1
   20 CONTINUE
      IF (MOD(NU,IFAC).NE.0) GO TO 30
      K=K+1
      JFAX(K)=IFAC
      IF (IFAC.NE.8) GO TO 25
      IF (K.EQ.1) GO TO 25
      JFAX(1)=8
      JFAX(K)=6
   25 CONTINUE
      NU=NU/IFAC
      IF (NU.EQ.1) GO TO 50
      IF (IFAC.NE.8) GO TO 20
   30 CONTINUE
      L=L+1
      IFAC=LFAX(L)
      IF (IFAC .GT. 1) GO TO 20

      ! Illegal factors:
      write (ErrorMessage, '(a,i4,a)') 'N = ', N,
     &                     ' contains illegal factors.'

      GOTO 9
!
!     NOW REVERSE ORDER OF FACTORS
   50 CONTINUE
      NFAX=K
      IFAX(1)=NFAX
      DO 60 I=1,NFAX
      IFAX(NFAX+2-I)=JFAX(I)
   60 CONTINUE

 9    continue

      END
