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
      SUBROUTINE HorizontalInterp(
     &           Scheme, LGlobal, LVector, 
     &           Len1In,  Len2In,
     &           Len1Out, Len2Out,
     &           LambdaOrigin,
     &           PhiOrigin,
     &           DLambdaIn,
     &           DPhiIn,
     &           LambdaOut,
     &           PhiOut,
     &           DataIn,
     &           DataOut,
     &           ErrorStatus,ErrorMessage)

! Description:
!   Interpolate by linear, cubic or quintic lagrange interpolation, a 2d
!   field to a 2d set of points.  Input can be on a sphere or within a
!   rectangular box. Requested output points must lie within the region
!   defined for the input data.
!     An option for Monotonicity can be enforced (currently disabled).
!
!     i=x=lambda, j=y=phi
!     assumption that nth point defined by:
!        (LambdaOrigin+(n-1)*DLambdaIn , PhiOrigin+(n-1)*DPhiIn)
!     assumption that input and output in radians in range
!        (lambda 0 to 2pi ; phi -pi/2 to +pi/2)                         
!     assumption that check for output positions outside input grid
!         domain is done externally
!
! Method: This is a modified version of the routine interpolation
!         written by Mark Mawson and described in:
!           The proposed semi-Lagrangian advection scheme for the
!              semi-Implicit Unified Model integration scheme.
!                    F.R. Division working paper No 162.
!                             Mark H. Mawson
!
!
! Owner: Stuart Bell
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   4.0   6/6/95   Equiv. to VAR code as at time of build:Stuart Bell
!   4.4   23/4/97  Relax checks so C90 dump OK on T3E:Stuart Bell
!   4.5   23/04/98 Add LVector arg. to enable Polar wind correction:SB
!
! Code Description:
!   Language:           Fortran 77 plus
!   Software Standards: "UM and Met O standards".
!
!
! Declarations:

      IMPLICIT NONE

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


!* Subroutine arguments
! Scalar arguments with :intent(in)
      INTEGER  Scheme   ! a code saying which order
!                                     ! of scheme to use:
!                                     ! 1 = linear;3 = cubic;5 = quintic
      LOGICAL  LGlobal      ! True, if global
      LOGICAL  LVector      ! True, if wind component

      INTEGER     Len1In    ! Dimension of DataIn in i direction.
      INTEGER     Len2In    ! Dimension of DataIn in j direction.
      INTEGER     Len1Out   ! Dimension of DataOut in i direction.
      INTEGER     Len2Out   ! Dimension of DataOut in j direction.
!  ErrorStatus
      INTEGER       ErrorStatus
      CHARACTER*256 ErrorMessage

      REAL  DLambdaIn    ! Holds spacing between points in
!                        ! the i direction for the input data field.
      REAL  DPhiIn       ! Holds spacing between points in
!                        ! the j direction for the input field.
      REAL  LambdaOrigin ! Position of 1st point in
!                        ! the i direction for the input data field.
      REAL  PhiOrigin    ! Position of 1st point in
!                        ! the j direction for the input field.

! Array  arguments with :intent(in)
      REAL  DataIn (Len1In,Len2In)      ! Data to be interpolated.

      REAL  LambdaOut (Len1Out,Len2Out)   ! Lambda coordinate of
!                                         ! output data on input grid
      REAL  PhiOut (Len1Out,Len2Out)      ! Phi coordinate of
!                                         ! output data on input grid

! Array  arguments with :intent(out)
      REAL  DataOut (Len1Out,Len2Out)     ! Data interpolated
!                                         ! to desired locations.

!* End of Subroutine arguments

! Local named constants:
      INTEGER Extend   ! extra points for DataExt
      LOGICAL L_mono   ! True, if interpolation required to be monotone
      PARAMETER  (Extend=3)
      PARAMETER  (L_mono=.false.)

! Local scalars:
      INTEGER   LowerBound   ! Lower Bound of DataExt in both directions
      INTEGER   i            !} Loop
      INTEGER   j            !} indices.

      INTEGER   HalfLen1In   ! Len1In/2
      INTEGER   VIorder      ! order of interpolation
      INTEGER   SwapPole     ! 1 if no sign swap, -1 otherwise

      REAL      HalfPi          ! pi/2
      REAL      TwoPi           ! pi*2
      REAL      RecipDLambdaIn  ! 1/DLambdaIn
      REAL      RecipDPhiIn     ! 1/DPhiIn
      REAL      HalfDPhiIn      ! DPhiIn/2

      REAL      Test1,Test2      !used in origin checking
      LOGICAL   FirstOnPole      ! "         "        "

! Local arrays:
      INTEGER IOut (Len1Out,Len2Out)
      INTEGER JOut (Len1Out,Len2Out)

      REAL   DataExt(1-Extend:Len1In+Extend,
     &          1-Extend:Len2In+Extend)  ! Data in extended
      REAL  DataHigh (Len1Out,Len2Out)
      REAL  DataMono (Len1Out,Len2Out)
      REAL  WtLambda (Len1Out,Len2Out)
      REAL  WtPhi (Len1Out,Len2Out)
      REAL P1,P2             ! used for relaxing real number checks
      LOGICAL LEQR           ! 'Logical EQual Real'


!- End of header -------------------------------------------------------

! ----------------------------------------------------------------------
!  Section 1.   Initialize
! ----------------------------------------------------------------------
      LEQR(P1,P2) = ((ABS(P1-P2)) .LT. (1.E-6*ABS(P1+P2)))
      ErrorStatus = 0
      ErrorMessage= "  "

!  0.0  Set derived variables
      HalfPi         = Pi * 0.5
      TwoPi          = Pi * 2.0
      RecipDLambdaIn = 1.0 / DLambdaIn
      RecipDPhiIn    = 1.0 / DPhiIn
      HalfDPhiIn     = DPhiIn * 0.5

!  1.0  Error trap unsupported options.

!  1.1  Check that a valid interpolation scheme has been specified.
      If (( Scheme .ne. 1 .and.
     &      Scheme .ne. 3 .and.
     &      Scheme .ne. 5       ) ) then

       ErrorStatus=-1
       ErrorMessage="Invalid value of 'Scheme':linear interp. done  "

       VIorder = 1

      Else

       VIorder = Scheme

      End If

!  1.2  Check code matches Extend parameter
      IF (Extend .ne. 3) THEN

       ErrorStatus=1
       ErrorMessage="Don't change 'Extend' without amending code  "
       GOTO 999

      END IF

!  1.3  For global case check grid offset from south pole (PhiOrigin)
!       using the local function LEQR to test near equal reals
      IF (LGlobal) THEN
       Test1=-HalfPi
       Test2=-HalfPi+HalfDPhiIn
       IF(LEQR(PhiOrigin,Test1)) THEN
! origin at pole
        FirstOnPole = .True.
       ELSEIF(LEQR(PhiOrigin,Test2)) THEN
! origin half a grid length from pole
        FirstOnPole = .False.
       ELSE
! unexpected origin
       ErrorStatus=1
        ErrorMessage="HORINT1A: Incorrect grid  "
       GOTO 999
       END IF
      END IF

!  1.4   Get lower bound for DataExt
       LowerBound=1-Extend
!       Len1In  = SIZE (DataIn,  1)
!       Len2In  = SIZE (DataIn,  2)
!       Len1Out = SIZE (DataOut, 1)
!       Len2Out = SIZE (DataOut, 2)

!-----------------------------------------------------------------------
!  2.0 Extend input data array to bigger area to allow interpolation to
!       be done without having to redo any end points.
!-----------------------------------------------------------------------

      HalfLen1In = Len1In/2

!  2.1  Extend in i direction.
      IF (LGlobal) THEN                 !    2.1.1  Global
      DO j = 1, Len2In
       DataExt (-2,j)       = DataIn (Len1In-2,j)
       DataExt (-1,j)       = DataIn (Len1In-1,j)
       DataExt (0,j)        = DataIn (Len1In,j)
       DataExt (Len1In+1,j) = DataIn (1,j)
       DataExt (Len1In+2,j) = DataIn (2,j)
       DataExt (Len1In+3,j) = DataIn (3,j)
      END DO

      ELSE                              !    2.1.2  Limited Area.
      DO j = 1, Len2In
       DataExt (-2,j)       = DataIn (1,j)
       DataExt (-1,j)       = DataIn (1,j)
       DataExt (0,j)        = DataIn (1,j)
       DataExt (Len1In+1,j) = DataIn (Len1In,j)
       DataExt (Len1In+2,j) = DataIn (Len1In,j)
       DataExt (Len1In+3,j) = DataIn (Len1In,j)
      END DO

      END IF

      DO j = 1, Len2In
       DO i = 1, Len1In
        DataExt (i,j) = DataIn (i,j)
       END DO
      END DO

!  2.2  Extend in j direction.

      IF (LGlobal) THEN                 !    2.2.1  Global.
       SwapPole = 1
       IF (LVector)SwapPole = -1

        IF (.not.FirstOnPole) THEN
!  2.2.1.1.3  v points.
        DO i = -2, HalfLen1In+2
          DataExt (i,-2)       =SwapPole*DataExt(i+HalfLen1In,3)
          DataExt (i,-1)       =SwapPole*DataExt(i+HalfLen1In,2)
          DataExt (i,0)        =SwapPole*DataExt(i+HalfLen1In,1)
          DataExt (i,Len2In+1) =SwapPole*DataExt(i+HalfLen1In,Len2In)
          DataExt (i,Len2In+2) =SwapPole*DataExt(i+HalfLen1In,Len2In-1)
          DataExt (i,Len2In+3) =SwapPole*DataExt(i+HalfLen1In,Len2In-2)
        END DO

        DO i = HalfLen1In+3, Len1In+3
          DataExt (i,-2)       =SwapPole*DataExt(i-HalfLen1In,3)
          DataExt (i,-1)       =SwapPole*DataExt(i-HalfLen1In,2)
          DataExt (i,0)        =SwapPole*DataExt(i-HalfLen1In,1)
          DataExt (i,Len2In+1) =SwapPole*DataExt(i-HalfLen1In,Len2In)
          DataExt (i,Len2In+2) =SwapPole*DataExt(i-HalfLen1In,Len2In-1)
          DataExt (i,Len2In+3) =SwapPole*DataExt(i-HalfLen1In,Len2In-2)
        END DO

        ELSE
!  2.2.1.1.4  Other points.

        DO i = -2, HalfLen1In+2
          DataExt (i,-2)       =SwapPole*DataExt(i+HalfLen1In,4)
          DataExt (i,-1)       =SwapPole*DataExt(i+HalfLen1In,3)
          DataExt (i,0)        =SwapPole*DataExt(i+HalfLen1In,2)
          DataExt (i,Len2In+1) =SwapPole*DataExt(i+HalfLen1In,Len2In-1)
          DataExt (i,Len2In+2) =SwapPole*DataExt(i+HalfLen1In,Len2In-2)
          DataExt (i,Len2In+3) =SwapPole*DataExt(i+HalfLen1In,Len2In-3)
        END DO

        DO i = HalfLen1In+3, Len1In+3
          DataExt (i,-2)       =SwapPole*DataExt(i-HalfLen1In,4)
          DataExt (i,-1)       =SwapPole*DataExt(i-HalfLen1In,3)
          DataExt (i,0)        =SwapPole*DataExt(i-HalfLen1In,2)
          DataExt (i,Len2In+1) =SwapPole*DataExt(i-HalfLen1In,Len2In-1)
          DataExt (i,Len2In+2) =SwapPole*DataExt(i-HalfLen1In,Len2In-2)
          DataExt (i,Len2In+3) =SwapPole*DataExt(i-HalfLen1In,Len2In-3)
        END DO

        END IF                    ! END IF .not.FirstOnPole

      ELSE
!  2.2.2  Limited Area code.

      DO i = -2, Len1In+3
       DataExt (i,-2)       = DataExt (i,1)
       DataExt (i,-1)       = DataExt (i,1)
       DataExt (i,0)        = DataExt (i,1)
       DataExt (i,Len2In+1) = DataExt (i,Len2In)
       DataExt (i,Len2In+2) = DataExt (i,Len2In)
       DataExt (i,Len2In+3) = DataExt (i,Len2In)
      END DO

      END IF                    ! END IF LGlobal

!-----------------------------------------------------------------------
!  3.0  For each output point find i and j so that the point on the
!       output grid lies between i and i+1 and j and j+1.
!       and compute the weights
!-----------------------------------------------------------------------

!  3.1  Get output coordinates and weights where LAM straddling meridian

       IF (.NOT.LGlobal. AND.
     &      LambdaOrigin+len1In*DLambdaIn.GT.TwoPi) THEN

       DO j = 1, Len2Out
        DO i = 1, Len1Out
          WtLambda(i,j) = (LambdaOut(i,j) - LambdaOrigin)
     &                   * RecipDLambdaIn

          IF (LambdaOut(i,j).LT.LambdaOrigin)
     &        WtLambda(i,j) = WtLambda(i,j) + (TwoPi * RecipDLambdaIn) 

        IOut(i,j)     = INT(WtLambda(i,j) + 1)
        WtLambda(i,j) = WtLambda(i,j) - (IOut(i,j) - 1.0)

        WtPhi(i,j) = (PhiOut(i,j) - PhiOrigin) * RecipDPhiIn
        JOut(i,j)     = INT(WtPhi(i,j) + 1)
        WtPhi(i,j) = WtPhi(i,j) - (JOut(i,j) - 1.0)

        END DO
       END DO

!  3.2  Get output coordinates and weights for normal case

       ELSE 
        DO j = 1, Len2Out                                               
         DO i = 1, Len1Out                                              
          WtLambda(i,j) = (LambdaOut(i,j) - LambdaOrigin)
     &                   * RecipDLambdaIn
          IOut(i,j)     = INT(WtLambda(i,j) + 1)                        
          WtLambda(i,j) = WtLambda(i,j) - (IOut(i,j) - 1.0)             
                                                                        
          WtPhi(i,j)    = (PhiOut(i,j) - PhiOrigin) * RecipDPhiIn       
          JOut(i,j)     = INT(WtPhi(i,j) + 1)                           
          WtPhi(i,j)    = WtPhi(i,j) - (JOut(i,j) - 1.0)                
                                                                        
         END DO                                                         
        END DO                                                          
       END IF       

!-----------------------------------------------------------------------
!  4.0    Perform required interpolations.
!-----------------------------------------------------------------------

!  4.1  Call the specified interpolation scheme:
      IF(VIorder.eq.1)Then
        CALL HorizontalInterpLinear (LowerBound,
     &                               Len1In,  Len2In,
     &                               Len1Out, Len2Out,
     &                               DataExt,
     &                               WtLambda,
     &                               WtPhi,
     &                               IOut,
     &                               JOut,
     &                               DataOut)

      ELSEIF(VIorder.eq.3)Then 
        CALL HorizontalInterpCubic (LowerBound,
     &                               Len1In,  Len2In,
     &                               Len1Out, Len2Out,
     &                               DataExt,
     &                               WtLambda,
     &                               WtPhi,
     &                               IOut,
     &                               JOut,
     &                               DataOut)


      ELSEIF(VIorder.eq.5)Then
        CALL HorizontalInterpQuintic (LowerBound,
     &                               Len1In,  Len2In,
     &                               Len1Out, Len2Out,
     &                               DataExt,
     &                               WtLambda,
     &                               WtPhi,
     &                               IOut,
     &                               JOut,
     &                               DataOut)
      EndIF

!-----------------------------------------------------------------------
!  5.0  Perform montonicity enforcement if high order scheme used and
!       monotonicity required.
! section commented out, montonicity code not in vertical interpolation
! and as such is currently disabled here (so adjoint not calculated)
!-----------------------------------------------------------------------

      IF(L_mono.AND.VIorder.gt.1)THEN
        CALL HorizontalInterpLinear (LowerBound,
     &                               Len1In,  Len2In,
     &                               Len1Out, Len2Out,
     &                               DataExt,
     &                               WtLambda,
     &                               WtPhi,
     &                               IOut,
     &                               JOut,
     &                               DataOut)


      DO j = 1, Len2Out
       DO i = 1, Len1Out
        DataHigh (i,j) = DataOut (i,j)
       END DO
      END DO

        CALL HorizontalInterpMonotone
     &                        (LowerBound,
     &                         Len1In,  Len2In,
     &                         Len1Out, Len2Out,
     &                         DataExt,
     &                         DataHigh,
     &                         DataMono,
     &                         IOut,
     &                         JOut,
     &                         DataOut)

      END IF


! End of routine.
999   CONTINUE
      RETURN
      END
