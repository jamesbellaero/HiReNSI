
! =========================== Modifications and Adjustments to my code ==============================================================================

! Allocating additional required arrays -------------------------------------------------------------------------------------------------------------
Allocate ( NGP(NEL), ID_BC(NEL,NDIM**2), LTRANS(100), PML_PARAM(NDIM * 2 , 4), STAT = ERR_Alloc) ;
  IF ( ERR_Alloc /= 0 ) Then ;
    Write (*, Fmt_ALLCT) ERR_Alloc ;  Write (UnInf, Fmt_ALLCT) ERR_Alloc ;
    !#Call BEEP_FAIL ;
    Write(*, Fmt_FL) ;  Write(UnInf, Fmt_FL) ; Write(*, Fmt_End) ; Read(*,*) ;  STOP ;
  End If ;


! We only use structured mesh; hence, NNode is constant throughout: ---------------------------------------------------------------------------------
NNode = MaxNNode


! Babak reads PMat in this order: 1.lambda 2.mu 3.rho 4.alpha_0 5.beta_0 6.m 7.? 8.ratio ------------------------------------------------------------
! We use 1.E 2.anu instead:

DO IMAT = 1 , NMAT
   temp_LAMBDA = PMAT ( IMAT, 1 )
   temp_MU     = PMAT ( IMAT, 2 )
   PMAT ( IMAT, 1 ) = temp_MU * ( 3.0d0 * temp_LAMBDA + 2.0d0 * temp_MU ) / ( temp_LAMBDA + temp_MU )
   PMAT ( IMAT, 2 ) = temp_LAMBDA / ( 2.0d0 * temp_LAMBDA + 2.0d0 * temp_MU )
   PMAT ( IMAT, 3 ) = 2.0d-3
   PMAT ( IMAT, 5 ) = 1.0d0
END DO


! Gauss points --------------------------------------------------------------------------------------------------------------------------------------
NGP = 3
NINT_Lobatto = 3


! PML parameters ------------------------------------------------------------------------------------------------------------------------------------
! various conventions:
! Babak: X LEFT	   L LEFT
!        X RIGHT	   L RIGHT
!        Y BOT	   L BOT
!        Y TOP	   L TOP

! Mine: [ alpha0 beta0 L s0]
!Allocate (PML_PARAM(NDIM * 2 , 4))

PML_PARAM = 0.0d0

Select Case (NDIM)

   Case (2)
      PML_PARAM (1,3) = PML_DIM(2,2)
      PML_PARAM (2,3) = PML_DIM(1,2)
      PML_PARAM (3,3) = PML_DIM(4,2)
      PML_PARAM (4,3) = PML_DIM(3,2)
      PML_PARAM (1,4) = PML_DIM(2,1)
      PML_PARAM (2,4) = PML_DIM(1,1)
      PML_PARAM (3,4) = PML_DIM(4,1)
      PML_PARAM (4,4) = PML_DIM(3,1)

      PML_PARAM (:,1) = 0.750d0
      PML_PARAM (:,2) = 100.0d0

   case (3)
      PML_PARAM (1,3) = PML_DIM(2,2)
      PML_PARAM (2,3) = PML_DIM(1,2)
      PML_PARAM (3,3) = PML_DIM(4,2)
      PML_PARAM (4,3) = PML_DIM(3,2)
      PML_PARAM (5,3) = PML_DIM(6,2)
      PML_PARAM (6,3) = PML_DIM(5,2)
      PML_PARAM (1,4) = PML_DIM(2,1)
      PML_PARAM (2,4) = PML_DIM(1,1)
      PML_PARAM (3,4) = PML_DIM(4,1)
      PML_PARAM (4,4) = PML_DIM(3,1)
      PML_PARAM (5,4) = PML_DIM(6,1)
      PML_PARAM (6,4) = PML_DIM(5,1)

      PML_PARAM (:,1) =   5.0d0
      PML_PARAM (:,2) = 400.0d0

End Select

! ID_BC ---------------------------------------------------------------------------------------------------------------------------------------------
ID_BC = 0

DO I = 1 , Param%IntM( 2, 4) ! NIDBC
   IEL = IDBC ( I , 1 )
   ID_BC ( IEL , 1 ) = IDBC ( I , 2 )
END DO

! time-stepping -------------------------------------------------------------------------------------------------------------------------------------
Select Case (NDim)

   Case (2)
      nstep = 200
      ALPHA = 0.2525062500d0
      DELTA = 0.505000d0
      dt    = 0.00010d0
      ltrans= 0
      ltrans(1) = 402

   Case (3)
      ltrans = 0
      ltrans(1) = 818
! inversion EX06
!      nstep  = 350
!      dt     = 6.d-3
! -----------------------------------------------
! inversion EX08
!      nstep  = 400
!      dt     = 0.001d0
! inversion EX08R
!      nstep  = 2000 * 2.0d0
!      dt     = 0.001d0 / 2.0d0
! -----------------------------------------------
! inversion EX11
!      nstep  = 390
!      dt     = 9.0d-4
! inversion EX11R
!      nstep  = 2224
!      dt     = 9.0d-4 / 2.0d0
! -----------------------------------------------
! inversion EX14
!      nstep  = 400
!      dt     = 1.0d-3
! -----------------------------------------------
! inversion EX16
      nstep  = 400     !450
      dt     = 1.0d-3  !1.0d-3

End Select



!====================================================================================================================================================
!  Inversion parameters
!====================================================================================================================================================

! 1- regularization factor continuation (adaptive regularization factor - See Sez p. 171)
  Regularization_factor_continuation = 0.3d0

! 2- maximum allowable update at each step
  Fac_Max_Update_Lambda = 5.0d0
  Fac_Max_Update_Mu     = 5.0d0

! 3- lower and upper bounds on Lambda & Mu
  bound_Lambda_upper    = 166.0d0
  bound_Mu_upper        = 166.0d0

  bound_Lambda_lower    = 50.0d0
  bound_Mu_lower        = 50.0d0

! 4- iterations for inversion
  Iter_begin            = 641
  Iter_end              = 10000
  MaxIter               = 10000

! 5- Restart Conjugate Gradient due to round-off error
  N_CG_Restart          = 10

! 6- Regularization method
  Regularization_Method = 'TV'                   ! TV | TN

! 7- Inversion (control) parameter. {Lambda} may mean {Lambda + 2 Mu}. See 9. ***********************************************************************
  Control_Parameter     = 'Lambda_Mu'                   ! Mu | Lambda | Lambda_Mu
!  Control_Parameter     = 'Lambda' !Output2                  ! Mu | Lambda | Lambda_Mu

! 8- Biasing of Lambda search direction based on Mu
  Bias_Lambda_search    = 'Yes'                   ! Yes | No
  NIterBias             = 50

! 9- Improve the insensitivity of the cost functional by using {Lambda + 2 Mu} and {Mu} instead of {Lambda} and {Mu}.
  Lambda_variable       = 'Lambda'               ! Lambda_2Mu | Lambda

! 10- Weighted-regularization scheme
  w_top                   = 50.0d0
  w_bot                   = 1.0d0
  w_Length                = 45.0d0
  Weighted_Regularization = 'No'                 ! Yes | No

! 11- L-BFGS number of vectors
  M_LBFGS               = 15

! 12- Search-direction method 
  Search_Direction_Method = 'LBFGS'              ! SD | CG | LBFGS

!====================================================================================================================================================
!
!====================================================================================================================================================

! modify material bounds if using {Lambda + 2 Mu} instead of {Lambda}
  If ( Lambda_variable == 'Lambda_2Mu' ) Then
     bound_Lambda_upper = bound_Lambda_upper + 2.0d0 * bound_Mu_upper
     bound_Lambda_lower = bound_Lambda_lower + 2.0d0 * bound_Mu_lower
  End If










! ===================================================================================================================================================
IF (RANK == 0 .AND. I_PRINT_SCREEN ) write(*,*) 'end of adjustments'



