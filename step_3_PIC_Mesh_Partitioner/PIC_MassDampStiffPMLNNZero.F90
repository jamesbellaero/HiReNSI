
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                                                                                                                  ++
! Last Update:  02 Sep 2011                                                                                                                        ++
!                                                                                                                                                  ++
! Description: This Module obtains non-zero entries of elemental stiffness, mass and damping matrices using in PETSc.                              ++
!                                                                                                                                                  ++
!                                                                                                                                                  ++
!              Subroutines                     Description                                                                                         ++
!              ==============                  ==============                                                                                      ++
!              MassDampStiffPML_2D_8N_Pattern       ELEMENT NUMBER IS 4 (PLANE STRAIN)                                                             ++
!              MassDampStiffSLD_2D_8N_Pattern       ELEMENT NUMBERS ARE 2 OR 4 (PLANE STRAIN)                                                      ++
!                                                                                                                                                  ++
!                                                                                                                                                  ++
!                                                                                                                                                  ++
!                                                                                                                                                  ++
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Module PML_Matrices_Pattern ;

Use Parameters ;
Use PML_PRMTRS ;
Use ShapeFunctions ;

Implicit None ;

  Interface ;
!    Module Procedure 
  End Interface    ;

Contains ;


! ***************************************************************************************************************************************************
!                                                                                                                                                  **
! Last Update:  02 Sep 2011                                                                                                                        **
! Description: THIS Subroutine CALCULATES MASS, DAMPING AND STIFFNESS MATRICES OF REGULAR DOMAIN FOR A PML TRUNCATED DOMAIN.                       **
!              ELENENT NUMBER IS 4 (PLANE STRAIN)                                                                                                  **
!                                                                                                                                                  **
! ***************************************************************************************************************************************************

!Subroutine MassDampStiffSLD_2D_8N_Pattern (    &
!IEL, NNode, NDim, NInt,                                                                                                         & ! Integer Variables
!Rho, Lambda, MU, BETTA, alpha,                                                                                                   & ! Real Variables
!INOD,                                                                                                                           & ! Integer Arrays
!XYZ,     KE, ME, CE,                                                                                                            & ! Real Arrays
!!                                                                                                                               & ! Characters
!GAUSS_PNT                                                                                                                       & ! Type 
!) ;

Subroutine MassDampStiffSLD_2D_8N_Pattern (                                                                                  &
IEl, NNode, NDim, NInt,                                                                                                              & ! Integer Variables
Rho, Lambda, MU, BETTA, alpha,                                                                                                   & ! Real Variables
!                                                                                                                               & ! Integer Arrays
KE, ME, CE,                                                                                                                     & ! Real Arrays
!                                                                                                                               & ! Characters
GAUSS_PNT                                                                                                                       & ! Type 
) ;

Implicit None ;

!?? <<<<<  fix the variables based on the new subroutine

! =========================== Global Variables ======================================================================================================

! - Integer Arrays ----------------------------------------------------------------------------------------------------------------------------------
Integer (Kind=Tiny), Intent(In)    :: NNode, NDim, NInt ;
Integer (Kind=Lng ), Intent(In)    :: IEL ;

! - Real Variables ----------------------------------------------------------------------------------------------------------------------------------
Real (Kind=DBL), Intent(In)    :: Rho, Lambda, MU, BETTA, alpha ;

! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
!#Complex, Intent(In)    ::  ;

! - Integer Arrays ----------------------------------------------------------------------------------------------------------------------------------
!Integer (Kind=Lng ), Intent(In), Dimension (:,:)  :: INOD ;

! - Real Arrays -------------------------------------------------------------------------------------------------------------------------------------
!Real (Kind=DBL), Intent(IN   ), Dimension (:,:)  :: XYZ ;
Real (Kind=DBL), Intent(InOut), Dimension (:,:)  :: KE, ME, CE ;

! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
!#Complex, Intent(In)    ::  ;
!#Complex, Intent(InOut) ::  ;
!#Complex, Intent(OUT)   ::  ;
! - Character Variables -----------------------------------------------------------------------------------------------------------------------------
!#Character   ::  ;
! - Logical Variables -------------------------------------------------------------------------------------------------------------------------------
!#Logical   ::  ;

! =========================== LOCAL Variables =======================================================================================================

! - Integer Arrays ------------------------------------------------------------------------------------------------------------------------------
Integer (Kind=Shrt)  :: I, J, K, LX, LY ;

! - Real Variables ----------------------------------------------------------------------------------------------------------------------------------
Real (Kind=DBL)     ::  FAC, WX, WY, WSTAR, DETJ ;

! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
!#Complex  ::  ;
! - Integer Arrays ----------------------------------------------------------------------------------------------------------------------------------
!#Integer (Kind=Shrt)  ::  ;

! - Real Arrays -------------------------------------------------------------------------------------------------------------------------------------
Real (Kind=DBL)  :: XT ( NDim, NNode ), DJ ( NDim, NDim ), DJI ( NDim, NDim ), DFX ( NNode, NDim ), PHI_PHI_T ( NNode, NNode ), PHIX_PHIX_T ( NNode, NNode ), PHIY_PHIY_T ( NNode, NNode ), PHIX_PHIY_T ( NNode, NNode ), PHIY_PHIX_T ( NNode, NNode ) ;
!Real (Kind=DBL)  :: DFX ( NNode, NDim ), PHI_PHI_T ( NNode, NNode ), PHIX_PHIX_T ( NNode, NNode ), PHIY_PHIY_T ( NNode, NNode ), PHIX_PHIY_T ( NNode, NNode ), PHIY_PHIX_T ( NNode, NNode ) ;
! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
!#Complex  ::  ;
! - Character Variables -----------------------------------------------------------------------------------------------------------------------------
!#Character   ::  ;
! - Logical Variables -------------------------------------------------------------------------------------------------------------------------------
!#Logical   ::  ;

! - Type DECLERATIONS -------------------------------------------------------------------------------------------------------------------------------
Type (  SF_2_8 ) ::  SF       ;  ! SHAPE FUNCTION
Type ( DSF_2_8 ) :: DSF       ;  ! DIFFERENTIALS OF SHAPE FUNCTION
Type ( GAUSS   ) :: GAUSS_PNT ;  ! GAUSS POINTS 

! =========================== Subroutine CODE =======================================================================================================

! COORDINATES OF THE ELEMENT
!ForAll ( I = 1:NDim, J = 1:NNode ) XT ( I, J ) = XYZ ( INOD(J,IEL), I ) ;

  ! INTEGRATE OVER EACH GAUSSIAN POINT
  DO LY = 1, NInt ;
    SF%X2  = GAUSS_PNT%XINT ( LY ) ;
    DSF%X2 = GAUSS_PNT%XINT ( LY ) ;
    WY     = GAUSS_PNT%WINT ( LY ) ;    !%%

      DO LX = 1, NInt ;
        SF%X1  = GAUSS_PNT%XINT ( LX ) ;
        DSF%X1 = GAUSS_PNT%XINT ( LX ) ;
        WX     = GAUSS_PNT%WINT ( LX ) ;    !%%

        WSTAR  = WX * WY ;    !%%

        ! SHAPE FUNCTIONS AND DIFFERENTIAL OF SHAPE FUNCTIONS
        Call     ShapeFuncSub (  SF ) ;
        Call DIF_ShapeFuncSub ( DSF ) ;  

        DJ   = MATMUL ( XT, DSF%DFXI ) ;    !%%
        DETJ = DJ ( 1, 1 ) * DJ ( 2, 2 ) - DJ ( 2, 1 ) * DJ ( 1, 2 ) ;    !%%
        FAC  = WSTAR * DETJ ;    !%%

!          IF ( DETJ <= 0.0_DBL ) Then ;
!            Write(*,"('ELEMENT NUMBER',3I6,'|J|<0  -ERROR-', 'DETJ = ',E15.10)" ) IEL, LX, LY, DETJ ;
!            Write(*, Fmt_End) ; Read(*,*) ;  STOP ;
!          End If ;

        ! CALCULATING THE INVERSE OF THE JACOBIAN
        !#Call DLINRG (NDim,DJ,NDim,DJI,NDim) ;
        DJI ( 1, 1 ) =   DJ ( 2, 2 ) ;    !%%
        DJI ( 2, 2 ) =   DJ ( 1, 1 ) ;    !%%
        DJI ( 1, 2 ) = - DJ ( 1, 2 ) ;    !%%
        DJI ( 2, 1 ) = - DJ ( 2, 1 ) ;    !%%

        !DJI = DJI  / DETJ ;    !%%

        DFX = MATMUL ( DSF%DFXI, DJI ) ;    !%%
        !DFX = DSF%DFXI ;

        !#PHI_PHI_T   = MATMUL ( SF%FN, TRANSPOSE (SF%FN) )                * FAC ;
        !#PHIX_PHIX_T = MATMUL ( DFX (:, 1), TRANSPOSE ( DFX ( : , 1 ) ) ) * FAC ;
        !#PHIY_PHIY_T = MATMUL ( DFX (:, 1), TRANSPOSE () )                * FAC ;
        !#PHIX_PHIY_T = MATMUL ( DFX (:, 1), TRANSPOSE (DFX ( : , 2 )) )   * FAC ;
        !#PHIY_PHIX_T = MATMUL ( DFX (:, 2), TRANSPOSE (DFX ( : , 1 )) )   * FAC ;

          DO I = 1, NNode   ;
            DO J = 1, NNode ;
              PHI_PHI_T   ( I, J ) = SF%FN ( I )     * SF%FN ( J )     * FAC ;       !%% *fac
              PHIX_PHIX_T ( I, J ) = DFX   ( I , 1 ) * DFX   ( J , 1 ) * FAC ;
              PHIY_PHIY_T ( I, J ) = DFX   ( I , 2 ) * DFX   ( J , 2 ) * FAC ;
              PHIX_PHIY_T ( I, J ) = DFX   ( I , 1 ) * DFX   ( J , 2 ) * FAC ;
              PHIY_PHIX_T ( I, J ) = DFX   ( I , 2 ) * DFX   ( J , 1 ) * FAC ;
            End Do ;
          End Do ;

        ! MASS MATRIX
        ForAll ( I = 1:NNode, J = 1:NNode, K = 1:NDim )   ME ( ( K -1 ) * NNode + I, ( K -1 ) * NNode + J ) = ME ( ( K -1 ) * NNode + I, ( K -1 ) * NNode + J ) + Rho * PHI_PHI_T ( I, J ) ;

        ! STIFFNESS MATRIX
          !ForAll ( I = 1:NNode, J = 1:NNode ) ;           !?? CHECK IT OUT
          !  S (         I,         J ) = S (         I,         J ) + ( Lambda + 2 * MU ) * PHIX_PHIX_T + MU * PHIY_PHIY_T ;
          !  S ( NNode + I, NNode + J ) = S ( NNode + I, NNode + J ) + ( Lambda + 2 * MU ) * PHIY_PHIY_T + MU * PHIX_PHIX_T ;
          !  S (         I, NNode + J ) = S (         I, NNode + J ) +   Lambda            * PHIX_PHIY_T + MU * PHIY_PHIX_T ;
          !  S ( NNode + I,         J ) = S ( NNode + I,         J ) +   Lambda            * PHIY_PHIX_T + MU * PHIX_PHIY_T ;
          !End ForAll ;

        KE ( 1:NNode          , 1:NNode           ) = KE ( 1:NNode, 1:NNode                     ) + ( Lambda + 2 * MU ) * PHIX_PHIX_T + MU * PHIY_PHIY_T  ;
        KE ( NNode + 1:2*NNode, NNode + 1:2*NNode ) = KE ( NNode + 1:2*NNode, NNode + 1:2*NNode ) + ( Lambda + 2 * MU ) * PHIY_PHIY_T + MU * PHIX_PHIX_T  ;
        KE ( 1:NNode          , NNode + 1:2*NNode ) = KE ( 1:NNode          , NNode + 1:2*NNode ) +   Lambda            * PHIX_PHIY_T + MU * PHIY_PHIX_T  ;  !?? fix it
        KE ( NNode + 1:2*NNode, 1:NNode           ) = KE ( NNode + 1:2*NNode, 1:NNode           ) +   Lambda            * PHIY_PHIX_T + MU * PHIX_PHIY_T  ;

      End Do ;
  End Do ;

! DAMPING MATRIX
CE = BETTA * KE + alpha * ME ; 

!#Write(*     ,*) 'End Subroutine < MassDampStiffSLD_2D_8N_Pattern >' ;
!#Write(UnInf,*) 'End Subroutine < MassDampStiffSLD_2D_8N_Pattern >' ;
Return ;
End Subroutine MassDampStiffSLD_2D_8N_Pattern ;


! ***************************************************************************************************************************************************
!                                                                                                                                                  **
! Last Update:  02 Sep 2011                                                                                                                        **
! Description: THIS Subroutine CALCULATES ALL REQUIRED MATRICES FOR PML TRUNCATED BOUNDARIES.                                                      **
!              ELENENT NUMBER IS 4 (PLANE STRAIN)                                                                                                  **
!                                                                                                                                                  **
! ***************************************************************************************************************************************************

!Subroutine MassDampStiffPML_2D_8N_Pattern (                                                                                &
!NDim, NNode, NInt, PARAM_Type,      IEL,                                                                                      & ! Integer Variables
!Lambda, MU, Rho, R_PML_alpha_0, C_REF_BETA_0, M_PML, B_PML,                                                                    & ! Real Variables
!INOD,                                                                                                                         & ! Integer Arrays
!XYZ,     KE, ME, CE, PML_DIM,                                                                                                 & ! Real Arrays
!!                                                                                                                             & ! Characters
!GAUSS_PNT                                                                                                                     & ! Type
!) ;

Subroutine MassDampStiffPML_2D_8N_Pattern (                                                                                &
NDim, NNode, NInt, PARAM_Type, IEL,                                                                                           & ! Integer Variables
Lambda, MU, Rho, R_PML_alpha_0, C_REF_BETA_0, M_PML, B_PML,                                                                    & ! Real Variables
INOD,                                                                                                                         & ! Integer Arrays
XYZ,     KE, ME, CE, PML_DIM,                                                                                                 & ! Real Arrays
!                                                                                                                             & ! Characters
GAUSS_PNT                                                                                                                     & ! Type
) ;

Implicit None ;

!?? <<<<<  fix the variables based on the new subroutine

! =========================== Global Variables ======================================================================================================

! - Integer Arrays ------------------------------------------------------------------------------------------------------------------------------
Integer (Kind=Tiny), Intent(In)    :: NDim, NNode, NInt, PARAM_Type ;
Integer (Kind=Lng ), Intent(In)    :: IEL ;
!#Integer (Kind=Shrt), Intent(InOut) ::  ;
!#Integer (Kind=Shrt), Intent(OUT)   ::  ;

! - Real Variables ----------------------------------------------------------------------------------------------------------------------------------
Real (Kind=DBL), Intent(In)    :: Lambda, MU, Rho, R_PML_alpha_0, C_REF_BETA_0, M_PML, B_PML ;
!#Real (Kind=DBL), Intent(InOut) ::  ;
!#Real (Kind=DBL), Intent(OUT)   ::  ;
! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
!#Complex, Intent(In)    ::  ;
!#Complex, Intent(InOut) ::  ;
!#Complex, Intent(OUT)   ::  ;

! - Integer Arrays ----------------------------------------------------------------------------------------------------------------------------------
!#Integer (Kind=Shrt), Intent(In), Dimension (:)  ::
Integer (Kind=Lng ), Intent(In), Dimension (:,:)  :: INOD ;
!#Integer (Kind=Shrt), Intent(InOut) ::  ;
!#Integer (Kind=Shrt), Intent(OUT)   ::  ;

! - Real Arrays -------------------------------------------------------------------------------------------------------------------------------------
!#Real (Kind=DBL), Intent(In), Dimension (:,:)  ::  
Real (Kind=DBL), Intent(IN   ), Dimension (:,:)  :: XYZ, PML_DIM ;
Real (Kind=DBL), Intent(InOut), Dimension (:,:)  :: KE, ME, CE ;
!#Real (Kind=DBL), Intent(OUT)   ::  ;

  ! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
!#Complex, Intent(In)    ::  ;
!#Complex, Intent(InOut) ::  ;
!#Complex, Intent(OUT)   ::  ;
! - Character Variables -----------------------------------------------------------------------------------------------------------------------------
!#Character   ::  ;
! - Logical Variables -------------------------------------------------------------------------------------------------------------------------------
!#Logical   ::  ;
! =========================== LOCAL Variables =======================================================================================================

! - Integer Arrays ------------------------------------------------------------------------------------------------------------------------------
Integer (Kind=Shrt)  :: I, J, K, LX, LY ;

! - Real Variables ----------------------------------------------------------------------------------------------------------------------------------
Real (Kind=DBL)     ::  FAC, WX, WY, WSTAR, DETJ ;
Real (Kind=DBL)     ::  FAC_A_Rho, FAC_B_Rho, FAC_C_Rho, FAC_A_1, FAC_A_2, FAC_A_3, FAC_B_1, FAC_B_2, FAC_B_3, FAC_C_1, FAC_C_2, FAC_C_3, Fac_A, Fac_B, Fac_C, alpha_X, alpha_Y, BETA_X, BETA_Y ;  ! PML Parameters
! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
!#Complex  ::  ;
! - Integer Arrays ----------------------------------------------------------------------------------------------------------------------------------
!#Integer (Kind=Shrt)  ::  ;

! - Real Arrays -------------------------------------------------------------------------------------------------------------------------------------
Real (Kind=DBL)  :: XT ( NDim, NNode ), DJ ( NDim, NDim ), DJI ( NDim, NDim ), DFX ( NNode, NDim ), PHI_PHI_T ( NNode, NNode ), PSI_PSI_T ( NNode, NNode ), PHIX_PSI_T ( NNode, NNode ), PHIY_PSI_T ( NNode, NNode ), PHIX_PHIX_T ( NNode, NNode ), PHIY_PHIY_T ( NNode, NNode ), PHIX_PHIY_T ( NNode, NNode ), PHIY_PHIX_T ( NNode, NNode ), XG ( NDim ) ;

! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
!#Complex  ::  ;
! - Character Variables -----------------------------------------------------------------------------------------------------------------------------
!#Character   ::  ;
! - Logical Variables -------------------------------------------------------------------------------------------------------------------------------
!#Logical   ::  ;

! - Type DECLERATIONS -------------------------------------------------------------------------------------------------------------------------------
Type (  SF_2_8 ) ::  SF       ;  ! SHAPE FUNCTION
Type ( DSF_2_8 ) :: DSF       ;  ! DIFFERENTIALS OF SHAPE FUNCTION
Type ( GAUSS   ) :: GAUSS_PNT ;  ! GAUSS POINTS 

! =========================== Subroutine CODE =======================================================================================================

! COORDINATES OF THE ELEMENT
ForAll ( I = 1:NDim, J = 1:NNode ) XT ( I, J ) = XYZ ( INOD( J, IEL ), I ) ;

  ! INTEGRATE OVER EACH GAUSSIAN POINT
  DO LY = 1, NInt ;
    SF%X2  = GAUSS_PNT%XINT ( LY ) ;
    DSF%X2 = GAUSS_PNT%XINT ( LY ) ;
   WY     = GAUSS_PNT%WINT ( LY ) ;         !%%

      DO LX = 1, NInt ;
        SF%X1  = GAUSS_PNT%XINT ( LX ) ;
        DSF%X1 = GAUSS_PNT%XINT ( LX ) ;
        WX     = GAUSS_PNT%WINT ( LX ) ;         !%%

        WSTAR  = WX * WY ;         !%%

        ! SHAPE FUNCTIONS AND DIFFERENTIAL OF SHAPE FUNCTIONS
        Call     ShapeFuncSub (  SF ) ;
        Call DIF_ShapeFuncSub ( DSF ) ;  

        DJ   = MATMUL ( XT, DSF%DFXI ) ;         !%%
        DETJ = DJ ( 1, 1 ) * DJ ( 2, 2 ) - DJ ( 2, 1 ) * DJ ( 1, 2 ) ;         !%%
        FAC  = WSTAR * DETJ ;         !%%

!          IF ( DETJ <= 0.0_DBL ) Then ;
!            Write(*,"('ELEMENT NUMBER',3I6,'|J|<0  -ERROR-', 'DETJ = ',E15.10)" ) IEL, LX, LY, DETJ ;
!            Write(*, Fmt_End) ; Read(*,*) ;  STOP ;
!          End If ;

        ! CALCULATING THE INVERSE OF THE JACOBIAN
        !#Call DLINRG (NDim,DJ,NDim,DJI,NDim) ;
        DJI ( 1, 1 ) =   DJ ( 2, 2 ) ;         !%%
        DJI ( 2, 2 ) =   DJ ( 1, 1 ) ;
        DJI ( 1, 2 ) = - DJ ( 1, 2 ) ;
        DJI ( 2, 1 ) = - DJ ( 2, 1 ) ;

        DJI = DJI  / DETJ ;         !%%

        DFX = MATMUL ( DSF%DFXI, DJI ) ;         !%%
!        DFX = DSF%DFXI

        !#PHI_PHI_T   = MATMUL ( SF%FN, TRANSPOSE (SF%FN) )                * FAC ;
        !#PSI_PSI_T   = MATMUL ( SF%FN, TRANSPOSE (SF%FN) )                * FAC ;
        !#PHIX_PSI_T  = MATMUL ( DFX (:, 1), TRANSPOSE (SF%FN) )           * FAC ;
        !#PHIY_PSI_T  = MATMUL ( DFX (:, 1), TRANSPOSE (SF%FN) )           * FAC;
        !#PHIX_PHIX_T = MATMUL ( DFX (:, 1), TRANSPOSE ( DFX ( : , 1 ) ) ) * FAC ;
        !#PHIY_PHIY_T = MATMUL ( DFX (:, 1), TRANSPOSE () )                * FAC ;
        !#PHIX_PHIY_T = MATMUL ( DFX (:, 1), TRANSPOSE (DFX ( : , 2 )) )   * FAC ;
        !#PHIY_PHIX_T = MATMUL ( DFX (:, 2), TRANSPOSE (DFX ( : , 1 )) )   * FAC ;

          DO I = 1, NNode   ;
            DO J = 1, NNode ;
              PHI_PHI_T   ( I, J ) = SF%FN ( I )     * SF%FN ( J )     * FAC ;
              PSI_PSI_T   ( I, J ) = SF%FN ( I )     * SF%FN ( J )     * FAC ;
              PHIX_PSI_T  ( I, J ) = DFX   ( I , 1 ) * SF%FN ( J )     * FAC ;
              PHIY_PSI_T  ( I, J ) = DFX   ( I , 2 ) * SF%FN ( J )     * FAC ;
              PHIX_PHIX_T ( I, J ) = DFX   ( I , 1 ) * DFX   ( J , 1 ) * FAC ;
              PHIY_PHIY_T ( I, J ) = DFX   ( I , 2 ) * DFX   ( J , 2 ) * FAC ;
              PHIX_PHIY_T ( I, J ) = DFX   ( I , 1 ) * DFX   ( J , 2 ) * FAC ;
              PHIY_PHIX_T ( I, J ) = DFX   ( I , 2 ) * DFX   ( J , 1 ) * FAC ;
            End Do ;
          End Do ;

        XG = MATMUL ( XT , SF%FN ) ;  !?? CHECK THIS ONE

        ! PML Parameters
        Call PML_Parameters_2D  ( NDim, PARAM_Type, IEL,     Lambda, MU, Rho,  R_PML_alpha_0, C_REF_BETA_0, M_PML, B_PML,     FAC_A_Rho, FAC_B_Rho, FAC_C_Rho, FAC_A_1, FAC_A_2, FAC_A_3, FAC_B_1, FAC_B_2, FAC_B_3, FAC_C_1, FAC_C_2, FAC_C_3, alpha_X, alpha_Y, BETA_X, BETA_Y, Fac_A, Fac_B, Fac_C,     XG, PML_DIM ) ;

          ForAll ( I = 1:NNode, J = 1:NNode ) ;

            ! MASS MATRIX
            ForAll ( K = 1:NDim             ) ME ( ( K - 1        ) * NNode + I, ( K - 1        ) * NNode + J ) = ME ( ( K - 1        ) * NNode + I, ( K - 1        ) * NNode + J ) +     PHI_PHI_T  ( I, J ) * FAC_A_Rho ;
            ForAll ( K = NDim + 1: 2 * NDim ) ME ( ( K - 1        ) * NNode + I, ( K - 1        ) * NNode + J ) = ME ( ( K - 1        ) * NNode + I, ( K - 1        ) * NNode + J ) +     PSI_PSI_T  ( I, J ) * FAC_A_1   ;
                                              ME (       NDim       * NNode + I, (     NDim + 1 ) * NNode + J ) = ME (       NDim       * NNode + I, (     NDim + 1 ) * NNode + J ) +     PSI_PSI_T  ( I, J ) * FAC_A_2   ;
                                              ME ( (     NDim + 1 ) * NNode + I,       NDim       * NNode + J ) = ME ( (     NDim + 1 ) * NNode + I,       NDim       * NNode + J ) +     PSI_PSI_T  ( I, J ) * FAC_A_2   ;
                                              ME ( ( 2 * NDim     ) * NNode + I, ( 2 * NDim     ) * NNode + J ) = ME ( ( 2 * NDim     ) * NNode + I, ( 2 * NDim     ) * NNode + J ) +     PSI_PSI_T  ( I, J ) * FAC_A_3   ;

            ! STIFFNESS MATRIX
            ForAll ( K = 1:NDim             ) KE ( ( K - 1        ) * NNode + I, ( K - 1        ) * NNode + J ) = KE ( ( K - 1        ) * NNode + I, ( K - 1        ) * NNode + J ) +     PHI_PHI_T  ( I, J ) * FAC_C_Rho ;
            ForAll ( K = NDim + 1: 2 * NDim ) KE ( ( K - 1        ) * NNode + I, ( K - 1        ) * NNode + J ) = KE ( ( K - 1        ) * NNode + I, ( K - 1        ) * NNode + J ) +     PSI_PSI_T  ( I, J ) * FAC_C_1   ;

                                              KE (       NDim       * NNode + I, (     NDim + 1 ) * NNode + J ) = KE (       NDim       * NNode + I, (     NDim + 1 ) * NNode + J ) +     PSI_PSI_T  ( I, J ) * FAC_C_2   ;
                                              KE ( (     NDim + 1 ) * NNode + I,       NDim       * NNode + J ) = KE ( (     NDim + 1 ) * NNode + I,       NDim       * NNode + J ) +     PSI_PSI_T  ( I, J ) * FAC_C_2   ;
                                              KE ( ( 2 * NDim     ) * NNode + I, ( 2 * NDim     ) * NNode + J ) = KE ( ( 2 * NDim     ) * NNode + I, ( 2 * NDim     ) * NNode + J ) +     PSI_PSI_T  ( I, J ) * FAC_C_3   ;

                                              KE (                            I,       NDim       * NNode + J ) = KE (                            I,       NDim       * NNode + J ) +     PHIX_PSI_T ( I, J ) * BETA_Y ;  ! 1
                                              KE (                            I,   2 * NDim       * NNode + J ) = KE (                            I,   2 * NDim       * NNode + J ) +     PHIY_PSI_T ( I, J ) * BETA_X ;  ! 2
                                              KE (     ( NDim - 1 ) * NNode + I, ( 2 * NDim - 1 ) * NNode + J ) = KE (     ( NDim - 1 ) * NNode + I, ( 2 * NDim - 1 ) * NNode + J ) +     PHIY_PSI_T ( I, J ) * BETA_X ;  ! 3
                                              KE (     ( NDim - 1 ) * NNode + I, ( 2 * NDim     ) * NNode + J ) = KE (     ( NDim - 1 ) * NNode + I, ( 2 * NDim     ) * NNode + J ) +     PHIX_PSI_T ( I, J ) * BETA_Y ;  ! 4

                                              KE (       NDim       * NNode + J,                            I ) = KE (       NDim       * NNode + J,                            I ) +     PHIX_PSI_T ( I, J ) * BETA_Y ;  ! 1
                                              KE (   2 * NDim       * NNode + J,                            I ) = KE (   2 * NDim       * NNode + J,                            I ) +     PHIY_PSI_T ( I, J ) * BETA_X ;  ! 2
                                              KE ( ( 2 * NDim - 1 ) * NNode + J,     ( NDim - 1 ) * NNode + I ) = KE ( ( 2 * NDim - 1 ) * NNode + J,     ( NDim - 1 ) * NNode + I ) +     PHIY_PSI_T ( I, J ) * BETA_X ;  ! 3
                                              KE ( ( 2 * NDim     ) * NNode + J,     ( NDim - 1 ) * NNode + I ) = KE ( ( 2 * NDim     ) * NNode + J,     ( NDim - 1 ) * NNode + I ) +     PHIX_PSI_T ( I, J ) * BETA_Y ;  ! 4

            ! DAMPING MATRIX
            ForAll ( K = 1:NDim             ) CE ( ( K - 1        ) * NNode + I, ( K - 1        ) * NNode + J ) = CE ( ( K - 1        ) * NNode + I, ( K - 1        ) * NNode + J ) +     PHI_PHI_T  ( I, J ) * FAC_B_Rho ;
            ForAll ( K = NDim + 1: 2 * NDim ) CE ( ( K - 1        ) * NNode + I, ( K - 1        ) * NNode + J ) = CE ( ( K - 1        ) * NNode + I, ( K - 1        ) * NNode + J ) +     PSI_PSI_T  ( I, J ) * FAC_B_1   ;

                                              CE (       NDim       * NNode + I, (     NDim + 1 ) * NNode + J ) = CE (       NDim       * NNode + I, (     NDim + 1 ) * NNode + J ) +     PSI_PSI_T  ( I, J ) * FAC_B_2   ;
                                              CE ( (     NDim + 1 ) * NNode + I,       NDim       * NNode + J ) = CE ( (     NDim + 1 ) * NNode + I,       NDim       * NNode + J ) +     PSI_PSI_T  ( I, J ) * FAC_B_2   ;
                                              CE ( ( 2 * NDim     ) * NNode + I, ( 2 * NDim     ) * NNode + J ) = CE ( ( 2 * NDim     ) * NNode + I, ( 2 * NDim     ) * NNode + J ) +     PSI_PSI_T  ( I, J ) * FAC_B_3   ;

                                              CE (                            I,       NDim       * NNode + J ) = CE (                            I,       NDim       * NNode + J ) +     PHIX_PSI_T ( I, J ) * alpha_Y ;  ! 1
                                              CE (                            I,   2 * NDim       * NNode + J ) = CE (                            I,   2 * NDim       * NNode + J ) +     PHIY_PSI_T ( I, J ) * alpha_X ;  ! 2
                                              CE (     ( NDim - 1 ) * NNode + I, ( 2 * NDim - 1 ) * NNode + J ) = CE (     ( NDim - 1 ) * NNode + I, ( 2 * NDim - 1 ) * NNode + J ) +     PHIY_PSI_T ( I, J ) * alpha_X ;  ! 3
                                              CE (     ( NDim - 1 ) * NNode + I, ( 2 * NDim     ) * NNode + J ) = CE (     ( NDim - 1 ) * NNode + I, ( 2 * NDim     ) * NNode + J ) +     PHIX_PSI_T ( I, J ) * alpha_Y ;  ! 4

                                              CE (       NDim       * NNode + J,                            I ) = CE (       NDim       * NNode + J,                            I ) +     PHIX_PSI_T ( I, J ) * alpha_Y ;  ! 1
                                              CE (   2 * NDim       * NNode + J,                            I ) = CE (   2 * NDim       * NNode + J,                            I ) +     PHIY_PSI_T ( I, J ) * alpha_X ;  ! 2
                                              CE ( ( 2 * NDim - 1 ) * NNode + J,     ( NDim - 1 ) * NNode + I ) = CE ( ( 2 * NDim - 1 ) * NNode + J,     ( NDim - 1 ) * NNode + I ) +     PHIY_PSI_T ( I, J ) * alpha_X ;  ! 3
                                              CE ( ( 2 * NDim     ) * NNode + J,     ( NDim - 1 ) * NNode + I ) = CE ( ( 2 * NDim     ) * NNode + J,     ( NDim - 1 ) * NNode + I ) +     PHIX_PSI_T ( I, J ) * alpha_Y ;  ! 4
          End ForAll ;

      End Do ;
  End Do ;


!#Write(*     ,*) 'End Subroutine < MassDampStiffPML_2D_8N_Pattern >' ;
!#Write(UnInf,*) 'End Subroutine < MassDampStiffPML_2D_8N_Pattern >' ;
Return ;
End Subroutine MassDampStiffPML_2D_8N_Pattern ;

End Module PML_Matrices_Pattern ;
