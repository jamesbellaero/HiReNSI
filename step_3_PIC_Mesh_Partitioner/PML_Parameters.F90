
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                                                                                                                  ++
! Last Update:  16 May 2012                                                                                                                        ++
!                                                                                                                                                  ++
! Description: THIS Module CALCULATES PML COEFFICIENTS AT EACH GAUSSIAN POINT.                                                                     ++
!                                                                                                                                                  ++
!              Subroutines                     Description                                                                                         ++
!              ==============                  ==============                                                                                      ++
!              PML_Parameters                                                                                                                      ++
!                                                                                                                                                  ++
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Module PML_PRMTRS ;

Use Parameters ;

Implicit None ;

  Interface
!    Module Procedure PML_Parameters_2D, PML_Parameters_3D ;
  End Interface    ;

Contains ;

! ***************************************************************************************************************************************************
!                                                                                                                                                  **
! Last Update:  12 JUNE 2011                                                                                                                       **
! Description: THIS Subroutine CALCULATES PML Parameters for a 2D model.                                                                           **
!                                                                                                                                                  **
! ***************************************************************************************************************************************************

Subroutine PML_Parameters_2D  (                                                                                                                    &
NDim, PARAM_Type, IEL,                                                                                                                             & ! Integer Variables
Lambda, MU, Rho, R_PML_Alpha_0, C_REF_Beta_0, M_PML, B_PML,                                                                                        & ! Real(IN ) Variables
FAC_A_Rho, FAC_B_Rho, FAC_C_Rho, FAC_A_1, FAC_A_2, FAC_A_3, FAC_B_1, FAC_B_2, FAC_B_3, FAC_C_1, FAC_C_2, FAC_C_3, Alpha_X, Alpha_Y, Beta_X, Beta_Y, A, B, C,  & ! Real(OUT) Variables
!                                                                                                                                                  & ! Integer Arrays
XG, PML_DIM                                                                                                                                        & ! Real Arrays
!                                                                                                                                                  & ! Characters
) ;

Implicit None ;

! =========================== Global Variables ======================================================================================================

! - Integer Arrays ------------------------------------------------------------------------------------------------------------------------------
Integer (Kind=Tiny), Intent(In)    :: NDim, PARAM_Type ;
Integer (Kind=Lng ), Intent(In)    :: IEL ;

! - Real Variables ----------------------------------------------------------------------------------------------------------------------------------
Real (Kind=Dbl), Intent(In)    :: Lambda, MU, Rho, R_PML_Alpha_0, C_REF_Beta_0, M_PML, B_PML ;
Real (Kind=Dbl), Intent(OUT)   :: FAC_A_Rho, FAC_B_Rho, FAC_C_Rho, FAC_A_1, FAC_A_2, FAC_A_3, FAC_B_1, FAC_B_2, FAC_B_3, FAC_C_1, FAC_C_2, FAC_C_3, Alpha_X, Alpha_Y, Beta_X, Beta_Y, A, B, C ;  ! PML Parameters

! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
! - Integer Arrays ----------------------------------------------------------------------------------------------------------------------------------

! - Real Arrays -------------------------------------------------------------------------------------------------------------------------------------
Real (Kind=Dbl), Intent(In), Dimension (:  )  :: XG ;
Real (Kind=Dbl), Intent(In), Dimension (:,:)  :: PML_DIM ;

! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
! - Character Variables -----------------------------------------------------------------------------------------------------------------------------
! - Logical Variables -------------------------------------------------------------------------------------------------------------------------------

! =========================== LOCAL Variables =======================================================================================================

! - Integer Arrays ------------------------------------------------------------------------------------------------------------------------------
Integer (Kind=Shrt)  :: I ;

! - Real Variables ----------------------------------------------------------------------------------------------------------------------------------
Real (Kind=Dbl)  :: S0, L_PML, N_S, Alpha_0, Beta_0 ;

! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
! - Integer Arrays ----------------------------------------------------------------------------------------------------------------------------------

! - Real Arrays -------------------------------------------------------------------------------------------------------------------------------------
Real (Kind=Dbl), Dimension ( NDim )  :: Alpha, Beta ;
! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
! - Character Variables -----------------------------------------------------------------------------------------------------------------------------
! - Logical Variables -------------------------------------------------------------------------------------------------------------------------------
! - Type DECLERATIONS -------------------------------------------------------------------------------------------------------------------------------

! =========================== Subroutine CODE =======================================================================================================

Alpha = 1.0_Dbl ;
Beta  = 0.0_Dbl ;

  ! FIND THE POSITION OF GAUSSIAN POINT WITH RESPECT TO THE REGULAR DOMAIN AND PML ZONE
  DO I = 1, NDim

      !IF     ( XG (I) <= PML_DIM ( (I - 1 ) * 2 + 1, 1) ) Then ;  ! GAUSSIAN POINT IS LOCATED BEFORE THE REGULAR DOMAIN ( LEFT OR BELOW )
      IF      ( XG (I) <  PML_DIM ( (I - 1 ) * 2 + 1, 1) ) Then ;  ! GAUSSIAN POINT IS LOCATED BEFORE THE REGULAR DOMAIN ( LEFT OR BELOW )
        S0    = PML_DIM ( (I - 1 ) * 2 + 1, 1) ;
        L_PML = PML_DIM ( (I - 1 ) * 2 + 1, 2) ;
        N_S   = - 1.0_Dbl                      ;
      !Else If ( XG (I) >= PML_DIM ( (I - 1 ) * 2 + 2, 1) ) Then ;  ! GAUSSIAN POINT IS LOCATED AFTER THE REGULAR DOMAIN ( RIGHT OR ABOVE )
      Else If ( XG (I) >  PML_DIM ( (I - 1 ) * 2 + 2, 1) ) Then ;  ! GAUSSIAN POINT IS LOCATED AFTER THE REGULAR DOMAIN ( RIGHT OR ABOVE )
        S0    = PML_DIM ( (I - 1 ) * 2 + 2, 1) ;
        L_PML = PML_DIM ( (I - 1 ) * 2 + 2, 2) ;
        N_S   = + 1.0_Dbl                      ;
      Else                                                         ! GAUSSIAN POINT IS LOCATED IN THE REGULAR DOMAIN, N0 NEED TO UPDATE Alpha AND Beta
        CYCLE ;
      End If ;

      IF      ( PARAM_Type == 0 ) Then ;
        Alpha_0 = ( M_PML + 1 ) * B_PML        * DLOG10 ( 1.0_Dbl / DABS (R_PML_Alpha_0) ) / ( 2 * L_PML ) ;
        Beta_0  = ( M_PML + 1 ) * C_REF_Beta_0 * DLOG10 ( 1.0_Dbl / DABS (R_PML_Alpha_0) ) / ( 2 * L_PML ) ;
      Else If ( PARAM_Type == 1 ) Then ;
        Alpha_0 = R_PML_Alpha_0 ;
        Beta_0  = C_REF_Beta_0 ;
      End If ;

    Alpha ( I ) = Alpha ( I ) + Alpha_0 * ( ( XG ( I ) - S0 ) * N_S / L_PML     ) ** M_PML ;
    Beta  ( I ) = Beta  ( I ) + Beta_0  * ( ( XG ( I ) - S0 ) * N_S / L_PML     ) ** M_PML ;

  End Do ;
 
! IF ( Alpha (1) == 1.0_Dbl .AND. Alpha (2) == 1.0_Dbl .AND. Beta (1) == 0.0_Dbl .AND. Beta (2) == 0.0_Dbl ) Then ;
!  Write (*     ,"(' THE GAUSSIAN POINT OF ELEMENET ', I12, '  IS LOCATED WITHIN THE REGULAR DOMAIN WHILE IT IS MARKED AS A PML DOMAIN.')") IEL ;
!  Write (UnInf,"(' THE GAUSSIAN POINT OF ELEMENET ', I12, '  IS LOCATED WITHIN THE REGULAR DOMAIN WHILE IT IS MARKED AS A PML DOMAIN.')") IEL ;
!  !#Call BEEP_FAIL ;
!  Write(*, Fmt_FL) ;  Write(UnInf, Fmt_FL) ; Write(*, Fmt_End) ; Read(*,*) ;  STOP ;
! End If ;
 
! PML REQUIRED COEFFICIENTS
Alpha_X = Alpha ( 1 ) ;
Alpha_Y = Alpha ( 2 ) ;
Beta_X  = Beta ( 1 ) ;
Beta_Y  = Beta ( 2 ) ;

A = Alpha_X * Alpha_Y ;
B = Alpha_X * Beta_Y  + Alpha_Y * Beta_X ;
C = Beta_X  * Beta_Y ;

FAC_A_Rho = A * Rho ;
FAC_B_Rho = B * Rho ;
FAC_C_Rho = C * Rho ;

FAC_A_1 = - A * ( Lambda + 2 * MU ) / ( 4 * MU * ( Lambda + MU ) ) ;
FAC_A_2 = + A *   Lambda            / ( 4 * MU * ( Lambda + MU ) ) ;
FAC_A_3 = - A                       / (     MU                 ) ;

FAC_B_1 = - B * ( Lambda + 2 * MU ) / ( 4 * MU * ( Lambda + MU ) ) ;
FAC_B_2 = + B *   Lambda            / ( 4 * MU * ( Lambda + MU ) ) ;
FAC_B_3 = - B                       / (     MU                 ) ;

FAC_C_1 = - C * ( Lambda + 2 * MU ) / ( 4 * MU * ( Lambda + MU ) ) ;
FAC_C_2 = + C *   Lambda            / ( 4 * MU * ( Lambda + MU ) ) ;
FAC_C_3 = - C                       / (     MU                 ) ;



!#Write(*    ,*) 'End Subroutine < PML_Parameters_2D >' ;
!#Write(UnInf,*) 'End Subroutine < PML_Parameters_2D >' ;
Return ;
End Subroutine PML_Parameters_2D ;

! ***************************************************************************************************************************************************
!                                                                                                                                                  **
! Last Update:  16 May 2012                                                                                                                        **
! Description: THIS Subroutine CALCULATES PML Parameters for a 3D model.                                                                           **
!                                                                                                                                                  **
! ***************************************************************************************************************************************************

Subroutine PML_Parameters_3D  (                                                                                                                    &
NDim, PARAM_Type, IEL,                                                                                                                             & ! Integer Variables
Lambda, MU, Rho, R_PML_Alpha_0, C_REF_Beta_0, M_PML, B_PML,                                                                                        & ! Real(IN ) Variables
FAC_A_Rho, FAC_B_Rho, FAC_C_Rho, FAC_D_Rho, FAC_A_1, FAC_A_2, FAC_A_3, FAC_B_1, FAC_B_2, FAC_B_3, FAC_C_1, FAC_C_2, FAC_C_3, FAC_D_1, FAC_D_2, FAC_D_3, Alpha_X, Alpha_Y, Alpha_Z, Beta_X, Beta_Y, Beta_Z,  & ! Real(OUT) Variables
!                                                                                                                                                  & ! Integer Arrays
XG, PML_DIM                                                                                                                                        & ! Real Arrays
!                                                                                                                                                  & ! Characters
) ;

Implicit None ;

! =========================== Global Variables ======================================================================================================

! - Integer Arrays ------------------------------------------------------------------------------------------------------------------------------
Integer (Kind=Tiny), Intent(In)    :: NDim, PARAM_Type ;
Integer (Kind=Lng ), Intent(In)    :: IEL ;
!#Integer (Kind=Shrt), Intent(InOut) ::  ;
!#Integer (Kind=Shrt), Intent(OUT)   ::  ;

! - Real Variables ----------------------------------------------------------------------------------------------------------------------------------
Real (Kind=Dbl), Intent(In)    :: Lambda, MU, Rho, R_PML_Alpha_0, C_REF_Beta_0, M_PML, B_PML ;
Real (Kind=Dbl), Intent(OUT)   :: FAC_A_Rho, FAC_B_Rho, FAC_C_Rho, FAC_D_Rho, FAC_A_1, FAC_A_2, FAC_A_3, FAC_B_1, FAC_B_2, FAC_B_3, FAC_C_1, FAC_C_2, FAC_C_3, FAC_D_1, FAC_D_2, FAC_D_3, Alpha_X, Alpha_Y, Alpha_Z, Beta_X, Beta_Y, Beta_Z ;  ! PML Parameters

! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
!#Complex, Intent(In)    ::  ;
!#Complex, Intent(InOut) ::  ;
!#Complex, Intent(OUT)   ::  ;
! - Integer Arrays ----------------------------------------------------------------------------------------------------------------------------------
!#Integer (Kind=Shrt), Intent(In), Dimension (:  )  ::
!#Integer (Kind=Shrt), Intent(In), Dimension (:,:)  ::
!#Integer (Kind=Shrt), Intent(In)    ::  ;
!#Integer (Kind=Shrt), Intent(InOut) ::  ;
!#Integer (Kind=Shrt), Intent(OUT)   ::  ;

! - Real Arrays -------------------------------------------------------------------------------------------------------------------------------------
Real (Kind=Dbl), Intent(In), Dimension (:  )  :: XG ;
Real (Kind=Dbl), Intent(In), Dimension (:,:)  :: PML_DIM ;
!#Real (Kind=Dbl), Intent(InOut) ::  ;
!#Real (Kind=Dbl), Intent(OUT)   ::  ;

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
Integer (Kind=Shrt)  :: I ;

! - Real Variables ----------------------------------------------------------------------------------------------------------------------------------
Real (Kind=Dbl)  :: S0, L_PML, N_S, Alpha_0, Beta_0, A, B, C, d ;

! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
!#Complex  ::  ;
! - Integer Arrays ----------------------------------------------------------------------------------------------------------------------------------
!#Integer (Kind=Shrt)  ::  ;

! - Real Arrays -------------------------------------------------------------------------------------------------------------------------------------
Real (Kind=Dbl), Dimension ( NDim )  :: Alpha, Beta ;
! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
!#Complex  ::  ;
! - Character Variables -----------------------------------------------------------------------------------------------------------------------------
!#Character   ::  ;
! - Logical Variables -------------------------------------------------------------------------------------------------------------------------------
!#Logical   ::  ;
! - Type DECLERATIONS -------------------------------------------------------------------------------------------------------------------------------
!#Type() :: ;

! =========================== Subroutine CODE =======================================================================================================

Alpha = 1.0_Dbl ;
Beta  = 0.0_Dbl ;

  ! FIND THE POSITION OF GAUSSIAN POINT WITH RESPECT TO THE REGULAR DOMAIN AND PML ZONE
  DO I = 1, NDim

      IF      ( XG (I) <= PML_DIM ( (I - 1 ) * 2 + 1, 1) ) Then ;  ! GAUSSIAN POINT IS LOCATED BEFORE THE REGULAR DOMAIN ( LEFT OR BELOW )
        S0    = PML_DIM ( (I - 1 ) * 2 + 1, 1) ;
        L_PML = PML_DIM ( (I - 1 ) * 2 + 1, 2) ;
        N_S   = - 1.0_Dbl                      ;
      Else If ( XG (I) >= PML_DIM ( (I - 1 ) * 2 + 2, 1) ) Then ;  ! GAUSSIAN POINT IS LOCATED AFTER THE REGULAR DOMAIN ( RIGHT OR ABOVE )
        S0    = PML_DIM ( (I - 1 ) * 2 + 2, 1) ;
        L_PML = PML_DIM ( (I - 1 ) * 2 + 2, 2) ;
        N_S   = + 1.0_Dbl                      ;
      Else                                                         ! GAUSSIAN POINT IS LOCATED IN THE REGULAR DOMAIN, N0 NEED TO UPDATE Alpha AND Beta
        CYCLE ;
      End If ;

      IF      ( PARAM_Type == 0 ) Then ;
        Alpha_0 = ( M_PML + 1 ) * B_PML        * DLOG10 ( 1.0_Dbl / DABS (R_PML_Alpha_0) ) / ( 2 * L_PML ) ;
        Beta_0  = ( M_PML + 1 ) * C_REF_Beta_0 * DLOG10 ( 1.0_Dbl / DABS (R_PML_Alpha_0) ) / ( 2 * L_PML ) ;
      Else If ( PARAM_Type == 1 ) Then ;
        Alpha_0 = R_PML_Alpha_0 ;
        Beta_0  = C_REF_Beta_0 ;
      End If ;

    Alpha ( I ) = Alpha ( I ) + Alpha_0 * (    ( XG ( I ) - S0 ) * N_S / L_PML     ) ** M_PML ;
    Beta  ( I ) = Beta  ( I ) + Beta_0  * (    ( XG ( I ) - S0 ) * N_S / L_PML     ) ** M_PML ;

  End Do ;

 IF ( Alpha (1) == 1.0_Dbl .AND. Alpha (2) == 1.0_Dbl .AND. Beta (1) == 0.0_Dbl .AND. Beta (2) == 0.0_Dbl ) Then ;
  Write (*     ,"(' THE GAUSSIAN POINT OF ELEMENET ', I12, '  IS LOCATED WITHIN THE REGULAR DOMAIN WHILE IT IS MARKED AS A PML DOMAIN.')") IEL ;
  Write (UnInf,"(' THE GAUSSIAN POINT OF ELEMENET ', I12, '  IS LOCATED WITHIN THE REGULAR DOMAIN WHILE IT IS MARKED AS A PML DOMAIN.')") IEL ;
  !#Call BEEP_FAIL ;
  Write(*, Fmt_FL) ;  Write(UnInf, Fmt_FL) ; Write(*, Fmt_End) ; Read(*,*) ;  STOP ;
 End If ;
 
! PML REQUIRED COEFFICIENTS
Alpha_X = Alpha ( 1 ) ;
Alpha_Y = Alpha ( 2 ) ;
Alpha_Z = Alpha ( 3 ) ;

Beta_X = Beta ( 1 ) ;
Beta_Y = Beta ( 2 ) ;
Beta_Z = Beta ( 3 ) ;

a = Alpha_X * Alpha_Y * Alpha_Z ;
b = Alpha_X * Alpha_Y * Beta_Z + Alpha_X * Alpha_Z * Beta_Y + Alpha_Y * Alpha_Z * Beta_X ;
c = Alpha_X * Beta_Y * Beta_Z + Alpha_Y * Beta_X * Beta_Z + Alpha_Z * Beta_X * Beta_Y ;
d = Beta_X * Beta_Y * Beta_Z ;

FAC_A_Rho = a * Rho ;
FAC_B_Rho = b * Rho ;
FAC_C_Rho = c * Rho ;
FAC_D_Rho = d * Rho ;

FAC_A_1 = - a * ( Lambda + MU ) / ( Mu *     ( 3 * Lambda + 2 * Mu ) ) ;
FAC_A_2 = + a *   Lambda        / ( 2 * Mu * ( 3 * Lambda + 2 * Mu ) ) ;
FAC_A_3 = - a                   / (     MU                           ) ;

FAC_b_1 = - b * ( Lambda + MU ) / ( Mu *     ( 3 * Lambda + 2 * Mu ) ) ;
FAC_b_2 = + b *   Lambda        / ( 2 * Mu * ( 3 * Lambda + 2 * Mu ) ) ;
FAC_b_3 = - b                   / (     MU                           ) ;

FAC_c_1 = - c * ( Lambda + MU ) / ( Mu *     ( 3 * Lambda + 2 * Mu ) ) ;
FAC_c_2 = + c *   Lambda        / ( 2 * Mu * ( 3 * Lambda + 2 * Mu ) ) ;
FAC_c_3 = - c                   / (     MU                           ) ;

FAC_d_1 = - d * ( Lambda + MU ) / ( Mu *     ( 3 * Lambda + 2 * Mu ) ) ;
FAC_d_2 = + d *   Lambda        / ( 2 * Mu * ( 3 * Lambda + 2 * Mu ) ) ;
FAC_d_3 = - d                   / (     MU                           ) ;

!#Write(*     ,*) 'End Subroutine < PML_Parameters_3D >' ;
!#Write(UnInf,*) 'End Subroutine < PML_Parameters_3D >' ;
Return ;
End Subroutine PML_Parameters_3D ;

! ***************************************************************************************************************************************************
!                                                                                                                                                  **
! Start:        17 August 2012                                                                                                                     **
! Last Update:  17 August 2012                                                                                                                     **
! Description: THIS Subroutine CALCULATES MPML Parameters for a 2D model.                                                                          **
!                                                                                                                                                  **
! ***************************************************************************************************************************************************

Subroutine MPML_Parameters_2D  (                                                                                                &
NDim, PARAM_Type,                                                                                                               & ! Integer (1) Variables
!                                                                                                                               & ! Integer (2) Variables
!                                                                                                                               & ! Integer (4) Variables
IEl,                                                                                                                            & ! Integer (8) Variables
Rho, R_PML_Alpha_0, C_REF_Beta_0, M_PML, B_PML, Ratio, FAC_A_Rho, FAC_B_Rho, FAC_C_Rho, Alpha_X, Alpha_Y, Beta_X, Beta_Y, DBetaX_Y, DBetaY_X, DAlphaX_Y, DAlphaY_X, A, B, C, & ! Real Variables
!                                                                                                                               & ! Integer Arrays
XG, PML_DIM                                                                                                                     & ! Real Arrays
!                                                                                                                               & ! Characters
!                                                                                                                               & ! Type
) ;

Implicit None ;

! =========================== Global Variables ======================================================================================================

! - Integer Arrays ------------------------------------------------------------------------------------------------------------------------------
Integer (Kind=Tiny), Intent(In)    :: NDim, PARAM_Type ;
Integer (Kind=Lng ), Intent(In)    :: IEL ;

! - Real Variables ----------------------------------------------------------------------------------------------------------------------------------
Real (Kind=Dbl), Intent(In)    :: Rho, R_PML_Alpha_0, C_REF_Beta_0, M_PML, B_PML, Ratio ;
Real (Kind=Dbl), Intent(OUT)   :: FAC_A_Rho, FAC_B_Rho, FAC_C_Rho, Alpha_X, Alpha_Y, Beta_X, Beta_Y, DBetaX_Y, DBetaY_X, DAlphaX_Y, DAlphaY_X, A, B, C ;  ! PML Parameters

! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
! - Integer Arrays ----------------------------------------------------------------------------------------------------------------------------------
! - Real Arrays -------------------------------------------------------------------------------------------------------------------------------------
Real (Kind=Dbl), Intent(In), Dimension (:  )  :: XG ;
Real (Kind=Dbl), Intent(In), Dimension (:,:)  :: PML_DIM ;

! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
! - Character Variables -----------------------------------------------------------------------------------------------------------------------------
! - Logical Variables -------------------------------------------------------------------------------------------------------------------------------

! =========================== LOCAL Variables =======================================================================================================

! - Integer Arrays ------------------------------------------------------------------------------------------------------------------------------
Integer (Kind=Shrt)  :: I,J ;

! - Real Variables ----------------------------------------------------------------------------------------------------------------------------------
Real (Kind=Dbl)  :: S0, L_PML, N_S, Alpha_0, Beta_0, Temp ;

! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
! - Integer Arrays ----------------------------------------------------------------------------------------------------------------------------------

! - Real Arrays -------------------------------------------------------------------------------------------------------------------------------------
Real (Kind=Dbl), Dimension ( NDim )        :: Alpha, Beta, Alpha_Temp, Beta_Temp ;
Real (Kind=Dbl), Dimension ( NDim, NDim )  :: DifAlpha, DifBeta ;
! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
! - Character Variables -----------------------------------------------------------------------------------------------------------------------------
! - Logical Variables -------------------------------------------------------------------------------------------------------------------------------
! - Type DECLERATIONS -------------------------------------------------------------------------------------------------------------------------------

! =========================== Subroutine CODE =======================================================================================================

! Initialize stretch function
Alpha (:) = 1.0_Dbl ;
Beta  (:) = 0.0_Dbl ; ! Attenuation function

! Initialize the differential of stretch function
DifAlpha (:,:) = 0.0_Dbl ;
DifBeta  (:,:) = 0.0_Dbl ;

  ! FIND THE POSITION OF GAUSSIAN POINT WITH RESPECT TO THE REGULAR DOMAIN AND PML ZONE
  DO I = 1, NDim

      !IF     ( XG (I) <= PML_DIM ( (I - 1 ) * 2 + 1, 1) ) Then ;  ! GAUSSIAN POINT IS LOCATED BEFORE THE REGULAR DOMAIN ( LEFT OR BELOW )
      IF      ( XG (I) <  PML_DIM ( (I - 1 ) * 2 + 1, 1) ) Then ;  ! GAUSSIAN POINT IS LOCATED BEFORE THE REGULAR DOMAIN ( LEFT OR BELOW )
        S0    = PML_DIM ( (I - 1 ) * 2 + 1, 1) ;
        L_PML = PML_DIM ( (I - 1 ) * 2 + 1, 2) ;
        N_S   = - 1.0_Dbl                      ;
      !Else If ( XG (I) >= PML_DIM ( (I - 1 ) * 2 + 2, 1) ) Then ;  ! GAUSSIAN POINT IS LOCATED AFTER THE REGULAR DOMAIN ( RIGHT OR ABOVE )
      Else If ( XG (I) >  PML_DIM ( (I - 1 ) * 2 + 2, 1) ) Then ;  ! GAUSSIAN POINT IS LOCATED AFTER THE REGULAR DOMAIN ( RIGHT OR ABOVE )
        S0    = PML_DIM ( (I - 1 ) * 2 + 2, 1) ;
        L_PML = PML_DIM ( (I - 1 ) * 2 + 2, 2) ;
        N_S   = + 1.0_Dbl                      ;
      Else                                                         ! GAUSSIAN POINT IS LOCATED IN THE REGULAR DOMAIN, N0 NEED TO UPDATE Alpha AND Beta
        CYCLE ;
      End If ;

      IF      ( PARAM_Type == 0 ) Then ;
        Alpha_0 = ( M_PML + 1 ) * B_PML        * DLOG10 ( 1.0_Dbl / DABS (R_PML_Alpha_0) ) / ( 2 * L_PML ) ;
        Beta_0  = ( M_PML + 1 ) * C_REF_Beta_0 * DLOG10 ( 1.0_Dbl / DABS (R_PML_Alpha_0) ) / ( 2 * L_PML ) ;
      Else If ( PARAM_Type == 1 ) Then ;
        Alpha_0 = R_PML_Alpha_0 ;
        Beta_0  = C_REF_Beta_0 ;
      End If ;

    ! Calculating stretch function
    Alpha ( I ) = Alpha ( I ) + Alpha_0 * ( ( XG ( I ) - S0 ) * N_S / L_PML     ) ** M_PML ;
    Beta  ( I ) = Beta  ( I ) + Beta_0  * ( ( XG ( I ) - S0 ) * N_S / L_PML     ) ** M_PML ;

      Do J = 1, NDim ;
        If ( I == J ) Cycle ;
        ! Calculating the derivatives of stretch function
        ! alpha
        DifAlpha ( J, I ) = Ratio * Alpha_0 * M_PML * ( N_S / L_PML ) *   ( ( ( XG ( I ) - S0 ) * N_S / L_PML     ) ** ( M_PML - 1 ) ) ;

        ! Beta
        DifBeta  ( J, I ) = Ratio * Beta_0  * M_PML * ( N_S / L_PML ) *   ( ( ( XG ( I ) - S0 ) * N_S / L_PML     ) ** ( M_PML - 1 ) ) ;
      End Do ;

  End Do ;

! Adding portion of damping in each direction to other direction
Alpha_Temp (:) = Alpha (:) ;
Beta_Temp (:)  = Beta (:) ;

!Where (Alpha_Temp == 1 ) Alpha_Temp = 0 ;
Alpha_Temp (:) = Alpha_Temp (:) - 1.0_Dbl ;

  Do I = 1, NDim ;


    Temp = Alpha_Temp ( I ) * Ratio ;
    Alpha (I) = Alpha (I) - Temp ;
    Alpha (:) = Alpha (:) + Temp ;


    Temp = Beta_Temp ( I ) * Ratio ;
    Beta (I) = Beta (I) - Temp ;
    Beta (:) = Beta (:) + Temp ;

  End Do ;

! IF ( Alpha (1) == 1.0_Dbl .AND. Alpha (2) == 1.0_Dbl .AND. Beta (1) == 0.0_Dbl .AND. Beta (2) == 0.0_Dbl ) Then ;
!  Write (*     ,"(' THE GAUSSIAN POINT OF ELEMENET ', I12, '  IS LOCATED WITHIN THE REGULAR DOMAIN WHILE IT IS MARKED AS A PML DOMAIN.')") IEL ;
!  Write (UnInf,"(' THE GAUSSIAN POINT OF ELEMENET ', I12, '  IS LOCATED WITHIN THE REGULAR DOMAIN WHILE IT IS MARKED AS A PML DOMAIN.')") IEL ;
!  !#Call BEEP_FAIL ;
!  Write(*, Fmt_FL) ;  Write(UnInf, Fmt_FL) ; Write(*, Fmt_End) ; Read(*,*) ;  STOP ;
! End If ;


! PML REQUIRED COEFFICIENTS
Alpha_X = Alpha ( 1 ) ;
Alpha_Y = Alpha ( 2 ) ;
Beta_X  = Beta  ( 1 ) ;
Beta_Y  = Beta  ( 2 ) ;

DBetaX_Y  = DifBeta  ( 1, 2 ) ;
DBetaY_X  = DifBeta  ( 2, 1 ) ;

DAlphaX_Y = DifAlpha ( 1, 2 ) ;
DAlphaY_X = DifAlpha ( 2, 1 ) ;

A = Alpha_X * Alpha_Y ;
B = Alpha_X * Beta_Y  + Alpha_Y * Beta_X ;
C = Beta_X  * Beta_Y ;

FAC_A_Rho = A * Rho ;
FAC_B_Rho = B * Rho ;
FAC_C_Rho = C * Rho ;

!#Write(*     ,*) 'End Subroutine < MPML_Parameters_2D >' ;
!#Write(UnInf,*) 'End Subroutine < MPML_Parameters_2D >' ;
Return ;
End Subroutine MPML_Parameters_2D ;

!###############################################

! ***************************************************************************************************************************************************
!                                                                                                                                                  **
! Start:        07 Oct 2012                                                                                                                        **
! Last Update:  07 Oct 2012                                                                                                                        **
! Description: THIS Subroutine CALCULATES MPML Parameters for a 2D model.                                                                          **
!                                                                                                                                                  **
! ***************************************************************************************************************************************************

Subroutine MPML_Parameters_3D  (                                                                                                &
NDim, PARAM_Type,                                                                                                               & ! Integer (1) Variables
!                                                                                                                               & ! Integer (2) Variables
!                                                                                                                               & ! Integer (4) Variables
IEl,                                                                                                                            & ! Integer (8) Variables
Rho, R_PML_Alpha_0, C_REF_Beta_0, M_PML, B_PML, Ratio, FAC_A_Rho, FAC_B_Rho, FAC_C_Rho, FAC_D_Rho, A, B, C, D, Coef_PML, D_Coef_PML,       & ! Real Variables
!                                                                                                                               & ! Integer Arrays
XG, PML_DIM                                                                                                                     & ! Real Arrays
!                                                                                                                               & ! Characters
!                                                                                                                               & ! Type
) ;

Implicit None ;

! =========================== Global Variables ======================================================================================================

! - Integer Arrays ------------------------------------------------------------------------------------------------------------------------------
Integer (Kind=Tiny), Intent(In)    :: NDim, PARAM_Type ;
Integer (Kind=Lng ), Intent(In)    :: IEL ;

! - Real Variables ----------------------------------------------------------------------------------------------------------------------------------
Real (Kind=Dbl), Intent(In)    :: Rho, R_PML_Alpha_0, C_REF_Beta_0, M_PML, B_PML, Ratio ;
Real (Kind=Dbl), Intent(OUT)   :: FAC_A_Rho, FAC_B_Rho, FAC_C_Rho, FAC_D_Rho, A, B, C, D ;

! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
! - Integer Arrays ----------------------------------------------------------------------------------------------------------------------------------
! - Real Arrays -------------------------------------------------------------------------------------------------------------------------------------
Real (Kind=Dbl), Intent(In),  Dimension (:  )  :: XG ;
Real (Kind=Dbl), Intent(In),  Dimension (:,:)  :: PML_DIM ;

Real (Kind=Dbl), Intent(Out), Dimension (:,:)  :: Coef_PML, D_Coef_PML ;

! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
! - Character Variables -----------------------------------------------------------------------------------------------------------------------------
! - Logical Variables -------------------------------------------------------------------------------------------------------------------------------

! =========================== LOCAL Variables =======================================================================================================

! - Integer Arrays ------------------------------------------------------------------------------------------------------------------------------
Integer (Kind=Shrt)  :: I,J ;

! - Real Variables ----------------------------------------------------------------------------------------------------------------------------------
Real (Kind=Dbl)  :: S0, L_PML, N_S, Alpha_0, Beta_0, Temp ;
Real (Kind=Dbl)  :: Alpha_X, Alpha_Y, Alpha_Z, Beta_X, Beta_Y, Beta_Z, DBetaX_Y, DBetaX_Z, DBetaY_X, DBetaY_Z, DBetaZ_X, DBetaZ_Y, DAlphaX_Y, DAlphaX_Z, DAlphaY_X, DAlphaY_Z, DAlphaZ_X, DAlphaZ_Y  ;  ! PML Parameters

! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
! - Integer Arrays ----------------------------------------------------------------------------------------------------------------------------------

! - Real Arrays -------------------------------------------------------------------------------------------------------------------------------------
Real (Kind=Dbl), Dimension ( NDim )        :: Alpha, Beta, Alpha_Temp, Beta_Temp ;
Real (Kind=Dbl), Dimension ( NDim, NDim )  :: DifAlpha, DifBeta ;
! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
! - Character Variables -----------------------------------------------------------------------------------------------------------------------------
! - Logical Variables -------------------------------------------------------------------------------------------------------------------------------
! - Type DECLERATIONS -------------------------------------------------------------------------------------------------------------------------------

! =========================== Subroutine CODE =======================================================================================================

! Initialize stretch function
Alpha (:) = 1.0_Dbl ;
Beta  (:) = 0.0_Dbl ; ! Attenuation function

! Initialize the differential of stretch function
DifAlpha (:,:) = 0.0_Dbl ;
DifBeta  (:,:) = 0.0_Dbl ;

  ! FIND THE POSITION OF GAUSSIAN POINT WITH RESPECT TO THE REGULAR DOMAIN AND PML ZONE
  DO I = 1, NDim

      IF      ( XG (I) <  PML_DIM ( (I - 1 ) * 2 + 1, 1) ) Then ;  ! GAUSSIAN POINT IS LOCATED BEFORE THE REGULAR DOMAIN ( LEFT OR BELOW )
        S0    = PML_DIM ( (I - 1 ) * 2 + 1, 1) ;
        L_PML = PML_DIM ( (I - 1 ) * 2 + 1, 2) ;
        N_S   = - 1.0_Dbl                      ;
      Else If ( XG (I) >  PML_DIM ( (I - 1 ) * 2 + 2, 1) ) Then ;  ! GAUSSIAN POINT IS LOCATED AFTER THE REGULAR DOMAIN ( RIGHT OR ABOVE )
        S0    = PML_DIM ( (I - 1 ) * 2 + 2, 1) ;
        L_PML = PML_DIM ( (I - 1 ) * 2 + 2, 2) ;
        N_S   = + 1.0_Dbl                      ;
      Else                                                         ! GAUSSIAN POINT IS LOCATED IN THE REGULAR DOMAIN, N0 NEED TO UPDATE Alpha AND Beta
        CYCLE ;
      End If ;

      IF      ( PARAM_Type == 0 ) Then ;
        Alpha_0 = ( M_PML + 1 ) * B_PML        * DLOG10 ( 1.0_Dbl / DABS (R_PML_Alpha_0) ) / ( 2 * L_PML ) ;
        Beta_0  = ( M_PML + 1 ) * C_REF_Beta_0 * DLOG10 ( 1.0_Dbl / DABS (R_PML_Alpha_0) ) / ( 2 * L_PML ) ;
      Else If ( PARAM_Type == 1 ) Then ;
        Alpha_0 = R_PML_Alpha_0 ;
        Beta_0  = C_REF_Beta_0 ;
      End If ;

    ! Calculating stretch function
    Alpha ( I ) = Alpha ( I ) + Alpha_0 * ( ( XG ( I ) - S0 ) * N_S / L_PML     ) ** M_PML ;
    Beta  ( I ) = Beta  ( I ) + Beta_0  * ( ( XG ( I ) - S0 ) * N_S / L_PML     ) ** M_PML ;

      Do J = 1, NDim ;
        If ( I == J ) Cycle ;
        ! Calculating the derivatives of stretch function
        ! alpha
        DifAlpha ( J, I ) = Ratio * Alpha_0 * M_PML * ( N_S / L_PML ) *   ( ( ( XG ( I ) - S0 ) * N_S / L_PML     ) ** ( M_PML - 1 ) ) ;

        ! Beta
        DifBeta  ( J, I ) = Ratio * Beta_0  * M_PML * ( N_S / L_PML ) *   ( ( ( XG ( I ) - S0 ) * N_S / L_PML     ) ** ( M_PML - 1 ) ) ;
      End Do ;

  End Do ;

! Adding portion of damping in each direction to other direction
Alpha_Temp (:) = Alpha (:) ;
Beta_Temp (:)  = Beta (:) ;

Alpha_Temp (:) = Alpha_Temp (:) - 1.0_Dbl ;

  Do I = 1, NDim ;

    Temp = Alpha_Temp ( I ) * Ratio ;
    Alpha (I) = Alpha (I) - Temp ;
    Alpha (:) = Alpha (:) + Temp ;

    Temp = Beta_Temp ( I ) * Ratio ;
    Beta (I) = Beta (I) - Temp ;
    Beta (:) = Beta (:) + Temp ;

  End Do ;

! IF ( Alpha (1) == 1.0_Dbl .AND. Alpha (2) == 1.0_Dbl .AND. Beta (1) == 0.0_Dbl .AND. Beta (2) == 0.0_Dbl ) Then ;
!  Write (*     ,"(' THE GAUSSIAN POINT OF ELEMENET ', I12, '  IS LOCATED WITHIN THE REGULAR DOMAIN WHILE IT IS MARKED AS A PML DOMAIN.')") IEL ;
!  Write (UnInf,"(' THE GAUSSIAN POINT OF ELEMENET ', I12, '  IS LOCATED WITHIN THE REGULAR DOMAIN WHILE IT IS MARKED AS A PML DOMAIN.')") IEL ;
!  !#Call BEEP_FAIL ;
!  Write(*, Fmt_FL) ;  Write(UnInf, Fmt_FL) ; Write(*, Fmt_End) ; Read(*,*) ;  STOP ;
! End If ;

! PML REQUIRED COEFFICIENTS
Alpha_X = Alpha ( 1 ) ;
Alpha_Y = Alpha ( 2 ) ;
Alpha_Z = Alpha ( 3 ) ;

Beta_X  = Beta  ( 1 ) ;
Beta_Y  = Beta  ( 2 ) ;
Beta_Z  = Beta  ( 3 ) ;


!if (iel == 1 ) write(*,"(6(e12.3))")Alpha ( 1 ),Alpha ( 2 ),Alpha ( 3 ),Beta  ( 1 ),Beta  ( 2 ),Beta  ( 3 )

DBetaX_Y  = DifBeta  ( 1, 2 ) ;
DBetaX_Z  = DifBeta  ( 1, 3 ) ;

DBetaY_X  = DifBeta  ( 2, 1 ) ;
DBetaY_Z  = DifBeta  ( 2, 3 ) ;

DBetaZ_X  = DifBeta  ( 3, 1 ) ;
DBetaZ_Y  = DifBeta  ( 3, 2 ) ;

DAlphaX_Y  = DifAlpha  ( 1, 2 ) ;
DAlphaX_Z  = DifAlpha  ( 1, 3 ) ;

DAlphaY_X  = DifAlpha  ( 2, 1 ) ;
DAlphaY_Z  = DifAlpha  ( 2, 3 ) ;

DAlphaZ_X  = DifAlpha  ( 3, 1 ) ;
DAlphaZ_Y  = DifAlpha  ( 3, 2 ) ;

A = Alpha_X * Alpha_Y * Alpha_Z ;

B = Alpha_X * Alpha_Y * Beta_Z + Alpha_X * Beta_Y * Alpha_Z + Beta_X * Alpha_Y * Alpha_Z ;

C = Alpha_X * Beta_Y * Beta_Z + Beta_X * Alpha_Y * Beta_Z + Beta_X * Beta_Y * Alpha_Z ;

D = Beta_X * Beta_Y * Beta_Z ;

FAC_A_Rho = A * Rho ;
FAC_B_Rho = B * Rho ;
FAC_C_Rho = C * Rho ;
FAC_D_Rho = D * Rho ;

Coef_PML ( 1, 1 ) = Alpha_Y * Alpha_Z ;
Coef_PML ( 1, 2 ) = Alpha_X * Alpha_Z ;
Coef_PML ( 1, 3 ) = Alpha_Y * Alpha_X ;

Coef_PML ( 3, 1 ) = Beta_Y * Beta_Z ;
Coef_PML ( 3, 2 ) = Beta_X * Beta_Z ;
Coef_PML ( 3, 3 ) = Beta_Y * Beta_X ;

Coef_PML ( 2, 1 ) = Alpha_Y * Beta_Z + Beta_Y * Alpha_Z ;
Coef_PML ( 2, 2 ) = Alpha_X * Beta_Z + Beta_X * Alpha_Z ;
Coef_PML ( 2, 3 ) = Alpha_Y * Beta_X + Beta_Y * Alpha_X ;


D_Coef_PML ( 1, 1 ) = DAlphaY_X * Alpha_Z + Alpha_y * DAlphaZ_X ;
D_Coef_PML ( 1, 2 ) = DAlphaX_Y * Alpha_Z + Alpha_X * DAlphaZ_Y ;
D_Coef_PML ( 1, 3 ) = DAlphaY_Z * Alpha_X + Alpha_y * DAlphaX_Z ;

D_Coef_PML ( 3, 1 ) = DBetaY_X * Beta_Z + Beta_y * DBetaZ_X ;
D_Coef_PML ( 3, 2 ) = DBetaX_Y * Beta_Z + Beta_X * DBetaZ_Y ;
D_Coef_PML ( 3, 3 ) = DBetaY_Z * Beta_X + Beta_y * DBetaX_Z ;

D_Coef_PML ( 2, 1 ) = DAlphaY_X * Beta_Z + Alpha_y * DBetaZ_X  +  DBetaY_X * Alpha_Z + Beta_y * DAlphaZ_X ;
D_Coef_PML ( 2, 2 ) = DAlphaX_Y * Beta_Z + Alpha_X * DBetaZ_Y  +  DBetaX_Y * Alpha_Z + Beta_X * DAlphaZ_Y ;
D_Coef_PML ( 2, 3 ) = DAlphaY_Z * Beta_X + Alpha_y * DBetaX_Z  +  DBetaY_Z * Alpha_X + Beta_y * DAlphaX_Z ;

!#Write(*    ,*) 'End Subroutine < MPML_Parameters_3D >' ;
!#Write(UnInf,*) 'End Subroutine < MPML_Parameters_3D >' ;
Return ;
End Subroutine MPML_Parameters_3D ;

End Module PML_PRMTRS ;
