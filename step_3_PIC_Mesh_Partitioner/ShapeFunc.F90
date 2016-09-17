
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                                                                                                                  ++
! Start:        01 August 2004                                                                                                                     ++
! Last Update:  16 Aug 2011                                                                                                                        ++
!                                                                                                                                                  ++
! Description: THIS Module Contains ALL SHAPE FUNCTIONS AND DIFFERENTIALS OF SHAPE FUNCTIONS.                                                      ++
!                                       -----------------------------------------------------                                                      ++
!                                                                                                                                                  ++
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Module ShapeFunctions ;

Use Parameters ;

Implicit None ;

! =========================== PUBLIC Variables ======================================================================================================

! General Shape Functions
Type, Public  :: ShapeF ;
  Real (Kind=DBL), Allocatable, Dimension(:) :: X ;   ! Coordinates
  Real (Kind=DBL), Allocatable, Dimension(:) :: FN ;  ! Shape Funcions
End Type ShapeF ;

! General Differential of Shape Functions
Type, Public  :: DShapeF ;
  Real (Kind=DBL), Allocatable, Dimension(:)   :: X ;   ! Coordinates
  Real (Kind=DBL), Allocatable, Dimension(:,:) :: DFXI ;! Shape Funcions
End Type DShapeF ;

! - 27-Node 3D Hexahedral Element -------------------------------------------------------------------------------------------------------------------
! Shape Functions
Type, Public ::  SF_3_27 ;
  Real (Kind=DBL) :: X1,  X2, X3 ;
  Real (Kind=DBL), Dimension ( 27 )  :: FN ;
End Type SF_3_27 ;

! Differentials
Type, Public ::  DSF_3_27 ;
  Real (Kind=DBL) :: X1, X2, X3 ;
  Real (Kind=DBL), Dimension ( 27, 3 )  :: DFXI ;
End Type DSF_3_27 ;

! - 20-Node 3D Hexahedral Element -------------------------------------------------------------------------------------------------------------------
! Shape Functions
Type, Public ::  SF_3_20 ;
  Real (Kind=DBL) :: X1,  X2, X3 ;
  Real (Kind=DBL), Dimension ( 20 )  :: FN ;
End Type SF_3_20 ;

! Differentials
Type, Public ::  DSF_3_20 ;
  Real (Kind=DBL) :: X1, X2, X3 ;
  Real (Kind=DBL), Dimension ( 20, 3 )  :: DFXI ;
End Type DSF_3_20 ;

! - 8-Node 3D Hexahedral Element --------------------------------------------------------------------------------------------------------------------
! Shape Functions
Type, Public ::  SF_3_8 ;
  Real (Kind=DBL) :: X1,  X2, X3 ;
  Real (Kind=DBL), Dimension ( 8 )  :: FN ;
End Type SF_3_8 ;

! Differentials
Type, Public ::  DSF_3_8 ;
  Real (Kind=DBL) :: X1, X2, X3 ;
  Real (Kind=DBL), Dimension ( 8 , 3 )  :: DFXI ;
End Type DSF_3_8 ;

! - 8-Node 2D Quadrilateral Element -----------------------------------------------------------------------------------------------------------------
! Shape Functions
Type, Public ::  SF_2_8 ;
  Real (Kind=DBL) :: X1, X2 ;
  Real (Kind=DBL), Dimension ( 8 )  :: FN ;
End Type SF_2_8 ;

! Differentials
Type, Public ::  DSF_2_8 ;
  Real (Kind=DBL) :: X1, X2 ;
  Real (Kind=DBL), Dimension ( 8 , 2 )  :: DFXI ;
End Type DSF_2_8 ;

! - 9-Node 2D Quadrilateral Spectal Element ---------------------------------------------------------------------------------------------------------
! Shape Functions
Type, Public ::  SF_2_9 ;
  Real (Kind=DBL) :: X1, X2 ;
  Real (Kind=DBL), Dimension ( 9 )  :: FN ;
End Type SF_2_9 ;

! Differentials
Type, Public ::  DSF_2_9 ;
  Real (Kind=DBL) :: X1, X2 ;
  Real (Kind=DBL), Dimension ( 9 , 2 )  :: DFXI ;
End Type DSF_2_9 ;

! - 4-Node 2D Quadrilateral Element -----------------------------------------------------------------------------------------------------------------
! Shape Functions
Type, Public ::  SF_2_4 ;
  Real (Kind=DBL) :: X1, X2 ;
  Real (Kind=DBL), Dimension ( 4 )  :: FN ;
End Type SF_2_4 ;

! Differentials
Type, Public ::  DSF_2_4 ;
  Real (Kind=DBL) :: X1, X2 ;
  Real (Kind=DBL), Dimension ( 4 , 2 )  :: DFXI ;
End Type DSF_2_4 ;

! - 6-Node 2D Triangle Element ----------------------------------------------------------------------------------------------------------------------
! Shape Functions 
Type, Public ::  SF_2_6 ;
  Real (Kind=DBL) :: X1, X2 ;
  Real (Kind=DBL), Dimension ( 6 )  :: FN ;
End Type SF_2_6 ;

! Differentials
Type, Public ::  DSF_2_6 ;
  Real (Kind=DBL) :: X1, X2 ;
  Real (Kind=DBL), Dimension ( 6 , 2 )  :: DFXI ;
End Type DSF_2_6 ;

! - 7-Node 2D Triangle Element ----------------------------------------------------------------------------------------------------------------------
! Shape Functions 
Type, Public ::  SF_2_7 ;
  Real (Kind=DBL) :: X1, X2 ;
  Real (Kind=DBL), Dimension ( 7 )  :: FN ;
End Type SF_2_7 ;

! Differentials
Type, Public ::  DSF_2_7 ;
  Real (Kind=DBL) :: X1, X2 ;
  Real (Kind=DBL), Dimension ( 7 , 2 )  :: DFXI ;
End Type DSF_2_7 ;

! - 3-Node 2D Triangle Element ----------------------------------------------------------------------------------------------------------------------
! Shape Functions 
Type, Public ::  SF_2_3 ;
  Real (Kind=DBL) :: r, s ;
  Real (Kind=DBL), Dimension ( 3 )  :: FN ;
End Type SF_2_3 ;

! Differentials
Type, Public ::  DSF_2_3 ;
  Real (Kind=DBL) :: r, s ;
  Real (Kind=DBL), Dimension ( 3 , 2 )  :: DFXI ;
End Type DSF_2_3 ;

! - 3-Node 1D Line Element --------------------------------------------------------------------------------------------------------------------------
! Shape Functions 
Type, Public ::  SF_1_3 ;
  Real (Kind=DBL) :: X1 ;
  Real (Kind=DBL), Dimension ( 3 )  :: FN ;
End Type SF_1_3 ;

! Differentials
Type, Public ::  DSF_1_3 ;
  Real (Kind=DBL) :: X1 ;
  Real (Kind=DBL), Dimension ( 3 , 1 )  :: DFXI ;
End Type DSF_1_3 ;


! =========================== PRIVATE Variables =====================================================================================================


! =========================== SHAPE FUNCTIONS Interface =============================================================================================
  Interface ShapeFuncSub ;
    Module Procedure  SH_FUNC_3_20, SH_FUNC_3_8, SH_FUNC_2_8, SH_FUNC_2_4, SH_FUNC_2_6 ; !, ShapeFunc_2_6, SH_FUNC_2_3 ;  ! SF_3_4 ,  SF_3_ 10 
  End Interface ShapeFuncSub ;

!  Interface ShapeFuncs ;
!    Module Procedure  ShapeFunc_2_6, ShapeFunc_2_8 ;   !, SH_FUNC_2_3 ;  ! SF_3_4 ,  SF_3_ 10 
!  End Interface ShapeFuncs ;
  

! =========================== DIFFERENTIALS OF SHAPE FUNCTIONS ======================================================================================
  Interface DIF_ShapeFuncSub ;
    Module Procedure  DIF_SH_FUNC_3_20, DIF_SH_FUNC_3_8, DIF_SH_FUNC_2_8, DIF_SH_FUNC_2_4, Dif_SH_FUNC_2_6 ;  !, DIF_ShapeFunc_2_6, DIF_SH_FUNC_2_3 ;  ! DSF_3_4 ,  DSF_3_ 10
  End Interface DIF_ShapeFuncSub ;

!  Interface DifShapeFuncs ;
!    Module Procedure  DIF_ShapeFunc_2_6, Dif_ShapeFunc_2_8  !, DIF_SH_FUNC_2_3 ;  ! DSF_3_4 ,  DSF_3_ 10
!  End Interface DifShapeFuncs ;

Contains ;

! ***************************************************************************************************************************************************
!                                                                                                                                                  **
! Last Update:  31 Jan 2012                                                                                                                        **
! Description: THIS Function calculates the shape functions of a                                                                                   **
!              { 2D SECOND ORDER  6-NODE SerEndipity ELEMENT }. ELEMENT NUMBER : { 10 OR 12 }.                                                     **
!                                                                                                                                                  **
! ***************************************************************************************************************************************************

Subroutine SH_FUNC_2_6 ( SF ) ;

Implicit None ;

Type ( SF_2_6 ) :: SF ;


! =========================== Function Code =========================================================================================================

SF%FN ( 1 ) = ( 1.0_DBL - SF%X1 - SF%X2 ) * ( 1.0_DBL - 2.0_DBL * SF%X1 - 2.0_DBL * SF%X2 ) ; 
SF%FN ( 2 ) = SF%X1 * ( 2.0_DBL * SF%X1 - 1.0_DBL ) ; 
SF%FN ( 3 ) = SF%X2 * ( 2.0_DBL * SF%X2 - 1.0_DBL ) ; 
SF%FN ( 4 ) = 4.0_DBL * SF%X1 * ( 1.0_DBL - SF%X1 - SF%X2 ) ;
SF%FN ( 5 ) = 4.0_DBL * SF%X1 * SF%X2 ;
SF%FN ( 6 ) = 4.0_DBL * SF%X2 * ( 1.0_DBL - SF%X1 - SF%X2 ) ;

End Subroutine SH_FUNC_2_6 ;


! ***************************************************************************************************************************************************
!                                                                                                                                                  **
! Last Update:  31 Jan 2011                                                                                                                        **
! Description: This Function calculates the differentials of the shape functions of a                                                              **
!              { 2D Second ORDER  3-NODE SerEndipity ELEMENT }. ELEMENT NUMBER : { 10 0R 12 }.                                                     **
!                                                                                                                                                  **
! ***************************************************************************************************************************************************

Subroutine DIF_SH_FUNC_2_6 ( DSF ) ;

Implicit None ;

Type ( DSF_2_6 ) :: DSF ;

! =========================== Subroutine CODE =======================================================================================================

DSF%DFXI ( 1 , 1 ) = - 3.0_DBL + 4.0_DBL * DSF%X1 + 4.0_DBL * DSF%X2 ; 
DSF%DFXI ( 1 , 2 ) = - 3.0_DBL + 4.0_DBL * DSF%X1 + 4.0_DBL * DSF%X2 ; 
DSF%DFXI ( 2 , 1 ) = + 4.0_DBL * DSF%X1 - 1.0_DBL ; 
DSF%DFXI ( 2 , 2 ) =   0.0_DBL ; 
DSF%DFXI ( 3 , 1 ) =   0.0_DBL ; 
DSF%DFXI ( 3 , 2 ) = + 4.0_DBL * DSF%X2 - 1.0_DBL ; 
DSF%DFXI ( 4 , 1 ) = + 4.0_DBL - 8.0_DBL * DSF%X1 - 4.0_DBL * DSF%X2 ; 
DSF%DFXI ( 4 , 2 ) = - 4.0_DBL * DSF%X1 ; 
DSF%DFXI ( 5 , 1 ) = + 4.0_DBL * DSF%X2 ; 
DSF%DFXI ( 5 , 2 ) = + 4.0_DBL * DSF%X1 ; 
DSF%DFXI ( 6 , 1 ) = - 4.0_DBL * DSF%X2 ; 
DSF%DFXI ( 6 , 2 ) = + 4.0_DBL - 4.0_DBL * DSF%X1 - 8.0_DBL * DSF%X2 ; 

End Subroutine DIF_SH_FUNC_2_6 ;

! ***************************************************************************************************************************************************
!                                                                                                                                                  **
! Last Update:  22 MAY 2011                                                                                                                        **
! Description: THIS Subroutines CALCULATES THE SHAPE FUNCTIONS OF A  { 3D SECOND ORDER 20-NODE SEREndEPITI ELEMENT }. ELEMENT NUMBER : { 6 }.      **
!                                                                     -----------------------------------------------------------------------      **
!                                                                                                                                                  **
! ***************************************************************************************************************************************************

Subroutine SH_FUNC_3_20 ( SF ) ;


Implicit None ;

Type ( SF_3_20 ) :: SF ;

! =========================== Subroutine CODE =======================================================================================================

SF%FN ( 1  ) = ( 1.0_DBL + SF%X1 ) * ( 1.0_DBL - SF%X2 ) * ( 1.0_DBL - SF%X3 ) * ( -2.0_DBL + SF%X1 - SF%X2 - SF%X3 ) * 0.125_DBL ;
SF%FN ( 2  ) = ( 1.0_DBL + SF%X1 ) * ( 1.0_DBL + SF%X2 ) * ( 1.0_DBL - SF%X3 ) * ( -2.0_DBL + SF%X1 + SF%X2 - SF%X3 ) * 0.125_DBL ;
SF%FN ( 3  ) = ( 1.0_DBL - SF%X1 ) * ( 1.0_DBL + SF%X2 ) * ( 1.0_DBL - SF%X3 ) * ( -2.0_DBL - SF%X1 + SF%X2 - SF%X3 ) * 0.125_DBL ;
SF%FN ( 4  ) = ( 1.0_DBL - SF%X1 ) * ( 1.0_DBL - SF%X2 ) * ( 1.0_DBL - SF%X3 ) * ( -2.0_DBL - SF%X1 - SF%X2 - SF%X3 ) * 0.125_DBL ;
SF%FN ( 5  ) = ( 1.0_DBL + SF%X1 ) * ( 1.0_DBL - SF%X2 ) * ( 1.0_DBL + SF%X3 ) * ( -2.0_DBL + SF%X1 - SF%X2 + SF%X3 ) * 0.125_DBL ;
SF%FN ( 6  ) = ( 1.0_DBL + SF%X1 ) * ( 1.0_DBL + SF%X2 ) * ( 1.0_DBL + SF%X3 ) * ( -2.0_DBL + SF%X1 + SF%X2 + SF%X3 ) * 0.125_DBL ;
SF%FN ( 7  ) = ( 1.0_DBL - SF%X1 ) * ( 1.0_DBL + SF%X2 ) * ( 1.0_DBL + SF%X3 ) * ( -2.0_DBL - SF%X1 + SF%X2 + SF%X3 ) * 0.125_DBL ;
SF%FN ( 8  ) = ( 1.0_DBL - SF%X1 ) * ( 1.0_DBL - SF%X2 ) * ( 1.0_DBL + SF%X3 ) * ( -2.0_DBL - SF%X1 - SF%X2 + SF%X3 ) * 0.125_DBL ;

SF%FN ( 9  ) = ( 1.0_DBL - SF%X2 * SF%X2 ) * ( 1.0_DBL + SF%X1 ) * ( 1.0_DBL - SF%X3 ) * 0.25_DBL ;
SF%FN ( 10 ) = ( 1.0_DBL - SF%X1 * SF%X1 ) * ( 1.0_DBL + SF%X2 ) * ( 1.0_DBL - SF%X3 ) * 0.25_DBL ;
SF%FN ( 11 ) = ( 1.0_DBL - SF%X2 * SF%X2 ) * ( 1.0_DBL - SF%X1 ) * ( 1.0_DBL - SF%X3 ) * 0.25_DBL ;
SF%FN ( 12 ) = ( 1.0_DBL - SF%X1 * SF%X1 ) * ( 1.0_DBL - SF%X2 ) * ( 1.0_DBL - SF%X3 ) * 0.25_DBL ;
SF%FN ( 13 ) = ( 1.0_DBL - SF%X3 * SF%X3 ) * ( 1.0_DBL + SF%X1 ) * ( 1.0_DBL - SF%X2 ) * 0.25_DBL ;
SF%FN ( 14 ) = ( 1.0_DBL - SF%X3 * SF%X3 ) * ( 1.0_DBL + SF%X1 ) * ( 1.0_DBL + SF%X2 ) * 0.25_DBL ;
SF%FN ( 15 ) = ( 1.0_DBL - SF%X3 * SF%X3 ) * ( 1.0_DBL - SF%X1 ) * ( 1.0_DBL + SF%X2 ) * 0.25_DBL ;
SF%FN ( 16 ) = ( 1.0_DBL - SF%X3 * SF%X3 ) * ( 1.0_DBL - SF%X1 ) * ( 1.0_DBL - SF%X2 ) * 0.25_DBL ;
SF%FN ( 17 ) = ( 1.0_DBL - SF%X2 * SF%X2 ) * ( 1.0_DBL + SF%X1 ) * ( 1.0_DBL + SF%X3 ) * 0.25_DBL ;
SF%FN ( 18 ) = ( 1.0_DBL - SF%X1 * SF%X1 ) * ( 1.0_DBL + SF%X2 ) * ( 1.0_DBL + SF%X3 ) * 0.25_DBL ;
SF%FN ( 19 ) = ( 1.0_DBL - SF%X2 * SF%X2 ) * ( 1.0_DBL - SF%X1 ) * ( 1.0_DBL + SF%X3 ) * 0.25_DBL ;
SF%FN ( 20 ) = ( 1.0_DBL - SF%X1 * SF%X1 ) * ( 1.0_DBL - SF%X2 ) * ( 1.0_DBL + SF%X3 ) * 0.25_DBL ;


!Write(*,"('End Subroutine < SH_FUNC_3_20 >')")  ;
!Write(UnInf,"('End Subroutine < SH_FUNC_3_20 >')") ;
Return ;
End Subroutine SH_FUNC_3_20 ;


! ***************************************************************************************************************************************************
!                                                                                                                                                  **
! Last Update:  22 MAY 2011                                                                                                                        **
! Description: THIS Subroutines CALCULATES THE DIFFERENTIALS OF SHAPE FUNCTIONS OF  A                                                              **
!              { 3D SECOND ORDER 20-NODE SEREndEPITI ELEMENT }. ELEMENT NUMBER : { 6 }.                                                            **
!                                                                                                                                                  **
! ***************************************************************************************************************************************************


Subroutine DIF_SH_FUNC_3_20 ( DSF ) ;


Implicit None ;

Type ( DSF_3_20 ) :: DSF ;


! =========================== Subroutine CODE =======================================================================================================

DSF%DFXI ( 1  , 1 ) = ( 1.0_DBL - DSF%X2 ) * ( 1.0_DBL - DSF%X3 ) * ( -1.0_DBL + 2 * DSF%X1 - DSF%X2 - DSF%X3 ) * 0.125_DBL ; ! dfn / dx
DSF%DFXI ( 1  , 2 ) = ( 1.0_DBL + DSF%X1 ) * ( 1.0_DBL - DSF%X3 ) * (  1.0_DBL - DSF%X1 + 2 * DSF%X2 + DSF%X3 ) * 0.125_DBL ; ! dfn / dy
DSF%DFXI ( 1  , 3 ) = ( 1.0_DBL + DSF%X1 ) * ( 1.0_DBL - DSF%X2 ) * (  1.0_DBL - DSF%X1 + DSF%X2 + 2 * DSF%X3 ) * 0.125_DBL ; ! dfn / dz

DSF%DFXI ( 2  , 1 ) = ( 1.0_DBL + DSF%X2 ) * ( 1.0_DBL - DSF%X3 ) * ( -1.0_DBL + 2 * DSF%X1 + DSF%X2 - DSF%X3 ) * 0.125_DBL ;
DSF%DFXI ( 2  , 2 ) = ( 1.0_DBL + DSF%X1 ) * ( 1.0_DBL - DSF%X3 ) * ( -1.0_DBL + DSF%X1 + 2 * DSF%X2 - DSF%X3 ) * 0.125_DBL ;
DSF%DFXI ( 2  , 3 ) = ( 1.0_DBL + DSF%X1 ) * ( 1.0_DBL + DSF%X2 ) * (  1.0_DBL - DSF%X1 - DSF%X2 + 2 * DSF%X3 ) * 0.125_DBL ;

DSF%DFXI ( 3  , 1 ) = ( 1.0_DBL + DSF%X2 ) * ( 1.0_DBL - DSF%X3 ) * (  1.0_DBL + 2 * DSF%X1 - DSF%X2 + DSF%X3 ) * 0.125_DBL ;
DSF%DFXI ( 3  , 2 ) = ( 1.0_DBL - DSF%X1 ) * ( 1.0_DBL - DSF%X3 ) * ( -1.0_DBL - DSF%X1 + 2 * DSF%X2 - DSF%X3 ) * 0.125_DBL ;
DSF%DFXI ( 3  , 3 ) = ( 1.0_DBL - DSF%X1 ) * ( 1.0_DBL + DSF%X2 ) * (  1.0_DBL + DSF%X1 - DSF%X2 + 2 * DSF%X3 ) * 0.125_DBL ;

DSF%DFXI ( 4  , 1 ) = ( 1.0_DBL - DSF%X2 ) * ( 1.0_DBL - DSF%X3 ) * (  1.0_DBL + 2 * DSF%X1 + DSF%X2 + DSF%X3 ) * 0.125_DBL ;
DSF%DFXI ( 4  , 2 ) = ( 1.0_DBL - DSF%X1 ) * ( 1.0_DBL - DSF%X3 ) * (  1.0_DBL + DSF%X1 + 2 * DSF%X2 + DSF%X3 ) * 0.125_DBL ;
DSF%DFXI ( 4  , 3 ) = ( 1.0_DBL - DSF%X1 ) * ( 1.0_DBL - DSF%X2 ) * (  1.0_DBL + DSF%X1 + DSF%X2 + 2 * DSF%X3 ) * 0.125_DBL ;

DSF%DFXI ( 5  , 1 ) = ( 1.0_DBL - DSF%X2 ) * ( 1.0_DBL + DSF%X3 ) * ( -1.0_DBL + 2 * DSF%X1 - DSF%X2 + DSF%X3 ) * 0.125_DBL ;
DSF%DFXI ( 5  , 2 ) = ( 1.0_DBL + DSF%X1 ) * ( 1.0_DBL + DSF%X3 ) * (  1.0_DBL - DSF%X1 + 2 * DSF%X2 - DSF%X3 ) * 0.125_DBL ;
DSF%DFXI ( 5  , 3 ) = ( 1.0_DBL + DSF%X1 ) * ( 1.0_DBL - DSF%X2 ) * ( -1.0_DBL + DSF%X1 - DSF%X2 + 2 * DSF%X3 ) * 0.125_DBL ;

DSF%DFXI ( 6  , 1 ) = ( 1.0_DBL + DSF%X2 ) * ( 1.0_DBL + DSF%X3 ) * ( -1.0_DBL + 2 * DSF%X1 + DSF%X2 + DSF%X3 ) * 0.125_DBL ;
DSF%DFXI ( 6  , 2 ) = ( 1.0_DBL + DSF%X1 ) * ( 1.0_DBL + DSF%X3 ) * ( -1.0_DBL + DSF%X1 + 2 * DSF%X2 + DSF%X3 ) * 0.125_DBL ;
DSF%DFXI ( 6  , 3 ) = ( 1.0_DBL + DSF%X1 ) * ( 1.0_DBL + DSF%X2 ) * ( -1.0_DBL + DSF%X1 + DSF%X2 + 2 * DSF%X3 ) * 0.125_DBL ;

DSF%DFXI ( 7  , 1 ) = ( 1.0_DBL + DSF%X2 ) * ( 1.0_DBL + DSF%X3 ) * (  1.0_DBL + 2 * DSF%X1 - DSF%X2 - DSF%X3 ) * 0.125_DBL ;
DSF%DFXI ( 7  , 2 ) = ( 1.0_DBL - DSF%X1 ) * ( 1.0_DBL + DSF%X3 ) * ( -1.0_DBL - DSF%X1 + 2 * DSF%X2 + DSF%X3 ) * 0.125_DBL ;
DSF%DFXI ( 7  , 3 ) = ( 1.0_DBL - DSF%X1 ) * ( 1.0_DBL + DSF%X2 ) * ( -1.0_DBL - DSF%X1 + DSF%X2 + 2 * DSF%X3 ) * 0.125_DBL ;

DSF%DFXI ( 8  , 1 ) = ( 1.0_DBL - DSF%X2 ) * ( 1.0_DBL + DSF%X3 ) * (  1.0_DBL + 2 * DSF%X1 + DSF%X2 - DSF%X3 ) * 0.125_DBL ;
DSF%DFXI ( 8  , 2 ) = ( 1.0_DBL - DSF%X1 ) * ( 1.0_DBL + DSF%X3 ) * (  1.0_DBL + DSF%X1 + 2 * DSF%X2 - DSF%X3 ) * 0.125_DBL ;
DSF%DFXI ( 8  , 3 ) = ( 1.0_DBL - DSF%X1 ) * ( 1.0_DBL - DSF%X2 ) * ( -1.0_DBL - DSF%X1 - DSF%X2 + 2 * DSF%X3 ) * 0.125_DBL ;


DSF%DFXI ( 9  , 1 ) =  ( 1.0_DBL -  DSF%X2 * DSF%X2 ) * ( 1.0_DBL - DSF%X3 ) * 0.25_DBL ; 
DSF%DFXI ( 9  , 2 ) = - DSF%X2 * ( 1.0_DBL + DSF%X1 ) * ( 1.0_DBL - DSF%X3 ) * 0.5_DBL  ; 
DSF%DFXI ( 9  , 3 ) = -( 1.0_DBL -  DSF%X2 * DSF%X2 ) * ( 1.0_DBL + DSF%X1 ) * 0.25_DBL ; 

DSF%DFXI ( 10 , 1 ) = - DSF%X1 * ( 1.0_DBL + DSF%X2 ) * ( 1.0_DBL - DSF%X3 ) * 0.5_DBL  ; 
DSF%DFXI ( 10 , 2 ) =  ( 1.0_DBL -  DSF%X1 * DSF%X1 ) * ( 1.0_DBL - DSF%X3 ) * 0.25_DBL ; 
DSF%DFXI ( 10 , 3 ) = -( 1.0_DBL -  DSF%X1 * DSF%X1 ) * ( 1.0_DBL + DSF%X2 ) * 0.25_DBL ; 

DSF%DFXI ( 11 , 1 ) = -( 1.0_DBL -  DSF%X2 * DSF%X2 ) * ( 1.0_DBL - DSF%X3 ) * 0.25_DBL ; 
DSF%DFXI ( 11 , 2 ) = - DSF%X2 * ( 1.0_DBL - DSF%X1 ) * ( 1.0_DBL - DSF%X3 ) * 0.5_DBL  ; 
DSF%DFXI ( 11 , 3 ) = -( 1.0_DBL -  DSF%X2 * DSF%X2 ) * ( 1.0_DBL - DSF%X1 ) * 0.25_DBL ; 

DSF%DFXI ( 12 , 1 ) = - DSF%X1 * ( 1.0_DBL - DSF%X2 ) * ( 1.0_DBL - DSF%X3 ) * 0.5_DBL  ; 
DSF%DFXI ( 12 , 2 ) = -( 1.0_DBL -  DSF%X1 * DSF%X1 ) * ( 1.0_DBL - DSF%X3 ) * 0.25_DBL ; 
DSF%DFXI ( 12 , 3 ) = -( 1.0_DBL -  DSF%X1 * DSF%X1 ) * ( 1.0_DBL - DSF%X2 ) * 0.25_DBL ; 

DSF%DFXI ( 13 , 1 ) =  ( 1.0_DBL -  DSF%X3 * DSF%X3 ) * ( 1.0_DBL - DSF%X2 ) * 0.25_DBL ; 
DSF%DFXI ( 13 , 2 ) = -( 1.0_DBL -  DSF%X3 * DSF%X3 ) * ( 1.0_DBL + DSF%X1 ) * 0.25_DBL ; 
DSF%DFXI ( 13 , 3 ) = - DSF%X3 * ( 1.0_DBL + DSF%X1 ) * ( 1.0_DBL - DSF%X2 ) * 0.5_DBL ; 

DSF%DFXI ( 14 , 1 ) =  ( 1.0_DBL -  DSF%X3 * DSF%X3 ) * ( 1.0_DBL + DSF%X2 ) * 0.25_DBL ; 
DSF%DFXI ( 14 , 2 ) =  ( 1.0_DBL -  DSF%X3 * DSF%X3 ) * ( 1.0_DBL + DSF%X1 ) * 0.25_DBL ; 
DSF%DFXI ( 14 , 3 ) = - DSF%X3 * ( 1.0_DBL + DSF%X1 ) * ( 1.0_DBL + DSF%X2 ) * 0.5_DBL ; 

DSF%DFXI ( 15 , 1 ) = -( 1.0_DBL -  DSF%X3 * DSF%X3 ) * ( 1.0_DBL + DSF%X2 ) * 0.25_DBL ; 
DSF%DFXI ( 15 , 2 ) =  ( 1.0_DBL -  DSF%X3 * DSF%X3 ) * ( 1.0_DBL - DSF%X1 ) * 0.25_DBL ; 
DSF%DFXI ( 15 , 3 ) = - DSF%X3 * ( 1.0_DBL - DSF%X1 ) * ( 1.0_DBL + DSF%X2 ) * 0.5_DBL  ; 

DSF%DFXI ( 16 , 1 ) = -( 1.0_DBL -  DSF%X3 * DSF%X3 ) * ( 1.0_DBL - DSF%X2 ) * 0.25_DBL ; 
DSF%DFXI ( 16 , 2 ) = -( 1.0_DBL -  DSF%X3 * DSF%X3 ) * ( 1.0_DBL - DSF%X1 ) * 0.25_DBL ; 
DSF%DFXI ( 16 , 3 ) = - DSF%X3 * ( 1.0_DBL - DSF%X1 ) * ( 1.0_DBL - DSF%X2 ) * 0.5_DBL  ; 

DSF%DFXI ( 17 , 1 ) =  ( 1.0_DBL -  DSF%X2 * DSF%X2 ) * ( 1.0_DBL + DSF%X3 ) * 0.25_DBL ; 
DSF%DFXI ( 17 , 2 ) = - DSF%X2 * ( 1.0_DBL + DSF%X1 ) * ( 1.0_DBL + DSF%X3 ) * 0.5_DBL  ; 
DSF%DFXI ( 17 , 3 ) =  ( 1.0_DBL -  DSF%X2 * DSF%X2 ) * ( 1.0_DBL + DSF%X1 ) * 0.25_DBL ; 

DSF%DFXI ( 18 , 1 ) = - DSF%X1 * ( 1.0_DBL + DSF%X2 ) * ( 1.0_DBL + DSF%X3 ) * 0.5_DBL  ; 
DSF%DFXI ( 18 , 2 ) =  ( 1.0_DBL -  DSF%X1 * DSF%X1 ) * ( 1.0_DBL + DSF%X3 ) * 0.25_DBL ; 
DSF%DFXI ( 18 , 3 ) =  ( 1.0_DBL -  DSF%X1 * DSF%X1 ) * ( 1.0_DBL + DSF%X2 ) * 0.25_DBL ; 

DSF%DFXI ( 19 , 1 ) = -( 1.0_DBL -  DSF%X2 * DSF%X2 ) * ( 1.0_DBL + DSF%X3 ) * 0.25_DBL ; 
DSF%DFXI ( 19 , 2 ) = - DSF%X2 * ( 1.0_DBL - DSF%X1 ) * ( 1.0_DBL + DSF%X3 ) * 0.5_DBL  ; 
DSF%DFXI ( 19 , 3 ) =  ( 1.0_DBL -  DSF%X2 * DSF%X2 ) * ( 1.0_DBL - DSF%X1 ) * 0.25_DBL ; 

DSF%DFXI ( 20 , 1 ) = - DSF%X1 * ( 1.0_DBL - DSF%X2 ) * ( 1.0_DBL + DSF%X3 ) * 0.5_DBL  ; 
DSF%DFXI ( 20 , 2 ) = -( 1.0_DBL -  DSF%X1 * DSF%X1 ) * ( 1.0_DBL + DSF%X3 ) * 0.25_DBL ; 
DSF%DFXI ( 20 , 3 ) =  ( 1.0_DBL -  DSF%X1 * DSF%X1 ) * ( 1.0_DBL - DSF%X2 ) * 0.25_DBL ; 


!Write(*,"('End Subroutine < DIF_SH_FUNC_3_20 >')")  ;
!Write(UnInf,"('End Subroutine < DIF_SH_FUNC_3_20 >')") ;
Return ;
End Subroutine DIF_SH_FUNC_3_20 ;


! ***************************************************************************************************************************************************
!                                                                                                                                                  **
! Last Update:  22 MAY 2011                                                                                                                        **
! Description: THIS Subroutines CALCULATES THE SHAPE FUNCTIONS OF A  { 3D FIRST  ORDER  8-NODE ELEMENT }. ELEMENT NUMBER : { 5 }.                  **
!                                                                    ------------------------------------------------------------                  **
!                                                                                                                                                  **
! ***************************************************************************************************************************************************

Subroutine SH_FUNC_3_8 ( SF ) ;


Implicit None ;

Type ( SF_3_8 ) :: SF ;

! =========================== Subroutine CODE =======================================================================================================

SF%FN ( 1 ) = ( 1.0_DBL + SF%X1 ) * ( 1.0_DBL - SF%X2 ) * ( 1.0_DBL - SF%X3 ) * 0.125_DBL ;
SF%FN ( 2 ) = ( 1.0_DBL + SF%X1 ) * ( 1.0_DBL + SF%X2 ) * ( 1.0_DBL - SF%X3 ) * 0.125_DBL ;
SF%FN ( 3 ) = ( 1.0_DBL - SF%X1 ) * ( 1.0_DBL + SF%X2 ) * ( 1.0_DBL - SF%X3 ) * 0.125_DBL ;
SF%FN ( 4 ) = ( 1.0_DBL - SF%X1 ) * ( 1.0_DBL - SF%X2 ) * ( 1.0_DBL - SF%X3 ) * 0.125_DBL ;
SF%FN ( 5 ) = ( 1.0_DBL + SF%X1 ) * ( 1.0_DBL - SF%X2 ) * ( 1.0_DBL + SF%X3 ) * 0.125_DBL ;
SF%FN ( 6 ) = ( 1.0_DBL + SF%X1 ) * ( 1.0_DBL + SF%X2 ) * ( 1.0_DBL + SF%X3 ) * 0.125_DBL ;
SF%FN ( 7 ) = ( 1.0_DBL - SF%X1 ) * ( 1.0_DBL + SF%X2 ) * ( 1.0_DBL + SF%X3 ) * 0.125_DBL ;
SF%FN ( 8 ) = ( 1.0_DBL - SF%X1 ) * ( 1.0_DBL - SF%X2 ) * ( 1.0_DBL + SF%X3 ) * 0.125_DBL ;

!Write(*,"('End Subroutine < SH_FUNC_3_8 >')")  ;
!Write(UnInf,"('End Subroutine < SH_FUNC_3_8 >')") ;
Return
End Subroutine SH_FUNC_3_8


! ***************************************************************************************************************************************************
!                                                                                                                                                  **
! Last Update:  22 MAY 2011                                                                                                                        **
! Description: THIS Subroutines CALCULATES THE DIFFERENTIALS OF SHAPE FUNCTIONS OF  A                                                              **
!              { 3D FIRST  ORDER  8-NODE ELEMENT }. ELEMENT NUMBER : { 5 }.                                                                        **
!                                                                                                                                                  **
! ***************************************************************************************************************************************************

Subroutine DIF_SH_FUNC_3_8 ( DSF ) ;


Implicit None ;

Type ( DSF_3_8 ) :: DSF ;

! =========================== Subroutine CODE =======================================================================================================

DSF%DFXI ( 8 , 1 ) = - ( 1.0_DBL - DSF%X2 ) * ( 1.0_DBL + DSF%X3 ) * 0.125_DBL ;
DSF%DFXI ( 8 , 2 ) = - ( 1.0_DBL - DSF%X1 ) * ( 1.0_DBL + DSF%X3 ) * 0.125_DBL ;
DSF%DFXI ( 8 , 3 ) = + ( 1.0_DBL - DSF%X1 ) * ( 1.0_DBL - DSF%X2 ) * 0.125_DBL ;
DSF%DFXI ( 7 , 1 ) = - ( 1.0_DBL + DSF%X2 ) * ( 1.0_DBL + DSF%X3 ) * 0.125_DBL ;
DSF%DFXI ( 7 , 2 ) = + ( 1.0_DBL - DSF%X1 ) * ( 1.0_DBL + DSF%X3 ) * 0.125_DBL ;
DSF%DFXI ( 7 , 3 ) = + ( 1.0_DBL - DSF%X1 ) * ( 1.0_DBL + DSF%X2 ) * 0.125_DBL ;
DSF%DFXI ( 6 , 1 ) = + ( 1.0_DBL + DSF%X2 ) * ( 1.0_DBL + DSF%X3 ) * 0.125_DBL ;
DSF%DFXI ( 6 , 2 ) = + ( 1.0_DBL + DSF%X1 ) * ( 1.0_DBL + DSF%X3 ) * 0.125_DBL ;
DSF%DFXI ( 6 , 3 ) = + ( 1.0_DBL + DSF%X1 ) * ( 1.0_DBL + DSF%X2 ) * 0.125_DBL ;
DSF%DFXI ( 5 , 1 ) = + ( 1.0_DBL - DSF%X2 ) * ( 1.0_DBL + DSF%X3 ) * 0.125_DBL ;
DSF%DFXI ( 5 , 2 ) = - ( 1.0_DBL + DSF%X1 ) * ( 1.0_DBL + DSF%X3 ) * 0.125_DBL ;
DSF%DFXI ( 5 , 3 ) = + ( 1.0_DBL + DSF%X1 ) * ( 1.0_DBL - DSF%X2 ) * 0.125_DBL ;
DSF%DFXI ( 4 , 1 ) = - ( 1.0_DBL - DSF%X2 ) * ( 1.0_DBL - DSF%X3 ) * 0.125_DBL ;
DSF%DFXI ( 4 , 2 ) = - ( 1.0_DBL - DSF%X1 ) * ( 1.0_DBL - DSF%X3 ) * 0.125_DBL ;
DSF%DFXI ( 4 , 3 ) = - ( 1.0_DBL - DSF%X1 ) * ( 1.0_DBL - DSF%X2 ) * 0.125_DBL ;
DSF%DFXI ( 3 , 1 ) = - ( 1.0_DBL + DSF%X2 ) * ( 1.0_DBL - DSF%X3 ) * 0.125_DBL ;
DSF%DFXI ( 3 , 2 ) = + ( 1.0_DBL - DSF%X1 ) * ( 1.0_DBL - DSF%X3 ) * 0.125_DBL ;
DSF%DFXI ( 3 , 3 ) = - ( 1.0_DBL - DSF%X1 ) * ( 1.0_DBL + DSF%X2 ) * 0.125_DBL ;
DSF%DFXI ( 2 , 1 ) = + ( 1.0_DBL + DSF%X2 ) * ( 1.0_DBL - DSF%X3 ) * 0.125_DBL ;
DSF%DFXI ( 2 , 2 ) = + ( 1.0_DBL + DSF%X1 ) * ( 1.0_DBL - DSF%X3 ) * 0.125_DBL ;
DSF%DFXI ( 2 , 3 ) = - ( 1.0_DBL + DSF%X1 ) * ( 1.0_DBL + DSF%X2 ) * 0.125_DBL ;
DSF%DFXI ( 1 , 1 ) = + ( 1.0_DBL - DSF%X2 ) * ( 1.0_DBL - DSF%X3 ) * 0.125_DBL ;
DSF%DFXI ( 1 , 2 ) = - ( 1.0_DBL + DSF%X1 ) * ( 1.0_DBL - DSF%X3 ) * 0.125_DBL ;
DSF%DFXI ( 1 , 3 ) = - ( 1.0_DBL + DSF%X1 ) * ( 1.0_DBL - DSF%X2 ) * 0.125_DBL ;

!Write(*,"('End Subroutine < SH_FUNC_3_20 >')")  ;
!Write(UnInf,"('End Subroutine < SH_FUNC_3_20 >')") ;
Return ;
End Subroutine DIF_SH_FUNC_3_8 ;


! ***************************************************************************************************************************************************
!                                                                                                                                                  **
! Last Update:  22 MAY 2011                                                                                                                        **
! Description: THIS Subroutines CALCULATES THE SHAPE FUNCTIONS OF A  { 2D SECOND ORDER  8-NODE SEREndEPITI ELEMENT }. ELEMENT NUMBER : { 2 0R 4 }  **
!                                                                     ---------------------------------------------------------------------------  **
!                                                                                                                                                  **
! ***************************************************************************************************************************************************

Subroutine SH_FUNC_2_8 ( SF ) ;


Implicit None ;
Type ( SF_2_8 ) :: SF ;

! =========================== Subroutine CODE =======================================================================================================

SF%FN ( 1 ) = ( 1.0_DBL + SF%X1 ) * ( 1.0_DBL - SF%X2 ) * ( -1.0_DBL + SF%X1 - SF%X2 ) * 0.25_DBL ; 
SF%FN ( 2 ) = ( 1.0_DBL + SF%X1 ) * ( 1.0_DBL + SF%X2 ) * ( -1.0_DBL + SF%X1 + SF%X2 ) * 0.25_DBL ; 
SF%FN ( 3 ) = ( 1.0_DBL - SF%X1 ) * ( 1.0_DBL + SF%X2 ) * ( -1.0_DBL - SF%X1 + SF%X2 ) * 0.25_DBL ; 
SF%FN ( 4 ) = ( 1.0_DBL - SF%X1 ) * ( 1.0_DBL - SF%X2 ) * ( -1.0_DBL - SF%X1 - SF%X2 ) * 0.25_DBL ; 
SF%FN ( 5 ) = ( 1.0_DBL - SF%X2 * SF%X2 ) * ( 1.0_DBL + SF%X1 ) * 0.5_DBL ; 
SF%FN ( 6 ) = ( 1.0_DBL - SF%X1 * SF%X1 ) * ( 1.0_DBL + SF%X2 ) * 0.5_DBL ; 
SF%FN ( 7 ) = ( 1.0_DBL - SF%X2 * SF%X2 ) * ( 1.0_DBL - SF%X1 ) * 0.5_DBL ; 
SF%FN ( 8 ) = ( 1.0_DBL - SF%X1 * SF%X1 ) * ( 1.0_DBL - SF%X2 ) * 0.5_DBL ; 


!Write(*,"('End Subroutine < SH_FUNC_3_20 >')")  ;
!Write(UnInf,"('End Subroutine < SH_FUNC_3_20 >')") ;
Return ;
End Subroutine SH_FUNC_2_8 ;


! ***************************************************************************************************************************************************
!                                                                                                                                                  **
! Last Update:  22 MAY 2011                                                                                                                        **
! Description: THIS Subroutines CALCULATES THE DIFFERENTIALS OF SHAPE FUNCTIONS OF  A                                                              **
!              { 2D SECOND ORDER  8-NODE SEREndEPITI ELEMENT }. ELEMENT NUMBER : { 2 OR 4 }.                                                       **
!                                                                                                                                                  **
! ***************************************************************************************************************************************************

Subroutine DIF_SH_FUNC_2_8 ( DSF ) ;


Implicit None ;

Type ( DSF_2_8 ) :: DSF ;

! =========================== Subroutine CODE =======================================================================================================

DSF%DFXI ( 1 , 1 ) =  ( 1.0_DBL - DSF%X2 ) * (   2.0_DBL * DSF%X1 - DSF%X2 ) * 0.25_DBL ; 
DSF%DFXI ( 1 , 2 ) =  ( 1.0_DBL + DSF%X1 ) * ( - DSF%X1 + 2.0_DBL * DSF%X2 ) * 0.25_DBL ; 
DSF%DFXI ( 2 , 1 ) =  ( 1.0_DBL + DSF%X2 ) * (   2.0_DBL * DSF%X1 + DSF%X2 ) * 0.25_DBL ; 
DSF%DFXI ( 2 , 2 ) =  ( 1.0_DBL + DSF%X1 ) * (   DSF%X1 + 2.0_DBL * DSF%X2 ) * 0.25_DBL ; 
DSF%DFXI ( 3 , 1 ) =  ( 1.0_DBL + DSF%X2 ) * (   2.0_DBL * DSF%X1 - DSF%X2 ) * 0.25_DBL ; 
DSF%DFXI ( 3 , 2 ) =  ( 1.0_DBL - DSF%X1 ) * ( - DSF%X1 + 2.0_DBL * DSF%X2 ) * 0.25_DBL ; 
DSF%DFXI ( 4 , 1 ) =  ( 1.0_DBL - DSF%X2 ) * (   2.0_DBL * DSF%X1 + DSF%X2 ) * 0.25_DBL ; 
DSF%DFXI ( 4 , 2 ) =  ( 1.0_DBL - DSF%X1 ) * (   DSF%X1 + 2.0_DBL * DSF%X2 ) * 0.25_DBL ; 

DSF%DFXI ( 5 , 1 ) =  ( 1.0_DBL - DSF%X2 * DSF%X2 ) * 0.5_DBL ; 
DSF%DFXI ( 5 , 2 ) = -( 1.0_DBL + DSF%X1 ) * DSF%X2 ;
DSF%DFXI ( 6 , 1 ) = -( 1.0_DBL + DSF%X2 ) * DSF%X1 ;
DSF%DFXI ( 6 , 2 ) =  ( 1.0_DBL - DSF%X1 * DSF%X1 ) * 0.5_DBL ; 
DSF%DFXI ( 7 , 1 ) = -( 1.0_DBL - DSF%X2 * DSF%X2 ) * 0.5_DBL ; 
DSF%DFXI ( 7 , 2 ) = -( 1.0_DBL - DSF%X1 ) * DSF%X2 ;
DSF%DFXI ( 8 , 1 ) = -( 1.0_DBL - DSF%X2 ) * DSF%X1 ;
DSF%DFXI ( 8 , 2 ) = -( 1.0_DBL - DSF%X1 * DSF%X1 ) * 0.5_DBL ; 


!Write(*,"('End Subroutine < SH_FUNC_3_20 >')")  ;
!Write(UnInf,"('End Subroutine < SH_FUNC_3_20 >')") ;
Return ;
End Subroutine DIF_SH_FUNC_2_8 ;


! ***************************************************************************************************************************************************
!                                                                                                                                                  **
! Last Update:  22 MAY 2011                                                                                                                        **
! Description: THIS Subroutines CALCULATES THE SHAPE FUNCTIONS OF A  { 2D FIRST  ORDER  4-NODE ELEMENT }. ELEMENT NUMBER : { 1 0R 3 }              **
!                                                                    ----------------------------------------------------------------              **
!                                                                                                                                                  **
! ***************************************************************************************************************************************************

Subroutine SH_FUNC_2_4 ( SF ) ;


Implicit None ;

Type ( SF_2_4 ) :: SF ;

! =========================== Subroutine CODE =======================================================================================================

SF%FN ( 1 ) = ( 1.0_DBL + SF%X1 ) * ( 1.0_DBL - SF%X2 ) * 0.25_DBL ; 
SF%FN ( 2 ) = ( 1.0_DBL + SF%X1 ) * ( 1.0_DBL + SF%X2 ) * 0.25_DBL ; 
SF%FN ( 3 ) = ( 1.0_DBL - SF%X1 ) * ( 1.0_DBL + SF%X2 ) * 0.25_DBL ; 
SF%FN ( 4 ) = ( 1.0_DBL - SF%X1 ) * ( 1.0_DBL - SF%X2 ) * 0.25_DBL ; 


!Write(*,"('End Subroutine < SH_FUNC_3_20 >')")  ;
!Write(UnInf,"('End Subroutine < SH_FUNC_3_20 >')") ;
Return ;
End Subroutine SH_FUNC_2_4 ;


! ***************************************************************************************************************************************************
!                                                                                                                                                  **
! Last Update:  22 MAY 2011                                                                                                                        **
! Description: THIS Subroutines CALCULATES THE DIFFERENTIALS OF SHAPE FUNCTIONS OF  A                                                              **
!              { 2D FIRTST ORDER  4-NODE ELEMENT }. ELEMENT NUMBER : { 1 0R 3 }.                                                                   **
!                                                                                                                                                  **
! ***************************************************************************************************************************************************

Subroutine DIF_SH_FUNC_2_4 ( DSF ) ;

Implicit None ;

Type ( DSF_2_4 ) :: DSF ;

! =========================== Subroutine CODE =======================================================================================================

DSF%DFXI ( 1 , 1 ) = + ( 1.0_DBL - DSF%X2 ) * 0.25_DBL ; 
DSF%DFXI ( 1 , 2 ) = - ( 1.0_DBL + DSF%X1 ) * 0.25_DBL ; 
DSF%DFXI ( 2 , 1 ) = + ( 1.0_DBL + DSF%X2 ) * 0.25_DBL ; 
DSF%DFXI ( 2 , 2 ) = + ( 1.0_DBL + DSF%X1 ) * 0.25_DBL ; 
DSF%DFXI ( 3 , 1 ) = - ( 1.0_DBL + DSF%X2 ) * 0.25_DBL ; 
DSF%DFXI ( 3 , 2 ) = + ( 1.0_DBL - DSF%X1 ) * 0.25_DBL ; 
DSF%DFXI ( 4 , 1 ) = - ( 1.0_DBL - DSF%X2 ) * 0.25_DBL ; 
DSF%DFXI ( 4 , 2 ) = - ( 1.0_DBL - DSF%X1 ) * 0.25_DBL ; 


!Write(*,"('End Subroutine < SH_FUNC_3_20 >')")  ;
!Write(UnInf,"('End Subroutine < SH_FUNC_3_20 >')") ;
Return ;
End Subroutine DIF_SH_FUNC_2_4 ;


! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

!  <><><><><><><><><><><><><><><><><><><><><><><><><><><> 1D 3N line element <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

! ***************************************************************************************************************************************************
!                                                                                                                                                  **
! Start:        20 July 2012                                                                                                                       **
! Last Update:  20 July 2012                                                                                                                       **
! Description: THIS function CALCULATES THE SHAPE FUNCTIONS OF A { 1D SECOND ORDER  3-NODE SEREndEPITI ELEMENT }. ELEMENT NUMBER : { }             **
!                                                                     ---------------------------------------------------------------------------  **
!                                                                                                                                                  **
! ***************************************************************************************************************************************************

Function ShapeFunc_1_3 ( X1 )   Result ( Funct ) ;

Implicit None ;

Real (Kind=DBL), Intent(In)    :: X1 ;
Real (Kind=DBL), Dimension(3)  :: Funct ;

! =========================== Function CODE =======================================================================================================

Funct ( 1 ) = -X1 * ( 1.0_DBL - X1 ) * 0.5_DBL ; 

Funct ( 2 ) = ( 1.0_DBL + X1 ) * ( 1.0_DBL - X1 ) ; 

Funct ( 3 ) = X1 * ( 1.0_DBL + X1 ) * 0.5_DBL ; 


End Function ShapeFunc_1_3 ;


! ***************************************************************************************************************************************************
!                                                                                                                                                  **
! Start:        20 July 2012                                                                                                                       **
! Last Update:  20 July 2012                                                                                                                       **
! Description: THIS Function CALCULATES THE DIFFERENTIALS OF SHAPE FUNCTIONS OF  A                                                                 **
!              { 1D SECOND ORDER  3-NODE SEREndEPITI ELEMENT }. ELEMENT NUMBER : {        }.                                                       **
!                                                                                                                                                  **
! ***************************************************************************************************************************************************

Function Dif_ShapeFunc_1_3 ( X1 )   Result ( DFXI ) ;


Implicit None ;

Real (Kind=DBL), Intent(In)    :: X1 ;
Real (Kind=DBL), Dimension(3,1):: DFXI ;

! =========================== Function CODE =======================================================================================================

DFXI ( 1 , 1 ) =  X1 - 0.5_Dbl ; 

DFXI ( 2 , 1 ) =  2.0_DBL * X1  ; 

DFXI ( 3 , 1 ) =  X1 + 0.5_DBL ; 

End Function Dif_ShapeFunc_1_3 ;


!  <><><><><><><><><><><><><><><><><><><><><><><><><><><> 2D 3N Linear Triangle  element <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

! ***************************************************************************************************************************************************
!                                                                                                                                                  **
! Start:        25 July 2012                                                                                                                       **
! Last Update:  25 July 2012                                                                                                                        **
! Description: THIS Function CALCULATES THE OF SHAPE FUNCTIONS OF  A                                                                               **
!              { 2D Fisrt ORDER  3-NODE ELEMENT }. ELEMENT NUMBER : { 9 OR 11 }.                                                                   **
!                                                                                                                                                  **
! ***************************************************************************************************************************************************

Function ShapeFunc_2_3 (r, s)  Result ( Funct ) ;

Implicit None ;

Real (Kind=DBL), Intent(In)    :: r, s ;
Real (Kind=DBL), Dimension(3)  :: Funct ;

! =========================== Function Code =========================================================================================================

Funct ( 1 ) = 1.0_DBL - r - s ; 
Funct ( 2 ) = r ; 
Funct ( 3 ) = s ; 

End Function ShapeFunc_2_3 ;


! ***************************************************************************************************************************************************
!                                                                                                                                                  **
! Start:        25 July 2012                                                                                                                       **
! Last Update:  25 July 2012                                                                                                                        **
! Description: This Function calculates the differentials of the shape functions of a                                                              **
!              { 2D FIRTST ORDER  3-NODE ELEMENT }. ELEMENT NUMBER : { 9 0R 11 }.                                                                  **
!                                                                                                                                                  **
! ***************************************************************************************************************************************************

Function DIF_ShapeFunc_2_3 ( r, s )   Result ( DFXI ) ;

Implicit None ;

Real (Kind=DBL), Intent(In)    :: r, s ;
Real (Kind=DBL), Dimension(3,2):: DFXI ;

! =========================== Subroutine CODE =======================================================================================================

DFXI ( 1 , 1 ) = -1.0_DBL ; 
DFXI ( 1 , 2 ) = -1.0_DBL ; 
DFXI ( 2 , 1 ) = +1.0_DBL ; 
DFXI ( 2 , 2 ) =  0.0_DBL ; 
DFXI ( 3 , 1 ) =  0.0_DBL ; 
DFXI ( 3 , 2 ) = +1.0_DBL ; 

End Function DIF_ShapeFunc_2_3 ;

!  <><><><><><><><><><><><><><><><><><><><><><><><><><><> 2D 6N Second order Triangle element <><><><><><><><><><><><><><><><><><><><><><><><><><><>

! ***************************************************************************************************************************************************
!                                                                                                                                                  **
! Last Update:  31 Jan 2012                                                                                                                        **
! Description: THIS Function calculates the shape functions of a                                                                                   **
!              { 2D SECOND ORDER  6-NODE SerEndipity ELEMENT }. ELEMENT NUMBER : { 10 OR 12 }.                                                     **
!                                                                                                                                                  **
! ***************************************************************************************************************************************************

Function ShapeFunc_2_6 (r, s)  Result ( Funct ) ;

Implicit None ;

Real (Kind=DBL), Intent(In)    :: r, s ;
Real (Kind=DBL), Dimension(6)  :: Funct ;


! =========================== Function Code =========================================================================================================

Funct ( 1 ) = ( 1.0_DBL - r - s ) * ( 1.0_DBL - 2.0_DBL * r - 2.0_DBL * s ) ; 
Funct ( 2 ) = r * ( 2.0_DBL * r - 1.0_DBL ) ; 
Funct ( 3 ) = s * ( 2.0_DBL * s - 1.0_DBL ) ; 
Funct ( 4 ) = 4.0_DBL * r * ( 1.0_DBL - r - s ) ;
Funct ( 5 ) = 4.0_DBL * r * s ;
Funct ( 6 ) = 4.0_DBL * s * ( 1.0_DBL - r - s ) ;

End Function ShapeFunc_2_6 ;


! ***************************************************************************************************************************************************
!                                                                                                                                                  **
! Last Update:  31 Jan 2011                                                                                                                        **
! Description: This Function calculates the differentials of the shape functions of a                                                              **
!              { 2D Second ORDER  3-NODE SerEndipity ELEMENT }. ELEMENT NUMBER : { 10 0R 12 }.                                                     **
!                                                                                                                                                  **
! ***************************************************************************************************************************************************

Function DIF_ShapeFunc_2_6 ( r, s )   Result ( DFXI ) ;

Implicit None ;

Real (Kind=DBL), Intent(In)    :: r, s ;
Real (Kind=DBL), Dimension(6,2):: DFXI ;

! =========================== Function CODE =======================================================================================================

DFXI ( 1 , 1 ) = - 3.0_DBL + 4.0_DBL * r + 4.0_DBL * s ; 
DFXI ( 1 , 2 ) = - 3.0_DBL + 4.0_DBL * r + 4.0_DBL * s ; 
DFXI ( 2 , 1 ) = + 4.0_DBL * r - 1.0_DBL ; 
DFXI ( 2 , 2 ) =   0.0_DBL ; 
DFXI ( 3 , 1 ) =   0.0_DBL ; 
DFXI ( 3 , 2 ) = + 4.0_DBL * s - 1.0_DBL ; 
DFXI ( 4 , 1 ) = + 4.0_DBL - 8.0_DBL * r - 4.0_DBL * s ; 
DFXI ( 4 , 2 ) = - 4.0_DBL * r ; 
DFXI ( 5 , 1 ) = + 4.0_DBL * s ; 
DFXI ( 5 , 2 ) = + 4.0_DBL * r ; 
DFXI ( 6 , 1 ) = - 4.0_DBL * s ; 
DFXI ( 6 , 2 ) = + 4.0_DBL - 4.0_DBL * r - 8.0_DBL * s ; 

End Function DIF_ShapeFunc_2_6 ;


!  <><><><><><><><><><><><><><><><><><><><><><><><><><><> 2D 7N Second order Triangle Spectral element <><><><><><><><><><><><><><><><><><><><><><><>

! ***************************************************************************************************************************************************
!                                                                                                                                                  **
! Last Update:  03 Jan 2013                                                                                                                        **
! Description: THIS Function calculates the shape functions of a                                                                                   **
!              { 2D SECOND ORDER  7-NODE SerEndipity ELEMENT }. ELEMENT NUMBER : { 36 OR 37 }.                                                     **
!                                                                                                                                                  **
! ***************************************************************************************************************************************************

Function ShapeFunc_2_7 (r, s)  Result ( Funct ) ;

Implicit None ;

Real (Kind=DBL), Intent(In)    :: r, s ;
Real (Kind=DBL), Dimension(7)  :: Funct ;


! =========================== Function Code =========================================================================================================

Funct ( 1 ) = 1.0_Dbl - 3.0_Dbl * (r+s) + 2.0_Dbl * (r * r +  s * s) + 7.0_Dbl * r * s- 3.0_Dbl * r* s * ( r + s ) ;
Funct ( 2 ) = r*(-1.0_Dbl+2.0_Dbl*r+3.0_Dbl*s-3.0_Dbl*s*(r+s)) ;
Funct ( 3 ) = s*(-1.0_Dbl+2.0_Dbl*s+3.0_Dbl*r-3.0_Dbl*r*(r+s)) ;
Funct ( 4 ) = 4.0_Dbl*r*(1.0_Dbl-r-4.0_Dbl*s+3.0_Dbl*s*(r+s) ) ;
Funct ( 5 ) = 4.0_Dbl*r*s*(-2.0_Dbl+3.0_Dbl*(r+s)) ;
Funct ( 6 ) = 4.0_Dbl*s*(1.0_Dbl-4.0_Dbl*r-s+3.0_Dbl*r*(r+s) ) ;
Funct ( 7 ) = 27.0_Dbl*r*s*(1.0_Dbl-r-s) ;

End Function ShapeFunc_2_7 ;


! ***************************************************************************************************************************************************
!                                                                                                                                                  **
! Last Update:  03 Jan 2013                                                                                                                        **
! Description: This Function calculates the derivatives of the shape functions of a                                                                **
!              { 2D Second ORDER  7-NODE SerEndipity ELEMENT }. ELEMENT NUMBER : { 36 0R 36 }.                                                     **
!                                                                                                                                                  **
! ***************************************************************************************************************************************************

Function DIF_ShapeFunc_2_7 ( r, s )   Result ( DFXI ) ;

Implicit None ;

Real (Kind=DBL), Intent(In)    :: r, s ;
Real (Kind=DBL), Dimension(7,2):: DFXI ;

! =========================== Function CODE =======================================================================================================

DFXI ( 1 , 1 ) = -3  +  4 * r  +  7 * s  -  6 * r * s  -  3 * s * s ;
DFXI ( 1 , 2 ) = -3  +  4 * s  +  7 * r  -  6 * r * s  -  3 * r * r ;

DFXI ( 2 , 1 ) = -1  +  4 * r  +  3 * s  -  6 * r * s  -  3 * s * s ;
DFXI ( 2 , 2 ) =  3 * r * ( 1  - r  - 2 * s ) ;

DFXI ( 3 , 1 ) =  3 * s * ( 1  - s  - 2 * r ) ;
DFXI ( 3 , 2 ) = -1  + 3 * r  + 4 * s  - 3 *r * r  - 6* r * s ;

DFXI ( 4 , 1 ) =  4 * ( 1  - 2 * r  - 4 * s  + 6 * r * s  + 3 * s * s ) ;
DFXI ( 4 , 2 ) =  4 * r * ( -4  + 3 * r  + 6 * s ) ;

DFXI ( 5 , 1 ) =  4 * s * ( -2  + 6 * r  + 3 * s ) ;
DFXI ( 5 , 2 ) =  4 * r * ( -2  + 6 * s  + 3 * r ) ;

DFXI ( 6 , 1 ) =  4 * s * ( -4  + 6 * r  + 3 * s ) ;
DFXI ( 6 , 2 ) =  4 * ( 1  - 4 * r  - 2 * s  + 3 * r *r  + 6 * r * s ) ;

DFXI ( 7 , 1 ) = 27 * s * ( 1  - 2 * r  - s ) ;
DFXI ( 7 , 2 ) = 27 * r * ( 1  - r - 2 * s ) ;

End Function DIF_ShapeFunc_2_7 ;

!  <><><><><><><><><><><><><><><><><><><><><><><><><><><> 2D 4N First order Quadrilateral element <><><><><><><><><><><><><><><><><><><><><><><><><>

! ***************************************************************************************************************************************************
!                                                                                                                                                  **
! Start:        25 July 2012                                                                                                                       **
! Last Update:  25 July 2012                                                                                                                        **
! Description: THIS function CALCULATES THE SHAPE FUNCTIONS OF A { 2D first ORDER  4-NODE SEREndEPITI ELEMENT }. ELEMENT NUMBER : { 2 0R 4 }      **
!                                                                     ---------------------------------------------------------------------------  **
!                                                                                                                                                  **
! ***************************************************************************************************************************************************

Function ShapeFunc_2_4 ( X1, X2 )   Result ( Funct ) ;

Implicit None ;

Real (Kind=DBL), Intent(In)    :: X1, X2 ;
Real (Kind=DBL), Dimension(4)  :: Funct ;

! =========================== Function CODE =======================================================================================================


Funct ( 1 ) = ( 1.0_DBL + X1 ) * ( 1.0_DBL - X2 ) * 0.25_DBL ; 
Funct ( 2 ) = ( 1.0_DBL + X1 ) * ( 1.0_DBL + X2 ) * 0.25_DBL ; 
Funct ( 3 ) = ( 1.0_DBL - X1 ) * ( 1.0_DBL + X2 ) * 0.25_DBL ; 
Funct ( 4 ) = ( 1.0_DBL - X1 ) * ( 1.0_DBL - X2 ) * 0.25_DBL ; 

End Function ShapeFunc_2_4 ;


! ***************************************************************************************************************************************************
!                                                                                                                                                  **
! Start:        25 July 2012                                                                                                                       **
! Last Update:  25 July 2012                                                                                                                       **
! Description: THIS Function CALCULATES THE DIFFERENTIALS OF SHAPE FUNCTIONS OF  A                                                                 **
!              { 2D First ORDER  4-NODE SEREndEPITI ELEMENT }. ELEMENT NUMBER : { 2 OR 4 }.                                                        **
!                                                                                                                                                  **
! ***************************************************************************************************************************************************

Function Dif_ShapeFunc_2_4 ( X1, X2 )   Result ( DFXI ) ;


Implicit None ;

Real (Kind=DBL), Intent(In)    :: X1, X2 ;
Real (Kind=DBL), Dimension(4,2):: DFXI ;

! =========================== Function CODE =======================================================================================================

DFXI ( 1 , 1 ) = + ( 1.0_DBL - X2 ) * 0.25_DBL ; 
DFXI ( 1 , 2 ) = - ( 1.0_DBL + X1 ) * 0.25_DBL ; 
DFXI ( 2 , 1 ) = + ( 1.0_DBL + X2 ) * 0.25_DBL ; 
DFXI ( 2 , 2 ) = + ( 1.0_DBL + X1 ) * 0.25_DBL ; 
DFXI ( 3 , 1 ) = - ( 1.0_DBL + X2 ) * 0.25_DBL ; 
DFXI ( 3 , 2 ) = + ( 1.0_DBL - X1 ) * 0.25_DBL ; 
DFXI ( 4 , 1 ) = - ( 1.0_DBL - X2 ) * 0.25_DBL ; 
DFXI ( 4 , 2 ) = - ( 1.0_DBL - X1 ) * 0.25_DBL ; 

End Function Dif_ShapeFunc_2_4 ;


!  <><><><><><><><><><><><><><><><><><><><><><><><><><><> 2D 8N Second order Quadrilateral element <><><><><><><><><><><><><><><><><><><><><><><><><>

! ***************************************************************************************************************************************************
!                                                                                                                                                  **
! Last Update:  03 Feb 2012                                                                                                                        **
! Description: THIS function CALCULATES THE SHAPE FUNCTIONS OF A { 2D SECOND ORDER  8-NODE SEREndEPITI ELEMENT }. ELEMENT NUMBER : { 2 0R 4 }      **
!                                                                     ---------------------------------------------------------------------------  **
!                                                                                                                                                  **
! ***************************************************************************************************************************************************

Function ShapeFunc_2_8 ( X1, X2 )   Result ( Funct ) ;

Implicit None ;

Real (Kind=DBL), Intent(In)    :: X1, X2 ;
Real (Kind=DBL), Dimension(8)  :: Funct ;

! =========================== Function CODE =======================================================================================================

Funct ( 1 ) = ( 1.0_DBL + X1 ) * ( 1.0_DBL - X2 ) * ( -1.0_DBL + X1 - X2 ) * 0.25_DBL ; 
Funct ( 2 ) = ( 1.0_DBL + X1 ) * ( 1.0_DBL + X2 ) * ( -1.0_DBL + X1 + X2 ) * 0.25_DBL ; 
Funct ( 3 ) = ( 1.0_DBL - X1 ) * ( 1.0_DBL + X2 ) * ( -1.0_DBL - X1 + X2 ) * 0.25_DBL ; 
Funct ( 4 ) = ( 1.0_DBL - X1 ) * ( 1.0_DBL - X2 ) * ( -1.0_DBL - X1 - X2 ) * 0.25_DBL ; 

Funct ( 5 ) = ( 1.0_DBL - X2 * X2 ) * ( 1.0_DBL + X1 ) * 0.5_DBL ; 
Funct ( 6 ) = ( 1.0_DBL - X1 * X1 ) * ( 1.0_DBL + X2 ) * 0.5_DBL ; 
Funct ( 7 ) = ( 1.0_DBL - X2 * X2 ) * ( 1.0_DBL - X1 ) * 0.5_DBL ; 
Funct ( 8 ) = ( 1.0_DBL - X1 * X1 ) * ( 1.0_DBL - X2 ) * 0.5_DBL ; 

End Function ShapeFunc_2_8 ;


! ***************************************************************************************************************************************************
!                                                                                                                                                  **
! Last Update:  03 Feb 2012                                                                                                                        **
! Description: THIS Function CALCULATES THE DIFFERENTIALS OF SHAPE FUNCTIONS OF  A                                                                 **
!              { 2D SECOND ORDER  8-NODE SEREndEPITI ELEMENT }. ELEMENT NUMBER : { 2 OR 4 }.                                                       **
!                                                                                                                                                  **
! ***************************************************************************************************************************************************


Function Dif_ShapeFunc_2_8 ( X1, X2 )   Result ( DFXI ) ;


Implicit None ;

Real (Kind=DBL), Intent(In)    :: X1, X2 ;
Real (Kind=DBL), Dimension(8,2):: DFXI ;

! =========================== Function CODE =======================================================================================================

DFXI ( 1 , 1 ) =  ( 1.0_DBL - X2 ) * (   2.0_DBL * X1 - X2 ) * 0.25_DBL ; 
DFXI ( 1 , 2 ) =  ( 1.0_DBL + X1 ) * ( - X1 + 2.0_DBL * X2 ) * 0.25_DBL ; 
DFXI ( 2 , 1 ) =  ( 1.0_DBL + X2 ) * (   2.0_DBL * X1 + X2 ) * 0.25_DBL ; 
DFXI ( 2 , 2 ) =  ( 1.0_DBL + X1 ) * (   X1 + 2.0_DBL * X2 ) * 0.25_DBL ; 
DFXI ( 3 , 1 ) =  ( 1.0_DBL + X2 ) * (   2.0_DBL * X1 - X2 ) * 0.25_DBL ; 
DFXI ( 3 , 2 ) =  ( 1.0_DBL - X1 ) * ( - X1 + 2.0_DBL * X2 ) * 0.25_DBL ; 
DFXI ( 4 , 1 ) =  ( 1.0_DBL - X2 ) * (   2.0_DBL * X1 + X2 ) * 0.25_DBL ; 
DFXI ( 4 , 2 ) =  ( 1.0_DBL - X1 ) * (   X1 + 2.0_DBL * X2 ) * 0.25_DBL ; 

DFXI ( 5 , 1 ) =  ( 1.0_DBL - X2 * X2 ) * 0.5_DBL ; 
DFXI ( 5 , 2 ) = -( 1.0_DBL + X1 ) * X2 ;
DFXI ( 6 , 1 ) = -( 1.0_DBL + X2 ) * X1 ;
DFXI ( 6 , 2 ) =  ( 1.0_DBL - X1 * X1 ) * 0.5_DBL ; 
DFXI ( 7 , 1 ) = -( 1.0_DBL - X2 * X2 ) * 0.5_DBL ; 
DFXI ( 7 , 2 ) = -( 1.0_DBL - X1 ) * X2 ;
DFXI ( 8 , 1 ) = -( 1.0_DBL - X2 ) * X1 ;
DFXI ( 8 , 2 ) = -( 1.0_DBL - X1 * X1 ) * 0.5_DBL ; 


End Function Dif_ShapeFunc_2_8 ;


!  <><><><><><><><><><><><><><><><><><><><><><><><><><><> 2D 9N Second order Quadrilateral element <><><><><><><><><><><><><><><><><><><><><><><><><>

! ***************************************************************************************************************************************************
!                                                                                                                                                  **
! Last Update:  03 Aug 2012                                                                                                                        **
! Description: THIS function CALCULATES THE SHAPE FUNCTIONS OF A { 2D SECOND ORDER  9-NODE SEREndEPITI ELEMENT }. ELEMENT NUMBER : { 2 0R 4 }      **
!                                                                     ---------------------------------------------------------------------------  **
!                                                                                                                                                  **
! ***************************************************************************************************************************************************

Function ShapeFunc_2_9 ( X1, X2 )   Result ( Funct ) ;

Implicit None ;

Real (Kind=DBL), Intent(In)    :: X1, X2 ;
Real (Kind=DBL), Dimension(9)  :: Funct ;

! =========================== Function CODE =======================================================================================================

Funct (1) = 0.25_Dbl *  X1 * -X2 * (1.0_Dbl+X1) * (1.0_Dbl-X2)
Funct (2) = 0.25_Dbl *  X1 *  X2 * (1.0_Dbl+X1) * (1.0_Dbl+X2)
Funct (3) = 0.25_Dbl * -X1 *  X2 * (1.0_Dbl-X1) * (1.0_Dbl+X2)
Funct (4) = 0.25_Dbl * -X1 * -X2 * (1.0_Dbl-X1) * (1.0_Dbl-X2)
Funct (5) = 0.50_Dbl *  X1 *       (1.0_Dbl+X1) *              (1.0_Dbl-X2*X2)
Funct (6) = 0.50_Dbl *        X2 *              (1.0_Dbl+X2) * (1.0_Dbl-X1*X1)
Funct (7) = 0.50_Dbl * -X1 *       (1.0_Dbl-X1) *              (1.0_Dbl-X2*X2)
Funct (8) = 0.50_Dbl *       -X2 *              (1.0_Dbl-X2) * (1.0_Dbl-X1*X1) 
Funct (9) = (1.0_Dbl-X1*X1) * (1.0_Dbl-X2*X2)

End Function ShapeFunc_2_9 ;

! ***************************************************************************************************************************************************
!                                                                                                                                                  **
! Last Update:  03 Aug 2012                                                                                                                        **
! Description: THIS Function CALCULATES THE DIFFERENTIALS OF SHAPE FUNCTIONS OF  A                                                                 **
!              { 2D SECOND ORDER  8-NODE SEREndEPITI ELEMENT }. ELEMENT NUMBER : { 2 OR 4 }.                                                       **
!                                                                                                                                                  **
! ***************************************************************************************************************************************************


Function Dif_ShapeFunc_2_9 ( X1, X2 )   Result ( DFXI ) ;


Implicit None ;

Real (Kind=DBL), Intent(In)    :: X1, X2 ;
Real (Kind=DBL), Dimension(9,2):: DFXI ;

! =========================== Function CODE =======================================================================================================

DFXI (1,1) = 0.25_Dbl * -X2 * (1.0_Dbl-X2) * ( 1.0_Dbl + 2.0_Dbl * X1)  
DFXI (2,1) = 0.25_Dbl *  X2 * (1.0_Dbl+X2) * ( 1.0_Dbl + 2.0_Dbl * X1)
DFXI (3,1) = 0.25_Dbl *  X2 * (1.0_Dbl+X2) * (-1.0_Dbl + 2.0_Dbl * X1)
DFXI (4,1) = 0.25_Dbl * -X2 * (1.0_Dbl-X2) * (-1.0_Dbl + 2.0_Dbl * X1)

DFXI (5,1) = 0.50_Dbl * (1.0_Dbl-X2*X2) * ( 1.0_Dbl + 2.0_Dbl * X1)
DFXI (7,1) = 0.50_Dbl * (1.0_Dbl-X2*X2) * (-1.0_Dbl + 2.0_Dbl * X1)

DFXI (6,1) =  X2 * (1.0_Dbl+X2) * -X1 
DFXI (8,1) = -X2 * (1.0_Dbl-X2) * -X1

DFXI (9,1) = -2.0_Dbl * X1 * (1.0_Dbl-X2*X2) 

DFXI (1,2) = 0.25_Dbl *  X1 * (1.0_Dbl+X1) * (-1.0_Dbl + 2.0_Dbl * X2)
DFXI (2,2) = 0.25_Dbl *  X1 * (1.0_Dbl+X1) * ( 1.0_Dbl + 2.0_Dbl * X2)
DFXI (3,2) = 0.25_Dbl * -X1 * (1.0_Dbl-X1) * ( 1.0_Dbl + 2.0_Dbl * X2)
DFXI (4,2) = 0.25_Dbl * -X1 * (1.0_Dbl-X1) * (-1.0_Dbl + 2.0_Dbl * X2)

DFXI (5,2) =  X1 * (1.0_Dbl+X1) * -X2
DFXI (7,2) = -X1 * (1.0_Dbl-X1) * -X2
DFXI (6,2) = 0.50_Dbl * (1.0_Dbl-X1*X1) * ( 1.0_Dbl + 2.0_Dbl * X2)
DFXI (8,2) = 0.50_Dbl * (1.0_Dbl-X1*X1) * (-1.0_Dbl + 2.0_Dbl * X2)

DFXI (9,2) = -2.0_Dbl * X2 * (1.0_Dbl-X1*X1)

End Function Dif_ShapeFunc_2_9 ;


!  <><><><><><><><><><><><><><><><><><><><><><><><><><><> 3D 4N First order tetralhedral element <><><><><><><><><><><><><><><><><><><><><><><><><><>

! ***************************************************************************************************************************************************
!                                                                                                                                                  **
! Start:        01 June 2012                                                                                                                       **
! Last Update:  01 June 2012                                                                                                                       **
! Description: This function computes the shape functions of a { 3d first order  4-node isoparametric element }. element number : { 13 }           **
!                                                                                                                                                  **
!                                                                                                                                                  **
! ***************************************************************************************************************************************************

Function ShapeFunc_3_4 ( r, s, t )   Result ( Funct ) ;

Implicit None ;

Real (Kind=DBL), Intent(In)     :: r, s, t ;
Real (Kind=DBL), Dimension(4)  :: Funct ;

! =========================== Function CODE =======================================================================================================

Funct ( 1  ) = 1.0_DBL - r - s - t ;
Funct ( 2  ) = r ;
Funct ( 3  ) = s ;
Funct ( 4  ) = t ;

End Function ShapeFunc_3_4 ;

! ***************************************************************************************************************************************************
!                                                                                                                                                  **
! Start:        01 June 2012                                                                                                                       **
! Last Update:  01 June 2012                                                                                                                       **
! Description: This Function calculates the differentials of the shape functions of a                                                              **
!              { 3D first order  4-node isoparametric element }. ELEMENT NUMBER : { 13 }.                                                          **
!                                                                                                                                                  **
! ***************************************************************************************************************************************************

Function DIF_ShapeFunc_3_4 ( r, s, t )   Result ( DFXI ) ;

Implicit None ;

Real (Kind=DBL), Intent(In)    :: r, s, t ;
Real (Kind=DBL), Dimension(4,3):: DFXI ;

! =========================== Function CODE =======================================================================================================

DFXI ( 1 , 1 ) = - 1.0_DBL ;
DFXI ( 1 , 2 ) = - 1.0_DBL ;
DFXI ( 1 , 3 ) = - 1.0_DBL ;

DFXI ( 2 , 1 ) = + 1.0_DBL ;
DFXI ( 2 , 2 ) =   0.0_DBL ;
DFXI ( 2 , 3 ) =   0.0_DBL ;

DFXI ( 3 , 1 ) =   0.0_DBL ;
DFXI ( 3 , 2 ) = + 1.0_DBL ;
DFXI ( 3 , 3 ) =   0.0_DBL ;

DFXI ( 4 , 1 ) =   0.0_DBL ;
DFXI ( 4 , 2 ) =   0.0_DBL ;
DFXI ( 4 , 3 ) = + 1.0_DBL ;

End Function DIF_ShapeFunc_3_4 ;

!  <><><><><><><><><><><><><><><><><><><><><><><><><><><> 3D 10N Second order tetralhedral element <><><><><><><><><><><><><><><><><><><><><><><><><>

! ***************************************************************************************************************************************************
!                                                                                                                                                  **
! Start:        01 June 2012                                                                                                                       **
! Last Update:  01 June 2012                                                                                                                       **
! Description: This function computes the shape functions of a { 3d second order  10-node isoparametric element }. element number : { 14 }         **
!                                                                                                                                                  **
!                                                                                                                                                  **
! ***************************************************************************************************************************************************

Function ShapeFunc_3_10 ( r, s, t )   Result ( Funct ) ;

Implicit None ;

Real (Kind=DBL), Intent(In)     :: r, s, t ;
Real (Kind=DBL), Dimension(10)  :: Funct ;

! =========================== Function CODE =======================================================================================================

Funct ( 1  ) = ( 1.0_DBL - r - s - t ) * ( 1.0_DBL - 2.0_Dbl * r - 2.0_Dbl * s - 2.0_Dbl * t );
Funct ( 2  ) = r * ( 2.0_Dbl * r - 1 ) ;
Funct ( 3  ) = s * ( 2.0_Dbl * s - 1 ) ;
Funct ( 4  ) = t * ( 2.0_Dbl * t - 1 ) ;

Funct ( 5  ) = 4.0_Dbl * r * ( 1.0_DBL - r - s - t ) ;
Funct ( 6  ) = 4.0_Dbl * s * ( 1.0_DBL - r - s - t ) ;
Funct ( 7  ) = 4.0_Dbl * t * ( 1.0_DBL - r - s - t ) ;

Funct ( 8  ) = 4.0_Dbl * r * s ;
Funct ( 9  ) = 4.0_Dbl * s * t ;
Funct ( 10 ) = 4.0_Dbl * t * r ;

End Function ShapeFunc_3_10 ;


! ***************************************************************************************************************************************************
!                                                                                                                                                  **
! Start:        01 June 2012                                                                                                                       **
! Last Update:  01 June 2012                                                                                                                       **
! Description: This Function calculates the differentials of the shape functions of a                                                              **
!              { 3D first order  4-node isoparametric element }. ELEMENT NUMBER : { 13 }.                                                          **
!                                                                                                                                                  **
! ***************************************************************************************************************************************************

Function DIF_ShapeFunc_3_10 ( r, s, t )   Result ( DFXI ) ;

Implicit None ;

Real (Kind=DBL), Intent(In)    :: r, s, t ;
Real (Kind=DBL), Dimension(10,3):: DFXI ;

! =========================== Function CODE =======================================================================================================

DFXI ( 1 , 1 ) = - 3.0_DBL + 4.0_DBL * r + 4.0_DBL * s + 4.0_DBL * t ; 
DFXI ( 1 , 2 ) = - 3.0_DBL + 4.0_DBL * r + 4.0_DBL * s + 4.0_DBL * t ; 
DFXI ( 1 , 3 ) = - 3.0_DBL + 4.0_DBL * r + 4.0_DBL * s + 4.0_DBL * t ; 

DFXI ( 2 , 1 ) = + 4.0_DBL * r - 1.0_DBL ; 
DFXI ( 2 , 2 ) =   0.0_DBL ; 
DFXI ( 2 , 3 ) =   0.0_DBL ; 

DFXI ( 3 , 1 ) =   0.0_DBL ; 
DFXI ( 3 , 2 ) = + 4.0_DBL * s - 1.0_DBL ; 
DFXI ( 3 , 3 ) =   0.0_DBL ; 

DFXI ( 4 , 1 ) =   0.0_DBL ; 
DFXI ( 4 , 2 ) =   0.0_DBL ; 
DFXI ( 4 , 3 ) = + 4.0_DBL * t - 1.0_DBL ; 

DFXI ( 5 , 1 ) = + 4.0_DBL - 8.0_DBL * r - 4.0_DBL * s - 4.0_DBL * t ; 
DFXI ( 5 , 2 ) = - 4.0_DBL * r ; 
DFXI ( 5 , 3 ) = - 4.0_DBL * r ; 

DFXI ( 6 , 1 ) = - 4.0_DBL * s ; 
DFXI ( 6 , 2 ) = + 4.0_DBL - 4.0_DBL * r - 8.0_DBL * s - 4.0_DBL * t ; 
DFXI ( 6 , 3 ) = - 4.0_DBL * s ; 

DFXI ( 7 , 1 ) = - 4.0_DBL * t ; 
DFXI ( 7 , 2 ) = - 4.0_DBL * t ; 
DFXI ( 7 , 3 ) = + 4.0_DBL - 4.0_DBL * r - 4.0_DBL * s - 8.0_DBL * t ; 

DFXI ( 8 , 1 ) = + 4.0_DBL * s ; 
DFXI ( 8 , 2 ) = + 4.0_DBL * r ; 
DFXI ( 8 , 3 ) =   0.0_DBL ; 

DFXI ( 9 , 1 ) =   0.0_DBL ; 
DFXI ( 9 , 2 ) = + 4.0_DBL * t ; 
DFXI ( 9 , 3 ) = + 4.0_DBL * s ; 

DFXI ( 10, 1 ) = + 4.0_DBL * t ; 
DFXI ( 10, 2 ) =   0.0_DBL ; 
DFXI ( 10, 3 ) = + 4.0_DBL * r ; 

End Function DIF_ShapeFunc_3_10 ;

!  <><><><><><><><><><><><><><><><><><><><><><><><><><><> 3D 8N Second order Hexahedral element  <><><><><><><><><><><><><><><><><><><><><><><><><><>

! ***************************************************************************************************************************************************
!                                                                                                                                                  **
! Last Update:  14 May 2013                                                                                                                        **
! Description: This function computes the shape functions of a { 3d second order  27-node Lagrange element }. element number : { ?? }              **
!                                                                                                                                                  **
!                                                                                                                                                  **
! ***************************************************************************************************************************************************

Function ShapeFunc_3_8 ( X1, X2, X3 )   Result ( Funct ) ;

Implicit None ;

Real (Kind=DBL), Intent(In)     :: X1, X2, X3 ;
Real (Kind=DBL), Dimension(8)  :: Funct ;

! =========================== Function CODE =======================================================================================================

Funct ( 1 ) = ( 1.0_DBL + X1 ) * ( 1.0_DBL - X2 ) * ( 1.0_DBL - X3 ) * 0.125_DBL ;
Funct ( 2 ) = ( 1.0_DBL + X1 ) * ( 1.0_DBL + X2 ) * ( 1.0_DBL - X3 ) * 0.125_DBL ;
Funct ( 3 ) = ( 1.0_DBL - X1 ) * ( 1.0_DBL + X2 ) * ( 1.0_DBL - X3 ) * 0.125_DBL ;
Funct ( 4 ) = ( 1.0_DBL - X1 ) * ( 1.0_DBL - X2 ) * ( 1.0_DBL - X3 ) * 0.125_DBL ;
Funct ( 5 ) = ( 1.0_DBL + X1 ) * ( 1.0_DBL - X2 ) * ( 1.0_DBL + X3 ) * 0.125_DBL ;
Funct ( 6 ) = ( 1.0_DBL + X1 ) * ( 1.0_DBL + X2 ) * ( 1.0_DBL + X3 ) * 0.125_DBL ;
Funct ( 7 ) = ( 1.0_DBL - X1 ) * ( 1.0_DBL + X2 ) * ( 1.0_DBL + X3 ) * 0.125_DBL ;
Funct ( 8 ) = ( 1.0_DBL - X1 ) * ( 1.0_DBL - X2 ) * ( 1.0_DBL + X3 ) * 0.125_DBL ;

End Function ShapeFunc_3_8 ;


! ***************************************************************************************************************************************************
!                                                                                                                                                  **
! Last Update:  14 May 2013                                                                                                                        **
! Description: This function computes the differentials of shape functions of  a                                                                   **
!              { 3d second order  8-node Lagrangian element }. element number : { 6 }.                                                            **
!                                                                                                                                                  **
! ***************************************************************************************************************************************************

Function Dif_ShapeFunc_3_8 ( X1, X2, X3 )   Result ( DFXI ) ;


Implicit None ;

Real (Kind=DBL), Intent(In)     :: X1, X2, X3 ;
Real (Kind=DBL), Dimension(8,3) :: DFXI ;

! =========================== Function CODE =======================================================================================================

DFXI ( 1 , 1 ) = + ( 1.0_DBL - X2 ) * ( 1.0_DBL - X3 ) * 0.125_DBL ;
DFXI ( 1 , 2 ) = - ( 1.0_DBL + X1 ) * ( 1.0_DBL - X3 ) * 0.125_DBL ;
DFXI ( 1 , 3 ) = - ( 1.0_DBL + X1 ) * ( 1.0_DBL - X2 ) * 0.125_DBL ;

DFXI ( 2 , 1 ) = + ( 1.0_DBL + X2 ) * ( 1.0_DBL - X3 ) * 0.125_DBL ;
DFXI ( 2 , 2 ) = + ( 1.0_DBL + X1 ) * ( 1.0_DBL - X3 ) * 0.125_DBL ;
DFXI ( 2 , 3 ) = - ( 1.0_DBL + X1 ) * ( 1.0_DBL + X2 ) * 0.125_DBL ;

DFXI ( 3 , 1 ) = - ( 1.0_DBL + X2 ) * ( 1.0_DBL - X3 ) * 0.125_DBL ;
DFXI ( 3 , 2 ) = + ( 1.0_DBL - X1 ) * ( 1.0_DBL - X3 ) * 0.125_DBL ;
DFXI ( 3 , 3 ) = - ( 1.0_DBL - X1 ) * ( 1.0_DBL + X2 ) * 0.125_DBL ;

DFXI ( 4 , 1 ) = - ( 1.0_DBL - X2 ) * ( 1.0_DBL - X3 ) * 0.125_DBL ;
DFXI ( 4 , 2 ) = - ( 1.0_DBL - X1 ) * ( 1.0_DBL - X3 ) * 0.125_DBL ;
DFXI ( 4 , 3 ) = - ( 1.0_DBL - X1 ) * ( 1.0_DBL - X2 ) * 0.125_DBL ;

DFXI ( 5 , 1 ) = + ( 1.0_DBL - X2 ) * ( 1.0_DBL + X3 ) * 0.125_DBL ;
DFXI ( 5 , 2 ) = - ( 1.0_DBL + X1 ) * ( 1.0_DBL + X3 ) * 0.125_DBL ;
DFXI ( 5 , 3 ) = + ( 1.0_DBL + X1 ) * ( 1.0_DBL - X2 ) * 0.125_DBL ;

DFXI ( 6 , 1 ) = + ( 1.0_DBL + X2 ) * ( 1.0_DBL + X3 ) * 0.125_DBL ;
DFXI ( 6 , 2 ) = + ( 1.0_DBL + X1 ) * ( 1.0_DBL + X3 ) * 0.125_DBL ;
DFXI ( 6 , 3 ) = + ( 1.0_DBL + X1 ) * ( 1.0_DBL + X2 ) * 0.125_DBL ;

DFXI ( 7 , 1 ) = - ( 1.0_DBL + X2 ) * ( 1.0_DBL + X3 ) * 0.125_DBL ;
DFXI ( 7 , 2 ) = + ( 1.0_DBL - X1 ) * ( 1.0_DBL + X3 ) * 0.125_DBL ;
DFXI ( 7 , 3 ) = + ( 1.0_DBL - X1 ) * ( 1.0_DBL + X2 ) * 0.125_DBL ;

DFXI ( 8 , 1 ) = - ( 1.0_DBL - X2 ) * ( 1.0_DBL + X3 ) * 0.125_DBL ;
DFXI ( 8 , 2 ) = - ( 1.0_DBL - X1 ) * ( 1.0_DBL + X3 ) * 0.125_DBL ;
DFXI ( 8 , 3 ) = + ( 1.0_DBL - X1 ) * ( 1.0_DBL - X2 ) * 0.125_DBL ;

End Function Dif_ShapeFunc_3_8 ;

!  <><><><><><><><><><><><><><><><><><><><><><><><><><><> 3D 20N Second order Hexahedral element <><><><><><><><><><><><><><><><><><><><><><><><><><>

! ***************************************************************************************************************************************************
!                                                                                                                                                  **
! Last Update:  14 May 2012                                                                                                                        **
! Description: This function computes the shape functions of a { 3d second order  20-node serendepiti element }. element number : { 6 }            **
!                                                                                                                                                  **
!                                                                                                                                                  **
! ***************************************************************************************************************************************************

Function ShapeFunc_3_20 ( X1, X2, X3 )   Result ( Funct ) ;

Implicit None ;

Real (Kind=DBL), Intent(In)     :: X1, X2, X3 ;
Real (Kind=DBL), Dimension(20)  :: Funct ;

! =========================== Function CODE =======================================================================================================

Funct ( 1  ) = ( 1.0_DBL + X1 ) * ( 1.0_DBL - X2 ) * ( 1.0_DBL - X3 ) * ( -2.0_DBL + X1 - X2 - X3 ) * 0.125_DBL ;
Funct ( 2  ) = ( 1.0_DBL + X1 ) * ( 1.0_DBL + X2 ) * ( 1.0_DBL - X3 ) * ( -2.0_DBL + X1 + X2 - X3 ) * 0.125_DBL ;
Funct ( 3  ) = ( 1.0_DBL - X1 ) * ( 1.0_DBL + X2 ) * ( 1.0_DBL - X3 ) * ( -2.0_DBL - X1 + X2 - X3 ) * 0.125_DBL ;
Funct ( 4  ) = ( 1.0_DBL - X1 ) * ( 1.0_DBL - X2 ) * ( 1.0_DBL - X3 ) * ( -2.0_DBL - X1 - X2 - X3 ) * 0.125_DBL ;
Funct ( 5  ) = ( 1.0_DBL + X1 ) * ( 1.0_DBL - X2 ) * ( 1.0_DBL + X3 ) * ( -2.0_DBL + X1 - X2 + X3 ) * 0.125_DBL ;
Funct ( 6  ) = ( 1.0_DBL + X1 ) * ( 1.0_DBL + X2 ) * ( 1.0_DBL + X3 ) * ( -2.0_DBL + X1 + X2 + X3 ) * 0.125_DBL ;
Funct ( 7  ) = ( 1.0_DBL - X1 ) * ( 1.0_DBL + X2 ) * ( 1.0_DBL + X3 ) * ( -2.0_DBL - X1 + X2 + X3 ) * 0.125_DBL ;
Funct ( 8  ) = ( 1.0_DBL - X1 ) * ( 1.0_DBL - X2 ) * ( 1.0_DBL + X3 ) * ( -2.0_DBL - X1 - X2 + X3 ) * 0.125_DBL ;

Funct ( 9  ) = ( 1.0_DBL - X2 * X2 ) * ( 1.0_DBL + X1 ) * ( 1.0_DBL - X3 ) * 0.25_DBL ;
Funct ( 10 ) = ( 1.0_DBL - X1 * X1 ) * ( 1.0_DBL + X2 ) * ( 1.0_DBL - X3 ) * 0.25_DBL ;
Funct ( 11 ) = ( 1.0_DBL - X2 * X2 ) * ( 1.0_DBL - X1 ) * ( 1.0_DBL - X3 ) * 0.25_DBL ;
Funct ( 12 ) = ( 1.0_DBL - X1 * X1 ) * ( 1.0_DBL - X2 ) * ( 1.0_DBL - X3 ) * 0.25_DBL ;
Funct ( 13 ) = ( 1.0_DBL - X3 * X3 ) * ( 1.0_DBL + X1 ) * ( 1.0_DBL - X2 ) * 0.25_DBL ;
Funct ( 14 ) = ( 1.0_DBL - X3 * X3 ) * ( 1.0_DBL + X1 ) * ( 1.0_DBL + X2 ) * 0.25_DBL ;
Funct ( 15 ) = ( 1.0_DBL - X3 * X3 ) * ( 1.0_DBL - X1 ) * ( 1.0_DBL + X2 ) * 0.25_DBL ;
Funct ( 16 ) = ( 1.0_DBL - X3 * X3 ) * ( 1.0_DBL - X1 ) * ( 1.0_DBL - X2 ) * 0.25_DBL ;
Funct ( 17 ) = ( 1.0_DBL - X2 * X2 ) * ( 1.0_DBL + X1 ) * ( 1.0_DBL + X3 ) * 0.25_DBL ;
Funct ( 18 ) = ( 1.0_DBL - X1 * X1 ) * ( 1.0_DBL + X2 ) * ( 1.0_DBL + X3 ) * 0.25_DBL ;
Funct ( 19 ) = ( 1.0_DBL - X2 * X2 ) * ( 1.0_DBL - X1 ) * ( 1.0_DBL + X3 ) * 0.25_DBL ;
Funct ( 20 ) = ( 1.0_DBL - X1 * X1 ) * ( 1.0_DBL - X2 ) * ( 1.0_DBL + X3 ) * 0.25_DBL ;

End Function ShapeFunc_3_20 ;

! ***************************************************************************************************************************************************
!                                                                                                                                                  **
! Last Update:  14 May 2012                                                                                                                        **
! Description: This function computes the differentials of shape functions of  a                                                                   **
!              { 3d second order  20-node serendepiti element }. element number : { 6 }.                                                           **
!                                                                                                                                                  **
! ***************************************************************************************************************************************************


Function Dif_ShapeFunc_3_20 ( X1, X2, X3 )   Result ( DFXI ) ;


Implicit None ;

Real (Kind=DBL), Intent(In)     :: X1, X2, X3 ;
Real (Kind=DBL), Dimension(20,3) :: DFXI ;

! =========================== Function CODE =======================================================================================================

DFXI ( 1  , 1 ) = ( 1.0_DBL - X2 ) * ( 1.0_DBL - X3 ) * ( -1.0_DBL + 2 * X1 - X2 - X3 ) * 0.125_DBL ; ! dfn / dx
DFXI ( 1  , 2 ) = ( 1.0_DBL + X1 ) * ( 1.0_DBL - X3 ) * (  1.0_DBL - X1 + 2 * X2 + X3 ) * 0.125_DBL ; ! dfn / dy
DFXI ( 1  , 3 ) = ( 1.0_DBL + X1 ) * ( 1.0_DBL - X2 ) * (  1.0_DBL - X1 + X2 + 2 * X3 ) * 0.125_DBL ; ! dfn / dz

DFXI ( 2  , 1 ) = ( 1.0_DBL + X2 ) * ( 1.0_DBL - X3 ) * ( -1.0_DBL + 2 * X1 + X2 - X3 ) * 0.125_DBL ;
DFXI ( 2  , 2 ) = ( 1.0_DBL + X1 ) * ( 1.0_DBL - X3 ) * ( -1.0_DBL + X1 + 2 * X2 - X3 ) * 0.125_DBL ;
DFXI ( 2  , 3 ) = ( 1.0_DBL + X1 ) * ( 1.0_DBL + X2 ) * (  1.0_DBL - X1 - X2 + 2 * X3 ) * 0.125_DBL ;

DFXI ( 3  , 1 ) = ( 1.0_DBL + X2 ) * ( 1.0_DBL - X3 ) * (  1.0_DBL + 2 * X1 - X2 + X3 ) * 0.125_DBL ;
DFXI ( 3  , 2 ) = ( 1.0_DBL - X1 ) * ( 1.0_DBL - X3 ) * ( -1.0_DBL - X1 + 2 * X2 - X3 ) * 0.125_DBL ;
DFXI ( 3  , 3 ) = ( 1.0_DBL - X1 ) * ( 1.0_DBL + X2 ) * (  1.0_DBL + X1 - X2 + 2 * X3 ) * 0.125_DBL ;

DFXI ( 4  , 1 ) = ( 1.0_DBL - X2 ) * ( 1.0_DBL - X3 ) * (  1.0_DBL + 2 * X1 + X2 + X3 ) * 0.125_DBL ;
DFXI ( 4  , 2 ) = ( 1.0_DBL - X1 ) * ( 1.0_DBL - X3 ) * (  1.0_DBL + X1 + 2 * X2 + X3 ) * 0.125_DBL ;
DFXI ( 4  , 3 ) = ( 1.0_DBL - X1 ) * ( 1.0_DBL - X2 ) * (  1.0_DBL + X1 + X2 + 2 * X3 ) * 0.125_DBL ;

DFXI ( 5  , 1 ) = ( 1.0_DBL - X2 ) * ( 1.0_DBL + X3 ) * ( -1.0_DBL + 2 * X1 - X2 + X3 ) * 0.125_DBL ;
DFXI ( 5  , 2 ) = ( 1.0_DBL + X1 ) * ( 1.0_DBL + X3 ) * (  1.0_DBL - X1 + 2 * X2 - X3 ) * 0.125_DBL ;
DFXI ( 5  , 3 ) = ( 1.0_DBL + X1 ) * ( 1.0_DBL - X2 ) * ( -1.0_DBL + X1 - X2 + 2 * X3 ) * 0.125_DBL ;

DFXI ( 6  , 1 ) = ( 1.0_DBL + X2 ) * ( 1.0_DBL + X3 ) * ( -1.0_DBL + 2 * X1 + X2 + X3 ) * 0.125_DBL ;
DFXI ( 6  , 2 ) = ( 1.0_DBL + X1 ) * ( 1.0_DBL + X3 ) * ( -1.0_DBL + X1 + 2 * X2 + X3 ) * 0.125_DBL ;
DFXI ( 6  , 3 ) = ( 1.0_DBL + X1 ) * ( 1.0_DBL + X2 ) * ( -1.0_DBL + X1 + X2 + 2 * X3 ) * 0.125_DBL ;

DFXI ( 7  , 1 ) = ( 1.0_DBL + X2 ) * ( 1.0_DBL + X3 ) * (  1.0_DBL + 2 * X1 - X2 - X3 ) * 0.125_DBL ;
DFXI ( 7  , 2 ) = ( 1.0_DBL - X1 ) * ( 1.0_DBL + X3 ) * ( -1.0_DBL - X1 + 2 * X2 + X3 ) * 0.125_DBL ;
DFXI ( 7  , 3 ) = ( 1.0_DBL - X1 ) * ( 1.0_DBL + X2 ) * ( -1.0_DBL - X1 + X2 + 2 * X3 ) * 0.125_DBL ;

DFXI ( 8  , 1 ) = ( 1.0_DBL - X2 ) * ( 1.0_DBL + X3 ) * (  1.0_DBL + 2 * X1 + X2 - X3 ) * 0.125_DBL ;
DFXI ( 8  , 2 ) = ( 1.0_DBL - X1 ) * ( 1.0_DBL + X3 ) * (  1.0_DBL + X1 + 2 * X2 - X3 ) * 0.125_DBL ;
DFXI ( 8  , 3 ) = ( 1.0_DBL - X1 ) * ( 1.0_DBL - X2 ) * ( -1.0_DBL - X1 - X2 + 2 * X3 ) * 0.125_DBL ;


DFXI ( 9  , 1 ) =  ( 1.0_DBL -  X2 * X2 ) * ( 1.0_DBL - X3 ) * 0.25_DBL ; 
DFXI ( 9  , 2 ) = - X2 * ( 1.0_DBL + X1 ) * ( 1.0_DBL - X3 ) * 0.5_DBL  ; 
DFXI ( 9  , 3 ) = -( 1.0_DBL -  X2 * X2 ) * ( 1.0_DBL + X1 ) * 0.25_DBL ; 

DFXI ( 10 , 1 ) = - X1 * ( 1.0_DBL + X2 ) * ( 1.0_DBL - X3 ) * 0.5_DBL  ; 
DFXI ( 10 , 2 ) =  ( 1.0_DBL -  X1 * X1 ) * ( 1.0_DBL - X3 ) * 0.25_DBL ; 
DFXI ( 10 , 3 ) = -( 1.0_DBL -  X1 * X1 ) * ( 1.0_DBL + X2 ) * 0.25_DBL ; 

DFXI ( 11 , 1 ) = -( 1.0_DBL -  X2 * X2 ) * ( 1.0_DBL - X3 ) * 0.25_DBL ; 
DFXI ( 11 , 2 ) = - X2 * ( 1.0_DBL - X1 ) * ( 1.0_DBL - X3 ) * 0.5_DBL  ; 
DFXI ( 11 , 3 ) = -( 1.0_DBL -  X2 * X2 ) * ( 1.0_DBL - X1 ) * 0.25_DBL ; 

DFXI ( 12 , 1 ) = - X1 * ( 1.0_DBL - X2 ) * ( 1.0_DBL - X3 ) * 0.5_DBL  ; 
DFXI ( 12 , 2 ) = -( 1.0_DBL -  X1 * X1 ) * ( 1.0_DBL - X3 ) * 0.25_DBL ; 
DFXI ( 12 , 3 ) = -( 1.0_DBL -  X1 * X1 ) * ( 1.0_DBL - X2 ) * 0.25_DBL ; 

DFXI ( 13 , 1 ) =  ( 1.0_DBL -  X3 * X3 ) * ( 1.0_DBL - X2 ) * 0.25_DBL ; 
DFXI ( 13 , 2 ) = -( 1.0_DBL -  X3 * X3 ) * ( 1.0_DBL + X1 ) * 0.25_DBL ; 
DFXI ( 13 , 3 ) = - X3 * ( 1.0_DBL + X1 ) * ( 1.0_DBL - X2 ) * 0.5_DBL ; 

DFXI ( 14 , 1 ) =  ( 1.0_DBL -  X3 * X3 ) * ( 1.0_DBL + X2 ) * 0.25_DBL ; 
DFXI ( 14 , 2 ) =  ( 1.0_DBL -  X3 * X3 ) * ( 1.0_DBL + X1 ) * 0.25_DBL ; 
DFXI ( 14 , 3 ) = - X3 * ( 1.0_DBL + X1 ) * ( 1.0_DBL + X2 ) * 0.5_DBL ; 

DFXI ( 15 , 1 ) = -( 1.0_DBL -  X3 * X3 ) * ( 1.0_DBL + X2 ) * 0.25_DBL ; 
DFXI ( 15 , 2 ) =  ( 1.0_DBL -  X3 * X3 ) * ( 1.0_DBL - X1 ) * 0.25_DBL ; 
DFXI ( 15 , 3 ) = - X3 * ( 1.0_DBL - X1 ) * ( 1.0_DBL + X2 ) * 0.5_DBL  ; 

DFXI ( 16 , 1 ) = -( 1.0_DBL -  X3 * X3 ) * ( 1.0_DBL - X2 ) * 0.25_DBL ; 
DFXI ( 16 , 2 ) = -( 1.0_DBL -  X3 * X3 ) * ( 1.0_DBL - X1 ) * 0.25_DBL ; 
DFXI ( 16 , 3 ) = - X3 * ( 1.0_DBL - X1 ) * ( 1.0_DBL - X2 ) * 0.5_DBL  ; 

DFXI ( 17 , 1 ) =  ( 1.0_DBL -  X2 * X2 ) * ( 1.0_DBL + X3 ) * 0.25_DBL ; 
DFXI ( 17 , 2 ) = - X2 * ( 1.0_DBL + X1 ) * ( 1.0_DBL + X3 ) * 0.5_DBL  ; 
DFXI ( 17 , 3 ) =  ( 1.0_DBL -  X2 * X2 ) * ( 1.0_DBL + X1 ) * 0.25_DBL ; 

DFXI ( 18 , 1 ) = - X1 * ( 1.0_DBL + X2 ) * ( 1.0_DBL + X3 ) * 0.5_DBL  ; 
DFXI ( 18 , 2 ) =  ( 1.0_DBL -  X1 * X1 ) * ( 1.0_DBL + X3 ) * 0.25_DBL ; 
DFXI ( 18 , 3 ) =  ( 1.0_DBL -  X1 * X1 ) * ( 1.0_DBL + X2 ) * 0.25_DBL ; 

DFXI ( 19 , 1 ) = -( 1.0_DBL -  X2 * X2 ) * ( 1.0_DBL + X3 ) * 0.25_DBL ; 
DFXI ( 19 , 2 ) = - X2 * ( 1.0_DBL - X1 ) * ( 1.0_DBL + X3 ) * 0.5_DBL  ; 
DFXI ( 19 , 3 ) =  ( 1.0_DBL -  X2 * X2 ) * ( 1.0_DBL - X1 ) * 0.25_DBL ; 

DFXI ( 20 , 1 ) = - X1 * ( 1.0_DBL - X2 ) * ( 1.0_DBL + X3 ) * 0.5_DBL  ; 
DFXI ( 20 , 2 ) = -( 1.0_DBL -  X1 * X1 ) * ( 1.0_DBL + X3 ) * 0.25_DBL ; 
DFXI ( 20 , 3 ) =  ( 1.0_DBL -  X1 * X1 ) * ( 1.0_DBL - X2 ) * 0.25_DBL ; 

End Function Dif_ShapeFunc_3_20 ;


!  <><><><><><><><><><><><><><><><><><><><><><><><><><><> 3D 27N Second order Hexahedral element <><><><><><><><><><><><><><><><><><><><><><><><><><>

! ***************************************************************************************************************************************************
!                                                                                                                                                  **
! Last Update:  07 Oct 2012                                                                                                                        **
! Description: This function computes the shape functions of a { 3d second order  27-node Lagrange element }. element number : { ?? }              **
!                                                                                                                                                  **
!                                                                                                                                                  **
! ***************************************************************************************************************************************************

Function ShapeFunc_3_27 ( X1, X2, X3 )   Result ( Funct ) ;

Implicit None ;

Real (Kind=DBL), Intent(In)     :: X1, X2, X3 ;
Real (Kind=DBL), Dimension(27)  :: Funct ;

! =========================== Function CODE =======================================================================================================

Funct ( 1  ) = 0.125_Dbl *  X1 * -X2 * -X3 * (1.0_Dbl+X1) * (1.0_Dbl-X2) * (1.0_Dbl-X3)
Funct ( 2  ) = 0.125_Dbl *  X1 *  X2 * -X3 * (1.0_Dbl+X1) * (1.0_Dbl+X2) * (1.0_Dbl-X3)
Funct ( 3  ) = 0.125_Dbl * -X1 *  X2 * -X3 * (1.0_Dbl-X1) * (1.0_Dbl+X2) * (1.0_Dbl-X3)
Funct ( 4  ) = 0.125_Dbl * -X1 * -X2 * -X3 * (1.0_Dbl-X1) * (1.0_Dbl-X2) * (1.0_Dbl-X3)
Funct ( 5  ) = 0.125_Dbl *  X1 * -X2 *  X3 * (1.0_Dbl+X1) * (1.0_Dbl-X2) * (1.0_Dbl+X3)
Funct ( 6  ) = 0.125_Dbl *  X1 *  X2 *  X3 * (1.0_Dbl+X1) * (1.0_Dbl+X2) * (1.0_Dbl+X3)
Funct ( 7  ) = 0.125_Dbl * -X1 *  X2 *  X3 * (1.0_Dbl-X1) * (1.0_Dbl+X2) * (1.0_Dbl+X3)
Funct ( 8  ) = 0.125_Dbl * -X1 * -X2 *  X3 * (1.0_Dbl-X1) * (1.0_Dbl-X2) * (1.0_Dbl+X3)

Funct ( 9  ) = 0.25_Dbl *  X1 * -X3 * (1.0_Dbl+X1) * (1.0_Dbl-X3) * (1.0_Dbl-X2**2)
Funct ( 11 ) = 0.25_Dbl * -X1 * -X3 * (1.0_Dbl-X1) * (1.0_Dbl-X3) * (1.0_Dbl-X2**2)
Funct ( 17 ) = 0.25_Dbl *  X1 *  X3 * (1.0_Dbl+X1) * (1.0_Dbl+X3) * (1.0_Dbl-X2**2)
Funct ( 19 ) = 0.25_Dbl * -X1 *  X3 * (1.0_Dbl-X1) * (1.0_Dbl+X3) * (1.0_Dbl-X2**2)

Funct ( 10 ) = 0.25_Dbl *  X2 * -X3 * (1.0_Dbl+X2) * (1.0_Dbl-X3) * (1.0_Dbl-X1**2)
Funct ( 12 ) = 0.25_Dbl * -X2 * -X3 * (1.0_Dbl-X2) * (1.0_Dbl-X3) * (1.0_Dbl-X1**2)
Funct ( 18 ) = 0.25_Dbl *  X2 *  X3 * (1.0_Dbl+X2) * (1.0_Dbl+X3) * (1.0_Dbl-X1**2)
Funct ( 20 ) = 0.25_Dbl * -X2 *  X3 * (1.0_Dbl-X2) * (1.0_Dbl+X3) * (1.0_Dbl-X1**2)

Funct ( 13 ) = 0.25_Dbl * -X2 *  X1 * (1.0_Dbl-X2) * (1.0_Dbl+X1) * (1.0_Dbl-X3**2)
Funct ( 14 ) = 0.25_Dbl *  X2 *  X1 * (1.0_Dbl+X2) * (1.0_Dbl+X1) * (1.0_Dbl-X3**2)
Funct ( 15 ) = 0.25_Dbl *  X2 * -X1 * (1.0_Dbl+X2) * (1.0_Dbl-X1) * (1.0_Dbl-X3**2)
Funct ( 16 ) = 0.25_Dbl * -X2 * -X1 * (1.0_Dbl-X2) * (1.0_Dbl-X1) * (1.0_Dbl-X3**2)

Funct ( 21 ) = 0.5_Dbl * -X3 * (1.0_Dbl-X3) * (1.0_Dbl-X1**2) * (1.0_Dbl-X2**2)
Funct ( 26 ) = 0.5_Dbl *  X3 * (1.0_Dbl+X3) * (1.0_Dbl-X1**2) * (1.0_Dbl-X2**2)

Funct ( 22 ) = 0.5_Dbl *  X1 * (1.0_Dbl+X1) * (1.0_Dbl-X2**2) * (1.0_Dbl-X3**2)
Funct ( 24 ) = 0.5_Dbl * -X1 * (1.0_Dbl-X1) * (1.0_Dbl-X2**2) * (1.0_Dbl-X3**2)

Funct ( 23 ) = 0.5_Dbl *  X2 * (1.0_Dbl+X2) * (1.0_Dbl-X1**2) * (1.0_Dbl-X3**2)
Funct ( 25 ) = 0.5_Dbl * -X2 * (1.0_Dbl-X2) * (1.0_Dbl-X1**2) * (1.0_Dbl-X3**2)

Funct ( 27 ) = (1.0_Dbl-X1**2) * (1.0_Dbl-X2**2) * (1.0_Dbl-X3**2)

End Function ShapeFunc_3_27 ;

! ***************************************************************************************************************************************************
!                                                                                                                                                  **
! Last Update:  07 Oct 2012                                                                                                                        **
! Description: This function computes the differentials of shape functions of  a                                                                   **
!              { 3d second order  27-node Lagrangian element }. element number : { 6 }.                                                            **
!                                                                                                                                                  **
! ***************************************************************************************************************************************************


Function Dif_ShapeFunc_3_27 ( X1, X2, X3 )   Result ( DFXI ) ;


Implicit None ;

Real (Kind=DBL), Intent(In)     :: X1, X2, X3 ;
Real (Kind=DBL), Dimension(27,3) :: DFXI ;

! =========================== Function CODE =======================================================================================================

DFXI(1,1) = 0.125_Dbl * -X2 * -X3 * (1.0_Dbl-X2) * (1.0_Dbl-X3) * ( 1.0_Dbl + 2.0_Dbl * X1 ) 
DFXI(2,1) = 0.125_Dbl *  X2 * -X3 * (1.0_Dbl+X2) * (1.0_Dbl-X3) * ( 1.0_Dbl + 2.0_Dbl * X1 )
DFXI(3,1) = 0.125_Dbl *  X2 * -X3 * (1.0_Dbl+X2) * (1.0_Dbl-X3) * (-1.0_Dbl + 2.0_Dbl * X1 )
DFXI(4,1) = 0.125_Dbl * -X2 * -X3 * (1.0_Dbl-X2) * (1.0_Dbl-X3) * (-1.0_Dbl + 2.0_Dbl * X1 )
DFXI(5,1) = 0.125_Dbl * -X2 *  X3 * (1.0_Dbl-X2) * (1.0_Dbl+X3) * ( 1.0_Dbl + 2.0_Dbl * X1 )
DFXI(6,1) = 0.125_Dbl *  X2 *  X3 * (1.0_Dbl+X2) * (1.0_Dbl+X3) * ( 1.0_Dbl + 2.0_Dbl * X1 )
DFXI(7,1) = 0.125_Dbl *  X2 *  X3 * (1.0_Dbl+X2) * (1.0_Dbl+X3) * (-1.0_Dbl + 2.0_Dbl * X1 )
DFXI(8,1) = 0.125_Dbl * -X2 *  X3 * (1.0_Dbl-X2) * (1.0_Dbl+X3) * (-1.0_Dbl + 2.0_Dbl * X1 )

DFXI(9 ,1) = 0.25_Dbl * -X3 * (1.0_Dbl-X3) * (1.0_Dbl-X2**2) * ( 1.0_Dbl + 2.0_Dbl * X1 )
DFXI(11,1) = 0.25_Dbl * -X3 * (1.0_Dbl-X3) * (1.0_Dbl-X2**2) * (-1.0_Dbl + 2.0_Dbl * X1 )
DFXI(17,1) = 0.25_Dbl *  X3 * (1.0_Dbl+X3) * (1.0_Dbl-X2**2) * ( 1.0_Dbl + 2.0_Dbl * X1 )
DFXI(19,1) = 0.25_Dbl *  X3 * (1.0_Dbl+X3) * (1.0_Dbl-X2**2) * (-1.0_Dbl + 2.0_Dbl * X1 )

DFXI(10,1) = 0.5_Dbl *  X2 * -X3 * (1.0_Dbl+X2) * (1.0_Dbl-X3) * (-X1)
DFXI(12,1) = 0.5_Dbl * -X2 * -X3 * (1.0_Dbl-X2) * (1.0_Dbl-X3) * (-X1)
DFXI(18,1) = 0.5_Dbl *  X2 *  X3 * (1.0_Dbl+X2) * (1.0_Dbl+X3) * (-X1)
DFXI(20,1) = 0.5_Dbl * -X2 *  X3 * (1.0_Dbl-X2) * (1.0_Dbl+X3) * (-X1)

DFXI(13,1) = 0.25_Dbl * -X2 * (1.0_Dbl-X2) * (1.0_Dbl-X3**2) * ( 1.0_Dbl + 2.0_Dbl * X1 )
DFXI(14,1) = 0.25_Dbl *  X2 * (1.0_Dbl+X2) * (1.0_Dbl-X3**2) * ( 1.0_Dbl + 2.0_Dbl * X1 )
DFXI(15,1) = 0.25_Dbl *  X2 * (1.0_Dbl+X2) * (1.0_Dbl-X3**2) * (-1.0_Dbl + 2.0_Dbl * X1 )
DFXI(16,1) = 0.25_Dbl * -X2 * (1.0_Dbl-X2) * (1.0_Dbl-X3**2) * (-1.0_Dbl + 2.0_Dbl * X1 )

DFXI(21,1) = -X3 * (1.0_Dbl-X3) * (-X1) * (1.0_Dbl-X2**2)
DFXI(26,1) =  X3 * (1.0_Dbl+X3) * (-X1) * (1.0_Dbl-X2**2)

DFXI(22,1) = 0.5_Dbl * (1.0_Dbl-X2**2) * (1.0_Dbl-X3**2) * ( 1.0_Dbl + 2.0_Dbl * X1 )
DFXI(24,1) = 0.5_Dbl * (1.0_Dbl-X2**2) * (1.0_Dbl-X3**2) * (-1.0_Dbl + 2.0_Dbl * X1 )

DFXI(23,1) =  X2 * (1.0_Dbl+X2) * (-X1) * (1.0_Dbl-X3**2)
DFXI(25,1) = -X2 * (1.0_Dbl-X2) * (-X1) * (1.0_Dbl-X3**2)

DFXI(27,1) = (-2.0_Dbl*X1) * (1.0_Dbl-X2**2) * (1.0_Dbl-X3**2)


DFXI(1,2) = 0.125_Dbl *  X1 * -X3 * (1.0_Dbl+X1) * (-1.0_Dbl+2.0_Dbl*X2) * (1.0_Dbl-X3)
DFXI(2,2) = 0.125_Dbl *  X1 * -X3 * (1.0_Dbl+X1) * (+1.0_Dbl+2.0_Dbl*X2) * (1.0_Dbl-X3)
DFXI(3,2) = 0.125_Dbl * -X1 * -X3 * (1.0_Dbl-X1) * (+1.0_Dbl+2.0_Dbl*X2) * (1.0_Dbl-X3)
DFXI(4,2) = 0.125_Dbl * -X1 * -X3 * (1.0_Dbl-X1) * (-1.0_Dbl+2.0_Dbl*X2) * (1.0_Dbl-X3)
DFXI(5,2) = 0.125_Dbl *  X1 *  X3 * (1.0_Dbl+X1) * (-1.0_Dbl+2.0_Dbl*X2) * (1.0_Dbl+X3)
DFXI(6,2) = 0.125_Dbl *  X1 *  X3 * (1.0_Dbl+X1) * (+1.0_Dbl+2.0_Dbl*X2) * (1.0_Dbl+X3)
DFXI(7,2) = 0.125_Dbl * -X1 *  X3 * (1.0_Dbl-X1) * (+1.0_Dbl+2.0_Dbl*X2) * (1.0_Dbl+X3)
DFXI(8,2) = 0.125_Dbl * -X1 *  X3 * (1.0_Dbl-X1) * (-1.0_Dbl+2.0_Dbl*X2) * (1.0_Dbl+X3)

DFXI(9 ,2) = 0.5_Dbl *  X1 * -X3 * (1.0_Dbl+X1) * (1.0_Dbl-X3) * (-X2)
DFXI(11,2) = 0.5_Dbl * -X1 * -X3 * (1.0_Dbl-X1) * (1.0_Dbl-X3) * (-X2)
DFXI(17,2) = 0.5_Dbl *  X1 *  X3 * (1.0_Dbl+X1) * (1.0_Dbl+X3) * (-X2)
DFXI(19,2) = 0.5_Dbl * -X1 *  X3 * (1.0_Dbl-X1) * (1.0_Dbl+X3) * (-X2)

DFXI(10,2) = 0.25_Dbl * -X3 * (+1.0_Dbl+2.0_Dbl*X2) * (1.0_Dbl-X3) * (1.0_Dbl-X1**2)
DFXI(12,2) = 0.25_Dbl * -X3 * (-1.0_Dbl+2.0_Dbl*X2) * (1.0_Dbl-X3) * (1.0_Dbl-X1**2)
DFXI(18,2) = 0.25_Dbl *  X3 * (+1.0_Dbl+2.0_Dbl*X2) * (1.0_Dbl+X3) * (1.0_Dbl-X1**2)
DFXI(20,2) = 0.25_Dbl *  X3 * (-1.0_Dbl+2.0_Dbl*X2) * (1.0_Dbl+X3) * (1.0_Dbl-X1**2)

DFXI(13,2) = 0.25_Dbl *  X1 * (-1.0_Dbl+2.0_Dbl*X2) * (1.0_Dbl+X1) * (1.0_Dbl-X3**2)
DFXI(14,2) = 0.25_Dbl *  X1 * (+1.0_Dbl+2.0_Dbl*X2) * (1.0_Dbl+X1) * (1.0_Dbl-X3**2)
DFXI(15,2) = 0.25_Dbl * -X1 * (+1.0_Dbl+2.0_Dbl*X2) * (1.0_Dbl-X1) * (1.0_Dbl-X3**2)
DFXI(16,2) = 0.25_Dbl * -X1 * (-1.0_Dbl+2.0_Dbl*X2) * (1.0_Dbl-X1) * (1.0_Dbl-X3**2)

DFXI(21,2) = -X3 * (1.0_Dbl-X3) * (1.0_Dbl-X1**2) * (-X2)
DFXI(26,2) =  X3 * (1.0_Dbl+X3) * (1.0_Dbl-X1**2) * (-X2)

DFXI(22,2) =  X1 * (1.0_Dbl+X1) * (-X2) * (1.0_Dbl-X3**2)
DFXI(24,2) = -X1 * (1.0_Dbl-X1) * (-X2) * (1.0_Dbl-X3**2)

DFXI(23,2) = 0.5_Dbl * (+1.0_Dbl+2.0_Dbl*X2) * (1.0_Dbl-X1**2) * (1.0_Dbl-X3**2)
DFXI(25,2) = 0.5_Dbl * (-1.0_Dbl+2.0_Dbl*X2) * (1.0_Dbl-X1**2) * (1.0_Dbl-X3**2)

DFXI(27,2) = (1.0_Dbl-X1**2) * (-2.0_Dbl*X2) * (1.0_Dbl-X3**2)


DFXI(1,3) = 0.125_Dbl *  X1 * -X2 * (1.0_Dbl+X1) * (1.0_Dbl-X2) * (-1.0_Dbl+2.0_Dbl*X3)
DFXI(2,3) = 0.125_Dbl *  X1 *  X2 * (1.0_Dbl+X1) * (1.0_Dbl+X2) * (-1.0_Dbl+2.0_Dbl*X3)
DFXI(3,3) = 0.125_Dbl * -X1 *  X2 * (1.0_Dbl-X1) * (1.0_Dbl+X2) * (-1.0_Dbl+2.0_Dbl*X3)
DFXI(4,3) = 0.125_Dbl * -X1 * -X2 * (1.0_Dbl-X1) * (1.0_Dbl-X2) * (-1.0_Dbl+2.0_Dbl*X3)
DFXI(5,3) = 0.125_Dbl *  X1 * -X2 * (1.0_Dbl+X1) * (1.0_Dbl-X2) * (+1.0_Dbl+2.0_Dbl*X3)
DFXI(6,3) = 0.125_Dbl *  X1 *  X2 * (1.0_Dbl+X1) * (1.0_Dbl+X2) * (+1.0_Dbl+2.0_Dbl*X3)
DFXI(7,3) = 0.125_Dbl * -X1 *  X2 * (1.0_Dbl-X1) * (1.0_Dbl+X2) * (+1.0_Dbl+2.0_Dbl*X3)
DFXI(8,3) = 0.125_Dbl * -X1 * -X2 * (1.0_Dbl-X1) * (1.0_Dbl-X2) * (+1.0_Dbl+2.0_Dbl*X3)

DFXI(9 ,3) = 0.25_Dbl *  X1 * (1.0_Dbl+X1) * (-1.0_Dbl+2.0_Dbl*X3) * (1.0_Dbl-X2**2)
DFXI(11,3) = 0.25_Dbl * -X1 * (1.0_Dbl-X1) * (-1.0_Dbl+2.0_Dbl*X3) * (1.0_Dbl-X2**2)
DFXI(17,3) = 0.25_Dbl *  X1 * (1.0_Dbl+X1) * (+1.0_Dbl+2.0_Dbl*X3) * (1.0_Dbl-X2**2)
DFXI(19,3) = 0.25_Dbl * -X1 * (1.0_Dbl-X1) * (+1.0_Dbl+2.0_Dbl*X3) * (1.0_Dbl-X2**2)

DFXI(10,3) = 0.25_Dbl *  X2 * (1.0_Dbl+X2) * (-1.0_Dbl+2.0_Dbl*X3) * (1.0_Dbl-X1**2)
DFXI(12,3) = 0.25_Dbl * -X2 * (1.0_Dbl-X2) * (-1.0_Dbl+2.0_Dbl*X3) * (1.0_Dbl-X1**2)
DFXI(18,3) = 0.25_Dbl *  X2 * (1.0_Dbl+X2) * (+1.0_Dbl+2.0_Dbl*X3) * (1.0_Dbl-X1**2)
DFXI(20,3) = 0.25_Dbl * -X2 * (1.0_Dbl-X2) * (+1.0_Dbl+2.0_Dbl*X3) * (1.0_Dbl-X1**2)

DFXI(13,3) = 0.5_Dbl * -X2 *  X1 * (1.0_Dbl-X2) * (1.0_Dbl+X1) * (-X3)
DFXI(14,3) = 0.5_Dbl *  X2 *  X1 * (1.0_Dbl+X2) * (1.0_Dbl+X1) * (-X3)
DFXI(15,3) = 0.5_Dbl *  X2 * -X1 * (1.0_Dbl+X2) * (1.0_Dbl-X1) * (-X3)
DFXI(16,3) = 0.5_Dbl * -X2 * -X1 * (1.0_Dbl-X2) * (1.0_Dbl-X1) * (-X3)

DFXI(21,3) = 0.5_Dbl * (-1.0_Dbl+2.0_Dbl*X3) * (1.0_Dbl-X1**2) * (1.0_Dbl-X2**2)
DFXI(26,3) = 0.5_Dbl * (+1.0_Dbl+2.0_Dbl*X3) * (1.0_Dbl-X1**2) * (1.0_Dbl-X2**2)

DFXI(22,3) =  X1 * (1.0_Dbl+X1) * (1.0_Dbl-X2**2) * (-X3)
DFXI(24,3) = -X1 * (1.0_Dbl-X1) * (1.0_Dbl-X2**2) * (-X3)

DFXI(23,3) =  X2 * (1.0_Dbl+X2) * (1.0_Dbl-X1**2) * (-X3)
DFXI(25,3) = -X2 * (1.0_Dbl-X2) * (1.0_Dbl-X1**2) * (-X3)

DFXI(27,3) = (1.0_Dbl-X1**2) * (1.0_Dbl-X2**2) * (-2.0_Dbl*X3)

End Function Dif_ShapeFunc_3_27 ;


End Module ShapeFunctions ;
