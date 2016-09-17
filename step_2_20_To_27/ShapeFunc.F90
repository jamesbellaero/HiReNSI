Module ShapeFunctions ;

Implicit none ;

contains;


! ***************************************************************************************************************************************************
!                                                                                                                                                  **
! Last Update:  14 May 2012                                                                                                                        **
! Description: This function computes the shape functions of a { 3d second order  20-node serendepiti element }. element number : { 6 }            **
!                                                                                                                                                  **
!                                                                                                                                                  **
! ***************************************************************************************************************************************************

Function ShapeFunc_3_20 ( X1, X2, X3 )   Result ( Funct ) ;

Implicit None ;

Integer(2), PARAMETER  :: DBL  = SELECTED_Real_Kind ( P = 13, R = 200 ) ;  ! EQUIVALENT TO Real (8) 

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

Function Dif_ShapeFunc_3_20 ( X1, X2, X3 )   Result ( DFXI ) ;


Implicit None ;

Integer(2), PARAMETER  :: DBL  = SELECTED_Real_Kind ( P = 13, R = 200 ) ;  ! EQUIVALENT TO Real (8) 

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


End Module ShapeFunctions