

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                                                                                                                  ++
! Start:        15 May 2011                                                                                                                        ++
! Last Update:  4 April 2013                                                                                                                       ++
!                                                                                                                                                  ++
! Description: THIS Module Contains ALL Input Subroutines.                                                                                         ++
!                                                                                                                                                  ++
!              Subroutines                     ELEMENT NUMBER                                                                                      ++
!              ==============                  ==============                                                                                      ++
!              Input_BASIC                     COLLECTS ALL BASIC Variables.                                                                       ++
!              Input_Arrays                    COLLECTS ALL Arrays                                                                                 ++
!                                                                                                                                                  ++
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Module Input_Subroutines ;

Use Parameters ;

Implicit None ;

  Interface Input ;
    Module Procedure Input_Basic, Input_Arrays, Input_Analysis, Input_Het_Mat, Input_Inv_DS_Basic, Input_Inv_DS_Array, Input_Inv_Measured_Data, Input_Mat_Vis_Basic, Input_Het_Mat_Resume ;
  End Interface ;

Contains ;

! ***************************************************************************************************************************************************
!                                                                                                                                                  **
! Start:        15 May 2011                                                                                                                        ++
! Last Update:  4 April 2013                                                                                                                       **
! Description: THIS Subroutine COLLECTS ALL BASIC DATA FOR EACH ANALYSIS.                                                                          **
!                                                                                                                                                  **
! ***************************************************************************************************************************************************

Subroutine Input_BASIC (                                                                                                        &
NDOF, MaxNNode, NDim,                                                                                                           & ! Integer (1) Variables
NGroup, NPM, NMat,                                                                                                              & ! Integer (2) Variables
NEL, NJ, NJTotal, NEQM, NEQMTotal, NEQM_Mapping,                                                                                & ! Integer (4) Variables
!                                                                                                                               & ! Integer (8) Variables
LoadC,                                                                                                                          & ! Real Variables
!                                                                                                                               & ! Integer Arrays
!                                                                                                                               & ! Real Arrays
!                                                                                                                               & ! Characters
Param                                                                                                                           & ! Type
) ;


Implicit None ;

! =========================== Global Variables ======================================================================================================

! - Integer Variables -------------------------------------------------------------------------------------------------------------------------------
Integer, Intent(Out)    :: NDOF, MaxNNode, NDim ;

Integer, Intent(Out)    :: NGroup, NPM, NMat ;

Integer, Intent(Out)    :: NEL, NJ, NJTotal, NEQM, NEQMTotal, NEQM_Mapping ;

! - Real Variables ----------------------------------------------------------------------------------------------------------------------------------
! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
! - Integer Arrays ----------------------------------------------------------------------------------------------------------------------------------
Integer, Intent(OUT)  , Dimension (:  )  :: LoadC ;

! - Real Arrays -------------------------------------------------------------------------------------------------------------------------------------
! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
! - Character Variables -----------------------------------------------------------------------------------------------------------------------------
! - Logical Variables -------------------------------------------------------------------------------------------------------------------------------
! - Type DECLERATIONS -------------------------------------------------------------------------------------------------------------------------------
Type ( BasicParam )    :: Param ;                ! Holds basic parameters of each load case

! =========================== Local Variables =======================================================================================================

! - Integer Variables -------------------------------------------------------------------------------------------------------------------------------
Integer  :: I, J, K ;               ! Loop Index

Integer  :: IO_Read ;               ! Holds error of read statements
Integer  :: UnFile ;                ! Holds Unit of a file for error message

! - Real Variables ----------------------------------------------------------------------------------------------------------------------------------
! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
! - Integer Arrays ----------------------------------------------------------------------------------------------------------------------------------
! - Real Arrays -------------------------------------------------------------------------------------------------------------------------------------
! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
! - Character Variables -----------------------------------------------------------------------------------------------------------------------------
! - Logical Variables -------------------------------------------------------------------------------------------------------------------------------
! - Type DECLERATIONS -------------------------------------------------------------------------------------------------------------------------------

! =========================== Subroutine CODE =======================================================================================================

Write(*,*)
Write(*,*)"<<<<<<<< Input Basic >>>>>>>>" ;

! Reading data from .dataModel

UnFile = UnInptMdl ;

! Basic Parameters for model
Write(*,*)"  Reading Basic Parameters for model ..." ;
Read (Unit = UnFile, FMT = "(3(I3,2X),1(I5,2X),6(I19,2X))", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Read, ERR = 1003, End = 1004, EOR = 1005 ) NDOF, NDim, MaxNNode,      NGroup,              NEL, NJ, NJTotal, NEQM, NEQMTotal, NEQM_Mapping ;

! Available Load Cases
Write(*,*)"  Reading Load Cases ..." ;
Read (Unit = UnFile, FMT = "(5(I3,2X))", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Read, ERR = 1003, End = 1004, EOR = 1005 ) ( LoadC ( I ), I = 1, 5 ) ;

NPM    = 8 ;
NMat   = NGroup ;

Allocate ( Param%IntM ( 7, 6),Param%RealM ( 7, 6)  ) ;   !?? modify this

! Reading Param Model
Write(*,*)"Reading Param Array ..." ;
Read (Unit = UnFile, FMT = "(<6*7>(I19,1X))", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Read, ERR = 1003, End = 1004, EOR = 1005 )  (( Param%IntM  ( I, J), I = 1, 7), J = 1, 6) ;
Read (Unit = UnFile, FMT = "(<6*7>(E31.23E3,1X))", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Read, ERR = 1003, End = 1004, EOR = 1005 ) (( Param%RealM  ( I, J), I = 1, 7), J = 1, 6) ;


Write(*    ,*) 'End Subroutine < Input_BASIC >' ;
Write(UnInf,*) 'End Subroutine < Input_BASIC >' ;
Return ;

! - ERROR IN READ STATEMENT -------------------------------------------------------------------------------------------------------------------------
1003    Write(*     , Fmt_READ1 ) UnFile, IO_READ ; Write( UnFile, Fmt_READ1 ) UnFile, IO_READ ;
        !#Call BEEP_FAIL ;
        Write(*, Fmt_FL) ;  Write(UnInf, Fmt_FL) ; Write(*, Fmt_End) ; Read(*,*) ;  Return ;

! - End-OF-FILE IN READ STATEMENT -------------------------------------------------------------------------------------------------------------------
1004    Write(*     , Fmt_READ2 ) UnFile, IO_READ ; Write( UnFile, Fmt_READ2 ) UnFile, IO_READ ;
        !#Call BEEP_FAIL ;
        Write(*, Fmt_FL) ;  Write(UnInf, Fmt_FL) ; Write(*, Fmt_End) ; Read(*,*) ;  Return ;

! - End-OF-FILE IN READ STATEMENT -------------------------------------------------------------------------------------------------------------------
1005    Write(*     , Fmt_READ3 ) UnFile, IO_READ ; Write( UnFile, Fmt_READ3 ) UnFile, IO_READ ;
        !#Call BEEP_FAIL ;
        Write(*, Fmt_FL) ;  Write(UnInf, Fmt_FL) ; Write(*, Fmt_End) ; Read(*,*) ;  Return ;

End Subroutine Input_BASIC ;

! ***************************************************************************************************************************************************
!                                                                                                                                                  **
! Last Update:  4 April 2013                                                                                                                       **
! Description: THIS Subroutine COLLECTS ALL Arrays.                                                                                                **
!                                                                                                                                                  **
! ***************************************************************************************************************************************************

Subroutine Input_Arrays (                                                                                                       &
NDOF, NDim, MaxNNode,                                                                                                           & ! Integer (1) Variables
NGroup, NPM, NMat,     NLCase,                                                                                                  & ! Integer (2) Variables
!                                                                                                                               & ! Integer (4) Variables
NEL, NJ, NEQM, NEQM_Mapping,                                                                                                    & ! Integer (8) Variables
!                                                                                                                               & ! Real Variables
LoadC,     IDBC,      MTEL, ELT,     ELGR,                                                                                      &
D_NNZ_Stiff, O_NNZ_Stiff, D_NNZ_Damp, O_NNZ_Damp, D_NNZ_Mass, O_NNZ_Mass,     App_Numbers, PETSc_Numbers,     Indices,          &
NDAN , NVAN, NAAN, NDPN , NVPN, NAPN,     NoBndry_DRM, NoLayer_DRM, JLoad,     EqDis, EqVel, EqAcc,      INOD, ID,              & ! Integer Arrays
PMat, PBLD, XYZ, UDis, PLoad, PML_DIM,                                                                                          & ! Real Arrays
!                                                                                                                               & ! Characters
Param                                                                                                                           & ! Type
) ;

Implicit None ;

#include "finclude/petscsys.h"
#include "finclude/petscmat.h"
#include "finclude/petscvec.h90"

PetscMPIInt    :: Rank ;                         ! Rank number

! =========================== Global Variables ======================================================================================================

! - Integer Arrays ------------------------------------------------------------------------------------------------------------------------------
Integer , Intent(In)    :: NDOF, NDim, MaxNNode ;

Integer , Intent(In)    :: NGroup, NPM, NMat ;
Integer , Intent(Out)   :: NLCase ;

Integer , Intent(In)    :: NEL, NJ, NEQM, NEQM_Mapping ;

! - Real Variables ----------------------------------------------------------------------------------------------------------------------------------
! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
! - Integer Arrays ----------------------------------------------------------------------------------------------------------------------------------
Integer , Intent(In),  Dimension (:  )  :: LoadC ;

Integer , Intent(Out), Dimension (:  )  :: MTEL, ELT ;

Integer , Intent(Out), Dimension (:  )  :: ELGR ;
Integer , Intent(Out), Dimension (:  )  :: D_NNZ_Stiff, O_NNZ_Stiff, D_NNZ_Damp, O_NNZ_Damp, D_NNZ_Mass, O_NNZ_Mass ;
Integer , Intent(Out), Dimension ( 0:NEQM_Mapping - 1  )  :: App_Numbers, PETSc_Numbers ;
Integer , Intent(Out), Dimension ( 0:NEQM - 1  )          :: Indices ;

Integer , Intent(Out), Dimension (:  )  :: NDAN , NVAN, NAAN, NDPN , NVPN, NAPN ;
Integer , Intent(Out), Dimension (:  )  :: NoBndry_DRM, NoLayer_DRM, JLoad ;
Integer , Intent(Out), Dimension (:  )  :: EqDis, EqVel, EqAcc ;
Integer , Intent(Out), Dimension (:,:)  :: INOD, ID ;
Integer , Intent(Out), Dimension (:,:)  :: IDBC ;

! - Real Arrays -------------------------------------------------------------------------------------------------------------------------------------
Real (8),     Intent(Out), Dimension (:,:)  :: PMat, PBLD, XYZ, UDis, PLoad, PML_DIM ;

! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
! - Character Variables -----------------------------------------------------------------------------------------------------------------------------
! - Logical Variables -------------------------------------------------------------------------------------------------------------------------------
! - Type DECLERATIONS -------------------------------------------------------------------------------------------------------------------------------
Type ( BasicParam )    :: Param ;                ! Holds basic parameters of each load case

! =========================== Local Variables =======================================================================================================
! - Integer Variables -------------------------------------------------------------------------------------------------------------------------------
Integer   :: IO_Read ;                ! Holds error of read statements.
Integer   :: UnFile ;                 ! Holds Unit of a file for error message.
Integer   :: INode;                   ! Loop counter for element node number.
Integer   :: IDOF ;                   ! Loop counter on the degree of freedom.

Integer   :: I1, I2 ;                 ! Loop indeces.
Integer   :: NODE ;                   ! Some variable for holding node number.

Integer   :: I, J, K ;                ! Loop indeces.
Integer   :: IEL ;                    ! Loop indeces over NEL - the number of elements.
Integer   :: IJ ;                     ! Loop index on NJ.

! - Real Variables ----------------------------------------------------------------------------------------------------------------------------------
! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
! - Integer Arrays ----------------------------------------------------------------------------------------------------------------------------------
! - Real Arrays -------------------------------------------------------------------------------------------------------------------------------------
! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
! - Character Variables -----------------------------------------------------------------------------------------------------------------------------
! - Logical Variables -------------------------------------------------------------------------------------------------------------------------------
! - Type DECLERATIONS -------------------------------------------------------------------------------------------------------------------------------

! =========================== Subroutine CODE =======================================================================================================

Write(*,*)
Write(*,*)"<<<<<<<< Input Arrays >>>>>>>>>" ;

NLCase = MaxVal ( LoadC ) ;

! - Geometry & Model -------------------------------------------------------------------------------------------------------------------------------

!! Number of elements of each group
!Write(*,*)"Reading Group numbers ..." ;
!  DO I = 1, NGroup ;
!    Read (Unit = UnFile, FMT = "(I19,2X,I5,2X)", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Read, ERR = 1003, End = 1004, EOR = 1005 ) NEG ( I, 1 ), NEG ( I, 2 ) ;
!  End Do ;

! Coordinates of nodes
UnFile = UnInptXYZ ;
Write(*,*)"Reading Coordinates ..." ;
  DO K = 1, NJ ;
    Read (Unit = UnFile, FMT = "(<NDim>(E31.23E3,2X))", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Read, ERR = 1003, End = 1004, EOR = 1005 )( XYZ ( K, J ), J = 1, NDim ) ;
  End Do ;

! Element connectivities
UnFile = UnInptCnn ;
Write(*,*)"Reading Element Connectivities ..." ;
  DO IEL = 1, NEL ;
    Read (Unit = UnFile, FMT = "(<MaxNNode>(I19,2X),I5,2X,I5,2X,I5)", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Read, ERR = 1003, End = 1004, EOR = 1005 ) ( INOD ( INode, IEL ), INode = 1, MaxNNode ), MTEL ( IEL ), ELT ( IEL ), ELGR ( IEL ) ;
  End Do ;

! Constraints
UnFile = UnInptCnt ;
Write(*,*)"Reading Constraints ..." ;
  DO IJ = 1, NJ ;
    Read (Unit = UnFile, FMT = "(<NDOF>(I19,2X))", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Read, ERR = 1003, End = 1004, EOR = 1005 ) ( ID ( IJ, IDOF ), IDOF = 1, NDOF ) ;
  End Do ;


UnFile = UnInptMdl ;

! PML teritory
  If ( LoadC (5) /= 0 ) Then ;
    Write(*,*)"Reading PML Teritory ..." ;
      DO J = 1, 2
        Read (Unit = UnFile, FMT = "(<2*NDim>(E31.23E3,2X))", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Read, ERR = 1003, End = 1004, EOR = 1005 ) (PML_DIM ( I, J ), I = 1, 2 * NDim ) ; 
      End Do ;
  End If ;

! Body force load
  If ( LoadC (1) /= 0 ) Then ;
    Write(*,*)"Reading Body Forces ..." ;
      DO I = 1, Param%IntM( 1, 1) ; ! Param%IntP( 1, 1)=NBLD: Number of body force load.   Param%IntP( 1, 2) = NPBL
        Read (Unit = UnFile, FMT = "(<Param%IntM( 1, 2)>(E31.23E3,2x))", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Read, ERR = 1003, End = 1004, EOR = 1005 ) ( PBLD ( I, J ) , J = 1, Param%IntM( 1, 2) ) ;  ! Param%IntM( 1, 2) = NDIm
      End Do ;

!    ! Load Types
!    LTEL = 0 ;
!    Read (Unit = UnFile, FMT = "(I10)", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Read, ERR = 1003, End = 1004, EOR = 1005 ) LDSET ;
!      DO I = 1, LDSET ;
!        Read (Unit = UnFile, FMT = "(<4>(I19,2X))", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Read, ERR = 1003, End = 1004, EOR = 1005 ) I1, I2, Inc, LType ;
!          DO J = I1, I2, Inc ;
!            LTEL ( J ) = LType ;
!          End Do ;
!      End Do ;

  End If ;

! Static Pressure
  If ( LoadC (2) /= 0 ) Then ;
    Write(*,*)"Reading Pressure Loads ..." ;
      DO J = 1,  Param%IntM( 2, 4 ) ; ! NIDBC
        Read (Unit = UnFile, FMT = "((I19,2X),<2*NDim>(I3,2X))", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Read, ERR = 1003, End = 1004, EOR = 1005 ) IDBC ( J, 1 ), ( IDBC ( J, I ), I = 2, 2 * NDim+1 ) ;
      End Do ;
  End If ;

  ! Joint load of every case - concentrated loads
  IF ( LoadC ( 3 ) /= 0 ) Then ;
    Write(*,*)"Reading Joint Loads ..." ;
      DO I = 1, Param%IntM( 3, 1) ; ! Param%IntP( 3, 1) = NLN
        Read (Unit = UnFile, FMT = "(I19,2X,<NDim>(E31.23E3,2X))", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Read, ERR = 1003, End = 1004, EOR = 1005 ) JLoad ( I ), ( PLoad( J, I ), J = 1, NDim ) ;
      End Do
  End If ;

! Supports Displacements
  IF ( LoadC ( 4 ) /= 0 ) Then ;
    Write(*,*)"Reading Support Displacements ..." ;
      DO I = 1, Param%IntM( 4, 1 ) ; ! Param%IntM( 4, 1 ) = NSND
        Read (Unit = UnFile, FMT = "(<NDOF+1>(E31.23E3,2X))", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Read, ERR = 1003, End = 1004, EOR = 1005 ) UDis( NDOF + 1, I ), ( UDis ( K, I ), K = 1, NDOF ) ;
      End Do ;
  End If ;

Write(*,*)"Reading Dynamic Loads ..." ;

! Dynamic loads
  If ( LoadC (5) /= 0 ) Then ;

    If ( Param%IntM( 5, 3 ) == 1 ) Then ;    ! RICKER PULSE or sine function for dynamic pressure
      Write(*,*)"Reading IDBC ..." ;
        DO J = 1, Param%IntM( 2, 4 ) ; ! NIDBC
          Read (Unit = UnFile, FMT = "((I19,2X),<2*NDim>(I3,2X))", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Read, ERR = 1003, End = 1004, EOR = 1005 ) IDBC ( J, 1 ), ( IDBC ( J, I ), I = 2, 2 * NDim+1 ) ;
        End Do ;

!        Else If ( Param%IntM( 5, 3 ) == 2_Tiny ) Then ;   ! Base acceleration
!
!          Write(*,*)"Reading Base Acceleration ..." ;
!            DO I = 1, Param%IntL( 7, 1) ; ! NAStep = Param%IntP( 7, 1)
!              Read (Unit = UnFile, FMT = "(<NDim>(E31.23E3,2X))", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Read, ERR = 1003, End = 1004, EOR = 1005 ) ( BACL ( I, K ), K = 1, NDim ) ;
!            End Do ;
!
    Else If ( Param%IntM( 5, 3 ) == 3 ) Then ;    ! Domain Reduction Method (DRM)

!          Write(*,*)"Reading wave information ..." ;
!          Read (Unit = UnFile, FMT = "(5(E28.20E3,1x))", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Read, ERR = 1003, End = 1004, EOR = 1005 ) InciWave(1), InciWave(2), InciWave(3), InciWave(4), InciWave(5) ;       ! Theta, Omega, amplitude, alpha1, alpha2

      ! Reading the node numbers on the DRM boundary
      Write(*,*)"Reading DRM boundary nodes ..."  ;
      Read (Unit = UnFile, FMT = "(<Param%IntM( 7, 1)>(I19,1X))", ADVANCE = 'Yes', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Read, ERR = 1003, End = 1004, EOR = 1005 ) ( NoBndry_DRM ( I ), I = 1, Param%IntM( 7, 1) ) ; ! NNBndry_DRM

      ! Reading the node numbers on the DRM neighbor
      Write(*,*)"Reading DRM neighbor nodes ..." ;
      Read (Unit = UnFile, FMT = "(<Param%IntM( 7, 2)>(I19,1X))", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Read, ERR = 1003, End = 1004, EOR = 1005 ) ( NoLayer_DRM ( I ), I = 1, Param%IntM( 7, 2)) ;  ! NNLayer_DRM

    End If ;

  End If ;

  ! Reads down node numbers and equation numbers of nodes in which history of Displacement is required
  If ( Param%IntM( 6, 1) /= 0 ) Then ; ! NNDH
    Write(*,*)"Reading Node Numbers for displacements ..." ;
    Read (Unit = UnFile, FMT = "(<Param%IntM( 6, 1)>(I19,2X),5x,<Param%IntM( 6, 1)>(I19,2X))", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Read, ERR = 1003, End = 1004, EOR = 1005 ) ( NDAN ( I ), I = 1, Param%IntM( 6, 1) ), ( NDPN ( I ), I = 1, Param%IntM( 6, 1) ) ;
    Read (Unit = UnFile, FMT = "(<Param%IntM( 6, 1)*NDIM>(I19,2X))"             , ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Read, ERR = 1003, End = 1004, EOR = 1005 ) ( EqDis ( I ), I = 1, Param%IntM( 6, 1)* NDIM ) ;
  End If ;

  ! Reads down node numbers and equation numbers of nodes in which history of velocity is required
  If ( Param%IntM( 6, 2) /= 0 ) Then ; ! NNVH
    Write(*,*)"Reading Node Numbers for velocity ..." ;
    Read (Unit = UnFile, FMT = "(<Param%IntM( 6, 2)>(I19,2X),5x,<Param%IntM( 6, 2)>(I19,2X))", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Read, ERR = 1003, End = 1004, EOR = 1005 ) ( NVAN ( I ), I = 1, Param%IntM( 6, 2) ), ( NVPN ( I ), I = 1, Param%IntM( 6, 2) ) ;
    Read (Unit = UnFile, FMT = "(<Param%IntM( 6, 2)*NDIM>(I19,2X))"             , ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Read, ERR = 1003, End = 1004, EOR = 1005 ) ( EqVel ( I ), I = 1, Param%IntM( 6, 2) * NDIM ) ;
  End If ;

  ! Reads down node numbers and equation numbers of nodes in which history of acceleration is required
  If ( Param%IntM( 6, 3) /= 0 ) Then ; ! NNAH
    Write(*,*)"Reading Node Numbers for acceleration ..." ;
    Read (Unit = UnFile, FMT = "(<Param%IntM( 6, 2)>(I19,2X),5x,<Param%IntM( 6, 2)>(I19,2X))", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Read, ERR = 1003, End = 1004, EOR = 1005 ) ( NAAN ( I ), I = 1, Param%IntM( 6, 3)), ( NAPN ( I ), I = 1, Param%IntM( 6, 3) ) ;
    Read (Unit = UnFile, FMT = "(<Param%IntM( 6, 2)*NDIM>(I19,2X))"             , ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Read, ERR = 1003, End = 1004, EOR = 1005 ) ( EqAcc ( I ), I = 1, Param%IntM( 6, 3) * NDIM ) ;
  End If ;

! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

! Material properties of each material
Write(*,*)"Reading Material Properties ..." ;
UnFile = UnInptMat ;
  DO I = 1, NMat ;
    Read (Unit = UnFile, FMT = "(<NPM>(E31.23E3,2X))", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Read, ERR = 1003, End = 1004, EOR = 1005 ) ( PMat ( I, J ), J = 1, NPM ) ; 
  End Do ;

! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

! Node numbering mapping
Write(*,*)"Reading Index Sets ..." ;
UnFile = UnInptApp ;
Read (Unit = UnFile, FMT = "(<NEQM_Mapping>(I19,2X))", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Read, ERR = 1003, End = 1004, EOR = 1005 ) ( App_Numbers   ( I ), I = 0, NEQM_Mapping - 1 ) ;
Read (Unit = UnFile, FMT = "(<NEQM_Mapping>(I19,2X))", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Read, ERR = 1003, End = 1004, EOR = 1005 ) ( PETSc_Numbers ( I ), I = 0, NEQM_Mapping - 1 ) ;
Read (Unit = UnFile, FMT = "(<NEQM>(I19,2X)        )", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Read, ERR = 1003, End = 1004, EOR = 1005 ) ( Indices       ( I ), I = 0, NEQM         - 1 ) ;

! Number of non-zero entris of PETSc objects   !!??  Check if they should start from 0 or 1 ;
UnFile = UnInptStiff ;
Write(*,*)"Reading Number of Non-Zeros of Stiffness ..." ;
!Read (Unit = UnFile, FMT = "(<NEqM_Mapping>(I19,2X))", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Read, ERR = 1003, End = 1004, EOR = 1005 ) ( D_NNZ_Stiff ( I ), I = 1, NEqM_Mapping ) ;  !!?? start if they should start from zero
!Read (Unit = UnFile, FMT = "(<NEqM_Mapping>(I19,2X))", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Read, ERR = 1003, End = 1004, EOR = 1005 ) ( O_NNZ_Stiff ( I ), I = 1, NEqM_Mapping ) ;

UnFile = UnInptDamp ;
Write(*,*)"Reading Number of Non-Zeros of Damping ..." ;
!Read (Unit = UnFile, FMT = "(<NEqM_Mapping>(I19,2X))", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Read, ERR = 1003, End = 1004, EOR = 1005 ) ( D_NNZ_Damp  ( I ), I = 1, NEqM_Mapping ) ;
!Read (Unit = UnFile, FMT = "(<NEqM_Mapping>(I19,2X))", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Read, ERR = 1003, End = 1004, EOR = 1005 ) ( O_NNZ_Damp  ( I ), I = 1, NEqM_Mapping ) ;

UnFile = UnInptMass ;
Write(*,*)"Reading Number of Non-Zeros of Mass ..." ;
!Read (Unit = UnFile, FMT = "(<NEqM_Mapping>(I19,2X))", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Read, ERR = 1003, End = 1004, EOR = 1005 ) ( D_NNZ_Mass  ( I ), I = 1, NEqM_Mapping ) ;
!Read (Unit = UnFile, FMT = "(<NEqM_Mapping>(I19,2X))", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Read, ERR = 1003, End = 1004, EOR = 1005 ) ( O_NNZ_Mass  ( I ), I = 1, NEqM_Mapping ) ;

Write(*    ,*) 'End Subroutine < Input_Arrays >' ;
Write(UnInf,*) 'End Subroutine < Input_Arrays >' ;
Return ;

! - ERROR IN READ STATEMENT -------------------------------------------------------------------------------------------------------------------------
1003    Write(*       , Fmt_READ1 ) UnFile, IO_READ ; Write( UnFile, Fmt_READ1 ) UnFile, IO_READ ;
        !#Call BEEP_FAIL ;
        Write(*, Fmt_FL) ;  Write(UnInf, Fmt_FL) ; Write(*, Fmt_End) ; Read(*,*) ;  Return ;

! - End-OF-FILE IN READ STATEMENT -------------------------------------------------------------------------------------------------------------------
1004    Write(*       , Fmt_READ2 ) UnFile, IO_READ ; Write( UnFile, Fmt_READ2 ) UnFile, IO_READ ;
        !#Call BEEP_FAIL ;
        Write(*, Fmt_FL) ;  Write(UnInf, Fmt_FL) ; Write(*, Fmt_End) ; Read(*,*) ;  Return ;

! - End-OF-FILE IN READ STATEMENT -------------------------------------------------------------------------------------------------------------------
1005    Write(*       , Fmt_READ3 ) UnFile, IO_READ ; Write( UnFile, Fmt_READ3 ) UnFile, IO_READ ;
        !#Call BEEP_FAIL ;
        Write(*, Fmt_FL) ;  Write(UnInf, Fmt_FL) ; Write(*, Fmt_End) ; Read(*,*) ;  Return ;


End Subroutine Input_Arrays ;


! ***************************************************************************************************************************************************
!                                                                                                                                                  **
! Start:        4 April 2013                                                                                                                       **
! Last Update:  4 April 2013                                                                                                                       **
! Description: THIS Subroutine COLLECTS ALL BASIC DATA FOR EACH ANALYSIS.                                                                          **
!                                                                                                                                                  **
! ***************************************************************************************************************************************************

Subroutine Input_Analysis (                                                                                                     &
NDim, Solver_Type, Int_Order,                                                                                                   & ! Integer (1) Variables
!NPM, NMat,                                                                                                                     & ! Integer (2) Variables
!                                                                                                                               & ! Integer (4) Variables
!                                                                                                                               & ! Integer (8) Variables
!                                                                                                                               & ! Real Variables
LoadC,     Step,                                                                                                                & ! Integer Arrays
InciWave,        BAcl,                                                                                                          & ! Real Arrays
!                                                                                                                               & ! Characters
Param                                                                                                                           & ! Type
) ;

Implicit None ;

! =========================== Global Variables ======================================================================================================

! - Integer Variables -------------------------------------------------------------------------------------------------------------------------------
Integer, Intent(In )    :: NDim ;

Integer, Intent(Out)    :: Solver_Type, Int_Order ;

! - Real Variables ----------------------------------------------------------------------------------------------------------------------------------
! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
! - Integer Arrays ----------------------------------------------------------------------------------------------------------------------------------
Integer , Intent(In), Dimension (:  )  :: LoadC ;
Integer , Allocatable, Intent(Out), Dimension (:  )  :: STEP ;


! - Real Arrays -------------------------------------------------------------------------------------------------------------------------------------
Real , Intent(Out), Dimension (:)    :: InciWave ;

Real , Allocatable, Intent(Out), Dimension (:,:)  :: BACL ;

! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
! - Character Variables -----------------------------------------------------------------------------------------------------------------------------
! - Logical Variables -------------------------------------------------------------------------------------------------------------------------------
! - Type DECLERATIONS -------------------------------------------------------------------------------------------------------------------------------
Type ( BasicParam )    :: Param ;                ! Holds basic parameters of each load case

! =========================== Local Variables =======================================================================================================

! - Integer Variables -------------------------------------------------------------------------------------------------------------------------------
Integer   :: I, J, K ;               ! Loop Index

Integer   :: IO_Read ;               ! Holds error of read statements
Integer   :: UnFile ;                ! Holds Unit of a file for error message

! - Real Variables ----------------------------------------------------------------------------------------------------------------------------------
! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
! - Integer Arrays ----------------------------------------------------------------------------------------------------------------------------------
! - Real Arrays -------------------------------------------------------------------------------------------------------------------------------------
! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
! - Character Variables -----------------------------------------------------------------------------------------------------------------------------
! - Logical Variables -------------------------------------------------------------------------------------------------------------------------------
! - Type DECLERATIONS -------------------------------------------------------------------------------------------------------------------------------

! =========================== Subroutine CODE =======================================================================================================

Write(*,*)
Write(*,*)"<<<<<<<< Input Analysis >>>>>>>>" ;


! Reading data from .dataAna
Allocate ( Param%IntL ( 7, 6),Param%RealL ( 7, 6)  ) ;   !?? modify the allocation numbers

UnFile = UnInptAna ;

! Basic Parameters
Write(*,*)"  Reading Basic Parameters for analysis..." ;

!Read  (UnFile, *) check ; Write(*,*)'1', check ! 1
!Read  (UnFile, *) check ; Write(*,*)'2', check ! 2
!Read  (UnFile, *) check ; Write(*,*)'3', check ! 3
!Read  (UnFile, *) check ; Write(*,*)'4', check ! 4
!Read  (UnFile, *) check ; Write(*,*)'5', check ! 5
!Read  (UnFile, *) check ; Write(*,*)'6', check ! 6
!Read  (UnFile, *) check ; Write(*,*)'7', check ! 7

Read  (UnFile, *) ; ! 1
Read  (UnFile, *) ; ! 2
Read  (UnFile, *) ; ! 3
Read  (UnFile, *) ; ! 4
Read  (UnFile, *) ; ! 5
Read  (UnFile, *) ; ! 6
Read  (UnFile, *) ; ! 7
Read  (UnFile, *) SOLVER_Type, Int_Order ;
Read  (UnFile, *) ; ! 10

! Reading Param Load
!Read (Unit = UnFile, FMT = "(<6*7>(I19,1X))", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Read, ERR = 1003, End = 1004, EOR = 1005 )  (( Param%IntL  ( I, J), I = 1, 7), J = 1, 6) ;
!Read (Unit = UnFile, FMT = "(<6*7>(E31.23E3,1X))", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Read, ERR = 1003, End = 1004, EOR = 1005 ) (( Param%RealL  ( I, J), I = 1, 7), J = 1, 6) ;

  If (LoadC (2) /= 0 ) Then ;  ! Pressure on elements
    Read  (UnInptAna, *);
    Read  (UnInptAna, *) Param%IntL( 2, 1), Param%IntL( 2, 2), Param%IntL( 2, 3),      Param%RealL( 2, 1), Param%RealL( 2, 2) ; ! PR_Type, LDIR, PRZERO,        HREF, PR
    Read  (UnInptAna, *);
  End If ;

  If ( LoadC (5) == 0 ) Then ;  ! Static Analysis
    Read  (UnInptAna, *) ;
    Read  (UnInptAna, *) Param%IntL( 5, 1 ) ; ! STStep
    Read  (UnInptAna, *) ;
  End If ;

  If ( LoadC (5) /= 0 ) Then ;  ! Dynamic Analysis
    Read  (UnInptAna, *) ;
    Read  (UnInptAna, *) ( Param%IntL( 5, I), I = 1, 2 ),Param%IntL( 5, 3), ( Param%IntL( 5, I), I = 4, 6 ),        ( Param%RealL( 5, I), I = 3, 4 ), Param%RealL( 5, 5), Param%RealL( 5, 6) ; ! NSTEP, NTST, NSParaview, NNodeDown, Solver(Implicit/Explicit 0/1), Analysis(Dynamics of structure/wave propagation 0/1),       DT, t0, TotalTime
    Read  (UnInptAna, *) ;

    Write(*,*)"Dynamic Loads" ;

      ! Loads
      If      ( Param%IntM( 5, 3) == 1 ) Then ;  ! Dynamic Pressure
        Read  (UnInptAna, *) ;
        Read  (UnInptAna, *) Param%IntL( 2, 1), Param%IntL( 2, 2), Param%IntL( 2, 3),     Param%RealL( 2, 1), Param%RealL( 2, 2) ; ! PR_Type, LDIR, PRZERO,        HREF, PR
        Read  (UnInptAna, *) ;

        Read  (UnInptAna, *) ;
        Read  (UnInptAna, *) Param%IntL( 7, 1) ; ! FUNC_Type  - SINE FUNCTION OR RICKER PULSE OR ... 
        Read  (UnInptAna, *) ;

        Write(*,*)"Dynamic pressure" ;

      Else If ( Param%IntM( 5, 3) == 2 ) Then ; ! Base Acceleration
        Read  (UnInptAna, *) ;
        Read  (UnInptAna, *) Param%IntL( 7, 1), Param%RealL( 7, 1) ; ! NASTEP, G   Number of base Accelerations to be read in the input file
        Read  (UnInptAna, *) ;

        Write(*,*)"Base Acceleration" ;

      Else If ( Param%IntM( 5, 3) == 3 ) Then ; ! Domain Reduction Method loading
        Read  (UnInptAna, *) ;
        Read  (UnInptAna, *) ( Param%IntL( 7, I), I = 3, 4 ), Param%RealL( 7, 1) ; ! Wave_Type, wave_func     AngleInci: angle of incident wave plane with the system coordinate ;
        Read  (UnInptAna, *) ;

        Write(*,*)"DRM" ;
      End If ;

      !  Implicit Solver
      If      ( Param%IntL( 5, 5) == 0 ) Then ;

        Read  (UnInptAna, *);
        Read  (UnInptAna, *) ( Param%RealL( 5, I), I = 1, 2 ) ; ! Delta, Gama        Newmark Coefficients
        Read  (UnInptAna, *);

      End If ;

  End If ;


!! Material properties of each material
!Write(*,*)"Reading Material Properties ..." ;
!  DO I = 1, NMat ;
!    Read (Unit = UnFile, FMT = "(<NPM>(E31.23E3,2X))", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Read, ERR = 1003, End = 1004, EOR = 1005 ) ( PMat ( I, J ), J = 1, NPM ) ; 
!  End Do ;


Read  (UnInptAna, *);
Read  (UnInptAna, *);

! Allocating required arrays
Allocate ( STEP  ( Param%IntL( 5, 2) ) ) ; !   Param%IntL( 5, 2) = NTST

! Required history of nodes
Write(*,*)"Reading Step Numbers ..." ;
Read  (UnInptAna, *);
Read  (UnInptAna, * ) ( STEP ( I ), I = 1, Param%IntL( 5, 2) ) ; 
Read  (UnInptAna, *);

    ! Dynamic loads
      If ( LoadC (5) /= 0 ) Then ;

        If ( Param%IntM( 5, 3 ) == 2 ) Then ;   ! Base acceleration

          ! Allocating required arrays
          Allocate ( BACL (Param%IntL( 7, 1 ), NDim) ) ; ! Param%IntL( 7, 1 ) = NASTEP

          Write(*,*)"Reading Base Acceleration ..." ;
            DO I = 1, Param%IntL( 7, 1) ; ! NAStep = Param%IntP( 7, 1)
              Read (UnInptAna, * ) ( BACL ( I, K ), K = 1, NDim ) ;
            End Do ;

        Else If ( Param%IntM( 5, 3 ) == 3 ) Then ;    ! Domain Reduction Method (DRM)

          Write(*,*)"Reading wave information ..." ;

          Read  (UnInptAna, *);
          Read  (UnInptAna, *);
          Read  (UnInptAna, * ) InciWave(1), InciWave(2), InciWave(3), InciWave(4), InciWave(5) ;       ! Theta, Omega, amplitude, alpha1, alpha2
          Read  (UnInptAna, *);


        End If ;

      End If ;

Write(*    ,*) 'End Subroutine < Input_Analysis >' ;
Write(UnInf,*) 'End Subroutine < Input_Analysis >' ;
Return ;

! - ERROR IN READ STATEMENT -------------------------------------------------------------------------------------------------------------------------
1003    Write(*     , Fmt_READ1 ) UnFile, IO_READ ; Write( UnFile, Fmt_READ1 ) UnFile, IO_READ ;
        !#Call BEEP_FAIL ;
        Write(*, Fmt_FL) ;  Write(UnInf, Fmt_FL) ; Write(*, Fmt_End) ; Read(*,*) ;  Return ;

! - End-OF-FILE IN READ STATEMENT -------------------------------------------------------------------------------------------------------------------
1004    Write(*     , Fmt_READ2 ) UnFile, IO_READ ; Write( UnFile, Fmt_READ2 ) UnFile, IO_READ ;
        !#Call BEEP_FAIL ;
        Write(*, Fmt_FL) ;  Write(UnInf, Fmt_FL) ; Write(*, Fmt_End) ; Read(*,*) ;  Return ;

! - End-OF-FILE IN READ STATEMENT -------------------------------------------------------------------------------------------------------------------
1005    Write(*     , Fmt_READ3 ) UnFile, IO_READ ; Write( UnFile, Fmt_READ3 ) UnFile, IO_READ ;
        !#Call BEEP_FAIL ;
        Write(*, Fmt_FL) ;  Write(UnInf, Fmt_FL) ; Write(*, Fmt_End) ; Read(*,*) ;  Return ;

End Subroutine Input_Analysis ;


! ***************************************************************************************************************************************************
!                                                                                                                                                  **
! Last Update:  10 July 2013                                                                                                                       **
! Description: THIS Subroutine reads heterogeneous material properties at each node.                                                               **
!                                                                                                                                                  **
! ***************************************************************************************************************************************************

Subroutine Input_Het_Mat (                                                                                                      &
!                                                                                                                               & ! Integer (1) Variables
!                                                                                                                               & ! Integer (2) Variables
!                                                                                                                               & ! Integer (4) Variables
NJ,                                                                                                                             & ! Integer (8) Variables
!                                                                                                                               & ! Real Variables
!                                                                                                                               &
!                                                                                                                               &
!                                                                                                                               & ! Integer Arrays
PMat_Lambda, PMat_Mu                                                                                                            & ! Real Arrays
!                                                                                                                               & ! Characters
!                                                                                                                               & ! Type
) ;

!Implicit None ;

#include "finclude/petscsys.h"
#include "finclude/petscmat.h"
#include "finclude/petscvec.h90"

PetscMPIInt    :: Rank ;                         ! Rank number

! =========================== Global Variables ======================================================================================================

! - Integer Arrays ------------------------------------------------------------------------------------------------------------------------------
Integer , Intent(In)    :: NJ;

! - Real Variables ----------------------------------------------------------------------------------------------------------------------------------
! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
! - Integer Arrays ----------------------------------------------------------------------------------------------------------------------------------
! - Real Arrays -------------------------------------------------------------------------------------------------------------------------------------
Real(8) ,     Intent(Out), Dimension (:)  :: PMat_Lambda, PMat_Mu ;

! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
! - Character Variables -----------------------------------------------------------------------------------------------------------------------------
! - Logical Variables -------------------------------------------------------------------------------------------------------------------------------
! - Type DECLERATIONS -------------------------------------------------------------------------------------------------------------------------------

! =========================== Local Variables =======================================================================================================
! - Integer Variables -------------------------------------------------------------------------------------------------------------------------------
Integer   :: IO_Read ;                ! Holds error of read statements.
Integer   :: UnFile ;                ! Holds Unit of a file for error message.
!Integer   :: INode;                   ! Loop counter for element node number.
!Integer   :: IDOF ;                   ! Loop counter on the degree of freedom.

!Integer   :: I1, I2 ;                 ! Loop indeces.
!Integer   :: NODE ;                   ! Some variable for holding node number.

Integer   :: I, J, K ;                ! Loop indeces.
!Integer   :: IEL ;                    ! Loop indeces over NEL - the number of elements.
Integer   :: IJ ;                     ! Loop index on NJ.

! - Real Variables ----------------------------------------------------------------------------------------------------------------------------------
! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
! - Integer Arrays ----------------------------------------------------------------------------------------------------------------------------------
! - Real Arrays -------------------------------------------------------------------------------------------------------------------------------------
! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
! - Character Variables -----------------------------------------------------------------------------------------------------------------------------
! - Logical Variables -------------------------------------------------------------------------------------------------------------------------------
! - Type DECLERATIONS -------------------------------------------------------------------------------------------------------------------------------

! =========================== Subroutine CODE =======================================================================================================

! Lambda
UnFile = UnInpt_Lambda ;
Write(*,*)"Reading Heterogeneous Lambda ..." ;
  DO K = 1, NJ ;
    Read (Unit = UnFile, FMT = "(E31.23E3)", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Read, ERR = 1003, End = 1004, EOR = 1005 ) PMat_Lambda ( K ) ;
  End Do ;

! Mu
UnFile = UnInpt_Mu ;
Write(*,*)"Reading Heterogeneous Mu ..." ;
  DO K = 1, NJ ;
    Read (Unit = UnFile, FMT = "(E31.23E3)", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Read, ERR = 1003, End = 1004, EOR = 1005 ) PMat_Mu ( K ) ;
  End Do ;

! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

Write(*    ,*) 'End Subroutine < Input_Het_Mat >' ;
Write(UnInf,*) 'End Subroutine < Input_Het_Mat >' ;
Return ;

! - ERROR IN READ STATEMENT -------------------------------------------------------------------------------------------------------------------------
1003    Write(*       , Fmt_READ1 ) UnFile, IO_READ ; Write( UnFile, Fmt_READ1 ) UnFile, IO_READ ;
        !#Call BEEP_FAIL ;
        Write(*, Fmt_FL) ;  Write(UnInf, Fmt_FL) ; Write(*, Fmt_End) ; Read(*,*) ;  Return ;

! - End-OF-FILE IN READ STATEMENT -------------------------------------------------------------------------------------------------------------------
1004    Write(*       , Fmt_READ2 ) UnFile, IO_READ ; Write( UnFile, Fmt_READ2 ) UnFile, IO_READ ;
        !#Call BEEP_FAIL ;
        Write(*, Fmt_FL) ;  Write(UnInf, Fmt_FL) ; Write(*, Fmt_End) ; Read(*,*) ;  Return ;

! - End-OF-FILE IN READ STATEMENT -------------------------------------------------------------------------------------------------------------------
1005    Write(*       , Fmt_READ3 ) UnFile, IO_READ ; Write( UnFile, Fmt_READ3 ) UnFile, IO_READ ;
        !#Call BEEP_FAIL ;
        Write(*, Fmt_FL) ;  Write(UnInf, Fmt_FL) ; Write(*, Fmt_End) ; Read(*,*) ;  Return ;


End Subroutine Input_Het_Mat ;


! ***************************************************************************************************************************************************
!                                                                                                                                                  **
! Last Update:  7 January 2014                                                                                                                     **
! Description: THIS Subroutine reads the basic data for Inversion Data Structure.                                                                  **
!                                                                                                                                                  **
! ***************************************************************************************************************************************************

Subroutine Input_Inv_DS_Basic (                                                                                                 &
!                                                                                                                               & ! Integer (1) Variables
!                                                                                                                               & ! Integer (2) Variables
!                                                                                                                               & ! Integer (4) Variables
NJ_Rank_IParts, NJ_Mapping, NStore_Mapping, NStore_Rank                                                                         & ! Integer (8) Variables
!                                                                                                                               & ! Real Variables
!                                                                                                                               &
!                                                                                                                               &
!                                                                                                                               & ! Integer Arrays
!                                                                                                                               & ! Real Arrays
!                                                                                                                               & ! Characters
!                                                                                                                               & ! Type
) ;

Implicit None ;

! =========================== Global Variables ======================================================================================================

! - Integer Arrays ------------------------------------------------------------------------------------------------------------------------------
Integer , Intent(Out)    :: NJ_Rank_IParts ; ! (=NJ)
Integer , Intent(Out)    :: NJ_Mapping ;
Integer , Intent(Out)    :: NStore_Mapping ;
Integer , Intent(Out)    :: NStore_Rank ;

! - Real Variables ----------------------------------------------------------------------------------------------------------------------------------
! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
! - Integer Arrays ----------------------------------------------------------------------------------------------------------------------------------
! - Real Arrays -------------------------------------------------------------------------------------------------------------------------------------
! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
! - Character Variables -----------------------------------------------------------------------------------------------------------------------------
! - Logical Variables -------------------------------------------------------------------------------------------------------------------------------
! - Type DECLERATIONS -------------------------------------------------------------------------------------------------------------------------------

! =========================== Local Variables =======================================================================================================
! - Integer Variables -------------------------------------------------------------------------------------------------------------------------------
Integer   :: IO_Read ;                ! Holds error of read statements.
Integer   :: UnFile ;                ! Holds Unit of a file for error message.

! - Real Variables ----------------------------------------------------------------------------------------------------------------------------------
! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
! - Integer Arrays ----------------------------------------------------------------------------------------------------------------------------------
! - Real Arrays -------------------------------------------------------------------------------------------------------------------------------------
! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
! - Character Variables -----------------------------------------------------------------------------------------------------------------------------
! - Logical Variables -------------------------------------------------------------------------------------------------------------------------------
! - Type DECLERATIONS -------------------------------------------------------------------------------------------------------------------------------

! =========================== Subroutine CODE =======================================================================================================

! Basic data
UnFile = Un_Inversion_DS ;
Write(*,*)"Reading Inversion Data Structure Basic Data ..." ;
Read (Unit = UnFile, FMT = "(4(I10,1X))", ADVANCE = 'YES', ASYNCHRONOUS = 'NO', IOSTAT = IO_Read, ERR = 1003, End = 1004, EOR = 1005 ) NJ_Rank_IParts, NJ_Mapping, NStore_Mapping, NStore_Rank ;


! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

Write(*    ,*) 'End Subroutine < Input_Inv_DS_Basic >' ;
Write(UnInf,*) 'End Subroutine < Input_Inv_DS_Basic >' ;
Return ;

! - ERROR IN READ STATEMENT -------------------------------------------------------------------------------------------------------------------------
1003    Write(*       , Fmt_READ1 ) UnFile, IO_READ ; Write( UnFile, Fmt_READ1 ) UnFile, IO_READ ;
        !#Call BEEP_FAIL ;
        Write(*, Fmt_FL) ;  Write(UnInf, Fmt_FL) ; Write(*, Fmt_End) ; Read(*,*) ;  Return ;

! - End-OF-FILE IN READ STATEMENT -------------------------------------------------------------------------------------------------------------------
1004    Write(*       , Fmt_READ2 ) UnFile, IO_READ ; Write( UnFile, Fmt_READ2 ) UnFile, IO_READ ;
        !#Call BEEP_FAIL ;
        Write(*, Fmt_FL) ;  Write(UnInf, Fmt_FL) ; Write(*, Fmt_End) ; Read(*,*) ;  Return ;


! - End-OF-FILE IN READ STATEMENT -------------------------------------------------------------------------------------------------------------------
1005    Write(*       , Fmt_READ3 ) UnFile, IO_READ ; Write( UnFile, Fmt_READ3 ) UnFile, IO_READ ;
        !#Call BEEP_FAIL ;
        Write(*, Fmt_FL) ;  Write(UnInf, Fmt_FL) ; Write(*, Fmt_End) ; Read(*,*) ;  Return ;


End Subroutine Input_Inv_DS_Basic ;


! ***************************************************************************************************************************************************
!                                                                                                                                                  **
! Last Update:  8 January 2014                                                                                                                     **
! Description: THIS Subroutine reads the array data for Inversion Data Structure.                                                                  **
!                                                                                                                                                  **
! ***************************************************************************************************************************************************

Subroutine Input_Inv_DS_Array (                                                                                                 &
!                                                                                                                               & ! Integer (1) Variables
!                                                                                                                               & ! Integer (2) Variables
!                                                                                                                               & ! Integer (4) Variables
NJ_Rank_IParts, NJ_Mapping, NStore_Mapping, NStore_Rank,                                                                        & ! Integer (8) Variables
!                                                                                                                               & ! Real Variables
!                                                                                                                               &
!                                                                                                                               &
idx_Mat_from, idx_Mat_to, idx_Mat_Extend_from, idx_Mat_Extend_to, U_Store_Numbers_Global, idx_u_from, idx_u_to                  & ! Integer Arrays
!                                                                                                                               & ! Real Arrays
!                                                                                                                               & ! Characters
!                                                                                                                               & ! Type
) ;

Implicit None ;

! =========================== Global Variables ======================================================================================================

! - Integer Arrays ------------------------------------------------------------------------------------------------------------------------------
Integer , Intent(In)    :: NJ_Rank_IParts ;
Integer , Intent(In)    :: NJ_Mapping ;
Integer , Intent(In)    :: NStore_Mapping ;
Integer , Intent(In)    :: NStore_Rank ;

! - Real Variables ----------------------------------------------------------------------------------------------------------------------------------
! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
! - Integer Arrays ----------------------------------------------------------------------------------------------------------------------------------
Integer , Intent(Out), Dimension ( : ) :: idx_Mat_from,        idx_Mat_to ;
Integer , Intent(Out), Dimension ( : ) :: idx_Mat_Extend_from, idx_Mat_Extend_to ;
Integer , Intent(Out), Dimension ( : ) :: U_Store_Numbers_Global ;
Integer , Intent(Out), Dimension ( : ) :: idx_u_from,          idx_u_to ;

! - Real Arrays -------------------------------------------------------------------------------------------------------------------------------------
! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
! - Character Variables -----------------------------------------------------------------------------------------------------------------------------
! - Logical Variables -------------------------------------------------------------------------------------------------------------------------------
! - Type DECLERATIONS -------------------------------------------------------------------------------------------------------------------------------

! =========================== Local Variables =======================================================================================================
! - Integer Variables -------------------------------------------------------------------------------------------------------------------------------
Integer   :: IO_Read ;                 ! Holds error of read statements.
Integer   :: UnFile ;                  ! Holds Unit of a file for error message.
Integer   :: I, J ;

! - Real Variables ----------------------------------------------------------------------------------------------------------------------------------
! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
! - Integer Arrays ----------------------------------------------------------------------------------------------------------------------------------
! - Real Arrays -------------------------------------------------------------------------------------------------------------------------------------
! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
! - Character Variables -----------------------------------------------------------------------------------------------------------------------------
! - Logical Variables -------------------------------------------------------------------------------------------------------------------------------
! - Type DECLERATIONS -------------------------------------------------------------------------------------------------------------------------------

! =========================== Subroutine CODE =======================================================================================================

! Index sets for inversion
UnFile = Un_Inversion_DS ;
Write(*,*)"Reading Inversion Data Structure Array Data ..." ;
Read (Unit = UnFile, FMT = "(<NJ_Rank_IParts>(I10,1X))", ADVANCE = 'YES', ASYNCHRONOUS = 'NO', IOSTAT = IO_Read, ERR = 1003, End = 1004, EOR = 1005 ) ( idx_Mat_from           ( I ), I = 1, NJ_Rank_IParts ) ;
Read (Unit = UnFile, FMT = "(<NJ_Mapping>(I10,1X))",     ADVANCE = 'YES', ASYNCHRONOUS = 'NO', IOSTAT = IO_Read, ERR = 1003, End = 1004, EOR = 1005 ) ( idx_Mat_Extend_from    ( I ), I = 1, NJ_Mapping ) ;
Read (Unit = UnFile, FMT = "(<NJ_Mapping>(I10,1X))",     ADVANCE = 'YES', ASYNCHRONOUS = 'NO', IOSTAT = IO_Read, ERR = 1003, End = 1004, EOR = 1005 ) ( idx_Mat_Extend_to      ( I ), I = 1, NJ_Mapping ) ;
Read (Unit = UnFile, FMT = "(<NStore_Mapping>(I10,1X))", ADVANCE = 'YES', ASYNCHRONOUS = 'NO', IOSTAT = IO_Read, ERR = 1003, End = 1004, EOR = 1005 ) ( U_Store_Numbers_Global ( I ), I = 1, NStore_Mapping ) ;
Read (Unit = UnFile, FMT = "(<NStore_Rank>(I10,1X))",    ADVANCE = 'YES', ASYNCHRONOUS = 'NO', IOSTAT = IO_Read, ERR = 1003, End = 1004, EOR = 1005 ) ( idx_u_from             ( I ), I = 1, NStore_Rank ) ;
Read (Unit = UnFile, FMT = "(<NStore_Rank>(I10,1X))",    ADVANCE = 'YES', ASYNCHRONOUS = 'NO', IOSTAT = IO_Read, ERR = 1003, End = 1004, EOR = 1005 ) ( idx_u_to               ( I ), I = 1, NStore_Rank ) ;

! Construction of idx_Mat_to:
ForAll ( I = 1:NJ_Rank_IParts ) idx_Mat_to (I) = I - 1;

! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

Write(*    ,*) 'End Subroutine < Input_Inv_DS_Array >' ;
Write(UnInf,*) 'End Subroutine < Input_Inv_DS_Array >' ;
Return ;

! - ERROR IN READ STATEMENT -------------------------------------------------------------------------------------------------------------------------
1003    Write(*       , Fmt_READ1 ) UnFile, IO_READ ; Write( UnFile, Fmt_READ1 ) UnFile, IO_READ ;
        !#Call BEEP_FAIL ;
        Write(*, Fmt_FL) ;  Write(UnInf, Fmt_FL) ; Write(*, Fmt_End) ; Read(*,*) ;  Return ;

! - End-OF-FILE IN READ STATEMENT -------------------------------------------------------------------------------------------------------------------
1004    Write(*       , Fmt_READ2 ) UnFile, IO_READ ; Write( UnFile, Fmt_READ2 ) UnFile, IO_READ ;
        !#Call BEEP_FAIL ;
        Write(*, Fmt_FL) ;  Write(UnInf, Fmt_FL) ; Write(*, Fmt_End) ; Read(*,*) ;  Return ;

! - End-OF-FILE IN READ STATEMENT -------------------------------------------------------------------------------------------------------------------
1005    Write(*       , Fmt_READ3 ) UnFile, IO_READ ; Write( UnFile, Fmt_READ3 ) UnFile, IO_READ ;
        !#Call BEEP_FAIL ;
        Write(*, Fmt_FL) ;  Write(UnInf, Fmt_FL) ; Write(*, Fmt_End) ; Read(*,*) ;  Return ;


End Subroutine Input_Inv_DS_Array ;


! ***************************************************************************************************************************************************
!                                                                                                                                                  **
! Last Update:  12 January 2014                                                                                                                    **
! Description: This Subroutine reads the measured data at select sensor locations to form misfit.                                                  **
!                                                                                                                                                  **
! ***************************************************************************************************************************************************

Subroutine Input_Inv_Measured_Data (                                                                                            &
!                                                                                                                               & ! Integer (1) Variables
!                                                                                                                               & ! Integer (2) Variables
!                                                                                                                               & ! Integer (4) Variables
NDim, NNDH,                                                                                                                     & ! Integer (8) Variables
!                                                                                                                               & ! Real Variables
!                                                                                                                               &
!                                                                                                                               &
!                                                                                                                               & ! Integer Arrays
Dis_meas                                                                                                                               & ! Real Arrays
!                                                                                                                               & ! Characters
!                                                                                                                               & ! Type
) ;

Implicit None ;

! =========================== Global Variables ======================================================================================================

! - Integer Arrays ------------------------------------------------------------------------------------------------------------------------------
Integer , Intent(In)    :: NDim ;
Integer , Intent(In)    :: NNDH ;

! - Real Variables ----------------------------------------------------------------------------------------------------------------------------------
! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
! - Integer Arrays ----------------------------------------------------------------------------------------------------------------------------------
Real(8) , Intent(Out), Dimension ( 0:NStep , NNDH * NDim ) :: Dis_meas ;

! - Real Arrays -------------------------------------------------------------------------------------------------------------------------------------
! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
! - Character Variables -----------------------------------------------------------------------------------------------------------------------------
! - Logical Variables -------------------------------------------------------------------------------------------------------------------------------
! - Type DECLERATIONS -------------------------------------------------------------------------------------------------------------------------------


! =========================== Local Variables =======================================================================================================
! - Integer Variables -------------------------------------------------------------------------------------------------------------------------------
Integer   :: IO_Read ;                 ! Holds error of read statements.
Integer   :: UnFile ;                  ! Holds Unit of a file for error message.
Integer   :: I, IStep ;

! - Real Variables ----------------------------------------------------------------------------------------------------------------------------------
Real(8)   :: Time ;
! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
! - Integer Arrays ----------------------------------------------------------------------------------------------------------------------------------
! - Real Arrays -------------------------------------------------------------------------------------------------------------------------------------
! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
! - Character Variables -----------------------------------------------------------------------------------------------------------------------------
Character ( Len = 250) :: temp_char ;   ! Temporary character (debugging)

! - Logical Variables -------------------------------------------------------------------------------------------------------------------------------
! - Type DECLERATIONS -------------------------------------------------------------------------------------------------------------------------------

! =========================== Subroutine CODE =======================================================================================================

UnFile = UN_Measured_HisD ;
If ( NNDH /= 0 ) Then
  Read  (UnFile, *) temp_char; ! 1
  Read  (UnFile, *) temp_char; ! 2
    Do IStep = 0, Nstep
      Read (Unit = UnFile, FMT = "(E14.6E3,2x, <NNDH>(<NDim>(E18.10E3,2x),3x) )", ADVANCE = 'YES', ASYNCHRONOUS = 'NO', IOSTAT = IO_Read, ERR = 1003, End = 1004, EOR = 1005) Time, ( Dis_meas ( IStep , I ), I = 1, NNDH * NDim ) ;
    End Do
End If

! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

Write(*    ,*) 'End Subroutine < Input_Inv_Measured_Data >' ;
Write(UnInf,*) 'End Subroutine < Input_Inv_Measured_Data >' ;
Return ;

! - ERROR IN READ STATEMENT -------------------------------------------------------------------------------------------------------------------------
1003    Write(*       , Fmt_READ1 ) UnFile, IO_READ ; Write( UnFile, Fmt_READ1 ) UnFile, IO_READ ;
        !#Call BEEP_FAIL ;
        Write(*, Fmt_FL) ;  Write(UnInf, Fmt_FL) ; Write(*, Fmt_End) ; Read(*,*) ;  Return ;

! - End-OF-FILE IN READ STATEMENT -------------------------------------------------------------------------------------------------------------------
1004    Write(*       , Fmt_READ2 ) UnFile, IO_READ ; Write( UnFile, Fmt_READ2 ) UnFile, IO_READ ;
        !#Call BEEP_FAIL ;
        Write(*, Fmt_FL) ;  Write(UnInf, Fmt_FL) ; Write(*, Fmt_End) ; Read(*,*) ;  Return ;

! - End-OF-FILE IN READ STATEMENT -------------------------------------------------------------------------------------------------------------------
1005    Write(*       , Fmt_READ3 ) UnFile, IO_READ ; Write( UnFile, Fmt_READ3 ) UnFile, IO_READ ;
        !#Call BEEP_FAIL ;
        Write(*, Fmt_FL) ;  Write(UnInf, Fmt_FL) ; Write(*, Fmt_End) ; Read(*,*) ;  Return ;


End Subroutine Input_Inv_Measured_Data ;


! ***************************************************************************************************************************************************
!                                                                                                                                                  **
! Last Update:  16 April 2014                                                                                                                      **
! Description: THIS Subroutine reads the number of material nodes per core for visualization.                                                      **
!                                                                                                                                                  **
! ***************************************************************************************************************************************************

Subroutine Input_Mat_Vis_Basic (                                                                                                &
!                                                                                                                               & ! Integer (1) Variables
!                                                                                                                               & ! Integer (2) Variables
!                                                                                                                               & ! Integer (4) Variables
NJ_Para,                                                                                                                        & ! Integer (8) Variables
!                                                                                                                               & ! Real Variables
!                                                                                                                               &
!                                                                                                                               &
!                                                                                                                               & ! Integer Arrays
!                                                                                                                               & ! Real Arrays
ModelName, IndexRank, IndexSize, Model_InDir                                                                                   & ! Characters
!                                                                                                                               & ! Type
) ;

Implicit None ;

! =========================== Global Variables ======================================================================================================

! - Integer Arrays ------------------------------------------------------------------------------------------------------------------------------

Integer , Intent(Out)    :: NJ_Para ;

! - Real Variables ----------------------------------------------------------------------------------------------------------------------------------
! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
! - Integer Arrays ----------------------------------------------------------------------------------------------------------------------------------
! - Real Arrays -------------------------------------------------------------------------------------------------------------------------------------
! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
! - Character Variables -----------------------------------------------------------------------------------------------------------------------------
Character (Kind = 1, Len = 30 ) :: ModelName ;   ! name of the model input file
Character (Kind = 1, Len = 20)  :: IndexRank ;   ! A variable for storing the Rank number in Character format for adding at the end of the input file Name.
Character (Kind = 1, Len = 20)  :: IndexSize ;   ! A variable for storing the Size number in Character format for adding at the end of the input file Name.
Character (Kind = 1, Len = 150) :: Model_InDir ; ! Directory of Model input file.

! - Logical Variables -------------------------------------------------------------------------------------------------------------------------------
! - Type DECLERATIONS -------------------------------------------------------------------------------------------------------------------------------

! =========================== Local Variables =======================================================================================================
! - Integer Variables -------------------------------------------------------------------------------------------------------------------------------
Integer, PARAMETER   :: Un_Mat_Vis     = 821  ;            ! the Unit number of Inversion Data Structure index sets
Integer   :: IO_File ;                ! Used for IOSTAT - Input Output Status - in the OPEN cammand.
Integer   :: IO_Read ;                ! Holds error of read statements.
Integer   :: UnFile ;                ! Holds Unit of a file for error message.

! - Real Variables ----------------------------------------------------------------------------------------------------------------------------------
! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
! - Integer Arrays ----------------------------------------------------------------------------------------------------------------------------------
! - Real Arrays -------------------------------------------------------------------------------------------------------------------------------------
! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
! - Character Variables -----------------------------------------------------------------------------------------------------------------------------
! - Logical Variables -------------------------------------------------------------------------------------------------------------------------------
! - Type DECLERATIONS -------------------------------------------------------------------------------------------------------------------------------

! =========================== Subroutine CODE =======================================================================================================

! Reading data from .Mat_Vis
    UnFile = Un_Mat_Vis ;
    Open ( Unit = UnFile, FILE = TRIM(ModelName)//'_'//Trim(AdjustL(IndexSize))//'_'//Trim(AdjustL(IndexRank))//'.Mat_Vis',  ERR =  1001, IOSTAT = IO_File, ACCESS = 'SEQUENTIAL', ACTION ='Read', ASYNCHRONOUS = 'NO', BLANK = 'NULL', BLOCKSIZE = 0, DEFAULTFILE = TRIM(Model_InDir), DisPOSE = 'KEEP', FORM = 'Formatted', POSITION = 'ASIS', STATUS ='Old') ;

    Write(*,*)"Reading Material Visualization Basic Data ..." ;
    Read (Unit = UnFile, FMT = "(1(I10,1X))", ADVANCE = 'YES', ASYNCHRONOUS = 'NO', IOSTAT = IO_Read, ERR = 1003, End = 1004, EOR = 1005 ) NJ_Para ;

    Close ( Unit = UnFile, STATUS = 'KEEP', ERR = 1002, IOSTAT = IO_File ) ;

! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

Write(*    ,*) 'End Subroutine < Input_Mat_Vis_Basic >' ;
Write(UnInf,*) 'End Subroutine < Input_Mat_Vis_Basic >' ;
Return ;

! =============================================== OPEN ERRORS =======================================================================================
1001  IF ( IO_File > 0 ) Then ;
        Write(*, Fmt_ERR1_OPEN) UnFile, IO_File  ;  Write(UnInf, Fmt_ERR1_OPEN) UnFile, IO_File  ; 
        Write(*, Fmt_FL) ; Write(UnInf, Fmt_FL) ;
        !#Call BEEP_FAIL ;
        Write(*, Fmt_End) ; Read(*,*) ;  Return ;
      Else If ( IO_File < 0 ) Then ;
        Write(*, Fmt_ERR1_OPEN) UnFile, IO_File  ; 
        Write(UnInf, Fmt_ERR1_OPEN) UnFile, IO_File  ;  Write(*, Fmt_FL) ; Write(UnInf, Fmt_FL) ;
        !#Call BEEP_FAIL ;
        Write(*, Fmt_End) ; Read(*,*) ;  Return ;
      End If ;

! =============================================== Close ERRORS ======================================================================================
1002  IF ( IO_File > 0 ) Then ;
        Write(*, Fmt_ERR1_Close) UnFile, IO_File  ;  Write(UnInf, Fmt_ERR1_Close) UnFile, IO_File  ; 
        Write(*, Fmt_FL) ; Write(UnInf, Fmt_FL) ;
        !Call BEEP_FAIL ;
        Write(*, Fmt_End) ; Read(*,*) ; ; Return ;
      End If ;

! - ERROR IN READ STATEMENT -------------------------------------------------------------------------------------------------------------------------
1003    Write(*       , Fmt_READ1 ) UnFile, IO_READ ; Write( UnFile, Fmt_READ1 ) UnFile, IO_READ ;
        !#Call BEEP_FAIL ;
        Write(*, Fmt_FL) ;  Write(UnInf, Fmt_FL) ; Write(*, Fmt_End) ; Read(*,*) ;  Return ;

! - End-OF-FILE IN READ STATEMENT -------------------------------------------------------------------------------------------------------------------
1004    Write(*       , Fmt_READ2 ) UnFile, IO_READ ; Write( UnFile, Fmt_READ2 ) UnFile, IO_READ ;
        !#Call BEEP_FAIL ;
        Write(*, Fmt_FL) ;  Write(UnInf, Fmt_FL) ; Write(*, Fmt_End) ; Read(*,*) ;  Return ;


! - End-OF-FILE IN READ STATEMENT -------------------------------------------------------------------------------------------------------------------
1005    Write(*       , Fmt_READ3 ) UnFile, IO_READ ; Write( UnFile, Fmt_READ3 ) UnFile, IO_READ ;
        !#Call BEEP_FAIL ;
        Write(*, Fmt_FL) ;  Write(UnInf, Fmt_FL) ; Write(*, Fmt_End) ; Read(*,*) ;  Return ;

End Subroutine Input_Mat_Vis_Basic ;


! ***************************************************************************************************************************************************
!                                                                                                                                                  **
! Last Update:  23 July 2014                                                                                                                       **
! Description: THIS Subroutine reads heterogeneous material properties at each node to be used if inversion is going to be resumed                 **
!                                                                                                                                                  **
! ***************************************************************************************************************************************************

Subroutine Input_Het_Mat_Resume (                                                                                               &
!                                                                                                                               & ! Integer (1) Variables
!                                                                                                                               & ! Integer (2) Variables
UnFile_L, UnFile_M,                                                                                                             & ! Integer (4) Variables
NJ,                                                                                                                             & ! Integer (8) Variables
!                                                                                                                               & ! Real Variables
!                                                                                                                               &
!                                                                                                                               &
!                                                                                                                               & ! Integer Arrays
PMat_Lambda, PMat_Mu                                                                                                            & ! Real Arrays
!                                                                                                                               & ! Characters
!                                                                                                                               & ! Type
) ;

!Implicit None ;

#include "finclude/petscsys.h"
#include "finclude/petscmat.h"
#include "finclude/petscvec.h90"

PetscMPIInt    :: Rank ;                         ! Rank number

! =========================== Global Variables ======================================================================================================

! - Integer Arrays ------------------------------------------------------------------------------------------------------------------------------
Integer , Intent(In)    :: NJ;

! - Real Variables ----------------------------------------------------------------------------------------------------------------------------------
! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
! - Integer Arrays ----------------------------------------------------------------------------------------------------------------------------------
! - Real Arrays -------------------------------------------------------------------------------------------------------------------------------------
Real(8) , Intent(Out), Dimension (:)  :: PMat_Lambda, PMat_Mu ;

! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
! - Character Variables -----------------------------------------------------------------------------------------------------------------------------
! - Logical Variables -------------------------------------------------------------------------------------------------------------------------------
! - Type DECLERATIONS -------------------------------------------------------------------------------------------------------------------------------

! =========================== Local Variables =======================================================================================================
! - Integer Variables -------------------------------------------------------------------------------------------------------------------------------
Integer   :: IO_Read ;                 ! Holds error of read statements.
Integer   :: UnFile ;                  ! Holds Unit of a file for error message.
Integer   :: UnFile_L, UnFile_M ;      ! Holds Unit file for Lambda and Mu corresponding to the last save in inversion (via Material_Save_Disk.F90).
!Integer   :: INode;                   ! Loop counter for element node number.
!Integer   :: IDOF ;                   ! Loop counter on the degree of freedom.

!Integer   :: I1, I2 ;                 ! Loop indeces.
!Integer   :: NODE ;                   ! Some variable for holding node number.

Integer   :: I, J, K ;                ! Loop indeces.
!Integer   :: IEL ;                    ! Loop indeces over NEL - the number of elements.
Integer   :: IJ ;                     ! Loop index on NJ.

! - Real Variables ----------------------------------------------------------------------------------------------------------------------------------
! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
! - Integer Arrays ----------------------------------------------------------------------------------------------------------------------------------
! - Real Arrays -------------------------------------------------------------------------------------------------------------------------------------
! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
! - Character Variables -----------------------------------------------------------------------------------------------------------------------------
! - Logical Variables -------------------------------------------------------------------------------------------------------------------------------
! - Type DECLERATIONS -------------------------------------------------------------------------------------------------------------------------------

! =========================== Subroutine CODE =======================================================================================================

! Lambda
UnFile = UnFile_L ;
Write(*,*)"Reading Heterogeneous Lambda ..." ;
  DO K = 1, NJ ;
    Read (Unit = UnFile, FMT = "(E31.23E3)", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Read, ERR = 1003, End = 1004, EOR = 1005 ) PMat_Lambda ( K ) ;
  End Do ;

! Mu
UnFile = UnFile_M ;
Write(*,*)"Reading Heterogeneous Mu ..." ;
  DO K = 1, NJ ;
    Read (Unit = UnFile, FMT = "(E31.23E3)", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Read, ERR = 1003, End = 1004, EOR = 1005 ) PMat_Mu ( K ) ;

  End Do ;

! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

Write(*    ,*) 'End Subroutine < Input_Het_Mat_Resume >' ;
Write(UnInf,*) 'End Subroutine < Input_Het_Mat_Resume >' ;
Return ;

! - ERROR IN READ STATEMENT -------------------------------------------------------------------------------------------------------------------------
1003    Write(*       , Fmt_READ1 ) UnFile, IO_READ ; Write( UnFile, Fmt_READ1 ) UnFile, IO_READ ;
        !#Call BEEP_FAIL ;
        Write(*, Fmt_FL) ;  Write(UnInf, Fmt_FL) ; Write(*, Fmt_End) ; Read(*,*) ;  Return ;

! - End-OF-FILE IN READ STATEMENT -------------------------------------------------------------------------------------------------------------------
1004    Write(*       , Fmt_READ2 ) UnFile, IO_READ ; Write( UnFile, Fmt_READ2 ) UnFile, IO_READ ;
        !#Call BEEP_FAIL ;
        Write(*, Fmt_FL) ;  Write(UnInf, Fmt_FL) ; Write(*, Fmt_End) ; Read(*,*) ;  Return ;

! - End-OF-FILE IN READ STATEMENT -------------------------------------------------------------------------------------------------------------------
1005    Write(*       , Fmt_READ3 ) UnFile, IO_READ ; Write( UnFile, Fmt_READ3 ) UnFile, IO_READ ;
        !#Call BEEP_FAIL ;
        Write(*, Fmt_FL) ;  Write(UnInf, Fmt_FL) ; Write(*, Fmt_End) ; Read(*,*) ;  Return ;


End Subroutine Input_Het_Mat_Resume ;


End Module Input_Subroutines ;




