
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                                                                                                                  ++
! Start:        29 June 2013                                                                                                                       ++
! Last Update:  29 June 2013                                                                                                                       ++
!                                                                                                                                                  ++
! Description: THIS Module Contains ALL Input Subroutines.                                                                                         ++
!                                                                                                                                                  ++
!              Subroutines                     ELEMENT NUMBER                                                                                      ++
!              ==============                  ==============                                                                                      ++
!              Input_BASIC                     COLLECTS ALL BASIC Variables.                                                                       ++
!              Input_Arrays                    COLLECTS ALL Arrays                                                                                 ++
!                                                                                                                                                  ++
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Module Input_Subroutines_HDF5 ;

Use HDF5 ;
Use Parameters ;

Implicit None ;

  Interface Input_HDF5 ;
    Module Procedure Input_Arrays_HDF5, Input_Analysis_HDF5 ;
  End Interface ;

Contains ;

! ***************************************************************************************************************************************************
!                                                                                                                                                  **
! Start:        29 June 2013                                                                                                                       **
! Last Update:  03 July 2013                                                                                                                       **
! Description: THIS Subroutine COLLECTS ALL Arrays in the HDF5 format                                                                              **
!                                                                                                                                                  **
! ***************************************************************************************************************************************************

Subroutine Input_Arrays_HDF5 (                                                                                                  &
NDOF, NDim, MaxNNode,                                                                                                           & ! Integer (1) Variables
NGroup, NPM, NMat,     NLCase,                                                                                                  & ! Integer (2) Variables
!                                                                                                                               & ! Integer (4) Variables
NEL, NJ, NEQM, NEQM_Mapping,                                                                                                    & ! Integer (8) Variables
!                                                                                                                               & ! Real Variables
LoadC,     IDBC,      MTEL, ELT,     ELGR,                                                                                      &
!D_NNZ_Stiff, O_NNZ_Stiff, D_NNZ_Damp, O_NNZ_Damp, D_NNZ_Mass, O_NNZ_Mass,     App_Numbers, PETSc_Numbers,     Indices,          &
App_Numbers, PETSc_Numbers,     Indices,                                                                                        &
NDAN , NVAN, NAAN, NDPN , NVPN, NAPN,     NoBndry_DRM, NoLayer_DRM, JLoad,     EqDis, EqVel, EqAcc,      INOD, ID,              & ! Integer Arrays
PMat, PBLD, XYZ, UDis, PLoad, PML_DIM,                                                                                          & ! Real Arrays
ModelName, IndexRank, IndexSize, Model_InDir,                                                                                   & ! Characters
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
!Integer (Kind=Shrt), Intent(Out), Dimension (:  )  :: D_NNZ_Stiff, O_NNZ_Stiff, D_NNZ_Damp, O_NNZ_Damp, D_NNZ_Mass, O_NNZ_Mass ;
Integer , Intent(Out), Dimension ( 0:NEQM_Mapping - 1  )  :: App_Numbers, PETSc_Numbers ;
Integer , Intent(Out), Dimension ( 0:NEQM - 1  )          :: Indices ;

Integer , Intent(Out), Dimension (:  )  :: NDAN , NVAN, NAAN, NDPN , NVPN, NAPN ;
Integer , Intent(Out), Dimension (:  )  :: NoBndry_DRM, NoLayer_DRM, JLoad ;
Integer , Intent(Out), Dimension (:  )  :: EqDis, EqVel, EqAcc ;
Integer , Intent(Out), Dimension (:,:)  :: INOD, ID ;
Integer , Intent(Out), Dimension (:,:)  :: IDBC ;

! - Real Arrays -------------------------------------------------------------------------------------------------------------------------------------
Real(8),     Intent(Out), Dimension (:,:)  :: PMat, PBLD, XYZ, UDis, PLoad, PML_DIM ;

! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
! - Character Variables -----------------------------------------------------------------------------------------------------------------------------
Character (Kind = 1, Len = 30 ) :: ModelName ;   ! name of the model input file
Character (Kind = 1, Len = 20)  :: IndexRank ;   ! A variable for storing the Rank number in Character format for adding at the end of the input file Name.
Character (Kind = 1, Len = 20)  :: IndexSize ;   ! A variable for storing the Size number in Character format for adding at the end of the input file Name.
Character (Kind = 1, Len = 150) :: Model_InDir ; ! Directory of Model input file.

! - Logical Variables -------------------------------------------------------------------------------------------------------------------------------
! - Type DECLERATIONS -------------------------------------------------------------------------------------------------------------------------------
Type ( BasicParam )    :: Param ;                ! Holds basic parameters of each load case

! =========================== Local Variables =======================================================================================================
PetscErrorCode :: ErrPTC ;                       ! Error 

! - Integer Variables -------------------------------------------------------------------------------------------------------------------------------
Integer   :: IO_File ;                ! Used for IOSTAT - Input Output Status - in the OPEN cammand.
Integer   :: IO_Read ;                ! Holds error of read statements.
Integer   :: UnFile ;                ! Holds Unit of a file for error message.
Integer   :: INode;                   ! Loop counter for element node number.
Integer   :: IDOF ;                   ! Loop counter on the degree of freedom.

Integer   :: I1, I2 ;                 ! Loop indeces.
Integer   :: NODE ;                   ! Some variable for holding node number.

Integer   :: I, J, K ;                ! Loop indeces.
Integer   :: IEL ;                    ! Loop indeces over NEL - the number of elements.
Integer   :: IJ ;                     ! Loop index on NJ.

! HDF5 variables
Integer(HID_T)       :: id_Geo ;                 ! File identifier for the Geometry in the main code

!Integer(HID_T)       :: dspace_id_XYZ ;          ! Dataspace identifier for coordinates
!Integer(HID_T)       :: dspace_id_Cnn ;          ! Dataspace identifier for connectivity
!Integer(HID_T)       :: dspace_id_Cnt ;          ! Dataspace identifier for constraints
!Integer(HID_T)       :: dspace_id_App ;          ! Dataspace identifier for Application Numbering
!Integer(HID_T)       :: dspace_id_PTC ;          ! Dataspace identifier for PETSc Numbering
!Integer(HID_T)       :: dspace_id_Ind ;          ! Dataspace identifier for Indices (Index Set)

Integer(HID_T)       :: dset_id_XYZ ;            ! Dataset identifier for coordinates
Integer(HID_T)       :: dset_id_Cnn ;            ! Dataset identifier for connectivity
Integer(HID_T)       :: dset_id_Cnt ;            ! Dataset identifier for constraints
Integer(HID_T)       :: dset_id_App ;            ! Dataset identifier for Application Numbering
Integer(HID_T)       :: dset_id_PTC ;            ! Dataset identifier for PETSc Numbering
Integer(HID_T)       :: dset_id_Ind ;            ! Dataset identifier for Indices (Index Set)

!Integer(HSIZE_T), DIMENSION(2) :: dims           ! Dataset dimensions
Integer(HSIZE_T), DIMENSION(2) :: data_dims

! - Real Variables ----------------------------------------------------------------------------------------------------------------------------------
! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
! - Integer Arrays ----------------------------------------------------------------------------------------------------------------------------------
Integer, Allocatable, Dimension(:,:)  :: dset_data_int ; ! Array to fill HDF5 file

! - Real Arrays -------------------------------------------------------------------------------------------------------------------------------------
! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
! - Character Variables -----------------------------------------------------------------------------------------------------------------------------
! - Logical Variables -------------------------------------------------------------------------------------------------------------------------------
! - Type DECLERATIONS -------------------------------------------------------------------------------------------------------------------------------

! =========================== Subroutine CODE =======================================================================================================

Write(*,*)
Write(*,*)"<<<<<<<< Input Arrays >>>>>>>>>" ;

! - Open Files --------------------------------------------------------------------------------------------------------------------------------------

! Material Properties (.Mat)
UnFile = UnInptMat ;
Open ( Unit = UnFile, FILE = TRIM(ModelName)//'_'//Trim(AdjustL(IndexSize))//'_'//Trim(AdjustL(IndexRank))//'.Mat', ERR =  1001, IOSTAT = IO_File, ACCESS = 'SEQUENTIAL', ACTION ='Read', ASYNCHRONOUS = 'NO', BLANK = 'NULL', BLOCKSIZE = 0, DEFAULTFILE = TRIM(Model_InDir), DisPOSE = 'KEEP', FORM = "Formatted", POSITION = 'ASIS', STATUS ='Old') ;

Call h5open_f(ErrPTC) ;

! Open and existing HDF5 file
Call h5fopen_f(TRIM(Model_InDir)//'/'//TRIM(ModelName)//'_'//'Geometry'//'_'//Trim(AdjustL(IndexSize))//'_'//Trim(AdjustL(IndexRank))//'.h5', H5F_ACC_RDONLY_F, id_Geo, ErrPTC)

! - Create the dataspaces
! Coordinate file for the main code(.XYZ)
!dims(1) = NJ ;
!dims(2) = NDim ;
!Call h5screate_simple_f(2, dims, dspace_id_XYZ, ErrPTC) ;

! Connectivity file for the main code (.Cnn)
!dims(1) = NEL ;
!dims(2) = MaxNNode + 3_Tiny ;
!Call h5screate_simple_f(2, dims, dspace_id_Cnn, ErrPTC) ;

! Coordinate file (.Cnt)
!dims(1) = NJ ;
!dims(2) = NDOF ;
!Call h5screate_simple_f(2, dims, dspace_id_Cnt, ErrPTC) ;

! Application Numbering file (.App)
!dims(1) = 1 ;
!dims(2) = NEqM_Mapping ;
!Call h5screate_simple_f(2, dims, dspace_id_App, ErrPTC) ;

! PETSc Numbering file (.PTC)
!dims(1) = 1 ;
!dims(2) = NEqM_Mapping ;
!Call h5screate_simple_f(2, dims, dspace_id_PTC, ErrPTC) ;

! Indices file (.Ind)
!dims(1) = 1 ;
!dims(2) = NEqM ;
!Call h5screate_simple_f(2, dims, dspace_id_Ind, ErrPTC) ;

! Create the dataset with default properties.
!Call h5dcreate_f(id_Geo, "XYZ",             H5T_NATIVE_DOUBLE,  dspace_id_XYZ, dset_id_XYZ, ErrPTC) ;
!Call h5dcreate_f(id_Geo, "XYZ",             H5T_NATIVE_REAL,    dspace_id_XYZ, dset_id_XYZ, ErrPTC) ;
!Call h5dcreate_f(id_Geo, "Connectivity",    H5T_NATIVE_INTEGER, dspace_id_Cnn, dset_id_Cnn, ErrPTC) ;
!Call h5dcreate_f(id_Geo, "Constraints",     H5T_NATIVE_INTEGER, dspace_id_Cnt, dset_id_Cnt, ErrPTC) ;   ! see type of integers for large numbers
!Call h5dcreate_f(id_Geo, "ApplicationNum",  H5T_NATIVE_INTEGER, dspace_id_App, dset_id_App, ErrPTC) ;   ! see type of integers for large numbers
!Call h5dcreate_f(id_Geo, "PETScNum",        H5T_NATIVE_INTEGER, dspace_id_App, dset_id_PTC, ErrPTC) ;   ! see type of integers for large numbers
!Call h5dcreate_f(id_Geo, "Indices",         H5T_NATIVE_INTEGER, dspace_id_Ind, dset_id_Ind, ErrPTC) ;   ! see type of integers for large numbers

Call h5dopen_f ( id_Geo, "XYZ",            dset_id_XYZ, ErrPTC ) ;
Call h5dopen_f ( id_Geo, "Connectivity",   dset_id_Cnn, ErrPTC ) ;
Call h5dopen_f ( id_Geo, "Constraints",    dset_id_Cnt, ErrPTC ) ;
Call h5dopen_f ( id_Geo, "ApplicationNum", dset_id_App, ErrPTC ) ;
Call h5dopen_f ( id_Geo, "PETScNum",       dset_id_PTC, ErrPTC ) ;
Call h5dopen_f ( id_Geo, "Indices",        dset_id_Ind, ErrPTC ) ;

NLCase = MaxVal ( LoadC ) ;

! - Geometry & Model -------------------------------------------------------------------------------------------------------------------------------

!! Number of elements of each group
!Write(*,*)"Reading Group numbers ..." ;
!  DO I = 1, NGroup ;
!    Read (Unit = UnFile, FMT = "(I19,2X,I5,2X)", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Read, ERR = 1003, End = 1004, EOR = 1005 ) NEG ( I, 1 ), NEG ( I, 2 ) ;
!  End Do ;

! Coordinates of nodes
Write(*,*)"Reading Coordinates ..." ;
!UnFile = UnInptXYZ ;
!  DO K = 1, NJ ;
!    Read (Unit = UnFile, FMT = "(<NDim>(E31.23E3,2X))", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Read, ERR = 1003, End = 1004, EOR = 1005 )( XYZ ( K, J ), J = 1, NDim ) ;
!  End Do ;

data_dims(1) = NJ ;
data_dims(2) = NDim ;
Call h5dread_f (dset_id_XYZ, H5T_NATIVE_DOUBLE, XYZ, data_dims, ErrPTC) ;

! Element connectivities
Write(*,*)"Reading Element Connectivities ..." ;
!UnFile = UnInptCnn ;
!  DO IEL = 1, NEL ;
!    Read (Unit = UnFile, FMT = "(<MaxNNode>(I19,2X),I5,2X,I5,2X,I5)", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Read, ERR = 1003, End = 1004, EOR = 1005 ) ( INOD ( INode, IEL ), INode = 1, MaxNNode ), MTEL ( IEL ), ELT ( IEL ), ELGR ( IEL ) ;
!  End Do ;
data_dims(1) = NEL ;
data_dims(2) = MaxNNode + 3 ;
Allocate ( dset_data_int ( NEl, MaxNNode + 3 ) ) ;

Call h5dread_f(dset_id_Cnn, H5T_NATIVE_INTEGER, dset_data_int, data_dims, ErrPTC) ;
ForAll (I = 1:MaxNNode, J = 1:NEl ) INOD ( I, J ) = dset_data_int ( J, I ) ;
MTEL ( : ) = dset_data_int ( :, MaxNNode+1 ) ;
ELT  ( : ) = dset_data_int ( :, MaxNNode+2 ) ;
ELGR ( : ) = dset_data_int ( :, MaxNNode+3 ) ;

DeAllocate ( dset_data_int ) ;

! Constraints
Write(*,*)"Reading Constraints ..." ;
!UnFile = UnInptCnt ;
!  DO IJ = 1, NJ ;
!    Read (Unit = UnFile, FMT = "(<NDOF>(I19,2X))", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Read, ERR = 1003, End = 1004, EOR = 1005 ) ( ID ( IJ, IDOF ), IDOF = 1, NDOF ) ;
!  End Do ;
Allocate ( dset_data_int ( NJ, NDOF ) ) ;

data_dims(1) = NJ ;
data_dims(2) = NDOF ;
Call h5dread_f( dset_id_Cnt, H5T_NATIVE_INTEGER, dset_data_int , data_dims, ErrPTC ) ;

ID (:,:) = dset_data_int (:,:) ;
DeAllocate ( dset_data_int ) ;


! PML teritory
UnFile = UnInptMdl ;
  If ( LoadC (5) /= 0 ) Then ;
    Write(*,*)"Reading PML Teritory ..." ;
      DO J = 1, 2
        Read (Unit = UnFile, FMT = "(<2*NDim>(E31.23E3,2X))", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Read, ERR = 1003, End = 1004 ) (PML_DIM ( I, J ), I = 1, 2 * NDim ) ; 
      End Do ;
  End If ;

! Body force load
  If ( LoadC (1) /= 0 ) Then ;
    Write(*,*)"Reading Body Forces ..." ;
      DO I = 1, Param%IntM( 1, 1) ; ! Param%IntP( 1, 1)=NBLD: Number of body force load.   Param%IntP( 1, 2) = NPBL
        Read (Unit = UnFile, FMT = "(<Param%IntM( 1, 2)>(E31.23E3,2x))", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Read, ERR = 1003, End = 1004 ) ( PBLD ( I, J ) , J = 1, Param%IntM( 1, 2) ) ;  ! Param%IntM( 1, 2) = NDIm
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
        Read (Unit = UnFile, FMT = "((I19,2X),<2*NDim>(I3,2X))", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Read, ERR = 1003, End = 1004 ) IDBC ( J, 1 ), ( IDBC ( J, I ), I = 2, 2 * NDim+1 ) ;
      End Do ;
  End If ;

  ! Joint load of every case - concentrated loads
  IF ( LoadC ( 3 ) /= 0 ) Then ;
    Write(*,*)"Reading Joint Loads ..." ;
      DO I = 1, Param%IntM( 3, 1) ; ! Param%IntP( 3, 1) = NLN
        Read (Unit = UnFile, FMT = "(I19,2X,<NDim>(E31.23E3,2X))", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Read, ERR = 1003, End = 1004 ) JLoad ( I ), ( PLoad( J, I ), J = 1, NDim ) ;
      End Do
  End If ;

! Supports Displacements
  IF ( LoadC ( 4 ) /= 0 ) Then ;
    Write(*,*)"Reading Support Displacements ..." ;
      DO I = 1, Param%IntM( 4, 1 ) ; ! Param%IntM( 4, 1 ) = NSND
        Read (Unit = UnFile, FMT = "(<NDOF+1>(E31.23E3,2X))", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Read, ERR = 1003, End = 1004 ) UDis( NDOF + 1, I ), ( UDis ( K, I ), K = 1, NDOF ) ;
      End Do ;
  End If ;

Write(*,*)"Reading Dynamic Loads ..." ;

! Dynamic loads
  If ( LoadC (5) /= 0 ) Then ;

    If ( Param%IntM( 5, 3 ) == 1 ) Then ;    ! RICKER PULSE or sine function for dynamic pressure
      Write(*,*)"Reading IDBC ..." ;
        DO J = 1, Param%IntM( 2, 4 ) ; ! NIDBC
          Read (Unit = UnFile, FMT = "((I19,2X),<2*NDim>(I3,2X))", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Read, ERR = 1003, End = 1004 ) IDBC ( J, 1 ), ( IDBC ( J, I ), I = 2, 2 * NDim+1 ) ;
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
      Read (Unit = UnFile, FMT = "(<Param%IntM( 7, 1)>(I19,1X))", ADVANCE = 'Yes', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Read, ERR = 1003, End = 1004 ) ( NoBndry_DRM ( I ), I = 1, Param%IntM( 7, 1) ) ; ! NNBndry_DRM

      ! Reading the node numbers on the DRM neighbor
      Write(*,*)"Reading DRM neighbor nodes ..."
      !Read (Unit = UnFile, FMT = "(<Param%IntM( 7, 2)>(I19,1X))", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Read, ERR = 1003, End = 1004, EOR = 1005 ) ( NoLayer_DRM ( I ), I = 1, Param%IntM( 7, 2)) ;  ! NNLayer_DRM
        Do I = 1, Param%IntM( 7, 2)
          Read (Unit = UnFile, FMT = "((I19,1X))", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Read, ERR = 1003, End = 1004 ) NoLayer_DRM ( I ) ;  ! NNLayer_DRM
        End Do ;
    End If ;

  End If ;
  ! Reads down node numbers and equation numbers of nodes in which history of Displacement is required
  If ( Param%IntM( 6, 1) /= 0 ) Then ; ! NNDH
    Write(*,*)"Reading Node Numbers for displacements ..." ;
    Read (Unit = UnFile, FMT = "(<Param%IntM( 6, 1)>(I19,2X),5x,<Param%IntM( 6, 1)>(I19,2X))", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Read, ERR = 1003, End = 1004 ) ( NDAN ( I ), I = 1, Param%IntM( 6, 1) ), ( NDPN ( I ), I = 1, Param%IntM( 6, 1) ) ;
    Read (Unit = UnFile, FMT = "(<Param%IntM( 6, 1)*NDIM>(I19,2X))"             , ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Read, ERR = 1003, End = 1004 ) ( EqDis ( I ), I = 1, Param%IntM( 6, 1)* NDIM ) ;
  End If ;

  ! Reads down node numbers and equation numbers of nodes in which history of velocity is required
  If ( Param%IntM( 6, 2) /= 0 ) Then ; ! NNVH
    Write(*,*)"Reading Node Numbers for velocity ..." ;
    Read (Unit = UnFile, FMT = "(<Param%IntM( 6, 2)>(I19,2X),5x,<Param%IntM( 6, 2)>(I19,2X))", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Read, ERR = 1003, End = 1004 ) ( NVAN ( I ), I = 1, Param%IntM( 6, 2) ), ( NVPN ( I ), I = 1, Param%IntM( 6, 2) ) ;
    Read (Unit = UnFile, FMT = "(<Param%IntM( 6, 2)*NDIM>(I19,2X))"             , ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Read, ERR = 1003, End = 1004) ( EqVel ( I ), I = 1, Param%IntM( 6, 2) * NDIM ) ;
  End If ;

  ! Reads down node numbers and equation numbers of nodes in which history of acceleration is required
  If ( Param%IntM( 6, 3) /= 0 ) Then ; ! NNAH
    Write(*,*)"Reading Node Numbers for acceleration ..." ;
    Read (Unit = UnFile, FMT = "(<Param%IntM( 6, 2)>(I19,2X),5x,<Param%IntM( 6, 2)>(I19,2X))", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Read, ERR = 1003, End = 1004 ) ( NAAN ( I ), I = 1, Param%IntM( 6, 3)), ( NAPN ( I ), I = 1, Param%IntM( 6, 3) ) ;
    Read (Unit = UnFile, FMT = "(<Param%IntM( 6, 2)*NDIM>(I19,2X))"             , ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Read, ERR = 1003, End = 1004 ) ( EqAcc ( I ), I = 1, Param%IntM( 6, 3) * NDIM ) ;
  End If ;

! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

! Material properties of each material
Write(*,*)"Reading Material Properties ..." ;
UnFile = UnInptMat ;
  DO I = 1, NMat ;
    Read (Unit = UnFile, FMT = "(<NPM>(E31.23E3,2X))", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Read, ERR = 1003, End = 1004 ) ( PMat ( I, J ), J = 1, NPM ) ; 
  End Do ;

! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

! Node numbering mapping
Write(*,*)"Reading Index Sets ..." ;
!UnFile = UnInptApp ;
!Read (Unit = UnFile, FMT = "(<NEQM_Mapping>(I19,2X))", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Read, ERR = 1003, End = 1004, EOR = 1005 ) ( App_Numbers   ( I ), I = 0, NEQM_Mapping - 1_Lng ) ;
!Read (Unit = UnFile, FMT = "(<NEQM_Mapping>(I19,2X))", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Read, ERR = 1003, End = 1004, EOR = 1005 ) ( PETSc_Numbers ( I ), I = 0, NEQM_Mapping - 1_Lng ) ;
!Read (Unit = UnFile, FMT = "(<NEQM>(I19,2X)        )", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Read, ERR = 1003, End = 1004, EOR = 1005 ) ( Indices       ( I ), I = 0, NEQM         - 1_Lng ) ;

! Application Numbering 
data_dims(1) = 1 ;
data_dims(2) = NEqM_Mapping;
Call h5dread_f( dset_id_App, H5T_NATIVE_INTEGER, App_Numbers, data_dims, ErrPTC ) ;

! PETSc Numbering 
data_dims(1) = 1 ;
data_dims(2) = NEqM_Mapping;
Call h5dread_f( dset_id_PTC, H5T_NATIVE_INTEGER, PETSc_Numbers, data_dims, ErrPTC ) ;

! Index set (Indices)
data_dims(1) = 1 ;
data_dims(2) = NEqM ;
Call h5dread_f( dset_id_Ind, H5T_NATIVE_INTEGER, Indices, data_dims, ErrPTC ) ;

!! Number of non-zero entris of PETSc objects   !!??  Check if they should start from 0 or 1 ;
!UnFile = UnInptStiff ;
!Write(*,*)"Reading Number of Non-Zeros of Stiffness ..." ;
!Read (Unit = UnFile, FMT = "(<NEqM_Mapping>(I19,2X))", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Read, ERR = 1003, End = 1004, EOR = 1005 ) ( D_NNZ_Stiff ( I ), I = 1_Lng, NEqM_Mapping ) ;  !!?? start if they should start from zero
!Read (Unit = UnFile, FMT = "(<NEqM_Mapping>(I19,2X))", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Read, ERR = 1003, End = 1004, EOR = 1005 ) ( O_NNZ_Stiff ( I ), I = 1_Lng, NEqM_Mapping ) ;
!
!UnFile = UnInptDamp ;
!Write(*,*)"Reading Number of Non-Zeros of Damping ..." ;
!Read (Unit = UnFile, FMT = "(<NEqM_Mapping>(I19,2X))", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Read, ERR = 1003, End = 1004, EOR = 1005 ) ( D_NNZ_Damp  ( I ), I = 1_Lng, NEqM_Mapping ) ;
!Read (Unit = UnFile, FMT = "(<NEqM_Mapping>(I19,2X))", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Read, ERR = 1003, End = 1004, EOR = 1005 ) ( O_NNZ_Damp  ( I ), I = 1_Lng, NEqM_Mapping ) ;
!
!UnFile = UnInptMass ;
!Write(*,*)"Reading Number of Non-Zeros of Mass ..." ;
!Read (Unit = UnFile, FMT = "(<NEqM_Mapping>(I19,2X))", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Read, ERR = 1003, End = 1004, EOR = 1005 ) ( D_NNZ_Mass  ( I ), I = 1_Lng, NEqM_Mapping ) ;
!Read (Unit = UnFile, FMT = "(<NEqM_Mapping>(I19,2X))", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Read, ERR = 1003, End = 1004, EOR = 1005 ) ( O_NNZ_Mass  ( I ), I = 1_Lng, NEqM_Mapping ) ;

! - Close Files -------------------------------------------------------------------------------------------------------------------------------------

! Material Properties (.Mat)
UnFile = UnInptMat ;
Close ( Unit = UnFile, STATUS = 'KEEP', ERR = 1002, IOSTAT = IO_File ) ;

! Close the dataset.
Call h5dclose_f( dset_id_XYZ, ErrPTC)
Call h5dclose_f( dset_id_Cnn, ErrPTC)
Call h5dclose_f( dset_id_Cnt, ErrPTC)
Call h5dclose_f( dset_id_App, ErrPTC)
Call h5dclose_f( dset_id_PTC, ErrPTC)
Call h5dclose_f( dset_id_Ind, ErrPTC)

! Close the HDF5 file.
CALL h5fclose_f( id_Geo, ErrPTC)

! Shut down HDF5.
CALL h5close_f(ErrPTC)

Write(*    ,*) 'End Subroutine < Input_Arrays_HDF5 >' ;
Write(UnInf,*) 'End Subroutine < Input_Arrays_HDF5 >' ;
Return ;

! =============================================== OPEN ERRORS =======================================================================================
1001  IF ( IO_File > 0 ) Then ;
        Write(*, Fmt_ERR1_OPEN) UnFile, IO_File  ;  Write(UnInf, Fmt_ERR1_OPEN) UnFile, IO_File  ; 
        Write(*, Fmt_FL) ; Write(UnInf, Fmt_FL) ;
        !#Call BEEP_FAIL ;
        Write(*, Fmt_End) ; Read(*,*) ;  STOP ;
      Else If ( IO_File < 0 ) Then ;
        Write(*, Fmt_ERR1_OPEN) UnFile, IO_File ; 
        Write(UnInf, Fmt_ERR1_OPEN) UnFile, IO_File  ;  Write(*, Fmt_FL) ; Write(UnInf, Fmt_FL) ;
        !#Call BEEP_FAIL ;
        Write(*, Fmt_End) ; Read(*,*) ;  STOP ;
      End If ;


! =============================================== Close ERRORS ======================================================================================
1002  IF ( IO_File > 0 ) Then ;
        Write(*, Fmt_ERR1_Close) UnFile, IO_File  ;  Write(UnInf, Fmt_ERR1_Close) UnFile, IO_File  ; 
        Write(*, Fmt_FL) ; Write(UnInf, Fmt_FL) ;
        !Call BEEP_FAIL ;
        Write(*, Fmt_End) ; Read(*,*) ; ; STOP ;
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


End Subroutine Input_Arrays_HDF5 ;


! ***************************************************************************************************************************************************
!                                                                                                                                                  **
! Start:        4 April 2013                                                                                                                       **
! Last Update:  4 April 2013                                                                                                                       **
! Description: THIS Subroutine COLLECTS ALL BASIC DATA FOR EACH ANALYSIS.                                                                          **
!                                                                                                                                                  **
! ***************************************************************************************************************************************************

Subroutine Input_Analysis_HDF5 (                                                                                                &
NDim, Solver_Type, Int_Order,                                                                                                   & ! Integer (1) Variables
!NPM, NMat,                                                                                                                     & ! Integer (2) Variables
!                                                                                                                               & ! Integer (4) Variables
NJ,                                                                                                                             & ! Integer (8) Variables
!                                                                                                                               & ! Real Variables
LoadC,     Step, idx_from,     ID_Para,                                                                                         & ! Integer Arrays
InciWave,        BAcl,                                                                                                          & ! Real Arrays
ModelName, IndexRank, IndexSize, Model_InDir,                                                                                   & ! Characters
Param                                                                                                                           & ! Type
) ;

Use HDF5 ;

Implicit None ;

! =========================== Global Variables ======================================================================================================
! - Integer Variables -------------------------------------------------------------------------------------------------------------------------------
Integer , Intent(In )    :: NDim ;

Integer , Intent(Out)    :: Solver_Type, Int_Order ;

Integer , Intent(In )    :: NJ ;

! - Real Variables ----------------------------------------------------------------------------------------------------------------------------------
! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
! - Integer Arrays ----------------------------------------------------------------------------------------------------------------------------------
Integer , Intent(In), Dimension (:  )  :: LoadC ;
Integer , Allocatable, Intent(Out), Dimension (:  )  :: STEP ;
Integer , Intent(Out), Dimension (:  )  :: idx_from ;
Integer , Intent(Out), Dimension (:,:)  :: ID_Para ;

! - Real Arrays -------------------------------------------------------------------------------------------------------------------------------------
Real(8) , Intent(Out), Dimension (:)    :: InciWave ;

Real(8) , Allocatable, Intent(Out), Dimension (:,:)  :: BACL ;

! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
! - Character Variables -----------------------------------------------------------------------------------------------------------------------------
Character (Kind = 1, Len = 30 ) :: ModelName ;   ! name of the model input file
Character (Kind = 1, Len = 20)  :: IndexRank ;   ! A variable for storing the Rank number in Character format for adding at the end of the input file Name.
Character (Kind = 1, Len = 20)  :: IndexSize ;   ! A variable for storing the Size number in Character format for adding at the end of the input file Name.
Character (Kind = 1, Len = 150) :: Model_InDir ; ! Directory of Model input file.

! - Logical Variables -------------------------------------------------------------------------------------------------------------------------------
! - Type DECLERATIONS -------------------------------------------------------------------------------------------------------------------------------
Type ( BasicParam )    :: Param ;                ! Holds basic parameters of each load case

! =========================== Local Variables =======================================================================================================
PetscErrorCode :: ErrPTC ;                       ! Error 

! - Integer Variables -------------------------------------------------------------------------------------------------------------------------------
Integer  :: I, J, K ;               ! Loop Index

!Integer (Kind=Smll)  :: IO_File ;                ! Used for IOSTAT - Input Output Status - in the OPEN cammand.
!Integer (Kind=Smll)  :: IO_Read ;               ! Holds error of read statements
Integer   :: UnFile ;                ! Holds Unit of a file for error message

! HDF5 variables
Integer(HID_T)       :: id_Data_PV ;              ! File identifier for the Geometry in Paraview

!Integer(HID_T)       :: dspace_id_Node_PV ;      ! Dataspace identifier for equation number of nodes on the rank for Paraview
!Integer(HID_T)       :: dspace_id_IS_PV ;        ! Dataspace identifier for the index set, Paraview purpose only

Integer(HID_T)       :: dset_id_Node_PV ;        ! Dataset identifier for global PETSc equation number of nodes on the rank for Paraview
Integer(HID_T)       :: dset_id_IS_PV ;          ! Dataset identifier for the index set, Paraview purposes only

!Integer(HSIZE_T), DIMENSION(2) :: dims           ! Dataset dimensions
Integer(HSIZE_T), DIMENSION(2) :: data_dims      ! Dataset dimensions

! - Real Variables ----------------------------------------------------------------------------------------------------------------------------------
! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
! - Integer Arrays ----------------------------------------------------------------------------------------------------------------------------------
! - Real Arrays -------------------------------------------------------------------------------------------------------------------------------------
! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
! - Character Variables -----------------------------------------------------------------------------------------------------------------------------
Character (Kind = 1, Len = 30 ) :: helpme ;     

! - Logical Variables -------------------------------------------------------------------------------------------------------------------------------
! - Type DECLERATIONS -------------------------------------------------------------------------------------------------------------------------------

! =========================== Subroutine CODE =======================================================================================================

Write(*,*)
Write(*,*)"<<<<<<<< Input Analysis HDF5 >>>>>>>>" ;

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
    Read  (UnInptAna, *) ( Param%IntL( 5, I), I = 1, 2 ),Param%IntL( 5, 3), ( Param%IntL( 5, I), I = 4, 6 ),        ( Param%RealL( 5, I), I = 3, 4 ), Param%RealL( 5, 5), Param%RealL( 5, 6) ; ! NSTEP, NTST, NSParaview, NNodeDown, Solver(Implicit/Explicit 0/1), Analysis(Dynamics of structure/wave propagation 0/1),       DT, t0, TotalTime, e_max
    Read  (UnInptAna, *) ;

    ! Reading data for visualization (Paraview) in the HDF5 format
!      If ( Param%IntL( 5, 3) < Param%IntL( 5, 1) ) Then ;    ! Param%IntL( 5, 3) = NSParaview -- Param%IntL( 5, 1 ) = NStep; if NSParaview< NStep, then it generates input files for Paraview, otherwise there would be no visualization

        Call h5open_f(ErrPTC) ;
        Call h5fopen_f( TRIM(Model_InDir)//'/'//TRIM(ModelName)//'_'//'Data_PV'//'_'//Trim(AdjustL(IndexSize))//'_'//Trim(AdjustL(IndexRank))//'.h5', H5F_ACC_RDONLY_F, id_Data_PV, ErrPTC) ;   ! Geometry file for input files

        ! Equation number of nodes on this rank used to extract data of each node for Paraview
        Call h5dopen_f  ( id_Data_PV, "EquationNum",  dset_id_Node_PV, ErrPTC ) ;
        Call h5dopen_f  ( id_Data_PV, "IndexSet",     dset_id_IS_PV,   ErrPTC ) ;

        data_dims(1) = Param%IntM ( 5, 2) ;
        data_dims(2) = NDim ;
        Call h5dread_f (dset_id_Node_PV, H5T_NATIVE_INTEGER, ID_Para, data_dims, ErrPTC) ;

        data_dims(1) = Param%IntM ( 5, 1) ; ! NE_Para ; ;
        data_dims(2) = 1 ;
        Call h5dread_f(dset_id_IS_PV, H5T_NATIVE_INTEGER,  idx_from, data_dims, ErrPTC) ;

        Call h5dclose_f( dset_id_Node_PV, ErrPTC) ;
        Call h5dclose_f( dset_id_IS_PV,   ErrPTC) ;

        ! Close the HDF5 file.
        CALL h5fclose_f( id_Data_PV, ErrPTC) ;

        ! Shut down HDF5.
        CALL h5close_f(ErrPTC) ;

!      End If ;

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

Write(*    ,*) 'End Subroutine < Input_Analysis_HDF5 >' ;
Write(UnInf,*) 'End Subroutine < Input_Analysis_HDF5 >' ;

Return ;

!! =============================================== OPEN ERRORS =======================================================================================
!1001  IF ( IO_File > 0 ) Then ;
!        Write(*, Fmt_ERR1_OPEN) UnFile, IO_File  ;  Write(UnInf, Fmt_ERR1_OPEN) UnFile, IO_File  ; 
!        Write(*, Fmt_FL) ; Write(UnInf, Fmt_FL) ;
!        !#Call BEEP_FAIL ;
!        Write(*, Fmt_End) ; Read(*,*) ;  STOP ;
!      Else If ( IO_File < 0 ) Then ;
!        Write(*, Fmt_ERR1_OPEN) UnFile, IO_File  ; 
!        Write(UnInf, Fmt_ERR1_OPEN) UnFile, IO_File  ;  Write(*, Fmt_FL) ; Write(UnInf, Fmt_FL) ;
!        !#Call BEEP_FAIL ;
!        Write(*, Fmt_End) ; Read(*,*) ;  STOP ;
!      End If ;
!
!
!! =============================================== Close ERRORS ======================================================================================
!1002  IF ( IO_File > 0 ) Then ;
!        Write(*, Fmt_ERR1_Close) UnFile, IO_File  ;  Write(UnInf, Fmt_ERR1_Close) UnFile, IO_File  ; 
!        Write(*, Fmt_FL) ; Write(UnInf, Fmt_FL) ;
!        !Call BEEP_FAIL ;
!        Write(*, Fmt_End) ; Read(*,*) ; ; STOP ;
!      End If ;
!
!! - ERROR IN READ STATEMENT -------------------------------------------------------------------------------------------------------------------------
!1003    Write(*     , Fmt_READ1 ) UnFile, IO_READ ; Write( UnFile, Fmt_READ1 ) UnFile, IO_READ ;
!        !#Call BEEP_FAIL ;
!        Write(*, Fmt_FL) ;  Write(UnInf, Fmt_FL) ; Write(*, Fmt_End) ; Read(*,*) ;  Return ;
!
!! - End-OF-FILE IN READ STATEMENT -------------------------------------------------------------------------------------------------------------------
!1004    Write(*     , Fmt_READ2 ) UnFile, IO_READ ; Write( UnFile, Fmt_READ2 ) UnFile, IO_READ ;
!        !#Call BEEP_FAIL ;
!        Write(*, Fmt_FL) ;  Write(UnInf, Fmt_FL) ; Write(*, Fmt_End) ; Read(*,*) ;  Return ;
!
!! - End-OF-FILE IN READ STATEMENT -------------------------------------------------------------------------------------------------------------------
!1005    Write(*     , Fmt_READ3 ) UnFile, IO_READ ; Write( UnFile, Fmt_READ3 ) UnFile, IO_READ ;
!        !#Call BEEP_FAIL ;
!        Write(*, Fmt_FL) ;  Write(UnInf, Fmt_FL) ; Write(*, Fmt_End) ; Read(*,*) ;  Return ;

End Subroutine Input_Analysis_HDF5 ;


! ***************************************************************************************************************************************************
!                                                                                                                                                  **
! Start:        16 April 2014                                                                                                                      **
! Last Update:  16 April 2014                                                                                                                      **
! Description: This Subroutine reads information for material visualization.                                                                       **
!                                                                                                                                                  **
! ***************************************************************************************************************************************************
! All truths are easy to understand once they are discovered; the point is to discover them.

Subroutine Input_Mat_Vis_HDF5 (                                                                                                &
!                                                                                                                              & ! Integer (1) Variables
!                                                                                                                              & ! Integer (2) Variables
!                                                                                                                              & ! Integer (4) Variables
NJ_Para,                                                                                                                       & ! Integer (8) Variables
!                                                                                                                              & ! Real Variables
idx_Mat_Vis_from,                                                                                                              & ! Integer Arrays
!                                                                                                                              & ! Real Arrays
ModelName, IndexRank, IndexSize, Model_InDir                                                                                   & ! Characters
!                                                                                                                              & ! Type
) ;

Use HDF5 ;

Implicit None ;

! =========================== Global Variables ======================================================================================================
! - Integer Variables -------------------------------------------------------------------------------------------------------------------------------
Integer , Intent(In)    :: NJ_Para ;

! - Real Variables ----------------------------------------------------------------------------------------------------------------------------------
! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
! - Integer Arrays ----------------------------------------------------------------------------------------------------------------------------------
Integer , Intent(Out), Dimension (:  )  :: idx_Mat_Vis_from ;

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
PetscErrorCode :: ErrPTC ;                       ! Error 

! - Integer Variables -------------------------------------------------------------------------------------------------------------------------------
!Integer  :: I, J, K ;               ! Loop Index

! HDF5 variables
Integer(HID_T)       :: id_Data_Mat_Vis_PV ;     ! File identifier for the material visualization in Paraview
Integer(HID_T)       :: dset_id_IS_Mat_Vis_PV ;  ! Dataset identifier for the index set, Paraview purposes only

Integer(HSIZE_T), DIMENSION(2) :: data_dims      ! Dataset dimensions

! - Unit NUMBERS OF EXTERNAL FILES ------------------------------------------------------------------------------------------------------------------
Integer, PARAMETER   :: Un_Mat_Vis     = 821  ;            ! the Unit number of Inversion Data Structure index sets
Integer   :: IO_File ;                ! Used for IOSTAT - Input Output Status - in the OPEN cammand.
Integer   :: IO_Read ;               ! Holds error of read statements
Integer   :: UnFile ;                ! Holds Unit of a file for error message


! - Real Variables ----------------------------------------------------------------------------------------------------------------------------------
! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
! - Integer Arrays ----------------------------------------------------------------------------------------------------------------------------------
! - Real Arrays -------------------------------------------------------------------------------------------------------------------------------------
! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------

! - Character Variables -----------------------------------------------------------------------------------------------------------------------------
Character (Kind = 1, Len = 30 ) :: helpme ;     

! - Logical Variables -------------------------------------------------------------------------------------------------------------------------------
! - Type DECLERATIONS -------------------------------------------------------------------------------------------------------------------------------

! =========================== Subroutine CODE =======================================================================================================

Write(*,*)
Write(*,*)"<<<<<<<< Input Analysis Mat-Vis HDF5 >>>>>>>>" ;

    ! Reading data for visualization (Paraview) in the HDF5 format

        Call h5open_f(ErrPTC) ;
        Call h5fopen_f( TRIM(Model_InDir)//'/'//TRIM(ModelName)//'_'//'Data_Mat_Vis_PV'//'_'//Trim(AdjustL(IndexSize))//'_'//Trim(AdjustL(IndexRank))//'.h5', H5F_ACC_RDONLY_F, id_Data_Mat_Vis_PV, ErrPTC) ;

        ! Equation number of nodes on this rank used to extract data of each node for Paraview
        Call h5dopen_f  ( id_Data_Mat_Vis_PV, "IndexSet_Mat_Vis",     dset_id_IS_Mat_Vis_PV,   ErrPTC ) ;

        data_dims(1) = NJ_Para ;
        data_dims(2) = 1 ;
        Call h5dread_f( dset_id_IS_Mat_Vis_PV, H5T_NATIVE_INTEGER,  idx_Mat_Vis_from, data_dims, ErrPTC) ;

        Call h5dclose_f( dset_id_IS_Mat_Vis_PV,   ErrPTC) ;

        ! Close the HDF5 file.
        CALL h5fclose_f( id_Data_Mat_Vis_PV, ErrPTC) ;

        ! Shut down HDF5.
        CALL h5close_f(ErrPTC) ;


Write(*    ,*) 'End Subroutine < Input_Mat_Vis_HDF5 >' ;
Write(UnInf,*) 'End Subroutine < Input_Mat_Vis_HDF5 >' ;

Return ;


End Subroutine Input_Mat_Vis_HDF5 ;


End Module Input_Subroutines_HDF5 ;
