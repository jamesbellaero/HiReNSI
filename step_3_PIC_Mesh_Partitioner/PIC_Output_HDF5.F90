
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                                                                                                                  ++
! Start:        18 June 2013                                                                                                                       ++
! Last Update:  05 August 2013                                                                                                                     ++
! Developed by: Babak Poursartip                                                                                                                   ++
!                                                                                                                                                  ++
! Description: THIS Module generates input files for the both the main code and Paraview based on HDF6 data format.                                ++
!                                                                                                                                                  ++
!              Subroutines                     Description                                                                                         ++
!              ==============                  ==============                                                                                      ++
!              Output                                                                                                                              ++
!                                                                                                                                                  ++
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Module Partition_Output_HDF5 ;

Use HDF5 ;
Use Parameters ;


Implicit None ;

  Interface
!    Module Procedure 
  End Interface    ;

Contains ;

! ***************************************************************************************************************************************************
!                                                                                                                                                  **
! Start:        18 June 2013                                                                                                                       **
! Last Update:  05 August 2013                                                                                                                     **
! Description: THIS Subroutine calculates application, PETSc and Local numbering of the elements and writes down all input data for each rank      **
! Called by:                                                                                                                                       **
!                                                                                                                                                  **
! ***************************************************************************************************************************************************

Subroutine OUTPUT_HDF5    (                                                                                                     &
NDim, MaxNNode, NDOF,                                                                                                           & ! Integer (1) Variables
NGroup, NMat, NPM,                                                                                                              & ! Integer (2) Variables
NParts,                                                                                                                         & ! Integer (4) Variables
NEL, NJ, NEQMTotal,                                                                                                             & ! Integer (8) Variables
!                                                                                                                               & ! Real Variables
LoadC,     EPart, MTEL, ELT, ELGR, JLoad, IDBC,                                                                                 &
D_NNZ_Stiff, O_NNZ_Stiff, D_NNZ_Damp, O_NNZ_Damp, D_NNZ_Mass, O_NNZ_Mass,     NDAN , NVAN, NAAN,  NEqRank, NNodeRank,           &
NoBndry_DRM, NoLayer_DRM,  NEL_Rank, NJ_Rank, Global_PETSc_Num, Local_PETSc_Num, INod, ID, ID_Application, ID_PETSc,            & ! Integer Arrays
PMat, PBLD, XYZ, UDis, PLoad, PML_DIM,                                                                                          & ! Real Arrays
ModelName, OutDir,                                                                                                              & ! Characters
Nodes, Param                                                                                                                    & ! Type
) ;


Implicit None ;

#include "finclude/petscsys.h"

! =========================== Global Variables ======================================================================================================

! - Integer Variables -------------------------------------------------------------------------------------------------------------------------------
Integer (Kind=Tiny), Intent(In)    :: NDim, MaxNNode, NDOF ;
Integer (Kind=Smll), Intent(In)    :: NGroup, NMat, NPM ;
Integer (Kind=Shrt), Intent(In)    :: NParts ;
Integer (Kind=Lng ), Intent(In)    :: NEL, NJ, NEQMTotal ;

! - Real Variables ----------------------------------------------------------------------------------------------------------------------------------
! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------

! - Integer Arrays ----------------------------------------------------------------------------------------------------------------------------------
Integer (Kind=Tiny), Intent(In),    Dimension (:  )  :: LoadC ;

Integer (Kind=Smll), Intent(In),    Dimension (:  )  :: MTEL, ELT ;  ! LTEL, 

Integer (Kind=Shrt), Intent(In),    Dimension (:  )  :: ELGR, EPart ;
Integer (Kind=Shrt), Intent(In),    Dimension (:  )  :: D_NNZ_Stiff, O_NNZ_Stiff, D_NNZ_Damp, O_NNZ_Damp, D_NNZ_Mass, O_NNZ_Mass ;   !?? see if required

Integer (Kind=Lng ), Intent(Inout), Dimension (:  )  :: NDAN , NVAN, NAAN ;
Integer (Kind=Lng ), Intent(In),    Dimension (:  )  :: NEqRank, NNodeRank, NoBndry_DRM, NoLayer_DRM, JLoad ;
Integer (Kind=Lng ), Intent(In),    Dimension (:  )  :: NEL_Rank, NJ_Rank, Global_PETSc_Num ;
Integer (Kind=Lng ), Intent(In),    Dimension (:,:)  :: Local_PETSc_Num, INod, ID, ID_Application, ID_PETSc ;
Integer (Kind=Lng ), Intent(In),    Dimension (:,:)  :: IDBC ;

! - Real Arrays -------------------------------------------------------------------------------------------------------------------------------------
Real (Kind=DBL),     Intent(In),    Dimension (:,:)  :: PMat, PBLD, XYZ, UDis, PLoad, PML_DIM ;

! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------

! - Character Variables -----------------------------------------------------------------------------------------------------------------------------
Character (Kind = 1, Len = 30 ) :: ModelName ;   ! Name of Input file
Character (Kind = 1, Len = 200) :: OutDir ;      ! Directory of output files (Results)

! - Logical Variables -------------------------------------------------------------------------------------------------------------------------------
! - Type DECLERATIONS -------------------------------------------------------------------------------------------------------------------------------
Type ( NodeID )      :: Nodes ;
Type ( BasicParam )  :: Param ;                  ! Holds basic parameters of each load case

! =========================== LOCAL Variables =======================================================================================================
PetscErrorCode :: ErrPTC ;                       ! Error 

! - Integer Variables -------------------------------------------------------------------------------------------------------------------------------
Integer (Kind=Tiny)  :: NNode ;                  ! Number of Nodes
Integer (Kind=Tiny)  :: X_CellType ;             ! XMDF cell type

Integer (Kind=Smll)  :: UnFile ;                 ! Unit Number of the Output file.
Integer (Kind=Smll)  :: IO_File ;                ! Used for IOSTAT - Input Output Status - in the Open cammand.
Integer (Kind=Smll)  :: IO_Write ;               ! Used for IOSTAT - Input Output Status - in the Write cammand.
Integer (Kind=Smll)  :: IDIM ;                   ! Loop index for NDim.
Integer (Kind=Smll)  :: IDOF ;                   ! Loop index for NDOF.
Integer (Kind=Smll)  :: INode ;                  ! Loop index for NNode.
Integer (Kind=Smll)  :: EType ;                  ! Element Type
Integer (Kind=Smll)  :: ERR_Alloc, ERR_DeAlloc ; ! Allocating and DeAllocating errors

Integer (Kind=Shrt)  :: IParts ;                 ! Loop index on the number of Partitions.
Integer (Kind=Shrt)  :: NLN_Rank ;               ! Number of Loaded Nodes of each Rank.
Integer (Kind=Shrt)  :: NSND_Rank ;              ! Number of Support Nodes with EqDisplacement in this rank
!Integer (Kind=Shrt)  :: LDSet ;                 ! Load Set.
Integer (Kind=Shrt)  :: ILType ;                 ! Load Type of Element.
!Integer (Kind=Shrt)  :: Inc ;                   ! Increment of elements for load sets.
Integer (Kind=Shrt)  :: NNDH_Rank ;              ! Holds Number of Nodes in which history of EqDisplacement is required in each partition
Integer (Kind=Shrt)  :: NNVH_Rank ;              ! Holds Number of Nodes in which history of Velocity is required in each partition
Integer (Kind=Shrt)  :: NNAH_Rank ;              ! Holds Number of Nodes in which history of Acceleration is required in each partition
Integer (Kind=Shrt)  :: NNBndry_DRM_Rank ;       ! Holds Number of nodes on the boudary of the internal and external domain on this rank for DRM analysis
Integer (Kind=Shrt)  :: NNLayer_DRM_Rank ;       ! Holds Number of nodes on the adjacent layer to the boudary of the internal and external domain on this rank for DRM analysis
Integer (Kind=Shrt)  :: NNDH ;                   ! Number of Nodes which History of Displacement is required in dynamic analysis.
Integer (Kind=Shrt)  :: NNVH ;                   ! Number of Nodes which History of Velocity is required in dynamic analysis.
Integer (Kind=Shrt)  :: NNAH ;                   ! Number of Nodes which History of Acceleration is required in dynamic analysis.
Integer (Kind=Shrt)  :: NNBndry_DRM ;            ! Number of nodes on the DRM boundary for the Domain Reduction Method
Integer (Kind=Shrt)  :: NNLayer_DRM ;            ! Number of nodes on the DRM layer for the Domain Reduction Method
Integer (Kind=Shrt)  :: NIDBC_Rank ;             ! Number of elements with pressure load on each rank 
Integer (Kind=Shrt)  :: LocalElN ;               ! counter on local element number

Integer (Kind=Lng )  :: NLN ;                    ! Number of Loaded Nodes in the original model
Integer (Kind=Lng )  :: NSND ;                   ! Number of Support Nodes with Displacements in the original model
Integer (Kind=Lng )  :: IEL ;                    ! Loop index on NEL.
Integer (Kind=Lng )  :: IJ ;                     ! Loop index on NJ.
Integer (Kind=Lng )  :: I, J, K, I1, I2 ;        ! Loop indeces.
Integer (Kind=Lng )  :: Lo_EL_Num ;              ! Load set of Element Number.
Integer (Kind=Lng )  :: NEqM_Mapping ;           ! Number of Modified EQuations for Mapping application numbering to PETSc numbering. !!?? at the end of the day see if we need this variable.
Integer (Kind=Lng )  :: LEqN ;                   ! Local Equation number, for ID_Local.
Integer (Kind=Lng )  :: NEQM ;                   ! Number of Modified Equations on this rank
Integer (Kind=Lng )  :: Counter ;                ! Counter
Integer (Kind=Lng )  :: LowerLimit ;             ! Lower Limit of node number based on the PETSc node numbering system, stored on a rank.
Integer (Kind=Lng )  :: UpperLimit ;             ! Upper Limit of node number based on the PETSc node numbering system, stored on a rank.
Integer (Kind=Lng )  :: PETScNodeNumber ;        ! A temporary varibale for saving node number based on the PETSc node numbering.
Integer (Kind=Lng )  :: Node ;                   ! Holds Node number
Integer (Kind=Lng )  :: NJ_Para ;                ! number joints in the model after downsampling, for Paraview
Integer (Kind=Lng )  :: NE_Para ;                ! number equations in the model after downsampling, for Paraview
Integer (Kind=Lng )  :: ConnSizePara ;           ! Size of connectivity vector for Paraview
Integer (Kind=Lng )  :: NEqM_Verifier ;          ! NEqM_Verifier

Integer(HID_T)       :: id_Geo ;                 ! File identifier for the Geometry in the main code
Integer(HID_T)       :: id_Geo_PV ;              ! File identifier for the Geometry in Paraview
Integer(HID_T)       :: id_Data_PV ;             ! File identifier for recording data for Paraview

Integer(HID_T)       :: dspace_id_XYZ ;          ! Dataspace identifier for coordinates
Integer(HID_T)       :: dspace_id_Cnn ;          ! Dataspace identifier for connectivity
Integer(HID_T)       :: dspace_id_Cnt ;          ! Dataspace identifier for constraints
Integer(HID_T)       :: dspace_id_App ;          ! Dataspace identifier for Application Numbering
Integer(HID_T)       :: dspace_id_PTC ;          ! Dataspace identifier for PETSc Numbering
Integer(HID_T)       :: dspace_id_Ind ;          ! Dataspace identifier for Indices (Index Set)

Integer(HID_T)       :: dspace_id_XYZ_PV ;       ! Dataspace identifier for coordinates for Paraview
Integer(HID_T)       :: dspace_id_Cnn_PV ;       ! Dataspace identifier for connectivity  for Paraview
Integer(HID_T)       :: dspace_id_Node_PV ;      ! Dataspace identifier for equation number of nodes on the rank for Paraview
Integer(HID_T)       :: dspace_id_IS_PV ;        ! Dataspace identifier for the index set, Paraview purpose only

Integer(HID_T)       :: dset_id_XYZ ;            ! Dataset identifier for coordinates
Integer(HID_T)       :: dset_id_Cnn ;            ! Dataset identifier for connectivity
Integer(HID_T)       :: dset_id_Cnt ;            ! Dataset identifier for constraints
Integer(HID_T)       :: dset_id_App ;            ! Dataset identifier for Application Numbering
Integer(HID_T)       :: dset_id_PTC ;            ! Dataset identifier for PETSc Numbering
Integer(HID_T)       :: dset_id_Ind ;            ! Dataset identifier for Indices (Index Set)

Integer(HID_T)       :: dset_id_XYZ_PV ;         ! Dataset identifier for coordinates for Paraview
Integer(HID_T)       :: dset_id_Cnn_PV ;         ! Dataset identifier for connectivity for Paraview
Integer(HID_T)       :: dset_id_Node_PV ;        ! Dataset identifier for global PETSc equation number of nodes on the rank for Paraview
Integer(HID_T)       :: dset_id_IS_PV ;          ! Dataset identifier for the index set, Paraview purposes only 

Integer(HSIZE_T), DIMENSION(2) :: dims           ! Dataset dimensions
Integer(HSIZE_T), DIMENSION(2) :: data_dims

! - Real Variables ----------------------------------------------------------------------------------------------------------------------------------
! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
! - Integer Arrays ----------------------------------------------------------------------------------------------------------------------------------
Integer (Kind=Shrt)  :: Local_ElN ( NEL ) ;                          ! Local Element Number
Integer (Kind=Shrt), Allocatable, Dimension(:  )  :: dset_data_ISet; ! Vector for Index Setting 
Integer (Kind=Shrt), Allocatable, Dimension(:,:)  :: dset_data_int ; ! Array to fill HDF5 file

! Watch out for large models, we may need to change the type to Lng. They are Shrt because HDF5 only accepts Shrt integer
Integer (Kind=Shrt), Allocatable, Dimension(:  )  :: App_Numbers ;   ! Holds Equation numbers based on the mesh generator node numbering 
Integer (Kind=Shrt), Allocatable, Dimension(:  )  :: PETSc_Numbers ; ! Holds Equation numbers based on PETSc node numbering 
Integer (Kind=Shrt), Allocatable, Dimension(:  )  :: Indices ;       ! 

Integer (Kind=Lng )  :: LTEL_Rank ( 4, Param%IntM( 1, 2) ) ;         ! Load Type Element - Param%IntM( 1, 2) = NPBL

Integer (Kind=Lng ), Allocatable, Dimension(:,:)  :: ID_Local ;      ! Holds Equation number of each rank
Integer (Kind=Lng ), Allocatable, Dimension(:  )  :: EqDis ;         ! Holds Equation numbers of all nodes' DOF in which the history of EqDisplacement is required based on PETSc node numbering 
Integer (Kind=Lng ), Allocatable, Dimension(:  )  :: EqVel ;         ! Holds Equation numbers of all nodes' DOF in which the history of velocity is required based on PETSc node numbering 
Integer (Kind=Lng ), Allocatable, Dimension(:  )  :: EqAcc ;         ! Holds Equation numbers of all nodes' DOF in which the history of acceleration is required based on PETSc node numbering 
Integer (Kind=Lng ), Allocatable, Dimension(:  )  :: NodeDis ;       ! Holds Node numbers in which the history of EqDisplacement is required based on Global node numbering 
Integer (Kind=Lng ), Allocatable, Dimension(:  )  :: NodeVel ;       ! Holds Node numbers in which the history of velocity is required based on Global node numbering 
Integer (Kind=Lng ), Allocatable, Dimension(:  )  :: NodeAcc ;       ! Holds Node numbers in which the history of acceleration is required based on Global node numbering 
Integer (Kind=Lng ), Allocatable, Dimension(:  )  :: ParaNodeNum ;   ! Temporary array used to determine the numnber of nodes in the model for Paraview, note that Paraview does not support spectral elements and we have to transform Spectral elements to Finite elements.
!Integer (Kind=Lng ), Allocatable, Dimension(:,:)  :: dset_data_int ; ! Array to fill HDF5 file

Integer (Kind=Lng ), Allocatable, Dimension(:,:)  :: List_NDAN_rank   ! af: this is for the prevention of "inverse crimes". The code map fine to coarse uses this information.


! - Real Arrays -------------------------------------------------------------------------------------------------------------------------------------
Real (Kind=DBL)   , Allocatable, Dimension(:,:)   :: dset_data_real ;! Array to fill HDF5 file  

! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
! - Character Variables -----------------------------------------------------------------------------------------------------------------------------
Character (Kind = 1, Len = 20)  :: IndexRank ;   ! A variable for storing the Rank number in Character format for adding at the end of the input file Name.
Character (Kind = 1, Len = 20)  :: IndexSize ;   ! A variable for storing the Size number in Character format for adding at the end of the input file Name.

! - Logical Variables -------------------------------------------------------------------------------------------------------------------------------
! - Type DECLERATIONS -------------------------------------------------------------------------------------------------------------------------------
Type ( BasicParam )    :: ParamPartition ;       ! Holds basic parameters of each load case in each partition

! =========================== Subroutine CODE =======================================================================================================

! Important notice: Due to the restriction of HDF5 format, we use real(4) and integer(4) in the output.

Write (*,*)"Output HDF5"

Allocate ( ParamPartition%IntM ( 7, 6), ParamPartition%RealM ( 7, 6) ) ;

  If ( LoadC (5) /= 0_Tiny ) Then ;  ! Dynamic Analysis

    NNDH = Param%IntM( 6, 1) ;
    NNVH = Param%IntM( 6, 2) ;
    NNAH = Param%IntM( 6, 3) ;

!/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_
    Allocate ( List_NDAN_rank ( NParts , NNDH ) )
    List_NDAN_rank = 0
!/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_


    ! Allocating Required Arrays Dynamic Analysis
    Allocate ( EqDis ( NNDH * NDim ), EqVel ( NNVH * NDim ), EqAcc ( NNAH * NDim ), NodeDis ( NNDH ), NodeVel ( NNVH ), NodeAcc ( NNAH ),     STAT = ERR_Alloc) ; 
      If ( ERR_Alloc /= 0 ) Then ;
        Write (*, Fmt_ALLCT) ERR_Alloc ;  Write (UnInf, Fmt_ALLCT) ERR_Alloc ;
        Write(*, Fmt_FL) ;  Write(UnInf, Fmt_FL) ; Read(*, Fmt_End) ;  Stop ;
      End If ;

  End If ;

NEqM_Verifier = 0_Lng ;

! - Writing down the input data for each process ----------------------------------------------------------------------------------------------------
  Do IParts = 1, NParts ;

    Write (IndexRank, *) IParts - 1_Shrt ; ! Converts Rank number to Character foramt for the file Name
    Write (IndexSize, *) NParts ;          ! Converts Size number to Character foramt for the file Name

    Write (*,*)"Writing files for partition: ", IParts;

    ! - Computations --------------------------------------------------------------------------------------------------------------------------------
    Write (*,*)"PreCalculations for Partition number:";

    ! Calculating Local Element Number for each rank - Used in IDBC
    Local_ElN = 0 ;
    LocalElN  = 0 ;
      Do IEl = 1, NEl ;
        If ( EPart ( IEl ) == IParts ) Then ;
          LocalElN = LocalElN + 1 ;
          Local_ElN ( IEl ) = LocalElN ;
        End If ;
      End Do ;

    ! Joint load of every case.
    NLN_Rank = 0_Shrt ;
      IF ( LoadC ( 3 ) /= 0 ) Then ;
        NLN = Param%IntM( 3, 1) ;
          DO I = 1, NLN ;
            If ( Local_PETSc_Num ( JLoad ( I ) , IParts ) /= 0_Lng ) Then ;
              If ( Nodes%Rep ( JLoad ( I ) ) == 1_Shrt ) Then ;
                NLN_Rank = NLN_Rank + 1_Shrt ;
              Else If ( Nodes%Rep ( JLoad ( I ) ) /= 1_Shrt .AND. IParts == MinVal ( Nodes%Locs( JLoad ( I ),1:Nodes%Rep ( JLoad ( I ) ) ) ) ) Then ;
                NLN_Rank = NLN_Rank + 1_Shrt ;
              End If ;
            End If ; 
          End Do ;
      EndIF ;

    ! Supports Displacements
    NSND_Rank = 0_Shrt ;
      IF ( LoadC ( 4 ) /= 0 ) Then ;
        NSND = Param%IntM( 4, 1) ;
          DO I = 1, NSND ;
            If ( Local_PETSc_Num ( UDis( NDOF + 1, I ), IParts ) /= 0_Lng ) Then ; 
              If ( Nodes%Rep ( UDis( NDOF + 1, I ) ) == 1_Shrt ) Then ;
                NSND_Rank = NSND_Rank + 1_Shrt ;
              Else If ( Nodes%Rep ( UDis( NDOF + 1, I ) ) /= 1_Shrt .AND. IParts == MinVal ( Nodes%Locs( UDis( NDOF + 1, I ),1:Nodes%Rep ( UDis( NDOF + 1, I ) ) ) ) ) Then ;
                NSND_Rank = NSND_Rank + 1_Shrt ;
              End If ;
            End If ; 
          End Do ;
      EndIF ;

    Allocate ( ID_Local ( NJ_Rank ( IParts ), NDOF ) ) ;

      ! Forming the Local ID for each Process
      DO J = 1, NJ ;
        If ( Local_PETSc_Num ( J, IParts ) /= 0_Lng ) Then ; !?? check the oreintation 
          ID_Local ( Local_PETSc_Num ( J, IParts ), : ) = ID ( J, : ) ;
        End If ;
      End Do ;

    ! Modified Number of Equations for each rank
    ! Obtaining the Local Equation Number
    LEqN = 0 ;
      Do IJ = 1, NJ_Rank ( IParts ) ;
        DO J = 1, NDOF ;
          IF ( ID_Local ( IJ, J ) == 1 ) Then ;
            ID_Local ( IJ, J ) = 0 ; 
          Else If ( ID_Local ( IJ, J ) == 0 ) Then ;
            LEqN = LEqN + 1 ;
            ID_Local ( IJ, J ) = LEqN ;
          End If ;
        End Do ;
      End Do ;
    NEQM = LEqN ; ! napp

    ! Obtanining the size of app_numbers and PETSc_numbers
      If ( IParts == 1_Shrt ) then ;
        NEqM_Mapping = NEqRank ( IParts ) ;
      Else ;
        NEqM_Mapping = NEqRank ( IParts ) - NEqRank ( IParts - 1_Shrt ) ;
      End If ;

NEqM_Verifier = NEqM_Verifier + NEqM_Mapping ;

    ! Pressure on elements
    NIDBC_Rank = 0_Shrt ;
      If ( LoadC (2) /= 0_Tiny .or. LoadC (5) /= 0_Tiny ) Then ;
        Do I = 1, Param%IntM( 2, 4) ;  ! NIDBC
          If ( EPart ( IDBC ( I, 1 ) ) == IParts ) Then ;
            NIDBC_Rank = NIDBC_Rank + 1 ;
          End If ;
        End Do ;
      End If ;

    ! Equation numbers of nodes in which history is required for dynamic analysis
    NNDH_Rank = 0_Lng ;
    NNVH_Rank = 0_Lng ;
    NNAH_Rank = 0_Lng ;

      ! Obtaining Equation number of nodes in which history of Displacement is required based on PETSc node numbering
      If ( LoadC (5) /= 0_Tiny ) Then ;
        NNDH = Param%IntM( 6, 1) ;
          Do I = 1, NNDH ;
            If ( NDAN ( I ) == 0_Shrt  ) Cycle ; ! why? because we do not want to apply the loads on the boundary of the partitions twice. see inside the loop
              If ( Local_PETSc_Num ( NDAN ( I ), IParts ) /= 0_Lng ) Then ;
                NNDH_Rank = NNDH_Rank + 1_Lng ;
                EqDis ( ( NNDH_Rank - 1_Lng )* NDim+1:( NNDH_Rank - 0_Lng )* NDim ) = ID_PETSc ( Global_Petsc_Num ( NDAN ( I ) ), 1:NDim )  ;
                NodeDis ( NNDH_Rank ) = NDAN ( I ) ;
                NDAN ( I ) = 0_Shrt ; ! In this way, the load on this node will be applied once in one rank, if the node is located at the border of two or more ranks 

!/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_
    
                List_NDAN_rank ( IParts , NNDH_Rank ) = NodeDis ( NNDH_Rank )

!/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_

              End If ;
          End Do ;
        EqDis (:) = EqDis (:)- 1_Lng ; ! PETSc numbering starts from 0

        ! Obtaining Equation number of nodes in which history of velocity is required based on PETSc node numbering
        NNVH = Param%IntM( 6, 2) ;
          Do I = 1, NNVH ;
            If ( NVAN ( I ) == 0_Shrt  ) Cycle ;
              If ( Local_PETSc_Num ( NVAN ( I ), IParts ) /= 0_Lng ) Then ;
                NNVH_Rank = NNVH_Rank + 1_Lng ;
                EqVel ( ( NNVH_Rank - 1_Lng )* NDim+1:( NNVH_Rank - 0_Lng )* NDim ) = ID_PETSc ( Global_Petsc_Num ( NVAN ( I ) ), 1:NDim )  ;
                NodeVel ( NNVH_Rank ) = NVAN ( I ) ;
                NVAN ( I ) = 0_Shrt ; ! In this way, the load on this node will be applied once in one rank, if the node is located at the border of two or more ranks 
              End If ;
          End Do ;
        EqVel (:) = EqVel (:)- 1_Lng ; ! PETSc numbering starts from 0

        ! Obtaining Equation number of nodes in which history of acceleration is required based on PETSc node numbering 
        NNAH = Param%IntM( 6, 3) ;
          Do I = 1, NNAH ;
            If ( NAAN ( I ) == 0_Shrt  ) Cycle ;
              If ( Local_PETSc_Num ( NAAN ( I ), IParts ) /= 0_Lng ) Then ;
                NNAH_Rank = NNAH_Rank + 1_Lng ;
                EqAcc ( ( NNAH_Rank - 1_Lng )* NDim+1:( NNAH_Rank - 0_Lng )* NDim ) = ID_PETSc ( Global_Petsc_Num ( NAAN ( I ) ), 1:NDim )  ;
                NodeAcc ( NNAH_Rank ) = NAAN ( I ) ;
                NAAN ( I ) = 0_Shrt ; ! In this way, the load on this node will be applied once in one rank, if the node is located at the border of two or more ranks 
              End If ;
          End Do ;
        EqAcc (:) = EqAcc (:)- 1_Lng ; ! PETSc numbering starts from 0

          ! required nodes on the boundary and adjacent layer for DRM analysis
          If ( Param%IntM( 5, 3) == 3_Tiny ) Then ;  ! LoadType

            ! Obtaining number of nodes on the boundary of the internal and external domain on this rank for DRM analysis
            NNBndry_DRM_Rank = 0_Shrt ;
            NNBndry_DRM = Param%IntM( 7, 1) ;
              Do I = 1, NNBndry_DRM ;
                Node = NoBndry_DRM ( I ) ;
                If ( ( Nodes%Rep (Node) == 1_Shrt .AND. Any ( Nodes%Locs ( Node, : ) == IParts ) ) .OR. ( Nodes%Rep (Node) /= 1_Shrt .AND. IParts == MinVal ( Nodes%Locs( Node, :), MASK = Nodes%Locs( Node, :) > 0_Shrt ) ) ) NNBndry_DRM_Rank = NNBndry_DRM_Rank + 1_Shrt ;
              End Do ;

            ! Obtaining Number of nodes on the adjacent layer to the boudary of the internal and external domain on this rank for DRM analysis
            NNLayer_DRM_Rank = 0_Shrt ;
            NNLayer_DRM = Param%IntM( 7, 2) ;
              Do I = 1, NNLayer_DRM ;
                Node = NoLayer_DRM ( I ) ;
                If ( ( Nodes%Rep (Node) == 1_Shrt .AND. Any ( Nodes%Locs ( Node, : ) == IParts ) ) .OR. ( Nodes%Rep (Node) /= 1_Shrt .AND. IParts == MinVal ( Nodes%Locs( Node, :), MASK = Nodes%Locs( Node, :) > 0_Shrt ) ) ) NNLayer_DRM_Rank = NNLayer_DRM_Rank + 1_Shrt ;
              End Do ;

          End If ;
      End If ;

    ParamPartition%IntM  (:,:) = 0_Lng ;
    ParamPartition%RealM (:,:) = 0._Dbl ;

      If (LoadC (1) /= 0_Tiny ) Then ;  ! Body Force
        ParamPartition%IntM( 1, 1) = Param%IntM( 1, 1) ;  ParamPartition%IntM( 1, 2) = Param%IntM( 1, 2) ; ! NBLD   = NGroup, NPBL   = NDim ;
      End If ;

      If (LoadC (2) /= 0_Tiny ) Then ;  ! Pressure on elements
        ParamPartition%IntM( 2, 4) = NIDBC_Rank ;       ! PR_Type, LDIR, PRZERO, HREF, PR
      End If ;

      If (LoadC (3) /= 0_Tiny ) Then ;  ! Concentrated force
        ParamPartition%IntM( 3, 1) = NLN_Rank ; ! NLN
      End If ;

      If (LoadC (4) /= 0_Tiny ) Then ;  ! Number of nodes (Supports) with predefined Displacements
        ParamPartition%IntM( 4, 1) = NSND_Rank; ! NSDN
      End If ;

!      If ( LoadC (5) == 0_Tiny ) Then ;  ! Static Analysis
!        ParamPartition%IntL( 5, 1 ) = Param%IntL( 5, 1 ); ! STStep
!      End If ;

      If ( LoadC (5) /= 0_Tiny ) Then ;  ! Dynamic Analysis
        ParamPartition%IntM( 5, 3) = Param%IntM( 5, 3) ; ! LoadType   

        ParamPartition%IntM( 6, 1) = NNDH_Rank ; ! NNDH
        ParamPartition%IntM( 6, 2) = NNVH_Rank ; ! NNVH
        ParamPartition%IntM( 6, 3) = NNAH_Rank ; ! NNAH
        ParamPartition%IntM( 6, 4) = Param%IntM( 6, 4) ; ! NDamp
        ParamPartition%IntM( 6, 5) = Param%IntM( 6, 5) ; ! NEnergy
        ParamPartition%IntM( 6, 6) = Param%IntM( 6, 6) ; ! PARAM_Type ! pml analysis

          If      ( Param%IntM( 5, 3) == 1_Tiny ) Then ;  ! Dynamic Pressure
            ParamPartition%IntM( 2, 4) = NIDBC_Rank ;  
          !Else If ( Param%IntM( 5, 3) == 2_Tiny ) Then ; ! Base Acceleration
          !  ParamPartition%IntL( 7, 1) = Param%IntL( 7, 1) ;  ParamPartition%RealL( 7, 1) = Param%RealL( 7, 1) ; ! NASTEP, G   Number of base Accelerations to be read in the input file
          Else If ( Param%IntM( 5, 3) == 3_Tiny ) Then ; ! Domain Reduction Method
            ParamPartition%IntM  ( 7, 1 ) = NNBndry_DRM_Rank ; ! NNBndry_DRM
            ParamPartition%IntM  ( 7, 2 ) = NNLayer_DRM_Rank ; ! NNLayer_DRM
          End If ;

      End If ;

    ! - Output FILEs --------------------------------------------------------------------------------------------------------------------------------
    Write (*,*)"Opening files ...";
    ! Data file
    UnFile = UN_Out ;
    Open ( Unit = UnFile, FILE = TRIM(ModelName)//'_'//Trim(AdjustL(IndexSize))//'_'//Trim(AdjustL(IndexRank))//'.dataModel', ERR =  1001, IOSTAT = IO_File, ACCESS = 'SEQUENTIAL', ACTION ='Write', ASYNCHRONOUS = 'NO', BLANK = 'NULL', BLOCKSIZE = 0, DEFAULTFILE = TRIM(OutDir), DisPOSE = 'KEEP', FORM = 'FORMATTED', POSITION = 'ASIS', STATUS ='REPLACE') ;

!!    ! Non-Zeros of Stiffness file (.NNZStiff)
!!    UnFile = Un_OutStiff ;
!!    Open ( Unit = UnFile, FILE = TRIM(ModelName)//'_'//Trim(AdjustL(IndexSize))//'_'//Trim(AdjustL(IndexRank))//'.NNZStiff', ERR =  1001, IOSTAT = IO_File, ACCESS = 'SEQUENTIAL', ACTION ='Write', ASYNCHRONOUS = 'NO', BLANK = 'NULL', BLOCKSIZE = 0, DEFAULTFILE = TRIM(OutDir), DisPOSE = 'KEEP', FORM = 'FORMATTED', POSITION = 'ASIS', STATUS ='REPLACE') ;

!!    ! Non-Zeros of Damping file (.NNZDamp)
!!    UnFile = Un_OutDamp ;
!!    Open ( Unit = UnFile, FILE = TRIM(ModelName)//'_'//Trim(AdjustL(IndexSize))//'_'//Trim(AdjustL(IndexRank))//'.NNZDamp', ERR =  1001, IOSTAT = IO_File, ACCESS = 'SEQUENTIAL', ACTION ='Write', ASYNCHRONOUS = 'NO', BLANK = 'NULL', BLOCKSIZE = 0, DEFAULTFILE = TRIM(OutDir), DisPOSE = 'KEEP', FORM = 'FORMATTED', POSITION = 'ASIS', STATUS ='REPLACE') ;

!!    ! Non-Zeros of Damping file (.NNZMass)
!!    UnFile = Un_OutMass ;
!!    Open ( Unit = UnFile, FILE = TRIM(ModelName)//'_'//Trim(AdjustL(IndexSize))//'_'//Trim(AdjustL(IndexRank))//'.NNZMass', ERR =  1001, IOSTAT = IO_File, ACCESS = 'SEQUENTIAL', ACTION ='Write', ASYNCHRONOUS = 'NO', BLANK = 'NULL', BLOCKSIZE = 0, DEFAULTFILE = TRIM(OutDir), DisPOSE = 'KEEP', FORM = 'FORMATTED', POSITION = 'ASIS', STATUS ='REPLACE') ;

    ! Material Properties (.Mat)
    UnFile = Un_OutMat ;
    Open ( Unit = UnFile, FILE = TRIM(ModelName)//'_'//Trim(AdjustL(IndexSize))//'_'//Trim(AdjustL(IndexRank))//'.Mat', ERR =  1001, IOSTAT = IO_File, ACCESS = 'SEQUENTIAL', ACTION ='Write', ASYNCHRONOUS = 'NO', BLANK = 'NULL', BLOCKSIZE = 0, DEFAULTFILE = TRIM(OutDir), DisPOSE = 'KEEP', FORM = 'FORMATTED', POSITION = 'ASIS', STATUS ='REPLACE') ;

    ! Create files based on HDF5 format for large files (connectivity, coordinates, constraints, ...   )
    Call h5open_f(ErrPTC)

    Call h5fcreate_f( TRIM(OutDir)//'/'//TRIM(ModelName)//'_'//'Geometry'//'_'//Trim(AdjustL(IndexSize))//'_'//Trim(AdjustL(IndexRank))//'.h5', H5F_ACC_TRUNC_F, id_Geo, ErrPTC) ;      ! Geometry file for input files

    ! - Create the dataspaces
    ! Coordinate file for the main code(.XYZ)
    dims(1) = NJ_Rank ( IParts ) ;
    dims(2) = NDim ;
    Call h5screate_simple_f(2, dims, dspace_id_XYZ, ErrPTC) ;

    ! Connectivity file for the main code (.Cnn)
    dims(1) = NEL_Rank ( IParts ) ;
    dims(2) = MaxNNode + 3_Tiny ;
    Call h5screate_simple_f(2, dims, dspace_id_Cnn, ErrPTC) ;

    ! Coordinate file (.Cnt)
    dims(1) = NJ_Rank ( IParts ) ;
    dims(2) = NDOF ;
    Call h5screate_simple_f(2, dims, dspace_id_Cnt, ErrPTC) ;

    ! Application Numbering file (.App)
    dims(1) = 1_Shrt ;
    dims(2) = NEqM_Mapping ;
    !dims(2) = 3000000 ;
    Call h5screate_simple_f(2, dims, dspace_id_App, ErrPTC) ;

    ! PETSc Numbering file (.PTC)
    dims(1) = 1_Shrt ;
    dims(2) = NEqM_Mapping ;
    Call h5screate_simple_f(2, dims, dspace_id_PTC, ErrPTC) ;

    ! Indices file (.Ind)
    dims(1) = 1_Shrt ;
    dims(2) = NEqM ;
    Call h5screate_simple_f(2, dims, dspace_id_Ind, ErrPTC) ;

    ! Create the dataset with default properties.
    Call h5dcreate_f(id_Geo, "XYZ",             H5T_NATIVE_DOUBLE,  dspace_id_XYZ, dset_id_XYZ, ErrPTC) ;
    Call h5dcreate_f(id_Geo, "Connectivity",    H5T_NATIVE_INTEGER, dspace_id_Cnn, dset_id_Cnn, ErrPTC) ;
    Call h5dcreate_f(id_Geo, "Constraints",     H5T_NATIVE_INTEGER, dspace_id_Cnt, dset_id_Cnt, ErrPTC) ;   ! see type of integers for large numbers
    Call h5dcreate_f(id_Geo, "ApplicationNum",  H5T_NATIVE_INTEGER, dspace_id_App, dset_id_App, ErrPTC) ;   ! see type of integers for large numbers
    Call h5dcreate_f(id_Geo, "PETScNum",        H5T_NATIVE_INTEGER, dspace_id_App, dset_id_PTC, ErrPTC) ;   ! see type of integers for large numbers
    Call h5dcreate_f(id_Geo, "Indices",         H5T_NATIVE_INTEGER, dspace_id_Ind, dset_id_Ind, ErrPTC) ;   ! see type of integers for large numbers

    ! - Main code geometry file  --------------------------------------------------------------------------------------------------------------------
    ! Coordinates of nodes for the main code - HDF5
    Write (*,*)"Write Coordinates ..." ;
    Allocate ( dset_data_real( NJ_Rank ( IParts ), NDim ) ) ;

    Counter = 0_Lng ;
      DO K = 1, NJ ;
        If ( Local_PETSc_Num ( K, IParts ) /= 0_Lng ) Then ; !?? check the oreintation 
          Counter = Counter + 1_Lng ;
          dset_data_real ( Counter, 1:NDim ) = XYZ ( K, 1:NDim ) ; 
        End If ;
      End Do ;

    ! Write the dataset.
    data_dims(1) = NJ_Rank ( IParts ) ;
    data_dims(2) = NDim ;
    Call h5dwrite_f(dset_id_XYZ, H5T_NATIVE_DOUBLE, dset_data_real, data_dims, ErrPTC)

    DeAllocate ( dset_data_real ) ;

    ! Element Connectivity for the main code
    Write (*,*)"Write Element Connectivity ..." ;
    Allocate ( dset_data_int ( NEL_Rank ( IParts ), MaxNNode + 3_Tiny) ) ;

    Counter = 0_Lng ;
    dset_data_int (:,:) = 0_Shrt ;

      DO IEL = 1, NEL ;

        If ( EPart ( IEL ) == IParts ) Then ;
          EType = ElT ( IEl ) ; ! Element Type
          Counter = Counter + 1_Lng ;

            ! determinig NNode and cell types based on http://www.xdmf.org/index.php/XDMF_Model_and_Format                ! determinig NNode for this element     !?? see if all elements included here
            IF      ( EType == El2d4NSldPN  ) Then ;  NNode = 4_Tiny ; ! 4 node  -2D - Solid - PLANE STRESS
            Else If ( EType == El2d4NPMLPN  ) Then ;  NNode = 4_Tiny ; ! 4 node  -2D - PML - PLANE STRESS
            Else If ( EType == El2d8NSldPN  ) Then ;  NNode = 8_Tiny ; ! 8 node  -2D - Solid - PLANE STRAIN
            Else If ( EType == El2d8NPMLPN  ) Then ;  NNode = 8_Tiny ; ! 8 node  -2D - PML - PLANE STRAIN
            Else If ( EType == El2d6NSldPN  ) Then ;  NNode = 6_Tiny ; ! 6 node  -2D - Solid - PLANE STRAIN
            Else If ( EType == El2d6NPMLPN  ) Then ;  NNode = 6_Tiny ; ! 6 node  -2D - PML - PLANE STRAIN
            Else If ( EType == El2d3NSldPN  ) Then ;  NNode = 3_Tiny ; ! 3 node  -2D - Solid - PLANE STRAIN
            Else If ( EType == El2d3NPMLPN  ) Then ;  NNode = 3_Tiny ; ! 3 node  -2D - PML - PLANE STRAIN
            Else If ( EType == SEl2d9NSldPN ) Then ;  NNode = 9_Tiny ; ! 9 node  -2D - Solid - PLANE STRAIN
            Else If ( EType == SEl2d9NPMLPN ) Then ;  NNode = 9_Tiny ; ! 9 node  -2D - NPML - PLANE STRAIN    !?? using 8 noded element for spectral 9 noded elements
            Else If ( EType == SEl2d9NMPMLPN) Then ;  NNode = 9_Tiny ; ! 9 node  -2D - MPML - PLANE STRAIN
            Else If ( EType == SEl3d27NSld  ) Then ;  NNode = 27_Tiny; ! 27 node -3D - Solid
            Else If ( EType == SEl3d27NMPML ) Then ;  NNode = 27_Tiny; ! 27 node -3D - MPML
            Else If ( EType == SEl2d7NSldPN ) Then ;  NNode = 7_Tiny ; ! 7 node -2D - solid - plane strain
            Else If ( EType == SEl2d7NMPMLPN) Then ;  NNode = 7_Tiny ; ! 7 node -2D - MPML - plane strain
            Else If ( EType == SEl2d9NSldSH ) Then ;  NNode = 9_Tiny ; ! 9 node - 2D - solid ( quadrilateral ) - Spectral Element - SH waves
            Else If ( EType == SEl2d9NPMLSH ) Then ;  NNode = 9_Tiny ; ! 9 node - 2D - PML   ( quadrilateral ) - Spectral Element - SH waves
            Else If ( EType == SEl3d8NSld   ) Then ;  NNode = 8_Tiny ; ! 8 node -3D - Solid
            Else If ( EType == SEl3d8NMPML  ) Then ;  NNode = 8_Tiny ; ! 8 node -3D - MPML
            End If ;

          dset_data_int ( Counter, 1:NNode)            = Local_PETSc_Num ( INOD ( 1:NNode, IEL ), IParts ) ;
          dset_data_int ( Counter, MaxNNode + 1_Tiny ) = MTEL ( IEL ) ;
          dset_data_int ( Counter, MaxNNode + 2_Tiny ) = ELT ( IEL ) ;
          dset_data_int ( Counter, MaxNNode + 3_Tiny ) = ELGR ( IEL ) ;
        End If ;

      End Do ;

      If ( Counter /= NEL_Rank ( IParts ) ) Then ; 
        Write (*,*)"There is a mistake in the connectivity portion of the output code."
        Stop ;
      End If ;

    data_dims(1) = NEL_Rank ( IParts ) ;
    data_dims(2) = MaxNNode + 3_Tiny ;
    Call h5dwrite_f(dset_id_Cnn, H5T_NATIVE_INTEGER, dset_data_int, data_dims, ErrPTC ) ;

    DeAllocate ( dset_data_int ) ;

    ! - Paraview geometry file  ---------------------------------------------------------------------------------------------------------------------
    ! determine the number of nodes in the model for Paraview model
    Write (*,*)"Write Data for Paraview ..." ;

    Allocate ( ParaNodeNum ( NJ ) ) ; ! 

    ParaNodeNum (:) = 0_Lng ;
    NJ_Para         = 0_Lng ;
    ConnSizePara    = 0_Lng ;

      ! This loop finds the nodes of Finite Elements, In other words, it eliminates spectral nodes.
      Do IEl = 1, NEl ;  !
        If ( EPart ( IEL ) == IParts ) Then ;

          EType = ElT ( IEl ) ; ! Element Type

          ! Determine the number of nodes corresponding to a FE.
            IF      ( EType == El2d4NSldPN  ) Then ;  NNode = 4_Tiny ; ! 4 node  -2D - Solid - PLANE STRESS
            Else If ( EType == El2d4NPMLPN  ) Then ;  NNode = 4_Tiny ; ! 4 node  -2D - PML - PLANE STRESS
            Else If ( EType == El2d8NSldPN  ) Then ;  NNode = 8_Tiny ; ! 8 node  -2D - Solid - PLANE STRAIN
            Else If ( EType == El2d8NPMLPN  ) Then ;  NNode = 8_Tiny ; ! 8 node  -2D - PML - PLANE STRAIN
            Else If ( EType == El2d6NSldPN  ) Then ;  NNode = 6_Tiny ; ! 6 node  -2D - Solid - PLANE STRAIN
            Else If ( EType == El2d6NPMLPN  ) Then ;  NNode = 6_Tiny ; ! 6 node  -2D - PML - PLANE STRAIN
            Else If ( EType == El2d3NSldPN  ) Then ;  NNode = 3_Tiny ; ! 3 node  -2D - Solid - PLANE STRAIN
            Else If ( EType == El2d3NPMLPN  ) Then ;  NNode = 3_Tiny ; ! 3 node  -2D - PML - PLANE STRAIN
            Else If ( EType == SEl2d9NSldPN ) Then ;  NNode = 8_Tiny ; ! 9 node  -2D - Solid - PLANE STRAIN
            Else If ( EType == SEl2d9NPMLPN ) Then ;  NNode = 8_Tiny ; ! 9 node  -2D - NPML - PLANE STRAIN    !?? using 8 noded element for spectral 9 noded elements
            Else If ( EType == SEl2d9NMPMLPN) Then ;  NNode = 8_Tiny ; ! 9 node  -2D - MPML - PLANE STRAIN
            Else If ( EType == SEl3d27NSld  ) Then ;  NNode = 20_Tiny; ! 27 node -3D - Solid
            Else If ( EType == SEl3d27NMPML ) Then ;  NNode = 20_Tiny; ! 27 node -3D - MPML
            Else If ( EType == SEl2d7NSldPN ) Then ;  NNode = 6_Tiny ; ! 7 node -2D - solid - plane strain
            Else If ( EType == SEl2d7NMPMLPN) Then ;  NNode = 6_Tiny ; ! 7 node -2D - MPML - plane strain
            Else If ( EType == SEl2d9NSldSH ) Then ;  NNode = 8_Tiny ; ! 9 node - 2D - solid ( quadrilateral ) - Spectral Element - SH waves
            Else If ( EType == SEl2d9NPMLSH ) Then ;  NNode = 8_Tiny ; ! 9 node - 2D - PML   ( quadrilateral ) - Spectral Element - SH waves
            Else If ( EType == SEl3d8NSld   ) Then ;  NNode = 8_Tiny ; ! 8 node -3D - Solid
            Else If ( EType == SEl3d8NMPML  ) Then ;  NNode = 8_Tiny ; ! 8 node -3D - MPML
            Else ; Write(*,*)"Element number is not in the list"; Stop ;
            End If ;

          ConnSizePara = ConnSizePara + Int4(NNode + 1_Tiny) ;

            Do INode = 1, NNode;

              Node = INod (INode, IEl) ;

                If ( ParaNodeNum ( Node ) == 0_Lng ) Then ;  ! This if condition is for determining the number of FEM nodes
                  NJ_Para = NJ_Para + 1_Lng ;
                  ParaNodeNum ( Node ) = NJ_Para ;           ! This array renumbers nodes for Paraview. THE ORDER OF NODE NUMBERING IN CONSISTENT EVERYWHERE.
                End If ;

            End Do ;

        End If ;

      End Do ;

    ParamPartition%IntM( 5, 4) = ConnSizePara ; ! Save this number for the main code 

    DeAllocate ( ParaNodeNum ) ;
    Print *, 'OutDir - Parview', OutDir ;
    Call h5fcreate_f( TRIM(OutDir)//'/'//TRIM(ModelName)//'_'//'Geometry_PV'//'_'//Trim(AdjustL(IndexSize))//'_'//Trim(AdjustL(IndexRank))//'.h5', H5F_ACC_TRUNC_F, id_Geo_PV, ErrPTC) ;      ! Geometry file for Paraview
    Call h5fcreate_f( TRIM(OutDir)//'/'//TRIM(ModelName)//'_'//'Data_PV'//'_'//Trim(AdjustL(IndexSize))//'_'//Trim(AdjustL(IndexRank))//'.h5', H5F_ACC_TRUNC_F, id_Data_PV, ErrPTC) ;         ! Index number of DOFs in Paraview to retrive data from the solution vector

    ! Coordinate file for Paraview (.XYZ)
    dims(2) = NJ_Para ;
    dims(1) = 3 ; !NDim ; ! This always should be 3. For 3D models set the third dimension to zero.

    Call h5screate_simple_f(2, dims, dspace_id_XYZ_PV, ErrPTC) ;

    ! Equation number of nodes (DOFs) on this rank used to retrieve solution for Paraview
    dims(1) = NJ_Para ;
    dims(2) = NDim ;
    Call h5screate_simple_f(2, dims, dspace_id_Node_PV, ErrPTC) ;

    ! Connectivity file for Paraview (.Cnn)
    dims(2) = ConnSizePara;
    dims(1) = 1_Smll ;
    Call h5screate_simple_f(2, dims, dspace_id_Cnn_PV, ErrPTC) ;

    !Call h5dcreate_f(id_Geo_PV,  "XYZ",          H5T_NATIVE_REAL,     dspace_id_XYZ_PV,  dset_id_XYZ_PV,  ErrPTC) ;
    Call h5dcreate_f(id_Geo_PV,  "XYZ",          H5T_NATIVE_DOUBLE,  dspace_id_XYZ_PV,  dset_id_XYZ_PV,  ErrPTC) ;
    Call h5dcreate_f(id_Geo_PV,  "Connectivity", H5T_NATIVE_INTEGER, dspace_id_Cnn_PV,  dset_id_Cnn_PV,  ErrPTC) ;
    Call h5dcreate_f(id_Data_PV, "EquationNum",  H5T_NATIVE_INTEGER, dspace_id_Node_PV, dset_id_Node_PV, ErrPTC) ;

    Write (*,*)"Write Data for Paraview - Coordinates..." ;

    ! Coordinates of nodes for Paraview
    Allocate ( dset_data_real( 3, NJ_Para  ) );  ! Holds coordinates     write xyz always
    Allocate ( dset_data_int ( NJ_Para, NDim ) );  ! Holds global PETSc equation numbers of nodes on the rank for Paraview

    ! Coordinates
    dset_data_real (3, :) = 0.0_Dbl ;

    Counter = 0_Lng ; 
      DO K = 1, NJ ;
        If ( Local_PETSc_Num ( K, IParts ) /= 0_Lng ) Then ; !?? check the oreintation 

          Counter = Counter + 1_Lng ;
          dset_data_real ( 1:NDim, Counter ) = XYZ ( K, 1:NDim ) ; 
          dset_data_int  ( Counter, 1:NDim ) = ID_PETSc ( Global_Petsc_Num ( K ), 1:NDim ) ; ! Equation number of nodes on Paraview to retrieve solution

          If ( Counter == NJ_Para ) Exit ;
        End If ;
      End Do ;

    ! Write the dataset.
    ParamPartition%IntM ( 5, 2) = NJ_Para ;
    data_dims(2) = NJ_Para ;
    data_dims(1) = 3 ; ! NDim ;
    Call h5dwrite_f(dset_id_XYZ_PV, H5T_NATIVE_DOUBLE, dset_data_real, data_dims, ErrPTC)

    ! Equation number of nodes on this rank used to extract data of each node for Paraview
    ! Remark: We already have equation numbers of this rank in Indices vector; However, we cannot use this vector to scatter solution for Paraview, because Paraview does not support Spectral Elements. Hence, we do not need to scatter all solution.
    NE_Para = Count ( dset_data_int ( :, : ) /= 0_Lng ) ;
    ParamPartition%IntM ( 5, 1) = NE_Para ;
    dims(1) = NE_Para ;
    dims(2) = 1 ;
    Call h5screate_simple_f(2, dims, dspace_id_IS_PV, ErrPTC) ;
    Call h5dcreate_f(id_Data_PV, "IndexSet",     H5T_NATIVE_INTEGER, dspace_id_IS_PV,   dset_id_IS_PV,   ErrPTC) ;

    dset_data_int = dset_data_int - 1 ; ! Equation numbers start from 0 in PETSc

    ! Write the dataset.
    data_dims(1) = NJ_Para ;
    data_dims(2) = NDim ;
    Call h5dwrite_f(dset_id_Node_PV, H5T_NATIVE_INTEGER, dset_data_int, data_dims, ErrPTC)

    Allocate ( dset_data_ISet ( NE_Para ) );

    Counter = 0 ;
      Do I = 1, NJ_Para ;
        Do J = 1, NDim ;
          If ( dset_data_int ( I, J) /= -1 ) Then ; 
            Counter = Counter + 1 ;
            dset_data_ISet ( Counter ) = dset_data_int ( I, J ) ;
          End If ; 
        End Do ;
      End Do ;

      If ( Counter /= NE_Para ) Then ; 
        Write (*,*)"There is a mistake in the Index setting portion of the output code.", Counter, NE_Para
        Stop ;
      End If ;

    ! Write the dataset for index setting of the equation numbers
    data_dims(1) = NE_Para ;
    data_dims(2) = 1 ;
    Call h5dwrite_f(dset_id_IS_PV, H5T_NATIVE_INTEGER, dset_data_ISet, data_dims, ErrPTC)

    DeAllocate ( dset_data_real, dset_data_int, dset_data_ISet ) ;

    Write (*,*)"Write Data for Paraview - Connectivity..." ;

    ! Element Connectivities for Paraview
    dims(2) = ConnSizePara ;
    dims(1) = 1 ;

    Allocate ( dset_data_int ( dims(1), dims(2)) );

    Counter = 1_Lng ;
    dset_data_int (:,:) = 0_Shrt ;

      DO IEL = 1, NEL ;

        If ( EPart ( IEL ) == IParts ) Then ;
          EType = ElT ( IEl ) ; ! Element Type

            ! determinig NNode and cell types based on 
!            // Topologies
!            #define XDMF_NOTOPOLOGY     0x0     - 0
!            #define XDMF_POLYVERTEX     0x1     - 1
!            #define XDMF_POLYLINE       0x2     - 2
!            #define XDMF_POLYGON        0x3     - 3
!            #define XDMF_TRI            0x4     - 4  !!
!            #define XDMF_QUAD           0x5     - 5  !!
!            #define XDMF_TET            0x6     - 6
!            #define XDMF_PYRAMID        0x7     - 7
!            #define XDMF_WEDGE          0x8     - 8
!            #define XDMF_HEX            0x9     - 9  !!
!            #define XDMF_EDGE_3         0x0022  - 
!            #define XDMF_TRI_6          0x0024  - 36  !!
!            #define XDMF_QUAD_8         0x0025  - 37  !!
!            #define XDMF_QUAD_9         0x0023  - 
!            #define XDMF_TET_10         0x0026  - 38  !!
!            #define XDMF_PYRAMID_13     0x0027  - 
!            #define XDMF_WEDGE_15       0x0028  - 
!            #define XDMF_WEDGE_18       0x0029  - 
!            #define XDMF_HEX_20         0x0030  - 48 !!
!            #define XDMF_HEX_24         0x0031  - 
!            #define XDMF_HEX_27         0x0032  - 
!            #define XDMF_MIXED          0x0070  - 
!            #define XDMF_2DSMESH        0x0100  - 
!            #define XDMF_2DRECTMESH     0x0101  - 
!            #define XDMF_2DCORECTMESH   0x0102  - 
!            #define XDMF_3DSMESH        0x1100  - 
!            #define XDMF_3DRECTMESH     0x1101  - 
!            #define XDMF_3DCORECTMESH   0x1102  - 

            IF      ( EType == El2d4NSldPN  ) Then ;  NNode = 4_Tiny ; X_CellType = 5_Tiny;     ! 4 node  -2D - Solid - PLANE STRESS
            Else If ( EType == El2d4NPMLPN  ) Then ;  NNode = 4_Tiny ; X_CellType = 5_Tiny;     ! 4 node  -2D - PML - PLANE STRESS
            Else If ( EType == El2d8NSldPN  ) Then ;  NNode = 8_Tiny ; X_CellType = 37_Tiny;    ! 8 node  -2D - Solid - PLANE STRAIN
            Else If ( EType == El2d8NPMLPN  ) Then ;  NNode = 8_Tiny ; X_CellType = 37_Tiny;    ! 8 node  -2D - PML - PLANE STRAIN
            Else If ( EType == El2d6NSldPN  ) Then ;  NNode = 6_Tiny ; X_CellType = 36_Tiny;    ! 6 node  -2D - Solid - PLANE STRAIN
            Else If ( EType == El2d6NPMLPN  ) Then ;  NNode = 6_Tiny ; X_CellType = 36_Tiny;    ! 6 node  -2D - PML - PLANE STRAIN
            Else If ( EType == El2d3NSldPN  ) Then ;  NNode = 3_Tiny ; X_CellType = 4_Tiny;     ! 3 node  -2D - Solid - PLANE STRAIN
            Else If ( EType == El2d3NPMLPN  ) Then ;  NNode = 3_Tiny ; X_CellType = 4_Tiny;     ! 3 node  -2D - PML - PLANE STRAIN
            Else If ( EType == SEl2d9NSldPN ) Then ;  NNode = 8_Tiny ; X_CellType = 37_Tiny;    ! 9 node  -2D - Solid - PLANE STRAIN
            Else If ( EType == SEl2d9NPMLPN ) Then ;  NNode = 8_Tiny ; X_CellType = 37_Tiny;    ! 9 node  -2D - NPML - PLANE STRAIN    !?? using 8 noded element for spectral 9 noded elements
            Else If ( EType == SEl2d9NMPMLPN) Then ;  NNode = 8_Tiny ; X_CellType = 37_Tiny;    ! 9 node  -2D - MPML - PLANE STRAIN
            Else If ( EType == SEl3d27NSld  ) Then ;  NNode = 20_Tiny; X_CellType = 48_Tiny;    ! 27 node -3D - Solid
            Else If ( EType == SEl3d27NMPML ) Then ;  NNode = 20_Tiny; X_CellType = 48_Tiny;    ! 27 node -3D - MPML
            Else If ( EType == SEl2d7NSldPN ) Then ;  NNode = 6_Tiny ; X_CellType = 36_Tiny;    ! 7 node -2D - solid - plane strain
            Else If ( EType == SEl2d7NMPMLPN) Then ;  NNode = 6_Tiny ; X_CellType = 36_Tiny;    ! 7 node -2D - MPML - plane strain
            Else If ( EType == SEl2d9NSldSH ) Then ;  NNode = 8_Tiny ; X_CellType = 37_Tiny;    ! 9 node - 2D - solid ( quadrilateral ) - Spectral Element - SH waves
            Else If ( EType == SEl2d9NPMLSH ) Then ;  NNode = 8_Tiny ; X_CellType = 37_Tiny;    ! 9 node - 2D - PML   ( quadrilateral ) - Spectral Element - SH waves
            Else If ( EType == SEl3d8NSld   ) Then ;  NNode = 8_Tiny ; X_CellType = 9_Tiny;     ! 8 node -3D - Solid
            Else If ( EType == SEl3d8NMPML  ) Then ;  NNode = 8_Tiny ; X_CellType = 9_Tiny;     ! 8 node -3D - MPML
            End If ;

          dset_data_int ( 1, Counter)                  = X_CellType ;

            If ( NNode /= 20_Tiny ) Then ;
              dset_data_int ( 1, Counter+1:Counter+NNode ) = Local_PETSc_Num ( INOD ( 1:NNode, IEL ), IParts ) - 1_Lng ;
            Else If ( NNode == 20_Tiny ) Then ;  ! Node numbering of the 20-node element in Paraview is similar to Ansys, so we have to change the node numbering
              ForAll( I = 1:12   ) dset_data_int ( 1, Counter+I ) = Local_PETSc_Num ( INOD ( I,   IEL ), IParts ) - 1_Lng ;
              ForAll( I = 13:16  ) dset_data_int ( 1, Counter+I ) = Local_PETSc_Num ( INOD ( I+4, IEL ), IParts ) - 1_Lng ;
              ForAll( I = 17:20  ) dset_data_int ( 1, Counter+I ) = Local_PETSc_Num ( INOD ( I-4, IEL ), IParts ) - 1_Lng ;
            End If ;

          Counter = Counter + Int4(NNode + 1_Lng) ;

        End If ;

      End Do ;

    If ( Counter-1_Lng /= ConnSizePara ) Then ; 
      Write (*,*)"There is a mistake in the connectivity portion of the output code.", Counter, ConnSizePara
      Stop ;
    End If ;

    data_dims(2) = ConnSizePara ;
    data_dims(1) = 1 ;
    Call h5dwrite_f(dset_id_Cnn_PV, H5T_NATIVE_INTEGER, dset_data_int, data_dims, ErrPTC ) ;

    DeAllocate ( dset_data_int ) ;

    ! Constraints
    Write (*,*)"Write Node constraints ..." ;
    Allocate ( dset_data_int( NJ_Rank ( IParts ), NDOF ) ) ;
    Counter = 0_Lng ;
    dset_data_int (:,:) = 0_Lng ;

    ! Writing the Local Equation Numbers (Local ID)
      Do IJ = 1, NJ_Rank ( IParts ) ;
        dset_data_int ( IJ, 1:NDOF ) = ID_Local ( IJ, 1:NDOF ) ;
      End Do ;

    data_dims(1) = NJ_Rank ( IParts ) ;
    data_dims(2) = NDOF ;
    Call h5dwrite_f( dset_id_Cnt, H5T_NATIVE_INTEGER, dset_data_int, data_dims, ErrPTC ) ;

    DeAllocate ( dset_data_int ) ;

! - Basic Data -------------------------------------------------------------------------------------------------------------------------------------
    ! Basic data
    UnFile = Un_Out ;
    Write (*,*)"Writing Basic Data ..." ;
    Write (Unit = UnFile, FMT = "(3(I3,2X),1(I5,2X),6(I19,2X))", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) NDOF, NDim, MaxNNode,      NGroup,         NEL_Rank ( IParts ), NJ_Rank ( IParts ), NJ, NEQM, NEQMTotal, NEqM_Mapping ;

    ! Available Load Cases
    Write(*,*)"writing Load Cases ..." ;
    Write (Unit = UnFile, FMT = "(5(I3,2X),'Load Case'        )", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) ( LoadC ( I ), I = 1, 5 ) ;

    ! Writing Param
    Write(*,*)"Writing Param Array ..." ;
    Write (Unit = UnFile, FMT = "(<6*7>(I19,1X),'Param%IntM')", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) (( ParamPartition%IntM  ( I, J), I = 1, 7), J = 1, 6) ;
    Write (Unit = UnFile, FMT = "(<6*7>(E31.23E3,1X),'Param%RealM')", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) (( ParamPartition%RealM  ( I, J), I = 1, 7), J = 1, 6) ;

! - Geometry & Model -------------------------------------------------------------------------------------------------------------------------------
    Write (*,*)"Writing Geometry ..." ;

!    ! Number of elements of each group
!    Write(*,*)"Writing Group numbers ..." ;
!      DO I = 1, NGroup ;
!        If ( I == 1_Lng ) then ;
!          Write (Unit = UnFile, FMT = "(I19,2X,I5,2X,'GROUP NUMBER, NUMBER OF ELEMENTS OF THIS GROUP , ELEMENT Type' )", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) NEG ( I, 1 ), NEG ( I, 2 ) ; ! NUMBER OF ELEMENTS OF THIS GROUP , ELEMENT Type
!        Else ;
!          Write (Unit = UnFile, FMT = "(I19,2X,I5,2X                                                                 )", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) NEG ( I, 1 ), NEG ( I, 2 ) ; ! NUMBER OF ELEMENTS OF THIS GROUP , ELEMENT Type
!        End If ;
!      End Do ;

    If ( LoadC (5) /= 0_Tiny ) Then ;
      Write (*,*)"Writing PML teritory ..." ;
      ! PML teritory
        DO J = 1, 2 ! See Related PML Subroutines
          Write (Unit = UnFile, FMT = "(<2*NDim>(E31.23E3,2X),'PML TERITORY')", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) (PML_DIM ( I, J ), I = 1, 2 * NDim ) ; 
        End Do ;
    End If ;

    ! Material properties of each material
    Write (*,*)"Writing Material Properties ..." ;
    UnFile = Un_OutMat ;
      DO I = 1, NMat ;
        Write (Unit = UnFile, FMT = "(<NPM>(E31.23E3,2X) )", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 )( PMat ( I, J ), J = 1, NPM ) ; 
      End Do ;

! - Loads -------------------------------------------------------------------------------------------------------------------------------------------
    Write (*,*)"Writing Load types ..." ;
    UnFile = Un_Out ;

    ! Body Force Load
      If ( LoadC (1) /= 0_Tiny ) Then ;
        Write(*,*)"Reading Body Forces ..." ;
          ! Body force load
          DO I = 1, ParamPartition%IntM( 1, 1) ; ! ParamPartition%IntP( 1, 1)=NBLD: Number of body force load.   ParamPartition%IntP( 1, 2) = NPBL
            If ( I == 1_Lng ) then ;
              Write (Unit = UnFile, FMT = "(<ParamPartition%IntM( 1, 2)>(E31.23E3,2x),'Body Force Load')", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) ( PBLD ( I, J ) , J = 1, ParamPartition%IntM( 1, 2) ) ;
            Else ;
              Write (Unit = UnFile, FMT = "(<ParamPartition%IntM( 1, 2)>(E31.23E3,2x)                  )", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) ( PBLD ( I, J ) , J = 1, ParamPartition%IntM( 1, 2) ) ;
            End If ;
          End Do ;

!        ! Load Types
!        LDSET = 0_Shrt ;
!        LTEL_Rank = 0_Lng ;
!        Lo_EL_Num = 0_Lng ;
!        I2 = 0_Lng ;
!
!          Do IEL = 1, NEL ;
!            If ( EPart ( IEL ) == IParts ) Lo_EL_Num = Lo_EL_Num + 1_Lng ;
!
!            If ( EPart ( IEL ) == IParts .AND. LDSET == 0_Shrt ) Then ;
!              ILType = LTEL ( IEL ) ;
!              I1 = Lo_EL_Num ;
!              LDSET = 1_Shrt ;
!            End If ;
!
!            If ( EPart ( IEL ) == IParts .AND. LTEL ( IEL ) /= ILType ) Then ;
!              I2 =  Lo_EL_Num - 1_Lng ;
!              INC = 1_Lng ;
!
!              LTEL_Rank ( 1, LDSET ) = I1  ;
!              LTEL_Rank ( 2, LDSET ) = I2  ;
!              LTEL_Rank ( 3, LDSET ) = INC ;
!              LTEL_Rank ( 4, LDSET ) = ILType ;
!              
!              ILType = LTEL ( IEL ) ;
!              LDSET = LDSET + 1_Shrt ;
!              I1 = Lo_EL_Num ;
!              I2 = 0_Lng ;
!            Else If ( I2 == 0_Lng .AND. IEL == NEL ) Then ;
!              I2 =  Lo_EL_Num ;
!              INC = 1_Lng ;
!
!              LTEL_Rank ( 1, LDSET ) = I1  ;
!              LTEL_Rank ( 2, LDSET ) = I2  ;
!              LTEL_Rank ( 3, LDSET ) = INC ;
!              LTEL_Rank ( 4, LDSET ) = ILType ;
!            End If ;
!          End Do ;
!
!        Write (Unit = UnFile, FMT = "(I10,2x,'LOAD SETS')", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) LDSET ;
!          DO I = 1, LDSET ;
!            Write (Unit = UnFile, FMT = "(<4>(I19,2X))", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) ( LTEL_Rank ( K, I ), K = 1, 4 ) ;
!          End Do ;

      End If ;

    ! Static Pressure
      If ( LoadC (2) /= 0_Tiny ) Then ;
        Write(*,*)"Writing Pressure Loads ..." ;
          Do I = 1, Param%IntM( 2, 4) ;  ! NIDBC
            If ( EPart ( IDBC ( I, 1 ) ) == IParts ) Then ;
              Write (Unit = UnFile, FMT = "((I19,2X),<2*NDim>(I3,2X))", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) Local_ElN ( IDBC ( I, 1 ) ),  ( IDBC ( I, J ), J = 2, 2 * NDim +1) ;
            End If ;
          End Do ;
      End If ;

    ! Joint load of every case - Concentrated loads
      If ( LoadC ( 3 ) /= 0 ) Then ;
        Write(*,*)"Writing Joint Loads ..." ;
          DO I = 1, NLN_Rank ;  ! ParamPartition%IntP( 3, 1) = NLN_Rank 
            If ( Local_PETSc_Num ( JLoad ( I ) , IParts ) /= 0_Lng ) Then ; 
              If ( Nodes%Rep ( JLoad ( I ) ) == 1_Shrt ) Then ;
                Write (Unit = UnFile, FMT = "(I19,2X,<NDim>(E31.23E3,2X),'JOINT LOADS')", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 )  Local_PETSc_Num ( JLoad ( I ), IParts ), ( PLoad( J, I ), J = 1, NDim ) ;
              Else If ( Nodes%Rep ( JLoad ( I ) ) /= 1_Shrt .AND. IParts == MinVal ( Nodes%Locs( JLoad ( I ),1:Nodes%Rep ( JLoad ( I ) ) ) ) ) Then ; ! Takes care of nodes on the boundary of the partitions. So, we do not have a repeatation of the loads
                Write (Unit = UnFile, FMT = "(I19,2X,<NDim>(E31.23E3,2X),'JOINT LOADS')", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 )  Local_PETSc_Num ( JLoad ( I ), IParts ), ( PLoad( J, I ), J = 1, NDim ) ;
              End If ;
            End If ; 
          End Do ;
      EndIF ;

    ! Supports Displacements
      If ( LoadC ( 4 ) /= 0 ) Then ;
        Write(*,*)"Writing Support Displacements ..." ;
          DO I = 1, NSND_Rank ;
            If ( Local_PETSc_Num ( UDis( NDOF + 1, I ), IParts ) /= 0_Lng ) Then ; 
              If ( Nodes%Rep ( UDis( NDOF + 1, I ) ) == 1_Shrt ) Then ;
                Write (Unit = UnFile, FMT = "(<NDOF+1>(E31.23E3,2X))", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 )  Local_PETSc_Num ( UDis( NDOF + 1, I ), IParts ), ( UDis ( K, I ), K = 1, NDOF ) ;
              Else If ( Nodes%Rep ( UDis( NDOF + 1, I ) ) /= 1_Shrt .AND. IParts == MinVal ( Nodes%Locs( UDis( NDOF + 1, I ),1:Nodes%Rep ( UDis( NDOF + 1, I ) ) ) ) ) Then ;
                Write (Unit = UnFile, FMT = "(<NDOF+1>(E31.23E3,2X))", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 )  Local_PETSc_Num ( UDis( NDOF + 1, I ), IParts ), ( UDis ( K, I ), K = 1, NDOF ) ;
              End If ;
            End If ; 
          End Do ;
      End If ;

    Write(*,*)"Writing Dynamic Loads ..." ;

    ! Dynamic loads
      If ( LoadC (5) /= 0_Tiny ) Then ;

        If ( ParamPartition%IntM( 5, 3 ) == 1 ) Then ;    ! RICKER PULSE or sine function for dynamic pressure
          Write(*,*)"Writing Function Type ..." ;
            Do I = 1, Param%IntM( 2, 4) ; ! NIDBC
              If ( EPart ( IDBC ( I, 1 ) ) == IParts ) Then ;
                Write (Unit = UnFile, FMT = "((I19,2X),<2*NDim>(I3,2X))", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) Local_ElN ( IDBC ( I, 1 ) ),  ( IDBC ( I, J ), J = 2, 2 * NDim +1) ;
              End If ;
            End Do ;

        Else If ( ParamPartition%IntM( 5, 3 ) == 2_Tiny ) Then ;   ! Base acceleration

!          Write(*,*)"Writing Base Acceleration ..." ;
!            DO I = 1, ParamPartition%IntL( 7, 1) ; ! NAStep = ParamPartition%IntP( 7, 1)
!              Write (Unit = UnFile, FMT = "(<NDim>(E31.23E3,2X),'ACCELERATION' )", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) ( BACL ( I, K ), K = 1, NDim ) ;
!            End Do ;

        Else If ( ParamPartition%IntM( 5, 3 ) == 3_Tiny ) Then ;    ! Domain Reduction Method (DRM)

!          Write(*,*)"Writing wave information ..." ;
!          Write (Unit = UnFile, FMT = "(5(E28.20E3,1x))", ADVANCE = 'Yes', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) InciWave(1), InciWave(2), InciWave(3), InciWave(4), InciWave(5) ;       ! Theta, Omega, amplitude, alpha1, alpha2

          ! writing the node numbers on the DRM boundary
          Write(*,*)"Writing DRM boundary nodes ..." ;
            Do I = 1, Param%IntM( 7, 1) ; ! NNBndry_DRM
              Node = NoBndry_DRM ( I ) ;
              If ( ( Nodes%Rep (Node) == 1_Shrt .AND. Any ( Nodes%Locs ( Node, : ) == IParts ) ) .OR. ( Nodes%Rep (Node) /= 1_Shrt .AND. IParts == MinVal ( Nodes%Locs( Node, :), MASK = Nodes%Locs( Node, :) > 0_Shrt ) ) ) Write (Unit = UnFile, FMT = "(1(I19,1X))", ADVANCE = 'No', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) Local_PETSc_Num ( Node, IParts ) ;
            End Do ;
          Write (Unit = UnFile, FMT = "()", ADVANCE = 'Yes', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) ;

          ! Writing the node numbers on the DRM neighbor
          Write(*,*)"Writing DRM neighbor nodes ..." ;
            Do I = 1, Param%IntM( 7, 2) ! NNLayer_DRM
              Node = NoLayer_DRM ( I ) ;
              If ( ( Nodes%Rep (Node) == 1_Shrt .AND. Any ( Nodes%Locs ( Node, : ) == IParts ) ) .OR. ( Nodes%Rep (Node) /= 1_Shrt .AND. IParts == MinVal ( Nodes%Locs( Node, :), MASK = Nodes%Locs( Node, :) > 0_Shrt ) ) ) Write (Unit = UnFile, FMT = "(1(I19,1X))", ADVANCE = 'Yes', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) Local_PETSc_Num ( Node, IParts ) ;
            End Do ;
          !Write (Unit = UnFile, FMT = "()", ADVANCE = 'Yes', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) ;

        End If ;

      End If ;

!    ! writes down step numbers for full results
!    Write(*,*)"Writing Step Numbers ..." ;
!    Write (Unit = UnFile, FMT = "(<ParamPartition%IntL( 5, 2)>(I19,2X),'Required STEPS')", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) ( STEP ( I ), I = 1, ParamPartition%IntL( 5, 2) ) ; ! NTST

      ! writes down node numbers and equation numbers of nodes in which history of Displacement is required
      If ( NNDH_Rank /= 0_Shrt ) Then ;
        Write(*,*)"Writing Node Numbers for displacements ..." ;
        Write (Unit = UnFile, FMT = "(<NNDH_Rank>(I19,2X),5x,<NNDH_Rank>(I19,2X),'Required NODES FOR DisPLACEMENTS, Global application num & Global PETSc num')", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) ( NodeDis ( I ), I = 1, NNDH_Rank ), ( Global_PETSc_Num ( NodeDis ( I ) ), I = 1, NNDH_Rank ) ;
        Write (Unit = UnFile, FMT = "(<NNDH_Rank*NDim>(I19,2X),'Required Equation numbers FOR DisPLACEMENTS')"                                                       , ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) ( EqDis ( I ), I = 1, NNDH_Rank * NDim ) ;
      End If ;

      ! writes down node numbers and equation numbers of nodes in which history of velocity is required
      If ( NNVH_Rank /= 0_Shrt ) Then ;
        Write(*,*)"Writing Node Numbers for velocity ..." ;
        Write (Unit = UnFile, FMT = "(<NNVH_Rank>(I19,2X),5x,<NNVH_Rank>(I19,2X),'Required NODES FOR velocity, Global application num & Global PETSc num')", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) ( NodeVel ( I ), I = 1, NNVH_Rank ), ( Global_PETSc_Num ( NodeVel ( I ) ), I = 1, NNVH_Rank ) ;
        Write (Unit = UnFile, FMT = "(<NNVH_Rank*NDim>(I19,2X),'Required Equation numbers FOR velocity')"                                                       , ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) ( EqVel ( I ), I = 1, NNVH_Rank * NDim ) ;
      End If ;

      ! writes down node numbers and equation numbers of nodes in which history of acceleration is required
      If ( NNAH_Rank /= 0_Shrt ) Then ;
        Write(*,*)"Writing Node Numbers for acceleration ..." ;
        Write (Unit = UnFile, FMT = "(<NNAH_Rank>(I19,2X),5x,<NNAH_Rank>(I19,2X),'Required NODES FOR acceleration, Global application num & Global PETSc num')", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) ( NodeAcc ( I ), I = 1, NNAH_Rank ), ( Global_PETSc_Num ( NodeAcc ( I ) ), I = 1, NNAH_Rank ) ;
        Write (Unit = UnFile, FMT = "(<NNAH_Rank*NDim>(I19,2X),'Required Equation numbers FOR acceleration')"                                                       , ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) ( EqAcc ( I ), I = 1, NNAH_Rank * NDim ) ;
      End If ;

! <><><><><><><><><><><><><><><><><><><><><>
    ! Application Numbering
    Write (*,*)"Writing application numberting ..." ;

    Allocate ( App_Numbers ( NEqM_Mapping ), PETSc_Numbers ( NEqM_Mapping ), Indices ( NEQM ) ) ;
    App_Numbers   (:) = 0_Lng ;
    PETSc_Numbers (:) = 0_Lng ;

    ! Obtaining the Application Numbering and PETSc Oredering for "Application Ordering (AO)" in PETSc.
    Counter = 0_Lng ;
      If ( IParts == 1_Shrt ) then ;
        LowerLimit = 1_Lng ;
      Else ;
        LowerLimit = NNodeRank ( IParts - 1_Shrt ) + 1_Lng ;
      End If ;

    UpperLimit = NNodeRank ( IParts ) ;

      Do IJ = 1_Lng, NJ ;

        PETScNodeNumber = Global_PETSc_Num ( IJ )

          If ( LowerLimit <= PETScNodeNumber .AND. PETScNodeNumber <= UpperLimit ) Then ;

            Do IDOF = 1, NDOF ;

              If ( ID_Application ( IJ, IDOF ) /= 0_Lng ) Then ;
                Counter = Counter + 1_Lng ;
                  If ( Counter > NEqM_Mapping ) Then ; !??
                    Write(*      , "('Number of NEqM_Mapping equations of this rank is greater than the estimated equation - Check Subroutine <OutputParMetis>', 3(I19,2x))" ) IParts, IJ, NEqM_Mapping ;
                    Write( UnFile, "('Number of NEqM_Mapping equations of this rank is greater than the estimated equation - Check Subroutine <OutputParMetis>', 3(I19,2x))" ) IParts, IJ, NEqM_Mapping ;
                    !#Call BEEP_FAIL ;
                    Write(*, Fmt_FL) ;  Write(UnInf, Fmt_FL) ; Write(*, Fmt_End) ; Read(*,*) ;  Return ;
                  End If ;
                App_Numbers   ( Counter ) = ID_Application ( IJ                     , IDOF ) ; ! myapp
                PETSc_Numbers ( Counter ) = ID_PETSc       ( Global_PETSc_Num ( IJ ), IDOF ) ; ! mypetsc
              End If ;
            End Do ;

          End If ;

      End Do ;


!    ! Obtaining the Application Numbering for "Application Ordering (AO)" in PETSc
!    Counter = 0_Lng ;
!      Do IJ = 1, NJ ;
!        If ( Local_PETSc_Num ( IJ, IParts ) /= 0_Lng .AND. Nodes%Rep ( IJ ) == 1_Shrt ) Then ;   !!??
!          Do IDOF = 1, NDOF ;
!            If ( ID_Application ( IJ, IDOF ) /= 0_Lng ) Then ;
!              Counter = Counter + 1_Lng ;
!              !Write(*,*)"debug-count",counter, ij
!                If ( Counter > NEqM_Mapping ) Then ; !!!!!!??
!                  Write(*       , "('Number of Non_Neighbor equations of this rank is greater than the estimated equation', I10)" ) ; ! Rank ; 
!                  Write( UnFile, "('Number of Non_Neighbor equations of this rank is greater than the estimated equation', I10)" ) ; ! Rank ; 
!                  !#Call BEEP_FAIL ;
!                  Write(*, Fmt_FL) ;  Write(UnInf, Fmt_FL) ; Write(*, Fmt_End) ; Read(*,*) ;  Return ;
!                End If ;
!              App_Numbers ( Counter ) = ID_Application ( IJ, IDOF ) ; ! myapp
!            End If ;
!          End Do ;
!        End If ;
!      End Do ;
!      
!    ! Obtaining the PETSc Numbering for "Application Ordering (AO)" in PETSc
!    Counter = 0_Lng ;
!      DO IJ = 1, NJ ;
!        If ( Local_PETSc_Num ( IJ, IParts ) /= 0_Lng .AND. Nodes%Rep ( IJ ) == 1_Shrt ) Then ;
!          Do IDOF = 1, NDOF ;
!            If ( ID_PETSc ( Global_PETSc_Num ( IJ ), IDOF ) /= 0_Lng ) Then ;
!              Counter = Counter + 1_Lng ;
!                If ( Counter > NEqM_Mapping ) Then ; !!!!!!??
!                  Write(*       , "('Number of Non_Neighbor equations of this rank is greater than the estimated equation', I10)" ) ; ! Rank ; 
!                  Write( UnFile, "('Number of Non_Neighbor equations of this rank is greater than the estimated equation', I10)" ) ; ! Rank ; 
!                  !#Call BEEP_FAIL ;
!                  Write(*, Fmt_FL) ;  Write(UnInf, Fmt_FL) ; Write(*, Fmt_End) ; Read(*,*) ;  Return ;
!                End If ;
!              PETSc_Numbers ( Counter ) = ID_PETSc ( Global_PETSc_Num ( IJ ), IDOF ) ;  ! mypetsc
!            End If ;
!          End Do ;
!        End If ;
!      End Do ;

    ! Mapping from the local numbering to Global PETSc numbering
      Do IJ = 1, NJ ;
        If ( Local_PETSc_Num ( IJ, IParts ) /= 0_Lng ) Then ;
          Do IDOF = 1, NDOF ;
            If ( ID_Local ( Local_PETSc_Num ( IJ, IParts ), IDOF ) /= 0_Lng ) Indices ( ID_Local ( Local_PETSc_Num ( IJ, IParts ), IDOF ) ) = ID_PETSc ( Global_PETSc_Num ( IJ ), IDOF)  ; 
          End Do ;
        End If ;
      End Do ;

    ! Equation numbers start form 0 in PETSc
    App_Numbers   (:) = App_Numbers   (:) - 1_Lng ;
    PETSc_Numbers (:) = PETSc_Numbers (:) - 1_Lng ;
    Indices       (:) = Indices       (:) - 1_Lng ;

    ! Writing down the data for node numbering mapping
    Write(*,*)"Writing Index Sets ..." ;
    ! Application Numbering 
    data_dims(1) = 1_Shrt ;
    data_dims(2) = NEqM_Mapping;

    Call h5dwrite_f( dset_id_App, H5T_NATIVE_INTEGER, App_Numbers, data_dims, ErrPTC ) ;

    ! PETSc Numbering 
    data_dims(1) = 1 ;
    data_dims(2) = NEqM_Mapping;
    Call h5dwrite_f( dset_id_PTC, H5T_NATIVE_INTEGER, PETSc_Numbers, data_dims, ErrPTC ) ;

    ! Index set (Indices)
    data_dims(1) = 1 ;
    data_dims(2) = NEqM ;
    Call h5dwrite_f( dset_id_Ind, H5T_NATIVE_INTEGER, Indices, data_dims, ErrPTC ) ;

! Since we are using a fixed number of nonzeros for all matrices, defined in the main, this part is not required. Uncomment and modify the related subroutines.
!    ! Number of non-zero entris of PETSc objects
!    Write (*,*)"Writing number of nonzero elements ..." ;
!
!      If ( IParts == 1_Shrt ) then ;
!        LowerLimit = 1_Lng ;
!      Else ;
!        LowerLimit = NEqRank ( IParts - 1_Shrt ) + 1_Lng ;
!      End If ;
!
!    UpperLimit = NEqRank ( IParts ) ;
!
!    UnFile = Un_OutStiff ;
!    Write(*,*)"Writing Number of Non-Zeros of Stiffness ..." ;
!    Write (Unit = UnFile, FMT = "(<NEqM_Mapping>(I19,2X),'NNZ_Stiff_Diagoanl    ')", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) ( D_NNZ_Stiff ( I ), I = LowerLimit, UpperLimit ) ;
!    Write (Unit = UnFile, FMT = "(<NEqM_Mapping>(I19,2X),'NNZ_Stiff Off diagonal')", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) ( O_NNZ_Stiff ( I ), I = LowerLimit, UpperLimit ) ;
!
!    UnFile = Un_OutDamp ;
!    Write(*,*)"Writing Number of Non-Zeros of Damping ..." ;
!    Write (Unit = UnFile, FMT = "(<NEqM_Mapping>(I19,2X),'NNZ_Damp  Diagonal    ')", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) ( D_NNZ_Damp  ( I ), I = LowerLimit, UpperLimit ) ;
!    Write (Unit = UnFile, FMT = "(<NEqM_Mapping>(I19,2X),'NNZ_Damp  Off diagonal')", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) ( O_NNZ_Damp  ( I ), I = LowerLimit, UpperLimit ) ;
!
!    UnFile = Un_OutMass ;
!    Write(*,*)"Writing Number of Non-Zeros of Mass ..." ;
!    Write (Unit = UnFile, FMT = "(<NEqM_Mapping>(I19,2X),'NNZ_Mass  Diagonal    ')", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) ( D_NNZ_Mass  ( I ), I = LowerLimit, UpperLimit ) ;
!    Write (Unit = UnFile, FMT = "(<NEqM_Mapping>(I19,2X),'NNZ_Mass  Off diagonal')", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) ( O_NNZ_Mass  ( I ), I = LowerLimit, UpperLimit ) ;

    DeAllocate ( ID_Local, App_Numbers, PETSc_Numbers, Indices ) ;

    ! - Closing the output file ---------------------------------------------------------------------------------------------------------------------
    UnFile =  UN_Out ;
    Close ( Unit = UnFile, Status = 'KEEP', ERR =  1002, IOSTAT = IO_File ) ;

!    UnFile =  UN_OutXYZ ;
!    Close ( Unit = UnFile, Status = 'KEEP', ERR =  1002, IOSTAT = IO_File ) ;

!    UnFile =  UN_OutCnn ;
!    Close ( Unit = UnFile, Status = 'KEEP', ERR =  1002, IOSTAT = IO_File ) ;

!    UnFile =  UN_OutCnt ;
!    Close ( Unit = UnFile, Status = 'KEEP', ERR =  1002, IOSTAT = IO_File ) ;

!    UnFile =  Un_OutStiff ;
!    Close ( Unit = UnFile, Status = 'KEEP', ERR =  1002, IOSTAT = IO_File ) ;

!    UnFile =  Un_OutDamp ;
!    Close ( Unit = UnFile, Status = 'KEEP', ERR =  1002, IOSTAT = IO_File ) ;

!    UnFile =  Un_OutMass ;
!    Close ( Unit = UnFile, Status = 'KEEP', ERR =  1002, IOSTAT = IO_File ) ;

!    UnFile =  Un_OutApp ;
!    Close ( Unit = UnFile, Status = 'KEEP', ERR =  1002, IOSTAT = IO_File ) ;

    UnFile = Un_OutMat ;
    Close ( Unit = UnFile, Status = 'KEEP', ERR =  1002, IOSTAT = IO_File ) ;

    ! Close the dataset.
    Call h5dclose_f( dset_id_XYZ, ErrPTC) ;
    Call h5dclose_f( dset_id_Cnn, ErrPTC) ;
    Call h5dclose_f( dset_id_Cnt, ErrPTC) ;
    Call h5dclose_f( dset_id_App, ErrPTC) ;
    Call h5dclose_f( dset_id_PTC, ErrPTC) ;
    Call h5dclose_f( dset_id_Ind, ErrPTC) ;

    ! Close the HDF5 file.
    CALL h5fclose_f( id_Geo,     ErrPTC) ;
    CALL h5fclose_f( id_Geo_PV,  ErrPTC) ;
    CALL h5fclose_f( id_Data_PV, ErrPTC) ;

    ! Shut down HDF5.
    CALL h5close_f(ErrPTC) ;

  End Do ;

  If ( NEqM_Verifier /= NEqMTotal ) Then ;
    Write(*,*)" MAJOR PROBLEM - CHECK NEQM_MAPPING", NEqMTotal, NEqM_Mapping ;
    Stop ;
  End If ;


! Deallocating 
DEAllocate( EqDis, EqVel, EqAcc, NodeDis, NodeVel, NodeAcc, STAT = ERR_DeAlloc ) ;
  IF ( ERR_DeAlloc /= 0 ) Then ;
    Write (*, Fmt_DEALLCT) ERR_DeAlloc ;  Write (UnInf, Fmt_DEALLCT) ERR_DeAlloc ;
    Write(*, Fmt_FL) ;  Write(UnInf, Fmt_FL) ; Write(*, Fmt_End) ; Read(*,*) ;  STOP ;
  End If ;


!/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_
! Write sensorList: List_NDAN_rank ( IParts , NNDH_Rank )

   UnFile = UN_Out_List_sensor ;

   Open ( Unit = UnFile, FILE = TRIM(ModelName)//'_'//Trim(AdjustL(IndexSize))//'.sensorList', ERR =  1001, IOSTAT = IO_File, ACCESS = 'SEQUENTIAL', ACTION ='Write', ASYNCHRONOUS = 'NO', BLANK = 'NULL', BLOCKSIZE = 0, DEFAULTFILE = TRIM(OutDir), DisPOSE = 'KEEP', FORM = 'FORMATTED', POSITION = 'ASIS', STATUS ='REPLACE') ;

   Write (Unit = UnFile, FMT = "(2(I10,1X))" ) NParts, NNDH ;

    Do IParts = 1 , NParts
       Write (Unit = UnFile, FMT = "(<NNDH>(I10,1X))" ) ( List_NDAN_rank ( IParts, I ) , I = 1 , NNDH ) ;
    End Do

!/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_



Write(*    ,*) 'End Subroutine < OUTPUT_HDF5 >' ;
Write(UnInf,*) 'End Subroutine < OUTPUT_HDF5 >' ;
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

! =============================================== ERROR IN Write STATEMENT ==========================================================================
1006    Write(*       , Fmt_Write1 ) UnFile, IO_Write ; Write( UnFile, Fmt_Write1 ) UnFile, IO_Write ;
        !#Call BEEP_FAIL ;
        Write(*, Fmt_FL) ;  Write(UnInf, Fmt_FL) ; Write(*, Fmt_End) ; Read(*,*) ;  STOP ;


End Subroutine OUTPUT_HDF5 ;

End Module Partition_Output_HDF5  ;
