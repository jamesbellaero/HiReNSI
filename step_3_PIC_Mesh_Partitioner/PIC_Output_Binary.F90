
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                                                                                                                  ++
! Start:        23 Feb 2013                                                                                                                        ++
! Last Update:  2 April 2013                                                                                                                       ++
!                                                                                                                                                  ++
! Description: THIS Module PRINTS THE OUTPUT OF MESH PARTIIONING for each rank in the BINARY format.                                               ++
!                                                                                                                                                  ++
!              Subroutines                     Description                                                                                         ++
!              ==============                  ==============                                                                                      ++
!              Output                                                                                                                              ++
!                                                                                                                                                  ++
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Module Partition_Output_Binary ;

Use Parameters ;

Implicit None ;

  Interface
!    Module Procedure 
  End Interface    ;

Contains ;

! ***************************************************************************************************************************************************
!                                                                                                                                                  **
! Start:        23 Feb 2013                                                                                                                        **
! Last Update:  2 April 2013                                                                                                                       **
! Description: THIS Subroutine calculates application, PETSc and Local numbering of the elements and writes down all input data for each rank      **
! Called by:                                                                                                                                       **
!                                                                                                                                                  **
! ***************************************************************************************************************************************************

Subroutine OUTPUT_Binary    (                                                                                                   &
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

Integer (Kind=Smll), Intent(In),    Dimension (:  )  :: MTEL, ELT ;

Integer (Kind=Shrt), Intent(In),    Dimension (:  )  :: ELGR, EPart ;
Integer (Kind=Shrt), Intent(In),    Dimension (:  )  :: D_NNZ_Stiff, O_NNZ_Stiff, D_NNZ_Damp, O_NNZ_Damp, D_NNZ_Mass, O_NNZ_Mass ;

Integer (Kind=Lng ), Intent(Inout), Dimension (:  )  :: NDAN , NVAN, NAAN ;
Integer (Kind=Lng ), Intent(In),    Dimension (:  )  :: NEqRank, NNodeRank, NoBndry_DRM, NoLayer_DRM, JLoad ;
Integer (Kind=Lng ), Intent(In),    Dimension (:  )  :: NEL_Rank, NJ_Rank, Global_PETSc_Num ;
Integer (Kind=Lng ), Intent(In),    Dimension (:,:)  :: Local_PETSc_Num, INod, ID, ID_Application, ID_PETSc ;
Integer (Kind=Lng ), Intent(In),    Dimension (:,:)  :: IDBC ;

! - Real Arrays -------------------------------------------------------------------------------------------------------------------------------------
Real (Kind=DBL),     Intent(In),    Dimension (:,:)  :: PMat, PBLD, XYZ, UDis, PLoad, PML_DIM ;

! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------

! - Character Variables -----------------------------------------------------------------------------------------------------------------------------
Character (Kind = 1, Len = 30 ) :: ModelName ;        ! Name of Input file
Character (Kind = 1, Len = 200) :: OutDir ;      ! Directory of output files (Results)

! - Logical Variables -------------------------------------------------------------------------------------------------------------------------------
! - Type DECLERATIONS -------------------------------------------------------------------------------------------------------------------------------
Type ( NodeID )        :: Nodes ;
Type ( BasicParam )    :: Param ;                ! Holds basic parameters of each load case

! =========================== LOCAL Variables =======================================================================================================
! - Integer Variables -------------------------------------------------------------------------------------------------------------------------------
Integer (Kind=Tiny)  :: NNode ;                  ! Number of Nodes

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
!Integer (Kind=Shrt)  :: LDSet ;                  ! Load Set.
Integer (Kind=Shrt)  :: ILType ;                 ! Load Type of Element.
!Integer (Kind=Shrt)  :: Inc ;                    ! Increment of elements for load sets.
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

Integer (Kind=Lng )  :: NLN ;                    ! NUmber of Loaded Nodes in the original model
Integer (Kind=Lng )  :: NSND ;                   ! Number of Support Nodes with Displacements in the original model
Integer (Kind=Lng )  :: IEL ;                    ! Loop index on NEL.
Integer (Kind=Lng )  :: IJ ;                     ! Loop index on NJ.
Integer (Kind=Lng )  :: I, J, K, I1, I2 ;        ! Loop indeces.
Integer (Kind=Lng )  :: Lo_EL_Num ;              ! Load set of Element Number.
Integer (Kind=Lng )  :: NEqM_Mapping ;           ! Number of Modified EQuations for Mapping application numbering to PETSc numbering. !!?? at the end of the day see if we need this variable.
Integer (Kind=Lng )  :: LEqN ;                   ! Local Equation number, for ID_Local.
Integer (Kind=Lng )  :: NEQM ;                   ! Number of Modified Equations
Integer (Kind=Lng )  :: Counter ;                ! Counter
Integer (Kind=Lng )  :: LowerLimit ;             ! Lower Limit of node number based on the PETSc node numbering system, stored on a rank.
Integer (Kind=Lng )  :: UpperLimit ;             ! Upper Limit of node number based on the PETSc node numbering system, stored on a rank.
Integer (Kind=Lng )  :: PETScNodeNumber ;        ! A temporary varibale for saving node number based on the PETSc node numbering.
Integer (Kind=Lng )  :: Node ;                   ! Holds Node number

! - Real Variables ----------------------------------------------------------------------------------------------------------------------------------
! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
! - Integer Arrays ----------------------------------------------------------------------------------------------------------------------------------
Integer (Kind=Shrt)  :: Local_ElN ( NEL ) ;                          ! Local Element Number

Integer (Kind=Lng )  :: LTEL_Rank ( 4, Param%IntM( 1, 2) ) ;         ! Load Type Element - Param%IntM( 1, 2) = NPBL

Integer (Kind=Lng ), Allocatable, Dimension(:  )  :: App_Numbers ;   ! 
Integer (Kind=Lng ), Allocatable, Dimension(:  )  :: PETSc_Numbers ; ! 
Integer (Kind=Lng ), Allocatable, Dimension(:  )  :: Indices ;       ! 
Integer (Kind=Lng ), Allocatable, Dimension(:,:)  :: ID_Local ;      ! 

Integer (Kind=Lng ), Allocatable, Dimension(:  )  :: EqDis ;         ! Holds Equation numbers of all nodes' DOF in which the history of EqDisplacement is required based on PETSc node numbering 
Integer (Kind=Lng ), Allocatable, Dimension(:  )  :: EqVel ;         ! Holds Equation numbers of all nodes' DOF in which the history of velocity is required based on PETSc node numbering 
Integer (Kind=Lng ), Allocatable, Dimension(:  )  :: EqAcc ;         ! Holds Equation numbers of all nodes' DOF in which the history of acceleration is required based on PETSc node numbering 
Integer (Kind=Lng ), Allocatable, Dimension(:  )  :: NodeDis ;       ! Holds Node numbers in which the history of EqDisplacement is required based on Global node numbering 
Integer (Kind=Lng ), Allocatable, Dimension(:  )  :: NodeVel ;       ! Holds Node numbers in which the history of velocity is required based on Global node numbering 
Integer (Kind=Lng ), Allocatable, Dimension(:  )  :: NodeAcc ;       ! Holds Node numbers in which the history of acceleration is required based on Global node numbering 

! - Real Arrays -------------------------------------------------------------------------------------------------------------------------------------
! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
! - Character Variables -----------------------------------------------------------------------------------------------------------------------------
Character (Kind = 1, Len = 20)  :: IndexRank ;   ! A variable for storing the Rank number in Character format for adding at the end of the input file Name.
Character (Kind = 1, Len = 20)  :: IndexSize ;   ! A variable for storing the Size number in Character format for adding at the end of the input file Name.

! - Logical Variables -------------------------------------------------------------------------------------------------------------------------------
! - Type DECLERATIONS -------------------------------------------------------------------------------------------------------------------------------
Type ( BasicParam )    :: ParamPartition ;       ! Holds basic parameters of each load case in each partition

! =========================== Subroutine CODE =======================================================================================================

Allocate ( ParamPartition%IntM ( 7, 6), ParamPartition%RealM ( 7, 6) ) ;

  If ( LoadC (5) /= 0_Tiny ) Then ;  ! Dynamic Analysis

    NNDH = Param%IntM( 6, 1) ;
    NNVH = Param%IntM( 6, 2) ;
    NNAH = Param%IntM( 6, 3) ;

    ! Allocating Required Arrays Dynamic Analysis
    Allocate ( EqDis ( NNDH * NDim ), EqVel ( NNVH * NDim ), EqAcc ( NNAH * NDim ), NodeDis ( NNDH ), NodeVel ( NNVH ), NodeAcc ( NNAH ),     STAT = ERR_Alloc) ; 
      If ( ERR_Alloc /= 0 ) Then ;
        Write (*, Fmt_ALLCT) ERR_Alloc ;  Write (UnInf, Fmt_ALLCT) ERR_Alloc ;
        Write(*, Fmt_FL) ;  Write(UnInf, Fmt_FL) ; Read(*, Fmt_End) ;  Stop ;
      End If ;

  End If ;

! - Writing down the input data for each process ----------------------------------------------------------------------------------------------------
  Do IParts = 1, NParts ; 

    Write (IndexRank, *) IParts - 1_Shrt ; ! Converts Rank number to Character foramt for the file Name
    Write (IndexSize, *) NParts ;          ! Converts Size number to Character foramt for the file Name

    Write (*,*)"Writing files for partition: ", IParts;

    ! - Output FILEs --------------------------------------------------------------------------------------------------------------------------------
    Write (*,*)"Opening files ...";

    ! Data file
    UnFile = UN_Out ;
    Open ( Unit = UnFile, FILE = TRIM(ModelName)//'_'//Trim(AdjustL(IndexSize))//'_'//Trim(AdjustL(IndexRank))//'.dataModel', ERR =  1001, IOSTAT = IO_File, ACCESS = 'SEQUENTIAL', ACTION ='Write', ASYNCHRONOUS = 'NO', BLANK = 'NULL', BLOCKSIZE = 0, DEFAULTFILE = TRIM(OutDir), DisPOSE = 'KEEP', FORM = 'Binary', POSITION = 'ASIS', STATUS ='REPLACE') ;

    ! Coordinate file (.XYZ)
    UnFile = Un_OutXYZ ;
    Open ( Unit = UnFile, FILE = TRIM(ModelName)//'_'//Trim(AdjustL(IndexSize))//'_'//Trim(AdjustL(IndexRank))//'.XYZ', ERR =  1001, IOSTAT = IO_File, ACCESS = 'SEQUENTIAL', ACTION ='Write', ASYNCHRONOUS = 'NO', BLANK = 'NULL', BLOCKSIZE = 0, DEFAULTFILE = TRIM(OutDir), DisPOSE = 'KEEP', FORM = 'Binary', POSITION = 'ASIS', STATUS ='REPLACE') ;

    ! Connectivity file (.Cnn)
    UnFile = Un_OutCnn ;
    Open ( Unit = UnFile, FILE = TRIM(ModelName)//'_'//Trim(AdjustL(IndexSize))//'_'//Trim(AdjustL(IndexRank))//'.Cnn', ERR =  1001, IOSTAT = IO_File, ACCESS = 'SEQUENTIAL', ACTION ='Write', ASYNCHRONOUS = 'NO', BLANK = 'NULL', BLOCKSIZE = 0, DEFAULTFILE = TRIM(OutDir), DisPOSE = 'KEEP', FORM = 'Binary', POSITION = 'ASIS', STATUS ='REPLACE') ;

    ! Constraint file (.Cnt)
    UnFile = Un_OutCnt ;
    Open ( Unit = UnFile, FILE = TRIM(ModelName)//'_'//Trim(AdjustL(IndexSize))//'_'//Trim(AdjustL(IndexRank))//'.Cnt', ERR =  1001, IOSTAT = IO_File, ACCESS = 'SEQUENTIAL', ACTION ='Write', ASYNCHRONOUS = 'NO', BLANK = 'NULL', BLOCKSIZE = 0, DEFAULTFILE = TRIM(OutDir), DisPOSE = 'KEEP', FORM = 'Binary', POSITION = 'ASIS', STATUS ='REPLACE') ;

    ! Non-Zeros of Stiffness file (.NNZStiff)
    UnFile = Un_OutStiff ;
    Open ( Unit = UnFile, FILE = TRIM(ModelName)//'_'//Trim(AdjustL(IndexSize))//'_'//Trim(AdjustL(IndexRank))//'.NNZStiff', ERR =  1001, IOSTAT = IO_File, ACCESS = 'SEQUENTIAL', ACTION ='Write', ASYNCHRONOUS = 'NO', BLANK = 'NULL', BLOCKSIZE = 0, DEFAULTFILE = TRIM(OutDir), DisPOSE = 'KEEP', FORM = 'Binary', POSITION = 'ASIS', STATUS ='REPLACE') ;

    ! Non-Zeros of Damping file (.NNZDamp)
    UnFile = Un_OutDamp ;
    Open ( Unit = UnFile, FILE = TRIM(ModelName)//'_'//Trim(AdjustL(IndexSize))//'_'//Trim(AdjustL(IndexRank))//'.NNZDamp', ERR =  1001, IOSTAT = IO_File, ACCESS = 'SEQUENTIAL', ACTION ='Write', ASYNCHRONOUS = 'NO', BLANK = 'NULL', BLOCKSIZE = 0, DEFAULTFILE = TRIM(OutDir), DisPOSE = 'KEEP', FORM = 'Binary', POSITION = 'ASIS', STATUS ='REPLACE') ;

    ! Non-Zeros of Damping file (.NNZMass)
    UnFile = Un_OutMass ;
    Open ( Unit = UnFile, FILE = TRIM(ModelName)//'_'//Trim(AdjustL(IndexSize))//'_'//Trim(AdjustL(IndexRank))//'.NNZMass', ERR =  1001, IOSTAT = IO_File, ACCESS = 'SEQUENTIAL', ACTION ='Write', ASYNCHRONOUS = 'NO', BLANK = 'NULL', BLOCKSIZE = 0, DEFAULTFILE = TRIM(OutDir), DisPOSE = 'KEEP', FORM = 'Binary', POSITION = 'ASIS', STATUS ='REPLACE') ;

    ! Application Ordering file (.App)
    UnFile = Un_OutApp ;
    Open ( Unit = UnFile, FILE = TRIM(ModelName)//'_'//Trim(AdjustL(IndexSize))//'_'//Trim(AdjustL(IndexRank))//'.App', ERR =  1001, IOSTAT = IO_File, ACCESS = 'SEQUENTIAL', ACTION ='Write', ASYNCHRONOUS = 'NO', BLANK = 'NULL', BLOCKSIZE = 0, DEFAULTFILE = TRIM(OutDir), DisPOSE = 'KEEP', FORM = 'Binary', POSITION = 'ASIS', STATUS ='REPLACE') ;

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

    ! Modified Number of Equations
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

    ! pressure on elements
    NIDBC_Rank = 0_Shrt ;
      If ( LoadC (2) /= 0_Tiny .or. LoadC (5) /= 0_Tiny ) Then ;
        Do I = 1, Param%IntM( 2, 4) ;
          If ( EPart ( IDBC ( I, 1 ) ) == IParts ) Then ;
            NIDBC_Rank = NIDBC_Rank + 1 ;
          End If ;
        End Do ;
      End If ;

    ! Equation numbers of nodes in which history is required for dynamic analysis
    NNDH_Rank = 0_Lng ;
    NNVH_Rank = 0_Lng ;
    NNAH_Rank = 0_Lng ;

      ! Obtaining Equation number of nodes in which history of EqDisplacement is required based on PETSc node numbering
      If ( LoadC (5) /= 0_Tiny ) Then ;
        NNDH = Param%IntM( 6, 1) ;
          Do I = 1, NNDH ;
            If ( NDAN ( I ) == 0_Shrt  ) Cycle ; ! why? because we do not want to apply the loads on the boundary of the partitions twice. see inside the loop
              If ( Local_PETSc_Num ( NDAN ( I ), IParts ) /= 0_Lng ) Then ;
                NNDH_Rank = NNDH_Rank + 1_Lng ;
                EqDis ( ( NNDH_Rank - 1_Lng )* NDim+1:( NNDH_Rank - 0_Lng )* NDim ) = ID_PETSc ( Global_Petsc_Num ( NDAN ( I ) ), 1:NDim )  ;
                NodeDis ( NNDH_Rank ) = NDAN ( I ) ;
                NDAN ( I ) = 0_Shrt ; ! In this way, the load on this node will be applied once in one rank, if the node is located at the border of two or more ranks 
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
        ParamPartition%IntM( 2, 4) = NIDBC_Rank ;  ! NIDBC
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
        ParamPartition%IntL( 5, 1:6) = Param%IntL( 5, 1:6) ; ParamPartition%IntM( 5, 3) = Param%IntM( 5, 3) ; ParamPartition%RealL( 5, 1:4 ) = Param%RealL( 5, 1:4 ) ; ! NSTEP, NTST, LoadType   Delta, Gama, DT, t0,            Newmark Coefficients, time

        ParamPartition%IntM( 6, 1) = NNDH_Rank ; ! NNDH
        ParamPartition%IntM( 6, 2) = NNVH_Rank ; ! NNVH
        ParamPartition%IntM( 6, 3) = NNAH_Rank ; ! NNAH
        ParamPartition%IntM( 6, 4) = Param%IntM( 6, 4) ; ! NDamp
        ParamPartition%IntM( 6, 5) = Param%IntM( 6, 5) ; ! NEnergy
        ParamPartition%IntM( 6, 6) = Param%IntM( 6, 6) ; ! PARAM_Type ! pml analysis


          If      ( Param%IntM( 5, 3) == 1_Tiny ) Then ;  ! Dynamic Pressure
            ParamPartition%IntM( 2, 4) = NIDBC_Rank ;   ! NIDBC
          !Else If ( Param%IntM( 5, 3) == 2_Tiny ) Then ; ! Base Acceleration
          !  ParamPartition%IntL( 7, 1) = Param%IntL( 7, 1) ;  ParamPartition%RealL( 7, 1) = Param%RealL( 7, 1) ; ! NASTEP, G   Number of base Accelerations to be read in the input file
          Else If ( Param%IntM( 5, 3) == 3_Tiny ) Then ; ! Domain Reduction Method
            ParamPartition%IntM  ( 7, 1 ) = NNBndry_DRM_Rank ; ! NNBndry_DRM
            ParamPartition%IntM  ( 7, 2 ) = NNLayer_DRM_Rank ; ! NNLayer_DRM
          End If ;

      End If ;

! - Basic Data -------------------------------------------------------------------------------------------------------------------------------------
    ! Basic data
    UnFile = Un_Out ;
    Write (*,*)"Writing Basic Data ..." ;
    !Write (Unit = UnFile, FMT = "(6(I3,2X),1(I5,2X),6(I19,2X))", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) NDOF, NDim, MaxNNode, NInt, NInt_Type, SOLVER_Type,      NGroup,         NEL_Rank ( IParts ), NJ_Rank ( IParts ), NJ, NEQM, NEQMTotal, NEqM_Mapping ;
    Write (Unit = UnFile, ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) NDOF, NDim, MaxNNode,        NGroup,         NEL_Rank ( IParts ), NJ_Rank ( IParts ), NJ, NEQM, NEQMTotal, NEqM_Mapping ;

    ! Available Load Cases
    Write(*,*)"writing Load Cases ..." ;
    Write (Unit = UnFile, ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) ( LoadC ( I ), I = 1, 5 ) ;

    ! Writing Param
    Write(*,*)"Writing Param Array ..." ;
    Write (Unit = UnFile, ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) (( ParamPartition%IntM  ( I, J), I = 1, 7), J = 1, 6) ;
    Write (Unit = UnFile, ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) (( ParamPartition%RealM  ( I, J), I = 1, 7), J = 1, 6) ;

! - Geometry & Model -------------------------------------------------------------------------------------------------------------------------------
    Write (*,*)"Writing Geometry ..." ;

!    ! Number of elements of each group
!    Write(*,*)"Writing Group numbers ..." ;
!      DO I = 1, NGroup ;
!        Write (Unit = UnFile, ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) NEG ( I, 1 ), NEG ( I, 2 ) ; ! NUMBER OF ELEMENTS OF THIS GROUP , ELEMENT Type
!      End Do ;

    If ( LoadC (5) /= 0_Tiny ) Then ;
      Write (*,*)"Writing PML teritory ..." ;
      ! PML teritory
        DO J = 1, 2 ! See Related PML Subroutines
          Write (Unit = UnFile, ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) (PML_DIM ( I, J ), I = 1, 2 * NDim ) ; 
        End Do ;
    End If ;

    ! Material properties of each material
    Write (*,*)"Writing Material Properties ..." ;
    UnFile = Un_OutMat ;
      DO I = 1, NMat ;
        Write (Unit = UnFile, ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 )( PMat ( I, J ), J = 1, NPM ) ; 
      End Do ;

    ! Coordinates of nodes
    Write (*,*)"Writing Coordinates ..." ;
    UnFile = Un_OutXYZ ;
      DO K = 1, NJ ;
        If ( Local_PETSc_Num ( K, IParts ) /= 0_Lng ) Then ; !?? check the oreintation 
          Write (Unit = UnFile, ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 )( XYZ ( K, J ), J = 1, NDim ) ; 
        End If ;
      End Do ;

    ! Element connectivities
    Write (*,*)"Writing Node connectivities ..." ;
    UnFile = Un_OutCnn ;
      DO IEL = 1, NEL ;

        If ( EPart ( IEL ) == IParts ) Then ;
          EType = ElT ( IEl ) ; ! Element Type

            ! determinig NNode for this element
            IF      ( EType == El2d4NSldPN  ) Then ;  NNode = 4 ;     ! 4 node  -2D - Solid - PLANE STRESS
            Else If ( EType == El2d4NPMLPN  ) Then ;  NNode = 4 ;     ! 4 node  -2D - PML - PLANE STRESS
            Else If ( EType == El2d8NSldPN  ) Then ;  NNode = 8 ;     ! 8 node  -2D - Solid - PLANE STRAIN
            Else If ( EType == El2d8NPMLPN  ) Then ;  NNode = 8 ;     ! 8 node  -2D - PML - PLANE STRAIN
            Else If ( EType == El2d6NSldPN  ) Then ;  NNode = 6 ;     ! 6 node  -2D - Solid - PLANE STRAIN
            Else If ( EType == El2d6NPMLPN  ) Then ;  NNode = 6 ;     ! 6 node  -2D - PML - PLANE STRAIN
            Else If ( EType == El2d3NSldPN  ) Then ;  NNode = 3 ;     ! 3 node  -2D - Solid - PLANE STRAIN
            Else If ( EType == El2d3NPMLPN  ) Then ;  NNode = 3 ;     ! 3 node  -2D - PML - PLANE STRAIN
            Else If ( EType == SEl2d9NSldPN ) Then ;  NNode = 9 ;     ! 9 node  -2D - Solid - PLANE STRAIN
            Else If ( EType == SEl2d9NPMLPN ) Then ;  NNode = 9 ;     ! 9 node  -2D - NPML - PLANE STRAIN
            Else If ( EType == SEl2d9NMPMLPN) Then ;  NNode = 9 ;     ! 9 node  -2D - MPML - PLANE STRAIN
            Else If ( EType == SEl3d27NSld  ) Then ;  NNode = 27;     ! 27 node -3D - MPML
            Else If ( EType == SEl3d27NMPML ) Then ;  NNode = 27;     ! 27 node -3D - MPML
            Else If ( EType == SEl2d7NSldPN ) Then ;  NNode = 7 ;     ! 7 node -2D - solid - plane strain
            Else If ( EType == SEl2d7NMPMLPN) Then ;  NNode = 7 ;     ! 7 node -2D - MPML - plane strain
            Else If ( EType == SEl2d9NSldSH ) Then ;  NNode = 9 ;     ! 9 node - 2D - solid ( quadrilateral ) - Spectral Element - SH waves
            Else If ( EType == SEl2d9NPMLSH ) Then ;  NNode = 9 ;     ! 9 node - 2D - PML   ( quadrilateral ) - Spectral Element - SH waves
            Else If ( EType == SEl3d8NSld   ) Then ;  NNode = 8 ;     ! 8 node -3D - Solid
            Else If ( EType == SEl3d8NMPML  ) Then ;  NNode = 8 ;     ! 8 node -3D - MPML
            End If ;

          Write (Unit = UnFile, ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) ( Local_PETSc_Num ( INOD ( I, IEL ), IParts ), I = 1, NNode ), (I*0, I = 1, MaxNNode - NNode ), MTEL ( IEL ), ELT ( IEL ), ELGR ( IEL ) ;
        End If ;

      End Do ;

    ! Constraints
    Write (*,*)"Writing Node constraints ..." ;
    UnFile = Un_OutCnt ;
    ! Writing the Local Equation Numbers (Local ID)
      Do IJ = 1, NJ_Rank ( IParts ) ;
        Write (Unit = UnFile, ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) ( ID_Local ( IJ, INode ), INode = 1, NDOF ) ;
      End Do ;

! - Loads -------------------------------------------------------------------------------------------------------------------------------------------
    Write (*,*)"Writing Load types ..." ;
    UnFile = Un_Out ;

    ! Body Force Load
      If ( LoadC (1) /= 0_Tiny ) Then ;
        Write(*,*)"Reading Body Forces ..." ;
          ! Body force load
          DO I = 1, ParamPartition%IntM( 1, 1) ; ! ParamPartition%IntP( 1, 1)=NBLD: Number of body force load.   ParamPartition%IntP( 1, 2) = NPBL
            Write (Unit = UnFile, ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) ( PBLD ( I, J ) , J = 1, ParamPartition%IntM( 1, 2) ) ;
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
!        Write (Unit = UnFile, ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) LDSET ;
!          DO I = 1, LDSET ;
!            Write (Unit = UnFile, ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) ( LTEL_Rank ( K, I ), K = 1, 4 ) ;
!          End Do ;

      End If ;

    ! Static Pressure
      If ( LoadC (2) /= 0_Tiny ) Then ;
        Write(*,*)"Writing Pressure Loads ..." ;
          Do I = 1, Param%IntM( 2, 4) ;  ! NIDBC
            If ( EPart ( IDBC ( I, 1 ) ) == IParts ) Then ;
              Write (Unit = UnFile, ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) Local_ElN ( IDBC ( I, 1 ) ),  ( IDBC ( I, J ), J = 2, 2 * NDim +1) ;
            End If ;
          End Do ;
      End If ;

    ! Joint load of every case - Concentrated loads
      If ( LoadC ( 3 ) /= 0 ) Then ;
        Write(*,*)"Writing Joint Loads ..." ;
          DO I = 1, NLN_Rank ;  ! ParamPartition%IntP( 3, 1) = NLN_Rank 
            If ( Local_PETSc_Num ( JLoad ( I ) , IParts ) /= 0_Lng ) Then ; 
              If ( Nodes%Rep ( JLoad ( I ) ) == 1_Shrt ) Then ;
                Write (Unit = UnFile, ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 )  Local_PETSc_Num ( JLoad ( I ), IParts ), ( PLoad( J, I ), J = 1, NDim ) ;
              Else If ( Nodes%Rep ( JLoad ( I ) ) /= 1_Shrt .AND. IParts == MinVal ( Nodes%Locs( JLoad ( I ),1:Nodes%Rep ( JLoad ( I ) ) ) ) ) Then ; ! Takes care of nodes on the boundary of the partitions. So, we do not have a repeatation of the loads
                Write (Unit = UnFile, ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 )  Local_PETSc_Num ( JLoad ( I ), IParts ), ( PLoad( J, I ), J = 1, NDim ) ;
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
                Write (Unit = UnFile, ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 )  Local_PETSc_Num ( UDis( NDOF + 1, I ), IParts ), ( UDis ( K, I ), K = 1, NDOF ) ;
              Else If ( Nodes%Rep ( UDis( NDOF + 1, I ) ) /= 1_Shrt .AND. IParts == MinVal ( Nodes%Locs( UDis( NDOF + 1, I ),1:Nodes%Rep ( UDis( NDOF + 1, I ) ) ) ) ) Then ;
                Write (Unit = UnFile, ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 )  Local_PETSc_Num ( UDis( NDOF + 1, I ), IParts ), ( UDis ( K, I ), K = 1, NDOF ) ;
              End If ;
            End If ; 
          End Do ;
      End If ;

    Write(*,*)"Writing Dynamic Loads ..." ;

    ! Dynamic loads
      If ( LoadC (5) /= 0_Tiny ) Then ;

        If ( ParamPartition%IntM( 5, 3 ) == 1 ) Then ;    ! RICKER PULSE or sine function for dynamic pressure

          Write(*,*)"Writing Function Type ..." ;
            Do I = 1, Param%IntM( 2, 4) ;
              If ( EPart ( IDBC ( I, 1 ) ) == IParts ) Then ;
                Write (Unit = UnFile, ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) Local_ElN ( IDBC ( I, 1 ) ),  ( IDBC ( I, J ), J = 2, 2 * NDim +1) ;
              End If ;
            End Do ;

        Else If ( ParamPartition%IntM( 5, 3 ) == 2_Tiny ) Then ;   ! Base acceleration

!          Write(*,*)"Writing Base Acceleration ..." ;
!            DO I = 1, ParamPartition%IntL( 7, 1) ; ! NAStep = ParamPartition%IntP( 7, 1)
!              Write (Unit = UnFile, ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) ( BACL ( I, K ), K = 1, NDim ) ;
!            End Do ;

        Else If ( ParamPartition%IntM( 5, 3 ) == 3_Tiny ) Then ;    ! Domain Reduction Method (DRM)

!          Write(*,*)"Writing wave information ..." ;
!          Write (Unit = UnFile, ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) InciWave(1), InciWave(2), InciWave(3), InciWave(4), InciWave(5) ;       ! Theta, Omega, amplitude, alpha1, alpha2

          ! writing the node numbers on the DRM boundary
          Write(*,*)"Writing DRM boundary nodes ..." ;
            Do I = 1, Param%IntM( 7, 1) ; ! NNBndry_DRM
              Node = NoBndry_DRM ( I ) ;
              If ( ( Nodes%Rep (Node) == 1_Shrt .AND. Any ( Nodes%Locs ( Node, : ) == IParts ) ) .OR. ( Nodes%Rep (Node) /= 1_Shrt .AND. IParts == MinVal ( Nodes%Locs( Node, :), MASK = Nodes%Locs( Node, :) > 0_Shrt ) ) ) Write (Unit = UnFile, ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) Local_PETSc_Num ( Node, IParts ) ;
            End Do ;

          ! Reading the node numbers on the DRM neighbor
          Write(*,*)"Writing DRM neighbor nodes ..." ;
            Do I = 1, Param%IntM( 7, 2) ! NNLayer_DRM
              Node = NoLayer_DRM ( I ) ;
              If ( ( Nodes%Rep (Node) == 1_Shrt .AND. Any ( Nodes%Locs ( Node, : ) == IParts ) ) .OR. ( Nodes%Rep (Node) /= 1_Shrt .AND. IParts == MinVal ( Nodes%Locs( Node, :), MASK = Nodes%Locs( Node, :) > 0_Shrt ) ) ) Write (Unit = UnFile, ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) Local_PETSc_Num ( Node, IParts ) ;
            End Do ;

        End If ;

      End If ;

!    ! writes down step numbers for full results
!    Write(*,*)"Writing Step Numbers ..." ;
!    Write (Unit = UnFile, ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) ( STEP ( I ), I = 1, ParamPartition%IntL( 5, 2) ) ; ! NTST

      ! writes down node numbers and equation numbers of nodes in which history of Displacement is required
      If ( NNDH_Rank /= 0_Shrt ) Then ;
        Write(*,*)"Writing Node Numbers for displacements ..." ;
        Write (Unit = UnFile, ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) ( NodeDis ( I ), I = 1, NNDH_Rank ), ( Global_PETSc_Num ( NodeDis ( I ) ), I = 1, NNDH_Rank ) ;
        Write (Unit = UnFile, ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) ( EqDis ( I ), I = 1, NNDH_Rank * NDim ) ;
      End If ;

      ! writes down node numbers and equation numbers of nodes in which history of velocity is required
      If ( NNVH_Rank /= 0_Shrt ) Then ;
        Write(*,*)"Writing Node Numbers for velocity ..." ;
        Write (Unit = UnFile, ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) ( NodeVel ( I ), I = 1, NNVH_Rank ), ( Global_PETSc_Num ( NodeVel ( I ) ), I = 1, NNVH_Rank ) ;
        Write (Unit = UnFile, ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) ( EqVel ( I ), I = 1, NNVH_Rank * NDim ) ;
      End If ;

      ! writes down node numbers and equation numbers of nodes in which history of acceleration is required
      If ( NNAH_Rank /= 0_Shrt ) Then ;
        Write(*,*)"Writing Node Numbers for acceleration ..." ;
        Write (Unit = UnFile, ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) ( NodeAcc ( I ), I = 1, NNAH_Rank ), ( Global_PETSc_Num ( NodeAcc ( I ) ), I = 1, NNAH_Rank ) ;
        Write (Unit = UnFile, ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) ( EqAcc ( I ), I = 1, NNAH_Rank * NDim ) ;
      End If ;

! <><><><><><><><><><><><><><><><><><><><><>
    ! Application Numbering
    Write (*,*)"Writing application numberting ..." ;

    Allocate ( App_Numbers ( NEqM_Mapping ), PETSc_Numbers ( NEqM_Mapping ), Indices ( NEQM ) ) ;
    App_Numbers   = 0_Lng ;
    PETSc_Numbers = 0_Lng ;

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
                  If ( Counter > NEqM_Mapping ) Then ; !!!!!!??
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
    
    App_Numbers = App_Numbers - 1_Lng ;
    PETSc_Numbers = PETSc_Numbers - 1_Lng ;
    Indices = Indices - 1_Lng ;
    
    ! Writing down the data for node numbering mapping
    UnFile = Un_OutApp ;
    Write(*,*)"Writing Index Sets ..." ;
    !Write (Unit = UnFile, ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) ( App_Numbers   ( I ), I = 1, NEqM_Mapping ) ;
    !Write (Unit = UnFile, ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) ( PETSc_Numbers ( I ), I = 1, NEqM_Mapping ) ;
    !Write (Unit = UnFile, ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) ( Indices       ( I ), I = 1, NEQM ) ;

    Write (Unit = UnFile, ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) int4(App_Numbers) ;
    Write (Unit = UnFile, ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) int4(PETSc_Numbers) ;
    Write (Unit = UnFile, ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) int4(Indices) ;

    ! Number of non-zero entris of PETSc objects
    Write (*,*)"Writing number of nonzero elements ..." ;

      If ( IParts == 1_Shrt ) then ;
        LowerLimit = 1_Lng ;
      Else ;
        LowerLimit = NEqRank ( IParts - 1_Shrt ) + 1_Lng ;
      End If ;

    UpperLimit = NEqRank ( IParts ) ;

    UnFile = Un_OutStiff ;
    Write(*,*)"Writing Number of Non-Zeros of Stiffness ..." ;
    Write (Unit = UnFile, ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) ( D_NNZ_Stiff ( I ), I = LowerLimit, UpperLimit ) ;
    Write (Unit = UnFile, ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) ( O_NNZ_Stiff ( I ), I = LowerLimit, UpperLimit ) ;

    UnFile = Un_OutDamp ;
    Write(*,*)"Writing Number of Non-Zeros of Damping ..." ;
    Write (Unit = UnFile, ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) ( D_NNZ_Damp  ( I ), I = LowerLimit, UpperLimit ) ;
    Write (Unit = UnFile, ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) ( O_NNZ_Damp  ( I ), I = LowerLimit, UpperLimit ) ;

    UnFile = Un_OutMass ;
    Write(*,*)"Writing Number of Non-Zeros of Mass ..." ;
    Write (Unit = UnFile, ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) ( D_NNZ_Mass  ( I ), I = LowerLimit, UpperLimit ) ;
    Write (Unit = UnFile, ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) ( O_NNZ_Mass  ( I ), I = LowerLimit, UpperLimit ) ;

    DeAllocate ( ID_Local, App_Numbers, PETSc_Numbers, Indices ) ;

    ! - Closing the output file ---------------------------------------------------------------------------------------------------------------------
    UnFile =  UN_Out ;
    Close ( Unit = UnFile, Status = 'KEEP', ERR =  1002, IOSTAT = IO_File ) ;

    UnFile =  UN_OutXYZ ;
    Close ( Unit = UnFile, Status = 'KEEP', ERR =  1002, IOSTAT = IO_File ) ;

    UnFile =  UN_OutCnn ;
    Close ( Unit = UnFile, Status = 'KEEP', ERR =  1002, IOSTAT = IO_File ) ;

    UnFile =  UN_OutCnt ;
    Close ( Unit = UnFile, Status = 'KEEP', ERR =  1002, IOSTAT = IO_File ) ;

    UnFile =  Un_OutStiff ;
    Close ( Unit = UnFile, Status = 'KEEP', ERR =  1002, IOSTAT = IO_File ) ;

    UnFile =  Un_OutDamp ;
    Close ( Unit = UnFile, Status = 'KEEP', ERR =  1002, IOSTAT = IO_File ) ;

    UnFile =  Un_OutMass ;
    Close ( Unit = UnFile, Status = 'KEEP', ERR =  1002, IOSTAT = IO_File ) ;

    UnFile =  Un_OutApp ;
    Close ( Unit = UnFile, Status = 'KEEP', ERR =  1002, IOSTAT = IO_File ) ;

  End Do ;


  ! This part is for debugging purposes. It generates input files for serial code based on the PETSc node numbering.
!  If (.False.) Then ;
!
!write(*,*)'debug 1001'
!
!    ! Generating coordinates files
!    UnFile = Un_OutXYZ ;
!    Open ( Unit = UnFile, FILE = TRIM(ModelName)//'_'//Trim(AdjustL(IndexSize))//'_Serial.XYZ', ERR =  1001, IOSTAT = IO_File, ACCESS = 'SEQUENTIAL', ACTION ='Write', ASYNCHRONOUS = 'NO', BLANK = 'NULL', BLOCKSIZE = 0, DEFAULTFILE = TRIM(OutDir), DisPOSE = 'KEEP', FORM = 'FORMATTED', POSITION = 'ASIS', STATUS ='REPLACE') ;
!      Do IJ = 1, NJ ;
!        Write (Unit = UnFile, FMT = "(I19,2x,<NDim>(E33.25E3,2X) )", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) Global_PETSc_Num ( IJ ), ( XYZ ( IJ, J ), J = 1, NDim ) ;
!      End Do ;
!    Close ( Unit = UnFile, Status = 'KEEP', ERR =  1002, IOSTAT = IO_File ) ;
!
!write(*,*)'debug 1002'
!
!    ! Connectivity file (.Cnn)
!    UnFile = Un_OutCnn ;
!    Open ( Unit = UnFile, FILE = TRIM(ModelName)//'_'//Trim(AdjustL(IndexSize))//'_Serial.Cnn', ERR =  1001, IOSTAT = IO_File, ACCESS = 'SEQUENTIAL', ACTION ='Write', ASYNCHRONOUS = 'NO', BLANK = 'NULL', BLOCKSIZE = 0, DEFAULTFILE = TRIM(OutDir), DisPOSE = 'KEEP', FORM = 'FORMATTED', POSITION = 'ASIS', STATUS ='REPLACE') ;
!
!      DO IEL = 1, NEL ;
!
!        If ( EPart ( IEL ) == IParts ) Then ;
!          EType = ElT ( IEl ) ; ! Element Type
!            ! determinig NNode for this element
!            IF      ( EType == 1  ) Then ;  NNode = 4 ;     ! 4 noded -2D - Solid - PLANE STRESS
!            Else If ( EType == 4  ) Then ;  NNode = 8 ;     ! 8 noded -2D - Solid - PLANE STRAIN
!            Else If ( EType == 104) Then ;  NNode = 8 ;     ! 8 noded -2D - PML - PLANE STRAIN
!            End If ;
!
!          Write (Unit = UnFile, FMT = "(<MaxNNode>(I19,2X),I5,2X,I5,2X,I5,2X,<2*NDim>(I5,2X))", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) ( Local_PETSc_Num ( INOD ( I, IEL ), IParts ), I = 1, NNode ), (0_Lng, I = 1, MaxNNode - NNode ), MTEL ( IEL ), ELT ( IEL ), ELGR ( IEL ) ;
!
!        End If ;
!      End Do ;
!    Close ( Unit = UnFile, Status = 'KEEP', ERR =  1002, IOSTAT = IO_File ) ;
!
!write(*,*)'debug 1003'
!
!    ! Constraint file (.Cnt)
!    UnFile = Un_OutCnt ;
!    Open ( Unit = UnFile, FILE = TRIM(ModelName)//'_'//Trim(AdjustL(IndexSize))//'_Serial.Cnt', ERR =  1001, IOSTAT = IO_File, ACCESS = 'SEQUENTIAL', ACTION ='Write', ASYNCHRONOUS = 'NO', BLANK = 'NULL', BLOCKSIZE = 0, DEFAULTFILE = TRIM(OutDir), DisPOSE = 'KEEP', FORM = 'FORMATTED', POSITION = 'ASIS', STATUS ='REPLACE') ;
!      Do IJ = 1, NJ ;
!        Write (Unit = UnFile, FMT = "(I19,2x,<NDOF>(I19,2X) )", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) Global_PETSc_Num ( IJ ), ( ID ( IJ, I ), I = 1, NDOF ) ;
!      End Do ;
!    Close ( Unit = UnFile, Status = 'KEEP', ERR =  1002, IOSTAT = IO_File ) ;
!
!write(*,*)'debug 1004'
!
!    UnFile = UN_Out ;
!    Open ( Unit = UnFile, FILE = TRIM(ModelName)//'_'//Trim(AdjustL(IndexSize))//'_Serial.data', ERR =  1001, IOSTAT = IO_File, ACCESS = 'SEQUENTIAL', ACTION ='Write', ASYNCHRONOUS = 'NO', BLANK = 'NULL', BLOCKSIZE = 0, DEFAULTFILE = TRIM(OutDir), DisPOSE = 'KEEP', FORM = 'FORMATTED', POSITION = 'ASIS', STATUS ='REPLACE') ;
!    ! node numbers on the DRM boundary
!    Write (Unit = UnFile, FMT = "('Boundary nodes')", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) ;
!!      Do I = 1, NNBndry_DRM ;
!!        Write (Unit = UnFile, FMT = "(I19)", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) Global_PETSc_Num ( NoBndry_DRM ( I ) ) ;
!!      End Do ;
!      Do IParts = 1, NParts ; 
!          Do I = 1, NNBndry_DRM ;
!            Node = NoBndry_DRM ( I ) ;
!            If ( ( Nodes%Rep (Node) == 1_Shrt .AND. Any ( Nodes%Locs ( Node, : ) == IParts ) ) .OR. ( Nodes%Rep (Node) /= 1_Shrt .AND. IParts == MinVal ( Nodes%Locs( Node, :), MASK = Nodes%Locs( Node, :) > 0_Shrt ) ) ) Write (Unit = UnFile, FMT = "(1(I19,1X))", ADVANCE = 'Yes', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) Global_PETSc_Num ( Node ) ;
!          End Do ;
!      End Do ;
!
!write(*,*)'debug 1005'
!
!    ! node numbers on the DRM layer
!    Write (Unit = UnFile, FMT = "()", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) ;
!    Write (Unit = UnFile, FMT = "('Layer nodes')", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) ;
!      Do IParts = 1, NParts ; 
!          Do I = 1, NNLayer_DRM ;
!            Node = NoLayer_DRM ( I ) ;
!            If ( ( Nodes%Rep (Node) == 1_Shrt .AND. Any ( Nodes%Locs ( Node, : ) == IParts ) ) .OR. ( Nodes%Rep (Node) /= 1_Shrt .AND. IParts == MinVal ( Nodes%Locs( Node, :), MASK = Nodes%Locs( Node, :) > 0_Shrt ) ) ) Write (Unit = UnFile, FMT = "(1(I19,1X))", ADVANCE = 'Yes', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) Global_PETSc_Num ( Node ) ;
!          End Do ;
!      End Do ;
!!      Do I = 1, NNLayer_DRM ;
!!        Write (Unit = UnFile, FMT = "(I19)", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) Global_PETSc_Num ( NoLayer_DRM ( I ) ) ;
!!      End Do ;
!    Close ( Unit = UnFile, Status = 'KEEP', ERR =  1002, IOSTAT = IO_File ) ;
!
!write(*,*)'debug 1006'
!
!  End If ;

! Deallocating 
DEAllocate( EqDis, EqVel, EqAcc, NodeDis, NodeVel, NodeAcc, STAT = ERR_DeAlloc ) ;
  IF ( ERR_DeAlloc /= 0 ) Then ;
    Write (*, Fmt_DEALLCT) ERR_DeAlloc ;  Write (UnInf, Fmt_DEALLCT) ERR_DeAlloc ;
    Write(*, Fmt_FL) ;  Write(UnInf, Fmt_FL) ; Write(*, Fmt_End) ; Read(*,*) ;  STOP ;
  End If ;

Write(*    ,*) 'End Subroutine < OUTPUT >' ;
Write(UnInf,*) 'End Subroutine < OUTPUT >' ;
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


End Subroutine OUTPUT_Binary ;

End Module Partition_Output_Binary  ;
