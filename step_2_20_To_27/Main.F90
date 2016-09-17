
! ***************************************************************************************************************************************************
! Name   : 20to27                                                                                                                                  **
!                                                                                                                                                  **
! DEVELOPER : Arash Fathi                                                                                                                          **
!             Babak Poursartip                                                                                                                     **
!                                                                                                                                                  **
! AMIRKABIR UNIVERSITY OF TECHNOLOGY (TEHRAN POLYTECHNIC), TEHRAN, IRAN                    ADVISOR : Vahid Lotfi                                   **
! THE UNIVERSITY OF TEXAS AT AUSTIN, TX, USA                                               ADVISOR : Loukas F. Kallivokas                          **
!                                                                                                                                                  **
! CIVIL, ARCHITECTRAL AND ENVIRONMENTAL ENGINEERNIG - STRUCTURAL ENGINEERING                                                                       **
! Institute for Computational and Engineering Science                                                                                              **
!                                                                                                                                                  **
! HISTORY : START    Aug 2012                                                                                                                      **
!                    Jan 2013                                                                                                                      **
!                    Aug 2013 Main algorithm modified                                                                                              **
!                                                                                                                                                  **
! Last Update: 07 August 2013                                                                                                                      **
!                                                                                                                                                  **
! Description: This code transforms the 20 node elements to 27 node elements for spectral element methods.                                         **
!              The algorithm used in the code is good for any 3D Hexahedral model, structured or unstructured. Number of adjacent elements         **
!              (NAdjacentElement) must be modified for unstrucutred models.                                                                        **
!                                                                                                                                                  **
! arash: 25 October 2013
! ***************************************************************************************************************************************************

Program El20ToEl27;

Use Copyfile ;
Use IFPORT ;
Use ShapeFunctions ;

Implicit None ;

! =========================== Global Variables ======================================================================================================
! =========================== Kind OF Integer AND Real ==============================================================================================
Integer(2), PARAMETER  :: SGL  = SELECTED_Real_Kind ( P = 6 , R = 37  ) ;  ! EQUIVALENT TO Real (4) 
Integer(2), PARAMETER  :: DBL  = SELECTED_Real_Kind ( P = 13, R = 200 ) ;  ! EQUIVALENT TO Real (8) 

Integer(2), PARAMETER  :: Tiny = SELECTED_INT_Kind  ( 1  ) ;               ! EQUIVALENT TO Integer (1)
Integer(2), PARAMETER  :: Smll = SELECTED_INT_Kind  ( 3  ) ;               ! EQUIVALENT TO Integer (2) 
Integer(2), PARAMETER  :: Shrt = SELECTED_INT_Kind  ( 8  ) ;               ! EQUIVALENT TO Integer (4) 
Integer(2), PARAMETER  :: Lng  = SELECTED_INT_Kind  ( 10 ) ;               ! EQUIVALENT TO Integer (8) 


! - Integer Variables -------------------------------------------------------------------------------------------------------------------------------

! Tiny Integers
Integer (Kind=Tiny)  :: NDOF ;                   ! Number of Degrees of Freedom
Integer (Kind=Tiny)  :: NDim ;                   ! Number of Dimensions
Integer (Kind=Tiny)  :: NNode ;                  ! MaxNumber of Nodes of elements in the model (in case we different element types in one model)
Integer (Kind=Tiny)  :: Analysis_Type ;          ! Used to determine the analysis type - 1:DRM 0:others
Integer (Kind=Tiny)  :: NAdjacentElement ;       ! Maximum numberof elements that can belong to

Integer (Kind=Tiny), PARAMETER  :: SEl3d27NSld   = 34 ;     ! Element: 3D - 27 noded -  Solid  ( quadrilateral ) - Spectral Element   - 6 for FEM
Integer (Kind=Tiny), PARAMETER  :: SEl3d27NMPML  = 35 ;     ! Element: 3D - 27 noded -  MPML   ( quadrilateral ) - Spectral Element   - 106 for SEM

! Smll Integers
Integer (Kind=Smll)  :: UnFile ;                 ! Holds Unit of a file for error message
Integer (Kind=Smll)  :: IO_File ;                ! Used for IOSTAT - Input Output Status - in the OPEN cammand.

Integer (Kind=Smll), PARAMETER  :: Un_ADR          = 500  ;            ! Unit number of address file to save the Name of input file and directories (.txt)
Integer (Kind=Smll), PARAMETER  :: UnInptXYZ       = 502  ;            ! the Unit number of the Input file for node coordinates (.XYZ)
Integer (Kind=Smll), PARAMETER  :: UnInptCnn       = 503  ;            ! the Unit number of the Input file for connectivities of elements (.Cnn)
Integer (Kind=Smll), PARAMETER  :: UnInptCnt       = 504  ;            ! the Unit number of the Input file for node constraints (.Cnt)
Integer (Kind=Smll), PARAMETER  :: UnInptDRM       = 509  ;            ! the Unit number of the input file for DRM nodes  (.DRM)
Integer (Kind=Smll), PARAMETER  :: UnInptDRMElm    = 510  ;            ! the Unit number of the input file for DRM nodes  (.DRMElm)

Integer (Kind=Smll), PARAMETER  :: UnInf           = 600  ;            ! the Unit number of the information file (.inf)
Integer (Kind=Smll), PARAMETER  :: UnOutXYZ        = 602  ;            ! the Unit number of the output file for coordinates of PIC  (.XYZ)
Integer (Kind=Smll), PARAMETER  :: UnOutCnn        = 603  ;            ! the Unit number of the output file for connectivities of PIC  (.Cnn)
Integer (Kind=Smll), PARAMETER  :: UnOutCnt        = 604  ;            ! the Unit number of the output file for constraints of PIC  (.Cnt)
Integer (Kind=Smll), PARAMETER  :: UnOutDRM        = 605  ;            ! the Unit number of the input file for DRM nodes  (.DRM)

! Shrt Integers
Integer (Kind=Shrt)  :: Counter ;                ! counter

! Lng Integers
Integer (Kind=Lng )  :: NEl ;                    ! Total Number of ELements of the model (on each rank)
Integer (Kind=Lng )  :: NJ ;                     ! Number of Joints (Nodes) of on this Rank
Integer (Kind=Lng )  :: I, J, K, L ;             ! Loop indices
Integer (Kind=Lng )  :: IEl, JEl, IJ, IDim ;     ! Loop indices
Integer (Kind=Lng )  :: NJN ;                    ! Number of new nodes (joints) need to be add to the available ones
Integer (Kind=Lng )  :: Node, Node_0 ;           ! Temporary variables to save node numbers
Integer (Kind=Lng )  :: Node_1, Node_2 ;         ! Temporary variables to save node numbers
Integer (Kind=Lng )  :: LocalNode_1, LocalNode_2;! Temporary variables to save node numbers
Integer (Kind=Lng )  :: EType ;                  ! Element Type
Integer (Kind=Lng )  :: NEL_DRM ;                ! Number of elements on the DRM boundary
Integer (Kind=Lng )  :: NNBndry ;                ! Number of nodes on the DRM boundary
Integer (Kind=Lng )  :: NNLayer ;                ! Number of nodes on the DRM Layer
Integer (Kind=Lng )  :: Temp ;                   ! Temporary variable
Integer (Kind=Lng )  :: NNewBound ;              ! Number of new (generated nodes due to the spectral element) nodes on the DRM boundary
Integer (Kind=Lng )  :: NNewLayer ;              ! Number of new (generated nodes due to the spectral element) nodes on the DRM layer
Integer (Kind=Lng )  :: INode ;                  ! Loop index on the node number
Integer (Kind=Lng )  :: NodeCounter ;            ! Counter on node
Integer (Kind=Lng )  :: NJCorner ;               ! Number of corner nodes (Comes from mesh generator-Ansys)
Integer (Kind=Lng )  :: Element ;                ! Temporary variable to save the element number
Integer (Kind=Lng )  :: NewNode ;                ! Temporary variable to save the node number

! - Real Variables ----------------------------------------------------------------------------------------------------------------------------------
Real (Kind=DBL)      :: X1, X2 ;                 ! TIME Variables for total runtime
Real (Kind=DBL)      :: Tol ;                    ! Tol for matching nodes
Real (Kind=DBL)      :: PMLX1, PMLX2, PMLY1, PMLY2, PMLZ1 ;  ! Tol for matching nodes
Real (Kind=DBL)      :: DetJ ;                    ! Determinant of Jacobian
Real (Kind=DBL)      :: Dist ;                    ! Distance between nodes
Real (Kind=DBL)      :: DRMxL, DRMxR, DRMyL, DRMyR, DRMzB ;  ! DRM boundaries
Real (Kind=DBL)      :: X, Y, Z

! - Logical Variable --------------------------------------------------------------------------------------------------------------------------------
! - Integer Arrays ----------------------------------------------------------------------------------------------------------------------------------
! - Real Arrays -------------------------------------------------------------------------------------------------------------------------------------

! - Integer ARRAY ALLOCATION ------------------------------------------------------------------------------------------------------------------------
! Tiny Integers
!Integer (Kind=Tiny), Allocatable, Dimension(:,:) :: DRM_Nodes ;      ! check matrix for nodes on the DRM layer
Integer (Kind=Tiny), Allocatable, Dimension(:  ) :: DRM_Nodes ;      ! check matrix for nodes on the DRM layer

! Smll Integers
Integer (Kind=Smll), Allocatable, Dimension(:)   :: LTEL ;           ! Load Type of each Element
Integer (Kind=Smll), Allocatable, Dimension(:)   :: MTEL ;           ! Material Type of each Element
Integer (Kind=Smll), Allocatable, Dimension(:)   :: ELT ;            ! Element Type of each Element

! Shrt Integers
Integer (Kind=Shrt), Allocatable, Dimension(:)   :: ELGR ;           ! Holds Group number of each element.
Integer (Kind=Shrt), Allocatable, Dimension(:)   :: DRM_NewBound ;   ! Holds new node numbers on the DRM boundary
Integer (Kind=Shrt), Allocatable, Dimension(:)   :: DRM_NewLayer ;   ! Holds new node numbers on the DRM layer

! Lng Integers
Integer (Kind=Lng ), Allocatable, Dimension(:  ) :: DRMElm ;         ! Element Connectivity

Integer (Kind=Lng ), Allocatable, Dimension(:,:) :: INOD ;           ! Node Connectivity.
Integer (Kind=Lng ), Allocatable, Dimension(:,:) :: ID ;             ! Identification matrix - Global constrains of the nodes.
Integer (Kind=Lng ), Allocatable, Dimension(:,:) :: IDF ;            ! Identification matrix for new generated elements- Global constrains of the nodes.
Integer (Kind=Lng ), Allocatable, Dimension(:  ) :: CornerNodeMapping;! Holds temporary node numbering for corner nodes of elements (We only corner nodes to identify sides of each element)
Integer (Kind=Lng ), Allocatable, Dimension(:,:) :: ElementCorner ;  ! Holds element numbers attached to each corner node

! - Real ARRAY ALLOCATION ---------------------------------------------------------------------------------------------------------------------------
Real (Kind=DBL)   :: DJ(3,3) ;            ! Node Coordinates.

Real (Kind=DBL)   , Allocatable, Dimension(:,:)  :: XYZ ;            ! Node Coordinates.
Real (Kind=DBL)   , Allocatable, Dimension(:,:)  :: XT ;             ! Coordinates of nodes of an element
Real (Kind=DBL)   , Allocatable, Dimension(:)    :: FN ;             ! Coordinates of nodes of an element
Real (Kind=DBL)   , Allocatable, Dimension(:)    :: XN ;             ! Coordinates of new nodes
Real (Kind=DBL)   , Allocatable, Dimension(:,:)  :: XYZN ;           ! Node Coordinates for new generated nodes for each element

! - Type ALLOCATION ---------------------------------------------------------------------------------------------------------------------------------
! - 20-Node 3D Hexahedral Element -------------------------------------------------------------------------------------------------------------------
! Shape Functions
Type ::  SF_3_20 ;
  Real (Kind=DBL) :: X1,  X2, X3 ;
  Real (Kind=DBL), Dimension ( 20 )  :: FN ;
End Type SF_3_20 ;

! Differentials
Type ::  DSF_3_20 ;
  Real (Kind=DBL) :: X1, X2, X3 ;
  Real (Kind=DBL), Dimension ( 20, 3 )  :: DFXI ;
End Type DSF_3_20 ;

Type (  SF_3_20 ) ::  SF       ;  ! shape function
Type ( DSF_3_20 ) :: DSF       ;  ! differentials of shape function

! - Character Variables -----------------------------------------------------------------------------------------------------------------------------
Character (Kind = 1, Len = 40 ) :: Name ;        ! Name of Input file
Character (Kind = 1, Len = 40 ) :: NameAns ;     ! Name of the output of Ansys
Character (Kind = 1, Len = 100) :: InDir ;       ! Directory of input file.
Character (Kind = 1, Len = 300) :: f1, f2;       ! temp files for copying .cnt

! =============================================== OPEN EXTERNAL FILES ===============================================================================

! - Address FILE ------------------------------------------------------------------------------------------------------------------------------------
UnFile = UN_ADR ;
Open ( Unit = UnFile, FILE = 'ADDRESS.TXT', ERR =  1001, IOSTAT = IO_File, ACCESS = 'SEQUENTIAL', ACTION = 'READ', ASYNCHRONOUS = 'NO', BLANK = 'NULL', BLOCKSIZE = 0                             , DisPOSE = 'KEEP', FORM = 'FORMATTED', POSITION = 'ASIS', STATUS = 'OLD' ) ;

! - Reading the input file name and directories from "address file" in the current directory -----------------------------------------------------------
Read(UN_ADR,*) ;
Read(UN_ADR,*) Name ;          ! File name
Read(UN_ADR,*) ;
Read(UN_ADR,*) InDir ;         ! Input files directories
Read(UN_ADR,*) ;
Read(UN_ADR,*) NDIM ;          ! Dimension of the model - For this case always 3
Read(UN_ADR,*) ;
Read(UN_ADR,*) NDOF ;          ! Number of degrees of freedom - For this case always 9
Read(UN_ADR,*) ;
Read(UN_ADR,*) NNode ;         ! Number of nodes of element -  For this case always 20
Read(UN_ADR,*) ;
Read(UN_ADR,*) Analysis_Type ; ! 0: PML ANALYSIS - 1: DRM ANALYSIS
Read(UN_ADR,*) ;
Read(UN_ADR,*) NEl ;           ! Number of Elements
Read(UN_ADR,*) ;
Read(UN_ADR,*) NJ ;            ! Number of Joints
Read(UN_ADR,*) ;
Read(UN_ADR,*) NJCorner ;      ! Number of Corner Nodes (Joints) - Comes from mesh generator Ansys
Read(UN_ADR,*) ;
Read(UN_ADR,*) Tol ;           ! Tolerance
Read(UN_ADR,*) ;
Read(UN_ADR,*) PMLX1 ;         ! PML boundaries
Read(UN_ADR,*) PMLX2 ;         ! PML boundaries
Read(UN_ADR,*) ;
Read(UN_ADR,*) PMLY1 ;         ! PML boundaries
Read(UN_ADR,*) PMLY2 ;         ! PML boundaries
Read(UN_ADR,*) ;
Read(UN_ADR,*) PMLZ1 ;         ! PML boundaries
  If ( Analysis_Type == 1_Tiny ) Then ;
    Read(UN_ADR,*) ;
    Read(UN_ADR,*) NNBndry ;       ! Number of nodes on the boundary of DRM layer - null for other analysis
    Read(UN_ADR,*) ;
    Read(UN_ADR,*) NNLayer ;
    Read(UN_ADR,*) ;
    Read(UN_ADR,*) NEL_DRM ;
    Read(UN_ADR,*) ;
    Read(UN_ADR,*) DRMxL ;
    Read(UN_ADR,*) DRMxR ;
    Read(UN_ADR,*) ;
    Read(UN_ADR,*) DRMyL ;
    Read(UN_ADR,*) DRMyR ;
    Read(UN_ADR,*) ;
    Read(UN_ADR,*) DRMzB ;
  End If ;

!InDir = TRIM(AdjustL (InDir))//'/'//TRIM(AdjustL (NAME))//"/Model" ;     ! This path is for Windows
InDir = TRIM(AdjustL (InDir))//'/'//TRIM(AdjustL (NAME))//"/Model" ;      ! This path is for Linux
NameAns = TRIM(AdjustL (Name))//'_Ans' ;


! =============================================== Input Files =======================================================================================
! - Input file for nodes' coordinates ---------------------------------------------------------------------------------------------------------------
UnFile = UnInptXYZ ;
print*,"open XYZ"
Open ( Unit = UnFile, FILE = TRIM(NameAns)//'.XYZ', ERR =  1001, IOSTAT = IO_File, ACCESS = 'SEQUENTIAL', ACTION = 'READ', ASYNCHRONOUS = 'NO', BLANK = 'NULL', BLOCKSIZE = 0, DEFAULTFILE = TRIM(InDir), DisPOSE = 'KEEP', FORM = 'FORMATTED', POSITION = 'ASIS', STATUS = 'OLD' ) ;

! - Input file for elements' connectivities ---------------------------------------------------------------------------------------------------------
UnFile = UnInptCnn ;
Open ( Unit = UnFile, FILE = TRIM(NameAns)//'.Cnn', ERR =  1001, IOSTAT = IO_File, ACCESS = 'SEQUENTIAL', ACTION = 'READ', ASYNCHRONOUS = 'NO', BLANK = 'NULL', BLOCKSIZE = 0, DEFAULTFILE = TRIM(InDir), DisPOSE = 'KEEP', FORM = 'FORMATTED', POSITION = 'ASIS', STATUS = 'OLD' ) ;

! - Input file for nodes' constraints ---------------------------------------------------------------------------------------------------------------
!UnFile = UnInptCnt ;
!Open ( Unit = UnFile, FILE = TRIM(NameAns)//'.Cnt', ERR =  1001, IOSTAT = IO_File, ACCESS = 'SEQUENTIAL', ACTION = 'READ', ASYNCHRONOUS = 'NO', BLANK = 'NULL', BLOCKSIZE = 0, DEFAULTFILE = TRIM(InDir), DisPOSE = 'KEEP', FORM = 'FORMATTED', POSITION = 'ASIS', STATUS = 'OLD' ) ;

! - Input file for DRM nodes ------------------------------------------------------------------------------------------------------------------------
  If ( Analysis_Type == 1 ) Then ;
    UnFile = UnInptDRMElm ;
    Open ( Unit = UnFile, FILE = TRIM(NameAns)//'.DRMElm', ERR =  1001, IOSTAT = IO_File, ACCESS = 'SEQUENTIAL', ACTION = 'READ', ASYNCHRONOUS = 'NO', BLANK = 'NULL', BLOCKSIZE = 0, DEFAULTFILE = TRIM(InDir), DisPOSE = 'KEEP', FORM = 'FORMATTED', POSITION = 'ASIS', STATUS = 'OLD' ) ;
    UnFile = UnInptDRM ;
    Open ( Unit = UnFile, FILE = TRIM(NameAns)//'.DRM',    ERR =  1001, IOSTAT = IO_File, ACCESS = 'SEQUENTIAL', ACTION = 'Read', ASYNCHRONOUS = 'NO', BLANK = 'NULL', BLOCKSIZE = 0, DEFAULTFILE = TRIM(InDir), DisPOSE = 'KEEP', FORM = 'FORMATTED', POSITION = 'ASIS', STATUS = 'OLD' ) ;
  End If ;

! =============================================== Reading Data ======================================================================================
!write(*,*) 'before' 
!write(*,*) nel, nj, njcorner, ndim, ndof
NAdjacentElement = 25_Tiny ;  ! For strucutred models, increase for unstructured
Allocate ( LTEL (NEL), MTEL (NEL), ELT (NEL), ELGR (NEL), INOD (NNode+7,NEL), XYZ (NJ,NDim), XT(NDim, NNode), XYZN(NEL*5,NDim), XN(NDim), ElementCorner ( NJCorner, NAdjacentElement ), CornerNodeMapping ( NJ ) ) ; ! DRMElm (NEL_DRM) , ID (NJ,NDOF)

CornerNodeMapping ( : ) = 0_Lng ; ! Initialize vector
ElementCorner     (:,:) = 0_Lng ; ! Initialize vector
NodeCounter             = 0_Lng ;
!write(*,*) 'after' 
! Coordinates of each node
Write (*,*)"Reading Coordinates ..." ;
  DO K = 1, NJ ;
    Read  (UnInptXYZ, *) I, ( XYZ ( I, J ), J = 1, NDim ) ; 
  End Do ;

! Element Connectivities
Write (*,*)"Reading Element Connectivities ..." ;
  ! Mind nodes orientation of elements

Temp = NEl /15 ;  ! A temporary variable to show how many elemenets are done.

  DO J = 1, NEl ;
    If ( Mod(J,Temp)==0) Write(*,*) "Element: ", J ;
    Read  (UnInptCnn, *) IEl, ( INOD ( I, IEl ), I = 1, NNode-8 ), MTEL ( IEl ), ELT ( IEl ), ELGR ( IEl ) ;
    Read  (UnInptCnn, *) ( INOD ( I, IEl ), I = 13, NNode ) ;
      Do INode = 1, 8 ;  ! finding element numbers attached to corner nodes
        Node = INod ( INode, IEl ) ;
          If ( CornerNodeMapping ( Node ) == 0_Lng ) Then ;
            NodeCounter = NodeCounter + 1_Lng ;
            CornerNodeMapping ( Node ) = NodeCounter ;
              Do K = 1, NAdjacentElement ;
                If ( ElementCorner ( NodeCounter, K ) == 0 ) Then ;
                  ElementCorner ( NodeCounter, K ) = IEl ;
                  Exit ;
                End If ;
              End Do ;
          Else ;
            Node_0 = CornerNodeMapping ( Node )
              Do K = 1, NAdjacentElement ;
                If ( ElementCorner ( Node_0, K ) == 0 ) Then ;
                  ElementCorner ( Node_0, K ) = IEl ;
                  Exit ;
                End If ;
              End Do ;
          End If ;
      End Do ;
  End Do ;

!  ! Constraints
!Write (*,*)"Reading Node constraints ..." ;
!  DO J = 1, NJ ;
!    Read  (UnInptCnt, *) NODE, ( ID ( NODE, I ), I = 1, NDOF ) ;
!  End Do ;

  ! DRM nodes
  If ( Analysis_Type == 1) Then ; 
    Write (*,*)"Reading DRM Elements ..." ;
    Read  (UnInptDRMElm, *) ;
      DO J = 1, NEL_DRM ;
        Read  (UnInptDRMElm, *) DRMElm(J) ;
      End Do ;
  End If ;

! =============================================== Transforming  Finite Elements to Spectral Elements ================================================

INod ( 21:27,:) = 0_Lng ; ! Initialize connectivity vector
NJN             = 0_Lng ;

! Adding nodes the connectivity vector
Write (*,*)"Finite Element to Spectral Element ..." ;
  Do IEl = 1, NEl ;

    If ( Mod(IEl,Temp)==0) print*,"Element:", IEl ;

    ! COORDINATES OF THE ELEMENT
    ForAll ( I = 1:NDim, J = 1:NNode ) XT ( I, J ) = XYZ ( INOD(J,IEl), I ) ;

    Do I = 1, 6 ! Loop on faces of Hexahedral element
      If ( INod ( 20+I, IEl ) == 0_Lng ) Then ;

          !  Local coordinates of the new node 
          SELECT CASE (I)
             CASE(1)
                SF%X1 =  0.0_Dbl ; SF%X2 =  0.0_Dbl ; SF%X3 = -1.0_Dbl
                DSF%X1 =  0.0_Dbl ; DSF%X2 =  0.0_Dbl ; DSF%X3 = -1.0_Dbl
                Node_1 = INod ( 1, IEl ) ;
                Node_2 = INod ( 3, IEl ) ;

             CASE(2)
                SF%X1 =  1.0_Dbl ; SF%X2 =  0.0_Dbl ; SF%X3 =  0.0_Dbl
                DSF%X1 =  1.0_Dbl ; DSF%X2 =  0.0_Dbl ; DSF%X3 =  0.0_Dbl
                Node_1 = INod ( 1, IEl ) ;
                Node_2 = INod ( 6, IEl ) ;

             CASE(3)
                SF%X1 =  0.0_Dbl ; SF%X2 =  1.0_Dbl ; SF%X3 =  0.0_Dbl
                DSF%X1 =  0.0_Dbl ; DSF%X2 =  1.0_Dbl ; DSF%X3 =  0.0_Dbl
                Node_1 = INod ( 2, IEl ) ;
                Node_2 = INod ( 7, IEl ) ;

             CASE(4)
                SF%X1 = -1.0_Dbl ; SF%X2 =  0.0_Dbl ; SF%X3 =  0.0_Dbl
                DSF%X1 = -1.0_Dbl ; DSF%X2 =  0.0_Dbl ; DSF%X3 =  0.0_Dbl
                Node_1 = INod ( 4, IEl );
                Node_2 = INod ( 7, IEl ) ;

             CASE(5)
                SF%X1 =  0.0_Dbl ; SF%X2 = -1.0_Dbl ; SF%X3 =  0.0_Dbl
                DSF%X1 =  0.0_Dbl ; DSF%X2 = -1.0_Dbl ; DSF%X3 =  0.0_Dbl
                Node_1 = INod ( 1, IEl ) ;
                Node_2 = INod ( 8, IEl ) ;

             CASE(6)
                SF%X1 =  0.0_Dbl ; SF%X2 =  0.0_Dbl ; SF%X3 =  1.0_Dbl
                DSF%X1 =  0.0_Dbl ; DSF%X2 =  0.0_Dbl ; DSF%X3 =  1.0_Dbl
                Node_1 = INod ( 5, IEl ) ;
                Node_2 = INod ( 7, IEl ) ;

          END SELECT

        SF%FN    =     ShapeFunc_3_20 (  SF%X1,  SF%X2,  SF%X3 ) ;
        DSF%DFXI = Dif_ShapeFunc_3_20 ( DSF%X1, DSF%X2, DSF%X3 ) ;


        DJ = MATMUL(XT,DSF%DFXI)

        DETJ = DJ(1,1) * DJ(2,2) * DJ(3,3) + DJ(2,1) * DJ(3,2) * DJ(1,3) + DJ(3,1) * DJ(1,2) * DJ(2,3) - &
               DJ(1,3) * DJ(2,2) * DJ(3,1) - DJ(2,3) * DJ(3,2) * DJ(1,1) - DJ(3,3) * DJ(1,2) * DJ(2,1)
          IF (DETJ.LE.0) THEN
             WRITE(*,*)'|J|<0 NEGATIVE DET. JACOBIAN',iel
             READ(*,*)
             STOP
          END IF

        ! Coordinate of the new node
        XN = MATMUL ( XT, SF%FN )
        NJN = NJN + 1_Lng ;
        XYZN( NJN,1:NDim) = XN(1:NDim) ;
        INod ( 20+I, IEl ) = NJ + NJN ;

        ! Seeking other elements that share this face 
          FirstLoop: Do J = 1, NAdjacentElement ;
            If ( ElementCorner ( CornerNodeMapping(Node_1), J ) == IEl ) Cycle ;
            Do K = 1, NAdjacentElement ; 
                If ( (ElementCorner ( CornerNodeMapping(Node_1), J ) == ElementCorner ( CornerNodeMapping(Node_2), K )) .And. ElementCorner ( CornerNodeMapping(Node_2), K ) /= 0) Then ; 
                  Element = ElementCorner ( CornerNodeMapping(Node_1), J ) ;
                  LocalNode_1 = 0 ;
                  LocalNode_2 = 0 ;
                    Do L = 1, 8 ; ! Find the local node number 
                      If      ( INod ( L, Element ) == Node_1 ) Then; LocalNode_1 = L ;
                      Else If ( INod ( L, Element ) == Node_2 ) Then; LocalNode_2 = L ;
                      End If ;
                    End Do ;
                    If ( LocalNode_1 == 0 .OR. LocalNode_2 == 0 ) Then ; Write (*,*)"Something is wrong. Cannot find Local nodes. Check it out."; Read(*,*); Stop ;
                    End If ;

                  ! Find the node number
                    If      ( (LocalNode_1 == 1 .And. LocalNode_2 == 3 ) .Or. (LocalNode_1 == 3 .And. LocalNode_2 == 1 ) .Or. (LocalNode_1 == 2 .And. LocalNode_2 == 4 ) .Or. (LocalNode_1 == 4 .And. LocalNode_2 == 2 ) ) Then ; NewNode = 21 ;
                    Else If ( (LocalNode_1 == 1 .And. LocalNode_2 == 6 ) .Or. (LocalNode_1 == 6 .And. LocalNode_2 == 1 ) .Or. (LocalNode_1 == 2 .And. LocalNode_2 == 5 ) .Or. (LocalNode_1 == 5 .And. LocalNode_2 == 2 ) ) Then ; NewNode = 22 ;
                    Else If ( (LocalNode_1 == 2 .And. LocalNode_2 == 7 ) .Or. (LocalNode_1 == 7 .And. LocalNode_2 == 2 ) .Or. (LocalNode_1 == 3 .And. LocalNode_2 == 6 ) .Or. (LocalNode_1 == 6 .And. LocalNode_2 == 3 ) ) Then ; NewNode = 23 ;
                    Else If ( (LocalNode_1 == 3 .And. LocalNode_2 == 8 ) .Or. (LocalNode_1 == 8 .And. LocalNode_2 == 3 ) .Or. (LocalNode_1 == 4 .And. LocalNode_2 == 7 ) .Or. (LocalNode_1 == 7 .And. LocalNode_2 == 4 ) ) Then ; NewNode = 24 ;
                    Else If ( (LocalNode_1 == 1 .And. LocalNode_2 == 8 ) .Or. (LocalNode_1 == 8 .And. LocalNode_2 == 1 ) .Or. (LocalNode_1 == 4 .And. LocalNode_2 == 5 ) .Or. (LocalNode_1 == 5 .And. LocalNode_2 == 4 ) ) Then ; NewNode = 25 ;
                    Else If ( (LocalNode_1 == 5 .And. LocalNode_2 == 7 ) .Or. (LocalNode_1 == 7 .And. LocalNode_2 == 5 ) .Or. (LocalNode_1 == 6 .And. LocalNode_2 == 8 ) .Or. (LocalNode_1 == 8 .And. LocalNode_2 == 6 ) ) Then ; NewNode = 26 ;
                    End If; 

                  INod ( NewNode, Element ) = NJ + NJN ;
                  Exit FirstLoop ;
                End If ;
            End Do ;
          End Do FirstLoop ;

      End If ;
    End Do ;

  End Do ;

Write (*,*)"Generating 27th node for elements..." ;

  DO IEL = 1, NEL

    If ( Mod(IEl,Temp)==0) print*,"Element:", IEl ;

     INOD(27, IEL) = NJ + NJN + IEL
    ! Lagrange center node (27)
    ! COORDINATES OF THE ELEMENT
    ForAll ( I = 1:NDim, J = 1:NNode ) XT ( I, J ) = XYZ ( INOD(J,IEl), I ) ;

    SF%X1 =  0.0_Dbl ; SF%X2 =  0.0_Dbl ; SF%X3 =  0.0_Dbl
    SF%FN    =     ShapeFunc_3_20 (  SF%X1,  SF%X2,  SF%X3 ) ;

    XN = MATMUL ( XT, SF%FN )
    I = NJN + IEL
    XYZN(I, 1:NDim) = XN(1:NDim)
  END DO

DeAllocate ( ElementCorner, CornerNodeMapping ) ;

! =============================================== Constraints of the added nodes ====================================================================
Allocate ( IDF ( NJN+NEL,NDOF ) ) 

IDF (:,:) = 0_Lng ;

Write (*,*)"Define constraints for the Regular Domain ..." ;

! Define constraints for added nodes in the Regular Domian:
! if we are in RD, activate displacement dofs and deactivate stress dofs ----------------
  DO IEL = 1, NEL

    If ( Mod(IEl,Temp)==0) print*,"Element:", IEl ;

    ! update (original, untouched) restraints ID (X-Y-Sxx-Syy-Sxy)
    EType = ELT ( IEl ) ;
      If ( EType == 14 ) Then ; ! 14 ==>6
        ELT ( IEl ) = 34 ;
          DO IDim =  1, 7 ;
            I = INOD ( 20+IDim, IEL)
            I = I - NJ ;
            IDF (I, 1:3) = 0 ; ! Ux Uy Uz
            IDF (I, 4:9) = 1 ; ! Sxx Syy Szz Sxy Syz Szx
          END DO
      End If ;
  End Do ;
! Remark: We should not combine these two loops, otherwise we may fix stresses on the boundary of regular domain and PML

Write (*,*)"Define constraints for PML ..." ;

! Define constraints for added nodes in PML:
  Do IEl = 1, NEl ;

    If ( Mod(IEl,Temp)==0) print*,"Element:", IEl ;

    EType = ELT ( IEl ) ;
      If ( EType == 114 ) Then ;  ! 114 ==> 106
        ELT ( IEl ) = 35 ;
          ! a. face nodes only need to be modified, 27 is already free
          DO IDim =  1, 7 !2 * NDIM
            I = INOD ( 20+IDim, IEL)
            I = I - NJ ;

!            IDF (I, 4:9) = 0 ; ! Sxx Syy Szz Sxy Syz Szx  ! Activating Stress DOFs on the Interface for new nodes"
            IDF (I, 1:9) = 0 ; ! Sxx Syy Szz Sxy Syz Szx  ! Activating Stress DOFs on the Interface for new nodes"
!              If (  DAbs  ( XYZN (I,3) ) < Tol ) Then ;  ! Traction free condition
!                IDF (I, 6) = 1 ; ! Szz
!                IDF (I, 8) = 1 ; ! Syz
!                IDF (I, 9) = 1 ; ! Szx
!              EndIf ;

              ! constraints on sides of the box
              !If (   ( XYZN(I, 1)-Tol < PMLX1 .or. XYZN(I, 1)+Tol > PMLX2 )   .or.   ( XYZN(I, 2)-Tol<PMLY1 .or. XYZN(I, 2)+Tol>PMLY2 )   .or.   ( XYZN(I, 3)-Tol<PMLZ1 ) ) Then ;  ! modify this
!              If (   ( XYZN(I, 1)-Tol < PMLX1 .or. XYZN(I, 1)+Tol > PMLX2 )   .or.   ( XYZN(I, 2)-Tol<PMLY1 .or. XYZN(I, 2)+Tol>PMLY2 )   .or.   ( Dabs(XYZN(I, 3))+Tol> Dabs( PMLZ1 ) ) ) Then ;
!                IDF (I, 1) = 1 ; ! Ux
!                IDF (I, 2) = 1 ; ! Uy
!                IDF (I, 3) = 1 ; ! Uz
!              End If ;

          END DO ;
      End If ;
  END DO ;


! apply rigid boundary conditions on the sides of the box ------------------------------------------
DO IEL = 1, NEL

   DO IDim =  1, 2 * NDIM
      I = INOD ( 20+IDim, IEL)
      I = I - NJ ;

         If ( ( XYZN(I, 1)-Tol < PMLX1 .or. XYZN(I, 1)+Tol > PMLX2 ) .or. ( XYZN(I, 2)-Tol<PMLY1 .or. XYZN(I, 2)+Tol>PMLY2 ) .or. ( XYZN(I, 3)-Tol<PMLZ1 ) ) Then ;
            IDF (I, 1) = 1 ; ! Ux
            IDF (I, 2) = 1 ; ! Uy
            IDF (I, 3) = 1 ; ! Uz
         End if
    End Do
End Do

! =============================================== DRM ===============================================================================================
! identifying new nodes located on the DRM boundary
!  If (Analysis_Type == 1) Then ;
!
!    Write (*,*)"Working on DRM ..." ;
!    Write (*,*)"Identifying the location of nodes in DRM layer ..." ;
!
!    ! Allocate matrix to identify nodes on the drm layer: 1 : nodes on the boundary - 2: noded adjacnet to the boundary
!    Allocate ( DRM_Nodes (NEL_DRM,7) ) ;
!
!    DRM_Nodes (:,:)  = 0_Tiny ;
!
!      Do I = 1, NEL_DRM ;
!        IEl = DRMElm ( I ) ;
!
!          Do J = 1, 7 ;
!            INode = NNode + J ; ! Local node number of the new node
!            Node  = INod ( INode, IEl ) ; ! Node number of the new node in the global system 
!
!            ! coordinates of the node
!            X = XYZN ( Node - NJ, 1 ) ;
!            Y = XYZN ( Node - NJ, 2 ) ;
!            Z = XYZN ( Node - NJ, 3 ) ;
!
!              If ( DRM_Nodes ( I, J ) == 0_Tiny ) Then ;
!                If ( ((X >= DRMxL-Tol .and. X <= DRMxL+Tol    .or.    X >= DRMxR-Tol .and. X <= DRMxR+Tol)    .and.    Y >= DRMyL-Tol .and. Y <= DRMyR+Tol    .and.    Z<+DRMzB+Tol)    .or.   (( Y >= DRMyL-Tol .and. Y <= DRMyL+Tol    .or.    Y >= DRMyR-Tol .and. Y <= DRMyR+Tol )    .and.    X >= DRMxL-Tol .and. X <= DRMxR+Tol    .and.    Z<+DRMzB+Tol)    .or.    ((Z >= DRMzB-Tol .and. Z <= DRMzB+Tol)    .and.    (X >= DRMxL-Tol .and. X <= DRMxR+Tol)    .and.    (Y >= DRMyL-Tol .and. Y <= DRMyR+Tol))) Then ; ! Node is located on the DRM boundary
!                  ! node on the boundary
!                  DRM_Nodes ( I, J ) = 1_Tiny ;
!
!                Else ;
!                  ! node on the layer
!                  DRM_Nodes ( I, J ) = 2_Tiny ;
!
!                    ! eliminating other nodes
!                    Do K = I+1, NEl_DRM; 
!                      JEl = DRMElm ( K ) ;
!                        Do L = 1, 7 ;
!                          Node_0 = INod ( NNode + L, JEl ) ;
!                          If (Node_0 == Node ) DRM_Nodes ( K, L ) = 3_Tiny ;
!                        End Do ;
!                    End Do ;
!                End If ;
!              End If ;
!
!          End Do ;
!
!      End Do ;
!
!    ! number of new nodes on the boundary
!    NNewBound = Count ( DRM_Nodes == 1_Lng ) ;
!    NNewLayer = Count ( DRM_Nodes == 2_Lng ) ;
!
!    ! Allocate a matrix to save the node numbers on the boundary of the DRM and on the layer
!    Allocate ( DRM_NewBound (NNewBound), DRM_NewLayer (NNewLayer) ) ;
!
!    Write (*,*)"Saving nodes on the DRM boundary ..." ;
!
!    ! saving node numbers on the DRM boundary in an array
!    Counter = 0_Shrt ;
!      Do I = 1, NEL_DRM ;
!        Do J = 1, 7 ;
!          If (DRM_Nodes ( I, J ) == 1_Tiny) Then ;
!            Counter = Counter + 1_Shrt ;
!            DRM_NewBound (Counter) = INod ( NNode + J, DRMElm ( I )) ;
!          End If ;
!        End Do ;
!      End Do
!
!    If ( Counter /= NNewBound ) Write (*,*)'There is a mistake in the number of node on the DRM boundary'
!
!    Write (*,*)"Saving nodes on the DRM Layer ..." ;
!
!    ! saving node numbers on the DRM layer in an array
!    Counter = 0_Shrt ;
!      Do I = 1, NEL_DRM ;
!        Do J = 1, 7 ;
!          If (DRM_Nodes ( I, J ) == 2_Tiny) Then; 
!            Counter = Counter + 1_Shrt ;
!            DRM_NewLayer (Counter) = INod ( NNode + J, DRMElm ( I )) ;
!          End If ;
!        End Do ;
!      End Do
!
!    If ( Counter /= NNewLayer ) Write (*,*)'There is a mistake in the number of node on the DRM layer'
!
!    Write (*,*)"Sorting DRM nodes ..." ;
!
!    ! sorting node numbers in array - notice that the PETSc library MatGetSubMatrix only works for sorted nodes
!    Call SORTQQ (LOC(DRM_NewBound), NNewBound, SRT$INTEGER4) ;
!    Call SORTQQ (LOC(DRM_NewLayer), NNewLayer, SRT$INTEGER4) ;
!
!  End If ;

  If (Analysis_Type == 1) Then ;

    Write (*,*)"Working on DRM ..." ;
    Write (*,*)"Identifying the location of nodes in DRM layer ..." ;

    ! Allocate matrix to identify nodes on the drm layer: 1 : nodes on the boundary - 2: noded adjacnet to the boundary
    Allocate ( DRM_Nodes (NJ+NJN+NEL) ) ;

    DRM_Nodes (:)  = 0_Tiny ;

      Do I = 1, NEL_DRM ;
        IEl = DRMElm ( I ) ;
          Do J = 1, 7 ;
            INode = NNode + J ; ! Local node number of the new node
            Node  = INod ( INode, IEl ) ; ! Node number of the new node in the global system 
              If ( DRM_Nodes ( Node ) == 0_Tiny ) Then ;
                ! coordinates of the node
                X = XYZN ( Node - NJ, 1 ) ;
                Y = XYZN ( Node - NJ, 2 ) ;
                Z = XYZN ( Node - NJ, 3 ) ;
                If ( ((X >= DRMxL-Tol .and. X <= DRMxL+Tol    .or.    X >= DRMxR-Tol .and. X <= DRMxR+Tol)    .and.    Y >= DRMyL-Tol .and. Y <= DRMyR+Tol    .and.    Z<+DRMzB+Tol)    .or.   (( Y >= DRMyL-Tol .and. Y <= DRMyL+Tol    .or.    Y >= DRMyR-Tol .and. Y <= DRMyR+Tol )    .and.    X >= DRMxL-Tol .and. X <= DRMxR+Tol    .and.    Z<+DRMzB+Tol)    .or.    ((Z >= DRMzB-Tol .and. Z <= DRMzB+Tol)    .and.    (X >= DRMxL-Tol .and. X <= DRMxR+Tol)    .and.    (Y >= DRMyL-Tol .and. Y <= DRMyR+Tol))) Then ; ! Node is located on the DRM boundary
                  ! node on the boundary
                  DRM_Nodes ( Node ) = 1_Tiny ;
                Else ;
                  ! node on the layer
                  DRM_Nodes ( Node ) = 2_Tiny ;
                End If ;
              End If ;
          End Do ;
      End Do ;

    ! number of new nodes on the boundary
    NNewBound = Count ( DRM_Nodes == 1_Lng ) ;
    NNewLayer = Count ( DRM_Nodes == 2_Lng ) ;

    ! Allocate a matrix to save the node numbers on the boundary of the DRM and on the layer
    Allocate ( DRM_NewBound (NNewBound), DRM_NewLayer (NNewLayer) ) ;

    Write (*,*)"Saving nodes on the DRM boundary ..." ;

    ! saving node numbers on the DRM boundary in an array
    Counter = 0_Shrt ;
      Do I = 1, NJ + NJN + NEl;
        If (DRM_Nodes ( I ) == 1_Tiny) Then ;
            Counter = Counter + 1_Shrt ;
            DRM_NewBound (Counter) = I ;
        End If ;
      End Do

    If ( Counter /= NNewBound ) Write (*,*)'There is a mistake in the number of node on the DRM boundary'

    Write (*,*)"Saving nodes on the DRM Layer ..." ;

    ! saving node numbers on the DRM layer in an array
    Counter = 0_Shrt ;
      Do I = 1, NJ + NJN + NEl ;
        If (DRM_Nodes ( I ) == 2_Tiny) Then ;
            Counter = Counter + 1_Shrt ;
            DRM_NewLayer (Counter) = I ;
        End If ;
      End Do

    If ( Counter /= NNewLayer ) Write (*,*)'There is a mistake in the number of node on the DRM layer'

  End If ;




! =============================================== Close Input Files =================================================================================
! - ADDRESS FILE ------------------------------------------------------------------------------------------------------------------------------------
UnFile =  UN_ADR ;
Close ( Unit = UnFile, STATUS = 'KEEP', ERR = 1002, IOSTAT = IO_File )

! - Input file for nodes' coordinates ---------------------------------------------------------------------------------------------------------------
UnFile = UnInptXYZ ;
Close ( Unit = UnFile, STATUS = 'KEEP', ERR = 1002, IOSTAT = IO_File )

! - Input file for elements' connectivities ---------------------------------------------------------------------------------------------------------
UnFile = UnInptCnn ;
Close ( Unit = UnFile, STATUS = 'KEEP', ERR = 1002, IOSTAT = IO_File )

! =============================================== Output Files ======================================================================================
! - Output file for Information ---------------------------------------------------------------------------------------------------------------------
UnFile = UnInf ;
Open ( Unit = UnFile, FILE = TRIM(Name)//'.Inf', ERR =  1001, IOSTAT = IO_File, ACCESS = 'SEQUENTIAL', ACTION = 'write', ASYNCHRONOUS = 'NO', BLANK = 'NULL', BLOCKSIZE = 0, DEFAULTFILE = TRIM(InDir), DisPOSE = 'KEEP', FORM = 'FORMATTED', POSITION = 'ASIS', STATUS = 'Replace' ) ;

! - Output file for Coordinates ---------------------------------------------------------------------------------------------------------------------
! Copy Ansys Coordinate file
!f1 = TRIM(InDir)//'\'//TRIM(NameAns)//'.XYZ' ;   ! This Path is for Windows
!f2 = TRIM(InDir)//'\'//TRIM(Name)//'.XYZ' ;
f1 = TRIM(InDir)//'/'//TRIM(NameAns)//'.XYZ' ;     ! This path is for Linux
f2 = TRIM(InDir)//'/'//TRIM(Name)//'.XYZ' ;

Call copy_file ( f1, f2) ;
UnFile = UnOutXYZ ;
Open ( Unit = UnFile, FILE = TRIM(Name)//'.XYZ', ERR =  1001, IOSTAT = IO_File, ACCESS = 'SEQUENTIAL', ACTION = 'write', ASYNCHRONOUS = 'NO', BLANK = 'NULL', BLOCKSIZE = 0, DEFAULTFILE = TRIM(InDir), DisPOSE = 'KEEP', FORM = 'FORMATTED', POSITION = 'APPEND', STATUS = 'Old' ) ;

! - Output file for elements' connectivities --------------------------------------------------------------------------------------------------------
UnFile = UnOutCnn ;
Open ( Unit = UnFile, FILE = TRIM(Name)//'.Cnn', ERR =  1001, IOSTAT = IO_File, ACCESS = 'SEQUENTIAL', ACTION = 'write', ASYNCHRONOUS = 'NO', BLANK = 'NULL', BLOCKSIZE = 0, DEFAULTFILE = TRIM(InDir), DisPOSE = 'KEEP', FORM = 'FORMATTED', POSITION = 'ASIS', STATUS = 'Replace' ) ;

! - Output file for nodes' constraints --------------------------------------------------------------------------------------------------------------
! Copy Ansys constraints file
!f1 = TRIM(InDir)//'\'//TRIM(NameAns)//'.Cnt' ;     ! This path is for Windows
!f2 = TRIM(InDir)//'\'//TRIM(Name)//'.Cnt' ;
f1 = TRIM(InDir)//'/'//TRIM(NameAns)//'.Cnt' ;      ! This path is for Linux
f2 = TRIM(InDir)//'/'//TRIM(Name)//'.Cnt' ;
Call copy_file ( f1, f2) ;
UnFile = UnOutCnt ;
Open  ( Unit = UnFile, FILE = TRIM(Name)//'.Cnt', ERR =  1001, IOSTAT = IO_File, ACCESS = 'SEQUENTIAL', ACTION = 'write', ASYNCHRONOUS = 'NO', BLANK = 'NULL', BLOCKSIZE = 0, DEFAULTFILE = TRIM(InDir), DisPOSE = 'KEEP', FORM = 'FORMATTED', POSITION = 'APPEND', STATUS = 'old' ) ;

! - Output file for DRM nodes ------------------------------------------------------------------------------------------------------------------------
  If ( Analysis_Type == 1_Tiny ) Then ;
    UnFile = UnOutDRM ;
    Open ( Unit = UnFile, FILE = TRIM(Name)//'.DRM', ERR =  1001, IOSTAT = IO_File, ACCESS = 'SEQUENTIAL', ACTION = 'Write', ASYNCHRONOUS = 'NO', BLANK = 'NULL', BLOCKSIZE = 0, DEFAULTFILE = TRIM(InDir), DisPOSE = 'KEEP', FORM = 'FORMATTED', POSITION = 'Append', STATUS = 'Replace' ) ;
  End If ;

! =============================================== Writing Modified Data =============================================================================

! Coordinates of each node
Write (*,*)"Write coordinates of new nodes ..." ;
!  DO I = 1, NJ ;
!    Write  (UnOutXYZ, "(I19, 2x, <NDim>(E31.23E3,2X) )") I, ( XYZ ( I, J ), J = 1, NDim ) ; 
!  End Do ;
  DO I = 1, NJN + NEL ;
    Write  (UnOutXYZ, "(I19, 2x, <NDim>(E18.10E3,2X) )") I+NJ, ( XYZN ( I, J ), J = 1, NDim ) ; 
  End Do ;

! Nodes connectivities
Write (*,*)"Write Element Connectivities ..." ;
  ! Mind nodes orientation of elements
  DO IEl = 1, NEl ;
    Write  (UnOutCnn, "(I19, 2x, <NNode+7>(I19,2X),I5,2X,I5,2X,I5)") IEl, ( INOD ( I, IEl ), I = 1, NNode+7 ), MTEL ( IEl ), ELT ( IEl ), ELGR ( IEl ) ;
  End Do ;

! Constraints
Write (*,*)"Write Nodes' Constraints of new nodes ..." ;
!  DO Node = 1, NJ ;
!    Write  (UnOutCnt, "(I19, 2x, <NDOF>(I19,2X) )") NODE, ( ID ( NODE, I ), I = 1, NDOF ) ;
!  End Do ;
  DO Node = 1, NJN+NEL ;
    Write  (UnOutCnt, "(I19, 2x, <NDOF>(I1,1X) )") NODE+NJ, ( IDF ( NODE, I ), I = 1, NDOF ) ;
  End Do ;

  ! DRM nodes on the boundary and on the layer adjacent to the boundary
  IF ( Analysis_Type == 1 ) Then ;
    Write (*,*)"DRM Output"
    Write  (UnOutDRM, "('NNbndry_nodes:',I19)") NNBndry + NNewBound ;
    Read (UnInptDRM,*)
      Do I = 1, NNBndry ;
        Read  (UnInptDRM,*) temp
        Write (UnOutDRM,*) temp
      End Do ;

      Do I = 1, NNewBound ;
        Write (UnOutDRM,*) DRM_NewBound ( I ) ;
      End Do

    Write  (UnOutDRM, "('NNLayer_nodes:',I19)") NNLayer + NNewLayer ;
    Read (UnInptDRM,*) 
      Do i = 1, NNLayer ;
        Read  (UnInptDRM,*) temp
        Write (UnOutDRM,*) temp
      End Do ;

      Do I = 1, NNewLayer ;
        Write (UnOutDRM,*) DRM_NewLayer ( I ) ;
      End Do

  End If ;


DeAllocate ( LTEL, MTEL, ELT, ELGR, INOD, XYZ, XYZN, XN ) !, DRMElm ) ;

! - Write Informtion --------------------------------------------------------------------------------------------------------------------------------
Write (*,*)"Writing Information file ..." ;
Write  (UnInf, "('Number of Elements:',I19 )") NEL ;
Write  (UnInf, "('Total Number of Joints:',I19 )") NJ+NJN+NEL ;
Write  (UnInf, "('Number of NOdes on the DRM boundary:',I19 )") NNBndry + NNewBound ;
Write  (UnInf, "('Number of NOdes on the DRM Layer:',I19 )") NNLayer + NNewLayer ;

! =============================================== Close Output Files ================================================================================
! - Output file for Information ---------------------------------------------------------------------------------------------------------------------
UnFile = UnInf ;
Close ( Unit = UnFile, STATUS = 'KEEP', ERR = 1002, IOSTAT = IO_File )

! - Output file for Coordinates ---------------------------------------------------------------------------------------------------------------------
UnFile = UnOutXYZ ;
Close ( Unit = UnFile, STATUS = 'KEEP', ERR = 1002, IOSTAT = IO_File )

! - Output file for elements' connectivities --------------------------------------------------------------------------------------------------------
UnFile = UnOutCnn ;
Close ( Unit = UnFile, STATUS = 'KEEP', ERR = 1002, IOSTAT = IO_File )

! - Output file for nodes' constraints --------------------------------------------------------------------------------------------------------------
UnFile = UnOutCnt ;
Close ( Unit = UnFile, STATUS = 'KEEP', ERR = 1002, IOSTAT = IO_File )

! - Output file for DRM nodes ------------------------------------------------------------------------------------------------------------------------
  If ( Analysis_Type == 1_Tiny ) Then ;
    UnFile = UnInptDRMElm ;
    Close ( Unit = UnFile, STATUS = 'KEEP', ERR = 1002, IOSTAT = IO_File )
    UnFile = UnInptDRM ;
    Close ( Unit = UnFile, STATUS = 'KEEP', ERR = 1002, IOSTAT = IO_File )

    UnFile = UnOutDRM ;
    Close ( Unit = UnFile, STATUS = 'KEEP', ERR = 1002, IOSTAT = IO_File )
  End If ;

! - Information on screen 
Write(*,*)' Data is ready!' ; 
Write(*,*)' -- End of the code --' ;
Write(*,*)' Press Enter ...' ;
!Read(*,*) ;
Stop ;

! =============================================== OPEN ERRORS =======================================================================================
1001  IF ( IO_File > 0 ) Then ;
        Write(*, "('error in open files. Unit: ', I4,'  -Code: ',I4)")UnFile, IO_File; 
        Write(*, *)"Input Dir:  ", InDir ; 
        Write(*, *)"File name:  ", Name ; 
!        Read(*,*) ;
        STOP ;
      End If ;

! =============================================== Close ERRORS ======================================================================================
1002  IF ( IO_File > 0 ) Then ;
        Write(*, *) 'error in close files' ;
!        Read(*,*) ;
        STOP ;
      End If ;

End Program El20ToEl27 ;
