! what really matters in life: kindness, love, and compassion.
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                                                                                                                  ++
! Start:         6  July 2013                                                                                                                      ++
! Last Update:   11 July 2013                                                                                                                      ++
!                                                                                                                                                  ++
! Description: This Module constructs appropriate data structures and mappings required for heterogeneous materials                                ++
!                                                                                                                                                  ++
!              Subroutines                     Description                                                                                         ++
!              ==============                  ==============                                                                                      ++
!                                                                                                                                                  ++
!                                                                                                                                                  ++
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Module Heterogeneous_Material ;

Use Parameters ;
Use ifport ; 

Implicit None ;

  Interface
!    Module Procedure 
  End Interface    ;

Contains ;

! ***************************************************************************************************************************************************
!                                                                                                                                                  **
! Last Update:  7 July 2013                                                                                                                        **
! Description: both 2D and 3D in the same subroutine                                                                                               **
!                                                                                                                                                  **
! ***************************************************************************************************************************************************

Subroutine Het_Mat_Numbering (                                                                                               &
NDim, MaxNNode, NDOF, NNeighbor, NParts, NEL, NJ, NEQMTotal,                                                                 & ! Integer Variables
!                                                                                                                            & ! Real Variables
MTEL, EPart, INod, ID, NEL_Rank, NJ_Rank, Global_PETSc_Num, NEqRank, NNodeRank, Local_PETSc_Num, ID_Application, ID_PETSc,   & ! Integer Arrays
Node_Mat_ID, Node_Mat_Mapping,                                                                                               &
PML_DIM, XYZ,                                                                                                                & ! Real Arrays
ModelName, OutDir                                                                                                            & ! Characters
!Nodes                                                                                                                       & ! Type 
) ;

Implicit None ;

! =========================== Global Variables ======================================================================================================

! - Integer Variables -------------------------------------------------------------------------------------------------------------------------------
Integer (Kind=Tiny), Intent(In)    :: NDim, MaxNNode, NDOF ;
Integer (Kind=Smll), Intent(In)    :: NNeighbor ;

Integer (Kind=Shrt), Intent(In)    :: NParts ;

Integer (Kind=Lng ), Intent(In)    :: NEL, NJ ;
Integer (Kind=Lng ), Intent(In)    :: NEQMTotal ;

! - Real Variables ----------------------------------------------------------------------------------------------------------------------------------
! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
! - Integer Arrays ----------------------------------------------------------------------------------------------------------------------------------
Integer (Kind=Shrt), Intent(IN), Dimension (:  )  :: EPart ;

Integer (Kind=Lng ), Intent(IN), Dimension (:,:)  :: INod, ID ;
Integer (Kind=Lng ), Intent(In), Dimension (:  )  :: NEL_Rank, NJ_Rank, Global_PETSc_Num, NEqRank, NNodeRank ;
Integer (Kind=Lng ), Intent(In), Dimension (:,:)  :: Local_PETSc_Num, ID_Application, ID_PETSc ;

Integer (Kind=Smll), Intent(In), Dimension (:  )  :: MTEL

Integer (Kind=Smll), Intent(Out), Dimension (:  )  :: Node_Mat_ID;
Integer (Kind=Lng),  Intent(Out), Dimension (:  )  :: Node_Mat_Mapping;

! - Real Arrays -------------------------------------------------------------------------------------------------------------------------------------
Real (Kind=DBL),     Intent(In),    Dimension (:,:)  :: PML_DIM ;
Real (Kind=DBL),     Intent(In),    Dimension (:,:)  :: XYZ

! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
! - Character Variables -----------------------------------------------------------------------------------------------------------------------------
! - Logical Variables -------------------------------------------------------------------------------------------------------------------------------

! - Type DECLERATIONS -------------------------------------------------------------------------------------------------------------------------------
!Type ( NodeID ) :: Nodes ;                                          ! Node%Locs :  Holds rank numbers of which this node belongs to.
                                                                     ! Node%Rep  : Holds the number of Repeatations on the ranks (How many times this node appears on ranks), Useful for neighboring

! =========================== LOCAL Variables =======================================================================================================
! - Integer Variables -------------------------------------------------------------------------------------------------------------------------------
Integer (Kind=Tiny)  :: INode ;                  ! Counter for Node Number
Integer (Kind=Tiny)  :: NNode ;                  ! Number of Nodes of element
Integer (Kind=Tiny)  :: IDOF ;                   ! Loop index over NDOF

Integer (Kind=Smll)  :: INeighbor ;              ! Counter for number of neighbors 

Integer (Kind=Shrt)  :: IParts ;                 ! Counter for Number of partiions
Integer (Kind=Shrt)  :: ILocalRow ;              ! Numner of local equation of each rank

Integer (Kind=Lng )  :: IEL ;                    ! Counter on the number of elements
Integer (Kind=Lng )  :: IJ ;                     ! Counter on the number of joints ( all nodes )
Integer (Kind=Lng )  :: Node ;                   ! A temporary variable to hold a node number

Integer (Kind=Lng )  :: GPETScNodeNum ;          ! A counter for Global PETSc Node Numbering
Integer (Kind=Lng )  :: LPETScNodeNum ;          ! Local PETSc Node Numbering
Integer (Kind=Lng )  :: ApplicationEqN ;
Integer (Kind=Lng )  :: PETScEqN ;
Integer (Kind=Lng )  :: I, J ;

Integer (Kind=Smll)  :: LOCx, LOCy, LOCz, Lreg ;                    
Integer (Kind=Lng )  :: Mat_Indep_Num_Total ;

Integer (Kind=Lng )  :: counter
! - Real Variables ----------------------------------------------------------------------------------------------------------------------------------
Real (Kind=DBL)      :: PML_X01, PML_X02, PML_X03, PML_X04, PML_X05, PML_X06
Real (Kind=DBL)      :: diffx, diffy, diffz
Real (Kind=DBL)      :: Tol
Real (Kind=DBL)      :: x0, y0, z0               ! center point of ellipsoidal inclusion
Real (Kind=DBL)      :: a0, b0, c0
Real (Kind=DBL)      :: radius
!
Real (Kind=DBL)      :: x02, y02, z02            
Real (Kind=DBL)      :: a02, b02, c02
!
Real (Kind=DBL)      :: x03, y03, z03            
Real (Kind=DBL)      :: a03, b03, c03
!
Real (Kind=DBL)      :: z_cp, z_cs
Real (Kind=DBL)      :: z_lambda, z_mu
!

! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
! - Integer Arrays ----------------------------------------------------------------------------------------------------------------------------------
Integer (Kind=Smll)  :: Element_Region ( NEL ) ; ! Region that element belongs to
Integer (Kind=Smll)  :: Node_Region ( NJ ) ;     ! Region that node belongs to
!Integer (Kind=Smll)  :: Node_Mat_ID ( NJ ) ;     ! 0 = RD , 1 = I , 2 = PML
Integer (Kind=Lng)   :: Mat_Indep_Num ( NJ ) ;   ! assigns a number to nodes with independent material property
!Integer (Kind=Lng)   :: Node_Mat_Mapping ( NJ ) ;! maps a (dependent) node to an independent node
Integer (Kind=Lng), Allocatable, Dimension ( : ) :: Node_Mat_Mapping_Indep;

! - Real Arrays -------------------------------------------------------------------------------------------------------------------------------------
Real (Kind=DBL)      :: X_Center ( NDim ) ;      ! Position of the center node of spectral element
Real (Kind=DBL)      :: X_Node ( NDim ) ;
Real (Kind=DBL)      :: X_Interface ( NDim ) ;
Real (Kind=DBL)      :: PMat_Lambda (NJ) ;       ! heterogeneous Lambda
Real (Kind=DBL)      :: PMat_Mu (NJ) ;           ! heterogeneous Mu
!Real (Kind=DBL)      :: PMat_beta_0 (NJ) ;       ! heterogeneous beta_0 (PML)


! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
! - Character Variables -----------------------------------------------------------------------------------------------------------------------------
! - Character Variables -----------------------------------------------------------------------------------------------------------------------------
Character (20)  :: IndexRank ;   ! A variable for storing the Rank number in Character format for adding at the end of the input file Name.
Character (Kind = 1, Len = 20)  :: IndexSize ;   ! A variable for storing the Size number in Character format for adding at the end of the input file Name.
Character (Kind = 1, Len = 30 ) :: ModelName ;   ! Name of Input file
Character (Kind = 1, Len = 100) :: OutDir ;      ! Directory of output files (Results)

! - Logical Variables -------------------------------------------------------------------------------------------------------------------------------
Logical   :: Check ;                             ! Used to find the end of while loop

! - Type DECLERATIONS -------------------------------------------------------------------------------------------------------------------------------

! - Unit NUMBERS OF EXTERNAL FILES ------------------------------------------------------------------------------------------------------------------
Integer (Kind=Smll)  :: UnFile ;                 ! Unit Number of the Output file.
Integer (Kind=Smll)  :: IO_File ;                ! Used for IOSTAT - Input Output Status - in the Open cammand.
Integer (Kind=Smll)  :: IO_Write ;               ! Used for IOSTAT - Input Output Status - in the Write cammand.
Integer (Kind=Smll), PARAMETER  :: Un_Lambda           = 801  ;            ! the Unit number of lambda
Integer (Kind=Smll), PARAMETER  :: Un_Mu               = 802  ;            ! the Unit number of mu
! =========================== Subroutine CODE =======================================================================================================
Tol = 0.000001_DBL
!goto 1947
Select Case (Ndim)
! ===================================================================================================================================================
!
!
!                                           T W O      D I M E N S I O N A L 
!
!
! ===================================================================================================================================================
 Case (2)
! ===================================================================================================================================================


! - determine region of elements --------------------------------------------------------------------------------------------------------------------
! we assume that we are using spectral elements (2d: 9 node; 3d: 27 node)
! -------------------------
! \       \       \       \
! \   5   \   RD  \   1   \
! \       \   6   \       \
! -------------------------
! \       \       \       \
! \   4   \   3   \   2   \
! \       \       \       \
! -------------------------

! make sure we are using spectral elements
Write(*,*) "Begin numbering for 2D heterogeneous materials"
IF ( MaxNNode /= 9 ) Then
   Write(*,*) 'Error in subroutine Het_Mat_Numbering'
   Stop
End If

! PML starting point location	  
PML_X01     = PML_Dim(2,1)
PML_X02     = PML_Dim(1,1)
PML_X03     = PML_Dim(4,1)
PML_X04     = PML_Dim(3,1)

! make sure we are reading PML boundaries correctly
IF ( PML_X02 > PML_X01 .OR. PML_X04 > PML_X03 ) Then
   Write(*,*) 'Error in subroutine Het_Mat_Numbering'
   Stop
End If

Do IEL = 1 , NEL

   X_Center (:) = XYZ ( INod ( MaxNNode , IEL ) , : )
   LOCx = 0
   LOCy = 0

Write(*,*) "Set PML boundaries"
! find the position of the center point (right or left zones)
!---------- ---------- ---------- ---------- ----------
! the point XG(1) can be either in right or left zone (or neither)
   IF (             X_Center(1) >= PML_X01 ) THEN         ! right PML zone
      LOCx        = 1

   ELSEIF (         X_Center(1) <= PML_X02 ) THEN         ! left PML zone
      LOCx        = 2

   END IF


! (top or bottom)
!---------- ---------- ---------- ---------- ----------
! the point XG(1) can be either in top or bottom zone (or neither)
! IF (              X_Center(2) > PML_X03 ) THEN       ! top PML zone
!     diffy       = X_Center(2) - PML_X03

   IF (             X_Center(2) <= PML_X04 ) THEN         ! bottom PML zone
      LOCy        = 4

   END IF


! find the element region that we are in it
!---------- ---------- ---------- ---------- ----------
! -------------------------
! \       \       \       \
! \   5   \   RD  \   1   \
! \       \   6   \       \
! -------------------------
! \       \       \       \
! \   4   \   3   \   2   \
! \       \       \       \
! -------------------------
   IF     ( LOCx == 1 .AND. LOCy == 0 ) THEN
      Lreg = 1
   ELSEIF ( LOCx == 1 .AND. LOCy == 4 ) THEN
      Lreg = 2
   ELSEIF ( LOCx == 0 .AND. LOCy == 4 ) THEN
      Lreg = 3
   ELSEIF ( LOCx == 2 .AND. LOCy == 4 ) THEN
      Lreg = 4
   ELSEIF ( LOCx == 2 .AND. LOCy == 0 ) THEN
      Lreg = 5
   ELSEIF ( LOCx == 0 .AND. LOCy == 0 ) THEN
      Lreg = 6
   ELSE
      write(*,*) 'Error in subroutine Het_Mat_Numbering'
      write(*,'(A10,2F10.5)') 'XG =', X_Center(1), X_Center(2)
      write(*,*) 'LOC-x,-y =', LOCx, LOCy
      STOP
   END IF

   Element_Region ( IEL ) = Lreg

End Do

! - determine region of nodes -----------------------------------------------------------------------------------------------------------------------
! ---------     -------     ---------
! \       \                 \       \
! \   5   \       RD        \   1   \
! \       \       6         \       \


! ---------     -------     ---------
! \       \                 \       \
! \   4   \        3        \   2   \
! \       \                 \       \
! ---------     -------     ---------
write(*,*) "Finding nodal regions"
Node_Region (:) = 6

! 1. non-corner regions: 1, 5, 3 (boundaries need to be modified)
Do IEL = 1 , NEL

   Lreg = Element_Region ( IEL )

   Select Case (Lreg)

      Case(6)
         Cycle

      Case(2)
         Cycle ! take care later

      Case(4)
         Cycle ! take care later

      Case(1)

         Do INode = 1 , MaxNNode
            Node_Region ( INod ( INode , IEL ) ) = 1
         End Do

      Case(5)

         Do INode = 1 , MaxNNode
            Node_Region ( INod ( INode , IEL ) ) = 5
         End Do

      Case(3)

         Do INode = 1 , MaxNNode
            Node_Region ( INod ( INode , IEL ) ) = 3
         End Do

   End Select

End Do

! 2. corner regions
Do IEL = 1 , NEL

   Lreg = Element_Region ( IEL )

   Select Case (Lreg)

      Case(6)
         Cycle

      Case(2)

         Do INode = 1 , MaxNNode
            Node_Region ( INod ( INode , IEL ) ) = 2
         End Do

      Case(4)

         Do INode = 1 , MaxNNode
            Node_Region ( INod ( INode , IEL ) ) = 4
         End Do

      Case(1)
         Cycle ! has already been taken care of

      Case(5)
         Cycle ! has already been taken care of

      Case(3)
         Cycle ! has already been taken care of

   End Select

End Do

! - determine node material id: 0 = RD , 1 = Interface , 2 = PML ------------------------------------------------------------------------------------

! 1. assume all nodes belong to RD
Node_Mat_ID (:) = 0

! 2. find PML nodes
Do IEL = 1 , NEL

   Lreg = Element_Region ( IEL )
   If ( Lreg == 6 ) Cycle

   Do J = 1 , MaxNNode
      Node_Mat_ID ( Inod ( J , IEL ) ) = 2
   End Do

End Do

! 3. find interface nodes (efficiently)

Do IEL = 1 , NEL

   Lreg = Element_Region ( IEL )

   Select Case (Lreg)

      Case(6)
         Cycle

      Case(2)
         Cycle

      Case(4)
         Cycle

      Case(1)

         Do INode = 1 , MaxNNode
            X_Node (:) = XYZ ( INod ( INode , IEL ) , : )
            diffx = Dabs ( X_Node (1) - PML_X01 )
            If ( diffx < Tol ) Then
               Node_Mat_ID ( Inod ( INode , IEL ) ) = 1
            End If
         End Do

      Case(5)

         Do INode = 1 , MaxNNode
            X_Node (:) = XYZ ( INod ( INode , IEL ) , : )
            diffx = Dabs ( X_Node (1) - PML_X02 )
            If ( diffx < Tol ) Then
               Node_Mat_ID ( Inod ( INode , IEL ) ) = 1
            End If
         End Do

      Case(3)

         Do INode = 1 , MaxNNode
            X_Node (:) = XYZ ( INod ( INode , IEL ) , : )
            diffy = Dabs ( X_Node (2) - PML_X04 )
            If ( diffy < Tol ) Then
               Node_Mat_ID ( Inod ( INode , IEL ) ) = 1
            End If
         End Do

   End Select

End Do

! - Independent Material Numbering ------------------------------------------------------------------------------------------------------------------

Mat_Indep_Num (:) = 0
J = 0
Do IJ = 1 , NJ
   If ( Node_Mat_ID ( IJ ) ==  2 ) Cycle
   J = J + 1
   Mat_Indep_Num ( IJ ) = J
End Do

Mat_Indep_Num_Total = MaxVal ( Mat_Indep_Num )

! - construct Node Material Mapping: maps a (dependent) node to an independent node -----------------------------------------------------------------

Node_Mat_Mapping = 0

! 1. construct the part with independent material (regular domain, including boundaries)
Do IJ = 1 , NJ
   If ( Mat_Indep_Num ( IJ ) /= 0 ) Then
      Node_Mat_Mapping ( IJ ) = IJ
   End IF
End Do

! 2. construct the part with dependent material ( PML, excluding boundaries): move on the interface and find associated nodes
Do IJ = 1 , NJ
   If ( Node_Mat_ID ( IJ ) /= 1 ) Cycle
   Lreg = Node_Region ( IJ )
   X_Interface (:) = XYZ ( IJ , : )

   Select Case (Lreg)

      Case (1)
         Do I = 1 , NJ
            If ( Node_Region ( I ) /= 1 ) Cycle
            diffy = Dabs (  X_Interface (2) - XYZ ( I , 2 )  )
            If ( diffy < Tol ) Then
               Node_Mat_Mapping ( I ) = IJ 
            End if
         End Do

      Case (2)
         Do I = 1 , NJ
            If ( Node_Region ( I ) /= 2 ) Cycle
            Node_Mat_Mapping ( I ) = IJ 
         End Do

      Case (3)
         Do I = 1 , NJ
            If ( Node_Region ( I ) /= 3 ) Cycle
            diffx = Dabs (  X_Interface (1) - XYZ ( I , 1 )  )
            If ( diffx < Tol ) Then
               Node_Mat_Mapping ( I ) = IJ 
            End if
         End Do

      Case (4)
         Do I = 1 , NJ
            If ( Node_Region ( I ) /= 4 ) Cycle
            Node_Mat_Mapping ( I ) = IJ 
         End Do

      Case (5)
         Do I = 1 , NJ
            If ( Node_Region ( I ) /= 5 ) Cycle
            diffy = Dabs (  X_Interface (2) - XYZ ( I , 2 )  )
            If ( diffy < Tol ) Then
               Node_Mat_Mapping ( I ) = IJ 
            End if
         End Do


   End Select

End Do

! - Construct global material work vectors mapping (used for material updates in inversion) ---------------------------------------------------------

! do it later


! - Material properties for the forward problem -----------------------------------------------------------------------------------------------------

PMat_Lambda (:) = 0.0_DBL
PMat_Mu     (:) = 0.0_DBL

! 1. assign material property to independent nodes
Do IJ = 1 , NJ

   If ( Mat_Indep_Num ( IJ ) /= 0 ) Then

      X_Node (:) = XYZ ( IJ , : )

! - Modifiable --------------------------------------------------------------------------------------------------------------------------------------------------
      PMat_Lambda ( IJ ) = 500.0_dbl !X_Node (1)
      PMat_Mu     ( IJ ) = 500.0_dbl !X_Node (2)
! ---------------------------------------------------------------------------------------------------------------------------------------------------

   End If

End Do

! 2. extend material property to PML zone
Do IJ = 1 , NJ

   If ( Mat_Indep_Num ( IJ ) == 0 ) Then

      PMat_Lambda ( IJ ) = PMat_Lambda ( Node_Mat_Mapping ( IJ ) )
      PMat_Mu     ( IJ ) = PMat_Mu     ( Node_Mat_Mapping ( IJ ) )

   End If

End Do

! ===================================================================================================================================================
!
!
!                                           T H R E E      D I M E N S I O N A L 
!
!
! ===================================================================================================================================================
Case (3)
! ===================================================================================================================================================

! - determine region of elements --------------------------------------------------------------------------------------------------------------------
! we assume that we are using spectral elements (2d: 9 node; 3d: 27 node)
!                     -------------------------
!                    /|       |       |       |
!                   / |   14  |   13  |   12  |
!                  /  |       |       |       |
!                 /   -------------------------
!                /    |       |       |       |
!               /     |   17  |   16  |   15  |
!              /      |       |       |       |
!             /       ------------------------- (rear zone: LOCx = -1)
!            /
!           -------------------------
!          /|       |       |       |
!         / |   8   |   RD  |   7   |
!        /  |       |       |       |
!       /   -------------------------
!      /    |       |       |       |
!     /     |   11  |   10  |   9   |
!    /      |       |       |       |
!   /       ------------------------- (middle zone: LOCx = 0)
!  /                       
! -------------------------
! |       |       |       |
! |   3   |   2   |   1   |
! |       |       |       |
! -------------------------
! |       |       |       |
! |   6   |   5   |   4   |
! |       |       |       |
! ------------------------- (front zone: LOCx = +1)

! make sure we are using spectral elements
Write(*,*) "3D Partitioning"
IF ( MaxNNode /= 27_Tiny ) Then
   Write(*,*) 'Error in subroutine Het_Mat_Numbering - this is only written for quadratic spectral elements'
   Stop
End If
Write(*,*) "Setting PML Dimensions"
! PML starting point location (positive x: 01 - negative x: 02 - positive y: 03 - negative y: 04 - surface z: 05 - bottom z: 06)
PML_X01 = PML_Dim(2,1)
PML_X02 = PML_Dim(1,1)
PML_X03 = PML_Dim(4,1)
PML_X04 = PML_Dim(3,1)
PML_X05 = PML_Dim(6,1)                                                                          ! ??????????  C H E C K  ??????????
PML_X06 = PML_Dim(5,1)                                                                          ! ??????????  C H E C K  ??????????
Write(*,*) "Dimensions Set"

! make sure we are reading PML boundaries correctly
IF ( PML_X02 > PML_X01 .OR. PML_X04 > PML_X03 .OR. PML_X06 > PML_X05) Then
   Write(*,*) 'Error in subroutine Het_Mat_Numbering (PML starting poing)'
   Stop
End If

Do IEL = 1_Lng , NEL

   X_Center (:) = XYZ ( INod ( MaxNNode , IEL ) , : )
   LOCx = 0_Smll
   LOCy = 0_Smll
   LOCz = 0_Smll


! find the position of the center point (right or left (or neither))
!---------- ---------- ---------- ---------- ----------
   IF (             X_Center(1) >= PML_X01 ) THEN         ! right PML zone
      LOCx        =  1_Smll
   ELSEIF (         X_Center(1) <= PML_X02 ) THEN         ! left PML zone
      LOCx        = -1_Smll
   END IF


! front or rear (or neither)
!---------- ---------- ---------- ---------- ----------
   IF (             X_Center(2) >= PML_X03 ) THEN         ! front PML zone
      LOCy        =  1_Smll
   ELSEIF (         X_Center(2) <= PML_X04 ) THEN         ! rear PML zone
      LOCy        = -1_Smll
   END IF

! top or at the bottom zone (or neither)
!---------- ---------- ---------- ---------- ----------
   IF (             X_Center(3) >  PML_X05 ) THEN         ! top PML zone
      LOCz        =  1_Smll
      write(*,*) 'something is wrong in het mat'
      stop
   ELSEIF (         X_Center(3) <= PML_X06 ) THEN         ! bottom PML zone
      LOCz        = -1_Smll
   END IF


!write(*,*) locx, locy, locz
!stop
! find the element region that we are in it
!---------- ---------- ---------- ---------- ----------
   IF     ( LOCx ==  1_Smll .AND. LOCy ==  1_Smll .AND. LOCz ==  0_Smll ) THEN
                                                                                 Lreg = 1_Smll
   ELSEIF ( LOCx ==  1_Smll .AND. LOCy ==  0_Smll .AND. LOCz ==  0_Smll ) THEN
                                                                                 Lreg = 2_Smll
   ELSEIF ( LOCx ==  1_Smll .AND. LOCy == -1_Smll .AND. LOCz ==  0_Smll ) THEN
                                                                                 Lreg = 3_Smll
   ELSEIF ( LOCx ==  1_Smll .AND. LOCy ==  1_Smll .AND. LOCz == -1_Smll ) THEN
                                                                                 Lreg = 4_Smll
   ELSEIF ( LOCx ==  1_Smll .AND. LOCy ==  0_Smll .AND. LOCz == -1_Smll ) THEN
                                                                                 Lreg = 5_Smll
   ELSEIF ( LOCx ==  1_Smll .AND. LOCy == -1_Smll .AND. LOCz == -1_Smll ) THEN
                                                                                 Lreg = 6_Smll
   ELSEIF ( LOCx ==  0_Smll .AND. LOCy ==  1_Smll .AND. LOCz ==  0_Smll ) THEN
                                                                                 Lreg = 7_Smll
   ELSEIF ( LOCx ==  0_Smll .AND. LOCy == -1_Smll .AND. LOCz ==  0_Smll ) THEN
                                                                                 Lreg = 8_Smll
   ELSEIF ( LOCx ==  0_Smll .AND. LOCy ==  1_Smll .AND. LOCz == -1_Smll ) THEN
                                                                                 Lreg = 9_Smll
   ELSEIF ( LOCx ==  0_Smll .AND. LOCy ==  0_Smll .AND. LOCz == -1_Smll ) THEN
                                                                                 Lreg = 10_Smll
   ELSEIF ( LOCx ==  0_Smll .AND. LOCy == -1_Smll .AND. LOCz == -1_Smll ) THEN
                                                                                 Lreg = 11_Smll
   ELSEIF ( LOCx == -1_Smll .AND. LOCy ==  1_Smll .AND. LOCz ==  0_Smll ) THEN
                                                                                 Lreg = 12_Smll
   ELSEIF ( LOCx == -1_Smll .AND. LOCy ==  0_Smll .AND. LOCz ==  0_Smll ) THEN
                                                                                 Lreg = 13_Smll
   ELSEIF ( LOCx == -1_Smll .AND. LOCy == -1_Smll .AND. LOCz ==  0_Smll ) THEN
                                                                                 Lreg = 14_Smll
   ELSEIF ( LOCx == -1_Smll .AND. LOCy ==  1_Smll .AND. LOCz == -1_Smll ) THEN
                                                                                 Lreg = 15_Smll
   ELSEIF ( LOCx == -1_Smll .AND. LOCy ==  0_Smll .AND. LOCz == -1_Smll ) THEN
                                                                                 Lreg = 16_Smll
   ELSEIF ( LOCx == -1_Smll .AND. LOCy == -1_Smll .AND. LOCz == -1_Smll ) THEN
                                                                                 Lreg = 17_Smll
   ELSEIF ( LOCx ==  0_Smll .AND. LOCy ==  0_Smll .AND. LOCz ==  0_Smll ) THEN
                                                                                 Lreg = 18_Smll
   ELSE
      write(*,*) 'Error in subroutine Het_Mat_Numbering'
      write(*,'(A10,2F10.5)') 'XG =', X_Center(1), X_Center(2)
      write(*,*) 'LOC-x,-y =', LOCx, LOCy
      STOP
   END IF

   Element_Region ( IEL ) = Lreg

End Do

!j = 0
!do i = 1, nel
!if ( element_region(i) == 4 ) then
!j = j +1
!write(*,*) j, element_region(i)
!end if
!end do

! - determine region of nodes (helps to speed up searching while constructing Node_Mat_Mapping) -----------------------------------------------------

!                     ---------  -------  ---------
!                    /|       |           |       |
!                   / |   14  |     13    |   12  |
!                  /  |       |           |       |


!                 /   ---------  -------  ---------
!                /    |       |           |       |
!               /     |   17  |     16    |   15  |
!              /      |       |           |       |
!             /       ---------  -------  --------- (rear zone: LOCx = -1)
!            /


!           ---------  -------  ---------
!          /|       |           |       |
!         / |   8   |     RD    |   7   |
!        /  |       |     18    |       |


!       /   ---------  -------  ---------
!      /    |       |           |       |
!     /     |   11  |     10    |   9   |
!    /      |       |           |       |
!   /       ---------  -------  --------- (middle zone: LOCx = 0)
!  /                       


! ---------  -------  ---------
! |       |           |       |
! |   3   |     2     |   1   |
! |       |           |       |


! ---------  -------  ---------
! |       |           |       |
! |   6   |     5     |   4   |
! |       |           |       |
! ---------  -------  --------- (front zone: LOCx = +1)

Write(*,*) "Determining node regions"
Node_Region (:) = 18_Smll

! 1. non-corner regions (type i: wall and mattress-like) (boundaries need to be modified)
Do IEL = 1_Lng , NEL

   Lreg = Element_Region ( IEL )

   Select Case (Lreg)

      Case ( 18_Smll )                                                              ! RD
         Cycle

      Case ( 1_Smll, 3_Smll, 5_Smll, 9_Smll, 11_Smll, 12_Smll, 14_Smll, 16_Smll )   ! tubes
         Cycle ! take care in the next part

      Case ( 4_Smll, 6_Smll, 15_Smll, 17_Smll )                                     ! boxes
         Cycle ! take care in the last part

      Case( 2_Smll )
         Do INode = 1_Tiny , MaxNNode
            Node_Region ( INod ( INode , IEL ) ) = 2_Smll
         End Do

      Case( 7_Smll )
         Do INode = 1_Tiny , MaxNNode
            Node_Region ( INod ( INode , IEL ) ) = 7_Smll
         End Do

      Case( 8_Smll )
         Do INode = 1_Tiny , MaxNNode
            Node_Region ( INod ( INode , IEL ) ) = 8_Smll
         End Do

      Case( 10_Smll )
         Do INode = 1_Tiny , MaxNNode
            Node_Region ( INod ( INode , IEL ) ) = 10_Smll
         End Do

      Case( 13_Smll )
         Do INode = 1_Tiny , MaxNNode
            Node_Region ( INod ( INode , IEL ) ) = 13_Smll
         End Do

   End Select

End Do

! 2. semi-corner regions (type ii: tube-like) (boundaries need to be modified)
Do IEL = 1_Lng , NEL

   Lreg = Element_Region ( IEL )

   Select Case (Lreg)

      Case ( 18_Smll )
         Cycle

      Case ( 2_Smll, 7_Smll, 8_Smll, 10_Smll, 13_Smll )
         Cycle     ! has already been taken care of (mattress)

      Case ( 4_Smll, 6_Smll, 15_Smll, 17_Smll )
         Cycle     ! take care in the last part (box)

      Case ( 1_Smll )
         Do INode = 1_Tiny , MaxNNode
            Node_Region ( INod ( INode , IEL ) ) = 1_Smll
         End Do

      Case ( 3_Smll )
         Do INode = 1_Tiny , MaxNNode
            Node_Region ( INod ( INode , IEL ) ) = 3_Smll
         End Do

      Case ( 5_Smll )
         Do INode = 1_Tiny , MaxNNode
            Node_Region ( INod ( INode , IEL ) ) = 5_Smll
         End Do

      Case ( 9_Smll )
         Do INode = 1_Tiny , MaxNNode
            Node_Region ( INod ( INode , IEL ) ) = 9_Smll
         End Do

      Case ( 11_Smll )
         Do INode = 1_Tiny , MaxNNode
            Node_Region ( INod ( INode , IEL ) ) = 11_Smll
         End Do
  
      Case ( 12_Smll )
         Do INode = 1_Tiny , MaxNNode
            Node_Region ( INod ( INode , IEL ) ) = 12_Smll
         End Do

      Case ( 14_Smll )
         Do INode = 1_Tiny , MaxNNode
            Node_Region ( INod ( INode , IEL ) ) = 14_Smll
         End Do

      Case ( 16_Smll )
         Do INode = 1_Tiny , MaxNNode
            Node_Region ( INod ( INode , IEL ) ) = 16_Smll
         End Do

   End Select

End Do

! 3. corner regions (type iii: box-like)
Do IEL = 1_Lng , NEL

   Lreg = Element_Region ( IEL )

   Select Case (Lreg)

      Case ( 18_Smll )
         Cycle

      Case ( 2_Smll, 7_Smll, 8_Smll, 10_Smll, 13_Smll )
         Cycle ! has already been taken care of (mattress)

      Case ( 1_Smll, 3_Smll, 5_Smll, 9_Smll, 11_Smll, 12_Smll, 14_Smll, 16_Smll )
         Cycle ! has already been taken care of (tubes)

      Case ( 4_Smll )
         Do INode = 1_Tiny , MaxNNode
            Node_Region ( INod ( INode , IEL ) ) = 4_Smll
         End Do

      Case ( 6_Smll )
         Do INode = 1_Tiny , MaxNNode
            Node_Region ( INod ( INode , IEL ) ) = 6_Smll
         End Do

      Case ( 15_Smll )
         Do INode = 1_Tiny , MaxNNode
            Node_Region ( INod ( INode , IEL ) ) = 15_Smll
         End Do

      Case ( 17_Smll )
         Do INode = 1_Tiny , MaxNNode
            Node_Region ( INod ( INode , IEL ) ) = 17_Smll
         End Do

   End Select

End Do
write(*,*) "Finished setting up node regions"
!j = 0
!do i = 1, nj
!if ( node_region(i) == 4 ) then
!j = j +1
!write(*,*) j, node_region(i)
!end if
!end do

! - determine node material id: 0 = RD , 1 = Interface , 2 = PML ------------------------------------------------------------------------------------

! 1. assume all nodes belong to RD
Node_Mat_ID (:) = 0_Smll

! 2. find PML nodes (interface needs to be modified in part 3)
Do IEL = 1_Lng , NEL

   Lreg = Element_Region ( IEL )
   If ( Lreg == 18_Smll ) Cycle

   Do INode = 1_Tiny , MaxNNode
      Node_Mat_ID ( INod ( INode , IEL ) ) = 2_Smll
   End Do

End Do

! 3. find interface nodes (efficiently)

Do IEL = 1_Lng , NEL

   Lreg = Element_Region ( IEL )

   Select Case (Lreg)

      Case ( 18_Smll, 1_Smll, 3_Smll, 4_Smll, 5_Smll, 6_Smll, 9_Smll, 11_Smll, 12_Smll, 14_Smll, 15_Smll, 16_Smll, 17_Smll )
         Cycle

      Case ( 2_Smll ) ! front x
         Do INode = 1_Tiny , MaxNNode
            X_Node (:) = XYZ ( INod ( INode , IEL ) , : )
            diffx = Dabs ( X_Node (1) - PML_X01 )
            If ( diffx < Tol ) Then
               Node_Mat_ID ( INod ( INode , IEL ) ) = 1_Smll
            End If
         End Do

      Case ( 13_Smll ) ! rear x
         Do INode = 1_Tiny , MaxNNode
            X_Node (:) = XYZ ( INod ( INode , IEL ) , : )
            diffx = Dabs ( X_Node (1) - PML_X02 )
            If ( diffx < Tol ) Then
               Node_Mat_ID ( INod ( INode , IEL ) ) = 1_Smll
            End If
         End Do

      Case ( 7_Smll ) ! front y
         Do INode = 1_Tiny , MaxNNode
            X_Node (:) = XYZ ( INod ( INode , IEL ) , : )
            diffy = Dabs ( X_Node (2) - PML_X03 )
            If ( diffy < Tol ) Then
               Node_Mat_ID ( INod ( INode , IEL ) ) = 1_Smll
            End If
         End Do

      Case ( 8_Smll ) ! rear y
         Do INode = 1_Tiny , MaxNNode
            X_Node (:) = XYZ ( INod ( INode , IEL ) , : )
            diffy = Dabs ( X_Node (2) - PML_X04 )
            If ( diffy < Tol ) Then
               Node_Mat_ID ( INod ( INode , IEL ) ) = 1_Smll
            End If
         End Do

      Case ( 10_Smll ) ! lower z
         Do INode = 1_Tiny , MaxNNode
            X_Node (:) = XYZ ( INod ( INode , IEL ) , : )
            diffz = Dabs ( X_Node (3) - PML_X06 )
            If ( diffz < Tol ) Then
               Node_Mat_ID ( INod ( INode , IEL ) ) = 1_Smll
            End If
         End Do

   End Select

End Do

!j = 0
!do i = 1, nj
!if ( node_mat_id(i) == 2 ) then
!j = j +1
!write(*,*) j, node_mat_id(i)
!end if
!end do

!write(*,*) 'nj=', nj

! - Independent Material Numbering ------------------------------------------------------------------------------------------------------------------
write(*,*) "Begin material numbering"
Mat_Indep_Num (:) = 0_Lng
J = 0_Lng
Do IJ = 1_Lng , NJ
   If ( Node_Mat_ID ( IJ ) ==  2_Smll ) Cycle
   J = J + 1_Lng
   Mat_Indep_Num ( IJ ) = J
End Do

Mat_Indep_Num_Total = MaxVal ( Mat_Indep_Num )

!write(*,*) Mat_Indep_Num_Total

! - construct Node Material Mapping: maps a (dependent) node to an independent node -----------------------------------------------------------------
! - We use this repeatedly in inversion to extend the material properties from the interface to the PML (Mat_Extend)

Node_Mat_Mapping = 0_Lng

! 1. construct the part with independent material (regular domain, including interfaces)
Do IJ = 1_Lng , NJ
   If ( Mat_Indep_Num ( IJ ) /= 0_Lng ) Then ! (Node_Mat_ID = 0 or 1)
      Node_Mat_Mapping ( IJ ) = IJ
   End IF
End Do

! 2. construct the part with dependent material ( PML, excluding interfaces): move on the interface and find associated nodes to be mapped to
Do IJ = 1_Lng , NJ
   If ( Node_Mat_ID ( IJ ) /= 1_Smll ) Cycle ! move only on the interface nodes
   Lreg = Node_Region ( IJ )
   X_Interface (:) = XYZ ( IJ , : )

   Select Case (Lreg)

!---------- ---------- ---------- ---------- ---------- side tubes (vertical)
      Case ( 1_Smll )
         Do I = 1_Lng , NJ
            If ( Node_Region ( I ) /= 1_Smll ) Cycle
            diffz = Dabs (  X_Interface (3) - XYZ ( I , 3 )  )
            If ( diffz < Tol ) Then
               Node_Mat_Mapping ( I ) = IJ 
            End if
         End Do

      Case ( 3_Smll )
         Do I = 1_Lng , NJ
            If ( Node_Region ( I ) /= 3_Smll ) Cycle
            diffz = Dabs (  X_Interface (3) - XYZ ( I , 3 )  )
            If ( diffz < Tol ) Then
               Node_Mat_Mapping ( I ) = IJ 
            End if
         End Do

      Case ( 12_Smll )
         Do I = 1_Lng , NJ
            If ( Node_Region ( I ) /= 12_Smll ) Cycle
            diffz = Dabs (  X_Interface (3) - XYZ ( I , 3 )  )
            If ( diffz < Tol ) Then
               Node_Mat_Mapping ( I ) = IJ 
            End if
         End Do

      Case ( 14_Smll )
         Do I = 1_Lng , NJ
            If ( Node_Region ( I ) /= 14_Smll ) Cycle
            diffz = Dabs (  X_Interface (3) - XYZ ( I , 3 )  )
            If ( diffz < Tol ) Then
               Node_Mat_Mapping ( I ) = IJ 
            End if
         End Do

!---------- ---------- ---------- ---------- ---------- corner cubes (bottom)
      Case ( 4_Smll )
         Do I = 1_Lng , NJ
            If ( Node_Region ( I ) /= 4_Smll ) Cycle
            Node_Mat_Mapping ( I ) = IJ 
         End Do

      Case ( 6_Smll )
         Do I = 1_Lng , NJ
            If ( Node_Region ( I ) /= 6_Smll ) Cycle
            Node_Mat_Mapping ( I ) = IJ 
         End Do

      Case ( 15_Smll )
         Do I = 1_Lng , NJ
            If ( Node_Region ( I ) /= 15_Smll ) Cycle
            Node_Mat_Mapping ( I ) = IJ 
         End Do

      Case ( 17_Smll )
         Do I = 1_Lng , NJ
            If ( Node_Region ( I ) /= 17_Smll ) Cycle
            Node_Mat_Mapping ( I ) = IJ 
         End Do

!---------- ---------- ---------- ---------- ---------- side walls
      Case ( 2_Smll ) ! front x
         Do I = 1_Lng , NJ
            If ( Node_Region ( I ) /= 2_Smll ) Cycle
            diffy = Dabs (  X_Interface (2) - XYZ ( I , 2 )  )
            diffz = Dabs (  X_Interface (3) - XYZ ( I , 3 )  )
            If ( (diffy + diffz) < Tol ) Then
               Node_Mat_Mapping ( I ) = IJ
            End if
         End Do

      Case ( 13_Smll ) ! rear x
         Do I = 1_Lng , NJ
            If ( Node_Region ( I ) /= 13_Smll ) Cycle
            diffy = Dabs (  X_Interface (2) - XYZ ( I , 2 )  )
            diffz = Dabs (  X_Interface (3) - XYZ ( I , 3 )  )
            If ( (diffy + diffz) < Tol ) Then
               Node_Mat_Mapping ( I ) = IJ 
            End if
         End Do

      Case ( 7_Smll ) ! front y
         Do I = 1_Lng , NJ
            If ( Node_Region ( I ) /= 7_Smll ) Cycle
            diffx = Dabs (  X_Interface (1) - XYZ ( I , 1 )  )
            diffz = Dabs (  X_Interface (3) - XYZ ( I , 3 )  )
            If ( (diffx + diffz) < Tol ) Then
               Node_Mat_Mapping ( I ) = IJ 
            End if
         End Do

      Case ( 8_Smll ) ! rear y
         Do I = 1_Lng , NJ
            If ( Node_Region ( I ) /= 8_Smll ) Cycle
            diffx = Dabs (  X_Interface (1) - XYZ ( I , 1 )  )
            diffz = Dabs (  X_Interface (3) - XYZ ( I , 3 )  )
            If ( (diffx + diffz) < Tol ) Then
               Node_Mat_Mapping ( I ) = IJ 
            End if
         End Do

!---------- ---------- ---------- ---------- ---------- bottom tubes (horizontal)
      Case ( 5_Smll ) ! front x
         Do I = 1_Lng , NJ
            If ( Node_Region ( I ) /= 5_Smll ) Cycle
            diffy = Dabs (  X_Interface (2) - XYZ ( I , 2 )  )
            If ( diffy < Tol ) Then
               Node_Mat_Mapping ( I ) = IJ 
            End if
         End Do

      Case ( 16_Smll ) ! rear x
         Do I = 1_Lng , NJ
            If ( Node_Region ( I ) /= 16_Smll ) Cycle
            diffy = Dabs (  X_Interface (2) - XYZ ( I , 2 )  )
            If ( diffy < Tol ) Then
               Node_Mat_Mapping ( I ) = IJ 
            End if
         End Do

      Case ( 9_Smll ) ! front y
         Do I = 1_Lng , NJ
            If ( Node_Region ( I ) /= 9_Smll ) Cycle
            diffx = Dabs (  X_Interface (1) - XYZ ( I , 1 )  )
            If ( diffx < Tol ) Then
               Node_Mat_Mapping ( I ) = IJ 
            End if
         End Do

      Case ( 11_Smll ) ! rear y
         Do I = 1_Lng , NJ
            If ( Node_Region ( I ) /= 11_Smll ) Cycle
            diffx = Dabs (  X_Interface (1) - XYZ ( I , 1 )  )
            If ( diffx < Tol ) Then
               Node_Mat_Mapping ( I ) = IJ 
            End if
         End Do

!---------- ---------- ---------- ---------- ---------- base mattress
      Case ( 10_Smll )
         Do I = 1_Lng , NJ
            If ( Node_Region ( I ) /= 10_Smll ) Cycle
            diffx = Dabs (  X_Interface (1) - XYZ ( I , 1 )  )
            diffy = Dabs (  X_Interface (2) - XYZ ( I , 2 )  )
            If ( (diffx + diffy) < Tol ) Then
               Node_Mat_Mapping ( I ) = IJ 
            End if
         End Do

   End Select

End Do
Write(*,*) "Finished material numbering"
! - Material properties for the forward problem -----------------------------------------------------------------------------------------------------

PMat_Lambda (:) = 0.0_DBL
PMat_Mu     (:) = 0.0_DBL

! 1. assign material property to independent nodes
Do IJ = 1_Lng , NJ

   If ( Mat_Indep_Num ( IJ ) /= 0_Lng ) Then

      X_Node (:) = XYZ ( IJ , : )

! - Modifiable --------------------------------------------------------------------------------------------------------------------------------------
! ---------------------------------------------------------------------------------------------------------------------------------------------------
!      PMat_Lambda ( IJ ) =  500.0_DBL
!      PMat_Mu     ( IJ ) =  500.0_DBL

!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
! Inverse first example: (smooth exponensial profile)
!      PMat_Lambda ( IJ ) =  80.0_DBL + 0.45_DBL * Dabs ( X_Node (3) ) + 35.0_DBL * exp ( -( Dabs ( X_Node (3) ) - 22.50_DBL )**2 / 150.0_DBL )
!      PMat_Mu     ( IJ ) =  80.0_DBL + 0.45_DBL * Dabs ( X_Node (3) ) + 35.0_DBL * exp ( -( Dabs ( X_Node (3) ) - 22.50_DBL )**2 / 150.0_DBL )

!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
! Inverse second example: (three layered media)
!   If ( X_Node (3) >= -12.0_DBL ) Then
!      PMat_Lambda ( IJ ) = 80.0_dbl
!      PMat_Mu     ( IJ ) = 80.0_dbl
!   ElseIf ( (X_Node (3) < -12.0_DBL) .AND. (X_Node (3) >= -27.0_DBL) ) Then
!      PMat_Lambda ( IJ ) = 101.250_dbl
!      PMat_Mu     ( IJ ) = 101.250_dbl
!   Else
!      PMat_Lambda ( IJ ) = 125.0_dbl
!      PMat_Mu     ( IJ ) = 125.0_dbl
!   End If

!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
! Inverse third example: (three layered media with ellipsoidal inclusion)

! ellipsoid parameters
x0 = 7.50_dbl ; y0 = 0.00_dbl ; z0 =-12.00_dbl 
a0 = 7.50_dbl ; b0 = 5.00_dbl ; c0 =  5.50_dbl ! radius

! a. generate a layered media
   If ( X_Node (3) >= -12.0_DBL ) Then
      PMat_Lambda ( IJ ) = 80.0_dbl
      PMat_Mu     ( IJ ) = 80.0_dbl
   ElseIf ( (X_Node (3) < -12.0_DBL) .AND. (X_Node (3) >= -27.0_DBL) ) Then
      PMat_Lambda ( IJ ) = 101.250_dbl
      PMat_Mu     ( IJ ) = 101.250_dbl
   Else
      PMat_Lambda ( IJ ) = 125.0_dbl
      PMat_Mu     ( IJ ) = 125.0_dbl
   End If

! b. carve elliptic inclusion
   radius = ( (XYZ( IJ, 1) - x0 ) / a0 )**2 + ( (XYZ( IJ, 2) - y0 ) / b0 )**2 + ( (XYZ( IJ, 3) - z0 ) / c0 )**2
   If ( radius <= 1.001_DBL ) Then
      PMat_Lambda ( IJ ) = 156.8_dbl
      PMat_Mu     ( IJ ) = 156.8_dbl
   End If

!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
! Inverse fourth example: (three layers and three inclusions)

! ellipsoid parameters
!x0 =-20.0_dbl  ; y0 = -20.0_dbl ; z0 =-8.75_dbl 
!a0 = 3.76_dbl  ; b0 = 19.99_dbl ; c0 =  3.76_dbl ! radius

! ellipsoid 2 parameters
!x02 =+20.0_dbl ; y02 = +20.0_dbl ; z02 =-30.00_dbl 
!a02 = 15.01_dbl ; b02 = 7.51_dbl ; c02 =  5.01_dbl 

! ellipsoid 3 parameters
!x03 =+20.0_dbl ; y03 = -20.0_dbl ; z03 =-35.00_dbl 
!a03 = 6.26_dbl ; b03 = 6.26_dbl ; c03 =  6.26_dbl 


! a. generate the layered media
!   If ( X_Node (3) >= -15.01_dbl ) Then
!      PMat_Lambda ( IJ ) = 80.0_dbl
!      PMat_Mu     ( IJ ) = 80.0_dbl
!   ElseIf ( (X_Node (3) < -15.01_dbl) .AND. (X_Node (3) >= -30.0_dbl) ) Then
!      PMat_Lambda ( IJ ) = 101.250_dbl
!      PMat_Mu     ( IJ ) = 101.250_dbl
!   Else
!      PMat_Lambda ( IJ ) = 125.0_dbl
!      PMat_Mu     ( IJ ) = 125.0_dbl
!   End If

! b. carve inclusion 1
!   radius = ( (XYZ( IJ, 1) - x0 ) / a0 )**2 + ( (XYZ( IJ, 2) - y0 ) / b0 )**2 + ( (XYZ( IJ, 3) - z0 ) / c0 )**2
!   If ( radius <= 1.000_DBL ) Then
!      PMat_Lambda ( IJ ) = 156.8_dbl
!      PMat_Mu     ( IJ ) = 156.8_dbl
!   End If

! b. carve inclusion 2
!   radius = ( (XYZ( IJ, 1) - x02 ) / a02 )**2 + ( (XYZ( IJ, 2) - y02 ) / b02 )**2 + ( (XYZ( IJ, 3) - z02 ) / c02 )**2
!   If ( radius <= 1.001_DBL ) Then
!      PMat_Lambda ( IJ ) = 156.8_dbl
!      PMat_Mu     ( IJ ) = 156.8_dbl
!   End If

! b. carve inclusion 3
!   radius = ( (XYZ( IJ, 1) - x03 ) / a03 )**2 + ( (XYZ( IJ, 2) - y03 ) / b03 )**2 + ( (XYZ( IJ, 3) - z03 ) / c03 )**2
!   If ( radius <= 1.001_DBL ) Then
!      PMat_Lambda ( IJ ) = 80.0_dbl
!      PMat_Mu     ( IJ ) = 80.0_dbl
!   End If

!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
! Garner Valley field experiment

!  a. generate the layered media based on SASW profiles (smooth profile)
!   PMat_Mu = 1.0_dbl

!   If ( X_Node (3) >= -2.01_dbl ) Then
!      PMat_Lambda ( IJ ) = 50.0_dbl

!   ElseIf ( (X_Node (3)<  -2.01_dbl) .AND. (X_Node (3) >= -4.01_dbl) ) Then
!      PMat_Lambda ( IJ ) = 57.9_dbl

!   ElseIf ( (X_Node (3) < -4.01_dbl) .AND. (X_Node (3) >= -7.61_dbl) ) Then
!      PMat_Lambda ( IJ ) = 88.4_dbl

!   ElseIf ( (X_Node (3) < -7.61_dbl) .AND. (X_Node (3) >= -12.21_dbl) ) Then
!      PMat_Lambda ( IJ ) = 137.2_dbl

!   ElseIf ( (X_Node (3) < -12.21_dbl) .AND. (X_Node (3) >= -16.81_dbl) ) Then
!      PMat_Lambda ( IJ ) = 167.6_dbl

!   ElseIf ( (X_Node (3) < -16.81_dbl) .AND. (X_Node (3) >= -22.51_dbl) ) Then
!      PMat_Lambda ( IJ ) = 182.9_dbl

!   ElseIf ( (X_Node (3) < -22.51_dbl) .AND. (X_Node (3) >= -30.51_dbl) ) Then
!      PMat_Lambda ( IJ ) = 240.0_dbl

!   ElseIf ( (X_Node (3) < -30.51_dbl) .AND. (X_Node (3) >= -40.01_dbl) ) Then
!      PMat_Lambda ( IJ ) = 300.0_dbl

!   Else
!      PMat_Lambda ( IJ ) = 300.0_dbl

!   End If

! for inversion with field data, use a linear initial guess:
!z_cs = 100.0_dbl - 6.50_dbl * X_Node (3)
!z_cp = 180.0_dbl - 14.0_dbl * X_Node (3)

!z_mu     = 0.0020_dbl * z_cs**2
!z_lambda = 0.0020_dbl * z_cp**2 - 2.0_dbl * z_mu

!PMat_Lambda ( IJ ) = z_lambda
!PMat_Mu     ( IJ ) = z_mu

! ---------------------------------------------------------------------------------------------------------------------------------------------------
! ---------------------------------------------------------------------------------------------------------------------------------------------------

   End If

End Do

write(*,*) "Extend material property to PML zone"
! 2. extend material property to PML zone
Do IJ = 1_Lng , NJ

   If ( Mat_Indep_Num ( IJ ) == 0_Lng ) Then

      PMat_Lambda ( IJ ) = PMat_Lambda ( Node_Mat_Mapping ( IJ ) )
      PMat_Mu     ( IJ ) = PMat_Mu     ( Node_Mat_Mapping ( IJ ) )

   End If

End Do

! ===================================================================================================================================================

End Select
write(*,*) "Finished extension"

!1947 continue


!do ij = 1, nj
!  pmat_mu(ij) = 2.0d0 + dabs(dsin(dble(ij)))
!  pmat_lambda(ij) = 2.0d0 + dabs(dcos(dble(ij)))
!end do

! a. generate two layered media
!   If ( XYZ( IJ, 1) >= 0.001_DBL ) Then
!      PMat_Lambda ( IJ ) = 1.0_dbl
!      PMat_Mu     ( IJ ) = 1.0_dbl
!   Else
!      PMat_Lambda ( IJ ) = 2.0_dbl
!      PMat_Mu     ( IJ ) = 2.0_dbl
!   End if


!write(*,*) '---------- 1947 is active ----------'
! - Modifiable --------------------------------------------------------------------------------------------------------------------------------------
! <><><><><><><><><><>()()()()()()()()()()**********##########@@@@@@@@@@!!!!!!!!!!%%%%%%%%%%<><><><><><><><><><>()()()()()()()()()()**********##########@@@@@@@@@@!!!!!!!!!!%%%%%%%%%%
! Example 1 (homogeneous case)
!      PMat_Lambda ( : ) = 500.0_dbl
!      PMat_Mu     ( : ) = 500.0_dbl
! Basu
!      PMat_Lambda ( : ) = 1.0_dbl
!      PMat_Mu     ( : ) = 1.0_dbl
! ---------------------------------------------------------------------------------------------------------------------------------------------------
! <><><><><><><><><><>()()()()()()()()()()**********##########@@@@@@@@@@!!!!!!!!!!%%%%%%%%%%<><><><><><><><><><>()()()()()()()()()()**********##########@@@@@@@@@@!!!!!!!!!!%%%%%%%%%%

! Example 2: layered media with inclusion
! ellipsoid parameters
!a0 = 15.0_dbl ; b0 =  5.0_dbl ; c0 =  5.0_dbl
!x0 = 25.0_dbl ; y0 = 25.0_dbl ; z0 =-20.0_dbl 
 
! a. generate two layered media
!Do IJ = 1 , NJ
!   If ( XYZ( IJ, 3) >= -20.001_DBL ) Then
!      PMat_Lambda ( IJ ) = 320.0_dbl
!      PMat_Mu     ( IJ ) = 320.0_dbl
!   Else
!      PMat_Lambda ( IJ ) = 500.0_dbl
!      PMat_Mu     ( IJ ) = 500.0_dbl
!   End if
!End Do

! b. carve elliptic inclusion
!counter = 0
!Do IJ = 1 , NJ
!
!   radius = ( (XYZ( IJ, 1) - x0 ) / a0 )**2 + ( (XYZ( IJ, 2) - y0 ) / b0 )**2 + ( (XYZ( IJ, 3) - z0 ) / c0 )**2
!
!   If ( radius <= 1.001_DBL ) Then
!      PMat_Lambda ( IJ ) = 720.0_dbl
!      PMat_Mu     ( IJ ) = 720.0_dbl
!      counter = counter + 1
!   End if
!End Do

!write(*,*) '--------------------------------------'
!write(*,*) 'Number of the points in the ellipsoid:'
!write(*,*) counter
!write(*,*) '--------------------------------------'

! =========================== Output ================================================================================================================

! - Writing down the input data for each process ----------------------------------------------------------------------------------------------------
write(*,*)"Begin writing files for partition"
  Do IParts = 1_Shrt, NParts ; 
    
    write(IndexRank,*) IParts - 1_Shrt ; ! Converts Rank number to Character foramt for the file Name;
    Write (IndexSize, *) NParts ;          ! Converts Size number to Character foramt for the file Name
    Write (*,*)"Writing files for partition: ", IParts;

    ! - Output FILEs --------------------------------------------------------------------------------------------------------------------------------
    Write (*,*)"Opening files for Heterogeneous materials ...";

    ! Heterogeneous Material Lambda
    UnFile = UN_Lambda ;
    Open ( Unit = UnFile, FILE = TRIM(ModelName)//'_'//Trim(AdjustL(IndexSize))//'_'//Trim(AdjustL(IndexRank))//'.Lambda', ERR =  1001, IOSTAT = IO_File, ACCESS = 'SEQUENTIAL', ACTION ='Write', ASYNCHRONOUS = 'NO', BLANK = 'NULL', BLOCKSIZE = 0, DEFAULTFILE = TRIM(OutDir), DisPOSE = 'KEEP', FORM = 'FORMATTED', POSITION = 'ASIS', STATUS ='REPLACE') ;

    ! Heterogeneous Material Mu
    UnFile = UN_Mu ;
    Open ( Unit = UnFile, FILE = TRIM(ModelName)//'_'//Trim(AdjustL(IndexSize))//'_'//Trim(AdjustL(IndexRank))//'.Mu', ERR =  1001, IOSTAT = IO_File, ACCESS = 'SEQUENTIAL', ACTION ='Write', ASYNCHRONOUS = 'NO', BLANK = 'NULL', BLOCKSIZE = 0, DEFAULTFILE = TRIM(OutDir), DisPOSE = 'KEEP', FORM = 'FORMATTED', POSITION = 'ASIS', STATUS ='REPLACE') ;

! - Write Lambda & Mu -------------------------------------------------------------------------------------------------------------------------------
! comment: this is consistent with the way we write XYZ.
    Write (*,*) "Writing Heterogeneous material properties ..." ;
    UnFile = Un_Lambda ;
      DO IJ = 1_Lng , NJ ;
        If ( Local_PETSc_Num ( IJ , IParts ) /= 0_Lng ) Then ; !?? check the oreintation
          Write (Unit = UnFile, FMT = "(<NDim>(E31.23E3,2X) )", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) PMat_Lambda ( IJ ) ;
        End If ;
      End Do ;

    UnFile = Un_Mu ;
      DO IJ = 1_Lng , NJ ;
        If ( Local_PETSc_Num ( IJ , IParts ) /= 0_Lng ) Then ;
          Write (Unit = UnFile, FMT = "(<NDim>(E31.23E3,2X) )", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) PMat_Mu ( IJ ) ; 
        End If ;
      End Do ;

    ! - Closing the output file ---------------------------------------------------------------------------------------------------------------------
    UnFile =  Un_Lambda ;
    Close ( Unit = UnFile, Status = 'KEEP', ERR =  1002, IOSTAT = IO_File ) ;

    UnFile =  Un_Mu ;
    Close ( Unit = UnFile, Status = 'KEEP', ERR =  1002, IOSTAT = IO_File ) ;

  End Do


Write(*     ,*) 'End Subroutine < Het_Mat_Numbering >' ;
Write(UnInf,*) 'End Subroutine < Het_Mat_Numbering >' ;
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


End Subroutine Het_Mat_Numbering ;

End Module Heterogeneous_Material ;
