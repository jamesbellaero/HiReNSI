
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                                                                                                                  ++
! Last Update:  21 June 2011                                                                                                                       ++
!                                                                                                                                                  ++
! Description: This Module evaluates number of non-zero (nnz) elements for major PETSc objects so that allocation time (malloc) reduces.           ++
!                                                                                                                                                  ++
!              Subroutines                     Description                                                                                         ++
!              ==============                  ==============                                                                                      ++
!              Object_Size_2D_8N               Evaluates nnz for 2D, 8-node elements                                                               ++
!                                                                                                                                                  ++
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Module PETScObject_Size ;

Use IFPORT ;
Use Parameters ;
Use AssembleSub ;
Use PML_Matrices_Pattern ;

Implicit None ;

  Interface Object ;
!    Module Procedure 
  End Interface    ;

Contains ;

! ***************************************************************************************************************************************************
!                                                                                                                                                  **
! Last Update:  21 June 2012                                                                                                                       **
! Description: THIS Subroutine EVALUATES THE SIZE OF MASS, DAMPING AND STIFFNESS MATRICES TO REDUCE THE ALLOCATION COST                            **
!                                                                                                                                                  **
! ***************************************************************************************************************************************************

Subroutine Object_Size  (                                                                                                       &
NDim, NLCase, Matrix_Type, NEqEl, El_Type, NNode, PARAM_Type, NInt, NDOF, NInt_Type,     NParts,     NJ, NEL,                   & ! Integer Variables
PR, HREF,                                                                                                                       & ! Real Variables
ELGR, LTEL, MTEL, ELT, NEqRank,     Global_PETSc_Num,     INOD, ID_PETSc,     D_NNZ_Stiff, O_NNZ_Stiff, D_NNZ_Damp, O_NNZ_Damp, D_NNZ_Mass, O_NNZ_Mass,   & ! Integer Arrays
PMat, XYZ, PML_DIM                                                                                                              & ! Real Arrays
!                                                                                                                               & ! Characters
!GAUSS_PNT                                                                                                                      & ! Type 
) ;

Implicit None ;

! =========================== Global Variables ======================================================================================================
! - Integer Arrays ----------------------------------------------------------------------------------------------------------------------------------
Integer (Kind=Tiny), Intent(In)    :: NDim, Matrix_Type, NNode, PARAM_Type, NInt, NDOF, NInt_Type  ;
Integer (Kind=Smll), Intent(In)    :: NLCase, NEqEl, El_Type ;
Integer (Kind=Shrt), Intent(In)    :: NParts ;
Integer (Kind=Lng ), Intent(In)    :: NJ, NEL ;

! - Real Variables ----------------------------------------------------------------------------------------------------------------------------------
Real (Kind=DBL), Intent(In)    :: PR, HREF ;

! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------

! - Integer Arrays ----------------------------------------------------------------------------------------------------------------------------------
Integer (Kind=Shrt), Intent(In ), Dimension (:  )  :: ELGR, LTEL, MTEL, ELT ;

Integer (Kind=Shrt), Intent(Out), Dimension (:  )  :: D_NNZ_Stiff, O_NNZ_Stiff, D_NNZ_Damp, O_NNZ_Damp, D_NNZ_Mass, O_NNZ_Mass;    ! , NNZ_MB, NNZ_P ;

Integer (Kind=Lng ), Intent(In ), Dimension (:  )  :: Global_PETSc_Num, NEqRank ;
Integer (Kind=Lng ), Intent(In ), Dimension (:,:)  :: INOD, ID_PETSc ;  !?? be careful with ID 

! - Real Arrays -------------------------------------------------------------------------------------------------------------------------------------
Real (Kind=DBL), Intent(In), Dimension (:,:)  :: PMat, XYZ, PML_DIM ;

! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------

! - Character Variables -----------------------------------------------------------------------------------------------------------------------------
! - Logical Variables -------------------------------------------------------------------------------------------------------------------------------
! - Type DECLERATIONS -------------------------------------------------------------------------------------------------------------------------------
!Type ( PML_Interface ) :: PML_INTER ;            ! STIFFNESS AND MASS OF PML ELEMENTS AT THE Interface NODES.

! =========================== LOCAL Variables =======================================================================================================
! - Integer Arrays ------------------------------------------------------------------------------------------------------------------------------
Integer (Kind=Smll)  :: I, J, K ;                ! Loop indices.
Integer (Kind=Smll)  :: MaxEl ;                  ! Maximum number of elements that might attach to a node in the model. By default it is set to 50. The code will automatically increases this number in case there are more elements and repeats calculations.
Integer (Kind=Smll)  :: ElementCounter ;         ! Counter on number of elements.
Integer (Kind=Smll)  :: IEL_Neigh ;              ! Loop counter for all elements having the main node
Integer (Kind=Smll)  :: NNeighbors ;             ! Holds total number of nodes connected to the main node.
Integer (Kind=Smll)  :: INode_Main ;             ! Holds local node number of the main node.
Integer (Kind=Smll)  :: INode;                   ! Loop counter for element node number.
Integer (Kind=Smll)  :: IDOF ;                   ! Loop counter on the degree of freedom.
Integer (Kind=Smll)  :: Location ;               ! Holds the location of the stiffness/mass/damping entirs of a particular neighbor of the main node in NNZ_Node_Stiff/NNZ_Node_Mass/NNZ_Node_Damp.
Integer (Kind=Smll)  :: ERR_Alloc, ERR_DeAlloc ; ! ERRORS.

Integer (Kind=Shrt)  :: LType ;                  ! Load Type of Element
Integer (Kind=Shrt)  :: MType ;                  ! Material Type of Element
Integer (Kind=Shrt)  :: IParts ;                 ! Loop counter for the number of partitions.

Integer (Kind=Lng )  :: LowerEq ;                ! Lower equation number of the diagonal on this rank  !?? clarify the Name
Integer (Kind=Lng )  :: UpperEq ;                ! Upper equation number of the diagonal on this rank  !?? clarify the Name
Integer (Kind=Lng )  :: Node ;                   ! Saves node number.
Integer (Kind=Lng )  :: IEL ;                    ! Loop counter for elements.
Integer (Kind=Lng )  :: NEq ;                    ! Equation number.
Integer (Kind=Lng )  :: IJ                       ! Loop index on NJ


! - Real Variables ----------------------------------------------------------------------------------------------------------------------------------
Real (Kind=DBL)  :: Lambda, MU ;                 ! Lame Coefficients
Real (Kind=DBL)  :: Rho ;                        ! Density
Real (Kind=DBL)  :: ALFA, BETTA ;                ! Rayleigh damping coefficients
Real (Kind=DBL)  :: R_PML_ALFA_0 ;               ! Power of decay or alfa_0
Real (Kind=DBL)  :: C_REF_BETA_0 ;               ! P-wave velocity or beta_0
Real (Kind=DBL)  :: M_PML ;                      ! Characteristic length of the domain
Real (Kind=DBL)  :: B_PML ;                      ! User tunable reflection coefficient

! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
! - Integer Arrays ----------------------------------------------------------------------------------------------------------------------------------
Integer (Kind=Tiny)  :: KE     ( NEqEl, NEqEl ) ;               ! Holds stiffness matrix of one element
Integer (Kind=Tiny)  :: ME     ( NEqEl, NEqEl ) ;               ! Holds mass matrix of one element
Integer (Kind=Tiny)  :: CE     ( NEqEl, NEqEl ) ;               ! Holds damping matrix of one element
Integer (Kind=Tiny)  :: KE_RD  ( NEqEl, NEqEl ) ;               ! Holds stiffness matrix of regular domain element
Integer (Kind=Tiny)  :: ME_RD  ( NEqEl, NEqEl ) ;               ! Holds mass matrix of regular domain element
Integer (Kind=Tiny)  :: CE_RD  ( NEqEl, NEqEl ) ;               ! Holds damping matrix of regular domain element
Integer (Kind=Tiny)  :: KE_PML ( NEqEl, NEqEl ) ;               ! Holds stiffness matrix of PML element
Integer (Kind=Tiny)  :: ME_PML ( NEqEl, NEqEl ) ;               ! Holds mass matrix of PML element
Integer (Kind=Tiny)  :: CE_PML ( NEqEl, NEqEl ) ;               ! Holds damping matrix of PML element

! - Integer Array Allocation ------------------------------------------------------------------------------------------------------------------------
Integer (Kind=Tiny), Allocatable, Dimension (:,:) :: NNZ_Node_Stiff ;                              ! Holds all non-zero entries of the main node equations in the total stiffness matrix
Integer (Kind=Tiny), Allocatable, Dimension (:,:) :: NNZ_Node_Mass ;                               ! Holds all non-zero entries of the main node equations in the total mass matrix
Integer (Kind=Tiny), Allocatable, Dimension (:,:) :: NNZ_Node_Damp ;                               ! Holds all non-zero entries of the main node equations in the total damping matrix

Integer (Kind=Lng ), Allocatable, Dimension (:)   :: El_Neighbors ;                                ! Holds element numbers having the main node as their node.
Integer (Kind=Lng ), Allocatable, Dimension (:)   :: Node_Connect_Total ;                          ! Holds all node numbers having a connection with the main node (Neighbor Nodes). Since we do not know how many nodes connected to the main node, we consider maximum possible number.
Integer (Kind=Lng ), Allocatable, Dimension (:)   :: Node_Connect ;                                ! Holds all node numbers having a connection with the main node (Neighbor Nodes). This time, we know exact number of neighbors.

! - Real Arrays -------------------------------------------------------------------------------------------------------------------------------------
Real (Kind=DBL)  :: KE_Exact ( NEqEl, NEqEl ) ;  ! Holds the exact stiffness matrix of any element
Real (Kind=DBL)  :: ME_Exact ( NEqEl, NEqEl ) ;  ! Holds the exact mass matrix of any element
Real (Kind=DBL)  :: CE_Exact ( NEqEl, NEqEl ) ;  ! Holds the exact damping matrix of any element

! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
! - Character Variables -----------------------------------------------------------------------------------------------------------------------------
! - Logical Variables -------------------------------------------------------------------------------------------------------------------------------
! - Type DECLERATIONS -------------------------------------------------------------------------------------------------------------------------------
Type ( GAUSS   )       :: GAUSS_PNT ;            ! Gauss Points

! =========================== Subroutine CODE =======================================================================================================

! <<<< IMPORTANT NOTICE: MIND the ordering of the equations in this code and the original code, make a notification >>> !?? 

! LOAD GAUSS INTEGRATION DATA POINTS
GAUSS_PNT = GAUSS_POINTS( NInt, NInt_Type ) ;

! Initializing object size arrays
D_NNZ_Stiff = 0_Shrt ;
O_NNZ_Stiff = 0_Shrt ;
D_NNZ_Damp  = 0_Shrt ;
O_NNZ_Damp  = 0_Shrt ;
D_NNZ_Mass  = 0_Shrt ;
O_NNZ_Mass  = 0_Shrt ;

  ! Approximate pattern of mass, damping and stiffness matrices. For approximate esttimate "Matrix_Type" is 0 and for exact estimation "Matrix_Type" is 1 .
  If ( Matrix_Type == 0_Smll ) Then ;

    Select Case ( El_Type ) ;

      Case ( 2, 4) ; ! 8-noded quadrilateral 2D element
        ! - Stiffness matrix ----------------------------------------------------------------------------------------------------------------------------
        ! - Mind the orientation of DOFs

        ! Regualar Domain
        KE_RD = 0_Tiny ;
        ForAll ( I =             1: 2 * NNode, J =             1: 2 * NNode ) KE_RD ( I, J ) = 1_Tiny ;

        ! PML
        KE_PML = 1_Tiny ;
        ForAll ( I =             1:     NNode, J = 3 * NNode + 1: 4 * NNode ) KE_PML ( I, J ) = 0_Tiny ;    ! 1
        ForAll ( I =     NNode + 1: 2 * NNode, J = 2 * NNode + 1: 3 * NNode ) KE_PML ( I, J ) = 0_Tiny ;    ! 2
        ForAll ( I = 2 * NNode + 1: 3 * NNode, J =     NNode + 1: 2 * NNode ) KE_PML ( I, J ) = 0_Tiny ;    ! 3
        ForAll ( I = 2 * NNode + 1: 3 * NNode, J = 4 * NNode + 1: 5 * NNode ) KE_PML ( I, J ) = 0_Tiny ;    ! 4
        ForAll ( I = 3 * NNode + 1: 4 * NNode, J =             1:     NNode ) KE_PML ( I, J ) = 0_Tiny ;    ! 5
        ForAll ( I = 3 * NNode + 1: 4 * NNode, J = 4 * NNode + 1: 5 * NNode ) KE_PML ( I, J ) = 0_Tiny ;    ! 6
        ForAll ( I = 4 * NNode + 1: 5 * NNode, J = 2 * NNode + 1: 4 * NNode ) KE_PML ( I, J ) = 0_Tiny ;    ! 7

        ! - Mass matrix ---------------------------------------------------------------------------------------------------------------------------------
        ! - Mind the orientation of DOFs

        ! Regualar Domain
        ME_RD = 0_Tiny ;
        ForAll ( I =             1:     NNode, J =             1:     NNode ) ME_RD ( I, J ) = 1_Tiny ;
        ForAll ( I =     NNode + 1: 2 * NNode, J =     NNode + 1: 2 * NNode ) ME_RD ( I, J ) = 1_Tiny ;

        ! PML
        ME_PML = 0_Tiny ;
        ForAll ( I =             1:     NNode, J =             1:     NNode ) ME_PML ( I, J ) = 1_Tiny ;
        ForAll ( I =     NNode + 1: 2 * NNode, J =     NNode + 1: 2 * NNode ) ME_PML ( I, J ) = 1_Tiny ;
        ForAll ( I = 2 * NNode + 1: 4 * NNode, J = 2 * NNode + 1: 4 * NNode ) ME_PML ( I, J ) = 1_Tiny ;
        ForAll ( I = 4 * NNode + 1: 5 * NNode, J = 4 * NNode + 1: 5 * NNode ) ME_PML ( I, J ) = 1_Tiny ;

        ! - Damping matrix ------------------------------------------------------------------------------------------------------------------------------
        ! - Mind the orientation of DOFs

        ! Regualar Domain
        CE_RD = 0_Tiny ;
        ForAll ( I =             1: 2 * NNode, J =             1: 2 * NNode ) CE_RD ( I, J ) = 1_Tiny ;

        ! PML
        CE_PML = 1_Tiny ;
        ForAll ( I =             1:     NNode, J = 3 * NNode + 1: 4 * NNode ) CE_PML ( I, J ) = 0_Tiny ;    ! 1
        ForAll ( I =     NNode + 1: 2 * NNode, J = 2 * NNode + 1: 3 * NNode ) CE_PML ( I, J ) = 0_Tiny ;    ! 2
        ForAll ( I = 2 * NNode + 1: 3 * NNode, J =     NNode + 1: 2 * NNode ) CE_PML ( I, J ) = 0_Tiny ;    ! 3
        ForAll ( I = 2 * NNode + 1: 3 * NNode, J = 4 * NNode + 1: 5 * NNode ) CE_PML ( I, J ) = 0_Tiny ;    ! 4
        ForAll ( I = 3 * NNode + 1: 4 * NNode, J =             1:     NNode ) CE_PML ( I, J ) = 0_Tiny ;    ! 5
        ForAll ( I = 3 * NNode + 1: 4 * NNode, J = 4 * NNode + 1: 5 * NNode ) CE_PML ( I, J ) = 0_Tiny ;    ! 6
        ForAll ( I = 4 * NNode + 1: 5 * NNode, J = 2 * NNode + 1: 4 * NNode ) CE_PML ( I, J ) = 0_Tiny ;    ! 7

      !Case () ;  ! 20-noded hexahedral 3D element
      !Case () ;  ! 20-noded hexahedral 3D element

    End Select ;


  End If ;


MaxEl = 50_Smll ; ! Maximum number of elements a node can attach to in an unstructured mesh (An initial guess).

Allocate ( El_Neighbors ( MaxEl ),     STAT = ERR_Alloc ) ; 
  If ( ERR_Alloc /= 0 ) Then ;
    Write (*, Fmt_ALLCT) ERR_Alloc ;  Write (UnInf, Fmt_ALLCT) ERR_Alloc ;
    Write(*, Fmt_FL) ;  Write(UnInf, Fmt_FL) ; Read(*, Fmt_End) ;  Stop ;
  End If ;

  ! Loop over all nodes (joints) BASED ON THE APPLICATION NODE NUMBERING to obtain number of non-zero entries of each equaion related to each Degree of freedom of a node. We Name this node "Main Node". 
  Do IJ = 1, NJ ;  ! Based on application numbering

    ! - Obtaning all elements having Main Node ----
    ElementCounter = 0_Smll ;
      Do IEL = 1, NEL ;                                         ! Loop over all elements to find out all elements connected to the Main Node.
          Do INode = 1, NNode ;                                 ! Loop over all nodes of each element to find out if the Main Node belongs to this element.
            If ( INod ( INode, IEL ) == IJ ) Then ;
!              check = .True. ;
              ElementCounter = ElementCounter + 1_Smll ;        ! Counting Number of elements having the main node.
              El_Neighbors ( ElementCounter ) = IEL ;           ! Saving the element number having the Main node
              Exit ;
            End If ;
          End Do ;
      End Do ;

    ! - Finding all node nieghbors of the main node ---
    ! Node_Connect_Total is a wild-sized guess of all nodes that might be in the nieghbor of the main node. Maximum possible nodes connected with the main node is equal to the NNode times the total number of elements that have the main node.
    Allocate ( Node_Connect_Total ( ElementCounter * NNode ),     STAT = ERR_Alloc) ;
      If ( ERR_Alloc /= 0 ) Then ;
        Write (*, Fmt_ALLCT) ERR_Alloc ;  Write (UnInf, Fmt_ALLCT) ERR_Alloc ;
        Write(*, Fmt_FL) ;  Write(UnInf, Fmt_FL) ; Read(*, Fmt_End) ;  Stop ;
      End If ;

    Node_Connect_Total = 0_Lng ;

      ! Storing all neighbor nodes of the main node. At the end of the day, we have collected all node numbers connected to the main node. Since we have saved node numbers based on the elements, there might be zeros in between. So, in the next loop we will eliminate zeros and gather only node numbers.
      Do IEL_Neigh = 1, ElementCounter ;

        Do INode = 1, NNode ;
          Node = INod ( INode, El_Neighbors ( IEL_Neigh ) ) ;  ! Node Number is based on the application node numbering.
          Node = Global_PETSc_Num ( Node ) ;                   ! Converting Application node numbering to PETSc node numbering
          If ( Any ( Node_Connect_Total == Node ) ) Cycle ;    ! We might save this node number before from a previous element.
          Node_Connect_Total ( ( IEL_Neigh - 1_Lng ) * NNode + INode ) = Node ;
        End Do ;
      End Do ;

    ! Total number of Neighbor nodes of the main node
    NNeighbors = Count ( Node_Connect_Total /= 0_Lng ) ;       ! This routine sorts the numbers in the array from large to small. Thus, zero entries go to the end.

    ! Node_Connect holds all nodes connected with the main node.
    Allocate ( Node_Connect ( NNeighbors ),     STAT = ERR_Alloc) ;
      If ( ERR_Alloc /= 0 ) Then ;
        Write (*, Fmt_ALLCT) ERR_Alloc ;  Write (UnInf, Fmt_ALLCT) ERR_Alloc ;
        Write(*, Fmt_FL) ;  Write(UnInf, Fmt_FL) ; Read(*, Fmt_End) ;  Stop ;
      End If ;

    ! Storing all neighbor node numbers from low to high
    Do I = NNeighbors, 1, -1 ;
      Node_Connect ( I ) = MaxVal ( Node_Connect_Total ) ;
      Node_Connect_Total ( MaxLoc ( Node_Connect_Total ) ) = 0_Lng ;
    End Do ;

    DeAllocate( Node_Connect_Total,    STAT = ERR_DeAlloc ) ; 
      IF ( ERR_DeAlloc /= 0 ) Then ;
        Write (*, Fmt_DEALLCT) ERR_DeAlloc ;  Write (UnInf, Fmt_DEALLCT) ERR_DeAlloc ;
        !#Call BEEP_FAIL ;
        Write(*, Fmt_FL) ;  Write(UnInf, Fmt_FL) ; Write(*, Fmt_End) ; Read(*,*) ;  STOP ;
      End If ;

    ! Definig matrices to store non-zero entries of all degrees of freedom of the main node
    Allocate ( NNZ_Node_Stiff ( NDOF, NDOF * NNeighbors ), NNZ_Node_Mass ( NDOF, NDOF * NNeighbors ), NNZ_Node_Damp ( NDOF, NDOF * NNeighbors ),     STAT = ERR_Alloc) ;  ! 
      If ( ERR_Alloc /= 0 ) Then ;
        Write (*, Fmt_ALLCT) ERR_Alloc ;  Write (UnInf, Fmt_ALLCT) ERR_Alloc ;
        Write(*, Fmt_FL) ;  Write(UnInf, Fmt_FL) ; Read(*, Fmt_End) ;  Stop ;
      End If ;

    ! Initialize matrices for each node.
    NNZ_Node_Stiff = 0_Tiny ;
    NNZ_Node_Mass  = 0_Tiny ;
    NNZ_Node_Damp  = 0_Tiny ;

      ! This is an approximation of number of non-zero entries of PETSc objects based on schematic form of matrices. !!!!!??
      Do IEL_Neigh = 1, ElementCounter ;

        IEL = El_Neighbors ( IEL_Neigh ) ; ! Element number 

          ! Update elemental matrices
          If      ( ELGR( IEL ) == 0 ) Then ;
          ! Elmenetal stifness, mass and damping matrices for regular domain

            If ( Matrix_Type == 0_Smll ) then ;

              ! Approximate estimation
              KE = KE_RD ;
              CE = CE_RD ;
              ME = ME_RD ;

            Else If ( Matrix_Type == 1_Smll ) Then ;
              ! Exact estimation

              ! Intialize the stiffness matrix
              KE_Exact = 0.0_DBL ;  ! Elemental stiffness matrix
              CE_Exact = 0.0_DBL ;  ! Elemental damping matrix
              ME_Exact = 0.0_DBL ;  ! Elemental mass matrix

              ! Material properties
              LType  = LTEL ( IEL ) ;        ! Group of Load Type of Element
              MType  = MTEL ( IEL ) ;        ! Group of Material Type of Element

              Lambda = PMat ( MType, 1 ) ;   ! Lame coe.
              MU     = PMat ( MType, 2 ) ;   ! Lame coe.
              Rho    = PMat ( MType, 3 ) ;   ! Density
              ALFA   = PMat ( MType, 4 ) ;   ! Reyleigh damping coe.
              BETTA  = PMat ( MType, 5 ) ;   ! Reyleigh damping coe.

              Call MassDampStiffSLD_2D_8N_Pattern ( IEl, NNode, NDim, NInt,     Rho, Lambda, MU, BETTA, ALFA,     KE_Exact, ME_Exact, CE_Exact,     GAUSS_PNT ) ;

              ! KE/CE/ME_Exact has real values, since only zero or nonzero patterns matters in this stage we transform real value to integer values to save memory.
              Where ( KE_Exact /= 0.0_Dbl ) KE_Exact = 1.0_Dbl ;
              Where ( CE_Exact /= 0.0_Dbl ) CE_Exact = 1.0_Dbl ;
              Where ( ME_Exact /= 0.0_Dbl ) ME_Exact = 1.0_Dbl ;

              KE = Int ( KE_Exact ) ;
              CE = Int ( CE_Exact ) ;
              ME = Int ( ME_Exact ) ;

            End If ; 


          Else If ( ELGR( IEL ) == 1 ) Then ;
          ! Elmenetal stifness, mass and damping matrices for PML domain

            If ( Matrix_Type == 0_Smll ) then ;
              ! Approximate estimation






              KE = KE_PML ;
              CE = CE_PML ;
              ME = ME_PML ;
            Else If ( Matrix_Type == 1_Smll ) Then ;  ! Exact
              ! Exact estimation

              ! Intialize the stiffness matrix
              KE_Exact = 0.0_DBL ;  ! Elemental stiffness matrix
              CE_Exact = 0.0_DBL ;  ! Elemental mass matrix
              ME_Exact = 0.0_DBL ;  ! Elemental mass matrix

              ! Material properties
              LType  = LTEL ( IEL ) ;               ! Group of Load Type of Element
              MType  = MTEL ( IEL ) ;               ! Group of Material Type of Element

              Lambda        = PMat ( MType, 1 ) ;   ! Lame coe.
              MU            = PMat ( MType, 2 ) ;   ! Lame coe.
              Rho           = PMat ( MType, 3 ) ;   ! Density
              R_PML_ALFA_0  = PMat ( MType, 4 ) ;   ! Power of decay or alfa_0
              C_REF_BETA_0  = PMat ( MType, 5 ) ;   ! P-wave velocity or beta_0
              M_PML         = PMat ( MType, 6 ) ;   ! Characteristic length of the domain
              B_PML         = PMat ( MType, 7 ) ;   ! User tunable reflection coefficient

!              Call MassDampStiffPML_2D_8N ( NDim, NNode, NInt, PARAM_Type, IEL,     Lambda, MU, Rho, R_PML_ALFA_0, C_REF_BETA_0, M_PML, B_PML,     INOD,     XYZ,     KE, ME, CE, PML_DIM,     GAUSS_PNT ) ;
              Call MassDampStiffPML_2D_8N_Pattern ( NDim, NNode, NInt, PARAM_Type,      IEL,     Lambda, MU, Rho, R_PML_ALFA_0, C_REF_BETA_0, M_PML, B_PML,     INOD,     XYZ,     KE_Exact, ME_Exact, CE_Exact, PML_DIM,     GAUSS_PNT ) ;


              ! KE/CE/ME_Exact has real values, since only zero or nonzero patterns matters in this stage we transform real value to integer values to save memory.
              Where ( KE_Exact /= 0.0_Dbl ) KE_Exact = 1.0_Dbl ;
              Where ( CE_Exact /= 0.0_Dbl ) CE_Exact = 1.0_Dbl ;
              Where ( ME_Exact /= 0.0_Dbl ) ME_Exact = 1.0_Dbl ;

              KE = Int ( KE_Exact ) ;
              CE = Int ( CE_Exact ) ;
              ME = Int ( ME_Exact ) ;

            End If ; 

          End If ;

          ! Find local node number of the main node in this specific element. Notice that IJ is based on the application numbering, thus 
          Do INode = 1, NNode ;
            If ( INod ( INode, El_Neighbors ( IEL_Neigh ) ) == IJ ) Then ;
              INode_Main = INode ;  ! Local node number of the main node in the current element.
              Exit ;
            End If ;
          End Do ;

          ! Finds the location of the !???????????????????????????????????????????????????????
          Do INode = 1, NNode ; ! Loop on nodes connected to the main node
            Node = INod ( INode, El_Neighbors ( IEL_Neigh ) ) ; ! Based on application node numbering
            Node = Global_PETSc_Num ( Node ) ;                  ! Converting Application node numbering to PETSc node numbering
              Do I = 1, NNeighbors ;
                If ( Node_Connect ( I ) == Node ) Then ;
                  Location = I ;
                  Exit ;
                End If ;
              End Do ;


!              Do IDOF = 1, NDOF ; ! Loop on DOFs of the main node
!                Do I = 1, NDOF ; ! Loop on DOFs of the connected nodes
!                  NNZ_Node_Stiff ( IDOF, ( Location - 1_Smll ) * NDOF + I ) = KE ( ( INode_Main - 1_Smll ) + IDOF, ( INode - 1_Smll ) * NDOF + I ) ;
!                  NNZ_Node_Damp  ( IDOF, ( Location - 1_Smll ) * NDOF + I ) = CE ( ( INode_Main - 1_Smll ) + IDOF, ( INode - 1_Smll ) * NDOF + I ) ;
!                  NNZ_Node_Mass  ( IDOF, ( Location - 1_Smll ) * NDOF + I ) = ME ( ( INode_Main - 1_Smll ) + IDOF, ( INode - 1_Smll ) * NDOF + I ) ;
!                End Do ;
!              End Do ;

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!              Do IDOF = 1, NDOF ;  ! Loop on DOFs of the main node
!                Do I = 1, NDOF ;  ! Loop on DOFs of the other nodes of the element
!                  If ( ID_PETSc ( Node, I ) /= 0_Lng ) Then ;
!                    NNZ_Node_Stiff ( IDOF, ( Location - 1_Smll ) * NDOF + I ) = NNZ_Node_Stiff ( IDOF, ( Location - 1_Smll ) * NDOF + I ) + KE ( ( IDOF - 1_Smll ) * NNode + INode_Main, ( I - 1_Smll ) * NNode + INode ) ;
!                    NNZ_Node_Damp  ( IDOF, ( Location - 1_Smll ) * NDOF + I ) = NNZ_Node_Damp  ( IDOF, ( Location - 1_Smll ) * NDOF + I ) + CE ( ( IDOF - 1_Smll ) * NNode + INode_Main, ( I - 1_Smll ) * NNode + INode ) ;
!                    NNZ_Node_Mass  ( IDOF, ( Location - 1_Smll ) * NDOF + I ) = NNZ_Node_Mass  ( IDOF, ( Location - 1_Smll ) * NDOF + I ) + ME ( ( IDOF - 1_Smll ) * NNode + INode_Main, ( I - 1_Smll ) * NNode + INode ) ;
!                  End If ;
!                End Do ;
!              End Do ;

              ForAll ( IDOF = 1:NDOF, I = 1:NDOF ) ; ! Loop on DOFs of the main node -  Loop on DOFs of the other nodes of the element.
                NNZ_Node_Stiff ( IDOF, ( Location - 1_Smll ) * NDOF + I ) = NNZ_Node_Stiff ( IDOF, ( Location - 1_Smll ) * NDOF + I ) + KE ( ( IDOF - 1_Smll ) * NNode + INode_Main, ( I - 1_Smll ) * NNode + INode ) ;
                NNZ_Node_Damp  ( IDOF, ( Location - 1_Smll ) * NDOF + I ) = NNZ_Node_Damp  ( IDOF, ( Location - 1_Smll ) * NDOF + I ) + CE ( ( IDOF - 1_Smll ) * NNode + INode_Main, ( I - 1_Smll ) * NNode + INode ) ;
                NNZ_Node_Mass  ( IDOF, ( Location - 1_Smll ) * NDOF + I ) = NNZ_Node_Mass  ( IDOF, ( Location - 1_Smll ) * NDOF + I ) + ME ( ( IDOF - 1_Smll ) * NNode + INode_Main, ( I - 1_Smll ) * NNode + INode ) ;
              End ForAll ;
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


            !Where ( NNZ_Node_Stiff /= 0_Tiny ) NNZ_Node_Stiff = 1_Tiny ;       ! UNcomment this later
            !Where ( NNZ_Node_Damp  /= 0_Tiny ) NNZ_Node_Damp  = 1_Tiny ;
            !Where ( NNZ_Node_Mass  /= 0_Tiny ) NNZ_Node_Mass  = 1_Tiny ;

          End Do ;
      End Do ;

    ! Find Number of NonZero entries for PETSc objects.
      Do IDOF = 1, NDOF ;
        NEq = ID_PETSc ( Global_PETSc_Num ( IJ ), IDOF ) ;  ! Equation number of the !???? in the global matrices based on PETSc node numbering
        If ( NEq == 0_Lng ) Cycle ;

        ! Searching for the rank number this equation belongs to. LowerLimit and Upper Limit determines the limits of the diagonal for PETSc objects. Refer to PETSc manual page 59.
          Do IParts = 1_Shrt, NParts ;
            If ( NEq <= NEqRank ( IParts ) ) then ;
                If ( IParts == 1_Shrt ) then ;
                  LowerEq = 0_Lng ;
                Else ;
                  LowerEq = NEqRank ( IParts - 1_Shrt ) + 1_Lng ;   ! Lower equation number of the diagonal on this rank  !?? clarify the Name
                End if ;
              UpperEq = NEqRank ( IParts ) ;   ! Upper equation number of the diagonal on this rank 
              Exit ;
            End If ;
          End Do ;


          ! Non Zero entris
          Do I = 1, NNeighbors ;
            Node = Node_Connect ( I ) ;
              Do J = 1, NDOF ;
                !K = ID_PETSc ( Global_PETSc_Num ( Node ), J ) ;
                K = ID_PETSc ( Node , J ) ;

                If ( K == 0_Smll ) Cycle ;

                  If ( LowerEq <= K .AND. K <= UpperEq ) then ; ! Non zero entries of diagonal part of the objects
                    If ( NNZ_Node_Stiff ( IDOF , ( I - 1_Smll ) * NDOF + J ) /= 0_Tiny )  D_NNZ_Stiff ( NEq ) = D_NNZ_Stiff ( NEq ) + 1_Shrt ;
                    If ( NNZ_Node_Damp  ( IDOF , ( I - 1_Smll ) * NDOF + J ) /= 0_Tiny )  D_NNZ_Damp  ( NEq ) = D_NNZ_Damp  ( NEq ) + 1_Shrt ;
                    If ( NNZ_Node_Mass  ( IDOF , ( I - 1_Smll ) * NDOF + J ) /= 0_Tiny )  D_NNZ_Mass  ( NEq ) = D_NNZ_Mass  ( NEq ) + 1_Shrt ;
                  Else ; ! Non zero entries of off-diagonal part of the objects
                    If ( NNZ_Node_Stiff ( IDOF , ( I - 1_Smll ) * NDOF + J ) /= 0_Tiny )  O_NNZ_Stiff ( NEq ) = O_NNZ_Stiff ( NEq ) + 1_Shrt ;
                    If ( NNZ_Node_Damp  ( IDOF , ( I - 1_Smll ) * NDOF + J ) /= 0_Tiny )  O_NNZ_Damp  ( NEq ) = O_NNZ_Damp  ( NEq ) + 1_Shrt ;
                    If ( NNZ_Node_Mass  ( IDOF , ( I - 1_Smll ) * NDOF + J ) /= 0_Tiny )  O_NNZ_Mass  ( NEq ) = O_NNZ_Mass  ( NEq ) + 1_Shrt ;
                  End If ;
              End Do ;
          End Do ;

      End Do ;

    ! Deallocating local arrays used in the loop
    Deallocate ( Node_Connect, NNZ_Node_Stiff, NNZ_Node_Mass, NNZ_Node_Damp,     STAT = ERR_Alloc) ; 
      If ( ERR_Alloc /= 0 ) Then ;
        Write (*, Fmt_ALLCT) ERR_Alloc ;  Write (UnInf, Fmt_ALLCT) ERR_Alloc ;
        Write(*, Fmt_FL) ;  Write(UnInf, Fmt_FL) ; Read(*, Fmt_End) ;  Stop ;
      End If ;


  End Do ;


! Deallocating all local arrays
    Deallocate ( El_Neighbors,     STAT = ERR_Alloc) ; 
      If ( ERR_Alloc /= 0 ) Then ;
        Write (*, Fmt_ALLCT) ERR_Alloc ;  Write (UnInf, Fmt_ALLCT) ERR_Alloc ;
        Write(*, Fmt_FL) ;  Write(UnInf, Fmt_FL) ; Read(*, Fmt_End) ;  Stop ;
      End If ;

Write(*     ,*) 'End Subroutine < Object_Size >' ;
!Write(UnInf,*) 'End Subroutine < Object_Size >' ;
Return ;
End Subroutine Object_Size ;


End Module PETScObject_Size ;
