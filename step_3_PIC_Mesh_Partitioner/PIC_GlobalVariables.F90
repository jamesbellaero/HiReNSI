! - Integer Variables -------------------------------------------------------------------------------------------------------------------------------

! Tiny Integers
Integer (Kind=Tiny)  :: Matrix_Type ;            ! Accuracy of calculations of Number of non-zero entries of the petsc objects. 0: approximate. 1: exact
Integer (Kind=Tiny)  :: MetisType ;              ! Metis Type 0: Metis. 1: ParMetis
Integer (Kind=Tiny)  :: Output_Type ;            ! Determines if the output files are formatted -0- or binary -1-

! Smll Integers
Integer (Kind=Smll)  :: NNeighbor ;              ! Maximum Number of Ranks that a node might belongs to - By default it is set to be 3 - if necessary, will be increased automatically
Integer (Kind=Shrt)  :: ETypeG ;                 ! Element type. Used only in Metis routines and not in the ParMetis routines

! Shrt Integers
Integer (Kind=Shrt)  :: NumFlag ;                ! NumFlag Used to indicate which numbering scheme is Used for the element node array. NumFlag can take the following two values:  0 C-style numbering is assumed that starts from 0 - 1 Fortran-style numbering is assumed that starts from 1
Integer (Kind=Shrt)  :: NParts ;                 ! Number of desired partitions. This is equal to the number of ranks in the main code.
Integer (Kind=Shrt)  :: EdgeCut ;                ! Dual Edgecut: stores the number of edges that are cut by the partition in the dual graph, upon successful completion
Integer (Kind=Shrt)  :: ObjVal ;                 ! Object value used in METIS_PartMeshDual in METIS v 5.1.0. see page 28 of the manual
Integer (Kind=Shrt)  :: WgtFlag ;                ! WeightFlag This is used to indicate if the elements of the mesh have weights associated with them. The wgtflag can take two values:
                                                 ! 0 No weights (elmwgt is NULL).
                                                 ! 2 Weights on the vertices only.
Integer (Kind=Shrt)  :: ncon = 1 ;               ! This is used to specify the number of weights that each vertex has. It is also the number of balance constraints that must be satisfied
Integer (Kind=Shrt)  :: ncommonnodes ;           ! see page 21 of ParMetis V3.2
Integer (Kind=Shrt)  :: Kway ;                   ! see page 21 of ParMetis V3.2


! Lng Integers
Integer (Kind=Lng )  :: NEQMTotal ;              ! Modified Total Number of Equatoins
Integer (Kind=Lng )  :: Neind ;                  ! 
Integer (Kind=Shrt)  :: NJG ;                    ! the number of joints of the graph, i.e., number of nodes of the model after Converting the model to first order elements

! - Real Variables ----------------------------------------------------------------------------------------------------------------------------------
Real (Kind=Dbl )     :: tpwgts2 ;

! - Logical Variable --------------------------------------------------------------------------------------------------------------------------------

! - Integer Arrays ----------------------------------------------------------------------------------------------------------------------------------
Integer (Kind=Shrt)  :: Poptions (0:3) ;
!Integer (Kind=Shrt)  :: MOptions ( METIS_NOPTIONS ) ;

! - Real Arrays -------------------------------------------------------------------------------------------------------------------------------------

! - Integer ARRAY ALLOCATION ------------------------------------------------------------------------------------------------------------------------


! Tiny Integers

! Smll Integers

! Shrt Integers
Integer (Kind=Shrt), Allocatable, Dimension(:)   :: ELMNTS ;         ! Element Connectivities in Metis
Integer (Kind=Shrt), Allocatable, Dimension(:)   :: NPart ;          ! Holds Node Partitioning, we do not use this array
Integer (Kind=Shrt), Allocatable, Dimension(:)   :: EPart ;          ! Holds Element Partitioning
Integer (Kind=Shrt), Allocatable, Dimension(:)   :: ElmDist ;        ! Element Distributions. see page 21 of ParMetis V3.2
Integer (Kind=Shrt), Allocatable, Dimension(:)   :: eptr ;           ! holds the starting point of node numbers in eind vector. size (nel +1) see page 21 of ParMetis V3.2
Integer (Kind=Shrt), Allocatable, Dimension(:)   :: eind ;           ! Connectivity. Holds nodes numbers of elemenets in a vector of size (nnode*nel)see page 21 of ParMetis V3.2
Integer (Kind=Shrt), Allocatable, Dimension(:)   :: eind2 ;          ! Connectivity. Holds nodes numbers of elemenets in a vector of size (nnode*nel)see page 21 of ParMetis V3.2
Integer (Kind=Shrt), Allocatable, Dimension(:)   :: ElmWgt ;         ! Element Weight. size (nel). Equal to NULL if unweighted.


Integer (Kind=Shrt), Allocatable, Dimension(:)   :: VWgt ;           ! Element Weight. size (nel). Equal to NULL if unweighted.
Integer (Kind=Shrt), Allocatable, Dimension(:)   :: VSize ;          ! Size of the elements 
Integer (Kind=Shrt), Allocatable, Dimension(:)   :: Moptions ;       ! Options for Mesh partitioning


Integer (Kind=Shrt), Allocatable, Dimension(:  ) :: D_NNZ_Stiff ;    ! Number of Non-Zero entries of Diagonal part of the stiffness matrix for each node
Integer (Kind=Shrt), Allocatable, Dimension(:  ) :: O_NNZ_Stiff ;    ! Number of Non-Zero entries of Off-Diagonal part of the stiffness matrix for each node
Integer (Kind=Shrt), Allocatable, Dimension(:  ) :: D_NNZ_Mass ;     ! Number of Non-Zero entries of Diagonal part of the mass matrix for each node
Integer (Kind=Shrt), Allocatable, Dimension(:  ) :: O_NNZ_Mass ;     ! Number of Non-Zero entries of Off-Diagonal part of the mass matrix for each node
Integer (Kind=Shrt), Allocatable, Dimension(:  ) :: D_NNZ_Damp ;     ! Number of Non-Zero entries of Diagonal part of the damp matrix for each node
Integer (Kind=Shrt), Allocatable, Dimension(:  ) :: O_NNZ_Damp ;     ! Number of Non-Zero entries of Off-Diagonal part of the damp matrix for each node

! Lng Integers
Integer (Kind=Lng ), Allocatable, Dimension(:)   :: NEqRank ;        ! Number of local rows on each rank that forms the global matrices.
Integer (Kind=Lng ), Allocatable, Dimension(:)   :: NNodeRank ;      ! Number of nodes belongs to each rank in terms of memory based on PETSc numbering. So, if NNodeRank (0)=x and NNodeRank (1)=y then entries of stiffness, ... matrices related to the first x nodes will be saved on rank 0 and the next (y-x) nodes will be saved on rank 1 and so on.
Integer (Kind=Lng ), Allocatable, Dimension(:  ) :: NEL_Rank ;       ! Number of elements in each partition ( rank )
Integer (Kind=Lng ), Allocatable, Dimension(:  ) :: NJ_Rank ;        ! Number of joints (nodes) in each partition ( rank )
Integer (Kind=Lng ), Allocatable, Dimension(:  ) :: Global_PETSc_Num;! Stores the global PETSc numbering of each node in its relevant application ordering
Integer (Kind=Lng ), Allocatable, Dimension(:,:) :: Local_PETSc_Num ;! Holds the local node numbers of each rank. Node numbering is based on Global application numbering - ( NJ, NParts )
Integer (Kind=Lng ), Allocatable, Dimension(:,:) :: ID_Application ; ! Holds global equation numbers based on the application node numbering.
Integer (Kind=Lng ), Allocatable, Dimension(:,:) :: ID_PETSc ;       ! Holds global equation numbers based on the PETSc numbering.


! ----------------- added for inversion DS ------
Integer (Kind=Smll), Allocatable, Dimension(:)   :: Node_Mat_ID;     ! 0 = RD , 1 = I , 2 = PML
Integer (Kind=Lng),  Allocatable, Dimension(:)   :: Node_Mat_Mapping ! maps a (dependent) node to an independent node
! -----------------------------------------------


! - Real ARRAY ALLOCATION ---------------------------------------------------------------------------------------------------------------------------
Real (Kind=DBL)   , Allocatable, Dimension(:)    :: tpwgts ;         ! see page 21 of ParMetis V3.2 - size (ncon * NParts)
Real (Kind=DBL)   , Allocatable, Dimension(:)    :: ubvec ;          ! see page 21 of ParMetis V3.2 - size (ncon)


! - Type ALLOCATION ---------------------------------------------------------------------------------------------------------------------------------
Type ( NodeID ) :: Nodes ;                                           ! Node%Locs :  Holds rank numbers of which this node belongs to.
                                                                     ! Node%Rep  : Holds the number of Repeatations on the ranks (How many times this node appears on ranks), Useful for neighboring
! - Character Variables -----------------------------------------------------------------------------------------------------------------------------


! - Pointers  ---------------------------------------------------------------------------------------------------------------------------------------
!Integer, Pointer :: PVWgt (:), PVSize (:), PMoptions (:)
!Real,    Pointer :: Ptpwgts (:)
Integer, Pointer :: PVWgt=>null(), PVSize=>null(), PMoptions=>null()
Real,    Pointer :: Ptpwgts=>null()
