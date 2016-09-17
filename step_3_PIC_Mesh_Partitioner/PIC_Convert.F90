
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                                                                                                                  ++
! Start:        01 July 2011                                                                                                                       ++
! Last Update:  23 August 2012                                                                                                                     ++
!                                                                                                                                                  ++
! Description: THIS Module Converts the second order elements to first order elements. METIS has only the first order elements.                    ++
!                                                                                                                                                  **
! Developed by: Babak Poursartip                                                                                                                   **
!                                                                                                                                                  ++
!              Subroutines                     Description                                                                                         ++
!              ==============                  ==============                                                                                      ++
!              Convert                                                                                                                             ++
!                                                                                                                                                  ++
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Module ConvertToMetis;

Use Parameters ;         ! Module comprising all Gauss Points, Type declarations, Constants, Definitions, ...

Implicit None ;

  Interface Convert
    Module Procedure MetisConvert_V4, MetisConvert_V51, ParConvert ;
  End Interface    ;

Contains ;

! ***************************************************************************************************************************************************
!                                                                                                                                                  **
! Start:        01 July 2011                                                                                                                       **
! Last Update:  31 May 2012                                                                                                                        **
! Description: THIS Module Converts the second order elements to first order elements. METIS has only the first order elements.                    **
!                                                                                                                                                  **
! Developed by: Babak Poursartip                                                                                                                   **
!                                                                                                                                                  **
! ***************************************************************************************************************************************************

Subroutine MetisConvert_V4 (                                                                                                 &
NNode, NEL, NJ, ETypeG, NJG,                                                                                                 & ! Integer Variables
!                                                                                                                            & ! Real Variables
ELT, ELMNTS,     INOD                                                                                                        & ! Integer Arrays
!                                                                                                                            & ! Real Arrays
!                                                                                                                            & ! Characters
!                                                                                                                            & ! Type 
) ;

Implicit None ;

! =========================== Global Variables ======================================================================================================
! - Integer Variables -------------------------------------------------------------------------------------------------------------------------------
Integer (Kind=Tiny), Intent(In)    :: NNode ;

Integer (Kind=Lng ), Intent(In)    :: NEL, NJ ;

Integer (Kind=Shrt), Intent(Out)   :: ETypeG ;
Integer (Kind=Shrt), Intent(OUT)   :: NJG ;

! - Real Variables ----------------------------------------------------------------------------------------------------------------------------------
! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
! - Integer Arrays ----------------------------------------------------------------------------------------------------------------------------------
Integer (Kind=Smll), Intent(IN), Dimension (:  )  :: ELT ;
Integer (Kind=Lng ), Intent(IN), Dimension (:,:)  :: INOD ;

Integer (Kind=Shrt), Intent(OUT), Allocatable, Dimension (:  )  :: ELMNTS ;

! - Real Arrays -------------------------------------------------------------------------------------------------------------------------------------
! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
! - Character Variables -----------------------------------------------------------------------------------------------------------------------------
! - Logical Variables -------------------------------------------------------------------------------------------------------------------------------
! - Type DECLERATIONS -------------------------------------------------------------------------------------------------------------------------------

! =========================== LOCAL Variables =======================================================================================================
! - Integer Variables -------------------------------------------------------------------------------------------------------------------------------
Integer (Kind=Tiny)  :: NNodeG ;                 ! Number of Nodes of Metis elements
Integer (Kind=Smll)  :: EType ;                  ! Element Type based on my code numbering
Integer (Kind=Shrt)  :: IEL ;                    ! Loop Index on the NEL - number of elements.
Integer (Kind=Shrt)  :: INode ;                  ! Loop index on the NNode - Number of nodes.
Integer (Kind=Shrt)  :: Counter ;                ! Counter.

! - Real Variables ----------------------------------------------------------------------------------------------------------------------------------
! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
! - Integer Arrays ----------------------------------------------------------------------------------------------------------------------------------
Integer (Kind=Shrt), Dimension ( NJ ) :: NodeReOrdering ;  ! Holds the reordered node numbers of elements. For second order elements, midnodes are eliminated and the rest are saved in this array.


! - Real Arrays -------------------------------------------------------------------------------------------------------------------------------------
! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
! - Character Variables -----------------------------------------------------------------------------------------------------------------------------
! - Logical Variables -------------------------------------------------------------------------------------------------------------------------------

! =========================== Subroutine CODE =======================================================================================================

! NOTICE : 1- Metis has only first order elements, thus in order to partition a model of second oreder elements we have to Convert them to first order ones.
!          2- This subroutine is for uni-type-element models, i.e., ELT (:) = is constant

Print *, "MetisConvert_V4:: Prepare model for partitioning"

EType = ELT(1) ;

! Assigning NNodeG, Number of nodes of Graph of ParMetis element and ETypeG, Element type of ParMetis. ParMetis has only first order elements AND filling element connectivities for Metis
  If     ( EType == El2d4NSldPS .OR. EType == El2d4NSldPN .OR. EType == El2d4NPMLPS .OR. EType == El2d4NPMLPN ) Then ; ! First order 2d 4-node elements
    NNodeG = 4 ;
    ETypeG = 4 ;

    Allocate ( ELMNTS ( NEl * NNodeG ) ) ;

      DO IEL = 1, NEL ;
        ELMNTS ( ( IEL - 1 ) * NNodeG + 1 ) = INOD ( 1, IEL ) ;
        ELMNTS ( ( IEL - 1 ) * NNodeG + 2 ) = INOD ( 4, IEL ) ;
        ELMNTS ( ( IEL - 1 ) * NNodeG + 3 ) = INOD ( 3, IEL ) ;
        ELMNTS ( ( IEL - 1 ) * NNodeG + 4 ) = INOD ( 2, IEL ) ;
      End Do ;

    NJG = NJ ;

  ElseIF ( EType == El2d8NSldPS .OR. EType == El2d8NSldPN .OR. EType == El2d8NPMLPS  .OR. EType == El2d8NPMLPN .OR. EType == El2d9NSldPS .OR. EType == El2d9NSldPN .OR. EType == El2d9NPMLPS .OR. EType == El2d9NPMLPN    .OR. EType == SEl2d9NSldPN .OR. EType == SEl2d9NPMLPN .OR. EType == SEl2d9NMPMLPN ) Then ;  ! Second order 2d 8-node elements OR Second order 2d 9-node
    NNodeG = 4 ;
    ETypeG = 4 ;

    Allocate ( ELMNTS ( NEl * NNodeG ) ) ;

    ! Converting second order elements to first order by reordering the node numbers 
    NodeReOrdering = 1_Shrt ;
    ForAll ( IEL = 1:NEL, INode = NNodeG + 1:NNode )  NodeReOrdering ( INOD ( INode, IEL ) ) = 0_Shrt ; ! Zeroing midside nodes 

    ! renumbering
    Counter = 0 ;
      DO INode = 1, NJ ;
        IF ( NodeReOrdering ( INode ) == 1 ) Then ;
          Counter = Counter + 1 ;
          NodeReOrdering ( INode ) = Counter ;
        End If ;
      End Do ;

    NJG = Counter ;

      DO IEL = 1, NEL ;
        ELMNTS ( ( IEL - 1 ) * NNodeG + 1 ) = NodeReOrdering ( INOD ( 1, IEL ) ) ;
        ELMNTS ( ( IEL - 1 ) * NNodeG + 2 ) = NodeReOrdering ( INOD ( 4, IEL ) ) ;
        ELMNTS ( ( IEL - 1 ) * NNodeG + 3 ) = NodeReOrdering ( INOD ( 3, IEL ) ) ;
        ELMNTS ( ( IEL - 1 ) * NNodeG + 4 ) = NodeReOrdering ( INOD ( 2, IEL ) ) ;
      End Do ;

  ElseIF ( EType == El3d8NSld .OR. EType == El2d9NPMLPS  .OR. EType == SEl3d8NSld .OR. EType == SEl3d8NMPML .OR. EType == El3d8NFld) Then ; ! First order 3d 8-node elements
    NNodeG = 8 ;
    ETypeG = 3 ;

    Allocate ( ELMNTS ( NEl * NNodeG ) ) ;

      DO IEL = 1, NEL ; !?? see if node ordering is correct 
        ELMNTS ( ( IEL - 1 ) * NNodeG + 1 ) = INOD ( 1, IEL ) ;
        ELMNTS ( ( IEL - 1 ) * NNodeG + 2 ) = INOD ( 4, IEL ) ;
        ELMNTS ( ( IEL - 1 ) * NNodeG + 3 ) = INOD ( 3, IEL ) ;
        ELMNTS ( ( IEL - 1 ) * NNodeG + 4 ) = INOD ( 2, IEL ) ;
        ELMNTS ( ( IEL - 1 ) * NNodeG + 5 ) = INOD ( 5, IEL ) ;
        ELMNTS ( ( IEL - 1 ) * NNodeG + 6 ) = INOD ( 8, IEL ) ;
        ELMNTS ( ( IEL - 1 ) * NNodeG + 7 ) = INOD ( 7, IEL ) ;
        ELMNTS ( ( IEL - 1 ) * NNodeG + 8 ) = INOD ( 6, IEL ) ;
      End Do ;

    NJG = NJ ;

  ElseIF ( EType == El3d20NSld .OR. EType == El3d20NPML .OR. EType == El3d27NSld .OR. EType == El3d27NPML .OR. EType == SEl3d27NMPML .OR. EType == SEl3d27NSld  ) Then ; ! Second order 3d 20-node elements OR Second order 3d 27-node elements OR Spectral elenement 27-node
    NNodeG = 8 ;
    ETypeG = 3 ;

    Allocate ( ELMNTS ( NEl * NNodeG ) ) ;

    ! Converting second order elements to first order by reordering the node numbers 
    NodeReOrdering = 1_Shrt ;
    ForAll ( IEL = 1:NEL, INode = NNodeG + 1:NNode )  NodeReOrdering ( INOD ( INode, IEL ) ) = 0_Shrt ; ! Zeroing midside nodes 

    Counter = 0 ;
      DO INode = 1, NJ ;
        IF ( NodeReOrdering ( INode ) == 1 ) Then ;
          Counter = Counter + 1 ;
          NodeReOrdering ( INode ) = Counter ;
        End If ;
      End Do ;

    NJG = Counter ;

      DO IEL = 1, NEL ;
        ELMNTS ( ( IEL - 1 ) * NNodeG + 1 ) = NodeReOrdering ( INOD ( 1, IEL ) ) ;
        ELMNTS ( ( IEL - 1 ) * NNodeG + 2 ) = NodeReOrdering ( INOD ( 4, IEL ) ) ;
        ELMNTS ( ( IEL - 1 ) * NNodeG + 3 ) = NodeReOrdering ( INOD ( 3, IEL ) ) ;
        ELMNTS ( ( IEL - 1 ) * NNodeG + 4 ) = NodeReOrdering ( INOD ( 2, IEL ) ) ;
        ELMNTS ( ( IEL - 1 ) * NNodeG + 5 ) = NodeReOrdering ( INOD ( 5, IEL ) ) ;
        ELMNTS ( ( IEL - 1 ) * NNodeG + 6 ) = NodeReOrdering ( INOD ( 8, IEL ) ) ;
        ELMNTS ( ( IEL - 1 ) * NNodeG + 7 ) = NodeReOrdering ( INOD ( 7, IEL ) ) ;
        ELMNTS ( ( IEL - 1 ) * NNodeG + 8 ) = NodeReOrdering ( INOD ( 6, IEL ) ) ;
      End Do ;

  ElseIF ( EType == El2d3NSldPS  .OR. EType == El2d3NSldPN ) Then ; ! First order 2d 3-node elements
    NNodeG = 3 ;
    ETypeG = 1 ;

    Allocate ( ELMNTS ( NEl * NNodeG ) ) ;

      DO IEL = 1, NEL ; !?? see if node ordering is correct 
        ELMNTS ( ( IEL - 1 ) * NNodeG + 1 ) = INOD ( 1, IEL ) ;
        ELMNTS ( ( IEL - 1 ) * NNodeG + 2 ) = INOD ( 2, IEL ) ;
        ELMNTS ( ( IEL - 1 ) * NNodeG + 3 ) = INOD ( 3, IEL ) ;
      End Do ;

    NJG = NJ ;

  ElseIF ( EType == El2d6NSldPS .OR. EType == El2d6NSldPN .OR. EType == El2d6NPMLPS .OR. EType == El2d6NPMLPN .OR. EType == El2d7NSldPS .OR. EType == El2d7NSldPN .OR. EType == El2d7NPMLPS .OR. EType == El2d7NPMLPN ) Then ; ! Second order 2d 6-node elements OR Second order 2d 7-node elements
    NNodeG = 3 ;
    ETypeG = 1 ;

    Allocate ( ELMNTS ( NEl * NNodeG ) ) ;

    ! Converting second order elements to first order by reordering the node numbers 
    NodeReOrdering = 1_Shrt ;
    ForAll ( IEL = 1:NEL, INode = NNodeG + 1:NNode )  NodeReOrdering ( INOD ( INode, IEL ) ) = 0_Shrt ; ! Zeroing midside nodes 

    Counter = 0 ;
      DO INode = 1, NJ ;
        IF ( NodeReOrdering ( INode ) == 1 ) Then ;
          Counter = Counter + 1 ;
          NodeReOrdering ( INode ) = Counter ;
        End If ;
      End Do ;

    NJG = Counter ;

      DO IEL = 1, NEL ;
        ELMNTS ( ( IEL - 1 ) * NNodeG + 1 ) = NodeReOrdering ( INOD ( 1, IEL ) ) ;
        ELMNTS ( ( IEL - 1 ) * NNodeG + 2 ) = NodeReOrdering ( INOD ( 2, IEL ) ) ;
        ELMNTS ( ( IEL - 1 ) * NNodeG + 3 ) = NodeReOrdering ( INOD ( 3, IEL ) ) ;
      End Do ;

  ElseIF ( EType == El3d4NSld           ) Then ; ! First order 3d 4-node elements
    NNodeG = 4 ;
    ETypeG = 2 ;

    Allocate ( ELMNTS ( NEl * NNodeG ) ) ;

      DO IEL = 1, NEL ; !?? see if node ordering is correct 
        ELMNTS ( ( IEL - 1 ) * NNodeG + 1 ) = INOD ( 1, IEL ) ;
        ELMNTS ( ( IEL - 1 ) * NNodeG + 2 ) = INOD ( 2, IEL ) ;
        ELMNTS ( ( IEL - 1 ) * NNodeG + 3 ) = INOD ( 3, IEL ) ;
        ELMNTS ( ( IEL - 1 ) * NNodeG + 4 ) = INOD ( 4, IEL ) ;
      End Do ;

    NJG = NJ ;

  ElseIF ( EType == El3d10NSld .OR. EType == El3d10NPML .OR. EType == El3d15NSld .OR. EType == El3d15NPML ) Then ; ! Second order 3d 10-node elements OR Second order 3d 15-node elements
    NNodeG = 4 ;
    ETypeG = 2 ;

    Allocate ( ELMNTS ( NEl * NNodeG ) ) ;

    ! Converting second order elements to first order by reordering the node numbers 
    NodeReOrdering = 1_Shrt ;
    ForAll ( IEL = 1:NEL, INode = NNodeG + 1:NNode )  NodeReOrdering ( INOD ( INode, IEL ) ) = 0_Shrt ; ! Zeroing midside nodes 

    Counter = 0 ;
      DO INode = 1, NJ ;
        IF ( NodeReOrdering ( INode ) == 1 ) Then ;
          Counter = Counter + 1 ;
          NodeReOrdering ( INode ) = Counter ;
        End If ;
      End Do ;

    NJG = Counter ;

      DO IEL = 1, NEL ;
        ELMNTS ( ( IEL - 1 ) * NNodeG + 1 ) = NodeReOrdering ( INOD ( 1, IEL ) ) ;
        ELMNTS ( ( IEL - 1 ) * NNodeG + 2 ) = NodeReOrdering ( INOD ( 2, IEL ) ) ;
        ELMNTS ( ( IEL - 1 ) * NNodeG + 3 ) = NodeReOrdering ( INOD ( 3, IEL ) ) ;
        ELMNTS ( ( IEL - 1 ) * NNodeG + 4 ) = NodeReOrdering ( INOD ( 4, IEL ) ) ;
      End Do ;

  End If ;

Write(*    ,*) 'End Subroutine < MetisConvert_V4 >' ;
Write(UnInf,*) 'End Subroutine < MetisConvert_V4 >' ;
Return ;
End Subroutine MetisConvert_V4 ;

! ***************************************************************************************************************************************************
!                                                                                                                                                  **
! Start:        12 July 2012                                                                                                                       **
! Last Update:  12 July 2012                                                                                                                       **
! Description: This subroutine prpares reguired information for the mesh partitioner, ParMetis 3.2.                                                **
! Called by:                                                                                                                                       **
!                                                                                                                                                  **
! Developed by: Babak Poursartip                                                                                                                   **
!                                                                                                                                                  **
! ***************************************************************************************************************************************************

Subroutine ParConvert    (                                                                                                      &
NDim,                                                                                                                           & ! Integer (1) Variables
KWay, NParts,                                                                                                                   & ! Integer (2) Variables
WgtFlag, NumFlag, NCommonNodes,                                                                                                 & ! Integer (4) Variables
NEl, NJ, NJG,                                                                                                                   & ! Integer (8) Variables
!                                                                                                                               & ! Real Variables
ELT,     ElmWgt, ElmDist, eptr, eind, Poptions,     INod,                                                                       & ! Integer Arrays
tpwgts, ubvec                                                                                                                   & ! Real Arrays
!                                                                                                                               & ! Characters
!                                                                                                                               & ! Type
) ;

Implicit None ;

! =========================== Global Variables ======================================================================================================

! - Integer Variables -------------------------------------------------------------------------------------------------------------------------------
Integer (Kind=Tiny), Intent(In)    :: NDim ;
Integer (Kind=Shrt), Intent(In)    :: KWay, NParts ;
Integer (Kind=Shrt), Intent(Out)   :: WgtFlag, NumFlag, NCommonNodes ;
Integer (Kind=Lng ), Intent(In)    :: NEl, NJ ;
Integer (Kind=Shrt ), Intent(Out)   :: NJG ;

! - Real Variables ----------------------------------------------------------------------------------------------------------------------------------
!#Real (Kind=Dbl),     Intent(In)    ::  ;
!#Real (Kind=Dbl),     Intent(InOut) ::  ;
!#Real (Kind=Dbl),     Intent(OUT)   ::  ;
! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
!#Complex,             Intent(In)    ::  ;
!#Complex,             Intent(InOut) ::  ;
!#Complex,             Intent(OUT)   ::  ;
! - Integer Arrays ----------------------------------------------------------------------------------------------------------------------------------
Integer (Kind=Smll), Intent(IN),     Dimension (:)  :: ELT ;
Integer (Kind=Shrt), Intent(Out),    Dimension (:)  :: ElmWgt, ElmDist, eptr, eind, Poptions ;
Integer (Kind=Lng ), Intent(IN),     Dimension (:,:):: INod ;

! - Real Arrays -------------------------------------------------------------------------------------------------------------------------------------
Real (Kind=Dbl),     Intent(OUT),   Dimension (:  )  :: tpwgts, ubvec ;
! - Character Variables -----------------------------------------------------------------------------------------------------------------------------
! - Logical Variables -------------------------------------------------------------------------------------------------------------------------------
! - Type DECLERATIONS -------------------------------------------------------------------------------------------------------------------------------

! - PETSC INTERNAL Variables ------------------------------------------------------------------------------------------------------------------------
!PetscMPIInt    :: Size ;                         ! Total number of ranks
!PetscMPIInt    :: Rank ;                         ! Rank number

! =========================== Local Variables =======================================================================================================
! - Integer Variables -------------------------------------------------------------------------------------------------------------------------------
Integer (Kind=Tiny)  :: NNode ;                  ! Number of Nodes of element in my code
Integer (Kind=Tiny)  :: NNodeG ;                 ! Number of Nodes of element in ParMetis
Integer (Kind=Tiny)  :: INode ;                  ! Loop index over NNode
Integer (Kind=Lng )  :: ETypeG ;                  ! Element Type
Integer (Kind=Lng )  :: NJCounter, Counter ;     ! Counter
Integer (Kind=Lng )  :: IEl ;                    ! Loop Index over NEl
Integer (Kind=Lng )  :: Node ;                   ! Node Number

! - Real Variables ----------------------------------------------------------------------------------------------------------------------------------
! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
! - Integer Arrays ----------------------------------------------------------------------------------------------------------------------------------
Integer (Kind=Shrt), Dimension ( NJ ) :: NodeReOrdering ;       ! Holds the reordered node numbers of elements. For second order elements, midnodes are eliminated and the rest are saved in this array.

! - Real Arrays -------------------------------------------------------------------------------------------------------------------------------------
! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
! - Character Variables -----------------------------------------------------------------------------------------------------------------------------
! - Logical Variables -------------------------------------------------------------------------------------------------------------------------------
! - Type DECLERATIONS -------------------------------------------------------------------------------------------------------------------------------

! =========================== Subroutine CODE =======================================================================================================

ElmDist ( 1  ) = 1 ;
ElmDist (KWay+1) = NEl+1 ;

! Converting second order elements to first order by reordering the node numbers 
NodeReOrdering = 0_Shrt ;
NJCounter = 0 ;
Counter = 0 ;

  Do IEl = 1, NEl ;

    ETypeG = ELT ( IEl ) ;

      ! Call corresponding GlobalMatrices subroutine
      IF      ( ETypeG == El2d4NSldPS ) Then ;  ! 4 noded -2D - Solid - PLANE STRESS
        NNode = 4 ; NNodeG = 4 ;

      Else If ( ETypeG == El2d8NSldPN ) Then ;  ! 8 noded -2D - Solid - PLANE STRAIN
        NNode = 8 ; NNodeG = 4 ;

      Else If ( ETypeG == El2d3NSldPS ) Then ;  ! 3 noded -2D - Solid - PLANE STRAIN
        NNode = 3 ; NNodeG = 3 ;

      Else If ( ETypeG == El2d8NPMLPN ) Then ;  ! 8 noded -2D - PML - PLANE STRAIN
        NNode = 8 ; NNodeG = 4 ;

      Else If ( ETypeG == El2d6NSldPN ) Then ;  ! 6 noded -2D - Solid - PLANE STRAIN
        NNode = 6 ; NNodeG = 3 ;

      Else If ( ETypeG == El2d6NPMLPN ) Then ;  ! 6 noded -2D - PML - PLANE STRAIN
        NNode = 6 ; NNodeG = 3 ;

      Else ; ! Element Type is not in the list
        Write(*, Fmt_Element2) ;  Write(UnInf, Fmt_Element2) ;
        Write(*, Fmt_FL) ; Write(UnInf, Fmt_FL) ;
        !#Call BEEP_FAIL ;
        Write(*, Fmt_End) ; Read(*,*) ;  STOP ;
      End If ;

    eptr ( IEl ) = Counter + 1 ;
      Do INode = 1, NNodeG ;
        Node = INod ( INode, IEl ) ;
          If ( NodeReOrdering ( Node ) == 0 ) Then ;
            NJCounter = NJCounter + 1;
            NodeReOrdering ( Node ) = NJCounter ;
          End If ;
        Counter = Counter + 1 ;
        eind ( Counter ) = NodeReOrdering ( Node )
      End Do ;

  End Do ;

NJG = Counter ;
eptr ( NEl + 1 ) = Counter + 1 ;
ElmWgt = 0 ;
WgtFlag = 0 ;
NumFlag = 1 ;

If      ( NDim  == 2 ) Then ; NCommonNodes = 2 ;
Else If ( NDim  == 3 ) Then ; NCommonNodes = 3 ;  ! check this 3D cases
End If ;

tpwgts (:)   = 1.0/NParts ;
ubvec  (:)   = 1.05 ;
Poptions (0) = 0 ;  ! 0: defined parameteres 1: user defined parameters
Poptions (1) = 0 ;  ! Timing 1
Poptions (2) = 0 ; ! ???
!Poptions (3) = PARMETIS_PSR_UNCOUPLED
Write(*    ,*) 'End Subroutine < ParConvert >' ;
Write(UnInf,*) 'End Subroutine < ParConvert >' ;
Return ;
End Subroutine ParConvert ;


! ***************************************************************************************************************************************************
!                                                                                                                                                  **
! Start:        11 August 2013                                                                                                                     **
! Last Update:  12 August 2013                                                                                                                     **
! Description: This subroutine prpares reguired information for the mesh partitioner, Metis 5.1. See page 28 of the manual for Metis 5.1.0.        **
!                                                                                                                                                  **
! Developed by: Babak Poursartip                                                                                                                   **
!                                                                                                                                                  **
! Called by:                                                                                                                                       **
!                                                                                                                                                  **
! ***************************************************************************************************************************************************

Subroutine MetisConvert_V51    (                                                                                                &
NDim,                                                                                                                           & ! Integer (1) Variables
NParts,                                                                                                                         & ! Integer (2) Variables
NCommon,                                                                                                                        & ! Integer (4) Variables
NEl, NJ, Neind, NJG,                                                                                                            & ! Integer (8) Variables
!                                                                                                                               & ! Real Variables
ELT,     eptr, eind, Vwgt, VSize,      INod,                                                                                    & ! Integer Arrays
tpwgts                                                                                                                          & ! Real Arrays
!                                                                                                                               & ! Characters
!                                                                                                                               & ! Type
) ;

!Metis Function:
!int METIS PartMeshDual(idx t *ne, idx t *nn, idx t *eptr, idx t *eind, idx t *vwgt, idx t *vsize,
!idx t *ncommon, idx t *nparts, real t *tpwgts, idx t *options, idx t *objval,
!idx t *epart, idx t *npart)

Implicit None ;

! =========================== Global Variables ======================================================================================================

! - Integer Variables -------------------------------------------------------------------------------------------------------------------------------
Integer (Kind=Tiny), Intent(In)    :: NDim ;
Integer (Kind=Shrt), Intent(In)    :: NParts ;
Integer (Kind=Shrt), Intent(Out)   :: NCommon ;
Integer (Kind=Lng ), Intent(In)    :: NEl, NJ ;
Integer (Kind=Lng ), Intent(Out)   :: Neind ;
Integer (Kind=Shrt), Intent(Out)   :: NJG ;


! - Real Variables ----------------------------------------------------------------------------------------------------------------------------------
!#Real (Kind=Dbl),     Intent(In)    ::  ;
!#Real (Kind=Dbl),     Intent(InOut) ::  ;
!#Real (Kind=Dbl),     Intent(OUT)   ::  ;
! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
!#Complex,             Intent(In)    ::  ;
!#Complex,             Intent(InOut) ::  ;
!#Complex,             Intent(OUT)   ::  ;
! - Integer Arrays ----------------------------------------------------------------------------------------------------------------------------------
Integer (Kind=Smll), Intent(IN),     Dimension (:)    :: ELT ;
Integer (Kind=Shrt), Intent(Out),    Dimension (:)    :: eptr, eind, Vwgt, VSize ;
Integer (Kind=Lng ), Intent(IN),     Dimension (:,:)  :: INod ;

! - Real Arrays -------------------------------------------------------------------------------------------------------------------------------------
Real (Kind=DBL),     Intent(Out),    Dimension (:)    :: tpwgts ;

! - Character Variables -----------------------------------------------------------------------------------------------------------------------------
! - Logical Variables -------------------------------------------------------------------------------------------------------------------------------
! - Type DECLERATIONS -------------------------------------------------------------------------------------------------------------------------------

! - PETSC INTERNAL Variables ------------------------------------------------------------------------------------------------------------------------
!PetscMPIInt    :: Size ;                         ! Total number of ranks
!PetscMPIInt    :: Rank ;                         ! Rank number

! =========================== Local Variables =======================================================================================================
! - Integer Variables -------------------------------------------------------------------------------------------------------------------------------
Integer (Kind=Tiny)  :: NNode ;                  ! Number of Nodes of element in my code
Integer (Kind=Tiny)  :: NNodeG ;                 ! Number of Nodes of element in ParMetis
Integer (Kind=Tiny)  :: INode ;                  ! Loop index over NNode
Integer (Kind=Lng )  :: ETypeG ;                  ! Element Type
Integer (Kind=Lng )  :: NJCounter, Counter ;     ! Counter
Integer (Kind=Lng )  :: IEl ;                    ! Loop Index over NEl
Integer (Kind=Lng )  :: Node ;                   ! Node Number

! - Real Variables ----------------------------------------------------------------------------------------------------------------------------------
! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
! - Integer Arrays ----------------------------------------------------------------------------------------------------------------------------------
Integer (Kind=Shrt), Dimension ( NJ ) :: NodeReOrdering ;       ! Holds the reordered node numbers of elements. For second order elements, midnodes are eliminated and the rest are saved in this array.

! - Real Arrays -------------------------------------------------------------------------------------------------------------------------------------
! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
! - Character Variables -----------------------------------------------------------------------------------------------------------------------------
! - Logical Variables -------------------------------------------------------------------------------------------------------------------------------
! - Type DECLERATIONS -------------------------------------------------------------------------------------------------------------------------------

! =========================== Subroutine CODE =======================================================================================================

! Converting second order elements to the first order ones by renumbering the corner nodes
NodeReOrdering (:) = 0_Shrt ;
NJCounter          = 0_Lng  ;
Counter            = 0_Lng  ;

! IMPORTANT NOTICE:  Not all the elements in the code considered here. Double check please!!!
  Do IEl = 1, NEl ;

    ETypeG = ELT ( IEl ) ;

      ! Call corresponding GlobalMatrices subroutine
      ! Ordinary Finite Elements
      IF      ( ETypeG == El2d4NSldPS ) Then ;  ! 4 noded - 2D - Solid - PLANE STRESS
        NNode = 4 ; NNodeG = 4 ;
        Vwgt ( IEl ) = 8_Shrt ;

      Else If ( ETypeG == El2d8NSldPN ) Then ;  ! 8 noded - 2D - Solid - PLANE STRAIN
        NNode = 8 ; NNodeG = 4 ;
        Vwgt ( IEl ) = 16_Shrt ;

      Else If ( ETypeG == El2d3NSldPS ) Then ;  ! 3 noded - 2D - Solid - PLANE STRAIN
        NNode = 3 ; NNodeG = 3 ;
        Vwgt ( IEl ) = 6_Shrt ;

      Else If ( ETypeG == El2d8NPMLPN ) Then ;  ! 8 noded - 2D - PML - PLANE STRAIN
        NNode = 8 ; NNodeG = 4 ;
        Vwgt ( IEl ) = 40_Shrt ;

      Else If ( ETypeG == El2d6NSldPN ) Then ;  ! 6 noded - 2D - Solid - PLANE STRAIN
        NNode = 6 ; NNodeG = 3 ;
        Vwgt ( IEl ) = 12_Shrt ;

      Else If ( ETypeG == El2d6NPMLPN ) Then ;  ! 6 noded - 2D - PML - PLANE STRAIN
        NNode = 6 ; NNodeG = 3 ;
        Vwgt ( IEl ) = 30_Shrt ;

      ! Spectral elements - Regular domain
      Else If ( ETypeG == SEl2d9NSldPN .OR.  ETypeG == SEl2d9NSldPS ) Then ;  ! 9 noded - 2D 
        NNode = 9 ; NNodeG = 4 ;
        Vwgt ( IEl ) = 18_Shrt ;

      Else If ( ETypeG == SEl2d7NSldPN .OR. ETypeG == SEl2d7NSldPS ) Then ;  ! 7 noded - 2D
        NNode = 7 ; NNodeG = 3 ;
        Vwgt ( IEl ) = 14_Shrt ;

      Else If ( ETypeG == SEl3d27NSld ) Then ;  ! 27 noded - 3D
        NNode = 27 ; NNodeG = 8 ;
        Vwgt ( IEl ) = 81_Shrt ;

      Else If ( ETypeG == SEl3d8NSld ) Then ;  ! 8 noded - 3D
        NNode = 8 ; NNodeG = 8 ;
        Vwgt ( IEl ) = 24_Shrt ;

      ! Spectral elements - PML
      Else If ( ETypeG == SEl2d9NPMLPN .OR. ETypeG == SEl2d9NMPMLPN .OR. ETypeG == SEl2d9NPMLPS .OR. ETypeG == SEl2d9NMPMLPS ) Then ;  ! 9 noded - 2D 
        NNode = 9 ; NNodeG = 4 ;
        Vwgt ( IEl ) = 45_Shrt ;

      Else If ( ETypeG == SEl2d7NMPMLPN .OR. ETypeG == SEl2d7NMPMLPS ) Then ;  ! 7 noded - 2D
        NNode = 7 ; NNodeG = 3 ;
        Vwgt ( IEl ) = 35_Shrt ;

      Else If ( ETypeG == SEl3d27NMPML ) Then ;  ! 27 noded - 3D
        NNode = 27 ; NNodeG = 8 ;
        Vwgt ( IEl ) = 243_Shrt ;

      Else If ( ETypeG == SEl3d8NMPML ) Then ;  ! 8 noded - 3D
        NNode = 8 ; NNodeG = 8 ;
        Vwgt ( IEl ) = 72_Shrt ;

      Else ; ! Element Type is not in the list
        Write(*, Fmt_Element2) ;  Write(UnInf, Fmt_Element2) ;
        Write(*, Fmt_FL) ; Write(UnInf, Fmt_FL) ;
        !#Call BEEP_FAIL ;
        Write(*, Fmt_End) ; Read(*,*) ;  STOP ;
      End If ;

    ! Remarks: 1- Node numbering starts from one. Fortran style. Check it.
    !          2- Vectors start 1. Fortran Style.
    eptr ( IEl ) = Counter + 1 ;
      Do INode = 1, NNodeG ;
        Node = INod ( INode, IEl ) ;
          If ( NodeReOrdering ( Node ) == 0 ) Then ;
            NJCounter = NJCounter + 1;
            NodeReOrdering ( Node ) = NJCounter ;
          End If ;
        Counter = Counter + 1 ;
        eind ( Counter ) = NodeReOrdering ( Node )
      End Do ;

  End Do ;

NJG = NJCounter ;
Neind = Counter ;
eptr ( NEl + 1 ) = Counter + 1 ;
!VWgt(:) = 1_Shrt ;         ! Vector of weights. Modify this for PML elements.
VSize(:) = 0_Shrt ;         ! Vector of the size of the elements -see page 28 of the manual.

If      ( NDim  == 2 ) Then ; NCommon = 2 ;
Else If ( NDim  == 3 ) Then ; NCommon = 4 ;  ! check this unstructred 3D models
End If ;

!tpwgts (:)   = 1.0/NParts ;  ! Set it equal to NULL ! of size nparts that specifies the desired weight for each partition. 
tpwgts (:)   = 0.0_Dbl ;  ! Set it equal to NULL ! of size nparts that specifies the desired weight for each partition. 

Write(*    ,*) 'End Subroutine < MetisConvert_V51 >' ;
Write(UnInf,*) 'End Subroutine < MetisConvert_V51 >' ;
Return ;
End Subroutine MetisConvert_V51 ;

End Module ConvertToMetis ;


