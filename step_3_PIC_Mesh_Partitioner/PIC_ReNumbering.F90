
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                                                                                                                  ++
! Start:        20 August 2011                                                                                                                     ++
! Last Update:  15 June 2012                                                                                                                       ++
!                                                                                                                                                  ++
! Description: THIS Module renumbers nodes and elements for each rank. It also finds PETSc numbering and its relation with application numbering.  ++
!                                                                                                                                                  ++
!              Subroutines                     Description                                                                                         ++
!              ==============                  ==============                                                                                      ++
!                                                                                                                                                  ++
!                                                                                                                                                  ++
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Module ReNumbering ;

Use Parameters ;

Implicit None ;

  Interface
!    Module Procedure 
  End Interface    ;

Contains ;

! ***************************************************************************************************************************************************
!                                                                                                                                                  **
! Last Update:  26 Oct 2011                                                                                                                        **
! Description: THIS Subroutine CALCULATES                                                                                                          **
!                                                                                                                                                  **
! ***************************************************************************************************************************************************

Subroutine Numbering    (                                                                                                    &
MaxNNode, NDOF, NNeighbor, NParts, NEL, NJ, NEQMTotal,                                                                       & ! Integer Variables
!                                                                                                                            & ! Real Variables
NPart, EPart, INod, ID, NEL_Rank, NJ_Rank, Global_PETSc_Num, NEqRank, NNodeRank, Local_PETSc_Num, ID_Application, ID_PETSc,  & ! Integer Arrays
!                                                                                                                            & ! Real Arrays
!                                                                                                                            & ! Characters
Nodes                                                                                                                        & ! Type 
) ;

Implicit None ;

! =========================== Global Variables ======================================================================================================

! - Integer Variables -------------------------------------------------------------------------------------------------------------------------------
Integer (Kind=Tiny), Intent(In)    :: MaxNNode, NDOF ;
Integer (Kind=Smll), Intent(InOut) :: NNeighbor ;

Integer (Kind=Shrt), Intent(In)    :: NParts ;

Integer (Kind=Lng ), Intent(In)    :: NEL, NJ ;
Integer (Kind=Lng ), Intent(Out)   :: NEQMTotal ;

! - Real Variables ----------------------------------------------------------------------------------------------------------------------------------
! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
! - Integer Arrays ----------------------------------------------------------------------------------------------------------------------------------
Integer (Kind=Shrt), Intent(IN ),   Dimension (:  )  :: EPart ;
Integer (Kind=Shrt), Intent(INOut), Dimension (:  )  :: NPart ;

Integer (Kind=Lng ), Intent(IN ), Dimension (:,:)  :: INod, ID ;
Integer (Kind=Lng ), Intent(OUT), Dimension (:  )  :: NEL_Rank, NJ_Rank, Global_PETSc_Num, NEqRank, NNodeRank ;
Integer (Kind=Lng ), Intent(OUT), Dimension (:,:)  :: Local_PETSc_Num, ID_Application, ID_PETSc ;

! - Real Arrays -------------------------------------------------------------------------------------------------------------------------------------
! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
! - Character Variables -----------------------------------------------------------------------------------------------------------------------------
! - Logical Variables -------------------------------------------------------------------------------------------------------------------------------

! - Type DECLERATIONS -------------------------------------------------------------------------------------------------------------------------------
Type ( NodeID ) :: Nodes ;                                           ! Node%Locs :  Holds rank numbers of which this node belongs to.
                                                                     ! Node%Rep  : Holds the number of Repeatations on the ranks (How many times this node appears on ranks), Useful for neighboring

! =========================== LOCAL Variables =======================================================================================================
! - Integer Variables -------------------------------------------------------------------------------------------------------------------------------
Integer (Kind=Tiny)  :: INode ;                  ! Counter for Node Number
Integer (Kind=Tiny)  :: NNode ;                  ! Number of Nodes of element

Integer (Kind=Smll)  :: INeighbor ;              ! Counter for number of neighbors 
Integer (Kind=Smll)  :: IDOF ;                   ! Loop index over NDOF

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

! - Real Variables ----------------------------------------------------------------------------------------------------------------------------------
! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
! - Integer Arrays ----------------------------------------------------------------------------------------------------------------------------------
! - Real Arrays -------------------------------------------------------------------------------------------------------------------------------------
! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
! - Character Variables -----------------------------------------------------------------------------------------------------------------------------
! - Logical Variables -------------------------------------------------------------------------------------------------------------------------------
Logical   :: Check ;                             ! Used to find the end of while loop
! - Type DECLERATIONS -------------------------------------------------------------------------------------------------------------------------------

! =========================== Subroutine CODE =======================================================================================================

Write (*,*)"ReNumbering:: Renumbering nodes for PETSc numbering ..."
Write (*,*)"ReNumbering:: Calculating number of nodes of each rank ..."
! - Calculating the total number of elements of each partition --------------------------------------------------------------------------------------
!?? Be careful if ranks start from 0
NEL_Rank = 0_Lng ;
  Do IEL = 1, NEL ;
    NEL_Rank ( EPart ( IEL )) = NEL_Rank ( EPart ( IEL )) + 1 ;
  End Do ; 

Write (*,*)"ReNumbering:: Locating nodes of each rank ..."

! - Indicating the location and repeats of each node in ranks ---------------------------------------------------------------------------------------
Check = .False. ;
  Do While ( Check == .False. ) ; ! The check varible ensures the adequacy of the number of defined neighbors ( NNeighbor )
    Allocate ( Nodes%Locs ( NJ, NNeighbor ), Nodes%Rep ( NJ ) ) ;
    Nodes%Locs(:,:) =  0_Shrt ; ! Holds rank numbers which a node belongs to 
    Nodes%Rep (:  ) =  0_Shrt ; ! Holds how many times a node is repeated in ranks
    Check = .True. ;

    MLoop:  Do IParts = 1, NParts ;

              Do IEL = 1, NEL ;
                IF ( EPart ( IEL ) == IParts ) Then ;

                    Do INode = 1, MaxNNode ;

                      Node = INod ( INode, IEL ) ; ! Node is the "Global Node NUmbering (Application Ordering)"
                      If (Node == -1) Exit ;   ! This case happens when different element types with different NNode exist in the same model. INod was initialized to -1 so that we can find null nodes.
                        Do INeighbor = 1, NNeighbor ;
                          If      ( Nodes%Locs ( Node, INeighbor ) == IParts ) Then ; 
                            Exit ; ! Node rank was indicated before
                          Else If ( Nodes%Locs ( Node, INeighbor ) == 0_Shrt ) Then ; 
                            Nodes%Rep ( Node ) = Nodes%Rep ( Node ) + 1_Smll ;
                            Nodes%Locs ( Node, INeighbor ) = IParts ;
                            Exit ;
                          End If ; 
                        End Do ;

                        ! If number of defined neighbors is less than the number of existing neighbors we have to increase NNeighbor. Notice the more neighbors you define here, the more memry the code needs. Be careful here.
                        If ( INeighbor == NNeighbor + 1 ) Then ;
                          Write (*,*) "Maximmum Number of PreDefined Neighbors is not enough. The Code will increasing this value automatiCally." ;
                          NNeighbor = NNeighbor + 5_Smll ; ! 5 is chosen for no special reason. Just a guess that might fit with the existing number of neighbors.
                          DEAllocate ( Nodes%Locs, Nodes%Rep ) ;
                          Check = .False.
                          Exit MLoop ;
                        End If ; 

                    End Do ;
                End If ;
              End Do ;
      End Do MLoop ;

  End Do ; ! End Do While

! - ReOrdering the node numbers for "Global PETSC Numbering" and number of equations of each rank ---------------------------------------------------
! Comments:
! 1- Refer to PETSC manual for Global PETSC Node Numbering (page 50). 
! 2- We have ghost nodes on each rank due to neighbor elements. In this code, the developer considered the ghost nodes belong to lower rank in PETSc numbering.
! 3- We first number the nodes that belong only to each rank, than, we number the ghost nodes. We use the variable repeatation for this aim.

Write (*,*)"ReNumbering:: PETSc global numbering ..."

GPETScNodeNum        = 0_Lng ; ! Global PETSc Node Numbering
Global_PETSc_Num (:) = 0_Lng ; ! Global PETSc Numbering

  Do IParts = 1_Lng, NParts ;  ! Numbering Nodes by rank

        ILocalRow = 0_Lng ;

          Do IJ = 1_Lng, NJ ;

            Do INeighbor = 1_Lng, NNeighbor ;
              If ( Nodes%Locs ( IJ, INeighbor ) == 0_Shrt ) Exit ; ! We already checked all ranks that this node belongs to.
              If ( Nodes%Locs ( IJ, INeighbor ) == IParts .AND. Nodes%Rep ( IJ ) == 1_Smll ) Then ;   ! Numbering the interior nodes

                  GPETScNodeNum = GPETScNodeNum + 1 ;
                  Global_PETSc_Num ( IJ ) = GPETScNodeNum ;
! ----------------- added for inversion DS ------
                  NPart(IJ) = IParts
! -----------------------------------------------

                  ! Finding number of equations related to this node and rank
                  ILocalRow = ILocalRow + Count ( ID ( IJ, : ) == 0_Smll ) ;
                  Exit ;

              End If ;

              ! Numbering the ghost nodes  
              If ( Nodes%Locs ( IJ, INeighbor ) == IParts .AND. Nodes%Rep ( IJ ) /= 1_Smll .AND. Global_PETSc_Num ( IJ ) == 0_Lng ) Then ;

                  GPETScNodeNum = GPETScNodeNum + 1 ;
                  Global_PETSc_Num ( IJ ) = GPETScNodeNum ;
! ----------------- added for inversion DS ------
                  NPart(IJ) = IParts
! -----------------------------------------------

                  ! Finding number of equations related to this node and rank
                  ILocalRow = ILocalRow + Count ( ID ( IJ, : ) == 0_Smll ) ;
                  Exit ;

              End If ;

            End Do ; ! End loop on NNeighbors
          End Do ; ! End loop on the Nparts
        NEqRank   ( IParts ) = ILocalRow ;      ! Number of Equations each rank stores in its own memory.
        NNodeRank ( IParts ) = GPETScNodeNum ;  ! Maximum PETSc node number saved on this rank.

  End Do ;

! Modifying Equation Numbers
  Do IParts = 2_Shrt, NParts ;  !
    NEqRank ( IParts ) = NEqRank ( IParts ) + NEqRank ( IParts -1_Shrt ) ; ! Now, NEqRank(i) indicates the maximum equation number that is saved on rank i.
  End Do ;

! - PETSc Local Node Numbering ----------------------------------------------------------------------------------------------------------------------
Local_PETSc_Num = 0_Lng ;

Write (*,*)"ReNumbering:: Renumbering nodes of each rank ..."

  Do IParts = 1_Shrt, NParts ;
    LPETScNodeNum = 0_Lng ; ! Local PETSc Node Numbering
      Do IJ = 1_Lng, NJ ;
        Do INeighbor = 1_Lng, NNeighbor ;
          If ( Nodes%Locs ( IJ, INeighbor ) == 0_Shrt ) Exit ;
          If ( Nodes%Locs ( IJ, INeighbor ) == IParts ) Then ;
            LPETScNodeNum = LPETScNodeNum + 1 ;
            Local_PETSc_Num ( IJ, IParts ) = LPETScNodeNum ;
            Exit ;
          End If ;
        End Do ; 
      End Do ;
  End Do ;

! - Number of nodes of each rank --------------------------------------------------------------------------------------------------------------------
NJ_Rank = MaxVal ( Local_PETSc_Num, Dim = 1 ) ; 

! - Application Equation Number ---------------------------------------------------------------------------------------------------------------------
! THIS ARRANGES X, Y, SXX, SYY, SXY degrees of freedom concessively(???) for each node
ApplicationEqN = 0_Lng ;
ID_Application = ID ;

Write (*,*)"ReNumbering:: Application equation numbering..."

  DO I = 1, NJ ;
    DO J = 1, NDOF ;
      IF ( ID_Application ( I, J ) == 1 ) Then ;
        ID_Application ( I, J ) = 0 ; 
      Else If ( ID_Application ( I, J ) == 0 ) Then ;
        ApplicationEqN = ApplicationEqN + 1 ;
        ID_Application ( I, J ) = ApplicationEqN ;
      Else ;
        Write (*,*)"A major mistake in constraints array - see numbering subroutine"
        Stop;
      End If ;
    End Do ;
  End Do ;

NEQMTotal = ApplicationEqN ;

! - PETSc Equation number ---------------------------------------------------------------------------------------------------------------------------
ForAll ( IJ = 1:NJ, IDOF = 1:NDOF ) ID_PETSc ( Global_PETSc_Num ( IJ ), IDOF ) = ID ( IJ, IDOF ) ;

PETScEqN = 0_Lng ;

Write (*,*)"ReNumbering:: PETSc equation numbering..."

  DO I = 1, NJ ;
    DO J = 1, NDOF ;
      IF ( ID_PETSc ( I, J ) == 1 ) Then ;
        ID_PETSc ( I, J ) = 0 ; 
      Else If ( ID_PETSc ( I, J ) == 0 ) Then ;
        PETScEqN = PETScEqN + 1 ;
        ID_PETSc ( I, J ) = PETScEqN ;
      End If ;
    End Do ;
  End Do ;

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ check

!Write (*,*)"nodes%locs" 
!do iel = 1, nj 
!Write (*,"(I4,3X,<nneighbor>(i5,2x),<1>(i5,2x),<1>(i5,2x),<2>(i5,2x))") IEL, (nodes%locs(iel, ij), ij = 1, nneighbor), nodes%rep(iel), Global_PETSc_Num (iel), (local_PETSc_Num (iel, ij), ij =1, NParts)
!End Do

!Write (*,*)"id" 
!do iel = 1, nj 
!Write (UN_chk,"(I4,3X,<ndof>(i5,2x),6X,<ndof>(i5,2x))") IEL, (id_application (iel,ij) -1 , ij=1,ndof),(id_petsc (iel,ij)-1, ij=1,ndof)
!Write (UN_chk,"(I4,3X,<ndof>(i5,2x),6X,<ndof>(i5,2x))") IEL, (id_application (iel,ij)  , ij=1,ndof),(id_petsc (iel,ij), ij=1,ndof)
!End Do

!Write (*,*)"nodes%rep"
!do iel = 1, nj 
!Write (*,"(I4,3X,<1>(i5,2x))") IEL, nodes%rep(iel)
!End Do

!Write (*,*)"Global NUMBERING"
!do iel = 1, nj 
!Write (*,"(I4,3X,<1>(i5,2x))") IEL, Global_PETSc_Num (iel)
!End Do

!Write (*,*)"LOCAL NUMBERING"
!do ij = 1, NParts
!Write (*,*)"LOCAL NUMBERING", ij
!do iel = 1, nj 
!Write (*,"(I4,3X,<1>(i5,2x))") IEL, local_PETSc_Num (iel, ij)
!End Do
!End Do 

Write(*    ,*) 'End Subroutine < Numbering >' ;
Write(UnInf,*) 'End Subroutine < Numbering >' ;
Return ;
End Subroutine Numbering ;

End Module ReNumbering ;
