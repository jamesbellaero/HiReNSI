
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                                                                                                                  ++
! Start:         6  Januardy 2014                                                                                                                  ++
! Last Update:   6  Januardy 2014                                                                                                                  ++
!                                                                                                                                                  ++
! Description: This Module constructs appropriate data structures and mappings required for inversion                                              ++
!              I. building material vectors and associated mappings and scatters                                                                   ++
!              II. material extension from interface into PML domain                                                                               ++
!              III. data structure for storing displacement, associated mappings and scatters                                                      ++
!                                                                                                                                                  ++
!              Subroutines                     Description                                                                                         ++
!              ==============                  ==============                                                                                      ++
!                                                                                                                                                  ++
!                                                                                                                                                  ++
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Module Inversion_Data_Structure ;

Use Parameters ;

Implicit None ;

  Interface
!    Module Procedure 
  End Interface    ;

Contains ;

! ***************************************************************************************************************************************************
!                                                                                                                                                  **
! Last Update:  6  Januardy 2014                                                                                                                   **
!                                                                                                                                                  **
! ***************************************************************************************************************************************************

Subroutine Inversion_DS (                                                                                                    &
NDim, MaxNNode, NDOF, NParts, NEL, NJ, NEQMTotal,                                                                            & ! Integer Variables
!                                                                                                                            & ! Real Variables
EPart, INod, ID, NEL_Rank, NJ_Rank, Global_PETSc_Num, NEqRank, NNodeRank, Local_PETSc_Num, ID_Application, ID_PETSc,         & ! Integer Arrays
NPart, Node_Mat_ID, Node_Mat_Mapping,                                                                                        &
!                                                                                                                            & ! Real Arrays
ModelName, OutDir                                                                                                            & ! Characters
!                                                                                                                            & ! Type 
) ;

Implicit None ;

! =========================== Global Variables ======================================================================================================

! - Integer Variables -------------------------------------------------------------------------------------------------------------------------------
Integer (Kind=Tiny), Intent(In)    :: NDim, MaxNNode, NDOF ;

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

Integer (Kind=Shrt), Intent(In), Dimension (:  )  :: NPart 
Integer (Kind=Smll), Intent(In), Dimension (:  )  :: Node_Mat_ID;
Integer (Kind=Lng),  Intent(In), Dimension (:  )  :: Node_Mat_Mapping;

! - Real Arrays -------------------------------------------------------------------------------------------------------------------------------------
! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
! - Character Variables -----------------------------------------------------------------------------------------------------------------------------
! - Logical Variables -------------------------------------------------------------------------------------------------------------------------------
! - Type DECLERATIONS -------------------------------------------------------------------------------------------------------------------------------

! =========================== LOCAL Variables =======================================================================================================
! - Integer Variables -------------------------------------------------------------------------------------------------------------------------------
Integer (Kind=Tiny)  :: INode ;                  ! Counter for Node Number
Integer (Kind=Tiny)  :: NNode ;                  ! Number of Nodes of element

Integer (Kind=Tiny)  :: IDOF ;                   ! Loop index over NDOF
!Integer (Kind=Smll)  :: IDim ;                   ! Loop index over NDim

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

Integer (Kind=Lng )  :: NJ_Mapping ;             ! Global nodes owned by each rank
Integer (Kind=Lng )  :: NStore_Mapping
Integer (Kind=Lng )  :: NStore_Rank

Integer (Kind=Lng )  :: Counter ; 
Integer (Kind=Lng )  :: LEqN ;                   ! Local Equation number, for ID_Local.

! - Real Variables ----------------------------------------------------------------------------------------------------------------------------------

! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
! - Integer Arrays ----------------------------------------------------------------------------------------------------------------------------------
Integer (Kind=Lng ), Allocatable, Dimension(:)    :: idx_Mat_from;
Integer (Kind=Lng ), Allocatable, Dimension(:)    :: idx_Mat_to;

Integer (Kind=Lng ), Allocatable, Dimension(:)    :: idx_Mat_Extend_from;
Integer (Kind=Lng ), Allocatable, Dimension(:)    :: idx_Mat_Extend_to;
 
Integer (Kind=Lng ), Allocatable, Dimension(:)    :: U_Store_Numbers_Global;
Integer (Kind=Lng ), Allocatable, Dimension(:)    :: idx_u_from;
Integer (Kind=Lng ), Allocatable, Dimension(:)    :: idx_u_to;

Integer (Kind=Lng ), Allocatable, Dimension(:,:)  :: ID_Local;      

! - Real Arrays -------------------------------------------------------------------------------------------------------------------------------------


! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
! - Character Variables -----------------------------------------------------------------------------------------------------------------------------
! - Character Variables -----------------------------------------------------------------------------------------------------------------------------
Character (Kind = 1, Len = 20)  :: IndexRank ;   ! A variable for storing the Rank number in Character format for adding at the end of the input file Name.
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
Integer (Kind=Smll), PARAMETER  :: Un_Inversion_DS     = 811  ;            ! the Unit number of Inversion Data Structure index sets

! =========================== Subroutine CODE =======================================================================================================
write(*,*) "Begin subroutine < InversionDS >"
  Do IParts = 1_Shrt, NParts ;

! I. building material vectors, associated mappings and scatters ------------------------------------------------------------------------------------ checked, works!


!   Global PETSc#           Local rank # 1     Local rank # 2

!  3 ---- 6   ---- 9           3 ---- 6           3 ---- 6
!  |      |        |           |      |           |      | 
!  |      |        |           |      |           |      |
!  2 ---- 5   ---- 8           2 ---- 5           2 ---- 5
!  |      |        |           |      |           |      |
!  |      |        |           |      |           |      |
!  1 ---- 4   ---- 7           1 ---- 4           1 ---- 4


! rank # 1: idx_Mat_from = [1 2 3 4 5 6] : Global numbers including ghost nodes (similar to Indices in PIC_Output.F90 (782)). 
! rank # 2: idx_Mat_from = [4 5 6 7 8 9] : Local numbering

! rank # 1: idx_Mat_to   = [1 2 3 4 5 6]
! rank # 2: idx_Mat_to   = [1 2 3 4 5 6]

! remark: NJ_Rank_IParts = NJ (NJ local in the main code, not here)

    Allocate ( idx_Mat_from ( NJ_Rank(IParts) ) )
!    Allocate ( idx_Mat_to   ( NJ_Rank(IParts) ) )

    ! Mapping between the local numbering and Global PETSc numbering - size of idx_Mat_from is NJ_Rank(IParts).
    Do IJ = 1_Lng, NJ ! (this NJ is global)
      If ( Local_PETSc_Num ( IJ, IParts ) /= 0_Lng ) Then
         idx_Mat_from ( Local_PETSc_Num ( IJ, IParts ) ) = Global_PETSc_Num ( IJ )
!         idx_Mat_to   ( Local_PETSc_Num ( IJ, IParts ) ) = Local_PETSc_Num ( IJ, IParts )
      End If 
    End Do
    ! idx_Mat_to can be constructed in the main program as follows:
    ! ForAll ( I = 1:NJ_Rank ( IParts ) ) idx_Mat_to (I) = I-1

! --------------------------------------------------------------------------------------------------------------------------------------------------- 

! II. material extension from interface into PML domain --------------------------------------------------------------------------------------------- checked, works!


! Q. Can we do vecscatter from a vector into itself?
! A. Yes. (See test_01.F90)


    ! Obtanining the size of Global nodes owned by each rank. -------------------------------------
    
    If ( IParts == 1_Shrt ) then ;
      NJ_Mapping = NNodeRank ( IParts ) ;
    Else ;
      NJ_Mapping = NNodeRank ( IParts ) - NNodeRank ( IParts - 1_Shrt ) ;
    End If ;
    Allocate ( idx_Mat_Extend_from (NJ_Mapping), idx_Mat_Extend_to (NJ_Mapping) )
    write(*,*) "NJ: ",NJ_Mapping
    ! Construct index sets - size of idx_Mat_Extend_from / idx_Mat_Extend_to is NJ_Mapping. ------- $ $ $ $ $ $ $ $ $ $
    Counter = 0_Lng
    Do IJ = 1_Lng, NJ
      If ( NPart ( IJ ) == IParts ) Then
       ! if(Counter>NJ_Mapping) then
        !  write(*,*) "shite",Counter,NJ_Mapping,IJ,size(Node_Mat_Mapping),size(Global_PETSc_Num),size(idx_Mat_Extend_from);
         ! write(*,*) Node_Mat_Mapping(IJ);
        !end if
        
        Counter = Counter + 1_Lng 
        idx_Mat_Extend_from ( Counter ) = Global_PETSc_Num ( Node_Mat_Mapping ( IJ ) )
        idx_Mat_Extend_to   ( Counter ) = Global_PETSc_Num (                  ( IJ ) )
      End If
    End Do
    
    ! Double check -------------------------------------------------------------------------------- $ $ $ $ $ $ $ $ $ $
    If ( Counter /= NJ_Mapping ) Then
       Write(*,*) ' Wrong in II Inv-DS '; Stop
    End If

! ---------------------------------------------------------------------------------------------------------------------------------------------------

! III. data structure for storing displacement, associated mappings and scatters --------------------------------------------------------------------


! Q. What nodes should be stored?
! A. Those with Node_Mat_ID = 0 or 1 should be stored (0 = RD, 1 = I).


!     Global #           Local rank # 1     Local rank # 2

!  3 ---- 6   ---- 9*          3 ---- 6           3 ---- 6*
!  |      |        |           |      |           |      | 
!  |      |        |           |      |           |      |
!  2*---- 5*  ---- 8           2*---- 5*          2*---- 5
!  |      |        |           |      |           |      |
!  |      |        |           |      |           |      |
!  1 ---- 4   ---- 7           1 ---- 4           1 ---- 4

! '*' indicates variable should be saved; i.e. Node_Mat_ID = 0 or 1 (only 1 dof per node in the above example for simplicity).

! rank # 1: U_Store_Numbers_Global = [2 5]  NStore_Mapping = 2                  
! rank # 2: U_Store_Numbers_Global = [9]    NStore_Mapping = 1

! Global:
! rank # 1: idx_u_from = [2 5]   NStore_Rank = 2
! rank # 2: idx_u_from = [5 9]   NStore_Rank = 2

! Local:
! rank # 1: idx_u_to   = [2 5]   NStore_Rank = 2
! rank # 2: idx_u_to   = [2 6]   NStore_Rank = 2


    ! Number of Equations each rank should store in its own memory. -------------------------------
    Counter = 0_Lng
    Do IJ = 1_Lng, NJ
      If ( NPart ( IJ ) == IParts .AND. ( Node_Mat_ID ( IJ ) == 0 .OR. Node_Mat_ID ( IJ ) == 1 ) ) Then
        Counter = Counter + NDim                 ! Need to store the three components of the displacement vector.
      End If
    End Do
    NStore_Mapping = Counter ;                   ! Number of Equations each rank stores in its own memory.


    ! Number of locally stored equations on each rank, including ghost nodes ----------------------
    ! (for assembly and gradient computation on each rank).
    Counter = 0_Lng
    Do IJ = 1_Lng, NJ
      If ( Local_PETSc_Num ( IJ, IParts ) /= 0_Lng .AND. ( Node_Mat_ID ( IJ ) == 0 .OR. Node_Mat_ID ( IJ ) == 1 ) ) Then
          Counter = Counter + NDim               ! Need to store the three components of the displacement vector.
      End If
    End Do
    NStore_Rank = Counter ;     


    ! Find which equation numbers need to be stored. ----------------------------------------------
    Allocate ( U_Store_Numbers_Global ( NStore_Mapping ), idx_u_from ( NStore_Rank ), idx_u_to ( NStore_Rank ) )
    
    ! Obtaining equation numbers of the nodes that need to be stored. -----------------------------
    Counter = 0_Lng ;
    Do IJ = 1_Lng, NJ ;
      If ( NPart ( IJ ) == IParts .AND. ( Node_Mat_ID ( IJ ) == 0 .OR. Node_Mat_ID ( IJ ) == 1 ) ) Then
        Do IDOF = 1, NDim ;                      ! Only displacement components are needed.

          If ( ID_Application ( IJ, IDOF ) == 0_Lng ) Then ;
            Write(*,*) 'This is wrong; because this dof must be active (we have PML surrounding RD).'
            Stop
          End If ;

          Counter = Counter + 1_Lng ;
          If ( Counter > NStore_Mapping ) Then ; 
            Write(*      , "('Number of equations for storage on this rank is greater than the estimated number - Check Subroutine < >', 3(I19,2x))" ) IParts, IJ, NStore_Mapping ;
            Write( UnFile, "('Number of equations for storage on this rank is greater than the estimated number - Check Subroutine < >', 3(I19,2x))" ) IParts, IJ, NStore_Mapping ;
            Write(*, Fmt_FL) ;  Write(UnInf, Fmt_FL) ; Write(*, Fmt_End) ; Read(*,*) ;  Return ;
          End If ;

          If ( ID_PETSc ( Global_PETSc_Num ( IJ ), IDOF ) == 0 ) Then
             Write(*,*) 'Inversion DS - Line 292 error.'; Stop
          End If

          U_Store_Numbers_Global ( Counter ) = ID_PETSc ( Global_PETSc_Num ( IJ ), IDOF ) ;

        End Do ;
      End If ;
    End Do ;


    ! idx_u_from ----------------------------------------------------------------------------------
    Counter = 0_Lng ;
    Do IJ = 1_Lng, NJ
      If ( Local_PETSc_Num ( IJ, IParts ) /= 0_Lng .AND. ( Node_Mat_ID ( IJ ) == 0 .OR. Node_Mat_ID ( IJ ) == 1 ) ) Then ;
          Do IDOF = 1, NDim
            Counter = Counter + 1_Lng ;
            idx_u_from ( Counter ) = ID_PETSc ( Global_PETSc_Num ( IJ ), IDOF ) ;  ! $ $ $ $ $
          End Do
      End if
    End Do


    ! construct ID_Local --------------------------------------------------------------------------
    Allocate ( ID_Local ( NJ_Rank ( IParts ), NDOF ) ) ;

      ! Forming the Local ID for each Process
      DO IJ = 1_Lng, NJ ;
        If ( Local_PETSc_Num ( IJ, IParts ) /= 0_Lng ) Then ;
          ID_Local ( Local_PETSc_Num ( IJ, IParts ), : ) = ID ( IJ, : ) ;
        End If ;
      End Do ;
    write(*,*) "Obtaining local equation numbers"
    ! Modified Number of Equations for each rank
    ! Obtaining the Local Equation Number
    LEqN = 0_Lng ;
      Do IJ = 1_Lng, NJ_Rank ( IParts ) ;
        DO IDof = 1, NDOF ;
          IF ( ID_Local ( IJ, IDof ) == 1 ) Then ;
            ID_Local ( IJ, IDof ) = 0 ; 
          Else If ( ID_Local ( IJ, IDof ) == 0 ) Then ;
            LEqN = LEqN + 1_Lng ;
            ID_Local ( IJ, IDof ) = LEqN ;
          End If ;
        End Do ;
      End Do ;


    ! idx_u_to ------------------------------------------------------------------------------------
    Counter = 0_Lng ;
    Do IJ = 1_Lng, NJ ;
      If ( Local_PETSc_Num ( IJ, IParts ) /= 0_Lng .AND. ( Node_Mat_ID ( IJ ) == 0 .OR. Node_Mat_ID ( IJ ) == 1 ) ) Then ;
          Do IDOF = 1, NDim ;                    ! Only displacement components are needed.
        
            If ( ID_Local ( Local_PETSc_Num ( IJ, IParts ), IDOF ) == 0_Lng ) Then
              Write(*,*) 'This is wrong; because we have PML surrounding RD.'
              Stop
            End If ;    

            Counter = Counter + 1_Lng ;
            idx_u_to ( Counter ) = ID_Local ( Local_PETSc_Num ( IJ, IParts ), IDOF )  ;  ! $ $ $ $ $

          End Do ;
      End If ;
    End Do ;

! ---------------------------------------------------------------------------------------------------------------------------------------------------

! C-indexing, as used in PETSc. -------------------------------------------------------------------
    idx_Mat_from           = idx_Mat_from           - 1_Lng          ! size: NJ_Rank(IParts)
    idx_Mat_Extend_from    = idx_Mat_Extend_from    - 1_Lng          ! size: NJ_Mapping
    idx_Mat_Extend_to      = idx_Mat_Extend_to      - 1_Lng          ! size: NJ_Mapping
    U_Store_Numbers_Global = U_Store_Numbers_Global - 1_Lng          ! size: NStore_Mapping
    idx_u_from             = idx_u_from             - 1_Lng          ! size: NStore_Rank
    idx_u_to               = idx_u_to               - 1_Lng          ! size: NStore_Rank

! =========================== Output ================================================================================================================

! - Writing down the input data for each process ----------------------------------------------------------------------------------------------------


    Write (IndexRank, *) IParts - 1_Shrt ; ! Converts Rank number to Character foramt for the file Name
    Write (IndexSize, *) NParts ;          ! Converts Size number to Character foramt for the file Name

    Write (*,*)"Writing files for partition: ", IParts;

    ! - Output FILEs --------------------------------------------------------------------------------------------------------------------------------
    Write (*,*)"Opening files for Inversion Data Structure ...";

    ! Inversion_DS - all index sets
    UnFile = Un_Inversion_DS ;
    Open ( Unit = UnFile, FILE = TRIM(ModelName)//'_'//Trim(AdjustL(IndexSize))//'_'//Trim(AdjustL(IndexRank))//'.Inv_DS', ERR =  1001, IOSTAT = IO_File, ACCESS = 'SEQUENTIAL', ACTION ='Write', ASYNCHRONOUS = 'NO', BLANK = 'NULL', BLOCKSIZE = 0, DEFAULTFILE = TRIM(OutDir), DisPOSE = 'KEEP', FORM = 'FORMATTED', POSITION = 'ASIS', STATUS ='REPLACE') ;

! - Write index sets --------------------------------------------------------------------------------------------------------------------------------
    UnFile = UN_Inversion_DS ;
    Write (*,*)"Writing Basic Data ..." ;
    Write (Unit = UnFile, FMT = "(4(I10,1X),'NJ_Rank(IParts) | ','NJ_Mapping | ','NStore_Mapping | ','NStore_Rank')", ADVANCE = 'YES', ASYNCHRONOUS = 'NO', IOSTAT = IO_Write, ERR = 1006 ) NJ_Rank(IParts), NJ_Mapping, NStore_Mapping, NStore_Rank ;
    Write (*,*) "Writing index sets for inversion ..." ;
    Write (Unit = UnFile, FMT = "(<NJ_Rank(IParts)>(I10,1X),'idx_Mat_from')",          ADVANCE = 'YES', ASYNCHRONOUS = 'NO',  IOSTAT = IO_Write, ERR = 1006 ) ( idx_Mat_from             ( I ), I = 1_Lng, NJ_Rank(IParts) ) ;
    Write (Unit = UnFile, FMT = "(<NJ_Mapping>(I10,1X),'idx_Mat_Extend_from')",        ADVANCE = 'YES', ASYNCHRONOUS = 'NO',  IOSTAT = IO_Write, ERR = 1006 ) ( idx_Mat_Extend_from      ( I ), I = 1_Lng, NJ_Mapping ) ;
    Write (Unit = UnFile, FMT = "(<NJ_Mapping>(I10,1X),'idx_Mat_Extend_to')",          ADVANCE = 'YES', ASYNCHRONOUS = 'NO',  IOSTAT = IO_Write, ERR = 1006 ) ( idx_Mat_Extend_to        ( I ), I = 1_Lng, NJ_Mapping ) ;
    Write (Unit = UnFile, FMT = "(<NStore_Mapping>(I10,1X),'U_Store_Numbers_Global')", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',  IOSTAT = IO_Write, ERR = 1006 ) ( U_Store_Numbers_Global   ( I ), I = 1_Lng, NStore_Mapping ) ;
    Write (Unit = UnFile, FMT = "(<NStore_Rank>(I10,1X),'idx_u_from')",                ADVANCE = 'YES', ASYNCHRONOUS = 'NO',  IOSTAT = IO_Write, ERR = 1006 ) ( idx_u_from               ( I ), I = 1_Lng, NStore_Rank ) ;
    Write (Unit = UnFile, FMT = "(<NStore_Rank>(I10,1X),'idx_u_to')",                  ADVANCE = 'YES', ASYNCHRONOUS = 'NO',  IOSTAT = IO_Write, ERR = 1006 ) ( idx_u_to                 ( I ), I = 1_Lng, NStore_Rank ) ;

    ! - Closing the output file ---------------------------------------------------------------------------------------------------------------------
    UnFile =  Un_Inversion_DS ;
    Close ( Unit = UnFile, Status = 'KEEP', ERR =  1002, IOSTAT = IO_File ) ;

    ! - Say goodbye to your love! -------------------------------------------------------------------------------------------------------------------
    DeAllocate ( idx_Mat_from, idx_Mat_Extend_from, idx_Mat_Extend_to, ID_Local, U_Store_Numbers_Global, idx_u_from, idx_u_to )
!    DeAllocate ( idx_Mat_to )


  End Do


Write(*    ,*) 'End Subroutine < Inversion_DS >' ;
Write(UnInf,*) 'End Subroutine < Inversion_DS >' ;
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
End Subroutine Inversion_DS ;

End Module Inversion_Data_Structure ;
