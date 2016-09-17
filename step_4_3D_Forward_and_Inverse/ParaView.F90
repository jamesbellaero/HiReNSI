
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                                                                                                                  ++
! Start:        01 May 2012                                                                                                                        ++
! Last Update:  26 July 2012                                                                                                                       ++
!                                                                                                                                                  ++
! Description: THIS Module generates the required input files for the visualizer software, ParaView.                                               ++
!              There are subroutines for serial and parallel code.                                                                                 ++
!                                                                                                                                                  ++
!              Subroutines                     Description                                                                                         ++
!              ==============                  ==============                                                                                      ++
!                                                                                                                                                  ++
!                                                                                                                                                  ++
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Module Visualizer ;

! Libraries and Modules
Use HDF5 ;
Use Parameters ;

Implicit None ;


  Interface
!    Module Procedure 
  End Interface    ;


Contains ;


! ***************************************************************************************************************************************************
!                                                                                                                                                  **
! Start:        01 May 2012                                                                                                                        ++
! Last Update:  26 July 2012                                                                                                                       **
! Description: This subroutine generates the input file for the visualizer software, ParaView, from the serial code.                               **
! Called from:                                                                                                                                     **
!                                                                                                                                                  **
! ***************************************************************************************************************************************************

Subroutine ParaView    (                                                                                                        &
NDim,     IStep,     NJ, NEl,                                                                                                   & ! Integer Variables
!                                                                                                                               & ! Real Variables
ELT,     INod, ID,                                                                                                              & ! Integer Arrays
U,     XYZ,                                                                                                                     & ! Real Arrays
Name, OutDir                                                                                                                    & ! Characters
!                                                                                                                               & ! Type
) ;

Implicit None ;

! =========================== Global Variables ======================================================================================================

! - Integer Variables -------------------------------------------------------------------------------------------------------------------------------
Integer , Intent(In)    :: NDim ;
Integer , Intent(In)    :: IStep ;
Integer , Intent(In)    :: NJ, NEL ;

! - Real Variables ----------------------------------------------------------------------------------------------------------------------------------
!#Real (Kind=Dbl), Intent(In)    ::  ;
!#Real (Kind=Dbl), Intent(InOut) ::  ;
!#Real (Kind=Dbl), Intent(OUT)   ::  ;
! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
!#Complex, Intent(In)    ::  ;
!#Complex, Intent(InOut) ::  ;
!#Complex, Intent(OUT)   ::  ;

! - Integer Arrays ----------------------------------------------------------------------------------------------------------------------------------
Integer , Intent(In),    Dimension (:  )  :: ELT ;

Integer , Intent(IN),    Dimension (:,:)  :: INOD, ID ;

! - Real Arrays -------------------------------------------------------------------------------------------------------------------------------------
Real(8), Intent(IN)  , Dimension (:  )  :: U ;
Real(8), Intent(IN)  , Dimension (:,:)  :: XYZ ;

! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
!#Complex, Intent(In)    ::  ;
!#Complex, Intent(InOut) ::  ;
!#Complex, Intent(OUT)   ::  ;
! - Character Variables -----------------------------------------------------------------------------------------------------------------------------
Character (Kind = 1, Len = 30 ) :: Name ;        ! Name of Input file
Character (Kind = 1, Len = 100) :: OutDir ;      ! Directory of output files (Results)

! - Logical Variables -------------------------------------------------------------------------------------------------------------------------------
!#Logical   ::  ;
! - Type DECLERATIONS -------------------------------------------------------------------------------------------------------------------------------
!#Type() :: ;

! =========================== Local Variables =======================================================================================================
! - Integer Variables -------------------------------------------------------------------------------------------------------------------------------
Integer   :: NNode ;                  ! Number of Nodes of element
Integer   :: ElType_PView ;           ! Stores the cell number based on the ParaView standard - see ParaView manual for the list of the numbers
Integer   :: IDim ;                   ! Loop index over the dimension of the model
Integer   :: INode ;                  ! Loop index over the node number of elements
Integer   :: EType ;                  ! Element Type
Integer   :: IO_Write ;               ! Used for IOSTAT - Input Output Status - in the Write cammand
Integer   :: IO_File ;                ! Used for IOSTAT - Input Output Status - in the OPEN cammand
Integer   :: UnFile ;                 ! Holds the Unit of a file for error message
Integer   :: IJ ;                     ! Loop index over the node numbers
Integer   :: IEl ;                    ! Loop index over the element numbers
Integer   :: OffSet ;                 ! OffSet for Cell node numbers 

! - Real Variables ----------------------------------------------------------------------------------------------------------------------------------
Real(8)      :: PointData ;              ! The value assined to a node in ParaView

! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
!#Complex              ::  ;
! - Integer Arrays ----------------------------------------------------------------------------------------------------------------------------------
!#Integer (Kind=Shrt), Dimension (?)  ::  ;
! - Real Arrays -------------------------------------------------------------------------------------------------------------------------------------
!#Real (Kind=Dbl), Dimension (?)      ::  ;
! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
!#Complex  ::  ;
! - Character Variables -----------------------------------------------------------------------------------------------------------------------------
Character (Kind = 1, Len = 22 ) :: Index ;       ! A variable for storing the step number in a Character format for adding at the end of the input file Name
! - Logical Variables -------------------------------------------------------------------------------------------------------------------------------
!#Logical   ::  ;

! =========================== Subroutine CODE =======================================================================================================

! <<<< I M P O R T A N T   N O T I C E >>>>: Node number in PARAVIEW starts from 0.

! open the output file for each time step
write(index, *)IStep ;

UnFile = UN_PView ;
Open ( Unit = UnFile, FILE = TRIM(Name)//'_'//Trim(AdjustL(Index))//'.vtu', ERR =  1001, IOSTAT = IO_File, ACCESS = 'SEQUENTIAL', ACTION = 'Write', ASYNCHRONOUS = 'NO', BLANK = 'NULL', BLOCKSIZE = 0, DEFAULTFILE = TRIM(OutDir), DISPOSE = 'KEEP', FORM = 'FORMATTED', POSITION = 'ASIS', STATUS = 'REPLACE' ) ;

! - header --------------------------------------
Write (Unit = UnFile, FMT = "(A73)", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">' ;
Write (Unit = UnFile, FMT = "('  <UnstructuredGrid>')", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 )  ;
Write (Unit = UnFile, FMT = "(A27,I19,A17,I19,A2)", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) '    <Piece NumberOfPoints="', NJ, '" NumberOfCells="', NEl, '">';
Write (Unit = UnFile, FMT = "('')", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 )  ;

! Points, i.e., nodes in our language -----------
Write (Unit = UnFile, FMT = "('      <Points>')", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 )  ;
Write (Unit = UnFile, FMT = "(A64)", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) '<DataArray type="Float32" NumberOfComponents="3" format="ascii">' ;

  ! coordinates
  Do IJ = 1, NJ ;
    If      (NDim == 2) Then ; Write (Unit = UnFile, FMT = "(3(E26.19,1x))", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) (XYZ (IJ,IDim),IDim = 1,NDim),0.0d0 ;
    Else If (NDim == 3) Then ; Write (Unit = UnFile, FMT = "(3(E26.19,1x))", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) (XYZ (IJ,IDim),IDim = 1,NDim) ;
    End If ;
  End Do ;

Write (Unit = UnFile, FMT = "('        </DataArray>')", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 )  ;
Write (Unit = UnFile, FMT = "('      </Points>')", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 )  ;

! - Cells, i.e., Elements in our language -------
Write (Unit = UnFile, FMT = "('      <Cells>')", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 )  ;
Write (Unit = UnFile, FMT = "(A67)", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) '        <DataArray type="Int32" Name="connectivity" format="ascii">' ;

  Do IEl = 1, NEl ;
    EType = ELT ( IEl ) ;

      IF      ( EType == El2d4NSldPN  ) Then ;      ! 4 noded -2D - Solid - PLANE STRESS
        NNode = 4   ; ElType_PView = 9 ;
      Else IF ( EType == El2d4NPMLPN  ) Then ;      ! 4 noded -2D - PML - PLANE STRESS
        NNode = 4   ; ElType_PView = 9 ;
      Else If ( EType == El2d8NSldPN  ) Then ;      ! 8 noded -2D - Solid - PLANE STRAIN
        NNode = 8   ; ElType_PView = 23 ;
      Else If ( EType == El2d8NPMLPN  ) Then ;      ! 8 noded -2D - PML - PLANE STRAIN
        NNode = 8   ; ElType_PView = 23 ;
      Else If ( EType == El2d6NSldPN  ) Then ;      ! 6 noded -2D - Solid - PLANE STRAIN
        NNode = 6   ; ElType_PView = 22 ;
      Else If ( EType == El2d6NPMLPN  ) Then ;      ! 6 noded -2D - PML - PLANE STRAIN
        NNode = 6   ; ElType_PView = 22 ;
      Else If ( EType == El2d3NSldPN ) Then ;      ! 3 noded -2D - Solid - PLANE STRAIN
        NNode = 3   ; ElType_PView = 5 ;
      Else If ( EType == El2d3NPMLPN ) Then ;      ! 3 noded -2D - PML - PLANE STRAIN
        NNode = 3   ; ElType_PView = 5 ;
      Else ; ! Element Type is not in the list
        Write(*, Fmt_Element2) ;  Write(UnInf, Fmt_Element2) ;
        Write(*, Fmt_FL) ; Write(UnInf, Fmt_FL) ;
        !#Call BEEP_FAIL ;
        Write(*, Fmt_End) ; Read(*,*) ;  STOP ;
      End If ;

    Write (Unit = UnFile, FMT = "(<NNode>(I19,1x))", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) ( INod (INode, IEl)-1, INode = 1, NNode) ; ! 1 is reduced because nodes in ParaView starts form 0.
  End Do ;

Write (Unit = UnFile, FMT = "('        </DataArray>')", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 )  ;
Write (Unit = UnFile, FMT = "(A62)", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) '        <DataArray type="Int32" Name="offsets" format="ascii">' ;


! - OffSet ---------------------------------------------------------
OffSet = 0 ;
  Do IEl = 1, NEl ;
    EType = ELT ( IEl ) ;

      IF      ( EType == El2d4NSldPN  ) Then ;      ! 4 noded -2D - Solid - PLANE STRESS
        NNode = 4   ; ElType_PView = 9 ;
      Else IF ( EType == El2d4NPMLPN  ) Then ;      ! 4 noded -2D - PML - PLANE STRESS
        NNode = 4   ; ElType_PView = 9 ;
      Else If ( EType == El2d8NSldPN  ) Then ;      ! 8 noded -2D - Solid - PLANE STRAIN
        NNode = 8   ; ElType_PView = 23 ;
      Else If ( EType == El2d8NPMLPN  ) Then ;      ! 8 noded -2D - PML - PLANE STRAIN
        NNode = 8   ; ElType_PView = 23 ;
      Else If ( EType == El2d6NSldPN  ) Then ;      ! 6 noded -2D - Solid - PLANE STRAIN
        NNode = 6   ; ElType_PView = 22 ;
      Else If ( EType == El2d6NPMLPN  ) Then ;      ! 6 noded -2D - PML - PLANE STRAIN
        NNode = 6   ; ElType_PView = 22 ;
      Else If ( EType == El2d3NSldPN ) Then ;      ! 3 noded -2D - Solid - PLANE STRAIN
        NNode = 3   ; ElType_PView = 5 ;
      Else If ( EType == El2d3NPMLPN ) Then ;      ! 3 noded -2D - PML - PLANE STRAIN
        NNode = 3   ; ElType_PView = 5 ;
      Else ; ! Element Type is not in the list
        Write(*, Fmt_Element2) ;  Write(UnInf, Fmt_Element2) ;
        Write(*, Fmt_FL) ; Write(UnInf, Fmt_FL) ;
        !#Call BEEP_FAIL ;
        Write(*, Fmt_End) ; Read(*,*) ;  STOP ;
      End If ;

    OffSet = OffSet + NNode ;

    !Write (Unit = UnFile, FMT = "(I19,1X)", ADVANCE = 'No', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) OffSet ;
    Write (Unit = UnFile, FMT = "(I19,1X)", ADVANCE = 'Yes', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) OffSet ;
  End Do ;

!Write (Unit = UnFile, FMT = "()", ADVANCE = 'Yes', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) ;
Write (Unit = UnFile, FMT = "('        </DataArray>')", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 )  ;
Write (Unit = UnFile, FMT = "(A60)", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) '        <DataArray type="UInt8" Name="types" format="ascii">' ;

! - Element Type
  Do IEl = 1, NEl ;
    EType = ELT ( IEl ) ;

      IF      ( EType == El2d4NSldPN  ) Then ;      ! 4 noded -2D - Solid - PLANE STRESS
        NNode = 4   ; ElType_PView = 9 ;
      Else IF ( EType == El2d4NPMLPN  ) Then ;      ! 4 noded -2D - PML - PLANE STRESS
        NNode = 4   ; ElType_PView = 9 ;
      Else If ( EType == El2d8NSldPN  ) Then ;      ! 8 noded -2D - Solid - PLANE STRAIN
        NNode = 8   ; ElType_PView = 23 ;
      Else If ( EType == El2d8NPMLPN  ) Then ;      ! 8 noded -2D - PML - PLANE STRAIN
        NNode = 8   ; ElType_PView = 23 ;
      Else If ( EType == El2d6NSldPN  ) Then ;      ! 6 noded -2D - Solid - PLANE STRAIN
        NNode = 6   ; ElType_PView = 22 ;
      Else If ( EType == El2d6NPMLPN  ) Then ;      ! 6 noded -2D - PML - PLANE STRAIN
        NNode = 6   ; ElType_PView = 22 ;
      Else If ( EType == El2d3NSldPN ) Then ;      ! 3 noded -2D - Solid - PLANE STRAIN
        NNode = 3   ; ElType_PView = 5 ;
      Else If ( EType == El2d3NPMLPN ) Then ;      ! 3 noded -2D - PML - PLANE STRAIN
        NNode = 3   ; ElType_PView = 5 ;
      Else ; ! Element Type is not in the list
        Write(*, Fmt_Element2) ;  Write(UnInf, Fmt_Element2) ;
        Write(*, Fmt_FL) ; Write(UnInf, Fmt_FL) ;
        !#Call BEEP_FAIL ;
        Write(*, Fmt_End) ; Read(*,*) ;  STOP ;
      End If ;

    Write (Unit = UnFile, FMT = "(I3,1X)", ADVANCE = 'No', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) ElType_PView ;
  End Do ;

Write (Unit = UnFile, FMT = "()", ADVANCE = 'Yes', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) ;
Write (Unit = UnFile, FMT = "('        </DataArray>')", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 )  ;
Write (Unit = UnFile, FMT = "('      </Cells>')", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 )  ;


! - Points data ---------------------------------
Write (Unit = UnFile, FMT = "('      <PointData>')", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 )  ;

! Total Displacement
Write (Unit = UnFile, FMT = "(A66)", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) '        <DataArray type="Float32" Name="total dis" format="ascii">' ;

  Do IJ = 1, NJ ;
    PointData = 0.0d0 ;

      Do IDim = 1, NDim ; 
        If ( ID ( IJ , IDim ) /= 0 ) PointData = PointData + U ( ID ( IJ , IDim ) ) * U ( ID ( IJ , IDim ) ) ;
      End Do ;

    PointData = Sqrt ( PointData ) ;
    !Write (Unit = UnFile, FMT = "(E31.23E3,1x)", ADVANCE = 'NO', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) PointData ;
    Write (Unit = UnFile, FMT = "(E31.23E3,1X)", ADVANCE = 'Yes', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) PointData ;
  End Do ;

!Write (Unit = UnFile, FMT = "()", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 )  ;
Write (Unit = UnFile, FMT = "('        </DataArray>')", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 )  ;

! Displacement in X direction
Write (Unit = UnFile, FMT = "(A66)", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) '        <DataArray type="Float32" Name="X dis" format="ascii">' ;

  Do IJ = 1, NJ ;
    PointData = 0.0d0 ;
    If ( ID ( IJ , 1 ) /= 0 ) PointData = U ( ID ( IJ , 1 ) ) ;

    !Write (Unit = UnFile, FMT = "(E31.23E3,1x)", ADVANCE = 'NO', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) PointData ;
    Write (Unit = UnFile, FMT = "(E31.23E3,1X)", ADVANCE = 'Yes', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) PointData ;
  End Do ;

!Write (Unit = UnFile, FMT = "()", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 )  ;
Write (Unit = UnFile, FMT = "('        </DataArray>')", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 )  ;

! Displacement in Y direction
Write (Unit = UnFile, FMT = "(A66)", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) '        <DataArray type="Float32" Name="Y dis" format="ascii">' ;

  Do IJ = 1, NJ ;
    PointData = 0.0d0 ;
    If ( ID ( IJ , 2 ) /= 0 ) PointData = U ( ID ( IJ , 2 ) ) ;

    !Write (Unit = UnFile, FMT = "(E31.23E3,1x)", ADVANCE = 'NO', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) PointData ;
    Write (Unit = UnFile, FMT = "(E31.23E3,1X)", ADVANCE = 'Yes', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) PointData ;
  End Do ;

!Write (Unit = UnFile, FMT = "()", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 )  ;
Write (Unit = UnFile, FMT = "('        </DataArray>')", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 )  ;

Write (Unit = UnFile, FMT = "('      </PointData>')", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 )  ;

! - Cell data -----------------------------------

! - Footer --------------------------------------
Write (Unit = UnFile, FMT = "('    </Piece>')", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 )  ;
Write (Unit = UnFile, FMT = "('  </UnstructuredGrid>')", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 )  ;
Write (Unit = UnFile, FMT = "('</VTKFile>')", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 )  ;

! - Closing the output file ---------------------------------------------------------------------------------------------------------------------
UnFile =  UN_PView ;
Close ( Unit = UnFile, Status = 'KEEP', ERR =  1002, IOSTAT = IO_File ) ;

!#Write(*    ,*) 'End Subroutine < ParaView >' ;
!#Write(UnInf,*) 'End Subroutine < ParaView >' ;
Return ;

! =============================================== OPEN ERRORS =======================================================================================
1001  IF ( IO_File > 0 ) Then ;
        Write(*, Fmt_ERR1_OPEN) UnFile, IO_File  ;  Write(UnInf, Fmt_ERR1_OPEN) UnFile, IO_File  ; 
        Write(*, Fmt_FL) ; Write(UnInf, Fmt_FL) ;
        !#Call BEEP_FAIL ;
        Write(*, Fmt_End) ; Read(*,*) ;  STOP ;
      Else If ( IO_File < 0 ) Then ;
        Write(*, Fmt_ERR1_OPEN) UnFile, IO_File  ; 
        Write(UnInf, Fmt_ERR1_OPEN) UnFile, IO_File  ;  Write(*, Fmt_FL) ; Write(UnInf, Fmt_FL) ;
        !#Call BEEP_FAIL ;
        Write(*, Fmt_End) ; Read(*,*) ;  STOP ;
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
        Write(*, Fmt_FL) ;  Write(UnInf, Fmt_FL) ; Write(*, Fmt_End) ; Read(*,*) ;  Return ;

End Subroutine ParaView ;


! ***************************************************************************************************************************************************
!                                                                                                                                                  **
! Start:        09 Sept 2012                                                                                                                       **
! Last Update:  14 March 2013                                                                                                                      **
! Description: This subroutine generates the input file for the visualizer software, ParaView, from the serial code.                               **
!              This subroutine exclusively designed for 9-noded elements.                                                                          **
! Called from:                                                                                                                                     **
!                                                                                                                                                  **
! ***************************************************************************************************************************************************

Subroutine ParaView_9noded    (                                                                                                 &
NDim,     IStep,     NJ, NEl,                                                                                                   & ! Integer Variables
!                                                                                                                               & ! Real Variables
ELT,     INod, ID,                                                                                                              & ! Integer Arrays
U, UD, UDD,     XYZ,                                                                                                            & ! Real Arrays
Name, OutDir                                                                                                                    & ! Characters
!                                                                                                                               & ! Type
) ;

Implicit None ;

! =========================== Global Variables ======================================================================================================

! - Integer Variables -------------------------------------------------------------------------------------------------------------------------------
Integer , Intent(In)    :: NDim ;
Integer , Intent(In)    :: IStep ;
Integer , Intent(In)    :: NJ, NEL ;

! - Real Variables ----------------------------------------------------------------------------------------------------------------------------------
! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------

! - Integer Arrays ----------------------------------------------------------------------------------------------------------------------------------
Integer , Intent(In),    Dimension (:  )  :: ELT ;

Integer , Intent(IN),    Dimension (:,:)  :: INOD, ID ;

! - Real Arrays -------------------------------------------------------------------------------------------------------------------------------------
Real(8), Intent(IN)  , Dimension (:  )  :: U, UD, UDD ;
Real(8), Intent(IN)  , Dimension (:,:)  :: XYZ ;

! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------

! - Character Variables -----------------------------------------------------------------------------------------------------------------------------
Character (Kind = 1, Len = 30 ) :: Name ;        ! Name of Input file
Character (Kind = 1, Len = 100) :: OutDir ;      ! Directory of output files (Results)

! - Logical Variables -------------------------------------------------------------------------------------------------------------------------------
! - Type DECLERATIONS -------------------------------------------------------------------------------------------------------------------------------

! =========================== Local Variables =======================================================================================================
! - Integer Variables -------------------------------------------------------------------------------------------------------------------------------
Integer   :: NNode ;                  ! Number of Nodes of element
Integer   :: ElType_PView ;           ! Stores the cell number based on the ParaView standard - see ParaView manual for the list of the numbers
Integer   :: IDim ;                   ! Loop index over the dimension of the model
Integer   :: INode ;                  ! Loop index over the node number of elements
Integer   :: EType ;                  ! Element Type
Integer   :: IO_Write ;               ! Used for IOSTAT - Input Output Status - in the Write cammand
Integer   :: IO_File ;                ! Used for IOSTAT - Input Output Status - in the OPEN cammand
Integer   :: UnFile ;                 ! Holds the Unit of a file for error message
Integer   :: IJ ;                     ! Loop index over the node numbers
Integer   :: IEl ;                    ! Loop index over the element numbers
Integer   :: OffSet ;                 ! OffSet for Cell node numbers 

! - Real Variables ----------------------------------------------------------------------------------------------------------------------------------
Real (8)      :: PointData ;              ! The value assined to a node in ParaView

! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
! - Integer Arrays ----------------------------------------------------------------------------------------------------------------------------------
! - Real Arrays -------------------------------------------------------------------------------------------------------------------------------------
! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------

! - Character Variables -----------------------------------------------------------------------------------------------------------------------------
Character (Kind = 1, Len = 22 ) :: Index ;       ! A variable for storing the step number in a Character format for adding at the end of the input file Name

! - Logical Variables -------------------------------------------------------------------------------------------------------------------------------

! =========================== Subroutine CODE =======================================================================================================

! <<<< I M P O R T A N T   N O T I C E >>>>: Node number in PARAVIEW starts from 0.

! open the output file for each time step
write(index, *)IStep ;

UnFile = UN_PView ;
Open ( Unit = UnFile, FILE = TRIM(Name)//'_'//Trim(AdjustL(Index))//'.vtu', ERR =  1001, IOSTAT = IO_File, ACCESS = 'SEQUENTIAL', ACTION = 'Write', ASYNCHRONOUS = 'NO', BLANK = 'NULL', BLOCKSIZE = 0, DEFAULTFILE = TRIM(OutDir), DISPOSE = 'KEEP', FORM = 'FORMATTED', POSITION = 'ASIS', STATUS = 'REPLACE' ) ;

! - header --------------------------------------
Write (Unit = UnFile, FMT = "(A73)", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">' ;
Write (Unit = UnFile, FMT = "('  <UnstructuredGrid>')", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 )  ;
Write (Unit = UnFile, FMT = "(A27,I19,A17,I19,A2)", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) '    <Piece NumberOfPoints="', NJ-NEl, '" NumberOfCells="', NEl, '">';
Write (Unit = UnFile, FMT = "('')", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 )  ;

! Points, i.e., nodes in our language -----------
Write (Unit = UnFile, FMT = "('      <Points>')", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 )  ;
Write (Unit = UnFile, FMT = "(A64)", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) '<DataArray type="Float32" NumberOfComponents="3" format="ascii">' ;

  ! coordinates
  Do IJ = 1, NJ - NEl ;
    If      (NDim == 2) Then ; Write (Unit = UnFile, FMT = "(3(E26.19,1x))", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) (XYZ (IJ,IDim),IDim = 1,NDim),0.0d0 ;
    Else If (NDim == 3) Then ; Write (Unit = UnFile, FMT = "(3(E26.19,1x))", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) (XYZ (IJ,IDim),IDim = 1,NDim) ;
    End If ;
  End Do ;

Write (Unit = UnFile, FMT = "('        </DataArray>')", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 )  ;
Write (Unit = UnFile, FMT = "('      </Points>')", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 )  ;

! - Cells, i.e., Elements in our language -------
Write (Unit = UnFile, FMT = "('      <Cells>')", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 )  ;
Write (Unit = UnFile, FMT = "(A67)", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) '        <DataArray type="Int32" Name="connectivity" format="ascii">' ;

  Do IEl = 1, NEl ;
    EType = ELT ( IEl ) ;

      IF      ( EType == El2d4NSldPN  ) Then ;      ! 4 noded -2D - Solid - PLANE STRESS
        NNode = 4   ; ElType_PView = 9 ;
      Else IF ( EType == El2d4NPMLPN  ) Then ;      ! 4 noded -2D - PML - PLANE STRESS
        NNode = 4   ; ElType_PView = 9 ;
      Else If ( EType == El2d8NSldPN  ) Then ;      ! 8 noded -2D - Solid - PLANE STRAIN
        NNode = 8   ; ElType_PView = 23 ;
      Else If ( EType == El2d8NPMLPN  ) Then ;      ! 8 noded -2D - PML - PLANE STRAIN
        NNode = 8   ; ElType_PView = 23 ;
      Else If ( EType == El2d6NSldPN  ) Then ;      ! 6 noded -2D - Solid - PLANE STRAIN
        NNode = 6   ; ElType_PView = 22 ;
      Else If ( EType == El2d6NPMLPN  ) Then ;      ! 6 noded -2D - PML - PLANE STRAIN
        NNode = 6   ; ElType_PView = 22 ;
      Else If ( EType == El2d3NSldPN  ) Then ;       ! 3 noded -2D - Solid - PLANE STRAIN
        NNode = 3   ; ElType_PView = 5 ;
      Else If ( EType == El2d3NPMLPN  ) Then ;       ! 3 noded -2D - PML - PLANE STRAIN
        NNode = 3   ; ElType_PView = 5 ;
      Else If ( EType == SEl2d9NSldPN ) Then ;      ! 9 noded -2D - Solid - PLANE STRAIN - Spectral Element
        NNode = 8   ; ElType_PView = 23 ;
      Else If ( EType == SEl2d9NPMLPN ) Then ;      ! 9 noded -2D - PML - PLANE STRAIN - Spectral Element
        NNode = 8   ; ElType_PView = 23 ;
      Else If ( EType == SEl2d9NMPMLPN) Then ;     ! 9 noded -2D - MPML - PLANE STRAIN - Spectral Element
        NNode = 8   ; ElType_PView = 23 ;
      Else If ( EType == SEl2d7NSldPN ) Then ;     ! 7 noded -2D - Solid - PLANE STRAIN - Spectral Element
        NNode = 6   ; ElType_PView = 22 ;
      Else If ( EType == SEl2d7NMPMLPN) Then ;     ! 7 noded -2D - MPML - PLANE STRAIN - Spectral Element
        NNode = 6   ; ElType_PView = 22 ;
      Else If ( EType == SEl2d9NSldSH ) Then ;     ! 9 noded -2D - Solid - PLANE STRAIN - Spectral Element - SH wave
        NNode = 8   ; ElType_PView = 23 ;
      Else If ( EType == SEl2d9NPMLSH ) Then ;     ! 9 noded -2D - MPML - PLANE STRAIN - Spectral Element - SH wave
        NNode = 8   ; ElType_PView = 23 ;
      Else ; ! Element Type is not in the list
        Write(*, Fmt_Element2) ;  Write(UnInf, Fmt_Element2) ;
        Write(*, Fmt_FL) ; Write(UnInf, Fmt_FL) ;
        !#Call BEEP_FAIL ;
        Write(*, Fmt_End) ; Read(*,*) ;  STOP ;
      End If ;

    Write (Unit = UnFile, FMT = "(<NNode>(I19,1x))", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) ( INod (INode, IEl)-1, INode = 1, NNode) ; ! 1 is reduced because nodes in ParaView starts form 0.
  End Do ;

Write (Unit = UnFile, FMT = "('        </DataArray>')", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 )  ;
Write (Unit = UnFile, FMT = "(A62)", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) '        <DataArray type="Int32" Name="offsets" format="ascii">' ;


! - OffSet ---------------------------------------------------------
OffSet = 0 ;
  Do IEl = 1, NEl ;
    EType = ELT ( IEl ) ;

      IF      ( EType == El2d4NSldPN  ) Then ;      ! 4 noded -2D - Solid - PLANE STRESS
        NNode = 4   ; ElType_PView = 9 ;
      Else IF ( EType == El2d4NPMLPN  ) Then ;      ! 4 noded -2D - PML - PLANE STRESS
        NNode = 4   ; ElType_PView = 9 ;
      Else If ( EType == El2d8NSldPN  ) Then ;      ! 8 noded -2D - Solid - PLANE STRAIN
        NNode = 8   ; ElType_PView = 23 ;
      Else If ( EType == El2d8NPMLPN  ) Then ;      ! 8 noded -2D - PML - PLANE STRAIN
        NNode = 8   ; ElType_PView = 23 ;
      Else If ( EType == El2d6NSldPN  ) Then ;      ! 6 noded -2D - Solid - PLANE STRAIN
        NNode = 6   ; ElType_PView = 22 ;
      Else If ( EType == El2d6NPMLPN  ) Then ;      ! 6 noded -2D - PML - PLANE STRAIN
        NNode = 6   ; ElType_PView = 22 ;
      Else If ( EType == El2d3NSldPN ) Then ;      ! 3 noded -2D - Solid - PLANE STRAIN
        NNode = 3   ; ElType_PView = 5 ;
      Else If ( EType == El2d3NPMLPN ) Then ;      ! 3 noded -2D - PML - PLANE STRAIN
        NNode = 3   ; ElType_PView = 5 ;
      Else If ( EType == SEl2d9NSldPN ) Then ;      ! 9 noded -2D - Solid - PLANE STRAIN - Spectral Element
        NNode = 8   ; ElType_PView = 23 ;
      Else If ( EType == SEl2d9NPMLPN ) Then ;      ! 9 noded -2D - PML - PLANE STRAIN - Spectral Element
        NNode = 8   ; ElType_PView = 23 ;
      Else If ( EType == SEl2d9NMPMLPN ) Then ;     ! 9 noded -2D - MPML - PLANE STRAIN - Spectral Element
        NNode = 8   ; ElType_PView = 23 ;
      Else If ( EType == SEl2d7NSldPN ) Then ;     ! 7 noded -2D - Solid - PLANE STRAIN - Spectral Element
        NNode = 6   ; ElType_PView = 22 ;
      Else If ( EType == SEl2d7NMPMLPN) Then ;     ! 7 noded -2D - MPML - PLANE STRAIN - Spectral Element
        NNode = 6   ; ElType_PView = 22 ;
      Else If ( EType == SEl2d9NSldSH ) Then ;     ! 9 noded -2D - Solid - PLANE STRAIN - Spectral Element - SH wave
        NNode = 8   ; ElType_PView = 23 ;
      Else If ( EType == SEl2d9NPMLSH ) Then ;     ! 9 noded -2D - MPML - PLANE STRAIN - Spectral Element - SH wave
        NNode = 8   ; ElType_PView = 23 ;
      Else ; ! Element Type is not in the list
        Write(*, Fmt_Element2) ;  Write(UnInf, Fmt_Element2) ;
        Write(*, Fmt_FL) ; Write(UnInf, Fmt_FL) ;
        !#Call BEEP_FAIL ;
        Write(*, Fmt_End) ; Read(*,*) ;  STOP ;
      End If ;

    OffSet = OffSet + NNode ;

    !Write (Unit = UnFile, FMT = "(I19,1X)", ADVANCE = 'No', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) OffSet ;
    Write (Unit = UnFile, FMT = "(I19,1X)", ADVANCE = 'Yes', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) OffSet ;
  End Do ;

!Write (Unit = UnFile, FMT = "()", ADVANCE = 'Yes', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) ;
Write (Unit = UnFile, FMT = "('        </DataArray>')", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 )  ;
Write (Unit = UnFile, FMT = "(A60)", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) '        <DataArray type="UInt8" Name="types" format="ascii">' ;

! - Element Type
  Do IEl = 1, NEl ;
    EType = ELT ( IEl ) ;

      IF      ( EType == El2d4NSldPN  ) Then ;      ! 4 noded -2D - Solid - PLANE STRESS
        NNode = 4   ; ElType_PView = 9 ;
      Else IF ( EType == El2d4NPMLPN  ) Then ;      ! 4 noded -2D - PML - PLANE STRESS
        NNode = 4   ; ElType_PView = 9 ;
      Else If ( EType == El2d8NSldPN  ) Then ;      ! 8 noded -2D - Solid - PLANE STRAIN
        NNode = 8   ; ElType_PView = 23 ;
      Else If ( EType == El2d8NPMLPN  ) Then ;      ! 8 noded -2D - PML - PLANE STRAIN
        NNode = 8   ; ElType_PView = 23 ;
      Else If ( EType == El2d6NSldPN  ) Then ;      ! 6 noded -2D - Solid - PLANE STRAIN
        NNode = 6   ; ElType_PView = 22 ;
      Else If ( EType == El2d6NPMLPN  ) Then ;      ! 6 noded -2D - PML - PLANE STRAIN
        NNode = 6   ; ElType_PView = 22 ;
      Else If ( EType == El2d3NSldPN ) Then ;      ! 3 noded -2D - Solid - PLANE STRAIN
        NNode = 3   ; ElType_PView = 5 ;
      Else If ( EType == El2d3NPMLPN ) Then ;      ! 3 noded -2D - PML - PLANE STRAIN
        NNode = 3   ; ElType_PView = 5 ;
      Else If ( EType == SEl2d9NSldPN ) Then ;      ! 9 noded -2D - Solid - PLANE STRAIN - Spectral Element
        NNode = 8   ; ElType_PView = 23 ;
      Else If ( EType == SEl2d9NPMLPN ) Then ;      ! 9 noded -2D - PML - PLANE STRAIN - Spectral Element
        NNode = 8   ; ElType_PView = 23 ;
      Else If ( EType == SEl2d9NMPMLPN ) Then ;     ! 9 noded -2D - MPML - PLANE STRAIN - Spectral Element
        NNode = 8   ; ElType_PView = 23 ;
      Else If ( EType == SEl2d7NSldPN ) Then ;     ! 7 noded -2D - Solid - PLANE STRAIN - Spectral Element
        NNode = 6   ; ElType_PView = 22 ;
      Else If ( EType == SEl2d7NMPMLPN) Then ;     ! 7 noded -2D - MPML - PLANE STRAIN - Spectral Element
        NNode = 6   ; ElType_PView = 22 ;
      Else If ( EType == SEl2d9NSldSH ) Then ;     ! 9 noded -2D - Solid - PLANE STRAIN - Spectral Element - SH wave
        NNode = 8   ; ElType_PView = 23 ;
      Else If ( EType == SEl2d9NPMLSH ) Then ;     ! 9 noded -2D - MPML - PLANE STRAIN - Spectral Element - SH wave
        NNode = 8   ; ElType_PView = 23 ;
      Else ; ! Element Type is not in the list
        Write(*, Fmt_Element2) ;  Write(UnInf, Fmt_Element2) ;
        Write(*, Fmt_FL) ; Write(UnInf, Fmt_FL) ;
        !#Call BEEP_FAIL ;
        Write(*, Fmt_End) ; Read(*,*) ;  STOP ;
      End If ;

    Write (Unit = UnFile, FMT = "(I3,1X)", ADVANCE = 'No', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) ElType_PView ;
  End Do ;

Write (Unit = UnFile, FMT = "()", ADVANCE = 'Yes', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) ;
Write (Unit = UnFile, FMT = "('        </DataArray>')", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 )  ;
Write (Unit = UnFile, FMT = "('      </Cells>')", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 )  ;

! - Points data ---------------------------------
Write (Unit = UnFile, FMT = "('      <PointData>')", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 )  ;

! Total Displacement
Write (Unit = UnFile, FMT = "(A66)", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) '        <DataArray type="Float32" Name="total dis" format="ascii">' ;

  Do IJ = 1, NJ - NEl ;
    PointData = 0.0d0 ;

      Do IDim = 1, NDim ; 
        If ( ID ( IJ , IDim ) /= 0 ) PointData = PointData + U ( ID ( IJ , IDim ) ) * U ( ID ( IJ , IDim ) ) ;
      End Do ;

    PointData = Sqrt ( PointData ) ;
    !Write (Unit = UnFile, FMT = "(E31.23E3,1x)", ADVANCE = 'NO', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) PointData ;
    Write (Unit = UnFile, FMT = "(E31.23E3,1X)", ADVANCE = 'Yes', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) PointData ;
  End Do ;

!Write (Unit = UnFile, FMT = "()", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 )  ;
Write (Unit = UnFile, FMT = "('        </DataArray>')", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 )  ;


! Displacement in X direction
Write (Unit = UnFile, FMT = "(A66)", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) '        <DataArray type="Float32" Name="X dir" format="ascii">' ;

  Do IJ = 1, NJ - NEl;
    PointData = 0.0d0 ;
    If ( ID ( IJ , 1 ) /= 0 ) PointData = U ( ID ( IJ , 1 ) ) ;

    !Write (Unit = UnFile, FMT = "(E31.23E3,1x)", ADVANCE = 'NO', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) PointData ;
    Write (Unit = UnFile, FMT = "(E18.10E3,1X)", ADVANCE = 'Yes', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) PointData ;
  End Do ;

!Write (Unit = UnFile, FMT = "()", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 )  ;
Write (Unit = UnFile, FMT = "('        </DataArray>')", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 )  ;

! Displacement in Y direction
Write (Unit = UnFile, FMT = "(A66)", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) '        <DataArray type="Float32" Name="Y dir" format="ascii">' ;

  Do IJ = 1, NJ - NEl ;
    PointData = 0.0d0 ;
    If ( ID ( IJ , 2 ) /= 0 ) PointData = U ( ID ( IJ , 2 ) ) ;

    !Write (Unit = UnFile, FMT = "(E31.23E3,1x)", ADVANCE = 'NO', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) PointData ;
    Write (Unit = UnFile, FMT = "(E18.10E3,1X)", ADVANCE = 'Yes', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) PointData ;
  End Do ;

!Write (Unit = UnFile, FMT = "()", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 )  ;
Write (Unit = UnFile, FMT = "('        </DataArray>')", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 )  ;

! Total Velocity
Write (Unit = UnFile, FMT = "(A66)", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) '        <DataArray type="Float32" Name="total Vel" format="ascii">' ;

  Do IJ = 1, NJ - NEl ;
    PointData = 0.0d0 ;

      Do IDim = 1, NDim ; 
        If ( ID ( IJ , IDim ) /= 0 ) PointData = PointData + UD ( ID ( IJ , IDim ) ) * UD ( ID ( IJ , IDim ) ) ;
      End Do ;

    PointData = Sqrt ( PointData ) ;
    !Write (Unit = UnFile, FMT = "(E31.23E3,1x)", ADVANCE = 'NO', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) PointData ;
    Write (Unit = UnFile, FMT = "(E18.10E3,1X)", ADVANCE = 'Yes', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) PointData ;
  End Do ;

!Write (Unit = UnFile, FMT = "()", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 )  ;
Write (Unit = UnFile, FMT = "('        </DataArray>')", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 )  ;

! Velocity in X direction
Write (Unit = UnFile, FMT = "(A66)", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) '        <DataArray type="Float32" Name="X vel" format="ascii">' ;

  Do IJ = 1, NJ - NEl ;
    PointData = 0.0d0 ;
    If ( ID ( IJ , 1 ) /= 0 ) PointData = UD ( ID ( IJ , 1 ) ) ;

    !Write (Unit = UnFile, FMT = "(E31.23E3,1x)", ADVANCE = 'NO', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) PointData ;
    Write (Unit = UnFile, FMT = "(E18.10E3,1X)", ADVANCE = 'Yes', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) PointData ;
  End Do ;

!Write (Unit = UnFile, FMT = "()", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 )  ;
Write (Unit = UnFile, FMT = "('        </DataArray>')", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 )  ;

! Velocity in Y direction
Write (Unit = UnFile, FMT = "(A66)", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) '        <DataArray type="Float32" Name="Y Vel" format="ascii">' ;

  Do IJ = 1, NJ - NEl ;
    PointData = 0.0d0 ;
    If ( ID ( IJ , 2 ) /= 0 ) PointData = UD ( ID ( IJ , 2 ) ) ;
    !Write (Unit = UnFile, FMT = "(E31.23E3,1x)", ADVANCE = 'NO', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) PointData ;
    Write (Unit = UnFile, FMT = "(E18.10E3,1X)", ADVANCE = 'Yes', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) PointData ;
  End Do ;
!Write (Unit = UnFile, FMT = "()", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 )  ;
Write (Unit = UnFile, FMT = "('        </DataArray>')", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 )  ;

! Total Acceleration
Write (Unit = UnFile, FMT = "(A66)", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) '        <DataArray type="Float32" Name="total Acc" format="ascii">' ;

  Do IJ = 1, NJ - NEl ;
    PointData = 0.0d0 ;

      Do IDim = 1, NDim ; 
        If ( ID ( IJ , IDim ) /= 0 ) PointData = PointData + UDD ( ID ( IJ , IDim ) ) * UDD ( ID ( IJ , IDim ) ) ;
      End Do ;

    PointData = Sqrt ( PointData ) ;
    !Write (Unit = UnFile, FMT = "(E31.23E3,1x)", ADVANCE = 'NO', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) PointData ;
    Write (Unit = UnFile, FMT = "(E18.10E3,1X)", ADVANCE = 'Yes', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) PointData ;
  End Do ;

!Write (Unit = UnFile, FMT = "()", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 )  ;
Write (Unit = UnFile, FMT = "('        </DataArray>')", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 )  ;

! Acceleration in the X direction
Write (Unit = UnFile, FMT = "(A66)", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) '        <DataArray type="Float32" Name="X Acc" format="ascii">' ;

  Do IJ = 1, NJ - NEl ;
    PointData = 0.0d0 ;
    If ( ID ( IJ , 1 ) /= 0 ) PointData = UDD ( ID ( IJ , 1 ) ) ;

    !Write (Unit = UnFile, FMT = "(E31.23E3,1x)", ADVANCE = 'NO', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) PointData ;
    Write (Unit = UnFile, FMT = "(E18.10E3,1X)", ADVANCE = 'Yes', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) PointData ;
  End Do ;

!Write (Unit = UnFile, FMT = "()", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 )  ;
Write (Unit = UnFile, FMT = "('        </DataArray>')", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 )  ;

! Acceleration in the Y direction
Write (Unit = UnFile, FMT = "(A66)", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) '        <DataArray type="Float32" Name="Y Acc" format="ascii">' ;

  Do IJ = 1, NJ - NEl ;
    PointData = 0.0d0 ;
    If ( ID ( IJ , 2 ) /= 0 ) PointData = UDD ( ID ( IJ , 2 ) ) ;
    !Write (Unit = UnFile, FMT = "(E31.23E3,1x)", ADVANCE = 'NO', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) PointData ;
    Write (Unit = UnFile, FMT = "(E18.10E3,1X)", ADVANCE = 'Yes', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) PointData ;
  End Do ;
!Write (Unit = UnFile, FMT = "()", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 )  ;
Write (Unit = UnFile, FMT = "('        </DataArray>')", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 )  ;

Write (Unit = UnFile, FMT = "('      </PointData>')", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 )  ;

! - Cell data -----------------------------------

! - Footer --------------------------------------
Write (Unit = UnFile, FMT = "('    </Piece>')", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 )  ;
Write (Unit = UnFile, FMT = "('  </UnstructuredGrid>')", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 )  ;
Write (Unit = UnFile, FMT = "('</VTKFile>')", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 )  ;

! - Closing the output file ---------------------------------------------------------------------------------------------------------------------
UnFile =  UN_PView ;
Close ( Unit = UnFile, Status = 'KEEP', ERR =  1002, IOSTAT = IO_File ) ;

!#Write(*    ,*) 'End Subroutine < ParaView_9noded >' ;
!#Write(UnInf,*) 'End Subroutine < ParaView_9noded >' ;
Return ;

! =============================================== OPEN ERRORS =======================================================================================
1001  IF ( IO_File > 0 ) Then ;
        Write(*, Fmt_ERR1_OPEN) UnFile, IO_File  ;  Write(UnInf, Fmt_ERR1_OPEN) UnFile, IO_File  ; 
        Write(*, Fmt_FL) ; Write(UnInf, Fmt_FL) ;
        !#Call BEEP_FAIL ;
        Write(*, Fmt_End) ; Read(*,*) ;  STOP ;
      Else If ( IO_File < 0 ) Then ;
        Write(*, Fmt_ERR1_OPEN) UnFile, IO_File  ; 
        Write(UnInf, Fmt_ERR1_OPEN) UnFile, IO_File  ;  Write(*, Fmt_FL) ; Write(UnInf, Fmt_FL) ;
        !#Call BEEP_FAIL ;
        Write(*, Fmt_End) ; Read(*,*) ;  STOP ;
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
        Write(*, Fmt_FL) ;  Write(UnInf, Fmt_FL) ; Write(*, Fmt_End) ; Read(*,*) ;  Return ;

End Subroutine ParaView_9noded ;


! ***************************************************************************************************************************************************
!                                                                                                                                                  **
! Start:        09 Nov 2012                                                                                                                        **
! Last Update:  09 Nov 2012                                                                                                                        **
! Description: This subroutine generates the input file for the visualizer software, ParaView, from the serial code.                               **
!              This subroutine exclusively designed for 27-noded elements.                                                                         **
! Called from:                                                                                                                                     **
!                                                                                                                                                  **
! ***************************************************************************************************************************************************

Subroutine ParaView_27noded   (                                                                                                 &
NDim, NNodeDown,    IStep,     NJ, NEl, NJ_Para,                                                                                & ! Integer Variables
!                                                                                                                               & ! Real Variables
ELT,     INod, ID, ParaNodeNum,                                                                                                 & ! Integer Arrays
U, UD,     XYZ,                                                                                                                 & ! Real Arrays
Name, OutDir                                                                                                                    & ! Characters
!                                                                                                                               & ! Type
) ;

Implicit None ;

! =========================== Global Variables ======================================================================================================

! - Integer Variables -------------------------------------------------------------------------------------------------------------------------------
Integer , Intent(In)    :: NDim ;
Integer , Intent(In)    :: IStep ;
Integer , Intent(In)    :: NJ, NEL, NJ_Para ;
Integer , Intent(In)    :: NNodeDown ;


! - Real Variables ----------------------------------------------------------------------------------------------------------------------------------
! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------

! - Integer Arrays ----------------------------------------------------------------------------------------------------------------------------------
Integer , Intent(In),    Dimension (:  )  :: ELT ;
Integer , Intent(In),    Dimension (:  )  :: ParaNodeNum ;
Integer , Intent(IN),    Dimension (:,:)  :: INOD, ID ;

! - Real Arrays -------------------------------------------------------------------------------------------------------------------------------------
Real(8), Intent(IN)  , Dimension (:  )  :: U, UD ;
Real(8), Intent(IN)  , Dimension (:,:)  :: XYZ ;

! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
! - Character Variables -----------------------------------------------------------------------------------------------------------------------------
Character (Kind = 1, Len = 30 ) :: Name ;        ! Name of Input file
Character (Kind = 1, Len = 100) :: OutDir ;      ! Directory of output files (Results)

! - Logical Variables -------------------------------------------------------------------------------------------------------------------------------
!#Logical   ::  ;
! - Type DECLERATIONS -------------------------------------------------------------------------------------------------------------------------------
!#Type() :: ;

! =========================== Local Variables =======================================================================================================
! - Integer Variables -------------------------------------------------------------------------------------------------------------------------------
Integer   :: NNode ;                  ! Number of Nodes of element
Integer   :: ElType_PView ;           ! Stores the cell number based on the ParaView standard - see ParaView manual for the list of the numbers
Integer   :: IDim ;                   ! Loop index over the dimension of the model
Integer   :: INode ;                  ! Loop index over the node number of elements
Integer   :: EType ;                  ! Element Type
Integer   :: IO_Write ;               ! Used for IOSTAT - Input Output Status - in the Write cammand
Integer   :: IO_File ;                ! Used for IOSTAT - Input Output Status - in the OPEN cammand
Integer   :: UnFile ;                 ! Holds the Unit of a file for error message
Integer   :: IJ ;                     ! Loop index over the node numbers
Integer   :: IEl ;                    ! Loop index over the element numbers
Integer   :: OffSet ;                 ! OffSet for Cell node numbers 

! - Real Variables ----------------------------------------------------------------------------------------------------------------------------------
Real(8)      :: PointData ;              ! The value assined to a node in ParaView

! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
! - Integer Arrays ----------------------------------------------------------------------------------------------------------------------------------
! - Real Arrays -------------------------------------------------------------------------------------------------------------------------------------
! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
! - Character Variables -----------------------------------------------------------------------------------------------------------------------------
Character (Kind = 1, Len = 23 ) :: Index ;       ! A variable for storing the step number in a Character format for adding at the end of the input file Name

! - Logical Variables -------------------------------------------------------------------------------------------------------------------------------

! =========================== Subroutine CODE =======================================================================================================

! <<<< I M P O R T A N T   N O T I C E >>>>: Node number in PARAVIEW starts from 0.

!Write(*,*)"Paraview 27node"

! open the output file for each time step
Write(Index, *) IStep ;

UnFile = UN_PView ;
Open ( Unit = UnFile, FILE = TRIM(Name)//'_'//Trim(AdjustL(Index))//'.vtu', ERR =  1001, IOSTAT = IO_File, ACCESS = 'SEQUENTIAL', ACTION = 'Write', ASYNCHRONOUS = 'NO', BLANK = 'NULL', BLOCKSIZE = 0, DEFAULTFILE = TRIM(OutDir), DISPOSE = 'KEEP', FORM = 'FORMATTED', POSITION = 'ASIS', STATUS = 'REPLACE' ) ;

! - header --------------------------------------
Write (Unit = UnFile, FMT = "(A73)", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">' ;
Write (Unit = UnFile, FMT = "('  <UnstructuredGrid>')", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 )  ;
Write (Unit = UnFile, FMT = "(A27,I19,A17,I19,A2)", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) '    <Piece NumberOfPoints="', NJ_Para, '" NumberOfCells="', NEl, '">';
Write (Unit = UnFile, FMT = "('')", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 )  ;

! Points, i.e., nodes in our language -----------
Write (Unit = UnFile, FMT = "('      <Points>')", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 )  ;
Write (Unit = UnFile, FMT = "(A64)", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) '<DataArray type="Float32" NumberOfComponents="3" format="ascii">' ;

  ! coordinates
  Do IJ = 1, NJ_Para ;
    !If ( ParaNodeNum ( IJ ) /= 0_Lng ) Write (Unit = UnFile, FMT = "(3(E25.16,1x))", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) (XYZ (IJ,IDim),IDim = 1,NDim) ;
    Write (Unit = UnFile, FMT = "(3(E25.16,1x))", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) (XYZ (IJ,IDim),IDim = 1,NDim) ;
  End Do ;

Write (Unit = UnFile, FMT = "('        </DataArray>')", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 )  ;
Write (Unit = UnFile, FMT = "('      </Points>')", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 )  ;

! - Cells, i.e., Elements in our language ------- (connectivity)
Write (Unit = UnFile, FMT = "('      <Cells>')", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 )  ;
Write (Unit = UnFile, FMT = "(A67)", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) '        <DataArray type="Int32" Name="connectivity" format="ascii">' ;

NNode = NNodeDown ; 
  If ( NNodeDown == 20 ) Then ;
    ElType_PView = 25 ;
  Else If ( NNodeDown == 8 ) Then ; 
    ElType_PView = 12 ;
  End If ;

  Do IEl = 1, NEl ;
    !Write (Unit = UnFile, FMT = "(<NNode>(I19,1x))", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) ( ParaNodeNum(INod (INode, IEl)) -1_Lng, INode = 1, NNode) ; ! 1 is reduced because nodes in ParaView starts form 0.
    Write (Unit = UnFile, FMT = "(<NNode>(I19,1x))", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) ( INod (INode, IEl) -1, INode = 1, 12), ( INod (INode, IEl) -1, INode = 17, 20), ( INod (INode, IEl) -1, INode = 13, 16); ! 1 is reduced because nodes in ParaView starts form 0.
  End Do ;

Write (Unit = UnFile, FMT = "('        </DataArray>')", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 )  ;
Write (Unit = UnFile, FMT = "(A62)", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) '        <DataArray type="Int32" Name="offsets" format="ascii">' ;

! - OffSet ---------------------------------------------------------
OffSet = 0 ;
  Do IEl = 1, NEl ;
    OffSet = OffSet + NNode ;
    !Write (Unit = UnFile, FMT = "(I19,1X)", ADVANCE = 'No', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) OffSet ;
    Write (Unit = UnFile, FMT = "(I19,1X)", ADVANCE = 'Yes', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) OffSet ;
  End Do ;

!Write (Unit = UnFile, FMT = "()", ADVANCE = 'Yes', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) ;
Write (Unit = UnFile, FMT = "('        </DataArray>')", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 )  ;
Write (Unit = UnFile, FMT = "(A60)", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) '        <DataArray type="UInt8" Name="types" format="ascii">' ;

! - Element Type
  Do IEl = 1, NEl ;
    Write (Unit = UnFile, FMT = "(I3,1X)", ADVANCE = 'No', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) ElType_PView ;
  End Do ;

Write (Unit = UnFile, FMT = "()", ADVANCE = 'Yes', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) ;
Write (Unit = UnFile, FMT = "('        </DataArray>')", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 )  ;
Write (Unit = UnFile, FMT = "('      </Cells>')", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 )  ;

! - Points data ---------------------------------
Write (Unit = UnFile, FMT = "('      <PointData>')", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 )  ;

! Total Displacement
Write (Unit = UnFile, FMT = "(A66)", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) '        <DataArray type="Float32" Name="total dis" format="ascii">' ;

  Do IJ = 1, NJ_para ;
    !If ( ParaNodeNum ( IJ ) == 0_Lng ) cycle ;
    PointData = 0.0d0 ;

      Do IDim = 1, NDim ; 
        If ( ID ( IJ , IDim ) /= 0 ) PointData = PointData + U ( ID ( IJ , IDim ) ) * U ( ID ( IJ , IDim ) ) ;
      End Do ;

    PointData = Sqrt ( PointData ) ;
    !Write (Unit = UnFile, FMT = "(E31.23E3,1x)", ADVANCE = 'NO', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) PointData ;
    Write (Unit = UnFile, FMT = "(E31.23E3,1X)", ADVANCE = 'Yes', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) PointData ;
  End Do ;

!Write (Unit = UnFile, FMT = "()", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 )  ;
Write (Unit = UnFile, FMT = "('        </DataArray>')", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 )  ;

!! Displacement in X direction
!Write (Unit = UnFile, FMT = "(A66)", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) '        <DataArray type="Float32" Name="X dis" format="ascii">' ;
!
!  Do IJ = 1, NJ_Para ;
!    !If ( ParaNodeNum ( IJ ) == 0_Lng ) cycle ;
!    PointData = 0.0_Dbl ;
!    If ( ID ( IJ , 1 ) /= 0_Lng ) PointData = U ( ID ( IJ , 1 ) ) ;
!
!    !Write (Unit = UnFile, FMT = "(E31.23E3,1x)", ADVANCE = 'NO', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) PointData ;
!    Write (Unit = UnFile, FMT = "(E31.23E3,1X)", ADVANCE = 'Yes', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) PointData ;
!  End Do ;
!
!!Write (Unit = UnFile, FMT = "()", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 )  ;
!Write (Unit = UnFile, FMT = "('        </DataArray>')", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 )  ;
!
!! Displacement in Y direction
!Write (Unit = UnFile, FMT = "(A66)", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) '        <DataArray type="Float32" Name="Y dis" format="ascii">' ;
!
!  Do IJ = 1, NJ_Para ;
!    !If ( ParaNodeNum ( IJ ) == 0_Lng ) cycle ;
!    PointData = 0.0_Dbl ;
!    If ( ID ( IJ , 2 ) /= 0_Lng ) PointData = U ( ID ( IJ , 2 ) ) ;
!
!    !Write (Unit = UnFile, FMT = "(E31.23E3,1x)", ADVANCE = 'NO', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) PointData ;
!    Write (Unit = UnFile, FMT = "(E31.23E3,1X)", ADVANCE = 'Yes', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) PointData ;
!  End Do ;
!
!!Write (Unit = UnFile, FMT = "()", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 )  ;
!Write (Unit = UnFile, FMT = "('        </DataArray>')", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 )  ;

!! Total Velocity
!Write (Unit = UnFile, FMT = "(A66)", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) '        <DataArray type="Float32" Name="total Vel" format="ascii">' ;
!
!  Do IJ = 1, NJ_Para ;
!    If ( ParaNodeNum ( IJ ) == 0_Lng ) cycle ;
!    PointData = 0.0_Dbl ;
!
!      Do IDim = 1, NDim ; 
!        If ( ID ( IJ , IDim ) /= 0_Lng ) PointData = PointData + UD ( ID ( IJ , IDim ) ) * UD ( ID ( IJ , IDim ) ) ;
!      End Do ;
!
!    PointData = Sqrt ( PointData ) ;
!    !Write (Unit = UnFile, FMT = "(E31.23E3,1x)", ADVANCE = 'NO', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) PointData ;
!    Write (Unit = UnFile, FMT = "(E31.23E3,1X)", ADVANCE = 'Yes', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) PointData ;
!  End Do ;
!
!!Write (Unit = UnFile, FMT = "()", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 )  ;
!Write (Unit = UnFile, FMT = "('        </DataArray>')", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 )  ;
!
!! Velocity in X direction
!Write (Unit = UnFile, FMT = "(A66)", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) '        <DataArray type="Float32" Name="X vel" format="ascii">' ;
!
!  Do IJ = 1, NJ_Para ;
!    If ( ParaNodeNum ( IJ ) == 0_Lng ) cycle ;
!    PointData = 0.0_Dbl ;
!    If ( ID ( IJ , 1 ) /= 0_Lng ) PointData = UD ( ID ( IJ , 1 ) ) ;
!
!    !Write (Unit = UnFile, FMT = "(E31.23E3,1x)", ADVANCE = 'NO', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) PointData ;
!    Write (Unit = UnFile, FMT = "(E31.23E3,1X)", ADVANCE = 'Yes', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) PointData ;
!  End Do ;
!
!!Write (Unit = UnFile, FMT = "()", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 )  ;
!Write (Unit = UnFile, FMT = "('        </DataArray>')", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 )  ;
!
!! Displacement in Y direction
!Write (Unit = UnFile, FMT = "(A66)", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) '        <DataArray type="Float32" Name="Y Vel" format="ascii">' ;
!
!  Do IJ = 1, NJ_Para ;
!    If ( ParaNodeNum ( IJ ) == 0_Lng ) cycle ;
!    PointData = 0.0_Dbl ;
!    If ( ID ( IJ , 2 ) /= 0_Lng ) PointData = UD ( ID ( IJ , 2 ) ) ;
!
!    !Write (Unit = UnFile, FMT = "(E31.23E3,1x)", ADVANCE = 'NO', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) PointData ;
!    Write (Unit = UnFile, FMT = "(E31.23E3,1X)", ADVANCE = 'Yes', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) PointData ;
!  End Do ;

!Write (Unit = UnFile, FMT = "()", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 )  ;
!Write (Unit = UnFile, FMT = "('        </DataArray>')", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 )  ;


Write (Unit = UnFile, FMT = "('      </PointData>')", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 )  ;

! - Cell data -----------------------------------

! - Footer --------------------------------------
Write (Unit = UnFile, FMT = "('    </Piece>')", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 )  ;
Write (Unit = UnFile, FMT = "('  </UnstructuredGrid>')", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 )  ;
Write (Unit = UnFile, FMT = "('</VTKFile>')", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 )  ;

! - Closing the output file ---------------------------------------------------------------------------------------------------------------------
UnFile =  UN_PView ;
Close ( Unit = UnFile, Status = 'KEEP', ERR =  1002, IOSTAT = IO_File ) ;

!#Write(*    ,*) 'End Subroutine < ParaView_27noded >' ;
!#Write(UnInf,*) 'End Subroutine < ParaView_27noded >' ;
Return ;

! =============================================== OPEN ERRORS =======================================================================================
1001  IF ( IO_File > 0 ) Then ;
        Write(*, Fmt_ERR1_OPEN) UnFile, IO_File  ;  Write(UnInf, Fmt_ERR1_OPEN) UnFile, IO_File  ; 
        Write(*, Fmt_FL) ; Write(UnInf, Fmt_FL) ;
        !#Call BEEP_FAIL ;
        Write(*, Fmt_End) ; Read(*,*) ;  STOP ;
      Else If ( IO_File < 0 ) Then ;
        Write(*, Fmt_ERR1_OPEN) UnFile, IO_File  ; 
        Write(UnInf, Fmt_ERR1_OPEN) UnFile, IO_File  ;  Write(*, Fmt_FL) ; Write(UnInf, Fmt_FL) ;
        !#Call BEEP_FAIL ;
        Write(*, Fmt_End) ; Read(*,*) ;  STOP ;
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
        Write(*, Fmt_FL) ;  Write(UnInf, Fmt_FL) ; Write(*, Fmt_End) ; Read(*,*) ;  Return ;

End Subroutine ParaView_27noded ;


! ***************************************************************************************************************************************************
!                                                                                                                                                  **
! Start:        09 Nov 2012                                                                                                                        **
! Last Update:  09 Nov 2012                                                                                                                        **
! Description: This subroutine renumbers the nodes for downsampling of data for paraview.                                                          **
! Called by: subrotines marching in time                                                                                                           **
!                                                                                                                                                  **
! ***************************************************************************************************************************************************

Subroutine downsampling (                                                                                                       &
NNodeDown,                                                                                                                      & ! Integer (1) Variables
!                                                                                                                               & ! Integer (2) Variables
!                                                                                                                               & ! Integer (4) Variables
NEl, NJ_Para,                                                                                                                   & ! Integer (8) Variables
!                                                                                                                               & ! Real Variables
INod, ParaNodeNum                                                                                                               & ! Integer Arrays
!                                                                                                                               & ! Real Arrays
!                                                                                                                               & ! Characters
!                                                                                                                               & ! Type
) ;

Implicit None ;

! =========================== Global Variables ======================================================================================================

! - Integer Variables -------------------------------------------------------------------------------------------------------------------------------
Integer , Intent(In)     :: NEl ;

Integer , Intent(In)     :: NNodeDown ;
Integer , Intent(Out)    :: NJ_Para ;

! - Real Variables ----------------------------------------------------------------------------------------------------------------------------------
! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------

! - Integer Arrays ----------------------------------------------------------------------------------------------------------------------------------
Integer , Intent(In),     Dimension (:,:)  :: INod ;

Integer , Intent(InOut),  Dimension (:)    :: ParaNodeNum ;

! - Real Arrays -------------------------------------------------------------------------------------------------------------------------------------
!#Real (Kind=Dbl),     Intent(In),    Dimension (:  )  ::  ;
!#Real (Kind=Dbl),     Intent(InOut), Dimension (:  )  ::  ;
!#Real (Kind=Dbl),     Intent(OUT),   Dimension (:  )  ::  ;
! - Character Variables -----------------------------------------------------------------------------------------------------------------------------
!Character (Kind = ?, Len = ? ) :: ;
! - Logical Variables -------------------------------------------------------------------------------------------------------------------------------
!#Logical   ::  ;
! - Type DECLERATIONS -------------------------------------------------------------------------------------------------------------------------------
!#Type() :: ;
! =========================== Local Variables =======================================================================================================
! - Integer Variables -------------------------------------------------------------------------------------------------------------------------------
Integer   :: INode ;                  ! Loop over node numbers of element 

Integer   :: IEl ;                    ! Loop over element number
Integer   :: Node ;                   ! temporary variable for node number

! - Real Variables ----------------------------------------------------------------------------------------------------------------------------------
!#Real (Kind=Dbl)      ::  ;
! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
!#Complex              ::  ;
! - Integer Arrays ----------------------------------------------------------------------------------------------------------------------------------
!#Integer (Kind=Shrt), Dimension (?)  ::  ;
! - Real Arrays -------------------------------------------------------------------------------------------------------------------------------------
!#Real (Kind=Dbl), Dimension (?)      ::  ;
! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
!#Complex  ::  ;
! - Character Variables -----------------------------------------------------------------------------------------------------------------------------
!Character (Kind = ?, Len = ? ) :: ;
! - Logical Variables -------------------------------------------------------------------------------------------------------------------------------
!#Logical   ::  ;
! - Type DECLERATIONS -------------------------------------------------------------------------------------------------------------------------------

! =========================== Subroutine CODE =======================================================================================================

! <<<<<<<< for now this subroutine is designed for serial codes >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

ParaNodeNum (:) = 0 ;
NJ_Para = 0 ;

  Do IEl = 1, NEl ;

    Do INode = 1, NNodeDown;

      Node = INod (INode, IEl) ;

        If ( ParaNodeNum ( Node ) == 0 ) Then ;
          NJ_Para = NJ_Para + 1 ;
          ParaNodeNum ( Node ) = NJ_Para ;
        End If ;

    End Do ;

  End Do ;

Write(*    ,*) 'End Subroutine < Downsampling >' ;
Write(UnInf,*) 'End Subroutine < Downsampling >' ;
Return ;
End Subroutine Downsampling ;


! ***************************************************************************************************************************************************
!                                                                                                                                                  **
! Start:        05 July 2013                                                                                                                       **
! Last Update:  05 July 2013                                                                                                                       **
!                                                                                                                                                  **
! Developed by: Babak Poursartip                                                                                                                   **
! Called from:  RK                                                                                                                                 **
!                                                                                                                                                  **
! Description: This subroutine writes down the solution , i.e. displacemenets, velocity, and acceleration for Paraview in HDF5 format              **
!                                                                                                                                                  **
! ***************************************************************************************************************************************************

Subroutine Paraview_HDF5    (                                                                                                   &
NDim, MaxNNode,                                                                                                                 & ! Integer (1) Variables
!                                                                                                                               & ! Integer (2) Variables
!                                                                                                                               & ! Integer (4) Variables
IStep, NEl, NJ_Para, ConnSizePara,                                                                                              & ! Integer (8) Variables
Time_n1,                                                                                                                        & ! Real Variables
ID_Para,                                                                                                                        & ! Integer Arrays
U_Rank, UD_Rank, UDD_Rank,                                                                                                      & ! Real Arrays
OutDir, ModelName, AnaName, IndexRank, IndexSize                                                                                & ! Characters
!                                                                                                                               & ! Type
) ;

Implicit None ;

! =========================== PETSC LIBRARIES =======================================================================================================
#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscmat.h"
#include "finclude/petscis.h"
#include "finclude/petscvec.h90"

! =========================== Global Variables ======================================================================================================

! - Integer Variables -------------------------------------------------------------------------------------------------------------------------------
Integer , Intent(In)    :: NDim, MaxNNode ;
Integer , Intent(In)    :: IStep, NEl, NJ_Para, ConnSizePara ;

! - Real Variables ----------------------------------------------------------------------------------------------------------------------------------
Real(8),     Intent(In)    :: Time_n1 ;

! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
!#Complex,             Intent(In)    ::  ;
!#Complex,             Intent(InOut) ::  ;
!#Complex,             Intent(OUT)   ::  ;
! - Integer Arrays ----------------------------------------------------------------------------------------------------------------------------------
Integer , Intent(In), Dimension (:,:)  :: ID_Para ;

! - Real Arrays -------------------------------------------------------------------------------------------------------------------------------------
!Real(8), Intent(IN)  , Dimension (:  )  :: U_Rank, UD_Rank, UDD_Rank ;
Real(8), Dimension (NJ_Para)  :: U_Rank, UD_Rank, UDD_Rank ;

! - Character Variables -----------------------------------------------------------------------------------------------------------------------------
Character (Kind = 1, Len = 150) :: OutDir ;      ! Directory of output files (Results)
Character (Kind = 1, Len = 20)  :: IndexRank ;   ! A variable for storing the Rank number in Character format for adding at the end of the input file Name.
Character (Kind = 1, Len = 20)  :: IndexSize ;   ! A variable for storing the Size number in Character format for adding at the end of the input file Name.
Character (Kind = 1, Len = 30 ) :: ModelName ;   ! name of the model input file
Character (Kind = 1, Len = 30 ) :: AnaName ;     ! Name of Input file

! - Logical Variables -------------------------------------------------------------------------------------------------------------------------------
!#Logical   ::  ;
! - Type DECLERATIONS -------------------------------------------------------------------------------------------------------------------------------
!#Type() :: ;

! ============================ PETSC Variables AND OBJECTS ==========================================================================================
! - PETSC INTERNAL Variables ------------------------------------------------------------------------------------------------------------------------
PetscErrorCode :: ErrPTC ;                       ! Error 
PetscMPIInt    :: Size ;                         ! Total number of ranks
PetscMPIInt    :: Rank ;                         ! Rank number.

! =========================== Local Variables =======================================================================================================
! - Integer Variables -------------------------------------------------------------------------------------------------------------------------------
Integer   :: UnFile ;                 ! Holds unit of a file for error message
Integer   :: IO_File ;                ! Used for IOSTAT - Input Output Status - in the OPEN cammand.
Integer   :: IO_Write ;               ! Used for IOSTAT - Input Output Status - in the Write cammand.

Integer   :: Counter ;                 ! Counter
Integer   :: I, J ;                    ! Loop Indices

Integer(HID_T)       :: id_Data ;                 ! File identifier for the Data file

Integer(HID_T)       :: dspace_id_Dis ;           ! Dataspace identifier for Displacement
Integer(HID_T)       :: dspace_id_Vel ;           ! Dataspace identifier for Velocity
Integer(HID_T)       :: dspace_id_Acc ;           ! Dataspace identifier for Acceleration

Integer(HID_T)       :: dset_id_Dis ;             ! Dataset identifier for Displacement 
Integer(HID_T)       :: dset_id_Vel ;             ! Dataset identifier for Velocity 
Integer(HID_T)       :: dset_id_Acc ;             ! Dataset identifier for Acceleration 

Integer(HSIZE_T), DIMENSION(2) :: dims           ! Dataset dimensions
Integer(HSIZE_T), DIMENSION(2) :: data_dims

! - Real Variables ----------------------------------------------------------------------------------------------------------------------------------
!#Real (Kind=Dbl)      ::  ;
! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
!#Complex              ::  ;
! - Integer Arrays ----------------------------------------------------------------------------------------------------------------------------------
!#Integer (Kind=Shrt), Dimension (?)  ::  ;
! - Real Arrays -------------------------------------------------------------------------------------------------------------------------------------
!Real (Kind=Dbl), Dimension ( 3, NJ_Para )      :: dset_data  ; ! NDim
Real(4), Dimension ( 3, NJ_Para )      :: dset_data  ; ! NDim

! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
!#Complex  ::  ;
! - Character Variables -----------------------------------------------------------------------------------------------------------------------------
Character (Kind = 1, Len = 25)  :: IndexStep ;   ! A variable for storing the Step number in Character format for adding at the end of the data file

! - Logical Variables -------------------------------------------------------------------------------------------------------------------------------
!#Logical   ::  ;
! - Type DECLERATIONS -------------------------------------------------------------------------------------------------------------------------------

! =========================== Subroutine CODE =======================================================================================================

! - Write data file on the rank ---------------------------------------------------------------------------------------------------------------------

Write(IndexStep, *)IStep;  ! Converts Step number to Character foramt for the data file

! Create files based on HDF5 format
Call h5open_f(ErrPTC)

Call h5fcreate_f( TRIM(OutDir)//'/'//'Data'//'_'//Trim(AdjustL(IndexSize))//'_'//Trim(AdjustL(IndexRank))//'_'//Trim(AdjustL(IndexStep))//'.h5', H5F_ACC_TRUNC_F, id_Data, ErrPTC) ;      ! Data file for the solution

! - Create the dataspaces
! Coordinate file for the main code(.XYZ)
dims(1) = 3 ; !NDim  ;
dims(2) = NJ_Para ;
Call h5screate_simple_f(2, dims, dspace_id_dis, ErrPTC) ;
Call h5screate_simple_f(2, dims, dspace_id_vel, ErrPTC) ;
Call h5screate_simple_f(2, dims, dspace_id_acc, ErrPTC) ;

!Call h5dcreate_f(id_Data, "Dis",             H5T_NATIVE_DOUBLE,  dspace_id_Dis, dset_id_Dis, ErrPTC) ;
!Call h5dcreate_f(id_Data, "Vel",             H5T_NATIVE_DOUBLE,  dspace_id_Vel, dset_id_Vel, ErrPTC) ;
!Call h5dcreate_f(id_Data, "Acc",             H5T_NATIVE_DOUBLE,  dspace_id_Acc, dset_id_Acc, ErrPTC) ;
Call h5dcreate_f(id_Data, "Dis",             H5T_NATIVE_REAL,  dspace_id_Dis, dset_id_Dis, ErrPTC) ;
Call h5dcreate_f(id_Data, "Vel",             H5T_NATIVE_REAL,  dspace_id_Vel, dset_id_Vel, ErrPTC) ;
Call h5dcreate_f(id_Data, "Acc",             H5T_NATIVE_REAL,  dspace_id_Acc, dset_id_Acc, ErrPTC) ;

! Write the dataset.
data_dims(1) = 3 ; !NDim ;
data_dims(2) = NJ_Para ;

!dset_data ( 3, : ) = 0.0_Dbl ;
dset_data ( 3, : ) = 0.0 ;

  ! Displacement
Counter = 0 ;
  Do I = 1, data_dims(2) ; ! NJ_Para
    Do J = 1, NDim ;
      If ( ID_Para ( I, J ) == -1 ) Then ;
        !dset_data ( J, I ) = 0.0_Dbl ;
        dset_data ( J, I ) = 0.0 ;
      Else 
        Counter = Counter + 1 ;
        !dset_data ( J, I ) = U_Rank ( Counter ) ;
        dset_data ( J, I ) = Sngl(U_Rank ( Counter )) ;
      End If ;
    End Do ;
  End Do ;


!Call h5dwrite_f(dset_id_Dis, H5T_NATIVE_DOUBLE, dset_data, data_dims, ErrPTC)
Call h5dwrite_f(dset_id_Dis, H5T_NATIVE_REAL, dset_data, data_dims, ErrPTC)

!dset_data ( 3, : ) = 0.0_Dbl ;
dset_data ( 3, : ) = 0.0 ;

! Velocity
Counter = 0 ;
  Do I = 1, data_dims(2) ; ! NJ_Para
    Do J = 1, NDim ;
      If ( ID_Para ( I, J ) == -1 ) Then ;
        !dset_data ( J, I ) = 0.0_Dbl ;
        dset_data ( J, I ) = 0.0 ;
      Else 
        Counter = Counter + 1 ;
        !dset_data ( J, I ) = UD_Rank ( Counter ) ;
        dset_data ( J, I ) = Sngl(UD_Rank ( Counter ) ) ;
      End If ;
    End Do ;
  End Do ;

!Call h5dwrite_f(dset_id_Vel, H5T_NATIVE_DOUBLE, dset_data, data_dims, ErrPTC)
Call h5dwrite_f(dset_id_Vel, H5T_NATIVE_REAL, dset_data, data_dims, ErrPTC)

!dset_data ( 3, : ) = 0.0_Dbl ;
dset_data ( 3, : ) = 0.0 ;

! Acceleration
Counter = 0 ;
  Do I = 1, data_dims(2) ; ! NJ_Para
    Do J = 1, NDim ;
      If ( ID_Para ( I, J ) == -1 ) Then ;
        !dset_data ( J, I ) = 0.0_Dbl ;
        dset_data ( J, I ) = 0.0 ;
      Else 
        Counter = Counter + 1 ;
        !dset_data ( J, I ) = UDD_Rank ( Counter ) ;
        dset_data ( J, I ) = Sngl (UDD_Rank ( Counter ) ) ;
      End If ;
    End Do ;
  End Do ;

!Call h5dwrite_f(dset_id_Acc, H5T_NATIVE_DOUBLE, dset_data, data_dims, ErrPTC)
Call h5dwrite_f(dset_id_Acc, H5T_NATIVE_REAL, dset_data, data_dims, ErrPTC)

! Close the dataset.
Call h5dclose_f( dset_id_dis, ErrPTC) ;
Call h5dclose_f( dset_id_vel, ErrPTC) ;
Call h5dclose_f( dset_id_acc, ErrPTC) ;

! Close the HDF5 file.
CALL h5fclose_f( id_Data,     ErrPTC) ;

! - Write rank wrap file  ---------------------------------------------------------------------------------------------------------------------------

! Open Rank wraper
UnFile = UN_OutPara ;
Open ( Unit = UnFile, FILE = TRIM(AnaNAME)//'_'//Trim(AdjustL(IndexSize))//'_'//Trim(AdjustL(IndexRank))//'_'//Trim(AdjustL(IndexStep))//'.xmf', ERR =  1001, IOSTAT = IO_File, ACCESS = 'SEQUENTIAL', ACTION ='Write', ASYNCHRONOUS = 'NO', BLANK = 'NULL', BLOCKSIZE = 0, DEFAULTFILE = TRIM(OutDir), DisPOSE = 'KEEP', FORM = 'FORMATTED', POSITION = 'ASIS', STATUS ='REPLACE') ;

Write (Unit = UnFile, FMT = "(A54)", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) "<Grid GridType='Collection' CollectionType='Temporal'>" ;
Write (Unit = UnFile, FMT = "(A27)", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) '  <Grid GridType="Uniform">' ;

Write (Unit = UnFile, FMT = "(A17 , F17.10, A4)", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) '    <Time Value="',Time_n1,'" />' ;
Write (Unit = UnFile, FMT = "(A53,I19,A9)", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) '    <Topology TopologyType="Mixed" NumberOfElements="',NEl,'" >' ;
!Write (Unit = UnFile, FMT = "(A51,I19,A9)", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) '  <Topology TopologyType="Quad_8" Dimensions="',NEl,'" >' ;   ! You can use this line instead of mixed if there is only one type of cell(element) in the model.
!Write (Unit = UnFile, FMT = "(A51,I19,A9)", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) '  <Topology TopologyType="HEX 20" Dimensions="',NEl,'" >' ;   ! You can use this line instead of mixed if there is only one type of cell(element) in the model.

Write (Unit = UnFile, FMT = "(A26,I19,A30,A184,A11)", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) '    <DataItem Dimensions="',ConnSizePara,'" DataType="Int" Format="HDF">',TRIM(ModelName)//'_'//'Geometry_PV'//'_'//Trim(AdjustL(IndexSize))//'_'//Trim(AdjustL(IndexRank))//'.h5:Connectivity','</DataItem>' ;

Write (Unit = UnFile, FMT = "(A15)", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) '    </Topology>' ;
Write (Unit = UnFile, FMT = "(A33)", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) '    <Geometry GeometryType="XYZ">' ;
Write (Unit = UnFile, FMT = "(A28,I19,A51,A175,A11)", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) '      <DataItem Dimensions="',NJ_Para,' 3" NumberType="Float" Precision="8" Format="HDF">',TRIM(ModelName)//'_'//'Geometry_PV'//'_'//Trim(AdjustL(IndexSize))//'_'//Trim(AdjustL(IndexRank))//'.h5:XYZ','</DataItem>' ;
Write (Unit = UnFile, FMT = "(A15)", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) '    </Geometry>' ;

Write (Unit = UnFile, FMT = "(A72)", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) '    <Attribute Name="Displacement" AttributeType="Vector" Center="Node">' ;
Write (Unit = UnFile, FMT = "(A28,I19,A50,A175,A11)", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) '      <DataItem Dimensions="',NJ_Para,' 3" NumberType="Float" Precision="8" Format="HDF">','Data'//'_'//Trim(AdjustL(IndexSize))//'_'//Trim(AdjustL(IndexRank))//'_'//Trim(AdjustL(IndexStep))//'.h5:Dis','</DataItem>' ;
Write (Unit = UnFile, FMT = "(A16)", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) '    </Attribute>' ;

Write (Unit = UnFile, FMT = "(A68)", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) '    <Attribute Name="Velocity" AttributeType="Vector" Center="Node">' ;
Write (Unit = UnFile, FMT = "(A28,I19,A50,A175,A11)", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) '      <DataItem Dimensions="',NJ_Para,' 3" NumberType="Float" Precision="8" Format="HDF">','Data'//'_'//Trim(AdjustL(IndexSize))//'_'//Trim(AdjustL(IndexRank))//'_'//Trim(AdjustL(IndexStep))//'.h5:Vel','</DataItem>' ;
Write (Unit = UnFile, FMT = "(A16)", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) '    </Attribute>' ;

Write (Unit = UnFile, FMT = "(A72)", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) '    <Attribute Name="Acceleration" AttributeType="Vector" Center="Node">' ;
Write (Unit = UnFile, FMT = "(A28,I19,A50,A175,A11)", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) '      <DataItem Dimensions="',NJ_Para,' 3" NumberType="Float" Precision="8" Format="HDF">','Data'//'_'//Trim(AdjustL(IndexSize))//'_'//Trim(AdjustL(IndexRank))//'_'//Trim(AdjustL(IndexStep))//'.h5:Acc','</DataItem>' ;
Write (Unit = UnFile, FMT = "(A16)", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) '    </Attribute>' ;

Write (Unit = UnFile, FMT = "(A9)", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) '  </Grid>' ;
Write (Unit = UnFile, FMT = "(A7)", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) '</Grid>' ;

!  <Attribute Name="scalar" AttributeType="Scalar" Center="Node">
!    <DataItem Dimensions="320 1" NumberType="Float" Precision="4" Format="HDF">example.h5:part-0-scalars</DataItem>
!  </Attribute>

! Coordinate file (.XYZ)
UnFile = UN_OutPara ;
Close ( Unit = UnFile, STATUS = 'KEEP', ERR = 1002, IOSTAT = IO_File ) ;

!#Write(*    ,*) 'End Subroutine < Paraview_HDF5 >' ;
!#Write(UnInf,*) 'End Subroutine < Paraview_HDF5 >' ;
Return ;

! =============================================== OPEN ERRORS =======================================================================================
1001  IF ( IO_File > 0 ) Then ;
        Write(*, Fmt_ERR1_OPEN) UnFile, IO_File  ;  Write(UnInf, Fmt_ERR1_OPEN) UnFile, IO_File  ; 
        Write(*, Fmt_FL) ; Write(UnInf, Fmt_FL) ;
        !#Call BEEP_FAIL ;
        Write(*, Fmt_End) ; Read(*,*) ;  STOP ;
      Else If ( IO_File < 0 ) Then ;
        Write(*, Fmt_ERR1_OPEN) UnFile, IO_File  ; 
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

! =============================================== ERROR IN Write STATEMENT ==========================================================================
1006    Write(*       , Fmt_Write1 ) UnFile, IO_Write ; Write( UnFile, Fmt_Write1 ) UnFile, IO_Write ;
        !#Call BEEP_FAIL ;
        Write(*, Fmt_FL) ;  Write(UnInf, Fmt_FL) ; Write(*, Fmt_End) ; Read(*,*) ;  Return ;

End Subroutine Paraview_HDF5 ;



! ***************************************************************************************************************************************************
!                                                                                                                                                  **
! Start:        03 April 2014                                                                                                                      **
! Last Update:  03 April 2014                                                                                                                      **
!                                                                                                                                                  **
! Developed by: Babak Poursartip, AF                                                                                                               **
! Called from:                                                                                                                                     **
!                                                                                                                                                  **
! Description: This subroutine writes down the material updates , i.e. Lambda & Mu for Paraview in HDF5 format                                     **
!                                                                                                                                                  **
! ***************************************************************************************************************************************************

Subroutine Paraview_Mat_HDF5 (                                                                                                  &
!                                                                                                                               & ! Integer (1) Variables
!                                                                                                                               & ! Integer (2) Variables
!                                                                                                                               & ! Integer (4) Variables
Iter, NJ_Para, ConnSizePara,                                                                                                    & ! Integer (8) Variables
Time_n1,                                                                                                                        & ! Real Variables
!                                                                                                                               & ! Integer Arrays
Lambda_Rank_Vis, Mu_Rank_Vis,                                                                                                   & ! Real Arrays
MatUpdateDir, ModelName, AnaName, IndexRank, IndexSize                                                                          & ! Characters
!                                                                                                                               & ! Type
) ;

Implicit None ;

! =========================== PETSC LIBRARIES =======================================================================================================
#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscmat.h"
#include "finclude/petscis.h"
#include "finclude/petscvec.h90"

! =========================== Global Variables ======================================================================================================

! - Integer Variables -------------------------------------------------------------------------------------------------------------------------------
Integer , Intent(In)    :: Iter, NJ_Para, ConnSizePara; ! ConnSizePara = Param%IntM( 5, 4)

! - Real Variables ----------------------------------------------------------------------------------------------------------------------------------
Real(8),  Intent(In)    :: Time_n1 ;

! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
! - Integer Arrays ----------------------------------------------------------------------------------------------------------------------------------
! - Real Arrays -------------------------------------------------------------------------------------------------------------------------------------
Real(8), Intent(In)  , Dimension (NJ_Para)  :: Lambda_Rank_Vis, Mu_Rank_Vis ;

! - Character Variables -----------------------------------------------------------------------------------------------------------------------------
Character (Kind = 1, Len = 150) :: MatUpdateDir ; ! Directory of output files (Results)
Character (Kind = 1, Len = 20)  :: IndexRank ;    ! A variable for storing the Rank number in Character format for adding at the end of the input file Name.
Character (Kind = 1, Len = 20)  :: IndexSize ;    ! A variable for storing the Size number in Character format for adding at the end of the input file Name.
Character (Kind = 1, Len = 30 ) :: ModelName ;    ! name of the model input file
Character (Kind = 1, Len = 30 ) :: AnaName ;      ! Name of Input file

! - Logical Variables -------------------------------------------------------------------------------------------------------------------------------
!#Logical   ::  ;
! - Type DECLERATIONS -------------------------------------------------------------------------------------------------------------------------------
!#Type() :: ;

! ============================ PETSC Variables AND OBJECTS ==========================================================================================
! - PETSC INTERNAL Variables ------------------------------------------------------------------------------------------------------------------------
PetscErrorCode :: ErrPTC ;                       ! Error 
PetscMPIInt    :: Size ;                         ! Total number of ranks
PetscMPIInt    :: Rank ;                         ! Rank number.

! =========================== Local Variables =======================================================================================================
! - Integer Variables -------------------------------------------------------------------------------------------------------------------------------
Integer   :: UnFile ;                 ! Holds unit of a file for error message
Integer   :: IO_File ;                ! Used for IOSTAT - Input Output Status - in the OPEN cammand.
Integer   :: IO_Write ;               ! Used for IOSTAT - Input Output Status - in the Write cammand.

Integer   :: Counter ;                ! Counter
Integer   :: I, J ;                   ! Loop Indices

Integer(HID_T)       :: id_Data ;                ! File identifier for the Data file

Integer(HID_T)       :: dspace_id_Lambda ;       ! Dataspace identifier for Lambda
Integer(HID_T)       :: dspace_id_Mu ;           ! Dataspace identifier for Mu

Integer(HID_T)       :: dset_id_Lambda ;         ! Dataset identifier for Lambda
Integer(HID_T)       :: dset_id_Mu ;             ! Dataset identifier for Mu

Integer(HSIZE_T), DIMENSION(2) :: dims           ! Dataset dimensions
Integer(HSIZE_T), DIMENSION(2) :: data_dims

! - Real Variables ----------------------------------------------------------------------------------------------------------------------------------
!#Real (Kind=Dbl)      ::  ;
! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
!#Complex              ::  ;
! - Integer Arrays ----------------------------------------------------------------------------------------------------------------------------------
!#Integer (Kind=Shrt), Dimension (?)  ::  ;
! - Real Arrays -------------------------------------------------------------------------------------------------------------------------------------
Real(4), Dimension ( 1, NJ_Para )      :: dset_data  ;

! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
!#Complex  ::  ;
! - Character Variables -----------------------------------------------------------------------------------------------------------------------------
Character (Kind = 1, Len = 25)  :: IndexStep ;   ! A variable for storing the Step number in Character format for adding at the end of the data file

! - Logical Variables -------------------------------------------------------------------------------------------------------------------------------
!#Logical   ::  ;
! - Type DECLERATIONS -------------------------------------------------------------------------------------------------------------------------------

! =========================== Subroutine CODE =======================================================================================================

! - Write data file on the rank ---------------------------------------------------------------------------------------------------------------------

Write(IndexStep, *) Iter;  ! Converts Step number to Character foramt for the data file

! Create files based on HDF5 format
Call h5open_f(ErrPTC)

Call h5fcreate_f(TRIM(MatUpdateDir)//'/'//'Mat'//'_'//Trim(AdjustL(IndexSize))//'_'//Trim(AdjustL(IndexRank))//'_'//Trim(AdjustL(IndexStep))//'.h5', H5F_ACC_TRUNC_F, id_Data, ErrPTC) ;      ! Data file for the solution

! - Create the dataspaces
dims(1) = 1 ; !NDim  ;
dims(2) = NJ_Para ;
Call h5screate_simple_f(2, dims, dspace_id_Lambda, ErrPTC) ;
Call h5screate_simple_f(2, dims, dspace_id_Mu,     ErrPTC) ;

Call h5dcreate_f(id_Data, "Lambda",         H5T_NATIVE_REAL,  dspace_id_Lambda, dset_id_Lambda, ErrPTC) ;
Call h5dcreate_f(id_Data, "Mu",             H5T_NATIVE_REAL,  dspace_id_Mu,     dset_id_Mu,     ErrPTC) ;

! Write the dataset.
data_dims(1) = 1 ;
data_dims(2) = NJ_Para ;


! Lambda
!  Do I = 1, data_dims(2) ; ! NJ_Para
!        dset_data ( 1, I ) = Sngl(Lambda_Rank_Vis ( I )) ;
!  End Do ;
  If ( Lambda_variable == 'Lambda_2Mu' ) Then
     Do I = 1, data_dims(2) ; ! NJ_Para
           dset_data ( 1, I ) = Sngl( Lambda_Rank_Vis ( I ) - 2.0d0 * Mu_Rank_Vis ( I ) ) ;
     End Do ;
  Else
     Do I = 1, data_dims(2) ; ! NJ_Para
           dset_data ( 1, I ) = Sngl( Lambda_Rank_Vis ( I ) ) ;
     End Do ;
  End If

Call h5dwrite_f(dset_id_Lambda, H5T_NATIVE_REAL, dset_data, data_dims, ErrPTC)


! Mu
  Do I = 1, data_dims(2) ; ! NJ_Para
     J = 1 ;
        dset_data ( J, I ) = Sngl( Mu_Rank_Vis ( I ) ) ;
  End Do ;

Call h5dwrite_f(dset_id_Mu, H5T_NATIVE_REAL, dset_data, data_dims, ErrPTC)


! Close the dataset.
Call h5dclose_f( dset_id_Lambda, ErrPTC) ;
Call h5dclose_f( dset_id_Mu,     ErrPTC) ;

! Close the HDF5 file.
CALL h5fclose_f( id_Data,     ErrPTC) ;

! - Write rank wrap file  ---------------------------------------------------------------------------------------------------------------------------

! Open Rank wraper
UnFile = UN_OutMatPara ;
Open ( Unit = UnFile, FILE = TRIM(AnaNAME)//'_'//Trim(AdjustL(IndexSize))//'_'//Trim(AdjustL(IndexRank))//'_'//Trim(AdjustL(IndexStep))//'.xmf', ERR =  1001, IOSTAT = IO_File, ACCESS = 'SEQUENTIAL', ACTION ='Write', ASYNCHRONOUS = 'NO', BLANK = 'NULL', BLOCKSIZE = 0, DEFAULTFILE = TRIM(MatUpdateDir), DisPOSE = 'KEEP', FORM = 'FORMATTED', POSITION = 'ASIS', STATUS ='REPLACE') ;

Write (Unit = UnFile, FMT = "(A54)", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) "<Grid GridType='Collection' CollectionType='Temporal'>" ;
Write (Unit = UnFile, FMT = "(A27)", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) '  <Grid GridType="Uniform">' ;

Write (Unit = UnFile, FMT = "(A17 , F17.10, A4)", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) '    <Time Value="',Time_n1,'" />' ;
Write (Unit = UnFile, FMT = "(A53,I19,A9)", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) '    <Topology TopologyType="Mixed" NumberOfElements="',NEl,'" >' ;
!Write (Unit = UnFile, FMT = "(A51,I19,A9)", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) '  <Topology TopologyType="Quad_8" Dimensions="',NEl,'" >' ;   ! You can use this line instead of mixed if there is only one type of cell(element) in the model.
!Write (Unit = UnFile, FMT = "(A51,I19,A9)", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) '  <Topology TopologyType="HEX 20" Dimensions="',NEl,'" >' ;   ! You can use this line instead of mixed if there is only one type of cell(element) in the model.

Write (Unit = UnFile, FMT = "(A26,I19,A30,A184,A11)", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) '    <DataItem Dimensions="',ConnSizePara,'" DataType="Int" Format="HDF">',TRIM(ModelName)//'_'//'Geometry_PV'//'_'//Trim(AdjustL(IndexSize))//'_'//Trim(AdjustL(IndexRank))//'.h5:Connectivity','</DataItem>' ;

Write (Unit = UnFile, FMT = "(A15)", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) '    </Topology>' ;
Write (Unit = UnFile, FMT = "(A33)", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) '    <Geometry GeometryType="XYZ">' ;
Write (Unit = UnFile, FMT = "(A28,I19,A51,A175,A11)", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) '      <DataItem Dimensions="',NJ_Para,' 3" NumberType="Float" Precision="8" Format="HDF">',TRIM(ModelName)//'_'//'Geometry_PV'//'_'//Trim(AdjustL(IndexSize))//'_'//Trim(AdjustL(IndexRank))//'.h5:XYZ','</DataItem>' ;
Write (Unit = UnFile, FMT = "(A15)", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) '    </Geometry>' ;

Write (Unit = UnFile, FMT = "(A72)", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) '    <Attribute Name="Lambda" AttributeType="Scalar" Center="Node">' ;
Write (Unit = UnFile, FMT = "(A28,I19,A50,A175,A11)", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) '      <DataItem Dimensions="',NJ_Para,' 1" NumberType="Float" Precision="8" Format="HDF">','Mat'//'_'//Trim(AdjustL(IndexSize))//'_'//Trim(AdjustL(IndexRank))//'_'//Trim(AdjustL(IndexStep))//'.h5:Lambda','</DataItem>' ;
Write (Unit = UnFile, FMT = "(A16)", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) '    </Attribute>' ;

Write (Unit = UnFile, FMT = "(A68)", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) '    <Attribute Name="Mu" AttributeType="Scalar" Center="Node">' ;
Write (Unit = UnFile, FMT = "(A28,I19,A50,A175,A11)", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) '      <DataItem Dimensions="',NJ_Para,' 1" NumberType="Float" Precision="8" Format="HDF">','Mat'//'_'//Trim(AdjustL(IndexSize))//'_'//Trim(AdjustL(IndexRank))//'_'//Trim(AdjustL(IndexStep))//'.h5:Mu','</DataItem>' ;
Write (Unit = UnFile, FMT = "(A16)", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) '    </Attribute>' ;


Write (Unit = UnFile, FMT = "(A9)", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) '  </Grid>' ;
Write (Unit = UnFile, FMT = "(A7)", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) '</Grid>' ;

!  <Attribute Name="scalar" AttributeType="Scalar" Center="Node">
!    <DataItem Dimensions="320 1" NumberType="Float" Precision="4" Format="HDF">example.h5:part-0-scalars</DataItem>
!  </Attribute>

! Coordinate file (.XYZ)
UnFile = UN_OutMatPara ;
Close ( Unit = UnFile, STATUS = 'KEEP', ERR = 1002, IOSTAT = IO_File ) ;

!#Write(*    ,*) 'End Subroutine < Paraview_HDF5 >' ;
!#Write(UnInf,*) 'End Subroutine < Paraview_HDF5 >' ;
Return ;

! =============================================== OPEN ERRORS =======================================================================================
1001  IF ( IO_File > 0 ) Then ;
        Write(*, Fmt_ERR1_OPEN) UnFile, IO_File  ;  Write(UnInf, Fmt_ERR1_OPEN) UnFile, IO_File  ; 
        Write(*, Fmt_FL) ; Write(UnInf, Fmt_FL) ;
        !#Call BEEP_FAIL ;
        Write(*, Fmt_End) ; Read(*,*) ;  STOP ;
      Else If ( IO_File < 0 ) Then ;
        Write(*, Fmt_ERR1_OPEN) UnFile, IO_File  ; 
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

! =============================================== ERROR IN Write STATEMENT ==========================================================================
1006    Write(*       , Fmt_Write1 ) UnFile, IO_Write ; Write( UnFile, Fmt_Write1 ) UnFile, IO_Write ;
        !#Call BEEP_FAIL ;
        Write(*, Fmt_FL) ;  Write(UnInf, Fmt_FL) ; Write(*, Fmt_End) ; Read(*,*) ;  Return ;


End Subroutine Paraview_Mat_HDF5 ;



End Module Visualizer ;
