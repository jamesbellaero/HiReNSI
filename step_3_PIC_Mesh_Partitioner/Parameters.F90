
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                                                                                                                  ++
! Start:        01 May 2011                                                                                                                        ++
! Last Update:  21 Jan 2013                                                                                                                        ++
!                                                                                                                                                  ++
! Description:  ALL Parameters HAS BEEN DEFINED IN THIS module                                                                                     ++
!                                                                                                                                                  ++
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Module Parameters ;

Implicit None ;


! =========================== Types OF Integer AND Real variables ===================================================================================
Integer(2), PARAMETER, PUBLIC :: SGL  = SELECTED_Real_Kind ( P = 6 , R = 37  ) ;  ! EQUIVALENT TO Real (4) ! CHECK IT in THE IEEE PAPER
Integer(2), PARAMETER, PUBLIC :: DBL  = SELECTED_Real_Kind ( P = 13, R = 200 ) ;  ! EQUIVALENT TO Real (8) 

Integer(2), PARAMETER, PUBLIC :: Tiny = SELECTED_INT_Kind  ( 1  ) ;               ! EQUIVALENT TO Integer (1)
Integer(2), PARAMETER, PUBLIC :: Smll = SELECTED_INT_Kind  ( 3  ) ;               ! EQUIVALENT TO Integer (2) 
Integer(2), PARAMETER, PUBLIC :: Shrt = SELECTED_INT_Kind  ( 8  ) ;               ! EQUIVALENT TO Integer (4) 
Integer(2), PARAMETER, PUBLIC :: Lng  = SELECTED_INT_Kind  ( 10 ) ;               ! EQUIVALENT TO Integer (8) 


! =========================== Type DECLARATOINS =====================================================================================================

! THIS DERIVED Type ARRAY HOLDS THE ABSCISSAE AND WEIGHTS OF THE GAUSS INTEGRATION 
Type, PUBLIC  :: GAUSS ;
  Real (Kind=DBL), Allocatable, Dimension(:) :: XINT ;  ! ABSCISSAE
  Real (Kind=DBL), Allocatable, Dimension(:) :: WINT ;  ! WEIGHTS
End Type GAUSS ;

! THIS DERIVED Type ARRAY HOLDS THE MASS AND STIFFNESS ENTRIES OF THE Interface NODES, SO THAT WE CAN SUBTARCT THEM FROM THE TOTAL STIFFNESS AND MASS MATRICES TO OBTAIN THE TOTAL EMERGY OF THE REGULAR DOMAIN
Type, PUBLIC  :: PML_Interface ;
  Integer (Kind=Tiny), Allocatable, Dimension(:)   :: NODE ;  ! HOLDS THE LOCATION OF A NODE - See subroutine NodeLocation for details
  Integer (Kind=Lng ), Allocatable, Dimension(:,:) :: Nghb ;  ! Holds all node neighbors of a node on the interface
  Real    (Kind=DBL ), Allocatable, Dimension(:,:) :: MASS ;  ! MASS ELEMENTS OF THE Interface NODES
  Real    (Kind=DBL ), Allocatable, Dimension(:,:) :: STIFF ; ! STIFFNESS ELEMENTS OF THE Interface NODES
End Type PML_Interface ;

! THIS DERIVED Type HOLDS NUMBER OF NON-ZERO ENTRIES OF MASS, STIFFNESS AND DAMPING MATRICES FOR EFFICEINT ALLOCATION ( MALLOC )
Type, PUBLIC  :: NON_ZERO_ENT ; 
  Integer (Kind=Lng ), Allocatable, Dimension(:) :: O_NNZ_STIFF, O_NNZ_MASS, O_NNZ_DAMP, O_NNZ_MB, O_NNZ_P ; ! NUMBER OF NON-ZERO ENTRIES OF DIAGONAL BLOCK
  Integer (Kind=Lng ), Allocatable, Dimension(:) :: D_NNZ_STIFF, D_NNZ_MASS, D_NNZ_DAMP, D_NNZ_MB, D_NNZ_P ; ! NUMBER OF NON-ZERO ENTRIES OF NON-DIAGONAL BLOCK
End Type NON_ZERO_ENT ;

Type, Public  :: NodeID ;
  Integer (Kind=Shrt), Allocatable, Dimension(:,:)  :: Locs  ;  ! Holds rank numbers of which this node belongs to.
  Integer (Kind=Shrt), Allocatable, Dimension(:)    :: Rep   ;  ! Holds the number of Repeats on the ranks (How many times this node appears on ranks), Useful for neighboring
End Type NodeID ;

! Holds basic parameters of each load case
Type, Public  :: BasicParam ;
  Integer (Kind=Lng ), Allocatable, Dimension(:,:)  :: IntM ;   ! Holds Integer Parameters for defining the Model
  Real    (Kind=Dbl ), Allocatable, Dimension(:,:)  :: RealM ;  ! Holds Real Parameters for defining the Model
  Integer (Kind=Lng ), Allocatable, Dimension(:,:)  :: IntL ;   ! Holds Integer Parameters for Loading
  Real    (Kind=Dbl ), Allocatable, Dimension(:,:)  :: RealL ;  ! Holds Real Parameters for Loading
End Type BasicParam ;

! =========================== MATHEMATICAL CONSTATNS ================================================================================================
!#Integer (Kind=Shrt), PARAMETER, PUBLIC  :: ;
Real (Kind=DBL), PARAMETER, PUBLIC  ::  PI = 3.141592653589793238_DBL ;


! =========================== FACE AND AREA ORIENTAION  =============================================================================================
Integer (Kind=Smll), PARAMETER, Dimension(4), PUBLIC  :: LFACE_2D = (/  2,  2,  1,  1                                       /) ;
Integer (Kind=Smll), PARAMETER, Dimension(2), PUBLIC  :: IPERM_2D = (/  2,  1                                               /) ;
Integer (Kind=Smll), PARAMETER, Dimension(6), PUBLIC  :: LFACE_3D = (/  2,  2,  3,  3,  1,  1                               /) ;
Integer (Kind=Smll), PARAMETER, Dimension(3), PUBLIC  :: IPERM_3D = (/  2,  3,  1                                           /) ;

Real    (Kind=DBL ), PARAMETER, Dimension(4), PUBLIC  :: FVAL_2D  = (/ +1._DBL, -1._DBL, +1._DBL, -1._DBL                   /) ;
Real    (Kind=DBL ), PARAMETER, Dimension(6), PUBLIC  :: FVAL_3D  = (/ +1._DBL, -1._DBL, +1._DBL, -1._DBL, +1._DBL, -1._DBL /) ;
!#Integer (Kind=Shrt), PARAMETER, PUBLIC  ::  ;

! =========================== FORMATS ===============================================================================================================
Character(87 ), PARAMETER, PUBLIC   :: Fmt_DATE        = "(' DATE :  ',I2.2,' - ',I2.2,' - ',I4,/,' TIME : ',I2.2,':',I2.2,':',I2.2,':',I2.2,/ )" 
Character(27 ), PARAMETER, PUBLIC   :: Fmt_End         = "('PRESS ENTER TO End ...')" ;
Character(80 ), PARAMETER, PUBLIC   :: Fmt_ERR1_OPEN   = "( 'ERROR IN OPEN STATEMENT. Unit NUMBER = ', I3, '   ERROR NUMBER IS = ', I4  )" ;
Character(91 ), PARAMETER, PUBLIC   :: Fmt_ERR2_OPEN   = "('End-OF-FILE ERROR IN OPEN STATEMENT. Unit NUMBER = ', I3, '   ERROR NUMBER IS = ', I4  )" ;
Character(81 ), PARAMETER, PUBLIC   :: Fmt_ERR1_Close  = "( 'ERROR IN Close STATEMENT. Unit NUMBER = ', I3, '   ERROR NUMBER IS = ', I4  )" ;
Character(163), PARAMETER, PUBLIC   :: Fmt_NM          = "(' FILE Name : ', A20,//,' Directories :',/, 'INPUT FILE DIRECTORY     : ', A100,/, 'OUTPUT FILES DIRECTORY   : ', A100,/, 'INTERNAL FILES DIRECTORY : ', A100,/ )" ;
Character(41 ), PARAMETER, PUBLIC   :: Fmt_SUC         = "('CONGRATULATIONS! DONE SUCCESSFULLY. ')" ;
Character(39 ), PARAMETER, PUBLIC   :: Fmt_FL          = "('OOPS!!!  FAIL TO OPERATE PROPERLY.')" ;
Character(78 ), PARAMETER, PUBLIC   :: Fmt_ALLCT       = "('ERROR IN ALLOCATING Arrays. ERROR NUMBER IS :', I4, '   LOCATION: ??????.')" ;
Character(80 ), PARAMETER, PUBLIC   :: Fmt_DEALLCT     = "('ERROR IN DEALLOCATING Arrays. ERROR NUMBER IS :', I4, '   LOCATION: ??????.')" ;
Character(23 ), PARAMETER, PUBLIC   :: Fmt_RUNTIME     = "(A,F50.2,'   SECONDS')" ;
Character(70 ), PARAMETER, PUBLIC   :: Fmt_READ1       = "('ERROR IN READ STATEMENT. Unit IS : ',I5,' ERROR NUMBER IS : ', I5 )" ;
Character(76 ), PARAMETER, PUBLIC   :: Fmt_READ2       = "('End OF FILE IN READ STATEMENT. Unit IS : ',I5,' ERROR NUMBER IS : ', I5 )" ;
Character(78 ), PARAMETER, PUBLIC   :: Fmt_READ3       = "('End OF RECORD IN READ STATEMENT. Unit IS : ',I5,' ERROR NUMBER IS : ', I5 )" ;
Character(71 ), PARAMETER, PUBLIC   :: Fmt_Write1      = "('ERROR IN Write STATEMENT. Unit IS : ',I5,' ERROR NUMBER IS : ', I5 )" ;
Character(143), PARAMETER, PUBLIC   :: Fmt_Element1    = "('Error in the element type. Either there is a mistake in the input file for element type or element type in not available in the code yet.')" ;
Character(144), PARAMETER, PUBLIC   :: Fmt_Element2    = "('Error in the element type. This element number',I3,'is not available in the list of this code. Check the input file for element number',I19)" ;

!#Character(??), PARAMETER, PUBLIC  :: Fmt_?  = "()" ;
!! " Type OF ANALYSIS IS NOT AVAILABLE IN THE ALLOCATION SELECT CASE - CHECK THE INPUT FILE" ;  MATCHING THE DEFINED ANALYSIS Type !??

! =========================== Unit NUMBERS OF EXTERNAL FILES ========================================================================================
! Address file
Integer (Kind=Smll), PARAMETER, PUBLIC  :: Un_ADR          = 500  ;            ! Unit number of address file to save the Name of input file and directories (.txt)

! Input files
Integer (Kind=Smll), PARAMETER, PUBLIC  :: UnInptMdl       = 501  ;            ! the Unit number of the Input file (.txt for serial and PIC code - .data for PTS code)
Integer (Kind=Smll), PARAMETER, PUBLIC  :: UnInptXYZ       = 502  ;            ! the Unit number of the Input file for node coordinates (.XYZ)
Integer (Kind=Smll), PARAMETER, PUBLIC  :: UnInptCnn       = 503  ;            ! the Unit number of the Input file for connectivities of elements (.Cnn)
Integer (Kind=Smll), PARAMETER, PUBLIC  :: UnInptCnt       = 504  ;            ! the Unit number of the Input file for node constraints (.Cnt)
Integer (Kind=Smll), PARAMETER, PUBLIC  :: UnInptStiff     = 505  ;            ! the Unit number of the input file for Number of Non-Zero entries of Stiffness matrix of PIC (.NNZStiff)
Integer (Kind=Smll), PARAMETER, PUBLIC  :: UnInptDamp      = 506  ;            ! the Unit number of the input file for Number of Non-Zero entries of Damp matrix of PIC (.NNZDamp)
Integer (Kind=Smll), PARAMETER, PUBLIC  :: UnInptMass      = 507  ;            ! the Unit number of the input file for Number of Non-Zero entries of Mass matrix of PIC  (.NNZMass)
Integer (Kind=Smll), PARAMETER, PUBLIC  :: UnInptApp       = 508  ;            ! the Unit number of the input file for Application Ordering of PIC  (.App)
Integer (Kind=Smll), PARAMETER, PUBLIC  :: UnInptDRM       = 509  ;            ! the Unit number of the input file for DRM nodes  (.DRM)
Integer (Kind=Smll), PARAMETER, PUBLIC  :: UnInptAna       = 510  ;            ! the Unit number of the Input file (.txt for serial and PIC code - .data for PTS code)
Integer (Kind=Smll), PARAMETER, PUBLIC  :: UnInptMat       = 511  ;            ! the Unit number of the Input file for material property (.Mat)

!Integer (Kind=Smll), PARAMETER, PUBLIC  :: UnInptDRMBnd    = 509  ;            ! the Unit number of the input file for DRM nodes  (.DRM)
!Integer (Kind=Smll), PARAMETER, PUBLIC  :: UnInptDRMLyr    = 510  ;            ! the Unit number of the input file for DRM nodes  (.DRM)
!Integer (Kind=Smll), PARAMETER, PUBLIC  :: UnInptElP       = 511  ;            ! the Unit number of the input file for element partitioning (.ElP)
!Integer (Kind=Smll), PARAMETER, PUBLIC  :: UnInptLPN       = 512  ;            ! the Unit number of the input file for Local PETSc Numbering (.LPN)
!Integer (Kind=Smll), PARAMETER, PUBLIC  :: UnInptNodes     = 513  ;            ! the Unit number of the input file for Nods (.Nodes)
!Integer (Kind=Smll), PARAMETER, PUBLIC  :: UnInptIDP       = 514  ;            ! the Unit number of the input file for ID PESTc (.IDP)
!Integer (Kind=Smll), PARAMETER, PUBLIC  :: UnInptGPN       = 515  ;            ! the Unit number of the input file for Golobal PESTc NUmbering(.GPN)
!Integer (Kind=Smll), PARAMETER, PUBLIC  :: UnInptRank      = 516  ;            ! the Unit number of the input file for Rank (.Rank)

! Debuggin files
Integer (Kind=Smll), PARAMETER, PUBLIC  :: Un_CHK          = 599  ;            ! Unit number of a scratch file for debuging (.Chk)

!Output files 
Integer (Kind=Smll), PARAMETER, PUBLIC  :: UnInf           = 600  ;            ! the Unit number of the information file (.inf)
Integer (Kind=Smll), PARAMETER, PUBLIC  :: Un_Out          = 601  ;            ! the Unit number of the output file for basic and general data of PIC (.data)
Integer (Kind=Smll), PARAMETER, PUBLIC  :: Un_OutXYZ       = 602  ;            ! the Unit number of the output file for coordinates of PIC  (.XYZ)
Integer (Kind=Smll), PARAMETER, PUBLIC  :: Un_OutCnn       = 603  ;            ! the Unit number of the output file for connectivities of PIC  (.Cnn)
Integer (Kind=Smll), PARAMETER, PUBLIC  :: Un_OutCnt       = 604  ;            ! the Unit number of the output file for constraints of PIC  (.Cnt)
Integer (Kind=Smll), PARAMETER, PUBLIC  :: Un_OutStiff     = 605  ;            ! the Unit number of the output file for Number of Non-Zero entries of Stiffness matrix of PIC (.NNZStiff)
Integer (Kind=Smll), PARAMETER, PUBLIC  :: Un_OutDamp      = 606  ;            ! the Unit number of the output file for Number of Non-Zero entries of Damp matrix of PIC (.NNZDamp)
Integer (Kind=Smll), PARAMETER, PUBLIC  :: Un_OutMass      = 607  ;            ! the Unit number of the output file for Number of Non-Zero entries of Mass matrix of PIC  (.NNZMass)
Integer (Kind=Smll), PARAMETER, PUBLIC  :: Un_OutApp       = 608  ;            ! the Unit number of the output file for Application Ordering of PIC  (.App)
Integer (Kind=Smll), PARAMETER, PUBLIC  :: Un_OutDRMBnd    = 609  ;            ! the Unit number of the output file for Application Ordering of PIC  (.App)
Integer (Kind=Smll), PARAMETER, PUBLIC  :: Un_OutDRMLyr    = 610  ;            ! the Unit number of the output file for Application Ordering of PIC  (.App)
Integer (Kind=Smll), PARAMETER, PUBLIC  :: Un_HisD         = 611  ;            ! the Unit number of the output file for writing down history of Displacements and stresses of nodes (.HisD)
Integer (Kind=Smll), PARAMETER, PUBLIC  :: Un_HisV         = 612  ;            ! the Unit number of the output file for writing down history of velocity of nodes (.HisV)
Integer (Kind=Smll), PARAMETER, PUBLIC  :: Un_HisA         = 613  ;            ! the Unit number of the output file for writing down history of acceleration of nodes (.HisA)
Integer (Kind=Smll), PARAMETER, PUBLIC  :: Un_Engy         = 614  ;            ! the Unit number of the output file for writing down total energy of the regular domain (.Enr)
Integer (Kind=Smll), PARAMETER, PUBLIC  :: Un_PView        = 615  ;            ! the Unit number of the output file for the ParaView visualizer, this is just for unstructured mesh (.vtu)
Integer (Kind=Smll), PARAMETER, PUBLIC  :: Un_OutMat       = 616  ;            ! the Unit number of the Input file for material property (.Mat)
Integer (Kind=Smll), PARAMETER, PUBLIC  :: Un_OutPara      = 617  ;            ! the Unit number of the Output file for Paraview (.xmf)
Integer (Kind=Smll), PARAMETER, PUBLIC  :: Un_OutWrp       = 618  ;            ! the Unit number of the wraper file for Paraview (.xmf)
Integer (Kind=Smll), PARAMETER, PUBLIC  :: Un_OutallWrp    = 618  ;            ! the Unit number of the wraper file for all Paraview data files(.xmf)
Integer (Kind=Smll), PARAMETER, PUBLIC  :: UN_Out_List_sensor = 701 ;          ! the Unit number of the file containing sensor nodes on all ranks

!Integer (Kind=Smll), PARAMETER, PUBLIC  :: Un_OutElP       = 616  ;            ! the Unit number of the output file for element partitioning (.ElP)
!Integer (Kind=Smll), PARAMETER, PUBLIC  :: Un_OutLPN       = 617  ;            ! the Unit number of the output file for Local PETSc Numbering (.LPN)
!Integer (Kind=Smll), PARAMETER, PUBLIC  :: Un_OutNodes     = 618  ;            ! the Unit number of the output file for Nods (.Nodes)
!Integer (Kind=Smll), PARAMETER, PUBLIC  :: Un_OutIDP       = 619  ;            ! the Unit number of the output file for ID PESTc (.IDP)
!Integer (Kind=Smll), PARAMETER, PUBLIC  :: Un_OutGPN       = 620  ;            ! the Unit number of the output file for Golobal PESTc NUmbering (.GPN)
!Integer (Kind=Smll), PARAMETER, PUBLIC  :: Un_OutRank      = 621  ;            ! the Unit number of the input file for Rank (.Rank)

! =========================== ANALYSIS CASE NUMBERS =================================================================================================
!#Integer (Kind=Smll), PARAMETER, PUBLIC  :: ACN_LI/NL_ST/DYN/STDYN_FREQ/TIME_DIR/MOD_SLD/FLD/SFLD_2D/3D/23D  = ?? ;                                   ! Analysis Code Number

! 1  - 20     MISC
! 21 - 40     LI ST
Integer (Kind=Smll), PARAMETER, PUBLIC  :: ACN_LI_ST                       = 21 ;      ! Linear Static Analysis - Solid Elements 2D and 3D
Integer (Kind=Smll), PARAMETER, PUBLIC  :: ACN_LI_ST_Spec                  = 22 ;      ! Linear Static Analysis - Solid Elements 2D and 3D - Spectral Element

! 41 - 60     NL ST 
! 61 - 80     LI ST DYN 
Integer (Kind=Smll), PARAMETER, PUBLIC  :: ACN_FEM_Implicit_2D             = 61 ;      ! FEM Implicit 2D
Integer (Kind=Smll), PARAMETER, PUBLIC  :: ACN_Spec_Implicit_2D            = 62 ;      ! Spec Implicit 2D
Integer (Kind=Smll), PARAMETER, PUBLIC  :: ACN_Spec_Explicit_2D            = 63 ;      ! Spec Explicit 2D

Integer (Kind=Smll), PARAMETER, PUBLIC  :: ACN_Spec_Implicit_2D_DRM        = 64 ;      ! Spec Implicit 2D DRM
Integer (Kind=Smll), PARAMETER, PUBLIC  :: ACN_FEM_Implicit_2D_DRM         = 65 ;      ! FEM Implicit 2D DRM
Integer (Kind=Smll), PARAMETER, PUBLIC  :: ACN_Spec_Explicit_2D_DRM        = 66 ;      ! Spec Explicit 2D DRM

Integer (Kind=Smll), PARAMETER, PUBLIC  :: ACN_FEM_Implicit_3D             = 67 ;      ! FEM Implicit 3D
Integer (Kind=Smll), PARAMETER, PUBLIC  :: ACN_Spec_Implicit_3D            = 68 ;      ! Spec Implicit 3D
Integer (Kind=Smll), PARAMETER, PUBLIC  :: ACN_Spec_Explicit_3D            = 69 ;      ! Spec Explicit 3D

Integer (Kind=Smll), PARAMETER, PUBLIC  :: ACN_Spec_Implicit_3D_DRM        = 70 ;      ! Spec Implicit 3D DRM
Integer (Kind=Smll), PARAMETER, PUBLIC  :: ACN_FEM_Implicit_3D_DRM         = 71 ;      ! FEM Implicit 3D DRM
Integer (Kind=Smll), PARAMETER, PUBLIC  :: ACN_Spec_Explicit_3D_DRM        = 72 ;      ! Spec Explicit 3D DRM


Integer (Kind=Smll), PARAMETER, PUBLIC  :: ACN_Analytical_DRM              = 73 ;      ! LINEAR DYNAMICS ANALYSIS - TIME DOMAIN - DIRECT APPROACH - SOLID ELEMENTS - PML TRUNCATED BOUNDARIES - FULL MATRICES - DRM

!___________________________________  Delete later !?? #####################################

Integer (Kind=Smll), PARAMETER, PUBLIC  :: ACN_LI_STDYN_TIME_DIR_SLD_23D   = 61 ;      ! LINEAR STATIC AND DYNAMICS ANALYSIS - TIME DOMAIN - DIRECT APPROACH - SOLID ELEMENTS (EMPTY RESERVOIR)
Integer (Kind=Smll), PARAMETER, PUBLIC  :: ACN_LI_STDYN_TIME_DIR_SFLD_23D  = 62 ;      ! LINEAR STATIC AND DYNAMICS ANALYSIS - TIME DOMAIN - DIRECT APPROACH (PSEUDO-SYMMETRIC)- SOLID AND FLUID ELEMENTS (FULL RESERVOIR)
Integer (Kind=Smll), PARAMETER, PUBLIC  :: ACN_LD_3                        = 63 ;      ! LINEAR DYNAMICS ANALYSIS - TIME DOMAIN - DIRECT APPROACH - SOLID ELEMENTS - PML TRUNCATED BOUNDARIES - FULL MATRICES
Integer (Kind=Smll), PARAMETER, PUBLIC  :: ACN_LD_4                        = 64 ;      ! LINEAR DYNAMICS ANALYSIS - TIME DOMAIN - DIRECT APPROACH - SOLID ELEMENTS - PML TRUNCATED BOUNDARIES - PETSC - SYMMETRIC APPROACH
Integer (Kind=Smll), PARAMETER, PUBLIC  :: ACN_LI_DYN_TIME_DIR_SLD_23D_DRM = 65 ;      ! LINEAR DYNAMICS ANALYSIS - TIME DOMAIN - DIRECT APPROACH - SOLID ELEMENTS - PML TRUNCATED BOUNDARIES - FULL MATRICES - DRM
Integer (Kind=Smll), PARAMETER, PUBLIC  :: ACN_LD_5                        = 67 ;      ! LINEAR DYNAMICS ANALYSIS - TIME DOMAIN - DIRECT APPROACH - SOLID ELEMENTS - MPML TRUNCATED BOUNDARIES - FULL MATRICES - Explicit Method
Integer (Kind=Smll), PARAMETER, PUBLIC  :: ACN_LD_6                        = 68 ;      ! LINEAR DYNAMICS ANALYSIS - TIME DOMAIN - DIRECT APPROACH - SOLID ELEMENTS - MPML TRUNCATED BOUNDARIES - FULL MATRICES - Implicit Method
Integer (Kind=Smll), PARAMETER, PUBLIC  :: ACN_LD_7                        = 69 ;      ! LINEAR DYNAMICS ANALYSIS - TIME DOMAIN - DIRECT APPROACH - SOLID ELEMENTS - MPML TRUNCATED BOUNDARIES - FULL MATRICES - Explicit - DRM
Integer (Kind=Smll), PARAMETER, PUBLIC  :: ACN_LD_8                        = 70 ;      ! LINEAR DYNAMICS ANALYSIS - TIME DOMAIN - DIRECT APPROACH - SOLID ELEMENTS - MPML TRUNCATED BOUNDARIES - Explicit Method - 3D
Integer (Kind=Smll), PARAMETER, PUBLIC  :: ACN_LD_9                        = 71 ;      ! LINEAR DYNAMICS ANALYSIS - TIME DOMAIN - DIRECT APPROACH - SOLID ELEMENTS - MPML TRUNCATED BOUNDARIES - Explicit Method - 3D -DRM
Integer (Kind=Smll), PARAMETER, PUBLIC  :: ACN_LD_10                       = 72 ;      ! LINEAR DYNAMICS ANALYSIS - TIME DOMAIN - DIRECT APPROACH - SOLID ELEMENTS - PML TRUNCATED BOUNDARIES - Explicit Method - 2D - SH waves


! 81 - 100    NL ST DYN 
! 101 - 120   LI DYN
! 121 - 140   NL DYN

! =========================== Element Types Number ==================================================================================================
! Solid Element
Integer (Kind=Tiny), PARAMETER, PUBLIC  :: El2d4NSldPS   = 1 ;      ! Element: 2D -  4 noded -  Solid - PLANE STress  ( quadrilateral )
Integer (Kind=Tiny), PARAMETER, PUBLIC  :: El2d8NSldPS   = 2 ;      ! Element: 2D -  8 noded -  Solid - PLANE STress  ( quadrilateral )
Integer (Kind=Tiny), PARAMETER, PUBLIC  :: El2d4NSldPN   = 3 ;      ! Element: 2D -  4 noded -  Solid - PLANE STrain  ( quadrilateral )
Integer (Kind=Tiny), PARAMETER, PUBLIC  :: El2d8NSldPN   = 4 ;      ! Element: 2D -  8 noded -  Solid - PLANE STrain  ( quadrilateral )
Integer (Kind=Tiny), PARAMETER, PUBLIC  :: El2d9NSldPS   = 15 ;     ! Element: 2D -  9 noded -  Solid - PLANE STress  ( quadrilateral )
Integer (Kind=Tiny), PARAMETER, PUBLIC  :: El2d9NSldPN   = 16 ;     ! Element: 2D -  9 noded -  Solid - PLANE STrain  ( quadrilateral )

Integer (Kind=Tiny), PARAMETER, PUBLIC  :: El2d3NSldPS   = 9 ;      ! Element: 2D -  3 noded -  Solid - PLANE STress  ( triangle )
Integer (Kind=Tiny), PARAMETER, PUBLIC  :: El2d6NSldPS   = 10 ;     ! Element: 2D -  6 noded -  Solid - PLANE STress  ( triangle )
Integer (Kind=Tiny), PARAMETER, PUBLIC  :: El2d3NSldPN   = 11 ;     ! Element: 2D -  3 noded -  Solid - PLANE STrain  ( triangle )
Integer (Kind=Tiny), PARAMETER, PUBLIC  :: El2d6NSldPN   = 12 ;     ! Element: 2D -  6 noded -  Solid - PLANE STrain  ( triangle )
Integer (Kind=Tiny), PARAMETER, PUBLIC  :: El2d7NSldPS   = 17 ;     ! Element: 2D -  7 noded -  Solid - PLANE STress  ( triangle )
Integer (Kind=Tiny), PARAMETER, PUBLIC  :: El2d7NSldPN   = 18 ;     ! Element: 2D -  7 noded -  Solid - PLANE STrain  ( triangle )

Integer (Kind=Tiny), PARAMETER, PUBLIC  :: El3d8NSld     = 5 ;      ! Element: 3D -  8 noded -  Solid  ( Hexahedral )
Integer (Kind=Tiny), PARAMETER, PUBLIC  :: El3d20NSld    = 6 ;      ! Element: 3D - 20 noded -  Solid  ( Hexahedral )
Integer (Kind=Tiny), PARAMETER, PUBLIC  :: El3d27NSld    = 19 ;     ! Element: 3D - 27 noded -  Solid  ( Hexahedral )

Integer (Kind=Tiny), PARAMETER, PUBLIC  :: El3d4NSld     = 13 ;     ! Element: 3D -  4 noded -  Solid  (  )
Integer (Kind=Tiny), PARAMETER, PUBLIC  :: El3d10NSld    = 14 ;     ! Element: 3D - 10 noded -  Solid  (  )
Integer (Kind=Tiny), PARAMETER, PUBLIC  :: El3d15NSld    = 20 ;     ! Element: 3D - 15 noded -  Solid  (  )


! PML Elements
Integer (Kind=Tiny), PARAMETER, PUBLIC  :: El2d4NPMLPS   = 101 ;    ! Element: 2D -  4 noded -  PML - PLANE STress  ( quadrilateral )
Integer (Kind=Tiny), PARAMETER, PUBLIC  :: El2d8NPMLPS   = 102 ;    ! Element: 2D -  8 noded -  PML - PLANE STress  ( quadrilateral )
Integer (Kind=Tiny), PARAMETER, PUBLIC  :: El2d4NPMLPN   = 103 ;    ! Element: 2D -  4 noded -  PML - PLANE STrain  ( quadrilateral )
Integer (Kind=Tiny), PARAMETER, PUBLIC  :: El2d8NPMLPN   = 104 ;    ! Element: 2D -  8 noded -  PML - PLANE STrain  ( quadrilateral )
Integer (Kind=Tiny), PARAMETER, PUBLIC  :: El2d9NPMLPS   = 115 ;    ! Element: 2D -  9 noded -  PML - PLANE STress  ( quadrilateral )
Integer (Kind=Tiny), PARAMETER, PUBLIC  :: El2d9NPMLPN   = 116 ;    ! Element: 2D -  9 noded -  PML - PLANE STrain  ( quadrilateral )

Integer (Kind=Tiny), PARAMETER, PUBLIC  :: El2d3NPMLPS   = 109 ;    ! Element: 2D -  3 noded -  PML - PLANE STress  ( triangle )
Integer (Kind=Tiny), PARAMETER, PUBLIC  :: El2d6NPMLPS   = 110 ;    ! Element: 2D -  6 noded -  PML - PLANE STress  ( triangle )
Integer (Kind=Tiny), PARAMETER, PUBLIC  :: El2d3NPMLPN   = 111 ;    ! Element: 2D -  3 noded -  PML - PLANE STrain  ( triangle )
Integer (Kind=Tiny), PARAMETER, PUBLIC  :: El2d6NPMLPN   = 112 ;    ! Element: 2D -  6 noded -  PML - PLANE STrain  ( triangle )
Integer (Kind=Tiny), PARAMETER, PUBLIC  :: El2d7NPMLPS   = 117 ;    ! Element: 2D -  7 noded -  PML - PLANE STress  ( triangle )
Integer (Kind=Tiny), PARAMETER, PUBLIC  :: El2d7NPMLPN   = 118 ;    ! Element: 2D -  7 noded -  PML - PLANE STrain  ( triangle )

Integer (Kind=Tiny), PARAMETER, PUBLIC  :: El3d8NPML     = 105 ;    ! Element: 3D -  8 noded -  PML  ( Hexahedral )
Integer (Kind=Tiny), PARAMETER, PUBLIC  :: El3d20NPML    = 106 ;    ! Element: 3D - 20 noded -  PML  ( Hexahedral )
Integer (Kind=Tiny), PARAMETER, PUBLIC  :: El3d27NPML    = 119 ;    ! Element: 3D - 27 noded -  PML  ( Hexahedral )

Integer (Kind=Tiny), PARAMETER, PUBLIC  :: El3d4NPML     = 113 ;    ! Element: 3D -  4 noded -  PML  (  )
Integer (Kind=Tiny), PARAMETER, PUBLIC  :: El3d10NPML    = 114 ;    ! Element: 3D - 10 noded -  PML  (  )
Integer (Kind=Tiny), PARAMETER, PUBLIC  :: El3d15NPML    = 120 ;    ! Element: 3D - 15 noded -  PML  (  )


! Interface Elements
Integer (Kind=Tiny), PARAMETER, PUBLIC  :: El2d6NInt     = 7 ;      ! Element: 2D -  6 noded -  Interface  ( Line )
Integer (Kind=Tiny), PARAMETER, PUBLIC  :: El3d16NInt    = 8 ;      ! Element: 3D - 16 noded -  Interface  ( Surface )


! Fluid Elements
Integer (Kind=Tiny), PARAMETER, PUBLIC  :: El2d4NFld     = 51 ;     ! Element: 2D -  4 noded -  Fluid  ( quadrilateral )
Integer (Kind=Tiny), PARAMETER, PUBLIC  :: El2d8NFld     = 52 ;     ! Element: 2D -  8 noded -  Fluid  ( quadrilateral )
Integer (Kind=Tiny), PARAMETER, PUBLIC  :: El3d8NFld     = 53 ;     ! Element: 3D -  8 noded -  Fluid  ( hexahedral )
Integer (Kind=Tiny), PARAMETER, PUBLIC  :: El3d20NFld    = 54 ;     ! Element: 3D - 20 noded -  Fluid  ( hexahedral )


! Spectral Element
! 2D - quadrilateral
Integer (Kind=Tiny), PARAMETER, PUBLIC  :: SEl2d9NSldPN  = 31 ;     ! Element: 2D -  9 node -  Solid  ( quadrilateral ) - Spectral Element - Plane Strain
Integer (Kind=Tiny), PARAMETER, PUBLIC  :: SEl2d9NPMLPN  = 32 ;     ! Element: 2D -  9 node -  PML    ( quadrilateral ) - Spectral Element - Plane Strain
Integer (Kind=Tiny), PARAMETER, PUBLIC  :: SEl2d9NMPMLPN = 33 ;     ! Element: 2D -  9 node -  MPML   ( quadrilateral ) - Spectral Element - Plane Strain
Integer (Kind=Tiny), PARAMETER, PUBLIC  :: SEl2d9NSldSH  = 38 ;     ! Element: 2D -  9 node -  Solid  ( quadrilateral ) - Spectral Element - SH waves
Integer (Kind=Tiny), PARAMETER, PUBLIC  :: SEl2d9NPMLSH  = 39 ;     ! Element: 2D -  9 node -  PML    ( quadrilateral ) - Spectral Element - SH waves
Integer (Kind=Tiny), PARAMETER, PUBLIC  :: SEl2d9NSldPS  = 42 ;     ! Element: 2D -  9 node -  Solid  ( quadrilateral ) - Spectral Element - Plane Stress
Integer (Kind=Tiny), PARAMETER, PUBLIC  :: SEl2d9NPMLPS  = 43 ;     ! Element: 2D -  9 node -  PML    ( quadrilateral ) - Spectral Element - Plane Stress
Integer (Kind=Tiny), PARAMETER, PUBLIC  :: SEl2d9NMPMLPS = 44 ;     ! Element: 2D -  9 node -  MPML   ( quadrilateral ) - Spectral Element - Plane Stress

! 2D - triangle
Integer (Kind=Tiny), PARAMETER, PUBLIC  :: SEl2d7NSldPN  = 36 ;     ! Element: 2D -  7 node -  Solid  ( triangle ) - Spectral Element -plane strain
Integer (Kind=Tiny), PARAMETER, PUBLIC  :: SEl2d7NMPMLPN = 37 ;     ! Element: 2D -  7 node -  MPML   ( triangle ) - Spectral Element -plane strain
Integer (Kind=Tiny), PARAMETER, PUBLIC  :: SEl2d7NSldPS  = 40 ;     ! Element: 2D -  7 node -  Solid  ( triangle ) - Spectral Element -plane stress
Integer (Kind=Tiny), PARAMETER, PUBLIC  :: SEl2d7NMPMLPS = 41 ;     ! Element: 2D -  7 node -  MPML   ( triangle ) - Spectral Element -plane stress

! 3D
Integer (Kind=Tiny), PARAMETER, PUBLIC  :: SEl3d27NSld   = 34 ;     ! Element: 3D - 27 node -  Solid  ( quadrilateral ) - Spectral Element
Integer (Kind=Tiny), PARAMETER, PUBLIC  :: SEl3d27NMPML  = 35 ;     ! Element: 3D - 27 node -  MPML   ( quadrilateral ) - Spectral Element
Integer (Kind=Tiny), PARAMETER, PUBLIC  :: SEl3d8NSld    = 45 ;     ! Element: 3D - 8  node -  Solid  ( quadrilateral ) - Spectral Element - Linear elements
Integer (Kind=Tiny), PARAMETER, PUBLIC  :: SEl3d8NMPML   = 46 ;     ! Element: 3D - 8  node -  MPML   ( quadrilateral ) - Spectral Element - Linear elements

! ===================================================================================================================================================

Contains ;

! THIS FUNCTION CALCULATES THE NUMERICAL VALUES FOR ABSCISSAE AND WEIGHTS OF GAUSSIAN QUADRATURE.
FUNCTION GAUSS_POINTS( NInt, NInt_Type ) ;

Implicit None ;

! =========================== LOCAL Variables =======================================================================================================
Integer (Kind=Tiny), Intent(In) :: NInt, NInt_Type ;

Type (GAUSS)  :: GAUSS_POINTS ;

! =========================== FUNCTION CODE =========================================================================================================

  SELECT CASE ( NInt_Type ) ;
    CASE ( 1_Smll ) ;    ! USING EXACT RELATIONS

      Allocate ( GAUSS_POINTS%XINT ( NInt ), GAUSS_POINTS%WINT ( NInt ) ) ;

      SELECT CASE ( NInt ) ;
        CASE (1_Smll) ;           ! polynomial degree 1
          GAUSS_POINTS%XINT = (/  0.0_DBL  /) ; ! ABSCISSAE
          GAUSS_POINTS%WINT = (/ +2.0_DBL  /) ; ! WEIGHTS

        CASE (2_Smll) ;           ! polynomial degree 3
          GAUSS_POINTS%XINT = (/  -DSQRT(1._DBL/3._DBL), +DSQRT(1._DBL/3._DBL)  /) ; ! ABSCISSAE
          GAUSS_POINTS%WINT = (/  +1.0_DBL             , +1.0_DBL               /) ; ! WEIGHTS
        
        CASE (3_Smll) ;           ! polynomial degree 5
          GAUSS_POINTS%XINT = (/  -DSQRT(3._DBL/5._DBL), 0._DBL       , +DSQRT(3._DBL/5._DBL)  /) ; ! ABSCISSAE
          GAUSS_POINTS%WINT = (/   5._DBL/9._DBL       , 8._DBL/9._DBL,  5._DBL/9._DBL         /) ; ! WEIGHTS

        CASE (4_Smll) ;           ! polynomial degree 7
          GAUSS_POINTS%XINT = (/  -DSQRT((3._DBL+2._DBL*DSQRT(6._DBL/5._DBL))/7._DBL), -DSQRT((3._DBL-2._DBL*DSQRT(6._DBL/5._DBL))/7._DBL), +DSQRT((3._DBL-2._DBL*DSQRT(6._DBL/5._DBL))/7._DBL), +DSQRT((3._DBL+2._DBL*DSQRT(6._DBL/5._DBL))/7._DBL)  /) ; ! ABSCISSAE
          GAUSS_POINTS%WINT = (/  (18._DBL-DSQRT(30._DBL))/36._DBL                   , (18._DBL+DSQRT(30._DBL))/36._DBL                    , (18._DBL+DSQRT(30._DBL))/36._DBL                    , (18._DBL-DSQRT(30._DBL))/36._DBL                      /) ; ! WEIGHTS

        CASE (5_Smll) ;           ! polynomial degree 9
          GAUSS_POINTS%XINT = (/  -DSQRT((5._DBL+2._DBL*DSQRT(10._DBL/7._DBL))/9._DBL), -DSQRT((5._DBL-2._DBL*DSQRT(10._DBL/7._DBL))/9._DBL),  0._DBL             , +DSQRT((5._DBL-2._DBL*DSQRT(10._DBL/7._DBL))/9._DBL), DSQRT((5._DBL+2._DBL*DSQRT(10._DBL/7._DBL))/9._DBL)  /) ; ! ABSCISSAE
          GAUSS_POINTS%WINT = (/  (322._DBL-13._DBL*DSQRT(70._DBL))/900._DBL          , (322._DBL+13._DBL*DSQRT(70._DBL))/900._DBL          ,  (128._DBL/225._DBL), (322._DBL+13._DBL*DSQRT(70._DBL))/900._DBL          , (322._DBL-13._DBL*DSQRT(70._DBL))/900._DBL            /) ; ! WEIGHTS

        CASE (6_Smll) ;           ! polynomial degree 11
          !GAUSS_POINTS%XINT = (/    /) ; ! ABSCISSAE
          !GAUSS_POINTS%WINT = (/    /) ; ! WEIGHTS

        CASE (7_Smll) ;           ! polynomial degree 13
          !GAUSS_POINTS%XINT = (/    /) ; ! ABSCISSAE
          !GAUSS_POINTS%WINT = (/    /) ; ! WEIGHTS

        CASE (8_Smll) ;           ! polynomial degree 15
          !GAUSS_POINTS%XINT = (/    /) ; ! ABSCISSAE
          !GAUSS_POINTS%WINT = (/    /) ; ! WEIGHTS

        CASE (10_Smll) ;          ! polynomial degree 19
          !GAUSS_POINTS%XINT = (/    /) ; ! ABSCISSAE
          !GAUSS_POINTS%WINT = (/    /) ; ! WEIGHTS

      End SELECT ;


    CASE ( 2_Smll ) ;   ! USING NUMERICAL REPRESENTATION

      Allocate ( GAUSS_POINTS%XINT ( NInt ), GAUSS_POINTS%WINT ( NInt ) ) ;

      SELECT CASE ( NInt ) ;
        CASE (1_Smll) ;           ! polynomial degree 1
          GAUSS_POINTS%XINT = (/  0.0_DBL  /) ; ! ABSCISSAE
          GAUSS_POINTS%WINT = (/ +2.0_DBL  /) ; ! WEIGHTS

        CASE (2_Smll) ;           ! polynomial degree 3
          GAUSS_POINTS%XINT = (/  -0.5773502691896257645091488_DBL, +0.5773502691896257645091488_DBL  /) ; ! ABSCISSAE
          GAUSS_POINTS%WINT = (/  +1.0_DBL                        , +1.0_DBL                          /) ; ! WEIGHTS
        
        CASE (3_Smll) ;           ! polynomial degree 5
          GAUSS_POINTS%XINT = (/  -0.7745966692414833770358531_DBL,  0.0_DBL                        , +0.7745966692414833770358531_DBL  /) ; ! ABSCISSAE
          GAUSS_POINTS%WINT = (/  +0.5555555555555555555555556_DBL, +0.8888888888888888888888889_DBL, +0.5555555555555555555555556_DBL  /) ; ! WEIGHTS

        CASE (4_Smll) ;           ! polynomial degree 7
          GAUSS_POINTS%XINT = (/  -0.8611363115940525752239465_DBL, -0.3399810435848562648026658_DBL, +0.3399810435848562648026658_DBL, +0.8611363115940525752239465_DBL  /) ; ! ABSCISSAE
          GAUSS_POINTS%WINT = (/  +0.3478548451374538573730639_DBL, +0.6521451548625461426269361_DBL, +0.6521451548625461426269361_DBL, +0.3478548451374538573730639_DBL  /) ; ! WEIGHTS

        CASE (5_Smll) ;           ! polynomial degree 9
          GAUSS_POINTS%XINT = (/  -0.9061798459386639927976269_DBL, -0.5384693101056830910363144_DBL,  0.0_DBL                        , +0.5384693101056830910363144_DBL, +0.9061798459386639927976269_DBL  /) ; ! ABSCISSAE
          GAUSS_POINTS%WINT = (/  +0.2369268850561890875142640_DBL, +0.4786286704993664680412915_DBL, +0.5688888888888888888888889_DBL, +0.4786286704993664680412915_DBL, +0.2369268850561890875142640_DBL  /) ; ! WEIGHTS

        CASE (6_Smll) ;           ! polynomial degree 11
          GAUSS_POINTS%XINT = (/  -0.9324695142031520278123016_DBL, -0.6612093864662645136613996_DBL, -0.2386191860831969086305017_DBL, +0.2386191860831969086305017_DBL, +0.6612093864662645136613996_DBL, +0.9324695142031520278123016_DBL  /) ; ! ABSCISSAE
          GAUSS_POINTS%WINT = (/  +0.1713244923791703450402961_DBL, +0.3607615730481386075698335_DBL, +0.4679139345726910473898703_DBL, +0.4679139345726910473898703_DBL, +0.3607615730481386075698335_DBL, +0.1713244923791703450402961_DBL  /) ; ! WEIGHTS

        CASE (7_Smll) ;           ! polynomial degree 13
          GAUSS_POINTS%XINT = (/  -0.9491079123427585245261897_DBL, -0.7415311855993944398638648_DBL, -0.4058451513773971669066064_DBL,  0.0_DBL                        , +0.4058451513773971669066064_DBL, +0.7415311855993944398638648_DBL, +0.9491079123427585245261897_DBL  /) ; ! ABSCISSAE
          GAUSS_POINTS%WINT = (/  +0.1294849661688696932706114_DBL, +0.2797053914892766679014678_DBL, +0.3818300505051189449503698_DBL, +0.4179591836734693877551020_DBL, +0.3818300505051189449503698_DBL, +0.2797053914892766679014678_DBL, +0.1294849661688696932706114_DBL  /) ; ! WEIGHTS

        CASE (8_Smll) ;           ! polynomial degree 15
          GAUSS_POINTS%XINT = (/  -0.9602898564975362316835609_DBL, -0.7966664774136267395915539_DBL, -0.5255324099163289858177390_DBL, -0.1834346424956498049394761_DBL, +0.1834346424956498049394761_DBL, +0.5255324099163289858177390_DBL, +0.7966664774136267395915539_DBL, +0.9602898564975362316835609_DBL  /) ; ! ABSCISSAE
          GAUSS_POINTS%WINT = (/  +0.1012285362903762591525314_DBL, +0.2223810344533744705443560_DBL, +0.3137066458778872873379622_DBL, +0.3626837833783619829651504_DBL, +0.3626837833783619829651504_DBL, +0.3137066458778872873379622_DBL, +0.2223810344533744705443560_DBL, +0.1012285362903762591525314_DBL  /) ; ! WEIGHTS

        CASE (10_Smll) ;          ! polynomial degree 19
          GAUSS_POINTS%XINT = (/  -0.9739065285171717200779640_DBL, -0.8650633666889845107320967_DBL, -0.6794095682990244062343274_DBL, -0.4333953941292471907992659_DBL, -0.1488743389816312108848260_DBL, +0.1488743389816312108848260_DBL, +0.4333953941292471907992659_DBL, +0.6794095682990244062343274_DBL, +0.8650633666889845107320967_DBL, +0.9739065285171717200779640_DBL  /) ; ! ABSCISSAE
          GAUSS_POINTS%WINT = (/  +0.0666713443086881375935688_DBL, +0.1494513491505805931457763_DBL, +0.2190863625159820439955349_DBL, +0.2692667193099963550912269_DBL, +0.2955242247147528701738930_DBL, +0.2955242247147528701738930_DBL, +0.2692667193099963550912269_DBL, +0.2190863625159820439955349_DBL, +0.1494513491505805931457763_DBL, +0.0666713443086881375935688_DBL  /) ; ! WEIGHTS

      End SELECT ;

    CASE ( 3_Smll ) ;   ! Guass Points for Triangle Elements - FEM
    ! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

      Allocate ( GAUSS_POINTS%XINT ( 2_Smll * NInt ), GAUSS_POINTS%WINT ( NInt ) ) ;

      SELECT CASE ( NInt ) ;
        CASE (1_Smll) ;           ! polynomial degree 1
          GAUSS_POINTS%XINT = (/ +1.0_DBL / +3.0_DBL,      +1.0_DBL / +3.0_DBL  /) ; ! ABSCISSAE
          GAUSS_POINTS%WINT = (/ +1.0_DBL                                       /) ; ! WEIGHTS

        CASE (3_Smll) ;           ! polynomial degree 5
          GAUSS_POINTS%XINT = (/ +2.0_DBL / +3.0_DBL, +1.0_DBL / +6.0_DBL, +1.0_DBL / +6.0_DBL,      +1.0_DBL / +6.0_DBL, +1.0_DBL / +6.0_DBL, +2.0_DBL / +3.0_DBL  /) ; ! ABSCISSAE
          GAUSS_POINTS%WINT = (/ +1.0_DBL / +3.0_DBL, +1.0_DBL / +3.0_DBL, +1.0_DBL / +3.0_DBL                                                                      /) ; ! WEIGHTS

        CASE (4_Smll) ;           ! polynomial degree 7
          GAUSS_POINTS%XINT = (/ + 1.0_DBL / + 3.0_DBL, + 3.0_DBL / + 5.0_DBL, + 1.0_DBL / + 5.0_DBL, + 1.0_DBL / + 5.0_DBL,     + 1.0_DBL / + 3.0_DBL, + 1.0_DBL / + 5.0_DBL, + 1.0_DBL / + 5.0_DBL, + 3.0_DBL / + 5.0_DBL /) ; ! ABSCISSAE
          GAUSS_POINTS%WINT = (/ -27.0_DBL / +48.0_DBL, +25.0_DBL / +48.0_DBL, +25.0_DBL / +48.0_DBL, +25.0_DBL / +48.0_DBL                                                                                                 /) ; ! WEIGHTS

        CASE (7_Smll) ;           ! polynomial degree 13

        CASE (8_Smll) ;           ! polynomial degree 15

        CASE (10_Smll) ;          ! polynomial degree 19

      End SELECT ;

    CASE ( 4_Smll ) ;   ! Gauss Points for Tetrahedra Elements

    CASE ( 5_Smll ) ;   ! Gauss-Lobatto quadrature  - Structured Mesh quadrilateral

      Allocate ( GAUSS_POINTS%XINT ( NInt ), GAUSS_POINTS%WINT ( NInt ) ) ;

      SELECT CASE ( NInt ) ;

        CASE (2_Smll) ; 
          GAUSS_POINTS%XINT = (/  -1.0_DBL,          +1.0_DBL          /) ; ! ABSCISSAE
          GAUSS_POINTS%WINT = (/  +1.0_DBL/2.0_DBL,  +1.0_DBL/2.0_DBL  /) ; ! WEIGHTS

        CASE (3_Smll) ; 
          GAUSS_POINTS%XINT = (/  -1.0_DBL,          0.0_DBL,         +1.0_DBL          /) ; ! ABSCISSAE
          GAUSS_POINTS%WINT = (/  +1.0_DBL/3.0_DBL, +4.0_DBL/3.0_DBL, +1.0_DBL/3.0_DBL  /) ; ! WEIGHTS
      End SELECT ;

    CASE ( 6_Smll ) ;   ! Gauss-Lobatto quadrature for triangle elements (spectral triangle quadrature) - 2D Unstructured mesh

      Allocate ( GAUSS_POINTS%XINT ( 2_Smll * NInt ), GAUSS_POINTS%WINT ( NInt ) ) ;

      SELECT CASE ( NInt ) ;
        CASE (7_Smll) ; 
          GAUSS_POINTS%XINT = (/ 0.0_Dbl, 1.0_Dbl, 0.0_Dbl, 0.5_Dbl, 0.5_Dbl, 0.0_Dbl, 1.0_Dbl/3.0_Dbl,      0.0_Dbl, 0.0_Dbl, 1.0_Dbl, 0.0_Dbl, 0.5_Dbl, 0.5_Dbl, 1.0_Dbl/3.0_Dbl  /) ; ! ABSCISSAE
          GAUSS_POINTS%WINT = (/ +1.0_DBL / +40.0_DBL, +1.0_DBL / +40.0_DBL, +1.0_DBL / +40.0_DBL,      +1.0_DBL / +15.0_DBL, +1.0_DBL / +15.0_DBL, +1.0_DBL / +15.0_DBL,       +9.0_DBL / +40.0_DBL    /) ; ! WEIGHTS
      End SELECT ;


    Case DeFault
      Write(*    ,*)" Type OF Integration point IS NOT AVAILABLE IN This code  - CHECK Parameteres subroutine." ;
      Write(UnInf,*)" Type OF Integration point IS NOT AVAILABLE IN This code  - CHECK Parameteres subroutine." ;
      Write(*,*) ;
      Write(UnInf,*) ;
      !#Call BEEP_FAIL ;
      Write(*, Fmt_End) ; Read(*,*) ;  STOP ;
      !Read(*," ('PRESS ENTER TO End ...') " ) ;

  End SELECT ;
!                                                                NInt = 1                                                                            NInt = 2                                                                                                                            NInt = 3                                                                                                                                                    NInt = 4                                                                                                                                                                             NInt = 5                                                                                                                                                                                                   NInt = 6                                                                                                                                                                                                                             NInt = 7                                                                                                                                                                                                                                                    NInt = 8                                                                                                                                                  
!                                                                ------------------------------------------------------------------------------      ------------------------------------------------------------------------------------------------------------------------------      ------------------------------------------------------------------------------------------------------------------------------------------------------      -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------     ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------      -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------     -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------     ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
!Real (Kind=DBL), PARAMETER, Dimension(8,8), PUBLIC  :: XINT = (/  0.0_DBL,  0.0_DBL,  0.0_DBL,  0.0_DBL,  0.0_DBL,  0.0_DBL,  0.0_DBL,  0.0_DBL,     -0.5773502691896257645091488_DBL, +0.5773502691896257645091488_DBL,  0.0_DBL,  0.0_DBL,  0.0_DBL,  0.0_DBL,  0.0_DBL,  0.0_DBL,     -0.7745966692414833770358531_DBL,  0.0_DBL                        , +0.7745966692414833770358531_DBL,  0.0_DBL,  0.0_DBL,  0.0_DBL,  0.0_DBL,  0.0_DBL,     -0.8611363115940525752239465_DBL, -0.3399810435848562648026658_DBL, +0.3399810435848562648026658_DBL, +0.8611363115940525752239465_DBL,  0.0_DBL,  0.0_DBL,  0.0_DBL,  0.0_DBL,     -0.9061798459386639927976269_DBL, -0.5384693101056830910363144_DBL,  0.0_DBL                        , +0.5384693101056830910363144_DBL, +0.9061798459386639927976269_DBL,  0.0_DBL,  0.0_DBL,  0.0_DBL,     -0.9324695142031520278123016_DBL, -0.6612093864662645136613996_DBL, -0.2386191860831969086305017_DBL, +0.2386191860831969086305017_DBL, +0.6612093864662645136613996_DBL, +0.9324695142031520278123016_DBL,  0.0_DBL,  0.0_DBL,     -0.9491079123427585245261897_DBL, -0.7415311855993944398638648_DBL, -0.4058451513773971669066064_DBL,  0.0_DBL                        , +0.4058451513773971669066064_DBL, +0.7415311855993944398638648_DBL, +0.9491079123427585245261897_DBL,  0.0_DBL,     -0.9602898564975362316835609_DBL, -0.7966664774136267395915539_DBL, -0.5255324099163289858177390_DBL, -0.1834346424956498049394761_DBL, +0.1834346424956498049394761_DBL, +0.5255324099163289858177390_DBL, +0.7966664774136267395915539_DBL, +0.9602898564975362316835609_DBL  /) ; 
!Real (Kind=DBL), PARAMETER, Dimension(8,8), PUBLIC  :: WINT = (/ +2.0_DBL,  0.0_DBL,  0.0_DBL,  0.0_DBL,  0.0_DBL,  0.0_DBL,  0.0_DBL,  0.0_DBL,     +1.0_DBL                        , +1.0_DBL                        ,  0.0_DBL,  0.0_DBL,  0.0_DBL,  0.0_DBL,  0.0_DBL,  0.0_DBL,     +0.5555555555555555555555556_DBL, +0.8888888888888888888888889_DBL, +0.5555555555555555555555556_DBL,  0.0_DBL,  0.0_DBL,  0.0_DBL,  0.0_DBL,  0.0_DBL,     +0.3478548451374538573730639_DBL, +0.6521451548625461426269361_DBL, +0.6521451548625461426269361_DBL, +0.3478548451374538573730639_DBL,  0.0_DBL,  0.0_DBL,  0.0_DBL,  0.0_DBL,     +0.2369268850561890875142640_DBL, +0.4786286704993664680412915_DBL, +0.5688888888888888888888889_DBL, +0.4786286704993664680412915_DBL, +0.2369268850561890875142640_DBL,  0.0_DBL,  0.0_DBL,  0.0_DBL,     +0.1713244923791703450402961_DBL, +0.3607615730481386075698335_DBL, +0.4679139345726910473898703_DBL, +0.4679139345726910473898703_DBL, +0.3607615730481386075698335_DBL, +0.1713244923791703450402961_DBL,  0.0_DBL,  0.0_DBL,     +0.1294849661688696932706114_DBL, +0.2797053914892766679014678_DBL, +0.3818300505051189449503698_DBL, +0.4179591836734693877551020_DBL, +0.3818300505051189449503698_DBL, +0.2797053914892766679014678_DBL, +0.1294849661688696932706114_DBL,  0.0_DBL,     +0.1012285362903762591525314_DBL, +0.2223810344533744705443560_DBL, +0.3137066458778872873379622_DBL, +0.3626837833783619829651504_DBL, +0.3626837833783619829651504_DBL, +0.3137066458778872873379622_DBL, +0.2223810344533744705443560_DBL, +0.1012285362903762591525314_DBL  /) ;
! SOURCE : http://www.holoborodko.com/pavel/numerical-methods/numerical-integration/


End FUNCTION GAUSS_POINTS ;

End Module Parameters ;

