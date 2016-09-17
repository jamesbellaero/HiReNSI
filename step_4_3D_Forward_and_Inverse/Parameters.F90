!************************************************
! PARAMETERS
!************************************************
! REVISION : M 22 Oct 2012

MODULE PARAMETERS

INTEGER         :: NJ, NDIM, NEL, NNODE, NDOF
INTEGER         :: NMAT, NPM
  
INTEGER         :: NEQ, NEQM, NEQ_DISP
  
INTEGER         :: MEMORY
  
INTEGER         :: NSTEP, NASTEP
  
REAL(8)         :: DT
REAL(8)         :: DELTA, ALPHA, ALPHA_DAMP, BETA_DAMP
REAL(8)         :: GRAVITY
  
INTEGER         :: NTRANS


! accomodate data structure for 9-noded elements
!---------- ---------- ---------- ---------- ----------
INTEGER         :: NJ_serendipity, NNODE_serendipity
INTEGER         :: NINT_Lobatto


! Weighted-regularization parameters
!---------- ---------- ---------- ---------- ----------
REAL(8)         :: w_top, w_bot, w_Length
Character (Kind = 1, Len = 3 )  :: Weighted_Regularization


! Visualization
!---------- ---------- ---------- ---------- ----------
Character (LEN=60) :: Vis_file_name

Integer   :: I_PRINT_SCREEN  = 0                 ! Screen printout


! parallel IO:
! =========================== Unit NUMBERS OF EXTERNAL FILES ========================================================================================
! Address file
Integer   :: Un_ADR          = 500  ;            ! Unit number of address file to save the Name of input file and directories (.txt)

! Input files
Integer   :: UnInptMdl       = 501  ;            ! the Unit number of the Input file (.txt for serial and PIC code - .data for PTS code)
Integer   :: UnInptXYZ       = 502  ;            ! the Unit number of the Input file for node coordinates (.XYZ)
Integer   :: UnInptCnn       = 503  ;            ! the Unit number of the Input file for connectivities of elements (.Cnn)
Integer   :: UnInptCnt       = 504  ;            ! the Unit number of the Input file for node constraints (.Cnt)
Integer   :: UnInptStiff     = 505  ;            ! the Unit number of the input file for Number of Non-Zero entries of Stiffness matrix of PIC (.NNZStiff)
Integer   :: UnInptDamp      = 506  ;            ! the Unit number of the input file for Number of Non-Zero entries of Damp matrix of PIC (.NNZDamp)
Integer   :: UnInptMass      = 507  ;            ! the Unit number of the input file for Number of Non-Zero entries of Mass matrix of PIC  (.NNZMass)
Integer   :: UnInptApp       = 508  ;            ! the Unit number of the input file for Application Ordering of PIC  (.App)
Integer   :: UnInptDRM       = 509  ;            ! the Unit number of the input file for DRM nodes  (.DRM)
Integer   :: UnInptAna       = 510  ;            ! the Unit number of the Input file (.txt for serial and PIC code - .data for PTS code)
Integer   :: UnInptMat       = 511  ;            ! the Unit number of the Input file for material property (.Mat)

! heterogeneous material properties files
Integer   :: UnInpt_Lambda   = 801  ;
Integer   :: UnInpt_Mu       = 802  ;
! heterogeneous material properties files - resume if inversion is interrupted
Integer   :: UnOut_Lambda_1 = 803  ;
Integer   :: UnOut_Mu_1     = 804  ;
Integer   :: UnOut_Lambda_2 = 805  ;
Integer   :: UnOut_Mu_2     = 806  ;

! Inversion files
Integer   :: Un_Inversion_DS = 811  ;            ! the Unit number of the Input file for Inversion Data Structure.
Integer   :: UN_Measured_HisD= 821  ;            ! the Unit number of the Input file for measured response (Displacement) at select sensor locations.
Integer   :: UN_Inversion_Iterations = 822       ! Inversion progress at each step

! Debuggin files
Integer   :: Un_CHK          = 599  ;            ! Unit number of a scratch file for debuging (.Chk)

!Output files 
Integer   :: UnInf           = 600  ;            ! the Unit number of the information file (.inf)
Integer   :: Un_Out          = 601  ;            ! the Unit number of the output file for basic and general data of PIC (.data)
Integer   :: Un_OutXYZ       = 602  ;            ! the Unit number of the output file for coordinates of PIC  (.XYZ)
Integer   :: Un_OutCnn       = 603  ;            ! the Unit number of the output file for connectivities of PIC  (.Cnn)
Integer   :: Un_OutCnt       = 604  ;            ! the Unit number of the output file for constraints of PIC  (.Cnt)
Integer   :: Un_OutStiff     = 605  ;            ! the Unit number of the output file for Number of Non-Zero entries of Stiffness matrix of PIC (.NNZStiff)
Integer   :: Un_OutDamp      = 606  ;            ! the Unit number of the output file for Number of Non-Zero entries of Damp matrix of PIC (.NNZDamp)
Integer   :: Un_OutMass      = 607  ;            ! the Unit number of the output file for Number of Non-Zero entries of Mass matrix of PIC  (.NNZMass)
Integer   :: Un_OutApp       = 608  ;            ! the Unit number of the output file for Application Ordering of PIC  (.App)
Integer   :: Un_OutDRMBnd    = 609  ;            ! the Unit number of the output file for Application Ordering of PIC  (.App)
Integer   :: Un_OutDRMLyr    = 610  ;            ! the Unit number of the output file for Application Ordering of PIC  (.App)
Integer   :: Un_HisD         = 611  ;            ! the Unit number of the output file for writing down history of Displacements and stresses of nodes (.HisD)
Integer   :: Un_HisV         = 612  ;            ! the Unit number of the output file for writing down history of velocity of nodes (.HisV)
Integer   :: Un_HisA         = 613  ;            ! the Unit number of the output file for writing down history of acceleration of nodes (.HisA)
Integer   :: Un_Engy         = 614  ;            ! the Unit number of the output file for writing down total energy of the regular domain (.Enr)
Integer   :: Un_PView        = 615  ;            ! the Unit number of the output file for the ParaView visualizer, this is just for unstructured mesh (.vtu)
Integer   :: Un_OutMat       = 616  ;            ! the Unit number of the Input file for material property (.Mat)
Integer   :: Un_OutPara      = 617  ;            ! the Unit number of the Output file for Paraview (.xmf)
Integer   :: Un_OutWrp       = 618  ;            ! the Unit number of the wraper file for Paraview (.xmf)
Integer   :: Un_OutallWrp    = 618  ;            ! the Unit number of the wraper file for all Paraview data files(.xmf)

! Material visualization
Integer   :: Un_OutMatPara   = 657  ;            ! the Unit number of the Output file for Paraview (.xmf)
Integer   :: UN_OutallWrp_Mat_Vis = 658 ;        ! the Unit number of the wraper file for all Paraview data files(.xmf)

! =========================== Type DECLARATOINS =====================================================================================================
! Holds basic parameters of each load case
Type, Public  :: BasicParam ;
  Integer , Allocatable, Dimension(:,:)  :: IntM ;   ! Holds Integer Parameters for defining the Model
  Real(8) , Allocatable, Dimension(:,:)  :: RealM ;  ! Holds Real Parameters for defining the Model
  Integer , Allocatable, Dimension(:,:)  :: IntL ;   ! Holds Integer Parameters for Loading
  Real(8) , Allocatable, Dimension(:,:)  :: RealL ;  ! Holds Real Parameters for Loading
End Type BasicParam ;

! =========================== Global Variables ======================================================================================================

! - Integer Variables -------------------------------------------------------------------------------------------------------------------------------
Integer   :: NJTotal ;                ! Total Number of Joints (Nodes) of model
Integer   :: NEQMTotal ;              ! Modified Total Number of Equatoins
Integer   :: NEQM_Mapping ;           ! Modified Number of Equations for Mapping between Application numbering and PETSc numbering. This also defines the Size of the PETSc objects on this Rank.

Integer, Dimension(5)   :: LoadC ;    ! Load Case, Indicates which kinds of static load should be consider in analysis (1 2 3 4) - 1: body force - 2: pressure - 3: concentrated force - 4: specified Displacements - 5: Static/Dynaimc analysis/both: 0/1/2

! - Logical Variable --------------------------------------------------------------------------------------------------------------------------------
Logical   :: Directory ;

! - Character Variables -----------------------------------------------------------------------------------------------------------------------------
Character ( Len = 50 ) :: ModelName ;   ! name of the model input file
Character ( Len = 50 ) :: AnaName ;     ! Name of the analysis input file
Character ( Len = 150) :: Model_InDir ; ! Directory of Model input file.
Character ( Len = 150) :: Ana_InDir ;   ! Directory of Analysis input file.
Character ( Len = 150) :: InlDir ;      ! Directory of internal files.
Character ( Len = 150) :: OutDir ;      ! Directory of output files (Results)
Character ( Len = 250) :: temp_char ;   ! Temporary character (debugging)

Character ( Len = 150) :: MatUpdateDir; ! Directory of output files (Results - Material Updates: 20-noded for Paraview).
Character ( Len = 150) :: MeasRespDir ; ! Directory of measured response at select sensor locations for computing the misfit.
Character ( Len = 150) :: MatResumeDir;
Character ( Len = 150) :: MatResumeDir_1; ! Directory of output files (Results - Material Updates in case inversion is interrupted: 27-noded). Iter: 10, 30, 50, ....
Character ( Len = 150) :: MatResumeDir_2; ! Directory of output files (Results - Material Updates in case inversion is interrupted: 27-noded). Iter: 20, 40, 60, ....
Character ( Len = 150) :: IterDir         ! Monitor iteration evolution

Character (Kind = 1, Len = 10)  :: Lambda_variable

! - Character Variables -----------------------------------------------------------------------------------------------------------------------------
Character ( Len = 20)  :: IndexRank ;   ! A variable for storing the Rank number in Character format for adding at the end of the input file Name.
Character ( Len = 20)  :: IndexSize ;   ! A variable for storing the Size number in Character format for adding at the end of the input file Name.
Character ( Len = 20)  :: Iter_begin_char ;  ! A variable for storing the iter_begin in Character format for adding at the end of the iteration monitoring file Name.
Character ( Len = 9)   :: Format_Type ;      ! determines if a file is formatted or binary

! =========================== FORMATS ===============================================================================================================
Character(87 )   :: Fmt_DATE        = "(' DATE :  ',I2.2,' - ',I2.2,' - ',I4,/,' TIME : ',I2.2,':',I2.2,':',I2.2,':',I2.2,/ )" ;
Character(27 )   :: Fmt_End         = "('PRESS ENTER TO End ...')" ;
Character(80 )   :: Fmt_ERR1_OPEN   = "( 'ERROR IN OPEN STATEMENT. Unit NUMBER = ', I3, '   ERROR NUMBER IS = ', I4  )" ;
Character(91 )   :: Fmt_ERR2_OPEN   = "('End-OF-FILE ERROR IN OPEN STATEMENT. Unit NUMBER = ', I3, '   ERROR NUMBER IS = ', I4  )" ;
Character(81 )   :: Fmt_ERR1_Close  = "( 'ERROR IN Close STATEMENT. Unit NUMBER = ', I3, '   ERROR NUMBER IS = ', I4  )" ;
Character(163)   :: Fmt_NM          = "(' FILE Name : ', A20,//,' Directories :',/, 'INPUT FILE DIRECTORY     : ', A100,/, 'OUTPUT FILES DIRECTORY   : ', A100,/, 'INTERNAL FILES DIRECTORY : ', A100,/ )" ;
Character(41 )   :: Fmt_SUC         = "('CONGRATULATIONS! DONE SUCCESSFULLY. ')" ;
Character(39 )   :: Fmt_FL          = "('OOPS!!!  FAIL TO OPERATE PROPERLY.')" ;
Character(78 )   :: Fmt_ALLCT       = "('ERROR IN ALLOCATING Arrays. ERROR NUMBER IS :', I4, '   LOCATION: ??????.')" ;
Character(80 )   :: Fmt_DEALLCT     = "('ERROR IN DEALLOCATING Arrays. ERROR NUMBER IS :', I4, '   LOCATION: ??????.')" ;
Character(23 )   :: Fmt_RUNTIME     = "(A,F50.2,'   SECONDS')" ;
Character(70 )   :: Fmt_READ1       = "('ERROR IN READ STATEMENT. Unit IS : ',I5,' ERROR NUMBER IS : ', I5 )" ;
Character(76 )   :: Fmt_READ2       = "('End OF FILE IN READ STATEMENT. Unit IS : ',I5,' ERROR NUMBER IS : ', I5 )" ;
Character(78 )   :: Fmt_READ3       = "('End OF RECORD IN READ STATEMENT. Unit IS : ',I5,' ERROR NUMBER IS : ', I5 )" ;
Character(71 )   :: Fmt_Write1      = "('ERROR IN Write STATEMENT. Unit IS : ',I5,' ERROR NUMBER IS : ', I5 )" ;
Character(143)   :: Fmt_Element1    = "('Error in the element type. Either there is a mistake in the input file for element type or element type in not available in the code yet.')" ;
Character(144)   :: Fmt_Element2    = "('Error in the element type. This element number',I3,'is not available in the list of this code. Check the input file for element number',I19)" ;

! =========================== Element Types Number ==================================================================================================
! Solid Element
Integer , PARAMETER, PUBLIC  :: El2d4NSldPS   = 1 ;      ! Element: 2D -  4 noded -  Solid - PLANE STress  ( quadrilateral )
Integer , PARAMETER, PUBLIC  :: El2d8NSldPS   = 2 ;      ! Element: 2D -  8 noded -  Solid - PLANE STress  ( quadrilateral )
Integer , PARAMETER, PUBLIC  :: El2d4NSldPN   = 3 ;      ! Element: 2D -  4 noded -  Solid - PLANE STrain  ( quadrilateral )
Integer , PARAMETER, PUBLIC  :: El2d8NSldPN   = 4 ;      ! Element: 2D -  8 noded -  Solid - PLANE STrain  ( quadrilateral )
Integer , PARAMETER, PUBLIC  :: El2d9NSldPS   = 15 ;     ! Element: 2D -  9 noded -  Solid - PLANE STress  ( quadrilateral )
Integer , PARAMETER, PUBLIC  :: El2d9NSldPN   = 16 ;     ! Element: 2D -  9 noded -  Solid - PLANE STrain  ( quadrilateral )

Integer , PARAMETER, PUBLIC  :: El2d3NSldPS   = 9 ;      ! Element: 2D -  3 noded -  Solid - PLANE STress  ( triangle )
Integer , PARAMETER, PUBLIC  :: El2d6NSldPS   = 10 ;     ! Element: 2D -  6 noded -  Solid - PLANE STress  ( triangle )
Integer , PARAMETER, PUBLIC  :: El2d3NSldPN   = 11 ;     ! Element: 2D -  3 noded -  Solid - PLANE STrain  ( triangle )
Integer , PARAMETER, PUBLIC  :: El2d6NSldPN   = 12 ;     ! Element: 2D -  6 noded -  Solid - PLANE STrain  ( triangle )
Integer , PARAMETER, PUBLIC  :: El2d7NSldPS   = 17 ;     ! Element: 2D -  7 noded -  Solid - PLANE STress  ( triangle )
Integer , PARAMETER, PUBLIC  :: El2d7NSldPN   = 18 ;     ! Element: 2D -  7 noded -  Solid - PLANE STrain  ( triangle )

Integer , PARAMETER, PUBLIC  :: El3d8NSld     = 5 ;      ! Element: 3D -  8 noded -  Solid  ( Hexahedral )
Integer , PARAMETER, PUBLIC  :: El3d20NSld    = 6 ;      ! Element: 3D - 20 noded -  Solid  ( Hexahedral )
Integer , PARAMETER, PUBLIC  :: El3d27NSld    = 19 ;     ! Element: 3D - 27 noded -  Solid  ( Hexahedral )

Integer , PARAMETER, PUBLIC  :: El3d4NSld     = 13 ;     ! Element: 3D -  4 noded -  Solid  (  )
Integer , PARAMETER, PUBLIC  :: El3d10NSld    = 14 ;     ! Element: 3D - 10 noded -  Solid  (  )
Integer , PARAMETER, PUBLIC  :: El3d15NSld    = 20 ;     ! Element: 3D - 15 noded -  Solid  (  )


! PML Elements
Integer , PARAMETER, PUBLIC  :: El2d4NPMLPS   = 101 ;    ! Element: 2D -  4 noded -  PML - PLANE STress  ( quadrilateral )
Integer , PARAMETER, PUBLIC  :: El2d8NPMLPS   = 102 ;    ! Element: 2D -  8 noded -  PML - PLANE STress  ( quadrilateral )
Integer , PARAMETER, PUBLIC  :: El2d4NPMLPN   = 103 ;    ! Element: 2D -  4 noded -  PML - PLANE STrain  ( quadrilateral )
Integer , PARAMETER, PUBLIC  :: El2d8NPMLPN   = 104 ;    ! Element: 2D -  8 noded -  PML - PLANE STrain  ( quadrilateral )
Integer , PARAMETER, PUBLIC  :: El2d9NPMLPS   = 115 ;    ! Element: 2D -  9 noded -  PML - PLANE STress  ( quadrilateral )
Integer , PARAMETER, PUBLIC  :: El2d9NPMLPN   = 116 ;    ! Element: 2D -  9 noded -  PML - PLANE STrain  ( quadrilateral )

Integer , PARAMETER, PUBLIC  :: El2d3NPMLPS   = 109 ;    ! Element: 2D -  3 noded -  PML - PLANE STress  ( triangle )
Integer , PARAMETER, PUBLIC  :: El2d6NPMLPS   = 110 ;    ! Element: 2D -  6 noded -  PML - PLANE STress  ( triangle )
Integer , PARAMETER, PUBLIC  :: El2d3NPMLPN   = 111 ;    ! Element: 2D -  3 noded -  PML - PLANE STrain  ( triangle )
Integer , PARAMETER, PUBLIC  :: El2d6NPMLPN   = 112 ;    ! Element: 2D -  6 noded -  PML - PLANE STrain  ( triangle )
Integer , PARAMETER, PUBLIC  :: El2d7NPMLPS   = 117 ;    ! Element: 2D -  7 noded -  PML - PLANE STress  ( triangle )
Integer , PARAMETER, PUBLIC  :: El2d7NPMLPN   = 118 ;    ! Element: 2D -  7 noded -  PML - PLANE STrain  ( triangle )

Integer , PARAMETER, PUBLIC  :: El3d8NPML     = 105 ;    ! Element: 3D -  8 noded -  PML  ( Hexahedral )
Integer , PARAMETER, PUBLIC  :: El3d20NPML    = 106 ;    ! Element: 3D - 20 noded -  PML  ( Hexahedral )
Integer , PARAMETER, PUBLIC  :: El3d27NPML    = 119 ;    ! Element: 3D - 27 noded -  PML  ( Hexahedral )

Integer , PARAMETER, PUBLIC  :: El3d4NPML     = 113 ;    ! Element: 3D -  4 noded -  PML  (  )
Integer , PARAMETER, PUBLIC  :: El3d10NPML    = 114 ;    ! Element: 3D - 10 noded -  PML  (  )
Integer , PARAMETER, PUBLIC  :: El3d15NPML    = 120 ;    ! Element: 3D - 15 noded -  PML  (  )


! Interface Elements
Integer , PARAMETER, PUBLIC  :: El2d6NInt     = 7 ;      ! Element: 2D -  6 noded -  Interface  ( Line )
Integer , PARAMETER, PUBLIC  :: El3d16NInt    = 8 ;      ! Element: 3D - 16 noded -  Interface  ( Surface )


! Fluid Elements
Integer , PARAMETER, PUBLIC  :: El2d4NFld     = 51 ;     ! Element: 2D -  4 noded -  Fluid  ( quadrilateral )
Integer , PARAMETER, PUBLIC  :: El2d8NFld     = 52 ;     ! Element: 2D -  8 noded -  Fluid  ( quadrilateral )
Integer , PARAMETER, PUBLIC  :: El3d8NFld     = 53 ;     ! Element: 3D -  8 noded -  Fluid  ( hexahedral )
Integer , PARAMETER, PUBLIC  :: El3d20NFld    = 54 ;     ! Element: 3D - 20 noded -  Fluid  ( hexahedral )


! Spectral Element
! 2D - quadrilateral
Integer , PARAMETER, PUBLIC  :: SEl2d9NSldPN  = 31 ;     ! Element: 2D -  9 node -  Solid  ( quadrilateral ) - Spectral Element - Plane Strain
Integer , PARAMETER, PUBLIC  :: SEl2d9NPMLPN  = 32 ;     ! Element: 2D -  9 node -  PML    ( quadrilateral ) - Spectral Element - Plane Strain
Integer , PARAMETER, PUBLIC  :: SEl2d9NMPMLPN = 33 ;     ! Element: 2D -  9 node -  MPML   ( quadrilateral ) - Spectral Element - Plane Strain
Integer , PARAMETER, PUBLIC  :: SEl2d9NSldSH  = 38 ;     ! Element: 2D -  9 node -  Solid  ( quadrilateral ) - Spectral Element - SH waves
Integer , PARAMETER, PUBLIC  :: SEl2d9NPMLSH  = 39 ;     ! Element: 2D -  9 node -  PML    ( quadrilateral ) - Spectral Element - SH waves
Integer , PARAMETER, PUBLIC  :: SEl2d9NSldPS  = 42 ;     ! Element: 2D -  9 node -  Solid  ( quadrilateral ) - Spectral Element - Plane Stress
Integer , PARAMETER, PUBLIC  :: SEl2d9NPMLPS  = 43 ;     ! Element: 2D -  9 node -  PML    ( quadrilateral ) - Spectral Element - Plane Stress
Integer , PARAMETER, PUBLIC  :: SEl2d9NMPMLPS = 44 ;     ! Element: 2D -  9 node -  MPML   ( quadrilateral ) - Spectral Element - Plane Stress

! 2D - triangle
Integer , PARAMETER, PUBLIC  :: SEl2d7NSldPN  = 36 ;     ! Element: 2D -  7 node -  Solid  ( triangle ) - Spectral Element -plane strain
Integer , PARAMETER, PUBLIC  :: SEl2d7NMPMLPN = 37 ;     ! Element: 2D -  7 node -  MPML   ( triangle ) - Spectral Element -plane strain
Integer , PARAMETER, PUBLIC  :: SEl2d7NSldPS  = 40 ;     ! Element: 2D -  7 node -  Solid  ( triangle ) - Spectral Element -plane stress
Integer , PARAMETER, PUBLIC  :: SEl2d7NMPMLPS = 41 ;     ! Element: 2D -  7 node -  MPML   ( triangle ) - Spectral Element -plane stress

! 3D
Integer , PARAMETER, PUBLIC  :: SEl3d27NSld   = 34 ;     ! Element: 3D - 27 node -  Solid  ( quadrilateral ) - Spectral Element
Integer , PARAMETER, PUBLIC  :: SEl3d27NMPML  = 35 ;     ! Element: 3D - 27 node -  MPML   ( quadrilateral ) - Spectral Element
Integer , PARAMETER, PUBLIC  :: SEl3d8NSld    = 45 ;     ! Element: 3D - 8  node -  Solid  ( quadrilateral ) - Spectral Element - Linear elements
Integer , PARAMETER, PUBLIC  :: SEl3d8NMPML   = 46 ;     ! Element: 3D - 8  node -  MPML   ( quadrilateral ) - Spectral Element - Linear elements

! ===================================================================================================================================================


END MODULE
