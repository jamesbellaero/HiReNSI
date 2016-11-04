
! ***************************************************************************************************************************************************
! NAME   : PIC  ( Partitioner and Input Creator )                                                                                                  **
!                                                                                                                                                  **
! DEVELOPER : Babak Poursartip, M.Sc.                                                                                                              **
!                                                                                                                                                  **
! The University of Texas at Austin, TX, USA                                               ADVISER : LOUKAS F. KALLIVOKAS                          **
!                                                                                                                                                  **
! CIVIL, ARCHITECTRAL AND ENVIRONMENTAL ENGINEERNIG - STRUCTURAL ENGINEERING                                                                       **
! Institue for Computational Engineering and Science                                                                                               **
!                                                                                                                                                  **
! History : START    JUL 2011                                                                                                                      **
!           V 1.0    JUL 2011                                                                                                                      **
!           V 1.1    Aug 2011                                                                                                                      **
!           V 1.2    Oct 2011     Number of non-zero entries of main PETSc objects to reduce malloc cost                                           **
!           V 1.3    May 2012     2D unstructured mesh - second oreder triangle elements                                                           **
!           V 1.4    Jun 2012     adding all elemets to the code                                                                                   **
!           V 1.5    Jun 2012     Domain Reduction Method                                                                                          **
!           V 2.0    July 2012    ParMetis                                                                                                         **
!           V 2.1    July 2012    Revise the structure of the code                                                                                 **
!           V 3.0    May  2014    Inversion data structure and material visualization
!                                                                                                                                                  **
! Last Update: 22-May-2014                                                                                                                         **
!                                                                                                                                                  **
! ***************************************************************************************************************************************************

Program PIC ;

! =========================== Libraries and Modules =================================================================================================
!Use Metis             ;
Use Parameters ;               ! Module comprising all Gauss Points, Type declarations, Constants, Definitions, ...
Use Information ;              ! Module for writing information down -  timing , model, etc.
Use ReNumbering ;              ! generates local node numbering and global PETSc numbering
Use Partition_Output ;         ! Writes down the output files
Use ConvertToMetis ;           ! Convert second oreder elements to first order for Metis
Use PETScObject_Size ;         ! Evaluates object sizes in PETSc.
Use Input_Subroutines ;        ! Contains all Input subroutines
Use Partition_Output_Binary ;  ! generates output files in binary format
Use Partition_Output_HDF5 ;
Use Heterogeneous_Material ;   ! mappings for heterogeneous materials
Use Inversion_Data_Structure ; ! data structure for inversion
Use Mat_Vis_HDF5 ;             ! material visualization

Implicit None ;

! =========================== PETSC LIBRARIES =======================================================================================================
#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscmat.h"
#include "finclude/petscvec.h90"
!#include "metis.h"

! ============================ PETSC Variables AND OBJECTS ==========================================================================================
! - PETSC OBJECTS -----------------------------------------------------------------------------------------------------------------------------------
!Vec            :: ;
!Mat            :: ;

! - PETSC INTERNAL Variables ------------------------------------------------------------------------------------------------------------------------
PetscErrorCode :: ErrPTC ;                       ! Error
!PetscTruth     :: FLG ;                         ! Flag
PetscMPIInt    :: Size ;                         ! Total number of ranks
PetscMPIInt    :: Rank ;                         ! Rank number

! - PETSC Variables ---------------------------------------------------------------------------------------------------------------------------------
!#PetscInt       ::  ;
!#PetscReal      ::  ;
!#PetscScalar    ::  ;

! - PETSC Arrays ------------------------------------------------------------------------------------------------------------------------------------
!#PetscInt       ::  ;
!#PetscReal      ::  ;
!#PetscScalar    ::  ;

! =========================== Global Variables ======================================================================================================
Include 'GlobalVariables.F90'              ! Global Variables used in all codes
Include 'PIC_GlobalVariables.F90'          ! Global Variables used sepecifically in PIC code

! ============================ START PETSC LIBRARY ===================================================================================================
Call PetscInitialize ( PETSC_NULL_Character  , ErrPTC ) ;
Call MPI_Comm_size   ( PETSC_COMM_WORLD, SIZE, ErrPTC ) ;
Call MPI_Comm_rank   ( PETSC_COMM_WORLD, RANK, ErrPTC ) ;

! =========================== TIME AND DATE =========================================================================================================
Call CPU_TIME( TimeS )  ;
Call GETDAT( Iyr, Imon, Iday ) ;
Call GETTIM( Ih, Im, Is, I100th ) ;

! =============================================== OPEN EXTERNAL FILES ===============================================================================
! - Input FILE --------------------------------------------------------------------------------------------------------------------------------------
UnFile = UN_ADR ;
Open ( Unit = UnFile, FILE = 'ADDRESS_PIC.TXT', ERR =  1001, IOSTAT = IO_File, ACCESS = 'SEQUENTIAL', ACTION = 'READ', ASYNCHRONOUS = 'NO', BLANK = 'NULL', BLOCKSIZE = 0                             , DisPOSE = 'KEEP', FORM = 'FORMATTED', POSITION = 'ASIS', STATUS = 'OLD' ) ;

! - READ THE Input FILE NAME AND DIRECTORIES FROM "ADDRESS FILE" IN THE CURRENT DIRECTORY -----------------------------------------------------------
Read(UN_ADR,*) ;
Read(UN_ADR,*) ModelNAME ; ! Name of the Input file
Read(UN_ADR,*) ;
Read(UN_ADR,*) NAT, Output_Type ;
Read(UN_ADR,*) ;
Read(UN_ADR,*) Model_InDir ; ! Direction of the Input file
Read(UN_ADR,*) ;
Read(UN_ADR,*) InlDir ; ! Direction of the internal files
Read(UN_ADR,*) ;
Read(UN_ADR,*) OutDir ; ! Direction of the output file

! - READ THE Input FILE NAME AND DIRECTORIES FROM THE MAIN Program ----------------------------------------------------------------------------------
!#NAME    ='RJ8' ;
!#Model_InDir  ='J:\Projects\Fortran\Thesis\Data\Input\RJ' ;
!#OutDir ='J:\Projects\Fortran\Thesis\Data\RESULT' ;
!#InlDir ='J:\Projects\Fortran\Thesis\Data\INTERNAL' ;

! - READ THE Input FILE NAME AND DIRECTORIES FROM THE TERMINAL (KEYBOARD) ---------------------------------------------------------------------------
!#Write(*,"(' ENTER THE NAME OF Input FILE : ',\)") ;
!#Read(*,*)NAME ;

!#Write(*,"(' ENTER THE ADDRESS OF Input FILE : ',\)") ;
!#Read(*,*)Model_InDir

!#Write(*,"(' ENTER THE ADDRESS OF OUTPUT FILES : ',\)")
!#Read(*,*)OutDir

Model_InDir = TRIM(AdjustL (Model_InDir)) ;

Write(*,*)"Model_InDir: ", Model_InDir ;
!Open ( 11, FILE = 'checkDir.txt', ERR =  1001, IOSTAT = IO_File, ACCESS = 'SEQUENTIAL', ACTION = 'READ', ASYNCHRONOUS = 'NO', BLANK = 'NULL', BLOCKSIZE = 0, DEFAULTFILE = TRIM(Model_InDir), DisPOSE = 'KEEP', FORM = 'FORMATTED', POSITION = 'ASIS', STATUS = 'unknown' ) ;

! - Input FILE --------------------------------------------------------------------------------------------------------------------------------------
UnFile = UnInptMdl ;
!Open ( Unit = UnFile, FILE = TRIM(ModelName)//'.txt', ERR =  1001, IOSTAT = IO_File, ACCESS = 'SEQUENTIAL', ACTION = 'READ', ASYNCHRONOUS = 'NO', BLANK = 'NULL', BLOCKSIZE = 0, DEFAULTFILE = TRIM(Model_InDir), DisPOSE = 'KEEP', FORM = 'FORMATTED', POSITION = 'ASIS', STATUS = 'OLD' ) ;

! - Input file for nodes' coordinates ---------------------------------------------------------------------------------------------------------------
UnFile = UnInptXYZ ;
Open ( Unit = UnFile, FILE = TRIM(ModelName)//'.XYZ', ERR =  1001, IOSTAT = IO_File, ACCESS = 'SEQUENTIAL', ACTION = 'READ', ASYNCHRONOUS = 'NO', BLANK = 'NULL', BLOCKSIZE = 0, DEFAULTFILE = TRIM(Model_InDir), DisPOSE = 'KEEP', FORM = 'FORMATTED', POSITION = 'ASIS', STATUS = 'OLD' ) ;

! - Input file for elements' connectivities ---------------------------------------------------------------------------------------------------------
UnFile = UnInptCnn ;
Open ( Unit = UnFile, FILE = TRIM(ModelName)//'.Cnn', ERR =  1001, IOSTAT = IO_File, ACCESS = 'SEQUENTIAL', ACTION = 'READ', ASYNCHRONOUS = 'NO', BLANK = 'NULL', BLOCKSIZE = 0, DEFAULTFILE = TRIM(Model_InDir), DisPOSE = 'KEEP', FORM = 'FORMATTED', POSITION = 'ASIS', STATUS = 'OLD' ) ;

! - Input file for nodes' constraints ---------------------------------------------------------------------------------------------------------------
UnFile = UnInptCnt ;
Open ( Unit = UnFile, FILE = TRIM(ModelName)//'.Cnt', ERR =  1001, IOSTAT = IO_File, ACCESS = 'SEQUENTIAL', ACTION = 'READ', ASYNCHRONOUS = 'NO', BLANK = 'NULL', BLOCKSIZE = 0, DEFAULTFILE = TRIM(Model_InDir), DisPOSE = 'KEEP', FORM = 'FORMATTED', POSITION = 'ASIS', STATUS = 'OLD' ) ;

! - Create the result folder ------------------------------------------------------------------------------------------------------------------------
Directory = MakeDirQQ (TRIM(AdjustL (OutDir))//'/'//TRIM(AdjustL (ModelName))) ;
  IF (Directory) THEN ;
     WRITE (*     ,*) 'New subdirectory successfully created' ;
     WRITE (UnInf,*) 'New subdirectory successfully created' ;
  ELSE ;
     WRITE (*     ,*) 'Failed to create subdirectory' ;
     WRITE (UnInf,*) 'Failed to create subdirectory' ;
  END IF ;

Directory = MakeDirQQ (TRIM(AdjustL (OutDir))//'/'//TRIM(AdjustL (ModelName))//'/'//'Model'   ) ;
  IF (Directory) THEN ;
     WRITE (*     ,*) 'New subdirectory successfully created' ;
     WRITE (UnInf,*) 'New subdirectory successfully created' ;
  ELSE ;
     WRITE (*     ,*) 'Failed to create subdirectory' ;
     WRITE (UnInf,*) 'Failed to create subdirectory' ;
  END IF ;

Directory = MakeDirQQ (TRIM(AdjustL (InlDir))//'/'//TRIM(AdjustL (ModelName))) ;
  IF (Directory) THEN ;
     WRITE (*     ,*) 'New subdirectory successfully created' ;
     WRITE (UnInf,*) 'New subdirectory successfully created' ;
  ELSE ;
     WRITE (*     ,*) 'Failed to create subdirectory' ;
     WRITE (UnInf,*) 'Failed to create subdirectory' ;
  END IF ;

OutDir = TRIM(AdjustL (OutDir))//'/'//'Model' ;
InlDir = TRIM(AdjustL (InlDir));

write (*,*)"Input Directory:     ",Model_InDir ;
write (*,*)"Output Directory:    ",OutDir ;
write (*,*)"Internal Directory:  ",InlDir ;

! - Informaton file ---------------------------------------------------------------------------------------------------------------------------------
UnFile = UnInf ;
Open ( Unit = UnFile, FILE = TRIM(ModelName)//'.INF', ERR =  1001, IOSTAT = IO_File, ACCESS = 'SEQUENTIAL', ACTION= 'Write', ASYNCHRONOUS = 'NO', BLANK = 'NULL', BLOCKSIZE = 0, DEFAULTFILE = TRIM(OutDir), DisPOSE = 'KEEP', FORM = 'FORMATTED', POSITION = 'ASIS', STATUS = 'REPLACE' ) ;

! - CHECK FILE --------------------------------------------------------------------------------------------------------------------------------------
UnFile = UN_CHK ;
Open ( Unit = UnFile, FILE = TRIM(ModelName)//'.Part', ERR =  1001, IOSTAT = IO_File, ACCESS = 'SEQUENTIAL', ACTION= 'Write', ASYNCHRONOUS = 'NO', BLANK = 'NULL', BLOCKSIZE = 0, DEFAULTFILE = TRIM(InlDir), DisPOSE = 'KEEP', FORM = 'FORMATTED', POSITION = 'ASIS', STATUS = 'REPLACE' ) ;

! =============================================== Write DOWN BASIC INFORMATION ======================================================================
Call INFO ( Iyr, Imon, Iday, Ih, Im, Is, I100th,     ModelName, Model_InDir, OutDir, InlDir ) ;

! =============================================== Output for different ANALYSIS =====================================================================

! Only for DRM
  If ( NAT == ACN_Spec_Implicit_2D_DRM .OR. NAT == ACN_FEM_Implicit_2D_DRM .OR. NAT == ACN_Spec_Explicit_2D_DRM .OR. NAT == ACN_Spec_Implicit_3D_DRM .OR. NAT == ACN_FEM_Implicit_3D_DRM .OR. NAT == ACN_Spec_Explicit_3D_DRM ) Then ;
    ! - Input file for DRM nodes  -----------------------------------------------------------------------------------------------------------------------
    UnFile = UnInptDRM ;
    Open ( Unit = UnFile, FILE = TRIM(ModelName)//'.DRM', ERR =  1001, IOSTAT = IO_File, ACCESS = 'SEQUENTIAL', ACTION = 'READ', ASYNCHRONOUS = 'NO', BLANK = 'NULL', BLOCKSIZE = 0, DEFAULTFILE = TRIM(Model_InDir), DisPOSE = 'KEEP', FORM = 'FORMATTED', POSITION = 'ASIS', STATUS = 'OLD' ) ;
  End If ;
! =============================================== Reading Data ======================================================================================
! BASIC DATA
Call Input (                                                                                                                    &
NDOF, MaxNNode, NDim, Matrix_Type, MetisType,                                                                                   & ! Integer (1) Variables
NGroup, NPM, NMat, NLCase,                                                                                                      & ! Integer (2) Variables
NumFlag, NParts, KWay,                                                                                                          & ! Integer (4) Variables
NEl, NJ,                                                                                                                        & ! Integer (8) Variables
!                                                                                                                               & ! Real Variables
LoadC,                                                                                                                          & ! Integer Arrays
!                                                                                                                               & ! Real Arrays
!                                                                                                                               & ! Characters
Param                                                                                                                           & ! Type
) ;

! Allocating required arrays
Allocate ( MTEL (NEl), ELT (NEl), ELGR (NEl), INOD (MaxNNode,NEl), PML_DIM (2*NDim, 2), XYZ (NJ,NDim), ID (NJ,NDOF), PMat (NMat,NPM), &
           IDBC  ( Param%IntM( 2, 4), 2*NDim + 1), &  ! Param%IntM( 2,4) = NIDBC
           JLoad ( Param%IntM( 3, 1) ), PLoad ( NDim, Param%IntM( 3,1) ), & ! Param%IntM(LoadC (3),1) = NLN
           NDAN  ( Param%IntM( 6, 1) ), NVAN ( Param%IntM( 6, 2) ), NAAN ( Param%IntM( 6, 3) ), & ! ( Param%IntM( NLCase + 2, I), I = 1, 3 ) ; ! NNDH, NNVH, NNAH
           PBLD  ( Param%IntM( 1, 1), Param%IntM( 1, 2) ), & ! Param%IntP( 1, 1), Param%IntP( 1, 2) ; ! NBLD   = NGroup, NPBL   = NDim ;
           UDis (NDOF+1, Param%IntM( 4, 1) ), & ! Param%IntM(LoadC (4),1) = NSDN
           STAT = ERR_Alloc) ;
  IF ( ERR_Alloc /= 0 ) Then ;
    Write (*, Fmt_ALLCT) ERR_Alloc ;  Write (UnInf, Fmt_ALLCT) ERR_Alloc ;
    !#Call BEEP_FAIL ;
    Write(*, Fmt_FL) ;  Write(UnInf, Fmt_FL) ; Write(*, Fmt_End) ; Read(*,*) ;  STOP ;
  End If ;

  If ( NAT == ACN_Spec_Implicit_2D_DRM .OR. NAT == ACN_FEM_Implicit_2D_DRM .OR. NAT == ACN_Spec_Explicit_2D_DRM .OR. NAT == ACN_Spec_Implicit_3D_DRM .OR. NAT == ACN_FEM_Implicit_3D_DRM .OR. NAT == ACN_Spec_Explicit_3D_DRM ) Then ;
    ! Allocating required arrays for drm analysis
    Allocate ( NoBndry_DRM( Param%IntM( 7, 1) ), NoLayer_DRM ( Param%IntM( 7, 2) ), InciWave (5),    STAT = ERR_Alloc) ;
      IF ( ERR_Alloc /= 0 ) Then ;
        Write (*, Fmt_ALLCT) ERR_Alloc ;  Write (UnInf, Fmt_ALLCT) ERR_Alloc ;
        !#Call BEEP_FAIL ;
        Write(*, Fmt_FL) ;  Write(UnInf, Fmt_FL) ; Write(*, Fmt_End) ; Read(*,*) ;  STOP ;
      End If ;
  End If ;

! READ Arrays
Call Input  (                                                                                                                   &
NDOF, NDim, MaxNNode,                                                                                                           & ! Integer (1) Variables
NPM, NMat,                                                                                                                      & ! Integer (2) Variables
!                                                                                                                               & ! Integer (4) Variables
NEl, NJ,                                                                                                                        & ! Integer (8) Variables
!                                                                                                                               & ! Real Variables
LoadC, IDBC,     MTEL, ELT, ELGR,         JLoad,  NDAN, NVAN, NAAN, NoBndry_DRM, NoLayer_DRM, INOD, ID,                         & ! Integer Arrays
Param,           PMat, PBLD, XYZ, UDis, PLoad,  PML_DIM                                                                         & ! Real Arrays
!                                                                                                                               & ! Characters
!                                                                                                                               & ! Type
) ;

! =============================================== MODIFYING InputS FOR METIS ========================================================================
! Allocating rquired arrays for mesh partitioning
Write (*,*)"Allocating arrays ..."
  If ( MetisType == 0_Tiny ) Then ;
    Allocate ( EPart ( NEl ),     STAT = ERR_Alloc ) ;
  Else If ( MetisType == 1_Tiny ) Then ;
    Allocate ( EPart ( NEl ), ElmDist ( KWay + 1 ), eptr ( NEl +1 ), eind ( MaxNNode * NEl ), ElmWgt ( NEl ), tpwgts ( 0:ncon * NParts-1 ), ubvec ( 0:ncon-1 ),      STAT = ERR_Alloc ) ;
  Else If ( MetisType == 2_Tiny ) Then ;
    !Allocate ( VWgt ( NEl ), VSize ( NEl ), EPart ( NEl ), eptr ( NEl +1 ), eind ( MaxNNode * NEl ),tpwgts ( 0:NParts-1 ), MOptions ( 1 ), NPart ( NJ ),      STAT = ERR_Alloc ) ; ! , MOptions ( METIS_NOPTIONS )
    !Allocate ( VWgt ( NEl ), VSize ( NEl ), EPart ( 0:NEl-1 ), eptr ( 0:NEl ), eind ( 0:MaxNNode * NEl-1 ),tpwgts ( 0:NParts-1 ), MOptions ( 1 ),      STAT = ERR_Alloc ) ; ! , MOptions ( METIS_NOPTIONS ) , NPart ( 0:NJ-1 )
    Allocate ( VWgt ( NEl ), VSize ( NEl ), EPart ( NEl ), eptr ( NEl+1 ), eind ( MaxNNode * NEl ),tpwgts ( 0:NParts-1 ),       STAT = ERR_Alloc ) ; ! , MOptions ( METIS_NOPTIONS ) , NPart ( 0:NJ-1 )
  End If ;
  IF ( ERR_Alloc /= 0 ) Then ;
    Write (*, Fmt_ALLCT) ERR_Alloc ;  Write (UnInf, Fmt_ALLCT) ERR_Alloc ;
    !#Call BEEP_FAIL ;
    Write(*, Fmt_FL) ;  Write(UnInf, Fmt_FL) ; Write(*, Fmt_End) ; Read(*,*) ;  STOP ;
  End If ;


! Converts the elements ordering to those that are compatible with Metis or ParMetis
  Write (*,*)"Prepare data for mesh partitioning ..."
  If ( NParts > 1_Shrt ) Then ;
    If ( MetisType == 0_Tiny ) Then ; ! Metis Version 4.0
      Call Convert         ( MaxNNode, NEl, NJ, ETypeG, NJG,     ELT, ELMNTS,     INOD ) ;
    Else If ( MetisType == 1_Tiny ) Then ; ! ParaMetis Version 3.2
      Call Convert  ( NDim, KWay, NParts, WgtFlag, NumFlag, NCommonNodes,     NEl, NJ, NJG,     ELT,      ElmWgt, ElmDist, eptr, eind, Poptions,     INod,     tpwgts, ubvec ) ;
    Else If ( MetisType == 2_Tiny ) Then ; ! Metis Version 5.1.0
      Call Convert  ( NDim,     NParts,     NCommonNodes,     NEl, NJ, Neind, NJG,     ELT,     eptr, eind, Vwgt, VSize,      INod, tpwgts ) ;
    End If ;
  End If ;

!Nullify (PVWgt, PVSize, PMoptions, Ptpwgts) ;

! =============================================== Partitionig ROUTINES (Metis or ParMetis) ==========================================================
! Mesh partitionig
  If ( NParts > 1_Shrt ) Then ;

    If ( MetisType == 0_Tiny ) Then ; ! Metis Version 4.0
      Write(*    ,*) "Partitioning in METIS version 4.0 ...", NEL, NJG, ETypeG, NParts ;
      Allocate (NPart ( NJG ))
        !Do IEl = 1,NEl ;
        !  Write(Un_CHK,"(8(I10))")(ELMNTS((IEl-1)*8+J),J=1,8)
        !End Do ;
write(*,*)'-----------------before Metis'
      Call METIS_PartMeshDual  ( NEl, NJG, ELMNTS, ETypeG, NumFlag, NParts, Edgecut, EPart , NPart  ) ;
write(*,*)'-----------------after Metis'
      Write(*    ,*) "End subroutine < Metis Version 4.0 >" ;
      Write(UnInf,*) "End subroutine < Metis Version 4.0 >" ;
    Else If ( MetisType == 1_Tiny ) Then ; ! ParMetis Version 3.2
      Call ParMetis_V3_PartMeshKway ( ElmDist, eptr, eind2, ElmWgt, Wgtflag, NumFlag, ncon, NCommonNodes, NParts, tpwgts2, ubvec, Poptions, EdgeCut, EPart, PETSC_COMM_WORLD ) ;
      write(*    ,*) "End subroutine < ParMetis Version 3.2 >"
      write(UnInf,*) "End subroutine < ParMetis Version 3.2 >"
    Else If ( MetisType == 2_Tiny ) Then ; ! Metis Version 5.1.0
      !MOptions(:) = 0 ;
      eptr (:)= eptr (:) - 1 ;
      eind (:)= eind (:) - 1 ;
!write(un_chk,*) "eind"
!!print *,4*nel, Neind
!do i = 0,Maxnnode*nel-1
!!if (eind (i) == 0 ) exit ;
!write(un_chk,*) i , eind (i);
!end do
!
!write(un_chk,*) "eptr"
!do i = 0,nel
!write(un_chk,*)i,eptr(i) ;
!end do ;

!allocate (eind2(0:Neind-1) , NPart ( 0:NJG-1 ))
allocate (eind2(Neind) , NPart ( NJG ))

eind2 = eind(1:Neind) ;

!      Call METIS_PartMeshDual  ( NEl, NJG, eptr, eind(1:Neind), Vwgt, VSize, NCommonNodes, NParts,  tpwgts, MOptions, ObjVal, EPart, NPart(1:NJG) ) ; ! Work on options, see pages 20 and 28 of the manual Metis version 5.1.0.
!      Call METIS_PartMeshDual  ( NEl, NJG, eptr, eind(0:Neind-1), 0, 0, NCommonNodes, NParts,  0, 0, ObjVal, EPart, NPart (0:NJG-1) ) ; ! Work on options, see pages 20 and 28 of the manual Metis version 5.1.0.
      Call METIS_PartMeshDual  ( NEl, NJG, eptr, eind2, 0, 0, NCommonNodes, NParts,  0, 0, ObjVal, EPart, NPart ) ; ! Work on options, see pages 20 and 28 of the manual Metis version 5.1.0.
      !Call METIS_PartMeshDual  ( NEl, NJG, eptr, eind2, PVWgt, PVSize, NCommonNodes, NParts,  Ptpwgts, PMoptions, ObjVal, EPart, NPart ) ; ! Work on options, see pages 20 and 28 of the manual Metis version 5.1.0.

      write(*    ,*) "End subroutine < Metis Version 5.1.0 >" ;
      Write(UnInf,*) "End subroutine < Metis Version 5.1.0 >" ;
    End If ;

  Else If ( NParts == 1_Shrt ) then ;
    EPart = 1_Shrt ;
  End If ;

! Deallocating Arrays used in Metis or ParMetis
  ERR_DeAlloc = 0 ! arash--------------------------------------------------------------

  If ( MetisType == 0_Tiny ) Then ;
       if ( NParts /= 1 ) then ! arash--------------------------------------------------------------
          DeAllocate ( NPart,     STAT = ERR_DeAlloc ) ;
       end if
write(*,*) 'ERR_DeAlloc', ERR_DeAlloc
  Else If ( MetisType == 1_Tiny ) Then ;
    DeAllocate ( ElmDist, eptr, eind, ElmWgt, tpwgts, ubvec,      STAT = ERR_DeAlloc ) ;
  Else If ( MetisType == 2_Tiny ) Then ;
    DeAllocate ( VWgt, EPart, eptr, eind,       STAT = ERR_Alloc ) ;
  End If ;
  IF ( ERR_DeAlloc /= 0 ) Then ;
    Write (*, Fmt_DEALLCT) ERR_DeAlloc ;  Write (UnInf, Fmt_DEALLCT) ERR_DeAlloc ;
    Write(*, Fmt_FL) ;  Write(UnInf, Fmt_FL) ; Write(*, Fmt_End) ; Read(*,*) ;  STOP ;
  End If ;

! =============================================== ReNumbering =======================================================================================
Write (*,*)"Main:: Allocating required matrices for ReNumbering ..."

! Allocating rquired arrays for Numbering and output
Allocate ( NEL_Rank ( NParts ), NJ_Rank ( NParts ), Global_PETSc_Num ( NJ ), Local_PETSc_Num ( NJ, NParts ), ID_Application ( NJ, NDOF ), ID_PETSc ( NJ, NDOF ), NEqRank ( NParts ), NNodeRank ( NParts ),     STAT = ERR_Alloc ) ;
  IF ( ERR_Alloc /= 0 ) Then ;
    Write (*, Fmt_ALLCT) ERR_Alloc ;  Write (UnInf, Fmt_ALLCT) ERR_Alloc ;
    !#Call BEEP_FAIL ;
    Write(*, Fmt_FL) ;  Write(UnInf, Fmt_FL) ; Write(*, Fmt_End) ; Read(*,*) ;  STOP ;
  End If ;


! ----------------- added for inversion DS ------
    Allocate ( NPart (NJ) )
    NPart = 1_Shrt                               ! if we are running on one processor
! -----------------------------------------------


! Obtaining local node and element numbers for each rank. Renumbering equation numbers to get PETSc numbering.
NNeighbor = 6 ;

Call Numbering    (                                                                                                             &
MaxNNode, NDOF, NNeighbor, NParts, NEL, NJ, NEQMTotal,                                                                          & ! Integer Variables
!                                                                                                                               & ! Real Variables
NPart, EPart, INod, ID, NEL_Rank, NJ_Rank, Global_PETSc_Num, NEqRank, NNodeRank, Local_PETSc_Num, ID_Application, ID_PETSc,     & ! Integer Arrays
!                                                                                                                               & ! Real Arrays
!                                                                                                                               & ! Characters
Nodes                                                                                                                           & ! Type 
) ;

!do i = 1,nj
!  write(*,*) i, npart(i)
!end do
!write(*,*) Count ( npart ( : ) == 1 ) , Count ( npart ( : ) == 2 )
!stop

! =============================================== Calculating number of non-zero entries of PETSc objects ===========================================

Allocate ( D_NNZ_Stiff ( NEqMTotal ), O_NNZ_Stiff ( NEqMTotal ), D_NNZ_Damp ( NEqMTotal ), O_NNZ_Damp ( NEqMTotal ), D_NNZ_Mass ( NEqMTotal ), O_NNZ_Mass ( NEqMTotal ),     STAT = ERR_Alloc ) ;
  IF ( ERR_Alloc /= 0 ) Then ;
    Write (*, Fmt_ALLCT) ERR_Alloc ;  Write (UnInf, Fmt_ALLCT) ERR_Alloc ;
    !#Call BEEP_FAIL ;
    Write(*, Fmt_FL) ;  Write(UnInf, Fmt_FL) ; Write(*, Fmt_End) ; Read(*,*) ;  STOP ;
  End If ;

! LOAD GAUSS INTEGRATION DATA POINTS
!GAUSS_PNT = GAUSS_POINTS( NINT, NINT_Type ) ;

!Call Object_Size  (                                                                                                             &
!NDim, NLCase, Matrix_Type, NEqEl, El_Type, NNode, PARAM_Type, NINT, NDOF, NINT_Type,     NParts,     NJ, NEL,                   & ! Integer Variables
!PR, HREF,                                                                                                                       & ! Real Variables
!ELGR, LTEL, MTEL, ELT, NEqRank,     Global_PETSc_Num,     INOD, ID_PETSc,     D_NNZ_Stiff, O_NNZ_Stiff, D_NNZ_Damp, O_NNZ_Damp, D_NNZ_Mass, O_NNZ_Mass,  & 
!PMat, XYZ, PML_DIM                                                                                                              & ! Real Arrays
!!                                                                                                                               & ! Characters
!!GAUSS_PNT                                                                                                                      & ! Type 
!) ;

D_NNZ_Stiff = 100 ;
O_NNZ_Stiff = 100 ;
D_NNZ_Damp  = 100 ;
O_NNZ_Damp  = 100 ;
D_NNZ_Mass  = 100 ;
O_NNZ_Mass  = 100 ;

! =============================================== Heterogeneous Materials ===========================================================================

Allocate ( Node_Mat_ID ( NJ ), Node_Mat_Mapping ( NJ ),     STAT = ERR_Alloc ) ;
  IF ( ERR_Alloc /= 0 ) Then ;
    Write (*, Fmt_ALLCT) ERR_Alloc ;  Write (UnInf, Fmt_ALLCT) ERR_Alloc ;
    Write(*, Fmt_FL) ;  Write(UnInf, Fmt_FL) ; Write(*, Fmt_End) ; Read(*,*) ;  STOP ;
  End If ;

!write(*,*) 'hey!-----------------Het_Mat_Numbering is next'
Call Het_Mat_Numbering    (                                                                                                  &
NDim, MaxNNode, NDOF, NNeighbor, NParts, NEL, NJ, NEQMTotal,                                                                 & ! Integer Variables
!                                                                                                                            & ! Real Variables
MTEL, EPart, INod, ID, NEL_Rank, NJ_Rank, Global_PETSc_Num, NEqRank, NNodeRank, Local_PETSc_Num, ID_Application, ID_PETSc,   & ! Integer Arrays
Node_Mat_ID, Node_Mat_Mapping,                                                                                               &
PML_DIM, XYZ,                                                                                                                & ! Real Arrays
ModelName, OutDir                                                                                                            & ! Characters
!Nodes                                                                                                                       & ! Type 
)
!write(*,*) 'hey!-----------------Het_Mat_Numbering is passed'

!do i = 1,nj
!  write(*,*) i, Node_Mat_ID(i) , Node_Mat_Mapping(i)
!end do
!stop

! =============================================== Inversion Data Structure ==========================================================================
! 1- data-structure

Call Inversion_DS (                                                                                                          &
NDim, MaxNNode, NDOF, NParts, NEL, NJ, NEQMTotal,                                                                            & ! Integer Variables
!                                                                                                                            & ! Real Variables
EPart, INod, ID, NEL_Rank, NJ_Rank, Global_PETSc_Num, NEqRank, NNodeRank, Local_PETSc_Num, ID_Application, ID_PETSc,         & ! Integer Arrays
NPart, Node_Mat_ID, Node_Mat_Mapping,                                                                                        &
!                                                                                                                            & ! Real Arrays
ModelName, OutDir                                                                                                            & ! Characters
!                                                                                                                            & ! Type 
) ;

! 2- material visualization

Call OUTPUT_Mat_Vis_HDF5 (                                                                                                      &
NDim, MaxNNode, NDOF,                                                                                                           & ! Integer (1) Variables
!                                                                                                                               & ! Integer (2) Variables
NParts,                                                                                                                         & ! Integer (4) Variables
NEL, NJ,                                                                                                                        & ! Integer (8) Variables
!                                                                                                                               & ! Real Variables
EPart, ELT,                                                                                                                     &
NNodeRank,                                                                                                                      &
NEL_Rank, NJ_Rank, Global_PETSc_Num, Local_PETSc_Num, INod,                                                                     & ! Integer Arrays
!                                                                                                                               & ! Real Arrays
ModelName, OutDir                                                                                                               & ! Characters
!                                                                                                                               & ! Type
) ;

! =============================================== PRINT THE RESULTS OUT =============================================================================
  If (Output_Type == 0_Tiny) Then ;
    Call OUTPUT    (                                                                                                                  &
      NDim, MaxNNode, NDOF,                                                                                                           & ! Integer (1) Variables
      NGroup, NMat, NPM,                                                                                                              & ! Integer (2) Variables
      NParts,                                                                                                                         & ! Integer (4) Variables
      NEL, NJ, NEQMTotal,                                                                                                             & ! Integer (8) Variables
      !                                                                                                                               & ! Real Variables
      LoadC,     EPart, MTEL, ELT, ELGR, JLoad, IDBC,                                                                                 &
      D_NNZ_Stiff, O_NNZ_Stiff, D_NNZ_Damp, O_NNZ_Damp, D_NNZ_Mass, O_NNZ_Mass,     NDAN , NVAN, NAAN,  NEqRank, NNodeRank,           &
      NoBndry_DRM, NoLayer_DRM,  NEL_Rank, NJ_Rank, Global_PETSc_Num, Local_PETSc_Num, INod, ID, ID_Application, ID_PETSc,            & ! Integer Arrays
      PMat, PBLD, XYZ, UDis, PLoad, PML_DIM,                                                                                          & ! Real Arrays
      ModelName, OutDir,                                                                                                              & ! Characters
      Nodes, Param                                                                                                                    & ! Type
      ) ;
  Else If (Output_Type == 1_Tiny) Then ;
    Call OUTPUT_Binary (                                                                                                              &
      NDim, MaxNNode, NDOF,                                                                                                           & ! Integer (1) Variables
      NGroup, NMat, NPM,                                                                                                              & ! Integer (2) Variables
      NParts,                                                                                                                         & ! Integer (4) Variables
      NEL, NJ, NEQMTotal,                                                                                                             & ! Integer (8) Variables
      !                                                                                                                               & ! Real Variables
      LoadC,     EPart, MTEL, ELT, ELGR, JLoad, IDBC,                                                                                 &
      D_NNZ_Stiff, O_NNZ_Stiff, D_NNZ_Damp, O_NNZ_Damp, D_NNZ_Mass, O_NNZ_Mass,     NDAN , NVAN, NAAN,  NEqRank, NNodeRank,           &
      NoBndry_DRM, NoLayer_DRM,  NEL_Rank, NJ_Rank, Global_PETSc_Num, Local_PETSc_Num, INod, ID, ID_Application, ID_PETSc,            & ! Integer Arrays
      PMat, PBLD, XYZ, UDis, PLoad, PML_DIM,                                                                                          & ! Real Arrays
      ModelName, OutDir,                                                                                                              & ! Characters
      Nodes, Param                                                                                                                    & ! Type
      ) ;
  Else If (Output_Type == 2_Tiny) Then ;
    Call OUTPUT_HDF5 (                                                                                                                &
      NDim, MaxNNode, NDOF,                                                                                                           & ! Integer (1) Variables
      NGroup, NMat, NPM,                                                                                                              & ! Integer (2) Variables
      NParts,                                                                                                                         & ! Integer (4) Variables
      NEL, NJ, NEQMTotal,                                                                                                             & ! Integer (8) Variables
      !                                                                                                                               & ! Real Variables
      LoadC,     EPart, MTEL, ELT, ELGR, JLoad, IDBC,                                                                                 &
      D_NNZ_Stiff, O_NNZ_Stiff, D_NNZ_Damp, O_NNZ_Damp, D_NNZ_Mass, O_NNZ_Mass,     NDAN , NVAN, NAAN,  NEqRank, NNodeRank,           &
      NoBndry_DRM, NoLayer_DRM,  NEL_Rank, NJ_Rank, Global_PETSc_Num, Local_PETSc_Num, INod, ID, ID_Application, ID_PETSc,            & ! Integer Arrays
      PMat, PBLD, XYZ, UDis, PLoad, PML_DIM,                                                                                          & ! Real Arrays
      ModelName, OutDir,                                                                                                              & ! Characters
      Nodes, Param                                                                                                                    & ! Type
      ) ;
  End If ;

Deallocate ( MTEL, ELT, ELGR, JLoad, NDAN, NVAN, NAAN, INOD, IDBC, ID, PBLD, XYZ, UDis, PLoad, PML_DIM, EPart, NEL_Rank, NJ_Rank, Global_PETSc_Num, Local_PETSc_Num, ID_Application, ID_PETSc, NEqRank, D_NNZ_Stiff, O_NNZ_Stiff, D_NNZ_Damp, O_NNZ_Damp, D_NNZ_Mass, O_NNZ_Mass,     STAT = ERR_Alloc ) ;   ! , NNZ_MB, NNZ_P
  IF ( ERR_Alloc /= 0 ) Then ;
    Write (*, Fmt_ALLCT) ERR_Alloc ;  Write (UnInf, Fmt_ALLCT) ERR_Alloc ;
    !#Call BEEP_FAIL ;
    Write(*, Fmt_FL) ;  Write(UnInf, Fmt_FL) ; Write(*, Fmt_End) ; Read(*,*) ;  STOP ;
  End If ;

!  Select Case ( NAT ) ;
!
!! ---------------------------------------------------------------------------------------------------------------------------------------------------
!! - 64: LINEAR DYNAMIC ANALYSIS - TIME DOMAIN - DIRECT APPROACH - SOLID ELEMENTS - PML TRUNCATED BOUNDARIES - PETSC - FULL MATRICES - ASSEMBLE TOTAL ELEMENT MATRICES WITH ZEORS
!! - 67: LINEAR DYNAMICS ANALYSIS - TIME DOMAIN - DIRECT APPROACH - SOLID ELEMENTS - PML TRUNCATED BOUNDARIES - FULL MATRICES - Explicit Method ----------
!! ---------------------------------------------------------------------------------------------------------------------------------------------------
!    CASE ( ACN_LD_4 , ACN_LD_5, ACN_LD_8 ) ; ! 64 and 67
!      Include 'PIC_CasePML.F90'
!
!! ---------------------------------------------------------------------------------------------------------------------------------------------------
!! - LINEAR DYNAMICS ANALYSIS - TIME DOMAIN - DIRECT APPROACH - SOLID ELEMENTS - PML TRUNCATED BOUNDARIES - FULL MATRICES - DRM ----------------------
!! ---------------------------------------------------------------------------------------------------------------------------------------------------
!    CASE ( ACN_LI_DYN_TIME_DIR_SLD_23D_DRM, ACN_LD_7, ACN_LD_9 ) ; !65 and 69 
!      Include 'PIC_CaseDRM.F90' 
!
!    CASE DEFAULT ;
!      Write(*,*)" Type OF ANALYSIS IS NOT AVAILABLE IN THE MAIN SELECT CASE - CHECK THE INPUT FILE" ;
!      Write(UnInf,*)" Type OF ANALYSIS IS NOT AVAILABLE IN THE MAIN SELECT CASE - CHECK THE INPUT FILE" ;
!      Write(*,*) ;
!      Write(UnInf,*) ;
!      Write(*,*)' Program TERMINATED DUE TO SOME TECHNICAL PROBLEM - CHECK THE .INF FILE FOR FURHTER INFORMATION' ;
!      !#Call BEEP_FAIL ;
!      Write(*, Fmt_End) ; Read(*,*) ; ! STOP ;
!      Read(*," ('PRESS ENTER TO End ...') " ) ;
!      STOP
!
!  End Select ;

! =========================== RUNNING TIME OF THE CODE ==============================================================================================
Call CPU_TIME ( TimeE ) ;

Write(*     ,Fmt_RUNTIME) "TOTAL"   ,TimeE - TimeS    ;
Write(UnInf,*)"---------- RUNNING TIME STATISTICS ----------" ;
Write(UnInf,Fmt_RUNTIME) "Input   ", TimeInputE      - TimeInputS    ;


! =========================== Close FILES ===========================================================================================================
! - Closing ADDRESS FILE ------------------------------------------------------------------------------------------------------------------------------
UnFile =  UN_ADR ;
Close ( Unit = UnFile, STATUS = 'KEEP', ERR = 1002, IOSTAT = IO_File ) ;

! - Closing Input FILE --------------------------------------------------------------------------------------------------------------------------------
!UnFile =  UnInptMdl ;
!Close ( Unit = UnFile, STATUS = 'KEEP', ERR = 1002, IOSTAT = IO_File ) ;

! - Closing file for nodes' coordinates ---------------------------------------------------------------------------------------------------------------
UnFile = UnInptXYZ ;
Close ( Unit = UnFile, STATUS = 'KEEP', ERR = 1002, IOSTAT = IO_File ) ;

! - Input file for elements' connectivities ---------------------------------------------------------------------------------------------------------
UnFile = UnInptCnn ;
Close ( Unit = UnFile, STATUS = 'KEEP', ERR = 1002, IOSTAT = IO_File ) ;

! - Input file for nodes' constraints ---------------------------------------------------------------------------------------------------------------
UnFile = UnInptCnt ;
Close ( Unit = UnFile, STATUS = 'KEEP', ERR = 1002, IOSTAT = IO_File ) ;

! - Close IFORMATION FILE ---------------------------------------------------------------------------------------------------------------------------
UnFile =  UnInf ;
Close ( Unit = UnFile, STATUS = 'KEEP', ERR = 1002, IOSTAT = IO_File ) ;

! - Close CHECK FILE --------------------------------------------------------------------------------------------------------------------------------
UnFile =  UN_CHK ;
Close ( Unit = UnFile, STATUS = 'KEEP', ERR = 1002, IOSTAT = IO_File ) ;

! =========================== FINISH THE CODE =======================================================================================================
!#Call BEEP_SUCC
Write (*, Fmt_SUC) ;  Write(UnInf, Fmt_SUC) ;
Write(*, Fmt_End) ;


! - SHUT DOWN PETSC ---------------------------------------------------------------------------------------------------------------------------------
Call PetscFinalize ( ErrPTC ) ;

!#Read(*,*);
STOP ;

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

End Program PIC ; 
