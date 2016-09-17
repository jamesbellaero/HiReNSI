
!****************************************************************************************************************************************************
! This component reads parallel input files created by Babak's PIC program
! June 20 2013
! November 8 2013 HDF5
!****************************************************************************************************************************************************

! =============================================== Opening external input files ======================================================================
! - Input FILE --------------------------------------------------------------------------------------------------------------------------------------
UnFile = UN_ADR ;
Open ( Unit = UnFile, FILE = 'ADDRESS_PTC.TXT', ERR =  1001, IOSTAT = IO_File, ACCESS = 'SEQUENTIAL', ACTION = 'READ', ASYNCHRONOUS = 'NO', BLANK = 'NULL', BLOCKSize = 0, DisPOSE = 'KEEP', FORM = 'FORMATTED', POSITION = 'ASIS', STATUS = 'OLD' ) ;

! - READ THE Input FILE Name AND DIRECTORIES FROM "ADDRESS FILE" IN THE CURRENT DIRECTORY -----------------------------------------------------------
Read(UN_ADR,*) ;
Read(UN_ADR,*) ModelName ;
Read(UN_ADR,*) ;
Read(UN_ADR,*) NAT, Output_Type ;
Read(UN_ADR,*) ;
Read(UN_ADR,*) Model_InDir ;
Read(UN_ADR,*) ;
Read(UN_ADR,*) InlDir ;
Read(UN_ADR,*) ;
Read(UN_ADR,*) OutDir ;
Read(UN_ADR,*) ;
Read(UN_ADR,*) AnaName ;

Ana_InDir   = TRIM(AdjustL (Model_InDir))//'/'//TRIM(AdjustL (ModelName))//'/'//'Analysis' ;
Model_InDir = TRIM(AdjustL (Model_InDir))//'/'//TRIM(AdjustL (ModelName))//'/'//'Model' ;

!write(*,*) Ana_InDir
!write(*,*) Model_InDir

! - Input FILE --------------------------------------------------------------------------------------------------------------------------------------

Write (IndexRank, *) Rank ; ! Converts Rank number to Character foramt for the file Name
Write (IndexSize, *) Size ; ! Converts Size number to Character foramt for the file Name

! Format of the input files ; Binary or Formatted
If      ( Output_Type == 0 ) Then ; Format_Type = "Formatted" ;
Else If ( Output_Type == 1 ) Then ; Format_Type = "Binary" ;
End If ;

!write(*,*) '---------- purplepoigirl ----------1'
! Model Data file 
UnFile = UnInptMdl ;
Open ( Unit = UnFile, FILE = TRIM(ModelName)//'_'//Trim(AdjustL(IndexSize))//'_'//Trim(AdjustL(IndexRank))//'.dataModel', ERR =  1001, IOSTAT = IO_File, ACCESS = 'SEQUENTIAL', ACTION = 'READ', ASYNCHRONOUS = 'NO', BLANK = 'NULL', BLOCKSize = 0, DEFAULTFILE = TRIM(Model_InDir), DisPOSE = 'KEEP', FORM = 'Formatted', POSITION = 'ASIS', STATUS = 'Old' ) ;! 


! - Create the result folder ------------------------------------------------------------------------------------------------------------------------
! results folder
Directory = MakeDirQQ (TRIM(AdjustL (OutDir))//'/'//TRIM(AdjustL (ModelName))  ) ;
IF (Directory) THEN ;
   IF ( I_PRINT_SCREEN == 1 )   WRITE (*     ,*) 'New subdirectory successfully created' ;
   WRITE (UnInf,*) 'New subdirectory successfully created' ;
ELSE ;
   IF ( I_PRINT_SCREEN == 1 )   WRITE (*     ,*) 'Failed to create subdirectory' ;
   WRITE (UnInf,*) 'Failed to create subdirectory' ;
END IF ;

Directory = MakeDirQQ (TRIM(AdjustL (OutDir))//'/'//TRIM(AdjustL (ModelName))//'/'//TRIM(AdjustL (AnaName))  ) ;
IF (Directory) THEN ;
   IF ( I_PRINT_SCREEN == 1 )   WRITE (*     ,*) 'New subdirectory successfully created' ;
   WRITE (UnInf,*) 'New subdirectory successfully created' ;
ELSE ;
      IF ( I_PRINT_SCREEN == 1 ) WRITE (*     ,*) 'Failed to create subdirectory' ;
      WRITE (UnInf,*) 'Failed to create subdirectory' ;
END IF ;

! internal folder
Directory = MakeDirQQ (TRIM(AdjustL (InlDir))//'/'//TRIM(AdjustL (ModelName))) ;
  IF (Directory) THEN ;
   IF ( I_PRINT_SCREEN == 1 )     WRITE (*    ,*) 'New subdirectory successfully created' ;
     WRITE (UnInf,*) 'New subdirectory successfully created' ;
  ELSE ;
   IF ( I_PRINT_SCREEN == 1 )     WRITE (*    ,*) 'Failed to create subdirectory' ;
     WRITE (UnInf,*) 'Failed to create subdirectory' ;
  END IF ;

Directory = MakeDirQQ (TRIM(AdjustL (InlDir))//'/'//TRIM(AdjustL (ModelName))//'/'//TRIM(AdjustL (AnaName))  ) ;
  IF (Directory) THEN ;
   IF ( I_PRINT_SCREEN == 1 )     WRITE (*    ,*) 'New subdirectory successfully created' ;
     WRITE (UnInf,*) 'New subdirectory successfully created' ;
  ELSE ;
   IF ( I_PRINT_SCREEN == 1 )     WRITE (*    ,*) 'Failed to create subdirectory' ;
     WRITE (UnInf,*) 'Failed to create subdirectory' ;
  END IF ;

OutDir = TRIM(AdjustL (OutDir))//'/'//TRIM(AdjustL (ModelName))//'/'//TRIM(AdjustL (AnaName)) ;
InlDir = TRIM(AdjustL (InlDir))//'/'//TRIM(AdjustL (ModelName))//'/'//TRIM(AdjustL (AnaName)) ;

! - INFORMATON FILE ---------------------------------------------------------------------------------------------------------------------------------
UnFile = UnInf ;
Open ( Unit = UnFile, FILE = TRIM(AnaName)//'_'//Trim(AdjustL(IndexSize))//'_'//Trim(AdjustL(IndexRank))//'.inf' , ERR =  1001, IOSTAT = IO_File, ACCESS = 'SEQUENTIAL', ACTION = 'Write', ASYNCHRONOUS = 'NO', BLANK = 'NULL', BLOCKSize = 0, DEFAULTFILE = TRIM(OutDir), DisPOSE = 'KEEP', FORM = 'FORMATTED', POSITION = 'ASIS', STATUS = 'REPLACE' ) ;



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

! =============================================== Reading model data ================================================================================

! - Reading Input FILE ------------------------------------------------------------------------------------------------------------------------
    Call Input (                                                                                                                    &
    NDOF, MaxNNode, NDim,                                                                                                           & ! Integer (1) Variables
    NGroup, NPM, NMat,                                                                                                              & ! Integer (2) Variables
    NEL, NJ, NJTotal, NEQM, NEQMTotal, NEQM_Mapping,                                                                                & ! Integer (4) Variables
    !                                                                                                                               & ! Integer (8) Variables
    LoadC,                                                                                                                          & ! Real Variables
    !                                                                                                                               & ! Integer Arrays
    !                                                                                                                               & ! Real Arrays
    !                                                                                                                               & ! Characters
    Param                                                                                                                           & ! Type
    ) ;


! Allocating required arrays
Allocate ( MTEL (NEl), ELT (NEl), ELGR (NEl), INOD (MaxNNode,NEl), PML_DIM (2*NDim, 2), XYZ (NJ,NDim), ID (NJ,NDOF), PMat (NMat,NPM),  &
           IDBC  ( Param%IntM( 2, 4), 2*NDim + 1), &  ! Param%IntP( 3,1) = NIDBC
           JLoad ( Param%IntM( 3, 1) ), PLoad ( NDim, Param%IntM( 3,1) ), & ! Param%IntP(LoadC (3),1) = NLN
           NDAN  ( Param%IntM( 6, 1) ), NVAN ( Param%IntM( 6, 2) ), NAAN ( Param%IntM( 6, 3) ), NDPN ( Param%IntM( 6, 1) ), NVPN ( Param%IntM( 6, 2) ), NAPN ( Param%IntM( 6, 3) ), EqDis (Param%IntM( 6, 1)*NDim), EqVel (Param%IntM( 6, 2)*NDim), EqAcc (Param%IntM( 6, 3)*NDim), & ! ( Param%IntP( NLCase + 2, I), I = 1, 3 ) ; ! NNDH, NNVH, NNAH
           PBLD  ( Param%IntM( 1, 1), Param%IntM( 1, 2) ), & ! Param%IntP( 1, 1), Param%IntP( 1, 2) ; ! NBLD   = NGroup, NPBL   = NDim ;
           UDis (NDOF+1, Param%IntM( 4, 1) ), & ! Param%IntP(LoadC (4),1) = NSDN
           STAT = ERR_Alloc) ;
  IF ( ERR_Alloc /= 0 ) Then ;
    Write (*, Fmt_ALLCT) ERR_Alloc ;  Write (UnInf, Fmt_ALLCT) ERR_Alloc ;
    !#Call BEEP_FAIL ;
    Write(*, Fmt_FL) ;  Write(UnInf, Fmt_FL) ; Write(*, Fmt_End) ; Read(*,*) ;  STOP ;
  End If ;

! Allocating required arrays for PETSc
Allocate ( App_Numbers ( 0:NEQM_Mapping - 1 ), PETSc_Numbers ( 0:NEQM_Mapping - 1 ), Indices ( 0:NEQM - 1 ), D_NNZ_Stiff ( NEQM_Mapping ), O_NNZ_Stiff ( NEQM_Mapping ), D_NNZ_Damp ( NEQM_Mapping ), O_NNZ_Damp ( NEQM_Mapping ), D_NNZ_Mass ( NEQM_Mapping ), O_NNZ_Mass ( NEQM_Mapping ),      STAT = ERR_Alloc ) ;
  IF ( ERR_Alloc /= 0 ) Then ;
    Write (*, Fmt_ALLCT) ERR_Alloc ;  Write (UnInf, Fmt_ALLCT) ERR_Alloc ;
    !#Call BEEP_FAIL ;
    Write(*, Fmt_FL) ;  Write(UnInf, Fmt_FL) ; Write(*, Fmt_End) ; Read(*,*) ;  STOP ;
  End If ;

  ! ALLOCATING REQUIRED Arrays for DRM analysis
  If ( Param%IntM( 5, 3 ) == 3 ) Then ; 
    Allocate ( NoBndry_DRM( Param%IntM( 7, 1) ), NoLayer_DRM ( Param%IntM( 7, 2) ), InciWave (5),  ND_b ( 0:NDim * Param%IntM( 7, 1) - 1 ), ND_e ( 0:NDim * Param%IntM( 7, 2) - 1 ),     STAT = ERR_Alloc) ;
      IF ( ERR_Alloc /= 0 ) Then ;
        Write (*, Fmt_ALLCT) ERR_Alloc ;  Write (UnInf, Fmt_ALLCT) ERR_Alloc ;
        !#Call BEEP_FAIL ;
        Write(*, Fmt_FL) ;  Write(UnInf, Fmt_FL) ; Write(*, Fmt_End) ; Read(*,*) ;  STOP ;
      End If ;
  End If ;

Include 'PTC_Open_Inc.F90'  ! Include file for opening files

! Reading input arrays
!    Call Input (                                                                                                                    &
!    NDOF, NDim, MaxNNode,                                                                                                           & ! Integer (1) Variables
!    NGroup, NPM, NMat,     NLCase,                                                                                                  & ! Integer (2) Variables
!    !                                                                                                                               & ! Integer (4) Variables
!    NEL, NJ, NEQM, NEQM_Mapping,                                                                                                    & ! Integer (8) Variables
!    !                                                                                                                               & ! Real Variables
!    LoadC,     IDBC,      MTEL, ELT,     ELGR,                                                                                      &
!    D_NNZ_Stiff, O_NNZ_Stiff, D_NNZ_Damp, O_NNZ_Damp, D_NNZ_Mass, O_NNZ_Mass,     App_Numbers, PETSc_Numbers,     Indices,          &
!    NDAN , NVAN, NAAN, NDPN , NVPN, NAPN,     NoBndry_DRM, NoLayer_DRM, JLoad,     EqDis, EqVel, EqAcc,      INOD, ID,              & ! Integer Arrays
!    PMat, PBLD, XYZ, UDis, PLoad, PML_DIM,                                                                                          & ! Real Arrays
!    !                                                                                                                               & ! Characters
!    Param                                                                                                                           & ! Type
!    ) ;

  If ( Output_Type == 0 ) Then ;
    Call Input (                                                                                                                    &
    NDOF, NDim, MaxNNode,                                                                                                           & ! Integer (1) Variables
    NGroup, NPM, NMat,     NLCase,                                                                                                  & ! Integer (2) Variables
    !                                                                                                                               & ! Integer (4) Variables
    NEL, NJ, NEQM, NEQM_Mapping,                                                                                                    & ! Integer (8) Variables
    !                                                                                                                               & ! Real Variables
    LoadC,     IDBC,      MTEL, ELT,     ELGR,                                                                                      &
    D_NNZ_Stiff, O_NNZ_Stiff, D_NNZ_Damp, O_NNZ_Damp, D_NNZ_Mass, O_NNZ_Mass,     App_Numbers, PETSc_Numbers,     Indices,          &
    NDAN , NVAN, NAAN, NDPN , NVPN, NAPN,     NoBndry_DRM, NoLayer_DRM, JLoad,     EqDis, EqVel, EqAcc,      INOD, ID,              & ! Integer Arrays
    PMat, PBLD, XYZ, UDis, PLoad, PML_DIM,                                                                                          & ! Real Arrays
    !                                                                                                                               & ! Characters
    Param                                                                                                                           & ! Type
    ) ;
!  Else If ( Output_Type == 1 ) Then ;
!    Call Input_Binary (                                                                                                             &
!    NDOF, NDim, MaxNNode,                                                                                                           & ! Integer (1) Variables
!    NGroup, NPM, NMat,     NLCase,                                                                                                  & ! Integer (2) Variables
!    !                                                                                                                               & ! Integer (4) Variables
!    NEL, NJ, NEQM, NEQM_Mapping,                                                                                                    & ! Integer (8) Variables
!    !                                                                                                                               & ! Real Variables
!    LoadC,     IDBC,      MTEL, ELT,     ELGR,                                                                                      &
!    D_NNZ_Stiff, O_NNZ_Stiff, D_NNZ_Damp, O_NNZ_Damp, D_NNZ_Mass, O_NNZ_Mass,     App_Numbers, PETSc_Numbers,     Indices,          &
!    NDAN , NVAN, NAAN, NDPN , NVPN, NAPN,     NoBndry_DRM, NoLayer_DRM, JLoad,     EqDis, EqVel, EqAcc,      INOD, ID,              & ! Integer Arrays
!    PMat, PBLD, XYZ, UDis, PLoad, PML_DIM,                                                                                          & ! Real Arrays
!    !                                                                                                                               & ! Characters
!    Param                                                                                                                           & ! Type
!    ) ;
  Else If ( Output_Type == 2 ) Then ;
    Call Input_HDF5 (                                                                                                               &
    NDOF, NDim, MaxNNode,                                                                                                           & ! Integer (1) Variables
    NGroup, NPM, NMat,     NLCase,                                                                                                  & ! Integer (2) Variables
    !                                                                                                                               & ! Integer (4) Variables
    NEL, NJ, NEQM, NEQM_Mapping,                                                                                                    & ! Integer (8) Variables
    !                                                                                                                               & ! Real Variables
    LoadC,     IDBC,      MTEL, ELT,     ELGR,                                                                                      &
!    D_NNZ_Stiff, O_NNZ_Stiff, D_NNZ_Damp, O_NNZ_Damp, D_NNZ_Mass, O_NNZ_Mass,     App_Numbers, PETSc_Numbers,     Indices,          &
    App_Numbers, PETSc_Numbers,     Indices,                                                                                        &
    NDAN , NVAN, NAAN, NDPN , NVPN, NAPN,     NoBndry_DRM, NoLayer_DRM, JLoad,     EqDis, EqVel, EqAcc,      INOD, ID,              & ! Integer Arrays
    PMat, PBLD, XYZ, UDis, PLoad, PML_DIM,                                                                                          & ! Real Arrays
    ModelName, IndexRank, IndexSize, Model_InDir,                                                                                   & ! Characters
    Param                                                                                                                           & ! Type
    ) ;
  End If ;

! ALLOCATING REQUIRED Arrays for Paraview
Allocate ( ID_Para ( Param%IntM ( 5, 2), NDim ), idx_from ( Param%IntM ( 5, 1) ),     STAT = ERR_Alloc) ;
  IF ( ERR_Alloc /= 0 ) Then ;
    Write (*, Fmt_ALLCT) ERR_Alloc ;  Write (UnInf, Fmt_ALLCT) ERR_Alloc ;
    !#Call BEEP_FAIL ;
    Write(*, Fmt_FL) ;  Write(UnInf, Fmt_FL) ; Write(*, Fmt_End) ; Read(*,*) ;  STOP ;
  End If ;

! =============================================== G matrix in 3D ====================================================================================
! Allocating required arrays for PETSc
Allocate ( O_NNZ_G ( NEQM_Mapping ), D_NNZ_G ( NEQM_Mapping ), STAT = ERR_Alloc ) ;
  IF ( ERR_Alloc /= 0 ) Then ;
    Write (*, Fmt_ALLCT) ERR_Alloc ;  Write (UnInf, Fmt_ALLCT) ERR_Alloc ;
    !#Call BEEP_FAIL ;
    Write(*, Fmt_FL) ;  Write(UnInf, Fmt_FL) ; Write(*, Fmt_End) ; Read(*,*) ;  STOP ;
  End If ;


! ===============================================  idx_from  ========================================================================================
      ! - Input data for analysis -------------------------------------------------------------------------------------------------------------------

        If ( Output_Type == 2 ) Then ;
          Call Input_HDF5 (                                                                                                               &
          NDim, Solver_Type, Int_Order,                                                                                                   & ! Integer (1) Variables
          !NPM, NMat,                                                                                                                     & ! Integer (2) Variables
          !                                                                                                                               & ! Integer (4) Variables
          NJ,                                                                                                                             & ! Integer (8) Variables
          !                                                                                                                               & ! Real Variables
          LoadC,    Step,    idx_from,     ID_Para,                                                                                       & ! Integer Arrays
          InciWave,        BAcl,                                                                                                          & ! Real Arrays
          ModelName, IndexRank, IndexSize, Model_InDir,                                                                                   & ! Characters
          Param                                                                                                                           & ! Type
          ) ;
        End If ;

! =============================================== Reading heterogeneous material properties =========================================================

! Lambda (.Lambda)
UnFile = UnInpt_Lambda ;
!Open ( Unit = UnFile, FILE = TRIM(ModelName)//'_'//Trim(AdjustL(IndexSize))//'_'//Trim(AdjustL(IndexRank))//'.Lambda',  ERR =  1001, IOSTAT = IO_File, ACCESS = 'SEQUENTIAL', ACTION ='Read', ASYNCHRONOUS = 'NO', BLANK = 'NULL', BLOCKSIZE = 0, DEFAULTFILE = TRIM(Model_InDir), DisPOSE = 'KEEP', FORM = Trim(AdjustL(Format_Type)), POSITION = 'ASIS', STATUS ='Old') ;
Open ( Unit = UnFile, FILE = TRIM(ModelName)//'_'//Trim(AdjustL(IndexSize))//'_'//Trim(AdjustL(IndexRank))//'.Lambda',  ERR =  1001, IOSTAT = IO_File, ACCESS = 'SEQUENTIAL', ACTION ='Read', ASYNCHRONOUS = 'NO', BLANK = 'NULL', BLOCKSIZE = 0, DEFAULTFILE = TRIM(Model_InDir), DisPOSE = 'KEEP', FORM = 'Formatted', POSITION = 'ASIS', STATUS ='Old') ;

! Mu (.Mu)
UnFile = UnInpt_Mu ;
Open ( Unit = UnFile, FILE = TRIM(ModelName)//'_'//Trim(AdjustL(IndexSize))//'_'//Trim(AdjustL(IndexRank))//'.Mu',  ERR =  1001, IOSTAT = IO_File, ACCESS = 'SEQUENTIAL', ACTION ='Read', ASYNCHRONOUS = 'NO', BLANK = 'NULL', BLOCKSIZE = 0, DEFAULTFILE = TRIM(Model_InDir), DisPOSE = 'KEEP', FORM = 'Formatted', POSITION = 'ASIS', STATUS ='Old') ;

! allocate arrays
Allocate ( PMat_Lambda (NJ) , PMat_Mu (NJ) )
Allocate ( PMat_Temp (NJ) )

    Call Input (                                                                                                                    &
    !                                                                                                                               & ! Integer (1) Variables
    !                                                                                                                               & ! Integer (2) Variables
    !                                                                                                                               & ! Integer (4) Variables
    NJ,                                                                                                                             & ! Integer (8) Variables
    !                                                                                                                               & ! Real Variables
    !                                                                                                                               &
    !                                                                                                                               &
    !                                                                                                                               & ! Integer Arrays
    PMat_Lambda, PMat_Mu                                                                                                            & ! Real Arrays
    !                                                                                                                               & ! Characters
    !                                                                                                                               & ! Type
    ) ;

! =========================== Material Visualization index-set from =================================================================================

     Call Input_Mat_Vis_Basic (                                                                                                     &
     !                                                                                                                              & ! Integer (1) Variables
     !                                                                                                                              & ! Integer (2) Variables
     !                                                                                                                              & ! Integer (4) Variables
     NJ_Para,                                                                                                                       & ! Integer (8) Variables
     !                                                                                                                              & ! Real Variables
     !                                                                                                                              &
     !                                                                                                                              &
     !                                                                                                                              & ! Integer Arrays
     !                                                                                                                              & ! Real Arrays
     ModelName, IndexRank, IndexSize, Model_InDir                                                                                   & ! Characters
     !                                                                                                                              & ! Type
     ) ;


! allocate arrays
Allocate ( idx_Mat_Vis_from (NJ_Para) )


     Call Input_Mat_Vis_HDF5 (                                                                                                      &
     !                                                                                                                              & ! Integer (1) Variables
     !                                                                                                                              & ! Integer (2) Variables
     !                                                                                                                              & ! Integer (4) Variables
     NJ_Para,                                                                                                                       & ! Integer (8) Variables
     !                                                                                                                              & ! Real Variables
     idx_Mat_Vis_from,                                                                                                              & ! Integer Arrays
     !                                                                                                                              & ! Real Arrays
     ModelName, IndexRank, IndexSize, Model_InDir                                                                                   & ! Characters
     !                                                                                                                              & ! Type
     ) ;

! =========================== Close FILES ===========================================================================================================
Include 'PTC_Close_Inc.F90' 


! =============================================== Opening external output files ======================================================================
   IF ( I_PRINT_SCREEN == 1 ) Write (*    ,"('Opening Output files')") ;
Write (UnInf,"('Opening Output files')") ;

Write (IndexRank, *) Rank ;  ! Converts Rank number to Character foramt for the file Name
Write (IndexSize, *) Size ;  ! Converts Size number to Character foramt for the file Name

! History of Displacements and stresses 
UnFile = UN_HisD ;
Open ( Unit = UnFile, FILE = TRIM(ModelNAME)//'_'//Trim(AdjustL(IndexSize))//'_'//Trim(AdjustL(IndexRank))//'.HisD.txt', ERR =  1001, IOSTAT = IO_File, ACCESS = 'SEQUENTIAL', ACTION ='Write', ASYNCHRONOUS = 'NO', BLANK = 'NULL', BLOCKSIZE = 0, DEFAULTFILE = TRIM(OutDir), DisPOSE = 'KEEP', FORM = 'FORMATTED', POSITION = 'ASIS', STATUS ='REPLACE') ;

! History of velocity 
UnFile = UN_HisV ;
Open ( Unit = UnFile, FILE = TRIM(ModelName)//'_'//Trim(AdjustL(IndexSize))//'_'//Trim(AdjustL(IndexRank))//'.HisV', ERR =  1001, IOSTAT = IO_File, ACCESS = 'SEQUENTIAL', ACTION ='Write', ASYNCHRONOUS = 'NO', BLANK = 'NULL', BLOCKSIZE = 0, DEFAULTFILE = TRIM(OutDir), DisPOSE = 'KEEP', FORM = 'FORMATTED', POSITION = 'ASIS', STATUS ='REPLACE') ;

! History of acceleration 
UnFile = UN_HisA ;
Open ( Unit = UnFile, FILE = TRIM(ModelName)//'_'//Trim(AdjustL(IndexSize))//'_'//Trim(AdjustL(IndexRank))//'.HisA', ERR =  1001, IOSTAT = IO_File, ACCESS = 'SEQUENTIAL', ACTION ='Write', ASYNCHRONOUS = 'NO', BLANK = 'NULL', BLOCKSIZE = 0, DEFAULTFILE = TRIM(OutDir), DisPOSE = 'KEEP', FORM = 'FORMATTED', POSITION = 'ASIS', STATUS ='REPLACE') ;

! History of energy
UnFile = UN_Engy ;
Open ( Unit = UnFile, FILE = TRIM(ModelName)//'_'//Trim(AdjustL(IndexSize))//'_'//Trim(AdjustL(IndexRank))//'.Enr.txt',  ERR =  1001, IOSTAT = IO_File, ACCESS = 'SEQUENTIAL', ACTION ='Write', ASYNCHRONOUS = 'NO', BLANK = 'NULL', BLOCKSIZE = 0, DEFAULTFILE = TRIM(OutDir), DisPOSE = 'KEEP', FORM = 'FORMATTED', POSITION = 'ASIS', STATUS ='REPLACE') ;

! Wraper file for Paraview
  If ( Rank == 0 ) Then ;  ! Top file wraper 
    UnFile = UN_OutallWrp ;
    Open ( Unit = UnFile, FILE = TRIM(AnaName)//'_'//'all'//'_'//Trim(AdjustL(IndexSize))//'.xmf',  ERR =  1001, IOSTAT = IO_File, ACCESS = 'SEQUENTIAL', ACTION ='Write', ASYNCHRONOUS = 'NO', BLANK = 'NULL', BLOCKSIZE = 0, DEFAULTFILE = TRIM(OutDir), DisPOSE = 'KEEP', FORM = 'FORMATTED', POSITION = 'ASIS', STATUS ='REPLACE') ;
  End If ;

NNDH = Param%IntM( 6, 1)
NNVH = Param%IntM( 6, 2)
NNAH = Param%IntM( 6, 3)

! - results time-history header file ----------------------------------------------------------------------------------------------------------------
    Call ResultsHeader (                                                                                                            &
    NNDH, NNVH, NNAH,                                                                                                               & ! Integer Variables
    !                                                                                                                               & ! Real Variables
    NDAN, NDPN, NVAN, NVPN, NAAN, NAPN                                                                                              & ! Integer Arrays
    !                                                                                                                               & ! Real Arrays
    !                                                                                                                               & ! Characters
    !                                                                                                                               & ! Type 
    ) ;

