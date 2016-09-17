
!****************************************************************************************************************************************************
! This component reads parallel input files created by Babak's PIC program
! June 20 2013
!****************************************************************************************************************************************************

! =============================================== Opening external input files ======================================================================
! - Input FILE --------------------------------------------------------------------------------------------------------------------------------------
UnFile = UN_ADR ;
Open ( Unit = UnFile, FILE = 'ADDRESS_PTC.TXT', ERR =  1001, IOSTAT = IO_File, ACCESS = 'SEQUENTIAL', ACTION = 'READ', ASYNCHRONOUS = 'NO', BLANK = 'NULL', BLOCKSize = 0                             , DisPOSE = 'KEEP', FORM = 'FORMATTED', POSITION = 'ASIS', STATUS = 'OLD' ) ;
 
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

! Model Data file 
UnFile = UnInptMdl ;
Open ( Unit = UnFile, FILE = TRIM(ModelName)//'_'//Trim(AdjustL(IndexSize))//'_'//Trim(AdjustL(IndexRank))//'.dataModel', ERR =  1001, IOSTAT = IO_File, ACCESS = 'SEQUENTIAL', ACTION = 'READ', ASYNCHRONOUS = 'NO', BLANK = 'NULL', BLOCKSize = 0, DEFAULTFILE = TRIM(Model_InDir), DisPOSE = 'KEEP', FORM = Trim(AdjustL(Format_Type)), POSITION = 'ASIS', STATUS = 'Old' ) ;! 

! Coordinate file (.XYZ)
UnFile = UnInptXYZ ;
Open ( Unit = UnFile, FILE = TRIM(ModelName)//'_'//Trim(AdjustL(IndexSize))//'_'//Trim(AdjustL(IndexRank))//'.XYZ',  ERR =  1001, IOSTAT = IO_File, ACCESS = 'SEQUENTIAL', ACTION ='Read', ASYNCHRONOUS = 'NO', BLANK = 'NULL', BLOCKSIZE = 0, DEFAULTFILE = TRIM(Model_InDir), DisPOSE = 'KEEP', FORM = Trim(AdjustL(Format_Type)), POSITION = 'ASIS', STATUS ='Old') ;

! Connectivity file (.Cnn)
UnFile = UnInptCnn ;
Open ( Unit = UnFile, FILE = TRIM(ModelName)//'_'//Trim(AdjustL(IndexSize))//'_'//Trim(AdjustL(IndexRank))//'.Cnn', ERR =  1001, IOSTAT = IO_File, ACCESS = 'SEQUENTIAL', ACTION ='Read', ASYNCHRONOUS = 'NO', BLANK = 'NULL', BLOCKSIZE = 0, DEFAULTFILE = TRIM(Model_InDir), DisPOSE = 'KEEP', FORM = Trim(AdjustL(Format_Type)), POSITION = 'ASIS', STATUS ='Old') ;

! Constraint file (.Cnt)
UnFile = UnInptCnt ;
Open ( Unit = UnFile, FILE = TRIM(ModelName)//'_'//Trim(AdjustL(IndexSize))//'_'//Trim(AdjustL(IndexRank))//'.Cnt', ERR =  1001, IOSTAT = IO_File, ACCESS = 'SEQUENTIAL', ACTION ='Read', ASYNCHRONOUS = 'NO', BLANK = 'NULL', BLOCKSIZE = 0, DEFAULTFILE = TRIM(Model_InDir), DisPOSE = 'KEEP', FORM = Trim(AdjustL(Format_Type)), POSITION = 'ASIS', STATUS ='Old') ;

! Material Properties (.Mat)
UnFile = UnInptMat ;
Open ( Unit = UnFile, FILE = TRIM(ModelName)//'_'//Trim(AdjustL(IndexSize))//'_'//Trim(AdjustL(IndexRank))//'.Mat', ERR =  1001, IOSTAT = IO_File, ACCESS = 'SEQUENTIAL', ACTION ='Read', ASYNCHRONOUS = 'NO', BLANK = 'NULL', BLOCKSIZE = 0, DEFAULTFILE = TRIM(Model_InDir), DisPOSE = 'KEEP', FORM = Trim(AdjustL(Format_Type)), POSITION = 'ASIS', STATUS ='Old') ;

! Non-Zeros of Stiffness file (.NNZStiff)
UnFile = UnInptStiff ;
Open ( Unit = UnFile, FILE = TRIM(ModelName)//'_'//Trim(AdjustL(IndexSize))//'_'//Trim(AdjustL(IndexRank))//'.NNZStiff', ERR =  1001, IOSTAT = IO_File, ACCESS = 'SEQUENTIAL', ACTION ='Read', ASYNCHRONOUS = 'NO', BLANK = 'NULL', BLOCKSIZE = 0, DEFAULTFILE = TRIM(Model_InDir), DisPOSE = 'KEEP', FORM = Trim(AdjustL(Format_Type)), POSITION = 'ASIS', STATUS ='Old') ;

! Non-Zeros of Damping file (.NNZDamp)
UnFile = UnInptDamp ;
Open ( Unit = UnFile, FILE = TRIM(ModelName)//'_'//Trim(AdjustL(IndexSize))//'_'//Trim(AdjustL(IndexRank))//'.NNZDamp', ERR =  1001, IOSTAT = IO_File, ACCESS = 'SEQUENTIAL', ACTION ='Read', ASYNCHRONOUS = 'NO', BLANK = 'NULL', BLOCKSIZE = 0, DEFAULTFILE = TRIM(Model_InDir), DisPOSE = 'KEEP', FORM = Trim(AdjustL(Format_Type)), POSITION = 'ASIS', STATUS ='Old') ;

! Non-Zeros of Damping file (.NNZMass)
UnFile = UnInptMass ;
Open ( Unit = UnFile, FILE = TRIM(ModelName)//'_'//Trim(AdjustL(IndexSize))//'_'//Trim(AdjustL(IndexRank))//'.NNZMass', ERR =  1001, IOSTAT = IO_File, ACCESS = 'SEQUENTIAL', ACTION ='Read', ASYNCHRONOUS = 'NO', BLANK = 'NULL', BLOCKSIZE = 0, DEFAULTFILE = TRIM(Model_InDir), DisPOSE = 'KEEP', FORM = Trim(AdjustL(Format_Type)), POSITION = 'ASIS', STATUS ='Old') ;

! Application Ordering file (.App)
UnFile = UnInptApp ;
Open ( Unit = UnFile, FILE = TRIM(ModelName)//'_'//Trim(AdjustL(IndexSize))//'_'//Trim(AdjustL(IndexRank))//'.App', ERR =  1001, IOSTAT = IO_File, ACCESS = 'SEQUENTIAL', ACTION ='Read', ASYNCHRONOUS = 'NO', BLANK = 'NULL', BLOCKSIZE = 0, DEFAULTFILE = TRIM(Model_InDir), DisPOSE = 'KEEP', FORM = Trim(AdjustL(Format_Type)), POSITION = 'ASIS', STATUS ='Old') ;

! - Create the result folder ------------------------------------------------------------------------------------------------------------------------
! results folder
Directory = MakeDirQQ (TRIM(AdjustL (OutDir))//'/'//TRIM(AdjustL (ModelName))  ) ;
IF (Directory) THEN ;
   WRITE (*     ,*) 'New subdirectory successfully created' ;
   WRITE (UnInf,*) 'New subdirectory successfully created' ;
ELSE ;
   WRITE (*     ,*) 'Failed to create subdirectory' ;
   WRITE (UnInf,*) 'Failed to create subdirectory' ;
END IF ;

Directory = MakeDirQQ (TRIM(AdjustL (OutDir))//'/'//TRIM(AdjustL (ModelName))//'/'//TRIM(AdjustL (AnaName))  ) ;
IF (Directory) THEN ;
   WRITE (*     ,*) 'New subdirectory successfully created' ;
   WRITE (UnInf,*) 'New subdirectory successfully created' ;
ELSE ;
   WRITE (*     ,*) 'Failed to create subdirectory' ;
   WRITE (UnInf,*) 'Failed to create subdirectory' ;
END IF ;

! internal folder
Directory = MakeDirQQ (TRIM(AdjustL (InlDir))//'/'//TRIM(AdjustL (ModelName))) ;
  IF (Directory) THEN ;
     WRITE (*    ,*) 'New subdirectory successfully created' ;
     WRITE (UnInf,*) 'New subdirectory successfully created' ;
  ELSE ;
     WRITE (*    ,*) 'Failed to create subdirectory' ;
     WRITE (UnInf,*) 'Failed to create subdirectory' ;
  END IF ;

Directory = MakeDirQQ (TRIM(AdjustL (InlDir))//'/'//TRIM(AdjustL (ModelName))//'/'//TRIM(AdjustL (AnaName))  ) ;
  IF (Directory) THEN ;
     WRITE (*    ,*) 'New subdirectory successfully created' ;
     WRITE (UnInf,*) 'New subdirectory successfully created' ;
  ELSE ;
     WRITE (*    ,*) 'Failed to create subdirectory' ;
     WRITE (UnInf,*) 'Failed to create subdirectory' ;
  END IF ;

OutDir = TRIM(AdjustL (OutDir))//'/'//TRIM(AdjustL (ModelName))//'/'//TRIM(AdjustL (AnaName)) ;
InlDir = TRIM(AdjustL (InlDir))//'/'//TRIM(AdjustL (ModelName))//'/'//TRIM(AdjustL (AnaName)) ;

! - INFORMATON FILE ---------------------------------------------------------------------------------------------------------------------------------
UnFile = UnInf ;
Open ( Unit = UnFile, FILE = TRIM(AnaName)//'_'//Trim(AdjustL(IndexSize))//'_'//Trim(AdjustL(IndexRank))//'.inf' , ERR =  1001, IOSTAT = IO_File, ACCESS = 'SEQUENTIAL', ACTION = 'Write', ASYNCHRONOUS = 'NO', BLANK = 'NULL', BLOCKSize = 0, DEFAULTFILE = TRIM(OutDir), DisPOSE = 'KEEP', FORM = 'FORMATTED', POSITION = 'ASIS', STATUS = 'REPLACE' ) ;

! - CHECK FILE --------------------------------------------------------------------------------------------------------------------------------------
UnFile = UN_CHK ;
Open ( Unit = UnFile, FILE = TRIM(AnaName)//'_'//Trim(AdjustL(IndexSize))//'_'//Trim(AdjustL(IndexRank))//'.Chk' , ERR =  1001, IOSTAT = IO_File, ACCESS = 'SEQUENTIAL', ACTION = 'Write', ASYNCHRONOUS = 'NO', BLANK = 'NULL', BLOCKSize = 0, DEFAULTFILE = TRIM(InlDir), DisPOSE = 'KEEP', FORM = 'FORMATTED', POSITION = 'ASIS', STATUS = 'REPLACE' ) ;

! Analysis Data file 
UnFile = UnInptAna ;
!Open ( Unit = UnFile, FILE = TRIM(AnaName)//'_'//Trim(AdjustL(IndexSize))//'_'//Trim(AdjustL(IndexRank))//'.txt', ERR =  1001, IOSTAT = IO_File, ACCESS = 'SEQUENTIAL', ACTION = 'READ', ASYNCHRONOUS = 'NO', BLANK = 'NULL', BLOCKSize = 0, DEFAULTFILE = TRIM(Ana_InDir), DisPOSE = 'KEEP', FORM = Trim(AdjustL(Format_Type)), POSITION = 'ASIS', STATUS = 'Old' ) ;
Open ( Unit = UnFile, FILE = TRIM(AnaName)//'.txt', ERR =  1001, IOSTAT = IO_File, ACCESS = 'SEQUENTIAL', ACTION = 'READ', ASYNCHRONOUS = 'NO', BLANK = 'NULL', BLOCKSize = 0, DEFAULTFILE = TRIM(Ana_InDir), DisPOSE = 'KEEP', FORM = Trim(AdjustL(Format_Type)), POSITION = 'ASIS', STATUS = 'Old' ) ;

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
    Call Input (                                                                                                              &
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

! Reading input arrays
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

! =============================================== G matrix in 3D ====================================================================================

! Allocating required arrays for PETSc
Allocate ( O_NNZ_G ( NEQM_Mapping ), D_NNZ_G ( NEQM_Mapping ), STAT = ERR_Alloc ) ;
  IF ( ERR_Alloc /= 0 ) Then ;
    Write (*, Fmt_ALLCT) ERR_Alloc ;  Write (UnInf, Fmt_ALLCT) ERR_Alloc ;
    !#Call BEEP_FAIL ;
    Write(*, Fmt_FL) ;  Write(UnInf, Fmt_FL) ; Write(*, Fmt_End) ; Read(*,*) ;  STOP ;
  End If ;

! =============================================== Reading heterogeneous material properties =========================================================

! Lambda (.Lambda)
UnFile = UnInpt_Lambda ;
Open ( Unit = UnFile, FILE = TRIM(ModelName)//'_'//Trim(AdjustL(IndexSize))//'_'//Trim(AdjustL(IndexRank))//'.Lambda',  ERR =  1001, IOSTAT = IO_File, ACCESS = 'SEQUENTIAL', ACTION ='Read', ASYNCHRONOUS = 'NO', BLANK = 'NULL', BLOCKSIZE = 0, DEFAULTFILE = TRIM(Model_InDir), DisPOSE = 'KEEP', FORM = Trim(AdjustL(Format_Type)), POSITION = 'ASIS', STATUS ='Old') ;

! Mu (.Mu)
UnFile = UnInpt_Mu ;
Open ( Unit = UnFile, FILE = TRIM(ModelName)//'_'//Trim(AdjustL(IndexSize))//'_'//Trim(AdjustL(IndexRank))//'.Mu',  ERR =  1001, IOSTAT = IO_File, ACCESS = 'SEQUENTIAL', ACTION ='Read', ASYNCHRONOUS = 'NO', BLANK = 'NULL', BLOCKSIZE = 0, DEFAULTFILE = TRIM(Model_InDir), DisPOSE = 'KEEP', FORM = Trim(AdjustL(Format_Type)), POSITION = 'ASIS', STATUS ='Old') ;

! allocate arrays
Allocate ( PMat_Lambda (NJ) , PMat_Mu (NJ) )

    Call Input               (                                                                                                      &
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


! =========================== Close FILES ===========================================================================================================
! - Close ADDRESS FILE ------------------------------------------------------------------------------------------------------------------------------
UnFile =  UN_ADR ;
Close ( Unit = UnFile, STATUS = 'KEEP', ERR = 1002, IOSTAT = IO_File )

! - Close Input FILE --------------------------------------------------------------------------------------------------------------------------------
! Data file of model
UnFile =  UnInptMdl ;
Close ( Unit = UnFile, STATUS = 'KEEP', ERR = 1002, IOSTAT = IO_File ) ;

! Coordinate file (.XYZ)
UnFile = UnInptXYZ ;
Close ( Unit = UnFile, STATUS = 'KEEP', ERR = 1002, IOSTAT = IO_File ) ;

! Material Properties (.Mat)
UnFile = UnInptMat ;
Close ( Unit = UnFile, STATUS = 'KEEP', ERR = 1002, IOSTAT = IO_File ) ;

! Connectivity file (.Cnn)
UnFile = UnInptCnn ;
Close ( Unit = UnFile, STATUS = 'KEEP', ERR = 1002, IOSTAT = IO_File ) ;

! Constraint file (.Cnt)
UnFile = UnInptCnt ;
Close ( Unit = UnFile, STATUS = 'KEEP', ERR = 1002, IOSTAT = IO_File ) ;

! Non-Zeros of Stiffness file (.NNZStiff)
UnFile = UnInptStiff ;
Close ( Unit = UnFile, STATUS = 'KEEP', ERR = 1002, IOSTAT = IO_File ) ;

! Non-Zeros of Damping file (.NNZDamp)
UnFile = UnInptDamp ;
Close ( Unit = UnFile, STATUS = 'KEEP', ERR = 1002, IOSTAT = IO_File ) ;

! Non-Zeros of Damping file (.NNZMass)
UnFile = UnInptMass ;
Close ( Unit = UnFile, STATUS = 'KEEP', ERR = 1002, IOSTAT = IO_File ) ;

! Application Ordering file (.App)
UnFile = UnInptApp ;
Close ( Unit = UnFile, STATUS = 'KEEP', ERR = 1002, IOSTAT = IO_File ) ;


! =============================================== Opening external output files ======================================================================
Write (*    ,"('Opening Output files')") ;
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

