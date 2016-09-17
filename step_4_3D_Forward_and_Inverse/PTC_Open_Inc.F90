
! - Input FILE --------------------------------------------------------------------------------------------------------------------------------------

Write (*,*) "Model_InDir: ", Model_InDir ;

! Format of the input files ; Binary, Formatted, HDF5
  If ( Output_Type == 2 ) Then ;  ! HDF5 format

  Else ;   ! Binary, Formatted
      If      ( Output_Type == 0 ) Then ; Format_Type = "Formatted" ;
      Else If ( Output_Type == 1 ) Then ; Format_Type = "Binary" ;
      End If ;

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

  End If ;

! - CHECK FILE --------------------------------------------------------------------------------------------------------------------------------------
UnFile = UN_CHK ;
Open ( Unit = UnFile, FILE = TRIM(AnaName)//'_'//Trim(AdjustL(IndexSize))//'_'//Trim(AdjustL(IndexRank))//'.Chk' , ERR =  1001, IOSTAT = IO_File, ACCESS = 'SEQUENTIAL', ACTION = 'Write', ASYNCHRONOUS = 'NO', BLANK = 'NULL', BLOCKSize = 0, DEFAULTFILE = TRIM(InlDir), DisPOSE = 'KEEP', FORM = 'FORMATTED', POSITION = 'ASIS', STATUS = 'REPLACE' ) ;

! Analysis Data file 
UnFile = UnInptAna ;
Open ( Unit = UnFile, FILE = TRIM(AnaName)//'.txt', ERR =  1001, IOSTAT = IO_File, ACCESS = 'SEQUENTIAL', ACTION = 'READ', ASYNCHRONOUS = 'NO', BLANK = 'NULL', BLOCKSize = 0, DEFAULTFILE = TRIM(Ana_InDir), DisPOSE = 'KEEP', FORM = 'FORMATTED', POSITION = 'ASIS', STATUS = 'Old' ) ;
