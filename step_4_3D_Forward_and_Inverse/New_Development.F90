
! =============================================== Read material updates to resume inversion when interrupted ========================================
! Comment: we record material properties at every 10 iterations, at two different files: Iter 10 (1) - 20 (2) - 30 (1) - 40 (2) - 50 (1) - 60 (2) ...

   Iter_File_10 = Mod ( Iter-1 , 10 )
   Iter_File_20 = Mod ( Iter-1 , 20 )

   If ( ( Iter_File_10 == 0 ) .AND. ( Iter_File_20 /= 0 ) ) Then
      UnFile_L     = UnOut_Lambda_1 ;
      UnFile_M     = UOut_Mu_1 ;
      MatResumeDir = MatResumeDir_1
   Else If ( Iter_File_20 == 0 ) Then
      UnFile_L     = UnOut_Lambda_2 ;
      UnFile_M     = UnOut_Mu_2 ;
      MatResumeDir = MatResumeDir_2
   End If


! - Read heterogeneous material properties to resume inversion --------------------------------------------------------------------------------------
   If ( Iter /= 1 ) Then

    ! - Input FILEs --------------------------------------------------------------------------------------------------------------------------------
         Open ( Unit = UnFile_L, FILE = TRIM(ModelName)//'_'//Trim(AdjustL(IndexSize))//'_'//Trim(AdjustL(IndexRank))//'.Lambda', ERR =  1001, IOSTAT = IO_File, ACCESS = 'SEQUENTIAL', ACTION ='Read', ASYNCHRONOUS = 'NO', BLANK = 'NULL', BLOCKSIZE = 0, DEFAULTFILE = TRIM(MatResumeDir), DisPOSE = 'KEEP', FORM = 'FORMATTED', POSITION = 'ASIS', STATUS ='Old') ;
         Open ( Unit = UnFile_M, FILE = TRIM(ModelName)//'_'//Trim(AdjustL(IndexSize))//'_'//Trim(AdjustL(IndexRank))//'.Mu',     ERR =  1001, IOSTAT = IO_File, ACCESS = 'SEQUENTIAL', ACTION ='Read', ASYNCHRONOUS = 'NO', BLANK = 'NULL', BLOCKSIZE = 0, DEFAULTFILE = TRIM(MatResumeDir), DisPOSE = 'KEEP', FORM = 'FORMATTED', POSITION = 'ASIS', STATUS ='Old') ;

    ! - Read Lambda & Mu ---------------------------------------------------------------------------------------------------------------------------
         If ( Rank == 0 .AND. I_PRINT_SCREEN == 1 ) Write (*,*) "Reading material properties ..." ;

         Call Input_Het_Mat_Resume (                                                                                                     &
         !                                                                                                                               & ! Integer (1) Variables
         !                                                                                                                               & ! Integer (2) Variables
         UnFile_L, UnFile_M,                                                                                                             & ! Integer (4) Variables
         NJ,                                                                                                                             & ! Integer (8) Variables
         !                                                                                                                               & ! Real Variables
         !                                                                                                                               &
         !                                                                                                                               &
         !                                                                                                                               & ! Integer Arrays
         PMat_Lambda, PMat_Mu                                                                                                            & ! Real Arrays
         !                                                                                                                               & ! Characters
         !                                                                                                                               & ! Type
         ) ;

    ! - Closing the files ---------------------------------------------------------------------------------------------------------------------------
         Close ( Unit = UnFile_L, Status = 'KEEP', ERR =  1002, IOSTAT = IO_File ) ;
         Close ( Unit = UnFile_M, Status = 'KEEP', ERR =  1002, IOSTAT = IO_File ) ;

   End If

