
! =============================================== Save material updates to resume inversion if interrupted ==========================================
! Comment: we record at every 10 iterations at two different files: Iter 10 (1) - 20 (2) - 30 (1) - 40 (2) - 50 (1) - 60 (2) ...

   Iter_File_10 = Mod ( Iter , 10 )
   Iter_File_20 = Mod ( Iter , 20 )
   I_Flag_Disk  = 0

   If ( ( Iter_File_10 == 0 ) .AND. ( Iter_File_20 /= 0 ) ) Then
      UnFile_L     = UnOut_Lambda_1 ;
      UnFile_M     = UOut_Mu_1 ;
      MatResumeDir = MatResumeDir_1
      I_Flag_Disk  = 1
   Else If ( Iter_File_20 == 0 ) Then
      UnFile_L     = UnOut_Lambda_2 ;
      UnFile_M     = UnOut_Mu_2 ;
      MatResumeDir = MatResumeDir_2
      I_Flag_Disk  = 1
   End If


! - Writing down the input data for each process ----------------------------------------------------------------------------------------------------
   If ( I_Flag_Disk  == 1 ) Then

    ! - Output FILEs --------------------------------------------------------------------------------------------------------------------------------
         Open ( Unit = UnFile_L, FILE = TRIM(ModelName)//'_'//Trim(AdjustL(IndexSize))//'_'//Trim(AdjustL(IndexRank))//'.Lambda', ERR =  1001, IOSTAT = IO_File, ACCESS = 'SEQUENTIAL', ACTION ='Write', ASYNCHRONOUS = 'NO', BLANK = 'NULL', BLOCKSIZE = 0, DEFAULTFILE = TRIM(MatResumeDir), DisPOSE = 'KEEP', FORM = 'FORMATTED', POSITION = 'ASIS', STATUS ='REPLACE') ;
         Open ( Unit = UnFile_M, FILE = TRIM(ModelName)//'_'//Trim(AdjustL(IndexSize))//'_'//Trim(AdjustL(IndexRank))//'.Mu',     ERR =  1001, IOSTAT = IO_File, ACCESS = 'SEQUENTIAL', ACTION ='Write', ASYNCHRONOUS = 'NO', BLANK = 'NULL', BLOCKSIZE = 0, DEFAULTFILE = TRIM(MatResumeDir), DisPOSE = 'KEEP', FORM = 'FORMATTED', POSITION = 'ASIS', STATUS ='REPLACE') ;

    ! - Write Lambda & Mu ---------------------------------------------------------------------------------------------------------------------------
         If ( Rank == 0 .AND. I_PRINT_SCREEN == 1 ) Write (*,*) "Storing material properties ..." ;

            If ( Lambda_variable == 'Lambda_2Mu' ) Then
               DO IJ = 1 , NJ ;
                  Write ( Unit = UnFile_L, FMT = "(<NDim>(E31.23E3,2X) )", ADVANCE = 'YES', ASYNCHRONOUS = 'NO', IOSTAT = IO_Write, ERR = 1006 ) PMat_Lambda ( IJ ) - 2.0d0 * PMat_Mu ( IJ );
               End Do ;
            Else
               DO IJ = 1 , NJ ;
                  Write ( Unit = UnFile_L, FMT = "(<NDim>(E31.23E3,2X) )", ADVANCE = 'YES', ASYNCHRONOUS = 'NO', IOSTAT = IO_Write, ERR = 1006 ) PMat_Lambda ( IJ ) ;
               End Do ;
            End If

            DO IJ = 1 , NJ ;
                Write ( Unit = UnFile_M, FMT = "(<NDim>(E31.23E3,2X) )", ADVANCE = 'YES', ASYNCHRONOUS = 'NO', IOSTAT = IO_Write, ERR = 1006 ) PMat_Mu ( IJ ) ; 
            End Do ;

    ! - Closing the files ---------------------------------------------------------------------------------------------------------------------------
         Close ( Unit = UnFile_L, Status = 'KEEP', ERR =  1002, IOSTAT = IO_File ) ;
         Close ( Unit = UnFile_M, Status = 'KEEP', ERR =  1002, IOSTAT = IO_File ) ;

   End If
