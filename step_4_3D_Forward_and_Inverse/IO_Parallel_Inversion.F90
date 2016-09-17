
!****************************************************************************************************************************************************
! This component reads parallel input files required for inversion
! January 12 2014
!****************************************************************************************************************************************************

! =============================================== Reading data structure for inversion ==============================================================

! - Reading Input FILE Basic ------------------------------------------------------------------------------------------------------------------------
UnFile =  Un_Inversion_DS ;
Open ( Unit = UnFile, FILE = TRIM(ModelName)//'_'//Trim(AdjustL(IndexSize))//'_'//Trim(AdjustL(IndexRank))//'.Inv_DS',  ERR =  1001, IOSTAT = IO_File, ACCESS = 'SEQUENTIAL', ACTION ='Read', ASYNCHRONOUS = 'NO', BLANK = 'NULL', BLOCKSIZE = 0, DEFAULTFILE = TRIM(Model_InDir), DisPOSE = 'KEEP', FORM = 'Formatted', POSITION = 'ASIS', STATUS ='Old') ;

    Call Input_Inv_DS_Basic (                                                                                                       &
    !                                                                                                                               & ! Integer (1) Variables
    !                                                                                                                               & ! Integer (2) Variables
    !                                                                                                                               & ! Integer (4) Variables
    NJ_Rank_IParts, NJ_Mapping, NStore_Mapping, NStore_Rank                                                                         & ! Integer (8) Variables
    !                                                                                                                               & ! Real Variables
    !                                                                                                                               &
    !                                                                                                                               &
    !                                                                                                                               & ! Integer Arrays
    !                                                                                                                               & ! Real Arrays
    !                                                                                                                               & ! Characters
    !                                                                                                                               & ! Type
    ) ;

! Allocating required arrays for PETSc --------------------------------------------------------------------------------------------------------------
Allocate ( idx_Mat_from ( NJ_Rank_IParts ), idx_Mat_to ( NJ_Rank_IParts ), idx_Mat_Extend_from ( NJ_Mapping ), idx_Mat_Extend_to ( NJ_Mapping ), U_Store_Numbers_Global ( NStore_Mapping ), idx_u_from ( NStore_Rank ), idx_u_to ( NStore_Rank ), STAT = ERR_Alloc ) ;
  IF ( ERR_Alloc /= 0 ) Then ;
    Write (*, Fmt_ALLCT) ERR_Alloc ;  Write (UnInf, Fmt_ALLCT) ERR_Alloc ;
    !#Call BEEP_FAIL ;
    Write(*, Fmt_FL) ;  Write(UnInf, Fmt_FL) ; Write(*, Fmt_End) ; Read(*,*) ;  STOP ;
  End If ;

! - Reading Input FILE Array ------------------------------------------------------------------------------------------------------------------------
    Call Input_Inv_DS_Array (                                                                                                       &
    !                                                                                                                               & ! Integer (1) Variables
    !                                                                                                                               & ! Integer (2) Variables
    !                                                                                                                               & ! Integer (4) Variables
    NJ_Rank_IParts, NJ_Mapping, NStore_Mapping, NStore_Rank,                                                                        & ! Integer (8) Variables
    !                                                                                                                               & ! Real Variables
    !                                                                                                                               &
    !                                                                                                                               &
    idx_Mat_from, idx_Mat_to, idx_Mat_Extend_from, idx_Mat_Extend_to, U_Store_Numbers_Global, idx_u_from, idx_u_to                  & ! Integer Arrays
    !                                                                                                                               & ! Real Arrays
    !                                                                                                                               & ! Characters
    !                                                                                                                               & ! Type
    ) ;

! checked: we are importing idx_Mat_Extend_from, idx_Mat_Extend_to correctly.

! Allocating required arrays for measured response at select sensor locations  ----------------------------------------------------------------------
NNDH = Param%IntM( 6, 1)

Allocate ( Dis_meas ( 0:NStep , NNDH * NDim ), STAT = ERR_Alloc ) ;
  IF ( ERR_Alloc /= 0 ) Then ;
    Write (*, Fmt_ALLCT) ERR_Alloc ;  Write (UnInf, Fmt_ALLCT) ERR_Alloc ;
    !#Call BEEP_FAIL ;
    Write(*, Fmt_FL) ;  Write(UnInf, Fmt_FL) ; Write(*, Fmt_End) ; Read(*,*) ;  STOP ;
  End If ;

! - Create the result folder for inversion ----------------------------------------------------------------------------------------------------------
! measured response
Directory = MakeDirQQ (TRIM(AdjustL (OutDir))//'/Measured_Response'  ) ;
IF (Directory) THEN ;
   IF ( I_PRINT_SCREEN == 1 )   WRITE (*     ,*) 'New subdirectory successfully created' ;
   WRITE (UnInf,*) 'New subdirectory successfully created' ;
ELSE ;
   IF ( I_PRINT_SCREEN == 1 )   WRITE (*     ,*) 'Failed to create subdirectory' ;
   WRITE (UnInf,*) 'Failed to create subdirectory' ;
END IF ;

! material updates
Directory = MakeDirQQ (TRIM(AdjustL (OutDir))//'/Material_Updates' ) ;
IF (Directory) THEN ;
   IF ( I_PRINT_SCREEN == 1 )   WRITE (*     ,*) 'New subdirectory successfully created' ;
   WRITE (UnInf,*) 'New subdirectory successfully created' ;
ELSE ;
   IF ( I_PRINT_SCREEN == 1 )   WRITE (*     ,*) 'Failed to create subdirectory' ;
   WRITE (UnInf,*) 'Failed to create subdirectory' ;
END IF ;

MatUpdateDir = TRIM(AdjustL (OutDir))//'/Material_Updates'
MeasRespDir  = TRIM(AdjustL (OutDir))//'/Measured_Response'

! Measured response at select sensor locations to form the misfit.
UnFile = UN_Measured_HisD ;
Open ( Unit = UnFile, FILE = TRIM(ModelNAME)//'_'//Trim(AdjustL(IndexSize))//'_'//Trim(AdjustL(IndexRank))//'.HisD.txt', ERR =  1001, IOSTAT = IO_File, ACCESS = 'SEQUENTIAL', ACTION ='Read', ASYNCHRONOUS = 'NO', BLANK = 'NULL', BLOCKSIZE = 0, DEFAULTFILE = TRIM(MeasRespDir), DisPOSE = 'KEEP', FORM = 'FORMATTED', POSITION = 'ASIS', STATUS ='OLD') ;
!goto 10010
! - Reading Measured response at select sensor locations to form the misfit -------------------------------------------------------------------------
    Call Input_Inv_Measured_Data (                                                                                                  &
    !                                                                                                                               & ! Integer (1) Variables
    !                                                                                                                               & ! Integer (2) Variables
    !                                                                                                                               & ! Integer (4) Variables
    NDim, NNDH,                                                                                                                     & ! Integer (8) Variables
    !                                                                                                                               & ! Real Variables
    !                                                                                                                               &
    !                                                                                                                               &
    !                                                                                                                               & ! Integer Arrays
    Dis_meas                                                                                                                        & ! Real Arrays
    !                                                                                                                               & ! Characters
    !                                                                                                                               & ! Type
    ) ;
!10010 continue
! Allocating required arrays for the regularization matrix (discrete Laplace operator) --------------------------------------------------------------
Allocate ( D_NNZ_Reg ( NJ_Mapping ), O_NNZ_Reg ( NJ_Mapping ), STAT = ERR_Alloc ) ;
  IF ( ERR_Alloc /= 0 ) Then ;
    Write (*, Fmt_ALLCT) ERR_Alloc ;  Write (UnInf, Fmt_ALLCT) ERR_Alloc ;
    !#Call BEEP_FAIL ;
    Write(*, Fmt_FL) ;  Write(UnInf, Fmt_FL) ; Write(*, Fmt_End) ; Read(*,*) ;  STOP ;
  End If ;


! =============================================== Process the inversion data structure ==============================================================

! - Construct a mapping between EqDis and U_Store_Numbers_Global needed for misfit computation ------------------------------------------------------

! Construct the mapping by exhaustive search.
! Mapping - left: equation index for misfit - right: equation index for displacement history storage
!  1  > 1
!  2  > 205
!  3  > 307

Allocate ( EqDis_MAP_U_Store_Numbers_Global (NNDH * NDim) )

Do IEq = 1 , NNDH * NDim
  Do IStore = 1 , NStore_Mapping
    If ( EqDis (IEq) == U_Store_Numbers_Global (IStore) ) Then
       EqDis_MAP_U_Store_Numbers_Global (IEq) = IStore
       Exit
    End If
  End Do
End Do


!================================================ Monitoring of the progress of inversion iterations ================================================
! - store iteration information
Directory = MakeDirQQ (TRIM(AdjustL (OutDir))//'/Iter_Monitor'   ) ; 
IterDir   = TRIM(AdjustL (OutDir))//'/Iter_Monitor'

Write (Iter_begin_char, *) Iter_begin

If ( Rank == 0 ) Then ;
   UnFile = UN_Inversion_Iterations ;
   Open ( Unit = UnFile, FILE = TRIM(AnaName)//'_'//'Iter'//'_'//Trim(AdjustL(Iter_begin_char))//'.txt',  ERR =  1001, IOSTAT = IO_File, ACCESS = 'SEQUENTIAL', ACTION ='Write', ASYNCHRONOUS = 'NO', BLANK = 'NULL', BLOCKSIZE = 0, DEFAULTFILE = TRIM(IterDir), DisPOSE = 'KEEP', FORM = 'FORMATTED', POSITION = 'ASIS', STATUS ='REPLACE') ;
End If


! =============================================== Save material updates to resume inversion if interrupted ==========================================

! material resume - create folder
Directory = MakeDirQQ (TRIM(AdjustL (OutDir))//'/Material_Resume'   ) ; 
Directory = MakeDirQQ (TRIM(AdjustL (OutDir))//'/Material_Resume/1' ) ; ! Iter: 10, 30, 50, ....
Directory = MakeDirQQ (TRIM(AdjustL (OutDir))//'/Material_Resume/2' ) ; ! Iter: 20, 40, 60, ....
MatResumeDir_1 = TRIM(AdjustL (OutDir))//'/Material_Resume/1'
MatResumeDir_2 = TRIM(AdjustL (OutDir))//'/Material_Resume/2'



















