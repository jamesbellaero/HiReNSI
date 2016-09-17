
! Try to leave the Earth a better place than when you arrived.
!====================================================================================================================================================

! - Output Lambda_PETSC, Mu_PETSC, Iter for Paraview ------------------------------------------------------------------------------------------------

! - 1) use repeatedly after each material update
   Call VecScatterBegin ( vscat_Mat_Vis, Lambda_PETSC, Lambda_Rank_Vis_PETSC, INSERT_VALUES, SCATTER_FORWARD, ErrPTC );
   Call VecScatterEnd   ( vscat_Mat_Vis, Lambda_PETSC, Lambda_Rank_Vis_PETSC, INSERT_VALUES, SCATTER_FORWARD, ErrPTC );

   Call VecScatterBegin ( vscat_Mat_Vis,     Mu_PETSC,     Mu_Rank_Vis_PETSC, INSERT_VALUES, SCATTER_FORWARD, ErrPTC );
   Call VecScatterEnd   ( vscat_Mat_Vis,     Mu_PETSC,     Mu_Rank_Vis_PETSC, INSERT_VALUES, SCATTER_FORWARD, ErrPTC );

! - 2) get
   Call VecGetArrayF90  ( Lambda_Rank_Vis_PETSC, Lambda_Rank_Vis, ErrPTC ) ;
   Call VecGetArrayF90  (     Mu_Rank_Vis_PETSC,     Mu_Rank_Vis, ErrPTC ) ;

! - 3) Write repeatedly before each call to Paraview
   If ( Rank == 0 ) Then ;  ! 
      Write (IndexIter, *) Iter ;
      Do I = 0, Size-1 ;
         Write (IndexRankTemp, *) I ;
         UnFile = UN_OutallWrp_Mat_Vis ;
         Write (UnFile, * ) "     <xi:include href='",TRIM(AnaNAME)//'_'//Trim(AdjustL(IndexSize))//'_'//Trim(AdjustL(IndexRankTemp))//'_'//Trim(AdjustL(IndexIter))//'.xmf',"'/>" ;
      End Do ;
   End If ;

! - 4) You may want to eliminate this later, since time is really for wavemotion visualization
   Time_n1 = Dble(Iter)

! - 5) Fasten your seat belt, and call Paraview for material visualization
! - -------------------------------------------------------------------------------------------------------------------------------------------------
   Call Paraview_Mat_HDF5 (                                                                                                        &
   !                                                                                                                               & ! Integer (1) Variables
   !                                                                                                                               & ! Integer (2) Variables
   !                                                                                                                               & ! Integer (4) Variables
   Iter, NJ_Para, Param%IntM( 5, 4),                                                                                               & ! Integer (8) Variables
   Time_n1,                                                                                                                        & ! Real Variables
   !                                                                                                                               & ! Integer Arrays
   Lambda_Rank_Vis, Mu_Rank_Vis,                                                                                                   & ! Real Arrays
   MatUpdateDir, ModelName, AnaName, IndexRank, IndexSize                                                                          & ! Characters
   !                                                                                                                               & ! Type
   ) ;
! - -------------------------------------------------------------------------------------------------------------------------------------------------

! - 6) restore
   Call VecRestoreArrayF90  ( Lambda_Rank_Vis_PETSC, Lambda_Rank_Vis, ErrPTC ) ;
   Call VecRestoreArrayF90  (     Mu_Rank_Vis_PETSC,     Mu_Rank_Vis, ErrPTC ) ;

!====================================================================================================================================================
