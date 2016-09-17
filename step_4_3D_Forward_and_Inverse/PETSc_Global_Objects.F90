
!==================================================================================================
! - DEFINE PETSC OBJECTS (Forward Problem)

      !Write (*    ,"('Creating Index Set ...')") ;
      !Write (UnInf,"('Creating Index Set ...')") ;

      !Call AOCreateBasic ( PETSC_COMM_WORLD, NEQM_Mapping, App_Numbers, PETSc_Numbers, AOOut, ErrPTC ) ;
      !Call AOView( AOOut, PETSC_VIEWER_STDOUT_WORLD , ErrPTC ) ;     ! comment out later

      ! - Local to Global Mapping -------------------------------------------------------------------------------------------------------------------
      Call ISLocalToGlobalMappingCreate ( PETSC_COMM_WORLD, NEqM, Indices, PETSC_COPY_VALUES, GMapping, ErrPTC ) ;
      !Call ISLocalToGlobalMappingView ( GMapping, PETSc_VIEWER_STDOUT_WORLD, ErrPTC ) ;

      IF ( I_PRINT_SCREEN == 1 )      Write (*    ,"('Index Set created successfully')") ;
      Write (UnInf,"('Index Set created successfully')") ;

      ! - Defining PETSc objects ----------------------------------------------------------------------------------------------------------------------
      ! Matrices
      IF ( I_PRINT_SCREEN == 1 )     Write (*    ,"('Creating PETSc objects - Matrices ...')") ;
      Write (UnInf,"('Creating PETSc objects - Matrices ...')") ;

Select Case (NDim)

   Case (2)
      D_NNZ_Stiff (:) = 200 ;   ! Delete this part after finalizing the ParaPIC code
      O_NNZ_Stiff (:) = 200 ;
      D_NNZ_Damp  (:) = 200 ;
      O_NNZ_Damp  (:) = 200 ;
      D_NNZ_Mass  (:) = 200 ;
      O_NNZ_Mass  (:) = 200 ;

   Case (3)
      D_NNZ_Stiff (:) = 200 ;   ! Delete this part after finalizing the ParaPIC code
      O_NNZ_Stiff (:) = 200 ;
      D_NNZ_Damp  (:) = 200 ;
      O_NNZ_Damp  (:) = 200 ;
      D_NNZ_Mass  (:) = 200 ;
      O_NNZ_Mass  (:) = 200 ;
      D_NNZ_G     (:) = 200 ;
      O_NNZ_G     (:) = 200 ;

End Select

      ! General Matrices
      CALL MatCreateAIJ ( PETSC_COMM_WORLD, NEQM_Mapping, NEQM_Mapping, NEQMTotal, NEQMTotal,      0, D_NNZ_Mass,         0,        O_NNZ_Mass,  AM_PETSC,    IERR ) ;
      CALL MatCreateAIJ ( PETSC_COMM_WORLD, NEQM_Mapping, NEQM_Mapping, NEQMTotal, NEQMTotal,      0, D_NNZ_Damp,         0,        O_NNZ_Damp,  AC_PETSC,    IERR ) ;
      CALL MatCreateAIJ ( PETSC_COMM_WORLD, NEQM_Mapping, NEQM_Mapping, NEQMTotal, NEQMTotal,      0, D_NNZ_Stiff,        0,        O_NNZ_Stiff, AK_PETSC,    IERR ) ;
      CALL MatCreateAIJ ( PETSC_COMM_WORLD, NEQM_Mapping, NEQM_Mapping, NEQMTotal, NEQMTotal,      0, D_NNZ_G,            0,        O_NNZ_G,     AG_PETSC,    IERR ) ;
      CALL MatCreateAIJ ( PETSC_COMM_WORLD, NEQM_Mapping, NEQM_Mapping, NEQMTotal, NEQMTotal,      0, D_NNZ_Mass,         0,        O_NNZ_Mass,  AM_RD_PETSC, IERR ) ;
      CALL MatCreateAIJ ( PETSC_COMM_WORLD, NEQM_Mapping, NEQM_Mapping, NEQMTotal, NEQMTotal,      0, D_NNZ_Stiff,        0,        O_NNZ_Stiff, AK_RD_PETSC, IERR ) ;

      ! MATRICES' OPTIONS
!      CALL MatSetOption    ( AM_PETSC,    MAT_STRUCTURALLY_SYMMETRIC, PETSC_TRUE, IERR ) ;
!      CALL MatSetOption    ( AC_PETSC,    MAT_STRUCTURALLY_SYMMETRIC, PETSC_TRUE, IERR ) ;
!      CALL MatSetOption    ( AK_PETSC,    MAT_STRUCTURALLY_SYMMETRIC, PETSC_TRUE, IERR ) ;
!      CALL MatSetOption    ( AG_PETSC,    MAT_STRUCTURALLY_SYMMETRIC, PETSC_TRUE, IERR ) ;
!      CALL MatSetOption    ( AM_RD_PETSC, MAT_STRUCTURALLY_SYMMETRIC, PETSC_TRUE, IERR ) ;
!      CALL MatSetOption    ( AK_RD_PETSC, MAT_STRUCTURALLY_SYMMETRIC, PETSC_TRUE, IERR ) ;

      CALL MatSetOption    ( AM_PETSC,    MAT_IGNORE_ZERO_ENTRIES,    PETSC_TRUE, IERR ) ;
      CALL MatSetOption    ( AC_PETSC,    MAT_IGNORE_ZERO_ENTRIES,    PETSC_TRUE, IERR ) ;
      CALL MatSetOption    ( AK_PETSC,    MAT_IGNORE_ZERO_ENTRIES,    PETSC_TRUE, IERR ) ;
      CALL MatSetOption    ( AG_PETSC,    MAT_IGNORE_ZERO_ENTRIES,    PETSC_TRUE, IERR ) ;
      CALL MatSetOption    ( AM_RD_PETSC, MAT_IGNORE_ZERO_ENTRIES,    PETSC_TRUE, IERR ) ;
      CALL MatSetOption    ( AK_RD_PETSC, MAT_IGNORE_ZERO_ENTRIES,    PETSC_TRUE, IERR ) ;

      ! MALLOC ERROR
      Call MatSetOption ( AC_PETSC, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE, IERR )
      Call MatSetOption ( AK_PETSC, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE, IERR )
      Call MatSetOption ( AG_PETSC, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE, IERR )

     ! Applying Local to Global Numbering for PETSc objects
      Call MatSetLocalToGlobalMapping   ( AM_PETSC,      GMapping,   GMapping, ErrPTC ) ; 
      Call MatSetLocalToGlobalMapping   ( AC_PETSC,      GMapping,   GMapping, ErrPTC ) ;
      Call MatSetLocalToGlobalMapping   ( AK_PETSC,      GMapping,   GMapping, ErrPTC ) ; 
      Call MatSetLocalToGlobalMapping   ( AG_PETSC,      GMapping,   GMapping, ErrPTC ) ;
      Call MatSetLocalToGlobalMapping   ( AM_RD_PETSC,   GMapping,   GMapping, ErrPTC ) ; 
      Call MatSetLocalToGlobalMapping   ( AK_RD_PETSC,   GMapping,   GMapping, ErrPTC ) ; 

      IF ( I_PRINT_SCREEN == 1 )   Write (*    ,"('PETSc matrices were created successfully')") ;
      Write (UnInf,"('PETSc matrices were created successfully')") ;

      ! VECTORS
      IF ( I_PRINT_SCREEN == 1 )   Write (*    ,"('Crreating PETSc objects - Vectors ...')") ;
      Write (UnInf,"('Crreating PETSc objects - Vectors ...')") ;

      ! Mass vector
      CALL VecCreateMPI    ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, DIAG_M_PETSC   , IERR ) ;
      CALL VecCreateMPI    ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, DIAG_iM_PETSC  , IERR ) ;
      CALL VecCreateMPI    ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, DIAG_M_RD_PETSC, IERR ) ;

      CALL VecSetOption    ( DIAG_M_PETSC   , VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE  , IERR ) ;
      CALL VecSetOption    ( DIAG_M_RD_PETSC, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE  , IERR ) ;

      Call VecSetLocalToGlobalMapping ( DIAG_M_PETSC,    GMapping, ErrPTC ) ;
      Call VecSetLocalToGlobalMapping ( DIAG_M_RD_PETSC, GMapping, ErrPTC ) ;

      CALL VecSet          ( DIAG_M_PETSC   , 0.0d0, IERR ) ;
      CALL VecSet          ( DIAG_M_RD_PETSC, 0.0d0, IERR ) ;

      ! Load vector
      CALL VecCreateMPI    ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, B_PETSC        , IERR ) ;
      CALL VecCreateMPI    ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, RMJ_PETSC      , IERR ) ;

      CALL VecSetOption    ( B_PETSC,   VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE  , IERR ) ;
      CALL VecSetOption    ( RMJ_PETSC, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE  , IERR ) ;

      Call VecSetLocalToGlobalMapping ( B_PETSC,   GMapping, ErrPTC ) ;
      Call VecSetLocalToGlobalMapping ( RMJ_PETSC, GMapping, ErrPTC ) ;

      CALL VecSet          ( B_PETSC        , 0.0d0, IERR ) ;
      CALL VecSet          ( RMJ_PETSC      , 0.0d0, IERR ) ;

      IF ( I_PRINT_SCREEN == 1 )   Write (*    ,"('PETSc Vectors were created successfully')") ;
      Write (UnInf,"('PETSc Vectors were created successfully')") ;

!==================================================================================================



