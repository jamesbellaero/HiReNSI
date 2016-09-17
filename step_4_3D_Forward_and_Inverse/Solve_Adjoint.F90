
! The past always looks better than it was. It's only pleasant because it isn't here.
! ---------------------------------------------------------------------------------------------------------------------------------------------------
! Solve the adjoint problem - Optimize-then-Discretize (because DtO was not going anywhere.)
! ---------------------------------------------------------------------------------------------------------------------------------------------------

! 1- initialize parallel sparce matrices
  Call MatZeroEntries ( AM_PETSC, IERR )
  Call MatZeroEntries ( AC_PETSC, IERR )
  Call MatZeroEntries ( AK_PETSC, IERR )
  Call MatZeroEntries ( AG_PETSC, IERR )
  Call MatZeroEntries ( AK_RD_PETSC, IERR )

! 2- initialize parallel vectors
  Call VecSet ( Diag_M_PETSC,    0.0d0, IERR )
  Call VecSet ( Diag_M_RD_PETSC, 0.0d0, IERR )

! 3- assemble regular domain and PML matrices
  IF (RANK == 0 .AND. I_PRINT_SCREEN == 1) WRITE(*,*) 'Assemble RD Adjoint'
  Call ASSEM_REGULAR_DOMAIN_PETSC_3D ( AK_PETSC, AM_PETSC, AC_PETSC, AK_RD_PETSC, AM_RD_PETSC, RMJ_PETSC, DIAG_M_PETSC, DIAG_M_RD_PETSC, XYZ, INOD, NGP, MTEL, PMAT, PMat_Lambda, PMat_Mu, PML_PARAM, ID_BC, ID, SIZE, RANK )

  IF (RANK == 0 .AND. I_PRINT_SCREEN == 1) WRITE(*,*) 'Assemble PML Adjoint'
  Call ASSEM_PML_PETSC_3D_Adjoint ( AK_PETSC, AM_PETSC, AC_PETSC, AG_PETSC, DIAG_M_PETSC, XYZ, INOD, NGP, MTEL, PMAT, PMat_Lambda, PMat_Mu, PML_PARAM, ID_BC, ID, SIZE, RANK )

  IF (RANK == 0 .AND. I_PRINT_SCREEN == 1) WRITE(*,*) 'MPI Assembly'
  CALL MatAssemblyBegin ( AK_PETSC, MAT_FINAL_ASSEMBLY, IERR )
  CALL MatAssemblyBegin ( AM_PETSC, MAT_FINAL_ASSEMBLY, IERR )
  CALL MatAssemblyBegin ( AC_PETSC, MAT_FINAL_ASSEMBLY, IERR )
  CALL MatAssemblyBegin ( AG_PETSC, MAT_FINAL_ASSEMBLY, IERR )
  CALL MatAssemblyEnd   ( AK_PETSC, MAT_FINAL_ASSEMBLY, IERR )
  CALL MatAssemblyEnd   ( AM_PETSC, MAT_FINAL_ASSEMBLY, IERR )
  CALL MatAssemblyEnd   ( AC_PETSC, MAT_FINAL_ASSEMBLY, IERR )
  CALL MatAssemblyEnd   ( AG_PETSC, MAT_FINAL_ASSEMBLY, IERR )

  CALL VecAssemblyBegin ( DIAG_M_PETSC   , IERR )
  CALL VecAssemblyBegin ( DIAG_M_RD_PETSC, IERR )
  CALL VecAssemblyEnd   ( DIAG_M_PETSC   , IERR )
  CALL VecAssemblyEnd   ( DIAG_M_RD_PETSC, IERR )

  IF (RANK == 0 .AND. I_PRINT_SCREEN == 1) WRITE(*,*) 'Adjoint: Optimize-then-Discretize'
  Call EXPLICIT_RK4_PETSC_PML_3D_Adjoint_OtD ( AK_PETSC, DIAG_M_PETSC, AC_PETSC, AG_PETSC, B_mis_PETSC, AK_RD_PETSC, DIAG_M_RD_PETSC, ID, &
                                               NNDH, EqDis, U_Store_Mapping, P_Store_Mapping, NStore_Mapping, U_Store_Numbers_Global, EqDis_MAP_U_Store_Numbers_Global, Dis_meas )

