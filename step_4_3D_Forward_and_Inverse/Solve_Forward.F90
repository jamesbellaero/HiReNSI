
! Computers are useless. They can only give you answers.

! ===================================================================================================================================================
! Solve the Forward Problem
! ===================================================================================================================================================

!---------- ---------- ---------- ---------- ----------  
CALL CPU_TIME(t1)
!IF (RANK == 0) WRITE(*,'(A40,F10.2,A10)')'pre-processing time:', t1-t0, 'seconds'
!IF (RANK == 0) WRITE(2,'(A40,F10.2,A10)')'pre-processing time:', t1-t0, 'seconds'
!IF (RANK == 0) WRITE(*,'(A88)')'---------- ---------- ---------- ---------- ---------- ---------- ---------- ---------- ---------- ----------'

SELECT CASE (NDIM)
   CASE(2)
      CALL ASSEM_REGULAR_DOMAIN_PETSC    ( AK_PETSC, AM_PETSC, AC_PETSC, AK_RD_PETSC, AM_RD_PETSC, RMJ_PETSC, DIAG_M_PETSC, DIAG_M_RD_PETSC, XYZ, INOD, NGP, MTEL, PMAT, PMat_Lambda, PMat_Mu, PML_PARAM, ID_BC, ID, SIZE, RANK )
   CASE(3)
      CALL ASSEM_REGULAR_DOMAIN_PETSC_3D ( AK_PETSC, AM_PETSC, AC_PETSC, AK_RD_PETSC, AM_RD_PETSC, RMJ_PETSC, DIAG_M_PETSC, DIAG_M_RD_PETSC, XYZ, INOD, NGP, MTEL, PMAT, PMat_Lambda, PMat_Mu, PML_PARAM, ID_BC, ID, SIZE, RANK )
END SELECT

CALL CPU_TIME(t2)
IF (RANK == 0 .AND. I_PRINT_SCREEN == 1) WRITE(*,'(A40,F10.2,A10)')'regular domain assemble time:', t2-t1, 'seconds'
IF (RANK == 0 .AND. I_PRINT_SCREEN == 1) WRITE(2,'(A40,F10.2,A10)')'regular domain assemble time:', t2-t1, 'seconds'
IF (RANK == 0 .AND. I_PRINT_SCREEN == 1) WRITE(*,'(A88)')'---------- ---------- ---------- ---------- ---------- ---------- ---------- ---------- ---------- ----------'

SELECT CASE (NDIM)
   CASE(2)
      CALL ASSEM_FORCE_PETSC    ( B_PETSC, XYZ, INOD, NGP, MTEL, PMAT, ID_BC, ID, SIZE, RANK )
   CASE(3)
      CALL ASSEM_FORCE_PETSC_3D ( B_PETSC, XYZ, INOD, NGP, MTEL, PMAT, ID_BC, ID, SIZE, RANK )
END SELECT

CALL CPU_TIME(t3)
IF (RANK == 0 .AND. I_PRINT_SCREEN == 1) WRITE(*,'(A40,F10.2,A10)')'force assembly time:', t3-t2, 'seconds'
IF (RANK == 0 .AND. I_PRINT_SCREEN == 1) WRITE(2,'(A40,F10.2,A10)')'force assembly time:', t3-t2, 'seconds'
IF (RANK == 0 .AND. I_PRINT_SCREEN == 1) WRITE(*,'(A88)')'---------- ---------- ---------- ---------- ---------- ---------- ---------- ---------- ---------- ----------'

!CALL CPU_TIME(t4)
!IF (RANK == 0) WRITE(*,'(A40,F10.2,A10)')'MPI time:', t4-t3, 'seconds'
!IF (RANK == 0) WRITE(2,'(A40,F10.2,A10)')'MPI time:', t4-t3, 'seconds'

SELECT CASE (NDIM)
   CASE(2)
      CALL ASSEM_PML_PETSC    ( AK_PETSC, AM_PETSC, AC_PETSC,           DIAG_M_PETSC, XYZ, INOD, NGP, MTEL, PMAT, PMat_Lambda, PMat_Mu, PML_PARAM, ID_BC, ID, SIZE, RANK )
!     add temporal damping for stability
!     CALL ASSEM_PML_PETSC_TD ( AK_PETSC, AM_PETSC, AC_PETSC, AG_PETSC, DIAG_M_PETSC, XYZ, INOD, NGP, MTEL, PMAT, PML_PARAM, ID_BC, ID, SIZE, RANK )
!     mixed FEM for elastodynamics (to check if mixed formulation is stable)
!     CALL ASSEM_MIXED_FEM_PETSC ( AK_PETSC, AM_PETSC, AC_PETSC, DIAG_M_PETSC, XYZ, INOD, NGP, MTEL, PMAT, PML_PARAM, ID_BC, ID, SIZE, RANK )
   CASE(3)
      CALL ASSEM_PML_PETSC_3D ( AK_PETSC, AM_PETSC, AC_PETSC, AG_PETSC, DIAG_M_PETSC, XYZ, INOD, NGP, MTEL, PMAT, PMat_Lambda, PMat_Mu, PML_PARAM, ID_BC, ID, SIZE, RANK )
!      CALL ASSEM_PML_PETSC_3D_Sezgin ( AK_PETSC, AM_PETSC, AC_PETSC, AG_PETSC, XYZ, INOD, NGP, MTEL, PMAT, PML_PARAM, ID_BC, ID, SIZE, RANK )
END SELECT


!DEALLOCATE ( XYZ, INOD, NGP, ID_BC, PMat_Lambda, PMat_Mu )


CALL CPU_TIME(t5)
IF (RANK == 0 .AND. I_PRINT_SCREEN == 1) WRITE(*,'(A40,F10.2,A10)')'PML assemble time:', t5-t3, 'seconds'
IF (RANK == 0 .AND. I_PRINT_SCREEN == 1) WRITE(2,'(A40,F10.2,A10)')'PML assemble time:', t5-t3, 'seconds'
IF (RANK == 0 .AND. I_PRINT_SCREEN == 1) WRITE(*,'(A88)')'---------- ---------- ---------- ---------- ---------- ---------- ---------- ---------- ---------- ----------'

CALL MatAssemblyBegin ( AK_PETSC,    MAT_FINAL_ASSEMBLY, IERR )
CALL MatAssemblyBegin ( AM_PETSC,    MAT_FINAL_ASSEMBLY, IERR )
CALL MatAssemblyBegin ( AC_PETSC,    MAT_FINAL_ASSEMBLY, IERR )
CALL MatAssemblyBegin ( AG_PETSC,    MAT_FINAL_ASSEMBLY, IERR )
CALL MatAssemblyBegin ( AK_RD_PETSC, MAT_FINAL_ASSEMBLY, IERR )
CALL MatAssemblyBegin ( AM_RD_PETSC, MAT_FINAL_ASSEMBLY, IERR )

CALL VecAssemblyBegin ( B_PETSC        , IERR )
CALL VecAssemblyBegin ( RMJ_PETSC      , IERR )
CALL VecAssemblyBegin ( DIAG_M_PETSC   , IERR )
CALL VecAssemblyBegin ( DIAG_M_RD_PETSC, IERR )

CALL MatAssemblyEnd   ( AK_PETSC, MAT_FINAL_ASSEMBLY, IERR )
CALL MatAssemblyEnd   ( AM_PETSC, MAT_FINAL_ASSEMBLY, IERR )
CALL MatAssemblyEnd   ( AC_PETSC, MAT_FINAL_ASSEMBLY, IERR )
CALL MatAssemblyEnd   ( AG_PETSC, MAT_FINAL_ASSEMBLY, IERR )
CALL MatAssemblyEnd   ( AK_RD_PETSC, MAT_FINAL_ASSEMBLY, IERR )
CALL MatAssemblyEnd   ( AM_RD_PETSC, MAT_FINAL_ASSEMBLY, IERR )

CALL VecAssemblyEnd   ( B_PETSC        , IERR )
CALL VecAssemblyEnd   ( RMJ_PETSC      , IERR )
CALL VecAssemblyEnd   ( DIAG_M_PETSC   , IERR )
CALL VecAssemblyEnd   ( DIAG_M_RD_PETSC, IERR )


CALL CPU_TIME(t6)
IF (RANK == 0 .AND. I_PRINT_SCREEN == 1) WRITE(*,'(A40,F10.2,A10)')'MPI time:', t6-t5, 'seconds'
IF (RANK == 0 .AND. I_PRINT_SCREEN == 1) WRITE(2,'(A40,F10.2,A10)')'MPI time:', t6-t5, 'seconds'
IF (RANK == 0 .AND. I_PRINT_SCREEN == 1) WRITE(*,'(A88)')'---------- ---------- ---------- ---------- ---------- ---------- ---------- ---------- ---------- ----------'


! Compute the inverse of the global system mass matrix:
!------------------------------------------------------------------------------------------------------------
Call VecCopy ( DIAG_M_PETSC, DIAG_iM_PETSC, IERR )
CALL VecReciprocal ( DIAG_iM_PETSC, IERR )
!------------------------------------------------------------------------------------------------------------


SELECT CASE (NDOF)
   CASE(2, 3, 5)
!     CALL NEWMARK_PETSC_0 ( AK_PETSC, AM_PETSC, AC_PETSC, B_PETSC, AK_RD_PETSC, AM_RD_PETSC, ID, LTRANS, RMJ_PETSC, BACL, RANK, SIZE, XYZ, INOD, NNDH, NNVH, NNAH, EqDis, EqVel, EqAcc )
!     CALL EXPLICIT_RK2_PETSC ( AK_PETSC, DIAG_M_PETSC, AC_PETSC, B_PETSC, AK_RD_PETSC, DIAG_M_RD_PETSC, ID, LTRANS, RMJ_PETSC, BACL, RANK )
!     CALL EXPLICIT_RK4_PETSC ( AK_PETSC, DIAG_M_PETSC, AC_PETSC, B_PETSC, AK_RD_PETSC, DIAG_M_RD_PETSC, ID, LTRANS, RMJ_PETSC, BACL, RANK, NNDH, NNVH, NNAH, EqDis, EqVel, EqAcc )
!     CALL Eigen ( AK_PETSC, AC_PETSC, DIAG_M_PETSC, RANK, SIZE )
!     CALL NEWMARK_PETSC_3 ( AK_PETSC, AM_PETSC, AC_PETSC, AG_PETSC, B_PETSC, AK_RD_PETSC, AM_RD_PETSC, ID, LTRANS, RANK, XYZ, INOD )

   CASE(9)
!     CALL NEWMARK_PETSC_3 ( AK_PETSC, AM_PETSC, AC_PETSC, AG_PETSC, B_PETSC, AK_RD_PETSC, AM_RD_PETSC, ID, LTRANS, RANK, NNDH, NNVH, NNAH, EqDis, EqVel, EqAcc )
!     CALL EXPLICIT_RK2_PETSC_PML_3D ( AK_PETSC, DIAG_M_PETSC, AC_PETSC, AG_PETSC, B_PETSC, AK_RD_PETSC, DIAG_M_RD_PETSC, ID, LTRANS, RANK )
!     CALL EXPLICIT_RK4_PETSC_PML_3D ( AK_PETSC, DIAG_M_PETSC, AC_PETSC, AG_PETSC, B_PETSC, AK_RD_PETSC, DIAG_M_RD_PETSC, ID, LTRANS, RANK, NNDH, NNVH, NNAH, EqDis, EqVel, EqAcc, idx_from, ID_Para, Param )
!     CALL NEWMARK_PETSC_0 ( AK_PETSC, AM_PETSC, AC_PETSC, B_PETSC, AK_RD_PETSC, AM_RD_PETSC, ID, LTRANS, RMJ_PETSC, BACL, RANK, SIZE, XYZ, INOD )
!     CALL EXPLICIT_RK4_PETSC ( AK_PETSC, DIAG_M_PETSC, AC_PETSC, B_PETSC, AK_RD_PETSC, DIAG_M_RD_PETSC, ID, LTRANS, RMJ_PETSC, BACL, RANK, NNDH, NNVH, NNAH, EqDis, EqVel, EqAcc )
!     CALL Eigen_3D ( AK_PETSC, AC_PETSC, AG_PETSC, DIAG_M_PETSC, RANK, SIZE )
!     CALL Eigen ( AK_PETSC, AC_PETSC, DIAG_M_PETSC, RANK, SIZE )
      Call EXPLICIT_RK4_PETSC_PML_3D ( AK_PETSC, DIAG_iM_PETSC, AC_PETSC, AG_PETSC, B_PETSC, AK_RD_PETSC, DIAG_M_RD_PETSC, ID, LTRANS, RANK, &
                                       NNDH, NNVH, NNAH, EqDis, EqVel, EqAcc, idx_from, ID_Para, Param, &
                                       U_Store_Mapping, NStore_Mapping, U_Store_Numbers_Global ) !, k1_m_Store_Mapping, k2_m_Store_Mapping, k3_m_Store_Mapping )

!      Call EXPLICIT_Central_Difference_2 ( AK_PETSC, DIAG_M_PETSC, AC_PETSC, B_PETSC, AK_RD_PETSC, DIAG_M_RD_PETSC, ID, RANK, NNDH, EqDis, U_Store_Mapping, NStore_Mapping, U_Store_Numbers_Global )
!      Call EXPLICIT_Central_Difference_3 ( AK_PETSC, DIAG_M_PETSC, AC_PETSC, AG_PETSC, B_PETSC, AK_RD_PETSC, DIAG_M_RD_PETSC, ID, RANK, NNDH, EqDis, U_Store_Mapping, NStore_Mapping, U_Store_Numbers_Global )

END SELECT

