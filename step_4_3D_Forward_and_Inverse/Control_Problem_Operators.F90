
! Control problem operators: Regularization matrix and inverse of the gradient mass vector
Call Assem_R_Mg_3D ( Reg_PETSC, iMg_PETSC, XYZ, INOD, NGP, MTEL, ID, idx_Mat_to, NJ_Rank_IParts )

CalL MatAssemblyBegin ( Reg_PETSC, MAT_FINAL_ASSEMBLY, IERR )
Call MatAssemblyEnd   ( Reg_PETSC, MAT_FINAL_ASSEMBLY, IERR )
Call VecAssemblyBegin ( iMg_PETSC, IERR )
Call VecAssemblyEnd   ( iMg_PETSC, IERR )


! Invert the gradient mass vector. This matrix contains zeros on the diagonal inside the PML; hence, we have to do this.
  Call VecGetArrayF90 ( iMg_PETSC, iMg_Rank, IERR )
! ----- -----
     Do I = 1, NJ_Mapping
        If ( iMg_Rank (I) > 1.0d-10 ) Then
           iMg_Rank (I) = 1.0d0 / iMg_Rank (I)
           restriction_RD ( I ) = 1.0d0
        Else
           iMg_Rank (I) = 0.0d0
           restriction_RD ( I ) = 0.0d0
        End If
     End Do
! ----- -----
  Call VecRestoreArrayF90 ( iMg_PETSC, iMg_Rank, IERR )


