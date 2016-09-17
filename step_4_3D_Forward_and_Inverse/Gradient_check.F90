! ===================================================================================================================================================
! Example 6 - CPU 12 - material point 12540: ( 0.625, 0.625, 0 )
! ===================================================================================================================================================
Call MPI_Barrier ( PETSC_COMM_WORLD, ierr )
If ( RANK == 12 ) Then
   Write ( *     ,'(A138)' ) '-----------------------------------------------------------------------------------------------------------------------------------------'
!   PMat_Mu (12540) = PMat_Mu (12540) * 1.0000100d0  ;
   Do i = 1, nj
!     dist = ( xyz(i,1) - 0.6250d0 )**2 + ( xyz(i,2) - 0.6250d0 )**2 + xyz(i,3)**2 
     dist = ( xyz(i,1) - 1.250d0 )**2 + ( xyz(i,2) - 1.250d0 )**2 + xyz(i,3)**2 
     If ( dist <= 1e-5 ) then
       Write(*,'(i10,3f10.5,i5)') i, xyz(i,1), xyz(i,2), xyz(i,3), rank
     End If
   End Do
  Write ( *     ,'(A138)' ) '-----------------------------------------------------------------------------------------------------------------------------------------'
End If
Call MPI_Barrier ( PETSC_COMM_WORLD, ierr ) ; stop
! ===================================================================================================================================================
Include 'Solve_Forward.F90'
Call Compute_Misfit ( B_mis_PETSC, NNDH, EqDis, U_Store_Mapping, NStore_Mapping, EqDis_MAP_U_Store_Numbers_Global, Dis_meas, Cost_Misfit )
! -----------------------------------------------------------------------------------------------------------------------------------------
! Example 6 - CPU 12 - material point 12540: ( 0.625, 0.625, 0 )
  If ( Rank == 0 ) Then       ! for 350 steps
    Write ( *     ,'(A138)' ) '-----------------------------------------------------------------------------------------------------------------------------------------'
    write(*,*) 'cost misfit', Cost_Misfit
    Write ( *     ,'(A138)' ) '-----------------------------------------------------------------------------------------------------------------------------------------' 
    write(*,'(a25,e30.16)') 'g_mu_1.1d0     =',  ( 2.925899639436918d-009 - 2.931579309405469d-009 ) / (400.0d0 * 0.1d0 )
    write(*,'(a25,e30.16)') 'g_mu_1.01d0    =',  ( 2.931003520150759d-009 - 2.931579309405469d-009 ) / (400.0d0 * 0.01d0 )
    write(*,'(a25,e30.16)') 'g_mu_1.001d0   =',  ( 2.931521651267220d-009 - 2.931579309405469d-009 ) / (400.0d0 * 0.001d0 ) ! = -0.1441453456222396E-12
    write(*,'(a25,e30.16)') 'g_mu_1.0001d0  =',  ( 2.931573542798511d-009 - 2.931579309405469d-009 ) / (400.0d0 * 0.0001d0 )
    write(*,'(a25,e30.16)') 'g_mu_1.00001d0 =',  ( 2.931578732736829d-009 - 2.931579309405469d-009 ) / (400.0d0 * 0.00001d0 )
    write(*,'(a25,e30.16)') 'OtD g_Mu (12540) =', -0.9961571941299456d-12 ! with m_inv
    write(*,'(a25,e30.16)') 'DtO g_Mu (12540) =',  0.2254769158439468d-10
! OtD - without the m_inv matrix:   g_Mu (12540) =       -0.1441199644285221E-12 
! DtO - only y_t approximation: g_Mu (12540) =       -0.1370887748275357E-12
    Write ( *     ,'(A138)') '-----------------------------------------------------------------------------------------------------------------------------------------'
  End If
! -----------------------------------------------------------------------------------------------------------------------------------------
! ===================================================================================================================================================
! ===================================================================================================================================================


! ===================================================================================================================================================
! Example 6 - CPU 12 - material point 363: ( 1.25, 1.25, 0 )
! ===================================================================================================================================================
Call MPI_Barrier ( PETSC_COMM_WORLD, IERR )
If ( RANK == 12 ) Then
   Write ( *     ,'(A138)' ) '-----------------------------------------------------------------------------------------------------------------------------------------'
   PMat_Mu (363) = PMat_Mu (363) * 1.000000d0  ;
   Do I = 1, NJ
     dist = ( XYZ(I,1) - 1.250d0 )**2 + ( XYZ(I,2) - 1.250d0 )**2 + XYZ(I,3)**2 
     If ( dist <= 1e-5 ) Then
       Write(*,'(I10,3F10.5,I5)') I, XYZ(I,1), XYZ(I,2), XYZ(I,3), RANK
     End If
   End Do
  Write ( *     ,'(A138)' ) '-----------------------------------------------------------------------------------------------------------------------------------------'
End If
Call MPI_Barrier ( PETSC_COMM_WORLD, IERR ) ; !stop
! ===================================================================================================================================================
Include 'Solve_Forward.F90'
Call Compute_Misfit ( B_mis_PETSC, NNDH, EqDis, U_Store_Mapping, NStore_Mapping, EqDis_MAP_U_Store_Numbers_Global, Dis_meas, Cost_Misfit )
! -----------------------------------------------------------------------------------------------------------------------------------------
! Example 6 - CPU 12 - material point 363: ( 1.25, 1.25, 0 )
  If ( Rank == 0 ) Then       ! for 350 steps
    Write ( *     ,'(A138)' ) '-----------------------------------------------------------------------------------------------------------------------------------------'
    write(*,*) 'cost misfit', Cost_Misfit
    Write ( *     ,'(A138)' ) '-----------------------------------------------------------------------------------------------------------------------------------------' 
    write(*,'(a25,e30.16)') 'g_mu_1.1d0     =',  ( 2.932539491353714d-009 - 2.931579309405469d-009 ) / (400.0d0 * 0.1d0 )
    write(*,'(a25,e30.16)') 'g_mu_1.01d0    =',  ( 2.931677098614636d-009 - 2.931579309405469d-009 ) / (400.0d0 * 0.01d0 )
    write(*,'(a25,e30.16)') 'g_mu_1.001d0   =',  ( 2.931589106642909d-009 - 2.931579309405469d-009 ) / (400.0d0 * 0.001d0 ) 
    write(*,'(a25,e30.16)') 'g_mu_1.0001d0  =',  ( 2.931580289313013d-009 - 2.931579309405469d-009 ) / (400.0d0 * 0.0001d0 )
    write(*,'(a25,e30.16)') 'g_mu_1.00001d0 =',  ( 2.931579407398058d-009 - 2.931579309405469d-009 ) / (400.0d0 * 0.00001d0 ) ! 0.2449814725558668E-13
    write(*,'(a25,e30.16)') 'OtD g_Mu (363) =', -0.0 ! with m_inv
    write(*,'(a25,e30.16)') 'DtO g_Mu (363) =',  0.0
! OtD - without the m_inv matrix:   g_Mu (363) =        0.2448795454141181E-13
! DtO - only y_t approximation: g_Mu (12540) =       -0.1370887748275357E-12
    Write ( *     ,'(A138)') '-----------------------------------------------------------------------------------------------------------------------------------------'
  End If
! -----------------------------------------------------------------------------------------------------------------------------------------
! ===================================================================================================================================================
! ===================================================================================================================================================

! Call MPI_Barrier ( PETSC_COMM_WORLD, ierr ) ; Write(*,*) 'end of misfit computation' ; Stop










! ===================================================================================================================================================
! gradient check (misfit)
! ===================================================================================================================================================


PMat_Lambda = 500.0d0
PMat_Mu     = 400.0d0


! ===================================================================================================================================================
! Example 6 - CPU 12 - material point 363: ( 1.25, 1.25, 0 )
! ===================================================================================================================================================
Call MPI_Barrier ( PETSC_COMM_WORLD, IERR )
If ( RANK == 12 ) Then
   Write ( *     ,'(A138)' ) '-----------------------------------------------------------------------------------------------------------------------------------------'
   PMat_Lambda (363) = PMat_Lambda (363) * 1.0000000000d0  ;
   Do I = 1, NJ
     dist = ( XYZ(I,1) - 1.250d0 )**2 + ( XYZ(I,2) - 1.250d0 )**2 + XYZ(I,3)**2 
     If ( dist <= 1e-5 ) Then
       Write(*,'(I10,3F10.5,I5)') I, XYZ(I,1), XYZ(I,2), XYZ(I,3), RANK
     End If
   End Do
  Write ( *     ,'(A138)' ) '-----------------------------------------------------------------------------------------------------------------------------------------'
End If
Call MPI_Barrier ( PETSC_COMM_WORLD, IERR ) ; !stop
! ===================================================================================================================================================
Include 'Solve_Forward.F90'
Call Compute_Misfit ( B_mis_PETSC, NNDH, EqDis, U_Store_Mapping, NStore_Mapping, EqDis_MAP_U_Store_Numbers_Global, Dis_meas, Cost_Misfit )
! -----------------------------------------------------------------------------------------------------------------------------------------
! Example 6 - CPU 12 - material point 363: ( 1.25, 1.25, 0 )
  If ( Rank == 0 ) Then       ! for 350 steps
    Write ( *     ,'(A138)' ) '-----------------------------------------------------------------------------------------------------------------------------------------'
    write(*,*) 'cost misfit', Cost_Misfit
    Write ( *     ,'(A138)' ) '-----------------------------------------------------------------------------------------------------------------------------------------' 
    write(*,'(a25,e30.16)') 'g_L_1.1d0     =',  ( 1.170668756822149d-010 - 1.170377194664848d-010 ) / (400.0d0 * 0.1d0 )
    write(*,'(a25,e30.16)') 'g_L_1.01d0    =',  ( 1.170406899876179d-010 - 1.170377194664848d-010 ) / (400.0d0 * 0.01d0 )
    write(*,'(a25,e30.16)') 'g_L_1.001d0   =',  ( 1.170380170792436d-010 - 1.170377194664848d-010 ) / (400.0d0 * 0.001d0 ) 
    write(*,'(a25,e30.16)') 'g_L_1.0001d0  =',  ( 1.170377492333794d-010 - 1.170377194664848d-010 ) / (400.0d0 * 0.0001d0 )
    write(*,'(a25,e30.16)') 'g_L_1.00001d0 =',  ( 1.170377224432302d-010 - 1.170377194664848d-010 ) / (400.0d0 * 0.00001d0 ) ! 0.7441863505213146E-15
    write(*,'(a25,e30.16)') 'OtD g_Mu (363) =', -0.0 ! with m_inv
    write(*,'(a25,e30.16)') 'DtO g_Mu (363) =',  0.0
! OtD - without the m_inv matrix:   g_Lambda (363) =        0.7442647494744367E-15
! DtO - only y_t approximation: g_Mu (12540) =       -0.1370887748275357E-12
    Write ( *     ,'(A138)') '-----------------------------------------------------------------------------------------------------------------------------------------'
  End If
! -----------------------------------------------------------------------------------------------------------------------------------------
! ===================================================================================================================================================
! ===================================================================================================================================================

! Call MPI_Barrier ( PETSC_COMM_WORLD, ierr ) ; Write(*,*) 'end of misfit computation' ; Stop


! ---------------------------------------------------------------------------------------------------------------------------------------------------
! compute gradient using the adjoint method

! 1- initialize the gradient vector
     Call VecSet ( g_Lambda_PETSC, 0.0d0, IERR )
     Call VecSet ( g_Mu_PETSC,     0.0d0, IERR )

! 2- Solve the Adjoint & Control Problems (misfit)

! ===============================================
! D-t-O
     Call EXPLICIT_RK4_PETSC_PML_3D_Adjoint_Control ( AK_PETSC, DIAG_M_PETSC, DIAG_iM_PETSC, AC_PETSC, AG_PETSC, B_mis_PETSC, ID, RANK,  NNDH, EqDis, U_Store_Mapping, NStore_Mapping, U_Store_Numbers_Global, EqDis_MAP_U_Store_Numbers_Global, Dis_meas, k1_m_Store_Mapping, k2_m_Store_Mapping, k3_m_Store_Mapping, Ux_big_PETSc, Ux_Bigger_Rank_PETSc, vscat_u_hist, P_Bigger_Rank_PETSc, Uk_big_PETSc, Uk_Bigger_Rank_PETSc, g_Lambda_PETSC, g_Mu_PETSC, idx_Mat_to, NJ_Rank_IParts, XYZ, INOD, NGP, MTEL, y_t_PETSC,  y_m_PETSC,  y_b_PETSC, p1_t_PETSC, p1_m_PETSC, p1_b_PETSC, p2_t_PETSC, p2_m_PETSC, p2_b_PETSC, p3_t_PETSC, p3_m_PETSC, p3_b_PETSC, p4_t_PETSC, p4_m_PETSC, p4_b_PETSC, help_PETSC, P_Store_Mapping )
! ===============================================
! O-t-D
!     Call EXPLICIT_RK4_PETSC_PML_3D_Adjoint_OtD ( AK_PETSC, DIAG_M_PETSC, DIAG_iM_PETSC, AC_PETSC, AG_PETSC, B_mis_PETSC, AK_RD_PETSC, DIAG_M_RD_PETSC, ID, NNDH, EqDis, U_Store_Mapping, P_Store_Mapping, NStore_Mapping, U_Store_Numbers_Global, EqDis_MAP_U_Store_Numbers_Global, Dis_meas )
!     Call EXPLICIT_RK4_PETSC_PML_3D_Control_OtD ( U_Store_Mapping, P_Store_Mapping, NStore_Mapping, U_Store_Numbers_Global, Ux_big_PETSc, Ux_Bigger_Rank_PETSc, Uk_big_PETSc, P_Bigger_Rank_PETSc, vscat_u_hist, g_Lambda_PETSC, g_Mu_PETSC, idx_Mat_to, NJ_Rank_IParts, XYZ, INOD, NGP, MTEL, ID, RANK )
! ===============================================

! Call MPI_Barrier ( PETSC_COMM_WORLD, ierr ) ; write(*,*)'adjoint gradient end' ; stop


! 3- assemble the distributed gradient vectors (only the part associated with the misfit; i.e. within the regular domain)
      Call VecAssemblyBegin ( g_Lambda_PETSC, IERR )  ;  Call VecAssemblyEnd   ( g_Lambda_PETSC, IERR )
      Call VecAssemblyBegin (     g_Mu_PETSC, IERR )  ;  Call VecAssemblyEnd   (     g_Mu_PETSC, IERR )

      Call VecScale ( g_Lambda_PETSC, -1.0d0, IERR )
      Call VecScale (     g_Mu_PETSC, -1.0d0, IERR )

! 6- apply the (inverse) mass matrix associated with the gradient vector
!      Call VecPointwiseMult ( g_Lambda_PETSC, g_Lambda_PETSC, iMg_PETSC, IERR )
!      Call VecPointwiseMult (     g_Mu_PETSC,     g_Mu_PETSC, iMg_PETSC, IERR )


! ---------------- print out ----------------------------------------------------------------------

      CALL VecCreateSeq ( PETSC_COMM_SELF, NJ, g_Lambda_Rank_PETSC, IERR ) 
      Call VecSet ( g_Lambda_Rank_PETSC, 0.0d0, IERR )

      Call VecScatterBegin ( vscat_Mat, g_Lambda_PETSC, g_Lambda_Rank_PETSC, INSERT_VALUES, SCATTER_FORWARD, IERR )
      Call VecScatterEnd   ( vscat_Mat, g_Lambda_PETSC, g_Lambda_Rank_PETSC, INSERT_VALUES, SCATTER_FORWARD, IERR )


     Call VecGetArrayF90 ( g_Lambda_Rank_PETSC, g_Lambda_Rank, IERR )


IF ( RANK == 12 ) then
   Write ( *     ,'(A138)' ) '-----------------------------------------------------------------------------------------------------------------------------------------'
!   write(*,'(a25,e30.16)') 'g_Mu (1833) =', g_Mu_Rank (1833)
!   write(*,'(a25,e30.16)') 'g_Mu (11373) =', g_Mu_Rank (11373)
   !write(*,'(a25,e30.16)') 'g_Mu (12540) =', g_Mu_Rank (12540)
   write(*,'(a25,e30.16)') 'g_Lambda (363) =', g_Lambda_Rank (363)
   Write ( *     ,'(A138)' ) '-----------------------------------------------------------------------------------------------------------------------------------------'
End If


!... end transfer and engineering manipulation
     Call VecRestoreArrayF90 ( g_Lambda_Rank_PETSC, g_Lambda_Rank, IERR )

Call MPI_Barrier ( PETSC_COMM_WORLD, ierr ) ; write(*,*)'officially end' ; stop













! check regularization gradient

pmat_lambda = 500.0d0
pmat_mu     = 400.0d0




! ===================================================================================================================================================
! Example 6 - CPU 12 - material point 363: ( 1.25, 1.25, 0 )
! ===================================================================================================================================================
Call MPI_Barrier ( PETSC_COMM_WORLD, IERR )
If ( RANK == 12 ) Then
   Write ( *     ,'(A138)' ) '-----------------------------------------------------------------------------------------------------------------------------------------'
!   PMat_Mu (363) = PMat_Mu (363) * 1.1000000000d0  ;
   PMat_lambda (363) = PMat_lambda (363) * 1.100000d0  ;
   Do I = 1, NJ
     dist = ( XYZ(I,1) - 1.250d0 )**2 + ( XYZ(I,2) - 1.250d0 )**2 + XYZ(I,3)**2 
     If ( dist <= 1e-5 ) Then
       Write(*,'(I10,3F10.5,I5)') I, XYZ(I,1), XYZ(I,2), XYZ(I,3), RANK
     End If
   End Do
  Write ( *     ,'(A138)' ) '-----------------------------------------------------------------------------------------------------------------------------------------'
End If
Call MPI_Barrier ( PETSC_COMM_WORLD, IERR ) ; !stop
! ===================================================================================================================================================

! 2.1- transfer to sequential PETSc vectors
  Call VecGetArrayF90 ( Lambda_Rank_PETSC, Lambda_Rank, IERR )
!... transfer ...
  ForAll ( I = 1:NJ ) Lambda_Rank (I) = PMat_Lambda (I)
  Call VecRestoreArrayF90 ( Lambda_Rank_PETSC, Lambda_Rank, IERR )

  Call VecGetArrayF90 ( Mu_Rank_PETSC, Mu_Rank, IERR )
!... transfer ...
  ForAll ( I = 1:NJ ) Mu_Rank (I) = PMat_Mu (I)
  Call VecRestoreArrayF90 ( Mu_Rank_PETSC, Mu_Rank, IERR )

! 2.2- scatter from Lambda_Rank_PETSC to Lambda_PETSC
  Call VecScatterBegin ( vscat_Mat, Lambda_Rank_PETSC, Lambda_PETSC, INSERT_VALUES, SCATTER_REVERSE_LOCAL, IERR )
  Call VecScatterEnd   ( vscat_Mat, Lambda_Rank_PETSC, Lambda_PETSC, INSERT_VALUES, SCATTER_REVERSE_LOCAL, IERR )
  Call VecScatterBegin ( vscat_Mat,     Mu_Rank_PETSC,     Mu_PETSC, INSERT_VALUES, SCATTER_REVERSE_LOCAL, IERR )
  Call VecScatterEnd   ( vscat_Mat,     Mu_Rank_PETSC,     Mu_PETSC, INSERT_VALUES, SCATTER_REVERSE_LOCAL, IERR )

! ---------------------------------------------------------------------------------------------------------------------------------------------------
! VIII- Extend the material properties from the interface into the PML domain -----------------------------------------------------------------------
! ---------------------------------------------------------------------------------------------------------------------------------------------------
     Call VecScatterBegin ( vscat_mat_extend, Lambda_PETSC, Lambda_PETSC, INSERT_VALUES, SCATTER_FORWARD, IERR )
     Call VecScatterEnd   ( vscat_mat_extend, Lambda_PETSC, Lambda_PETSC, INSERT_VALUES, SCATTER_FORWARD, IERR )
     Call VecScatterBegin ( vscat_mat_extend,     Mu_PETSC,     Mu_PETSC, INSERT_VALUES, SCATTER_FORWARD, IERR )
     Call VecScatterEnd   ( vscat_mat_extend,     Mu_PETSC,     Mu_PETSC, INSERT_VALUES, SCATTER_FORWARD, IERR )


! TN
      Call Assem_TV_TN_Regularization ( Reg_Lambda_PETSC, Reg_Mu_PETSC, PMat_Lambda, PMat_Mu, idx_Mat_to, NJ_Rank_IParts, XYZ, INOD, NGP, MTEL, ID, cost_reg_Lambda_rank, cost_reg_Mu_rank )
      Call VecAssemblyBegin ( Reg_Lambda_PETSC, IERR )  ;  Call VecAssemblyEnd   ( Reg_Lambda_PETSC, IERR )
      Call VecAssemblyBegin (     Reg_Mu_PETSC, IERR )  ;  Call VecAssemblyEnd   (     Reg_Mu_PETSC, IERR )

! 2- compute the cost of the regularization component of the discrete objective functional
      Call MPI_Reduce ( cost_reg_Lambda_rank, Cost_Regularization_Lambda, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, PETSC_COMM_WORLD, IERR)
      Call MPI_Reduce ( cost_reg_Mu_rank,     Cost_Regularization_Mu,     1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, PETSC_COMM_WORLD, IERR)

!      Call VecDot ( Reg_Lambda_PETSC, Lambda_PETSC, Cost_Regularization_Lambda_2, IERR )
!      Call VecDot (     Reg_Mu_PETSC,     Mu_PETSC, Cost_Regularization_Mu_2,     IERR )

! -----------------------------------------------------------------------------------------------------------------------------------------
! Example 6 - CPU 12 - material point 363: ( 1.25, 1.25, 0 )
  If ( Rank == 0 ) Then       ! for 350 steps
    Write ( *     ,'(A138)' ) '-----------------------------------------------------------------------------------------------------------------------------------------'
    write(*,*) 'cost reg', Cost_Regularization_Mu
    write(*,*) 'cost reg', Cost_Regularization_lambda
    Write ( *     ,'(A138)' ) '-----------------------------------------------------------------------------------------------------------------------------------------' 
!    write(*,'(a25,e30.16)') 'g_mu_1.1d0     =',  ( 2.932539491353714d-009 - 23.3750905571258d0 ) / (400.0d0 * 0.1d0 )
    write(*,'(a25,e30.16)') 'g_mu_1.01d0    =',  ( 25.7119836697561d0 - 23.3750905571258d0 ) / (400.0d0 * 0.01d0 )
    write(*,'(a25,e30.16)') 'g_mu_1.001d0   =',  ( 23.5989107505086d0 - 23.3750905571258d0 ) / (400.0d0 * 0.001d0 ) 
    write(*,'(a25,e30.16)') 'g_mu_1.0001d0  =',  ( 23.3975651600810d0 - 23.3750905571258d0 ) / (400.0d0 * 0.0001d0 )
    write(*,'(a25,e30.16)') 'g_mu_1.00001d0 =',  ( 23.3738408031019d0 - 23.3750905571258d0 ) / (400.0d0 * 0.00001d0 ) ! 0.2449814725558668E-13
    write(*,'(a25,e30.16)') 'g_mu_1.00001d0 =',  ( 23.3748568102750d0 - 23.3750905571258d0 ) / (400.0d0 * 0.00001d0 ) ! 0.2449814725558668E-13
    write(*,'(a25,e30.16)') 'OtD g_Mu (363) =', -0.0 ! with m_inv
    write(*,'(a25,e30.16)') 'DtO g_Mu (363) =',  0.0
! OtD - without the m_inv matrix:   g_Mu (363) =        0.2448795454141181E-13
! DtO - only y_t approximation: g_Mu (12540) =       -0.1370887748275357E-12
    Write ( *     ,'(A138)') '-----------------------------------------------------------------------------------------------------------------------------------------'
  End If
! -----------------------------------------------------------------------------------------------------------------------------------------
! ===================================================================================================================================================
! ===================================================================================================================================================


! ---------------- print out ----------------------------------------------------------------------

      CALL VecCreateSeq ( PETSC_COMM_SELF, NJ, g_Lambda_Rank_PETSC, IERR ) 
      Call VecSet ( g_Lambda_Rank_PETSC, 0.0d0, IERR )

      Call VecScatterBegin ( vscat_Mat, Reg_lambda_PETSC, g_Lambda_Rank_PETSC, INSERT_VALUES, SCATTER_FORWARD, IERR )
      Call VecScatterEnd   ( vscat_Mat, Reg_lambda_PETSC, g_Lambda_Rank_PETSC, INSERT_VALUES, SCATTER_FORWARD, IERR )


     Call VecGetArrayF90 ( g_Lambda_Rank_PETSC, g_Lambda_Rank, IERR )


IF ( RANK == 12 ) then
   Write ( *     ,'(A138)' ) '-----------------------------------------------------------------------------------------------------------------------------------------'
!   write(*,'(a25,e30.16)') 'g_Mu (1833) =', g_Mu_Rank (1833)
!   write(*,'(a25,e30.16)') 'g_Mu (11373) =', g_Mu_Rank (11373)
   !write(*,'(a25,e30.16)') 'g_Mu (12540) =', g_Mu_Rank (12540)
   write(*,'(a25,e30.16)') 'g_Lambda (363) =', g_Lambda_Rank (363)
   Write ( *     ,'(A138)' ) '-----------------------------------------------------------------------------------------------------------------------------------------'
End If


!... end transfer and engineering manipulation
     Call VecRestoreArrayF90 ( g_Lambda_Rank_PETSC, g_Lambda_Rank, IERR )

Call MPI_Barrier ( PETSC_COMM_WORLD, ierr ) ; write(*,*)'officially end' ; stop







 Call MPI_Barrier ( PETSc_COMM_WORLD, IERR ) ; Stop !===















