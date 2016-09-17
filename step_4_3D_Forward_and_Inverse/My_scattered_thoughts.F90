! ===================================================================================================================================================
! gradient check (misfit)
! ===================================================================================================================================================


! for 8 processors, (rank 1) I = 190 is the XYZ = 0.
!                   (rank 7) I = 187.
! rank 7: I = 1833. X = 1, Y = 1, Z = 0.



! ===================================================================
Call MPI_Barrier ( PETSC_COMM_WORLD, ierr )
IF ( RANK == 7 ) then
   Write ( *     ,'(A138)' ) '-----------------------------------------------------------------------------------------------------------------------------------------'
   PMat_Mu (1833) = PMat_Mu (1833) !* 1.0001d0

   write(*,'(a20,f10.5)') 'PMat_Mu (1833) =', PMat_Mu (1833)

   do i = 1, nj
     dist = ( xyz(i,1) - 1.0d0 )**2 + ( xyz(i,2) - 1.0d0 )**2 + xyz(i,3)**2
     if ( dist <= 1e-3 ) then
       write(*,'(i5,3f10.5)') i, xyz(i,1), xyz(i,2), xyz(i,3)
     end if
   end do

End If
Call MPI_Barrier ( PETSC_COMM_WORLD, ierr )
! ===================================================================

Include 'Solve_Forward.F90'
call Compute_Misfit ( DIAG_M_PETSC, B_mis_PETSC, NNDH, EqDis, U_Store_Mapping, NStore_Mapping, EqDis_MAP_U_Store_Numbers_Global, Dis_meas, Cost_Misfit )


  If ( Rank == 0 ) Then
    Write ( *     ,'(A138)' ) '-----------------------------------------------------------------------------------------------------------------------------------------'
    write(*,*) 'cost misfit', Cost_Misfit
    Write ( *     ,'(A138)' ) '-----------------------------------------------------------------------------------------------------------------------------------------'
    write(*,'(a25,e30.16)') 'g_mu_1.1d0     =',  (8.714894578478584d-003 - 8.716161292096132d-003) / (5.0d0 * 0.1d0 )
    write(*,'(a25,e30.16)') 'g_mu_1.01d0    =',  (8.716036544278719d-003 - 8.716161292096132d-003) / (5.0d0 * 0.01d0 )
    write(*,'(a25,e30.16)') 'g_mu_1.001d0   =',  (8.716148837682876d-003 - 8.716161292096132d-003) / (5.0d0 * 0.001d0 )
    write(*,'(a25,e30.16)') 'g_mu_1.0001d0  =',  (8.716160046859663d-003 - 8.716161292096132d-003) / (5.0d0 * 0.0001d0 )
    write(*,'(a25,e30.16)') 'g_mu_1.00001d0 =',  (8.716161167574539d-003 - 8.716161292096132d-003) / (5.0d0 * 0.00001d0 )
    Write ( *     ,'(A138)') '-----------------------------------------------------------------------------------------------------------------------------------------'
  End If

! -----------------------------------------------------------------------------------------------------------------------------------------
!         cost misfit    =        8.716161292096134E-003
!         g_mu_1.1d0     =       -0.2533427235097002E-05
!         g_mu_1.01d0    =       -0.2494956348256250E-05
!         g_mu_1.001d0   =       -0.2490882651182136E-05
!         g_mu_1.0001d0  =       -0.2490472938987409E-05
!         g_mu_1.00001d0 =       -0.2490431853796604E-05
!
!            g_Mu (1833) =        0.6259354340137708E-03
! -----------------------------------------------------------------------------------------------------------------------------------------


Call MPI_Barrier ( PETSC_COMM_WORLD, ierr )


! ---------------------------------------------------------------------------------------------------------------------------------------------------
! compute gradient using the adjoint method

! 1- initialize the gradient vector
     Call VecSet ( g_Lambda_PETSC, 0.0d0, IERR )
     Call VecSet ( g_Mu_PETSC,     0.0d0, IERR )

! 2- Solve the Adjoint & Control Problems (misfit)
     Call EXPLICIT_RK4_PETSC_PML_3D_Adjoint_Control ( AK_PETSC, DIAG_M_PETSC, AC_PETSC, AG_PETSC, B_mis_PETSC, AK_RD_PETSC, DIAG_M_RD_PETSC, ID, RANK, &
                                                      NNDH, EqDis, U_Store_Mapping, NStore_Mapping, U_Store_Numbers_Global, EqDis_MAP_U_Store_Numbers_Global, Dis_meas, &
                                                      k1_m_Store_Mapping, k2_m_Store_Mapping, k3_m_Store_Mapping, &
                                                      Ux_big_PETSc, Ux_Bigger_Rank_PETSc, vscat_u_hist, P_Bigger_Rank_PETSc, Uk_big_PETSc, Uk_Bigger_Rank_PETSc, &
                                                      g_Lambda_PETSC, g_Mu_PETSC, idx_Mat_to, NJ_Rank_IParts, &
                                                      XYZ, INOD, NGP, MTEL )

! 3- assemble the distributed gradient vectors (only the part associated with the misfit; i.e. within the regular domain)
      Call VecAssemblyBegin ( g_Lambda_PETSC, IERR )  ;  Call VecAssemblyEnd   ( g_Lambda_PETSC, IERR )
      Call VecAssemblyBegin (     g_Mu_PETSC, IERR )  ;  Call VecAssemblyEnd   (     g_Mu_PETSC, IERR )

      Call VecScale ( g_Lambda_PETSC, -1.0d0, IERR )
      Call VecScale (     g_Mu_PETSC, -1.0d0, IERR )

! 6- apply the (inverse) mass matrix associated with the gradient vector
      Call VecPointwiseMult ( g_Lambda_PETSC, g_Lambda_PETSC, iMg_PETSC, IERR )
      Call VecPointwiseMult (     g_Mu_PETSC,     g_Mu_PETSC, iMg_PETSC, IERR )


! ---------------- print out ----------------------------------------------------------------------

      CALL VecCreateSeq ( PETSC_COMM_SELF, NJ, g_Mu_Rank_PETSC, IERR ) 
      Call VecSet ( g_Mu_Rank_PETSC, 0.0d0, IERR )

      Call VecScatterBegin ( vscat_Mat, g_Mu_PETSC, g_Mu_Rank_PETSC, INSERT_VALUES, SCATTER_FORWARD, IERR )
      Call VecScatterEnd   ( vscat_Mat, g_Mu_PETSC, g_Mu_Rank_PETSC, INSERT_VALUES, SCATTER_FORWARD, IERR )


     Call VecGetArrayF90 ( g_Mu_Rank_PETSC, g_Mu_Rank, IERR )
!... transfer to appropriate Fortran matrices & ascertain Lambda is within pre-specified limits


IF ( RANK == 7 ) then
   Write ( *     ,'(A138)' ) '-----------------------------------------------------------------------------------------------------------------------------------------'
   write(*,'(a25,e30.16)') 'g_Mu (1833) =', g_Mu_Rank (1833)
   Write ( *     ,'(A138)' ) '-----------------------------------------------------------------------------------------------------------------------------------------'
End If


!... end transfer and engineering manipulation
     Call VecRestoreArrayF90 ( g_Mu_Rank_PETSC, g_Mu_Rank, IERR )







! ===================================================================================================================================================
! Visualize the initial guess without solving the inverse problem
! ===================================================================================================================================================

! 1- material properties in distributed PETSc vectors

! 1.1-
  Call VecGetArrayF90 ( Lambda_Rank_PETSC, Lambda_Rank, IERR )
!... transfer ...
  ForAll ( I = 1:NJ ) Lambda_Rank (I) = PMat_Lambda (I)
  Call VecRestoreArrayF90 ( Lambda_Rank_PETSC, Lambda_Rank, IERR )

  Call VecGetArrayF90         ( Mu_Rank_PETSC,     Mu_Rank, IERR )
!... transfer ...
  ForAll ( I = 1:NJ ) Mu_Rank (I) = PMat_Mu (I)
  Call VecRestoreArrayF90     ( Mu_Rank_PETSC,     Mu_Rank, IERR )

! 1.2- scatter from {}_Rank_PETSC to {}_PETSC
  Call VecScatterBegin ( vscat_Mat, Lambda_Rank_PETSC, Lambda_PETSC, INSERT_VALUES, SCATTER_REVERSE_LOCAL, IERR )
  Call VecScatterEnd   ( vscat_Mat, Lambda_Rank_PETSC, Lambda_PETSC, INSERT_VALUES, SCATTER_REVERSE_LOCAL, IERR )
  Call VecScatterBegin ( vscat_Mat,     Mu_Rank_PETSC,     Mu_PETSC, INSERT_VALUES, SCATTER_REVERSE_LOCAL, IERR )
  Call VecScatterEnd   ( vscat_Mat,     Mu_Rank_PETSC,     Mu_PETSC, INSERT_VALUES, SCATTER_REVERSE_LOCAL, IERR )

! 1.3- Extend the material properties from the interface into the PML domain -------------------------------------------------------------------------
  Call VecScatterBegin ( vscat_mat_extend, Lambda_PETSC, Lambda_PETSC, INSERT_VALUES, SCATTER_FORWARD, IERR )
  Call VecScatterEnd   ( vscat_mat_extend, Lambda_PETSC, Lambda_PETSC, INSERT_VALUES, SCATTER_FORWARD, IERR )
  Call VecScatterBegin ( vscat_mat_extend,     Mu_PETSC,     Mu_PETSC, INSERT_VALUES, SCATTER_FORWARD, IERR )
  Call VecScatterEnd   ( vscat_mat_extend,     Mu_PETSC,     Mu_PETSC, INSERT_VALUES, SCATTER_FORWARD, IERR )

!  Call VecScatterView ( vscat_mat_extend, PETSC_VIEWER_STDOUT, IERR )

! 1.4- monitor
  Call VecScatterBegin ( vscat_Mat, Lambda_PETSC, Lambda_Rank_PETSC, INSERT_VALUES, SCATTER_FORWARD, IERR )
  Call VecScatterEnd   ( vscat_Mat, Lambda_PETSC, Lambda_Rank_PETSC, INSERT_VALUES, SCATTER_FORWARD, IERR )
  Call VecScatterBegin ( vscat_Mat,     Mu_PETSC,     Mu_Rank_PETSC, INSERT_VALUES, SCATTER_FORWARD, IERR )
  Call VecScatterEnd   ( vscat_Mat,     Mu_PETSC,     Mu_Rank_PETSC, INSERT_VALUES, SCATTER_FORWARD, IERR )

  Call VecGetArrayF90 ( Lambda_Rank_PETSC, Lambda_Rank, IERR )
  Call VecGetArrayF90     ( Mu_Rank_PETSC,     Mu_Rank, IERR )

icounter = 0
if ( rank == 0 ) then
  do ij = 1, nj
    res_val = Lambda_Rank (ij) - Mu_Rank (ij)
    if ( dabs(res_val) > 1.0d-3 ) then
       icounter = icounter +1
       write(*,'(2i5,2f5.0)') icounter, ij, Lambda_Rank (ij), Mu_Rank (ij)
    end if
  end do
end if

  Call VecRestoreArrayF90 ( Lambda_Rank_PETSC, Lambda_Rank, IERR )
  Call VecRestoreArrayF90     ( Mu_Rank_PETSC,     Mu_Rank, IERR )

! - Output Lambda_PETSC, Mu_PETSC, Iter for Paraview ------------------------------------------------------------------------------------------------
Iter = 0
Include 'Material_Visualization.F90'                                                                ! R E P E A T E D
! - --------------------------------*----------------------------------------------------------------------------------------------------------------

! - Call only once at the end of all material updates to shut down top file wrapper -----------------------------------------------------------------
Include 'Material_Visualization_End.F90'                                                            ! O N L Y   O N C E
! - --------------------------------*----------------------------------------------------------------------------------------------------------------

Call MPI_Barrier( PETSC_COMM_WORLD, IERR )
write(*,*) 'wtf'
Stop
! ===================================================================================================================================================

!  Call VecMin ( Lambda_PETSC, Petsc_NULL_INTEGER, val_min_Lambda, IERR )
!  Call VecMin (     Mu_PETSC, Petsc_NULL_INTEGER, val_min_Mu,     IERR )

!if (rank == 0) then
!  write(*,*) 'val_min_Lambda =', val_min_Lambda, 'val_min_Mu =', val_min_Mu
!end if


!Call VecView ( g_Mu_PETSc, PETSc_VIEWER_STDOUT_World, ierr ) 

! ===================================================================================================================================================
! Testing Zone
! ===================================================================================================================================================

goto 10
! 2- get it in Fortran format for system matrix assembly
  Call VecGetArrayF90         ( Lambda_Rank_PETSC,     Lambda_Rank, IERR )
!... transfer to appropriate Fortran matrices ...
  ForAll ( I = 1:NJ ) Lambda_Rank (I) = PMat_Lambda (I)
  Call VecRestoreArrayF90     ( Lambda_Rank_PETSC,     Lambda_Rank, IERR )

  Call VecScatterBegin ( vscat_Mat,     Lambda_Rank_PETSC,     Lambda_PETSC, INSERT_VALUES, SCATTER_REVERSE_LOCAL, IERR )
  Call VecScatterEnd   ( vscat_Mat,     Lambda_Rank_PETSC,     Lambda_PETSC, INSERT_VALUES, SCATTER_REVERSE_LOCAL, IERR )

CALL VecSetValue ( Lambda_PETSC, 9,   2.01d0, INSERT_VALUES, IERR )
Call VecAssemblyBegin ( Lambda_PETSC, IERR )
Call VecAssemblyEnd   ( Lambda_PETSC, IERR )

  Call VecScatterBegin ( vscat_Mat,     Lambda_PETSC,     Lambda_Rank_PETSC, INSERT_VALUES, SCATTER_FORWARD, IERR )
  Call VecScatterEnd   ( vscat_Mat,     Lambda_PETSC,     Lambda_Rank_PETSC, INSERT_VALUES, SCATTER_FORWARD, IERR )

  Call VecGetArrayF90         ( Lambda_Rank_PETSC,     Lambda_Rank, IERR )
!... transfer to appropriate Fortran matrices ...
  ForAll ( I = 1:NJ ) PMat_Lambda (I) = Lambda_Rank (I)
  Call VecRestoreArrayF90     ( Lambda_Rank_PETSC,     Lambda_Rank, IERR )


!write(*,*) pmat_mu
!Call VecView ( Mu_PETSc, PETSc_VIEWER_STDOUT_World, ierr )

!write(*,*) pmat_mu
!stop
10 continue
! ===================================================================================================================================================
! End of Testing Zone
! ===================================================================================================================================================


Call MPI_Barrier ( PETSC_COMM_WORLD, ierr )

FD_gradient = (2.951925090849184d0 - 2.951637778506873d0) / 0.01d0
!adj_grad09 = -0.0500032

!write(*,'(a20,f20.10)') 'FD_gradient', FD_gradient


!write(*,*) PMat_Mu


! Kids playground -----------------------------------------------------------------------------------------------------------------------------------
!  *                                                                                                        *
! ---------------------------------------------------------------------------------------------------------------------------------------------------
Call MPI_Barrier ( PETSC_COMM_WORLD, ierr )
IF ( RANK == 0 ) then

do i = 1, nj

dist = xyz(i,1)**2 + xyz(i,2)**2 + xyz(i,3)**2

  if ( dist <= 1e-3 ) then
!    write(*,*) i, xyz(i,1), xyz(i,2), xyz(i,3)
  end if

end do



End If
Call MPI_Barrier ( PETSC_COMM_WORLD, ierr )
!stop
! Kids playground -----------------------------------------------------------------------------------------------------------------------------------
!  *                                                                                                        *
! ---------------------------------------------------------------------------------------------------------------------------------------------------




!write(*,*) pmat_lambda(48), pmat_mu(48)




  Include 'Solve_Forward.F90'

stop


! Kids playground -----------------------------------------------------------------------------------------------------------------------------------
!  *                                                                                                        *
! ---------------------------------------------------------------------------------------------------------------------------------------------------

! for 8 processors, (rank 1) I = 190 is the XYZ = 0.
!                   (rank 7) I = 187.
! rank 7: I = 1833. X = 1, Y = 1, Z = 0.

! ===================================================================
Call MPI_Barrier ( PETSC_COMM_WORLD, ierr )
IF ( RANK == 7 ) then
   Write ( *     ,'(A138)' ) '-----------------------------------------------------------------------------------------------------------------------------------------'
   PMat_Mu (1833) = PMat_Mu (1833) * 1.0100d0

   write(*,'(a20,f10.5)') 'PMat_Mu (1833) =', PMat_Mu (1833)

   do i = 1, nj
     dist = ( xyz(i,1) - 1.0d0 )**2 + ( xyz(i,2) - 1.0d0 )**2 + xyz(i,3)**2
     if ( dist <= 1e-3 ) then
       write(*,'(i5,3f10.5)') i, xyz(i,1), xyz(i,2), xyz(i,3)
     end if
   end do

End If
Call MPI_Barrier ( PETSC_COMM_WORLD, ierr )
! ===================================================================


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

! 4- compute regularization vectors
      Call MatMult ( Reg_PETSC, Lambda_PETSC, Reg_Lambda_PETSC, IERR )
      Call MatMult ( Reg_PETSC,     Mu_PETSC,     Reg_Mu_PETSC, IERR )

! 4.2- compute the cost of the regularization component of the discrete objective functional
      Call VecDot ( Reg_Lambda_PETSC, Lambda_PETSC, Cost_Regularization_Lambda, IERR )
      Call VecDot (     Reg_Mu_PETSC,     Mu_PETSC, Cost_Regularization_Mu,     IERR )
      Cost_Regularization_Lambda = 0.50d0 * Cost_Regularization_Lambda 
      Cost_Regularization_Mu     = 0.50d0 * Cost_Regularization_Mu

! 6- apply the (inverse) mass matrix associated with the gradient vector
!      Call VecPointwiseMult ( Reg_Mu_PETSC, Reg_Mu_PETSC, iMg_PETSC, IERR )








Call MPI_Barrier ( PETSC_COMM_WORLD, ierr )

IF ( RANK == 0 ) then
   Write ( *     ,'(A138)' ) '-----------------------------------------------------------------------------------------------------------------------------------------'
   write(*,'(a25,e30.16)') 'cost reg mu =', Cost_Regularization_Mu
   Write ( *     ,'(A138)' ) '-----------------------------------------------------------------------------------------------------------------------------------------'
End If





  If ( Rank == 0 ) Then
    Write ( *     ,'(A138)' ) '-----------------------------------------------------------------------------------------------------------------------------------------'
    write(*,*) 'cost misfit', Cost_Reg_Mu
    Write ( *     ,'(A138)' ) '-----------------------------------------------------------------------------------------------------------------------------------------'
    write(*,'(a25,e30.16)') 'g_mu_1.1d0     =',  (0.5555555555573649E+00 ) / (5.0d0 * 0.1d0 )
    write(*,'(a25,e30.16)') 'g_mu_1.01d0    =',  (0.5555555557383282E-02 ) / (5.0d0 * 0.01d0 )
    write(*,'(a25,e30.16)') 'g_mu_1.001d0   =',  (0.5555555739205492E-04 ) / (5.0d0 * 0.001d0 )
    write(*,'(a25,e30.16)') 'g_mu_1.0001d0  =',  (0.5555573920727752E-06 ) / (5.0d0 * 0.0001d0 )
    write(*,'(a25,e30.16)') 'g_mu_1.00001d0 =',  (0.5557400954464727E-08 ) / (5.0d0 * 0.00001d0 )
    Write ( *     ,'(A138)') '-----------------------------------------------------------------------------------------------------------------------------------------'
  End If



call VecCopy ( Reg_Mu_PETSC, g_Mu_PETSC, IERR )

! ---------------- print out ----------------------------------------------------------------------

      CALL VecCreateSeq ( PETSC_COMM_SELF, NJ, g_Mu_Rank_PETSC, IERR ) 
      Call VecSet ( g_Mu_Rank_PETSC, 0.0d0, IERR )

      Call VecScatterBegin ( vscat_Mat, g_Mu_PETSC, g_Mu_Rank_PETSC, INSERT_VALUES, SCATTER_FORWARD, IERR )
      Call VecScatterEnd   ( vscat_Mat, g_Mu_PETSC, g_Mu_Rank_PETSC, INSERT_VALUES, SCATTER_FORWARD, IERR )


     Call VecGetArrayF90 ( g_Mu_Rank_PETSC, g_Mu_Rank, IERR )
!... transfer to appropriate Fortran matrices & ascertain Lambda is within pre-specified limits


IF ( RANK == 7 ) then
   Write ( *     ,'(A138)' ) '-----------------------------------------------------------------------------------------------------------------------------------------'
   write(*,'(a25,e30.16)') 'g_Mu (1833) =', g_Mu_Rank (1833)
   Write ( *     ,'(A138)' ) '-----------------------------------------------------------------------------------------------------------------------------------------'
End If


!... end transfer and engineering manipulation
     Call VecRestoreArrayF90 ( g_Mu_Rank_PETSC, g_Mu_Rank, IERR )







Call MPI_Barrier ( PETSC_COMM_WORLD, ierr )
stop
! Kids playground -----------------------------------------------------------------------------------------------------------------------------------
!  *                                                                                                        *
! ---------------------------------------------------------------------------------------------------------------------------------------------------



