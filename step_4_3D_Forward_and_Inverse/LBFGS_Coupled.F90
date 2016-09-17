! *************************************************************************************************
!                                                L  B  F  G  S - Coupled
! *************************************************************************************************
           If ( Iter == Iter_begin ) Then
              Call VecAXPY ( p_Lambda_PETSC, -1.0d0, g_Lambda_PETSC, IERR )
              Call VecAXPY (     p_Mu_PETSC, -1.0d0,     g_Mu_PETSC, IERR )
           Else
! 2.1- LBFGS k
             K_LBFGS = Iter - Iter_begin

! 2.2- compute y, s, rho; form q
! comment: Lambda & Mu encompass the entire domain; however, g belongs only to the regular domain. So, the dot-product belongs to the regular domain.
! Lambda
             Call VecSet   ( y_Lambda_LBFGS_PETSC, 0.0d0, IERR )
             Call VecMAXPY ( y_Lambda_LBFGS_PETSC, 2, [ 1.0d0, -1.0d0 ], [ g_Lambda_PETSC, g_Lambda_PETSC_pre], IERR )
             Call VecSet   ( s_Lambda_LBFGS_PETSC, 0.0d0, IERR )
             Call VecMAXPY ( s_Lambda_LBFGS_PETSC, 2, [ 1.0d0, -1.0d0 ], [   Lambda_PETSC,   Lambda_PETSC_pre], IERR )
             Call VecDot   ( y_Lambda_LBFGS_PETSC, s_Lambda_LBFGS_PETSC, rho_Lambda_LBFGS, IERR )
! Mu
             Call VecSet   (     y_Mu_LBFGS_PETSC, 0.0d0, IERR )
             Call VecMAXPY (     y_Mu_LBFGS_PETSC, 2, [ 1.0d0, -1.0d0 ], [     g_Mu_PETSC,     g_Mu_PETSC_pre], IERR )
             Call VecSet   (     s_Mu_LBFGS_PETSC, 0.0d0, IERR )
             Call VecMAXPY (     s_Mu_LBFGS_PETSC, 2, [ 1.0d0, -1.0d0 ], [       Mu_PETSC,       Mu_PETSC_pre], IERR )
             Call VecDot   ( y_Mu_LBFGS_PETSC, s_Mu_LBFGS_PETSC, rho_Mu_LBFGS, IERR )
! rho
             rho_LM_LBFGS = 1.0d0 / ( rho_Lambda_LBFGS + rho_Mu_LBFGS )

!2.3- warn user if rho is large
             Tol_rho_LBFGS = 1.0d30
             If ( Dabs (rho_LM_LBFGS) > Tol_rho_LBFGS ) Then
                If ( rank == 0 ) Then
                   Write(*    ,'(A33,E33.10E3)') 'rho_LM_LBFGS = ', rho_LM_LBFGS
                   Write(UnInf,'(A33,E33.10E3)') 'rho_LM_LBFGS = ', rho_LM_LBFGS
                End If
             End If

! 2.2.1- form the q vector: This is okay, we access the local values of a global vector (test_06).
! Lambda
             Call VecGetArrayF90 ( g_Lambda_PETSC, g_Lambda_Rank, IERR )
             ForAll ( I = 1:NJ_Mapping ) q_Lambda ( I ) = g_Lambda_Rank ( I )
             Call VecRestoreArrayF90 ( g_Lambda_PETSC, g_Lambda_Rank, IERR )
! Mu
             Call VecGetArrayF90 ( g_Mu_PETSC, g_Mu_Rank, IERR )
             ForAll ( I = 1:NJ_Mapping ) q_Mu ( I ) = g_Mu_Rank ( I )
             Call VecRestoreArrayF90 ( g_Mu_PETSC, g_Mu_Rank, IERR )

! 2.3- the algorithm
             If ( (K_LBFGS > 0) .AND. (K_LBFGS <= M_LBFGS) ) Then ! early stages of LBFGS                                   --- E A R L Y --- S T A G E S --- O F --- L B F G S
                N_use = K_LBFGS ! number of vectors to be used for inverse Hessian approximation
                Lcounter = 0
! 2.3.1- vector sequence
                Do I = K_LBFGS, 1, -1
                   Lcounter = Lcounter + 1
                   LBFGS_vector_sequence ( Lcounter ) = I
                End Do
! 2.3.2- where to store the latest vectors                
                If ( (mod ( K_LBFGS , M_LBFGS )) /= 0 ) Then
                   N_store = mod ( K_LBFGS , M_LBFGS )
                Else
                   N_store = M_LBFGS
                End If
! 2.3.3- store scalars
                Rhos_LM_LBFGS ( N_store ) = rho_LM_LBFGS
! 2.3.4- store vectors: y_Lambda : access local values of a global vector
                Call VecGetArrayF90 ( y_Lambda_LBFGS_PETSC, y_Lambda_LBFGS_Rank, IERR )
                ForAll ( I = 1:NJ_Mapping ) Ys_Lambda_LBFGS ( I , N_store ) = y_Lambda_LBFGS_Rank (I)
                Call VecRestoreArrayF90 ( y_Lambda_LBFGS_PETSC, y_Lambda_LBFGS_Rank, IERR )
! 2.3.5- store vectors: y_Mu
                Call VecGetArrayF90 ( y_Mu_LBFGS_PETSC, y_Mu_LBFGS_Rank, IERR )
                ForAll ( I = 1:NJ_Mapping ) Ys_Mu_LBFGS ( I , N_store ) = y_Mu_LBFGS_Rank (I)
                Call VecRestoreArrayF90 ( y_Mu_LBFGS_PETSC, y_Mu_LBFGS_Rank, IERR )
! 2.3.6- store vectors: s_Lambda
                Call VecGetArrayF90 ( s_Lambda_LBFGS_PETSC, s_Lambda_LBFGS_Rank, IERR )
                ForAll ( I = 1:NJ_Mapping ) Ss_Lambda_LBFGS ( I , N_store ) = s_Lambda_LBFGS_Rank (I) * restriction_RD (I)
                Call VecRestoreArrayF90 ( s_Lambda_LBFGS_PETSC, s_Lambda_LBFGS_Rank, IERR )
! 2.3.7- store vectors: s_Mu
                Call VecGetArrayF90 ( s_Mu_LBFGS_PETSC, s_Mu_LBFGS_Rank, IERR )
                ForAll ( I = 1:NJ_Mapping ) Ss_Mu_LBFGS ( I , N_store ) = s_Mu_LBFGS_Rank (I) * restriction_RD (I)
                Call VecRestoreArrayF90 ( s_Mu_LBFGS_PETSC, s_Mu_LBFGS_Rank, IERR )

! 2.3.8- loop 1 =================================================================
                Do I_vec = 1 , N_use
                   I = LBFGS_vector_sequence ( I_vec )
                   Alphas_LM_LBFGS_Rank ( I ) = Rhos_LM_LBFGS ( I ) * Dot_Product ( q_Lambda (:) , Ss_Lambda_LBFGS ( : , I ) ) &
                                              + Rhos_LM_LBFGS ( I ) * Dot_Product ( q_Mu     (:) , Ss_Mu_LBFGS     ( : , I ) )
                   Call MPI_Barrier ( PETSc_COMM_WORLD, IERR )
                   Call MPI_Allreduce ( Alphas_LM_LBFGS_Rank ( I ), Alphas_LM_LBFGS ( I ), 1, MPI_DOUBLE_PRECISION, MPI_SUM, PETSC_COMM_WORLD, IERR )
                   Call MPI_Barrier ( PETSc_COMM_WORLD, IERR )
                   ! Lambda
                   q_Lambda (:) = q_Lambda (:) - Alphas_LM_LBFGS ( I ) * Ys_Lambda_LBFGS ( : , I )
                   ! Mu
                   q_Mu (:)     = q_Mu (:)     - Alphas_LM_LBFGS ( I ) * Ys_Mu_LBFGS ( : , I )
                End Do
! 2.3.9- H^0 ====================================================================
                Call VecDot ( y_Lambda_LBFGS_PETSC, s_Lambda_LBFGS_PETSC, H0_Lambda_numerator,   IERR )
                Call VecDot ( y_Lambda_LBFGS_PETSC, y_Lambda_LBFGS_PETSC, H0_Lambda_denominator, IERR )
                Call VecDot (     y_Mu_LBFGS_PETSC,     s_Mu_LBFGS_PETSC, H0_Mu_numerator,       IERR )
                Call VecDot (     y_Mu_LBFGS_PETSC,     y_Mu_LBFGS_PETSC, H0_Mu_denominator,     IERR )
                H0_LM_numerator   = H0_Lambda_numerator   + H0_Mu_numerator
                H0_LM_denominator = H0_Lambda_denominator + H0_Mu_denominator
                ! Lambda
                r_Lambda = ( H0_LM_numerator / H0_LM_denominator ) * q_Lambda
                ! Mu
                r_Mu     = ( H0_LM_numerator / H0_LM_denominator ) * q_Mu
! 2.3.10- loop 2 ================================================================
                Do I_vec = N_use, 1, -1
                   I = LBFGS_vector_sequence ( I_vec )
                   beta_LM_LBFGS_Rank = Rhos_LM_LBFGS ( I ) * Dot_Product ( r_Lambda (:) , Ys_Lambda_LBFGS ( : , I ) ) &
                                      + Rhos_LM_LBFGS ( I ) * Dot_Product ( r_Mu     (:) , Ys_Mu_LBFGS     ( : , I ) )
                   Call MPI_Barrier ( PETSc_COMM_WORLD, IERR )
                   Call MPI_Allreduce ( beta_LM_LBFGS_Rank, beta_LM_LBFGS, 1, MPI_DOUBLE_PRECISION, MPI_SUM, PETSC_COMM_WORLD, IERR )
                   Call MPI_Barrier ( PETSc_COMM_WORLD, IERR )
                   ! Lambda
                   r_Lambda (:) = r_Lambda (:) + Ss_Lambda_LBFGS ( : , I ) * ( Alphas_LM_LBFGS ( I ) - beta_LM_LBFGS )
                   ! Mu
                   r_Mu (:) = r_Mu (:) + Ss_Mu_LBFGS ( : , I ) * ( Alphas_LM_LBFGS ( I ) - beta_LM_LBFGS )
                End Do
!                ================================================================

! 2.3.11- end of early stages of LBFGS

             Else ! middle stages of LBFGS                                                                                  --- M I D D L E --- S T A G E S --- O F --- L B F G S
                N_use = M_LBFGS ! number of vectors to be used for inverse Hessian approximation
                LEnd_sequence = mod ( K_LBFGS , M_LBFGS )                                                                                            ! 
                Lcounter = 0
! 2.3.12- vector sequence: part 1
                Do I = LEnd_sequence, 1, -1
                   Lcounter = Lcounter + 1
                   LBFGS_vector_sequence ( Lcounter ) = I
                End Do
! 2.3.13- vector sequence: part 2
                J = 0
                Do I = Lcounter + 1 , M_LBFGS
                   LBFGS_vector_sequence ( I ) = M_LBFGS - J
                   J = J + 1
                End Do
! 2.3.14- where to store the latest vectors
                If ( (mod ( K_LBFGS , M_LBFGS )) /= 0 ) Then
                   N_store = mod ( K_LBFGS , M_LBFGS )
                Else
                   N_store = M_LBFGS
                End If
! 2.3.15- store scalars
                Rhos_LM_LBFGS ( N_store ) = rho_LM_LBFGS
! 2.3.16- store vectors: y_Lambda
                Call VecGetArrayF90 ( y_Lambda_LBFGS_PETSC, y_Lambda_LBFGS_Rank, IERR )
                ForAll ( I = 1:NJ_Mapping ) Ys_Lambda_LBFGS ( I , N_store ) = y_Lambda_LBFGS_Rank (I)
                Call VecRestoreArrayF90 ( y_Lambda_LBFGS_PETSC, y_Lambda_LBFGS_Rank, IERR )
! 2.3.17- store vectors: y_Mu
                Call VecGetArrayF90 ( y_Mu_LBFGS_PETSC, y_Mu_LBFGS_Rank, IERR )
                ForAll ( I = 1:NJ_Mapping ) Ys_Mu_LBFGS ( I , N_store ) = y_Mu_LBFGS_Rank (I)
                Call VecRestoreArrayF90 ( y_Mu_LBFGS_PETSC, y_Mu_LBFGS_Rank, IERR )
! 2.3.18- store vectors: s_Lambda
                Call VecGetArrayF90 ( s_Lambda_LBFGS_PETSC, s_Lambda_LBFGS_Rank, IERR )
                ForAll ( I = 1:NJ_Mapping ) Ss_Lambda_LBFGS ( I , N_store ) = s_Lambda_LBFGS_Rank (I) * restriction_RD (I)
                Call VecRestoreArrayF90 ( s_Lambda_LBFGS_PETSC, s_Lambda_LBFGS_Rank, IERR )
! 2.3.19- store vectors: s_Mu
                Call VecGetArrayF90 ( s_Mu_LBFGS_PETSC, s_Mu_LBFGS_Rank, IERR )
                ForAll ( I = 1:NJ_Mapping ) Ss_Mu_LBFGS ( I , N_store ) = s_Mu_LBFGS_Rank (I) * restriction_RD (I)
                Call VecRestoreArrayF90 ( s_Mu_LBFGS_PETSC, s_Mu_LBFGS_Rank, IERR )

! 2.3.20- loop 1 =================================================================
                Do I_vec = 1 , N_use
                   I = LBFGS_vector_sequence ( I_vec )
                   Alphas_LM_LBFGS_Rank ( I ) = Rhos_LM_LBFGS ( I ) * Dot_Product ( q_Lambda (:) , Ss_Lambda_LBFGS ( : , I ) ) &
                                              + Rhos_LM_LBFGS ( I ) * Dot_Product ( q_Mu     (:) , Ss_Mu_LBFGS     ( : , I ) )
                   Call MPI_Barrier ( PETSc_COMM_WORLD, IERR )
                   Call MPI_Allreduce ( Alphas_LM_LBFGS_Rank ( I ), Alphas_LM_LBFGS ( I ), 1, MPI_DOUBLE_PRECISION, MPI_SUM, PETSC_COMM_WORLD, IERR )
                   Call MPI_Barrier ( PETSc_COMM_WORLD, IERR )
                   ! Lambda
                   q_Lambda (:) = q_Lambda (:) - Alphas_LM_LBFGS ( I ) * Ys_Lambda_LBFGS ( : , I )
                   ! Mu
                   q_Mu (:)     = q_Mu (:)     - Alphas_LM_LBFGS ( I ) * Ys_Mu_LBFGS ( : , I )
                End Do
! 2.3.21- H^0 ====================================================================
                Call VecDot ( y_Lambda_LBFGS_PETSC, s_Lambda_LBFGS_PETSC, H0_Lambda_numerator,   IERR )
                Call VecDot ( y_Lambda_LBFGS_PETSC, y_Lambda_LBFGS_PETSC, H0_Lambda_denominator, IERR )
                Call VecDot ( y_Mu_LBFGS_PETSC, s_Mu_LBFGS_PETSC, H0_Mu_numerator,   IERR )
                Call VecDot ( y_Mu_LBFGS_PETSC, y_Mu_LBFGS_PETSC, H0_Mu_denominator, IERR )
                H0_LM_numerator   = H0_Lambda_numerator   + H0_Mu_numerator
                H0_LM_denominator = H0_Lambda_denominator + H0_Mu_denominator
                ! Lambda
                r_Lambda = ( H0_LM_numerator / H0_LM_denominator ) * q_Lambda
                ! Mu
                r_Mu     = ( H0_LM_numerator / H0_LM_denominator ) * q_Mu
! 2.3.22- loop 2 ================================================================
                Do I_vec = N_use, 1, -1
                   I = LBFGS_vector_sequence ( I_vec )
                   beta_LM_LBFGS_Rank = Rhos_LM_LBFGS ( I ) * Dot_Product ( r_Lambda (:) , Ys_Lambda_LBFGS ( : , I ) ) &
                                      + Rhos_LM_LBFGS ( I ) * Dot_Product ( r_Mu     (:) , Ys_Mu_LBFGS     ( : , I ) )
                   Call MPI_Barrier ( PETSc_COMM_WORLD, IERR )
                   Call MPI_Allreduce ( beta_LM_LBFGS_Rank, beta_LM_LBFGS, 1, MPI_DOUBLE_PRECISION, MPI_SUM, PETSC_COMM_WORLD, IERR )
                   Call MPI_Barrier ( PETSc_COMM_WORLD, IERR )
                   ! Lambda
                   r_Lambda (:) = r_Lambda (:) + Ss_Lambda_LBFGS ( : , I ) * ( Alphas_LM_LBFGS ( I ) - beta_LM_LBFGS )
                   ! Mu
                   r_Mu (:) = r_Mu (:) + Ss_Mu_LBFGS ( : , I ) * ( Alphas_LM_LBFGS ( I ) - beta_LM_LBFGS )
                End Do
!                ================================================================

! 2.3.23- end if
             End If

! 2.3.24- assemble in PETSc distributed vectors
             ! Lambda
             Call VecGetArrayF90 ( p_Lambda_PETSC, p_Lambda_Rank, IERR )
             ForAll ( I = 1:NJ_Mapping ) p_Lambda_Rank ( I ) = r_Lambda ( I )
             Call VecRestoreArrayF90 ( p_Lambda_PETSC, p_Lambda_Rank, IERR )
             ! Mu
             Call VecGetArrayF90 ( p_Mu_PETSC, p_Mu_Rank, IERR )
             ForAll ( I = 1:NJ_Mapping ) p_Mu_Rank ( I ) = r_Mu ( I )
             Call VecRestoreArrayF90 ( p_Mu_PETSC, p_Mu_Rank, IERR )

! 2.4- spit out the BFGS search direction
             Call VecScale ( p_Lambda_PETSC, -1.0d0, IERR )
             Call VecScale (     p_Mu_PETSC, -1.0d0, IERR )

           End If
! *************************************************************************************************
!                                                L  B  F  G  S    (end)
! *************************************************************************************************
