
! If the path be beautiful, let us not ask where it leads.

!====================================================================================================================================================
! - DEFINE PETSC OBJECTS (Inverse Problem)
!====================================================================================================================================================

! I. building material vectors, associated mappings and scatters ------------------------------------------------------------------------------------
                                                                                                                                                    !
! 1- Lambda(big) , Mu(big)                                                                                                                          !
      CALL VecCreateMPI ( PETSC_COMM_WORLD, NJ_Mapping, NJTotal, Lambda_PETSC, IERR )                                                               !
      CALL VecCreateMPI ( PETSC_COMM_WORLD, NJ_Mapping, NJTotal,     Mu_PETSC, IERR )                                                               !
      CALL VecSet ( Lambda_PETSC, 0.0d0, IERR )                                                                                                     !
      CALL VecSet (     Mu_PETSC, 0.0d0, IERR )                                                                                                     !
                                                                                                                                                    !
      CALL VecCreateMPI ( PETSC_COMM_WORLD, NJ_Mapping, NJTotal, Lambda_PETSC_pre, IERR )                                                           !
      CALL VecCreateMPI ( PETSC_COMM_WORLD, NJ_Mapping, NJTotal,     Mu_PETSC_pre, IERR )                                                           !
      CALL VecSet ( Lambda_PETSC_pre, 0.0d0, IERR )                                                                                                 !
      CALL VecSet (     Mu_PETSC_pre, 0.0d0, IERR )                                                                                                 !
                                                                                                                                                    !
! 2- Lambda(Bigger) , Mu(Bigger): obtained via scattering from Lambda(big) and Mu(big)                                                              !
! remark: NJ_Rank_IParts = NJ                                                                                                                       !
      CALL VecCreateSeq ( PETSC_COMM_SELF, NJ_Rank_IParts, Lambda_Rank_PETSC, IERR )                                                                !
      CALL VecCreateSeq ( PETSC_COMM_SELF, NJ_Rank_IParts,     Mu_Rank_PETSC, IERR )                                                                !
                                                                                                                                                    !
! 3- Index sets                                                                                                                                     !
      Call ISCreateGeneral ( PETSC_COMM_WORLD, NJ_Rank_IParts, idx_Mat_from, PETSC_COPY_VALUES, IS_Mat_from, IERR )                                 !
      Call ISCreateGeneral ( PETSC_COMM_WORLD, NJ_Rank_IParts, idx_Mat_to,   PETSC_COPY_VALUES, IS_Mat_to,   IERR )                                 !
!     DeAllocate ( idx_Mat_from, idx_Mat_to )                                                                                                       !
!     DeAllocate ( idx_Mat_to ) ! need idx_Mat_from later for control problem operators                                                             !
                                                                                                                                                    !
! 4- scatter context                                                                                                                                !
      Call VecScatterCreate ( Lambda_PETSC, IS_Mat_from, Lambda_Rank_PETSC, IS_Mat_to, vscat_Mat, IERR )                                            !
                                                                                                                                                    !
! 5- scattering in action                                                                                                                           !
!      Call VecScatterBegin ( vscat_Mat, Lambda_PETSC, Lambda_Rank_PETSC, INSERT_VALUES, SCATTER_FORWARD, IERR )                                    !
!      Call VecScatterEnd   ( vscat_Mat, Lambda_PETSC, Lambda_Rank_PETSC, INSERT_VALUES, SCATTER_FORWARD, IERR )                                    !
!      Call VecScatterBegin ( vscat_Mat,     Mu_PETSC,     Mu_Rank_PETSC, INSERT_VALUES, SCATTER_FORWARD, IERR )                                    !
!      Call VecScatterEnd   ( vscat_Mat,     Mu_PETSC,     Mu_Rank_PETSC, INSERT_VALUES, SCATTER_FORWARD, IERR )                                    !
                                                                                                                                                    !
! 6- get it in Fortran format for system matrix assembly                                                                                            !
!      Call VecGetArrayF90     ( Lambda_Rank_PETSC, Lambda_Rank, IERR )                                                                             !
!... assemble system matrices ...                                                                                                                   !
!      Call VecRestoreArrayF90 ( Lambda_Rank_PETSC, Lambda_Rank, IERR )                                                                             !
!      Call VecGetArrayF90         ( Mu_Rank_PETSC,     Mu_Rank, IERR )                                                                             !
!... assemble system matrices ...                                                                                                                   !
!      Call VecRestoreArrayF90     ( Mu_Rank_PETSC,     Mu_Rank, IERR )                                                                             !
                                                                                                                                                    !
! II. material extension from interface into PML domain ---------------------------------------------------------------------------------------------
                                                                                                                                                    !
! 1- Index sets                                                                                                                                     !
      Call ISCreateGeneral ( PETSC_COMM_WORLD, NJ_Mapping, idx_Mat_Extend_from, PETSC_COPY_VALUES, is_mat_extend_from, IERR )                       !
      Call ISCreateGeneral ( PETSC_COMM_WORLD, NJ_Mapping, idx_Mat_Extend_to,   PETSC_COPY_VALUES, is_mat_extend_to,   IERR )                       !
      Call MPI_Barrier ( PETSC_COMM_WORLD, IERR )                                                                                                   !
      DeAllocate ( idx_Mat_Extend_from, idx_Mat_Extend_to )                                                                                         !
!     Call ISView ( IS_Mat_Extend_to, PETSC_VIEWER_STDOUT, IERR )                                                                                   !
                                                                                                                                                    !
! 2- scatter context                                                                                                                                !
      Call VecScatterCreate ( Lambda_PETSC, is_mat_extend_from, Lambda_PETSC, is_mat_extend_to, vscat_mat_extend, IERR )                            !
!     Call VecScatterView   ( vscat_mat_extend, PETSC_VIEWER_STDOUT, IERR )                                                                         !
                                                                                                                                                    !
! 3- scattering in action                                                                                                                           !
!      Call VecScatterBegin ( vscat_Mat_Extend, Lambda_PETSC, Lambda_PETSC, INSERT_VALUES, SCATTER_FORWARD, IERR )                                  !
!      Call VecScatterEnd   ( vscat_Mat_Extend, Lambda_PETSC, Lambda_PETSC, INSERT_VALUES, SCATTER_FORWARD, IERR )                                  !
!      Call VecScatterBegin ( vscat_Mat_Extend,     Mu_PETSC,     Mu_PETSC, INSERT_VALUES, SCATTER_FORWARD, IERR )                                  !
!      Call VecScatterEnd   ( vscat_Mat_Extend,     Mu_PETSC,     Mu_PETSC, INSERT_VALUES, SCATTER_FORWARD, IERR )                                  !
                                                                                                                                                    !
! III. data structure for storing displacement, associated mappings and scatters --------------------------------------------------------------------

! 1- create a dense matrix to store U in the forward problem: storing at time steps 0:N. We assume U at t=0 is zero.
      Allocate (    U_Store_Mapping ( NStore_Mapping , 0:NStep ) )
      Allocate (    P_Store_Mapping ( NStore_Mapping , 0:NStep ) )
!      Allocate ( k1_m_Store_Mapping ( NStore_Mapping , 0:NStep ) )
!      Allocate ( k2_m_Store_Mapping ( NStore_Mapping , 0:NStep ) )
!      Allocate ( k3_m_Store_Mapping ( NStore_Mapping , 0:NStep ) )
!     Allocate ( k4_m_Store_Mapping ( NStore_Mapping , 0:NStep ) )

         U_Store_Mapping ( : , : ) = 0.0d0
         P_Store_Mapping ( : , : ) = 0.0d0
!      k1_m_Store_Mapping ( : , : ) = 0.0d0
!      k2_m_Store_Mapping ( : , : ) = 0.0d0
!      k3_m_Store_Mapping ( : , : ) = 0.0d0
!      k4_m_Store_Mapping ( : , : ) = 0.0d0

! 2- Index sets
      Call ISCreateGeneral ( PETSC_COMM_WORLD, NStore_Rank, idx_u_from, PETSC_COPY_VALUES, IS_u_from, IERR ) ! from Global long-thin vector
      Call ISCreateGeneral ( PETSC_COMM_WORLD, NStore_Rank, idx_u_to,   PETSC_COPY_VALUES, IS_u_to,   IERR ) ! to Local sequential vectors
      DeAllocate ( idx_u_from, idx_u_to )

! 3- U_big_PETSc: contains state variable solution only within the regural domain. this is the long-thin vector that we will scatter from, into local sequential vectors for gradient computation.
      CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, Ux_big_PETSc, IERR )
      CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, Uk_big_PETSc, IERR )
      CALL VecSet ( Ux_big_PETSc, 0.0d0, IERR )
      CALL VecSet ( Uk_big_PETSc, 0.0d0, IERR )

! 4- U_Bigger_Rank_PETSc: obtained via scattering from U_big_PETSc
      CALL VecCreateSeq ( PETSC_COMM_SELF, NEQM, Ux_Bigger_Rank_PETSc, IERR )
      CALL VecCreateSeq ( PETSC_COMM_SELF, NEQM, Uk_Bigger_Rank_PETSc, IERR )
      CALL VecCreateSeq ( PETSC_COMM_SELF, NEQM,  P_Bigger_Rank_PETSc, IERR )
      CALL VecSet ( Ux_Bigger_Rank_PETSc, 0.0d0, IERR )
      CALL VecSet ( Uk_Bigger_Rank_PETSc, 0.0d0, IERR )
      CALL VecSet (  P_Bigger_Rank_PETSc, 0.0d0, IERR )

! 5- scatter context
      Call VecScatterCreate ( Ux_big_PETSc, IS_u_from, Ux_Bigger_Rank_PETSc, IS_u_to, vscat_u_hist, IERR )
!     Call VecScatterView   ( vscat_u_hist, PETSC_VIEWER_STDOUT, IERR ) 

! 6- scattering in action
!     Call VecScatterBegin ( vscat_u_hist, Ux_big_PETSc, U_Bigger_Rank_PETSc, INSERT_VALUES, SCATTER_FORWARD, IERR )
!     Call VecScatterEnd   ( vscat_u_hist, Ux_big_PETSc, U_Bigger_Rank_PETSc, INSERT_VALUES, SCATTER_FORWARD, IERR )

! IV. vectors for gradient and search direction -----------------------------------------------------------------------------------------------------
      CALL VecCreateMPI ( PETSC_COMM_WORLD, NJ_Mapping, NJTotal, g_Lambda_PETSC, IERR )
      CALL VecCreateMPI ( PETSC_COMM_WORLD, NJ_Mapping, NJTotal,     g_Mu_PETSC, IERR )
      CALL VecCreateMPI ( PETSC_COMM_WORLD, NJ_Mapping, NJTotal, p_Lambda_PETSC, IERR )
      CALL VecCreateMPI ( PETSC_COMM_WORLD, NJ_Mapping, NJTotal,     p_Mu_PETSC, IERR )

      CALL VecCreateMPI ( PETSC_COMM_WORLD, NJ_Mapping, NJTotal, g_Lambda_PETSC_pre, IERR )
      CALL VecCreateMPI ( PETSC_COMM_WORLD, NJ_Mapping, NJTotal,     g_Mu_PETSC_pre, IERR )
      CALL VecCreateMPI ( PETSC_COMM_WORLD, NJ_Mapping, NJTotal, p_Lambda_PETSC_pre, IERR )
      CALL VecCreateMPI ( PETSC_COMM_WORLD, NJ_Mapping, NJTotal,     p_Mu_PETSC_pre, IERR )

      Call VecSet ( g_Lambda_PETSC, 0.0d0, IERR )
      Call VecSet ( g_Mu_PETSC,     0.0d0, IERR )

! V. vectors for computing the adjoint problem ------------------------------------------------------------------------------------------------------

! 1- misfit vector
      CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, B_mis_PETSC, IERR )
      CALL VecSetOption ( B_mis_PETSC, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE  , IERR )
! comment: this is not required since we use Global Equation Numbers for constructing the misfit vector.
!     Call VecSetLocalToGlobalMapping ( B_misfit_PETSC,   GMapping, ErrPTC ) ;

! VI. vectors and matrices for computing the control problem ----------------------------------------------------------------------------------------

! 1- regularization matrix R (discrete laplace operator)
!      D_NNZ_Reg (:) = 80 !125
!      O_NNZ_Reg (:) = 80 !125
      D_NNZ_Reg (:) = 125
      O_NNZ_Reg (:) = 125

      CALL MatCreateAIJ ( PETSC_COMM_WORLD, NJ_Mapping, NJ_Mapping, NJTotal, NJTotal, 0, D_NNZ_Reg, 0, O_NNZ_Reg, Reg_PETSC, IERR ) ! local set values
      CALL MatSetOption ( Reg_PETSC, MAT_IGNORE_ZERO_ENTRIES, PETSC_TRUE, IERR ) ;

! 2- inverse mass of the gradient vector
      CALL VecCreateMPI ( PETSC_COMM_WORLD, NJ_Mapping, NJTotal, iMg_PETSC, IERR ) ! local set values
      CALL VecSet ( iMg_PETSC, 0.0d0, IERR )

! 3- regularization vectors
      CALL VecCreateMPI ( PETSC_COMM_WORLD, NJ_Mapping, NJTotal, Reg_Lambda_PETSC, IERR )     ! Reg_Lambda = Reg * Lambda
      CALL VecCreateMPI ( PETSC_COMM_WORLD, NJ_Mapping, NJTotal,     Reg_Mu_PETSC, IERR )
      CALL VecSet ( Reg_Lambda_PETSC, 0.0d0, IERR )
      CALL VecSet (     Reg_Mu_PETSC, 0.0d0, IERR )
! comment: TV regularization is implemented in a matrix-free fasion. Therefore, we assemble directly into regularization vectors. We need to set values locally, which is done in VIII-3.

! VII. Local to Global Mapping for gradient and Laplacian assembly (material nodes) -----------------------------------------------------------------

! 1- create the mapping. Note: See PIC_Inversion_Data_Structure.F90 (160)
      Call ISLocalToGlobalMappingCreate ( PETSC_COMM_WORLD, NJ_Rank_IParts, idx_Mat_from, PETSC_COPY_VALUES, GMapping_Mat, IERR )

! 2- applying Local to Global Numbering for PETSc objects
      Call MatSetLocalToGlobalMapping ( Reg_PETSC, GMapping_Mat, GMapping_Mat, IERR )

      Call VecSetLocalToGlobalMapping (      iMg_PETSC, GMapping_Mat, IERR ) ;
      Call VecSetLocalToGlobalMapping (     g_Mu_PETSC, GMapping_Mat, IERR ) ;
      Call VecSetLocalToGlobalMapping ( g_Lambda_PETSC, GMapping_Mat, IERR ) ;

! 3- Local to Global Numbering for regularization (from VI)
      Call VecSetLocalToGlobalMapping (     Reg_Mu_PETSC, GMapping_Mat, IERR ) ;
      Call VecSetLocalToGlobalMapping ( Reg_Lambda_PETSC, GMapping_Mat, IERR ) ;

!====================================================================================================================================================
! - DEFINE PETSC OBJECTS (Material Visualization using Paraview)
!====================================================================================================================================================

! VII. material visualization -----------------------------------------------------------------------------------------------------------------------

! 1- define Index Set for Paraview ------------------------------------------------------------------------------------------------------------------
      Allocate ( idx_Mat_Vis_to ( NJ_Para ) ) ;                                          ! (NJ 20-noded FEM for visualization)
      ForAll ( I = 1 : NJ_Para ) idx_Mat_Vis_to (I) = I-1 ;

      Call VecCreateSeq ( PETSC_COMM_SELF, NJ_Para, Lambda_Rank_Vis_PETSC, ErrPTC );     ! Temporary vector to hold Lambda in 20-nodes for Paraview
      Call VecCreateSeq ( PETSC_COMM_SELF, NJ_Para,     Mu_Rank_Vis_PETSC, ErrPTC );     ! Temporary vector to hold Mu in 20-nodes for Paraview

      Call ISCreateGeneral ( PETSC_COMM_SELF, NJ_Para, idx_Mat_Vis_from, PETSC_COPY_VALUES, IS_Mat_Vis_from, ErrPTC );
      Call ISCreateGeneral ( PETSC_COMM_SELF, NJ_Para, idx_Mat_Vis_to,   PETSC_COPY_VALUES, IS_Mat_Vis_to,   ErrPTC );

      DeAllocate ( idx_Mat_Vis_to ) ;

      Call VecScatterCreate ( Lambda_PETSC, IS_Mat_Vis_from, Lambda_Rank_Vis_PETSC, IS_Mat_Vis_to, vscat_Mat_Vis, IERR )

! 2- Write total wraper file for Paraview (only once at the beginning) ------------------------------------------------------------------------------
      If ( Rank == 0 ) Then ;  ! Top file wraper 
         UnFile = UN_OutallWrp_Mat_Vis ;
         Open ( Unit = UnFile, FILE = TRIM(AnaName)//'_'//'all'//'_'//Trim(AdjustL(IndexSize))//'.xmf',  ERR =  1001, IOSTAT = IO_File, ACCESS = 'SEQUENTIAL', ACTION ='Write', ASYNCHRONOUS = 'NO', BLANK = 'NULL', BLOCKSIZE = 0, DEFAULTFILE = TRIM(MatUpdateDir), DisPOSE = 'KEEP', FORM = 'FORMATTED', POSITION = 'ASIS', STATUS ='REPLACE') ;

         Write (Unit = UnFile, FMT = "(A63)", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) '<Xdmf xmlns:xi="http://www.w3.org/2003/XInclude" Version="2.2">' ;
         Write (Unit = UnFile, FMT = "(A10)", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) '  <Domain>' ;
         Write (Unit = UnFile, FMT = "(A57)", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) "    <Grid GridType='Collection' CollectionType='Spatial'>" ;
      End If ;
! - -------------------------------------------------------------------------------------------------------------------------------------------------



!====================================================================================================================================================
! - Vectors for adjoint solution
!====================================================================================================================================================


CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, y_t_PETSC, IERR )
CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, y_m_PETSC, IERR )
CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, y_b_PETSC, IERR )

CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, p1_t_PETSC, IERR )
CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, p1_m_PETSC, IERR )
CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, p1_b_PETSC, IERR )

CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, p2_t_PETSC, IERR )
CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, p2_m_PETSC, IERR )
CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, p2_b_PETSC, IERR )

CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, p3_t_PETSC, IERR )
CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, p3_m_PETSC, IERR )
CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, p3_b_PETSC, IERR )

CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, p4_t_PETSC, IERR )
CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, p4_m_PETSC, IERR )
CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, p4_b_PETSC, IERR )

CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, help_PETSC, IERR )


!====================================================================================================================================================
! - L-BFGS
!====================================================================================================================================================

! 1- vectors
!      CALL VecCreateMPI ( PETSC_COMM_WORLD, NJ_Mapping, NJTotal, q_Lambda_PETSC, IERR )
!      CALL VecCreateMPI ( PETSC_COMM_WORLD, NJ_Mapping, NJTotal,     q_Mu_PETSC, IERR )
      Allocate ( q_Lambda ( NJ_Mapping ) )
      Allocate (     q_Mu ( NJ_Mapping ) )
      Allocate ( r_Lambda ( NJ_Mapping ) )
      Allocate (     r_Mu ( NJ_Mapping ) )

      q_Lambda ( : ) = 0.0d0
      q_Mu     ( : ) = 0.0d0
      r_Lambda ( : ) = 0.0d0
      r_Mu     ( : ) = 0.0d0

      CALL VecCreateMPI ( PETSC_COMM_WORLD, NJ_Mapping, NJTotal, y_Lambda_LBFGS_PETSC, IERR )
      CALL VecCreateMPI ( PETSC_COMM_WORLD, NJ_Mapping, NJTotal,     y_Mu_LBFGS_PETSC, IERR )
      CALL VecCreateMPI ( PETSC_COMM_WORLD, NJ_Mapping, NJTotal, s_Lambda_LBFGS_PETSC, IERR )
      CALL VecCreateMPI ( PETSC_COMM_WORLD, NJ_Mapping, NJTotal,     s_Mu_LBFGS_PETSC, IERR )

! 2- dense matrices, comprising the m vectors corresponding to previous iterations.
!      Call MatCreateDense ( PETSC_COMM_WORLD, NJ_Mapping, PETSC_DECIDE, NJTotal, M_LBFGS, NULL_SCALAR, Ys_Lambda_LBFGS_PETSC, IERR )
!      Call MatCreateDense ( PETSC_COMM_WORLD, NJ_Mapping, PETSC_DECIDE, NJTotal, M_LBFGS, NULL_SCALAR,     Ys_Mu_LBFGS_PETSC, IERR )
!      Call MatCreateDense ( PETSC_COMM_WORLD, NJ_Mapping, PETSC_DECIDE, NJTotal, M_LBFGS, NULL_SCALAR, Ss_Lambda_LBFGS_PETSC, IERR )
!      Call MatCreateDense ( PETSC_COMM_WORLD, NJ_Mapping, PETSC_DECIDE, NJTotal, M_LBFGS, NULL_SCALAR,     Ss_Mu_LBFGS_PETSC, IERR )

! 3- allocate memory for scalar parameters
      Allocate (   Rhos_Lambda_LBFGS (M_LBFGS) )
      Allocate (       Rhos_Mu_LBFGS (M_LBFGS) )
      Allocate (       Rhos_LM_LBFGS (M_LBFGS) )
      Allocate ( Alphas_Lambda_LBFGS (M_LBFGS) )
      Allocate (     Alphas_Mu_LBFGS (M_LBFGS) )
      Allocate (     Alphas_LM_LBFGS (M_LBFGS) )
      Allocate ( Alphas_Lambda_LBFGS_Rank (M_LBFGS) )     
      Allocate (     Alphas_Mu_LBFGS_Rank (M_LBFGS) )
      Allocate (     Alphas_LM_LBFGS_Rank (M_LBFGS) )
      Allocate (    LBFGS_vector_sequence (M_LBFGS) )

! 4- allocate memory for m column vectors, corresponding to previous iterations
      Allocate ( Ys_Lambda_LBFGS ( NJ_Mapping , M_LBFGS ) )
      Allocate (     Ys_Mu_LBFGS ( NJ_Mapping , M_LBFGS ) )
      Allocate ( Ss_Lambda_LBFGS ( NJ_Mapping , M_LBFGS ) )
      Allocate (     Ss_Mu_LBFGS ( NJ_Mapping , M_LBFGS ) )

      Ys_Lambda_LBFGS ( : , : ) = 0.0d0
          Ys_Mu_LBFGS ( : , : ) = 0.0d0
      Ss_Lambda_LBFGS ( : , : ) = 0.0d0
          Ss_Mu_LBFGS ( : , : ) = 0.0d0

! 5- restriction of a vector to the regular domain (RD = 1, PML = 0)
      Allocate ( restriction_RD ( NJ_Mapping ) )












































