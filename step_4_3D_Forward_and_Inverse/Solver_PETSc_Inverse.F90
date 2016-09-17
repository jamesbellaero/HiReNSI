
!************************************************
! Explicit time stepping of the adjoint & control problem, based on Discretize-Then-Optimized Runge-Kutta method, 4 stage-4th order, PML-3D
!************************************************
! REVISION : Sunday 16 January 2014: A bee is never as busy as it seems; it's just that it can't buzz any slower.
! REVISION : Tuesday 5 February 2014: (77-18-77)

SUBROUTINE EXPLICIT_RK4_PETSC_PML_3D_Adjoint_Control ( AK_PETSC, DIAG_M_PETSC, DIAG_iM_PETSC, AC_PETSC, AG_PETSC, B_mis_PETSC, ID, RANK, &
                                                       NNDH, EqDis, U_Store_Mapping, NStore_Mapping, U_Store_Numbers_Global, EqDis_MAP_U_Store_Numbers_Global, Dis_meas, &
                                                       k1_m_Store_Mapping, k2_m_Store_Mapping, k3_m_Store_Mapping, &
                                                       Ux_big_PETSc, Ux_Bigger_Rank_PETSc, vscat_u_hist, P_Bigger_Rank_PETSc, Uk_big_PETSc, Uk_Bigger_Rank_PETSc, &
                                                       g_Lambda_PETSC, g_Mu_PETSC, idx_Mat_to, NJ_Rank_IParts, &
                                                       XYZ, INOD, NGP, MTEL, &
y_t_PETSC,  y_m_PETSC,  y_b_PETSC, p1_t_PETSC, p1_m_PETSC, p1_b_PETSC, p2_t_PETSC, p2_m_PETSC, p2_b_PETSC, p3_t_PETSC, p3_m_PETSC, p3_b_PETSC, p4_t_PETSC, p4_m_PETSC, p4_b_PETSC, help_PETSC, P_Store_Mapping )

USE PARAMETERS
Use Results
Use Visualizer
USE IFPORT
IMPLICIT DOUBLE PRECISION (A-H,O-Z)

! - Type DECLERATIONS -------------------------------------------------------------------------------------------------------------------------------
Integer, Intent(In), Dimension ( NNDH * NDim )              :: EqDis
Integer, Intent(In), Dimension ( NNDH * NDim )              :: EqDis_MAP_U_Store_Numbers_Global
Integer, Intent(In), Dimension ( NStore_Mapping )           ::           U_Store_Numbers_Global
Real(8), Intent(In), Dimension ( NStore_Mapping , 0:NStep ) ::    U_Store_Mapping
Real(8), Intent(Out),Dimension ( NStore_Mapping , 0:NStep ) ::    P_Store_Mapping
Real(8), Intent(In), Dimension ( NStore_Mapping , 0:NStep ) :: k1_m_Store_Mapping
Real(8), Intent(In), Dimension ( NStore_Mapping , 0:NStep ) :: k2_m_Store_Mapping
Real(8), Intent(In), Dimension ( NStore_Mapping , 0:NStep ) :: k3_m_Store_Mapping
!==================================================================================================
! - PETSC INCLUDE FILES
#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"
#include "finclude/petscis.h"

! - PETSC VARIABLES AND MATRICES
Mat            :: AK_PETSC, AC_PETSC, AG_PETSC
Vec            :: B_mis_PETSC
Vec            :: DIAG_M_PETSC, DIAG_iM_PETSC

PetscErrorCode :: IERR, ErrPTC
PetscMPIInt    :: SIZE, RANK

! - PETSc local vectors
Vec            ::  y_t_PETSC,  y_m_PETSC,  y_b_PETSC
Vec            :: p1_t_PETSC, p1_m_PETSC, p1_b_PETSC
Vec            :: p2_t_PETSC, p2_m_PETSC, p2_b_PETSC
Vec            :: p3_t_PETSC, p3_m_PETSC, p3_b_PETSC
Vec            :: p4_t_PETSC, p4_m_PETSC, p4_b_PETSC
Vec            :: help_PETSC

PetscScalar    :: u_misfit_PETSC ( NNDH * NDim )

! - Control Problem ----------------------------------------------------------------------------------------*-----------------------------------------
Vec            :: Ux_big_PETSc,         Uk_big_PETSc
Vec            :: Ux_Bigger_Rank_PETSc, Uk_Bigger_Rank_PETSc
Vec            ::  P_Bigger_Rank_PETSc
Vec            :: g_Lambda_PETSC, g_Mu_PETSC

PetscScalar    :: u_small_PETSC ( NStore_Mapping )

PetscScalar, pointer :: U_Bigger_Rank(:), P_Bigger_Rank(:)

VecScatter     :: vscat_u_hist

PetscScalar, pointer :: Yt(:), Ym(:), Yb(:), P(:)
!==================================================================================================


! in
!---------- ---------- ---------- ---------- ----------
DIMENSION XYZ ( NJ, NDIM ), INOD ( NNODE, NEL ), NGP ( NEL ), MTEL ( NEL ), ID ( NJ, NDOF )
!---------- ---------- ---------- ---------- ----------


! in
!---------- ---------- ---------- ---------- ----------
Dimension Dis_meas ( 0:NStep , NNDH * NDim )               ! Measured displacements at select sensor locations
Dimension idx_Mat_to ( NJ_Rank_IParts )                    ! Local equation numbers for materials, C-indexing
!---------- ---------- ---------- ---------- ----------


! local
!---------- ---------- ---------- ---------- ----------
Dimension u_misfit ( NNDH * NDim )
Dimension u_small  ( NStore_Mapping )
Dimension p_small  ( NStore_Mapping )
!---------- ---------- ---------- ---------- ----------


Call VecGetOwnershipRange ( y_b_PETSC, IStart, IEnd, IERR )


I_PRINT_STEP = 1
IF ( RANK == 0 .AND. I_PRINT_SCREEN == 1 ) WRITE(*,*) '3D RK4 Adjoint'
!==================================================================================================
! - DEFINE PETSC OBJECTS
!CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, y_t_PETSC, IERR )
!CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, y_m_PETSC, IERR )
!CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, y_b_PETSC, IERR )

!CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, p1_t_PETSC, IERR )
!CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, p1_m_PETSC, IERR )
!CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, p1_b_PETSC, IERR )

!CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, p2_t_PETSC, IERR )
!CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, p2_m_PETSC, IERR )
!CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, p2_b_PETSC, IERR )

!CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, p3_t_PETSC, IERR )
!CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, p3_m_PETSC, IERR )
!CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, p3_b_PETSC, IERR )

!CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, p4_t_PETSC, IERR )
!CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, p4_m_PETSC, IERR )
!CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, p4_b_PETSC, IERR )

!CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, help_PETSC, IERR )

CALL VecSet ( y_t_PETSC, 0.0d0, ierr)
CALL VecSet ( y_m_PETSC, 0.0d0, ierr)
CALL VecSet ( y_b_PETSC, 0.0d0, ierr)
!==================================================================================================

! Comment: System matrices K, C, G have been left-scaled with Mass inverse, and also by a minus sign. See Forward solver.
! form system matrices (cancel left-scaling, do right-scaling)
!---------- ---------- ---------- ---------- ----------
CALL MatDiagonalScale ( AK_PETSC, DIAG_M_PETSC, DIAG_iM_PETSC, IERR )
CALL MatDiagonalScale ( AC_PETSC, DIAG_M_PETSC, DIAG_iM_PETSC, IERR )
CALL MatDiagonalScale ( AG_PETSC, DIAG_M_PETSC, DIAG_iM_PETSC, IERR )


! Explicit Integration ------------------------------------------------------------------------------------------------------------------------------
DO ISTEP = NSTEP , 1 , -1

! -------------------------------------------------------------------------------------------------
! this is to see if adjoint diverges or not
! CALL TOTAL_ENERGY_RD_PETSC_EXPLICIT ( AK_RD_PETSC, DIAG_M_RD_PETSC, y_m_PETSC, y_b_PETSC, help_PETSC, ENERGY )
! if (rank ==0 ) write(*,*) energy
! -------------------------------------------------------------------------------------------------


! storage of the adjoint variable
!---------- ---------- ---------- ---------- ---------- ---------------------------------------------------------------------------------------------
! 1- selection: during time integration, access p(big), pick p(small) entries.
   Call VecGetArrayF90 ( y_t_PETSC, P, IERR )
   Do I = 1, NStore_Mapping
     p_small ( I ) = P ( U_Store_Numbers_Global ( I ) - IStart + 1 )
   End Do
   Call VecRestoreArrayF90 ( y_t_PETSC, P, IERR )

! 2- store p(small) entries in the dense matrix: We are storing 0:N-1.
   P_Store_Mapping ( : , ISTEP ) = p_small (:)
!---------- ---------- ---------- ---------- ---------- ---------------------------------------------------------------------------------------------


! 1- screen printout for nonbelivers; may they be enlightened
   T = DBLE( ISTEP ) * DT
   IF ( RANK == 0 .AND. MOD ( ISTEP , I_PRINT_STEP ) == 0 ) THEN
      WRITE (*,'(A30,F10.5,a30,e30.16)') 'TIME=', T !, 'Energy =', energy
   END IF

! 3- Last time-step requires special care
   If ( ISTEP == NSTEP ) Then

! 3.1- form the misfit (only for the last time step: NStep)
      Do IEq = 1, NNDH * NDim
        u_misfit ( IEq ) = U_Store_Mapping ( EqDis_MAP_U_Store_Numbers_Global ( IEq ) , IStep ) - Dis_meas ( IStep , IEq )
      End Do

! 3.2- initialize the misfit vector
      Call VecSet ( B_mis_PETSC, 0.0d0, IERR )

! 3.3- assemble the distributed misfit vector
      u_misfit_PETSC = u_misfit
      If ( NNDH /= 0 ) CALL VecSetValues ( B_mis_PETSC, NNDH * NDim, EqDis, u_misfit_PETSC, ADD_VALUES, IERR ) ! insert_values should be the same here.
      Call VecAssemblyBegin ( B_mis_PETSC, IERR )
      Call VecAssemblyEnd   ( B_mis_PETSC, IERR )

      Call VecSet  ( y_m_PETSC, 0.0d0, IERR )
      Call VecAXPY ( y_m_PETSC, DT, B_mis_PETSC, IERR )

! 3.4- Q is scaled with the inverse of the mass matrix
      CALL VecPointwiseMult ( B_mis_PETSC, B_mis_PETSC, DIAG_iM_PETSC, IERR )

   End If

! 4- construction of the misfit at select sensor locations (all time-steps except NStep)
! comment: this does not need an if-statement, since we are going to update y_m^{n-1}, which needs misfit^{n-1}
   Do IEq = 1, NNDH * NDim
     u_misfit ( IEq ) = U_Store_Mapping ( EqDis_MAP_U_Store_Numbers_Global ( IEq ) , IStep-1 ) - Dis_meas ( IStep-1 , IEq )
   End Do

! 5.1- initialize the misfit vector
   Call VecSet ( B_mis_PETSC, 0.0d0, IERR )

! 5.2- assemble the distributed misfit vector
   u_misfit_PETSC = u_misfit
   If ( NNDH /= 0 ) CALL VecSetValues ( B_mis_PETSC, NNDH * NDim, EqDis, u_misfit_PETSC, ADD_VALUES, IERR )
   Call VecAssemblyBegin ( B_mis_PETSC, IERR )
   Call VecAssemblyEnd   ( B_mis_PETSC, IERR )

! 5.3- Q is scaled with the inverse of the mass matrix
   CALL VecPointwiseMult ( B_mis_PETSC, B_mis_PETSC, DIAG_iM_PETSC, IERR )

! 6- middle steps: compute p4_t, p4_m, p4_b
   Call VecSet  ( p4_t_PETSC, 0.0d0, IERR )
   Call VecAXPY ( p4_t_PETSC, DT/6.0d0, y_t_PETSC, IERR )
   Call VecSet  ( p4_m_PETSC, 0.0d0, IERR )
   Call VecAXPY ( p4_m_PETSC, DT/6.0d0, y_m_PETSC, IERR )
   Call VecSet  ( p4_b_PETSC, 0.0d0, IERR )
   Call VecAXPY ( p4_b_PETSC, DT/6.0d0, y_b_PETSC, IERR )

GO TO 7718
! ggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggg
! -------------------------------------------------------------------------------------------------
! a- retrieve from the dense matrix: X^{n-1}
   u_small (:) = U_Store_Mapping ( : , ISTEP-1 )
! b- initlalize the parallel vector
   Call VecSet ( Ux_big_PETSc, 0.0d0, IERR )
! c- move from the dense matrix to parallel vector (long-thin vector)
   u_small_PETSC = u_small
   If ( NStore_Mapping /= 0 ) Call VecSetValues ( Ux_big_PETSc, NStore_Mapping, U_Store_Numbers_Global, u_small_PETSC, ADD_VALUES, IERR ) ! INSERT_VALUES should give the same result.
! d- this is not needed since everything is on the right processor, but let's do it.
   Call VecAssemblyBegin ( Ux_big_PETSc, IERR )
   Call VecAssemblyEnd   ( Ux_big_PETSc, IERR )
! e- scatter state variables
   Call VecScatterBegin ( vscat_u_hist, Ux_big_PETSc, Ux_Bigger_Rank_PETSc, INSERT_VALUES, SCATTER_FORWARD, IERR )
   Call VecScatterEnd   ( vscat_u_hist, Ux_big_PETSc, Ux_Bigger_Rank_PETSc, INSERT_VALUES, SCATTER_FORWARD, IERR )
! f- scatter adjoint variables
   Call VecScatterBegin ( vscat_u_hist,   p4_m_PETSC,  P_Bigger_Rank_PETSc, INSERT_VALUES, SCATTER_FORWARD, IERR )
   Call VecScatterEnd   ( vscat_u_hist,   p4_m_PETSC,  P_Bigger_Rank_PETSc, INSERT_VALUES, SCATTER_FORWARD, IERR )
! f- access data for assembly
   Call VecGetArrayF90 ( Ux_Bigger_Rank_PETSc, U_Bigger_Rank, IERR )
   Call VecGetArrayF90 (  P_Bigger_Rank_PETSc, P_Bigger_Rank, IERR )
! ... update the gradients ...
      Call Update_g ( g_Lambda_PETSC, g_Mu_PETSC, U_Bigger_Rank, P_Bigger_Rank, idx_Mat_to, NJ_Rank_IParts, XYZ, INOD, NGP, MTEL, ID, 1.0d0 )
! ... restore arrays
   Call VecRestoreArrayF90 (  P_Bigger_Rank_PETSc, P_Bigger_Rank, IERR )
   Call VecRestoreArrayF90 ( Ux_Bigger_Rank_PETSc, U_Bigger_Rank, IERR )
! -------------------------------------------------------------------------------------------------
! ggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggg
! ggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggg
! -------------------------------------------------------------------------------------------------
! a- retrieve from the dense matrix
   u_small (:) = k3_m_Store_Mapping ( : , ISTEP )
! b- initlalize the parallel vector
   Call VecSet ( Uk_big_PETSc, 0.0d0, IERR )
! c- move from the dense matrix to parallel vector
   u_small_PETSC = u_small
   If ( NStore_Mapping /= 0 ) Call VecSetValues ( Uk_big_PETSc, NStore_Mapping, U_Store_Numbers_Global, u_small_PETSC, ADD_VALUES, IERR ) ! one-to-one map: ADD = INSERT
! d- this is not needed since everything is on the right processor, but let's do it.
   Call VecAssemblyBegin ( Uk_big_PETSc, IERR )
   Call VecAssemblyEnd   ( Uk_big_PETSc, IERR )
! e- scatter state variables
   Call VecScatterBegin ( vscat_u_hist, Uk_big_PETSc, Uk_Bigger_Rank_PETSc, INSERT_VALUES, SCATTER_FORWARD, IERR )
   Call VecScatterEnd   ( vscat_u_hist, Uk_big_PETSc, Uk_Bigger_Rank_PETSc, INSERT_VALUES, SCATTER_FORWARD, IERR )
! f- scatter adjoint variables
!   Call VecScatterBegin ( vscat_u_hist, p4_m_PETSC, P_Bigger_Rank_PETSc, INSERT_VALUES, SCATTER_FORWARD, IERR )
!   Call VecScatterEnd   ( vscat_u_hist, p4_m_PETSC, P_Bigger_Rank_PETSc, INSERT_VALUES, SCATTER_FORWARD, IERR )
! f- access data for assembly
   Call VecGetArrayF90 ( Uk_Bigger_Rank_PETSc, U_Bigger_Rank, IERR )
   Call VecGetArrayF90 (  P_Bigger_Rank_PETSc, P_Bigger_Rank, IERR )
! ... update the gradients ...
      Call Update_g ( g_Lambda_PETSC, g_Mu_PETSC, U_Bigger_Rank, P_Bigger_Rank, idx_Mat_to, NJ_Rank_IParts, XYZ, INOD, NGP, MTEL, ID, DT )
! ... restore arrays
   Call VecRestoreArrayF90 (  P_Bigger_Rank_PETSc, P_Bigger_Rank, IERR )
   Call VecRestoreArrayF90 ( Uk_Bigger_Rank_PETSc, U_Bigger_Rank, IERR )
! -------------------------------------------------------------------------------------------------
! ggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggg
7718 CONTINUE

! 7- middle steps: compute p3_t, p3_m, p3_b
   Call VecSet   ( p3_t_PETSC, 0.0d0, IERR )
   Call VecMAXPY ( p3_t_PETSC, 2, [DT, DT/3.0d0], [p4_m_PETSC, y_t_PETSC], IERR )
   Call VecSet   ( p3_m_PETSC, 0.0d0, IERR )
   Call VecMAXPY ( p3_m_PETSC, 2, [DT, DT/3.0d0], [p4_b_PETSC, y_m_PETSC], IERR )
   Call VecSet   ( p3_b_PETSC, 0.0d0, IERR )
   Call MatMultTransposeAdd ( AG_PETSC, p4_t_PETSC, p3_b_PETSC, p3_b_PETSC, IERR )
   Call MatMultTransposeAdd ( AK_PETSC, p4_m_PETSC, p3_b_PETSC, p3_b_PETSC, IERR )
   Call MatMultTransposeAdd ( AC_PETSC, p4_b_PETSC, p3_b_PETSC, p3_b_PETSC, IERR )
   Call VecScale ( p3_b_PETSC, DT, IERR )
   Call VecAXPY  ( p3_b_PETSC, DT/3.0d0, y_b_PETSC, IERR )

GO TO 7719
! ggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggg
! -------------------------------------------------------------------------------------------------
! a- retrieve from the dense matrix
!   u_small (:) = U_Store_Mapping ( : , ISTEP-1 )
! b- initlalize the parallel vector
!   Call VecSet ( Ux_big_PETSc, 0.0d0, IERR )
! c- move from the dense matrix to parallel vector
!   u_small_PETSC = u_small
!   Call VecSetValues ( Ux_big_PETSc, NStore_Mapping, U_Store_Numbers_Global, u_small_PETSC, INSERT_VALUES, IERR )
! d- this is not needed since everything is on the right processor, but let's do it.
!   Call VecAssemblyBegin ( Ux_big_PETSc, IERR )
!   Call VecAssemblyEnd   ( Ux_big_PETSc, IERR )
! e- scatter state variables
!   Call VecScatterBegin ( vscat_u_hist, Ux_big_PETSc, Ux_Bigger_Rank_PETSc, INSERT_VALUES, SCATTER_FORWARD, IERR )
!   Call VecScatterEnd   ( vscat_u_hist, Ux_big_PETSc, Ux_Bigger_Rank_PETSc, INSERT_VALUES, SCATTER_FORWARD, IERR )
! f- scatter adjoint variables
   Call VecScatterBegin ( vscat_u_hist, p3_m_PETSC, P_Bigger_Rank_PETSc, INSERT_VALUES, SCATTER_FORWARD, IERR )
   Call VecScatterEnd   ( vscat_u_hist, p3_m_PETSC, P_Bigger_Rank_PETSc, INSERT_VALUES, SCATTER_FORWARD, IERR )
! f- access data for assembly
   Call VecGetArrayF90 ( Ux_Bigger_Rank_PETSc, U_Bigger_Rank, IERR )
   Call VecGetArrayF90 (  P_Bigger_Rank_PETSc, P_Bigger_Rank, IERR )
! ... update the gradients ...
      Call Update_g ( g_Lambda_PETSC, g_Mu_PETSC, U_Bigger_Rank, P_Bigger_Rank, idx_Mat_to, NJ_Rank_IParts, XYZ, INOD, NGP, MTEL, ID, 1.0d0 )
! ... restore arrays
   Call VecRestoreArrayF90 (  P_Bigger_Rank_PETSc, P_Bigger_Rank, IERR )
   Call VecRestoreArrayF90 ( Ux_Bigger_Rank_PETSc, U_Bigger_Rank, IERR )
! -------------------------------------------------------------------------------------------------
! ggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggg
! ggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggg
! -------------------------------------------------------------------------------------------------
! a- retrieve from the dense matrix
   u_small (:) = k2_m_Store_Mapping ( : , ISTEP )
! b- initlalize the parallel vector
   Call VecSet ( Uk_big_PETSc, 0.0d0, IERR )
! c- move from the dense matrix to parallel vector
   u_small_PETSC = u_small
   If ( NStore_Mapping /= 0 ) Call VecSetValues ( Uk_big_PETSc, NStore_Mapping, U_Store_Numbers_Global, u_small_PETSC, INSERT_VALUES, IERR )
! d- this is not needed since everything is on the right processor, but let's do it.
   Call VecAssemblyBegin ( Uk_big_PETSc, IERR )
   Call VecAssemblyEnd   ( Uk_big_PETSc, IERR )
! e- scatter state variables
   Call VecScatterBegin ( vscat_u_hist, Uk_big_PETSc, Uk_Bigger_Rank_PETSc, INSERT_VALUES, SCATTER_FORWARD, IERR )
   Call VecScatterEnd   ( vscat_u_hist, Uk_big_PETSc, Uk_Bigger_Rank_PETSc, INSERT_VALUES, SCATTER_FORWARD, IERR )
! f- scatter adjoint variables
!   Call VecScatterBegin ( vscat_u_hist, p3_m_PETSC, P_Bigger_Rank_PETSc, INSERT_VALUES, SCATTER_FORWARD, IERR )
!   Call VecScatterEnd   ( vscat_u_hist, p3_m_PETSC, P_Bigger_Rank_PETSc, INSERT_VALUES, SCATTER_FORWARD, IERR )
! f- access data for assembly
   Call VecGetArrayF90 ( Uk_Bigger_Rank_PETSc, U_Bigger_Rank, IERR )
   Call VecGetArrayF90 (  P_Bigger_Rank_PETSc, P_Bigger_Rank, IERR )
! ... update the gradients ...
      Call Update_g ( g_Lambda_PETSC, g_Mu_PETSC, U_Bigger_Rank, P_Bigger_Rank, idx_Mat_to, NJ_Rank_IParts, XYZ, INOD, NGP, MTEL, ID, DT/2.0d0 )
! ... restore arrays
   Call VecRestoreArrayF90 (  P_Bigger_Rank_PETSc, P_Bigger_Rank, IERR )
   Call VecRestoreArrayF90 ( Uk_Bigger_Rank_PETSc, U_Bigger_Rank, IERR )
! -------------------------------------------------------------------------------------------------
! ggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggg
7719 CONTINUE

! 8- middle steps: compute p2_t, p2_m, p2_b
   Call VecSet   ( p2_t_PETSC, 0.0d0, IERR )
   Call VecMAXPY ( p2_t_PETSC, 2, [DT/2.0d0, DT/3.0d0], [p3_m_PETSC, y_t_PETSC], IERR )
   Call VecSet   ( p2_m_PETSC, 0.0d0, IERR )
   Call VecMAXPY ( p2_m_PETSC, 2, [DT/2.0d0, DT/3.0d0], [p3_b_PETSC, y_m_PETSC], IERR )
   Call VecSet   ( p2_b_PETSC, 0.0d0, IERR )
   Call MatMultTransposeAdd ( AG_PETSC, p3_t_PETSC, p2_b_PETSC, p2_b_PETSC, IERR )
   Call MatMultTransposeAdd ( AK_PETSC, p3_m_PETSC, p2_b_PETSC, p2_b_PETSC, IERR )
   Call MatMultTransposeAdd ( AC_PETSC, p3_b_PETSC, p2_b_PETSC, p2_b_PETSC, IERR )
   Call VecScale ( p2_b_PETSC, DT/2.0d0, IERR )
   Call VecAXPY  ( p2_b_PETSC, DT/3.0d0, y_b_PETSC, IERR )

GO TO 7720
! ggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggg
! -------------------------------------------------------------------------------------------------
! a- retrieve from the dense matrix
!   u_small (:) = U_Store_Mapping ( : , ISTEP-1 )
! b- initlalize the parallel vector
!   Call VecSet ( Ux_big_PETSc, 0.0d0, IERR )
! c- move from the dense matrix to parallel vector
!   u_small_PETSC = u_small
!   Call VecSetValues ( Ux_big_PETSc, NStore_Mapping, U_Store_Numbers_Global, u_small_PETSC, INSERT_VALUES, IERR )
! d- this is not needed since everything is on the right processor, but let's do it.
!   Call VecAssemblyBegin ( Ux_big_PETSc, IERR )
!   Call VecAssemblyEnd   ( Ux_big_PETSc, IERR )
! e- scatter state variables
!   Call VecScatterBegin ( vscat_u_hist, Ux_big_PETSc, Ux_Bigger_Rank_PETSc, INSERT_VALUES, SCATTER_FORWARD, IERR )
!   Call VecScatterEnd   ( vscat_u_hist, Ux_big_PETSc, Ux_Bigger_Rank_PETSc, INSERT_VALUES, SCATTER_FORWARD, IERR )
! f- scatter adjoint variables
   Call VecScatterBegin ( vscat_u_hist, p2_m_PETSC, P_Bigger_Rank_PETSc, INSERT_VALUES, SCATTER_FORWARD, IERR )
   Call VecScatterEnd   ( vscat_u_hist, p2_m_PETSC, P_Bigger_Rank_PETSc, INSERT_VALUES, SCATTER_FORWARD, IERR )
! f- access data for assembly
   Call VecGetArrayF90 ( Ux_Bigger_Rank_PETSc, U_Bigger_Rank, IERR )
   Call VecGetArrayF90 (  P_Bigger_Rank_PETSc, P_Bigger_Rank, IERR )
! ... update the gradients ...
      Call Update_g ( g_Lambda_PETSC, g_Mu_PETSC, U_Bigger_Rank, P_Bigger_Rank, idx_Mat_to, NJ_Rank_IParts, XYZ, INOD, NGP, MTEL, ID, 1.0d0 )
! ... restore arrays
   Call VecRestoreArrayF90 (  P_Bigger_Rank_PETSc, P_Bigger_Rank, IERR )
   Call VecRestoreArrayF90 ( Ux_Bigger_Rank_PETSc, U_Bigger_Rank, IERR )
! -------------------------------------------------------------------------------------------------
! ggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggg
! ggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggg
! -------------------------------------------------------------------------------------------------
! a- retrieve from the dense matrix
   u_small (:) = k1_m_Store_Mapping ( : , ISTEP )
! b- initlalize the parallel vector
   Call VecSet ( Uk_big_PETSc, 0.0d0, IERR )
! c- move from the dense matrix to parallel vector
   u_small_PETSC = u_small
   If ( NStore_Mapping /= 0 ) Call VecSetValues ( Uk_big_PETSc, NStore_Mapping, U_Store_Numbers_Global, u_small_PETSC, INSERT_VALUES, IERR )
! d- this is not needed since everything is on the right processor, but let's do it.
   Call VecAssemblyBegin ( Uk_big_PETSc, IERR )
   Call VecAssemblyEnd   ( Uk_big_PETSc, IERR )
! e- scatter state variables
   Call VecScatterBegin ( vscat_u_hist, Uk_big_PETSc, Uk_Bigger_Rank_PETSc, INSERT_VALUES, SCATTER_FORWARD, IERR )
   Call VecScatterEnd   ( vscat_u_hist, Uk_big_PETSc, Uk_Bigger_Rank_PETSc, INSERT_VALUES, SCATTER_FORWARD, IERR )
! f- scatter adjoint variables
!   Call VecScatterBegin ( vscat_u_hist, p3_m_PETSC, P_Bigger_Rank_PETSc, INSERT_VALUES, SCATTER_FORWARD, IERR )
!   Call VecScatterEnd   ( vscat_u_hist, p3_m_PETSC, P_Bigger_Rank_PETSc, INSERT_VALUES, SCATTER_FORWARD, IERR )
! f- access data for assembly
   Call VecGetArrayF90 ( Uk_Bigger_Rank_PETSc, U_Bigger_Rank, IERR )
   Call VecGetArrayF90 (  P_Bigger_Rank_PETSc, P_Bigger_Rank, IERR )
! ... update the gradients ...
      Call Update_g ( g_Lambda_PETSC, g_Mu_PETSC, U_Bigger_Rank, P_Bigger_Rank, idx_Mat_to, NJ_Rank_IParts, XYZ, INOD, NGP, MTEL, ID, DT/2.0d0 )
! ... restore arrays
   Call VecRestoreArrayF90 (  P_Bigger_Rank_PETSc, P_Bigger_Rank, IERR )
   Call VecRestoreArrayF90 ( Uk_Bigger_Rank_PETSc, U_Bigger_Rank, IERR )
! -------------------------------------------------------------------------------------------------
! ggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggg
7720 CONTINUE

! 9- middle steps: compute p1_t, p1_m, p1_b
   Call VecSet   ( p1_t_PETSC, 0.0d0, IERR )
   Call VecMAXPY ( p1_t_PETSC, 2, [DT/2.0d0, DT/6.0d0], [p2_m_PETSC, y_t_PETSC], IERR )
   Call VecSet   ( p1_m_PETSC, 0.0d0, IERR )
   Call VecMAXPY ( p1_m_PETSC, 2, [DT/2.0d0, DT/6.0d0], [p2_b_PETSC, y_m_PETSC], IERR )
   Call VecSet ( p1_b_PETSC, 0.0d0, IERR )
   Call MatMultTransposeAdd ( AG_PETSC, p2_t_PETSC, p1_b_PETSC, p1_b_PETSC, IERR )
   Call MatMultTransposeAdd ( AK_PETSC, p2_m_PETSC, p1_b_PETSC, p1_b_PETSC, IERR )
   Call MatMultTransposeAdd ( AC_PETSC, p2_b_PETSC, p1_b_PETSC, p1_b_PETSC, IERR )
   Call VecScale ( p1_b_PETSC, DT/2.0d0, IERR )
   Call VecAXPY  ( p1_b_PETSC, DT/6.0d0, y_b_PETSC, IERR )

GO TO 7721
! ggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggg
! -------------------------------------------------------------------------------------------------
! a- retrieve from the dense matrix
!   u_small (:) = U_Store_Mapping ( : , ISTEP-1 )
! b- initlalize the parallel vector
!   Call VecSet ( Ux_big_PETSc, 0.0d0, IERR )
! c- move from the dense matrix to parallel vector
!   u_small_PETSC = u_small
!   Call VecSetValues ( Ux_big_PETSc, NStore_Mapping, U_Store_Numbers_Global, u_small_PETSC, INSERT_VALUES, IERR )
! d- this is not needed since everything is on the right processor, but let's do it.
!   Call VecAssemblyBegin ( Ux_big_PETSc, IERR )
!   Call VecAssemblyEnd   ( Ux_big_PETSc, IERR )
! e- scatter state variables
!   Call VecScatterBegin ( vscat_u_hist, Ux_big_PETSc, Ux_Bigger_Rank_PETSc, INSERT_VALUES, SCATTER_FORWARD, IERR )
!   Call VecScatterEnd   ( vscat_u_hist, Ux_big_PETSc, Ux_Bigger_Rank_PETSc, INSERT_VALUES, SCATTER_FORWARD, IERR )
! f- scatter adjoint variables
   Call VecScatterBegin ( vscat_u_hist, p1_m_PETSC, P_Bigger_Rank_PETSc, INSERT_VALUES, SCATTER_FORWARD, IERR )
   Call VecScatterEnd   ( vscat_u_hist, p1_m_PETSC, P_Bigger_Rank_PETSc, INSERT_VALUES, SCATTER_FORWARD, IERR )
! f- access data for assembly
   Call VecGetArrayF90 ( Ux_Bigger_Rank_PETSc, U_Bigger_Rank, IERR )
   Call VecGetArrayF90 (  P_Bigger_Rank_PETSc, P_Bigger_Rank, IERR )
! ... update the gradients ...
      Call Update_g ( g_Lambda_PETSC, g_Mu_PETSC, U_Bigger_Rank, P_Bigger_Rank, idx_Mat_to, NJ_Rank_IParts, XYZ, INOD, NGP, MTEL, ID, 1.0d0 )
! ... restore arrays
   Call VecRestoreArrayF90 (  P_Bigger_Rank_PETSc, P_Bigger_Rank, IERR )
   Call VecRestoreArrayF90 ( Ux_Bigger_Rank_PETSc, U_Bigger_Rank, IERR )
! -------------------------------------------------------------------------------------------------
! ggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggg
7721 CONTINUE

! 10- middle steps: y_t^{n-1}
   ! compute h_m
   Call VecSet  ( help_PETSC, 0.0d0, IERR )
   Call VecMAXPY( help_PETSC, 4, [1.0d0, 1.0d0, 1.0d0, 1.0d0], [p1_m_PETSC, p2_m_PETSC, p3_m_PETSC, p4_m_PETSC], IERR )
   ! compute (update) y_t^{n-1}
   Call VecAXPY ( y_t_PETSC, 1.0d0, help_PETSC, IERR )
   ! update y_b^{n-1}
   Call MatMultTransposeAdd ( AK_PETSC, help_PETSC, y_b_PETSC, y_b_PETSC, IERR )

! 11- middle steps: y_m^{n-1}
   ! compute h_b
   Call VecSet  ( help_PETSC, 0.0d0, IERR )
   Call VecMAXPY( help_PETSC, 4, [1.0d0, 1.0d0, 1.0d0, 1.0d0], [p1_b_PETSC, p2_b_PETSC, p3_b_PETSC, p4_b_PETSC], IERR )
   ! compute (update) y_m^{n-1}
   Call VecMAXPY ( y_m_PETSC, 2, [1.0d0, DT], [help_PETSC, B_mis_PETSC], IERR )
   ! update y_b^{n-1}
   Call MatMultTransposeAdd ( AC_PETSC, help_PETSC, y_b_PETSC, y_b_PETSC, IERR )

! 12- middle steps: y_b^{n-1}
   ! compute h_t
   Call VecSet  ( help_PETSC, 0.0d0, IERR )
   Call VecMAXPY( help_PETSC, 4, [1.0d0, 1.0d0, 1.0d0, 1.0d0], [p1_t_PETSC, p2_t_PETSC, p3_t_PETSC, p4_t_PETSC], IERR )
   ! compute y_b^{n-1}
   Call MatMultTransposeAdd ( AG_PETSC, help_PETSC, y_b_PETSC, y_b_PETSC, IERR )






! =================================================================================================
! =================================================================================================
! Yt (D) = Ym (C). wtf is going on?
!   Call VecGetArrayF90 ( y_t_PETSC, Yt, IERR ) ; Call VecGetArrayF90 ( y_m_PETSC, Ym, IERR ) ; Call VecGetArrayF90 ( y_b_PETSC, Yb, IERR )
!      If ( rank == 12 ) then
!         Write ( 73 , '(3e20.10)' ) Yt ( 1 ) , Ym ( 1 ) , Yb ( 1 ) 
!      End If
!   Call VecRestoreArrayF90 ( y_t_PETSC, Yt, IERR ) ; Call VecRestoreArrayF90 ( y_m_PETSC, Ym, IERR ) ; Call VecRestoreArrayF90 ( y_b_PETSC, Yb, IERR )
! =================================================================================================
! =================================================================================================



!GO TO 7722
! ggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggg
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! -------------------------------------------------------------------------------------------------
! a- retrieve from the dense matrix: X^{n-1}
   u_small (:) = U_Store_Mapping ( : , ISTEP-1 )
!   u_small (:) = U_Store_Mapping ( : , ISTEP-1 )
! b- initlalize the parallel vector
   Call VecSet ( Ux_big_PETSc, 0.0d0, IERR )
! c- move from the dense matrix to parallel vector (long-thin vector)
   u_small_PETSC = u_small
   If ( NStore_Mapping /= 0 ) Call VecSetValues ( Ux_big_PETSc, NStore_Mapping, U_Store_Numbers_Global, u_small_PETSC, ADD_VALUES, IERR )
! d- this is not needed since everything is on the right processor, but let's do it.
   Call VecAssemblyBegin ( Ux_big_PETSc, IERR )
   Call VecAssemblyEnd   ( Ux_big_PETSc, IERR )
! e- scatter state variables
   Call VecScatterBegin ( vscat_u_hist, Ux_big_PETSc, Ux_Bigger_Rank_PETSc, INSERT_VALUES, SCATTER_FORWARD, IERR )
   Call VecScatterEnd   ( vscat_u_hist, Ux_big_PETSc, Ux_Bigger_Rank_PETSc, INSERT_VALUES, SCATTER_FORWARD, IERR )
! f- scatter adjoint variables
   Call VecScatterBegin ( vscat_u_hist,    y_t_PETSC,  P_Bigger_Rank_PETSc, INSERT_VALUES, SCATTER_FORWARD, IERR )
   Call VecScatterEnd   ( vscat_u_hist,    y_t_PETSC,  P_Bigger_Rank_PETSc, INSERT_VALUES, SCATTER_FORWARD, IERR )
! f- access data for assembly
   Call VecGetArrayF90 ( Ux_Bigger_Rank_PETSc, U_Bigger_Rank, IERR )
   Call VecGetArrayF90 (  P_Bigger_Rank_PETSc, P_Bigger_Rank, IERR )
! ... update the gradients ...
      Call Update_g ( g_Lambda_PETSC, g_Mu_PETSC, U_Bigger_Rank, P_Bigger_Rank, idx_Mat_to, NJ_Rank_IParts, XYZ, INOD, NGP, MTEL, ID, DT )
! ... restore arrays
   Call VecRestoreArrayF90 (  P_Bigger_Rank_PETSc, P_Bigger_Rank, IERR )
   Call VecRestoreArrayF90 ( Ux_Bigger_Rank_PETSc, U_Bigger_Rank, IERR )
! -------------------------------------------------------------------------------------------------
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! ggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggg
7722 CONTINUE




END DO
! End of Explicit Integration loop ------------------------------------------------------------------------------------------------------------------


!CALL VecDestroy ( y_t_PETSC, IERR )
!CALL VecDestroy ( y_m_PETSC, IERR )
!CALL VecDestroy ( y_b_PETSC, IERR )

!CALL VecDestroy ( p1_t_PETSC, IERR )
!CALL VecDestroy ( p1_m_PETSC, IERR )
!CALL VecDestroy ( p1_b_PETSC, IERR )

!CALL VecDestroy ( p2_t_PETSC, IERR )
!CALL VecDestroy ( p2_m_PETSC, IERR )
!CALL VecDestroy ( p2_b_PETSC, IERR )

!CALL VecDestroy ( p3_t_PETSC, IERR )
!CALL VecDestroy ( p3_m_PETSC, IERR )
!CALL VecDestroy ( p3_b_PETSC, IERR )

!CALL VecDestroy ( p4_t_PETSC, IERR )
!CALL VecDestroy ( p4_m_PETSC, IERR )
!CALL VecDestroy ( p4_b_PETSC, IERR )

!CALL VecDestroy ( help_PETSC, IERR )


Return ;

! =============================================== Close ERRORS ======================================================================================
1002  IF ( IO_File > 0 ) Then ;
        Write(*, Fmt_ERR1_Close) UnFile, IO_File  ;  Write(UnInf, Fmt_ERR1_Close) UnFile, IO_File  ; 
        Write(*, Fmt_FL) ; Write(UnInf, Fmt_FL) ;
        !Call BEEP_FAIL ;
        Write(*, Fmt_End) ; Read(*,*) ; ; STOP ;
      End If ;

! =============================================== ERROR IN Write STATEMENT ==========================================================================
1006    Write(*       , Fmt_Write1 ) UnFile, IO_Write ; Write( UnFile, Fmt_Write1 ) UnFile, IO_Write ;
        !#Call BEEP_FAIL ;
        Write(*, Fmt_FL) ;  Write(UnInf, Fmt_FL) ; Write(*, Fmt_End) ; Read(*,*) ; Return ;

END


!************************************************
! Compute the misfit cost
!************************************************
! REVISION : Monday 28 April 2014

SUBROUTINE Compute_Misfit ( B_mis_PETSC, NNDH, EqDis, U_Store_Mapping, NStore_Mapping, EqDis_MAP_U_Store_Numbers_Global, Dis_meas, Cost_Misfit )

USE PARAMETERS
USE IFPORT
IMPLICIT DOUBLE PRECISION (A-H,O-Z)

! - Type DECLERATIONS -------------------------------------------------------------------------------------------------------------------------------
Integer, Intent(In), Dimension ( NNDH * NDim )              :: EqDis
Integer, Intent(In), Dimension ( NNDH * NDim )              :: EqDis_MAP_U_Store_Numbers_Global
Real(8), Intent(In), Dimension ( NStore_Mapping , 0:NStep ) :: U_Store_Mapping
!==================================================================================================
! - PETSC INCLUDE FILES
#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

! - PETSC VARIABLES AND MATRICES
Vec            :: B_mis_PETSC

PetscErrorCode :: IERR, ErrPTC
PetscMPIInt    :: SIZE, RANK

PetscScalar    :: u_misfit_PETSC ( NNDH * NDim )
PetscScalar    :: Cost_Misfit_IStep
!==================================================================================================


! in
!---------- ---------- ---------- ---------- ----------
Dimension Dis_meas ( 0:NStep , NNDH * NDim )               ! Measured displacements at select sensor locations
!---------- ---------- ---------- ---------- ----------


! local
!---------- ---------- ---------- ---------- ----------
Dimension u_misfit ( NNDH * NDim )
!---------- ---------- ---------- ---------- ----------

Cost_Misfit  = 0.0d0


! Start Integration ------------------------------------------------------------------------------------------------------------------------------
DO ISTEP = NSTEP , 1 , -1

! 2- construction of the misfit at select sensor locations
   Do IEq = 1, NNDH * NDim
     u_misfit ( IEq ) = U_Store_Mapping ( EqDis_MAP_U_Store_Numbers_Global ( IEq ) , IStep ) - Dis_meas ( IStep , IEq )
   End Do

! 3- initialize the misfit vector
   Call VecSet ( B_mis_PETSC, 0.0d0, IERR )

! 4- assemble the distributed misfit vector
   u_misfit_PETSC = u_misfit
   If ( NNDH /= 0 ) CALL VecSetValues ( B_mis_PETSC, NNDH * NDim, EqDis, u_misfit_PETSC, ADD_VALUES, IERR )
   Call VecAssemblyBegin ( B_mis_PETSC, IERR ) ;  Call VecAssemblyEnd   ( B_mis_PETSC, IERR )

! 5- compute the misfit component of the discrete objective functional
   Call VecDot ( B_mis_PETSC, B_mis_PETSC, Cost_Misfit_IStep, IERR )
   Cost_Misfit = Cost_Misfit +  0.5d0 * Cost_Misfit_IStep * DT

END DO
! End of Integration loop ------------------------------------------------------------------------------------------------------------------


Return
END





