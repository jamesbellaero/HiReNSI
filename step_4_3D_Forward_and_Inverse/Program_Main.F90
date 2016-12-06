!---------- ---------- ---------- ---------- ----------
! Questions may be directed to Arash Fathi at arash.fathi@utexas.edu
!---------- ---------- ---------- ---------- ----------

! REVISION    : M 21 April  2014 - Material visualization.
! REVISION    : T 29 April  2014 - Armijo condition.
! REVISION    : W 23 July   2014 - Save material updates to resume inversion if interrupted.
! Comment     : M 28 July   2014 - successful reconstruction for single parameter inversion (Lambda / Mu).
! REVISION    : W 30 July   2014 - TV regularization.
! REVISION    : F 1  August 2014 - Biased search direction for Lambda.
! REVISION    : S 2  August 2014 - New control parameters: {Lambda + 2 Mu} , {Mu}
! REVISION    : T 9  Sep    2014 - Variable regularization scheme: w(z)
! REVISION    : T 30 Sep    2014 - LBFGS-Decoupled
! REVISION    : M 6  Oct    2014 - LBFGS-Coupled

!---------- ---------- ---------- ---------- ----------
!     I don't tell the murky world
!     to turn pure.
!     I purify myself
!     and check my reflection
!     in the water of the valley brook.
!---------- ---------- ---------- ---------- ----------

PROGRAM MAIN
USE PARAMETERS
Use Input_Subroutines
Use Input_Subroutines_HDF5
USE IFPORT
Use Results
Use Visualizer

IMPLICIT DOUBLE PRECISION (A-H,O-Z)

!==================================================================================================
! - PETSC INCLUDE FILES
#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"
#include "finclude/petscao.h"
#include "finclude/petscis.h"

! - PETSC VARIABLES AND MATRICES (Forward Problem)
AO                      :: AOOut ;               ! Application Ordering Object.
ISLocalToGlobalMapping  :: GMapping ;            ! Global Mapping (displacement and stress DOFs)

Mat                     :: AM_PETSC
Mat                     :: AC_PETSC
Mat                     :: AK_PETSC
Mat                     :: AG_PETSC
Mat                     :: AM_RD_PETSC
Mat                     :: AK_RD_PETSC

Vec                     :: B_PETSC
Vec                     :: RMJ_PETSC
Vec                     :: DIAG_M_PETSC, DIAG_iM_PETSC                                             ! Diagonal system mass matrix, and its inverse.
Vec                     :: DIAG_M_RD_PETSC

PetscErrorCode          :: IERR, ErrPTC
PetscMPIInt             :: SIZE, RANK

! - PETSC VARIABLES AND MATRICES (Inverse Problem)
ISLocalToGlobalMapping  :: GMapping_Mat ;                                                          ! Global Mapping (material nodes-inverse)
Vec                     :: Lambda_PETSC,       Mu_PETSC                                            ! Distributed material vectors
Vec                     :: Lambda_PETSC_pre,   Mu_PETSC_pre                                        ! Distributed material vectors (previous iteraton)
Vec                     :: Lambda_Rank_PETSC,  Mu_Rank_PETSC                                       ! Sequential material vectors on each rank, including ghost nodes for computation and matrix assembly
Vec                     :: g_Lambda_PETSC,     g_Mu_PETSC                                          ! gradient vector
Vec                     :: g_Lambda_PETSC_pre, g_Mu_PETSC_pre                                      ! gradient vector (previous iteraton)
Vec                     :: p_Lambda_PETSC,     p_Mu_PETSC                                          ! search direction
Vec                     :: p_Lambda_PETSC_pre, p_Mu_PETSC_pre                                      ! search direction (previous iteraton)
IS                      :: IS_Mat_from, IS_Mat_to                                                  ! material index set (from Global vector to Local sequential vectors)
VecScatter              :: vscat_Mat                                                               ! material scatter (for the entire problem)
PetscScalar, pointer    :: Lambda_Rank(:), Mu_Rank(:)                                              ! pointer for sequential material vectors Lambda_Rank_PETSC, Mu_Rank_PETSC

IS                      :: is_mat_extend_from, is_mat_extend_to                                    ! extension of material properties from the interface into the PML
IS                      :: IS_u_from, IS_u_to                                                      ! index set for displacements (for computing the control problem)
VecScatter              :: vscat_mat_extend                                                        ! material scatterer (into the PML)

Vec                     :: B_mis_PETSC                                                             ! misfit vector for the adjoint problem
Vec                     :: Ux_big_PETSc, Uk_big_PETSc                                              ! long-thin vector containing solution only within the regular domain
Vec                     :: Ux_Bigger_Rank_PETSc, Uk_Bigger_Rank_PETSc, P_Bigger_Rank_PETSc
VecScatter              :: vscat_u_hist                                                            ! scatter displacement history (for computing the control problem)
! Control problem:
Mat                     :: Reg_PETSC                                                               ! Laplacian operator (H1-regularization)
Vec                     :: iMg_PETSC                                                               ! inverse of the gradient mass matrix
Vec                     :: Reg_Lambda_PETSC, Reg_Mu_PETSC                                          ! regularization vectors
PetscScalar             :: Cost_Regularization_Lambda, Cost_Regularization_Mu
PetscScalar, pointer    :: iMg_Rank(:)

! Material Visualization
Vec                     :: Lambda_Rank_Vis_PETSC, Mu_Rank_Vis_PETSC                                ! Material vectors for visualization (20-noded sampling)
IS                      :: IS_Mat_Vis_from, IS_Mat_Vis_to                                          ! index sets for material visualization
VecScatter              :: vscat_Mat_Vis                                                           ! scatter for material visualization
PetscScalar, pointer    :: Lambda_Rank_Vis(:), Mu_Rank_Vis(:)                                      ! pointer to Lambda_Rank_Vis_PETSC, Mu_Rank_Vis_PETSC, to transfer to Fortran vectors, and later to HDF5

! Finite Differencing for checking the gradient
Vec                     :: g_Lambda_Rank_PETSC
PetscScalar, pointer    :: g_Lambda_Rank(:), g_Mu_Rank(:)

Vec                     ::  y_t_PETSC,  y_m_PETSC,  y_b_PETSC
Vec                     :: p1_t_PETSC, p1_m_PETSC, p1_b_PETSC
Vec                     :: p2_t_PETSC, p2_m_PETSC, p2_b_PETSC
Vec                     :: p3_t_PETSC, p3_m_PETSC, p3_b_PETSC
Vec                     :: p4_t_PETSC, p4_m_PETSC, p4_b_PETSC
Vec                     :: help_PETSC

! L-BFGS
Vec                     :: y_Lambda_LBFGS_PETSC, y_Mu_LBFGS_PETSC
Vec                     :: s_Lambda_LBFGS_PETSC, s_Mu_LBFGS_PETSC

PetscScalar, pointer    :: y_Lambda_LBFGS_Rank(:)
PetscScalar, pointer    ::     y_Mu_LBFGS_Rank(:)
PetscScalar, pointer    :: s_Lambda_LBFGS_Rank(:)
PetscScalar, pointer    ::     s_Mu_LBFGS_Rank(:)

!Vec                     :: p_Lambda_Rank_PETSC, p_Mu_Rank_PETSC
PetscScalar, pointer    :: p_Lambda_Rank(:), p_Mu_Rank(:)

!==================================================================================================
! variable declarations
Type ( BasicParam )                  :: Param
Integer                              :: ERR_Alloc, ERR_DeAlloc ;                                   ! Allocating and DeAllocating errors
Integer                              :: Output_Type ;                                              ! Determines if the output files are formatted -0- or binary -1-
Integer                              :: SOLVER_Type ;                                              ! Solver Type, 1: Iterative solver - 2: super LU direct solver
Integer                              :: Int_Order ;                                                ! Order of Quadrature - INtegration

! for inversion:
Integer                              :: NJ_Rank_IParts                                             ! Number of material nodes on each rank, including ghost nodes. (=NJ)
Integer                              :: NJ_Mapping                                                 ! Number of material nodes on each rank, excluding ghost nodes (node ownership).
Integer                              :: NStore_Mapping                                             ! Number of equations on each rank that should be stored, excluding ghost nodes (ownership).
Integer                              :: NStore_Rank                                                ! Number of equations on each rank that should be stored, including ghost nodes (for computing).
Integer                              :: Armijo_Counter
Integer                              :: IJ                                                         ! Loop counter for debugging purposes
Integer                              :: UnFile_L, UnFile_M
! - Material Visualization --------------------------------------------------------------------------------------------------------------------------
Integer                              :: NJ_Para                                                    ! Number of material nodes to be visualized (based on 20-noded elements)
Character (Kind = 1, Len = 20)  :: IndexRankTemp;! A variable for storing the Rank number in Character format for adding at the end of the input file Name in the do loop 
Character (Kind = 1, Len = 22)  :: IndexIter ;   ! A variable for storing the Material Iteration number in Character format for adding at the end of the input file Name in the do loop 
Character (Kind = 1, Len = 2 )  :: Regularization_Method
Character (Kind = 1, Len = 9 )  :: Control_Parameter
Character (Kind = 1, Len = 3 )  :: Bias_Lambda_search
Character (Kind = 1, Len = 5 )  :: Search_Direction_Method

Integer, Allocatable, Dimension(:,:) :: ID ;
Integer, Allocatable, Dimension(:,:) :: INOD ;
Integer, Allocatable, Dimension(:)   :: NGP ;
Integer, Allocatable, Dimension(:)   :: MTEL ;
Integer, Allocatable, Dimension(:,:) :: ID_BC ;
Real(8), Allocatable, Dimension(:,:) :: BACL ;
Real(8), Allocatable, Dimension(:,:) :: XYZ_serendipity ;
Integer, Allocatable, Dimension(:,:) :: ID_serendipity ;
Real(8), Allocatable, Dimension(:,:) :: PML_Param ;
Integer, Allocatable, Dimension(:)   :: LTRANS ;
Real(8), Allocatable, Dimension(:,:) :: XYZ_Lagrange ;
Integer, Allocatable, Dimension(:,:) :: ID_Lagrange ;
Real(8), Allocatable, Dimension(:,:) :: XYZ ;
Real(8), Allocatable, Dimension(:,:) :: PMAT ;
! heterogeneous material properties:
Real(8), Allocatable, Dimension(:)   :: PMat_Lambda ;
Real(8), Allocatable, Dimension(:)   :: PMat_Mu ;
Real(8), Allocatable, Dimension(:)   :: PMat_Temp ;

! additional arrays from PIC
Integer, Allocatable, Dimension(:)   :: ELT ;                        ! Element Type of each Element
Integer, Allocatable, Dimension(:)   :: ELGR ;                       ! Holds Group number of each element.
Real(8), Allocatable, Dimension(:,:) :: PML_DIM ;                    ! PML Dimensions.
Integer, Allocatable, Dimension(:,:) :: IDBC ;                       ! Identification of Boundary Condition
Integer, Allocatable, Dimension(:)   :: JLoad ;                      ! Joint Load, Holds node numbers with concentrated force.
Real(8), Allocatable, Dimension(:,:) :: PLoad ;                      ! Holds concentrated forces on nodes (See JLoad).
Integer, Allocatable, Dimension(:)   :: NDAN ;                       ! Holds Node numbers in which history of Displacements are required based on Application Numbering
Integer, Allocatable, Dimension(:)   :: NVAN ;                       ! Holds Node numbers in which history of Velocity are required based on Application Numbering
Integer, Allocatable, Dimension(:)   :: NAAN ;                       ! Holds Node numbers in which history of Acceleration are required based on Application Numbering
Integer, Allocatable, Dimension(:)   :: NDPN ;                       ! Holds Node numbers in which history of Displacements are required based on PETSc Numbering
Integer, Allocatable, Dimension(:)   :: NVPN ;                       ! Holds Node numbers in which history of Velocity are required based on PETSc Numbering
Integer, Allocatable, Dimension(:)   :: NAPN ;                       ! Holds Node numbers in which history of Acceleration are required based on PETSc Numbering
Integer, Allocatable, Dimension(:)   :: EqDis ;                      ! Holds equation numbers of nodes in which history of Displacement is required based on global PETSc numbering
Integer, Allocatable, Dimension(:)   :: EqVel ;                      ! Holds equation numbers of nodes in which history of Velocity is required based on global PETSc numbering
Integer, Allocatable, Dimension(:)   :: EqAcc ;                      ! Holds equation numbers of nodes in which history of Acceleration is required based on global PETSc numbering
Real(8), Allocatable, Dimension(:,:) :: PBLD ;                       ! Properties of the Body force LoaD.
Real(8), Allocatable, Dimension(:,:) :: UDis ;                       ! Holds predefined Displacements in nodes (like settlement).
Integer, Allocatable, Dimension(:)   :: App_Numbers ;                ! Holds application equation numbers used in defining application ordering - This matrix must be of Shrt type.    !{{{{
Integer, Allocatable, Dimension(:)   :: PETSc_Numbers ;              ! Holds Petsc equation numbers used in defining application ordering - This matrix must be of Shrt type.    !{{{{
Integer, Allocatable, Dimension(:)   :: Indices ;                    ! Holds equivalent Global equatin numbers for local numbering - THis array must of Shrt type    !{{{{
Integer, Allocatable, Dimension(:)   :: D_NNZ_Stiff ;                ! Number of Non-Zero entries of Diagonal part of the stiffness matrix for each node
Integer, Allocatable, Dimension(:)   :: O_NNZ_Stiff ;                ! Number of Non-Zero entries of Off-Diagonal part of the stiffness matrix for each node
Integer, Allocatable, Dimension(:)   :: D_NNZ_Mass ;                 ! Number of Non-Zero entries of Diagonal part of the mass matrix for each node
Integer, Allocatable, Dimension(:)   :: O_NNZ_Mass ;                 ! Number of Non-Zero entries of Off-Diagonal part of the mass matrix for each node
Integer, Allocatable, Dimension(:)   :: D_NNZ_Damp ;                 ! Number of Non-Zero entries of Diagonal part of the damp matrix for each node
Integer, Allocatable, Dimension(:)   :: O_NNZ_Damp ;                 ! Number of Non-Zero entries of Off-Diagonal part of the damp matrix for each node
Integer, Allocatable, Dimension(:)   :: D_NNZ_G ;                    ! Number of Non-Zero entries of Diagonal part of the G matrix for each node
Integer, Allocatable, Dimension(:)   :: O_NNZ_G ;                    ! Number of Non-Zero entries of Off-Diagonal part of the G matrix for each node
Integer, Allocatable, Dimension(:)   :: NoBndry_DRM ;                ! Holds Node numbers on the DRM boundary
Integer, Allocatable, Dimension(:)   :: NoLayer_DRM ;                ! Holds Node numbers on the DRM layer
Real(8), Allocatable, Dimension(:)   :: InciWave ;                   ! Properties of the incident wave for DRM. i.e. 1: theta, 2: omega 3: amplitude(U)
Integer, Allocatable, Dimension(:)   :: ND_b ;                       ! Holds the equation number of DRM boundary nodes. This array cannot be a LNG, because it has a conflict with PETSc lib.
Integer, Allocatable, Dimension(:)   :: ND_e ;                       ! Holds the equation number of DRM layer nodes
Integer, Allocatable, Dimension(:)   :: STEP ;                       ! Holds steps in which full information in dynamic analysis is required.
! for response visualization:
Integer, Allocatable, Dimension(:)   :: idx_from ;                   ! Holds all equation numbers on each rank for FE.
Integer, Allocatable, Dimension(:,:) :: ID_Para ;                    ! Holds all equation numbers on each rank for Paraview.
! for material visualization:
Integer, Allocatable, Dimension(:)   :: idx_Mat_Vis_from
Integer, Allocatable, Dimension(:)   :: idx_Mat_Vis_to
! for inversion:
Integer, Allocatable, Dimension(:)   :: idx_Mat_from                 ! Holdes node numbers (PETSc Numbering) of material properties - Global. See illustration in Subroutine Inversion_DS.
Integer, Allocatable, Dimension(:)   :: idx_Mat_to                   ! Local node numbers.
Integer, Allocatable, Dimension(:)   :: idx_Mat_Extend_from          ! Index set from: extension of material property from the interface into the PML.
Integer, Allocatable, Dimension(:)   :: idx_Mat_Extend_to            ! Index set to: extension of material property from the interface into the PML.
Integer, Allocatable, Dimension(:)   :: U_Store_Numbers_Global       ! Equation numbers that their response should be stored for inversion - Global PETSc.
Integer, Allocatable, Dimension(:)   :: idx_u_from                   ! Index set from: Equation numbers needed on each rank for gradient computation. (from the global long-thin vector)
Integer, Allocatable, Dimension(:)   :: idx_u_to                     ! Index set to: Equation numbers needed on each rank for gradient computation. (to the local sequential vectors)
! storage of solution history for the forward problem:
Real(8), Allocatable, Dimension(:,:) ::    U_Store_Mapping           ! Dense matrix for storing displacement time-history of the forward problem
Real(8), Allocatable, Dimension(:,:) ::    P_Store_Mapping           ! Dense matrix for storing adjoint variable time-history of the adjoint problem
Real(8), Allocatable, Dimension(:,:) :: k1_m_Store_Mapping           ! Dense matrix for storing k1_m (RK-4) time-history of the forward problem
Real(8), Allocatable, Dimension(:,:) :: k2_m_Store_Mapping           ! Dense matrix for storing k2_m (RK-4) time-history of the forward problem
Real(8), Allocatable, Dimension(:,:) :: k3_m_Store_Mapping           ! Dense matrix for storing k3_m (RK-4) time-history of the forward problem
!Real(8), Allocatable, Dimension(:,:) :: k4_m_Store_Mapping          ! Dense matrix for storing k4_m (RK-4) time-history of the forward problem
Real(8), Allocatable, Dimension(:,:) :: Dis_meas                     ! Measured response at select sensor locations to form the misfit.

Integer, Allocatable, Dimension(:)   :: EqDis_MAP_U_Store_Numbers_Global  ! Mapping between these two arrays.

Integer, Allocatable, Dimension(:)   :: D_NNZ_Reg ;                  ! Number of Non-Zero entries of Diagonal part of the Regularization matrix for each node
Integer, Allocatable, Dimension(:)   :: O_NNZ_Reg ;                  ! Number of Non-Zero entries of Off-Diagonal part of the Regularization matrix for each node

! L-BFGS
Real(8), Allocatable, Dimension(:)   ::   Rhos_Lambda_LBFGS ;        !
Real(8), Allocatable, Dimension(:)   ::       Rhos_Mu_LBFGS ;        !
Real(8), Allocatable, Dimension(:)   ::       Rhos_LM_LBFGS ;        !
Real(8), Allocatable, Dimension(:)   :: Alphas_Lambda_LBFGS ;        !
Real(8), Allocatable, Dimension(:)   ::     Alphas_Mu_LBFGS ;        !
Real(8), Allocatable, Dimension(:)   ::     Alphas_LM_LBFGS ;        !
Real(8), Allocatable, Dimension(:)   :: Alphas_Lambda_LBFGS_Rank ;   !
Real(8), Allocatable, Dimension(:)   ::     Alphas_Mu_LBFGS_Rank ;   !
Real(8), Allocatable, Dimension(:)   ::     Alphas_LM_LBFGS_Rank ;   !
Integer, Allocatable, Dimension(:)   :: LBFGS_vector_sequence ;      !

Real(8), Allocatable, Dimension(:)   :: q_Lambda
Real(8), Allocatable, Dimension(:)   :: q_Mu
Real(8), Allocatable, Dimension(:)   :: r_Lambda
Real(8), Allocatable, Dimension(:)   :: r_Mu

Real(8), Allocatable, Dimension(:,:) :: Ys_Lambda_LBFGS              ! dense matrices, comprising the m vectors corresponding to previous iterations
Real(8), Allocatable, Dimension(:,:) ::     Ys_Mu_LBFGS              ! dense matrices, comprising the m vectors corresponding to previous iterations
Real(8), Allocatable, Dimension(:,:) :: Ss_Lambda_LBFGS              ! dense matrices, comprising the m vectors corresponding to previous iterations
Real(8), Allocatable, Dimension(:,:) ::     Ss_Mu_LBFGS              ! dense matrices, comprising the m vectors corresponding to previous iterations

Real(8), Allocatable, Dimension(:)   :: restriction_RD               ! restriction of a vector to the regular domain (RD = 1, PML = 0)
!==================================================================================================
! - START PETSC
CALL PetscInitialize ( PETSC_NULL_CHARACTER  , IERR )
CALL MPI_Comm_size   ( PETSC_COMM_WORLD, SIZE, IERR )
CALL MPI_Comm_rank   ( PETSC_COMM_WORLD, RANK, IERR )

write(*,*) 'size =', size, 'rank =', rank
!==================================================================================================
CALL CPU_TIME(t0)
! Read input files in parallel from PIC -------------------------------------------------------------------------------------------------------------
!Include 'IO_Parallel_hdf5.F90'
Include 'IO_Parallel_hdf5.F90'
write(*,*) 'hello'
Include 'Adjustments.F90'
Include 'IO_Parallel_Inversion.F90'
write(*,*) 'hello2'
Include 'PETSc_Global_Objects.F90'
Include 'PETSc_Global_Objects_Inversion.F90'
Include 'Control_Problem_Operators.F90'
write(*,*) 'hello3'


!!========================================================
!  Include 'Solve_Forward.F90'                        !===
!  Call MPI_Barrier ( PETSc_COMM_WORLD, IERR ) ; Stop !===
!!========================================================



! ===================================================================================================================================================
! Solve the Inverse Problem
! ===================================================================================================================================================

! ---------------------------------------------------------------------------------------------------------------------------------------------------
! I- Initialize -------------------------------------------------------------------------------------------------------------------------------------<><><><><><><><><><><><><><><><><><><><>
! ---------------------------------------------------------------------------------------------------------------------------------------------------
! 1- material properties in Fortran vectors (for solving the forward problem): also, already in place in PMat_Lambda & PMat_Mu
 write(*,*) 'abcde' 
 If ( Iter_begin == 1 ) Then
    ! initial guess
    Select Case (Control_Parameter)
       Case ('Lambda')
          PMat_Lambda = 80.01d0
          ! read Mu from file (target)
       Case ('Mu')
          PMat_Mu     = 80.02d0
          ! read Lambda from file (target)
       Case ('Lambda_Mu')
          ! simultaneous inversion
          PMat_Lambda = 80.01d0
          PMat_Mu     = 80.02d0
    End Select
  Else
    ! Resume from the last save
    Include'Resume_Inversion_Read_Material.F90'
  End If
write(*,*) 'abcdef'
! 1.1 We use mo = {Lambda + 2 Mu} instead of {Lambda} as the control variable. Moreover, we store {Lambda + 2 Mu} in variable {Lambda} to avoid creating a new variable. 
!     for visualization though, we output {Lambda}.
  If ( Lambda_variable == 'Lambda_2Mu' ) Then
    PMat_Lambda = PMat_Lambda + 2.0d0 * PMat_Mu
  End If
write(*,*) 'abc'
! 2- material properties in distributed PETSc vectors (for computing the gradients)

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

! 3- read (synthetic or field) measured data at select sensor locations: already done in IO_Parallel_Inversion.F90

! 4- Solve State
! 4.1- use actual {Lambda} and {Mu} for solving state
  If ( Lambda_variable == 'Lambda_2Mu' ) Then
    PMat_Lambda = PMat_Lambda - 2.0d0 * PMat_Mu
  End If
! 4.2- solve
  Include 'Solve_Forward.F90'
! 4.3- use {Lambda + 2 Mu} instead of {Lambd} for inversion.
  If ( Lambda_variable == 'Lambda_2Mu' ) Then
    PMat_Lambda = PMat_Lambda + 2.0d0 * PMat_Mu
  End If
write(*,*) 'abcd'
! 5- Compute misfit
  Call Compute_Misfit ( B_mis_PETSC, NNDH, EqDis, U_Store_Mapping, NStore_Mapping, EqDis_MAP_U_Store_Numbers_Global, Dis_meas, Cost_Misfit )

! 6- Printout initial misfit and cost
  If ( Iter_begin == 1 ) Then
    ! compute the cost functional due to the initial guess
    Iter = 0
    Cost_Regularization_Lambda = 0.0d0
    Cost_Regularization_Mu     = 0.0d0
    g_Lambda_Norm_2            = 0.0d0
    g_Mu_Norm_2                = 0.0d0
    Cost_Functional = Cost_Misfit + Cost_Regularization_Mu + Cost_Regularization_Lambda
    If ( Rank == 0 ) Then
      UnFile = UN_Inversion_Iterations ;
      Write ( *     ,'(A146)') '--------------------------------------------------------------------------------------------------------------------------------------------------'
      Write ( *     ,'(A6,7A20)' ) 'Iter', 'Cost_Functional', 'Misfit', 'Reg_Lambda', 'Reg_Mu', '||g_Lambda||', '||g_Mu||', 'Reg_Fac_Cont'
      Write ( UnFile,'(A6,7A20)' ) 'Iter', 'Cost_Functional', 'Misfit', 'Reg_Lambda', 'Reg_Mu', '||g_Lambda||', '||g_Mu||', 'Reg_Fac_Cont'
      Write ( *     ,'(A146)') '--------------------------------------------------------------------------------------------------------------------------------------------------'
      Write ( *     ,'(I6,6E20.10E3,F20.4)' ) Iter, Cost_Functional, Cost_Misfit, Cost_Regularization_Lambda, Cost_Regularization_Mu, g_Lambda_Norm_2, g_Mu_Norm_2, Regularization_factor_continuation
      Write ( UnFile,'(I6,6E20.10E3,F20.4)' ) Iter, Cost_Functional, Cost_Misfit, Cost_Regularization_Lambda, Cost_Regularization_Mu, g_Lambda_Norm_2, g_Mu_Norm_2, Regularization_factor_continuation
    End If
  Else
    ! resume inversion - no need to compute the misfit, but print out the header
    If ( Rank == 0 ) Then
      Write ( * ,'(A146)') '--------------------------------------------------------------------------------------------------------------------------------------------------'
      Write ( * ,'(A6,7A20)' ) 'Iter', 'Cost_Functional', 'Misfit', 'Reg_Lambda', 'Reg_Mu', '||g_Lambda||', '||g_Mu||', 'Reg_Fac_Cont'
      Write ( * ,'(A146)') '--------------------------------------------------------------------------------------------------------------------------------------------------'
    End If
  End If


! ===================================================================================================================================================
! *** iterate until convergence: *** 
! ===================================================================================================================================================
   Do Iter = Iter_begin , Iter_end

! ---------------------------------------------------------------------------------------------------------------------------------------------------
! II- Store old results -----------------------------------------------------------------------------------------------------------------------------
! ---------------------------------------------------------------------------------------------------------------------------------------------------
! 1- misfit
     Cost_Misfit_pre = Cost_Misfit

! 2- gradients, search directions, Lambda & Mu
     Call VecCopy ( g_Lambda_PETSC , g_Lambda_PETSC_pre, IERR )
     Call VecCopy (     g_Mu_PETSC ,     g_Mu_PETSC_pre, IERR )
     Call VecCopy ( p_Lambda_PETSC , p_Lambda_PETSC_pre, IERR )
     Call VecCopy (     p_Mu_PETSC ,     p_Mu_PETSC_pre, IERR )
! I moved these guys down, after computing LBFGS 's'.
!     Call VecCopy (   Lambda_PETSC ,   Lambda_PETSC_pre, IERR )
!     Call VecCopy (       Mu_PETSC ,       Mu_PETSC_pre, IERR )

! ---------------------------------------------------------------------------------------------------------------------------------------------------
! III- Solve the Adjoint problem, and compute the gradient on the fly -------------------------------------------------------------------------------
! ---------------------------------------------------------------------------------------------------------------------------------------------------
! 1- initialize the gradient vector
     Call VecSet ( g_Lambda_PETSC, 0.0d0, IERR )
     Call VecSet ( g_Mu_PETSC,     0.0d0, IERR )

! 2- Solve the Adjoint & Control Problems (Optimize-then-Discretize)
     Call EXPLICIT_RK4_PETSC_PML_3D_Adjoint_OtD ( AK_PETSC, DIAG_M_PETSC, DIAG_iM_PETSC, AC_PETSC, AG_PETSC, B_mis_PETSC, AK_RD_PETSC, DIAG_M_RD_PETSC, ID, NNDH, EqDis, U_Store_Mapping, P_Store_Mapping, NStore_Mapping, U_Store_Numbers_Global, EqDis_MAP_U_Store_Numbers_Global, Dis_meas )
     Call EXPLICIT_RK4_PETSC_PML_3D_Control_OtD ( U_Store_Mapping, P_Store_Mapping, NStore_Mapping, U_Store_Numbers_Global, Ux_big_PETSc, Ux_Bigger_Rank_PETSc, Uk_big_PETSc, P_Bigger_Rank_PETSc, vscat_u_hist, g_Lambda_PETSC, g_Mu_PETSC, idx_Mat_to, NJ_Rank_IParts, XYZ, INOD, NGP, MTEL, ID, RANK )

! 3- assemble the distributed gradient vectors (only the part associated with the misfit and only within the regular domain)
      Call VecAssemblyBegin ( g_Lambda_PETSC, IERR )  ;  Call VecAssemblyEnd   ( g_Lambda_PETSC, IERR )
      Call VecAssemblyBegin (     g_Mu_PETSC, IERR )  ;  Call VecAssemblyEnd   (     g_Mu_PETSC, IERR )

! 4- correct the sign of the gradient vector
      Call VecScale ( g_Lambda_PETSC, -1.0d0, IERR )
      Call VecScale (     g_Mu_PETSC, -1.0d0, IERR )

! ---------------------------------------------------------------------------------------------------------------------------------------------------
! IV- Regularization --------------------------------------------------------------------------------------------------------------------------------
! ---------------------------------------------------------------------------------------------------------------------------------------------------
! 1- compute regularization vectors -- needs to be scaled later
      Select Case (Regularization_Method)
         Case ('TN')
            Call MatMult ( Reg_PETSC, Lambda_PETSC, Reg_Lambda_PETSC, IERR )
            Call MatMult ( Reg_PETSC,     Mu_PETSC,     Reg_Mu_PETSC, IERR )
         Case ('TV')
            Call VecSet ( Reg_Lambda_PETSC, 0.0d0, IERR )
            Call VecSet ( Reg_Mu_PETSC,     0.0d0, IERR )
            Call Assem_TV_TN_Regularization ( Reg_Lambda_PETSC, Reg_Mu_PETSC, PMat_Lambda, PMat_Mu, idx_Mat_to, NJ_Rank_IParts, XYZ, INOD, NGP, MTEL, ID, cost_reg_Lambda_rank, cost_reg_Mu_rank )
            Call VecAssemblyBegin ( Reg_Lambda_PETSC, IERR )  ;  Call VecAssemblyEnd   ( Reg_Lambda_PETSC, IERR )
            Call VecAssemblyBegin (     Reg_Mu_PETSC, IERR )  ;  Call VecAssemblyEnd   (     Reg_Mu_PETSC, IERR )
      End Select

! 2- scale with regularization factor (adaptive regularization factor - See Sez p. 171)
   ! Lambda ------------------------
      Call VecNorm (   g_Lambda_PETSC, Norm_2, Fac_Misfit_Lambda_Norm_2, IERR )
      Call VecNorm ( Reg_Lambda_PETSC, Norm_2,    Fac_Reg_Lambda_Norm_2, IERR )
      Reg_Lambda_Factor = Regularization_factor_continuation * Fac_Misfit_Lambda_Norm_2 / Fac_Reg_Lambda_Norm_2
      Call VecScale ( Reg_Lambda_PETSC, Reg_Lambda_Factor, IERR )

   ! Mu ----------------------------
      Call VecNorm (   g_Mu_PETSC, Norm_2, Fac_Misfit_Mu_Norm_2, IERR )
      Call VecNorm ( Reg_Mu_PETSC, Norm_2,    Fac_Reg_Mu_Norm_2, IERR )
      Reg_Mu_Factor     = Regularization_factor_continuation * Fac_Misfit_Mu_Norm_2 / Fac_Reg_Mu_Norm_2
      Call VecScale ( Reg_Mu_PETSC, Reg_Mu_Factor, IERR )

! 3- add regularization to the misfit-component of the gradient
      ! Note: Do not add regularization for Iter = 1, because we start with a homogeneous guess for material distribution which results in zero regularization.
      !       Whatever we compute for Iter = 1 is due to numerical round-off error.
      !       If you are starting with a non-homogeneous initial guess, remove the if-statement.
      If ( Iter /= 1 ) Then
         Call VecAXPY ( g_Lambda_PETSC, +1.0d0, Reg_Lambda_PETSC, IERR )
         Call VecAXPY (     g_Mu_PETSC, +1.0d0,     Reg_Mu_PETSC, IERR )
      End If

! 4- apply the (inverse) mass matrix associated with the gradient vector
      Call VecPointwiseMult ( g_Lambda_PETSC, g_Lambda_PETSC, iMg_PETSC, IERR )
      Call VecPointwiseMult (     g_Mu_PETSC,     g_Mu_PETSC, iMg_PETSC, IERR )

! ---------------------------------------------------------------------------------------------------------------------------------------------------
! V- Cost of the objective functional ---------------------------------------------------------------------------------------------------------------
! ---------------------------------------------------------------------------------------------------------------------------------------------------
! 1- Compute misfit
      Call Compute_Misfit ( B_mis_PETSC, NNDH, EqDis, U_Store_Mapping, NStore_Mapping, EqDis_MAP_U_Store_Numbers_Global, Dis_meas, Cost_Misfit )

! 2- compute the cost of the regularization component of the discrete objective functional
      Select Case (Regularization_Method)
         Case ('TN')
            Call VecDot ( Reg_Lambda_PETSC, Lambda_PETSC, Cost_Regularization_Lambda, IERR )
            Call VecDot (     Reg_Mu_PETSC,     Mu_PETSC, Cost_Regularization_Mu,     IERR )
            Cost_Regularization_Lambda = 0.50d0 * Cost_Regularization_Lambda 
            Cost_Regularization_Mu     = 0.50d0 * Cost_Regularization_Mu
         Case ('TV')
            Call MPI_Reduce ( cost_reg_Lambda_rank, Cost_Regularization_Lambda, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, PETSC_COMM_WORLD, IERR )
            Call MPI_Reduce ( cost_reg_Mu_rank,     Cost_Regularization_Mu,     1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, PETSC_COMM_WORLD, IERR )
            Cost_Regularization_Lambda = Cost_Regularization_Lambda * Reg_Lambda_Factor
            Cost_Regularization_Mu     = Cost_Regularization_Mu     * Reg_Mu_Factor
      End Select

! 3- we probably have no regularization at the beginning
      If ( Iter == 1 ) Then
         Cost_Regularization_Lambda = 0.0d0
         Cost_Regularization_Mu     = 0.0d0
      End If

! 4- value of the objective functional
      Cost_Functional = Cost_Misfit + Cost_Regularization_Mu + Cost_Regularization_Lambda

! ---------------------------------------------------------------------------------------------------------------------------------------------------
! VI- Find a search direction -----------------------------------------------------------------------------------------------------------------------
! ---------------------------------------------------------------------------------------------------------------------------------------------------
! 1- initialize the search direction
     Call VecSet ( p_Lambda_PETSC, 0.0d0, IERR )
     Call VecSet ( p_Mu_PETSC,     0.0d0, IERR )

! 2- compute a search direction 
     Select Case ( Search_Direction_Method )
        Case ('SD') 
           ! using the steepest descent method
           Call VecAXPY ( p_Lambda_PETSC, -1.0d0, g_Lambda_PETSC, IERR )
           Call VecAXPY (     p_Mu_PETSC, -1.0d0,     g_Mu_PETSC, IERR )
        Case ('CG')
           ! using the Fletcher-Revees (FR) Conjugate Gradient (CG) method
           If ( (Iter == Iter_begin) .OR. (MOD ( Iter, N_CG_Restart ) == 0) ) Then
              Call VecAXPY ( p_Lambda_PETSC, -1.0d0, g_Lambda_PETSC, IERR )
              Call VecAXPY (     p_Mu_PETSC, -1.0d0,     g_Mu_PETSC, IERR )
           Else
              Call VecNorm (     g_Mu_PETSC,     Norm_2, g_Mu_Norm_2,         IERR )
              Call VecNorm (     g_Mu_PETSC_pre, Norm_2, g_Mu_Norm_2_pre,     IERR )
              Call VecNorm ( g_Lambda_PETSC,     Norm_2, g_Lambda_Norm_2,     IERR )
              Call VecNorm ( g_Lambda_PETSC_pre, Norm_2, g_Lambda_Norm_2_pre, IERR )

              beta_Mu     =     g_Mu_Norm_2 /     g_Mu_Norm_2_pre
              beta_Lambda = g_Lambda_Norm_2 / g_Lambda_Norm_2_pre

              Call VecMAXPY (     p_Mu_PETSC, 2, [-1.0d0,     beta_Mu], [    g_Mu_PETSC,     p_Mu_PETSC_pre], IERR )
              Call VecMAXPY ( p_Lambda_PETSC, 2, [-1.0d0, beta_Lambda], [g_Lambda_PETSC, p_Lambda_PETSC_pre], IERR )
           End If
        Case ('LBFGS')
!====================================================================================================================================================
!====================================================================================================================================================
!           Include 'LBFGS_Decoupled.F90'
           Include 'LBFGS_Coupled.F90'
!====================================================================================================================================================
!====================================================================================================================================================
     End Select

! 3- from part II- Store old results:
     Call VecCopy ( Lambda_PETSC , Lambda_PETSC_pre, IERR )
     Call VecCopy (     Mu_PETSC ,     Mu_PETSC_pre, IERR )

! 4- make sure we have a descent direction
! 4.1- Lambda
     check_grad_search_Lambda = -1.0d0
     Call VecDot ( g_Lambda_PETSC, p_Lambda_PETSC, grad_search_Lambda, IERR )
     If ( grad_search_Lambda > 0 ) Then
        check_grad_search_Lambda = +1.0d0
        If ( rank == 0 ) Then
           Write(*,'(A33,E33.10E3)') 'g.p Lambda = ', grad_search_Lambda
           Write(*,'(A55)') 'Not a descent search direction. Using p = -g instead.'
        End If
        Call VecSet  ( p_Lambda_PETSC,  0.0d0, IERR )
        Call VecAXPY ( p_Lambda_PETSC, -1.0d0, g_Lambda_PETSC, IERR )
        Call VecDot  ( g_Lambda_PETSC, p_Lambda_PETSC, grad_search_Lambda, IERR )
     End If

! 4.2- Mu
     check_grad_search_Mu = -1.0d0
     Call VecDot ( g_Mu_PETSC, p_Mu_PETSC, grad_search_Mu, IERR )
     If ( grad_search_Mu > 0 ) Then
       check_grad_search_Mu = +1.0d0
        If ( rank == 0 ) Then
           Write(*,'(A33,E33.10E3)') 'g.p Mu = ', grad_search_Mu
           Write(*,'(A55)') 'Not a descent search direction. Using p = -g instead.'
        End If
        Call VecSet  ( p_Mu_PETSC,  0.0d0, IERR )
        Call VecAXPY ( p_Mu_PETSC, -1.0d0, g_Mu_PETSC, IERR )
        Call VecDot  ( g_Mu_PETSC, p_Mu_PETSC, grad_search_Mu, IERR )
     End If

! 5- compute a reasonable step size
! 5.1- for SD and CG
     Call VecNorm ( p_Lambda_PETSC, NORM_INFINITY, p_Lambda_max_val, IERR )
     Call VecNorm (     p_Mu_PETSC, NORM_INFINITY,     p_Mu_max_val, IERR )
     alpha_Lambda = Fac_Max_Update_Lambda / p_Lambda_max_val
     alpha_Mu     = Fac_Max_Update_Mu     / p_Mu_max_val
!    this is needed if LBFGS gives ascent direction
     alpha_Lambda_SD = alpha_Lambda
     alpha_Mu_SD     = alpha_Mu

! 5.2- for LBFGS
     If ( (Search_Direction_Method == 'LBFGS') .AND. (Iter /= Iter_begin) ) Then
        alpha_Lambda = 1.0d0
        alpha_Mu     = 1.0d0

       ! engineering watchdog: limit maximum allowable update, even in Qusi-Newton methods
       ! Lambda
        If ( p_Lambda_max_val > Fac_Max_Update_Lambda ) Then
           alpha_Lambda = Fac_Max_Update_Lambda / p_Lambda_max_val
        End If
        ! Mu
        If ( p_Mu_max_val > Fac_Max_Update_Mu ) Then
           alpha_Mu = Fac_Max_Update_Mu / p_Mu_max_val
        End If

     End If

! 5.3- correction in step-length if LBFGS gives ascent direction
! Lambda
     If ( (check_grad_search_Lambda > 0.0d0) .AND. (Search_Direction_Method == 'LBFGS') ) Then
        alpha_Lambda = alpha_Lambda_SD
     End If
! Mu
     If ( (check_grad_search_Mu > 0.0d0) .AND. (Search_Direction_Method == 'LBFGS') ) Then
        alpha_Mu = alpha_Mu_SD
     End If

! ---------------------------------------------------------------------------------------------------------------------------------------------------
! VII- Update the material properties within the Regular Domain (check for bounds later) ------------------------------------------------------------<><><><><><><><><><><><><><><><><><><><>
! ---------------------------------------------------------------------------------------------------------------------------------------------------
     Select Case (Control_Parameter)
        Case ('Lambda')
           Call VecWAXPY ( Lambda_PETSC, alpha_Lambda, p_Lambda_PETSC, Lambda_PETSC_pre, IERR )
        Case ('Mu')
           Call VecWAXPY (     Mu_PETSC, alpha_Mu,         p_Mu_PETSC,     Mu_PETSC_pre, IERR )
        Case ('Lambda_Mu')
!          Call VecWAXPY ( Lambda_PETSC, alpha_Lambda, p_Lambda_PETSC, Lambda_PETSC_pre, IERR )
!          Call VecWAXPY (     Mu_PETSC, alpha_Mu,         p_Mu_PETSC,     Mu_PETSC_pre, IERR )
           Select Case (Bias_Lambda_search)
              Case ('Yes')
                 If ( Iter < NIterBias ) Then
                    ! Bias Lambda.
                    Call VecNorm (     p_Mu_PETSC, Norm_2,     p_Mu_Norm_2, IERR )
                    Call VecNorm ( p_Lambda_PETSC, Norm_2, p_Lambda_Norm_2, IERR )
                    W = 1.0d0 - Dble ( Iter / NIterBias )
                    FAC_Bias_Mu     = W * p_Lambda_Norm_2 / p_Mu_Norm_2 * alpha_Lambda
                    FAC_Bias_Lambda = ( 1.0d0 - W ) * p_Lambda_Norm_2   * alpha_Lambda
                    Call VecCopy  ( Lambda_PETSC_pre , Lambda_PETSC, IERR )
                    Call VecMAXPY ( Lambda_PETSC, 2, [FAC_Bias_Mu, FAC_Bias_Lambda], [p_Mu_PETSC, p_Lambda_PETSC], IERR )
                    ! Mu is fine 
                    Call VecWAXPY (     Mu_PETSC, alpha_Mu,         p_Mu_PETSC,     Mu_PETSC_pre, IERR )
                 Else
                    ! Bias is enough.
                    Call VecWAXPY ( Lambda_PETSC, alpha_Lambda, p_Lambda_PETSC, Lambda_PETSC_pre, IERR )
                    Call VecWAXPY (     Mu_PETSC, alpha_Mu,         p_Mu_PETSC,     Mu_PETSC_pre, IERR )
                 End If
              Case ('No')
                 Call VecWAXPY ( Lambda_PETSC, alpha_Lambda, p_Lambda_PETSC, Lambda_PETSC_pre, IERR )
                 Call VecWAXPY (     Mu_PETSC, alpha_Mu,         p_Mu_PETSC,     Mu_PETSC_pre, IERR )
           End Select !Bias_Lambda_search   
     End Select !Control_Parameter

! ---------------------------------------------------------------------------------------------------------------------------------------------------
! VIII- Extend the material properties from the interface into the PML domain -----------------------------------------------------------------------
! ---------------------------------------------------------------------------------------------------------------------------------------------------
     Call VecScatterBegin ( vscat_mat_extend, Lambda_PETSC, Lambda_PETSC, INSERT_VALUES, SCATTER_FORWARD, IERR )
     Call VecScatterEnd   ( vscat_mat_extend, Lambda_PETSC, Lambda_PETSC, INSERT_VALUES, SCATTER_FORWARD, IERR )
     Call VecScatterBegin ( vscat_mat_extend,     Mu_PETSC,     Mu_PETSC, INSERT_VALUES, SCATTER_FORWARD, IERR )
     Call VecScatterEnd   ( vscat_mat_extend,     Mu_PETSC,     Mu_PETSC, INSERT_VALUES, SCATTER_FORWARD, IERR )

! ---------------------------------------------------------------------------------------------------------------------------------------------------
! IX- Prepare Fortran vectors for matrix computation ------------------------------------------------------------------------------------------------
! ---------------------------------------------------------------------------------------------------------------------------------------------------
! 1- at each rank, all nodes (including ghost nodes) should have material properties, which is needed for assembly of system matrices.
     Call VecScatterBegin ( vscat_Mat, Lambda_PETSC, Lambda_Rank_PETSC, INSERT_VALUES, SCATTER_FORWARD, IERR )
     Call VecScatterEnd   ( vscat_Mat, Lambda_PETSC, Lambda_Rank_PETSC, INSERT_VALUES, SCATTER_FORWARD, IERR )
     Call VecScatterBegin ( vscat_Mat,     Mu_PETSC,     Mu_Rank_PETSC, INSERT_VALUES, SCATTER_FORWARD, IERR )
     Call VecScatterEnd   ( vscat_Mat,     Mu_PETSC,     Mu_Rank_PETSC, INSERT_VALUES, SCATTER_FORWARD, IERR )

! 2- get it in Fortran format for system matrix assembly
     ! Lambda ----------------------
     Call VecGetArrayF90 ( Lambda_Rank_PETSC, Lambda_Rank, IERR )
!... transfer to appropriate Fortran matrices & ascertain Lambda is within pre-specified limits
     Do I = 1 , NJ
        If ( Lambda_Rank (I) < bound_Lambda_lower ) Then
           PMat_Lambda (I) = bound_Lambda_lower
           Lambda_Rank (I) = bound_Lambda_lower
        Else If ( Lambda_Rank (I) > bound_Lambda_upper ) Then
           PMat_Lambda (I) = bound_Lambda_upper
           Lambda_Rank (I) = bound_Lambda_upper
        Else
           PMat_Lambda (I) = Lambda_Rank (I)
        End If
     End Do
!... end transfer and engineering manipulation
     Call VecRestoreArrayF90 ( Lambda_Rank_PETSC, Lambda_Rank, IERR )

     ! Mu --------------------------
     Call VecGetArrayF90 ( Mu_Rank_PETSC, Mu_Rank, IERR )
!... transfer to appropriate Fortran matrices & ascertain Mu is within pre-specified limits
     Do I = 1 , NJ
        If ( Mu_Rank (I) < bound_Mu_lower ) Then
           PMat_Mu (I) = bound_Mu_lower
           Mu_Rank (I) = bound_Mu_lower
        Else If ( Mu_Rank (I) > bound_Mu_upper ) Then
           PMat_Mu (I) = bound_Mu_upper
           Mu_Rank (I) = bound_Mu_upper
        Else
           PMat_Mu (I) = Mu_Rank (I)
        End If
     End Do
!... end transfer and engineering manipulation     
     Call VecRestoreArrayF90 ( Mu_Rank_PETSC, Mu_Rank, IERR )

! 3- scatter from Lambda_Rank_PETSC to Lambda_PETSC (update Global vectors)
     Call VecScatterBegin ( vscat_Mat, Lambda_Rank_PETSC, Lambda_PETSC, INSERT_VALUES, SCATTER_REVERSE_LOCAL, IERR )
     Call VecScatterEnd   ( vscat_Mat, Lambda_Rank_PETSC, Lambda_PETSC, INSERT_VALUES, SCATTER_REVERSE_LOCAL, IERR )
     Call VecScatterBegin ( vscat_Mat,     Mu_Rank_PETSC,     Mu_PETSC, INSERT_VALUES, SCATTER_REVERSE_LOCAL, IERR )
     Call VecScatterEnd   ( vscat_Mat,     Mu_Rank_PETSC,     Mu_PETSC, INSERT_VALUES, SCATTER_REVERSE_LOCAL, IERR )

! ---------------------------------------------------------------------------------------------------------------------------------------------------
! X- Ensure sufficient decrease of the cost functional ----------------------------------------------------------------------------------------------
! ---------------------------------------------------------------------------------------------------------------------------------------------------
! 1- initialize parallel sparce matrices
     Call MatZeroEntries ( AM_PETSC, IERR )
     Call MatZeroEntries ( AC_PETSC, IERR )
     Call MatZeroEntries ( AK_PETSC, IERR )
     Call MatZeroEntries ( AG_PETSC, IERR )
     Call MatZeroEntries ( AK_RD_PETSC, IERR )
     Call MatZeroEntries ( AM_RD_PETSC, IERR )

! 2- initialize parallel vectors
     Call VecSet ( B_PETSC,         0.0d0, IERR )
     Call VecSet ( Diag_M_PETSC,    0.0d0, IERR )
     Call VecSet ( Diag_iM_PETSC,   0.0d0, IERR )
     Call VecSet ( RMJ_PETSC,       0.0d0, IERR )
     Call VecSet ( Diag_M_RD_PETSC, 0.0d0, IERR )

! 3- Solve State
! 3.1- use actual {Lambda} and {Mu} for solving state
     If ( Lambda_variable == 'Lambda_2Mu' ) Then
       PMat_Lambda = PMat_Lambda - 2.0d0 * PMat_Mu
     End If
! 3.2- solve
     Include 'Solve_Forward.F90'
! 3.3- use {Lambda + 2 Mu} instead of {Lambd} for inversion.
     If ( Lambda_variable == 'Lambda_2Mu' ) Then
       PMat_Lambda = PMat_Lambda + 2.0d0 * PMat_Mu
     End If

! 4- Compute misfit
     Call Compute_Misfit ( B_mis_PETSC, NNDH, EqDis, U_Store_Mapping, NStore_Mapping, EqDis_MAP_U_Store_Numbers_Global, Dis_meas, Cost_Misfit )
     Cost_Functional = Cost_Misfit + Cost_Regularization_Mu + Cost_Regularization_Lambda

! 5- Check Armijo <><><><><><><><><><><><>                                                                                                          
     Armijo_Counter = 0
     Fac_Armijo_rho = 0.50d0 ! suggested to be 0.5, but Sez uses 0.1.
     Fac_Armijo_c   = 1.0d-4 ! should be 4

     Select Case (Control_Parameter)
        Case ('Lambda')
           Armijo = Fac_Armijo_c * DAbs(alpha_Lambda) * grad_search_Lambda
        Case ('Mu')
           Armijo = Fac_Armijo_c * DAbs(alpha_Mu) * grad_search_Mu
        Case ('Lambda_Mu')
           Armijo = Fac_Armijo_c * DAbs(alpha_Mu) * grad_search_Mu + Fac_Armijo_c * DAbs(alpha_Lambda) * grad_search_Lambda
     End Select

! - $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
     Do While ( Cost_Misfit > Cost_Misfit_pre + Armijo )   ! ----------------------------$$$$$$$$$$$

        Armijo_Counter = Armijo_Counter + 1
        If (rank == 0) Write(*,'(I6,E40.10E3,A20)') Armijo_Counter, Cost_Misfit, 'Armijo triggered'

! 5.1- adaptive step-length: depending on Armijo behaviour, adjust alpha_Mu, alpha_Lambda for future iterations.
        ! comment: Our heuristic is that if we get into Armijo twice, consecutively, alpha needs adjustment.
        If ( Armijo_Counter == 3 ) Then
           Fac_Max_Update_Lambda = Fac_Max_Update_Lambda / 6.0d0
           Fac_Max_Update_Mu     = Fac_Max_Update_Mu / 6.0d0
        End If

! 5.2- compute new step-size length
        alpha_Lambda = alpha_Lambda * Fac_Armijo_rho
        alpha_Mu     = alpha_Mu     * Fac_Armijo_rho

! ** This is wrong, and shows that the gradient is inconsistent
        If ( Armijo_Counter == 9 ) Then
           alpha_Lambda = -alpha_Lambda
           alpha_Mu     = -alpha_Mu
        If (rank == 0) write(*,*) 'gradient re-signed'
        End If

! 5.3- update Lambda & Mu within the RD                                                                                                              !<><><><><><><><><><><><><><><><><><><><>
        Select Case (Control_Parameter)
           Case ('Lambda')
              Call VecWAXPY ( Lambda_PETSC, alpha_Lambda, p_Lambda_PETSC, Lambda_PETSC_pre, IERR )
           Case ('Mu')
              Call VecWAXPY (     Mu_PETSC, alpha_Mu,         p_Mu_PETSC,     Mu_PETSC_pre, IERR )
           Case ('Lambda_Mu')
!             Call VecWAXPY ( Lambda_PETSC, alpha_Lambda, p_Lambda_PETSC, Lambda_PETSC_pre, IERR )
!             Call VecWAXPY (     Mu_PETSC, alpha_Mu,         p_Mu_PETSC,     Mu_PETSC_pre, IERR )
              Select Case (Bias_Lambda_search)
                 Case ('Yes')
                    If ( Iter < NIterBias ) Then
                       ! Bias Lambda.
                       Call VecNorm (     p_Mu_PETSC, Norm_2,     p_Mu_Norm_2, IERR )
                       Call VecNorm ( p_Lambda_PETSC, Norm_2, p_Lambda_Norm_2, IERR )
                       W = 1.0d0 - Dble ( Iter / NIterBias )
                       FAC_Bias_Mu     = W * p_Lambda_Norm_2 / p_Mu_Norm_2 * alpha_Lambda
                       FAC_Bias_Lambda = ( 1.0d0 - W ) * p_Lambda_Norm_2   * alpha_Lambda
                       Call VecCopy  ( Lambda_PETSC_pre , Lambda_PETSC, IERR )
                       Call VecMAXPY ( Lambda_PETSC, 2, [FAC_Bias_Mu, FAC_Bias_Lambda], [p_Mu_PETSC, p_Lambda_PETSC], IERR )
                       ! Mu is fine
                       Call VecWAXPY (     Mu_PETSC, alpha_Mu,         p_Mu_PETSC,     Mu_PETSC_pre, IERR )
                    Else
                       ! Bias is enough.
                       Call VecWAXPY ( Lambda_PETSC, alpha_Lambda, p_Lambda_PETSC, Lambda_PETSC_pre, IERR )
                       Call VecWAXPY (     Mu_PETSC, alpha_Mu,         p_Mu_PETSC,     Mu_PETSC_pre, IERR )
                    End If
                 Case ('No')
                    Call VecWAXPY ( Lambda_PETSC, alpha_Lambda, p_Lambda_PETSC, Lambda_PETSC_pre, IERR )
                    Call VecWAXPY (     Mu_PETSC, alpha_Mu,         p_Mu_PETSC,     Mu_PETSC_pre, IERR )
              End Select !Bias_Lambda_search   
        End Select !Control_Parameter

! 5.4- Extend the material properties from the interface into the PML domain
        Call VecScatterBegin ( vscat_mat_extend, Lambda_PETSC, Lambda_PETSC, INSERT_VALUES, SCATTER_FORWARD, IERR )
        Call VecScatterEnd   ( vscat_mat_extend, Lambda_PETSC, Lambda_PETSC, INSERT_VALUES, SCATTER_FORWARD, IERR )
        Call VecScatterBegin ( vscat_mat_extend,     Mu_PETSC,     Mu_PETSC, INSERT_VALUES, SCATTER_FORWARD, IERR )
        Call VecScatterEnd   ( vscat_mat_extend,     Mu_PETSC,     Mu_PETSC, INSERT_VALUES, SCATTER_FORWARD, IERR )

! 5.5- at each rank, all nodes (including ghost nodes) should have material properties, which is needed for assembly of system matrices.
        Call VecScatterBegin ( vscat_Mat, Lambda_PETSC, Lambda_Rank_PETSC, INSERT_VALUES, SCATTER_FORWARD, IERR )
        Call VecScatterEnd   ( vscat_Mat, Lambda_PETSC, Lambda_Rank_PETSC, INSERT_VALUES, SCATTER_FORWARD, IERR )
        Call VecScatterBegin ( vscat_Mat,     Mu_PETSC,     Mu_Rank_PETSC, INSERT_VALUES, SCATTER_FORWARD, IERR )
        Call VecScatterEnd   ( vscat_Mat,     Mu_PETSC,     Mu_Rank_PETSC, INSERT_VALUES, SCATTER_FORWARD, IERR )

! 5.6- get it in Fortran format for system matrix assembly
     ! Lambda ----------------------
        Call VecGetArrayF90 ( Lambda_Rank_PETSC, Lambda_Rank, IERR )
!... transfer to appropriate Fortran matrices & ascertain Lambda is within pre-specified limits
        Do I = 1 , NJ
           If ( Lambda_Rank (I) < bound_Lambda_lower ) Then
              PMat_Lambda (I) = bound_Lambda_lower
              Lambda_Rank (I) = bound_Lambda_lower
           Else If ( Lambda_Rank (I) > bound_Lambda_upper ) Then
              PMat_Lambda (I) = bound_Lambda_upper
              Lambda_Rank (I) = bound_Lambda_upper
           Else
              PMat_Lambda (I) = Lambda_Rank (I)
           End If
        End Do
!... end transfer and engineering manipulation
        Call VecRestoreArrayF90 ( Lambda_Rank_PETSC, Lambda_Rank, IERR )

     ! Mu --------------------------
        Call VecGetArrayF90 ( Mu_Rank_PETSC, Mu_Rank, IERR )
!... transfer to appropriate Fortran matrices & ascertain Mu is within pre-specified limits
        Do I = 1 , NJ
           If ( Mu_Rank (I) < bound_Mu_lower ) Then
              PMat_Mu (I) = bound_Mu_lower
              Mu_Rank (I) = bound_Mu_lower
           Else If ( Mu_Rank (I) > bound_Mu_upper ) Then
              PMat_Mu (I) = bound_Mu_upper
              Mu_Rank (I) = bound_Mu_upper
           Else
              PMat_Mu (I) = Mu_Rank (I)
           End If
        End Do
!... end transfer and engineering manipulation
        Call VecRestoreArrayF90 ( Mu_Rank_PETSC, Mu_Rank, IERR )

! 5.7- scatter from Lambda_Rank_PETSC to Lambda_PETSC (update Global vectors)
        Call VecScatterBegin ( vscat_Mat, Lambda_Rank_PETSC, Lambda_PETSC, INSERT_VALUES, SCATTER_REVERSE_LOCAL, IERR )
        Call VecScatterEnd   ( vscat_Mat, Lambda_Rank_PETSC, Lambda_PETSC, INSERT_VALUES, SCATTER_REVERSE_LOCAL, IERR )
        Call VecScatterBegin ( vscat_Mat,     Mu_Rank_PETSC,     Mu_PETSC, INSERT_VALUES, SCATTER_REVERSE_LOCAL, IERR )
        Call VecScatterEnd   ( vscat_Mat,     Mu_Rank_PETSC,     Mu_PETSC, INSERT_VALUES, SCATTER_REVERSE_LOCAL, IERR )

! 5.8- initialize parallel sparce matrices
        Call MatZeroEntries ( AM_PETSC, IERR )
        Call MatZeroEntries ( AC_PETSC, IERR )
        Call MatZeroEntries ( AK_PETSC, IERR )
        Call MatZeroEntries ( AG_PETSC, IERR )
        Call MatZeroEntries ( AK_RD_PETSC, IERR )
        Call MatZeroEntries ( AM_RD_PETSC, IERR )

! 5.9- initialize parallel vectors
        Call VecSet ( B_PETSC,         0.0d0, IERR )
        Call VecSet ( Diag_M_PETSC,    0.0d0, IERR )
        Call VecSet ( Diag_iM_PETSC,   0.0d0, IERR )
        Call VecSet ( RMJ_PETSC,       0.0d0, IERR )
        Call VecSet ( Diag_M_RD_PETSC, 0.0d0, IERR )

! 5.10- Solve State
! 5.10.1- use actual {Lambda} and {Mu} for solving state
        If ( Lambda_variable == 'Lambda_2Mu' ) Then
          PMat_Lambda = PMat_Lambda - 2.0d0 * PMat_Mu
        End If
! 5.10.2- solve
        Include 'Solve_Forward.F90'
! 5.10.3- use {Lambda + 2 Mu} instead of {Lambd} for inversion.
        If ( Lambda_variable == 'Lambda_2Mu' ) Then
          PMat_Lambda = PMat_Lambda + 2.0d0 * PMat_Mu
        End If

! 5.11- Compute misfit
        Call Compute_Misfit ( B_mis_PETSC, NNDH, EqDis, U_Store_Mapping, NStore_Mapping, EqDis_MAP_U_Store_Numbers_Global, Dis_meas, Cost_Misfit )
        Cost_Functional = Cost_Misfit + Cost_Regularization_Mu + Cost_Regularization_Lambda

! 5.12- Check Armijo again
! comment: The absolute value is not necessary; it is a temporary check due to gradient re-sign, which indicates something is indeed wrong.          !<><><><><><><><><><><><><><><><><><><><>
        Select Case (Control_Parameter)
           Case ('Lambda')
              Armijo = Fac_Armijo_c * DAbs(alpha_Lambda) * grad_search_Lambda
           Case ('Mu')
              Armijo = Fac_Armijo_c * DAbs(alpha_Mu) * grad_search_Mu
           Case ('Lambda_Mu')
              Armijo = Fac_Armijo_c * DAbs(alpha_Mu) * grad_search_Mu + Fac_Armijo_c * DAbs(alpha_Lambda) * grad_search_Lambda
        End Select

! 5.13- Express grief and pray to the Lord
        If ( Armijo_Counter == 10 ) Then
           If (rank == 0) Write(*,*) ' Armijo condition failure. '
           Include 'Material_Visualization_End.F90'
           Call MPI_Barrier ( PETSc_COMM_WORLD, IERR )
           Stop
        End If

     End Do        ! --------------- End Armijo loop ------------------------------------$$$$$$$$$$$
! - $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


! ---------------------------------------------------------------------------------------------------------------------------------------------------
! XI- Printout misfit, cost, visualization, and prepare for the unexpected --------------------------------------------------------------------------
! ---------------------------------------------------------------------------------------------------------------------------------------------------
! 0- compute gradient norm
     Call VecNorm (     g_Mu_PETSC, Norm_2, g_Mu_Norm_2,     IERR )
     Call VecNorm ( g_Lambda_PETSC, Norm_2, g_Lambda_Norm_2, IERR )

! 1- print screen
     If ( Rank == 0 ) Then
       UnFile = UN_Inversion_Iterations ;
       Write ( *     ,'(I6,6E20.10E3,F20.4)' ) Iter, Cost_Functional, Cost_Misfit, Cost_Regularization_Lambda, Cost_Regularization_Mu, g_Lambda_Norm_2, g_Mu_Norm_2, Regularization_factor_continuation
       Write ( UnFile,'(I6,6E20.10E3,F20.4)' ) Iter, Cost_Functional, Cost_Misfit, Cost_Regularization_Lambda, Cost_Regularization_Mu, g_Lambda_Norm_2, g_Mu_Norm_2, Regularization_factor_continuation
     End If

! 2- store material updates to resume inversion if interrupted
     Include 'Material_Save_Disk.F90'

! 3- Output Lambda_PETSC, Mu_PETSC, Iter for Paraview -----------------------------------------------------------------------------------------------
     Include 'Material_Visualization.F90'                                                                ! R E P E A T E D
! - --------------------------------*----------------------------------------------------------------------------------------------------------------

   End Do
! ===================================================================================================================================================
! *** end iteration ***
! ===================================================================================================================================================

! - Call only once at the end of all material updates to shut down top file wrapper -----------------------------------------------------------------
Include 'Material_Visualization_End.F90'                                                            ! O N L Y   O N C E
! - --------------------------------*----------------------------------------------------------------------------------------------------------------

!==================================================================================================
! - SHUT DOWN PETSc
IF (RANK == 0) WRITE(*,'(A100)')'Destroying PETSc matrices'
IF (RANK == 0) WRITE(*,'(A88)')'---------- ---------- ---------- ---------- ---------- ---------- ---------- ---------- ---------- ----------'
CALL MatDestroy    ( AK_PETSC        , IERR )
CALL MatDestroy    ( AM_PETSC        , IERR )
CALL MatDestroy    ( AC_PETSC        , IERR )
CALL MatDestroy    ( AG_PETSC        , IERR )
CALL MatDestroy    ( AK_RD_PETSC     , IERR )
CALL MatDestroy    ( AM_RD_PETSC     , IERR )

CALL VecDestroy    ( B_PETSC         , IERR )
CALL VecDestroy    ( RMJ_PETSC       , IERR )
CALL VecDestroy    ( DIAG_M_PETSC    , IERR )
CALL VecDestroy    ( DIAG_iM_PETSC   , IERR )
CALL VecDestroy    ( DIAG_M_RD_PETSC , IERR )

!Call AODestroy ( AOOut, ErrPTC ) ; 
Call ISLocalToGlobalMappingDestroy ( GMapping, ErrPTC ) ;
Call ISLocalToGlobalMappingDestroy ( GMapping_Mat, ErrPTC )

! Objects associated to the inverse problem.
Call VecDestroy (      Lambda_PETSC , IERR )
Call VecDestroy (          Mu_PETSC , IERR )
Call VecDestroy ( Lambda_Rank_PETSC , IERR )
Call VecDestroy (     Mu_Rank_PETSC , IERR )
Call VecDestroy (    g_Lambda_PETSC , IERR )
Call VecDestroy (        g_Mu_PETSC , IERR )
Call VecDestroy (    p_Lambda_PETSC , IERR )
Call VecDestroy (        p_Mu_PETSC , IERR )

CALL ISDestroy ( IS_Mat_from,        IERR )
CALL ISDestroy ( IS_Mat_to,          IERR )
CALL ISDestroy ( IS_Mat_Extend_from, IERR )
CALL ISDestroy ( IS_Mat_Extend_to,   IERR )
CALL ISDestroy ( IS_u_from,          IERR )
CALL ISDestroy ( IS_u_to,            IERR )

Call VecScatterDestroy ( vscat_Mat,        IERR )
Call VecScatterDestroy ( vscat_Mat_Extend, IERR )

CALL VecDestroy (  B_mis_PETSC, IERR )
CALL VecDestroy ( Ux_big_PETSc, IERR )
CALL VecDestroy ( Uk_big_PETSc, IERR )

CALL VecDestroy ( Ux_Bigger_Rank_PETSc, IERR )
CALL VecDestroy ( Uk_Bigger_Rank_PETSc, IERR )
CALL VecDestroy (  P_Bigger_Rank_PETSc, IERR )
! Control problem
CALL MatDestroy ( Reg_PETSC, IERR )
CALL VecDestroy ( iMg_PETSC, IERR )
CALL VecDestroy ( Reg_Lambda_PETSC, IERR )
CALL VecDestroy (     Reg_Mu_PETSC, IERR )

CALL PetscFinalize ( IERR )
!==================================================================================================


CALL CPU_TIME (t7)
IF (RANK == 0) WRITE(*,'(A40,F10.2,A10)')'total operation time:', t7-t0, 'seconds'
IF (RANK == 0) WRITE(2,'(A40,F10.2,A10)')'total operation time:', t7-t0, 'seconds'
IF (RANK == 0) WRITE(*,'(A88)')'---------- ---------- ---------- ---------- ---------- ---------- ---------- ---------- ---------- ----------'
IF (RANK == 0) WRITE(*,'(A10)')'((!))'

Stop

! =============================================== ERROR IN Write STATEMENT ==========================================================================
1006    Write(*       , Fmt_Write1 ) UnFile, IO_Write ; Write( UnFile, Fmt_Write1 ) UnFile, IO_Write ;
        !#Call BEEP_FAIL ;
        Write(*, Fmt_FL) ;  Write(UnInf, Fmt_FL) ; Write(*, Fmt_End) ; Read(*,*) ; Stop ;

! ===================================================================================================================================================


END
