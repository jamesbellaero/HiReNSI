!************************************************
! NEWMARK: implicit time-stepping - advanced preconditioning
!************************************************
! REVISION : Wed, 18 January 2012

SUBROUTINE NEWMARK_PETSC_1 ( AK_PETSC, AM_PETSC, AC_PETSC, B_PETSC, AK_RD_PETSC, AM_RD_PETSC, ID, LTRANS, RMJ_PETSC, BACL, RANK, SIZE )
USE PARAMETERS
USE IFPORT
IMPLICIT DOUBLE PRECISION (A-H,O-Z)

character(30)     solver_package

!==================================================================================================
! - PETSC INCLUDE FILES
#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscmat.h"
#include "finclude/petscksp.h"
#include "finclude/petscpc.h"
#include "finclude/petscis.h"
#include "finclude/petscis.h90"
#include "finclude/petscdmmesh.h"
!#include "finclude/petscdmcomposite.h"
#include "finclude/ftn-custom/petscdmcomposite.h90"
#include "finclude/petscdmda.h"

! - PETSC VARIABLES AND MATRICES
Mat            :: AK_PETSC, AM_PETSC, AC_PETSC, AK_RD_PETSC, AM_RD_PETSC
Vec            :: B_PETSC
Vec            :: RMJ_PETSC

KSP            :: ksp
PC             :: pc
KSPConvergedReason :: reason

IS             :: perm
MatFactorInfo  :: mat_factor_info(MAT_FACTORINFO_SIZE)

PetscReal      :: rtol
PetscInt       :: its, m

PetscErrorCode :: IERR
PetscMPIInt    :: SIZE, RANK
PetscScalar    :: ENERGY

PetscScalar    :: zero, one, minus_one
PetscScalar    :: A0, A1, A2, A3, A4, A5
PetscScalar    :: AG_X

! - preconditioner
IS             :: is_u, is_s
PetscInt       :: NEQ_U, NEQ_S
! - preconditioner BJacobi
PetscInt       :: local_size, blocks, lens(SIZE)

! - PETSc local vectors
Vec            :: R_PETSC, U_PETSC, UD_PETSC, UDD_PETSC, DU_PETSC, UDD_OLD_PETSC
Vec            :: HELP_PETSC
!==================================================================================================


! in 
!---------- ---------- -------- -- ---------- ----------
DIMENSION BACL(NSTEP,NDIM)
DIMENSION ID(NJ, NDOF), LTRANS(NTRANS)
!---------- ---------- ---------- ---------- ----------


! local
!---------- ---------- ---------- ---------- ----------
DIMENSION AG(NDIM)
!---------- ---------- ---------- ---------- ----------


zero      =  0.0d0
one       =  1.0d0
minus_one = -1.0d0
rtol      =  1.d-7
its       =  0


! (1) earthquake analysis        vs  (2) specified external load
! (1) SuperLU_Dist direct solver vs  (2) Krylov iterative solver
! (1) direct solver factor LU    vs  (2) Cholesky <symmetric>
!---------- ---------- ---------- ---------- ----------
IFLAG_LOAD   = 2
IFLAG_SOLVER = 2
IFLAG_FACTOR = 1


!==================================================================================================
! - DEFINE PETSC OBJECTS
CALL VecCreateMPI ( PETSC_COMM_WORLD, PETSC_DECIDE, NEQ, R_PETSC      , IERR)
CALL VecCreateMPI ( PETSC_COMM_WORLD, PETSC_DECIDE, NEQ, U_PETSC      , IERR)
CALL VecCreateMPI ( PETSC_COMM_WORLD, PETSC_DECIDE, NEQ, UD_PETSC     , IERR)
CALL VecCreateMPI ( PETSC_COMM_WORLD, PETSC_DECIDE, NEQ, UDD_PETSC    , IERR)
CALL VecCreateMPI ( PETSC_COMM_WORLD, PETSC_DECIDE, NEQ, DU_PETSC     , IERR)
CALL VecCreateMPI ( PETSC_COMM_WORLD, PETSC_DECIDE, NEQ, UDD_OLD_PETSC, IERR)
CALL VecCreateMPI ( PETSC_COMM_WORLD, PETSC_DECIDE, NEQ, HELP_PETSC   , IERR)

CALL VecSet       ( R_PETSC  , zero, ierr)
CALL VecSet       ( U_PETSC  , zero, ierr)
CALL VecSet       ( UD_PETSC , zero, ierr)
CALL VecSet       ( UDD_PETSC, zero, ierr)
CALL VecSet       ( DU_PETSC , zero, ierr)

! factored matrix for direct solver; see if you can factor in place and get rid of this matrix
CALL MatCreateMPIAIJ ( PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, NEQ, NEQ, 100, PETSC_NULL_INTEGER, 100, PETSC_NULL_INTEGER, AF_PETSC, IERR )
!==================================================================================================


! Newmark parameters, initialize
!---------- ---------- ---------- ---------- ----------  
A0 = 1.0d0 / (ALPHA * DT**2)
A1 = DELTA / (ALPHA * DT)
A2 = 1.0d0 / (ALPHA * DT)
A3 = 1.0d0 / (2.0d0 * ALPHA) - 1.0d0
A4 = DELTA / ALPHA - 1.0d0
A5 = (DELTA / (2.0d0 * ALPHA) - 1.0d0) * DT
!---------- ---------- ---------- ---------- ---------- 
T  = 0.0d0
!---------- ---------- ---------- ---------- ---------- 


! form the effective "stiffness matrix"
!---------- ---------- ---------- ---------- ----------
! AK_f  = A0 * AM_f  +  A1 * AC_f  +  AK_f
CALL MatAXPY ( AK_PETSC, A0, AM_PETSC, DIFFERENT_NONZERO_PATTERN, IERR )
CALL MatAXPY ( AK_PETSC, A1, AC_PETSC, DIFFERENT_NONZERO_PATTERN, IERR )


! SuperLU_Dist direct solver(1) vs Krylov iterative solver(2)
! for direct solver, perform matrix factorization, for iterative solver, setup appropriate solver and preconditioner
!---------- ---------- ---------- ---------- ---------- ---------- ---------- ---------- ----------
CALL CPU_TIME(t0)
SELECT CASE (IFLAG_SOLVER)

   CASE(1)
! SuperLU_Dist solver 
!---------- ---------- ---------- ---------- ----------
      SELECT CASE (IFLAG_FACTOR)
         
! LU factorization
         CASE(1)

            WRITE(*,*) ' LU factorization begins ...'
            CALL MatGetFactor        ( AK_PETSC, "superlu_dist", MAT_FACTOR_LU, AF_PETSC, IERR )
            CALL MatGetSize          ( AK_PETSC, m, PETSC_NULL, IERR )
            CALL ISCreateStride      ( PETSC_COMM_WORLD, m, 0, 1, perm, IERR )
            CALL MatLUFactorSymbolic ( AF_PETSC, AK_PETSC, perm, perm, mat_factor_info, IERR )
            CALL MatLUFactorNumeric  ( AF_PETSC, AK_PETSC, mat_factor_info, IERR )


! Cholesky factorization
         CASE(2)
             write(*,*) ' Cholesky factorization: coming soon!'
             stop

      END SELECT


! spit out what direct solver we used
      CALL MatFactorGetSolverPackage ( AF_PETSC, solver_package, IERR )
      WRITE(*,*) 'Solver Package: ', solver_package
      WRITE(2,*) 'Solver Package: ', solver_package


   CASE(2)
! KSP solver 
!---------- ---------- ---------- ---------- ----------
      CALL KSPCreate         ( PETSC_COMM_WORLD, ksp, IERR )
      CALL KSPSetOperators   ( ksp, AK_PETSC, AK_PETSC, SAME_PRECONDITIONER, IERR )
      CALL KSPSetType        ( ksp, KSPGMRES  , IERR)
!      CALL KSPSetInitialGuessNonzero ( ksp, PETSC_FALSE, IERR )
! Question: when I make this PETSC_TRUE, the code cannot pass the KSPSolve routine on more than 1 processor; why?!


! spit out iterative solver
      CALL PetscPrintf ( PETSC_COMM_WORLD, "Solver Package: Krylov iterative solver \n", IERR )
      WRITE(2,*) "Solver Package: Krylov iterative solver"


! Pre-conditioner
IFLAG_PRECONDITIONER = 1
!---------- ---------- ---------- ---------- ---------- ----------
      CALL KSPGetPC ( ksp, pc, IERR )


    SELECT CASE (IFLAG_PRECONDITIONER)

       CASE(1)
! 1. block Jacobi preconditioning 
!    comment - currently only works for the small example on 2 processors
!---------- ---------- ---------- ----------
!            CALL PCSetType ( pc, PCBJACOBI, IERR )
            CALL PCSetType ( pc, PCNONE, IERR )
 !           CALL VecGetLocalSize ( U_PETSC, local_size, IERR )

            blocks  = 2 
            lens(1) = 95
            lens(2) = 94

!            CALL PCBJacobiSetTotalBlocks ( pc, blocks, lens, IERR )
!            CALL KSPSetFromOptions ( ksp, IERR )
!            CALL PCSetFromOptions  (  pc, IERR )
!Comment: use -sub_pc_type ilu -sub_ksp_type gmres -sub_pc_factor_levels 1, at the front end            
!---------- ---------- ---------- ----------


       CASE(2)
! block preconditioning
!---------- ---------- ---------- ----------
! 1. construct the index sets corresponding to each block
            NEQ_U = NEQ_DISP
            NEQ_S = NEQ - NEQ_DISP
            CALL ISCreateStride ( PETSC_COMM_WORLD, NEQ_U, 0    , 1, is_u, IERR )
            CALL ISCreateStride ( PETSC_COMM_WORLD, NEQ_S, NEQ_U, 1, is_s, IERR )

!call ISView(is_u,PETSC_VIEWER_STDOUT_WORLD,ierr)
!write(*,*)'**************************'
!call ISView(is_s,PETSC_VIEWER_STDOUT_WORLD,ierr)



! 2. specify the general type of preconditioner
            CALL PCSetType         ( pc, PCFIELDSPLIT, IERR )

! 3. indicate exactly which rows/columns of the matrix belong to a particular block
            CALL PCFieldSplitSetIS ( pc , PETSC_NULL_CHARACTER, is_u, IERR )
print*,'isu', ierr

            CALL PCFieldSplitSetIS ( pc , PETSC_NULL_CHARACTER, is_s, IERR )

print*,'is_s',ierr
! 4. type
!            CALL PCFieldSplitSetType ( pc, PC_COMPOSITE_SCHUR, IERR )
!            CALL PCFieldSplitSchurPrecondition ( pc, PC_FIELDSPLIT_SCHUR_PRE_DIAG, PETSC_NULL_OBJECT, IERR )

! 5. from options
            CALL KSPSetFromOptions ( ksp, IERR )
            CALL PCSetFromOptions  (  pc, IERR )

call PCSetUp( pc, IERR)

    END SELECT
!---------- ---------- ---------- ---------- ---------- ----------


      CALL KSPSetTolerances ( ksp, PETSC_DEFAULT_DOUBLE_PRECISION, PETSC_DEFAULT_DOUBLE_PRECISION, PETSC_DEFAULT_DOUBLE_PRECISION, PETSC_DEFAULT_INTEGER, IERR)


END SELECT
CALL CPU_TIME(t1)
IF ( RANK == 0 ) WRITE(*,'(A40,F10.2,A10)') 'Factorization / Preconditioning time:', t1-t0, 'seconds'
IF ( RANK == 0 ) WRITE(2,'(A40,F10.2,A10)') 'Factorization / Preconditioning time:', t1-t0, 'seconds'
!---------- ---------- ---------- ---------- ---------- ---------- ---------- ---------- ----------


! NEWMARK INTEGRATION
!---------- ---------- ---------- ---------- ----------
          DO ISTEP = 1 , NSTEP


! print out
!---------- ---------- ---------- ---------- ----------
CALL OUT_NEWMARK_PETSC ( U_PETSC, ID, LTRANS, T, ENERGY )
!---------- ---------- ---------- ---------- ---------- 


             T = DBLE(ISTEP) * DT
             IF ( RANK == 0 ) WRITE (*,'(A20,F10.5,A20,I10)') 'TIME=', T, 'KSP Iter=', its


! earthquake analysis (1) vs specified external load (2)
!---------- ---------- ---------- ---------- ---------- ---------- ----------
SELECT CASE (IFLAG_LOAD)


   CASE(1)
! earthquake analysis case
! only x-component of acceleration has been considered. see Assemble_PETSc.F90
!---------- ---------- ---------- ---------- ----------
      AG_X = BACL(ISTEP, 1)
      CALL VecSet   ( R_PETSC, zero             , IERR )
      CALL VecAXPY  ( R_PETSC, AG_X, RMJ_PETSC  , IERR )
!---------- ---------- ---------- ---------- ---------- 


   CASE(2)
!  external load case ( forward & inverse problem )
!---------- ---------- ---------- ---------- ---------- 
      CALL TIME_FUNCTION ( T, FAC1 )
      CALL VecSet   ( R_PETSC, zero             , IERR )
      CALL VecAXPY  ( R_PETSC, FAC1, B_PETSC    , IERR )
!---------- ---------- ---------- ---------- ----------


END SELECT
!---------- ---------- ---------- ---------- ---------- ---------- ----------


! Newmark's method
!---------- ---------- ---------- ---------- ----------        
! R = R + MATMUL(AM_f, (A0 * U + A2 * UD + A3 * UDD)) + MATMUL(AC_f, (A1 * U + A4 * UD + A5 * UDD))
CALL VecSet     ( HELP_PETSC, zero         , IERR )
CALL VecAXPY    ( HELP_PETSC, A0, U_PETSC  , IERR )
CALL VecAXPY    ( HELP_PETSC, A2, UD_PETSC , IERR )
CALL VecAXPY    ( HELP_PETSC, A3, UDD_PETSC, IERR )
CALL MatMultAdd ( AM_PETSC, HELP_PETSC, R_PETSC, R_PETSC, IERR )

CALL VecSet     ( HELP_PETSC, zero         , IERR )
CALL VecAXPY    ( HELP_PETSC, A1, U_PETSC  , IERR )
CALL VecAXPY    ( HELP_PETSC, A4, UD_PETSC , IERR )
CALL VecAXPY    ( HELP_PETSC, A5, UDD_PETSC, IERR )
CALL MatMultAdd ( AC_PETSC, HELP_PETSC, R_PETSC, R_PETSC, IERR )


! SuperLU_Dist direct solver(1) vs Krylov iterative solver(2)
! solve for the same coefficient matrix but different right hand sides
!---------- ---------- ---------- ---------- ---------- ---------- ----------
SELECT CASE (IFLAG_SOLVER)


   CASE(1)
! SuperLu_Dist solver
!---------- ---------- ---------- ---------- ----------
! solve
      CALL MatSolve   ( AF_PETSC, R_PETSC, HELP_PETSC, IERR )
! transfer
      CALL VecCopy    ( HELP_PETSC, R_PETSC, IERR )


   CASE(2)
! Krylov solver
!---------- ---------- ---------- ---------- ----------
! initial guess from previous iteration
      CALL VecCopy    ( U_PETSC, HELP_PETSC, IERR )
! solve
      write(*,*) 'before calling KSPSolve'
      CALL KSPSolve   ( ksp, R_PETSC, HELP_PETSC, IERR )
      write(*,*) 'after  calling KSPSolve'
! transfer
      CALL VecCopy    ( HELP_PETSC, R_PETSC, IERR )

! convergence reason
      CALL KSPGetConvergedReason ( ksp, reason, IERR )
      IF ( reason <= 0 ) THEN
         CALL PetscPrintf ( PETSC_COMM_WORLD, "Failure to converge by arash's strict policy!\n", IERR )
         STOP
      END IF
! note: this is off by one step in the printout
      CALL KSPGetIterationNumber ( ksp, its, IERR )
!     CALL KSPView    ( ksp, PETSC_VIEWER_STDOUT_WORLD, IERR )


END SELECT
!---------- ---------- ---------- ---------- ---------- ---------- ----------


! Update Newmark parameters for the next iteration
!---------- ---------- ---------- ---------- ----------
! DU      = R - U
CALL VecCopy    ( R_PETSC,        DU_PETSC      , IERR )
CALL VecAXPY    ( DU_PETSC, minus_one, U_PETSC  , IERR )
! U       = R
CALL VecCopy    ( R_PETSC,        U_PETSC       , IERR )
! UDD_OLD = UDD
CALL VecCopy    ( UDD_PETSC,      UDD_OLD_PETSC , IERR )
! UDD     = A0 * DU - A2 * UD - A3 * UDD_OLD
CALL VecSet     ( UDD_PETSC, zero               , IERR )
CALL VecAXPY    ( UDD_PETSC,  A0, DU_PETSC      , IERR )
CALL VecAXPY    ( UDD_PETSC, -A2, UD_PETSC      , IERR )
CALL VecAXPY    ( UDD_PETSC, -A3, UDD_OLD_PETSC , IERR )
! UD      = A1 * DU - A4 * UD - A5 * UDD_OLD
CALL VecScale   ( UD_PETSC,  -A4                , IERR )
CALL VecAXPY    ( UD_PETSC,   A1, DU_PETSC      , IERR )
CALL VecAXPY    ( UD_PETSC,  -A5, UDD_OLD_PETSC , IERR )
!---------- ---------- ---------- ---------- ----------


          END DO


CALL VecDestroy    ( R_PETSC      , IERR )
CALL VecDestroy    ( U_PETSC      , IERR )
CALL VecDestroy    ( UD_PETSC     , IERR )
CALL VecDestroy    ( UDD_PETSC    , IERR )
CALL VecDestroy    ( DU_PETSC     , IERR )
CALL VecDestroy    ( UDD_OLD_PETSC, IERR )
CALL VecDestroy    ( HELP_PETSC   , IERR )


! SuperLU_Dist direct solver (1) vs Krylov iterative solver (2)
!---------- ---------- ---------- ---------- ---------- ---------- ----------
SELECT CASE (IFLAG_SOLVER)


   CASE(1)
! SuperLu_Dist solver 
!---------- ---------- ---------- ---------- ----------
      CALL MatDestroy ( AF_PETSC, IERR )
      CALL ISDestroy  ( perm    , IERR )


   CASE(2)
! KSP solver 
!---------- ---------- ---------- ---------- ----------
      CALL KSPDestroy ( ksp, IERR )

      IF (IFLAG_PRECONDITIONER == 2) THEN
         CALL ISDestroy ( is_u, IERR )
         CALL ISDestroy ( is_s, IERR )
      END IF


END SELECT
!---------- ---------- ---------- ---------- ---------- ---------- ----------


RETURN
END


!************************************************
! NEWMARK: implicit time-stepping - diagonal mass (to check if diagonalization leads to instability; Newmark is stable)
!************************************************
! REVISION : F 10 Feb 2012

SUBROUTINE NEWMARK_PETSC_2 ( AK_PETSC, DIAG_M_PETSC, AC_PETSC, B_PETSC, AK_RD_PETSC, AM_RD_PETSC, ID, LTRANS, RMJ_PETSC, BACL, RANK )
USE PARAMETERS
USE IFPORT
IMPLICIT DOUBLE PRECISION (A-H,O-Z)

character(30)     solver_package

!==================================================================================================
! - PETSC INCLUDE FILES
#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscmat.h"
#include "finclude/petscksp.h"
#include "finclude/petscpc.h"
#include "finclude/petscis.h"

! - PETSC VARIABLES AND MATRICES
Mat            :: AK_PETSC, AC_PETSC, AK_RD_PETSC, AM_RD_PETSC
Vec            :: B_PETSC
Vec            :: RMJ_PETSC
Vec            :: DIAG_M_PETSC

KSP            :: ksp
PC             :: pc
KSPConvergedReason :: reason

IS             :: perm
MatFactorInfo  :: mat_factor_info(MAT_FACTORINFO_SIZE)

PetscReal      :: rtol
PetscInt       :: its, m

PetscErrorCode :: IERR
!PetscTruth     :: FLAG
PetscMPIInt    :: SIZE, RANK
PetscScalar    :: ENERGY

PetscScalar    :: zero, one, minus_one
PetscScalar    :: A0, A1, A2, A3, A4, A5
PetscScalar    :: AG_X

! - PETSc local vectors
Vec            :: R_PETSC, U_PETSC, UD_PETSC, UDD_PETSC, DU_PETSC, UDD_OLD_PETSC
Vec            :: HELP_PETSC
!==================================================================================================


! in 
!---------- ---------- -------- -- ---------- ----------
DIMENSION BACL(NSTEP,NDIM)
DIMENSION ID(NJ, NDOF), LTRANS(NTRANS)
!---------- ---------- ---------- ---------- ----------


! local
!---------- ---------- ---------- ---------- ----------
DIMENSION AG(NDIM)
!---------- ---------- ---------- ---------- ----------


zero      =  0.0d0
one       =  1.0d0
minus_one = -1.0d0
rtol      =  1.d-7
its       =  0


! (1) earthquake analysis        vs  (2) specified external load
! (1) SuperLU_Dist direct solver vs  (2) Krylov iterative solver
! (1) direct solver factor LU    vs  (2) Cholesky <symmetric>
!---------- ---------- ---------- ---------- ----------
IFLAG_LOAD   = 2
IFLAG_SOLVER = 1
IFLAG_FACTOR = 1


!==================================================================================================
! - DEFINE PETSC OBJECTS
CALL VecCreateMPI ( PETSC_COMM_WORLD, PETSC_DECIDE, NEQ, R_PETSC      , IERR)
CALL VecCreateMPI ( PETSC_COMM_WORLD, PETSC_DECIDE, NEQ, U_PETSC      , IERR)
CALL VecCreateMPI ( PETSC_COMM_WORLD, PETSC_DECIDE, NEQ, UD_PETSC     , IERR)
CALL VecCreateMPI ( PETSC_COMM_WORLD, PETSC_DECIDE, NEQ, UDD_PETSC    , IERR)
CALL VecCreateMPI ( PETSC_COMM_WORLD, PETSC_DECIDE, NEQ, DU_PETSC     , IERR)
CALL VecCreateMPI ( PETSC_COMM_WORLD, PETSC_DECIDE, NEQ, UDD_OLD_PETSC, IERR)
CALL VecCreateMPI ( PETSC_COMM_WORLD, PETSC_DECIDE, NEQ, HELP_PETSC   , IERR)

CALL VecSet       ( R_PETSC  , zero, ierr)
CALL VecSet       ( U_PETSC  , zero, ierr)
CALL VecSet       ( UD_PETSC , zero, ierr)
CALL VecSet       ( UDD_PETSC, zero, ierr)
CALL VecSet       ( DU_PETSC , zero, ierr)
CALL VecSet       ( HELP_PETSC,zero, ierr)

! factored matrix for direct solver; see if you can factor in place and get rid of this matrix
CALL MatCreateMPIAIJ ( PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, NEQ, NEQ, 100, PETSC_NULL_INTEGER, 100, PETSC_NULL_INTEGER, AF_PETSC, IERR )
!==================================================================================================


! Newmark parameters, initialize
!---------- ---------- ---------- ---------- ----------  
A0 = 1.0d0 / (ALPHA * DT**2)
A1 = DELTA / (ALPHA * DT)
A2 = 1.0d0 / (ALPHA * DT)
A3 = 1.0d0 / (2.0d0 * ALPHA) - 1.0d0
A4 = DELTA / ALPHA - 1.0d0
A5 = (DELTA / (2.0d0 * ALPHA) - 1.0d0) * DT
!---------- ---------- ---------- ---------- ---------- 
T  = 0.0d0
!---------- ---------- ---------- ---------- ---------- 


! form the effective "stiffness matrix"
!---------- ---------- ---------- ---------- ----------
! AK_f  = A0 * AM_f  +  A1 * AC_f  +  AK_f
!CALL MatAXPY ( AK_PETSC, A0, AM_PETSC, SAME_NONZERO_PATTERN, IERR )
!CALL MatAXPY ( AK_PETSC, A1, AC_PETSC, SAME_NONZERO_PATTERN, IERR )
! this was the global mass matrix:
!CALL MatAXPY ( AK_PETSC, A0, AM_PETSC, DIFFERENT_NONZERO_PATTERN, IERR )
! we replace it with a diagonal mass:
CALL VecAXPY ( HELP_PETSC, A0, DIAG_M_PETSC, IERR )
CALL MatDiagonalSet ( AK_PETSC, HELP_PETSC, ADD_VALUES, IERR )
CALL MatAXPY ( AK_PETSC, A1, AC_PETSC, DIFFERENT_NONZERO_PATTERN, IERR )


! SuperLU_Dist direct solver(1) vs Krylov iterative solver(2)
! for direct solver, perform matrix factorization, for iterative solver, setup appropriate solver and preconditioner
!---------- ---------- ---------- ---------- ---------- ---------- ---------- ---------- ----------
CALL CPU_TIME(t0)
SELECT CASE (IFLAG_SOLVER)

   CASE(1)
! SuperLu_Dist solver 
!---------- ---------- ---------- ---------- ----------
! LU factorization
      SELECT CASE (IFLAG_FACTOR)

         CASE(1)

            WRITE(*,*) ' LU factorization begins ...'
            CALL MatGetFactor        ( AK_PETSC, "superlu_dist", MAT_FACTOR_LU, AF_PETSC, IERR )
            CALL MatGetSize          ( AK_PETSC, m, PETSC_NULL, IERR )
            CALL ISCreateStride      ( PETSC_COMM_WORLD, m, 0, 1, perm, IERR )
            CALL MatLUFactorSymbolic ( AF_PETSC, AK_PETSC, perm, perm, mat_factor_info, IERR )
            CALL MatLUFactorNumeric  ( AF_PETSC, AK_PETSC, mat_factor_info, IERR )

! Cholesky factorization
         CASE(2)
write(*,*) ' Cholesky factorization: coming soon!'
stop
!            CALL MatGetFactor              ( AK_PETSC, "superlu_dist", MAT_FACTOR_CHOLESKY, AF_PETSC, IERR )
!            CALL MatGetSize                ( AK_PETSC, m, PETSC_NULL, IERR )
!            CALL ISCreateStride            ( PETSC_COMM_WORLD, m, 0, 1, perm, IERR )
!            CALL MatCholeskyFactorSymbolic ( AF_PETSC, AK_PETSC, perm, mat_factor_info, IERR )
!            CALL MatCholeskyFactorNumeric  ( AF_PETSC, AK_PETSC, mat_factor_info, IERR )

      END SELECT


! spit out what direct solver we used
      CALL MatFactorGetSolverPackage ( AF_PETSC, solver_package, IERR )
      WRITE(*,*) 'Solver Package: ', solver_package
      WRITE(2,*) 'Solver Package: ', solver_package


   CASE(2)
! KSP solver 
!---------- ---------- ---------- ---------- ----------
      CALL KSPCreate         ( PETSC_COMM_WORLD, ksp, IERR )
      CALL KSPSetOperators   ( ksp, AK_PETSC, AK_PETSC, SAME_PRECONDITIONER, IERR )
!     CALL KSPSetType        ( ksp, KSPCG     , IERR)
!     CALL KSPSetType        ( ksp, KSPMINRES , IERR)
!      CALL KSPSetType        ( ksp, KSPGMRES  , IERR)
!      CALL KSPSetInitialGuessNonzero ( ksp, PETSC_TRUE, IERR )

      WRITE(*,*) "Solver Package: Krylov iterative solver"
      WRITE(2,*) "Solver Package: Krylov iterative solver"

! Pre-conditioner
!---------- ---------- ---------- ---------- ----------
      CALL KSPGetPC         ( ksp, pc       , IERR )
!     CALL PCSetType        ( pc,  PCNONE   , IERR )
      CALL PCSetType        ( pc,  PCJACOBI , IERR )
!     CALL PCSetType        ( pc,  PCILU    , IERR )
      CALL PCFactorSetLevels( pc,  3        , IERR )
      CALL KSPSetTolerances ( ksp, PETSC_DEFAULT_DOUBLE_PRECISION, PETSC_DEFAULT_DOUBLE_PRECISION, PETSC_DEFAULT_DOUBLE_PRECISION, PETSC_DEFAULT_INTEGER, IERR)

END SELECT


CALL CPU_TIME(t1)
IF ( RANK == 0 ) WRITE(*,'(A40,F10.2,A10)') 'Factorization / Preconditioning time:', t1-t0, 'seconds'
IF ( RANK == 0 ) WRITE(2,'(A40,F10.2,A10)') 'Factorization / Preconditioning time:', t1-t0, 'seconds'
!---------- ---------- ---------- ---------- ---------- ---------- ---------- ---------- ----------


! NEWMARK INTEGRATION
!---------- ---------- ---------- ---------- ----------
          DO ISTEP = 1 , NSTEP


! Energy in regular domain
!---------- ---------- ---------- ---------- ----------      
!CALL TOTAL_ENERGY_RD_PETSC ( AK_RD_PETSC, AM_RD_PETSC, U_PETSC, UD_PETSC, ENERGY)
!---------- ---------- ---------- ---------- ----------              


! print out
!---------- ---------- ---------- ---------- ----------
CALL OUT_NEWMARK_PETSC ( U_PETSC, ID, LTRANS, T, ENERGY )
!---------- ---------- ---------- ---------- ---------- 


             T = DBLE(ISTEP) * DT
             IF ( RANK == 0 ) WRITE (*,'(A20,F10.5,A20,I10)') 'TIME=', T, 'KSP Iter=', its


! earthquake analysis (1) vs specified external load (2)
!---------- ---------- ---------- ---------- ---------- ---------- ----------
SELECT CASE (IFLAG_LOAD)

   CASE(1)
! earthquake analysis case
! only x-component of acceleration has been considered. see Assemble_PETSc.F90
!---------- ---------- ---------- ---------- ----------
      AG_X = BACL(ISTEP, 1)
      CALL VecSet   ( R_PETSC, zero             , IERR )
      CALL VecAXPY  ( R_PETSC, AG_X, RMJ_PETSC  , IERR )
!---------- ---------- ---------- ---------- ---------- 


   CASE(2)
!  external load case ( forward & inverse problem )
!---------- ---------- ---------- ---------- ---------- 
      CALL TIME_FUNCTION ( T, FAC1 )
      CALL VecSet   ( R_PETSC, zero             , IERR )
      CALL VecAXPY  ( R_PETSC, FAC1, B_PETSC    , IERR )
!---------- ---------- ---------- ---------- ----------


END SELECT
!---------- ---------- ---------- ---------- ---------- ---------- ----------


! Newmark's method
!---------- ---------- ---------- ---------- ----------        
! R = R + MATMUL(AM_f, (A0 * U + A2 * UD + A3 * UDD)) + MATMUL(AC_f, (A1 * U + A4 * UD + A5 * UDD))
CALL VecSet     ( HELP_PETSC, zero         , IERR )
CALL VecAXPY    ( HELP_PETSC, A0, U_PETSC  , IERR )
CALL VecAXPY    ( HELP_PETSC, A2, UD_PETSC , IERR )
CALL VecAXPY    ( HELP_PETSC, A3, UDD_PETSC, IERR )
! this was the global mass matrix:
!CALL MatMultAdd ( AM_PETSC, HELP_PETSC, R_PETSC, R_PETSC, IERR )
! we replace it with a diagonal mass:
CALL VecPointwiseMult ( HELP_PETSC, DIAG_M_PETSC, HELP_PETSC, IERR )
CALL VecAXPY ( R_PETSC, 1.0d0, HELP_PETSC, IERR )

CALL VecSet     ( HELP_PETSC, zero         , IERR )
CALL VecAXPY    ( HELP_PETSC, A1, U_PETSC  , IERR )
CALL VecAXPY    ( HELP_PETSC, A4, UD_PETSC , IERR )
CALL VecAXPY    ( HELP_PETSC, A5, UDD_PETSC, IERR )
CALL MatMultAdd ( AC_PETSC, HELP_PETSC, R_PETSC, R_PETSC, IERR )


! SuperLU_Dist direct solver(1) vs Krylov iterative solver(2)
! solve for the same coefficient matrix but different right hand sides
!---------- ---------- ---------- ---------- ---------- ---------- ----------
SELECT CASE (IFLAG_SOLVER)


   CASE(1)
! SuperLu_Dist solver
!---------- ---------- ---------- ---------- ----------
! solve
      CALL MatSolve   ( AF_PETSC, R_PETSC, HELP_PETSC, IERR )
! transfer
      CALL VecCopy    ( HELP_PETSC, R_PETSC, IERR )


   CASE(2)
! Krylov solver
!---------- ---------- ---------- ---------- ----------
! initial guess from previous iteration
      CALL VecCopy    ( U_PETSC, HELP_PETSC, IERR )
! solve
      CALL KSPSolve   ( ksp, R_PETSC, HELP_PETSC, IERR )
! transfer
      CALL VecCopy    ( HELP_PETSC, R_PETSC, IERR )

! convergence reason
      CALL KSPGetConvergedReason ( ksp, reason, IERR )
      IF ( reason < 0 ) THEN
         CALL PetscPrintf ( PETSC_COMM_WORLD, "Failure to converge \n", IERR )
         STOP
      END IF
! note: this is off by one step in the printout
      CALL KSPGetIterationNumber ( ksp, its, IERR )
!     CALL KSPView    ( ksp, PETSC_VIEWER_STDOUT_WORLD, IERR )


END SELECT
!---------- ---------- ---------- ---------- ---------- ---------- ----------


! Update Newmark parameters for the next iteration
!---------- ---------- ---------- ---------- ----------
! DU      = R - U
CALL VecCopy    ( R_PETSC,        DU_PETSC      , IERR )
CALL VecAXPY    ( DU_PETSC, minus_one, U_PETSC  , IERR )
! U       = R
CALL VecCopy    ( R_PETSC,        U_PETSC       , IERR )
! UDD_OLD = UDD
CALL VecCopy    ( UDD_PETSC,      UDD_OLD_PETSC , IERR )
! UDD     = A0 * DU - A2 * UD - A3 * UDD_OLD
CALL VecSet     ( UDD_PETSC, zero               , IERR )
CALL VecAXPY    ( UDD_PETSC,  A0, DU_PETSC      , IERR )
CALL VecAXPY    ( UDD_PETSC, -A2, UD_PETSC      , IERR )
CALL VecAXPY    ( UDD_PETSC, -A3, UDD_OLD_PETSC , IERR )
! UD      = A1 * DU - A4 * UD - A5 * UDD_OLD
CALL VecScale   ( UD_PETSC,  -A4                , IERR )
CALL VecAXPY    ( UD_PETSC,   A1, DU_PETSC      , IERR )
CALL VecAXPY    ( UD_PETSC,  -A5, UDD_OLD_PETSC , IERR )
!---------- ---------- ---------- ---------- ----------


          END DO


CALL VecDestroy    ( R_PETSC      , IERR )
CALL VecDestroy    ( U_PETSC      , IERR )
CALL VecDestroy    ( UD_PETSC     , IERR )
CALL VecDestroy    ( UDD_PETSC    , IERR )
CALL VecDestroy    ( DU_PETSC     , IERR )
CALL VecDestroy    ( UDD_OLD_PETSC, IERR )
CALL VecDestroy    ( HELP_PETSC   , IERR )
CALL VecDestroy    ( DIAG_M_PETSC , IERR )

! SuperLU_Dist direct solver (1) vs Krylov iterative solver (2)
!---------- ---------- ---------- ---------- ---------- ---------- ----------
SELECT CASE (IFLAG_SOLVER)


   CASE(1)
! SuperLu_Dist solver 
!---------- ---------- ---------- ---------- ----------
      CALL MatDestroy ( AF_PETSC, IERR )
      CALL ISDestroy  ( perm    , IERR )


   CASE(2)
! KSP solver 
!---------- ---------- ---------- ---------- ----------
      CALL KSPDestroy ( ksp     , IERR )


END SELECT
!---------- ---------- ---------- ---------- ---------- ---------- ----------


RETURN
END
