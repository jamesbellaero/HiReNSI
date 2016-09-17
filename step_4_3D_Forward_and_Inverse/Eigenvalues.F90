!************************************************
! eigenvalue computation: fully functional for small problems
!************************************************
! REVISION : M 6 Feb 2012

SUBROUTINE Eigen ( AK_PETSC, AC_PETSC, DIAG_M_PETSC, RANK, SIZE )
USE PARAMETERS
USE IFPORT
IMPLICIT DOUBLE PRECISION (A-H,O-Z)

!==================================================================================================
! - PETSC INCLUDE FILES
#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscmat.h"
#include "finclude/petscksp.h"
#include "finclude/petscpc.h"
#include "finclude/petscis.h"
#include "finclude/petscis.h90"


! - PETSC VARIABLES AND MATRICES
Mat            :: AK_PETSC, AC_PETSC
Vec            :: DIAG_M_PETSC, v11_PETSC, v12_PETSC, v21_PETSC, v22_PETSC

Mat            :: Ab_PETSC, I_PETSC
Vec            :: xb_PETSC, rb_PETSC

KSP            :: ksp
PC             :: pc

PetscErrorCode :: IERR
PetscMPIInt    :: SIZE, RANK

PetscReal      :: r( 2 * NEQ ), c( 2 * NEQ )

IS             :: is_1, is_2
!==================================================================================================


zero      =  0.0d0
one       =  1.0d0
minus_one = -1.0d0


!==================================================================================================
! - DEFINE PETSC OBJECTS
CALL MatCreateAIJ  ( PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, NEQ, NEQ, 1, PETSC_NULL_INTEGER, 1, PETSC_NULL_INTEGER, I_PETSC, IERR )
     DO I = 1, NEQ
        CALL MatSetValue ( I_PETSC, I-1, I-1, one, ADD_VALUES, IERR )
     END DO
CALL MatAssemblyBegin ( I_PETSC, MAT_FINAL_ASSEMBLY, IERR )
CALL MatAssemblyEnd   ( I_PETSC, MAT_FINAL_ASSEMBLY, IERR )
CALL MatAssemblyBegin ( I_PETSC, MAT_FINAL_ASSEMBLY, IERR )
CALL MatAssemblyEnd   ( I_PETSC, MAT_FINAL_ASSEMBLY, IERR )

CALL VecCreateMPI ( PETSC_COMM_WORLD, PETSC_DECIDE, NEQ, v11_PETSC, IERR )
CALL VecCreateMPI ( PETSC_COMM_WORLD, PETSC_DECIDE, NEQ, v12_PETSC, IERR )
CALL VecCreateMPI ( PETSC_COMM_WORLD, PETSC_DECIDE, NEQ, v21_PETSC, IERR )
CALL VecCreateMPI ( PETSC_COMM_WORLD, PETSC_DECIDE, NEQ, v22_PETSC, IERR )
CALL VecSet ( v11_PETSC, zero, IERR )
CALL VecSet ( v12_PETSC, zero, IERR )
CALL VecSet ( v21_PETSC, one , IERR )
CALL VecSet ( v22_PETSC, one , IERR )
!==================================================================================================
! block objects
CALL ISCreateStride ( PETSC_COMM_WORLD, NEQ, 0  , 1, is_1, IERR )
CALL ISCreateStride ( PETSC_COMM_WORLD, NEQ, NEQ, 1, is_2, IERR )
CALL VecCreateNest  ( PETSC_COMM_WORLD, 2, [is_1, is_2], [v11_PETSC, v12_PETSC], xb_PETSC, IERR )
CALL VecCreateNest  ( PETSC_COMM_WORLD, 2, [is_1, is_2], [v21_PETSC, v22_PETSC], rb_PETSC, IERR )
!==================================================================================================
! since vec-nest is not working for xb, rb:
CALL VecCreateMPI ( PETSC_COMM_WORLD, PETSC_DECIDE, 2*NEQ, x_B_PETSC, IERR )
CALL VecCreateMPI ( PETSC_COMM_WORLD, PETSC_DECIDE, 2*NEQ, r_B_PETSC, IERR )
CALL VecSet ( x_B_PETSC, one , IERR )
CALL VecSet ( r_B_PETSC, one , IERR )
!==================================================================================================


! inverse of the diagonal mass matrix
!---------- ---------- ---------- ---------- ----------  
CALL VecReciprocal ( DIAG_M_PETSC, IERR )
CALL MatDiagonalScale ( AK_PETSC, DIAG_M_PETSC, PETSC_NULL, IERR )
CALL MatDiagonalScale ( AC_PETSC, DIAG_M_PETSC, PETSC_NULL, IERR )
CALL MatScale ( AK_PETSC, -1.0D0, IERR )
CALL MatScale ( AC_PETSC, -1.0D0, IERR )
CALL MatCreateNest ( PETSC_COMM_WORLD, 2, [is_1, is_2], 2, [is_1, is_2], [PETSC_NULL_OBJECT, I_PETSC, AK_PETSC, AC_PETSC], Ab_PETSC, IERR )


! assembly
!---------- ---------- ---------- ---------- ---------- 
CALL MatAssemblyBegin ( Ab_PETSC, MAT_FINAL_ASSEMBLY, IERR )
CALL MatAssemblyEnd   ( Ab_PETSC, MAT_FINAL_ASSEMBLY, IERR )
CALL VecAssemblyBegin ( xb_PETSC, IERR )
CALL VecAssemblyEnd   ( xb_PETSC, IERR ) 
CALL VecAssemblyBegin ( rb_PETSC, IERR )
CALL VecAssemblyEnd   ( rb_PETSC, IERR )


!---------- ---------- ---------- ---------- ----------
CALL KSPCreate       ( PETSC_COMM_WORLD, ksp, IERR )
CALL KSPSetOperators ( ksp, Ab_PETSC, Ab_PETSC, SAME_PRECONDITIONER, IERR )
CALL KSPSetType      ( ksp, KSPGMRES, IERR)
CALL PetscPrintf     ( PETSC_COMM_WORLD, "2D Eigenvalue computation ... \n", IERR )
CALL KSPGetPC        ( ksp, pc, IERR )
CALL PCSetType       ( pc, PCNONE, IERR )
CALL KSPSetTolerances( ksp, PETSC_DEFAULT_DOUBLE_PRECISION, PETSC_DEFAULT_DOUBLE_PRECISION, PETSC_DEFAULT_DOUBLE_PRECISION, 20, IERR)

!CALL KSPSetComputeEigenvalues( ksp, PETSC_TRUE, IERR )

CALL KSPSetUp        ( ksp, IERR )
CALL PetscPrintf     ( PETSC_COMM_WORLD, "Before solve ... \n", IERR )
CALL KSPSolve        ( ksp, x_B_PETSC, r_B_PETSC, IERR )
CALL PetscPrintf     ( PETSC_COMM_WORLD, "After solve ... \n", IERR )

CALL MatGetSize      ( Ab_PETSC, m, PETSC_NULL, IERR )
CALL KSPComputeEigenvaluesExplicitly( ksp, m, r, c, IERR ) 
!---------- ---------- ---------- ---------- ----------

! this is not working. use commandline option: -ksp_compute_eigenvalues
!CALL KSPComputeEigenvalues ( ksp, m, r, c, 10, IERR )


CALL VecDestroy ( v11_PETSC, IERR )
CALL VecDestroy ( v12_PETSC, IERR )
CALL VecDestroy ( v21_PETSC, IERR )
CALL VecDestroy ( v22_PETSC, IERR )
CALL VecDestroy (  xb_PETSC, IERR )
CALL VecDestroy (  rb_PETSC, IERR )

CALL VecDestroy (  x_B_PETSC, IERR )
CALL VecDestroy (  r_B_PETSC, IERR )

CALL MatDestroy (   I_PETSC, IERR )
CALL MatDestroy (  Ab_PETSC, IERR )

CALL KSPDestroy ( ksp, IERR )

J = 0
DO I = 1, 2 * NEQ
   IF ( r(I) >= 1.0D-12 ) THEN
      J = J + 1
      WRITE(* ,'(2I10,2E15.5)') J, I, r(I), c(I)
      WRITE(95,'(2I10,2E15.5)') J, I, r(I), c(I)
   END IF
END DO

RETURN
END


!************************************************
! eigenvalue computation: fully functional for small problems
!************************************************
! REVISION : Sat 20 April 2013: after the tea party!

SUBROUTINE Eigen_3D ( AK_PETSC, AC_PETSC, AG_PETSC, DIAG_M_PETSC, RANK, SIZE )
USE PARAMETERS
USE IFPORT
IMPLICIT DOUBLE PRECISION (A-H,O-Z)

!==================================================================================================
! - PETSC INCLUDE FILES
#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscmat.h"
#include "finclude/petscksp.h"
#include "finclude/petscpc.h"
#include "finclude/petscis.h"
#include "finclude/petscis.h90"


! - PETSC VARIABLES AND MATRICES
Mat            :: AK_PETSC, AC_PETSC, AG_PETSC
Vec            :: DIAG_M_PETSC, v11_PETSC, v12_PETSC, v13_PETSC, v21_PETSC, v22_PETSC, v23_PETSC

Mat            :: Ab_PETSC,  I_PETSC
Vec            :: xb_PETSC, rb_PETSC

KSP            :: ksp
PC             :: pc

PetscErrorCode :: IERR
PetscMPIInt    :: SIZE, RANK

PetscReal      :: r( 3 * NEQ ), c( 3 * NEQ )

IS             :: is_1, is_2, is_3
!==================================================================================================

!==================================================================================================
! - DEFINE PETSC OBJECTS
CALL MatCreateAIJ  ( PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, NEQ, NEQ, 1, PETSC_NULL_INTEGER, 1, PETSC_NULL_INTEGER, I_PETSC, IERR )
     DO I = 1, NEQ
        CALL MatSetValue ( I_PETSC, I-1, I-1, 1.0d0, ADD_VALUES, IERR )
     END DO
CALL MatAssemblyBegin ( I_PETSC, MAT_FINAL_ASSEMBLY, IERR )
CALL MatAssemblyEnd   ( I_PETSC, MAT_FINAL_ASSEMBLY, IERR )

CALL VecCreateMPI ( PETSC_COMM_WORLD, PETSC_DECIDE, NEQ, v11_PETSC, IERR )
CALL VecCreateMPI ( PETSC_COMM_WORLD, PETSC_DECIDE, NEQ, v12_PETSC, IERR )
CALL VecCreateMPI ( PETSC_COMM_WORLD, PETSC_DECIDE, NEQ, v13_PETSC, IERR )
CALL VecCreateMPI ( PETSC_COMM_WORLD, PETSC_DECIDE, NEQ, v21_PETSC, IERR )
CALL VecCreateMPI ( PETSC_COMM_WORLD, PETSC_DECIDE, NEQ, v22_PETSC, IERR )
CALL VecCreateMPI ( PETSC_COMM_WORLD, PETSC_DECIDE, NEQ, v23_PETSC, IERR )
CALL VecSet ( v11_PETSC, 0.0d0, IERR )
CALL VecSet ( v12_PETSC, 0.0d0, IERR )
CALL VecSet ( v13_PETSC, 0.0d0, IERR )

CALL VecSet ( v21_PETSC, 1.0d0, IERR )
CALL VecSet ( v22_PETSC, 1.0d0, IERR )
CALL VecSet ( v23_PETSC, 1.0d0, IERR )
!==================================================================================================
! block objects
CALL ISCreateStride ( PETSC_COMM_WORLD, NEQ,   0  , 1, is_1, IERR )
CALL ISCreateStride ( PETSC_COMM_WORLD, NEQ,   NEQ, 1, is_2, IERR )
CALL ISCreateStride ( PETSC_COMM_WORLD, NEQ, 2*NEQ, 1, is_3, IERR )
CALL VecCreateNest  ( PETSC_COMM_WORLD, 3, [is_1, is_2, is_3], [v11_PETSC, v12_PETSC, v13_PETSC], xb_PETSC, IERR )
CALL VecCreateNest  ( PETSC_COMM_WORLD, 3, [is_1, is_2, is_3], [v21_PETSC, v22_PETSC, v23_PETSC], rb_PETSC, IERR )
!==================================================================================================
! since vec-nest is not working for xb, rb:
CALL VecCreateMPI ( PETSC_COMM_WORLD, PETSC_DECIDE, 3*NEQ, x_B_PETSC, IERR )
CALL VecCreateMPI ( PETSC_COMM_WORLD, PETSC_DECIDE, 3*NEQ, r_B_PETSC, IERR )
CALL VecSet ( x_B_PETSC, 1.0d0 , IERR )
CALL VecSet ( r_B_PETSC, 1.0d0 , IERR )
!==================================================================================================


! inverse of the diagonal mass matrix
!---------- ---------- ---------- ---------- ----------  
CALL VecReciprocal ( DIAG_M_PETSC, IERR )
CALL MatDiagonalScale ( AK_PETSC, DIAG_M_PETSC, PETSC_NULL, IERR )
CALL MatDiagonalScale ( AC_PETSC, DIAG_M_PETSC, PETSC_NULL, IERR )
CALL MatDiagonalScale ( AG_PETSC, DIAG_M_PETSC, PETSC_NULL, IERR )
CALL MatScale ( AK_PETSC, -1.0D0, IERR )
CALL MatScale ( AC_PETSC, -1.0D0, IERR )
CALL MatScale ( AG_PETSC, -1.0D0, IERR )
CALL MatCreateNest ( PETSC_COMM_WORLD, 3, [is_1, is_2, is_3], 3, [is_1, is_2, is_3], [PETSC_NULL_OBJECT, I_PETSC          , PETSC_NULL_OBJECT, &
                                                                                      PETSC_NULL_OBJECT, PETSC_NULL_OBJECT, I_PETSC          , &                                                                                                     AG_PETSC         , AK_PETSC         , AC_PETSC]        , Ab_PETSC, IERR )

! assembly
!---------- ---------- ---------- ---------- ---------- 
CALL MatAssemblyBegin ( Ab_PETSC, MAT_FINAL_ASSEMBLY, IERR )
CALL MatAssemblyEnd   ( Ab_PETSC, MAT_FINAL_ASSEMBLY, IERR )
CALL VecAssemblyBegin ( xb_PETSC, IERR )
CALL VecAssemblyEnd   ( xb_PETSC, IERR )
CALL VecAssemblyBegin ( rb_PETSC, IERR )
CALL VecAssemblyEnd   ( rb_PETSC, IERR )


!---------- ---------- ---------- ---------- ----------
CALL KSPCreate       ( PETSC_COMM_WORLD, ksp, IERR )
CALL KSPSetOperators ( ksp, Ab_PETSC, Ab_PETSC, SAME_PRECONDITIONER, IERR )
CALL KSPSetType      ( ksp, KSPGMRES, IERR)
CALL PetscPrintf     ( PETSC_COMM_WORLD, "3D Eigenvalue computation ... \n", IERR )
CALL KSPGetPC        ( ksp, pc, IERR )
CALL PCSetType       ( pc, PCNONE, IERR )
CALL KSPSetTolerances( ksp, PETSC_DEFAULT_DOUBLE_PRECISION, PETSC_DEFAULT_DOUBLE_PRECISION, PETSC_DEFAULT_DOUBLE_PRECISION, 20, IERR)

CALL KSPSetUp        ( ksp, IERR )
CALL PetscPrintf     ( PETSC_COMM_WORLD, "Before solve ... \n", IERR )
CALL KSPSolve        ( ksp, x_B_PETSC, r_B_PETSC, IERR )
CALL PetscPrintf     ( PETSC_COMM_WORLD, "After solve ... \n", IERR )

CALL MatGetSize      ( Ab_PETSC, m, PETSC_NULL, IERR )
CALL KSPComputeEigenvaluesExplicitly( ksp, m, r, c, IERR )
!---------- ---------- ---------- ---------- ----------

! this is not working. use commandline option: -ksp_compute_eigenvalues
!CALL KSPComputeEigenvalues ( ksp, m, r, c, 10, IERR )


CALL VecDestroy ( v11_PETSC, IERR )
CALL VecDestroy ( v12_PETSC, IERR )
CALL VecDestroy ( v13_PETSC, IERR )
CALL VecDestroy ( v21_PETSC, IERR )
CALL VecDestroy ( v22_PETSC, IERR )
CALL VecDestroy ( v23_PETSC, IERR )
CALL VecDestroy (  xb_PETSC, IERR )
CALL VecDestroy (  rb_PETSC, IERR )

CALL VecDestroy (  x_B_PETSC, IERR )
CALL VecDestroy (  r_B_PETSC, IERR )

CALL MatDestroy (   I_PETSC, IERR )
CALL MatDestroy (  Ab_PETSC, IERR )

CALL KSPDestroy ( ksp, IERR )

J = 0
DO I = 1, 3 * NEQ
   IF ( r(I) >= 1.0D-12 ) THEN
      J = J + 1
      WRITE(* ,'(2I10,2E15.5)') J, I, r(I), c(I)
      WRITE(95,'(2I10,2E15.5)') J, I, r(I), c(I)
   END IF
END DO

RETURN
END
