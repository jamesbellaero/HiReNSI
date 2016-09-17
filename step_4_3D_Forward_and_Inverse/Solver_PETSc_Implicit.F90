!************************************************
! NEWMARK: implicit time-stepping - basic preconditioning
!************************************************
! REVISION : M 22 Oct 2012

SUBROUTINE NEWMARK_PETSC_0 ( AK_PETSC, AM_PETSC, AC_PETSC, B_PETSC, AK_RD_PETSC, AM_RD_PETSC, ID, LTRANS, RMJ_PETSC, BACL, RANK, SIZE, XYZ, INOD, NNDH, NNVH, NNAH, EqDis, EqVel, EqAcc )

USE PARAMETERS
Use Results
USE IFPORT
IMPLICIT DOUBLE PRECISION (A-H,O-Z)

Integer, Intent(In), Dimension ( NNDH * NDim )  :: EqDis
Integer, Intent(In), Dimension ( NNVH * NDim )  :: EqVel
Integer, Intent(In), Dimension ( NNAH * NDim )  :: EqAcc
character(30)     solver_package

!==================================================================================================
! - PETSC INCLUDE FILES
#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"
#include "finclude/petscksp.h"
#include "finclude/petscpc.h"
#include "finclude/petscis.h"

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

PetscErrorCode :: IERR, ErrPTC
PetscMPIInt    :: SIZE, RANK
PetscScalar    :: ENERGY

PetscScalar    :: zero, one, minus_one
PetscScalar    :: A0, A1, A2, A3, A4, A5
PetscScalar    :: AG_X

! - PETSc local vectors
Vec            :: R_PETSC, U_PETSC, UD_PETSC, UDD_PETSC, DU_PETSC, UDD_OLD_PETSC
Vec            :: HELP_PETSC

PetscScalar, pointer    :: U(:), UD(:), UDD(:)
!==================================================================================================


! in 
!---------- ---------- -------- -- ---------- ----------
DIMENSION BACL(NSTEP,NDIM)
DIMENSION ID(NJ, NDOF), LTRANS(NTRANS), XYZ(NJ, NDIM), INOD(NNODE, NEL)
!---------- ---------- ---------- ---------- ----------


! local
!---------- ---------- ---------- ---------- ----------
DIMENSION AG(NDIM)
!---------- ---------- ---------- ---------- ----------


I_PRINT_STEP  =  10
zero          =  0.0d0
one           =  1.0d0
minus_one     = -1.0d0
rtol          =  1.d-7
its           =  0
T             =  0.0d0


! (1) earthquake analysis        vs  (2) specified external load
! (1) SuperLU_Dist direct solver vs  (2) Krylov iterative solver
! (1) direct solver factor LU    vs  (2) Cholesky <symmetric>
!---------- ---------- ---------- ---------- ----------
IFLAG_LOAD    = 2
IFLAG_SOLVER  = 1
IFLAG_FACTOR  = 1


!==================================================================================================
! - DEFINE PETSC OBJECTS
CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, R_PETSC      , IERR)
CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, U_PETSC      , IERR)
CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, UD_PETSC     , IERR)
CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, UDD_PETSC    , IERR)
CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, DU_PETSC     , IERR)
CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, UDD_OLD_PETSC, IERR)
CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, HELP_PETSC   , IERR)

CALL VecSet       ( R_PETSC  , 0.0d0, ierr)
CALL VecSet       ( U_PETSC  , 0.0d0, ierr)
CALL VecSet       ( UD_PETSC , 0.0d0, ierr)
CALL VecSet       ( UDD_PETSC, 0.0d0, ierr)
CALL VecSet       ( DU_PETSC , 0.0d0, ierr)

! factored matrix for direct solver; see if you can factor in place and get rid of this matrix
CALL MatCreateAIJ ( PETSC_COMM_WORLD, NEQM_Mapping, NEQM_Mapping, NEQMTotal, NEQMTotal, 200, PETSC_NULL_INTEGER, 200, PETSC_NULL_INTEGER,  AF_PETSC,    IERR ) ;
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
! SuperLU_dist solver 
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
            WRITE(*,*) ' Cholesky factorization: coming soon!'
            STOP

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
!     CALL KSPSetType        ( ksp, KSPGMRES  , IERR)
!     CALL KSPSetInitialGuessNonzero ( ksp, PETSC_TRUE, IERR )

      WRITE(*,*) "Solver Package: Krylov iterative solver"
      WRITE(2,*) "Solver Package: Krylov iterative solver"


! Pre-conditioner
!---------- ---------- ---------- ---------- ----------
      CALL KSPGetPC         ( ksp, pc       , IERR )
!     CALL PCSetType        ( pc,  PCNONE   , IERR )
!     CALL PCSetType        ( pc,  PCJACOBI , IERR )
      CALL PCSetType        ( pc,  PCILU    , IERR )
      CALL PCFactorSetLevels( pc,  3        , IERR )
      CALL KSPSetTolerances ( ksp, PETSC_DEFAULT_DOUBLE_PRECISION, PETSC_DEFAULT_DOUBLE_PRECISION, PETSC_DEFAULT_DOUBLE_PRECISION, PETSC_DEFAULT_INTEGER, IERR)

 END SELECT
!---------- ---------- ---------- ---------- ----------


CALL CPU_TIME(t1)
IF ( RANK == 0 ) WRITE(*,'(A40,F10.2,A10)') 'Factorization / Preconditioning time:', t1-t0, 'seconds'
IF ( RANK == 0 ) WRITE(2,'(A40,F10.2,A10)') 'Factorization / Preconditioning time:', t1-t0, 'seconds'
!---------- ---------- ---------- ---------- ---------- ---------- ---------- ---------- ----------


! print out & Visualization
!---------- ---------- ---------- ---------- ---------- ---------------------------------------------------------------------------------------------
Call VecGetOwnershipRange ( U_PETSC, IStart0, IEnd0, ErrPTC ) ;
IStart = IStart0
IEnd   = IEnd0
!---------- ---------- ---------- ---------- ---------- ---------------------------------------------------------------------------------------------


! NEWMARK INTEGRATION
!---------- ---------- ---------- ---------- ----------
DO ISTEP = 1 , NSTEP


! Energy in regular domain
!---------- ---------- ---------- ---------- ----------      
   CALL TOTAL_ENERGY_RD_PETSC ( AK_RD_PETSC, AM_RD_PETSC, U_PETSC, UD_PETSC, ENERGY )
!---------- ---------- ---------- ---------- ----------              


! print out & Visualization
!---------- ---------- ---------- ---------- ---------- ---------------------------------------------------------------------------------------------

    If ( NNDH /= 0 ) Call VecGetArrayF90     (    U_PETSC,   U,  ErrPTC ) ; ! NNDH
    If ( NNVH /= 0 ) Call VecGetArrayF90     (   UD_PETSC,  UD,  ErrPTC ) ; ! NNVH
    If ( NNAH /= 0 ) Call VecGetArrayF90     (  UDD_PETSC, UDD,  ErrPTC ) ; ! NNAH

    ! Writing down the history of displacement, velocity and acceleration
    Call ResDyn   ( NDIM,     NNDH, NNVH, NNAH, IStep,     IStart, IEnd,     T, Energy,     EqDis, EqVel, EqAcc,    U, UD, UDD ) ;  !?? see why interface does not working
    !Call ParaView ( NDim,     IStep,            NJ, NEl,     ELT,     INod, ID,      U,     XYZ,     Name, OutDir   ) ;
!    If ( Any( U /= 0.0_dbl ) .and. Mod(IStep,500)==0 ) Call ParaView_9noded ( NDim,     IStep,            NJ, NEl,     ELT,     INod, ID,      U, UD, UDD,    XYZ,     Name, OutDir   ) ;

!### this must be changed to paraview_9node

    If ( NNDH /= 0 ) Call VecRestoreArrayF90 (    U_PETSC,   U,  ErrPTC ) ; ! NNDH
    If ( NNVH /= 0 ) Call VecRestoreArrayF90 (   UD_PETSC,  UD,  ErrPTC ) ; ! NNVH
    If ( NNAH /= 0 ) Call VecRestoreArrayF90 (  UDD_PETSC, UDD,  ErrPTC ) ; ! NNAH


   CALL OUT_NEWMARK_PETSC ( U_PETSC, ID, LTRANS, T, ENERGY )

   IF ( MOD ( ISTEP - 1 , 1 * I_PRINT_STEP ) == 0 ) THEN
!      CALL ParaView ( XYZ, INOD, ID, U_PETSC, ISTEP )
   END IF
!---------- ---------- ---------- ---------- ---------- ---------------------------------------------------------------------------------------------



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


idebug = 9
if (idebug == 95) then
    write(95,'(2f10.5)') T, FAC1
end if


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
!       CALL KSPView    ( ksp, PETSC_VIEWER_STDOUT_WORLD, IERR )


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
!---------- ---------- ---------- ---------- ----------

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
      CALL KSPDestroy ( ksp     , IERR )

END SELECT
!---------- ---------- ---------- ---------- ---------- ---------- ----------


RETURN
END


!************************************************
! extended NEWMARK: implicit time-stepping (for 3rd order ODEs)
! M u_ddd + C u_dd + K u_d + G u = f
!************************************************
! REVISION : Th 2 May 2013

SUBROUTINE NEWMARK_PETSC_3 ( AK_PETSC, AM_PETSC, AC_PETSC, AG_PETSC, B_PETSC, AK_RD_PETSC, AM_RD_PETSC, ID, LTRANS, RANK, XYZ, INOD, NNDH, NNVH, NNAH, EqDis, EqVel, EqAcc )
USE PARAMETERS
Use Results
USE IFPORT
IMPLICIT DOUBLE PRECISION (A-H,O-Z)

Integer, Intent(In), Dimension ( NNDH * NDim )  :: EqDis
Integer, Intent(In), Dimension ( NNVH * NDim )  :: EqVel
Integer, Intent(In), Dimension ( NNAH * NDim )  :: EqAcc
character(30)     solver_package

!==================================================================================================
! - PETSC INCLUDE FILES
#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"
#include "finclude/petscis.h"

! - PETSC VARIABLES AND MATRICES
Mat            :: AK_PETSC, AM_PETSC, AC_PETSC, AG_PETSC
Mat            :: AK_RD_PETSC, AM_RD_PETSC
Vec            :: B_PETSC

IS             :: perm
MatFactorInfo  :: mat_factor_info(MAT_FACTORINFO_SIZE)
PetscInt       :: m

PetscErrorCode :: IERR, ErrPTC
PetscMPIInt    :: SIZE, RANK

PetscScalar    :: a0, a1, a2, a3, a4, a5, a6, FAC1

! - PETSc local vectors
Vec            :: R_PETSC, U_PETSC, UD_PETSC, UDD_PETSC, UDDD_PETSC, HELP_PETSC

PetscScalar, pointer    :: U(:), UD(:), UDD(:)
!==================================================================================================


! in 
!---------- ---------- -------- -- ---------- ----------
DIMENSION ID(NJ, NDOF), LTRANS(NTRANS)
DIMENSION XYZ(NJ, NDIM), INOD(NNODE, NEL)
!---------- ---------- ---------- ---------- ----------


!==================================================================================================
! - DEFINE PETSC OBJECTS
CALL VecCreateMPI ( PETSC_COMM_WORLD, PETSC_DECIDE, NEQ,    R_PETSC, IERR )
CALL VecCreateMPI ( PETSC_COMM_WORLD, PETSC_DECIDE, NEQ,    U_PETSC, IERR )
CALL VecCreateMPI ( PETSC_COMM_WORLD, PETSC_DECIDE, NEQ,   UD_PETSC, IERR )
CALL VecCreateMPI ( PETSC_COMM_WORLD, PETSC_DECIDE, NEQ,  UDD_PETSC, IERR )
CALL VecCreateMPI ( PETSC_COMM_WORLD, PETSC_DECIDE, NEQ, UDDD_PETSC, IERR )
CALL VecCreateMPI ( PETSC_COMM_WORLD, PETSC_DECIDE, NEQ, HELP_PETSC, IERR )

CALL VecSet (    R_PETSC, 0.0d0, IERR )
CALL VecSet (    U_PETSC, 0.0d0, IERR )
CALL VecSet (   UD_PETSC, 0.0d0, IERR )
CALL VecSet (  UDD_PETSC, 0.0d0, IERR )
CALL VecSet ( UDDD_PETSC, 0.0d0, IERR )

! factored matrix for direct solver; see if you can factor in place and get rid of this matrix
CALL MatCreateAIJ ( PETSC_COMM_WORLD, NEQM_Mapping, NEQM_Mapping, NEQMTotal, NEQMTotal, 200, PETSC_NULL_INTEGER, 200, PETSC_NULL_INTEGER,  AF_PETSC,    IERR ) ;
! I know, I know, this is stupid and can be fixed (L stands for linear operator).
CALL MatCreateAIJ ( PETSC_COMM_WORLD, NEQM_Mapping, NEQM_Mapping, NEQMTotal, NEQMTotal, 200, PETSC_NULL_INTEGER, 200, PETSC_NULL_INTEGER,  AL_PETSC,    IERR ) ;
!==================================================================================================


! Newmark parameters, initialize
!---------- ---------- ---------- ---------- ----------  
alpha_NM_3 = 1.0d0 / 12.0d0                      ! constant average jerk (Sezgin's dissertation, p. 92)
 beta_NM_3 = 1.0d0 /  4.0d0
gamma_NM_3 = 1.0d0 /  2.0d0

a0 = gamma_NM_3 * dt
a1 =  beta_NM_3 * dt**2
a2 = alpha_NM_3 * dt**3
a3 = ( 1.0d0 - gamma_NM_3 ) * dt
a4 = ( 0.5d0 -  beta_NM_3 ) * dt**2
a5 =   0.5d0 * dt**2
a6 = ( 1.0d0 / 6.0d0 - alpha_NM_3 ) * dt**3
!---------- ---------- ---------- ---------- ---------- 
T  = 0.0d0
!---------- ---------- ---------- ---------- ---------- 


! form the effective "stiffness matrix"
!---------- ---------- ---------- ---------- ----------
CALL MatCopy ( AM_PETSC,     AL_PETSC, DIFFERENT_NONZERO_PATTERN, IERR )
CALL MatAXPY ( AL_PETSC, a0, AC_PETSC, DIFFERENT_NONZERO_PATTERN, IERR )
CALL MatAXPY ( AL_PETSC, a1, AK_PETSC, DIFFERENT_NONZERO_PATTERN, IERR )
CALL MatAXPY ( AL_PETSC, a2, AG_PETSC, DIFFERENT_NONZERO_PATTERN, IERR )


! SuperLU_dist solver 
!---------- ---------- ---------- ---------- ----------
CALL CPU_TIME(t0)
WRITE(*,*) ' LU factorization begins (Extended Newmark)... fasten your seat belts ...'

CALL MatGetFactor        ( AL_PETSC, "superlu_dist", MAT_FACTOR_LU, AF_PETSC, IERR )
CALL MatGetSize          ( AL_PETSC, m, PETSC_NULL_INTEGER, IERR )
CALL ISCreateStride      ( PETSC_COMM_WORLD, m, 0, 1, perm, IERR )
CALL MatLUFactorSymbolic ( AF_PETSC, AL_PETSC, perm, perm, mat_factor_info, IERR )
CALL MatLUFactorNumeric  ( AF_PETSC, AL_PETSC, mat_factor_info, IERR )

CALL MatFactorGetSolverPackage ( AF_PETSC, solver_package, IERR )
WRITE(*,*) 'Solver Package: ', solver_package
WRITE(2,*) 'Solver Package: ', solver_package


CALL CPU_TIME(t1)
IF ( RANK == 0 ) WRITE(*,'(A40,F10.2,A10)') 'Factorization time:', t1-t0, 'seconds'
IF ( RANK == 0 ) WRITE(2,'(A40,F10.2,A10)') 'Factorization time:', t1-t0, 'seconds'
!---------- ---------- ---------- ---------- ---------- ---------- ---------- ---------- ----------


! print out & Visualization
!---------- ---------- ---------- ---------- ---------- ---------------------------------------------------------------------------------------------
Call VecGetOwnershipRange ( U_PETSC, IStart0, IEnd0, ErrPTC ) ;
IStart = IStart0
IEnd   = IEnd0
!---------- ---------- ---------- ---------- ---------- ---------------------------------------------------------------------------------------------


! NEWMARK INTEGRATION
!---------- ---------- ---------- ---------- ----------
DO ISTEP = 1 , NSTEP           


! Energy in regular domain
!---------- ---------- ---------- ---------- ----------      
!   CALL TOTAL_ENERGY_RD_PETSC ( AK_RD_PETSC, AM_RD_PETSC, U_PETSC, UD_PETSC, ENERGY )
   CALL TOTAL_ENERGY_RD_PETSC ( AK_RD_PETSC, AM_RD_PETSC, UD_PETSC, UDD_PETSC, ENERGY )
!---------- ---------- ---------- ---------- ----------              


! print out & Visualization
!---------- ---------- ---------- ---------- ---------- ---------------------------------------------------------------------------------------------

    If ( NNDH /= 0 ) Call VecGetArrayF90     (    U_PETSC,   U,  ErrPTC ) ; ! NNDH
    If ( NNVH /= 0 ) Call VecGetArrayF90     (   UD_PETSC,  UD,  ErrPTC ) ; ! NNVH
    If ( NNAH /= 0 ) Call VecGetArrayF90     (  UDD_PETSC, UDD,  ErrPTC ) ; ! NNAH

    ! Writing down the history of displacement, velocity and acceleration
    Call ResDyn   ( NDIM,     NNDH, NNVH, NNAH, IStep,     IStart, IEnd,     T, Energy,     EqDis, EqVel, EqAcc,    U, UD, UDD ) ;  !?? see why interface does not working
    !Call ParaView ( NDim,     IStep,            NJ, NEl,     ELT,     INod, ID,      U,     XYZ,     Name, OutDir   ) ;
!    If ( Any( U /= 0.0_dbl ) .and. Mod(IStep,500)==0 ) Call ParaView_9noded ( NDim,     IStep,            NJ, NEl,     ELT,     INod, ID,      U, UD, UDD,    XYZ,     Name, OutDir   ) ;

!### this must be changed to paraview_9node

    If ( NNDH /= 0 ) Call VecRestoreArrayF90 (    U_PETSC,   U,  ErrPTC ) ; ! NNDH
    If ( NNVH /= 0 ) Call VecRestoreArrayF90 (   UD_PETSC,  UD,  ErrPTC ) ; ! NNVH
    If ( NNAH /= 0 ) Call VecRestoreArrayF90 (  UDD_PETSC, UDD,  ErrPTC ) ; ! NNAH


   CALL OUT_NEWMARK_PETSC ( U_PETSC, ID, LTRANS, T, ENERGY )

   IF ( MOD ( ISTEP - 1 , 1 * I_PRINT_STEP ) == 0 ) THEN
!      CALL ParaView ( XYZ, INOD, ID, U_PETSC, ISTEP )
   END IF
!---------- ---------- ---------- ---------- ---------- ---------------------------------------------------------------------------------------------



   T = DBLE(ISTEP) * DT
   IF ( RANK == 0 ) WRITE (*,'(A20,F10.5,A20,I10)') 'TIME=', T, 'KSP Iter=', its


!  external load case ( forward & inverse problem )
!---------- ---------- ---------- ---------- ---------- 
   CALL TIME_FUNCTION ( T, FAC1 )
   CALL VecSet  ( R_PETSC, 0.0d0        , IERR )
   CALL VecAXPY ( R_PETSC, FAC1, B_PETSC, IERR )
!---------- ---------- ---------- ---------- ----------


! extended Newmark's method (rhs vector)
!---------- ---------- ---------- ---------- ----------        
   CALL VecSet     ( HELP_PETSC, 0.0d0, IERR )
   CALL VecMAXPY   ( HELP_PETSC, 2, [-1.0d0, -a3], [UDD_PETSC, UDDD_PETSC], IERR )
   CALL MatMultAdd ( AC_PETSC, HELP_PETSC, R_PETSC, R_PETSC, IERR )

   CALL VecSet     ( HELP_PETSC, 0.0d0, IERR )
   CALL VecMAXPY   ( HELP_PETSC, 3, [-1.0d0, -dt, -a4], [UD_PETSC, UDD_PETSC, UDDD_PETSC], IERR )
   CALL MatMultAdd ( AK_PETSC, HELP_PETSC, R_PETSC, R_PETSC, IERR )

   CALL VecSet     ( HELP_PETSC, 0.0d0, IERR )
   CALL VecMAXPY   ( HELP_PETSC, 4, [-1.0d0, -dt, -a5, -a6], [U_PETSC, UD_PETSC, UDD_PETSC, UDDD_PETSC], IERR )
   CALL MatMultAdd ( AG_PETSC, HELP_PETSC, R_PETSC, R_PETSC, IERR )


! Triangular solves ( k_eff u_ddd = rhs )
!---------- ---------- ---------- ---------- ----------
   CALL MatSolve ( AF_PETSC, R_PETSC, HELP_PETSC, IERR )


! Newmark Updates
!---------- ---------- ---------- ---------- ----------
! update u_n+1
   CALL VecMAXPY (   U_PETSC, 4, [dt, a5, a6, a2], [UD_PETSC, UDD_PETSC, UDDD_PETSC, HELP_PETSC], IERR )
! update u_d_n+1
   CALL VecMAXPY (  UD_PETSC, 3,     [dt, a4, a1],           [UDD_PETSC, UDDD_PETSC, HELP_PETSC], IERR )
! update u_dd_n+1
   CALL VecMAXPY ( UDD_PETSC, 2,         [a3, a0],                      [UDDD_PETSC, HELP_PETSC], IERR )
! transfer u_ddd_n+1
   CALL VecCopy  ( HELP_PETSC, UDDD_PETSC, IERR )


END DO
!---------- ---------- ---------- ---------- ----------


CALL VecDestroy (    R_PETSC, IERR )
CALL VecDestroy (    U_PETSC, IERR )
CALL VecDestroy (   UD_PETSC, IERR )
CALL VecDestroy (  UDD_PETSC, IERR )
CALL VecDestroy ( UDDD_PETSC, IERR )
CALL VecDestroy ( HELP_PETSC, IERR )
CALL MatDestroy (   AF_PETSC, IERR )
CALL MatDestroy (   AL_PETSC, IERR )
CALL ISDestroy  (       perm, IERR )


RETURN
END


!************************************************
! 3D-PML NEWMARK: implicit time-stepping (for 3rd order ODEs)
! We re-write the 3rd order system as a 2nd order system and apply the standard Newmark method.
! Should be called from SUBROUTINE Solve_PML_3D_2nd_order_ODE
!************************************************
! REVISION : Th 26 July 2012

SUBROUTINE NEWMARK_PETSC_2nd_order_ODE ( AM_b_PETSC, AC_b_PETSC, AK_b_PETSC, B_b_PETSC, ID, LTRANS, RANK )
USE PARAMETERS
USE IFPORT
IMPLICIT DOUBLE PRECISION (A-H,O-Z)

character(30)     solver_package

!==================================================================================================
! - PETSC INCLUDE FILES
#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscmat.h"
#include "finclude/petscis.h"

! - PETSC VARIABLES AND MATRICES
Mat            :: AM_b_PETSC, AC_b_PETSC, AK_b_PETSC
Vec            ::  B_b_PETSC

Mat            :: AF_PETSC

Vec            :: Xdd_b_PETSC, Xd_b_PETSC, X_b_PETSC, R_b_PETSC
Vec            :: HELP_PETSC, Xdd_OLD_b_PETSC, U_PETSC

IS             :: is_1, is_2, perm
MatFactorInfo  :: mat_factor_info(MAT_FACTORINFO_SIZE)
PetscInt       :: m

PetscErrorCode :: IERR
PetscMPIInt    :: SIZE, RANK

PetscScalar    :: a0, a1, a2, a3, a4, a5, FAC1
!==================================================================================================


! in 
!---------- ---------- -------- -- ---------- ----------
DIMENSION ID(NJ, NDOF), LTRANS(NTRANS)
!---------- ---------- ---------- ---------- ----------


write(*,*) 'NEWMARK_PETSC_2nd_order_ODE should be fixed'
stop

! matrices
!---------- ---------- ---------- ---------- ----------
! factored matrix for direct solver; see if you can factor in place and get rid of this matrix
CALL MatCreateAIJ ( PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, 2*NEQ, 2*NEQ, 100, PETSC_NULL_INTEGER, 100, PETSC_NULL_INTEGER, AF_PETSC, IERR )

!---------- ---------- ---------- ---------- ----------


! vectors
!---------- ---------- ---------- ---------- ----------
CALL VecCreateMPI ( PETSC_COMM_WORLD, PETSC_DECIDE, 2*NEQ,     Xdd_b_PETSC, IERR )
CALL VecCreateMPI ( PETSC_COMM_WORLD, PETSC_DECIDE, 2*NEQ,      Xd_b_PETSC, IERR )
CALL VecCreateMPI ( PETSC_COMM_WORLD, PETSC_DECIDE, 2*NEQ,       X_b_PETSC, IERR )
CALL VecCreateMPI ( PETSC_COMM_WORLD, PETSC_DECIDE, 2*NEQ,       R_b_PETSC, IERR )
CALL VecCreateMPI ( PETSC_COMM_WORLD, PETSC_DECIDE, 2*NEQ, Xdd_OLD_b_PETSC, IERR )
CALL VecCreateMPI ( PETSC_COMM_WORLD, PETSC_DECIDE, 2*NEQ,      HELP_PETSC, IERR )


! system is at rest
!---------- ---------- ---------- ---------- ----------
CALL VecSet ( Xdd_b_PETSC, 0.0d0, IERR )
CALL VecSet (  Xd_b_PETSC, 0.0d0, IERR )
CALL VecSet (   X_b_PETSC, 0.0d0, IERR )


! index sets
!---------- ---------- ---------- ---------- ----------
CALL ISCreateStride ( PETSC_COMM_WORLD, NEQ, 0  , 1, is_1, IERR )
CALL ISCreateStride ( PETSC_COMM_WORLD, NEQ, NEQ, 1, is_2, IERR )


! time integration
!---------- ---------- ---------- ---------- ----------


! Newmark parameters, initialize
!---------- ---------- ---------- ---------- ----------  
a0 = 1.0d0 / (ALPHA * DT**2)
a1 = DELTA / (ALPHA * DT)
a2 = 1.0d0 / (ALPHA * DT)
a3 = 1.0d0 / (2.0d0 * ALPHA) - 1.0d0
a4 = DELTA / ALPHA - 1.0d0
a5 = (DELTA / (2.0d0 * ALPHA) - 1.0d0) * DT
!---------- ---------- ---------- ---------- ---------- 
T  = 0.0d0


! form the effective "stiffness matrix"
!---------- ---------- ---------- ---------- ----------
CALL MatAXPY ( AK_b_PETSC, A0, AM_b_PETSC, DIFFERENT_NONZERO_PATTERN, IERR )
CALL MatAXPY ( AK_b_PETSC, A1, AC_b_PETSC, DIFFERENT_NONZERO_PATTERN, IERR )


! SuperLU_dist: matrix factorization
!---------- ---------- ---------- ---------- ----------
CALL CPU_TIME(t0)
WRITE(*,*) ' LU factorization begins ... fasten your seat belts ...'

CALL MatGetFactor        ( AK_b_PETSC, "superlu_dist", MAT_FACTOR_LU, AF_PETSC, IERR )
CALL MatGetSize          ( AK_b_PETSC, m, PETSC_NULL_INTEGER, IERR )
CALL ISCreateStride      ( PETSC_COMM_WORLD, m, 0, 1, perm, IERR )
CALL MatLUFactorSymbolic ( AF_PETSC, AK_b_PETSC, perm, perm, mat_factor_info, IERR )
CALL MatLUFactorNumeric  ( AF_PETSC, AK_b_PETSC, mat_factor_info, IERR )
!---------- ---------- ---------- ---------- ----------

CALL MatFactorGetSolverPackage ( AF_PETSC, solver_package, IERR )
WRITE(*,*) 'Solver Package: ', solver_package
WRITE(2,*) 'Solver Package: ', solver_package


CALL CPU_TIME(t1)
IF ( RANK == 0 ) WRITE(*,'(A40,F10.2,A10)') 'Factorization time:', t1-t0, 'seconds'
IF ( RANK == 0 ) WRITE(2,'(A40,F10.2,A10)') 'Factorization time:', t1-t0, 'seconds'
!---------- ---------- ---------- ---------- ----------


! NEWMARK INTEGRATION
!---------- ---------- ---------- ---------- ----------
DO ISTEP = 1 , NSTEP


! print out
!---------- ---------- ---------- ---------- ----------
! get the second sub-vector
   CALL VecGetSubVector ( X_b_PETSC, is_1, U_PETSC, IERR )
! output
   CALL OUT_NEWMARK_PETSC ( U_PETSC, ID, LTRANS, T, ENERGY )
! re-store the second sub-vector
   CALL VecRestoreSubVector ( X_b_PETSC, is_1, U_PETSC, IERR )
!---------- ---------- ---------- ---------- ---------- 


   T = DBLE(ISTEP) * DT
   IF ( RANK == 0 ) WRITE (*,'(A20,F10.5,A20,I10)') 'TIME=', T, 'KSP Iter=', its


!  external load case ( forward & inverse problem )
!---------- ---------- ---------- ---------- ---------- 
   CALL TIME_FUNCTION ( T, FAC1 )
   CALL VecSet  ( R_b_PETSC, 0.0d0          , IERR )
   CALL VecAXPY ( R_b_PETSC, FAC1, B_b_PETSC, IERR )
!---------- ---------- ---------- ---------- ----------


! Newmark's method (rhs vector)
!---------- ---------- ---------- ---------- ----------        
   CALL VecSet     ( HELP_PETSC, 0.0d0, IERR )
   CALL VecMAXPY   ( HELP_PETSC, 3, [a0, a2, a3], [X_b_PETSC, Xd_b_PETSC, Xdd_b_PETSC], IERR )
   CALL MatMultAdd ( AM_b_PETSC, HELP_PETSC, R_b_PETSC, R_b_PETSC, IERR )

   CALL VecSet     ( HELP_PETSC, 0.0d0, IERR )
   CALL VecMAXPY   ( HELP_PETSC, 3, [a1, a4, a5], [X_b_PETSC, Xd_b_PETSC, Xdd_b_PETSC], IERR )
   CALL MatMultAdd ( AC_b_PETSC, HELP_PETSC, R_b_PETSC, R_b_PETSC, IERR )


! Triangular solves ( k_eff u = rhs )
!---------- ---------- ---------- ---------- ----------
   CALL MatSolve ( AF_PETSC, R_b_PETSC, HELP_PETSC, IERR )


! Newmark Updates
!---------- ---------- ---------- ---------- ----------
! keep the old Xdd
   CALL VecCopy  ( Xdd_b_PETSC, Xdd_OLD_b_PETSC, IERR )
! update Xdd_n+1
   CALL VecScale ( Xdd_b_PETSC, -a3, IERR )
   CALL VecMAXPY ( Xdd_b_PETSC, 3, [a0, -a0, -a2], [HELP_PETSC, X_b_PETSC, Xd_b_PETSC], IERR )
! update Xd_n+1
   CALL VecScale (  Xd_b_PETSC, -a4, IERR )
   CALL VecMAXPY (  Xd_b_PETSC, 3, [a1, -a1, -a5], [HELP_PETSC, X_b_PETSC, Xdd_OLD_b_PETSC], IERR )
! transfer X_n+1
   CALL VecCopy  (  HELP_PETSC, X_b_PETSC, IERR )


END DO
!---------- ---------- ---------- ---------- ----------


! destroy PETSc objects
!---------- ---------- ---------- ---------- ----------
CALL VecDestroy (Xdd_b_PETSC, IERR )
CALL VecDestroy ( Xd_b_PETSC, IERR )
CALL VecDestroy (  X_b_PETSC, IERR )
CALL VecDestroy (  R_b_PETSC, IERR )
CALL VecDestroy (    U_PETSC, IERR )

CALL VecDestroy (     Help_PETSC, IERR )
CALL VecDestroy (Xdd_OLD_b_PETSC, IERR )

CALL MatDestroy ( AF_PETSC, IERR )
CALL ISDestroy  ( perm, IERR )
CALL ISDestroy  ( is_1, IERR )
CALL ISDestroy  ( is_2, IERR )
!---------- ---------- ---------- ---------- ----------


RETURN
END
