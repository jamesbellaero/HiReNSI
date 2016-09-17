!************************************************
! Explicit time stepping based on Runge-Kutta method, 2 stage-2nd order
!************************************************
! REVISION : MONDAY, 21 May 2012

SUBROUTINE EXPLICIT_RK2_PETSC ( AK_PETSC, DIAG_M_PETSC, AC_PETSC, B_PETSC, AK_RD_PETSC, DIAG_M_RD_PETSC, ID, LTRANS, RMJ_PETSC, BACL, RANK )
USE PARAMETERS
IMPLICIT DOUBLE PRECISION (A-H,O-Z)

!==================================================================================================
! - PETSC INCLUDE FILES
#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscmat.h"

! - PETSC VARIABLES AND MATRICES
Mat            :: AK_PETSC, AC_PETSC, AK_RD_PETSC
Vec            :: B_PETSC
Vec            :: RMJ_PETSC, DIAG_M_PETSC, DIAG_M_RD_PETSC

PetscErrorCode :: IERR
PetscMPIInt    :: SIZE, RANK
PetscScalar    :: ENERGY

PetscScalar    :: zero, mass_shift
PetscScalar    :: fac7, fac8
PetscScalar    :: AG_X
PetscScalar    :: NORM_RMJ

! - PETSc local vectors
Vec            :: x1_PETSC, x2_PETSC
Vec            :: k1_t_PETSC, k1_b_PETSC, k2_t_PETSC, k2_b_PETSC, help_PETSC
!==================================================================================================


! in 
!---------- ---------- ---------- ---------- ----------
DIMENSION BACL(NSTEP,NDIM)
DIMENSION ID(NJ, NDOF), LTRANS(NTRANS)
!---------- ---------- ---------- ---------- ----------


! local
!---------- ---------- ---------- ---------- ----------
DIMENSION AG(NDIM)
!---------- ---------- ---------- ---------- ----------


I_PRINT_STEP =  10
zero         =  0.0d0
mass_shift   =  1.d-7


! (1) earthquake analysis        vs  (2) specified external load
!---------- ---------- ---------- ---------- ----------
IFLAG_LOAD   = 2


!==================================================================================================
! - DEFINE PETSC OBJECTS
CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, x1_PETSC,   IERR )
CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, x2_PETSC,   IERR )
CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, k1_t_PETSC, IERR )
CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, k1_b_PETSC, IERR )
CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, k2_t_PETSC, IERR )
CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, k2_b_PETSC, IERR )
CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, help_PETSC, IERR )

CALL VecSet       ( x1_PETSC, 0.0d0, ierr)
CALL VecSet       ( x2_PETSC, 0.0d0, ierr)


! inverse mass for explicit methods
!CALL VecCreateMPI ( PETSC_COMM_WORLD, PETSC_DECIDE, NEQ, DIAG_M_INVERSE_PETSC, IERR )
!==================================================================================================


! inverse of the diagonal mass matrix
!---------- ---------- ---------- ---------- ----------  
CALL VecReciprocal ( DIAG_M_PETSC, IERR )


! modify RMJ since inv(M) is premultiplied by this vector
! note that we are only considering x-component of ground acceleration
! k1_t_PETSC acts as Help vector here
!---------- ---------- ---------- ---------- ----------  
CALL VecNorm ( RMJ_PETSC, NORM_INFINITY, NORM_RMJ, IERR )
mass_shift = mass_shift * NORM_RMJ

CALL VecCopy            (  RMJ_PETSC, k1_t_PETSC           , IERR )
CALL VecShift           ( k1_t_PETSC, mass_shift           , IERR )
CALL VecPointwiseDivide (  RMJ_PETSC, RMJ_PETSC, k1_t_PETSC, IERR )


! form right hand side matrices
!---------- ---------- ---------- ---------- ----------
CALL MatDiagonalScale   ( AK_PETSC, DIAG_M_PETSC, PETSC_NULL, IERR )
CALL MatDiagonalScale   ( AC_PETSC, DIAG_M_PETSC, PETSC_NULL, IERR )
CALL VecPointwiseMult   ( B_PETSC,  B_PETSC, DIAG_M_PETSC   , IERR )

CALL MatScale ( AK_PETSC, -1.0D0, IERR )
CALL MatScale ( AC_PETSC, -1.0D0, IERR )


! Explicit Integration
!---------- ---------- ---------- ---------- ----------
DO ISTEP = 1 , NSTEP

   T = DBLE( ISTEP - 1 ) * DT
   IF ( RANK == 0 .AND. MOD ( ISTEP - 1 , I_PRINT_STEP ) == 0 ) THEN
      WRITE (*,'(A30,F10.5)') 'TIME=', T
   END IF


! earthquake analysis (1) vs specified external load (2)
!---------- ---------- ---------- ---------- ---------- ---------- ----------
   SELECT CASE (IFLAG_LOAD)

      CASE(1)
! earthquake analysis case
! only x-component of acceleration has been considered. see Assemble_PETSc.F90
!---------- ---------- ---------- ---------- ----------              
         AG_X = BACL(ISTEP, 1)
         CALL VecSet   ( k1_b_PETSC, zero           , IERR )
         CALL VecAXPY  ( k1_b_PETSC, AG_X, RMJ_PETSC, IERR )

         AG_X = BACL(ISTEP+1, 1)
         CALL VecSet   ( k2_b_PETSC, zero           , IERR )
         CALL VecAXPY  ( k2_b_PETSC, AG_X, RMJ_PETSC, IERR )
!---------- ---------- ---------- ---------- ---------- 


      CASE(2)
!  external load case ( forward & inverse problem )
!---------- ---------- ---------- ---------- ----------
         CALL TIME_FUNCTION ( T, fac7 )
         CALL VecSet   ( k1_b_PETSC, zero         , IERR )
         CALL VecAXPY  ( k1_b_PETSC, fac7, B_PETSC, IERR )

         CALL TIME_FUNCTION ( T+dt, fac8 )
         CALL VecSet   ( k2_b_PETSC, zero         , IERR )
         CALL VecAXPY  ( k2_b_PETSC, fac8, B_PETSC, IERR )
!---------- ---------- ---------- ---------- ----------


   END SELECT
!---------- ---------- ---------- ---------- ---------- ---------- ----------


! RK-2 Explicit time-stepping:
!---------- ---------- ---------- ---------- ----------        
! 1. compute k1_t
   CALL VecCopy ( x2_PETSC, k1_t_PETSC, IERR )
! 2. compute k1_b
   CALL MatMultAdd ( AC_PETSC, x2_PETSC, k1_b_PETSC, k1_b_PETSC, IERR )
   CALL MatMultAdd ( AK_PETSC, x1_PETSC, k1_b_PETSC, k1_b_PETSC, IERR )
! 3. compute k2_t
   CALL VecWAXPY   ( k2_t_PETSC, dt, k1_b_PETSC, x2_PETSC, IERR )
! 4. compute k2_b
   CALL VecWAXPY   ( help_PETSC, dt, k1_t_PETSC, x1_PETSC, IERR )
   CALL MatMultAdd ( AC_PETSC, k2_t_PETSC, k2_b_PETSC, k2_b_PETSC, IERR )
   CALL MatMultAdd ( AK_PETSC, help_PETSC, k2_b_PETSC, k2_b_PETSC, IERR )
! 5. compute disp_n+1
   CALL VecAXPY ( x1_PETSC, 0.5d0*dt, k1_t_PETSC, IERR )
   CALL VecAXPY ( x1_PETSC, 0.5d0*dt, k2_t_PETSC, IERR )
! 6. compute velocity_n+1
   CALL VecAXPY ( x2_PETSC, 0.5d0*dt, k1_b_PETSC, IERR )
   CALL VecAXPY ( x2_PETSC, 0.5d0*dt, k2_b_PETSC, IERR )

! print out
!---------- ---------- ---------- ---------- ---------- 
   IF ( MOD ( ISTEP - 1 , I_PRINT_STEP ) == 0 ) THEN
      CALL OUT_NEWMARK_PETSC ( x1_PETSC, ID, LTRANS, T, ENERGY )
   END IF


END DO
!---------- ---------- ---------- ---------- ----------


CALL VecDestroy ( x1_PETSC  , IERR )
CALL VecDestroy ( x2_PETSC  , IERR )
CALL VecDestroy ( k1_t_PETSC  , IERR )
CALL VecDestroy ( k1_b_PETSC  , IERR )
CALL VecDestroy ( k2_t_PETSC  , IERR )
CALL VecDestroy ( k2_b_PETSC  , IERR )
CALL VecDestroy ( help_PETSC  , IERR )


RETURN
END


!************************************************
! Explicit time stepping based on Runge-Kutta method, 4 stage-4th order
!************************************************
! REVISION : W 26 June 2013

SUBROUTINE EXPLICIT_RK4_PETSC ( AK_PETSC, DIAG_M_PETSC, AC_PETSC, B_PETSC, AK_RD_PETSC, DIAG_M_RD_PETSC, ID, LTRANS, RMJ_PETSC, BACL, RANK, NNDH, NNVH, NNAH, EqDis, EqVel, EqAcc )

USE PARAMETERS
Use Results
USE IFPORT
IMPLICIT DOUBLE PRECISION (A-H,O-Z)

Integer, Intent(In), Dimension ( NNDH * NDim )  :: EqDis
Integer, Intent(In), Dimension ( NNVH * NDim )  :: EqVel
Integer, Intent(In), Dimension ( NNAH * NDim )  :: EqAcc

!==================================================================================================
! - PETSC INCLUDE FILES
#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"
#include "finclude/petscis.h"

! - PETSC VARIABLES AND MATRICES
Mat            :: AK_PETSC, AC_PETSC, AK_RD_PETSC
Vec            :: B_PETSC
Vec            :: RMJ_PETSC, DIAG_M_PETSC, DIAG_M_RD_PETSC

PetscErrorCode :: IERR, ErrPTC
PetscMPIInt    :: SIZE, RANK
PetscScalar    :: ENERGY

PetscScalar    :: zero, mass_shift
PetscScalar    :: fac7, fac8, fac9
PetscScalar    :: AG_X
PetscScalar    :: NORM_RMJ

! - PETSc local vectors
Vec            :: x1_PETSC, x2_PETSC
Vec            :: k1_t_PETSC, k1_b_PETSC, k2_t_PETSC, k2_b_PETSC, help_PETSC
Vec            :: k3_t_PETSC, k3_b_PETSC, k4_t_PETSC, k4_b_PETSC
Vec            :: HELP_1_PETSC

PetscScalar, pointer    :: U(:), UD(:), UDD(:)
!==================================================================================================


! in 
!---------- ---------- ---------- ---------- ----------
DIMENSION BACL(NSTEP,NDIM)
DIMENSION ID(NJ, NDOF), LTRANS(NTRANS)
!---------- ---------- ---------- ---------- ----------


! local
!---------- ---------- ---------- ---------- ----------
DIMENSION AG(NDIM)
!---------- ---------- ---------- ---------- ----------


IF (RANK == 0) WRITE(*,*)'2D RK4 Explicit time integration'


I_PRINT_STEP =  1
zero         =  0.0d0
mass_shift   =  1.d-7


! (1) earthquake analysis        vs  (2) specified external load
!     earthquake case INCOMPLETE
!---------- ---------- ---------- ---------- ----------
IFLAG_LOAD   = 2


!==================================================================================================
! - DEFINE PETSC OBJECTS
CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, x1_PETSC  , IERR )
CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, x2_PETSC  , IERR )
CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, k1_t_PETSC, IERR )
CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, k1_b_PETSC, IERR )
CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, k2_t_PETSC, IERR )
CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, k2_b_PETSC, IERR )
CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, help_PETSC, IERR )

CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, k3_t_PETSC, IERR )
CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, k3_b_PETSC, IERR )
CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, k4_t_PETSC, IERR )
CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, k4_b_PETSC, IERR )

CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, HELP_1_PETSC, IERR )

CALL VecSet       ( x1_PETSC, 0.0d0, ierr)
CALL VecSet       ( x2_PETSC, 0.0d0, ierr)


! inverse mass for explicit methods
!CALL VecCreateMPI ( PETSC_COMM_WORLD, PETSC_DECIDE, NEQ, DIAG_M_INVERSE_PETSC, IERR )
!==================================================================================================




! inverse of the diagonal mass matrix
!---------- ---------- ---------- ---------- ----------  
CALL VecReciprocal ( DIAG_M_PETSC, IERR )


! modify RMJ since inv(M) is premultiplied by this vector
! note that we are only considering x-component of ground acceleration
! k1_t_PETSC acts as Help vector here
!---------- ---------- ---------- ---------- ----------  
CALL VecNorm ( RMJ_PETSC, NORM_INFINITY, NORM_RMJ, IERR )
mass_shift = mass_shift * NORM_RMJ


CALL VecCopy            (  RMJ_PETSC, k1_t_PETSC           , IERR )
CALL VecShift           ( k1_t_PETSC, mass_shift           , IERR )
CALL VecPointwiseDivide (  RMJ_PETSC, RMJ_PETSC, k1_t_PETSC, IERR )


IF (RANK == 0) WRITE(*,*)'Rhs is complete'


! form right hand side matrices
!---------- ---------- ---------- ---------- ----------
CALL MatDiagonalScale   ( AK_PETSC, DIAG_M_PETSC, PETSC_NULL_OBJECT, IERR )
CALL MatDiagonalScale   ( AC_PETSC, DIAG_M_PETSC, PETSC_NULL_OBJECT, IERR )
CALL VecPointwiseMult   ( B_PETSC,  B_PETSC, DIAG_M_PETSC   , IERR )

CALL MatScale ( AK_PETSC, -1.0D0, IERR )
CALL MatScale ( AC_PETSC, -1.0D0, IERR )


IF (RANK == 0) WRITE(*,*)'Lhs is complete'


! print out & Visualization
!---------- ---------- ---------- ---------- ---------- ---------------------------------------------------------------------------------------------
Call VecGetOwnershipRange ( x1_PETSC, IStart0, IEnd0, ErrPTC ) ;
IStart = IStart0
IEnd   = IEnd0
!---------- ---------- ---------- ---------- ---------- ---------------------------------------------------------------------------------------------


! Explicit Integration
!---------- ---------- ---------- ---------- ----------
DO ISTEP = 1 , NSTEP

   T = DBLE( ISTEP - 1 ) * DT
   IF ( RANK == 0 .AND. MOD ( ISTEP - 1 , I_PRINT_STEP ) == 0 ) THEN
      WRITE (*,'(A30,F10.5)') 'TIME=', T
   END IF


! earthquake analysis (1) vs specified external load (2)
!---------- ---------- ---------- ---------- ---------- ---------- ----------
   SELECT CASE (IFLAG_LOAD)


      CASE(1)
! earthquake analysis case
! only x-component of acceleration has been considered. see Assemble_PETSc.F90
!---------- ---------- ---------- ---------- ----------              
         AG_X = BACL(ISTEP, 1)
         CALL VecSet   ( k1_b_PETSC, 0.0d0          , IERR )
         CALL VecAXPY  ( k1_b_PETSC, AG_X, RMJ_PETSC, IERR )

         AG_X = BACL(ISTEP+1, 1)
         CALL VecSet   ( k2_b_PETSC, 0.0d0          , IERR )
         CALL VecAXPY  ( k2_b_PETSC, AG_X, RMJ_PETSC, IERR )
!---------- ---------- ---------- ---------- ---------- 


      CASE(2)
!  external load case ( forward & inverse problem )
!---------- ---------- ---------- ---------- ----------
         CALL TIME_FUNCTION ( T, fac7 )
         CALL VecSet   ( k1_b_PETSC, 0.0d0        , IERR )
         CALL VecAXPY  ( k1_b_PETSC, fac7, B_PETSC, IERR )

         CALL TIME_FUNCTION ( T+dt, fac8 )
         CALL VecSet   ( k4_b_PETSC, 0.0d0        , IERR )
         CALL VecAXPY  ( k4_b_PETSC, fac8, B_PETSC, IERR )

         CALL TIME_FUNCTION ( T + 0.50d0 * dt, fac9 )
         CALL VecSet   ( k2_b_PETSC, 0.0d0        , IERR )
         CALL VecAXPY  ( k2_b_PETSC, fac9, B_PETSC, IERR )
         CALL VecSet   ( k3_b_PETSC, 0.0d0        , IERR )
         CALL VecAXPY  ( k3_b_PETSC, fac9, B_PETSC, IERR )
!---------- ---------- ---------- ---------- ----------


   END SELECT
!---------- ---------- ---------- ---------- ---------- ---------- ----------


! RK-4 Explicit time-stepping:
!---------- ---------- ---------- ---------- ----------        
! 1. compute k1_t
   CALL VecCopy ( x2_PETSC, k1_t_PETSC, IERR )
! 2. compute k1_b
   CALL MatMultAdd ( AC_PETSC, x2_PETSC, k1_b_PETSC, k1_b_PETSC, IERR )
   CALL MatMultAdd ( AK_PETSC, x1_PETSC, k1_b_PETSC, k1_b_PETSC, IERR )
! 3. compute k2_t
   CALL VecWAXPY   ( k2_t_PETSC, 0.5d0*dt, k1_b_PETSC, x2_PETSC, IERR )
! 4. compute k2_b
   CALL VecWAXPY   ( help_PETSC, 0.5d0*dt, k1_t_PETSC, x1_PETSC, IERR )
   CALL MatMultAdd ( AC_PETSC, k2_t_PETSC, k2_b_PETSC, k2_b_PETSC, IERR )
   CALL MatMultAdd ( AK_PETSC, help_PETSC, k2_b_PETSC, k2_b_PETSC, IERR )

! 5. compute k3_t
   CALL VecWAXPY   ( k3_t_PETSC, 0.5d0*dt, k2_b_PETSC, x2_PETSC, IERR )
! 6. compute k3_b
   CALL VecWAXPY   ( help_PETSC, 0.5d0*dt, k2_t_PETSC, x1_PETSC, IERR )
   CALL MatMultAdd ( AC_PETSC, k3_t_PETSC, k3_b_PETSC, k3_b_PETSC, IERR )
   CALL MatMultAdd ( AK_PETSC, help_PETSC, k3_b_PETSC, k3_b_PETSC, IERR )

! 7. compute k4_t
   CALL VecWAXPY   ( k4_t_PETSC, dt, k3_b_PETSC, x2_PETSC, IERR )
! 8. compute k4_b
   CALL VecWAXPY   ( help_PETSC, dt, k3_t_PETSC, x1_PETSC, IERR )
   CALL MatMultAdd ( AC_PETSC, k4_t_PETSC, k4_b_PETSC, k4_b_PETSC, IERR )
   CALL MatMultAdd ( AK_PETSC, help_PETSC, k4_b_PETSC, k4_b_PETSC, IERR )

! 9. compute disp_n+1
   CALL VecSet  ( help_PETSC, 0.0d0, ierr)
   CALL VecMAXPY( help_PETSC, 4, [1.0d0, 2.0d0, 2.0d0, 1.0d0], [k1_t_PETSC, k2_t_PETSC, k3_t_PETSC, k4_t_PETSC], IERR )
   CALL VecAXPY ( x1_PETSC, dt/6.0d0, help_PETSC, IERR )
!10. compute velocity_n+1
   CALL VecSet  ( help_PETSC, 0.0d0, ierr)
   CALL VecMAXPY( help_PETSC, 4, [1.0d0, 2.0d0, 2.0d0, 1.0d0], [k1_b_PETSC, k2_b_PETSC, k3_b_PETSC, k4_b_PETSC], IERR )
   CALL VecAXPY ( x2_PETSC, dt/6.0d0, help_PETSC, IERR )



! print out
!---------- ---------- ---------- ---------- ---------- 
   IF ( MOD ( ISTEP - 1 , I_PRINT_STEP ) == 0 ) THEN

      CALL TOTAL_ENERGY_RD_PETSC_EXPLICIT ( AK_RD_PETSC, DIAG_M_RD_PETSC, x1_PETSC, x2_PETSC, HELP_1_PETSC, ENERGY )
!      CALL OUT_NEWMARK_PETSC ( x1_PETSC, ID, LTRANS, T, ENERGY )


! print out & Visualization
!---------- ---------- ---------- ---------- ---------- ---------------------------------------------------------------------------------------------

    If ( NNDH /= 0 ) Call VecGetArrayF90     (   x1_PETSC,   U,  ErrPTC ) ; ! NNDH
    If ( NNVH /= 0 ) Call VecGetArrayF90     (   x2_PETSC,  UD,  ErrPTC ) ; ! NNVH
    If ( NNAH /= 0 ) Call VecGetArrayF90     ( k1_b_PETSC, UDD,  ErrPTC ) ; ! NNAH

    ! Writing down the history of displacement, velocity and acceleration
    Call ResDyn   ( NDIM,     NNDH, NNVH, NNAH, IStep,     IStart, IEnd,     T, Energy,     EqDis, EqVel, EqAcc,    U, UD, UDD ) ;
    !Call ParaView ( NDim,     IStep,            NJ, NEl,     ELT,     INod, ID,      U,     XYZ,     Name, OutDir   ) ;
!    If ( Any( U /= 0.0_dbl ) .and. Mod(IStep,500)==0 ) Call ParaView_9noded ( NDim,     IStep,            NJ, NEl,     ELT,     INod, ID,      U, UD, UDD,    XYZ,     Name, OutDir   ) ;

!### this must be changed to paraview_9node

    If ( NNDH /= 0 ) Call VecRestoreArrayF90 (   x1_PETSC,   U,  ErrPTC ) ; ! NNDH
    If ( NNVH /= 0 ) Call VecRestoreArrayF90 (   x2_PETSC,  UD,  ErrPTC ) ; ! NNVH
    If ( NNAH /= 0 ) Call VecRestoreArrayF90 ( k1_b_PETSC, UDD,  ErrPTC ) ; ! NNAH


!   CALL OUT_NEWMARK_PETSC ( x1_PETSC, ID, LTRANS, T, ENERGY )

   IF ( MOD ( ISTEP - 1 , 1 * I_PRINT_STEP ) == 0 ) THEN
!      CALL ParaView ( XYZ, INOD, ID, U_PETSC, ISTEP )
   END IF
!---------- ---------- ---------- ---------- ---------- ---------------------------------------------------------------------------------------------


   END IF


END DO
!---------- ---------- ---------- ---------- ----------


CALL VecDestroy (   x1_PETSC, IERR )
CALL VecDestroy (   x2_PETSC, IERR )
CALL VecDestroy ( k1_t_PETSC, IERR )
CALL VecDestroy ( k1_b_PETSC, IERR )
CALL VecDestroy ( k2_t_PETSC, IERR )
CALL VecDestroy ( k2_b_PETSC, IERR )
CALL VecDestroy ( help_PETSC, IERR )

CALL VecDestroy ( k3_t_PETSC, IERR )
CALL VecDestroy ( k3_b_PETSC, IERR )
CALL VecDestroy ( k4_t_PETSC, IERR )
CALL VecDestroy ( k4_b_PETSC, IERR )

CALL VecDestroy ( HELP_1_PETSC, IERR )

RETURN
END                                      


!************************************************
! Explicit time stepping based on Runge-Kutta method, 2 stage-2nd order, PML-3D
! Comment:
! for test09.txt, dt = 0.008 and it becomes unstable at t = 25 sec. RK4 works fine with dt = 0.04 until t = 60 sec where instability infuses the solution.
!************************************************
! REVISION : F, 22 Jun 2012

SUBROUTINE EXPLICIT_RK2_PETSC_PML_3D ( AK_PETSC, DIAG_M_PETSC, AC_PETSC, AG_PETSC, B_PETSC, AK_RD_PETSC, DIAG_M_RD_PETSC, ID, LTRANS, RANK )
USE PARAMETERS
IMPLICIT DOUBLE PRECISION (A-H,O-Z)

!==================================================================================================
! - PETSC INCLUDE FILES
#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscmat.h"

! - PETSC VARIABLES AND MATRICES
Mat            :: AK_PETSC, AC_PETSC, AK_RD_PETSC
Mat            :: AG_PETSC
Vec            :: B_PETSC
Vec            :: DIAG_M_PETSC, DIAG_M_RD_PETSC

PetscErrorCode :: IERR
PetscMPIInt    :: SIZE, RANK
PetscScalar    :: ENERGY

PetscScalar    :: zero, mass_shift
PetscScalar    :: fac7, fac8

! - PETSc local vectors
Vec            :: x1_PETSC, x2_PETSC, x3_PETSC
Vec            :: k1_t_PETSC, k1_m_PETSC, k1_b_PETSC, k2_t_PETSC, k2_m_PETSC, k2_b_PETSC, help_PETSC
!==================================================================================================


! in 
!---------- ---------- ---------- ---------- ----------
DIMENSION ID(NJ, NDOF), LTRANS(NTRANS)
!---------- ---------- ---------- ---------- ----------


I_PRINT_STEP =  1
zero         =  0.0d0


!==================================================================================================
! - DEFINE PETSC OBJECTS
CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, x1_PETSC, IERR )
CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, x2_PETSC, IERR )
CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, x3_PETSC, IERR )

CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, k1_t_PETSC, IERR )
CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, k1_m_PETSC, IERR )
CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, k1_b_PETSC, IERR )

CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, k2_t_PETSC, IERR )
CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, k2_m_PETSC, IERR )
CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, k2_b_PETSC, IERR )

CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, help_PETSC, IERR )

CALL VecSet       ( x1_PETSC, 0.0d0, ierr)
CALL VecSet       ( x2_PETSC, 0.0d0, ierr)
CALL VecSet       ( x3_PETSC, 0.0d0, ierr)
!==================================================================================================


! inverse of the diagonal mass matrix
!---------- ---------- ---------- ---------- ----------  
CALL VecReciprocal ( DIAG_M_PETSC, IERR )


! form right hand side matrices
!---------- ---------- ---------- ---------- ----------
CALL MatDiagonalScale ( AK_PETSC, DIAG_M_PETSC, PETSC_NULL_OBJECT, IERR )
CALL MatDiagonalScale ( AC_PETSC, DIAG_M_PETSC, PETSC_NULL_OBJECT, IERR )
CALL MatDiagonalScale ( AG_PETSC, DIAG_M_PETSC, PETSC_NULL_OBJECT, IERR )
CALL VecPointwiseMult (  B_PETSC, B_PETSC, DIAG_M_PETSC   , IERR )

CALL MatScale ( AK_PETSC, -1.0D0, IERR )
CALL MatScale ( AC_PETSC, -1.0D0, IERR )
CALL MatScale ( AG_PETSC, -1.0D0, IERR )


! Explicit Integration
!---------- ---------- ---------- ---------- ----------
DO ISTEP = 1 , NSTEP

   T = DBLE( ISTEP - 1 ) * DT
   IF ( RANK == 0 .AND. MOD ( ISTEP - 1 , I_PRINT_STEP ) == 0 ) THEN
      WRITE (*,'(A30,F10.5)') 'TIME=', T
   END IF

   CALL TIME_FUNCTION ( T, fac7 )
   CALL VecSet   ( k1_b_PETSC, zero         , IERR )
   CALL VecAXPY  ( k1_b_PETSC, fac7, B_PETSC, IERR )

   CALL TIME_FUNCTION ( T+dt, fac8 )
   CALL VecSet   ( k2_b_PETSC, zero         , IERR )
   CALL VecAXPY  ( k2_b_PETSC, fac8, B_PETSC, IERR )


! RK-2 Explicit time-stepping:
!---------- ---------- ---------- ---------- ----------        
! 1. compute k1_t
   CALL VecCopy ( x2_PETSC, k1_t_PETSC, IERR )
! 2. compute k1_m
   CALL VecCopy ( x3_PETSC, k1_m_PETSC, IERR )
! 3. compute k1_b
   CALL MatMultAdd ( AC_PETSC, x3_PETSC, k1_b_PETSC, k1_b_PETSC, IERR )
   CALL MatMultAdd ( AK_PETSC, x2_PETSC, k1_b_PETSC, k1_b_PETSC, IERR )
   CALL MatMultAdd ( AG_PETSC, x1_PETSC, k1_b_PETSC, k1_b_PETSC, IERR )
! 4. compute k2_t
   CALL VecWAXPY   ( k2_t_PETSC, dt, k1_m_PETSC, x2_PETSC, IERR )
! 5. compute k2_m
   CALL VecWAXPY   ( k2_m_PETSC, dt, k1_b_PETSC, x3_PETSC, IERR )
! 6. compute k2_b
   CALL VecWAXPY   ( help_PETSC, dt, k1_t_PETSC, x1_PETSC, IERR )
   CALL MatMultAdd ( AC_PETSC, k2_m_PETSC, k2_b_PETSC, k2_b_PETSC, IERR )
   CALL MatMultAdd ( AK_PETSC, k2_t_PETSC, k2_b_PETSC, k2_b_PETSC, IERR )
   CALL MatMultAdd ( AG_PETSC, help_PETSC, k2_b_PETSC, k2_b_PETSC, IERR )
! 7. compute disp_n+1
   CALL VecAXPY ( x1_PETSC, 0.5d0*dt, k1_t_PETSC, IERR )
   CALL VecAXPY ( x1_PETSC, 0.5d0*dt, k2_t_PETSC, IERR )
! 8. compute velocity_n+1
   CALL VecAXPY ( x2_PETSC, 0.5d0*dt, k1_m_PETSC, IERR )
   CALL VecAXPY ( x2_PETSC, 0.5d0*dt, k2_m_PETSC, IERR )
! 9. compute acceleration_n+1
   CALL VecAXPY ( x3_PETSC, 0.5d0*dt, k1_b_PETSC, IERR )
   CALL VecAXPY ( x3_PETSC, 0.5d0*dt, k2_b_PETSC, IERR )

! print out
!---------- ---------- ---------- ---------- ---------- 
   IF ( MOD ( ISTEP - 1 , I_PRINT_STEP ) == 0 ) THEN
      CALL OUT_NEWMARK_PETSC ( x1_PETSC, ID, LTRANS, T, ENERGY )
   END IF


END DO
!---------- ---------- ---------- ---------- ----------


CALL VecDestroy ( x1_PETSC  , IERR )
CALL VecDestroy ( x2_PETSC  , IERR )
CALL VecDestroy ( x3_PETSC  , IERR )

CALL VecDestroy ( k1_t_PETSC, IERR )
CALL VecDestroy ( k1_m_PETSC, IERR )
CALL VecDestroy ( k1_b_PETSC, IERR )

CALL VecDestroy ( k2_t_PETSC, IERR )
CALL VecDestroy ( k2_m_PETSC, IERR )
CALL VecDestroy ( k2_b_PETSC, IERR )

CALL VecDestroy ( help_PETSC, IERR )


RETURN
END


!************************************************
! Explicit time stepping based on Runge-Kutta method, 4 stage-4th order, PML-3D
!************************************************
! REVISION : M 11 November 2013: be a man!
! REVISION : Th 9 January 2014: a man of principle!
! You're only here for a short visit. Don't hurry, don't worry. And be sure to smell the flowers along the way.

SUBROUTINE EXPLICIT_RK4_PETSC_PML_3D ( AK_PETSC, DIAG_iM_PETSC, AC_PETSC, AG_PETSC, B_PETSC, AK_RD_PETSC, DIAG_M_RD_PETSC, ID, LTRANS, RANK, &
                                       NNDH, NNVH, NNAH, EqDis, EqVel, EqAcc, idx_from, ID_Para, Param, &
                                       U_Store_Mapping, NStore_Mapping, U_Store_Numbers_Global ) !, k1_m_Store_Mapping, k2_m_Store_Mapping, k3_m_Store_Mapping )

USE PARAMETERS
Use Results
Use Visualizer
USE IFPORT
IMPLICIT DOUBLE PRECISION (A-H,O-Z)

! - Type DECLERATIONS -------------------------------------------------------------------------------------------------------------------------------
Type ( BasicParam )    :: Param ;                ! Holds basic parameters of each load case

Integer, Intent(In), Dimension ( NNDH * NDim )  :: EqDis
Integer, Intent(In), Dimension ( NNVH * NDim )  :: EqVel
Integer, Intent(In), Dimension ( NNAH * NDim )  :: EqAcc

!Integer, Intent(In), Dimension ( Param%IntL ( 5, 2) )        :: STEP
Integer, Intent(In), Dimension ( Param%IntM ( 5, 1) )        :: idx_from 
Integer, Intent(In), Dimension ( Param%IntM ( 5, 2), NDim )  :: ID_Para 
Integer, Allocatable, Dimension (:)                          :: idx_to 

Integer, Intent(In), Dimension ( NStore_Mapping )  ::  U_Store_Numbers_Global

Real(8), Intent(Out), Dimension ( NStore_Mapping , 0:NStep ) ::    U_Store_Mapping
!Real(8), Intent(Out), Dimension ( NStore_Mapping , 0:NStep ) :: k1_m_Store_Mapping
!Real(8), Intent(Out), Dimension ( NStore_Mapping , 0:NStep ) :: k2_m_Store_Mapping
!Real(8), Intent(Out), Dimension ( NStore_Mapping , 0:NStep ) :: k3_m_Store_Mapping
!Real(8), Intent(Out), Dimension ( NStore_Mapping , 0:NStep ) :: k4_m_Store_Mapping

Character (Kind = 1, Len = 20)  :: IndexRankTemp;! A variable for storing the Rank number in Character format for adding at the end of the input file Name in the do loop 
Character (Kind = 1, Len = 22)  :: IndexStep ;   ! A variable for storing the Step number in Character format for adding at the end of the input file Name in the do loop 
!==================================================================================================
! - PETSC INCLUDE FILES
#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"
#include "finclude/petscis.h"

! - PETSC VARIABLES AND MATRICES
Mat            :: AK_PETSC, AC_PETSC, AK_RD_PETSC
Mat            :: AG_PETSC
Vec            :: B_PETSC
Vec            :: DIAG_iM_PETSC, DIAG_M_RD_PETSC

PetscErrorCode :: IERR, ErrPTC
PetscMPIInt    :: SIZE, RANK
PetscScalar    :: ENERGY

PetscScalar    :: zero, mass_shift
PetscScalar    :: fac7, fac8, fac9

! - PETSc local vectors
Vec            :: x1_PETSC, x2_PETSC, x3_PETSC
Vec            :: k1_t_PETSC, k1_m_PETSC, k1_b_PETSC, k2_t_PETSC, k2_m_PETSC, k2_b_PETSC, help_PETSC
Vec            :: k3_t_PETSC, k3_m_PETSC, k3_b_PETSC, k4_t_PETSC, k4_m_PETSC, k4_b_PETSC
Vec            :: HELP_1_PETSC

PetscScalar, pointer    :: U(:),      UD(:),      UDD(:)
PetscScalar, pointer    :: U_Rank(:), UD_Rank(:), UDD_Rank(:) ;

Vec            ::   U_Rank_PETSC ;                 ! 
Vec            ::  UD_Rank_PETSC ;                 ! 
Vec            :: UDD_Rank_PETSC ;                 !

! - Index setting for visualization -----------------------------------------------------------------------------------------------------------------
VecScatter     :: vscat_U, vscat_UD, vscat_UDD ; ! Used to scatter ghost node solution, this is just for Paraview
IS             :: IS_from, IS_to  ;              ! Index sets for retrieving the ghost node solution

!==================================================================================================


! in
!---------- ---------- ---------- ---------- ----------
DIMENSION ID ( NJ, NDOF ), LTRANS ( NTRANS )
!---------- ---------- ---------- ---------- ----------


! local
!---------- ---------- ---------- ---------- ----------
Dimension u_small ( NStore_Mapping )
!---------- ---------- ---------- ---------- ----------


IF (RANK == 0 .AND. I_PRINT_SCREEN == 1) WRITE(*,*) '3D RK4 State'


I_PRINT_STEP = 1000
zero         = 0.0d0


! define Index Set for Paraview ---------------------------------------------------------------------------------------------------------------------
Allocate ( idx_to ( Param%IntM ( 5, 1) ) ) ;                                        ! ( 5, 1) corresponds to NE_Para (NEq 20-noded FEM for visualization)

Call VecCreateSeq ( PETSC_COMM_SELF, Param%IntM ( 5, 1),   U_Rank_PETSC, ErrPTC );  ! Temporary vector to hold displacements in nodes for Paraview
Call VecCreateSeq ( PETSC_COMM_SELF, Param%IntM ( 5, 1),  UD_Rank_PETSC, ErrPTC );  ! Temporary vector to hold velocity in nodes for Paraview
Call VecCreateSeq ( PETSC_COMM_SELF, Param%IntM ( 5, 1), UDD_Rank_PETSC, ErrPTC );  ! Temporary vector to hold acceleration in nodes for Paraview

ForAll ( I = 1:Param%IntM ( 5, 1) ) idx_to (I) = I-1 ;

Call ISCreateGeneral ( PETSC_COMM_SELF, int4(Param%IntM ( 5, 1)), idx_from, PETSC_COPY_VALUES, IS_from, ErrPTC );
Call ISCreateGeneral ( PETSC_COMM_SELF, int4(Param%IntM ( 5, 1)), idx_to,   PETSC_COPY_VALUES, IS_to,   ErrPTC );

DeAllocate ( idx_to ) ;

NSParaview = Param%IntL ( 5, 3) ;  ! Step numbers for Paraview
! ---------------------------------------------------------------------------------------------------------------------------------------------------
!==================================================================================================
! - DEFINE PETSC OBJECTS
CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, x1_PETSC, IERR )
CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, x2_PETSC, IERR )
CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, x3_PETSC, IERR )

CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, k1_t_PETSC, IERR )
CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, k1_m_PETSC, IERR )
CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, k1_b_PETSC, IERR )

CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, k2_t_PETSC, IERR )
CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, k2_m_PETSC, IERR )
CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, k2_b_PETSC, IERR )

CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, k3_t_PETSC, IERR )
CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, k3_m_PETSC, IERR )
CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, k3_b_PETSC, IERR )

CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, k4_t_PETSC, IERR )
CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, k4_m_PETSC, IERR )
CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, k4_b_PETSC, IERR )

CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, help_PETSC, IERR )
CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, HELP_1_PETSC, IERR )

CALL VecSet       ( x1_PETSC, 0.0d0, ierr)
CALL VecSet       ( x2_PETSC, 0.0d0, ierr)
CALL VecSet       ( x3_PETSC, 0.0d0, ierr)
! just for hygenic concerns:
CALL VecSet       ( k1_m_PETSC, 0.0d0, ierr)
CALL VecSet       ( k2_m_PETSC, 0.0d0, ierr)
CALL VecSet       ( k3_m_PETSC, 0.0d0, ierr)
CALL VecSet       ( k4_m_PETSC, 0.0d0, ierr)
!==================================================================================================

Call VecScatterCreate ( x2_PETSC,   IS_from, U_Rank_PETSC,   IS_to, vscat_U, ErrPTC );
!Call VecScatterCreate ( x3_PETSC,   IS_from, UD_Rank_PETSC,  IS_to, vscat_UD, ErrPTC );
!Call VecScatterCreate ( k1_b_PETSC, IS_from, UDD_Rank_PETSC, IS_to, vscat_UDD, ErrPTC );

! inverse of the diagonal mass matrix
!---------- ---------- ---------- ---------- ----------  
! Not needed. We did that in Solve_Forward.F90
!CALL VecReciprocal ( DIAG_M_PETSC, IERR )


! form right-hand-side matrices (left-scaling)
!---------- ---------- ---------- ---------- ----------
CALL MatDiagonalScale ( AK_PETSC, DIAG_iM_PETSC, PETSC_NULL_OBJECT, IERR )
CALL MatDiagonalScale ( AC_PETSC, DIAG_iM_PETSC, PETSC_NULL_OBJECT, IERR )
CALL MatDiagonalScale ( AG_PETSC, DIAG_iM_PETSC, PETSC_NULL_OBJECT, IERR )
CALL VecPointwiseMult ( B_PETSC,  B_PETSC, DIAG_iM_PETSC, IERR )

!IF (RANK == 0) WRITE(*,*)'Rhs is complete'

CALL MatScale ( AK_PETSC, -1.0D0, IERR )
CALL MatScale ( AC_PETSC, -1.0D0, IERR )
CALL MatScale ( AG_PETSC, -1.0D0, IERR )

!IF (RANK == 0) WRITE(*,*)'Lhs is complete'




! print out & Visualization
!---------- ---------- ---------- ---------- ---------- ---------------------------------------------------------------------------------------------
Call VecGetOwnershipRange ( x1_PETSC, IStart0, IEnd0, ErrPTC ) ;
IStart = IStart0
IEnd   = IEnd0
!---------- ---------- ---------- ---------- ---------- ---------------------------------------------------------------------------------------------


CALL MPI_Comm_size   ( PETSC_COMM_WORLD, SIZE, IERR )
CALL MPI_Comm_rank   ( PETSC_COMM_WORLD, RANK, IERR )

! - Write total wraper file  ------------------------------------------------------------------------------------------------------------------------
  If ( Rank == 0 ) Then ;
    UnFile = UN_OutallWrp ;
    Write (Unit = UnFile, FMT = "(A63)", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) '<Xdmf xmlns:xi="http://www.w3.org/2003/XInclude" Version="2.2">' ;
    Write (Unit = UnFile, FMT = "(A10)", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) '  <Domain>' ;
    Write (Unit = UnFile, FMT = "(A57)", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) "    <Grid GridType='Collection' CollectionType='Spatial'>" ;
  End If ;

! Explicit Integration
!---------- ---------- ---------- ---------- ----------
DO ISTEP = 0 , NSTEP

   T = DBLE( ISTEP ) * DT
   IF ( RANK == 0 .AND. MOD ( ISTEP+1 , I_PRINT_STEP ) == 0 ) THEN
      WRITE (*,'(A30,F10.5)') 'TIME=', T
   END IF


! storage of displacement for the inverse problem - (see test03.F90 for a simple example)
!---------- ---------- ---------- ---------- ---------- ---------------------------------------------------------------------------------------------
! 1- selection: during time integration, access u(big), pick u(small) entries.
   Call VecGetArrayF90 ( x2_PETSC, U, IERR )
   Do I = 1, NStore_Mapping
     u_small ( I ) = U ( U_Store_Numbers_Global ( I ) - IStart + 1 )
   End Do
   Call VecRestoreArrayF90 ( x2_PETSC, U, IERR )

! 2- store u(small) entries in the dense matrix: We are storing 0:N-1.
   U_Store_Mapping ( : , ISTEP ) = u_small (:)
!---------- ---------- ---------- ---------- ---------- ---------------------------------------------------------------------------------------------


! storage of k1_m, k2_m, k3_m, k4_m for the inverse problem
!---------- ---------- ---------- ---------- ---------- ---------------------------------------------------------------------------------------------
! 1.1- selection: access k1_m(big) and pick k1_m(small) entries.
!   Call VecGetArrayF90 ( k1_m_PETSC, U, IERR )
!   Do I = 1, NStore_Mapping
!     u_small ( I ) = U ( U_Store_Numbers_Global ( I ) - IStart + 1 )
!   End Do
!   Call VecRestoreArrayF90 ( k1_m_PETSC, U, IERR )
! 1.2- store k1_m(small) entries in the dense matrix: We are storing 0:N-1.
!   k1_m_Store_Mapping ( : , ISTEP ) = u_small (:)
! 2.1- selection: access k2_m(big) and pick k2_m(small) entries.
!   Call VecGetArrayF90 ( k2_m_PETSC, U, IERR )
!   Do I = 1, NStore_Mapping
!     u_small ( I ) = U ( U_Store_Numbers_Global ( I ) - IStart + 1 )
!   End Do
!   Call VecRestoreArrayF90 ( k2_m_PETSC, U, IERR )
! 2.2- store k2_m(small) entries in the dense matrix: We are storing 0:N-1.
!   k2_m_Store_Mapping ( : , ISTEP ) = u_small (:)
! 3.1- selection: access k3_m(big) and pick k3_m(small) entries.
!   Call VecGetArrayF90 ( k3_m_PETSC, U, IERR )
!   Do I = 1, NStore_Mapping
!     u_small ( I ) = U ( U_Store_Numbers_Global ( I ) - IStart + 1 )
!   End Do
!   Call VecRestoreArrayF90 ( k3_m_PETSC, U, IERR )
! 3.2- store k2_m(small) entries in the dense matrix: We are storing 0:N-1.
!   k3_m_Store_Mapping ( : , ISTEP ) = u_small (:)
! 4.1- selection: access k4_m(big) and pick k4_m(small) entries.
!   Call VecGetArrayF90 ( k4_m_PETSC, U, IERR )
!   Do I = 1, NStore_Mapping
!     u_small ( I ) = U ( U_Store_Numbers_Global ( I ) - IStart + 1 )
!   End Do
!   Call VecRestoreArrayF90 ( k4_m_PETSC, U, IERR )
! 4.2- store k2_m(small) entries in the dense matrix: We are storing 0:N-1.
!   k4_m_Store_Mapping ( : , ISTEP ) = u_small (:)
!---------- ---------- ---------- ---------- ---------- ---------------------------------------------------------------------------------------------


! print out
!---------- ---------- ---------- ---------- ---------- ---------------------------------------------------------------------------------------------
!---------- ---------- ---------- ---------- ---------- ---------------------------------------------------------------------------------------------
! if RHS is f(t), then x1: displacement history, x2: displacement, x3: velocity
   IF ( MOD ( ISTEP , I_PRINT_STEP ) == 0 ) THEN
      CALL TOTAL_ENERGY_RD_PETSC_EXPLICIT ( AK_RD_PETSC, DIAG_M_RD_PETSC, x2_PETSC, x3_PETSC, HELP_1_PETSC, ENERGY )


! print out & Visualization
!---------- ---------- ---------- ---------- ---------- ---------------------------------------------------------------------------------------------

    If ( NNDH /= 0 ) Call VecGetArrayF90     (   x2_PETSC,   U,  ErrPTC ) ; ! NNDH
    If ( NNVH /= 0 ) Call VecGetArrayF90     (   x3_PETSC,  UD,  ErrPTC ) ; ! NNVH
    If ( NNAH /= 0 ) Call VecGetArrayF90     ( k1_b_PETSC, UDD,  ErrPTC ) ; ! NNAH ! we are passing the wrong vector for acceleration. kill this later

    ! Writing down the history of displacement, velocity and acceleration
    Call ResDyn   ( NDIM,     NNDH, NNVH, NNAH, IStep,     IStart, IEnd,     T, Energy,     EqDis, EqVel, EqAcc,    U, UD, UDD ) ;

    If ( NNDH /= 0 ) Call VecRestoreArrayF90 (   x2_PETSC,   U,  ErrPTC ) ; ! NNDH
    If ( NNVH /= 0 ) Call VecRestoreArrayF90 (   x3_PETSC,  UD,  ErrPTC ) ; ! NNVH
    If ( NNAH /= 0 ) Call VecRestoreArrayF90 ( k1_b_PETSC, UDD,  ErrPTC ) ; ! NNAH


   END IF
!---------- ---------- ---------- ---------- ---------- ---------------------------------------------------------------------------------------------
    ! Paraview 
    NSParaview = 20000

    If ( Mod(IStep,NSParaview)==0 ) Then ;
        If ( Rank == 0 ) Then ;  ! 
          Write (IndexStep, *) IStep ;  ! Converts Step number to Character foramt for the file Name
            Do I = 0, Size-1 ;
              Write (IndexRankTemp, *) I ;  ! Converts Rank number to Character foramt for the file Name
              UnFile = UN_OutallWrp ;
              Write (UnFile, * ) "     <xi:include href='",TRIM(AnaNAME)//'_'//Trim(AdjustL(IndexSize))//'_'//Trim(AdjustL(IndexRankTemp))//'_'//Trim(AdjustL(IndexStep))//'.xmf',"'/>" ;
            End Do ;
        End If ;

      Call VecScatterBegin ( vscat_U,   x2_PETSC,   U_Rank_PETSC,   INSERT_VALUES, SCATTER_FORWARD, ErrPTC );
      Call VecScatterEnd   ( vscat_U,   x2_PETSC,   U_Rank_PETSC,   INSERT_VALUES, SCATTER_FORWARD, ErrPTC );

!      Call VecScatterBegin ( vscat_UD,  x3_PETSC,   UD_Rank_PETSC,  INSERT_VALUES, SCATTER_FORWARD, ErrPTC );
!      Call VecScatterEnd   ( vscat_UD,  x3_PETSC,   UD_Rank_PETSC,  INSERT_VALUES, SCATTER_FORWARD, ErrPTC );

!      Call VecScatterBegin ( vscat_UDD, k1_b_PETSC, UDD_Rank_PETSC, INSERT_VALUES, SCATTER_FORWARD, ErrPTC );
!      Call VecScatterEnd   ( vscat_UDD, k1_b_PETSC, UDD_Rank_PETSC, INSERT_VALUES, SCATTER_FORWARD, ErrPTC );

      Call VecGetArrayF90 (   U_Rank_PETSC,   U_Rank,  ErrPTC ) ; ! NNDH Displacement
      Call VecGetArrayF90 (  UD_Rank_PETSC,  UD_Rank,  ErrPTC ) ; ! NNDH Displacement
      Call VecGetArrayF90 ( UDD_Rank_PETSC, UDD_Rank,  ErrPTC ) ; ! NNDH Displacement

        Call Paraview_HDF5    (                                                                                                         &
        NDim, NNode,                                                                                                                    & ! Integer (1) Variables
        !                                                                                                                               & ! Integer (2) Variables
        !                                                                                                                               & ! Integer (4) Variables
        IStep, NEl, Param%IntM( 5, 2), Param%IntM( 5, 4),                                                                               & ! Integer (8) Variables
        T,                                                                                                                              & ! Real Variables
        ID_Para,                                                                                                                        & ! Integer Arrays
        U_Rank, UD_Rank, UDD_Rank,                                                                                                      & ! Real Arrays
        OutDir, ModelName, AnaName, IndexRank, IndexSize                                                                                & ! Characters
        !                                                                                                                               & ! Type
        ) ;

      Call VecRestoreArrayF90 (   U_Rank_PETSC,   U_Rank,  ErrPTC ) ; ! NNDH
      Call VecRestoreArrayF90 (  UD_Rank_PETSC,  UD_Rank,  ErrPTC ) ; ! NNDH
      Call VecRestoreArrayF90 ( UDD_Rank_PETSC, UDD_Rank,  ErrPTC ) ; ! NNDH

    End If ;
!---------- ---------- ---------- ---------- ---------- ---------------------------------------------------------------------------------------------
!---------- ---------- ---------- ---------- ---------- ---------------------------------------------------------------------------------------------


   CALL TIME_FUNCTION ( T, fac7 )
   CALL VecSet   ( k1_b_PETSC, 0.0d0        , IERR )
   CALL VecAXPY  ( k1_b_PETSC, fac7, B_PETSC, IERR )

   CALL TIME_FUNCTION ( T+dt, fac8 )
   CALL VecSet   ( k4_b_PETSC, 0.0d0        , IERR )
   CALL VecAXPY  ( k4_b_PETSC, fac8, B_PETSC, IERR )

   CALL TIME_FUNCTION ( T + 0.5d0 * dt, fac9 )
   CALL VecSet   ( k2_b_PETSC, 0.0d0        , IERR )
   CALL VecAXPY  ( k2_b_PETSC, fac9, B_PETSC, IERR )
   CALL VecSet   ( k3_b_PETSC, 0.0d0        , IERR )
   CALL VecAXPY  ( k3_b_PETSC, fac9, B_PETSC, IERR )


! RK-4 Explicit time-stepping:
!---------- ---------- ---------- ---------- ----------
! 1. compute k1_t
   CALL VecCopy ( x2_PETSC, k1_t_PETSC, IERR )
! 2. compute k1_m
   CALL VecCopy ( x3_PETSC, k1_m_PETSC, IERR )
! 3. compute k1_b
   CALL MatMultAdd ( AC_PETSC, x3_PETSC, k1_b_PETSC, k1_b_PETSC, IERR )
   CALL MatMultAdd ( AK_PETSC, x2_PETSC, k1_b_PETSC, k1_b_PETSC, IERR )
   CALL MatMultAdd ( AG_PETSC, x1_PETSC, k1_b_PETSC, k1_b_PETSC, IERR )
! 4. compute k2_t
   CALL VecWAXPY   ( k2_t_PETSC, 0.5d0*dt, k1_m_PETSC, x2_PETSC, IERR )
! 5. compute k2_m
   CALL VecWAXPY   ( k2_m_PETSC, 0.5d0*dt, k1_b_PETSC, x3_PETSC, IERR )
! 6. compute k2_b
   CALL VecWAXPY   ( help_PETSC, 0.5d0*dt, k1_t_PETSC, x1_PETSC, IERR )
   CALL MatMultAdd ( AC_PETSC, k2_m_PETSC, k2_b_PETSC, k2_b_PETSC, IERR )
   CALL MatMultAdd ( AK_PETSC, k2_t_PETSC, k2_b_PETSC, k2_b_PETSC, IERR )
   CALL MatMultAdd ( AG_PETSC, help_PETSC, k2_b_PETSC, k2_b_PETSC, IERR )
! 7. compute k3_t
   CALL VecWAXPY   ( k3_t_PETSC, 0.5d0*dt, k2_m_PETSC, x2_PETSC, IERR )
! 8. compute k3_m
   CALL VecWAXPY   ( k3_m_PETSC, 0.5d0*dt, k2_b_PETSC, x3_PETSC, IERR )
! 9. compute k3_b
   CALL VecWAXPY   ( help_PETSC, 0.5d0*dt, k2_t_PETSC, x1_PETSC, IERR )
   CALL MatMultAdd ( AC_PETSC, k3_m_PETSC, k3_b_PETSC, k3_b_PETSC, IERR )
   CALL MatMultAdd ( AK_PETSC, k3_t_PETSC, k3_b_PETSC, k3_b_PETSC, IERR )
   CALL MatMultAdd ( AG_PETSC, help_PETSC, k3_b_PETSC, k3_b_PETSC, IERR )
! 10. compute k4_t
   CALL VecWAXPY   ( k4_t_PETSC, dt, k3_m_PETSC, x2_PETSC, IERR )
! 11. compute k4_m
   CALL VecWAXPY   ( k4_m_PETSC, dt, k3_b_PETSC, x3_PETSC, IERR )
! 12. compute k4_b
   CALL VecWAXPY   ( help_PETSC, dt, k3_t_PETSC, x1_PETSC, IERR )
   CALL MatMultAdd ( AC_PETSC, k4_m_PETSC, k4_b_PETSC, k4_b_PETSC, IERR )
   CALL MatMultAdd ( AK_PETSC, k4_t_PETSC, k4_b_PETSC, k4_b_PETSC, IERR )
   CALL MatMultAdd ( AG_PETSC, help_PETSC, k4_b_PETSC, k4_b_PETSC, IERR )
! 13. compute x1_n+1
   CALL VecSet  ( help_PETSC, 0.0d0, ierr)
   CALL VecMAXPY( help_PETSC, 4, [1.0d0, 2.0d0, 2.0d0, 1.0d0], [k1_t_PETSC, k2_t_PETSC, k3_t_PETSC, k4_t_PETSC], IERR )
   CALL VecAXPY ( x1_PETSC, dt/6.0d0, help_PETSC, IERR )
! 14. compute x2_n+1
   CALL VecSet  ( help_PETSC, 0.0d0, ierr)
   CALL VecMAXPY( help_PETSC, 4, [1.0d0, 2.0d0, 2.0d0, 1.0d0], [k1_m_PETSC, k2_m_PETSC, k3_m_PETSC, k4_m_PETSC], IERR )
   CALL VecAXPY ( x2_PETSC, dt/6.0d0, help_PETSC, IERR )
! 15. compute x3_n+1
   CALL VecSet  ( help_PETSC, 0.0d0, ierr)
   CALL VecMAXPY( help_PETSC, 4, [1.0d0, 2.0d0, 2.0d0, 1.0d0], [k1_b_PETSC, k2_b_PETSC, k3_b_PETSC, k4_b_PETSC], IERR )
   CALL VecAXPY ( x3_PETSC, dt/6.0d0, help_PETSC, IERR )

END DO
!---------- ---------- ---------- ---------- ----------

  If ( Rank == 0 ) Then ; ! End top file wraper

    UnFile = UN_OutallWrp ;
    Write (Unit = UnFile, FMT = "(A11)", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) '    </Grid>' ;
    Write (Unit = UnFile, FMT = "(A11)", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) '  </Domain>' ;
    Write (Unit = UnFile, FMT = "(A7 )", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) '</Xdmf>' ;

  End If ;



CALL VecDestroy ( x1_PETSC  , IERR )
CALL VecDestroy ( x2_PETSC  , IERR )
CALL VecDestroy ( x3_PETSC  , IERR )

CALL VecDestroy ( k1_t_PETSC, IERR )
CALL VecDestroy ( k1_m_PETSC, IERR )
CALL VecDestroy ( k1_b_PETSC, IERR )

CALL VecDestroy ( k2_t_PETSC, IERR )
CALL VecDestroy ( k2_m_PETSC, IERR )
CALL VecDestroy ( k2_b_PETSC, IERR )

CALL VecDestroy ( k3_t_PETSC, IERR )
CALL VecDestroy ( k3_m_PETSC, IERR )
CALL VecDestroy ( k3_b_PETSC, IERR )

CALL VecDestroy ( k4_t_PETSC, IERR )
CALL VecDestroy ( k4_m_PETSC, IERR )
CALL VecDestroy ( k4_b_PETSC, IERR )

CALL VecDestroy ( help_PETSC, IERR )
CALL VecDestroy ( HELP_1_PETSC, IERR )

Call VecDestroy   ( U_Rank_PETSC,   ErrPTC ) ;
Call VecDestroy   ( UD_Rank_PETSC,  ErrPTC ) ;
Call VecDestroy   ( UDD_Rank_PETSC, ErrPTC ) ;

! - Destroying PETSc Index Set ----------------------------------------------------------------------------------------------------------------------
Call ISDestroy(IS_from, ErrPTC );
Call ISDestroy(IS_to,   ErrPTC );

Call VecScatterDestroy(vscat_U, ErrPTC);
!Call VecScatterDestroy(vscat_UD, ErrPTC);
!Call VecScatterDestroy(vscat_UDD, ErrPTC);

  If ( Rank == 0 ) Then ;
    UnFile = UN_OutallWrp ;
    Close ( Unit = UnFile, STATUS = 'KEEP', ERR = 1002, IOSTAT = IO_File ) ;
  End If ;

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
! Explicit time stepping based on Runge-Kutta method, 4 stage-4th order, PML-3D: adjoint problem, Optimize-then-Discretize
!************************************************
! REVISION : T 8 July 2014

SUBROUTINE EXPLICIT_RK4_PETSC_PML_3D_Adjoint_OtD ( AK_PETSC, DIAG_M_PETSC, DIAG_iM_PETSC, AC_PETSC, AG_PETSC, B_mis_PETSC, AK_RD_PETSC, DIAG_M_RD_PETSC, ID, &
                                                   NNDH, EqDis, U_Store_Mapping, P_Store_Mapping, NStore_Mapping, U_Store_Numbers_Global, EqDis_MAP_U_Store_Numbers_Global, Dis_meas )

USE PARAMETERS
USE IFPORT
IMPLICIT DOUBLE PRECISION (A-H,O-Z)

! - Type DECLERATIONS -------------------------------------------------------------------------------------------------------------------------------
Integer, Intent(In),  Dimension ( NNDH * NDim )              :: EqDis
Integer, Intent(In),  Dimension ( NNDH * NDim )              :: EqDis_MAP_U_Store_Numbers_Global
Integer, Intent(In),  Dimension ( NStore_Mapping )           ::           U_Store_Numbers_Global
Real(8), Intent(In),  Dimension ( NStore_Mapping , 0:NStep ) :: U_Store_Mapping
Real(8), Intent(Out), Dimension ( NStore_Mapping , 0:NStep ) :: P_Store_Mapping
!==================================================================================================
! - PETSC INCLUDE FILES
#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"
#include "finclude/petscis.h"

! - PETSC VARIABLES AND MATRICES
Mat            :: AK_PETSC, AC_PETSC, AG_PETSC, AK_RD_PETSC
Vec            :: B_mis_PETSC
Vec            :: DIAG_M_PETSC, DIAG_iM_PETSC, DIAG_M_RD_PETSC

PetscErrorCode :: IERR, ErrPTC
PetscMPIInt    :: SIZE, RANK

! - PETSc local vectors
Vec            ::   x1_PETSC,   x2_PETSC,   x3_PETSC
Vec            :: k1_t_PETSC, k1_m_PETSC, k1_b_PETSC, k2_t_PETSC, k2_m_PETSC, k2_b_PETSC, help_PETSC
Vec            :: k3_t_PETSC, k3_m_PETSC, k3_b_PETSC, k4_t_PETSC, k4_m_PETSC, k4_b_PETSC
Vec            :: HELP_1_PETSC

PetscScalar    :: u_misfit_PETSC ( NNDH * NDim )
PetscScalar, pointer :: P(:)
PetscScalar, pointer :: X(:)
!==================================================================================================


! in
!---------- ---------- ---------- ---------- ----------
DIMENSION ID ( NJ, NDOF )
Dimension Dis_meas ( 0:NStep , NNDH * NDim )               ! Measured displacements at select sensor locations
!---------- ---------- ---------- ---------- ----------


! local
!---------- ---------- ---------- ---------- ----------
Dimension u_misfit ( NNDH * NDim )
Dimension p_small ( NStore_Mapping )
!---------- ---------- ---------- ---------- ----------


I_PRINT_STEP = 1000
CALL MPI_Comm_rank   ( PETSC_COMM_WORLD, RANK, IERR )
IF (RANK == 0 .AND. I_PRINT_SCREEN == 1) WRITE(*,*) '3D RK4 Adjoint'


!==================================================================================================
! - DEFINE PETSC OBJECTS
CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, x1_PETSC, IERR )
CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, x2_PETSC, IERR )
CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, x3_PETSC, IERR )

CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, k1_t_PETSC, IERR )
CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, k1_m_PETSC, IERR )
CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, k1_b_PETSC, IERR )

CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, k2_t_PETSC, IERR )
CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, k2_m_PETSC, IERR )
CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, k2_b_PETSC, IERR )

CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, k3_t_PETSC, IERR )
CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, k3_m_PETSC, IERR )
CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, k3_b_PETSC, IERR )

CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, k4_t_PETSC, IERR )
CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, k4_m_PETSC, IERR )
CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, k4_b_PETSC, IERR )

CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, help_PETSC, IERR )
CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, HELP_1_PETSC, IERR )

CALL VecSet       ( x1_PETSC, 0.0d0, ierr)
CALL VecSet       ( x2_PETSC, 0.0d0, ierr)
CALL VecSet       ( x3_PETSC, 0.0d0, ierr)
! just for hygenic concerns:
CALL VecSet       ( k1_m_PETSC, 0.0d0, ierr)
CALL VecSet       ( k2_m_PETSC, 0.0d0, ierr)
CALL VecSet       ( k3_m_PETSC, 0.0d0, ierr)
CALL VecSet       ( k4_m_PETSC, 0.0d0, ierr)
!==================================================================================================
! Comment: M, C, K, G, have already been formed for the solution of the forward problem; left-multiplied by mass inverse, and except M, multiplied by a minus sign.
!          for the adjoint problem, only K should retain its minus sign. C and G are positive when moved to the right-hand-side. Hence, we hit them with a minus sign to fix them.

CALL MatScale ( AC_PETSC, -1.0d0, IERR )
CALL MatScale ( AG_PETSC, -1.0d0, IERR )

! form system matrices (cancel left-scaling, do right-scaling)
!---------- ---------- ---------- ---------- ----------
CALL MatDiagonalScale ( AK_PETSC, DIAG_M_PETSC, DIAG_iM_PETSC, IERR )
CALL MatDiagonalScale ( AC_PETSC, DIAG_M_PETSC, DIAG_iM_PETSC, IERR )
CALL MatDiagonalScale ( AG_PETSC, DIAG_M_PETSC, DIAG_iM_PETSC, IERR )


CALL MPI_Comm_size   ( PETSC_COMM_WORLD, SIZE, IERR )
CALL MPI_Comm_rank   ( PETSC_COMM_WORLD, RANK, IERR )
Call VecGetOwnershipRange ( x2_PETSC, IStart, IEnd, IERR )


! Explicit Integration
!---------- ---------- ---------- ---------- ----------
DO ISTEP = NSTEP, 1, -1

   T = DBLE( ISTEP ) * DT
   IF ( RANK == 0 .AND. MOD ( ISTEP+1 , I_PRINT_STEP ) == 0 ) THEN
      WRITE (*,'(A30,F10.5)') 'TIME=', T
   END IF

! -------------------------------------------------------------------------------------------------
! this is to see if adjoint diverges or not
!  Call VecDot ( x2_PETSC, x2_PETSC, energy, IERR ); if (rank ==0 ) write(*,*) energy
! -------------------------------------------------------------------------------------------------


! storage of the adjoint variable
!---------- ---------- ---------- ---------- ---------- ---------------------------------------------------------------------------------------------
! 1- selection: during time integration, access p(big), pick p(small) entries.
   Call VecGetArrayF90 ( x2_PETSC, P, IERR )
   Do I = 1, NStore_Mapping
     p_small ( I ) = P ( U_Store_Numbers_Global ( I ) - IStart + 1 )
   End Do
   Call VecRestoreArrayF90 ( x2_PETSC, P, IERR )

! 2- store p(small) entries in the dense matrix: We are storing 0:N-1.
   P_Store_Mapping ( : , ISTEP ) = p_small (:)
!---------- ---------- ---------- ---------- ---------- ---------------------------------------------------------------------------------------------


! right-hand-side force at 4 stages (misfit vector)
!---------- ---------- ---------- ---------- ---------- ---------------------------------------------------------------------------------------------
! 3.1- form the misfit for the time step n
   Do IEq = 1, NNDH * NDim
      u_misfit ( IEq ) = U_Store_Mapping ( EqDis_MAP_U_Store_Numbers_Global ( IEq ) , IStep ) - Dis_meas ( IStep , IEq )
   End Do

! 3.2- initialize the misfit vector
   Call VecSet ( B_mis_PETSC, 0.0d0, IERR )

! 3.3- assemble the distributed misfit vector
   u_misfit_PETSC = u_misfit
   If ( NNDH /= 0 ) CALL VecSetValues ( B_mis_PETSC, NNDH * NDim, EqDis, u_misfit_PETSC, ADD_VALUES, IERR )
   Call VecAssemblyBegin ( B_mis_PETSC, IERR ) ; Call VecAssemblyEnd ( B_mis_PETSC, IERR )

! 3.4- Q is scaled with the inverse of the mass matrix
   CALL VecPointwiseMult ( B_mis_PETSC, B_mis_PETSC, DIAG_iM_PETSC, IERR )

! 3.5- Send it to the right vector
   CALL VecSet  ( k1_b_PETSC, 0.0d0             , IERR )
   CALL VecAXPY ( k1_b_PETSC, 1.0d0, B_mis_PETSC, IERR )

! 4.1- form the misfit for the time step n-1
   Do IEq = 1, NNDH * NDim
      u_misfit ( IEq ) = U_Store_Mapping ( EqDis_MAP_U_Store_Numbers_Global ( IEq ) , IStep-1 ) - Dis_meas ( IStep-1 , IEq )
   End Do

! 4.2- initialize the misfit vector
   Call VecSet ( B_mis_PETSC, 0.0d0, IERR )

! 4.3- assemble the distributed misfit vector
   u_misfit_PETSC = u_misfit
   If ( NNDH /= 0 ) CALL VecSetValues ( B_mis_PETSC, NNDH * NDim, EqDis, u_misfit_PETSC, ADD_VALUES, IERR )
   Call VecAssemblyBegin ( B_mis_PETSC, IERR ) ; Call VecAssemblyEnd ( B_mis_PETSC, IERR )

! 4.4- Q is scaled with the inverse of the mass matrix
   CALL VecPointwiseMult ( B_mis_PETSC, B_mis_PETSC, DIAG_iM_PETSC, IERR )

! 4.5- Send it to the right vector
   CALL VecSet  ( k4_b_PETSC, 0.0d0             , IERR )
   CALL VecAXPY ( k4_b_PETSC, 1.0d0, B_mis_PETSC, IERR )

! 5.1- use average for the middle steps
   CALL VecSet  ( k2_b_PETSC, 0.0d0, IERR )
   CALL VecSet  ( k3_b_PETSC, 0.0d0, IERR )
   Call VecMAXPY( k2_b_PETSC, 2, [0.5d0, 0.5d0], [k1_b_PETSC, k4_b_PETSC], IERR )
   Call VecMAXPY( k3_b_PETSC, 2, [0.5d0, 0.5d0], [k1_b_PETSC, k4_b_PETSC], IERR )
!---------- ---------- ---------- ---------- ---------- ---------------------------------------------------------------------------------------------


! RK-4 Explicit time-stepping:
!---------- ---------- ---------- ---------- ----------
! 1. compute k1_t
   CALL VecCopy ( x2_PETSC, k1_t_PETSC, IERR )
! 2. compute k1_m
   CALL VecCopy ( x3_PETSC, k1_m_PETSC, IERR )
! 3. compute k1_b
   CALL MatMultTransposeAdd ( AC_PETSC, x3_PETSC, k1_b_PETSC, k1_b_PETSC, IERR )
   CALL MatMultTransposeAdd ( AK_PETSC, x2_PETSC, k1_b_PETSC, k1_b_PETSC, IERR )
   CALL MatMultTransposeAdd ( AG_PETSC, x1_PETSC, k1_b_PETSC, k1_b_PETSC, IERR )
! 4. compute k2_t
   CALL VecWAXPY   ( k2_t_PETSC,-0.5d0*dt, k1_m_PETSC, x2_PETSC, IERR )
! 5. compute k2_m
   CALL VecWAXPY   ( k2_m_PETSC,-0.5d0*dt, k1_b_PETSC, x3_PETSC, IERR )
! 6. compute k2_b
   CALL VecWAXPY   ( help_PETSC,-0.5d0*dt, k1_t_PETSC, x1_PETSC, IERR )
   CALL MatMultTransposeAdd ( AC_PETSC, k2_m_PETSC, k2_b_PETSC, k2_b_PETSC, IERR )
   CALL MatMultTransposeAdd ( AK_PETSC, k2_t_PETSC, k2_b_PETSC, k2_b_PETSC, IERR )
   CALL MatMultTransposeAdd ( AG_PETSC, help_PETSC, k2_b_PETSC, k2_b_PETSC, IERR )
! 7. compute k3_t
   CALL VecWAXPY   ( k3_t_PETSC,-0.5d0*dt, k2_m_PETSC, x2_PETSC, IERR )
! 8. compute k3_m
   CALL VecWAXPY   ( k3_m_PETSC,-0.5d0*dt, k2_b_PETSC, x3_PETSC, IERR )
! 9. compute k3_b
   CALL VecWAXPY   ( help_PETSC,-0.5d0*dt, k2_t_PETSC, x1_PETSC, IERR )
   CALL MatMultTransposeAdd ( AC_PETSC, k3_m_PETSC, k3_b_PETSC, k3_b_PETSC, IERR )
   CALL MatMultTransposeAdd ( AK_PETSC, k3_t_PETSC, k3_b_PETSC, k3_b_PETSC, IERR )
   CALL MatMultTransposeAdd ( AG_PETSC, help_PETSC, k3_b_PETSC, k3_b_PETSC, IERR )
! 10. compute k4_t
   CALL VecWAXPY   ( k4_t_PETSC, -dt, k3_m_PETSC, x2_PETSC, IERR )
! 11. compute k4_m
   CALL VecWAXPY   ( k4_m_PETSC, -dt, k3_b_PETSC, x3_PETSC, IERR )
! 12. compute k4_b
   CALL VecWAXPY   ( help_PETSC, -dt, k3_t_PETSC, x1_PETSC, IERR )
   CALL MatMultTransposeAdd ( AC_PETSC, k4_m_PETSC, k4_b_PETSC, k4_b_PETSC, IERR )
   CALL MatMultTransposeAdd ( AK_PETSC, k4_t_PETSC, k4_b_PETSC, k4_b_PETSC, IERR )
   CALL MatMultTransposeAdd ( AG_PETSC, help_PETSC, k4_b_PETSC, k4_b_PETSC, IERR )
! 13. compute x1_n-1
   CALL VecSet  ( help_PETSC, 0.0d0, ierr)
   CALL VecMAXPY( help_PETSC, 4, [1.0d0, 2.0d0, 2.0d0, 1.0d0], [k1_t_PETSC, k2_t_PETSC, k3_t_PETSC, k4_t_PETSC], IERR )
   CALL VecAXPY ( x1_PETSC, -dt/6.0d0, help_PETSC, IERR )
! 14. compute x2_n-1
   CALL VecSet  ( help_PETSC, 0.0d0, ierr)
   CALL VecMAXPY( help_PETSC, 4, [1.0d0, 2.0d0, 2.0d0, 1.0d0], [k1_m_PETSC, k2_m_PETSC, k3_m_PETSC, k4_m_PETSC], IERR )
   CALL VecAXPY ( x2_PETSC, -dt/6.0d0, help_PETSC, IERR )
! 15. compute x3_n-1
   CALL VecSet  ( help_PETSC, 0.0d0, ierr)
   CALL VecMAXPY( help_PETSC, 4, [1.0d0, 2.0d0, 2.0d0, 1.0d0], [k1_b_PETSC, k2_b_PETSC, k3_b_PETSC, k4_b_PETSC], IERR )
   CALL VecAXPY ( x3_PETSC, -dt/6.0d0, help_PETSC, IERR )


! =================================================================================================
! =================================================================================================
!   Call VecGetArrayF90 ( x2_PETSC, X, IERR )
!      if ( rank == 12 ) then
!         write ( 72 , '(2e20.10)' ) X ( 1 ) !, P_Bigger_Rank ( 1 )
!      end if
!   Call VecRestoreArrayF90 ( x2_PETSC, X, IERR )
! =================================================================================================
! =================================================================================================


END DO
!---------- ---------- ---------- ---------- ----------


CALL VecDestroy ( x1_PETSC  , IERR )
CALL VecDestroy ( x2_PETSC  , IERR )
CALL VecDestroy ( x3_PETSC  , IERR )

CALL VecDestroy ( k1_t_PETSC, IERR )
CALL VecDestroy ( k1_m_PETSC, IERR )
CALL VecDestroy ( k1_b_PETSC, IERR )

CALL VecDestroy ( k2_t_PETSC, IERR )
CALL VecDestroy ( k2_m_PETSC, IERR )
CALL VecDestroy ( k2_b_PETSC, IERR )

CALL VecDestroy ( k3_t_PETSC, IERR )
CALL VecDestroy ( k3_m_PETSC, IERR )
CALL VecDestroy ( k3_b_PETSC, IERR )

CALL VecDestroy ( k4_t_PETSC, IERR )
CALL VecDestroy ( k4_m_PETSC, IERR )
CALL VecDestroy ( k4_b_PETSC, IERR )

CALL VecDestroy ( help_PETSC, IERR )
!CALL VecDestroy ( HELP_1_PETSC, IERR )

END


!************************************************
! The control problem: Optimize-then-Discretize
!************************************************
! REVISION : Th 3 July 2014

SUBROUTINE EXPLICIT_RK4_PETSC_PML_3D_Control_OtD ( U_Store_Mapping, P_Store_Mapping, NStore_Mapping, U_Store_Numbers_Global, Ux_big_PETSc, Ux_Bigger_Rank_PETSc, Uk_big_PETSc, P_Bigger_Rank_PETSc, &
                                                   vscat_u_hist, g_Lambda_PETSC, g_Mu_PETSC, idx_Mat_to, NJ_Rank_IParts, XYZ, INOD, NGP, MTEL, ID, RANK )

USE PARAMETERS
USE IFPORT
IMPLICIT DOUBLE PRECISION (A-H,O-Z)

! - Type DECLERATIONS -------------------------------------------------------------------------------------------------------------------------------
Integer, Intent(In), Dimension ( NStore_Mapping )           ::    U_Store_Numbers_Global
Real(8), Intent(In), Dimension ( NStore_Mapping , 0:NStep ) ::    U_Store_Mapping
Real(8), Intent(In), Dimension ( NStore_Mapping , 0:NStep ) ::    P_Store_Mapping
!==================================================================================================
! - PETSC INCLUDE FILES
#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"
#include "finclude/petscis.h"

! - PETSC VARIABLES AND MATRICES
PetscErrorCode :: IERR, ErrPTC
PetscMPIInt    :: SIZE, RANK

! - Control Problem ----------------------------------------------------------------------------------------*-----------------------------------------
Vec            :: Ux_big_PETSc,         Uk_big_PETSc
Vec            :: Ux_Bigger_Rank_PETSc, P_Bigger_Rank_PETSc
Vec            :: g_Lambda_PETSC, g_Mu_PETSC

PetscScalar    :: u_small_PETSC ( NStore_Mapping )
PetscScalar    :: p_small_PETSC ( NStore_Mapping )

PetscScalar, pointer :: U_Bigger_Rank(:), P_Bigger_Rank(:)

VecScatter     :: vscat_u_hist
!==================================================================================================


! in
!---------- ---------- ---------- ---------- ----------
DIMENSION XYZ ( NJ, NDIM ), INOD ( NNODE, NEL ), NGP ( NEL ), MTEL ( NEL ), ID ( NJ, NDOF )
!---------- ---------- ---------- ---------- ----------


! in
!---------- ---------- ---------- ---------- ----------
Dimension idx_Mat_to ( NJ_Rank_IParts )                    ! Local equation numbers for materials, C-indexing
!---------- ---------- ---------- ---------- ----------


! local
!---------- ---------- ---------- ---------- ----------
Dimension u_small ( NStore_Mapping )
Dimension p_small ( NStore_Mapping )
!---------- ---------- ---------- ---------- ----------
I_PRINT_STEP = 1000


IF ( RANK == 0 .AND. I_PRINT_SCREEN == 1 ) WRITE(*,*) '3D Control'


! Integration ---------------------------------------------------------------------------------------------------------------------------------------
DO ISTEP = 1 , NSTEP-1 ! at t_{0}, u=0, and at t_{NSTEP}, p = 0. Trapezoidal rule.

   T = DBLE( ISTEP ) * DT
   IF ( RANK == 0 .AND. MOD ( ISTEP+1 , I_PRINT_STEP ) == 0 ) THEN
      WRITE (*,'(A30,F10.5)') 'TIME=', T
   END IF

! ggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggggg
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! -------------------------------------------------------------------------------------------------
! a- retrieve from the dense matrix
   u_small (:) = U_Store_Mapping ( : , ISTEP )
   p_small (:) = P_Store_Mapping ( : , ISTEP )
! b- initlalize the parallel vector
   Call VecSet ( Ux_big_PETSc, 0.0d0, IERR )
   Call VecSet ( Uk_big_PETSc, 0.0d0, IERR )
! c- move from the dense matrix to parallel vector (long-thin vector)
   u_small_PETSC = u_small
   p_small_PETSC = p_small
   If ( NStore_Mapping /= 0 ) Then
      Call VecSetValues ( Ux_big_PETSc, NStore_Mapping, U_Store_Numbers_Global, u_small_PETSC, ADD_VALUES, IERR )
      Call VecSetValues ( Uk_big_PETSc, NStore_Mapping, U_Store_Numbers_Global, p_small_PETSC, ADD_VALUES, IERR )
   End If
! d- this is not needed since everything is on the right processor, but let's do it.
   Call VecAssemblyBegin ( Ux_big_PETSc, IERR ) ; Call VecAssemblyEnd ( Ux_big_PETSc, IERR )
   Call VecAssemblyBegin ( Uk_big_PETSc, IERR ) ; Call VecAssemblyEnd ( Uk_big_PETSc, IERR )
! e- scatter state variables
   Call VecScatterBegin ( vscat_u_hist, Ux_big_PETSc, Ux_Bigger_Rank_PETSc, INSERT_VALUES, SCATTER_FORWARD, IERR )
   Call VecScatterEnd   ( vscat_u_hist, Ux_big_PETSc, Ux_Bigger_Rank_PETSc, INSERT_VALUES, SCATTER_FORWARD, IERR )
! f- scatter adjoint variables
   Call VecScatterBegin ( vscat_u_hist, Uk_big_PETSc,  P_Bigger_Rank_PETSc, INSERT_VALUES, SCATTER_FORWARD, IERR )
   Call VecScatterEnd   ( vscat_u_hist, Uk_big_PETSc,  P_Bigger_Rank_PETSc, INSERT_VALUES, SCATTER_FORWARD, IERR )
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

END DO
! End Integration loop ------------------------------------------------------------------------------------------------------------------------------

Return
END


