!************************************************
! ENERGY DECAY IN REGULAR DOMAIN
!************************************************
! REVISION : T, 26 June 2012

SUBROUTINE TOTAL_ENERGY_RD_PETSC ( AK_RD_PETSC, AM_RD_PETSC, U_PETSC, UD_PETSC, ENERGY )
USE PARAMETERS
IMPLICIT DOUBLE PRECISION (A-H,O-Z)


!==================================================================================================
! - PETSC INCLUDE FILES
#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscmat.h"

! - PETSC VARIABLES AND MATRICES
Mat            :: AK_RD_PETSC, AM_RD_PETSC

PetscErrorCode :: IERR
PetscMPIInt    :: SIZE, RANK
PetscScalar    :: ENERGY_1, ENERGY_2, ENERGY

! - PETSc local vectors
Vec            :: U_PETSC, UD_PETSC
Vec            :: HELP_PETSC
!==================================================================================================


!==================================================================================================
! - DEFINE PETSC OBJECTS
CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, HELP_PETSC   , IERR )
!==================================================================================================


! AK_RD_U  = MATMUL(AK_RD_f, U)
CALL MatMult ( AK_RD_PETSC, U_PETSC, HELP_PETSC, IERR )
! ENERGY_1 = DOT_PRODUCT(U , AK_RD_U)
CALL VecDot  ( U_PETSC, HELP_PETSC, ENERGY_1, IERR )

! AM_RD_UD = MATMUL(AM_RD_f, UD)
CALL MatMult ( AM_RD_PETSC, UD_PETSC, HELP_PETSC, IERR )
! ENERGY_2 = DOT_PRODUCT(UD, AM_RD_UD)
CALL VecDot  ( UD_PETSC, HELP_PETSC, ENERGY_2, IERR )

ENERGY   = ( ENERGY_1 + ENERGY_2 ) / 2.0d0

CALL VecDestroy ( HELP_PETSC, IERR )

RETURN
END


!************************************************
! ENERGY DECAY IN REGULAR DOMAIN: Explicit method
!************************************************
! REVISION : M, 16 June 2014
! My father considered a walk among the mountains as the equivalent of churchgoing.

SUBROUTINE TOTAL_ENERGY_RD_EXPLICIT_Central_Difference ( AK_RD_PETSC, DIAG_M_RD_PETSC, U_nn_PETSC, U_nm_PETSC, HELP_1_PETSC, HELP_2_PETSC, ENERGY )
USE PARAMETERS
IMPLICIT DOUBLE PRECISION (A-H,O-Z)

!==================================================================================================
! - PETSC INCLUDE FILES
#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscmat.h"

! - PETSC VARIABLES AND MATRICES
Mat            :: AK_RD_PETSC
Vec            :: DIAG_M_RD_PETSC
Vec            :: U_nn_PETSC, U_nm_PETSC

PetscErrorCode :: IERR
PetscMPIInt    :: SIZE, RANK
PetscScalar    :: ENERGY_1, ENERGY_2, ENERGY

! - PETSc local vectors
Vec            :: HELP_1_PETSC, HELP_2_PETSC
!==================================================================================================

! stiffness
CALL MatMult ( AK_RD_PETSC, U_nn_PETSC, HELP_1_PETSC, IERR )
CALL VecDot  ( U_nn_PETSC, HELP_1_PETSC, ENERGY_1, IERR )

! mass
! 1- compute the velocity vector
CALL VecCopy ( U_nn_PETSC, HELP_1_PETSC, IERR )
Call VecAXPY ( HELP_1_PETSC, -1.0d0, U_nm_PETSC, IERR )
Call VecScale ( HELP_1_PETSC, 1.0d0 / DT, IERR )

CALL VecPointwiseMult ( HELP_2_PETSC, DIAG_M_RD_PETSC, HELP_1_PETSC, IERR )
CALL VecDot ( HELP_2_PETSC, HELP_1_PETSC, ENERGY_2, IERR )

ENERGY = 0.50D0 * ( ENERGY_1 + ENERGY_2 )


RETURN
END


!************************************************
! ENERGY DECAY IN REGULAR DOMAIN: Explicit method
!************************************************
! REVISION : W, 17 April 2012

SUBROUTINE TOTAL_ENERGY_RD_PETSC_EXPLICIT ( AK_RD_PETSC, DIAG_M_RD_PETSC, U_PETSC, UD_PETSC, HELP_1_PETSC, ENERGY )
USE PARAMETERS
IMPLICIT DOUBLE PRECISION (A-H,O-Z)


!==================================================================================================
! - PETSC INCLUDE FILES
#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscmat.h"

! - PETSC VARIABLES AND MATRICES
Mat            :: AK_RD_PETSC
Vec            :: DIAG_M_RD_PETSC
Vec            :: U_PETSC, UD_PETSC

PetscErrorCode :: IERR
PetscMPIInt    :: SIZE, RANK
PetscScalar    :: ENERGY_1, ENERGY_2, ENERGY

! - PETSc local vectors
Vec            :: HELP_1_PETSC
!==================================================================================================

! stiffness
CALL MatMult ( AK_RD_PETSC, U_PETSC, HELP_1_PETSC, IERR )
CALL VecDot  ( U_PETSC, HELP_1_PETSC, ENERGY_1, IERR )

! mass
CALL VecPointwiseMult ( HELP_1_PETSC, DIAG_M_RD_PETSC, UD_PETSC, IERR )
CALL VecDot ( HELP_1_PETSC, UD_PETSC, ENERGY_2, IERR )

ENERGY = 0.50D0 * ( ENERGY_1 + ENERGY_2 )


RETURN
END


!************************************************
! 3D-PML NEWMARK: implicit time-stepping (for 3rd order ODEs)
! We re-write the 3rd order system as a 2nd order system and apply the standard Newmark method.
! implementation is not complete; NEST matrix needs to be (converted and) factorized which is not implemented in PETSc yet.
!************************************************
! REVISION : W 25 July 2012

SUBROUTINE NEWMARK_PETSC_4_does_not_work ( AK_PETSC, AM_PETSC, AC_PETSC, AG_PETSC, B_PETSC, ID, LTRANS, RANK )
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
Mat            :: AK_PETSC, AM_PETSC, AC_PETSC, AG_PETSC
Vec            :: B_PETSC

Mat            :: I_p_PETSC, I_m_PETSC
Mat            :: AK_b_PETSC, AM_b_PETSC, AC_b_PETSC
Mat            :: AF_PETSC, AL_PETSC

Vec            :: Vdd_PETSC, Vd_PETSC, V_PETSC, Udd_PETSC, Ud_PETSC, U_PETSC
Vec            :: Xdd_b_PETSC, Xd_b_PETSC, X_b_PETSC, B_b_PETSC, R_b_PETSC
Vec            :: HELP_PETSC, Xdd_OLD_b_PETSC, zero_PETSC

IS             :: is_1, is_2, perm
MatFactorInfo  :: mat_factor_info(MAT_FACTORINFO_SIZE)
PetscInt       :: m

PetscErrorCode :: IERR
PetscMPIInt    :: SIZE, RANK

PetscScalar    :: A0, A1, A2, A3, A4, A5, FAC1
!==================================================================================================


! in 
!---------- ---------- -------- -- ---------- ----------
DIMENSION ID(NJ, NDOF), LTRANS(NTRANS)
!---------- ---------- ---------- ---------- ----------


! matrices
!---------- ---------- ---------- ---------- ----------
CALL MatCreateAIJ ( PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, NEQ, NEQ, 1, PETSC_NULL_INTEGER, 1, PETSC_NULL_INTEGER, I_p_PETSC, IERR )
CALL MatCreateAIJ ( PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, NEQ, NEQ, 1, PETSC_NULL_INTEGER, 1, PETSC_NULL_INTEGER, I_m_PETSC, IERR )

DO I = 1, NEQ
   CALL MatSetValue ( I_p_PETSC, I-1, I-1, +1.0d0, ADD_VALUES, IERR )
   CALL MatSetValue ( I_m_PETSC, I-1, I-1, -1.0d0, ADD_VALUES, IERR )
END DO

CALL MatAssemblyBegin ( I_p_PETSC, MAT_FINAL_ASSEMBLY, IERR )
CALL MatAssemblyBegin ( I_m_PETSC, MAT_FINAL_ASSEMBLY, IERR )
CALL MatAssemblyEnd   ( I_p_PETSC, MAT_FINAL_ASSEMBLY, IERR )
CALL MatAssemblyEnd   ( I_m_PETSC, MAT_FINAL_ASSEMBLY, IERR )

! factored matrix for direct solver; see if you can factor in place and get rid of this matrix
CALL MatCreateAIJ ( PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, 2*NEQ, 2*NEQ, 100, PETSC_NULL_INTEGER, 100, PETSC_NULL_INTEGER, AF_PETSC, IERR )


! vectors
!---------- ---------- ---------- ---------- ----------
CALL VecCreateMPI ( PETSC_COMM_WORLD, PETSC_DECIDE, NEQ, Vdd_PETSC, IERR )
CALL VecCreateMPI ( PETSC_COMM_WORLD, PETSC_DECIDE, NEQ,  Vd_PETSC, IERR )
CALL VecCreateMPI ( PETSC_COMM_WORLD, PETSC_DECIDE, NEQ,   V_PETSC, IERR )
CALL VecCreateMPI ( PETSC_COMM_WORLD, PETSC_DECIDE, NEQ, Udd_PETSC, IERR )
CALL VecCreateMPI ( PETSC_COMM_WORLD, PETSC_DECIDE, NEQ,  Ud_PETSC, IERR )
CALL VecCreateMPI ( PETSC_COMM_WORLD, PETSC_DECIDE, NEQ,   U_PETSC, IERR )
CALL VecCreateMPI ( PETSC_COMM_WORLD, PETSC_DECIDE, NEQ,zero_PETSC, IERR )

CALL VecCreateMPI ( PETSC_COMM_WORLD, PETSC_DECIDE, 2*NEQ,       R_b_PETSC, IERR )
CALL VecCreateMPI ( PETSC_COMM_WORLD, PETSC_DECIDE, 2*NEQ,      HELP_PETSC, IERR )
CALL VecCreateMPI ( PETSC_COMM_WORLD, PETSC_DECIDE, 2*NEQ,      DX_b_PETSC, IERR )
CALL VecCreateMPI ( PETSC_COMM_WORLD, PETSC_DECIDE, 2*NEQ, Xdd_OLD_b_PETSC, IERR )


! system is at rest
!---------- ---------- ---------- ---------- ----------
CALL VecSet ( Vdd_PETSC, 0.0d0, IERR )
CALL VecSet (  Vd_PETSC, 0.0d0, IERR )
CALL VecSet (   V_PETSC, 0.0d0, IERR )
CALL VecSet ( Udd_PETSC, 0.0d0, IERR )
CALL VecSet (  Ud_PETSC, 0.0d0, IERR )
CALL VecSet (   U_PETSC, 0.0d0, IERR )
CALL VecSet (zero_PETSC, 0.0d0, IERR )


! index sets
!---------- ---------- ---------- ---------- ----------
CALL ISCreateStride ( PETSC_COMM_WORLD, NEQ, 0  , 1, is_1, IERR )
CALL ISCreateStride ( PETSC_COMM_WORLD, NEQ, NEQ, 1, is_2, IERR )


! re-write the 3rd order system as a 2nd order system
!---------- ---------- ---------- ---------- ----------
CALL MatCreateNest ( PETSC_COMM_WORLD, 2, [is_1, is_2], 2, [is_1, is_2], [AM_PETSC, PETSC_NULL_OBJECT, PETSC_NULL_OBJECT, I_p_PETSC], AM_b_PETSC, IERR )
CALL MatCreateNest ( PETSC_COMM_WORLD, 2, [is_1, is_2], 2, [is_1, is_2], [AC_PETSC, PETSC_NULL_OBJECT, I_m_PETSC, PETSC_NULL_OBJECT], AC_b_PETSC, IERR )
CALL MatCreateNest ( PETSC_COMM_WORLD, 2, [is_1, is_2], 2, [is_1, is_2], [AK_PETSC, AG_PETSC,  PETSC_NULL_OBJECT, PETSC_NULL_OBJECT], AK_b_PETSC, IERR )

CALL VecCreateNest ( PETSC_COMM_WORLD, 2, [is_1, is_2], [Vdd_PETSC, Udd_PETSC], Xdd_b_PETSC, IERR )
CALL VecCreateNest ( PETSC_COMM_WORLD, 2, [is_1, is_2], [ Vd_PETSC,  Ud_PETSC],  Xd_b_PETSC, IERR )
CALL VecCreateNest ( PETSC_COMM_WORLD, 2, [is_1, is_2], [  V_PETSC,   U_PETSC],   X_b_PETSC, IERR )
CALL VecCreateNest ( PETSC_COMM_WORLD, 2, [is_1, is_2], [  B_PETSC,zero_PETSC],   B_b_PETSC, IERR )


! assembly
!---------- ---------- ---------- ---------- ---------- 
CALL MatAssemblyBegin ( AM_b_PETSC, MAT_FINAL_ASSEMBLY, IERR )
CALL MatAssemblyBegin ( AC_b_PETSC, MAT_FINAL_ASSEMBLY, IERR )
CALL MatAssemblyBegin ( AK_b_PETSC, MAT_FINAL_ASSEMBLY, IERR )
CALL MatAssemblyEnd   ( AM_b_PETSC, MAT_FINAL_ASSEMBLY, IERR )
CALL MatAssemblyEnd   ( AC_b_PETSC, MAT_FINAL_ASSEMBLY, IERR )
CALL MatAssemblyEnd   ( AK_b_PETSC, MAT_FINAL_ASSEMBLY, IERR )

CALL VecAssemblyBegin ( Xdd_b_PETSC, IERR )
CALL VecAssemblyBegin (  Xd_b_PETSC, IERR )
CALL VecAssemblyBegin (   X_b_PETSC, IERR )
CALL VecAssemblyBegin (   B_b_PETSC, IERR )
CALL VecAssemblyBegin (   R_b_PETSC, IERR )
CALL VecAssemblyEnd   ( Xdd_b_PETSC, IERR )
CALL VecAssemblyEnd   (  Xd_b_PETSC, IERR )
CALL VecAssemblyEnd   (   X_b_PETSC, IERR )
CALL VecAssemblyEnd   (   B_b_PETSC, IERR )
CALL VecAssemblyEnd   (   R_b_PETSC, IERR )


! time integration
!---------- ---------- ---------- ---------- ----------


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


! form the effective "stiffness matrix" (not supported)
!---------- ---------- ---------- ---------- ----------
WRITE(*,*) 'Operation not supported in NEWMARK_PETSC_4'
PAUSE

CALL MatAXPY ( AK_b_PETSC, A0, AM_b_PETSC, DIFFERENT_NONZERO_PATTERN, IERR )
CALL MatAXPY ( AK_b_PETSC, A1, AC_b_PETSC, DIFFERENT_NONZERO_PATTERN, IERR )
!CALL MatConvert ( AK_b_PETSC, MATMPIAIJ, MAT_INITIAL_MATRIX, AL_PETSC, IERR )


! SuperLU_dist: matrix factorization
!---------- ---------- ---------- ---------- ----------
CALL CPU_TIME(t0)
WRITE(*,*) ' LU factorization begins ... fasten your seat belts ...'

CALL MatGetFactor        ( AK_b_PETSC, "petsc", MAT_FACTOR_LU, AF_PETSC, IERR )
CALL MatGetSize          ( AK_b_PETSC, m, PETSC_NULL_INTEGER, IERR )
CALL ISCreateStride      ( PETSC_COMM_WORLD, m, 0, 1, perm, IERR )
CALL MatLUFactorSymbolic ( AF_PETSC, AK_b_PETSC, perm, perm, mat_factor_info, IERR )
CALL MatLUFactorNumeric  ( AF_PETSC, AK_b_PETSC, mat_factor_info, IERR )

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
   CALL OUT_NEWMARK_PETSC ( U_PETSC, ID, LTRANS, T, ENERGY )
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
! keep old acceleration
   CALL VecCopy  ( Xdd_b_PETSC, Xdd_OLD_b_PETSC, IERR )
! update Xdd_n+1
   CALL VecScale ( Xdd_b_PETSC, -a3, IERR )
   CALL VecMAXPY ( Xdd_b_PETSC, 4, [a0, -a0, -a2], [HELP_PETSC, X_b_PETSC, Xd_b_PETSC], IERR )
! update Xd_n+1
   CALL VecScale (  Xd_b_PETSC, -a4, IERR )
   CALL VecMAXPY (  Xd_b_PETSC, 4, [a1, -a1, -a5], [HELP_PETSC, X_b_PETSC, Xdd_OLD_b_PETSC], IERR )
! transfer X_n+1
   CALL VecCopy  (  HELP_PETSC, X_b_PETSC, IERR )


END DO
!---------- ---------- ---------- ---------- ----------


! destroy PETSc objects
!---------- ---------- ---------- ---------- ----------
CALL MatDestroy ( AM_b_PETSC, IERR )
CALL MatDestroy ( AC_b_PETSC, IERR )
CALL MatDestroy ( AK_b_PETSC, IERR )

CALL MatDestroy (  I_p_PETSC, IERR )
CALL MatDestroy (  I_m_PETSC, IERR )

CALL VecDestroy (Xdd_b_PETSC, IERR )
CALL VecDestroy ( Xd_b_PETSC, IERR )
CALL VecDestroy (  X_b_PETSC, IERR )
CALL VecDestroy (  B_b_PETSC, IERR )
CALL VecDestroy (  R_b_PETSC, IERR )

CALL VecDestroy (  Vdd_PETSC, IERR )
CALL VecDestroy (   Vd_PETSC, IERR )
CALL VecDestroy (    V_PETSC, IERR )
CALL VecDestroy (  Udd_PETSC, IERR )
CALL VecDestroy (   Ud_PETSC, IERR )
CALL VecDestroy (    U_PETSC, IERR )

CALL VecDestroy (     zero_PETSC, IERR )
CALL VecDestroy (     Help_PETSC, IERR )
CALL VecDestroy (Xdd_OLD_b_PETSC, IERR )


RETURN
END

!************************************************
! gradient update (Mu & Lambda)
!************************************************
! REVISION : W, 15 January 2014 : A life spent making mistakes is not only more honorable, but more useful than a life spent doing nothing.

SUBROUTINE ELEMENT_g_3D_old ( XYZ, INOD, NGP, IEL, E_U, E_P, E_g_Mu, E_g_Lambda )
USE PARAMETERS
IMPLICIT DOUBLE PRECISION (A-H,O-Z)


! in 
!---------- ---------- ---------- ---------- ----------
DIMENSION XYZ ( NJ, NDIM ), INOD ( NNODE, NEL ), NGP ( NEL )
Dimension E_U ( NNODE * NDIM )
Dimension E_P ( NNODE * NDIM )
!---------- ---------- ---------- ---------- ----------


! out
!---------- ---------- ---------- ---------- ----------
Dimension E_g_Mu     ( NNODE )
Dimension E_g_Lambda ( NNODE )
!---------- ---------- ---------- ---------- ----------


! local
!---------- ---------- ---------- ---------- ----------
Dimension E_core_Mu     ( NNODE * NDIM , NNODE * NDIM )
Dimension E_core_Lambda ( NNODE * NDIM , NNODE * NDIM )
!---------- ---------- ---------- ---------- ----------


! local
!---------- ---------- ---------- ---------- ----------
DIMENSION XT(NDIM,NNODE), DJ(NDIM,NDIM), DJI(NDIM,NDIM), FN(NNODE), DFX(NNODE,NDIM), DFXI(NNODE,NDIM)
DIMENSION XINT(9,9), WINT(9,9)

DIMENSION H1(NNODE, NNODE), H2(NNODE, NNODE), H3(NNODE, NNODE), H4(NNODE, NNODE), H5(NNODE, NNODE), H6(NNODE, NNODE), H7(NNODE, NNODE), H8(NNODE, NNODE), H9(NNODE, NNODE)
DIMENSION G1(NNODE, NNODE), G2(NNODE, NNODE), G3(NNODE, NNODE), G4(NNODE, NNODE), G5(NNODE, NNODE), G6(NNODE, NNODE), G7(NNODE, NNODE), G8(NNODE, NNODE), G9(NNODE, NNODE)
!---------- ---------- ---------- ---------- ----------


CALL GAUSS ( XINT, WINT )


! load properties
!---------- ---------- ---------- ---------- ----------
NINT = NGP(IEL)


! define factors and initialize
!---------- ---------- ---------- ---------- ---------- 
E_core_Mu     = 0.0d0
E_core_Lambda = 0.0d0
E_g_Mu        = 0.0d0
E_g_Lambda    = 0.0d0


DO I = 1, NNODE
   K = INOD(I,IEL)
   DO J = 1, NDIM
      XT(J,I) = XYZ(K,J)
   END DO
END DO


! Gauss integration        
!---------- ---------- ---------- ---------- ----------
DO LZ = 1, NINT
   X3 = XINT(LZ,NINT)
   WZ = WINT(LZ,NINT)

   DO LY = 1, NINT
      X2 = XINT(LY,NINT)
      WY = WINT(LY,NINT)

      DO LX = 1, NINT
         X1 = XINT(LX,NINT)
         WX = WINT(LX,NINT)

         WSTAR = WX * WY * WZ


! Jacobian and shape functions    
!---------- ---------- ---------- ---------- ----------         
         SELECT CASE (NNODE)

            CASE (20)                            ! 20-noded elements (serendipity)
               CALL SHAPE_3_20 (X1, X2, X3, FN, DFXI)

            CASE (27)                            ! 27-noded elements (Lagrange)
               CALL SHAPE_3_27 (X1, X2, X3, FN, DFXI)

            CASE DEFAULT
               WRITE(*,'("FLAG ERROR IN EK_3D - NNODE")')
               STOP

         END SELECT


         DJ = MATMUL(XT,DFXI)

         DETJ = DJ(1,1) * DJ(2,2) * DJ(3,3) + DJ(2,1) * DJ(3,2) * DJ(1,3) + DJ(3,1) * DJ(1,2) * DJ(2,3) - &
                DJ(1,3) * DJ(2,2) * DJ(3,1) - DJ(2,3) * DJ(3,2) * DJ(1,1) - DJ(3,3) * DJ(1,2) * DJ(2,1)

         IF (DETJ.LE.0) THEN
            WRITE(*,*)'|J|<0 NEGATIVE DET. JACOBIAN'
            READ(*,*)
            STOP
         END IF

         DJI(1,1) = + DJ(3,3) * DJ(2,2) - DJ(3,2) * DJ(2,3)
         DJI(1,2) = - DJ(3,3) * DJ(1,2) + DJ(3,2) * DJ(1,3)
         DJI(1,3) = + DJ(2,3) * DJ(1,2) - DJ(2,2) * DJ(1,3)
 
         DJI(2,1) = - DJ(3,3) * DJ(2,1) + DJ(3,1) * DJ(2,3)
         DJI(2,2) = + DJ(3,3) * DJ(1,1) - DJ(3,1) * DJ(1,3)
         DJI(2,3) = - DJ(2,3) * DJ(1,1) + DJ(2,1) * DJ(1,3)

         DJI(3,1) = + DJ(3,2) * DJ(2,1) - DJ(3,1) * DJ(2,2)
         DJI(3,2) = - DJ(3,2) * DJ(1,1) + DJ(3,1) * DJ(1,2)
         DJI(3,3) = + DJ(2,2) * DJ(1,1) - DJ(2,1) * DJ(1,2)

         DJI = DJI / DETJ


         DFX = MATMUL(DFXI,DJI)
!---------- ---------- ---------- ---------- ---------- 
         FAC0 = WSTAR * DETJ


! form the matrix
!---------- ---------- ---------- ---------- ----------
         DO I = 1 , NNODE
            DO J = 1 , NNODE
               H1(I,J) = DFX(I,1) * DFX(J,1)
               H2(I,J) = DFX(I,2) * DFX(J,2)
               H3(I,J) = DFX(I,3) * DFX(J,3)
               H4(I,J) = DFX(I,1) * DFX(J,2)
               H5(I,J) = DFX(I,2) * DFX(J,1)
               H6(I,J) = DFX(I,1) * DFX(J,3)
               H7(I,J) = DFX(I,3) * DFX(J,1)
               H8(I,J) = DFX(I,2) * DFX(J,3)
               H9(I,J) = DFX(I,3) * DFX(J,2)
            END DO
         END DO

         G1 = 2.0d0 * H1 + H2 + H3
         G2 = H5
         G3 = H7
         G4 = H4
         G5 = H1 + 2.0d0 * H2 + H3
         G6 = H9
         G7 = H6
         G8 = H8
         G9 = H1 + H2 + 2.0d0 * H3

         DO I = 1 , NNODE
            I1 = I
            I2 = I1 + NNODE
            I3 = I2 + NNODE

            DO J = 1 , NNODE
               J1 = J
               J2 = J1 + NNODE
               J3 = J2 + NNODE

! Mu ----------------------------------------------------------------------------------------------
               E_core_Mu(I1,J1) = G1(I,J) ! 2xx + yy + zz
               E_core_Mu(I1,J2) = G2(I,J) ! yx
               E_core_Mu(I1,J3) = G3(I,J) ! zx

               E_core_Mu(I2,J1) = G4(I,J) ! xy
               E_core_Mu(I2,J2) = G5(I,J) ! xx + 2yy + zz
               E_core_Mu(I2,J3) = G6(I,J) ! zy

               E_core_Mu(I3,J1) = G7(I,J) ! xz          
               E_core_Mu(I3,J2) = G8(I,J) ! yz
               E_core_Mu(I3,J3) = G9(I,J) ! xx + yy + 2zz

! Lambda ------------------------------------------------------------------------------------------
               E_core_Lambda(I1,J1) = H1(I,J) 
               E_core_Lambda(I1,J2) = H4(I,J)
               E_core_Lambda(I1,J3) = H6(I,J)

               E_core_Lambda(I2,J1) = H5(I,J)
               E_core_Lambda(I2,J2) = H2(I,J)
               E_core_Lambda(I2,J3) = H8(I,J)

               E_core_Lambda(I3,J1) = H7(I,J)
               E_core_Lambda(I3,J2) = H9(I,J)
               E_core_Lambda(I3,J3) = H3(I,J)

            END DO
         END DO
!---------- ---------- ---------- ---------- ----------

         E_g_Mu     = E_g_Mu     + FN * Dot_Product ( E_P , MatMul( E_core_Mu     , E_U ) ) * FAC0
         E_g_Lambda = E_g_Lambda + FN * Dot_Product ( E_P , MatMul( E_core_Lambda , E_U ) ) * FAC0

      END DO
   END DO
END DO
!---------- ---------- ---------- ---------- ----------


RETURN
END



!************************************************
! Explicit time stepping based on lagged Central Difference Method: (only 2 system matrices: C, K)
!************************************************
! REVISION : F, 13 June 2014
! An intellectual says a simple thing in a hard way. An artist says a hard thing in a simple way.

SUBROUTINE EXPLICIT_Central_Difference_2 ( AK_PETSC, DIAG_M_PETSC, AC_PETSC, B_PETSC, AK_RD_PETSC, DIAG_M_RD_PETSC, ID, RANK, NNDH, EqDis, U_Store_Mapping, NStore_Mapping, U_Store_Numbers_Global )
USE PARAMETERS
Use Results
USE IFPORT
IMPLICIT DOUBLE PRECISION (A-H,O-Z)

! - Type DECLERATIONS -------------------------------------------------------------------------------------------------------------------------------
Integer, Intent(In),  Dimension ( NNDH * NDim )              :: EqDis
Integer, Intent(In),  Dimension ( NStore_Mapping )           :: U_Store_Numbers_Global
Real(8), Intent(Out), Dimension ( NStore_Mapping , 0:NStep ) :: U_Store_Mapping

Integer,              Dimension ( NNDH * NDim )              :: EqVel
Integer,              Dimension ( NNDH * NDim )              :: EqAcc
!==================================================================================================


!==================================================================================================
! - PETSC INCLUDE FILES
#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"

! - PETSC VARIABLES AND MATRICES
Mat            :: AK_PETSC, AC_PETSC, AK_RD_PETSC
Vec            :: B_PETSC
Vec            :: DIAG_M_PETSC, DIAG_M_RD_PETSC

PetscErrorCode :: IERR
PetscMPIInt    :: SIZE, RANK

PetscScalar    :: FAC

! - PETSc local vectors
Vec            :: U_nm_PETSC, U_nn_PETSC, U_np_PETSC
Vec            :: HELP_1_PETSC, HELP_2_PETSC

PetscScalar, pointer :: U(:), UD(:), UDD(:)
!==================================================================================================

! in 
!---------- ---------- ---------- ---------- ----------
DIMENSION ID ( NJ, NDOF)
!---------- ---------- ---------- ---------- ----------

! local
!---------- ---------- ---------- ---------- ----------
Dimension u_small ( NStore_Mapping )

EqVel = EqDis ! garbage
EqAcc = EqDis
!---------- ---------- ---------- ---------- ----------

!==================================================================================================
! - DEFINE PETSC OBJECTS
CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, U_nm_PETSC, IERR )
CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, U_nn_PETSC, IERR )
CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, U_np_PETSC, IERR )
! for energy computation:
CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, HELP_1_PETSC, IERR )
CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, HELP_2_PETSC, IERR )

CALL VecSet ( U_nm_PETSC, 0.0d0, IERR)
CALL VecSet ( U_nn_PETSC, 0.0d0, IERR)
CALL VecSet ( U_np_PETSC, 0.0d0, IERR)
!==================================================================================================


! inverse of the diagonal mass matrix
!---------- ---------- ---------- ---------- ----------  
CALL VecReciprocal ( DIAG_M_PETSC, IERR )


! form right hand side matrices: inv(M) * A * U(i) + inv(M) * B * U(i-1)
!---------- ---------- ---------- ---------- ----------
CALL MatDiagonalScale ( AK_PETSC, DIAG_M_PETSC, PETSC_NULL_OBJECT, IERR )
CALL MatDiagonalScale ( AC_PETSC, DIAG_M_PETSC, PETSC_NULL_OBJECT, IERR )
CALL VecPointwiseMult (  B_PETSC, DIAG_M_PETSC, B_PETSC, IERR )

CALL MatScale ( AK_PETSC, -DT**2, IERR )
CALL MatAXPY  ( AK_PETSC, -DT, AC_PETSC, DIFFERENT_NONZERO_PATTERN, IERR )
CALL MatScale ( AC_PETSC,     DT, IERR )
Call VecScale (  B_PETSC,  DT**2, IERR )

! print out & Visualization
!---------- ---------- ---------- ---------- ---------- ---------------------------------------------------------------------------------------------
I_PRINT_STEP = 1
Call VecGetOwnershipRange ( U_nn_PETSC, IStart, IEnd, IERR )
!---------- ---------- ---------- ---------- ---------- ---------------------------------------------------------------------------------------------


! Explicit Integration
!---------- ---------- ---------- ---------- ----------
DO ISTEP = 1 , NSTEP

   T = DBLE( ISTEP ) * DT
   IF ( RANK == 0 ) then !.AND. MOD ( ISTEP+1 , I_PRINT_STEP ) == 0 ) THEN
      WRITE (*,'(A30,F10.5)') 'TIME=', T
   END IF


! storage of displacement for the inverse problem - (see test03.F90 for a simple example)
!---------- ---------- ---------- ---------- ---------- ---------------------------------------------------------------------------------------------
! 1- selection: during time integration, access u(big), pick u(small) entries.
   Call VecGetArrayF90 ( U_nn_PETSC, U, IERR )
   Do I = 1, NStore_Mapping
     u_small ( I ) = U ( U_Store_Numbers_Global ( I ) - IStart + 1 )
   End Do
   Call VecRestoreArrayF90 ( U_nn_PETSC, U, IERR )

! 2- store u(small) entries in the dense matrix: We are storing 0:N-1.
   U_Store_Mapping ( : , ISTEP ) = u_small (:)
!---------- ---------- ---------- ---------- ---------- ---------------------------------------------------------------------------------------------


! print out
!---------- ---------- ---------- ---------- ---------- ---------------------------------------------------------------------------------------------
!---------- ---------- ---------- ---------- ---------- ---------------------------------------------------------------------------------------------
   IF ( MOD ( ISTEP , I_PRINT_STEP ) == 0 ) THEN
      CALL TOTAL_ENERGY_RD_EXPLICIT_Central_Difference ( AK_RD_PETSC, DIAG_M_RD_PETSC, U_nn_PETSC, U_nm_PETSC, HELP_1_PETSC, HELP_2_PETSC, ENERGY )

! print out
!---------- ---------- ---------- ---------- ---------- ---------------------------------------------------------------------------------------------
    If ( NNDH /= 0 ) Call VecGetArrayF90 ( U_nn_PETSC,   U,  IErr ) ; ! NNDH
    If ( NNDH /= 0 ) Call VecGetArrayF90 ( U_np_PETSC,  UD,  IErr ) ; ! NNVH ! we are passing the wrong vector
    If ( NNDH /= 0 ) Call VecGetArrayF90 ( U_nm_PETSC, UDD,  IErr ) ; ! NNAH ! we are passing the wrong vector for acceleration. kill this later

    ! Writing down the history of displacement, velocity and acceleration
    Call ResDyn   ( NDIM, NNDH, NNDH, NNDH, IStep, IStart, IEnd, T, Energy, EqDis, EqVel, EqAcc, U, UD, UDD ) ;

    If ( NNDH /= 0 ) Call VecRestoreArrayF90 ( U_nn_PETSC,   U,  IErr ) ; ! NNDH
    If ( NNDH /= 0 ) Call VecRestoreArrayF90 ( U_np_PETSC,  UD,  IErr ) ; ! NNVH
    If ( NNDH /= 0 ) Call VecRestoreArrayF90 ( U_nm_PETSC, UDD,  IErr ) ; ! NNAH


   END IF
!---------- ---------- ---------- ---------- ---------- ---------------------------------------------------------------------------------------------
!---------- ---------- ---------- ---------- ---------- ---------------------------------------------------------------------------------------------


!  external load case ( forward & inverse problem )
!---------- ---------- ---------- ---------- ----------
   CALL TIME_FUNCTION ( T, FAC )
   CALL VecSet ( U_np_PETSC, 0.0d0, IERR )
   CALL VecAXPY( U_np_PETSC, FAC, B_PETSC, IERR )
!---------- ---------- ---------- ---------- ----------


! Explicit time-stepping: central difference for UDD and backward difference for UD
! after this step, we obtain u(n+1) from u(n) and u(n-1)
!---------- ---------- ---------- ---------- ----------
   Call VecMAXPY   ( U_np_PETSC, 2, [ 2.0d0, -1.0d0 ], [ U_nn_PETSC, U_nm_PETSC], IERR )
   CALL MatMultAdd ( AK_PETSC, U_nn_PETSC, U_np_PETSC, U_np_PETSC, IERR )
   CALL MatMultAdd ( AC_PETSC, U_nm_PETSC, U_np_PETSC, U_np_PETSC, IERR )

! update for next iteration
!---------- ---------- ---------- ---------- ----------
   CALL VecCopy ( U_nn_PETSC, U_nm_PETSC, IERR )
   CALL VecCopy ( U_np_PETSC, U_nn_PETSC, IERR )

!---------- ---------- ---------- ---------- ----------
END DO


CALL VecDestroy ( U_nm_PETSC  , IERR )
CALL VecDestroy ( U_nn_PETSC  , IERR )
CALL VecDestroy ( U_np_PETSC  , IERR )
CALL VecDestroy ( HELP_1_PETSC, IERR )
CALL VecDestroy ( HELP_2_PETSC, IERR )


RETURN
END


!************************************************
! Explicit time stepping based on lagged Central Difference Method: (3 system matrices: C, K, G, as in 3D PML)
!************************************************
! REVISION : Th, 19 June 2014
! Important principles may, and must, be inflexible.

SUBROUTINE EXPLICIT_Central_Difference_3 ( AK_PETSC, DIAG_M_PETSC, AC_PETSC, AG_PETSC, B_PETSC, AK_RD_PETSC, DIAG_M_RD_PETSC, ID, RANK, NNDH, EqDis, U_Store_Mapping, NStore_Mapping, U_Store_Numbers_Global )
USE PARAMETERS
Use Results
USE IFPORT
IMPLICIT DOUBLE PRECISION (A-H,O-Z)

! - Type DECLERATIONS -------------------------------------------------------------------------------------------------------------------------------
Integer, Intent(In),  Dimension ( NNDH * NDim )              :: EqDis
Integer, Intent(In),  Dimension ( NStore_Mapping )           :: U_Store_Numbers_Global
Real(8), Intent(Out), Dimension ( NStore_Mapping , 0:NStep ) :: U_Store_Mapping

Integer,              Dimension ( NNDH * NDim )              :: EqVel
Integer,              Dimension ( NNDH * NDim )              :: EqAcc
!==================================================================================================


!==================================================================================================
! - PETSC INCLUDE FILES
#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"

! - PETSC VARIABLES AND MATRICES
Mat            :: AK_PETSC, AC_PETSC, AG_PETSC, AK_RD_PETSC
Vec            :: B_PETSC
Vec            :: DIAG_M_PETSC, DIAG_M_RD_PETSC

PetscErrorCode :: IERR
PetscMPIInt    :: SIZE, RANK

PetscScalar    :: FAC

! - PETSc local vectors
Vec            :: U_nm_PETSC, U_nn_PETSC, U_np_PETSC, U_his_PETSC
Vec            :: HELP_1_PETSC, HELP_2_PETSC

PetscScalar, pointer :: U(:), UD(:), UDD(:)
!==================================================================================================

! in 
!---------- ---------- ---------- ---------- ----------
DIMENSION ID ( NJ, NDOF)
!---------- ---------- ---------- ---------- ----------

! local
!---------- ---------- ---------- ---------- ----------
Dimension u_small ( NStore_Mapping )

EqVel = EqDis ! garbage
EqAcc = EqDis
!---------- ---------- ---------- ---------- ----------

!==================================================================================================
! - DEFINE PETSC OBJECTS
CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, U_his_PETSC, IERR )
CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, U_nm_PETSC, IERR )
CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, U_nn_PETSC, IERR )
CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, U_np_PETSC, IERR )
! for energy computation:
CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, HELP_1_PETSC, IERR )
CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, HELP_2_PETSC, IERR )

CALL VecSet ( U_his_PETSC, 0.0d0, IERR)
CALL VecSet ( U_nm_PETSC, 0.0d0, IERR)
CALL VecSet ( U_nn_PETSC, 0.0d0, IERR)
CALL VecSet ( U_np_PETSC, 0.0d0, IERR)
!==================================================================================================


! inverse of the diagonal mass matrix
!---------- ---------- ---------- ---------- ----------  
CALL VecReciprocal ( DIAG_M_PETSC, IERR )


! form right hand side matrices:
!---------- ---------- ---------- ---------- ----------
CALL MatDiagonalScale ( AK_PETSC, DIAG_M_PETSC, PETSC_NULL_OBJECT, IERR )
CALL MatDiagonalScale ( AC_PETSC, DIAG_M_PETSC, PETSC_NULL_OBJECT, IERR )
CALL MatDiagonalScale ( AG_PETSC, DIAG_M_PETSC, PETSC_NULL_OBJECT, IERR )
CALL VecPointwiseMult (  B_PETSC, DIAG_M_PETSC, B_PETSC, IERR )

CALL MatScale ( AK_PETSC, -DT**2, IERR )
CALL MatAXPY  ( AK_PETSC, -DT   , AC_PETSC, DIFFERENT_NONZERO_PATTERN, IERR )
CALL MatAXPY  ( AK_PETSC, -DT**3, AG_PETSC, DIFFERENT_NONZERO_PATTERN, IERR )
CALL MatScale ( AC_PETSC,     DT, IERR )
CALL MatScale ( AG_PETSC, -DT**2, IERR )
Call VecScale (  B_PETSC,  DT**2, IERR )

! print out & Visualization
!---------- ---------- ---------- ---------- ---------- ---------------------------------------------------------------------------------------------
I_PRINT_STEP = 1
Call VecGetOwnershipRange ( U_nn_PETSC, IStart, IEnd, IERR )
!---------- ---------- ---------- ---------- ---------- ---------------------------------------------------------------------------------------------


! Explicit Integration
!---------- ---------- ---------- ---------- ----------
DO ISTEP = 1 , NSTEP

   T = DBLE( ISTEP ) * DT
   IF ( RANK == 0 ) then !.AND. MOD ( ISTEP+1 , I_PRINT_STEP ) == 0 ) THEN
      WRITE (*,'(A30,F10.5)') 'TIME=', T
   END IF


! storage of displacement for the inverse problem - (see test03.F90 for a simple example)
!---------- ---------- ---------- ---------- ---------- ---------------------------------------------------------------------------------------------
! 1- selection: during time integration, access u(big), pick u(small) entries.
   Call VecGetArrayF90 ( U_nn_PETSC, U, IERR )
   Do I = 1, NStore_Mapping
     u_small ( I ) = U ( U_Store_Numbers_Global ( I ) - IStart + 1 )
   End Do
   Call VecRestoreArrayF90 ( U_nn_PETSC, U, IERR )

! 2- store u(small) entries in the dense matrix: We are storing 0:N-1.
   U_Store_Mapping ( : , ISTEP ) = u_small (:)
!---------- ---------- ---------- ---------- ---------- ---------------------------------------------------------------------------------------------


! print out
!---------- ---------- ---------- ---------- ---------- ---------------------------------------------------------------------------------------------
!---------- ---------- ---------- ---------- ---------- ---------------------------------------------------------------------------------------------
   IF ( MOD ( ISTEP , I_PRINT_STEP ) == 0 ) THEN
      CALL TOTAL_ENERGY_RD_EXPLICIT_Central_Difference ( AK_RD_PETSC, DIAG_M_RD_PETSC, U_nn_PETSC, U_nm_PETSC, HELP_1_PETSC, HELP_2_PETSC, ENERGY )

! print out
!---------- ---------- ---------- ---------- ---------- ---------------------------------------------------------------------------------------------
    If ( NNDH /= 0 ) Call VecGetArrayF90 ( U_nn_PETSC,   U,  IErr ) ; ! NNDH
    If ( NNDH /= 0 ) Call VecGetArrayF90 ( U_np_PETSC,  UD,  IErr ) ; ! NNVH ! we are passing the wrong vector
    If ( NNDH /= 0 ) Call VecGetArrayF90 ( U_nm_PETSC, UDD,  IErr ) ; ! NNAH ! we are passing the wrong vector for acceleration. kill this later

    ! Writing down the history of displacement, velocity and acceleration
    Call ResDyn   ( NDIM, NNDH, NNDH, NNDH, IStep, IStart, IEnd, T, Energy, EqDis, EqVel, EqAcc, U, UD, UDD ) ;

    If ( NNDH /= 0 ) Call VecRestoreArrayF90 ( U_nn_PETSC,   U,  IErr ) ; ! NNDH
    If ( NNDH /= 0 ) Call VecRestoreArrayF90 ( U_np_PETSC,  UD,  IErr ) ; ! NNVH
    If ( NNDH /= 0 ) Call VecRestoreArrayF90 ( U_nm_PETSC, UDD,  IErr ) ; ! NNAH

   END IF
!---------- ---------- ---------- ---------- ---------- ---------------------------------------------------------------------------------------------
!---------- ---------- ---------- ---------- ---------- ---------------------------------------------------------------------------------------------


!  external load case ( forward & inverse problem )
!---------- ---------- ---------- ---------- ----------
   CALL TIME_FUNCTION ( T, FAC )
   CALL VecSet ( U_np_PETSC, 0.0d0, IERR )
   CALL VecAXPY( U_np_PETSC, FAC, B_PETSC, IERR )
!---------- ---------- ---------- ---------- ----------


! Explicit time-stepping: central difference for UDD and backward difference for UD
! after this step, we obtain u(n+1) from u(n) and u(n-1)
!---------- ---------- ---------- ---------- ----------
   Call VecMAXPY   ( U_np_PETSC, 2, [ 2.0d0, -1.0d0 ], [ U_nn_PETSC, U_nm_PETSC], IERR )
   CALL MatMultAdd ( AK_PETSC, U_nn_PETSC, U_np_PETSC, U_np_PETSC, IERR )
   CALL MatMultAdd ( AC_PETSC, U_nm_PETSC, U_np_PETSC, U_np_PETSC, IERR )
   CALL MatMultAdd ( AG_PETSC, U_his_PETSC,U_np_PETSC, U_np_PETSC, IERR )


! update for next iteration
!---------- ---------- ---------- ---------- ----------
   CALL VecAXPY ( U_his_PETSC, DT, U_np_PETSC, IERR )
   CALL VecCopy ( U_nn_PETSC, U_nm_PETSC, IERR )
   CALL VecCopy ( U_np_PETSC, U_nn_PETSC, IERR )

!---------- ---------- ---------- ---------- ----------
END DO


CALL VecDestroy ( U_his_PETSC  , IERR )
CALL VecDestroy ( U_nm_PETSC  , IERR )
CALL VecDestroy ( U_nn_PETSC  , IERR )
CALL VecDestroy ( U_np_PETSC  , IERR )
CALL VecDestroy ( HELP_1_PETSC, IERR )
CALL VecDestroy ( HELP_2_PETSC, IERR )


RETURN
END

