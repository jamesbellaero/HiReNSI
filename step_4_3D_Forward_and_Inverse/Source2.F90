! This file is aimed to transform the 3rd order PML probrem into a 2nd order system, and use a classical Newmark method to march in time.
! A better way of doing this is to utilize NEST matrices (NEWMARK_PETSC_4) however, PETSc does not support NEST matrix factorization now.
! Hence, we have this mess here.

!************************************************
! main subroutine that calls other subroutines 
!************************************************
! REVISION : Th 2 Aug 2012

SUBROUTINE Solve_PML_3D_2nd_order_ODE ( XYZ, INOD, NGP, MTEL, PMAT, PML_PARAM, ID, ID_BC, LTRANS, SIZE, RANK )
USE PARAMETERS
IMPLICIT DOUBLE PRECISION (A-H,O-Z)

!==================================================================================================
! - PETSC INCLUDE FILES
#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscmat.h"

! - PETSC VARIABLES AND MATRICES
Mat            ::     AM_b_PETSC, AC_b_PETSC, AK_b_PETSC
Vec            :: DIAG_M_b_PETSC
Vec            ::      B_b_PETSC

PetscErrorCode :: IERR
PetscMPIInt    :: SIZE, RANK
!==================================================================================================


! in 
!---------- ---------- ---------- ---------- ----------
DIMENSION XYZ(NJ,NDIM), INOD(NNODE,NEL), NGP(NEL), MTEL(NEL), PMAT(NMAT,NPM), PML_PARAM(NDIM * 2 , 4), ID(NJ,NDOF), ID_BC(NEL,NDIM**2), LTRANS(NTRANS)
!---------- ---------- ---------- ---------- ----------


! Define PETSc onjects
!---------- ---------- ---------- ---------- ----------
CALL MatCreateAIJ ( PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, 2*NEQ, 2*NEQ, 1000, PETSC_NULL_INTEGER,  10, PETSC_NULL_INTEGER, AM_b_PETSC, IERR )
CALL MatCreateAIJ ( PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, 2*NEQ, 2*NEQ, 1000, PETSC_NULL_INTEGER,  10, PETSC_NULL_INTEGER, AC_b_PETSC, IERR )
CALL MatCreateAIJ ( PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, 2*NEQ, 2*NEQ, 3000, PETSC_NULL_INTEGER,  10, PETSC_NULL_INTEGER, AK_b_PETSC, IERR )

CALL MatSetOption    ( AM_b_PETSC, MAT_IGNORE_ZERO_ENTRIES, PETSC_TRUE, IERR )
CALL MatSetOption    ( AC_b_PETSC, MAT_IGNORE_ZERO_ENTRIES, PETSC_TRUE, IERR )
CALL MatSetOption    ( AK_b_PETSC, MAT_IGNORE_ZERO_ENTRIES, PETSC_TRUE, IERR )

CALL MatSetOption    ( AM_b_PETSC, MAT_ROW_ORIENTED, PETSC_FALSE, IERR )
CALL MatSetOption    ( AC_b_PETSC, MAT_ROW_ORIENTED, PETSC_FALSE, IERR )
CALL MatSetOption    ( AK_b_PETSC, MAT_ROW_ORIENTED, PETSC_FALSE, IERR )

CALL VecCreateMPI    ( PETSC_COMM_WORLD, PETSC_DECIDE, 2*NEQ, DIAG_M_b_PETSC, IERR )
CALL VecCreateMPI    ( PETSC_COMM_WORLD, PETSC_DECIDE, 2*NEQ,      B_b_PETSC, IERR )

CALL VecSetOption    ( DIAG_M_b_PETSC, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE, IERR )
CALL VecSetOption    (      B_b_PETSC, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE, IERR )

CALL VecSet          ( DIAG_M_b_PETSC, 0.0d0, IERR )
CALL VecSet          (      B_b_PETSC, 0.0d0, IERR )
!---------- ---------- ---------- ---------- ----------


! assemble system matrices
!---------- ---------- ---------- ---------- ----------
! regular domain
CALL ASSEM_REGULAR_DOMAIN_PETSC_3D_2nd_order_ODE  ( AM_b_PETSC, AK_b_PETSC, DIAG_M_b_PETSC, XYZ, INOD, NGP, MTEL, PMAT, PML_PARAM, ID_BC, ID, SIZE, RANK )
! force
CALL ASSEM_FORCE_PETSC_3D                                      ( B_b_PETSC,                 XYZ, INOD, NGP, MTEL, PMAT,            ID_BC, ID, SIZE, RANK )
! PML
CALL ASSEM_PML_PETSC_3D_2nd_order_ODE ( AM_b_PETSC, AC_b_PETSC, AK_b_PETSC, DIAG_M_b_PETSC, XYZ, INOD, NGP, MTEL, PMAT, PML_PARAM, ID_BC, ID, SIZE, RANK )
!---------- ---------- ---------- ---------- ----------


! PETSc assembly
!---------- ---------- ---------- ---------- ----------
CALL MatAssemblyBegin ( AK_b_PETSC, MAT_FINAL_ASSEMBLY, IERR )
CALL MatAssemblyBegin ( AC_b_PETSC, MAT_FINAL_ASSEMBLY, IERR )
CALL MatAssemblyBegin ( AM_b_PETSC, MAT_FINAL_ASSEMBLY, IERR )
CALL MatAssemblyEnd   ( AM_b_PETSC, MAT_FINAL_ASSEMBLY, IERR )
CALL MatAssemblyEnd   ( AC_b_PETSC, MAT_FINAL_ASSEMBLY, IERR )
CALL MatAssemblyEnd   ( AK_b_PETSC, MAT_FINAL_ASSEMBLY, IERR )

CALL VecAssemblyBegin ( DIAG_M_b_PETSC, IERR )
CALL VecAssemblyBegin (      B_b_PETSC, IERR )
CALL VecAssemblyEnd   ( DIAG_M_b_PETSC, IERR )
CALL VecAssemblyEnd   (      B_b_PETSC, IERR )
!---------- ---------- ---------- ---------- ----------


! Solve the 2nd order system of ODEs
!---------- ---------- ---------- ---------- ----------
!  CALL      NEWMARK_PETSC_2nd_order_ODE (     AM_b_PETSC, AC_b_PETSC, AK_b_PETSC, B_b_PETSC, ID, LTRANS, RANK )
!   CALL EXPLICIT_RK4_PETSC_2nd_order_ODE ( DIAG_M_b_PETSC, AC_b_PETSC, AK_b_PETSC, B_b_PETSC, ID, LTRANS, RANK )
write(*,*)'un-comment source2 line 91'
!---------- ---------- ---------- ---------- ----------


! Destroy PETSc objects
!---------- ---------- ---------- ---------- ----------
CALL MatDestroy (     AM_b_PETSC, IERR )
CALL MatDestroy (     AC_b_PETSC, IERR )
CALL MatDestroy (     AK_b_PETSC, IERR )
CALL VecDestroy ( DIAG_M_b_PETSC, IERR )
CALL VecDestroy (      B_b_PETSC, IERR )
!---------- ---------- ---------- ---------- ----------


RETURN
END


!************************************************
! ASSEMBLE REGULAR DOMAIN MASS AND STIFFNESS MATRIX
!************************************************
! REVISION : Th 2 Aug 2012

SUBROUTINE ASSEM_REGULAR_DOMAIN_PETSC_3D_2nd_order_ODE ( AM_b_PETSC, AK_b_PETSC, DIAG_M_b_PETSC, XYZ, INOD, NGP, MTEL, PMAT, PML_PARAM, ID_BC, ID, SIZE, RANK )
USE PARAMETERS
IMPLICIT DOUBLE PRECISION (A-H,O-Z)

!==================================================================================================
! - PETSC INCLUDE FILES
#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscmat.h"

! - PETSC VARIABLES AND MATRICES
Mat            ::     AM_b_PETSC, AK_b_PETSC
Vec            :: DIAG_M_b_PETSC

PetscErrorCode :: IERR
PetscMPIInt    :: SIZE, RANK

PetscScalar    :: EM_PETSC(NNODE * NDIM, NNODE * NDIM), EK_PETSC(NNODE * NDIM, NNODE * NDIM)
PetscScalar    :: EDIAG_M_PETSC(NNODE * NDIM)
!==================================================================================================


! in 
!---------- ---------- ---------- ---------- ---------- 
DIMENSION XYZ(NJ,NDIM), INOD(NNODE,NEL), NGP(NEL), MTEL(NEL), PMAT(NMAT,NPM), PML_PARAM(NDIM * 2 , 4), ID(NJ,NDOF), ID_BC(NEL,NDIM**2)
!---------- ---------- ---------- ---------- ----------


! local
!---------- ---------- ---------- ---------- ---------- 
DIMENSION EM(NNODE * NDIM, NNODE * NDIM), EK(NNODE * NDIM, NNODE * NDIM), ND1(NNODE * NDIM)
DIMENSION EDIAG_M(NNODE * NDIM)
!---------- ---------- ---------- ---------- ----------


DO IEL = 1 + RANK , NEL , SIZE

   MTYPE = MTEL(IEL)
!  MTYPE 1 : REGULAR DOMAIN
!  MTYPE 2 : PML
   IF ( MTYPE /= 1 ) CYCLE


   WRITE(*,'(A30, I10)')'REGULAR DOMAIN - ELEMENT NO.', IEL


! regular domain EM, EK
! flag 201 : regular domain mass matrix
!---------- ---------- ---------- ---------- ----------           
   CALL ELEMENT_M_3D ( XYZ, INOD, NGP, MTEL, PMAT, PML_PARAM, IEL, EM, 201 )
   CALL ELEMENT_K_3D ( XYZ, INOD, NGP, MTEL, PMAT,            IEL, EK )


! Diagonal mass matrix for Explicit time-stepping
!---------- ---------- ---------- ---------- ----------
   CALL Extract_Diagonal ( EM, NNODE * NDIM, EDIAG_M )
!  CALL HRZ_Mass_Diagonalization ( EM, NDIM, EDIAG_M )


! Send Fortran matrix to PETSc data structure
!---------- ---------- ---------- ---------- ---------- 
   EK_PETSC      = EK
   EM_PETSC      = EM

   EDIAG_M_PETSC = EDIAG_M


! Local to Global mapping
!---------- ---------- ---------- ---------- ----------   
   DO I = 1 , NNODE
      I1 = I
      I2 = I1 + NNODE
      I3 = I2 + NNODE

      ND1(I1) = ID(INOD(I,IEL) , 1)
      ND1(I2) = ID(INOD(I,IEL) , 2)
      ND1(I3) = ID(INOD(I,IEL) , 3)
   END DO
!---------- ---------- ---------- ---------- ----------


! Modify ND according to C indexing
!---------- ---------- ---------- ---------- ----------
   DO I = 1 , NNODE * NDIM
      ND1(I) = ND1(I) - 1
   END DO


! assemble mass and stiffness matrices (MAT_ROW_ORIENTED)
!---------- ---------- ---------- ---------- ----------   
   NEQEL = NNODE * NDIM
   CALL MatSetValues ( AM_b_PETSC, NEQEL, ND1, NEQEL, ND1, EM_PETSC, ADD_VALUES, IERR )
   CALL MatSetValues ( AK_b_PETSC, NEQEL, ND1, NEQEL, ND1, EK_PETSC, ADD_VALUES, IERR )
!---------- ---------- ---------- ---------- ----------


! assemble Diagonal mass for Explicit time stepping
!---------- ---------- ---------- ---------- ----------
   CALL VecSetValues ( DIAG_M_b_PETSC, NEQEL, ND1, EDIAG_M_PETSC, ADD_VALUES, IERR )
!---------- ---------- ---------- ---------- ----------


END DO
!---------- ---------- ---------- ---------- ----------


RETURN
END


!************************************************
! ASSEMBLE PML elements
!************************************************
! REVISION : F 3 Aug 2012

SUBROUTINE ASSEM_PML_PETSC_3D_2nd_order_ODE ( AM_b_PETSC, AC_b_PETSC, AK_b_PETSC, DIAG_M_b_PETSC, XYZ, INOD, NGP, MTEL, PMAT, PML_PARAM, ID_BC, ID, SIZE, RANK)
USE PARAMETERS
IMPLICIT DOUBLE PRECISION (A-H,O-Z)

!==================================================================================================
! - PETSC INCLUDE FILES
#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscmat.h"

! - PETSC VARIABLES AND MATRICES
Mat            :: AM_b_PETSC, AC_b_PETSC, AK_b_PETSC
Vec            :: DIAG_M_b_PETSC

PetscErrorCode :: IERR
PetscMPIInt    :: SIZE, RANK

PetscScalar    :: EM_a_PETSC(NNODE * NDIM, NNODE * NDIM), EM_b_PETSC(NNODE * NDIM, NNODE * NDIM), EM_c_PETSC(NNODE * NDIM, NNODE * NDIM), &
                  EM_d_PETSC(NNODE * NDIM, NNODE * NDIM)

PetscScalar    :: EA_e_PETSC(NNODE * 3, NNODE * 6), EA_p_PETSC(NNODE * 3, NNODE * 6), EA_w_PETSC(NNODE * 3, NNODE * 6)
PetscScalar    :: EN_a_PETSC(NNODE * 6, NNODE * 6), EN_b_PETSC(NNODE * 6, NNODE * 6), EN_c_PETSC(NNODE * 6, NNODE * 6), EN_d_PETSC(NNODE * 6, NNODE * 6)

PetscScalar    :: EA_e_PETSC_T(NNODE * 6, NNODE * 3), EA_p_PETSC_T(NNODE * 6, NNODE * 3), EA_w_PETSC_T(NNODE * 6, NNODE * 3)

PetscScalar    :: EM_a_DIAG_PETSC(NNODE * NDIM), EN_a_DIAG_PETSC(NNODE * 6)
!==================================================================================================


! in 
!---------- ---------- ---------- ---------- ---------- 
DIMENSION XYZ(NJ,NDIM), INOD(NNODE,NEL), NGP(NEL), MTEL(NEL), PMAT(NMAT,NPM), PML_PARAM(NDIM * 2 , 4), ID(NJ,NDOF), ID_BC(NEL,NDIM**2)
!---------- ---------- ---------- ---------- ----------


! local
!---------- ---------- ---------- ---------- ----------
DIMENSION EM_a(NNODE * NDIM, NNODE * NDIM), EM_b(NNODE * NDIM, NNODE * NDIM), EM_c(NNODE * NDIM, NNODE * NDIM), EM_d(NNODE * NDIM, NNODE * NDIM)
DIMENSION EA_e_upper(NNODE * 3, NNODE * 6), EA_p_upper(NNODE * 3, NNODE * 6), EA_w_upper(NNODE * 3, NNODE * 6)
DIMENSION EA_e_lower(NNODE * 3, NNODE * 6), EA_p_lower(NNODE * 3, NNODE * 6), EA_w_lower(NNODE * 3, NNODE * 6)
DIMENSION EN_a(NNODE * 6, NNODE * 6), EN_b(NNODE * 6, NNODE * 6), EN_c(NNODE * 6, NNODE * 6), EN_d(NNODE * 6, NNODE * 6)

DIMENSION ND1(NNODE * NDIM), ND2(NNODE * 6), ND1_b(NNODE * NDIM), ND2_b(NNODE * 6)

DIMENSION EM_a_DIAG(NNODE * NDIM), EN_a_DIAG(NNODE * 6)
DIMENSION EM_x_DIAG(NNODE * NDIM), EN_x_DIAG(NNODE * 6)

DIMENSION I_PML(NEQ)
!---------- ---------- ---------- ---------- ----------


DO IEL  = 1 + RANK , NEL , SIZE

   MTYPE = MTEL(IEL)
!  MTYPE 1 : REGULAR DOMAIN
!  MTYPE 2 : PML
   IF ( MTYPE /= 2 ) CYCLE


   WRITE(*,'(A30, I10)')'PML - ELEMENT NO.', IEL


! PML matrices
!---------- ---------- ---------- ---------- ----------   
   CALL ELEMENT_M_3D  ( XYZ, INOD, NGP, MTEL, PMAT, PML_PARAM, IEL, EM_a,       202 )    ! PML M_a
   CALL ELEMENT_M_3D  ( XYZ, INOD, NGP, MTEL, PMAT, PML_PARAM, IEL, EM_b,       203 )    ! PML M_b
   CALL ELEMENT_M_3D  ( XYZ, INOD, NGP, MTEL, PMAT, PML_PARAM, IEL, EM_c,       204 )    ! PML M_c
   CALL ELEMENT_M_3D  ( XYZ, INOD, NGP, MTEL, PMAT, PML_PARAM, IEL, EM_d,       205 )    ! PML M_d


!---------- ---------- ---------- ---------- ----------
I_Constitutive = 0

   SELECT CASE (I_Constitutive)


      CASE (0)
! arash's style constitutive
!---------- ---------- ---------- ---------- ----------
         CALL ELEMENT_EA_3D ( XYZ, INOD, NGP, MTEL, PMAT, PML_PARAM, IEL, EA_e_upper, 101 )    ! PML EA_e_upper
         CALL ELEMENT_EA_3D ( XYZ, INOD, NGP, MTEL, PMAT, PML_PARAM, IEL, EA_p_upper, 102 )    ! PML EA_p_upper
         CALL ELEMENT_EA_3D ( XYZ, INOD, NGP, MTEL, PMAT, PML_PARAM, IEL, EA_w_upper, 103 )    ! PML EA_w_upper
         CALL ELEMENT_EA_3D ( XYZ, INOD, NGP, MTEL, PMAT, PML_PARAM, IEL, EA_e_lower, 104 )    ! PML EA_e_lower
         CALL ELEMENT_EA_3D ( XYZ, INOD, NGP, MTEL, PMAT, PML_PARAM, IEL, EA_p_lower, 105 )    ! PML EA_p_lower
         CALL ELEMENT_EA_3D ( XYZ, INOD, NGP, MTEL, PMAT, PML_PARAM, IEL, EA_w_lower, 106 )    ! PML EA_w_lower

         CALL ELEMENT_EN_3D ( XYZ, INOD, NGP, MTEL, PMAT, PML_PARAM, IEL, EN_a,       301 )    ! PML EN_a
         CALL ELEMENT_EN_3D ( XYZ, INOD, NGP, MTEL, PMAT, PML_PARAM, IEL, EN_b,       302 )    ! PML EN_b
         CALL ELEMENT_EN_3D ( XYZ, INOD, NGP, MTEL, PMAT, PML_PARAM, IEL, EN_c,       303 )    ! PML EN_c
         CALL ELEMENT_EN_3D ( XYZ, INOD, NGP, MTEL, PMAT, PML_PARAM, IEL, EN_d,       304 )    ! PML EN_d


      CASE (1)
! sezgin's style constitutive
!---------- ---------- ---------- ---------- ----------
         CALL ELEMENT_EA_3D_Sezgin ( XYZ, INOD, NGP, MTEL, PMAT, PML_PARAM, IEL, EA_e_upper, 101 )    ! PML EA_e_upper
         CALL ELEMENT_EA_3D_Sezgin ( XYZ, INOD, NGP, MTEL, PMAT, PML_PARAM, IEL, EA_p_upper, 102 )    ! PML EA_p_upper
         CALL ELEMENT_EA_3D_Sezgin ( XYZ, INOD, NGP, MTEL, PMAT, PML_PARAM, IEL, EA_w_upper, 103 )    ! PML EA_w_upper
         CALL ELEMENT_EA_3D_Sezgin ( XYZ, INOD, NGP, MTEL, PMAT, PML_PARAM, IEL, EA_e_lower, 104 )    ! PML EA_e_lower
         CALL ELEMENT_EA_3D_Sezgin ( XYZ, INOD, NGP, MTEL, PMAT, PML_PARAM, IEL, EA_p_lower, 105 )    ! PML EA_p_lower
         CALL ELEMENT_EA_3D_Sezgin ( XYZ, INOD, NGP, MTEL, PMAT, PML_PARAM, IEL, EA_w_lower, 106 )    ! PML EA_w_lower

         CALL ELEMENT_EN_3D_Sezgin ( XYZ, INOD, NGP, MTEL, PMAT, PML_PARAM, IEL, EN_a,       301 )    ! PML EN_a
         CALL ELEMENT_EN_3D_Sezgin ( XYZ, INOD, NGP, MTEL, PMAT, PML_PARAM, IEL, EN_b,       302 )    ! PML EN_b
         CALL ELEMENT_EN_3D_Sezgin ( XYZ, INOD, NGP, MTEL, PMAT, PML_PARAM, IEL, EN_c,       303 )    ! PML EN_c
         CALL ELEMENT_EN_3D_Sezgin ( XYZ, INOD, NGP, MTEL, PMAT, PML_PARAM, IEL, EN_d,       304 )    ! PML EN_d


   END SELECT


! Diagonal mass matrix for Explicit time-stepping
!---------- ---------- ---------- ---------- ----------
! 1. EM matrices
   CALL Extract_Diagonal ( EM_a, NNODE * NDIM, EM_a_DIAG )
!  CALL HRZ_Mass_Diagonalization ( EM_a, NDIM, EM_a_DIAG )
!  CALL HRZ_Mass_Diagonalization ( EM_b, NDIM, EM_x_DIAG )
!  CALL HRZ_Mass_Diagonalization ( EM_c, NDIM, EM_x_DIAG )
!  CALL HRZ_Mass_Diagonalization ( EM_d, NDIM, EM_x_DIAG )
! 2. EN matrices
   CALL Extract_Diagonal ( EN_a, NNODE *  6  , EN_a_DIAG )
!  CALL HRZ_Mass_Diagonalization ( EN_a,  6  , EN_a_DIAG )
!  CALL HRZ_Mass_Diagonalization ( EN_b,  6  , EN_x_DIAG )
!  CALL HRZ_Mass_Diagonalization ( EN_c,  6  , EN_x_DIAG )
!  CALL HRZ_Mass_Diagonalization ( EN_d,  6  , EN_x_DIAG )


! Send Fortran matrix to PETSc data structure
!---------- ---------- ---------- ---------- ---------- 
   EM_a_PETSC   = EM_a
   EM_b_PETSC   = EM_b
   EM_c_PETSC   = EM_c
   EM_d_PETSC   = EM_d

   EA_e_PETSC   = EA_e_upper
   EA_p_PETSC   = EA_p_upper
   EA_w_PETSC   = EA_w_upper

   EA_e_PETSC_T = transpose(EA_e_lower)
   EA_p_PETSC_T = transpose(EA_p_lower)
   EA_w_PETSC_T = transpose(EA_w_lower)

   EN_a_PETSC   = EN_a
   EN_b_PETSC   = EN_b
   EN_c_PETSC   = EN_c
   EN_d_PETSC   = EN_d

   EM_a_DIAG_PETSC = EM_a_DIAG
   EN_a_DIAG_PETSC = EN_a_DIAG


! Local to Global mapping
!---------- ---------- ---------- ---------- ----------   
   DO I = 1 , NNODE

      I1 = I
      I2 = I1 + NNODE
      I3 = I2 + NNODE
      I4 = I3 + NNODE
      I5 = I4 + NNODE
      I6 = I5 + NNODE

      ND1(I1) = ID(INOD(I,IEL) , 1)
      ND1(I2) = ID(INOD(I,IEL) , 2)
      ND1(I3) = ID(INOD(I,IEL) , 3)

      ND2(I1) = ID(INOD(I,IEL) , 4)
      ND2(I2) = ID(INOD(I,IEL) , 5)
      ND2(I3) = ID(INOD(I,IEL) , 6)
      ND2(I4) = ID(INOD(I,IEL) , 7)
      ND2(I5) = ID(INOD(I,IEL) , 8)
      ND2(I6) = ID(INOD(I,IEL) , 9)

   END DO
!---------- ---------- ---------- ---------- ----------


! Modify ND according to C indexing
!---------- ---------- ---------- ---------- ----------
   DO I = 1 , NNODE * NDIM
      ND1(I) = ND1(I) - 1
   END DO

   DO I = 1 , NNODE * 6
      ND2(I) = ND2(I) - 1
   END DO


! Local to Global mapping for the AG block
!---------- ---------- ---------- ---------- ----------
   ND1_b = ND1 + NEQ
   ND2_b = ND2 + NEQ


! assemble EM matrices (MAT_ROW_ORIENTED) 
!---------- ---------- ---------- ---------- ----------   
   NEQEL_U = NNODE * NDIM
   CALL MatSetValues ( AM_b_PETSC, NEQEL_U, ND1, NEQEL_U, ND1,   EM_a_PETSC, ADD_VALUES, IERR )
   CALL MatSetValues ( AC_b_PETSC, NEQEL_U, ND1, NEQEL_U, ND1,   EM_b_PETSC, ADD_VALUES, IERR )
   CALL MatSetValues ( AK_b_PETSC, NEQEL_U, ND1, NEQEL_U, ND1,   EM_c_PETSC, ADD_VALUES, IERR )
   CALL MatSetValues ( AK_b_PETSC, NEQEL_U, ND1, NEQEL_U, ND1_b, EM_d_PETSC, ADD_VALUES, IERR )
!---------- ---------- ---------- ---------- ----------


! assemble EN matrices (MAT_ROW_ORIENTED)
!---------- ---------- ---------- ---------- ----------   
   NEQEL_S = NNODE * 6
   CALL MatSetValues ( AM_b_PETSC, NEQEL_S, ND2, NEQEL_S, ND2,   EN_a_PETSC, ADD_VALUES, IERR )
   CALL MatSetValues ( AC_b_PETSC, NEQEL_S, ND2, NEQEL_S, ND2,   EN_b_PETSC, ADD_VALUES, IERR )
   CALL MatSetValues ( AK_b_PETSC, NEQEL_S, ND2, NEQEL_S, ND2,   EN_c_PETSC, ADD_VALUES, IERR )
   CALL MatSetValues ( AK_b_PETSC, NEQEL_S, ND2, NEQEL_S, ND2_b, EN_d_PETSC, ADD_VALUES, IERR )
!---------- ---------- ---------- ---------- ----------


! assemble EA matrices (MAT_ROW_ORIENTED)
!---------- ---------- ---------- ---------- ----------   
   NEQEL_U = NNODE * NDIM
   NEQEL_S = NNODE * 6
   CALL MatSetValues ( AC_b_PETSC, NEQEL_u, ND1, NEQEL_s, ND2,   EA_e_PETSC, ADD_VALUES, IERR )
   CALL MatSetValues ( AK_b_PETSC, NEQEL_u, ND1, NEQEL_s, ND2,   EA_p_PETSC, ADD_VALUES, IERR )
   CALL MatSetValues ( AK_b_PETSC, NEQEL_u, ND1, NEQEL_s, ND2_b, EA_w_PETSC, ADD_VALUES, IERR )
! assemble transpose
   CALL MatSetValues ( AC_b_PETSC, NEQEL_s, ND2, NEQEL_u, ND1,   -EA_e_PETSC_T, ADD_VALUES, IERR )
   CALL MatSetValues ( AK_b_PETSC, NEQEL_s, ND2, NEQEL_u, ND1,   -EA_p_PETSC_T, ADD_VALUES, IERR )
   CALL MatSetValues ( AK_b_PETSC, NEQEL_s, ND2, NEQEL_u, ND1_b, -EA_w_PETSC_T, ADD_VALUES, IERR )
!---------- ---------- ---------- ---------- ----------


! Assemble Diagonal mass for Explicit time stepping
!---------- ---------- ---------- ---------- ----------
   CALL VecSetValues ( DIAG_M_b_PETSC, NEQEL_U, ND1, EM_a_DIAG_PETSC, ADD_VALUES, IERR )
   CALL VecSetValues ( DIAG_M_b_PETSC, NEQEL_S, ND2, EN_a_DIAG_PETSC, ADD_VALUES, IERR )
!---------- ---------- ---------- ---------- ----------


END DO
!---------- ---------- ---------- ---------- ----------
 

! assemble identity blocks (two possibilities)
!---------- ---------- ---------- ---------- ----------
CALL PML_Identity_Operator ( I_PML, INOD, MTEL, ID )


!---------- ---------- ---------- ---------- ----------
I_Identity_Operator = 0


SELECT CASE (I_Identity_Operator)


   CASE (0)
!---------- ---------- ---------- ---------- ----------
! in M and C
      DO I = 1, NEQ

         IF ( I_PML(I) == 1 ) THEN 
            CALL MatSetValue (     AM_b_PETSC, NEQ+I-1, NEQ+I-1, +1.0d0, ADD_VALUES, IERR )
            CALL MatSetValue (     AC_b_PETSC, NEQ+I-1,     I-1, -1.0d0, ADD_VALUES, IERR )
            CALL VecSetValue ( DIAG_M_b_PETSC, NEQ+I-1,          +1.0d0, ADD_VALUES, IERR )
         END IF

      END DO


   CASE (1)
!---------- ---------- ---------- ---------- ----------
! in C and K
      DO I = 1, NEQ

         IF ( I_PML(I) == 1 ) THEN
            CALL MatSetValue ( AC_b_PETSC, NEQ+I-1, NEQ+I-1, +1.0d0, ADD_VALUES, IERR )
            CALL MatSetValue ( AK_b_PETSC, NEQ+I-1,     I-1, -1.0d0, ADD_VALUES, IERR )
         END IF

      END DO
!---------- ---------- ---------- ---------- ----------


END SELECT


RETURN
END


!************************************************
! PML Identity operator (only PML dofs are active) - sequential
!************************************************
! REVISION : F 10 August 2012

SUBROUTINE PML_Identity_Operator ( I_PML, INOD, MTEL, ID )
USE PARAMETERS
IMPLICIT DOUBLE PRECISION (A-H,O-Z)


! in 
!---------- ---------- ---------- ---------- ----------
DIMENSION I_PML(NEQ)
DIMENSION INOD(NNODE,NEL), MTEL(NEL), ID(NJ,NDOF)
!---------- ---------- ---------- ---------- ----------


! local
!---------- ---------- ---------- ---------- ----------
DIMENSION ND1(NNODE * NDIM), ND2(NNODE * 6)
!---------- ---------- ---------- ---------- ----------


! initialize
!---------- ---------- ---------- ---------- ---------- 
I_PML = 0
!---------- ---------- ---------- ---------- ----------


DO IEL = 1, NEL

   MTYPE = MTEL(IEL)
!  MTYPE 1 : REGULAR DOMAIN
!  MTYPE 2 : PML
   IF ( MTYPE /= 2 ) CYCLE


! Local to Global mapping
!---------- ---------- ---------- ---------- ----------   
   DO I = 1 , NNODE

      I1 = I
      I2 = I1 + NNODE
      I3 = I2 + NNODE
      I4 = I3 + NNODE
      I5 = I4 + NNODE
      I6 = I5 + NNODE

      ND1(I1) = ID(INOD(I,IEL) , 1)
      ND1(I2) = ID(INOD(I,IEL) , 2)
      ND1(I3) = ID(INOD(I,IEL) , 3)

      ND2(I1) = ID(INOD(I,IEL) , 4)
      ND2(I2) = ID(INOD(I,IEL) , 5)
      ND2(I3) = ID(INOD(I,IEL) , 6)
      ND2(I4) = ID(INOD(I,IEL) , 7)
      ND2(I5) = ID(INOD(I,IEL) , 8)
      ND2(I6) = ID(INOD(I,IEL) , 9)

   END DO
!---------- ---------- ---------- ---------- ----------


! construct the PML identity operator
!---------- ---------- ---------- ---------- ----------
! 1. deformation equations
   DO I = 1, NNODE * NDIM
      IF ( ND1(I) > 0 ) I_PML( ND1(I) ) = 1
   END DO

! 2. stress equations
   DO I = 1, NNODE * 6
      IF ( ND2(I) > 0 ) I_PML( ND2(I) ) = 1
   END DO
!---------- ---------- ---------- ---------- ----------


END DO


RETURN
END


!************************************************
! 3D-PML RK4: explicit time-stepping (for 3rd order ODEs)
! We re-write the 3rd order system as a 2nd order system and apply the standard Newmark method.
! Should be called from SUBROUTINE Solve_PML_3D_2nd_order_ODE
!************************************************
! REVISION : Th 2 Aug 2012

SUBROUTINE EXPLICIT_RK4_PETSC_2nd_order_ODE ( DIAG_M_b_PETSC, AC_b_PETSC, AK_b_PETSC, B_b_PETSC, ID, LTRANS, RANK )
USE PARAMETERS
IMPLICIT DOUBLE PRECISION (A-H,O-Z)

!==================================================================================================
! - PETSC INCLUDE FILES
#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscmat.h"

! - PETSC VARIABLES AND MATRICES
Mat            :: AC_b_PETSC, AK_b_PETSC
Vec            ::  B_b_PETSC, DIAG_M_b_PETSC, U_PETSC

IS             :: is_1, is_2

PetscErrorCode :: IERR
PetscMPIInt    :: SIZE, RANK

PetscScalar    :: fac7, fac8, fac9

! - PETSc local vectors
Vec            :: x1_PETSC, x2_PETSC
Vec            :: k1_t_PETSC, k1_b_PETSC, k2_t_PETSC, k2_b_PETSC, help_PETSC
Vec            :: k3_t_PETSC, k3_b_PETSC, k4_t_PETSC, k4_b_PETSC
!==================================================================================================


! in 
!---------- ---------- ---------- ---------- ----------
DIMENSION ID(NJ, NDOF), LTRANS(NTRANS)
!---------- ---------- ---------- ---------- ----------

write(*,*) 'see EXPLICIT_RK4_PETSC_2nd_order_ODE'
stop

! vectors
!---------- ---------- ---------- ---------- ----------
CALL VecCreateMPI ( PETSC_COMM_WORLD, PETSC_DECIDE, 2*NEQ,   x1_PETSC, IERR )
CALL VecCreateMPI ( PETSC_COMM_WORLD, PETSC_DECIDE, 2*NEQ,   x2_PETSC, IERR )
CALL VecCreateMPI ( PETSC_COMM_WORLD, PETSC_DECIDE, 2*NEQ, k1_t_PETSC, IERR )
CALL VecCreateMPI ( PETSC_COMM_WORLD, PETSC_DECIDE, 2*NEQ, k1_b_PETSC, IERR )
CALL VecCreateMPI ( PETSC_COMM_WORLD, PETSC_DECIDE, 2*NEQ, k2_t_PETSC, IERR )
CALL VecCreateMPI ( PETSC_COMM_WORLD, PETSC_DECIDE, 2*NEQ, k2_b_PETSC, IERR )
CALL VecCreateMPI ( PETSC_COMM_WORLD, PETSC_DECIDE, 2*NEQ, help_PETSC, IERR )

CALL VecCreateMPI ( PETSC_COMM_WORLD, PETSC_DECIDE, 2*NEQ, k3_t_PETSC, IERR )
CALL VecCreateMPI ( PETSC_COMM_WORLD, PETSC_DECIDE, 2*NEQ, k3_b_PETSC, IERR )
CALL VecCreateMPI ( PETSC_COMM_WORLD, PETSC_DECIDE, 2*NEQ, k4_t_PETSC, IERR )
CALL VecCreateMPI ( PETSC_COMM_WORLD, PETSC_DECIDE, 2*NEQ, k4_b_PETSC, IERR )


! system is at rest
!---------- ---------- ---------- ---------- ----------
CALL VecSet ( x1_PETSC, 0.0d0, ierr)
CALL VecSet ( x2_PETSC, 0.0d0, ierr)


! index sets
!---------- ---------- ---------- ---------- ----------
CALL ISCreateStride ( PETSC_COMM_WORLD, NEQ, 0  , 1, is_1, IERR )
CALL ISCreateStride ( PETSC_COMM_WORLD, NEQ, NEQ, 1, is_2, IERR )


! inverse of the diagonal mass matrix
!---------- ---------- ---------- ---------- ----------  
CALL VecReciprocal ( DIAG_M_b_PETSC, IERR )


! form right hand side matrices
!---------- ---------- ---------- ---------- ----------
CALL MatDiagonalScale   ( AK_b_PETSC, DIAG_M_b_PETSC, PETSC_NULL_OBJECT, IERR )
CALL MatDiagonalScale   ( AC_b_PETSC, DIAG_M_b_PETSC, PETSC_NULL_OBJECT, IERR )
CALL VecPointwiseMult   (  B_b_PETSC,      B_b_PETSC,    DIAG_M_b_PETSC, IERR )

CALL MatScale ( AK_b_PETSC, -1.0D0, IERR )
CALL MatScale ( AC_b_PETSC, -1.0D0, IERR )


! Explicit Integration
!---------- ---------- ---------- ---------- ----------
DO ISTEP = 1 , NSTEP

   T = DBLE( ISTEP - 1 ) * DT
!   IF ( RANK == 0 .AND. MOD ( ISTEP - 1 , I_PRINT_STEP ) == 0 ) THEN
      WRITE (*,'(A30,F10.5)') 'TIME=', T
!   END IF

   CALL TIME_FUNCTION ( T, fac7 )
   CALL VecSet   ( k1_b_PETSC, 0.0d0          , IERR )
   CALL VecAXPY  ( k1_b_PETSC, fac7, B_b_PETSC, IERR )

   CALL TIME_FUNCTION ( T+dt, fac8 )
   CALL VecSet   ( k4_b_PETSC, 0.0d0          , IERR )
   CALL VecAXPY  ( k4_b_PETSC, fac8, B_b_PETSC, IERR )

   CALL TIME_FUNCTION ( T + 0.5d0 * dt, fac9 )
   CALL VecSet   ( k2_b_PETSC, 0.0d0          , IERR )
   CALL VecAXPY  ( k2_b_PETSC, fac9, B_b_PETSC, IERR )
   CALL VecSet   ( k3_b_PETSC, 0.0d0          , IERR )
   CALL VecAXPY  ( k3_b_PETSC, fac9, B_b_PETSC, IERR )


! RK-4 Explicit time-stepping:
!---------- ---------- ---------- ---------- ----------        
! 1. compute k1_t
   CALL VecCopy ( x2_PETSC, k1_t_PETSC, IERR )
! 2. compute k1_b
   CALL MatMultAdd ( AC_b_PETSC, x2_PETSC, k1_b_PETSC, k1_b_PETSC, IERR )
   CALL MatMultAdd ( AK_b_PETSC, x1_PETSC, k1_b_PETSC, k1_b_PETSC, IERR )
! 3. compute k2_t
   CALL VecWAXPY   ( k2_t_PETSC, 0.5d0*dt, k1_b_PETSC,   x2_PETSC, IERR )
! 4. compute k2_b
   CALL VecWAXPY   ( help_PETSC, 0.5d0*dt, k1_t_PETSC,   x1_PETSC, IERR )
   CALL MatMultAdd ( AC_b_PETSC, k2_t_PETSC, k2_b_PETSC, k2_b_PETSC, IERR )
   CALL MatMultAdd ( AK_b_PETSC, help_PETSC, k2_b_PETSC, k2_b_PETSC, IERR )

! 5. compute k3_t
   CALL VecWAXPY   ( k3_t_PETSC, 0.5d0*dt, k2_b_PETSC, x2_PETSC, IERR )
! 6. compute k3_b
   CALL VecWAXPY   ( help_PETSC, 0.5d0*dt,   k2_t_PETSC,   x1_PETSC, IERR )
   CALL MatMultAdd ( AC_b_PETSC, k3_t_PETSC, k3_b_PETSC, k3_b_PETSC, IERR )
   CALL MatMultAdd ( AK_b_PETSC, help_PETSC, k3_b_PETSC, k3_b_PETSC, IERR )

! 7. compute k4_t
   CALL VecWAXPY   ( k4_t_PETSC, dt, k3_b_PETSC, x2_PETSC, IERR )
! 8. compute k4_b
   CALL VecWAXPY   ( help_PETSC, dt, k3_t_PETSC, x1_PETSC, IERR )
   CALL MatMultAdd ( AC_b_PETSC, k4_t_PETSC, k4_b_PETSC, k4_b_PETSC, IERR )
   CALL MatMultAdd ( AK_b_PETSC, help_PETSC, k4_b_PETSC, k4_b_PETSC, IERR )

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
! get the second sub-vector
   CALL VecGetSubVector ( x1_PETSC, is_1, U_PETSC, IERR )
! output
   CALL OUT_NEWMARK_PETSC ( U_PETSC, ID, LTRANS, T, ENERGY )
! re-store the second sub-vector
   CALL VecRestoreSubVector ( x1_PETSC, is_1, U_PETSC, IERR )
!---------- ---------- ---------- ---------- ---------- 


END DO
!---------- ---------- ---------- ---------- ----------


! destroy PETSc objects
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

CALL VecDestroy (    U_PETSC, IERR )

CALL ISDestroy  ( is_1, IERR )
CALL ISDestroy  ( is_2, IERR )
!---------- ---------- ---------- ---------- ----------


RETURN
END

