!************************************************
! ASSEMBLE PML MATRICES for the adjoint problem (optimize-then-discretize)
!************************************************
! REVISION : W 2 July 2014: If you can't convince them, confuse them.

SUBROUTINE ASSEM_PML_PETSC_3D_Adjoint ( AK_PETSC, AM_PETSC, AC_PETSC, AG_PETSC, DIAG_M_PETSC, XYZ, INOD, NGP, MTEL, PMAT, PMat_Lambda, PMat_Mu, PML_PARAM, ID_BC, ID, SIZE, RANK )
USE PARAMETERS
IMPLICIT DOUBLE PRECISION (A-H,O-Z)

!==================================================================================================
! - PETSC INCLUDE FILES
#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscmat.h"
#include "finclude/petscis.h"

! - PETSC VARIABLES AND MATRICES
ISLocalToGlobalMapping :: GMapping ;             ! Global Mapping.
Mat            :: AK_PETSC, AM_PETSC, AC_PETSC
Mat            :: AG_PETSC
Vec            :: DIAG_M_PETSC

PetscErrorCode :: IERR
PetscMPIInt    :: SIZE, RANK

PetscScalar    :: zero

PetscScalar    :: EM_a_PETSC(NNODE * NDIM, NNODE * NDIM), EM_b_PETSC(NNODE * NDIM, NNODE * NDIM), EM_c_PETSC(NNODE * NDIM, NNODE * NDIM), &
                  EM_d_PETSC(NNODE * NDIM, NNODE * NDIM)

PetscScalar    :: EA_e_PETSC(NNODE * 3, NNODE * 6), EA_p_PETSC(NNODE * 3, NNODE * 6), EA_w_PETSC(NNODE * 3, NNODE * 6)
PetscScalar    :: EN_a_PETSC(NNODE * 6, NNODE * 6), EN_b_PETSC(NNODE * 6, NNODE * 6), EN_c_PETSC(NNODE * 6, NNODE * 6), EN_d_PETSC(NNODE * 6, NNODE * 6)

PetscScalar    :: EA_e_PETSC_T(NNODE * 6, NNODE * 3), EA_p_PETSC_T(NNODE * 6, NNODE * 3), EA_w_PETSC_T(NNODE * 6, NNODE * 3)

PetscScalar    :: EM_a_DIAG_PETSC(NNODE * NDIM), EN_a_DIAG_PETSC(NNODE * 6)
!PetscScalar    :: EM_b_DIAG_PETSC(NNODE * NDIM), EN_b_DIAG_PETSC(NNODE * 6)
!PetscScalar    :: EM_c_DIAG_PETSC(NNODE * NDIM), EN_c_DIAG_PETSC(NNODE * 6)
!PetscScalar    :: EM_d_DIAG_PETSC(NNODE * NDIM), EN_d_DIAG_PETSC(NNODE * 6)
!==================================================================================================


! in 
!---------- ---------- ---------- ---------- ----------
DIMENSION XYZ(NJ,NDIM), INOD(NNODE,NEL), NGP(NEL), MTEL(NEL), PMAT(NMAT,NPM), ID(NJ,NDOF), ID_BC(NEL,NDIM**2)
DIMENSION PML_PARAM(NDIM * 2 , 4)
Dimension PMat_Lambda ( NJ ), PMat_Mu ( NJ )
!---------- ---------- ---------- ---------- ----------


! local
!---------- ---------- ---------- ---------- ----------
DIMENSION EM_a(NNODE * NDIM, NNODE * NDIM), EM_b(NNODE * NDIM, NNODE * NDIM), EM_c(NNODE * NDIM, NNODE * NDIM), EM_d(NNODE * NDIM, NNODE * NDIM)

DIMENSION EA_e_upper(NNODE * 3, NNODE * 6), EA_p_upper(NNODE * 3, NNODE * 6), EA_w_upper(NNODE * 3, NNODE * 6)
DIMENSION EA_e_lower(NNODE * 3, NNODE * 6), EA_p_lower(NNODE * 3, NNODE * 6), EA_w_lower(NNODE * 3, NNODE * 6)
DIMENSION EN_a(NNODE * 6, NNODE * 6), EN_b(NNODE * 6, NNODE * 6), EN_c(NNODE * 6, NNODE * 6), EN_d(NNODE * 6, NNODE * 6)

DIMENSION ND1(NNODE * NDIM), ND2(NNODE * 6)

DIMENSION EM_a_DIAG(NNODE * NDIM), EN_a_DIAG(NNODE * 6)
DIMENSION EM_b_DIAG(NNODE * NDIM), EN_b_DIAG(NNODE * 6)
DIMENSION EM_c_DIAG(NNODE * NDIM), EN_c_DIAG(NNODE * 6)
DIMENSION EM_d_DIAG(NNODE * NDIM), EN_d_DIAG(NNODE * 6)
!---------- ---------- ---------- ---------- ----------  


!DO IEL  = 1 + RANK , NEL , SIZE
DO IEL = 1 , NEL


   MTYPE = MTEL(IEL)
!  MTYPE 1 : REGULAR DOMAIN
!  MTYPE 2 : PML
   IF ( MTYPE /= 2 ) CYCLE


   IF ( MOD ( IEL , 1000 ) == 0 ) THEN
      WRITE(*,'(A30, I10)')'PML DOMAIN - ELEMENT NO.', IEL
!     WRITE(2,'(A30, I10)')'PML DOMAIN - ELEMENT NO.', IEL
   END IF


! PML matrices
!---------- ---------- ---------- ---------- ---------- 
   CALL ELEMENT_M_3D  ( XYZ, INOD, NGP, MTEL, PMAT, PML_PARAM, IEL, EM_a,       202 )    ! PML M_a
   CALL ELEMENT_M_3D  ( XYZ, INOD, NGP, MTEL, PMAT, PML_PARAM, IEL, EM_b,       203 )    ! PML M_b
   CALL ELEMENT_M_3D  ( XYZ, INOD, NGP, MTEL, PMAT, PML_PARAM, IEL, EM_c,       204 )    ! PML M_c
   CALL ELEMENT_M_3D  ( XYZ, INOD, NGP, MTEL, PMAT, PML_PARAM, IEL, EM_d,       205 )    ! PML M_d

   CALL ELEMENT_EA_3D ( XYZ, INOD, NGP, MTEL, PMAT, PMat_Lambda, PMat_Mu, PML_PARAM, IEL, EA_e_upper, 101 )    ! PML EA_e_upper
   CALL ELEMENT_EA_3D ( XYZ, INOD, NGP, MTEL, PMAT, PMat_Lambda, PMat_Mu, PML_PARAM, IEL, EA_p_upper, 102 )    ! PML EA_p_upper
   CALL ELEMENT_EA_3D ( XYZ, INOD, NGP, MTEL, PMAT, PMat_Lambda, PMat_Mu, PML_PARAM, IEL, EA_w_upper, 103 )    ! PML EA_w_upper
   CALL ELEMENT_EA_3D ( XYZ, INOD, NGP, MTEL, PMAT, PMat_Lambda, PMat_Mu, PML_PARAM, IEL, EA_e_lower, 104 )    ! PML EA_e_lower
   CALL ELEMENT_EA_3D ( XYZ, INOD, NGP, MTEL, PMAT, PMat_Lambda, PMat_Mu, PML_PARAM, IEL, EA_p_lower, 105 )    ! PML EA_p_lower
   CALL ELEMENT_EA_3D ( XYZ, INOD, NGP, MTEL, PMAT, PMat_Lambda, PMat_Mu, PML_PARAM, IEL, EA_w_lower, 106 )    ! PML EA_w_lower
  
   CALL ELEMENT_EN_3D ( XYZ, INOD, NGP, MTEL, PMAT, PML_PARAM, IEL, EN_a,       301 )    ! PML EN_a
   CALL ELEMENT_EN_3D ( XYZ, INOD, NGP, MTEL, PMAT, PML_PARAM, IEL, EN_b,       302 )    ! PML EN_b
   CALL ELEMENT_EN_3D ( XYZ, INOD, NGP, MTEL, PMAT, PML_PARAM, IEL, EN_c,       303 )    ! PML EN_c
   CALL ELEMENT_EN_3D ( XYZ, INOD, NGP, MTEL, PMAT, PML_PARAM, IEL, EN_d,       304 )    ! PML EN_d


idebug = 9 
if (idebug == 98) then
   write(97,*) iel
   do i = 1 , nnode * ndim
      write(97,'(16e12.3)') (em_a(i,j) , j = 1 , nnode * ndim)

   end do

   write(98,*) iel
   do i = 1 , nnode * ndim
      write(98,'(24e12.3)') (ea_e(i,j) , j = 1 , nnode * 3)
   end do
  
   write(99,*) iel
   do i = 1 , nnode * 3
      write(99,'(24e12.3)') (en_a(i,j) , j = 1 , nnode * 3)
   end do
end if


! Diagonal mass matrix for Explicit time-stepping: HRZ lumping
! if you use spectral elements, mass matrix is automatically diagonal and these lines basically extract the diagonal part
! W A R N I N G : if you uncomment b,c,d, you'll get NaN for PML; i don't know why.
!---------- ---------- ---------- ---------- ----------
! 1. EM matrices
    CALL Extract_Diagonal ( EM_a, NNODE * NDIM, EM_a_DIAG )
!   CALL HRZ_Mass_Diagonalization ( EM_a, NDIM, EM_a_DIAG )
!   CALL HRZ_Mass_Diagonalization ( EM_b, NDIM, EM_b_DIAG )
!   CALL HRZ_Mass_Diagonalization ( EM_c, NDIM, EM_c_DIAG )
!   CALL HRZ_Mass_Diagonalization ( EM_d, NDIM, EM_d_DIAG )
! 2. EN matrices
    CALL Extract_Diagonal ( EN_a, NNODE *  6  , EN_a_DIAG )
!   CALL HRZ_Mass_Diagonalization ( EN_a,  6  , EN_a_DIAG )
!   CALL HRZ_Mass_Diagonalization ( EN_b,  6  , EN_b_DIAG )
!   CALL HRZ_Mass_Diagonalization ( EN_c,  6  , EN_c_DIAG )
!   CALL HRZ_Mass_Diagonalization ( EN_d,  6  , EN_d_DIAG )


! Send Fortran matrix to PETSc data structure
!---------- ---------- ---------- ---------- ---------- 
   EM_a_PETSC   = EM_a
   EM_b_PETSC   = EM_b
   EM_c_PETSC   = EM_c
   EM_d_PETSC   = EM_d

   EA_e_PETSC   = EA_e_lower
   EA_p_PETSC   = EA_p_lower
   EA_w_PETSC   = EA_w_lower

   EA_e_PETSC_T = transpose(EA_e_upper)
   EA_p_PETSC_T = transpose(EA_p_upper)
   EA_w_PETSC_T = transpose(EA_w_upper)

   EN_a_PETSC   = EN_a
   EN_b_PETSC   = EN_b
   EN_c_PETSC   = EN_c
   EN_d_PETSC   = EN_d

   EM_a_DIAG_PETSC = EM_a_DIAG
   EN_a_DIAG_PETSC = EN_a_DIAG


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


! ASSEMBLE PML MATRICES
! EM : (U , U); EA : (U , S); EN : (S , S)

! assemble EM matrices 
!---------- ---------- ---------- ---------- ----------   
   NEQEL_U = NNODE * NDIM
   CALL MatSetValuesLocal ( AM_PETSC, NEQEL_U, ND1, NEQEL_U, ND1, +EM_a_PETSC, ADD_VALUES, IERR )
   CALL MatSetValuesLocal ( AC_PETSC, NEQEL_U, ND1, NEQEL_U, ND1, -EM_b_PETSC, ADD_VALUES, IERR )
   CALL MatSetValuesLocal ( AK_PETSC, NEQEL_U, ND1, NEQEL_U, ND1, +EM_c_PETSC, ADD_VALUES, IERR )
   CALL MatSetValuesLocal ( AG_PETSC, NEQEL_U, ND1, NEQEL_U, ND1, -EM_d_PETSC, ADD_VALUES, IERR )
!---------- ---------- ---------- ---------- ----------


! assemble EN matrices 
!---------- ---------- ---------- ---------- ----------   
   NEQEL_S = NNODE * 6
   CALL MatSetValuesLocal ( AM_PETSC, NEQEL_S, ND2, NEQEL_S, ND2, +EN_a_PETSC, ADD_VALUES, IERR )
   CALL MatSetValuesLocal ( AC_PETSC, NEQEL_S, ND2, NEQEL_S, ND2, -EN_b_PETSC, ADD_VALUES, IERR )
   CALL MatSetValuesLocal ( AK_PETSC, NEQEL_S, ND2, NEQEL_S, ND2, +EN_c_PETSC, ADD_VALUES, IERR )
   CALL MatSetValuesLocal ( AG_PETSC, NEQEL_S, ND2, NEQEL_S, ND2, -EN_d_PETSC, ADD_VALUES, IERR )
!---------- ---------- ---------- ---------- ----------


! assemble EA matrices
! note: check this part later; it is really tricky (and dirty)! 
! Petsc have a row-oriented matrix strategy!
!---------- ---------- ---------- ---------- ----------   
   NEQEL_U = NNODE * NDIM
   NEQEL_S = NNODE * 6
   CALL MatSetValuesLocal ( AC_PETSC, NEQEL_s, ND2, NEQEL_u, ND1, -EA_e_PETSC, ADD_VALUES, IERR )
   CALL MatSetValuesLocal ( AK_PETSC, NEQEL_s, ND2, NEQEL_u, ND1, +EA_p_PETSC, ADD_VALUES, IERR )
   CALL MatSetValuesLocal ( AG_PETSC, NEQEL_s, ND2, NEQEL_u, ND1, -EA_w_PETSC, ADD_VALUES, IERR )
! assemble transpose
   CALL MatSetValuesLocal ( AC_PETSC, NEQEL_u, ND1, NEQEL_s, ND2, +EA_e_PETSC_T, ADD_VALUES, IERR )
   CALL MatSetValuesLocal ( AK_PETSC, NEQEL_u, ND1, NEQEL_s, ND2, -EA_p_PETSC_T, ADD_VALUES, IERR )
   CALL MatSetValuesLocal ( AG_PETSC, NEQEL_u, ND1, NEQEL_s, ND2, +EA_w_PETSC_T, ADD_VALUES, IERR )
!---------- ---------- ---------- ---------- ----------


! Assemble Diagonal mass for Explicit time stepping
!---------- ---------- ---------- ---------- ----------
   CALL VecSetValuesLocal ( DIAG_M_PETSC, NEQEL_U, ND1, EM_a_DIAG_PETSC, ADD_VALUES, IERR )
   CALL VecSetValuesLocal ( DIAG_M_PETSC, NEQEL_S, ND2, EN_a_DIAG_PETSC, ADD_VALUES, IERR )
!---------- ---------- ---------- ---------- ----------


END DO
!---------- ---------- ---------- ---------- ----------


RETURN 
END
