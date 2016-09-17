!************************************************
! ASSEMBLE REGULAR DOMAIN MASS AND STIFFNESS MATRIX
!************************************************
! REVISION : F, 29 March 2013

SUBROUTINE ASSEM_REGULAR_DOMAIN_PETSC ( AK_PETSC, AM_PETSC, AC_PETSC, AK_RD_PETSC, AM_RD_PETSC, RMJ_PETSC, DIAG_M_PETSC, DIAG_M_RD_PETSC, XYZ, INOD, NGP, MTEL, PMAT, PMat_Lambda, PMat_Mu, PML_PARAM, ID_BC, ID, SIZE, RANK )
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
Mat            :: AK_RD_PETSC, AM_RD_PETSC
Vec            :: RMJ_PETSC
Vec            :: DIAG_M_PETSC
Vec            :: DIAG_M_RD_PETSC

PetscErrorCode :: IERR
PetscMPIInt    :: SIZE, RANK

PetscScalar    :: zero

PetscScalar    :: EK_PETSC(NNODE * NDIM, NNODE * NDIM), EM_PETSC(NNODE * NDIM, NNODE * NDIM)
PetscScalar    :: EC_PETSC(NNODE * NDIM, NNODE * NDIM)
PetscScalar    :: ERMJ_PETSC(NNODE * NDIM)
PetscScalar    :: EDIAG_M_PETSC(NNODE * NDIM)
!==================================================================================================


! in 
!---------- ---------- ---------- ---------- ----------	
DIMENSION XYZ(NJ,NDIM), INOD(NNODE,NEL), NGP(NEL), MTEL(NEL), PMAT(NMAT,NPM), ID(NJ,NDOF), ID_BC(NEL,NDIM**2)
DIMENSION PML_PARAM(NDIM * 2 , 4)
Dimension PMat_Lambda ( NJ ), PMat_Mu ( NJ )
!---------- ---------- ---------- ---------- ----------


! out
!---------- ---------- ---------- ---------- ----------	
!
!---------- ---------- ---------- ---------- ----------


! local
!---------- ---------- ---------- ---------- ----------	
DIMENSION EK(NNODE * NDIM, NNODE * NDIM), EM(NNODE * NDIM, NNODE * NDIM), ND1(NNODE * NDIM)
DIMENSION EC(NNODE * NDIM, NNODE * NDIM)
DIMENSION ERMJ(NNODE * NDIM)
DIMENSION EDIAG_M(NNODE * NDIM)
!---------- ---------- ---------- ---------- ----------

DO IEL = 1 , NEL

   MTYPE = MTEL(IEL)
!  MTYPE 1 : REGULAR DOMAIN
!  MTYPE 2 : PML
   IF ( MTYPE /= 1 ) CYCLE


   IF ( MOD ( IEL , 1000 ) == 0 ) THEN
      WRITE(*,'(A30, I10)')'REGULAR DOMAIN - ELEMENT NO.', IEL
!     WRITE(2,'(A30, I10)')'REGULAR DOMAIN - ELEMENT NO.', IEL
   END IF

! initialize the right hand side "force" vector
!---------- ---------- ---------- ---------- ---------- 
   ERMJ = 0.0d0


! initialize this since we are assembling that; we don't want to insert garbage into C, don't touch this; modify in "Rayleigh damping" if necessary
!---------- ---------- ---------- ---------- ----------
   EC   = 0.0d0


! regular domain EM, EK
! flag 201 : regular domain mass matrix
!---------- ---------- ---------- ---------- ----------
   CALL ELEMENT_M ( XYZ, INOD, NGP, MTEL, PMAT, PML_PARAM, IEL, EM, 201 )
   CALL ELEMENT_K ( XYZ, INOD, NGP, MTEL, PMAT, PMat_Lambda, PMat_Mu, IEL, EK )


idebug = 9
if (idebug == 96) then
   write(*,*) iel
   do i = 1,nnode*ndim
      write(*,'(16e12.3)') (ek(i,j),j=1,nnode*ndim)
   end do
   write(*,*) iel
   do i = 1,nnode*ndim
      write(*,'(16e12.3)') (ek(i,j),j=1,nnode*ndim)
   end do
end if


! Rayleigh damping
!---------- ---------- ---------- ---------- ----------	
!   EC = ALPHA_DAMP * EM + BETA_DAMP * EK


! form right hand side "force" vector
! assumption: only x-component of acceleration is considered i.e. ag = ( ax , 0 )
! notice the minus sign: - M J ag
!---------- ---------- ---------- ---------- ----------
   DO I = 1 , NNODE
      DO J = 1 , NNODE
         ERMJ(I) = ERMJ(I) - EM(I,J)
      END DO
   END DO


! Diagonal mass matrix for Explicit time-stepping: HRZ lumping
! if you use spectral elements, mass matrix is diagonal and these lines basically extract the diagonal part
!---------- ---------- ---------- ---------- ----------
   CALL Extract_Diagonal ( EM, NNODE * NDIM, EDIAG_M )
!  CALL HRZ_Mass_Diagonalization ( EM, NDIM, EDIAG_M )


! Send Fortran matrix to PETSc data structure
!---------- ---------- ---------- ---------- ---------- 
   EK_PETSC      = EK
   EM_PETSC      = EM
   EC_PETSC      = EC

   ERMJ_PETSC    = ERMJ
   EDIAG_M_PETSC = EDIAG_M


!---------- ---------- ---------- ---------- ----------	  
   DO I = 1 , NNODE
      I1 = I
      I2 = I1 + NNODE

      ND1(I1) = ID(INOD(I,IEL) , 1)
      ND1(I2) = ID(INOD(I,IEL) , 2)
   END DO
!---------- ---------- ---------- ---------- ----------


! Modify ND according to C indexing
!---------- ---------- ---------- ---------- ----------
   DO I = 1 , NNODE * NDIM 
      ND1(I) = ND1(I) - 1
   END DO


! ASSEMBLE MASS, DAMPING AND STIFFNESS MATRIX
!---------- ---------- ---------- ---------- ----------	  
   NEQEL = NNODE * NDIM
   CALL MatSetValuesLocal ( AK_PETSC, NEQEL, ND1, NEQEL, ND1, EK_PETSC, ADD_VALUES, IERR )
   CALL MatSetValuesLocal ( AM_PETSC, NEQEL, ND1, NEQEL, ND1, EM_PETSC, ADD_VALUES, IERR )
!  CALL MatSetValues ( AC_PETSC, NEQEL, ND1, NEQEL, ND1, EC_PETSC, ADD_VALUES, IERR )

! Matrices that correspond to the Regular Domain only. This is for energy decay assessment within the RD.
   CALL MatSetValuesLocal ( AK_RD_PETSC, NEQEL, ND1, NEQEL, ND1, EK_PETSC, ADD_VALUES, IERR )
   CALL MatSetValuesLocal ( AM_RD_PETSC, NEQEL, ND1, NEQEL, ND1, EM_PETSC, ADD_VALUES, IERR )
!---------- ---------- ---------- ---------- ----------


! ASSEMBLE EXACT DYNAMIC FORCE vector (=-MJU") and Diagonal mass for Explicit time stepping
! only x-component of acceleration is considered in force 
!---------- ---------- ---------- ---------- ----------
   CALL VecSetValuesLocal ( RMJ_PETSC   , NEQEL, ND1, ERMJ_PETSC   , ADD_VALUES, IERR )
   CALL VecSetValuesLocal ( DIAG_M_PETSC, NEQEL, ND1, EDIAG_M_PETSC, ADD_VALUES, IERR )

! Matrices that correspond to the Regular Domain only. This is for energy decay assessment within the RD.
   CALL VecSetValuesLocal ( DIAG_M_RD_PETSC, NEQEL, ND1, EDIAG_M_PETSC, ADD_VALUES, IERR )
!---------- ---------- ---------- ---------- ----------	


END DO
!---------- ---------- ---------- ---------- ----------	


RETURN
END


!************************************************
! ASSEMBLE PML MATRICES
!************************************************
! REVISION : Tue, 24 July 2012

SUBROUTINE ASSEM_PML_PETSC ( AK_PETSC, AM_PETSC, AC_PETSC, DIAG_M_PETSC, XYZ, INOD, NGP, MTEL, PMAT, PMat_Lambda, PMat_Mu, PML_PARAM, ID_BC, ID, SIZE, RANK )
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
Vec            :: DIAG_M_PETSC

PetscErrorCode :: IERR
PetscMPIInt    :: SIZE, RANK

PetscScalar    :: zero

PetscScalar    :: EM_a_PETSC(NNODE * NDIM, NNODE * NDIM), EM_b_PETSC(NNODE * NDIM, NNODE * NDIM), EM_c_PETSC(NNODE * NDIM, NNODE * NDIM)
PetscScalar    :: EA_e_PETSC(NNODE * 2, NNODE * 3), EA_p_PETSC(NNODE * 2, NNODE * 3)
PetscScalar    :: EN_a_PETSC(NNODE * 3, NNODE * 3), EN_b_PETSC(NNODE * 3, NNODE * 3), EN_c_PETSC(NNODE * 3, NNODE * 3)

PetscScalar    :: EA_e_PETSC_T(NNODE * 3, NNODE * 2), EA_p_PETSC_T(NNODE * 3, NNODE * 2)

PetscScalar    :: EM_a_DIAG_PETSC(NNODE * NDIM), EN_a_DIAG_PETSC(NNODE * 3)
!PetscScalar    :: EM_b_DIAG_PETSC(NNODE * NDIM), EN_b_DIAG_PETSC(NNODE * 3)
!PetscScalar    :: EM_c_DIAG_PETSC(NNODE * NDIM), EN_c_DIAG_PETSC(NNODE * 3)
!==================================================================================================


! in 
!---------- ---------- ---------- ---------- ----------
DIMENSION XYZ(NJ,NDIM), INOD(NNODE,NEL), NGP(NEL), MTEL(NEL), PMAT(NMAT,NPM), ID(NJ,NDOF), ID_BC(NEL,NDIM**2)
DIMENSION PML_PARAM(NDIM * 2 , 4)
Dimension PMat_Lambda (NJ), PMat_Mu(NJ)
!---------- ---------- ---------- ---------- ----------


! local
!---------- ---------- ---------- ---------- ----------
DIMENSION EM_a(NNODE * NDIM, NNODE * NDIM), EM_b(NNODE * NDIM, NNODE * NDIM), EM_c(NNODE * NDIM, NNODE * NDIM)
DIMENSION EA_e_upper(NNODE * 2, NNODE * 3), EA_p_upper(NNODE * 2, NNODE * 3)
DIMENSION EA_e_lower(NNODE * 2, NNODE * 3), EA_p_lower(NNODE * 2, NNODE * 3)
DIMENSION EN_a(NNODE * 3, NNODE * 3), EN_b(NNODE * 3, NNODE * 3), EN_c(NNODE * 3, NNODE * 3)

DIMENSION ND1(NNODE * NDIM), ND2(NNODE * 3)

DIMENSION EM_a_DIAG(NNODE * NDIM), EN_a_DIAG(NNODE * 3)
DIMENSION EM_b_DIAG(NNODE * NDIM), EN_b_DIAG(NNODE * 3)
DIMENSION EM_c_DIAG(NNODE * NDIM), EN_c_DIAG(NNODE * 3)
!---------- ---------- ---------- ---------- ----------  
  

DO IEL = 1 , NEL


   MTYPE = MTEL(IEL)
!  MTYPE 1 : REGULAR DOMAIN
!  MTYPE 2 : PML
   IF ( MTYPE /= 2 ) CYCLE

   IF ( MOD ( IEL , 1000 ) == 0 ) THEN
      WRITE(*,'(A30, I10)')'PML - ELEMENT NO.', IEL
!     WRITE(2,'(A30, I10)')'PML - ELEMENT NO.', IEL
   END IF

! PML matrices
!---------- ---------- ---------- ---------- ---------- 
   CALL ELEMENT_M  (XYZ, INOD, NGP, MTEL, PMAT, PML_PARAM, IEL, EM_a,       202)         ! PML M_a
   CALL ELEMENT_M  (XYZ, INOD, NGP, MTEL, PMAT, PML_PARAM, IEL, EM_b,       203)         ! PML M_b
   CALL ELEMENT_M  (XYZ, INOD, NGP, MTEL, PMAT, PML_PARAM, IEL, EM_c,       204)         ! PML M_c
  
   CALL ELEMENT_EA (XYZ, INOD, NGP, MTEL, PMAT, PMat_Lambda, PMat_Mu, PML_PARAM, IEL, EA_e_upper, 101)         ! PML EA_e_upper
   CALL ELEMENT_EA (XYZ, INOD, NGP, MTEL, PMAT, PMat_Lambda, PMat_Mu, PML_PARAM, IEL, EA_p_upper, 102)         ! PML EA_p_upper
   CALL ELEMENT_EA (XYZ, INOD, NGP, MTEL, PMAT, PMat_Lambda, PMat_Mu, PML_PARAM, IEL, EA_e_lower, 103)         ! PML EA_e_lower
   CALL ELEMENT_EA (XYZ, INOD, NGP, MTEL, PMAT, PMat_Lambda, PMat_Mu, PML_PARAM, IEL, EA_p_lower, 104)         ! PML EA_p_lower
  
   CALL ELEMENT_EN (XYZ, INOD, NGP, MTEL, PMAT, PML_PARAM, IEL, EN_a,       301)         ! PML EN_a
   CALL ELEMENT_EN (XYZ, INOD, NGP, MTEL, PMAT, PML_PARAM, IEL, EN_b,       302)         ! PML EN_b
   CALL ELEMENT_EN (XYZ, INOD, NGP, MTEL, PMAT, PML_PARAM, IEL, EN_c,       303)         ! PML EN_c


idebug = 9
if (idebug == 98) then
   write(97,*) iel
   do i = 1 , nnode*ndim
      write(97,'(16e12.3)') (em_a(i,j) , j = 1 , nnode*ndim)
   end do

   write(98,*) iel
   do i = 1 , nnode*ndim
      write(98,'(24e12.3)') (ea_e_upper(i,j) , j = 1 , nnode * 3)
   end do
  
   write(99,*) iel
   do i = 1 , nnode * 3
      write(99,'(24e12.3)') (en_a(i,j) , j = 1 , nnode * 3)
   end do
end if


! Diagonal mass matrix for Explicit time-stepping: HRZ lumping
! if you use spectral elements, mass matrix is diagonal and these lines basically extract the diagonal part
!---------- ---------- ---------- ---------- ----------
! 1. EM matrices
    CALL Extract_Diagonal ( EM_a, NNODE * NDIM, EM_a_DIAG )
!   CALL HRZ_Mass_Diagonalization ( EM_a, NDIM, EM_a_DIAG )
!   CALL HRZ_Mass_Diagonalization ( EM_b, NDIM, EM_b_DIAG )
!   CALL HRZ_Mass_Diagonalization ( EM_c, NDIM, EM_c_DIAG )

! 2. EN matrices
    CALL Extract_Diagonal ( EN_a, NNODE *  3  , EN_a_DIAG )
!   CALL HRZ_Mass_Diagonalization ( EN_a, 3   , EN_a_DIAG )
!   CALL HRZ_Mass_Diagonalization ( EN_b, 3   , EN_b_DIAG )
!   CALL HRZ_Mass_Diagonalization ( EN_c, 3   , EN_c_DIAG )


! Send Fortran matrix to PETSc data structure
!---------- ---------- ---------- ---------- ---------- 
   EM_a_PETSC   = EM_a
   EM_b_PETSC   = EM_b
   EM_c_PETSC   = EM_c

   EA_e_PETSC   = EA_e_upper
   EA_p_PETSC   = EA_p_upper

   EA_e_PETSC_T = transpose(EA_e_lower)
   EA_p_PETSC_T = transpose(EA_p_lower)

   EN_a_PETSC   = EN_a
   EN_b_PETSC   = EN_b
   EN_c_PETSC   = EN_c

   EM_a_DIAG_PETSC = EM_a_DIAG
   EN_a_DIAG_PETSC = EN_a_DIAG


!---------- ---------- ---------- ---------- ----------	  
   DO I = 1 , NNODE

      I1 = I
      I2 = I1 + NNODE
      I3 = I2 + NNODE
      
      ND1(I1) = ID(INOD(I,IEL) , 1)
      ND1(I2) = ID(INOD(I,IEL) , 2)
      
      ND2(I1) = ID(INOD(I,IEL) , 3)
      ND2(I2) = ID(INOD(I,IEL) , 4)
      ND2(I3) = ID(INOD(I,IEL) , 5)
     
   END DO
!---------- ---------- ---------- ---------- ----------


! Modify ND according to C indexing
!---------- ---------- ---------- ---------- ----------
   DO I = 1 , NNODE * NDIM
      ND1(I) = ND1(I) - 1 
   END DO
   DO I = 1 , NNODE * 3
      ND2(I) = ND2(I) - 1
   END DO


! ASSEMBLE PML MATRICES
! EM : (U , U); EA : (U , S); EN : (S , S)

! assemble EM matrices
!---------- ---------- ---------- ---------- ----------
   NEQEL_U = NNODE * NDIM
   CALL MatSetValuesLocal ( AM_PETSC, NEQEL_U, ND1, NEQEL_U, ND1, EM_a_PETSC, ADD_VALUES, IERR )
   CALL MatSetValuesLocal ( AC_PETSC, NEQEL_U, ND1, NEQEL_U, ND1, EM_b_PETSC, ADD_VALUES, IERR )
   CALL MatSetValuesLocal ( AK_PETSC, NEQEL_U, ND1, NEQEL_U, ND1, EM_c_PETSC, ADD_VALUES, IERR )
!---------- ---------- ---------- ---------- ----------


! assemble EN matrices
!---------- ---------- ---------- ---------- ----------
   NEQEL_S = NNODE * 3
   CALL MatSetValuesLocal ( AM_PETSC, NEQEL_S, ND2, NEQEL_S, ND2, EN_a_PETSC, ADD_VALUES, IERR )
   CALL MatSetValuesLocal ( AC_PETSC, NEQEL_S, ND2, NEQEL_S, ND2, EN_b_PETSC, ADD_VALUES, IERR )
   CALL MatSetValuesLocal ( AK_PETSC, NEQEL_S, ND2, NEQEL_S, ND2, EN_c_PETSC, ADD_VALUES, IERR )
!---------- ---------- ---------- ---------- ----------


! assemble EA matrices
! note: check this part later; it is really tricky (and dirty)! 
! Petsc have a row-oriented matrix strategy!
!---------- ---------- ---------- ---------- ----------   
   NEQEL_U = NNODE * NDIM
   NEQEL_S = NNODE * 3
   CALL MatSetValuesLocal ( AC_PETSC, NEQEL_S, ND2, NEQEL_U, ND1, EA_e_PETSC, ADD_VALUES, IERR )
   CALL MatSetValuesLocal ( AK_PETSC, NEQEL_S, ND2, NEQEL_U, ND1, EA_p_PETSC, ADD_VALUES, IERR )
! assemble transpose
   CALL MatSetValuesLocal ( AC_PETSC, NEQEL_U, ND1, NEQEL_S, ND2, -EA_e_PETSC_T, ADD_VALUES, IERR )
   CALL MatSetValuesLocal ( AK_PETSC, NEQEL_U, ND1, NEQEL_S, ND2, -EA_p_PETSC_T, ADD_VALUES, IERR )
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


!************************************************
! ASSEMBLE REGULAR DOMAIN FORCE
!************************************************
! REVISION : M, 18 March 2013

SUBROUTINE ASSEM_FORCE_PETSC ( B_PETSC, XYZ, INOD, NGP, MTEL, PMAT, ID_BC, ID, SIZE, RANK )
USE PARAMETERS
IMPLICIT DOUBLE PRECISION (A-H,O-Z)


!==================================================================================================
! - PETSC INCLUDE FILES
#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscis.h"

! - PETSC VARIABLES AND MATRICES
ISLocalToGlobalMapping :: GMapping ;             ! Global Mapping.
Vec            :: B_PETSC

PetscErrorCode :: IERR
PetscMPIInt    :: SIZE, RANK

PetscScalar    :: EFS_PETSC ( NNODE * NDIM )
!==================================================================================================


! in 
!---------- ---------- ---------- ---------- ---------- 
DIMENSION XYZ(NJ,NDIM), INOD(NNODE,NEL), NGP(NEL), MTEL(NEL), PMAT(NMAT,NPM), ID(NJ,NDOF), ID_BC(NEL,NDIM**2)
!---------- ---------- ---------- ---------- ----------


! local
!---------- ---------- ---------- ---------- ---------- 
DIMENSION EFS(NNODE * NDIM), ND1(NNODE * NDIM)
DIMENSION EFS_surface_traction(NNODE * NDIM), EFS_disc(NNODE * NDIM)
!---------- ---------- ---------- ---------- ----------


DO IEL  = 1 , NEL

   MTYPE = MTEL(IEL)
!  MTYPE 1 : REGULAR DOMAIN
!  MTYPE 2 : PML
   IF ( MTYPE /= 1 ) CYCLE

!  proceed if there is boundary traction         
   NFACE = ID_BC(IEL, 1)
   IF (NFACE == 0) CYCLE


   WRITE(*,'(A30, I10)')'SURFACE TRACTION - ELEMENT NO.', IEL
!  WRITE(2,'(A30, I10)')'SURFACE TRACTION - ELEMENT NO.', IEL
        
        
! surface traction or disc body force (just one of them may be attributed to an element)
!---------- ---------- ---------- ---------- ----------           
   CALL ELEMENT_FS      ( XYZ, INOD, NGP, MTEL, PMAT, ID_BC, IEL, EFS_surface_traction )
   CALL ELEMENT_FB_disc ( XYZ, INOD, NGP, MTEL, PMAT, ID_BC, IEL, EFS_disc )


! Send Fortran matrix to PETSc data structure
!---------- ---------- ---------- ---------- ---------- 
   EFS_PETSC   = EFS_surface_traction + EFS_disc


!---------- ---------- ---------- ---------- ----------   
   DO I = 1 , NNODE
      I1 = I
      I2 = I1 + NNODE

      ND1(I1) = ID(INOD(I,IEL) , 1)
      ND1(I2) = ID(INOD(I,IEL) , 2)
   END DO
!---------- ---------- ---------- ---------- ----------


! Modify ND according to C indexing
!---------- ---------- ---------- ---------- ----------
   DO I = 1 , NNODE * NDIM
      ND1(I) = ND1(I) - 1
   END DO


idebug = 9
if (idebug == 96) then
   write(96,*) iel
   do i = 1,nnode*ndim
      write(*,'(i5,f12.3)') nd1(i), efs_petsc(i)
   end do
end if


! ASSEMBLE surface traction
!---------- ---------- ---------- ---------- ----------   
   NEQEL = NNODE * NDIM
   CALL VecSetValuesLocal ( B_PETSC, NEQEL, ND1, EFS_PETSC, ADD_VALUES, IERR )
!---------- ---------- ---------- ---------- ----------


END DO
!---------- ---------- ---------- ---------- ----------


RETURN
END

!************************************************
! ASSEMBLE PML MATRICES - add temporal damping for stability (not helpful)
!************************************************
! REVISION : TH, 2 May 2013

SUBROUTINE ASSEM_PML_PETSC_TD ( AK_PETSC, AM_PETSC, AC_PETSC, AG_PETSC, DIAG_M_PETSC, XYZ, INOD, NGP, MTEL, PMAT, PML_PARAM, ID_BC, ID, SIZE, RANK )
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
Mat            :: AK_PETSC, AM_PETSC, AC_PETSC, AG_PETSC
Vec            :: DIAG_M_PETSC

PetscErrorCode :: IERR
PetscMPIInt    :: SIZE, RANK

PetscScalar    :: zero

PetscScalar    :: EM_a_PETSC(NNODE * NDIM, NNODE * NDIM), EM_b_PETSC(NNODE * NDIM, NNODE * NDIM), EM_c_PETSC(NNODE * NDIM, NNODE * NDIM)
PetscScalar    :: EM_d_PETSC(NNODE * NDIM, NNODE * NDIM)             ! extra term that involves displacement history  

PetscScalar    :: EA_e_PETSC(NNODE * 2, NNODE * 3), EA_p_PETSC(NNODE * 2, NNODE * 3)
PetscScalar    :: EN_a_PETSC(NNODE * 3, NNODE * 3), EN_b_PETSC(NNODE * 3, NNODE * 3), EN_c_PETSC(NNODE * 3, NNODE * 3)

PetscScalar    :: EA_e_PETSC_T(NNODE * 3, NNODE * 2), EA_p_PETSC_T(NNODE * 3, NNODE * 2)

PetscScalar    :: EM_a_DIAG_PETSC(NNODE * NDIM), EN_a_DIAG_PETSC(NNODE * 3)
!PetscScalar    :: EM_b_DIAG_PETSC(NNODE * NDIM), EN_b_DIAG_PETSC(NNODE * 3)
!PetscScalar    :: EM_c_DIAG_PETSC(NNODE * NDIM), EN_c_DIAG_PETSC(NNODE * 3)
!==================================================================================================


! in 
!---------- ---------- ---------- ---------- ----------
DIMENSION XYZ(NJ,NDIM), INOD(NNODE,NEL), NGP(NEL), MTEL(NEL), PMAT(NMAT,NPM), ID(NJ,NDOF), ID_BC(NEL,NDIM**2)
DIMENSION PML_PARAM(NDIM * 2 , 4)
!---------- ---------- ---------- ---------- ----------


! local
!---------- ---------- ---------- ---------- ----------
DIMENSION EM_a(NNODE * NDIM, NNODE * NDIM), EM_b(NNODE * NDIM, NNODE * NDIM), EM_c(NNODE * NDIM, NNODE * NDIM)
DIMENSION EA_e_upper(NNODE * 2, NNODE * 3), EA_p_upper(NNODE * 2, NNODE * 3)
DIMENSION EA_e_lower(NNODE * 2, NNODE * 3), EA_p_lower(NNODE * 2, NNODE * 3)
DIMENSION EN_a(NNODE * 3, NNODE * 3), EN_b(NNODE * 3, NNODE * 3), EN_c(NNODE * 3, NNODE * 3)

DIMENSION ND1(NNODE * NDIM), ND2(NNODE * 3)

DIMENSION EM_a_DIAG(NNODE * NDIM), EN_a_DIAG(NNODE * 3)
DIMENSION EM_b_DIAG(NNODE * NDIM), EN_b_DIAG(NNODE * 3)
DIMENSION EM_c_DIAG(NNODE * NDIM), EN_c_DIAG(NNODE * 3)
!---------- ---------- ---------- ---------- ----------  
  

DO IEL  = 1 , NEL

   MTYPE = MTEL(IEL)
!  MTYPE 1 : REGULAR DOMAIN
!  MTYPE 2 : PML
   IF ( MTYPE /= 2 ) CYCLE

   WRITE(*,'(A30, I10)')'PML - ELEMENT NO.', IEL
!  WRITE(2,'(A30, I10)')'PML - ELEMENT NO.', IEL


! PML matrices
!---------- ---------- ---------- ---------- ---------- 
   CALL ELEMENT_M  (XYZ, INOD, NGP, MTEL, PMAT, PML_PARAM, IEL, EM_a,       202)         ! PML M_a
   CALL ELEMENT_M  (XYZ, INOD, NGP, MTEL, PMAT, PML_PARAM, IEL, EM_b,       203)         ! PML M_b
   CALL ELEMENT_M  (XYZ, INOD, NGP, MTEL, PMAT, PML_PARAM, IEL, EM_c,       204)         ! PML M_c
  
   CALL ELEMENT_EA (XYZ, INOD, NGP, MTEL, PMAT, PML_PARAM, IEL, EA_e_upper, 101)         ! PML EA_e_upper
   CALL ELEMENT_EA (XYZ, INOD, NGP, MTEL, PMAT, PML_PARAM, IEL, EA_p_upper, 102)         ! PML EA_p_upper
   CALL ELEMENT_EA (XYZ, INOD, NGP, MTEL, PMAT, PML_PARAM, IEL, EA_e_lower, 103)         ! PML EA_e_lower
   CALL ELEMENT_EA (XYZ, INOD, NGP, MTEL, PMAT, PML_PARAM, IEL, EA_p_lower, 104)         ! PML EA_p_lower
  
   CALL ELEMENT_EN (XYZ, INOD, NGP, MTEL, PMAT, PML_PARAM, IEL, EN_a,       301)         ! PML EN_a
   CALL ELEMENT_EN (XYZ, INOD, NGP, MTEL, PMAT, PML_PARAM, IEL, EN_b,       302)         ! PML EN_b
   CALL ELEMENT_EN (XYZ, INOD, NGP, MTEL, PMAT, PML_PARAM, IEL, EN_c,       303)         ! PML EN_c


idebug = 9 
if (idebug == 98) then
   write(97,*) iel
   do i = 1 , nnode*ndim
      write(97,'(16e12.3)') (em_a(i,j) , j = 1 , nnode*ndim)
   end do

   write(98,*) iel
   do i = 1 , nnode*ndim
      write(98,'(24e12.3)') (ea_e(i,j) , j = 1 , nnode * 3)

   end do
  
   write(99,*) iel
   do i = 1 , nnode * 3
      write(99,'(24e12.3)') (en_a(i,j) , j = 1 , nnode * 3)
   end do
end if


! Diagonal mass matrix for Explicit time-stepping: HRZ lumping
! if you use spectral elements, mass matrix is diagonal and these lines basically extract the diagonal part
!---------- ---------- ---------- ---------- ----------
! 1. EM matrices
    CALL Extract_Diagonal ( EM_a, NNODE * NDIM, EM_a_DIAG )
!   CALL HRZ_Mass_Diagonalization ( EM_a, NDIM, EM_a_DIAG )
!   CALL HRZ_Mass_Diagonalization ( EM_b, NDIM, EM_b_DIAG )
!   CALL HRZ_Mass_Diagonalization ( EM_c, NDIM, EM_c_DIAG )

! 2. EN matrices
    CALL Extract_Diagonal ( EN_a, NNODE *  3  , EN_a_DIAG )
!   CALL HRZ_Mass_Diagonalization ( EN_a, 3   , EN_a_DIAG )
!   CALL HRZ_Mass_Diagonalization ( EN_b, 3   , EN_b_DIAG )
!   CALL HRZ_Mass_Diagonalization ( EN_c, 3   , EN_c_DIAG )


! Send Fortran matrix to PETSc data structure
!---------- ---------- ---------- ---------- ---------- 
   EM_a_PETSC   = EM_a
   EM_b_PETSC   = EM_b + ALPHA_DAMP * EM_a
   EM_c_PETSC   = EM_c + ALPHA_DAMP * EM_b
   EM_d_PETSC   =        ALPHA_DAMP * EM_c

   EA_e_PETSC   = EA_e_upper
   EA_p_PETSC   = EA_p_upper

   EA_e_PETSC_T = transpose(EA_e_lower)
   EA_p_PETSC_T = transpose(EA_p_lower)

   EN_a_PETSC   = EN_a
   EN_b_PETSC   = EN_b
   EN_c_PETSC   = EN_c

   EM_a_DIAG_PETSC = EM_a_DIAG
   EN_a_DIAG_PETSC = EN_a_DIAG


!---------- ---------- ---------- ---------- ----------	  
   DO I = 1 , NNODE

      I1 = I
      I2 = I1 + NNODE
      I3 = I2 + NNODE
      
      ND1(I1) = ID(INOD(I,IEL) , 1)
      ND1(I2) = ID(INOD(I,IEL) , 2)
      
      ND2(I1) = ID(INOD(I,IEL) , 3)
      ND2(I2) = ID(INOD(I,IEL) , 4)
      ND2(I3) = ID(INOD(I,IEL) , 5)
     
   END DO
!---------- ---------- ---------- ---------- ----------


! Modify ND according to C indexing
!---------- ---------- ---------- ---------- ----------
   DO I = 1 , NNODE * NDIM
      ND1(I) = ND1(I) - 1 
   END DO
   DO I = 1 , NNODE * 3
      ND2(I) = ND2(I) - 1
   END DO


! ASSEMBLE PML MATRICES
! EM : (U , U); EA : (U , S); EN : (S , S)

! assemble EM matrices 
!---------- ---------- ---------- ---------- ----------   
   NEQEL_U = NNODE * NDIM
   CALL MatSetValuesLocal ( AM_PETSC, NEQEL_U, ND1, NEQEL_U, ND1, EM_a_PETSC, ADD_VALUES, IERR )
   CALL MatSetValuesLocal ( AC_PETSC, NEQEL_U, ND1, NEQEL_U, ND1, EM_b_PETSC, ADD_VALUES, IERR )
   CALL MatSetValuesLocal ( AK_PETSC, NEQEL_U, ND1, NEQEL_U, ND1, EM_c_PETSC, ADD_VALUES, IERR )
! the extra term
   CALL MatSetValuesLocal ( AG_PETSC, NEQEL_U, ND1, NEQEL_U, ND1, EM_d_PETSC, ADD_VALUES, IERR )
!---------- ---------- ---------- ---------- ----------


! assemble EN matrices 
!---------- ---------- ---------- ---------- ----------   
   NEQEL_S = NNODE * 3
   CALL MatSetValuesLocal ( AM_PETSC, NEQEL_S, ND2, NEQEL_S, ND2, EN_a_PETSC, ADD_VALUES, IERR )
   CALL MatSetValuesLocal ( AC_PETSC, NEQEL_S, ND2, NEQEL_S, ND2, EN_b_PETSC, ADD_VALUES, IERR )
   CALL MatSetValuesLocal ( AK_PETSC, NEQEL_S, ND2, NEQEL_S, ND2, EN_c_PETSC, ADD_VALUES, IERR )
!---------- ---------- ---------- ---------- ----------


! assemble EA matrices
! note: check this part later; it is really tricky (and dirty)! 
! Petsc have a row-oriented matrix strategy!
!---------- ---------- ---------- ---------- ----------   
   NEQEL_U = NNODE * NDIM
   NEQEL_S = NNODE * 3
   CALL MatSetValuesLocal ( AC_PETSC, NEQEL_S, ND2, NEQEL_U, ND1, EA_e_PETSC, ADD_VALUES, IERR )
   CALL MatSetValuesLocal ( AK_PETSC, NEQEL_S, ND2, NEQEL_U, ND1, EA_p_PETSC, ADD_VALUES, IERR )
! assemble transpose
   CALL MatSetValuesLocal ( AC_PETSC, NEQEL_U, ND1, NEQEL_S, ND2, -EA_e_PETSC_T, ADD_VALUES, IERR )
   CALL MatSetValuesLocal ( AK_PETSC, NEQEL_U, ND1, NEQEL_S, ND2, -EA_p_PETSC_T, ADD_VALUES, IERR )
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


!************************************************
! ASSEMBLE Elastodynamic matrices in mixed formulation
! This helps us understand if our mixed formulation is stable.
!************************************************
! REVISION : W, 8 May 2013

SUBROUTINE ASSEM_MIXED_FEM_PETSC ( AK_PETSC, AM_PETSC, AC_PETSC, DIAG_M_PETSC, XYZ, INOD, NGP, MTEL, PMAT, PML_PARAM, ID_BC, ID, SIZE, RANK )
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
Vec            :: DIAG_M_PETSC

PetscErrorCode :: IERR
PetscMPIInt    :: SIZE, RANK

PetscScalar    :: EM_PETSC(NNODE * NDIM, NNODE * NDIM)
PetscScalar    :: EA_PETSC(NNODE * 2, NNODE * 3), EA_PETSC_T(NNODE * 3, NNODE * 2)
PetscScalar    :: EN_PETSC(NNODE * 3, NNODE * 3)


PetscScalar    :: EM_DIAG_PETSC(NNODE * NDIM), EN_DIAG_PETSC(NNODE * 3)
!==================================================================================================


! in 
!---------- ---------- ---------- ---------- ----------
DIMENSION XYZ(NJ,NDIM), INOD(NNODE,NEL), NGP(NEL), MTEL(NEL), PMAT(NMAT,NPM), ID(NJ,NDOF), ID_BC(NEL,NDIM**2)
DIMENSION PML_PARAM(NDIM * 2 , 4)
!---------- ---------- ---------- ---------- ----------


! local
!---------- ---------- ---------- ---------- ----------
DIMENSION EM(NNODE * NDIM, NNODE * NDIM)
DIMENSION EA_upper(NNODE * 2, NNODE * 3)
DIMENSION EA_lower(NNODE * 2, NNODE * 3)
DIMENSION EN(NNODE * 3, NNODE * 3)

DIMENSION ND1(NNODE * NDIM), ND2(NNODE * 3)

DIMENSION EM_DIAG(NNODE * NDIM), EN_DIAG(NNODE * 3)
!---------- ---------- ---------- ---------- ----------  
  

DO IEL  = 1 , NEL

   MTYPE = MTEL(IEL)
!  MTYPE 1 : REGULAR DOMAIN
!  MTYPE 2 : PML
   IF ( MTYPE /= 2 ) CYCLE


   WRITE(*,'(A30, I10)')'Mixed FE Elastodynamics - ELEMENT NO.', IEL
!  WRITE(2,'(A30, I10)')'Mixed FE Elastodynamics - ELEMENT NO.', IEL


! Mixed elastodynamic matrices
!---------- ---------- ---------- ---------- ---------- 
   CALL ELEMENT_M ( XYZ, INOD, NGP, MTEL, PMAT, PML_PARAM, IEL, EM,       201)           ! M mass matrix

   CALL ELEMENT_EA (XYZ, INOD, NGP, MTEL, PMAT, PML_PARAM, IEL, EA_upper, 110)           ! EA_upper
   CALL ELEMENT_EA (XYZ, INOD, NGP, MTEL, PMAT, PML_PARAM, IEL, EA_lower, 111)           ! EA_lower
  
   CALL ELEMENT_EN (XYZ, INOD, NGP, MTEL, PMAT, PML_PARAM, IEL, EN,       300)           ! EN


! Diagonal mass matrix for Explicit time-stepping: HRZ lumping
! if you use spectral elements, mass matrix is diagonal and these lines basically extract the diagonal part
!---------- ---------- ---------- ---------- ----------
! 1. EM
   CALL Extract_Diagonal ( EM, NNODE * NDIM, EM_DIAG )

! 2. EN matrices
   CALL Extract_Diagonal ( EN, NNODE *  3  , EN_DIAG )

! Send Fortran matrix to PETSc data structure
!---------- ---------- ---------- ---------- ---------- 
   EM_PETSC      = EM

   EA_PETSC      = EA_upper
   EA_PETSC_T    = transpose(EA_lower)

   EN_PETSC      = EN

   EM_DIAG_PETSC = EM_DIAG
   EN_DIAG_PETSC = EN_DIAG


!---------- ---------- ---------- ---------- ----------	  
   DO I = 1 , NNODE

      I1 = I
      I2 = I1 + NNODE
      I3 = I2 + NNODE
      
      ND1(I1) = ID(INOD(I,IEL) , 1)
      ND1(I2) = ID(INOD(I,IEL) , 2)
      
      ND2(I1) = ID(INOD(I,IEL) , 3)
      ND2(I2) = ID(INOD(I,IEL) , 4)
      ND2(I3) = ID(INOD(I,IEL) , 5)
     
   END DO
!---------- ---------- ---------- ---------- ----------


! Modify ND according to C indexing
!---------- ---------- ---------- ---------- ----------
   ND1 = ND1 - 1 
   ND2 = ND2 - 1


! ASSEMBLE PML MATRICES
! EM : (U , U); EA : (U , S); EN : (S , S)

! assemble EM
!---------- ---------- ---------- ---------- ----------   
   NEQEL_U = NNODE * NDIM
   CALL MatSetValuesLocal ( AM_PETSC, NEQEL_U, ND1, NEQEL_U, ND1, EM_PETSC, ADD_VALUES, IERR )
!---------- ---------- ---------- ---------- ----------


! assemble EN
!---------- ---------- ---------- ---------- ----------   
   NEQEL_S = NNODE * 3
   CALL MatSetValuesLocal ( AK_PETSC, NEQEL_S, ND2, NEQEL_S, ND2, EN_PETSC, ADD_VALUES, IERR )
!---------- ---------- ---------- ---------- ----------


! assemble EA matrices
! note: check this part later; it is really tricky (and dirty)! 
! Petsc have a row-oriented matrix strategy!
!---------- ---------- ---------- ---------- ----------   
   NEQEL_U = NNODE * NDIM
   NEQEL_S = NNODE * 3
   CALL MatSetValuesLocal ( AK_PETSC, NEQEL_S, ND2, NEQEL_U, ND1, EA_PETSC,   ADD_VALUES, IERR )
! assemble transpose
   CALL MatSetValuesLocal ( AK_PETSC, NEQEL_U, ND1, NEQEL_S, ND2,-EA_PETSC_T, ADD_VALUES, IERR )
!---------- ---------- ---------- ---------- ----------


! Assemble Diagonal mass for Explicit time stepping
!---------- ---------- ---------- ---------- ----------
   CALL VecSetValuesLocal ( DIAG_M_PETSC, NEQEL_U, ND1, EM_DIAG_PETSC, ADD_VALUES, IERR )
   CALL VecSetValuesLocal ( DIAG_M_PETSC, NEQEL_S, ND2, EN_DIAG_PETSC, ADD_VALUES, IERR )
!---------- ---------- ---------- ---------- ----------


END DO
!---------- ---------- ---------- ---------- ----------


RETURN 
END
