! These subroutines provide input data, pre/post -processing services, etc. for some tailor-made examples, considered for code verification.

!************************************************
! specify active dofs for Lagrange nodes in Pine Flat-3D example
! this is a plane-stress model.
!************************************************
! REVISION : F, 25 May 2012

SUBROUTINE Pine_Flat_3D_serendipity_to_Lagrange ( XYZ, ID )
USE PARAMETERS
IMPLICIT DOUBLE PRECISION (A-H,O-Z)


! in 
!---------- ---------- ---------- ---------- ----------
DIMENSION XYZ(NJ,NDIM), ID(NJ,NDOF)
!---------- ---------- ---------- ---------- ----------


TOL = 1.0D-3

DO I = 1, NJ

   X_1 = XYZ(I,1)
   X_2 = XYZ(I,2)
   X_3 = XYZ(I,3)

! dam base is fixed
   IF ( DABS(X_3) < TOL ) THEN
      ID(I,1) = 1                                ! translation in x-dir is restricted
      ID(I,2) = 0                                ! translation in y-dir is permitted (plane-stress)
      ID(I,3) = 1                                ! translation in z-dir is restricted      
      IF ( DABS(X_2 - 5.0D0) < TOL ) THEN
         ID(I,2) = 1                             ! avoid rigid body motion
      END IF
   END IF

! mid plane does not move along y-dir
   IF ( DABS(X_2 - 5.0D0) < TOL ) THEN
      ID(I,2) = 1
   END IF

END DO


RETURN
END


!************************************************
! ASSEMBLE PML MATRICES
! this is based on Sezgin's formulation, which leads to an implicit scheme
!************************************************
! REVISION : Th 19 July 2012

SUBROUTINE ASSEM_PML_PETSC_3D_Sezgin ( AK_PETSC, AM_PETSC, AC_PETSC, AG_PETSC, XYZ, INOD, NGP, MTEL, PMAT, PML_PARAM, ID_BC, ID, SIZE, RANK )
USE PARAMETERS
IMPLICIT DOUBLE PRECISION (A-H,O-Z)

!==================================================================================================
! - PETSC INCLUDE FILES
#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscmat.h"

! - PETSC VARIABLES AND MATRICES
Mat            :: AK_PETSC, AM_PETSC, AC_PETSC, AG_PETSC

PetscErrorCode :: IERR
PetscMPIInt    :: SIZE, RANK

PetscScalar    :: EM_a_PETSC(NNODE * NDIM, NNODE * NDIM), EM_b_PETSC(NNODE * NDIM, NNODE * NDIM), EM_c_PETSC(NNODE * NDIM, NNODE * NDIM), &
                  EM_d_PETSC(NNODE * NDIM, NNODE * NDIM)

PetscScalar    :: EA_e_PETSC(NNODE * 3, NNODE * 6), EA_p_PETSC(NNODE * 3, NNODE * 6), EA_w_PETSC(NNODE * 3, NNODE * 6)
PetscScalar    :: EN_a_PETSC(NNODE * 6, NNODE * 6), EN_b_PETSC(NNODE * 6, NNODE * 6), EN_c_PETSC(NNODE * 6, NNODE * 6), EN_d_PETSC(NNODE * 6, NNODE * 6)

PetscScalar    :: EA_e_PETSC_T(NNODE * 6, NNODE * 3), EA_p_PETSC_T(NNODE * 6, NNODE * 3), EA_w_PETSC_T(NNODE * 6, NNODE * 3)
!==================================================================================================


! in 
!---------- ---------- ---------- ---------- ----------
DIMENSION XYZ(NJ,NDIM), INOD(NNODE,NEL), NGP(NEL), MTEL(NEL), PMAT(NMAT,NPM), ID(NJ,NDOF), ID_BC(NEL,NDIM**2)
DIMENSION PML_PARAM(NDIM * 2 , 4)
!---------- ---------- ---------- ---------- ----------


! local
!---------- ---------- ---------- ---------- ----------
DIMENSION EM_a(NNODE * NDIM, NNODE * NDIM), EM_b(NNODE * NDIM, NNODE * NDIM), EM_c(NNODE * NDIM, NNODE * NDIM), EM_d(NNODE * NDIM, NNODE * NDIM)
DIMENSION EA_e_upper(NNODE * 3, NNODE * 6), EA_p_upper(NNODE * 3, NNODE * 6), EA_w_upper(NNODE * 3, NNODE * 6)
DIMENSION EA_e_lower(NNODE * 3, NNODE * 6), EA_p_lower(NNODE * 3, NNODE * 6), EA_w_lower(NNODE * 3, NNODE * 6)
DIMENSION EN_a(NNODE * 6, NNODE * 6), EN_b(NNODE * 6, NNODE * 6), EN_c(NNODE * 6, NNODE * 6), EN_d(NNODE * 6, NNODE * 6)

DIMENSION ND1(NNODE * NDIM), ND2(NNODE * 6)
!---------- ---------- ---------- ---------- ----------  
 

DO IEL  = 1 + RANK , NEL , SIZE

   MTYPE = MTEL(IEL)
!  MTYPE 1 : REGULAR DOMAIN
!  MTYPE 2 : PML
   IF ( MTYPE /= 2 ) CYCLE


   WRITE(*,'(A30, I10)')'PML - ELEMENT NO.', IEL
!  WRITE(2,'(A30, I10)')'PML - ELEMENT NO.', IEL


! PML matrices
!---------- ---------- ---------- ---------- ---------- 
   CALL ELEMENT_M_3D  ( XYZ, INOD, NGP, MTEL, PMAT, PML_PARAM, IEL, EM_a,       202 )    ! PML M_a
   CALL ELEMENT_M_3D  ( XYZ, INOD, NGP, MTEL, PMAT, PML_PARAM, IEL, EM_b,       203 )    ! PML M_b
   CALL ELEMENT_M_3D  ( XYZ, INOD, NGP, MTEL, PMAT, PML_PARAM, IEL, EM_c,       204 )    ! PML M_c
   CALL ELEMENT_M_3D  ( XYZ, INOD, NGP, MTEL, PMAT, PML_PARAM, IEL, EM_d,       205 )    ! PML M_d

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
   CALL MatSetValues ( AM_PETSC, NEQEL_U, ND1, NEQEL_U, ND1, EM_a_PETSC, ADD_VALUES, IERR )
   CALL MatSetValues ( AC_PETSC, NEQEL_U, ND1, NEQEL_U, ND1, EM_b_PETSC, ADD_VALUES, IERR )
   CALL MatSetValues ( AK_PETSC, NEQEL_U, ND1, NEQEL_U, ND1, EM_c_PETSC, ADD_VALUES, IERR )
   CALL MatSetValues ( AG_PETSC, NEQEL_U, ND1, NEQEL_U, ND1, EM_d_PETSC, ADD_VALUES, IERR )
!---------- ---------- ---------- ---------- ----------


! assemble EN matrices 
!---------- ---------- ---------- ---------- ----------   
   NEQEL_S = NNODE * 6
   CALL MatSetValues ( AM_PETSC, NEQEL_S, ND2, NEQEL_S, ND2, -EN_a_PETSC, ADD_VALUES, IERR )
   CALL MatSetValues ( AC_PETSC, NEQEL_S, ND2, NEQEL_S, ND2, -EN_b_PETSC, ADD_VALUES, IERR )
   CALL MatSetValues ( AK_PETSC, NEQEL_S, ND2, NEQEL_S, ND2, -EN_c_PETSC, ADD_VALUES, IERR )
   CALL MatSetValues ( AG_PETSC, NEQEL_S, ND2, NEQEL_S, ND2, -EN_d_PETSC, ADD_VALUES, IERR )
!---------- ---------- ---------- ---------- ----------


! assemble EA matrices
! note: check this part later; it is really tricky (and dirty)! 
! Petsc have a row-oriented matrix strategy!
!---------- ---------- ---------- ---------- ----------   
   NEQEL_U = NNODE * NDIM
   NEQEL_S = NNODE * 6
   CALL MatSetValues ( AC_PETSC, NEQEL_s, ND2, NEQEL_u, ND1, EA_e_PETSC, ADD_VALUES, IERR )
   CALL MatSetValues ( AK_PETSC, NEQEL_s, ND2, NEQEL_u, ND1, EA_p_PETSC, ADD_VALUES, IERR )
   CALL MatSetValues ( AG_PETSC, NEQEL_s, ND2, NEQEL_u, ND1, EA_w_PETSC, ADD_VALUES, IERR )
! assemble transpose
   CALL MatSetValues ( AC_PETSC, NEQEL_u, ND1, NEQEL_s, ND2, EA_e_PETSC_T, ADD_VALUES, IERR )
   CALL MatSetValues ( AK_PETSC, NEQEL_u, ND1, NEQEL_s, ND2, EA_p_PETSC_T, ADD_VALUES, IERR )
   CALL MatSetValues ( AG_PETSC, NEQEL_u, ND1, NEQEL_s, ND2, EA_w_PETSC_T, ADD_VALUES, IERR )
!---------- ---------- ---------- ---------- ----------


END DO
!---------- ---------- ---------- ---------- ----------


RETURN 
END


!************************************************
! Assemble the delta source
! Greenberg's delta function (p.12) is used to verify with Kausel's analytical solution (p.78)
!************************************************
! REVISION : F 27 July 2012

SUBROUTINE ASSEM_delta_verification_PETSC_3D ( B_PETSC, XYZ, INOD, NGP, MTEL, PMAT, ID_BC, ID, SIZE, RANK )
USE PARAMETERS
IMPLICIT DOUBLE PRECISION (A-H,O-Z)


!==================================================================================================
! - PETSC INCLUDE FILES
#include "finclude/petscsys.h"
#include "finclude/petscvec.h"

! - PETSC VARIABLES AND MATRICES
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
!DIMENSION EFS(NNODE * NDIM), ND1(NNODE * NDIM)
!---------- ---------- ---------- ---------- ----------


WRITE(*,'(A30)') 'delta function verification'


! center node that the delta source is acting on it
!---------- ---------- ---------- ---------- ----------
!NODE_delta = 3042                               ! IN_Ansys_3D_01
 NODE_delta = 6098                               ! IN_Ansys_3D_02

! Local to Global mapping
!---------- ---------- ---------- ---------- ----------
ND = ID ( NODE_delta, 3)


! Modify ND according to C indexing
!---------- ---------- ---------- ---------- ----------
ND = ND - 1


! Assemble the vertical source
!---------- ---------- ---------- ---------- ----------
CALL VecSetValue ( B_PETSC, ND, 1.0d0, ADD_VALUES, IERR )


RETURN
END


!************************************************
! READS ANSYS.OUT AND PREPARES INPUT-3D
! assumption: coordinate system is located on the surface, in the middle; +z is upward.
! this is a 3D imitation of a 1D bar, truncated by PML at the bottom boundary.
!************************************************
! REVISION : Sat 20 April 2013: after the tea party!

SUBROUTINE ANSYS_3D_Import_1D ( XYZ, ID, INOD, NGP, MTEL, ID_BC, PMAT, PML_PARAM, XYZ_Lagrange, ID_Lagrange )
USE PARAMETERS
IMPLICIT DOUBLE PRECISION (A-H,O-Z)


! out     
!---------- ---------- ---------- ---------- ----------   
DIMENSION XYZ(NJ,NDIM), ID(NJ,NDOF), INOD(NNODE,NEL), NGP(NEL), MTEL(NEL), ID_BC(NEL, NDIM**2), PMAT(NMAT,NPM), PML_PARAM(NDIM * 2 , 4)
DIMENSION XYZ_Lagrange(NEL, 6 * NDIM), ID_Lagrange(NEL, 6)
!---------- ---------- ---------- ---------- ----------


! local   
!---------- ---------- ---------- ---------- ----------   
CHARACTER*100 HELP
!---------- ---------- ---------- ---------- ----------


IF ( NDIM /= 3 ) RETURN
ID_BC      = 0


! number of integration points. for diagonal mass matrix, 9-noded elements should be used with 3 Lobatto points (defined in INPUT1).
!---------- ---------- ---------- ---------- ----------
NGP        = 3


! read geometry
READ (21,'(1X,D20.8)') d_p_Lx
READ (21,'(1X,D20.8)') d_m_Lx
READ (21,'(1X,D20.8)') d_p_Ly
READ (21,'(1X,D20.8)') d_m_Ly
READ (21,'(1X,D20.8)') d_p_Lz
READ (21,'(1X,D20.8)') d_m_Lz
READ (21,'(1X,D20.8)') pml_thickness

WRITE(22,'(F10.2)') d_p_Lx
WRITE(22,'(F10.2)') d_m_Lx
WRITE(22,'(F10.2)') d_p_Ly
WRITE(22,'(F10.2)') d_m_Ly
WRITE(22,'(F10.2)') d_p_Lz
WRITE(22,'(F10.2)') d_m_Lz
WRITE(22,'(F10.2)') pml_thickness


! read number of nodes and elements
READ (21, '(2(1X,I12))') NJ_no_use, NEL_no_use


! read joints
READ (21, '(1X,A5)') HELP
DO I = 1, NJ_serendipity
   READ (21, '(1X,I12,1X,3D20.8)') II, (XYZ(II,J), J = 1, NDIM)
END DO


! read elements (same ordering)
READ (21, '(1X,A5)') HELP
DO I = 1, NEL
   READ (21, '(9(1X,I9))') II,   (INOD(J,I), J = 1, 8)
   READ (21, '(10X, 12(1X,I9))') (INOD(J,I), J = 9, 20)
END DO


! read surface load info
READ (21, '(1X,A5)') HELP
READ (21, '(1X,I12)') NLEL

DO I = 1, NLEL
   READ (21, '(2(1X,I12))') IEL , ID_BC(IEL,1)
END DO


! complete data:
! 0. material properties
!---------- ---------- ---------- ---------- ----------  
PMAT(1,1) = 1.250d3
PMAT(1,2) = 0.250d0
PMAT(1,3) = 2000.0d-6
PMAT(1,4) = 0.0d0
PMAT(1,5) = 1.0d0

PMAT(2,1) = 1.250d3
PMAT(2,2) = 0.250d0
PMAT(2,3) = 2000.0d-6
PMAT(2,4) = 0.0d0
PMAT(2,5) = 1.0d0


! Chopra's explicit example
PMAT(1,1) = 2.50d0
PMAT(1,2) = 0.250d0
PMAT(1,3) = 1.0d0

PMAT(2,1) = 2.50d0
PMAT(2,2) = 0.250d0
PMAT(2,3) = 1.0d0

! 1. complete PML data
!---------- ---------- ---------- ---------- ----------
! PML length
PML_PARAM( 1 , 3 ) = 0.0d0
PML_PARAM( 2 , 3 ) = 0.0d0
PML_PARAM( 3 , 3 ) = 0.0d0
PML_PARAM( 4 , 3 ) = 0.0d0
PML_PARAM( 5 , 3 ) = 0.0d0
PML_PARAM( 6 , 3 ) = pml_thickness

! PML starting point location
PML_PARAM( 1 , 4 ) = d_p_Lx
PML_PARAM( 2 , 4 ) = d_m_Lx
PML_PARAM( 3 , 4 ) = d_p_Ly
PML_PARAM( 4 , 4 ) = d_m_Ly
PML_PARAM( 5 , 4 ) = 0.0d0
PML_PARAM( 6 , 4 ) = d_m_Lz + pml_thickness


! 2. update geometry for Lagrange elements
!---------- ---------- ---------- ---------- ----------
IF ( NDIM == 3 .AND. NNODE == 27 ) THEN
   CALL GEOMETRY_3D_20_to_27 ( XYZ, ID, INOD ,XYZ_Lagrange, ID_Lagrange )
END IF


RETURN
END


!************************************************
! Build ID array: identify constraints 
! this is a 3D imitation of a 1D bar, truncated by PML at the bottom boundary.
!************************************************
! REVISION : Sat 20 April 2013: after the tea party!

SUBROUTINE ANSYS_3D_Process_1D ( XYZ, ID, INOD, NGP, MTEL, ID_BC, PMAT, PML_PARAM )
USE PARAMETERS
IMPLICIT DOUBLE PRECISION (A-H,O-Z)


! in
!---------- ---------- ---------- ---------- ----------   
DIMENSION XYZ(NJ,NDIM), INOD(NNODE,NEL), NGP(NEL), MTEL(NEL), ID_BC(NEL, NDIM**2), PMAT(NMAT,NPM), PML_PARAM(NDIM * 2 , 4)
!---------- ---------- ---------- ---------- ----------


! out
!---------- ---------- ---------- ---------- ----------
DIMENSION ID(NJ,NDOF)
!---------- ---------- ---------- ---------- ----------


! local: for Lagrange elements, the name is clear; for serendipity elements, it denotes the center of the element (i.e. the node does not exist).
!---------- ---------- ---------- ---------- ----------
DIMENSION X27(NDIM)
DIMENSION XT(NDIM, NNODE_serendipity), FN(NNODE_serendipity), DFXI(NNODE_serendipity, NDIM), DJ(NDIM, NDIM)
!---------- ---------- ---------- ---------- ----------


IF ( NDIM /= 3 ) RETURN
TOL = 1.0d-6


SELECT CASE (NDOF)


   CASE (3)                                      ! elastodynamic medium
!---------- ---------- ---------- ---------- ----------
! remove PML
      pml_thickness = PML_PARAM( 6 , 3 )
! edges of the computational domain
      d_p_Lx = PML_PARAM(1,4)
      d_m_Lx = PML_PARAM(2,4)
      d_p_Ly = PML_PARAM(3,4)
      d_m_Ly = PML_PARAM(4,4)
      d_p_Lz = PML_PARAM(5,4)
      d_m_Lz = PML_PARAM(6,4) - pml_thickness

      MTEL = 1
      ID   = 0

      DO I = 1, NJ

         X_1 = XYZ(I,1)
         X_2 = XYZ(I,2)
         X_3 = XYZ(I,3)

! base is fixed
         IF ( DABS(X_3 - d_m_Lz) < TOL ) THEN
            ID(I,1) = 1                          ! translation in x-dir is restricted
            ID(I,2) = 1                          ! translation in y-dir is restricted
            ID(I,3) = 1                          ! translation in z-dir is restricted   
         END IF

! x = d_p_Lx is fixed (domain plus Lx)
         IF ( DABS(X_1 - d_p_Lx) < TOL ) THEN
!            ID(I,1) = 1                          ! translation in x-dir is restricted
!           ID(I,2) = 1                          ! translation in y-dir is restricted
!           ID(I,3) = 1                          ! translation in z-dir is restricted   
         END IF

! x = d_m_Lx is fixed
         IF ( DABS(X_1 - d_m_Lx) < TOL ) THEN
!            ID(I,1) = 1                          ! translation in x-dir is restricted
!           ID(I,2) = 1                          ! translation in y-dir is restricted
!           ID(I,3) = 1                          ! translation in z-dir is restricted   
         END IF

! y = d_p_Ly is fixed
         IF ( DABS(X_2 - d_p_Ly) < TOL ) THEN
!           ID(I,1) = 1                          ! translation in x-dir is restricted
!            ID(I,2) = 1                          ! translation in y-dir is restricted
!           ID(I,3) = 1                          ! translation in z-dir is restricted   
         END IF

! y = d_m_Ly is fixed
         IF ( DABS(X_2 - d_m_Ly) < TOL ) THEN
!           ID(I,1) = 1                          ! translation in x-dir is restricted
!            ID(I,2) = 1                          ! translation in y-dir is restricted
!           ID(I,3) = 1                          ! translation in z-dir is restricted   
         END IF

      END DO


   CASE (9)                                      ! elastic medium, surrounded by PML
!---------- ---------- ---------- ---------- ---------- ---------- ----------
! 1. identify element type (RD vs PML)
      MTEL = 0

! PML starting point locations
      PML_X01 = PML_PARAM(1,4)
      PML_X02 = PML_PARAM(2,4)
      PML_X03 = PML_PARAM(3,4)
      PML_X04 = PML_PARAM(4,4)
      PML_X05 = PML_PARAM(5,4)
      PML_X06 = PML_PARAM(6,4)

! edges of the computational domain
      d_p_Lx = PML_PARAM(1,4)
      d_m_Lx = PML_PARAM(2,4)
      d_p_Ly = PML_PARAM(3,4)
      d_m_Ly = PML_PARAM(4,4)
      d_p_Lz = PML_PARAM(5,4)
      d_m_Lz = PML_PARAM(6,4) - PML_PARAM(6,3)


      DO IEL = 1, NEL

         SELECT CASE (NNODE)

            CASE (20)                            ! compute the center of the serendipity element

               DO I = 1, NNODE_serendipity
                  K = INOD(I,IEL)
                  DO J = 1, NDIM
                     XT(J,I) = XYZ(K,J)
                  END DO
               END DO

               X1 =  0.0D0 ; X2 =  0.0D0 ; X3 =  0.0D0

               CALL SHAPE_3_20 ( X1, X2, X3, FN, DFXI )

               DJ = MATMUL(XT,DFXI)

               DETJ = DJ(1,1) * DJ(2,2) * DJ(3,3) + DJ(2,1) * DJ(3,2) * DJ(1,3) + DJ(3,1) * DJ(1,2) * DJ(2,3) - &
                      DJ(1,3) * DJ(2,2) * DJ(3,1) - DJ(2,3) * DJ(3,2) * DJ(1,1) - DJ(3,3) * DJ(1,2) * DJ(2,1)

               IF (DETJ.LE.0) THEN
                  WRITE(*,*)'|J|<0 NEGATIVE DET. JACOBIAN'
                  READ(*,*)
                  STOP
               END IF

               X27 = MATMUL ( XT, FN )

            CASE (27)                            ! load the center node (27) of Lagrange element

               K = INOD(27,IEL)
               DO J = 1, NDIM
                  X27(J) = XYZ(K,J)
               END DO

            CASE DEFAULT
               WRITE(*,'("FLAG ERROR IN ANSYS_3D_Process")')
               STOP

         END SELECT

! check if the node resides in Omega_RD or Omega_PML
         IF ( X27(3) <= PML_X06 )  THEN
            MTEL(IEL) = 2
         ELSE
            MTEL(IEL) = 1
         END IF

      END DO


! 2. a. identify the degrees of freedom for each node and set them free (ID that contains all potentially active dofs)
      ID = 1

      DO IEL = 1 , NEL
         MTYPE = MTEL(IEL)

         SELECT CASE (MTYPE)

            CASE (1)                             ! regular domain
               DO I = 1 , NNODE
                  DO J = 1 , NDIM
                     ID(INOD(I, IEL), J) = 0
                  END DO
               END DO

            CASE (2)                             ! PML
               DO I = 1 , NNODE
                  DO J = 1 , NDOF
                     ID(INOD(I, IEL), J) = 0
                  END DO
               END DO

            CASE DEFAULT
               WRITE(*,'("FLAG ERROR IN ANSYS_3D_Process")')
               STOP

         END SELECT

      END DO

! 2. b. apply displacement boundary conditions
      DO I = 1, NJ

         X_1 = XYZ(I,1)
         X_2 = XYZ(I,2)
         X_3 = XYZ(I,3)

! base is fixed
         IF ( DABS(X_3 - d_m_Lz) < TOL ) THEN
            ID(I,1) = 1                          ! translation in x-dir is restricted
            ID(I,2) = 1                          ! translation in y-dir is restricted
            ID(I,3) = 1                          ! translation in z-dir is restricted   
         END IF

! x = d_p_Lx is fixed
         IF ( DABS(X_1 - d_p_Lx) < TOL ) THEN
!            ID(I,1) = 1                          ! translation in x-dir is restricted
!           ID(I,2) = 1                          ! translation in y-dir is restricted
!           ID(I,3) = 1                          ! translation in z-dir is restricted   
         END IF

! x = d_m_Lx is fixed
         IF ( DABS(X_1 - d_m_Lx) < TOL ) THEN
!            ID(I,1) = 1                          ! translation in x-dir is restricted
!           ID(I,2) = 1                          ! translation in y-dir is restricted
!           ID(I,3) = 1                          ! translation in z-dir is restricted   
         END IF

! y = d_p_Ly is fixed
         IF ( DABS(X_2 - d_p_Ly) < TOL ) THEN
!           ID(I,1) = 1                          ! translation in x-dir is restricted
!            ID(I,2) = 1                          ! translation in y-dir is restricted
!           ID(I,3) = 1                          ! translation in z-dir is restricted   
         END IF

! y = d_m_Ly is fixed
         IF ( DABS(X_2 - d_m_Ly) < TOL ) THEN
!           ID(I,1) = 1                          ! translation in x-dir is restricted
!            ID(I,2) = 1                          ! translation in y-dir is restricted
!           ID(I,3) = 1                          ! translation in z-dir is restricted   
         END IF

      END DO


   END SELECT


RETURN
END


!************************************************
! READS ANSYS.OUT AND PREPARES INPUT-3D
! assumption: coordinate system is located on the surface, in the middle; +z is upward.
! this is a 3D imitation of a 2D problem in the XZ plane, truncated by PML at surrounding boundaries.
!************************************************
! REVISION : T 23 April 2013

SUBROUTINE ANSYS_3D_Import_2D ( XYZ, ID, INOD, NGP, MTEL, ID_BC, PMAT, PML_PARAM, XYZ_Lagrange, ID_Lagrange )
USE PARAMETERS
IMPLICIT DOUBLE PRECISION (A-H,O-Z)


! out     
!---------- ---------- ---------- ---------- ----------   
DIMENSION XYZ(NJ,NDIM), ID(NJ,NDOF), INOD(NNODE,NEL), NGP(NEL), MTEL(NEL), ID_BC(NEL, NDIM**2), PMAT(NMAT,NPM), PML_PARAM(NDIM * 2 , 4)
DIMENSION XYZ_Lagrange(NEL, 6 * NDIM), ID_Lagrange(NEL, 6)
!---------- ---------- ---------- ---------- ----------


! local   
!---------- ---------- ---------- ---------- ----------   
CHARACTER*100 HELP
!---------- ---------- ---------- ---------- ----------


IF ( NDIM /= 3 ) RETURN
ID_BC      = 0


! number of integration points. for diagonal mass matrix, 9-noded elements should be used with 3 Lobatto points (defined in INPUT1).
!---------- ---------- ---------- ---------- ----------
NGP        = 3


! read geometry
READ (21,'(1X,D20.8)') d_p_Lx
READ (21,'(1X,D20.8)') d_m_Lx
READ (21,'(1X,D20.8)') d_p_Ly
READ (21,'(1X,D20.8)') d_m_Ly
READ (21,'(1X,D20.8)') d_p_Lz
READ (21,'(1X,D20.8)') d_m_Lz
READ (21,'(1X,D20.8)') pml_thickness

WRITE(22,'(F10.2)') d_p_Lx
WRITE(22,'(F10.2)') d_m_Lx
WRITE(22,'(F10.2)') d_p_Ly
WRITE(22,'(F10.2)') d_m_Ly
WRITE(22,'(F10.2)') d_p_Lz
WRITE(22,'(F10.2)') d_m_Lz
WRITE(22,'(F10.2)') pml_thickness


! read number of nodes and elements
READ (21, '(2(1X,I12))') NJ_no_use, NEL_no_use


! read joints
READ (21, '(1X,A5)') HELP
DO I = 1, NJ_serendipity
   READ (21, '(1X,I12,1X,3D20.8)') II, (XYZ(II,J), J = 1, NDIM)
END DO


! read elements (same ordering)
READ (21, '(1X,A5)') HELP
DO I = 1, NEL
   READ (21, '(9(1X,I9))') II,   (INOD(J,I), J = 1, 8)
   READ (21, '(10X, 12(1X,I9))') (INOD(J,I), J = 9, 20)
END DO


! read surface load info
READ (21, '(1X,A5)') HELP
READ (21, '(1X,I12)') NLEL

DO I = 1, NLEL
   READ (21, '(2(1X,I12))') IEL , ID_BC(IEL,1)
END DO


! complete data:
! 0. material properties
!---------- ---------- ---------- ---------- ----------  
PMAT(1,1) = 1.250d3
PMAT(1,2) = 0.250d0
PMAT(1,3) = 2000.0d-6
PMAT(1,4) = 0.0d0
PMAT(1,5) = 1.0d0

PMAT(2,1) = 1.250d3
PMAT(2,2) = 0.250d0
PMAT(2,3) = 2000.0d-6
PMAT(2,4) = 0.0d0
PMAT(2,5) = 1.0d0


! Chopra's explicit example
PMAT(1,1) = 2.50d0
PMAT(1,2) = 0.250d0
PMAT(1,3) = 1.0d0

PMAT(2,1) = 2.50d0
PMAT(2,2) = 0.250d0
PMAT(2,3) = 1.0d0

! 1. complete PML data
!---------- ---------- ---------- ---------- ----------
! PML length
PML_PARAM( 1 , 3 ) = pml_thickness
PML_PARAM( 2 , 3 ) = pml_thickness
PML_PARAM( 3 , 3 ) = 0.0d0
PML_PARAM( 4 , 3 ) = 0.0d0
PML_PARAM( 5 , 3 ) = 0.0d0
PML_PARAM( 6 , 3 ) = pml_thickness

! PML starting point location
PML_PARAM( 1 , 4 ) = d_p_Lx - pml_thickness
PML_PARAM( 2 , 4 ) = d_m_Lx + pml_thickness
PML_PARAM( 3 , 4 ) = d_p_Ly
PML_PARAM( 4 , 4 ) = d_m_Ly
PML_PARAM( 5 , 4 ) = 0.0d0
PML_PARAM( 6 , 4 ) = d_m_Lz + pml_thickness


! 2. update geometry for Lagrange elements
!---------- ---------- ---------- ---------- ----------
IF ( NDIM == 3 .AND. NNODE == 27 ) THEN
   CALL GEOMETRY_3D_20_to_27 ( XYZ, ID, INOD ,XYZ_Lagrange, ID_Lagrange )
END IF


RETURN
END


!************************************************
! Build ID array: identify constraints 
! this is a 3D imitation of a 2D problem in the XZ plane, truncated by PML at surrounding boundaries.
!************************************************
! REVISION : T 23 April 2013

SUBROUTINE ANSYS_3D_Process_2D ( XYZ, ID, INOD, NGP, MTEL, ID_BC, PMAT, PML_PARAM )
USE PARAMETERS
IMPLICIT DOUBLE PRECISION (A-H,O-Z)


! in
!---------- ---------- ---------- ---------- ----------   
DIMENSION XYZ(NJ,NDIM), INOD(NNODE,NEL), NGP(NEL), MTEL(NEL), ID_BC(NEL, NDIM**2), PMAT(NMAT,NPM), PML_PARAM(NDIM * 2 , 4)
!---------- ---------- ---------- ---------- ----------


! out
!---------- ---------- ---------- ---------- ----------
DIMENSION ID(NJ,NDOF)
!---------- ---------- ---------- ---------- ----------


! local: for Lagrange elements, the name is clear; for serendipity elements, it denotes the center of the element (i.e. the node does not exist).
!---------- ---------- ---------- ---------- ----------
DIMENSION X27(NDIM)
DIMENSION XT(NDIM, NNODE_serendipity), FN(NNODE_serendipity), DFXI(NNODE_serendipity, NDIM), DJ(NDIM, NDIM)
!---------- ---------- ---------- ---------- ----------


IF ( NDIM /= 3 ) RETURN
TOL = 1.0d-6


SELECT CASE (NDOF)


   CASE (3)                                      ! elastodynamic medium
!---------- ---------- ---------- ---------- ----------
      write(*,*) 'incomplete: ANSYS_3D_Process_2D'
      stop
! remove PML
      pml_thickness = PML_PARAM( 6 , 3 )
! edges of the computational domain
      d_p_Lx = PML_PARAM(1,4)
      d_m_Lx = PML_PARAM(2,4)
      d_p_Ly = PML_PARAM(3,4)
      d_m_Ly = PML_PARAM(4,4)
      d_p_Lz = PML_PARAM(5,4)
      d_m_Lz = PML_PARAM(6,4) - pml_thickness

      MTEL = 1
      ID   = 0
      DO I = 1, NJ

         X_1 = XYZ(I,1)
         X_2 = XYZ(I,2)
         X_3 = XYZ(I,3)

! base is fixed
         IF ( DABS(X_3 - d_m_Lz) < TOL ) THEN
            ID(I,1) = 1                          ! translation in x-dir is restricted
            ID(I,2) = 1                          ! translation in y-dir is restricted
            ID(I,3) = 1                          ! translation in z-dir is restricted   
         END IF

! x = d_p_Lx is fixed (domain plus Lx)
         IF ( DABS(X_1 - d_p_Lx) < TOL ) THEN
!            ID(I,1) = 1                          ! translation in x-dir is restricted
!           ID(I,2) = 1                          ! translation in y-dir is restricted
!           ID(I,3) = 1                          ! translation in z-dir is restricted   
         END IF

! x = d_m_Lx is fixed
         IF ( DABS(X_1 - d_m_Lx) < TOL ) THEN
!            ID(I,1) = 1                          ! translation in x-dir is restricted
!           ID(I,2) = 1                          ! translation in y-dir is restricted
!           ID(I,3) = 1                          ! translation in z-dir is restricted   
         END IF

! y = d_p_Ly is fixed
         IF ( DABS(X_2 - d_p_Ly) < TOL ) THEN
!           ID(I,1) = 1                          ! translation in x-dir is restricted
!            ID(I,2) = 1                          ! translation in y-dir is restricted
!           ID(I,3) = 1                          ! translation in z-dir is restricted   
         END IF

! y = d_m_Ly is fixed
         IF ( DABS(X_2 - d_m_Ly) < TOL ) THEN
!           ID(I,1) = 1                          ! translation in x-dir is restricted
!            ID(I,2) = 1                          ! translation in y-dir is restricted
!           ID(I,3) = 1                          ! translation in z-dir is restricted   
         END IF

      END DO


   CASE (9)                                      ! elastic medium, surrounded by PML
!---------- ---------- ---------- ---------- ---------- ---------- ----------
! 1. identify element type (RD vs PML)
      MTEL = 0

! PML starting point locations
      PML_X01 = PML_PARAM(1,4)
      PML_X02 = PML_PARAM(2,4)
      PML_X03 = PML_PARAM(3,4)
      PML_X04 = PML_PARAM(4,4)
      PML_X05 = PML_PARAM(5,4)
      PML_X06 = PML_PARAM(6,4)

! edges of the computational domain
      d_p_Lx = PML_PARAM(1,4) + PML_PARAM(1,3)
      d_m_Lx = PML_PARAM(2,4) - PML_PARAM(2,3)
      d_p_Ly = PML_PARAM(3,4)
      d_m_Ly = PML_PARAM(4,4)
      d_p_Lz = PML_PARAM(5,4)
      d_m_Lz = PML_PARAM(6,4) - PML_PARAM(6,3)


      DO IEL = 1, NEL

         SELECT CASE (NNODE)

            CASE (20)                            ! compute the center of the serendipity element

               DO I = 1, NNODE_serendipity
                  K = INOD(I,IEL)
                  DO J = 1, NDIM
                     XT(J,I) = XYZ(K,J)
                  END DO
               END DO

               X1 =  0.0D0 ; X2 =  0.0D0 ; X3 =  0.0D0

               CALL SHAPE_3_20 ( X1, X2, X3, FN, DFXI )

               DJ = MATMUL(XT,DFXI)

               DETJ = DJ(1,1) * DJ(2,2) * DJ(3,3) + DJ(2,1) * DJ(3,2) * DJ(1,3) + DJ(3,1) * DJ(1,2) * DJ(2,3) - &
                      DJ(1,3) * DJ(2,2) * DJ(3,1) - DJ(2,3) * DJ(3,2) * DJ(1,1) - DJ(3,3) * DJ(1,2) * DJ(2,1)

               IF (DETJ.LE.0) THEN
                  WRITE(*,*)'|J|<0 NEGATIVE DET. JACOBIAN'
                  READ(*,*)
                  STOP
               END IF

               X27 = MATMUL ( XT, FN )

            CASE (27)                            ! load the center node (27) of Lagrange element

               K = INOD(27,IEL)
               DO J = 1, NDIM
                  X27(J) = XYZ(K,J)
               END DO

            CASE DEFAULT
               WRITE(*,'("FLAG ERROR IN ANSYS_3D_Process")')
               STOP

         END SELECT

! check if the node resides in Omega_RD or Omega_PML
         IF ( X27(1) >= PML_X01 .OR. X27(1) <= PML_X02 .OR. X27(3) <= PML_X06 )  THEN
            MTEL(IEL) = 2
         ELSE
            MTEL(IEL) = 1
         END IF

      END DO


! 2. a. identify the degrees of freedom for each node and set them free (ID that contains all potentially active dofs)
      ID = 1

      DO IEL = 1 , NEL
         MTYPE = MTEL(IEL)

         SELECT CASE (MTYPE)

            CASE (1)                             ! regular domain
               DO I = 1 , NNODE
                  DO J = 1 , NDIM
                     ID(INOD(I, IEL), J) = 0
                  END DO
               END DO

            CASE (2)                             ! PML
               DO I = 1 , NNODE
                  DO J = 1 , NDOF
                     ID(INOD(I, IEL), J) = 0
                  END DO
               END DO

            CASE DEFAULT
               WRITE(*,'("FLAG ERROR IN ANSYS_3D_Process")')
               STOP
                                           
         END SELECT

      END DO

! 2. b. apply displacement boundary conditions
      DO I = 1, NJ

         X_1 = XYZ(I,1)
         X_2 = XYZ(I,2)
         X_3 = XYZ(I,3)

! base is fixed
         IF ( DABS(X_3 - d_m_Lz) < TOL ) THEN
            ID(I,1) = 1                          ! translation in x-dir is restricted
            ID(I,2) = 1                          ! translation in y-dir is restricted
            ID(I,3) = 1                          ! translation in z-dir is restricted   
         END IF

! x = d_p_Lx is fixed
         IF ( DABS(X_1 - d_p_Lx) < TOL ) THEN
            ID(I,1) = 1                          ! translation in x-dir is restricted
            ID(I,2) = 1                          ! translation in y-dir is restricted
            ID(I,3) = 1                          ! translation in z-dir is restricted   
         END IF

! x = d_m_Lx is fixed
         IF ( DABS(X_1 - d_m_Lx) < TOL ) THEN
            ID(I,1) = 1                          ! translation in x-dir is restricted
            ID(I,2) = 1                          ! translation in y-dir is restricted
            ID(I,3) = 1                          ! translation in z-dir is restricted   
         END IF

! y = d_p_Ly is fixed
         IF ( DABS(X_2 - d_p_Ly) < TOL ) THEN
!           ID(I,1) = 1                          ! translation in x-dir is restricted
            ID(I,2) = 1                          ! translation in y-dir is restricted
!           ID(I,3) = 1                          ! translation in z-dir is restricted   
         END IF

! y = d_m_Ly is fixed
         IF ( DABS(X_2 - d_m_Ly) < TOL ) THEN
!           ID(I,1) = 1                          ! translation in x-dir is restricted
            ID(I,2) = 1                          ! translation in y-dir is restricted
!           ID(I,3) = 1                          ! translation in z-dir is restricted   
         END IF

      END DO


   END SELECT


RETURN
END
