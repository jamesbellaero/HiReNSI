!************************************************
! Explicit time stepping of the adjoint problem, based on Discretized-Then-Optimized Runge-Kutta method, 4 stage-4th order, PML-3D
!************************************************
! REVISION : Sunday 14 January 2014: There is much pleasure to be gained from useless knowledge.

SUBROUTINE EXPLICIT_RK4_PETSC_PML_3D_Adjoint ( AK_PETSC, DIAG_M_PETSC, AC_PETSC, AG_PETSC, B_mis_PETSC, AK_RD_PETSC, DIAG_M_RD_PETSC, ID, LTRANS, RANK, &
                                       NNDH, NNVH, NNAH, EqDis, EqVel, EqAcc, idx_from, ID_Para, Param, &
                                       U_Store_Mapping, NStore_Mapping, IS_u_from, IS_u_to, U_Store_Numbers_Global, &
                                       EqDis_MAP_U_Store_Numbers_Global, Dis_meas )

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

Integer, Intent(In), Dimension ( Param%IntM ( 5, 1) )        :: idx_from 
Integer, Intent(In), Dimension ( Param%IntM ( 5, 2), NDim )  :: ID_Para 
Integer, Allocatable, Dimension (:)                          :: idx_to 

Integer, Intent(In), Dimension ( NStore_Mapping )  :: U_Store_Numbers_Global
Integer, Intent(In), Dimension ( NNDH * NDim )     :: EqDis_MAP_U_Store_Numbers_Global

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
Vec            :: B_mis_PETSC
Vec            :: DIAG_M_PETSC, DIAG_M_RD_PETSC

PetscErrorCode :: IERR, ErrPTC
PetscMPIInt    :: SIZE, RANK
PetscScalar    :: ENERGY

PetscScalar    :: fac7, fac8, fac9

! - PETSc local vectors
Vec            :: y_t_PETSC, y_m_PETSC, y_b_PETSC
Vec            :: p1_t_PETSC, p1_m_PETSC, p1_b_PETSC, p2_t_PETSC, p2_m_PETSC, p2_b_PETSC, help_PETSC
Vec            :: p3_t_PETSC, p3_m_PETSC, p3_b_PETSC, p4_t_PETSC, p4_m_PETSC, p4_b_PETSC
Vec            :: HELP_1_PETSC

PetscScalar, pointer    :: U(:),      UD(:),      UDD(:)
PetscScalar, pointer    :: U_Rank(:), UD_Rank(:), UDD_Rank(:) ;

Vec            :: U_Rank_PETSC ;                   ! 
Vec            :: UD_Rank_PETSC ;                  ! 
Vec            :: UDD_Rank_PETSC ;                 ! 


PetscScalar    :: u_misfit_PETSC ( NNDH * NDim )
! - Index setting for visualization -----------------------------------------------------------------------------------------------------------------
VecScatter     :: vscat_U, vscat_UD, vscat_UDD ; ! Used to scatter ghost node solution, this is just for Paraview
IS             :: IS_from, IS_to  ;              ! Index sets for retrieving the ghost node solution

! - Storage -----------------------------------------------------------------------------------------------------------------------------------------
IS             :: IS_u_from, IS_u_to

!==================================================================================================


! in
!---------- ---------- ---------- ---------- ----------
DIMENSION ID ( NJ, NDOF ), LTRANS ( NTRANS )
!---------- ---------- ---------- ---------- ----------


! in
!---------- ---------- ---------- ---------- ----------
Dimension U_Store_Mapping ( NStore_Mapping , 0:NStep )     ! Displacement components of the forward problem
Dimension Dis_meas ( 0:NStep , NNDH * NDim )               ! Measured displacements at select sensor locations
!---------- ---------- ---------- ---------- ----------


! local
!---------- ---------- ---------- ---------- ----------
Dimension u_misfit ( NNDH * NDim )
!---------- ---------- ---------- ---------- ----------

IF (RANK == 0) WRITE(*,*)'3D RK4 Adjoint'












! Kids playground -----------------------------------------------------------------------------------------------------------------------------------
!  *                                                                                                        *
! ---------------------------------------------------------------------------------------------------------------------------------------------------
Call MPI_Barrier ( PETSC_COMM_WORLD, ierr )





IF ( RANK == 0 ) then


!write(*,*) 'rank =', rank
!write(*,*) 'eqdis',eqdis

!write(*,*) 

!write(*,*) U_Store_Numbers_Global ( EqDis_MAP_U_Store_Numbers_Global ( 1 ) )

do I = 1, NStore_Mapping
!  write(*,*) I, U_Store_Numbers_Global(I)
end do


End If

Call MPI_Barrier ( PETSC_COMM_WORLD, ierr )
!stop
! Kids playground -----------------------------------------------------------------------------------------------------------------------------------
!  *                                                                                                        *
! ---------------------------------------------------------------------------------------------------------------------------------------------------











I_PRINT_STEP =  1


! define Index Set for Paraview ---------------------------------------------------------------------------------------------------------------------
Allocate ( idx_to ( Param%IntM ( 5, 1) ) ) ;

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
CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, y_t_PETSC, IERR )
CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, y_m_PETSC, IERR )
CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, y_b_PETSC, IERR )

CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, p1_t_PETSC, IERR )
CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, p1_m_PETSC, IERR )
CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, p1_b_PETSC, IERR )

CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, p2_t_PETSC, IERR )
CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, p2_m_PETSC, IERR )
CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, p2_b_PETSC, IERR )

CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, p3_t_PETSC, IERR )
CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, p3_m_PETSC, IERR )
CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, p3_b_PETSC, IERR )

CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, p4_t_PETSC, IERR )
CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, p4_m_PETSC, IERR )
CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, p4_b_PETSC, IERR )

CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, help_PETSC, IERR )
CALL VecCreateMPI ( PETSC_COMM_WORLD, NEQM_Mapping, NEQMTotal, HELP_1_PETSC, IERR )

CALL VecSet       ( y_t_PETSC, 0.0d0, ierr)
CALL VecSet       ( y_m_PETSC, 0.0d0, ierr)
CALL VecSet       ( y_b_PETSC, 0.0d0, ierr)
!==================================================================================================

Call VecScatterCreate ( y_m_PETSC,   IS_from, U_Rank_PETSC,   IS_to, vscat_U, ErrPTC );
Call VecScatterCreate ( y_b_PETSC,   IS_from, UD_Rank_PETSC,  IS_to, vscat_UD, ErrPTC );
Call VecScatterCreate ( p1_b_PETSC, IS_from, UDD_Rank_PETSC, IS_to, vscat_UDD, ErrPTC );

! inverse of the diagonal mass matrix
!---------- ---------- ---------- ---------- ----------  
! CALL VecReciprocal ( DIAG_M_PETSC, IERR )
! This has already been accomplished when solving the Forward problem.


! form right hand side matrices
!---------- ---------- ---------- ---------- ----------
!CALL MatDiagonalScale ( AK_PETSC, DIAG_M_PETSC, PETSC_NULL_OBJECT, IERR )
!CALL MatDiagonalScale ( AC_PETSC, DIAG_M_PETSC, PETSC_NULL_OBJECT, IERR )
!CALL MatDiagonalScale ( AG_PETSC, DIAG_M_PETSC, PETSC_NULL_OBJECT, IERR )
! This has already been accomplished when solving the Forward problem.

!CALL VecPointwiseMult ( B_PETSC,  B_PETSC, DIAG_M_PETSC   , IERR ) !----------$$$$$$$$$$ take care of this later

!CALL MatScale ( AK_PETSC, -1.0D0, IERR )
!CALL MatScale ( AC_PETSC, -1.0D0, IERR )
!CALL MatScale ( AG_PETSC, -1.0D0, IERR )
! This has already been accomplished when solving the Forward problem.




! print out & Visualization
!---------- ---------- ---------- ---------- ---------- ---------------------------------------------------------------------------------------------
Call VecGetOwnershipRange ( y_t_PETSC, IStart0, IEnd0, ErrPTC ) ;
IStart = IStart0
IEnd   = IEnd0
!---------- ---------- ---------- ---------- ---------- ---------------------------------------------------------------------------------------------


CALL MPI_Comm_size   ( PETSC_COMM_WORLD, SIZE, IERR )
CALL MPI_Comm_rank   ( PETSC_COMM_WORLD, RANK, IERR )

! - Write total wraper file  ------------------------------------------------------------------------------------------------------------------------
!  If ( Rank == 0 ) Then ;
!    UnFile = UN_OutallWrp ;
!    Write (Unit = UnFile, FMT = "(A63)", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) '<Xdmf xmlns:xi="http://www.w3.org/2003/XInclude" Version="2.2">' ;
!    Write (Unit = UnFile, FMT = "(A10)", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) '  <Domain>' ;
!    Write (Unit = UnFile, FMT = "(A57)", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) "    <Grid GridType='Collection' CollectionType='Spatial'>" ;
!  End If ;


! Explicit Integration ------------------------------------------------------------------------------------------------------------------------------
DO ISTEP = NSTEP , 0 , -1

! 1- screen printout for nonbelivers
   T = DBLE( ISTEP ) * DT
   IF ( RANK == 0 .AND. MOD ( ISTEP , I_PRINT_STEP ) == 0 ) THEN
      WRITE (*,'(A30,F10.5)') 'TIME=', T
   END IF

! 2- construction of the misfit at select sensor locations
   Do IEq = 1, NNDH * NDim
     u_misfit ( IEq ) = U_Store_Mapping ( EqDis_MAP_U_Store_Numbers_Global ( IEq ) , IStep ) - Dis_meas ( IStep , IEq )
   End Do

! 3- initialize the misfit vector
   Call VecSet ( B_mis_PETSC, 0.0d0, IERR )

! 4- assemble the distributed misfit vector
   u_misfit_PETSC = u_misfit
   If ( NNDH /= 0 ) CALL VecSetValues ( B_mis_PETSC, NNDH * NDim, EqDis, u_misfit_PETSC, ADD_VALUES, IERR )
   Call VecAssemblyBegin ( B_mis_PETSC, IERR )
   Call VecAssemblyEnd   ( B_mis_PETSC, IERR )
  !Call VecView ( B_mis_PETSC, PETSc_VIEWER_STDOUT_World, IERR )

! 5- Last time-step requires special care
   If ( ISTEP == NSTEP ) Then
     Call VecSet   ( y_m_PETSC, 0.0d0, IERR )
     Call VecAXPY  ( y_m_PETSC, DT, B_mis_PETSC, IERR )
   End If

! 6- middle steps: compute p4_t, p4_m, p4_b
   Call VecSet  ( p4_t_PETSC, 0.0d0, IERR )
   Call VecAXPY ( p4_t_PETSC, DT/6.0d0, y_t_PETSC, IERR )
   Call VecSet  ( p4_m_PETSC, 0.0d0, IERR )
   Call VecAXPY ( p4_m_PETSC, DT/6.0d0, y_m_PETSC, IERR )
   Call VecSet  ( p4_b_PETSC, 0.0d0, IERR )
   Call VecAXPY ( p4_b_PETSC, DT/6.0d0, y_b_PETSC, IERR )

! 7- middle steps: compute p3_t, p3_m, p3_b
   Call VecSet   ( p3_t_PETSC, 0.0d0, IERR )
   Call VecMAXPY ( p3_t_PETSC, 2, [DT, DT/3.0d0], [p4_m_PETSC, y_t_PETSC], IERR )
   Call VecSet   ( p3_m_PETSC, 0.0d0, IERR )
   Call VecMAXPY ( p3_m_PETSC, 2, [DT, DT/3.0d0], [p4_b_PETSC, y_m_PETSC], IERR )
   Call VecSet ( p3_b_PETSC, 0.0d0, IERR )
   Call MatMultTransposeAdd ( AG_PETSC, p4_t_PETSC, p3_b_PETSC, p3_b_PETSC, IERR )
   Call MatMultTransposeAdd ( AK_PETSC, p4_m_PETSC, p3_b_PETSC, p3_b_PETSC, IERR )
   Call MatMultTransposeAdd ( AC_PETSC, p4_b_PETSC, p3_b_PETSC, p3_b_PETSC, IERR )
   Call VecScale ( p3_b_PETSC, DT, IERR )
   Call VecAXPY ( p3_b_PETSC, DT/3.0d0, y_b_PETSC, IERR )

! 8- middle steps: compute p2_t, p2_m, p2_b
   Call VecSet   ( p2_t_PETSC, 0.0d0, IERR )
   Call VecMAXPY ( p2_t_PETSC, 2, [DT/2.0d0, DT/3.0d0], [p3_m_PETSC, y_t_PETSC], IERR )
   Call VecSet   ( p2_m_PETSC, 0.0d0, IERR )
   Call VecMAXPY ( p2_m_PETSC, 2, [DT/2.0d0, DT/3.0d0], [p3_b_PETSC, y_m_PETSC], IERR )
   Call VecSet ( p2_b_PETSC, 0.0d0, IERR )
   Call MatMultTransposeAdd ( AG_PETSC, p3_t_PETSC, p2_b_PETSC, p2_b_PETSC, IERR )
   Call MatMultTransposeAdd ( AK_PETSC, p3_m_PETSC, p2_b_PETSC, p2_b_PETSC, IERR )
   Call MatMultTransposeAdd ( AC_PETSC, p3_b_PETSC, p2_b_PETSC, p2_b_PETSC, IERR )
   Call VecScale ( p2_b_PETSC, DT/2.0d0, IERR )
   Call VecAXPY ( p2_b_PETSC, DT/3.0d0, y_b_PETSC, IERR )

! 9- middle steps: compute p1_t, p1_m, p1_b
   Call VecSet   ( p1_t_PETSC, 0.0d0, IERR )
   Call VecMAXPY ( p1_t_PETSC, 2, [DT/2.0d0, DT/6.0d0], [p2_m_PETSC, y_t_PETSC], IERR )
   Call VecSet   ( p1_m_PETSC, 0.0d0, IERR )
   Call VecMAXPY ( p1_m_PETSC, 2, [DT/2.0d0, DT/6.0d0], [p2_b_PETSC, y_m_PETSC], IERR )
   Call VecSet ( p1_b_PETSC, 0.0d0, IERR )
   Call MatMultTransposeAdd ( AG_PETSC, p2_t_PETSC, p1_b_PETSC, p1_b_PETSC, IERR )
   Call MatMultTransposeAdd ( AK_PETSC, p2_m_PETSC, p1_b_PETSC, p1_b_PETSC, IERR )
   Call MatMultTransposeAdd ( AC_PETSC, p2_b_PETSC, p1_b_PETSC, p1_b_PETSC, IERR )
   Call VecScale ( p1_b_PETSC, DT/2.0d0, IERR )
   Call VecAXPY ( p1_b_PETSC, DT/6.0d0, y_b_PETSC, IERR )

! 10- middle steps: y_t^{n-1}
   ! compute h_m
   Call VecSet  ( help_PETSC, 0.0d0, IERR )
   Call VecMAXPY( help_PETSC, 4, [1.0d0, 1.0d0, 1.0d0, 1.0d0], [p1_m_PETSC, p2_m_PETSC, p3_m_PETSC, p4_m_PETSC], IERR )
   ! compute y_t^{n-1}
   Call VecAXPY ( y_t_PETSC, 1.0d0, help_PETSC, IERR )
   ! update y_b^{n-1}
   Call MatMultTransposeAdd ( AK_PETSC, help_PETSC, y_b_PETSC, y_b_PETSC, IERR )

! 11- middle steps: y_m^{n-1}
   ! compute h_b
   Call VecSet  ( help_PETSC, 0.0d0, IERR )
   Call VecMAXPY( help_PETSC, 4, [1.0d0, 1.0d0, 1.0d0, 1.0d0], [p1_b_PETSC, p2_b_PETSC, p3_b_PETSC, p4_b_PETSC], IERR )
   ! compute y_m^{n-1}
   Call VecMAXPY ( y_m_PETSC, 2, [1.0d0, DT], [help_PETSC, B_mis_PETSC], IERR )
   ! update y_b^{n-1}
   Call MatMultTransposeAdd ( AC_PETSC, help_PETSC, y_b_PETSC, y_b_PETSC, IERR )

! 12- middle steps: y_b^{n-1}
   ! compute h_t
   Call VecSet  ( help_PETSC, 0.0d0, IERR )
   Call VecMAXPY( help_PETSC, 4, [1.0d0, 1.0d0, 1.0d0, 1.0d0], [p1_t_PETSC, p2_t_PETSC, p3_t_PETSC, p4_t_PETSC], IERR )
   ! compute y_b^{n-1}
   Call MatMultTransposeAdd ( AG_PETSC, help_PETSC, y_b_PETSC, y_b_PETSC, IERR )


END DO
! End of Explicit Integration loop ------------------------------------------------------------------------------------------------------------------


CALL VecDestroy ( y_t_PETSC, IERR )
CALL VecDestroy ( y_m_PETSC, IERR )
CALL VecDestroy ( y_b_PETSC, IERR )

CALL VecDestroy ( p1_t_PETSC, IERR )
CALL VecDestroy ( p1_m_PETSC, IERR )
CALL VecDestroy ( p1_b_PETSC, IERR )

CALL VecDestroy ( p2_t_PETSC, IERR )
CALL VecDestroy ( p2_m_PETSC, IERR )
CALL VecDestroy ( p2_b_PETSC, IERR )

CALL VecDestroy ( p3_t_PETSC, IERR )
CALL VecDestroy ( p3_m_PETSC, IERR )
CALL VecDestroy ( p3_b_PETSC, IERR )

CALL VecDestroy ( p4_t_PETSC, IERR )
CALL VecDestroy ( p4_m_PETSC, IERR )
CALL VecDestroy ( p4_b_PETSC, IERR )

CALL VecDestroy ( help_PETSC, IERR )
CALL VecDestroy ( HELP_1_PETSC, IERR )

Call VecDestroy   ( U_Rank_PETSC,   ErrPTC ) ;
Call VecDestroy   ( UD_Rank_PETSC,  ErrPTC ) ;
Call VecDestroy   ( UDD_Rank_PETSC, ErrPTC ) ;

! - Destroying PETSc Index Set ----------------------------------------------------------------------------------------------------------------------
Call ISDestroy(IS_from, ErrPTC );
Call ISDestroy(IS_to,   ErrPTC );

Call VecScatterDestroy(vscat_U, ErrPTC);
Call VecScatterDestroy(vscat_UD, ErrPTC);
Call VecScatterDestroy(vscat_UDD, ErrPTC);

!  If ( Rank == 0 ) Then ;
!    UnFile = UN_OutallWrp ;
!    Close ( Unit = UnFile, STATUS = 'KEEP', ERR = 1002, IOSTAT = IO_File ) ;
!  End If ;

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


