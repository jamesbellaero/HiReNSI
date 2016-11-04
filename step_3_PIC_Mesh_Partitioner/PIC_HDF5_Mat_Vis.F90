
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                                                                                                                  ++
! Start:        6 April 2014                                                                                                                       ++
! Last Update:  15 April 2014                                                                                                                      ++
! Developed by: Babak Poursartip, AF                                                                                                               ++
!                                                                                                                                                  ++
! Description: THIS Module generates Material Visualization data structure for Paraview based on HDF5 format.                                      ++
!                                                                                                                                                  ++
!              Subroutines                     Description                                                                                         ++
!              ==============                  ==============                                                                                      ++
!              Output                                                                                                                              ++
!                                                                                                                                                  ++
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Module Mat_Vis_HDF5 ;

Use HDF5 ;
Use Parameters ;


Implicit None ;

  Interface
!    Module Procedure 
  End Interface    ;

Contains ;

! ***************************************************************************************************************************************************
!                                                                                                                                                  **
! Start:        6 April 2014                                                                                                                       **
! Last Update:  15 April 2014                                                                                                                      **
! Description: Material Visualization                                                                                                              **
! Called by:                                                                                                                                       **
!                                                                                                                                                  **
! ***************************************************************************************************************************************************

Subroutine OUTPUT_Mat_Vis_HDF5    (                                                                                             &
NDim, MaxNNode, NDOF,                                                                                                           & ! Integer (1) Variables
!                                                                                                                               & ! Integer (2) Variables
NParts,                                                                                                                         & ! Integer (4) Variables
NEL, NJ,                                                                                                                        & ! Integer (8) Variables
!                                                                                                                               & ! Real Variables
EPart, ELT,                                                                                                                     &
NNodeRank,                                                                                                                      &
NEL_Rank, NJ_Rank, Global_PETSc_Num, Local_PETSc_Num, INod,                                                                     & ! Integer Arrays
!                                                                                                                               & ! Real Arrays
ModelName, OutDir                                                                                                               & ! Characters
!                                                                                                                               & ! Type
) ;


Implicit None ;
#include "finclude/petscsys.h"
!#include "metis.h"

! =========================== Global Variables ======================================================================================================

! - Integer Variables -------------------------------------------------------------------------------------------------------------------------------
Integer (Kind=Tiny), Intent(In)    :: NDim, MaxNNode, NDOF ;
Integer (Kind=Shrt), Intent(In)    :: NParts ;
Integer (Kind=Lng ), Intent(In)    :: NEL, NJ;

! - Real Variables ----------------------------------------------------------------------------------------------------------------------------------
! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------

! - Integer Arrays ----------------------------------------------------------------------------------------------------------------------------------
Integer (Kind=Smll), Intent(In),    Dimension (:  )  :: ELT ;

Integer (Kind=Shrt), Intent(In),    Dimension (:  )  :: EPart ;

Integer (Kind=Lng ), Intent(In),    Dimension (:  )  :: NNodeRank ;
Integer (Kind=Lng ), Intent(In),    Dimension (:  )  :: NEL_Rank, NJ_Rank, Global_PETSc_Num ;
Integer (Kind=Lng ), Intent(In),    Dimension (:,:)  :: Local_PETSc_Num, INod ;

! - Real Arrays -------------------------------------------------------------------------------------------------------------------------------------

! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------

! - Character Variables -----------------------------------------------------------------------------------------------------------------------------
Character (Kind = 1, Len = 30 ) :: ModelName ;   ! Name of Input file
Character (Kind = 1, Len = 200) :: OutDir ;      ! Directory of output files (Results)

! - Logical Variables -------------------------------------------------------------------------------------------------------------------------------
! - Type DECLERATIONS -------------------------------------------------------------------------------------------------------------------------------
Type ( NodeID )      :: Nodes ;
Type ( BasicParam )  :: Param ;                  ! Holds basic parameters of each load case

! =========================== LOCAL Variables =======================================================================================================
PetscErrorCode :: ErrPTC ;                       ! Error 

! - Integer Variables -------------------------------------------------------------------------------------------------------------------------------
Integer (Kind=Tiny)  :: NNode ;                  ! Number of Nodes
Integer (Kind=Tiny)  :: X_CellType ;             ! XMDF cell type

Integer (Kind=Smll)  :: UnFile ;                 ! Unit Number of the Output file.
Integer (Kind=Smll)  :: IO_File ;                ! Used for IOSTAT - Input Output Status - in the Open cammand.
Integer (Kind=Smll)  :: IO_Write ;               ! Used for IOSTAT - Input Output Status - in the Write cammand.
Integer (Kind=Smll)  :: IDIM ;                   ! Loop index for NDim.
Integer (Kind=Smll)  :: IDOF ;                   ! Loop index for NDOF.
Integer (Kind=Smll)  :: INode ;                  ! Loop index for NNode.
Integer (Kind=Smll)  :: EType ;                  ! Element Type
Integer (Kind=Smll)  :: ERR_Alloc, ERR_DeAlloc ; ! Allocating and DeAllocating errors

Integer (Kind=Shrt)  :: IParts ;                 ! Loop index on the number of Partitions.
Integer (Kind=Shrt)  :: LocalElN ;               ! counter on local element number

Integer (Kind=Lng )  :: IEL ;                    ! Loop index on NEL.
Integer (Kind=Lng )  :: IJ ;                     ! Loop index on NJ.
Integer (Kind=Lng )  :: I, J, K, I1, I2 ;        ! Loop indeces.
Integer (Kind=Lng )  :: NEqM_Mapping ;           ! Number of Modified EQuations for Mapping application numbering to PETSc numbering. !!?? at the end of the day see if we need this variable.
Integer (Kind=Lng )  :: LEqN ;                   ! Local Equation number, for ID_Local.
Integer (Kind=Lng )  :: NEQM ;                   ! Number of Modified Equations on this rank
Integer (Kind=Lng )  :: Counter ;                ! Counter
Integer (Kind=Lng )  :: PETScNodeNumber ;        ! A temporary varibale for saving node number based on the PETSc node numbering.
Integer (Kind=Lng )  :: Node ;                   ! Holds Node number
Integer (Kind=Lng )  :: NJ_Para ;                ! number joints in the model after downsampling, for Paraview
Integer (Kind=Lng )  :: NE_Mat_Vis_Para ;        ! number equations in the model after downsampling, for Paraview
Integer (Kind=Lng )  :: ConnSizePara ;           ! Size of connectivity vector for Paraview

Integer(HID_T)       :: id_Data_Mat_Vis_PV ;     ! File identifier for recording data for Paraview (material visualization)
Integer(HID_T)       :: dspace_id_IS_Mat_Vis_PV; ! Dataspace identifier for the index set, Paraview purpose only
Integer(HID_T)       ::   dset_id_IS_Mat_Vis_PV; ! Dataset identifier for the index set, Paraview purposes only 

Integer(HSIZE_T), DIMENSION(2) :: dims           ! Dataset dimensions
Integer(HSIZE_T), DIMENSION(2) :: data_dims

! - Real Variables ----------------------------------------------------------------------------------------------------------------------------------
! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
! - Integer Arrays ----------------------------------------------------------------------------------------------------------------------------------
Integer (Kind=Shrt)  :: Local_ElN ( NEL ) ;                          ! Local Element Number
Integer (Kind=Shrt), Allocatable, Dimension(:  )  :: dset_data_ISet; ! Vector for Index Setting 
Integer (Kind=Shrt), Allocatable, Dimension(:  )  :: dset_data_int ; ! Array to fill HDF5 file


Integer (Kind=Lng ), Allocatable, Dimension(:  )  :: ParaNodeNum ;   ! Temporary array used to determine the numnber of nodes in the model for Paraview, note that Paraview does not support spectral elements and we have to transform Spectral elements to Finite elements.

! - Real Arrays -------------------------------------------------------------------------------------------------------------------------------------
!Real (Kind=DBL)   , Allocatable, Dimension(:,:)   :: dset_data_real ;! Array to fill HDF5 file  

! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
! - Character Variables -----------------------------------------------------------------------------------------------------------------------------
Character (Kind = 1, Len = 20)  :: IndexRank ;   ! A variable for storing the Rank number in Character format for adding at the end of the input file Name.
Character (Kind = 1, Len = 20)  :: IndexSize ;   ! A variable for storing the Size number in Character format for adding at the end of the input file Name.

! - Logical Variables -------------------------------------------------------------------------------------------------------------------------------
! - Type DECLERATIONS -------------------------------------------------------------------------------------------------------------------------------

! - Unit NUMBERS OF EXTERNAL FILES ------------------------------------------------------------------------------------------------------------------
Integer (Kind=Smll), PARAMETER  :: Un_Mat_Vis     = 821  ;            ! the Unit number of Inversion Data Structure index sets

! =========================== Subroutine CODE =======================================================================================================

! Important notice: Due to the restriction of HDF5 format, we use real(4) and integer(4) in the output.

Write (*,*)"Output Material-Visualization HDF5"


! - Writing down the input data for each process ----------------------------------------------------------------------------------------------------
  Do IParts = 1, NParts ;

    Write (IndexRank, *) IParts - 1_Shrt ; ! Converts Rank number to Character foramt for the file Name
    Write (IndexSize, *) NParts ;          ! Converts Size number to Character foramt for the file Name

    Write (*,*)"Writing files for partition: ", IParts;

    ! - Computations --------------------------------------------------------------------------------------------------------------------------------
    ! Create files based on HDF5 format for large files
    Call h5open_f(ErrPTC)

    ! - Paraview geometry file  ---------------------------------------------------------------------------------------------------------------------
    ! determine the number of nodes in the model for Paraview model

    Allocate ( ParaNodeNum ( NJ ) ) ; ! 

    ParaNodeNum (:) = 0_Lng ;
    NJ_Para         = 0_Lng ;
    ConnSizePara    = 0_Lng ;

      ! This loop finds the nodes of Finite Elements, In other words, it eliminates spectral nodes.
      Do IEl = 1, NEl ;  !
        If ( EPart ( IEL ) == IParts ) Then ;

          EType = ElT ( IEl ) ; ! Element Type

          ! Determine the number of nodes corresponding to a FE.
            IF      ( EType == El2d4NSldPN  ) Then ;  NNode = 4_Tiny ; ! 4 node  -2D - Solid - PLANE STRESS
            Else If ( EType == El2d4NPMLPN  ) Then ;  NNode = 4_Tiny ; ! 4 node  -2D - PML - PLANE STRESS
            Else If ( EType == El2d8NSldPN  ) Then ;  NNode = 8_Tiny ; ! 8 node  -2D - Solid - PLANE STRAIN
            Else If ( EType == El2d8NPMLPN  ) Then ;  NNode = 8_Tiny ; ! 8 node  -2D - PML - PLANE STRAIN
            Else If ( EType == El2d6NSldPN  ) Then ;  NNode = 6_Tiny ; ! 6 node  -2D - Solid - PLANE STRAIN
            Else If ( EType == El2d6NPMLPN  ) Then ;  NNode = 6_Tiny ; ! 6 node  -2D - PML - PLANE STRAIN
            Else If ( EType == El2d3NSldPN  ) Then ;  NNode = 3_Tiny ; ! 3 node  -2D - Solid - PLANE STRAIN
            Else If ( EType == El2d3NPMLPN  ) Then ;  NNode = 3_Tiny ; ! 3 node  -2D - PML - PLANE STRAIN
            Else If ( EType == SEl2d9NSldPN ) Then ;  NNode = 8_Tiny ; ! 9 node  -2D - Solid - PLANE STRAIN
            Else If ( EType == SEl2d9NPMLPN ) Then ;  NNode = 8_Tiny ; ! 9 node  -2D - NPML - PLANE STRAIN    !?? using 8 noded element for spectral 9 noded elements
            Else If ( EType == SEl2d9NMPMLPN) Then ;  NNode = 8_Tiny ; ! 9 node  -2D - MPML - PLANE STRAIN
            Else If ( EType == SEl3d27NSld  ) Then ;  NNode = 20_Tiny; ! 27 node -3D - Solid
            Else If ( EType == SEl3d27NMPML ) Then ;  NNode = 20_Tiny; ! 27 node -3D - MPML
            Else If ( EType == SEl2d7NSldPN ) Then ;  NNode = 6_Tiny ; ! 7 node -2D - solid - plane strain
            Else If ( EType == SEl2d7NMPMLPN) Then ;  NNode = 6_Tiny ; ! 7 node -2D - MPML - plane strain
            Else If ( EType == SEl2d9NSldSH ) Then ;  NNode = 8_Tiny ; ! 9 node - 2D - solid ( quadrilateral ) - Spectral Element - SH waves
            Else If ( EType == SEl2d9NPMLSH ) Then ;  NNode = 8_Tiny ; ! 9 node - 2D - PML   ( quadrilateral ) - Spectral Element - SH waves
            Else If ( EType == SEl3d8NSld   ) Then ;  NNode = 8_Tiny ; ! 8 node -3D - Solid
            Else If ( EType == SEl3d8NMPML  ) Then ;  NNode = 8_Tiny ; ! 8 node -3D - MPML
            Else ; Write(*,*)"Element number is not in the list"; Stop ;
            End If ;

          ConnSizePara = ConnSizePara + Int4(NNode + 1_Tiny) ;

            Do INode = 1, NNode;

              Node = INod (INode, IEl) ;

                If ( ParaNodeNum ( Node ) == 0_Lng ) Then ;  ! This if condition is for determining the number of FEM nodes
                  NJ_Para = NJ_Para + 1_Lng ;
                  ParaNodeNum ( Node ) = NJ_Para ;           ! This array renumbers nodes for Paraview. THE ORDER OF NODE NUMBERING IN CONSISTENT EVERYWHERE.
                End If ;

            End Do ;

        End If ;

      End Do ;

    DeAllocate ( ParaNodeNum ) ;

    Print *, 'OutDir - Parview', OutDir ;
    Call h5fcreate_f( TRIM(OutDir)//'/'//TRIM(ModelName)//'_'//'Data_Mat_Vis_PV'//'_'//Trim(AdjustL(IndexSize))//'_'//Trim(AdjustL(IndexRank))//'.h5', H5F_ACC_TRUNC_F, id_Data_Mat_Vis_PV, ErrPTC) ;         ! Index number of DOFs in Paraview to retrive data from the solution vector

    Allocate ( dset_data_int( NJ_Para ) ) ;

    Counter = 0_Lng ;
      DO K = 1, NJ ;
        If ( Local_PETSc_Num ( K, IParts ) /= 0_Lng ) Then ;

          Counter = Counter + 1_Lng ;
          dset_data_int  ( Counter ) = Global_Petsc_Num ( K ) ; ! Equation number of nodes on Paraview to retrieve solution

          If ( Counter == NJ_Para ) Exit ;
        End If ;
      End Do ;

    ! Equation number of nodes on this rank used to extract data of each node for Paraview
    ! Remark: We already have equation numbers of this rank in Indices vector; However, we cannot use this vector to scatter solution for Paraview, because Paraview does not support Spectral Elements. Hence, we do not need to scatter all solution.
    NE_Mat_Vis_Para = Count ( dset_data_int ( : ) /= 0_Lng ) ;

     if ( NE_Mat_Vis_Para /= NJ_Para ) then
        write(*,*) 'NE_Mat_Vis_Para, NJ_Para', NE_Mat_Vis_Para, NJ_Para
        write(*,*) 'material visualization - error'
        stop
    end if

    dims(1) = NE_Mat_Vis_Para ;
    dims(2) = 1 ;
    Call h5screate_simple_f(2, dims, dspace_id_IS_Mat_Vis_PV, ErrPTC) ;
    Call h5dcreate_f(id_Data_Mat_Vis_PV, "IndexSet_Mat_Vis",     H5T_NATIVE_INTEGER, dspace_id_IS_Mat_Vis_PV, dset_id_IS_Mat_Vis_PV,   ErrPTC) ;

    dset_data_int = dset_data_int - 1 ; ! Equation numbers start from 0 in PETSc >>> see if this has already been done.

!    if (iparts == 2)  then 
!       write(*,*)  dset_data_int; stop
!    end if

    ! Write the dataset for index setting of the equation numbers
    data_dims(1) = NE_Mat_Vis_Para ;
    data_dims(2) = 1 ;
    Call h5dwrite_f(dset_id_IS_Mat_Vis_PV, H5T_NATIVE_INTEGER, dset_data_int, data_dims, ErrPTC)

    DeAllocate ( dset_data_int ) ;


! - Basic Data -------------------------------------------------------------------------------------------------------------------------------------
! - Geometry & Model -------------------------------------------------------------------------------------------------------------------------------
! <><><><><><><><><><><><><><><><><><><><><>
    ! Application Numbering


    ! - Closing the output file ---------------------------------------------------------------------------------------------------------------------


    ! Close the HDF5 file.
    CALL h5fclose_f( id_Data_Mat_Vis_PV, ErrPTC) ;

    ! Shut down HDF5.
    CALL h5close_f(ErrPTC) ;

    ! - ---------------------------------------------------------------------------------------------------------------------------------------------
    Write (*,*)"Opening files for Mat_Vis ...";

    ! 
    UnFile = Un_Mat_Vis ;
    Open ( Unit = UnFile, FILE = TRIM(ModelName)//'_'//Trim(AdjustL(IndexSize))//'_'//Trim(AdjustL(IndexRank))//'.Mat_Vis', ERR =  1001, IOSTAT = IO_File, ACCESS = 'SEQUENTIAL', ACTION ='Write', ASYNCHRONOUS = 'NO', BLANK = 'NULL', BLOCKSIZE = 0, DEFAULTFILE = TRIM(OutDir), DisPOSE = 'KEEP', FORM = 'FORMATTED', POSITION = 'ASIS', STATUS ='REPLACE') ;

! - Write index sets --------------------------------------------------------------------------------------------------------------------------------
    UnFile = UN_Mat_Vis ;
!    Write (*,*)"Writing Basic Data ..." ;
    Write (Unit = UnFile, FMT = "(1(I10,1X),'NE_Mat_Vis_Para')", ADVANCE = 'YES', ASYNCHRONOUS = 'NO', IOSTAT = IO_Write, ERR = 1006 ) NE_Mat_Vis_Para ;

    ! - Closing the output file ---------------------------------------------------------------------------------------------------------------------
    UnFile =  Un_Mat_Vis ;
    Close ( Unit = UnFile, Status = 'KEEP', ERR =  1002, IOSTAT = IO_File ) ;
    ! - ---------------------------------------------------------------------------------------------------------------------------------------------


  End Do ;


Write(*    ,*) 'End Subroutine < OUTPUT_Mat_Vis_HDF5 >' ;
Write(UnInf,*) 'End Subroutine < OUTPUT_Mat_Vis_HDF5 >' ;
Return ;

! =============================================== OPEN ERRORS =======================================================================================
1001  IF ( IO_File > 0 ) Then ;
        Write(*, Fmt_ERR1_OPEN) UnFile, IO_File  ;  Write(UnInf, Fmt_ERR1_OPEN) UnFile, IO_File  ; 
        Write(*, Fmt_FL) ; Write(UnInf, Fmt_FL) ;
        !#Call BEEP_FAIL ;
        Write(*, Fmt_End) ; Read(*,*) ;  Return ;
      Else If ( IO_File < 0 ) Then ;
        Write(*, Fmt_ERR1_OPEN) UnFile, IO_File  ; 
        Write(UnInf, Fmt_ERR1_OPEN) UnFile, IO_File  ;  Write(*, Fmt_FL) ; Write(UnInf, Fmt_FL) ;
        !#Call BEEP_FAIL ;
        Write(*, Fmt_End) ; Read(*,*) ;  Return ;
      End If ;


! =============================================== Close ERRORS ======================================================================================
1002  IF ( IO_File > 0 ) Then ;
        Write(*, Fmt_ERR1_Close) UnFile, IO_File  ;  Write(UnInf, Fmt_ERR1_Close) UnFile, IO_File  ; 
        Write(*, Fmt_FL) ; Write(UnInf, Fmt_FL) ;
        !Call BEEP_FAIL ;
        Write(*, Fmt_End) ; Read(*,*) ; ; Return ;
      End If ;

! =============================================== ERROR IN Write STATEMENT ==========================================================================
1006    Write(*       , Fmt_Write1 ) UnFile, IO_Write ; Write( UnFile, Fmt_Write1 ) UnFile, IO_Write ;
        !#Call BEEP_FAIL ;
        Write(*, Fmt_FL) ;  Write(UnInf, Fmt_FL) ; Write(*, Fmt_End) ; Read(*,*) ;  STOP ;


End Subroutine OUTPUT_Mat_Vis_HDF5 ;

End Module Mat_Vis_HDF5  ;
