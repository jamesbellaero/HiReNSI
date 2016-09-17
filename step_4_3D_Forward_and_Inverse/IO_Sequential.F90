
!****************************************************************************************************************************************************
! This component reads my old input files, sequentially
! June 24 2013
!****************************************************************************************************************************************************

CALL IO_FILES
CALL INPUT1()

ALLOCATE ( XYZ_serendipity(NJ,NDIM), ID_serendipity(NJ,NDOF), INOD(NNODE,NEL), NGP(NEL), MTEL(NEL), ID_BC(NEL,NDIM**2), PMAT(NMAT,NPM) )
ALLOCATE ( LTRANS(NTRANS) )
ALLOCATE ( BACL(NSTEP,NDIM) )
ALLOCATE ( XYZ_Lagrange(NEL, 6 * NDIM), ID_Lagrange(NEL, 6) )
ALLOCATE ( PML_PARAM(NDIM * 2 , 4) )

CALL INPUT_LTRANS ( LTRANS )
CALL INPUT_PML ( PML_PARAM )
!CALL INPUT2         ( XYZ_serendipity, ID_serendipity, INOD, NGP, MTEL, ID_BC, PMAT, BACL, XYZ_Lagrange, ID_Lagrange )
CALL ANSYS_SEZGIN    ( XYZ_serendipity, ID_serendipity, INOD, NGP, MTEL, ID_BC, PMAT, PML_PARAM )
!CALL ANSYS_3D_Import ( XYZ_serendipity, ID_serendipity, INOD, NGP, MTEL, ID_BC, PMAT, PML_PARAM, XYZ_Lagrange, ID_Lagrange )
!CALL ANSYS_3D_Import_2D ( XYZ_serendipity, ID_serendipity, INOD, NGP, MTEL, ID_BC, PMAT, PML_PARAM, XYZ_Lagrange, ID_Lagrange )

! this part deals with Lagrange element generation from serendipity data
!---------- ---------- ---------- ---------- ----------
ALLOCATE ( XYZ(NJ,NDIM), ID(NJ,NDOF) )
DO I = 1, NJ
   XYZ(I,:) = XYZ_serendipity(I,:)
    ID(I,:) =  ID_serendipity(I,:)
END DO
DEALLOCATE ( XYZ_serendipity, ID_serendipity )
DEALLOCATE ( XYZ_Lagrange   , ID_Lagrange )
!---------- ---------- ---------- ---------- ----------
!CALL Pine_Flat_3D_serendipity_to_Lagrange ( XYZ, ID )

!CALL ANSYS_3D_Process ( XYZ, ID, INOD, NGP, MTEL, ID_BC, PMAT, PML_PARAM )
!CALL ANSYS_3D_Process_2D ( XYZ, ID, INOD, NGP, MTEL, ID_BC, PMAT, PML_PARAM )
CALL Assign_Equation_Numbers ( ID )
CALL Print_Out_Input_Data ( XYZ, ID, INOD, NGP, MTEL, ID_BC, PMAT )


!---------- ---------- ---------- ---------- ----------
!CALL Solve_PML_3D_2nd_order_ODE ( XYZ, INOD, NGP, MTEL, PMAT, PML_PARAM, ID, ID_BC, LTRANS, SIZE, RANK )
!GOTO 1031
!---------- ---------- ---------- ---------- ----------
