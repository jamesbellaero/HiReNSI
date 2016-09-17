!************************************************
! Visualization with Paraview
!************************************************
! REVISION : M 22 Oct 2012

SUBROUTINE ParaView ( XYZ, INOD, ID, U_PETSC, ISTEP )
USE PARAMETERS
IMPLICIT DOUBLE PRECISION (A-H,O-Z)


!==================================================================================================
! - PETSC INCLUDE FILES
#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"


Vec                     :: U_PETSC
PetscScalar, pointer    :: U(:)
!==================================================================================================


! in
!---------- ---------- ---------- ---------- ---------- 
DIMENSION XYZ(NJ, NDIM), INOD(NNODE, NEL), ID(NJ,NDOF)
!---------- ---------- ---------- ---------- ----------


! local
!---------- ---------- ---------- ---------- ----------
!Character (LEN=60)      :: OutDir = '/home/arash/2D_Forward_Unix/vis'
!Character (LEN=60)      :: OutDir = '/terra1/arash/vis'
Character (LEN=60)      :: Index
Integer                 :: IO_File
Integer                 :: IO_Write
Integer                 :: UnFile
!---------- ---------- ---------- ---------- ----------


CALL VecGetArrayF90 ( U_PETSC, U, IERR )


! open the output file for each time step -------
write(index, *) IStep
UnFile = UN_PView

Open ( Unit = UnFile, FILE = TRIM(ModelName)//'_'//Trim(AdjustL(Index))//'.vtu', ACTION = 'Write', DEFAULTFILE = TRIM(OutDir), DISPOSE = 'KEEP', STATUS = 'REPLACE' )


! - header --------------------------------------
Write (Unit = UnFile, FMT = "(A73)",                    ADVANCE = 'YES', ASYNCHRONOUS = 'NO', IOSTAT = IO_Write) '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">'
Write (Unit = UnFile, FMT = "('  <UnstructuredGrid>')", ADVANCE = 'YES', ASYNCHRONOUS = 'NO', IOSTAT = IO_Write)
Write (Unit = UnFile, FMT = "(A27,I19,A17,I19,A2)",     ADVANCE = 'YES', ASYNCHRONOUS = 'NO', IOSTAT = IO_Write) '    <Piece NumberOfPoints="', NJ, '" NumberOfCells="', NEl, '">'
Write (Unit = UnFile, FMT = "('')",                     ADVANCE = 'YES', ASYNCHRONOUS = 'NO', IOSTAT = IO_Write)


! Points, i.e., nodes in our language -----------
Write (Unit = UnFile, FMT = "('      <Points>')",       ADVANCE = 'YES', ASYNCHRONOUS = 'NO', IOSTAT = IO_Write)
Write (Unit = UnFile, FMT = "(A64)",                    ADVANCE = 'YES', ASYNCHRONOUS = 'NO', IOSTAT = IO_Write) '<DataArray type="Float32" NumberOfComponents="3" format="ascii">'


! coordinates -----------------------------------
Do IJ = 1, NJ
   If (NDim == 2) Then
      Write (Unit = UnFile, FMT = "(3(E15.09,1x))",     ADVANCE = 'YES', ASYNCHRONOUS = 'NO', IOSTAT = IO_Write) (XYZ (IJ,IDim), IDim = 1,NDim), 0.0d0
   Else If (NDim == 3) Then
      Write (Unit = UnFile, FMT = "(3(E15.09,1x))",     ADVANCE = 'YES', ASYNCHRONOUS = 'NO', IOSTAT = IO_Write) (XYZ (IJ,IDim), IDim = 1,NDim)
   End If
End Do


Write (Unit = UnFile, FMT = "('        </DataArray>')", ADVANCE = 'YES', ASYNCHRONOUS = 'NO', IOSTAT = IO_Write)
Write (Unit = UnFile, FMT = "('      </Points>')",      ADVANCE = 'YES', ASYNCHRONOUS = 'NO', IOSTAT = IO_Write)


! - Cells, i.e., Elements in our language -------
Write (Unit = UnFile, FMT = "('      <Cells>')",        ADVANCE = 'YES', ASYNCHRONOUS = 'NO', IOSTAT = IO_Write)
Write (Unit = UnFile, FMT = "(A67)",                    ADVANCE = 'YES', ASYNCHRONOUS = 'NO', IOSTAT = IO_Write) '        <DataArray type="Int32" Name="connectivity" format="ascii">'


SELECT CASE (NNODE)

   CASE (8)
      ElType_PView = 23

   CASE (9) 
      ElType_PView = 23

   CASE DEFAULT
      WRITE(*,'("ERROR IN Visualization")')
      STOP

END SELECT


Do IEL = 1, NEL
! subtract 1 because nodes in ParaView start form 0.
   Write (Unit = UnFile, FMT = "(<NNode>(I19,1x))",     ADVANCE = 'YES', ASYNCHRONOUS = 'NO', IOSTAT = IO_Write) ( INod (INode, IEl)-1, INode = 1, NNode)
END DO


Write (Unit = UnFile, FMT = "('        </DataArray>')", ADVANCE = 'YES', ASYNCHRONOUS = 'NO', IOSTAT = IO_Write)
Write (Unit = UnFile, FMT = "(A62)",                    ADVANCE = 'YES', ASYNCHRONOUS = 'NO', IOSTAT = IO_Write) '        <DataArray type="Int32" Name="offsets" format="ascii">' ;


! - OffSet --------------------------------------
OffSet = 0
Do IEl = 1, NEl
   OffSet = OffSet + NNode
   Write (Unit = UnFile, FMT = "(I19,1X)",              ADVANCE = 'Yes', ASYNCHRONOUS = 'NO', IOSTAT = IO_Write) INT(OffSet)
End Do


Write (Unit = UnFile, FMT = "('        </DataArray>')", ADVANCE = 'YES', ASYNCHRONOUS = 'NO', IOSTAT = IO_Write)
Write (Unit = UnFile, FMT = "(A60)",                    ADVANCE = 'YES', ASYNCHRONOUS = 'NO', IOSTAT = IO_Write) '        <DataArray type="UInt8" Name="types" format="ascii">'


! - Element Type --------------------------------
Do IEl = 1, NEl
   Write (Unit = UnFile, FMT = "(I3,1X)",               ADVANCE = 'No',  ASYNCHRONOUS = 'NO', IOSTAT = IO_Write) INT(ElType_PView)
End Do


Write (Unit = UnFile, FMT = "()",                       ADVANCE = 'Yes', ASYNCHRONOUS = 'NO', IOSTAT = IO_Write)
Write (Unit = UnFile, FMT = "('        </DataArray>')", ADVANCE = 'YES', ASYNCHRONOUS = 'NO', IOSTAT = IO_Write)
Write (Unit = UnFile, FMT = "('      </Cells>')",       ADVANCE = 'YES', ASYNCHRONOUS = 'NO', IOSTAT = IO_Write)


! - Points data ---------------------------------
Write (Unit = UnFile, FMT = "('      <PointData>')",    ADVANCE = 'YES', ASYNCHRONOUS = 'NO', IOSTAT = IO_Write)


! Total Displacement ----------------------------
Write (Unit = UnFile, FMT = "(A66)",                    ADVANCE = 'YES', ASYNCHRONOUS = 'NO', IOSTAT = IO_Write) '        <DataArray type="Float32" Name="total dis" format="ascii">'


Do IJ = 1, NJ
   PointData = 0.0d0
      Do IDim = 1, NDim
         If ( ID ( IJ , IDim ) /= 0 ) PointData = PointData + U ( ID ( IJ , IDim ) ) * U ( ID ( IJ , IDim ) )
      End Do
   PointData = Sqrt ( PointData )
   Write (Unit = UnFile, FMT = "(E14.6E3,1X)",          ADVANCE = 'Yes', ASYNCHRONOUS = 'NO', IOSTAT = IO_Write) PointData
End Do
Write (Unit = UnFile, FMT = "('        </DataArray>')", ADVANCE = 'YES', ASYNCHRONOUS = 'NO', IOSTAT = IO_Write)


I_X_Y_Visualization = 0
!---------- ---------- ---------- ---------- ----------
IF ( I_X_Y_Visualization == 1 ) THEN

! Displacement in X direction -------------------
   Write (Unit = UnFile, FMT = "(A66)",                    ADVANCE = 'YES', ASYNCHRONOUS = 'NO', IOSTAT = IO_Write) '        <DataArray type="Float32" Name="X dis" format="ascii">'

   Do IJ = 1, NJ
      PointData = 0.0d0
      If ( ID ( IJ , 1 ) /= 0 ) PointData = U ( ID ( IJ , 1 ) )
      Write (Unit = UnFile, FMT = "(E21.13E3,1X)",         ADVANCE = 'Yes', ASYNCHRONOUS = 'NO', IOSTAT = IO_Write) PointData
   End Do
   Write (Unit = UnFile, FMT = "('        </DataArray>')", ADVANCE = 'YES', ASYNCHRONOUS = 'NO', IOSTAT = IO_Write)


! Displacement in Y direction -------------------
   Write (Unit = UnFile, FMT = "(A66)",                    ADVANCE = 'YES', ASYNCHRONOUS = 'NO', IOSTAT = IO_Write) '        <DataArray type="Float32" Name="Y dis" format="ascii">'

   Do IJ = 1, NJ 
      PointData = 0.0d0
      If ( ID ( IJ , 2 ) /= 0 ) PointData = U ( ID ( IJ , 2 ) )
      Write (Unit = UnFile, FMT = "(E21.13E3,1X)",         ADVANCE = 'Yes', ASYNCHRONOUS = 'NO', IOSTAT = IO_Write) PointData
   End Do
   Write (Unit = UnFile, FMT = "('        </DataArray>')", ADVANCE = 'YES', ASYNCHRONOUS = 'NO', IOSTAT = IO_Write)


END IF
!---------- ---------- ---------- ---------- ----------


Write (Unit = UnFile, FMT = "('      </PointData>')",   ADVANCE = 'YES', ASYNCHRONOUS = 'NO', IOSTAT = IO_Write)


! - Cell data -----------------------------------
! - Footer --------------------------------------
Write (Unit = UnFile, FMT = "('    </Piece>')",         ADVANCE = 'YES', ASYNCHRONOUS = 'NO', IOSTAT = IO_Write)
Write (Unit = UnFile, FMT = "('  </UnstructuredGrid>')",ADVANCE = 'YES', ASYNCHRONOUS = 'NO', IOSTAT = IO_Write)
Write (Unit = UnFile, FMT = "('</VTKFile>')",           ADVANCE = 'YES', ASYNCHRONOUS = 'NO', IOSTAT = IO_Write)


! - Closing the output file ---------------------
UnFile = UN_PView
Close ( Unit = UnFile, Status = 'KEEP', IOSTAT = IO_File )


CALL VecRestoreArrayF90 ( U_PETSC, U, IERR )


RETURN
END
