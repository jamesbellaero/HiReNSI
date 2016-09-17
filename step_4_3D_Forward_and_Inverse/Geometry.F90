!************************************************
! Receives geometry based on serendipity-noded elements and constructs data structure for 9-noded spectral elements
!************************************************
! REVISION : F, 18 May 2012

SUBROUTINE GEOMETRY_2D_8_to_9 ( XYZ, ID, INOD, MTEL )
USE PARAMETERS
IMPLICIT DOUBLE PRECISION (A-H,O-Z)


! in 
!---------- ---------- ---------- ---------- ----------
DIMENSION XYZ(NJ,NDIM), ID(NJ,NDOF), INOD(NNODE,NEL), MTEL(NEL)
!---------- ---------- ---------- ---------- ----------


! local
!---------- ---------- ---------- ---------- ----------
DIMENSION XT(NDIM, NNODE_serendipity), FN(NNODE_serendipity), DFXI(NNODE_serendipity, NDIM), DJ(NDIM, NDIM)
DIMENSION X9(NDIM)
!---------- ---------- ---------- ---------- ----------


! 1. compute the coordinates of the 9th node
!---------- ---------- ---------- ---------- ----------
DO IEL = 1, NEL


! read existing geometry
   DO I = 1, NNODE_serendipity
      K = INOD(I,IEL)
      DO J = 1, NDIM
          XT(J,I) = XYZ(K,J)
      END DO
   END DO


! compute the coordinates of the 9th node
   X2 = 0.0D0
   X1 = 0.0D0
   
   CALL SHAPE8 ( X1, X2, FN, DFXI )

   DJ = MATMUL(XT,DFXI)
   DETJ = DJ(1,1) * DJ(2,2) - DJ(1,2) * DJ(2,1)
   IF (DETJ.LE.0) THEN
      WRITE(*,*)'|J|<0 NEGATIVE DET. JACOBIAN'
      READ(*,*)
      STOP
   END IF

   X9 = MATMUL ( XT, FN )


! store the new node in coordinate list XYZ
   I = NJ_serendipity + IEL
   DO J = 1, NDIM
      XYZ (I, J) = X9(J)
   END DO


! add the new node to element connectivity matrix INOD
   INOD(9, IEL) = I


! update (original, untouched) restraints ID (X-Y-Sxx-Syy-Sxy)
   MTYPE   = MTEL(IEL)
   
   IF ( MTYPE == 1 ) THEN

      SELECT CASE (NDOF)

         CASE (2)                                ! regular domain
            ID (I, 1) = 0                        ! X
            ID (I, 2) = 0                        ! Y
                     
         CASE (5)                                ! PML
            ID (I, 1) = 0                        ! X
            ID (I, 2) = 0                        ! Y
            ID (I, 3) = 1                        ! Sxx
            ID (I, 4) = 1                        ! Syy
            ID (I, 5) = 1                        ! Sxy
       
         CASE DEFAULT
            WRITE(*,'("FLAG ERROR IN Geometry")')
            STOP

      END SELECT

   ELSE

      SELECT CASE (NDOF)

         CASE (2)                                ! regular domain
            ID (I, 1) = 0                        ! X
            ID (I, 2) = 0                        ! Y
            WRITE(*,*) 'Geometry: this should be wrong! if there are PML elements, NDOF should be 5, not 2'
            STOP

         CASE (5)                                ! PML
            ID (I, 1) = 0                        ! X
            ID (I, 2) = 0                        ! Y
            ID (I, 3) = 0                        ! Sxx
            ID (I, 4) = 0                        ! Syy
            ID (I, 5) = 0                        ! Sxy

         CASE DEFAULT
            WRITE(*,'("FLAG ERROR IN Geometry")')
            STOP

      END SELECT
      
   END IF


END DO
!---------- ---------- ---------- ---------- ---------- 1


RETURN
END


!************************************************
! Receives geometry based on 20-noded serendipity elements and constructs data structure for 27-noded Lagrange elements
!************************************************
! REVISION : Th, 24 May 2012 (Funeral of Dr. A. Noorzad)
! REVISION : M 22 April 2012

SUBROUTINE GEOMETRY_3D_20_to_27 ( XYZ, ID, INOD, XYZ_Lagrange, ID_Lagrange )
USE PARAMETERS
IMPLICIT DOUBLE PRECISION (A-H,O-Z)


! in 
!---------- ---------- ---------- ---------- ----------
DIMENSION XYZ(NJ,NDIM), ID(NJ,NDOF), INOD(NNODE,NEL)
!---------- ---------- ---------- ---------- ----------


! local
!---------- ---------- ---------- ---------- ----------
DIMENSION XT(NDIM, NNODE_serendipity), FN(NNODE_serendipity), DFXI(NNODE_serendipity, NDIM), DJ(NDIM, NDIM)
DIMENSION X_int(NDIM), X_1(NDIM), X_2(NDIM)
!---------- ---------- ---------- ---------- ----------


! out
!---------- ---------- ---------- ---------- ----------
DIMENSION XYZ_Lagrange(NEL, 6 * NDIM), ID_Lagrange(NEL, 6)
!---------- ---------- ---------- ---------- ----------


! initialize
!---------- ---------- ---------- ---------- ----------
ID_Lagrange = 0
NJ_Lagrange = 0


! part i: interpolate unknown points
!---------- ---------- ---------- ---------- ----------
DO IEL = 1, NEL


! read existing geometry
   DO I = 1, NNODE_serendipity
      K = INOD(I,IEL)
      DO J = 1, NDIM
         XT(J,I) = XYZ(K,J)
      END DO
   END DO


! compute the coordinates of nodes located in the middle of element faces (1>21  2>22  3>23  4>24  5>25  6>26)
   DO I = 1, 2 * NDIM


      SELECT CASE (I)

         CASE(1)
            X1 =  0.0D0 ; X2 =  0.0D0 ; X3 = -1.0D0

         CASE(2)
            X1 =  1.0D0 ; X2 =  0.0D0 ; X3 =  0.0D0

         CASE(3)
            X1 =  0.0D0 ; X2 =  1.0D0 ; X3 =  0.0D0

         CASE(4)
            X1 = -1.0D0 ; X2 =  0.0D0 ; X3 =  0.0D0

         CASE(5)
            X1 =  0.0D0 ; X2 = -1.0D0 ; X3 =  0.0D0

         CASE(6)
            X1 =  0.0D0 ; X2 =  0.0D0 ; X3 =  1.0D0

      END SELECT


      CALL SHAPE_3_20 ( X1, X2, X3, FN, DFXI )

      DJ = MATMUL(XT,DFXI)

      DETJ = DJ(1,1) * DJ(2,2) * DJ(3,3) + DJ(2,1) * DJ(3,2) * DJ(1,3) + DJ(3,1) * DJ(1,2) * DJ(2,3) - &
             DJ(1,3) * DJ(2,2) * DJ(3,1) - DJ(2,3) * DJ(3,2) * DJ(1,1) - DJ(3,3) * DJ(1,2) * DJ(2,1)
 
      IF (DETJ.LE.0) THEN
         WRITE(*,*)'|J|<0 NEGATIVE DET. JACOBIAN'
         READ(*,*)
         STOP
      END IF

      X_int = MATMUL ( XT, FN )

      I3 = 3 * I
      I2 = I3 - 1
      I1 = I2 - 1

      XYZ_Lagrange (IEL, I1) = X_int(1)
      XYZ_Lagrange (IEL, I2) = X_int(2)
      XYZ_Lagrange (IEL, I3) = X_int(3)

   END DO


END DO
!---------- ---------- ---------- ---------- ----------


! part ii: search for similar nodes and identify them
!---------- ---------- ---------- ---------- ----------
TOL = 1.0d-6

DO IEL = 1, NEL-1

   DO L_node =  1, 2 * NDIM

      IF ( ID_Lagrange(IEL, L_node) /= 0 ) THEN 
         CYCLE
      ELSE
         NJ_Lagrange = NJ_Lagrange + 1
         ID_Lagrange(IEL, L_node) = NJ_Lagrange

         I3 = 3 * L_node
         I2 = I3 - 1
         I1 = I2 - 1
         X_1(1) = XYZ_Lagrange (IEL, I1)
         X_1(2) = XYZ_Lagrange (IEL, I2)
         X_1(3) = XYZ_Lagrange (IEL, I3)

! inner loop
!---------- ---------- ---------- ---------- ----------
         ID_Flag = 0

         DO JEL = IEL + 1, NEL

            IF ( ID_Flag == 1 ) EXIT

            DO M_node =  1, 2 * NDIM

               IF ( ID_Lagrange(JEL, M_node) /= 0 ) THEN
                  CYCLE
               ELSE
                  I3 = 3 * M_node
                  I2 = I3 - 1
                  I1 = I2 - 1
                  X_2(1) = XYZ_Lagrange (JEL, I1)
                  X_2(2) = XYZ_Lagrange (JEL, I2)
                  X_2(3) = XYZ_Lagrange (JEL, I3)
                  dist = ( X_2(1) - X_1(1) )**2 + ( X_2(2) - X_1(2) )**2 + ( X_2(3) - X_1(3) )**2 
                  IF ( dist < TOL ) THEN 
                     ID_Lagrange(JEL, M_node) = NJ_Lagrange
                     ID_Flag = 1
                  END IF
              
                END IF

            END DO

         END DO
!---------- ---------- ---------- ---------- ----------

      END IF

   END DO

END DO

! last row needs special care
IEL = NEL
DO L_node =  1, 2 * NDIM

   IF ( ID_Lagrange(IEL, L_node) /= 0 ) THEN
      CYCLE
   ELSE
      NJ_Lagrange = NJ_Lagrange + 1
      ID_Lagrange(IEL, L_node) = NJ_Lagrange
   END IF

END DO


! part iii: store Lagrange data in appropriate data structure
!---------- ---------- ---------- ---------- ----------
! 0. total number of joints
NJ = NJ_serendipity + NJ_Lagrange + NEL

! 1. update ID_Lagrange according to the total number of joints
ID_Lagrange = ID_Lagrange + NJ_serendipity

! 2. add Lagrange face nodes (21-26) and Lagrange center node (27) to the element connectivity array
DO IEL = 1, NEL
   INOD(21, IEL) = ID_Lagrange(IEL, 1)
   INOD(22, IEL) = ID_Lagrange(IEL, 2)
   INOD(23, IEL) = ID_Lagrange(IEL, 3)
   INOD(24, IEL) = ID_Lagrange(IEL, 4)
   INOD(25, IEL) = ID_Lagrange(IEL, 5)
   INOD(26, IEL) = ID_Lagrange(IEL, 6)
   INOD(27, IEL) = NJ_serendipity + NJ_Lagrange + IEL
END DO

! 3. store Lagrange nodes in coordinate list XYZ
DO IEL = 1, NEL

! Lagrange face nodes (21-26)
   DO L_node =  1, 2 * NDIM
      I = ID_Lagrange(IEL, L_node)

      I3 = 3 * L_node
      I2 = I3 - 1
      I1 = I2 - 1
      XYZ(I,1) = XYZ_Lagrange (IEL, I1)
      XYZ(I,2) = XYZ_Lagrange (IEL, I2)
      XYZ(I,3) = XYZ_Lagrange (IEL, I3)
   END DO

! Lagrange center node (27)
! read serendipity geometry
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

   X_int = MATMUL ( XT, FN )

   I = NJ_serendipity + NJ_Lagrange + IEL

   DO J = 1, NDIM
      XYZ(I,J) = X_int(J)
   END DO

END DO


RETURN
END
