!************************************************
! READS ANSYS.OUT AND PREPARES INPUT-2D
! assumption: coordinate system is located on the surface, in the middle.
!************************************************
! REVISION : W 24 April 2013

SUBROUTINE ANSYS_SEZGIN ( XYZ, ID, INOD, NGP, MTEL, ID_BC, PMAT, PML_PARAM )
USE PARAMETERS
IMPLICIT DOUBLE PRECISION (A-H,O-Z)


! out	  
!---------- ---------- ---------- ---------- ----------	  
DIMENSION XYZ(NJ,NDIM), ID(NJ,NDOF), INOD(NNODE,NEL), NGP(NEL), MTEL(NEL), ID_BC(NEL, NDIM**2)
DIMENSION PMAT(NMAT,NPM)
!---------- ---------- ---------- ---------- ----------	 


! in-out (update)	  
!---------- ---------- ---------- ---------- ----------
DIMENSION PML_PARAM(NDIM * 2 , 4)
!---------- ---------- ---------- ---------- ----------


! local	  
!---------- ---------- ---------- ---------- ----------	  
CHARACTER*100 HELP
DIMENSION ID_HELP(NJ_serendipity, NDIM)
DIMENSION XT(NDIM, NNODE_serendipity), FN(NNODE_serendipity), DFXI(NNODE_serendipity, NDIM), DJ(NDIM, NDIM)
DIMENSION X9(NDIM)
!---------- ---------- ---------- ---------- ----------


IF ( NDIM /= 2 ) RETURN


ID_HELP    = 0
ID_BC      = 0


! number of integration points. for diagonal mass matrix, 9-noded elements should be used with 3 Lobatto points (defined in INPUT1).
!---------- ---------- ---------- ---------- ----------
NGP        = 4


! read FEM geometry
READ (21,'(1X,D20.8)') rdwd
READ (21,'(1X,D20.8)') rdht
READ (21,'(1X,D20.8)') PML_PARAM( 2 , 3 )        ! L left
READ (21,'(1X,D20.8)') PML_PARAM( 1 , 3 )        ! L right
READ (21,'(1X,D20.8)') PML_PARAM( 4 , 3 )        ! L bottom
     
WRITE(22,'(F10.2)') rdwd
WRITE(22,'(F10.2)') rdht
WRITE(22,'(F10.2)') PML_PARAM( 2 , 3 )     
WRITE(22,'(F10.2)') PML_PARAM( 1 , 3 )     
WRITE(22,'(F10.2)') PML_PARAM( 4 , 3 )         


! read number of nodes and elements
READ (21, '(2(1X,I12))') NJ_no_use, NEL_no_use
READ (21, '(3(1X,I12))') n_nodesRD, n_nodesINT, n_nodesPML
READ (21, '(2(1X,I12))') n_elemsRD, n_elemsPML


! read joints
READ (21, '(1X,A5)') HELP
DO I = 1, NJ_serendipity
   READ (21, '(1X,I12,1X,2D20.8)') II, (XYZ(II,J), J = 1, NDIM)         
END DO
!DO I = 1, NJ_serendipity
!   WRITE(22, '(I10,3F10.2)') I,  (XYZ(I ,J), J = 1, NDIM)         
!END DO   
   
   
! read elements
! ANSYS uses a different ordering which is modified to be consistent with Dr. Lotfi's method
! Now, local and global coordinate systems coinside w.r.t. their directions. this is helpful when appnying boundary conditions on element's faces
READ (21, '(1X,A5)') HELP
DO J = 1, NEL
   READ (21, '(10(1X,I12))') I, INOD(4,I), INOD(1,I), INOD(2,I), INOD(3,I), &
                                INOD(8,I), INOD(5,I), INOD(6,I), INOD(7,I), MTEL(I)         
END DO


! read boundary data
READ (21, '(1X,A5)') HELP
DO
   READ (21, '(1X,A12)') HELP
   IF (TRIM(HELP) /= 'ENDS') THEN
      READ (HELP, '(I12)') II
      ID_HELP (II,1) = 1
   ELSE
      EXIT
   END IF
END DO
READ (21, '(1X,A5)') HELP
DO
   READ (21, '(1X,A12)') HELP
   IF (TRIM(HELP) /= 'ENDS') THEN
      READ (HELP, '(I12)') II
      ID_HELP (II,2) = 1
   ELSE
      EXIT
   END IF
END DO


! read surface load data (currently supports surface traction (PLOAD) on element edges and disc load (DLOAD))
READ (21, '(1X,A5)') HELP
READ (21, '(1X,I12)') NLEL
   
SELECT CASE ( TRIM(HELP) )

   CASE ('PLOAD')
      DO I = 1, NLEL
         READ (21, '(3(1X,I12))') IEL , ID_BC(IEL,1) , II
      END DO

   CASE ('DLOAD')
      DO I = 1, NLEL
         READ (21, '(1X,I12)') IEL
         ID_BC(IEL,1) = 11
      END DO
      READ(21,'(1X,D20.8)') x_c                  ! x-coordinate of disk center (I don't store disc geometry automatically)
      READ(21,'(1X,D20.8)') y_c                  ! y-coordinate of disk center
      READ(21,'(1X,D20.8)') r_d                  ! disk radius

   CASE DEFAULT
      WRITE(*,'("FLAG ERROR IN EM")')
      STOP

END SELECT


! complete the file  
! PML (update)
!---------- ---------- ---------- ---------- ----------
PML_PARAM( 1 , 4 ) =  rdwd
PML_PARAM( 2 , 4 ) = -rdwd
PML_PARAM( 3 , 4 ) =  0.0d0
PML_PARAM( 4 , 4 ) = -rdht


! if the element is not solid, it should be PML     
!---------- ---------- ---------- ---------- ----------
! PML starting point locations
PML_X01 = PML_PARAM(1,4)
PML_X02 = PML_PARAM(2,4)
PML_X03 = PML_PARAM(3,4)
PML_X04 = PML_PARAM(4,4)

DO IEL = 1, NEL

   DO I = 1, NNODE_serendipity
      K = INOD(I,IEL)
      DO J = 1, NDIM
         XT(J,I) = XYZ(K,J)
      END DO
   END DO

   X1 =  0.0D0 ; X2 =  0.0D0 

   CALL SHAPE8 ( X1, X2, FN, DFXI )

   DJ = MATMUL(XT,DFXI)
   DETJ = DJ(1,1) * DJ(2,2) - DJ(1,2) * DJ(2,1)
   IF (DETJ.LE.0) THEN
      WRITE(*,*)'|J|<0 NEGATIVE DET. JACOBIAN'
      READ(*,*)
      STOP
   END IF

   X9 = MATMUL ( XT, FN )

! check if the node resides in Omega_RD or Omega_PML
         IF ( X9(1) >= PML_X01 .OR. X9(1) <= PML_X02 .OR. X9(2) >= PML_X03 .OR. X9(2) <= PML_X04 )  THEN
            MTEL(IEL) = 2
         ELSE
            MTEL(IEL) = 1
         END IF

END DO


! restraints
!---------- ---------- ---------- ---------- ----------
ID = 1

! first, identify the degrees of freedom for each node and set them free (basic ID that contains all potentially active dofs)
DO IEL = 1 , NEL
   MTYPE = MTEL(IEL)
         
   SELECT CASE (MTYPE)
             
      CASE (1)                                   ! regular domain
         DO I = 1 , NNODE_serendipity
            DO J = 1 , NDIM
               ID(INOD(I, IEL), J) = 0
            END DO
         END DO    
             
      CASE (2, 3, 4, 5, 6)                       ! PML
         DO I = 1 , NNODE_serendipity
            DO J = 1 , NDOF
               ID(INOD(I, IEL), J) = 0
            END DO
         END DO            
             
      CASE DEFAULT
         WRITE(*,'("FLAG ERROR IN ANSYS_SEZGIN")')
         STOP

   END SELECT
         
END DO
      
! second, apply displacement boundary conditions identified by ANSYS
DO I = 1, NJ_serendipity
   DO J = 1, NDIM
      IF (ID_HELP(I,J) == 1) ID(I,J) = 1 
   END DO
END DO


! third:  surface boundary condition
!---------- ---------- ---------- ---------- ----------
I_surface_bc = 0

SELECT CASE (I_surface_bc)


   CASE(0)
!     do nothing (Sezgin)


   CASE(1)
! apply traction free bc on the surface for PML elements i.e. SYY = 0, SXY = 0
! note: coordinate system is on the surface
! note: Sezgin is not doing this
      TOL = 1.E-6
      DO IEL = 1, NEL
         MTYPE = MTEL(IEL)
         
          SELECT CASE (MTYPE)   
             
             CASE (2, 3, 4, 5, 6)                ! PML
                 DO I = 1 , NNODE
                     Y = XYZ( INOD(I, IEL) , 2 )
                     IF ( ABS(Y) < TOL ) THEN 
                         ID( INOD(I, IEL) , 4 ) = 1 
                         ID( INOD(I, IEL) , 5 ) = 1                                              
                     END IF 
                 END DO            
     
          END SELECT    
        
      END DO


   CASE(2)
! zero displacement bc on the surface for PML elements i.e. ux = 0, uy = 0
! note: coordinate system is on the surface
! note: I added this part when I saw instability initiates from this region in mixed FE elastodynamics. 
! This removes instability in mixed FE formulation; same as case(0); but instability is still present in the PML formulation.
      TOL = 1.E-6
      DO IEL = 1, NEL
         MTYPE = MTEL(IEL)

          SELECT CASE (MTYPE)

             CASE (2, 3, 4, 5, 6)                ! PML
                 DO I = 1 , NNODE
                     Y = XYZ( INOD(I, IEL) , 2 )
                     IF ( ABS(Y) < TOL ) THEN
                         ID( INOD(I, IEL) , 1 ) = 1
                         ID( INOD(I, IEL) , 2 ) = 1
                     END IF
                 END DO

          END SELECT

      END DO


   END SELECT

! print out restraints
!          DO I = 1, NJ
!              WRITE(22,'(7I10)') I, (ID(I,J), J = 1, NDOF)
!          END DO


! update geometry for 9-noded elements
!---------- ---------- ---------- ---------- ----------
IF ( NDIM == 2 .AND. NNODE == 9 ) THEN
   CALL GEOMETRY_2D_8_to_9 ( XYZ, ID, INOD, MTEL )
END IF


! assign equation numbers
!---------- ---------- ---------- ---------- ----------
!CALL Assign_Equation_Numbers ( ID )  


! print out stuff
!---------- ---------- ---------- ---------- ----------
!CALL Print_Out_Input_Data ( XYZ, ID, INOD, NGP, MTEL, ID_BC, PMAT )


! print out ONLY THE LAST LINEs of the final dof numbering
!          DO I = NJ - 5 , NJ
!              WRITE(95,'(7I10)') I, (ID(I,J), J = 1, NDOF)
!          END DO


! material properties
!---------- ---------- ---------- ---------- ----------	 
I_material_property = 0

SELECT CASE (I_material_property)

   CASE(0)                                       ! disc load - wave motion paper
      PMAT(1,1) = 1.250d3   !* 1.0d6
      PMAT(1,2) = 0.250d0    
      PMAT(1,3) = 2000.0d-6 !* 1.0d6
      PMAT(1,4) = 0.0d0    
      PMAT(1,5) = 1.0d0

      PMAT(2,1) = 1.250d3   !* 1.0d6
      PMAT(2,2) = 0.250d0    
      PMAT(2,3) = 2000.0d-6 !* 1.0d6 
      PMAT(2,4) = 0.0d0    
      PMAT(2,5) = 1.0d0

   CASE(1)                                       ! disc load - grazing waves - CMA paper example 2
      PMAT(1,1) = 178.80d3 
      PMAT(1,2) = 0.20d0
      PMAT(1,3) = 2200.0d0
      PMAT(1,4) = 0.0d0
      PMAT(1,5) = 1.0d0

      PMAT(2,1) = 178.80d3   
      PMAT(2,2) = 0.20d0
      PMAT(2,3) = 2200.0d0
      PMAT(2,4) = 0.0d0
      PMAT(2,5) = 1.0d0

   CASE(2)                                       ! Chopra's explicit example
      PMAT(1,1) = 2.50d0 * 100.0
      PMAT(1,2) = 0.250d0
      PMAT(1,3) = 1.0d0
      PMAT(1,4) = 0.0d0
      PMAT(1,5) = 1.0d0

      PMAT(2,1) = 2.50d0 * 100.0
      PMAT(2,2) = 0.250d0
      PMAT(2,3) = 1.0d0
      PMAT(2,4) = 0.0d0
      PMAT(2,5) = 1.0d0

END SELECT


! print out
!---------- ---------- ---------- ---------- ----------	 
DO I = 1, NMAT
   WRITE(22,'(I10, E10.5, 4F10.2)') I, (PMAT(I,J), J = 1, NPM)
END DO
  
DO I = 1 , 2 * NDIM
   WRITE(22,'(I10,4F10.2)') I , (PML_PARAM( I , J ) , J = 1 , 4)
END DO
     
     
RETURN
END  


!************************************************
! READS ANSYS.OUT AND PREPARES INPUT-3D
! assumption: coordinate system is located on the surface, in the middle; +z is upward.
!************************************************
! REVISION : T 7 August 2012

SUBROUTINE ANSYS_3D_Import ( XYZ, ID, INOD, NGP, MTEL, ID_BC, PMAT, PML_PARAM, XYZ_Lagrange, ID_Lagrange )
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
PML_PARAM( 3 , 3 ) = pml_thickness
PML_PARAM( 4 , 3 ) = pml_thickness
PML_PARAM( 5 , 3 ) = 0.0d0
PML_PARAM( 6 , 3 ) = pml_thickness

! PML starting point location
PML_PARAM( 1 , 4 ) = d_p_Lx - pml_thickness
PML_PARAM( 2 , 4 ) = d_m_Lx + pml_thickness
PML_PARAM( 3 , 4 ) = d_p_Ly - pml_thickness
PML_PARAM( 4 , 4 ) = d_m_Ly + pml_thickness
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
!************************************************
! REVISION : M 30 July 2012

SUBROUTINE ANSYS_3D_Process ( XYZ, ID, INOD, NGP, MTEL, ID_BC, PMAT, PML_PARAM )
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
! edges of the computational domain
      d_p_Lx = PML_PARAM(1,4)
      d_m_Lx = PML_PARAM(2,4)
      d_p_Ly = PML_PARAM(3,4)
      d_m_Ly = PML_PARAM(4,4)
      d_p_Lz = PML_PARAM(5,4)
      d_m_Lz = PML_PARAM(6,4)

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
            ID(I,1) = 1                          ! translation in x-dir is restricted
            ID(I,2) = 1                          ! translation in y-dir is restricted
            ID(I,3) = 1                          ! translation in z-dir is restricted   
         END IF

! y = d_m_Ly is fixed
         IF ( DABS(X_2 - d_m_Ly) < TOL ) THEN
            ID(I,1) = 1                          ! translation in x-dir is restricted
            ID(I,2) = 1                          ! translation in y-dir is restricted
            ID(I,3) = 1                          ! translation in z-dir is restricted   
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
      d_p_Ly = PML_PARAM(3,4) + PML_PARAM(3,3) 
      d_m_Ly = PML_PARAM(4,4) - PML_PARAM(4,3)
      d_p_Lz = PML_PARAM(5,4) + PML_PARAM(5,3)
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
         IF ( X27(1) >= PML_X01 .OR. X27(1) <= PML_X02 .OR. X27(2) >= PML_X03 .OR. X27(2) <= PML_X04 .OR. X27(3) >= PML_X05 .OR. X27(3) <= PML_X06 )  THEN
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
            ID(I,1) = 1                          ! translation in x-dir is restricted
            ID(I,2) = 1                          ! translation in y-dir is restricted
            ID(I,3) = 1                          ! translation in z-dir is restricted   
         END IF

! y = d_m_Ly is fixed
         IF ( DABS(X_2 - d_m_Ly) < TOL ) THEN
            ID(I,1) = 1                          ! translation in x-dir is restricted
            ID(I,2) = 1                          ! translation in y-dir is restricted
            ID(I,3) = 1                          ! translation in z-dir is restricted   
         END IF

      END DO

! 2. c. apply traction free bc on the surface of PML elements i.e. Szx = 0, Szy = 0, Szz = 0
! dofs in ID are sorted as follows: Ux, Uy, Uz, Sxx, Syy, Szz, Sxy, Sxz, Syz; therefore we are interested in column 6, 8, 9.
      DO IEL = 1, NEL
         MTYPE = MTEL(IEL)

          SELECT CASE (MTYPE)   

             CASE (2)                            ! PML
                DO I = 1 ,NNODE
                   Z = XYZ( INOD(I, IEL) , 3 )
                   IF ( ABS(Z) < TOL ) THEN 
                      ID( INOD(I, IEL) , 6 ) = 1 
                      ID( INOD(I, IEL) , 8 ) = 1                                              
                      ID( INOD(I, IEL) , 9 ) = 1
                   END IF 
                END DO            

          END SELECT    

      END DO


   END SELECT


RETURN
END
