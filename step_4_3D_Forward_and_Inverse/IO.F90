!************************************************
! IO
!************************************************
! REVISION : Th, 5 July 2012

SUBROUTINE IO_FILES
USE PARAMETERS
IMPLICIT DOUBLE PRECISION (A-H,O-Z)

OPEN(01, FILE = 'test03.txt',     STATUS='OLD', ACTION='READ')
!OPEN(01, FILE = 'test61.txt',     STATUS='OLD', ACTION='READ')


! 2D
!---------- ---------- ---------- ---------- ----------
!OPEN(02, FILE = 'TEST03inf.txt',       STATUS='UNKNOWN')
!OPEN(03, FILE = '03_Spectral_PML.txt', STATUS='UNKNOWN') 
!OPEN(21, FILE = 'IN_Ansys03.txt', STATUS='OLD', ACTION='READ')
!OPEN(22, FILE = 'OT_Ansys03.txt',      STATUS='UNKNOWN')


! Visualization
!---------- ---------- ---------- ---------- ----------
Vis_file_name = 'PML03_16'


! 3D
!---------- ---------- ---------- ---------- ----------
OPEN(02, FILE = 'TEST03inf.txt', STATUS='UNKNOWN')
OPEN(03, FILE = 'TEST03his_16.txt', STATUS='UNKNOWN')
!OPEN(03, FILE = 'TEST08his_node_20_order_2_Gauss_4_Im.txt', STATUS='UNKNOWN') 
!OPEN(03, FILE = 'TEST10his_delta_verification_node_27_Gauss_3_Ex_alpha_5_beta_5.txt', STATUS='UNKNOWN')
OPEN(21, FILE = 'IN_Ansys03.txt', STATUS='OLD', ACTION='READ')
OPEN(22, FILE = 'OT_Ansys.txt', STATUS='UNKNOWN')


! debugging files
!---------- ---------- ---------- ---------- ----------
 OPEN(95, FILE = '95_2dCheck.txt', STATUS='UNKNOWN')
 OPEN(96, FILE = '96_Check.txt', STATUS='UNKNOWN')
! OPEN(97, FILE = '97_Check.txt', STATUS='UNKNOWN')
! OPEN(98, FILE = '98_Check.txt', STATUS='UNKNOWN')
! OPEN(99, FILE = '99_Check.txt', STATUS='UNKNOWN')


RETURN 
END
     

!************************************************
! IMPORT INITIAL DATA FROM DATA FILE
!************************************************
! REVISION : F, 18 May 2012

SUBROUTINE INPUT1 ()
USE PARAMETERS
IMPLICIT DOUBLE PRECISION (A-H,O-Z)


! local 
!---------- ---------- ---------- ---------- ----------	 
CHARACTER*100 HELP
!---------- ---------- ---------- ---------- ----------


READ (1,*) HELP
WRITE(2,*) HELP
READ (1,*) HELP
WRITE(2,*) HELP
READ (1,*) HELP
WRITE(2,*) HELP
READ (1,*) HELP
WRITE(2,*) HELP
READ (1,*) HELP
WRITE(2,*) HELP
READ (1,*) HELP
WRITE(2,*) HELP
READ (1,*) HELP
WRITE(2,*) HELP
READ (1,*)        NJ, NDIM
WRITE(2,'(2I10)') NJ, NDIM

READ (1,*) HELP
WRITE(2,*) HELP
READ (1,*)        NEL, NNODE, NDOF
WRITE(2,'(3I10)') NEL, NNODE, NDOF

READ (1,*) HELP
WRITE(2,*) HELP
READ (1,*)        NMAT, NPM
WRITE(2,'(2I10)') NMAT, NPM

GRAVITY = 9.810d0

READ (1,*) HELP
WRITE(2,*) HELP
READ (1,*)              NSTEP, NASTEP, DT
WRITE(2,'(2I10,F10.5)') NSTEP, NASTEP, DT

READ (1,*) HELP
WRITE(2,*) HELP
READ (1,*)          DELTA, ALPHA, ALPHA_DAMP, BETA_DAMP
WRITE(2,'(4F10.5)') DELTA, ALPHA, ALPHA_DAMP, BETA_DAMP
  
READ (1,*) HELP
WRITE(2,*) HELP
READ (1,*)       NTRANS
WRITE(2,'(I10)') NTRANS


! accomodate data structure for spectral elements, and
! for diagonal mass matrix, quadratic spectral elements should be used with 3 Lobatto points.
!---------- ---------- ---------- ---------- ----------
NINT_Lobatto      = 4
NJ_serendipity    = NJ
NNODE_serendipity = NNODE

IF ( NDIM == 2 .AND. NNODE == 9 ) THEN
   NJ = NJ_serendipity + NEL
   NNODE_serendipity = NNODE - 1
END IF

IF ( NDIM == 3 .AND. NNODE == 27 ) THEN
   NJ = NJ_serendipity + NEL * 7                 ! this is a (brutal) rough estimation
   NNODE_serendipity = NNODE - 7
END IF
!---------- ---------- ---------- ---------- ----------


RETURN
END


!************************************************
! IMPORT SUPPLEMENTARY DATA FROM DATA FILE
!************************************************
! REVISION : M, 22 April 2013

SUBROUTINE INPUT2 ( XYZ_serendipity, ID_serendipity, INOD, NGP, MTEL, ID_BC, PMAT, BACL, XYZ_Lagrange, ID_Lagrange )
USE PARAMETERS
IMPLICIT DOUBLE PRECISION (A-H,O-Z)


! out	  
!---------- ---------- ---------- ---------- ----------	  
DIMENSION XYZ_serendipity(NJ, NDIM), ID_serendipity(NJ, NDOF), INOD(NNODE, NEL), NGP(NEL), MTEL(NEL), ID_BC(NEL, NDIM**2), PMAT(NMAT, NPM)
DIMENSION BACL(NSTEP,NDIM)
DIMENSION XYZ_Lagrange(NEL, 6 * NDIM), ID_Lagrange(NEL, 6)
!---------- ---------- ---------- ---------- ---------- 


! local	  
!---------- ---------- ---------- ---------- ----------	
CHARACTER*100 HELP
!---------- ---------- ---------- ---------- ----------

! geometry
!---------- ---------- ---------- ---------- ---------- 	  
READ (1,*) HELP
WRITE(2,*) HELP

! read geometry for serendipity-noded elements; construct geometry for spectral elements later if necessary
DO I = 1, NJ_serendipity
   READ (1,*  ) II, (XYZ_serendipity(I,J), J = 1, NDIM)
!  WRITE(2,'(I10, 3F10.2)') II, (XYZ(I,J), J = 1, NDIM)
END DO


! restraints (just read, do not assign equation numbers)
!---------- ---------- ---------- ---------- ---------- 
ID_serendipity = 0
READ (1,*) HELP
WRITE(2,*) HELP
DO I = 1, NJ_serendipity
   READ(1,*) II, (ID_serendipity(I,J), J = 1, NDOF)
END DO


! element connectivity
!---------- ---------- ---------- ---------- ----------	    
READ (1,*) HELP
WRITE(2,*) HELP

SELECT CASE (NDIM)
   
   CASE(2)
      DO I = 1, NEL
         READ (1,*) II, (INOD(J,I), J = 1, NNODE_serendipity), NGP(I), MTEL(I), (ID_BC(I,K), K = 1, NDIM**2)
!        WRITE(2,'(30I10)') II, (INOD(J,I), J = 1, NNODE_serendipity), NGP(I), MTEL(I), (ID_BC(I,K), K = 1, NDIM**2)
      END DO

   CASE(3)
      DO I = 1, NEL
         READ (1,*) II, (INOD(J,I), J = 1, 8), NGP(I), MTEL(I), (ID_BC(I,K), K = 1, 4)
         IF ( NNODE == 20 .OR. NNODE == 27 ) READ (1,*) (INOD(J,I), J = 9, 20)
!        WRITE(2,'(30I10)') II, (INOD(J,I), J = 1, 8), NGP(I), MTEL(I), (ID_BC(I,K), K = 1, 4)
!        IF ( NNODE == 20 ) WRITE (2,'(30I10)') (INOD(J,I), J = 9, 20)
      END DO

END SELECT


! material properties
!---------- ---------- ---------- ---------- ----------		  			  
READ (1,*) HELP
WRITE(2,*) HELP
DO I = 1, NMAT
   READ (1,*) II, (PMAT(I,J), J = 1, NPM)
!  WRITE(2,'(I10, E10.5, 4F10.2)') II, (PMAT(I,J), J = 1, NPM)
END DO


! make appropriate corrections for plane stress problems 
! (2D problem is formulated as a plane strain problem)
!---------- ---------- ---------- ---------- ----------
I_plane_stress = 0
IF ( NDIM == 2 .AND. I_plane_stress == 1 ) THEN

   DO I = 1, NMAT
      E   = PMAT(I,1)
      ANU = PMAT(I,2)

      PMAT(I,1) = E * ( 1.0d0 + 2.0d0 * ANU ) / ( 1.0d0 + ANU )**2 
      PMAT(I,2) = ANU / ( 1.0d0 + ANU )
   END DO

END IF 


! earthquake acceleration time-history
!---------- ---------- ---------- ---------- ----------	
CALL INPUT4(BACL)


! update geometry for spectral elements
!---------- ---------- ---------- ---------- ----------
IF ( NDIM == 2 .AND. NNODE == 9 ) THEN
   CALL GEOMETRY_2D_8_to_9   ( XYZ_serendipity, ID_serendipity, INOD, MTEL )
END IF

IF ( NDIM == 3 .AND. NNODE == 27 ) THEN
   CALL GEOMETRY_3D_20_to_27 ( XYZ_serendipity, ID_serendipity, INOD, XYZ_Lagrange, ID_Lagrange )
END IF


! assign equation numbers
!---------- ---------- ---------- ---------- ----------
!CALL Assign_Equation_Numbers ( ID )


! print out stuff
!---------- ---------- ---------- ---------- ----------
!CALL Print_Out_Input_Data ( XYZ_serendipity, ID, INOD, NGP, MTEL, ID_BC, PMAT )


RETURN
END


!************************************************
! Assign equation numbers according to the restraints ID
!************************************************
! REVISION : T 7 Aug 2012

SUBROUTINE Print_Out_Input_Data ( XYZ, ID, INOD, NGP, MTEL, ID_BC, PMAT )
USE PARAMETERS
IMPLICIT DOUBLE PRECISION (A-H,O-Z)


! in
!---------- ---------- ---------- ---------- ----------   
DIMENSION XYZ(NJ, NDIM), ID(NJ, NDOF), INOD(NNODE, NEL), NGP(NEL), MTEL(NEL), ID_BC(NEL, NDIM**2), PMAT(NMAT, NPM)
!---------- ---------- ---------- ---------- ---------- 


! 1. geometry
DO I = 1, NJ
   WRITE(2,'(I10, 3F10.2)') I, (XYZ(I,J), J = 1, NDIM)
END DO

! 2. equation numbers
DO I = 1, NJ
   WRITE(2,'(10I10)') I, (ID(I,J), J = 1, NDOF)
END DO

! 3. element connectivity
SELECT CASE (NDIM)

   CASE(2)
      DO I = 1, NEL
         WRITE(2,'(32I10)') I, (INOD(J,I), J = 1, NNODE), NGP(I), MTEL(I), (ID_BC(I,K), K = 1, NDIM)
      END DO

   CASE(3)
      DO I = 1, NEL
         WRITE(2,'(I10, 32I10)') I, (INOD(J,I), J = 1, 8), NGP(I), MTEL(I), (ID_BC(I,K), K = 1, 4)
         IF ( NNODE == 20 .OR. NNODE == 27 ) WRITE(2,'(10X, 32I10)') (INOD(J,I), J = 9, 20)
         IF ( NNODE == 27 )                  WRITE(2,'(10X, 32I10)') (INOD(J,I), J = 21,27)
      END DO

END SELECT

! 4. material properties
DO I = 1, NMAT
   WRITE(2,'(I10, E10.5, 4F10.2)') I, (PMAT(I,J), J = 1, NPM)
END DO


RETURN
END


!************************************************
! Assign equation numbers according to the restraints ID
!************************************************
! REVISION : Tue, 7 August 2012

SUBROUTINE Assign_Equation_Numbers ( ID )
USE PARAMETERS
IMPLICIT DOUBLE PRECISION (A-H,O-Z)


! in
!---------- ---------- ---------- ---------- ----------
DIMENSION ID(NJ,NDOF)
!---------- ---------- ---------- ---------- ----------


IFLAG_NUMBERING = 1
NEQ = 0


SELECT CASE (IFLAG_NUMBERING)

   CASE(1)
! old way of numbering: displacement and stress dofs are not separated; j goes from 1 to ndof
          DO I = 1, NJ
             DO J = 1, NDOF
                 IF (ID(I,J).EQ.1) THEN
                     ID(I,J) = 0
                 ELSE
                     NEQ = NEQ + 1
                     ID(I,J) = NEQ
                 END IF
             END DO
          END DO


   CASE(2)
! alternative way of numbering: displacement dofs are numbered first, followed by stress dofs
          DO I = 1, NJ
             DO J = 1, NDIM
                 IF (ID(I,J).EQ.1) THEN
                     ID(I,J) = 0
                 ELSE
                     NEQ = NEQ + 1
                     ID(I,J) = NEQ
                 END IF
             END DO
          END DO

! size of displacement dofs
          NEQ_DISP = NEQ


          DO I = 1, NJ
             DO J = NDIM + 1, NDOF
                 IF (ID(I,J).EQ.1) THEN
                     ID(I,J) = 0
                 ELSE
                     NEQ = NEQ + 1
                     ID(I,J) = NEQ
                 END IF
             END DO
          END DO

   CASE(3)
! all seperated: first x, then y, sxx, syy, sxy
          DO J = 1, NDOF
             DO I = 1, NJ
                 IF (ID(I,J).EQ.1) THEN
                     ID(I,J) = 0
                 ELSE
                     NEQ = NEQ + 1
                     ID(I,J) = NEQ
                 END IF
             END DO
          END DO


END SELECT


! print out the final dof numbering
          DO I = 1, NJ
              WRITE(96,'(7I10)') I, (ID(I,J), J = 1, NDOF)
          END DO


RETURN
END

  
!************************************************
! IMPORT EARTHQUAKE RECORD FROM DATA FILE
!************************************************
! REVISION : Th, 17 May 2011

SUBROUTINE INPUT4 ( BACL )
USE PARAMETERS
IMPLICIT DOUBLE PRECISION (A-H,O-Z)
  

! out	  
!---------- ---------- ---------- ---------- ----------	
DIMENSION BACL(NSTEP, NDIM)
!---------- ---------- ---------- ---------- ----------


! local	  
!---------- ---------- ---------- ---------- ----------	
DIMENSION BACL0(NASTEP, NDIM+1) 
CHARACTER*100 HELP
!---------- ---------- ---------- ---------- ---------- 


! read acceleration time-history. first column is time.	  
!---------- ---------- ---------- ---------- ----------	
DO ISTEP = 1, NASTEP
   READ (1,*  ) (BACL0(ISTEP, J), J = 1, NDIM+1)
!  WRITE(2,'(4F10.5)') (BACL0(ISTEP, J), J = 1, NDIM+1)
END DO


! interpolate the missing parts  
!---------- ---------- ---------- ---------- ----------	
ICOL = 0
DO ISTEP = 1, NASTEP-1
   T1 = BACL0(ISTEP  , 1)
   T2 = BACL0(ISTEP+1, 1)
   DELTAT = T2 - T1
   N = DELTAT / DT + 0.10d0 * DT
   TBAR = (-1.0d0) * DT

   DO I = 1, N
      TBAR = TBAR + DT
      F1 = 1.0d0 - TBAR / DELTAT
      F2 =      TBAR / DELTAT
      ICOL = ICOL + 1
  
      DO J = 2, NDIM+1
         IF (ICOL.LE.NSTEP) THEN
            BACL(ICOL, J-1) = F1 * BACL0(ISTEP, J) + F2 * BACL0(ISTEP+1, J)
         END IF
      END DO
   END DO
END DO

!DO ISTEP = 1, NSTEP
!   WRITE(2,'(4F10.5)') (ISTEP-1)*DT, (BACL(ISTEP,J), J = 1, NDIM)
!END DO

! THIS IS DUE TO DR LOTFI'S PROGRAMMING STYLE (to avoid matrix factorization for the initial conditions).
! CHECK NEWMARK SUBROUTINE.
DO J = 1, NDIM
   BACL(1,J) = 0.0d0
END DO

BACL = BACL * GRAVITY

RETURN
END


!************************************************
! IMPORT PML PARAMETERS
!************************************************
! REVISION : Th, 17 May 2012

SUBROUTINE INPUT_PML ( PML_PARAM )
USE PARAMETERS
IMPLICIT DOUBLE PRECISION (A-H,O-Z)


! out	  
!---------- ---------- ---------- ---------- ----------	  
DIMENSION PML_PARAM(NDIM * 2 , 4)
!---------- ---------- ---------- ---------- ----------	 


! local	  
!---------- ---------- ---------- ---------- ----------	  
CHARACTER*100 HELP
!---------- ---------- ---------- ---------- ----------	

  
READ (1,*) HELP
WRITE(2,*) HELP
DO I = 1 , 2 * NDIM
   READ (1,*)              II , (PML_PARAM( I , J ) , J = 1 , 4)
   WRITE(2,'(I10,4F10.2)') II , (PML_PARAM( I , J ) , J = 1 , 4)
END DO
  
RETURN
END


!************************************************
! IMPORT OUTPUT JOINTS
!************************************************
! REVISION : Th, 17 May 2012

SUBROUTINE INPUT_LTRANS ( LTRANS )
USE PARAMETERS
IMPLICIT DOUBLE PRECISION (A-H,O-Z)


! out	  
!---------- ---------- ---------- ---------- ----------	  
DIMENSION LTRANS(NTRANS)
!---------- ---------- ---------- ---------- ----------	


DO I = 1, NTRANS
   READ(1,*) II, LTRANS(I)
   WRITE(2,'(2I10)') I, LTRANS (I)
END DO
      
WRITE(3,'(2A10,20I15)') 'Time', 'ENERGY', ( LTRANS(I), LTRANS(I), I = 1, NTRANS )
      
RETURN
END


!************************************************
! LOAD TIME FUNCTION
! receives t and returns f(t)
!************************************************
! REVISION : Sun 5 May 2013

SUBROUTINE TIME_FUNCTION ( T, F )
IMPLICIT DOUBLE PRECISION (A-H,O-Z)


! define parameters
!---------- ---------- ---------- ---------- ----------	 
pi     = 3.14159265358979323846264338327950288419716939937510d0

! Ricker pulse central frequency       
!fr     = 15.0d0                                  ! Sezgin uses 15 Hz in WM paper and 4 Hz in CMA paper
fr     = 0.120d0 ! Ricker pulse for inversion 

! - --------------------------------*----------------------------------------------------------------------------------------------------------------
! Gaussian pulse mean value (mu) and standard deviation (sigma)
! - --------------------------------*----------------------------------------------------------------------------------------------------------------
! 1- specify the maximum frequency of the Gaussian here
L_Gaussian = 40

! 2- selection of mu and sigma
Select Case ( L_Gaussian )
   Case ( 3  ) ; G_MEAN = 0.70d0 ; G_SDEV = 0.07000d0 ; ! t_end = 1.4
   Case ( 5  ) ; G_MEAN = 0.40d0 ; G_SDEV = 0.02600d0 ; ! t_end = 0.8
   Case ( 6  ) ; G_MEAN = 0.36d0 ; G_SDEV = 0.01700d0 ; ! t_end = 0.7
   Case ( 10 ) ; G_MEAN = 0.22d0 ; G_SDEV = 0.00580d0 ; ! t_end = 0.4
   Case ( 20 ) ; G_MEAN = 0.11d0 ; G_SDEV = 0.00140d0 ; ! t_end = 0.2
   Case ( 30 ) ; G_MEAN = 0.08d0 ; G_SDEV = 0.00070d0 ; ! t_end = 0.15
   Case ( 40 ) ; G_MEAN = 0.06d0 ; G_SDEV = 0.00040d0 ; ! t_end = 0.12
End Select

!G_MEAN = 0.150d0
!G_SDEV = 0.00060d0
! - --------------------------------*----------------------------------------------------------------------------------------------------------------
      
! Chopra's cosine wave
t_d    = 10.0d0
w_f    = 3.00d0

! delta function
k      = 10
t0     = 0.0005d0

!---------- ---------- ---------- ---------- ----------
! temporal loading case options
!        1   simple sine
!        2   Ricker pulse
!        3   Gaussian pulse
!        4   PML 3D: Chopra's cosine wave      f(t)    
!        5   PML 3D: Chopra's cosine wave d_dt f(t)
!        6   delta function (Greenberg, p. 12)
!        7   derivative of Ricker pulse
!        8   f(t) = 1.0
!        9   f(t) = 0.0
!       10   f(t) = Smooth Heaviside
LCASE  = 3
!---------- ---------- ---------- ---------- ----------


SELECT CASE (LCASE)

      
   CASE(1)
! simple sine	  
!---------- ---------- ---------- ---------- ----------	  
      F = DSIN( 2.0d0 * pi * fr * T )
!---------- ---------- ---------- ---------- ----------
      
      
   CASE(2)
! Ricker pulse	  
!---------- ---------- ---------- ---------- ----------
      wr    = 2.0d0 * pi * fr
      t_max = 6.0d0 * DSQRT(6.0d0) / wr
      u     = wr * T - 3.0d0 * DSQRT(6.0d0)
     
      IF ( T < t_max ) THEN
         FAC1  = ( 0.25d0 * u**2 - 0.5d0 ) * DEXP( -0.25d0 * u**2 ) - 13.0d0 * DEXP( -13.5d0 )
         FAC2  = 0.50d0 + 13.0d0 * DEXP( -13.5d0 )
         F     = FAC1 / FAC2
      ELSE
         F     = 0.0d0
      END IF
!---------- ---------- ---------- ---------- ----------

      
   CASE(3)
! Gaussian pulse	  
!---------- ---------- ---------- ---------- ----------
      F = DEXP( -( T - G_MEAN ) * ( T - G_MEAN ) / G_SDEV )
!---------- ---------- ---------- ---------- ----------


   CASE(4)
! Chopra's cosine wave f(t)
!---------- ---------- ---------- ---------- ----------
      T_f   = 2.0d0 * pi / w_f                   ! dominant forcing period
      nc    = t_d / T_f - 0.50d0                 ! number of full cycles
      T_f   = t_d / ( DBLE(nc) + 0.50d0 )        ! adjust for consistency

      IF ( T < T_f / 2.0d0 ) THEN
         F = 0.5d0 * ( 1.0d0 - dcos( 2.0d0 * pi * T / T_f ) )
      ELSEIF ( T >= T_f / 2.0d0 .AND. T < DBLE(nc) * T_f ) THEN
         F = dcos( 2.0d0 * pi * ( T - T_f / 2.0d0 ) / T_f )
      ELSEIF ( T >= DBLE(nc) * T_f .AND. T <= t_d ) THEN
         F = 0.50d0 * ( 1.0d0 - dcos( 2.0d0 * pi * ( T - DBLE(nc) * T_f ) / T_f ) ) - 1.0d0
      ELSE
         F = 0.0d0
      END IF
!---------- ---------- ---------- ---------- ----------


   CASE(5)
! Chopra's cosine wave d_dt f(t)
!---------- ---------- ---------- ---------- ----------
      T_f   = 2.0d0 * pi / w_f                   ! dominant forcing period
      nc    = t_d / T_f - 0.50d0                 ! number of full cycles
      T_f   = t_d / ( DBLE(nc) + 0.50d0 )        ! adjust for consistency

      IF ( T < T_f / 2.0d0 ) THEN
         F = ( pi / T_f ) * dsin( 2.0d0 * pi * T / T_f )
      ELSEIF ( T >= T_f / 2.0d0 .AND. T < DBLE(nc) * T_f ) THEN
         F = ( -2.0d0 * pi / T_f ) * dsin( 2.0d0 * pi * ( T - T_f / 2.0d0 ) / T_f )
      ELSEIF ( T >= DBLE(nc) * T_f .AND. T <= t_d ) THEN
         F = ( pi / T_f ) * dsin( 2.0d0 * pi * ( T - DBLE(nc) * T_f ) / T_f )
      ELSE
         F = 0.0d0
      END IF
!---------- ---------- ---------- ---------- ----------


   CASE(6)
! Greenberg's delta function d(t-t0)
!---------- ---------- ---------- ---------- ----------
      F = DBLE(k) / ( pi * ( 1.0d0 + ( DBLE(k) * (T - t0) )**2 ) )
!---------- ---------- ---------- ---------- ----------


   CASE(7)
! derivative of Ricker pulse    
!---------- ---------- ---------- ---------- ----------
      wr    = 2.0d0 * pi * fr
      t_max = 6.0d0 * DSQRT(6.0d0) / wr
      u     = wr * T - 3.0d0 * DSQRT(6.0d0)
      up    = wr                                 ! d/dt u(t)

      IF ( T < t_max ) THEN
         FAC1  = 0.25d0 * 2.0d0 * up * u * DEXP( -0.25d0 * u**2 ) + ( 0.25d0 * u**2 - 0.5d0 ) * ( -0.25d0 * 2.0d0 * up * u) * DEXP( -0.25d0 * u**2 )
         FAC2  = 0.50d0 + 13.0d0 * DEXP( -13.5d0 )
         F     = FAC1 / FAC2
      ELSE
         F     = 0.0d0
      END IF
!---------- ---------- ---------- ---------- ----------


   CASE(8)
! f(t) = 1.0
!---------- ---------- ---------- ---------- ----------   
      F = 1.0d0
!---------- ---------- ---------- ---------- ----------


   CASE(9)
! f(t) = 0.0
!---------- ---------- ---------- ---------- ----------   
      F = 0.0d0
!---------- ---------- ---------- ---------- ----------


   CASE(10)
! smooth Heaviside
!---------- ---------- ---------- ---------- ----------   
      if ( T <= t0 ) then
         F = T/t0
      else
         F = 1.0d0
      end if
!---------- ---------- ---------- ---------- ----------


   CASE DEFAULT
      WRITE(*,*) 'ERROR IN TIME_FUNCTION'
      STOP
      
END SELECT
      

RETURN
END


!************************************************
! PRINT OUT TIME-HISTORY FOR SELECTED JOINTS
!************************************************
! REVISION : F, 18 May 2012      

SUBROUTINE OUT_NEWMARK_PETSC ( U_PETSC, ID, LTRANS, T, ENERGY )
USE PARAMETERS
IMPLICIT DOUBLE PRECISION (A-H,O-Z)


!==================================================================================================
! - PETSC INCLUDE FILES
#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"


Vec                     :: U_PETSC
PetscScalar             :: ENERGY

PetscScalar, pointer    :: U(:)
PetscInt                :: Istart, Iend
!==================================================================================================


! in  
!---------- ---------- ---------- ---------- ----------   
DIMENSION ID(NJ,NDOF), LTRANS(NTRANS)
!---------- ---------- ---------- ---------- ---------- 


! prints out to screen
!---------- ---------- ---------- ---------- ----------
! CALL VecView ( U_PETSC, PETSC_VIEWER_STDOUT_SELF, IERR )


! start: first local element
! end  : one more than the last local element
!---------- ---------- ---------- ---------- ----------
CALL VecGetOwnershipRange ( U_PETSC, Istart, Iend, IERR )


! transform start / end to Fortran way. Now, Iend is the last one
Istart = Istart + 1 


! print out data related to the first entry of LTRANS
!---------- ---------- ---------- ---------- ----------
ID_X_GLOBAL = ID( LTRANS(1) , 1 )
ID_Y_GLOBAL = ID( LTRANS(1) , 2 )

IF ( NDIM == 3 ) ID_Z_GLOBAL = ID( LTRANS(1) , 3 )

IF ( ID_X_GLOBAL >= Istart .AND. ID_X_GLOBAL <= Iend ) THEN

   CALL VecGetArrayF90     ( U_PETSC, U, IERR )

   ID_X_LOCAL = ID_X_GLOBAL - Istart + 1
   ID_Y_LOCAL = ID_Y_GLOBAL - Istart + 1         ! a dirty implementation, fix it later
   IF ( NDIM == 3 ) ID_Z_LOCAL = ID_Z_GLOBAL - Istart + 1

! this format reserves 15 spaces, and prints 3 digits of the exponent
   SELECT CASE (NDIM)

      CASE(2)
         WRITE(3,'(F10.5,20E15.5E3)') T, ENERGY, U(ID_X_LOCAL), U(ID_Y_LOCAL)

      CASE(3)
         WRITE(3,'(F10.5,20E15.5E3)') T, ENERGY, U(ID_X_LOCAL), U(ID_Y_LOCAL), U(ID_Z_LOCAL)

   END SELECT

   CALL VecRestoreArrayF90 ( U_PETSC, U, IERR )

END IF


RETURN
END


!************************************************
! Print out Petsc objects for debugging purposes
!************************************************
! REVISION : MONDAY, 16 JULY 2012      

SUBROUTINE PETSC_Vec_Debug ( Vec_PETSC, T )
USE PARAMETERS
IMPLICIT DOUBLE PRECISION (A-H,O-Z)


!==================================================================================================
! - PETSC INCLUDE FILES
#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"


Vec                     :: Vec_PETSC

PetscScalar, pointer    :: Vec_Fortran(:)
PetscInt                :: Istart, Iend
!==================================================================================================


CALL VecGetOwnershipRange ( Vec_PETSC, Istart, Iend, IERR )

!transform start / end to Fortran way. Now, Iend is the last one
!Istart = Istart + 1
!write(*,'(A10,I10,A10,I10)') "Istart =", Istart, "Iend =", Iend


CALL VecGetArrayF90     ( Vec_PETSC, Vec_Fortran, IERR )

!WRITE(*,'(F10.5,20E10.2)') T, Vec_Fortran(1)
WRITE(95,'(15F10.0,F9.0)') Vec_Fortran

CALL VecRestoreArrayF90 ( Vec_PETSC, Vec_Fortran, IERR )


RETURN
END 
