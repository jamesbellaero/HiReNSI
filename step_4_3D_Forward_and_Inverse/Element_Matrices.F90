!************************************************
! SOLID ELEMENT EK MATRIX
!************************************************
! REVISION : Th, 17 May 2012

SUBROUTINE ELEMENT_K ( XYZ, INOD, NGP, MTEL, PMAT, PMat_Lambda, PMat_Mu, IEL, EK )
USE PARAMETERS
IMPLICIT DOUBLE PRECISION (A-H,O-Z)
  

! in 
!---------- ---------- ---------- ---------- ----------
DIMENSION XYZ(NJ,NDIM), INOD(NNODE,NEL), NGP(NEL), MTEL(NEL), PMAT(NMAT,NPM)
Dimension PMat_Lambda ( NJ ), PMat_Mu ( NJ )
!---------- ---------- ---------- ---------- ----------


! out
!---------- ---------- ---------- ---------- ----------
DIMENSION EK(NNODE * NDIM, NNODE * NDIM)
!---------- ---------- ---------- ---------- ----------


! local
!---------- ---------- ---------- ---------- ----------
DIMENSION XT(NDIM,NNODE), DJ(NDIM,NDIM), DJI(NDIM,NDIM), FN(NNODE), DFX(NNODE,NDIM), DFXI(NNODE,NDIM)
DIMENSION XINT(9,9), WINT(9,9)
       
DIMENSION H1(NNODE, NNODE), H2(NNODE, NNODE), H3(NNODE, NNODE), H4(NNODE, NNODE)
DIMENSION G1(NNODE, NNODE), G2(NNODE, NNODE), G3(NNODE, NNODE), G4(NNODE, NNODE)           

Dimension zLambda_e(NNODE), zMu_e(NNODE)
!---------- ---------- ---------- ---------- ----------


CALL GAUSS(XINT,WINT)


! load material properties
!---------- ---------- ---------- ---------- ----------
MTYPE   = MTEL(IEL)
NINT    = NGP(IEL)
E       = PMAT(MTYPE,1)
ANU     = PMAT(MTYPE,2)
RO      = PMAT(MTYPE,3)
THICK   = PMAT(MTYPE,5)
! homogeneous material property
Zlambda = E * ANU / ( (1.0d0 + ANU) * (1.0d0 - 2.0d0 * ANU) )
Zmu     = E / ( 2.0d0*(1.0d0 + ANU) )


! define factors and initialize (Homogeneous material property)
!---------- ---------- ---------- ---------- ----------	
FAC1 = zLAMBDA + 2.0d0 * zMU
FAC2 = zLAMBDA
FAC3 = zMU
!---------- ---------- ---------- ---------- ----------
EK = 0.0d0


DO I = 1, NNODE
   K = INOD(I,IEL)
   DO J = 1, NDIM
      XT(J,I) = XYZ(K,J)
   END DO
   zLambda_e (I) = PMat_Lambda ( K )
   zMu_e     (I) = PMat_Mu ( K )
END DO
  
  
! Gauss integration	   
!---------- ---------- ---------- ---------- ----------
DO LY = 1, NINT
   X2 = XINT(LY,NINT)
   WY = WINT(LY,NINT)

   DO LX = 1, NINT
      X1 = XINT(LX,NINT)
      WX = WINT(LX,NINT)

      WSTAR = WX * WY
         
         
! Jacobian and shape functions	  
!---------- ---------- ---------- ---------- ----------
      SELECT CASE (NNODE)

         CASE (8)                                ! 8-noded elements
            CALL SHAPE8(X1, X2, FN, DFXI)

         CASE (9)                                ! 9-noded elements
            CALL SHAPE9(X1, X2, FN, DFXI)

         CASE DEFAULT
            WRITE(*,'("FLAG ERROR IN EM - NNODE")')
            STOP

      END SELECT


      DJ = MATMUL(XT,DFXI)

      DETJ = DJ(1,1) * DJ(2,2) - DJ(1,2) * DJ(2,1)

      IF (DETJ.LT.0) THEN
         WRITE(*,*)'|J|<0 NEGATIVE DET. JACOBIAN'
         READ(*,*)
         STOP
      END IF

      DJI(1,1) =  DJ(2,2) / DETJ
      DJI(2,2) =  DJ(1,1) / DETJ
      DJI(1,2) = -DJ(1,2) / DETJ
      DJI(2,1) = -DJ(2,1) / DETJ

      DFX  = MATMUL(DFXI,DJI)
!---------- ---------- ---------- ---------- ----------	
      FAC0 = WSTAR * DETJ


! define factors (Heterogeneous material property)
!---------- ---------- ---------- ---------- ----------	
      zLAMBDA = Dot_Product ( FN , zLambda_e )
      zMU     = Dot_Product ( FN , zMu_e )
      FAC1 = zLAMBDA + 2.0d0 * zMU
      FAC2 = zLAMBDA
      FAC3 = zMU
!---------- ---------- ---------- ---------- ----------


! form the matrix
!---------- ---------- ---------- ---------- ----------
      DO I = 1 , NNODE
         DO J = 1 , NNODE
            H1(I,J) = DFX(I,1) * DFX(J,1)
            H2(I,J) = DFX(I,2) * DFX(J,2)
            H3(I,J) = DFX(I,1) * DFX(J,2)
            H4(I,J) = DFX(I,2) * DFX(J,1)
         END DO
      END DO
      
      G1 = FAC1 * H1 + FAC3 * H2
      G2 = FAC2 * H3 + FAC3 * H4
      G3 = FAC2 * H4 + FAC3 * H3
      G4 = FAC1 * H2 + FAC3 * H1
      
      DO I = 1 , NNODE
         I1 = I
         I2 = I1 + NNODE
      
         DO J = 1 , NNODE
            J1 = J
            J2 = J1 + NNODE
        
            EK(I1,J1) = EK(I1,J1) + G1(I,J) * FAC0
            EK(I1,J2) = EK(I1,J2) + G2(I,J) * FAC0
            EK(I2,J1) = EK(I2,J1) + G3(I,J) * FAC0
            EK(I2,J2) = EK(I2,J2) + G4(I,J) * FAC0
         END DO
      END DO
!---------- ---------- ---------- ---------- ----------


   END DO
END DO
!---------- ---------- ---------- ---------- ----------


RETURN
END


!************************************************
! SOLID ELEMENT EM MATRIX
!************************************************
! REVISION : M 24 Sep 2012

SUBROUTINE ELEMENT_M ( XYZ, INOD, NGP, MTEL, PMAT, PML_PARAM, IEL, EM, IFLAG )
USE PARAMETERS
IMPLICIT DOUBLE PRECISION (A-H,O-Z)
  

! in 
!---------- ---------- ---------- ---------- ----------
DIMENSION XYZ(NJ,NDIM), INOD(NNODE,NEL), NGP(NEL), MTEL(NEL), PMAT(NMAT,NPM)
DIMENSION PML_PARAM(NDIM * 2 , 4)
!---------- ---------- ---------- ---------- ----------


! out
!---------- ---------- ---------- ---------- ----------
DIMENSION EM(NNODE * NDIM, NNODE * NDIM)
!---------- ---------- ---------- ---------- ----------


! local
!---------- ---------- ---------- ---------- ----------
DIMENSION XT(NDIM,NNODE), DJ(NDIM,NDIM), FN(NNODE), DFXI(NNODE,NDIM)
DIMENSION XINT(9,9), WINT(9,9)
  
DIMENSION H1(NNODE, NNODE)
DIMENSION G1(NNODE, NNODE)
DIMENSION XG(NDIM), PML_ALPHA_BETA(19)
!---------- ---------- ---------- ---------- ----------	  


SELECT CASE (IFLAG)

   CASE (201, 202, 203, 204)                     ! M1, Ma, Mb, Mc
      CALL GAUSS_Lobatto ( XINT,WINT )
      NINT  = NINT_Lobatto

   CASE DEFAULT
      WRITE(*,'("FLAG ERROR IN EM")')
      STOP

END SELECT


! load material properties
!---------- ---------- ---------- ---------- ----------	  
MTYPE = MTEL(IEL)
RO    = PMAT(MTYPE,3)
THICK = PMAT(MTYPE,5)
  
  
! define factors and initialize	 	  
!---------- ---------- ---------- ---------- ----------	  
EM   = 0.0d0


DO I = 1, NNODE
   K = INOD(I,IEL)
   DO J = 1, NDIM
      XT(J,I) = XYZ(K,J)
   END DO
END DO


! Gauss integration	   
!---------- ---------- ---------- ---------- ----------	  
DO LY = 1, NINT
   X2 = XINT(LY,NINT)
   WY = WINT(LY,NINT)
  
   DO LX = 1, NINT
      X1 = XINT(LX,NINT)
      WX = WINT(LX,NINT)
  
      WSTAR = WX * WY
  
  
! Jacobian and shape functions	  
!---------- ---------- ---------- ---------- ----------		  
      SELECT CASE (NNODE)

         CASE (8)                                ! 8-noded elements
            CALL SHAPE8(X1, X2, FN, DFXI)

         CASE (9)                                ! 9-noded elements
            CALL SHAPE9(X1, X2, FN, DFXI)

         CASE DEFAULT
            WRITE(*,'("FLAG ERROR IN EM - NNODE")')
            STOP

      END SELECT


      DJ = MATMUL(XT,DFXI)
   
      DETJ = DJ(1,1) * DJ(2,2) - DJ(1,2) * DJ(2,1)
      IF (DETJ.LT.0) THEN
         WRITE(*,*)'|J|<0 NEGATIVE DET. JACOBIAN'
         READ (*,*)
         STOP
      END IF
!---------- ---------- ---------- ---------- ----------	
      FAC0 = WSTAR * DETJ


! define factors	 	  
!---------- ---------- ---------- ---------- ----------	
      SELECT CASE (IFLAG)
             
         CASE (202, 203, 204)
            XG = MATMUL(XT, FN)
            CALL PML_FINDER (XG , PML_PARAM , PML_ALPHA_BETA)         
     
      END SELECT          
      
     
      SELECT CASE (IFLAG)
      
         CASE (201)                              ! M1: mass matrix, regular domain
            FAC1 = RO
    
         CASE (202)                              ! Ma: mass matrix, PML
            FAC1 = RO * PML_ALPHA_BETA(5)                  

         CASE (203)                              ! Mb: mass matrix, PML
            FAC1 = RO * PML_ALPHA_BETA(6)
      
         CASE (204)                              ! Mc: mass matrix, PML
            FAC1 = RO * PML_ALPHA_BETA(7)
             
         CASE DEFAULT
            WRITE(*,'("FLAG ERROR IN EM")')
            STOP
            
      END SELECT  
!---------- ---------- ---------- ---------- ----------


idebug = 9
if (idebug == 95 .and. iflag == 202) then
    write(95,'(i10,7f10.5)') iel, pml_alpha_beta
end if

  
! form element matrix	  
!---------- ---------- ---------- ---------- ----------
      DO I = 1 , NNODE
         DO J = 1 , NNODE
            H1(I,J) = FN(I) * FN(J)
         END DO
      END DO
      
      G1 = FAC1 * H1
      
      DO I = 1 , NNODE
         I1 = I
         I2 = I1 + NNODE
      
         DO J = 1 , NNODE
            J1 = J
            J2 = J1 + NNODE
          
            EM(I1,J1) = EM(I1,J1) + G1(I,J) * FAC0
            EM(I2,J2) = EM(I2,J2) + G1(I,J) * FAC0
         END DO
      END DO
!---------- ---------- ---------- ---------- ----------		 

 
   END DO
END DO
!---------- ---------- ---------- ---------- ----------

  
RETURN
END

  
!************************************************
! PML ELEMENT EN MATRIX
! flag   case
! 300    EN (elastodynamics with mixed FEM formulation)
! 301    ENa
! 302    ENb
! 303    ENc
!************************************************
! REVISION : W 8 May 2013

SUBROUTINE ELEMENT_EN ( XYZ, INOD, NGP, MTEL, PMAT, PML_PARAM, IEL, EN, IFLAG )
USE PARAMETERS
IMPLICIT DOUBLE PRECISION (A-H,O-Z)


! in 
!---------- ---------- ---------- ---------- ----------
DIMENSION XYZ(NJ,NDIM), INOD(NNODE,NEL), NGP(NEL), MTEL(NEL), PMAT(NMAT,NPM)
DIMENSION PML_PARAM(NDIM * 2 , 4)
!---------- ---------- ---------- ---------- ----------


! out
!---------- ---------- ---------- ---------- ----------
DIMENSION EN(NNODE * 3, NNODE * 3)
!---------- ---------- ---------- ---------- ----------


! local
!---------- ---------- ---------- ---------- ----------
DIMENSION XT(NDIM,NNODE), FN(NNODE), DJ(NDIM,NDIM), DFXI(NNODE,NDIM)
DIMENSION XINT(9,9), WINT(9,9)
          
DIMENSION H1(NNODE, NNODE)
DIMENSION G1(NNODE, NNODE), G2(NNODE, NNODE)
    
DIMENSION XG(NDIM), PML_ALPHA_BETA(19)
!---------- ---------- ---------- ---------- ----------

  
SELECT CASE (IFLAG)

   CASE (300, 301, 302, 303)                     ! EN, ENa, ENb, ENc
      CALL GAUSS_Lobatto ( XINT,WINT )
      NINT  = NINT_Lobatto

   CASE DEFAULT
      WRITE(*,'("FLAG ERROR IN EN")')
      STOP

END SELECT


! load material properties
!---------- ---------- ---------- ---------- ----------	
MTYPE   = MTEL(IEL)
E       = PMAT(MTYPE,1)
ANU     = PMAT(MTYPE,2)
RO      = PMAT(MTYPE,3)
THICK   = PMAT(MTYPE,5)


! define factors and initialize. move them down later as they become distributed parameters
!---------- ---------- ---------- ---------- ----------
FAC1 = 2.0D0
!---------- ---------- ---------- ---------- ----------
EN = 0.0d0


DO I = 1, NNODE
   K = INOD(I,IEL)
   DO J = 1, NDIM
      XT(J,I) = XYZ(K,J)
   END DO
END DO


! Gauss quadrature
!---------- ---------- ---------- ---------- ----------
DO LY = 1, NINT
   X2 = XINT(LY,NINT)
   WY = WINT(LY,NINT)

   DO LX = 1, NINT
      X1 = XINT(LX,NINT)
      WX = WINT(LX,NINT)

      WSTAR = WX * WY
  
  
! Jacobian and shape functions	  
!---------- ---------- ---------- ---------- ----------		
      SELECT CASE (NNODE)

         CASE (8)                                ! 8-noded elements
            CALL SHAPE8 ( X1, X2, FN, DFXI )

         CASE (9)                                ! 9-noded elements
            CALL SHAPE9 ( X1, X2, FN, DFXI )

         CASE DEFAULT
            WRITE(*,'("FLAG ERROR IN EM - NNODE")')
            STOP

      END SELECT


      DJ = MATMUL(XT,DFXI)

      DETJ = DJ(1,1) * DJ(2,2) - DJ(1,2) * DJ(2,1)
      IF (DETJ.LT.0) THEN
         WRITE(*,*)'|J|<0 NEGATIVE DET. JACOBIAN'
         READ(*,*)
         STOP
      END IF
!---------- ---------- ---------- ---------- ----------	
      FAC0 = WSTAR * DETJ
         

! define factors	 	  
!---------- ---------- ---------- ---------- ----------	
      XG = MATMUL(XT, FN)
      CALL PML_FINDER ( XG , PML_PARAM , PML_ALPHA_BETA )
      
      SELECT CASE (IFLAG)

         CASE (300)
            FAC2 = 1.0d0                         ! EN
   
         CASE (301)                              ! ENa
            FAC2 = PML_ALPHA_BETA(5)         
    
         CASE (302)                              ! ENb
            FAC2 = PML_ALPHA_BETA(6)         

         CASE (303)                              ! ENc
            FAC2 = PML_ALPHA_BETA(7)         
              
         CASE DEFAULT
            WRITE(*,'("FLAG ERROR IN EA")')
            STOP
      
      END SELECT  
!---------- ---------- ---------- ---------- ----------	


! form element matrix
!---------- ---------- ---------- ---------- ----------
      DO I = 1 , NNODE
         DO J = 1 , NNODE
            H1(I,J) = FN(I) * FN(J)
         END DO
      END DO
      
      G1 =        FAC2 * H1
      G2 = FAC1 * FAC2 * H1
      
      DO I = 1 , NNODE
         I1 = I
         I2 = I1 + NNODE
         I3 = I2 + NNODE
      
         DO J = 1 , NNODE
            J1 = J
            J2 = J1 + NNODE
            J3 = J2 + NNODE
          
            EN(I1,J1) = EN(I1,J1) + G1(I,J) * FAC0
            EN(I2,J2) = EN(I2,J2) + G1(I,J) * FAC0
            EN(I3,J3) = EN(I3,J3) + G2(I,J) * FAC0
         END DO
      END DO
!---------- ---------- ---------- ---------- ----------


   END DO
END DO
!---------- ---------- ---------- ---------- ----------


RETURN
END

  
!************************************************
! PML ELEMENT EA MATRIX

! IFLAG == 101 evanescent  wave => EAe_upper 
! IFLAG == 102 propagaring wave => EAp_upper
! IFLAG == 103 evanescent  wave => EAe_lower
! IFLAG == 104 propagaring wave => EAp_lower
!          110 elastodynamics, mixed FEM upper
!          111 elastodynamics, mixed FEM lower
!************************************************
! REVISION : W 8 May 2013

SUBROUTINE ELEMENT_EA ( XYZ, INOD, NGP, MTEL, PMAT, PMat_Lambda, PMat_Mu, PML_PARAM, IEL, EA, IFLAG )
USE PARAMETERS
IMPLICIT DOUBLE PRECISION (A-H,O-Z)


! in 
!---------- ---------- ---------- ---------- ----------
DIMENSION XYZ(NJ,NDIM), INOD(NNODE,NEL), NGP(NEL), MTEL(NEL), PMAT(NMAT,NPM)
DIMENSION PML_PARAM(NDIM * 2 , 4)
Dimension PMat_Lambda (NJ), PMat_Mu(NJ)
!---------- ---------- ---------- ---------- ----------


! out 
!---------- ---------- ---------- ---------- ----------
DIMENSION EA(NNODE * 2, NNODE * 3) 
!---------- ---------- ---------- ---------- ----------


! local
!---------- ---------- ---------- ---------- ----------
DIMENSION XT(NDIM,NNODE), DJ(NDIM,NDIM), DJI(NDIM,NDIM), FN(NNODE), DFX(NNODE,NDIM), DFXI(NNODE,NDIM)
DIMENSION XINT(9,9), WINT(9,9)
       
DIMENSION H1(NNODE, NNODE), H2(NNODE, NNODE), H3(NNODE, NNODE)
DIMENSION G1(NNODE, NNODE), G2(NNODE, NNODE), G3(NNODE, NNODE), G4(NNODE, NNODE), G5(NNODE, NNODE), G6(NNODE, NNODE)

DIMENSION XG(NDIM), PML_ALPHA_BETA(19)    

Dimension zLambda_e(NNODE), zMu_e(NNODE) 
!---------- ---------- ---------- ---------- ----------


CALL GAUSS(XINT,WINT)


! load material properties
!---------- ---------- ---------- ---------- ----------	
MTYPE = MTEL(IEL)
NINT  = NGP(IEL)
E     = PMAT(MTYPE,1)
ANU   = PMAT(MTYPE,2)
RO    = PMAT(MTYPE,3)
THICK = PMAT(MTYPE,5)
! homogeneous material property
Zlambda = E * ANU / ( (1.0d0 + ANU)*(1.0d0 - 2.0d0 * ANU) )
Zmu     = E / ( 2.0d0 * (1.0d0 + ANU) )

        
! define factors and initialize. move them down later as they become distributed parameters
!---------- ---------- ---------- ---------- ---------- 
SELECT CASE (IFLAG)

   CASE (101, 102, 110)                          ! upper block
      FAC5 = 1.0d0
      FAC6 = 0.0d0
      FAC7 = 1.0d0
      FAC8 = 1.0d0
      FAC9 = 0.0d0

   CASE (103, 104, 111)                          ! lower block
      FAC5 = Zlambda + 2.0d0 * Zmu
      FAC6 = Zlambda
      FAC7 =           2.0d0 * Zmu
      FAC8 = 0.0d0
      FAC9 = 1.0d0
  
   CASE DEFAULT
      WRITE(*,'("FLAG ERROR IN EA")')
      STOP

END SELECT
!---------- ---------- ---------- ---------- ----------
EA = 0.0d0


DO I = 1, NNODE
   K = INOD(I,IEL)
   DO J = 1, NDIM
      XT(J,I) = XYZ(K,J)
   END DO
   zLambda_e (I) = PMat_Lambda ( K )
   zMu_e     (I) = PMat_Mu ( K )
END DO
  

! Gauss quadrature
!---------- ---------- ---------- ---------- ----------
DO LY = 1, NINT
   X2 = XINT(LY,NINT)
   WY = WINT(LY,NINT)

   DO LX = 1, NINT
      X1 = XINT(LX,NINT)
      WX = WINT(LX,NINT)

      WSTAR = WX * WY


! Jacobian and shape functions	  
!---------- ---------- ---------- ---------- ----------	
      SELECT CASE (NNODE)

         CASE (8)                                ! 8-noded elements
            CALL SHAPE8(X1, X2, FN, DFXI)

         CASE (9)                                ! 9-noded elements
            CALL SHAPE9(X1, X2, FN, DFXI)

         CASE DEFAULT
            WRITE(*,'("FLAG ERROR IN EM - NNODE")')
            STOP

      END SELECT


      DJ = MATMUL(XT,DFXI)

      DETJ = DJ(1,1) * DJ(2,2) - DJ(1,2) * DJ(2,1)
      IF (DETJ.LT.0) THEN
         WRITE(*,*)'|J|<0 NEGATIVE DET. JACOBIAN'
         READ(*,*)
         STOP
      END IF

      DJI(1,1)= DJ(2,2)/DETJ
      DJI(2,2)= DJ(1,1)/DETJ
      DJI(1,2)=-DJ(1,2)/DETJ
      DJI(2,1)=-DJ(2,1)/DETJ

      DFX = MATMUL(DFXI,DJI)
!---------- ---------- ---------- ---------- ----------
      FAC0 = WSTAR * DETJ


! define factors (Heterogeneous material property)
!---------- ---------- ---------- ---------- ----------	
      zLAMBDA = Dot_Product ( FN , zLambda_e )
      zMU     = Dot_Product ( FN , zMu_e )
      SELECT CASE (IFLAG)
         CASE (103, 104, 111)                          ! lower block
            FAC5 = Zlambda + 2.0d0 * Zmu
            FAC6 = Zlambda
            FAC7 =           2.0d0 * Zmu
      End Select
!---------- ---------- ---------- ---------- ----------


! define factors	 	  
!---------- ---------- ---------- ---------- ----------	
! [ alpha_x , alpha_y , beta_x , beta_y , a , b , c, d/dx alpha_y, d/dy alpha_x, d/dx beta_y, d/dy beta_x ]
      XG = MATMUL(XT, FN)
      CALL PML_FINDER ( XG , PML_PARAM , PML_ALPHA_BETA )
      
      SELECT CASE (IFLAG)

         CASE (101, 103)                         ! evanescent wave
            FAC1  = PML_ALPHA_BETA(1)            ! alpha_x
            FAC2  = PML_ALPHA_BETA(2)            ! alpha_y
            FAC3  = PML_ALPHA_BETA(8)            ! d/dx alpha_y    (interpretation: use ratio for computation of the other  axis)
            FAC4  = PML_ALPHA_BETA(9)            ! d/dy alpha_x

            FAC11 = PML_ALPHA_BETA(12)           ! alpha_x_yt (honest mpml terms)
            FAC12 = PML_ALPHA_BETA(13)           ! alpha_y_xt
            FAC13 = PML_ALPHA_BETA(16)           ! d/dx alpha_x_yt
            FAC14 = PML_ALPHA_BETA(17)           ! d/dy alpha_y_xt
    
         CASE (102, 104)                         ! propagating wave
            FAC1  = PML_ALPHA_BETA(3)            ! beta_x
            FAC2  = PML_ALPHA_BETA(4)            ! beta_y
            FAC3  = PML_ALPHA_BETA(10)           ! d/dx beta_y
            FAC4  = PML_ALPHA_BETA(11)           ! d/dy beta_x
  
            FAC11 = PML_ALPHA_BETA(14)           ! beta_x_yt (honest mpml)
            FAC12 = PML_ALPHA_BETA(15)           ! beta_y_xt
            FAC13 = PML_ALPHA_BETA(18)           ! d/dx beta_x_yt
            FAC14 = PML_ALPHA_BETA(19)           ! d/dy beta_y_xt

         CASE (110, 111)                         ! elastodynamics, mixed FEM
            FAC1 = 1.0d0
            FAC2 = 1.0d0
            FAC3 = 0.0d0
            FAC4 = 0.0d0

            FAC11 = 0.0d0
            FAC12 = 0.0d0
            FAC13 = 0.0d0
            FAC14 = 0.0d0

         CASE DEFAULT
            WRITE(*,'("FLAG ERROR IN EA")')
            STOP
      
      END SELECT
!---------- ---------- ---------- ---------- ----------


! form element matrix
!---------- ---------- ---------- ---------- ----------
      DO I = 1 , NNODE
         DO J = 1 , NNODE
            H1(I,J) = DFX(I,1) * FN(J)
            H2(I,J) = DFX(I,2) * FN(J)
            H3(I,J) =    FN(I) * FN(J)
         END DO
      END DO
      
      honesty = 0.0d0                            ! use 1.0 for mpml with corrected derivatives, or 0.0 for Papageorgiou's mpml

      G1 = FAC5 * FAC2 * H1 + FAC3 * FAC8 * H3   +  ( FAC5 * FAC12 * H2 + FAC8 * FAC14 * H3 ) * honesty
      G2 = FAC6 * FAC2 * H1 * FAC9               +  ( FAC6 * FAC12 * H2                     ) * honesty
      G3 = FAC7 * FAC1 * H2 + FAC4 * FAC8 * H3   +  ( FAC7 * FAC11 * H1 + FAC8 * FAC13 * H3 ) * honesty
      G4 = FAC6 * FAC1 * H2 * FAC9               +  ( FAC6 * FAC11 * H1                     ) * honesty
      G5 = FAC5 * FAC1 * H2 + FAC4 * FAC8 * H3   +  ( FAC5 * FAC11 * H1 + FAC8 * FAC13 * H3 ) * honesty
      G6 = FAC7 * FAC2 * H1 + FAC3 * FAC8 * H3   +  ( FAC7 * FAC12 * H2 + FAC8 * FAC14 * H3 ) * honesty
      
      DO I = 1 , NNODE
         I1 = I
         I2 = I1 + NNODE
      
         DO J = 1 , NNODE
            J1 = J
            J2 = J1 + NNODE
            J3 = J2 + NNODE
           
            EA(I1,J1) = EA(I1,J1) + G1(I,J) * FAC0
            EA(I1,J2) = EA(I1,J2) + G2(I,J) * FAC0
            EA(I1,J3) = EA(I1,J3) + G3(I,J) * FAC0
            EA(I2,J1) = EA(I2,J1) + G4(I,J) * FAC0
            EA(I2,J2) = EA(I2,J2) + G5(I,J) * FAC0
            EA(I2,J3) = EA(I2,J3) + G6(I,J) * FAC0
         END DO
      END DO
!---------- ---------- ---------- ---------- ----------


   END DO
END DO
!---------- ---------- ---------- ---------- ----------


RETURN
END


!************************************************
! REGULAR DOMAIN surface tractions
!************************************************
! REVISION : F 26 Oct 2012

SUBROUTINE ELEMENT_FS ( XYZ, INOD, NGP, MTEL, PMAT, ID_BC, IEL, EFS )
USE PARAMETERS
IMPLICIT DOUBLE PRECISION (A-H,O-Z)


! in 
!---------- ---------- ---------- ---------- ----------
DIMENSION XYZ(NJ,NDIM), INOD(NNODE,NEL), NGP(NEL), MTEL(NEL), PMAT(NMAT,NPM), ID_BC(NEL,NDIM**2)
!---------- ---------- ---------- ---------- ----------


! out
!---------- ---------- ---------- ---------- ----------
DIMENSION EFS(NNODE * NDIM)
!---------- ---------- ---------- ---------- ----------


! local
!---------- ---------- ---------- ---------- ----------
DIMENSION XT(NDIM,NNODE), DJ(NDIM,NDIM), FN(NNODE), DFXI(NNODE,NDIM)
DIMENSION LFACE(4), FVAL(4), IPERM(2), XI(NDIM), VN(NDIM)
DIMENSION XINT(9,9), WINT(9,9)
  
DIMENSION H1(NNODE)
DIMENSION G1(NNODE), G2(NNODE)

! for face 3, 4 we are integrating along xi; face 1, 2 along eta       
DATA LFACE / 2, 2, 1, 1/
DATA FVAL  /+1.0d0, -1.0d0, +1.0d0, -1.0d0/
DATA IPERM / 2, 1/
!---------- ---------- ---------- ---------- ----------	  
  
 
CALL GAUSS(XINT,WINT)


! load material properties
!---------- ---------- ---------- ---------- ----------	  
MTYPE = MTEL(IEL)
NINT  = NGP(IEL)
THICK = PMAT(MTYPE,5)
  
  
! define factors and initialize: only Ty component
!---------- ---------- ---------- ---------- ----------	  
FAC1 = 0.0d0                                     ! Tx = 0 
FAC2 = 1.0d0                                     ! Ty = 1
!---------- ---------- ---------- ---------- ----------	  
EFS  = 0.0d0
      
      
NFACE = ID_BC(IEL, 1)
IF ( NFACE == 0 .OR. NFACE == 11 ) RETURN         ! 11 stands for disc load
!WRITE(*,'(A30, I10)')'SURFACE TRACTION'


! correct normal direction according to the face number
!---------- ---------- ---------- ---------- ----------	
SELECT CASE (NFACE)
      
   CASE (1)
      FACLAN = -1.0d0
         
   CASE (2)
      FACLAN =  1.0d0
                  
   CASE (3)
      FACLAN =  1.0d0
                  
   CASE (4)
      FACLAN = -1.0d0
                 
   CASE DEFAULT
      WRITE(*,*) 'FACE ERROR IN ELEMENT_FN'
      STOP
      
END SELECT
!---------- ---------- ---------- ---------- ----------	     

      
DO I = 1, NNODE
   K = INOD(I,IEL)
   DO J = 1, NDIM
      XT(J,I) = XYZ(K,J)
   END DO
END DO
  
  
L     = LFACE (NFACE)                            ! principal direction of integration
M     = IPERM (L)                                ! perpendicular to principal direction, hence constant
XI(M) = FVAL  (NFACE)                            ! constant Gauss value at boundary
  
  
! Gauss integration on boundary	   
!---------- ---------- ---------- ---------- ----------	  
DO LX = 1, NINT
   XI(L) = XINT(LX,NINT)
   WX    = WINT(LX,NINT)


! Jacobian and shape functions
! VN is the outward normal vector to the boundary	  
!---------- ---------- ---------- ---------- ----------		  
   SELECT CASE (NNODE)

      CASE (8)                                   ! 8-noded elements
         CALL SHAPE8(XI(1), XI(2), FN, DFXI)

      CASE (9)                                   ! 9-noded elements
         CALL SHAPE9(XI(1), XI(2), FN, DFXI)

      CASE DEFAULT
         WRITE(*,'("FLAG ERROR IN EM - NNODE")')
         STOP

   END SELECT


   DJ = MATMUL(XT,DFXI)
      
   VN(1) = -DJ(2,L) * FACLAN
   VN(2) =  DJ(1,L) * FACLAN
     
   DETJ_SURF = DSQRT( VN(1)**2 + VN(2)**2 )
!---------- ---------- ---------- ---------- ----------	
   FAC0 = WX * DETJ_SURF * THICK
    
    
idebug = 9
if (idebug == 96) then
    write(*,'(2i5,2f10.3)') iel, nface, vn(1)/DETJ_SURF, vn(2)/DETJ_SURF
end if
         
  
! form surface force vector	  
!---------- ---------- ---------- ---------- ----------
   DO I = 1 , NNODE
      H1(I) = FN(I)
   END DO
      
   G1 = FAC1 * H1
   G2 = FAC2 * H1
      
   DO I = 1 , NNODE
      I1 = I
      I2 = I1 + NNODE
          
      EFS(I1) = EFS(I1) + G1(I) * FAC0
      EFS(I2) = EFS(I2) + G2(I) * FAC0
   END DO
!---------- ---------- ---------- ---------- ----------		 


END DO
!---------- ---------- ---------- ---------- ----------


RETURN
END


!************************************************
! REGULAR DOMAIN body force (disc spatial component)
!************************************************
! REVISION : W 6 March 2013

SUBROUTINE ELEMENT_FB_disc ( XYZ, INOD, NGP, MTEL, PMAT, ID_BC, IEL, EFB )
USE PARAMETERS
IMPLICIT DOUBLE PRECISION (A-H,O-Z)


! in
!---------- ---------- ---------- ---------- ----------
DIMENSION XYZ(NJ,NDIM), INOD(NNODE,NEL), NGP(NEL), MTEL(NEL), PMAT(NMAT,NPM), ID_BC(NEL,NDIM**2)
!---------- ---------- ---------- ---------- ----------


! out
!---------- ---------- ---------- ---------- ----------
DIMENSION EFB(NNODE * NDIM)
!---------- ---------- ---------- ---------- ----------


! local
!---------- ---------- ---------- ---------- ----------
DIMENSION XT(NDIM,NNODE), DJ(NDIM,NDIM), FN(NNODE), DFXI(NNODE,NDIM)
DIMENSION XINT(9,9), WINT(9,9)

DIMENSION H1(NNODE)
DIMENSION G1(NNODE), G2(NNODE)
DIMENSION XG(NDIM)
!---------- ---------- ---------- ---------- ----------


CALL GAUSS ( XINT, WINT )


! load material properties
!---------- ---------- ---------- ---------- ----------
MTYPE = MTEL(IEL)
NINT  =  NGP(IEL)
THICK = PMAT(MTYPE,5)


! define disc parameters and initialize
!---------- ---------- ---------- ---------- ----------
! WM paper
x_c  =    0.0d0                                  ! disc center
y_c  = -125.0d0
r_d  =    2.5d0                                  ! disc radius

! CMA paper
!x_c  =  -10.0d0                                  ! disc center
!y_c  =   -4.0d0
!r_d  =    0.4d0                                  ! disc radius

EFB  =    0.0d0


NFACE = ID_BC(IEL, 1)
IF ( NFACE /= 11 ) RETURN                        ! NFACE = 11 stands for disc load
WRITE(*,'(A30, I10)')'DISC LOAD'


DO I = 1, NNODE
   K = INOD(I,IEL)
   DO J = 1, NDIM
      XT(J,I) = XYZ(K,J)
   END DO
END DO


! Gauss integration        
!---------- ---------- ---------- ---------- ----------   
DO LY = 1, NINT
   X2 = XINT(LY,NINT)
   WY = WINT(LY,NINT)

   DO LX = 1, NINT
      X1 = XINT(LX,NINT)
      WX = WINT(LX,NINT)

      WSTAR = WX * WY


! Jacobian and shape functions    
!---------- ---------- ---------- ---------- ----------           
      SELECT CASE (NNODE)

         CASE (8)                                ! 8-noded elements
            CALL SHAPE8(X1, X2, FN, DFXI)

         CASE (9)                                ! 9-noded elements
            CALL SHAPE9(X1, X2, FN, DFXI)

         CASE DEFAULT
            WRITE(*,'("FLAG ERROR IN FB_disc - NNODE")')
            STOP

      END SELECT


      DJ = MATMUL(XT,DFXI)

      DETJ = DJ(1,1) * DJ(2,2) - DJ(1,2) * DJ(2,1)
      IF (DETJ.LE.0) THEN
         WRITE(*,*)'|J|<0 NEGATIVE DET. JACOBIAN'
         READ (*,*)
         STOP
      END IF
!---------- ---------- ---------- ---------- ---------- 
      FAC0 = WSTAR * DETJ * THICK


! define and integrate body force
!---------- ---------- ---------- ---------- ----------
      XG   = MATMUL(XT, FN)
      r    = DSQRT(  ( XG(1) - x_c )**2 + ( XG(2) - y_c )**2  )

      IF ( r > 1.001d0 * r_d  ) THEN
         WRITE (*,*) 'Error in disc load'
         STOP
      END IF

      FAC1 = ( 1.0d0 - ( r/r_d )**2   )**3 * ( XG(1) - x_c ) / r
      FAC2 = ( 1.0d0 - ( r/r_d )**2   )**3 * ( XG(2) - y_c ) / r


! form disc force vector       
!---------- ---------- ---------- ---------- ----------
      DO I = 1 , NNODE
         H1(I) = FN(I)
      END DO

      G1 = FAC1 * H1
      G2 = FAC2 * H1

      DO I = 1 , NNODE
         I1 = I
         I2 = I1 + NNODE

         EFB(I1) = EFB(I1) + G1(I) * FAC0
         EFB(I2) = EFB(I2) + G2(I) * FAC0
      END DO
!---------- ---------- ---------- ---------- ----------       


   END DO
END DO
!---------- ---------- ---------- ---------- ----------


RETURN
END
