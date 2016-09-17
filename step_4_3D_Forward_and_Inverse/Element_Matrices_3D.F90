!************************************************
! SOLID ELEMENT EK MATRIX
!************************************************
! REVISION : Th, 24 May 2012
! REVISION : Tu, 31 March 2015

SUBROUTINE ELEMENT_K_3D ( XYZ, INOD, NGP, MTEL, PMAT, PMat_Lambda, PMat_Mu, IEL, EK )
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

DIMENSION FN_8(8)



DIMENSION H1(NNODE, NNODE), H2(NNODE, NNODE), H3(NNODE, NNODE), H4(NNODE, NNODE), H5(NNODE, NNODE), H6(NNODE, NNODE), H7(NNODE, NNODE), H8(NNODE, NNODE), &
          H9(NNODE, NNODE)
DIMENSION G1(NNODE, NNODE), G2(NNODE, NNODE), G3(NNODE, NNODE), G4(NNODE, NNODE), G5(NNODE, NNODE), G6(NNODE, NNODE), G7(NNODE, NNODE), G8(NNODE, NNODE), &
          G9(NNODE, NNODE)

Dimension zLambda_e(NNODE), zMu_e(NNODE)
!---------- ---------- ---------- ---------- ----------


CALL GAUSS ( XINT, WINT )


! load material properties
!---------- ---------- ---------- ---------- ----------
MTYPE   = MTEL(IEL)
NINT    = NGP(IEL)
E       = PMAT(MTYPE,1)
ANU     = PMAT(MTYPE,2)
RO      = PMAT(MTYPE,3)
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
DO LZ = 1, NINT
   X3 = XINT(LZ,NINT)
   WZ = WINT(LZ,NINT)

   DO LY = 1, NINT
      X2 = XINT(LY,NINT)
      WY = WINT(LY,NINT)

      DO LX = 1, NINT
         X1 = XINT(LX,NINT)
         WX = WINT(LX,NINT)

         WSTAR = WX * WY * WZ


! Jacobian and shape functions    
!---------- ---------- ---------- ---------- ----------         
         SELECT CASE (NNODE)

            CASE (20)                            ! 20-noded elements (serendipity)
               CALL SHAPE_3_20 (X1, X2, X3, FN, DFXI)

            CASE (27)                            ! 27-noded elements (Lagrange)
               CALL SHAPE_3_27 (X1, X2, X3, FN, DFXI)

            CASE DEFAULT
               WRITE(*,'("FLAG ERROR IN EK_3D - NNODE")')
               STOP

         END SELECT


         DJ = MATMUL(XT,DFXI)

         DETJ = DJ(1,1) * DJ(2,2) * DJ(3,3) + DJ(2,1) * DJ(3,2) * DJ(1,3) + DJ(3,1) * DJ(1,2) * DJ(2,3) - &
                DJ(1,3) * DJ(2,2) * DJ(3,1) - DJ(2,3) * DJ(3,2) * DJ(1,1) - DJ(3,3) * DJ(1,2) * DJ(2,1)

         IF (DETJ.LE.0) THEN
            WRITE(*,*)'|J|<0 NEGATIVE DET. JACOBIAN'
            READ(*,*)
            STOP
         END IF

         DJI(1,1) = + DJ(3,3) * DJ(2,2) - DJ(3,2) * DJ(2,3)
         DJI(1,2) = - DJ(3,3) * DJ(1,2) + DJ(3,2) * DJ(1,3)
         DJI(1,3) = + DJ(2,3) * DJ(1,2) - DJ(2,2) * DJ(1,3)
 
         DJI(2,1) = - DJ(3,3) * DJ(2,1) + DJ(3,1) * DJ(2,3)
         DJI(2,2) = + DJ(3,3) * DJ(1,1) - DJ(3,1) * DJ(1,3)
         DJI(2,3) = - DJ(2,3) * DJ(1,1) + DJ(2,1) * DJ(1,3)

         DJI(3,1) = + DJ(3,2) * DJ(2,1) - DJ(3,1) * DJ(2,2)
         DJI(3,2) = - DJ(3,2) * DJ(1,1) + DJ(3,1) * DJ(1,2)
         DJI(3,3) = + DJ(2,2) * DJ(1,1) - DJ(2,1) * DJ(1,2)

         DJI = DJI / DETJ
              

         DFX = MATMUL(DFXI,DJI)
!---------- ---------- ---------- ---------- ---------- 
         FAC0 = WSTAR * DETJ


! define factors (Heterogeneous material property)
!---------- ---------- ---------- ---------- ----------	
! option 1: use 27 noded elements
      zLAMBDA = Dot_Product ( FN , zLambda_e )
      zMU     = Dot_Product ( FN , zMu_e )
!/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-
! option 2: use 8 noded elements for material interpolation -- inversion crime run to address the reviewers
! for unknown reasons, this dosn't work. The energy plot looks unreasonable.
!      Call SHAPE_3_8 ( X1, X2, X3, FN_8 )
!      zLAMBDA = Dot_Product ( FN_8 , zLambda_e(1:8) )
!      zMU     = Dot_Product ( FN_8 , zMu_e(1:8) )
!/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-/-
! compute factors
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
               H3(I,J) = DFX(I,3) * DFX(J,3)
               H4(I,J) = DFX(I,1) * DFX(J,2)
               H5(I,J) = DFX(I,2) * DFX(J,1)
               H6(I,J) = DFX(I,1) * DFX(J,3)
               H7(I,J) = DFX(I,3) * DFX(J,1)
               H8(I,J) = DFX(I,2) * DFX(J,3)
               H9(I,J) = DFX(I,3) * DFX(J,2)
            END DO
         END DO

         G1 = FAC1 * H1 + FAC3 * ( H2 + H3 )
         G2 = FAC2 * H4 + FAC3 * H5
         G3 = FAC2 * H6 + FAC3 * H7
         G4 = FAC2 * H5 + FAC3 * H4
         G5 = FAC1 * H2 + FAC3 * ( H1 + H3 )
         G6 = FAC2 * H8 + FAC3 * H9
         G7 = FAC2 * H7 + FAC3 * H6
         G8 = FAC2 * H9 + FAC3 * H8
         G9 = FAC1 * H3 + FAC3 * ( H1 + H2 )

         DO I = 1 , NNODE
            I1 = I
            I2 = I1 + NNODE
            I3 = I2 + NNODE

            DO J = 1 , NNODE
               J1 = J
               J2 = J1 + NNODE
               J3 = J2 + NNODE
        
               EK(I1,J1) = EK(I1,J1) + G1(I,J) * FAC0
               EK(I1,J2) = EK(I1,J2) + G2(I,J) * FAC0
               EK(I1,J3) = EK(I1,J3) + G3(I,J) * FAC0

               EK(I2,J1) = EK(I2,J1) + G4(I,J) * FAC0
               EK(I2,J2) = EK(I2,J2) + G5(I,J) * FAC0
               EK(I2,J3) = EK(I2,J3) + G6(I,J) * FAC0

               EK(I3,J1) = EK(I3,J1) + G7(I,J) * FAC0                
               EK(I3,J2) = EK(I3,J2) + G8(I,J) * FAC0
               EK(I3,J3) = EK(I3,J3) + G9(I,J) * FAC0

            END DO
         END DO
!---------- ---------- ---------- ---------- ----------


      END DO
   END DO
END DO
!---------- ---------- ---------- ---------- ----------

RETURN
END


!************************************************
! SOLID ELEMENT EM MATRIX
!************************************************
! REVISION : Th, 24 May 2012

SUBROUTINE ELEMENT_M_3D ( XYZ, INOD, NGP, MTEL, PMAT, PML_PARAM, IEL, EM, IFLAG )
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

   CASE (201, 202, 203, 204, 205)                ! M1, Ma, Mb, Mc, Md
      CALL GAUSS_Lobatto ( XINT, WINT )
      NINT = NINT_Lobatto

   CASE DEFAULT
      WRITE(*,'("FLAG ERROR IN EM_3D")')
      STOP

END SELECT


! load material properties
!---------- ---------- ---------- ---------- ----------   
MTYPE = MTEL(IEL)
RO    = PMAT(MTYPE,3)


! define factors and initialize           
!---------- ---------- ---------- ---------- ----------   
EM    = 0.0d0


DO I = 1, NNODE
   K = INOD(I,IEL)
   DO J = 1, NDIM
      XT(J,I) = XYZ(K,J)
   END DO
END DO


! Gauss integration        
!---------- ---------- ---------- ---------- ----------   
DO LZ = 1, NINT
   X3 = XINT(LZ,NINT)
   WZ = WINT(LZ,NINT)

   DO LY = 1, NINT
      X2 = XINT(LY,NINT)
      WY = WINT(LY,NINT)

      DO LX = 1, NINT
         X1 = XINT(LX,NINT)
         WX = WINT(LX,NINT)

         WSTAR = WX * WY * WZ


! Jacobian and shape functions    
!---------- ---------- ---------- ---------- ----------           
         SELECT CASE (NNODE)

            CASE (20)                            ! 20-noded elements (serendipity)
               CALL SHAPE_3_20 (X1, X2, X3, FN, DFXI)

            CASE (27)                            ! 27-noded elements (Lagrange)
               CALL SHAPE_3_27 (X1, X2, X3, FN, DFXI)

            CASE DEFAULT
               WRITE(*,'("FLAG ERROR IN EM_3D - NNODE")')
               STOP

         END SELECT


         DJ = MATMUL(XT,DFXI)

         DETJ = DJ(1,1) * DJ(2,2) * DJ(3,3) + DJ(2,1) * DJ(3,2) * DJ(1,3) + DJ(3,1) * DJ(1,2) * DJ(2,3) - &
                DJ(1,3) * DJ(2,2) * DJ(3,1) - DJ(2,3) * DJ(3,2) * DJ(1,1) - DJ(3,3) * DJ(1,2) * DJ(2,1)

         IF (DETJ.LE.0) THEN
            WRITE(*,*)'|J|<0 NEGATIVE DET. JACOBIAN'
            READ(*,*)
            STOP
         END IF
!---------- ---------- ---------- ---------- ---------- 
         FAC0 = WSTAR * DETJ


! define factors                  
!---------- ---------- ---------- ---------- ---------- 
         SELECT CASE (IFLAG)

            CASE (202, 203, 204, 205)
               XG = MATMUL(XT, FN)
               CALL PML_FINDER_3D ( XG, PML_PARAM, PML_ALPHA_BETA )

         END SELECT


         SELECT CASE (IFLAG)

            CASE (201)                           ! M1: mass matrix, regular domain
               FAC1 = RO

            CASE (202)                           ! Ma: mass matrix, PML
               FAC1 = RO * PML_ALPHA_BETA(7)

            CASE (203)                           ! Mb: mass matrix, PML
               FAC1 = RO * PML_ALPHA_BETA(8)

            CASE (204)                           ! Mc: mass matrix, PML
               FAC1 = RO * PML_ALPHA_BETA(9)

            CASE (205)                           ! Md: mass matrix, PML
               FAC1 = RO * PML_ALPHA_BETA(10)

            CASE DEFAULT
               WRITE(*,'("FLAG ERROR IN EM_3D")')
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
            I3 = I2 + NNODE

            DO J = 1 , NNODE
               J1 = J
               J2 = J1 + NNODE
               J3 = J2 + NNODE

               EM(I1,J1) = EM(I1,J1) + G1(I,J) * FAC0
               EM(I2,J2) = EM(I2,J2) + G1(I,J) * FAC0
               EM(I3,J3) = EM(I3,J3) + G1(I,J) * FAC0
        
            END DO
         END DO
!---------- ---------- ---------- ---------- ----------          


      END DO
   END DO
END DO
!---------- ---------- ---------- ---------- ----------


RETURN
END


!************************************************
! PML ELEMENT EN MATRIX
!************************************************
! REVISION : Th, 24 May 2012

SUBROUTINE ELEMENT_EN_3D ( XYZ, INOD, NGP, MTEL, PMAT, PML_PARAM, IEL, EN, IFLAG )
USE PARAMETERS
IMPLICIT DOUBLE PRECISION (A-H,O-Z)


! in 
!---------- ---------- ---------- ---------- ----------
DIMENSION XYZ(NJ,NDIM), INOD(NNODE,NEL), NGP(NEL), MTEL(NEL), PMAT(NMAT,NPM)
DIMENSION PML_PARAM(NDIM * 2 , 4)
!---------- ---------- ---------- ---------- ----------


! out
!---------- ---------- ---------- ---------- ----------
DIMENSION EN(NNODE * 6, NNODE * 6)
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

   CASE (301, 302, 303, 304)                     ! ENa, ENb, ENc, ENd
      CALL GAUSS_Lobatto ( XINT, WINT )
      NINT = NINT_Lobatto

   CASE DEFAULT
      WRITE(*,'("FLAG ERROR IN EN_3D")')
      STOP

END SELECT


! define factors and initialize.
!---------- ---------- ---------- ---------- ----------
FAC1 = 2.0d0
!---------- ---------- ---------- ---------- ----------
EN   = 0.0d0


DO I = 1, NNODE
   K = INOD(I,IEL)
   DO J = 1, NDIM
      XT(J,I) = XYZ(K,J)
   END DO
END DO


! Gauss integration        
!---------- ---------- ---------- ---------- ----------   
DO LZ = 1, NINT
   X3 = XINT(LZ,NINT)
   WZ = WINT(LZ,NINT)

   DO LY = 1, NINT
      X2 = XINT(LY,NINT)
      WY = WINT(LY,NINT)

      DO LX = 1, NINT
         X1 = XINT(LX,NINT)
         WX = WINT(LX,NINT)

         WSTAR = WX * WY * WZ


! Jacobian and shape functions    
!---------- ---------- ---------- ---------- ----------           
         SELECT CASE (NNODE)

            CASE (20)                            ! 20-noded elements (serendipity)
               CALL SHAPE_3_20 (X1, X2, X3, FN, DFXI)

            CASE (27)                            ! 27-noded elements (Lagrange)
               CALL SHAPE_3_27 (X1, X2, X3, FN, DFXI)

            CASE DEFAULT
               WRITE(*,'("FLAG ERROR IN EN_3D - NNODE")')
               STOP

         END SELECT


         DJ = MATMUL(XT,DFXI)

         DETJ = DJ(1,1) * DJ(2,2) * DJ(3,3) + DJ(2,1) * DJ(3,2) * DJ(1,3) + DJ(3,1) * DJ(1,2) * DJ(2,3) - &
                DJ(1,3) * DJ(2,2) * DJ(3,1) - DJ(2,3) * DJ(3,2) * DJ(1,1) - DJ(3,3) * DJ(1,2) * DJ(2,1)

         IF (DETJ.LE.0) THEN
            WRITE(*,*)'|J|<0 NEGATIVE DET. JACOBIAN'
            READ(*,*)
            STOP
         END IF
!---------- ---------- ---------- ---------- ---------- 
         FAC0 = WSTAR * DETJ


! define factors                  
!---------- ---------- ---------- ---------- ---------- 
         XG = MATMUL(XT, FN)
         CALL PML_FINDER_3D ( XG, PML_PARAM, PML_ALPHA_BETA )

         SELECT CASE (IFLAG)

            CASE (301)                           ! ENa
               FAC2 = PML_ALPHA_BETA(7)

            CASE (302)                           ! ENb
               FAC2 = PML_ALPHA_BETA(8)

            CASE (303)                           ! ENc
               FAC2 = PML_ALPHA_BETA(9)

            CASE (304)                           ! ENd
               FAC2 = PML_ALPHA_BETA(10)

            CASE DEFAULT
               WRITE(*,'("FLAG ERROR IN EN_3D")')
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
            I4 = I3 + NNODE
            I5 = I4 + NNODE
            I6 = I5 + NNODE

            DO J = 1 , NNODE
               J1 = J
               J2 = J1 + NNODE
               J3 = J2 + NNODE
               J4 = J3 + NNODE
               J5 = J4 + NNODE
               J6 = J5 + NNODE
        
               EN(I1,J1) = EN(I1,J1) + G1(I,J) * FAC0
               EN(I2,J2) = EN(I2,J2) + G1(I,J) * FAC0
               EN(I3,J3) = EN(I3,J3) + G1(I,J) * FAC0

               EN(I4,J4) = EN(I4,J4) + G2(I,J) * FAC0
               EN(I5,J5) = EN(I5,J5) + G2(I,J) * FAC0
               EN(I6,J6) = EN(I6,J6) + G2(I,J) * FAC0

            END DO
         END DO
!---------- ---------- ---------- ---------- ----------


      END DO
   END DO
END DO
!---------- ---------- ---------- ---------- ----------


RETURN
END


!************************************************
! PML ELEMENT EA MATRIX

! IFLAG == 101 evanescent  wave => EAe_upper 
! IFLAG == 102 propagaring wave => EAp_upper
! IFLAG == 103                  => EAw_upper
! IFLAG == 104 evanescent  wave => EAe_lower
! IFLAG == 105 propagaring wave => EAp_lower
! IFLAG == 106                  => EAw_lower
!************************************************
! REVISION : M, 11 Jun 2012

SUBROUTINE ELEMENT_EA_3D ( XYZ, INOD, NGP, MTEL, PMAT, PMat_Lambda, PMat_Mu, PML_PARAM, IEL, EA, IFLAG )
USE PARAMETERS
IMPLICIT DOUBLE PRECISION (A-H,O-Z)


! in 
!---------- ---------- ---------- ---------- ----------
DIMENSION XYZ(NJ,NDIM), INOD(NNODE,NEL), NGP(NEL), MTEL(NEL), PMAT(NMAT,NPM)
DIMENSION PML_PARAM(NDIM * 2 , 4)
Dimension PMat_Lambda ( NJ ), PMat_Mu ( NJ )
!---------- ---------- ---------- ---------- ----------


! out 
!---------- ---------- ---------- ---------- ----------
DIMENSION EA(NNODE * NDIM, NNODE * 6)
!---------- ---------- ---------- ---------- ----------


! local
!---------- ---------- ---------- ---------- ----------
DIMENSION XT(NDIM,NNODE), DJ(NDIM,NDIM), DJI(NDIM,NDIM), FN(NNODE), DFX(NNODE,NDIM), DFXI(NNODE,NDIM)
DIMENSION XINT(9,9), WINT(9,9)

DIMENSION H1(NNODE, NNODE), H2(NNODE, NNODE), H3(NNODE, NNODE), H4(NNODE, NNODE)
DIMENSION G1(NNODE, NNODE), G2(NNODE, NNODE), G3(NNODE, NNODE), G4(NNODE, NNODE), G5(NNODE, NNODE), G6(NNODE, NNODE), G7(NNODE, NNODE), G8(NNODE, NNODE), &
          G9(NNODE, NNODE)

DIMENSION XG(NDIM), PML_ALPHA_BETA(19)

Dimension zLambda_e(NNODE), zMu_e(NNODE)
!---------- ---------- ---------- ---------- ----------


CALL GAUSS ( XINT, WINT )


! load material properties
!---------- ---------- ---------- ---------- ---------- 
MTYPE   = MTEL(IEL)
NINT    = NGP(IEL)
E       = PMAT(MTYPE,1)
ANU     = PMAT(MTYPE,2)
! homogeneous material property
Zlambda = E * ANU / ( (1.0d0 + ANU) * (1.0d0 - 2.0d0 * ANU) )
Zmu     = E / ( 2.0d0 * (1.0d0 + ANU) )


! define factors and initialize. move them down later as they become distributed parameters
!---------- ---------- ---------- ---------- ---------- 
SELECT CASE (IFLAG)

   CASE (101, 102, 103)                          ! upper block
      FAC7  = 1.0d0
      FAC8  = 0.0d0
      FAC9  = 1.0d0
      FAC10 = 1.0d0

   CASE (104, 105, 106)                          ! lower block
      FAC7  = Zlambda + 2.0d0 * Zmu
      FAC8  = Zlambda
      FAC9  =           2.0d0 * Zmu
      FAC10 = 0.0d0

   CASE DEFAULT
      WRITE(*,'("FLAG ERROR IN EA_3D")')
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
        

! Gauss integration        
!---------- ---------- ---------- ---------- ----------
DO LZ = 1, NINT
   X3 = XINT(LZ,NINT)
   WZ = WINT(LZ,NINT)

   DO LY = 1, NINT
      X2 = XINT(LY,NINT)
      WY = WINT(LY,NINT)

      DO LX = 1, NINT
         X1 = XINT(LX,NINT)
         WX = WINT(LX,NINT)

         WSTAR = WX * WY * WZ


! Jacobian and shape functions    
!---------- ---------- ---------- ---------- ----------         
         SELECT CASE (NNODE)

            CASE (20)                            ! 20-noded elements (serendipity)
               CALL SHAPE_3_20 (X1, X2, X3, FN, DFXI)

            CASE (27)                            ! 27-noded elements (Lagrange)
               CALL SHAPE_3_27 (X1, X2, X3, FN, DFXI)

            CASE DEFAULT
               WRITE(*,'("FLAG ERROR IN EA_3D - NNODE")')
               STOP

         END SELECT


         DJ = MATMUL(XT,DFXI)

         DETJ = DJ(1,1) * DJ(2,2) * DJ(3,3) + DJ(2,1) * DJ(3,2) * DJ(1,3) + DJ(3,1) * DJ(1,2) * DJ(2,3) - &
                DJ(1,3) * DJ(2,2) * DJ(3,1) - DJ(2,3) * DJ(3,2) * DJ(1,1) - DJ(3,3) * DJ(1,2) * DJ(2,1)

         IF (DETJ.LE.0) THEN
            WRITE(*,*)'|J|<0 NEGATIVE DET. JACOBIAN'
            READ(*,*)
            STOP
         END IF

         DJI(1,1) = + DJ(3,3) * DJ(2,2) - DJ(3,2) * DJ(2,3)
         DJI(1,2) = - DJ(3,3) * DJ(1,2) + DJ(3,2) * DJ(1,3)
         DJI(1,3) = + DJ(2,3) * DJ(1,2) - DJ(2,2) * DJ(1,3)

         DJI(2,1) = - DJ(3,3) * DJ(2,1) + DJ(3,1) * DJ(2,3)
         DJI(2,2) = + DJ(3,3) * DJ(1,1) - DJ(3,1) * DJ(1,3)
         DJI(2,3) = - DJ(2,3) * DJ(1,1) + DJ(2,1) * DJ(1,3)

         DJI(3,1) = + DJ(3,2) * DJ(2,1) - DJ(3,1) * DJ(2,2)
         DJI(3,2) = - DJ(3,2) * DJ(1,1) + DJ(3,1) * DJ(1,2)
         DJI(3,3) = + DJ(2,2) * DJ(1,1) - DJ(2,1) * DJ(1,2)

         DJI = DJI / DETJ


         DFX  = MATMUL(DFXI,DJI)
!---------- ---------- ---------- ---------- ---------- 
         FAC0 = WSTAR * DETJ


! define factors (Heterogeneous material property)
!---------- ---------- ---------- ---------- ----------	
      zLAMBDA = Dot_Product ( FN , zLambda_e )
      zMU     = Dot_Product ( FN , zMu_e )
      SELECT CASE (IFLAG)
         CASE (104, 105, 106)                          ! lower block
            FAC7  = Zlambda + 2.0d0 * Zmu
            FAC8  = Zlambda
            FAC9  =           2.0d0 * Zmu
      End Select
!---------- ---------- ---------- ---------- ----------


! define factors                  
!---------- ---------- ---------- ---------- ---------- 
         XG = MATMUL(XT, FN)
         CALL PML_FINDER_3D ( XG , PML_PARAM , PML_ALPHA_BETA )

         SELECT CASE (IFLAG)

            CASE (101, 104)                                                                           ! e-wave
               FAC1 = PML_ALPHA_BETA(1) * PML_ALPHA_BETA(2)                                           !    alpha_x * alpha_y
               FAC2 = PML_ALPHA_BETA(1) * PML_ALPHA_BETA(3)                                           !    alpha_x * alpha_z
               FAC3 = PML_ALPHA_BETA(2) * PML_ALPHA_BETA(3)                                           !    alpha_y * alpha_z
               FAC4 = PML_ALPHA_BETA(11)                                                              ! d (alpha_x * alpha_y) 
               FAC5 = PML_ALPHA_BETA(12)                                                              ! d (alpha_x * alpha_z)
               FAC6 = PML_ALPHA_BETA(13)                                                              ! d (alpha_y * alpha_z)

            CASE (102, 105)                                                                           ! p-wave  
               FAC1 = PML_ALPHA_BETA(1) * PML_ALPHA_BETA(5) + PML_ALPHA_BETA(2) * PML_ALPHA_BETA(4)   !    alpha_x * beta_y + alpha_y * beta_x
               FAC2 = PML_ALPHA_BETA(1) * PML_ALPHA_BETA(6) + PML_ALPHA_BETA(3) * PML_ALPHA_BETA(4)   !    alpha_x * beta_z + alpha_z * beta_x        
               FAC3 = PML_ALPHA_BETA(2) * PML_ALPHA_BETA(6) + PML_ALPHA_BETA(3) * PML_ALPHA_BETA(5)   !    alpha_y * beta_z + alpha_z * beta_y         
               FAC4 = PML_ALPHA_BETA(14)                                                              ! d (alpha_x * beta_y + alpha_y * beta_x) 
               FAC5 = PML_ALPHA_BETA(15)                                                              ! d (alpha_x * beta_z + alpha_z * beta_x)
               FAC6 = PML_ALPHA_BETA(16)                                                              ! d (alpha_y * beta_z + alpha_z * beta_y)

            CASE (103, 106)                                                                           ! w-wave    
               FAC1 = PML_ALPHA_BETA(4) * PML_ALPHA_BETA(5)                                           !    beta_x * beta_y
               FAC2 = PML_ALPHA_BETA(4) * PML_ALPHA_BETA(6)                                           !    beta_x * beta_z
               FAC3 = PML_ALPHA_BETA(5) * PML_ALPHA_BETA(6)                                           !    beta_y * beta_z
               FAC4 = PML_ALPHA_BETA(17)                                                              ! d (beta_x * beta_y)
               FAC5 = PML_ALPHA_BETA(18)                                                              ! d (beta_x * beta_z)
               FAC6 = PML_ALPHA_BETA(19)                                                              ! d (beta_y * beta_z)

            CASE DEFAULT
               WRITE(*,'("FLAG ERROR IN EA_3D")')
               STOP

         END SELECT
!---------- ---------- ---------- ---------- ----------


! form element matrix
!---------- ---------- ---------- ---------- ----------
         DO I = 1 , NNODE
            DO J = 1 , NNODE
               H1(I,J) = DFX(I,1) * FN(J)
               H2(I,J) = DFX(I,2) * FN(J)
               H3(I,J) = DFX(I,3) * FN(J)
               H4(I,J) =  FN(I)   * FN(J)
            END DO
         END DO

         G1 = FAC7 * FAC3 * H1 + FAC10 * FAC6 * H4
         G2 = FAC8 * FAC3 * H1
         G3 = FAC9 * FAC2 * H2 + FAC10 * FAC5 * H4
         G4 = FAC9 * FAC1 * H3 + FAC10 * FAC4 * H4
         G5 = FAC8 * FAC2 * H2
         G6 = FAC7 * FAC2 * H2 + FAC10 * FAC5 * H4
         G7 = FAC9 * FAC3 * H1 + FAC10 * FAC6 * H4
         G8 = FAC8 * FAC1 * H3
         G9 = FAC7 * FAC1 * H3 + FAC10 * FAC4 * H4

         DO I = 1 , NNODE
            I1 = I
            I2 = I1 + NNODE
            I3 = I2 + NNODE

            DO J = 1 , NNODE
               J1 = J
               J2 = J1 + NNODE
               J3 = J2 + NNODE
               J4 = J3 + NNODE
               J5 = J4 + NNODE
               J6 = J5 + NNODE

               EA(I1,J1) = EA(I1,J1) + G1(I,J) * FAC0
               EA(I1,J2) = EA(I1,J2) + G2(I,J) * FAC0
               EA(I1,J3) = EA(I1,J3) + G2(I,J) * FAC0
               EA(I1,J4) = EA(I1,J4) + G3(I,J) * FAC0
               EA(I1,J5) = EA(I1,J5) + G4(I,J) * FAC0

               EA(I2,J1) = EA(I2,J1) + G5(I,J) * FAC0
               EA(I2,J2) = EA(I2,J2) + G6(I,J) * FAC0
               EA(I2,J3) = EA(I2,J3) + G5(I,J) * FAC0
               EA(I2,J4) = EA(I2,J4) + G7(I,J) * FAC0
               EA(I2,J6) = EA(I2,J6) + G4(I,J) * FAC0

               EA(I3,J1) = EA(I3,J1) + G8(I,J) * FAC0
               EA(I3,J2) = EA(I3,J2) + G8(I,J) * FAC0
               EA(I3,J3) = EA(I3,J3) + G9(I,J) * FAC0
               EA(I3,J5) = EA(I3,J5) + G7(I,J) * FAC0
               EA(I3,J6) = EA(I3,J6) + G3(I,J) * FAC0

            END DO
         END DO
!---------- ---------- ---------- ---------- ----------


      END DO
   END DO
END DO
!---------- ---------- ---------- ---------- ----------


RETURN
END


!************************************************
! REGULAR DOMAIN surface tractions
!************************************************
! REVISION : Th, 24 May 2012

SUBROUTINE ELEMENT_FS_3D ( XYZ, INOD, NGP, MTEL, PMAT, ID_BC, IEL, EFS )
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
DIMENSION LFACE(6), FVAL(6), IPERM(3), XI(NDIM), VR(NDIM), VS(NDIM), VN(NDIM)
DIMENSION XINT(9,9), WINT(9,9)

DIMENSION H1(NNODE)
DIMENSION G1(NNODE), G2(NNODE), G3(NNODE)

DATA LFACE / 2, 2, 3, 3, 1, 1/
DATA FVAL  /+1.0d0, -1.0d0, +1.0d0, -1.0d0, +1.0d0, -1.0d0/
DATA IPERM / 2, 3, 1/

dimension x_check(ndim)
!---------- ---------- ---------- ---------- ----------   


CALL GAUSS ( XINT, WINT )


! load element properties
!---------- ---------- ---------- ---------- ----------
NINT = NGP (IEL)


! define factors and initialize: only Tz component
!---------- ---------- ---------- ---------- ----------   
FAC1 = 0.0d0                                     ! Tx = 0 
FAC2 = 0.0d0                                     ! Ty = 0
FAC3 = 1.0d0                                     ! Tz = 1
!---------- ---------- ---------- ---------- ----------   
EFS  = 0.0d0


NFACE = ID_BC(IEL, 1)
IF (NFACE == 0) RETURN

! modify this later
if ( nface == 6 ) then
   nface = 5
!   write(*,*) 'temporary adjustment for nface in fs_3d'
end if


DO I = 1, NNODE
   K = INOD(I,IEL)
   DO J = 1, NDIM
      XT(J,I) = XYZ(K,J)
   END DO
END DO

!write(*,'(8f8.3)') xt(1,1:8)
!write(*,'(8f8.3)') xt(2,1:8)
!write(*,'(8f8.3)') xt(3,1:8)


L     = LFACE (NFACE)                            ! principal direction of integration: R
M     = IPERM (L)                                ! principal direction of integration: S
N     = IPERM (M)                                ! perpendicular to the principal directions, and constant
XI(N) = FVAL  (NFACE)                            ! constant Gauss value at boundary (surface)


!write(*,*) l,m,n
!write(*,*) xi(n)
!stop

! Gauss-Legendre quadrature on boundary (surface)
!---------- ---------- ---------- ---------- ----------   
DO LJ = 1, NINT
   XI(M) = XINT(LJ,NINT)
   WJ    = WINT(LJ,NINT)

   DO LI = 1, NINT
      XI(L) = XINT(LI,NINT)
      WI    = WINT(LI,NINT)


!xi(1) = 0.0d0
!xi(2) = 0.0d0
!xi(3) = 1.0d0
!wi = 1.0d0
!wj = 1.0d0



      WSTAR = WI * WJ


! Jacobian and shape functions
! VN is the outward normal vector to the boundary         
!---------- ---------- ---------- ---------- ----------           
      SELECT CASE (NNODE)

         CASE (20)                               ! 20-noded elements (serendipity)
            CALL SHAPE_3_20 ( XI(1), XI(2), XI(3), FN, DFXI )

         CASE (27)                               ! 27-noded elements (Lagrange)
            CALL SHAPE_3_27 ( XI(1), XI(2), XI(3), FN, DFXI )

         CASE DEFAULT
            WRITE(*,'("ERROR IN FS - NNODE")')
            STOP

      END SELECT

      DJ = MATMUL(XT,DFXI)

      DO I = 1, NDIM
         VR(I) = DJ(I,L)
         VS(I) = DJ(I,M)
      END DO

      CALL CROSS ( VR, VS, VN )
      CALL NORM  ( VN, NDIM, DETJ_SURF )

! correct the normal direction
      FACLAN = FVAL(NFACE)
      VN     = VN * FACLAN
!---------- ---------- ---------- ---------- ---------- 
      FAC0 = WSTAR * DETJ_SURF


! smooth loading function on the surface
!---------- ---------- ---------- ---------- ---------- ========== ========== ========== ========== ========== M-W-M-W-M-W-M-W-M-W-M-W-M-W-M-W-M-W-M-W-
XG = Dot_Product ( FN, XT( 1 , : ) )
YG = Dot_Product ( FN, XT( 2 , : ) )
ZG = Dot_Product ( FN, XT( 3 , : ) )
f  = ( 1.0d0 - ( XG/2.0d0 )**2 ) * ( 1.0d0 - ( YG/2.0d0 )**2 )
!write(*,*) XG, YG, ZG
!FAC3 = 1.0d0 * f
!dist = DSQRT ( XG**2 + YG**2 + ZG**2 )
!if ( dist >= 0.0001 ) then
!   cycle
!else
!   write(*,*) XG, YG, ZG
!   fac0 = 1.0d0
!   fac3 = 0.250d0
!end if

!---------- ---------- ---------- ---------- ---------- ========== ========== ========== ========== ========== M-W-M-W-M-W-M-W-M-W-M-W-M-W-M-W-M-W-M-W-


! two things are infinite: the universe and ... -----------------------------------------------------------------------------------------------------
! additional control: fix nface and remove this redundant mess later
! make sure load is being applied on the surface only
x_check = matmul(xt, fn)
if ( dabs(x_check(3)) > 1.e-3  ) then
   write(*,'(3f8.3)') x_check
   write(*,*) 'loaded face nont on the surface. see element matrices 3d'
   stop
end if
! ---------------------------------------------------------------------------------------------------------------------------------------------------


idebug = 9
if (idebug == 96) then
    write(*,'(2i5,3f10.3)') iel, nface, vn(1)/DETJ_SURF, vn(2)/DETJ_SURF, vn(3)/DETJ_SURF
end if



! form surface force vector       
!---------- ---------- ---------- ---------- ----------
      DO I = 1 , NNODE
         H1(I) = FN(I)
      END DO

      G1 = FAC1 * H1
      G2 = FAC2 * H1
      G3 = FAC3 * H1

      DO I = 1 , NNODE
         I1 = I
         I2 = I1 + NNODE
         I3 = I2 + NNODE

         EFS(I1) = EFS(I1) + G1(I) * FAC0
         EFS(I2) = EFS(I2) + G2(I) * FAC0
         EFS(I3) = EFS(I3) + G3(I) * FAC0
      END DO
!---------- ---------- ---------- ---------- ----------          


   END DO
END DO
!---------- ---------- ---------- ---------- ----------
!write(*,*) efs

RETURN
END


!************************************************
! PML ELEMENT EN MATRIX - Sezgin's way (implicit)
!************************************************
! REVISION : W 18 July 2012

SUBROUTINE ELEMENT_EN_3D_Sezgin ( XYZ, INOD, NGP, MTEL, PMAT, PML_PARAM, IEL, EN, IFLAG )
USE PARAMETERS
IMPLICIT DOUBLE PRECISION (A-H,O-Z)


! in 
!---------- ---------- ---------- ---------- ----------
DIMENSION XYZ(NJ,NDIM), INOD(NNODE,NEL), NGP(NEL), MTEL(NEL), PMAT(NMAT,NPM)
DIMENSION PML_PARAM(NDIM * 2 , 4)
!---------- ---------- ---------- ---------- ----------


! out
!---------- ---------- ---------- ---------- ----------
DIMENSION EN(NNODE * 6, NNODE * 6)
!---------- ---------- ---------- ---------- ----------


! local
!---------- ---------- ---------- ---------- ----------
DIMENSION XT(NDIM,NNODE), FN(NNODE), DJ(NDIM,NDIM), DFXI(NNODE,NDIM)
DIMENSION XINT(9,9), WINT(9,9)

DIMENSION H1(NNODE, NNODE)
DIMENSION G1(NNODE, NNODE), G2(NNODE, NNODE), G3(NNODE, NNODE)

DIMENSION XG(NDIM), PML_ALPHA_BETA(19)
!---------- ---------- ---------- ---------- ----------


SELECT CASE (IFLAG)

   CASE (301, 302, 303, 304)                     ! ENa, ENb, ENc, ENd
      CALL GAUSS_Lobatto ( XINT, WINT )
      NINT = NINT_Lobatto

   CASE DEFAULT
      WRITE(*,'("FLAG ERROR IN EN_3D")')
      STOP

END SELECT


! load material properties
!---------- ---------- ---------- ---------- ----------
MTYPE   = MTEL(IEL)
E       = PMAT(MTYPE,1)
ANU     = PMAT(MTYPE,2)
Zlambda = E * ANU / ( (1.0d0 + ANU) * (1.0d0 - 2.0d0 * ANU) )
Zmu     = E / ( 2.0d0*(1.0d0 + ANU) )


! define factors and initialize.
!---------- ---------- ---------- ---------- ----------
FAC1 = ( zLAMBDA + zMU ) / ( zMU * ( 3.0d0 * zLAMBDA + 2.0d0 * zMU ) ) 
FAC2 =  -zLAMBDA / ( 2.0d0 * zMU * ( 3.0d0 * zLAMBDA + 2.0d0 * zMU ) )
FAC3 = 1.0d0 / zMU
!---------- ---------- ---------- ---------- ----------
EN = 0.0d0


DO I = 1, NNODE
   K = INOD(I,IEL)
   DO J = 1, NDIM
      XT(J,I) = XYZ(K,J)
   END DO
END DO


! Gauss integration        
!---------- ---------- ---------- ---------- ----------   
DO LZ = 1, NINT
   X3 = XINT(LZ,NINT)
   WZ = WINT(LZ,NINT)

   DO LY = 1, NINT
      X2 = XINT(LY,NINT)
      WY = WINT(LY,NINT)

      DO LX = 1, NINT
         X1 = XINT(LX,NINT)
         WX = WINT(LX,NINT)

         WSTAR = WX * WY * WZ


! Jacobian and shape functions    
!---------- ---------- ---------- ---------- ----------           
         SELECT CASE (NNODE)

            CASE (20)                            ! 20-noded elements (serendipity)
               CALL SHAPE_3_20 (X1, X2, X3, FN, DFXI)

            CASE (27)                            ! 27-noded elements (Lagrange)
               CALL SHAPE_3_27 (X1, X2, X3, FN, DFXI)

            CASE DEFAULT
               WRITE(*,'("FLAG ERROR IN EN_3D - NNODE")')
               STOP

         END SELECT


         DJ = MATMUL(XT,DFXI)

         DETJ = DJ(1,1) * DJ(2,2) * DJ(3,3) + DJ(2,1) * DJ(3,2) * DJ(1,3) + DJ(3,1) * DJ(1,2) * DJ(2,3) - &
                DJ(1,3) * DJ(2,2) * DJ(3,1) - DJ(2,3) * DJ(3,2) * DJ(1,1) - DJ(3,3) * DJ(1,2) * DJ(2,1)

         IF (DETJ.LE.0) THEN
            WRITE(*,*)'|J|<0 NEGATIVE DET. JACOBIAN'
            READ(*,*)
            STOP
         END IF
!---------- ---------- ---------- ---------- ---------- 
         FAC0 = WSTAR * DETJ


! define factors                  
!---------- ---------- ---------- ---------- ---------- 
         XG = MATMUL(XT, FN)
         CALL PML_FINDER_3D ( XG, PML_PARAM, PML_ALPHA_BETA )

         SELECT CASE (IFLAG)

            CASE (301)                           ! ENa
               FAC4 = PML_ALPHA_BETA(7)

            CASE (302)                           ! ENb
               FAC4 = PML_ALPHA_BETA(8)

            CASE (303)                           ! ENc
               FAC4 = PML_ALPHA_BETA(9)

            CASE (304)                           ! ENd
               FAC4 = PML_ALPHA_BETA(10)

            CASE DEFAULT
               WRITE(*,'("FLAG ERROR IN EN_3D")')
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

         G1 = FAC1 * FAC4 * H1
         G2 = FAC2 * FAC4 * H1
         G3 = FAC3 * FAC4 * H1

         DO I = 1 , NNODE
            I1 = I
            I2 = I1 + NNODE
            I3 = I2 + NNODE
            I4 = I3 + NNODE
            I5 = I4 + NNODE
            I6 = I5 + NNODE

            DO J = 1 , NNODE
               J1 = J
               J2 = J1 + NNODE
               J3 = J2 + NNODE
               J4 = J3 + NNODE
               J5 = J4 + NNODE
               J6 = J5 + NNODE

               EN(I1,J1) = EN(I1,J1) + G1(I,J) * FAC0
               EN(I1,J2) = EN(I1,J2) + G2(I,J) * FAC0
               EN(I1,J3) = EN(I1,J3) + G2(I,J) * FAC0

               EN(I2,J1) = EN(I2,J1) + G2(I,J) * FAC0
               EN(I2,J2) = EN(I2,J2) + G1(I,J) * FAC0
               EN(I2,J3) = EN(I2,J3) + G2(I,J) * FAC0

               EN(I3,J1) = EN(I3,J1) + G2(I,J) * FAC0
               EN(I3,J2) = EN(I3,J2) + G2(I,J) * FAC0
               EN(I3,J3) = EN(I3,J3) + G1(I,J) * FAC0

               EN(I4,J4) = EN(I4,J4) + G3(I,J) * FAC0
               EN(I5,J5) = EN(I5,J5) + G3(I,J) * FAC0
               EN(I6,J6) = EN(I6,J6) + G3(I,J) * FAC0

            END DO
         END DO
!---------- ---------- ---------- ---------- ----------


      END DO
   END DO
END DO
!---------- ---------- ---------- ---------- ----------


RETURN
END


!************************************************
! PML ELEMENT EA MATRIX - Sezgin's way

! IFLAG == 101 evanescent  wave => EAe_upper 
! IFLAG == 102 propagaring wave => EAp_upper
! IFLAG == 103                  => EAw_upper
! IFLAG == 104 evanescent  wave => EAe_lower
! IFLAG == 105 propagaring wave => EAp_lower
! IFLAG == 106                  => EAw_lower
!************************************************
! REVISION : Sun, 22 July 2012

SUBROUTINE ELEMENT_EA_3D_Sezgin ( XYZ, INOD, NGP, MTEL, PMAT, PML_PARAM, IEL, EA, IFLAG )
USE PARAMETERS
IMPLICIT DOUBLE PRECISION (A-H,O-Z)


! in 
!---------- ---------- ---------- ---------- ----------
DIMENSION XYZ(NJ,NDIM), INOD(NNODE,NEL), NGP(NEL), MTEL(NEL), PMAT(NMAT,NPM)
DIMENSION PML_PARAM(NDIM * 2 , 4)
!---------- ---------- ---------- ---------- ----------


! out 
!---------- ---------- ---------- ---------- ----------
DIMENSION EA(NNODE * NDIM, NNODE * 6)
!---------- ---------- ---------- ---------- ----------


! local
!---------- ---------- ---------- ---------- ----------
DIMENSION XT(NDIM,NNODE), DJ(NDIM,NDIM), DJI(NDIM,NDIM), FN(NNODE), DFX(NNODE,NDIM), DFXI(NNODE,NDIM)
DIMENSION XINT(9,9), WINT(9,9)

DIMENSION H1(NNODE, NNODE), H2(NNODE, NNODE), H3(NNODE, NNODE), H4(NNODE, NNODE)
DIMENSION G1(NNODE, NNODE), G3(NNODE, NNODE), G4(NNODE, NNODE), G6(NNODE, NNODE), G7(NNODE, NNODE), G9(NNODE, NNODE)

DIMENSION XG(NDIM), PML_ALPHA_BETA(19)
!---------- ---------- ---------- ---------- ----------


CALL GAUSS ( XINT, WINT )


! load element properties
!---------- ---------- ---------- ---------- ----------
NINT = NGP (IEL)


! define factors and initialize.
!---------- ---------- ---------- ---------- ---------- 
SELECT CASE (IFLAG)

   CASE (101, 102, 103)                          ! upper block
      FAC10 = 1.0d0

   CASE (104, 105, 106)                          ! lower block
      FAC10 = 0.0d0

   CASE DEFAULT
      WRITE(*,'("FLAG ERROR IN EA_3D")')
      STOP

END SELECT
!---------- ---------- ---------- ---------- ----------
EA = 0.0d0


DO I = 1, NNODE
   K = INOD(I,IEL)
   DO J = 1, NDIM
      XT(J,I) = XYZ(K,J)
   END DO
END DO
        

! Gauss integration        
!---------- ---------- ---------- ---------- ----------
DO LZ = 1, NINT
   X3 = XINT(LZ,NINT)
   WZ = WINT(LZ,NINT)

   DO LY = 1, NINT
      X2 = XINT(LY,NINT)
      WY = WINT(LY,NINT)

      DO LX = 1, NINT
         X1 = XINT(LX,NINT)
         WX = WINT(LX,NINT)

         WSTAR = WX * WY * WZ


! Jacobian and shape functions    
!---------- ---------- ---------- ---------- ----------         
         SELECT CASE (NNODE)

            CASE (20)                            ! 20-noded elements (serendipity)
               CALL SHAPE_3_20 (X1, X2, X3, FN, DFXI)

            CASE (27)                            ! 27-noded elements (Lagrange)
               CALL SHAPE_3_27 (X1, X2, X3, FN, DFXI)

            CASE DEFAULT
               WRITE(*,'("FLAG ERROR IN EA_3D - NNODE")')
               STOP

         END SELECT


         DJ = MATMUL(XT,DFXI)

         DETJ = DJ(1,1) * DJ(2,2) * DJ(3,3) + DJ(2,1) * DJ(3,2) * DJ(1,3) + DJ(3,1) * DJ(1,2) * DJ(2,3) - &
                DJ(1,3) * DJ(2,2) * DJ(3,1) - DJ(2,3) * DJ(3,2) * DJ(1,1) - DJ(3,3) * DJ(1,2) * DJ(2,1)

         IF (DETJ.LE.0) THEN
            WRITE(*,*)'|J|<0 NEGATIVE DET. JACOBIAN'
            READ(*,*)
            STOP
         END IF

         DJI(1,1) = + DJ(3,3) * DJ(2,2) - DJ(3,2) * DJ(2,3)
         DJI(1,2) = - DJ(3,3) * DJ(1,2) + DJ(3,2) * DJ(1,3)
         DJI(1,3) = + DJ(2,3) * DJ(1,2) - DJ(2,2) * DJ(1,3)

         DJI(2,1) = - DJ(3,3) * DJ(2,1) + DJ(3,1) * DJ(2,3)
         DJI(2,2) = + DJ(3,3) * DJ(1,1) - DJ(3,1) * DJ(1,3)
         DJI(2,3) = - DJ(2,3) * DJ(1,1) + DJ(2,1) * DJ(1,3)

         DJI(3,1) = + DJ(3,2) * DJ(2,1) - DJ(3,1) * DJ(2,2)
         DJI(3,2) = - DJ(3,2) * DJ(1,1) + DJ(3,1) * DJ(1,2)
         DJI(3,3) = + DJ(2,2) * DJ(1,1) - DJ(2,1) * DJ(1,2)

         DJI = DJI / DETJ


         DFX  = MATMUL(DFXI,DJI)
!---------- ---------- ---------- ---------- ---------- 
         FAC0 = WSTAR * DETJ


! define factors                  
!---------- ---------- ---------- ---------- ---------- 
         XG = MATMUL(XT, FN)
         CALL PML_FINDER_3D ( XG , PML_PARAM , PML_ALPHA_BETA )

         SELECT CASE (IFLAG)

            CASE (101, 104)                                                                           ! e-wave
               FAC1 = PML_ALPHA_BETA(1) * PML_ALPHA_BETA(2)                                           !    alpha_x * alpha_y
               FAC2 = PML_ALPHA_BETA(1) * PML_ALPHA_BETA(3)                                           !    alpha_x * alpha_z
               FAC3 = PML_ALPHA_BETA(2) * PML_ALPHA_BETA(3)                                           !    alpha_y * alpha_z
               FAC4 = PML_ALPHA_BETA(11)                                                              ! d (alpha_x * alpha_y) 
               FAC5 = PML_ALPHA_BETA(12)                                                              ! d (alpha_x * alpha_z)
               FAC6 = PML_ALPHA_BETA(13)                                                              ! d (alpha_y * alpha_z)

            CASE (102, 105)                                                                           ! p-wave  
               FAC1 = PML_ALPHA_BETA(1) * PML_ALPHA_BETA(5) + PML_ALPHA_BETA(2) * PML_ALPHA_BETA(4)   !    alpha_x * beta_y + alpha_y * beta_x
               FAC2 = PML_ALPHA_BETA(1) * PML_ALPHA_BETA(6) + PML_ALPHA_BETA(3) * PML_ALPHA_BETA(4)   !    alpha_x * beta_z + alpha_z * beta_x        
               FAC3 = PML_ALPHA_BETA(2) * PML_ALPHA_BETA(6) + PML_ALPHA_BETA(3) * PML_ALPHA_BETA(5)   !    alpha_y * beta_z + alpha_z * beta_y         
               FAC4 = PML_ALPHA_BETA(14)                                                              ! d (alpha_x * beta_y + alpha_y * beta_x) 
               FAC5 = PML_ALPHA_BETA(15)                                                              ! d (alpha_x * beta_z + alpha_z * beta_x)
               FAC6 = PML_ALPHA_BETA(16)                                                              ! d (alpha_y * beta_z + alpha_z * beta_y)

            CASE (103, 106)                                                                           ! w-wave    
               FAC1 = PML_ALPHA_BETA(4) * PML_ALPHA_BETA(5)                                           !    beta_x * beta_y
               FAC2 = PML_ALPHA_BETA(4) * PML_ALPHA_BETA(6)                                           !    beta_x * beta_z
               FAC3 = PML_ALPHA_BETA(5) * PML_ALPHA_BETA(6)                                           !    beta_y * beta_z
               FAC4 = PML_ALPHA_BETA(17)                                                              ! d (beta_x * beta_y)
               FAC5 = PML_ALPHA_BETA(18)                                                              ! d (beta_x * beta_z)
               FAC6 = PML_ALPHA_BETA(19)                                                              ! d (beta_y * beta_z)

            CASE DEFAULT
               WRITE(*,'("FLAG ERROR IN EA_3D")')
               STOP

         END SELECT
!---------- ---------- ---------- ---------- ----------


! form element matrix
!---------- ---------- ---------- ---------- ----------
         DO I = 1 , NNODE
            DO J = 1 , NNODE
               H1(I,J) = DFX(I,1) * FN(J)
               H2(I,J) = DFX(I,2) * FN(J)
               H3(I,J) = DFX(I,3) * FN(J)
               H4(I,J) =  FN(I)   * FN(J)
            END DO
         END DO

         G1 =  FAC3 * H1 + FAC10 * FAC6 * H4
         G3 =  FAC2 * H2 + FAC10 * FAC5 * H4
         G4 =  FAC1 * H3 + FAC10 * FAC4 * H4
         G6 =  FAC2 * H2 + FAC10 * FAC5 * H4
         G7 =  FAC3 * H1 + FAC10 * FAC6 * H4
         G9 =  FAC1 * H3 + FAC10 * FAC4 * H4

         DO I = 1 , NNODE
            I1 = I
            I2 = I1 + NNODE
            I3 = I2 + NNODE

            DO J = 1 , NNODE
               J1 = J
               J2 = J1 + NNODE
               J3 = J2 + NNODE
               J4 = J3 + NNODE
               J5 = J4 + NNODE
               J6 = J5 + NNODE

               EA(I1,J1) = EA(I1,J1) + G1(I,J) * FAC0
               EA(I1,J4) = EA(I1,J4) + G3(I,J) * FAC0
               EA(I1,J5) = EA(I1,J5) + G4(I,J) * FAC0

               EA(I2,J2) = EA(I2,J2) + G6(I,J) * FAC0
               EA(I2,J4) = EA(I2,J4) + G7(I,J) * FAC0
               EA(I2,J6) = EA(I2,J6) + G4(I,J) * FAC0

               EA(I3,J3) = EA(I3,J3) + G9(I,J) * FAC0
               EA(I3,J5) = EA(I3,J5) + G7(I,J) * FAC0
               EA(I3,J6) = EA(I3,J6) + G3(I,J) * FAC0

            END DO
         END DO
!---------- ---------- ---------- ---------- ----------


      END DO
   END DO
END DO
!---------- ---------- ---------- ---------- ----------


RETURN
END


!************************************************
! Regularization matrix R (weak-form of the Laplacian operator)
! inverse of the gradient Mass matrix (due to Spectral-Elements, it is diagonal)
!************************************************
! REVISION : Th, 16 January 2014: You are remembered for the rules you break.
! It is not enough to have a good mind; the main thing is to use it well.
! REVISION : T, 9 September 2014: Weighted-regularization scheme


SUBROUTINE ELEMENT_R_Mg_3D ( XYZ, INOD, NGP, MTEL, IEL, ER, EM_g_diagonal )
USE PARAMETERS
IMPLICIT DOUBLE PRECISION (A-H,O-Z)


! in
!---------- ---------- ---------- ---------- ----------
DIMENSION XYZ ( NJ, NDIM ), INOD ( NNODE , NEL ), NGP ( NEL ), MTEL ( NEL )
!---------- ---------- ---------- ---------- ----------


! out
!---------- ---------- ---------- ---------- ----------
DIMENSION ER   ( NNODE , NNODE )
Dimension EM_g ( NNODE , NNODE ), EM_g_diagonal ( NNODE )
!---------- ---------- ---------- ---------- ----------


! local
!---------- ---------- ---------- ---------- ----------
DIMENSION XT ( NDIM, NNODE ), DJ ( NDIM, NDIM ), DJI ( NDIM, NDIM ), FN ( NNODE ), DFX ( NNODE, NDIM ), DFXI ( NNODE, NDIM )
DIMENSION XINT(9,9), WINT(9,9)

DIMENSION H1 ( NNODE, NNODE ), H2 ( NNODE, NNODE ), H3 ( NNODE, NNODE )
DIMENSION H0 ( NNODE, NNODE )
!---------- ---------- ---------- ---------- ----------


CALL GAUSS ( XINT, WINT )


! load material properties
!---------- ---------- ---------- ---------- ----------
NINT  = NGP(IEL)


! define factors and initialize
!---------- ---------- ---------- ---------- ---------- 
ER    = 0.0d0
EM_g  = 0.0d0
!EiM_g = 0.0d0


DO I = 1, NNODE
   K = INOD(I,IEL)
   DO J = 1, NDIM
      XT(J,I) = XYZ(K,J)
   END DO
END DO


! Gauss integration        
!---------- ---------- ---------- ---------- ----------
DO LZ = 1, NINT
   X3 = XINT(LZ,NINT)
   WZ = WINT(LZ,NINT)

   DO LY = 1, NINT
      X2 = XINT(LY,NINT)
      WY = WINT(LY,NINT)

      DO LX = 1, NINT
         X1 = XINT(LX,NINT)
         WX = WINT(LX,NINT)

         WSTAR = WX * WY * WZ


! Jacobian and shape functions    
!---------- ---------- ---------- ---------- ----------         
         SELECT CASE (NNODE)

            CASE (20)                            ! 20-noded elements (serendipity)
               CALL SHAPE_3_20 (X1, X2, X3, FN, DFXI)

            CASE (27)                            ! 27-noded elements (Lagrange)
               CALL SHAPE_3_27 (X1, X2, X3, FN, DFXI)

            CASE DEFAULT
               WRITE(*,'("FLAG ERROR IN EK_3D - NNODE")')
               STOP

         END SELECT


         DJ = MATMUL(XT,DFXI)

         DETJ = DJ(1,1) * DJ(2,2) * DJ(3,3) + DJ(2,1) * DJ(3,2) * DJ(1,3) + DJ(3,1) * DJ(1,2) * DJ(2,3) - &
                DJ(1,3) * DJ(2,2) * DJ(3,1) - DJ(2,3) * DJ(3,2) * DJ(1,1) - DJ(3,3) * DJ(1,2) * DJ(2,1)

         IF (DETJ.LE.0) THEN
            WRITE(*,*)'|J|<0 NEGATIVE DET. JACOBIAN'
            READ(*,*)
            STOP
         END IF

         DJI(1,1) = + DJ(3,3) * DJ(2,2) - DJ(3,2) * DJ(2,3)
         DJI(1,2) = - DJ(3,3) * DJ(1,2) + DJ(3,2) * DJ(1,3)
         DJI(1,3) = + DJ(2,3) * DJ(1,2) - DJ(2,2) * DJ(1,3)
 
         DJI(2,1) = - DJ(3,3) * DJ(2,1) + DJ(3,1) * DJ(2,3)
         DJI(2,2) = + DJ(3,3) * DJ(1,1) - DJ(3,1) * DJ(1,3)
         DJI(2,3) = - DJ(2,3) * DJ(1,1) + DJ(2,1) * DJ(1,3)

         DJI(3,1) = + DJ(3,2) * DJ(2,1) - DJ(3,1) * DJ(2,2)
         DJI(3,2) = - DJ(3,2) * DJ(1,1) + DJ(3,1) * DJ(1,2)
         DJI(3,3) = + DJ(2,2) * DJ(1,1) - DJ(2,1) * DJ(1,2)

         DJI = DJI / DETJ
              

         DFX = MATMUL(DFXI,DJI)
!---------- ---------- ---------- ---------- ---------- 
         FAC0 = WSTAR * DETJ


! Weighted-regularization scheme
!---------- ---------- ---------- ---------- ---------- $$$$$$$
         Z_Gauss = Dot_Product ( FN , XT ( 3 , : ) )
         f1      = 1.0d0 - Dabs ( Z_Gauss / w_Length )
         f2      =         Dabs ( Z_Gauss / w_Length )
         weight  = f1 * w_top + f2 * w_bot
         Select Case ( Weighted_Regularization )
            Case ('Yes')
               Fac1 = weight
            Case ('No')
               Fac1 = 1.0d0
         End Select
!---------- ---------- ---------- ---------- ---------- $$$$$$$


! form the matrix
!---------- ---------- ---------- ---------- ----------
         DO I = 1 , NNODE
            DO J = 1 , NNODE
               H1(I,J) = DFX(I,1) * DFX(J,1)
               H2(I,J) = DFX(I,2) * DFX(J,2)
               H3(I,J) = DFX(I,3) * DFX(J,3)
               H0(I,J) = FN(I) * FN(J)
            END DO
         END DO

         ! regularization matrix
         ER = ER + ( H1 + H2 + H3 ) * FAC0 * Fac1

         ! gradient mass matrix (diagonal)
         EM_g = EM_g + H0 * FAC0
!---------- ---------- ---------- ---------- ----------


      END DO
   END DO
END DO
!---------- ---------- ---------- ---------- ----------

         ! extract diagonal
         Do I = 1 , NNode
           EM_g_diagonal (I) = EM_g (I,I)
         End Do


RETURN
END


!************************************************
! gradient update (Mu & Lambda)
!************************************************
! REVISION : W, 15 January 2014 : A life spent making mistakes is not only more honorable, but more useful than a life spent doing nothing.

SUBROUTINE ELEMENT_g_3D ( XYZ, INOD, NGP, IEL, E_U, E_P, E_g_Mu, E_g_Lambda )
USE PARAMETERS
IMPLICIT DOUBLE PRECISION (A-H,O-Z)


! in 
!---------- ---------- ---------- ---------- ----------
DIMENSION XYZ ( NJ, NDIM ), INOD ( NNODE, NEL ), NGP ( NEL )
Dimension E_U ( NNODE * NDIM )
Dimension E_P ( NNODE * NDIM )
!---------- ---------- ---------- ---------- ----------


! out
!---------- ---------- ---------- ---------- ----------
Dimension E_g_Mu     ( NNODE )
Dimension E_g_Lambda ( NNODE )
!---------- ---------- ---------- ---------- ----------


! local
!---------- ---------- ---------- ---------- ----------
!Dimension E_core_Mu     ( NNODE * NDIM , NNODE * NDIM )
!Dimension E_core_Lambda ( NNODE * NDIM , NNODE * NDIM )
!---------- ---------- ---------- ---------- ----------


! local
!---------- ---------- ---------- ---------- ----------
DIMENSION XT(NDIM,NNODE), DJ(NDIM,NDIM), DJI(NDIM,NDIM), FN(NNODE), DFX(NNODE,NDIM), DFXI(NNODE,NDIM)
DIMENSION XINT(9,9), WINT(9,9)

!DIMENSION H1(NNODE, NNODE), H2(NNODE, NNODE), H3(NNODE, NNODE), H4(NNODE, NNODE), H5(NNODE, NNODE), H6(NNODE, NNODE), H7(NNODE, NNODE), H8(NNODE, NNODE), H9(NNODE, NNODE)
!DIMENSION G1(NNODE, NNODE), G2(NNODE, NNODE), G3(NNODE, NNODE), G4(NNODE, NNODE), G5(NNODE, NNODE), G6(NNODE, NNODE), G7(NNODE, NNODE), G8(NNODE, NNODE), G9(NNODE, NNODE)
!---------- ---------- ---------- ---------- ----------


CALL GAUSS ( XINT, WINT )


! load properties
!---------- ---------- ---------- ---------- ----------
NINT = NGP(IEL)


! define factors and initialize
!---------- ---------- ---------- ---------- ---------- 
E_core_Mu     = 0.0d0
E_core_Lambda = 0.0d0
E_g_Mu        = 0.0d0
E_g_Lambda    = 0.0d0


DO I = 1, NNODE
   K = INOD(I,IEL)
   DO J = 1, NDIM
      XT(J,I) = XYZ(K,J)
   END DO
END DO


! Gauss integration        
!---------- ---------- ---------- ---------- ----------
DO LZ = 1, NINT
   X3 = XINT(LZ,NINT)
   WZ = WINT(LZ,NINT)

   DO LY = 1, NINT
      X2 = XINT(LY,NINT)
      WY = WINT(LY,NINT)

      DO LX = 1, NINT
         X1 = XINT(LX,NINT)
         WX = WINT(LX,NINT)

         WSTAR = WX * WY * WZ


! Jacobian and shape functions    
!---------- ---------- ---------- ---------- ----------         
         SELECT CASE (NNODE)

            CASE (20)                            ! 20-noded elements (serendipity)
               CALL SHAPE_3_20 (X1, X2, X3, FN, DFXI)

            CASE (27)                            ! 27-noded elements (Lagrange)
               CALL SHAPE_3_27 (X1, X2, X3, FN, DFXI)

            CASE DEFAULT
               WRITE(*,'("FLAG ERROR IN EK_3D - NNODE")')
               STOP

         END SELECT


         DJ = MATMUL(XT,DFXI)

         DETJ = DJ(1,1) * DJ(2,2) * DJ(3,3) + DJ(2,1) * DJ(3,2) * DJ(1,3) + DJ(3,1) * DJ(1,2) * DJ(2,3) - &
                DJ(1,3) * DJ(2,2) * DJ(3,1) - DJ(2,3) * DJ(3,2) * DJ(1,1) - DJ(3,3) * DJ(1,2) * DJ(2,1)

         IF (DETJ.LE.0) THEN
            WRITE(*,*)'|J|<0 NEGATIVE DET. JACOBIAN'
            READ(*,*)
            STOP
         END IF

         DJI(1,1) = + DJ(3,3) * DJ(2,2) - DJ(3,2) * DJ(2,3)
         DJI(1,2) = - DJ(3,3) * DJ(1,2) + DJ(3,2) * DJ(1,3)
         DJI(1,3) = + DJ(2,3) * DJ(1,2) - DJ(2,2) * DJ(1,3)
 
         DJI(2,1) = - DJ(3,3) * DJ(2,1) + DJ(3,1) * DJ(2,3)
         DJI(2,2) = + DJ(3,3) * DJ(1,1) - DJ(3,1) * DJ(1,3)
         DJI(2,3) = - DJ(2,3) * DJ(1,1) + DJ(2,1) * DJ(1,3)

         DJI(3,1) = + DJ(3,2) * DJ(2,1) - DJ(3,1) * DJ(2,2)
         DJI(3,2) = - DJ(3,2) * DJ(1,1) + DJ(3,1) * DJ(1,2)
         DJI(3,3) = + DJ(2,2) * DJ(1,1) - DJ(2,1) * DJ(1,2)

         DJI = DJI / DETJ


         DFX = MATMUL(DFXI,DJI)
!---------- ---------- ---------- ---------- ---------- 
         FAC0 = WSTAR * DETJ


! form the components
!---------- ---------- ---------- ---------- ----------
         u1_x = Dot_Product ( E_U (            1 :     NNODE ) , DFX(:,1) )
         u2_x = Dot_Product ( E_U (     NNODE +1 : 2 * NNODE ) , DFX(:,1) )
         u3_x = Dot_Product ( E_U ( 2 * NNODE +1 : 3 * NNODE ) , DFX(:,1) )
         u1_y = Dot_Product ( E_U (            1 :     NNODE ) , DFX(:,2) )
         u2_y = Dot_Product ( E_U (     NNODE +1 : 2 * NNODE ) , DFX(:,2) )
         u3_y = Dot_Product ( E_U ( 2 * NNODE +1 : 3 * NNODE ) , DFX(:,2) )
         u1_z = Dot_Product ( E_U (            1 :     NNODE ) , DFX(:,3) )
         u2_z = Dot_Product ( E_U (     NNODE +1 : 2 * NNODE ) , DFX(:,3) )
         u3_z = Dot_Product ( E_U ( 2 * NNODE +1 : 3 * NNODE ) , DFX(:,3) )

         p1_x = Dot_Product ( E_P (            1 :     NNODE ) , DFX(:,1) )
         p2_x = Dot_Product ( E_P (     NNODE +1 : 2 * NNODE ) , DFX(:,1) )
         p3_x = Dot_Product ( E_P ( 2 * NNODE +1 : 3 * NNODE ) , DFX(:,1) )
         p1_y = Dot_Product ( E_P (            1 :     NNODE ) , DFX(:,2) )
         p2_y = Dot_Product ( E_P (     NNODE +1 : 2 * NNODE ) , DFX(:,2) )
         p3_y = Dot_Product ( E_P ( 2 * NNODE +1 : 3 * NNODE ) , DFX(:,2) )
         p1_z = Dot_Product ( E_P (            1 :     NNODE ) , DFX(:,3) )
         p2_z = Dot_Product ( E_P (     NNODE +1 : 2 * NNODE ) , DFX(:,3) )
         p3_z = Dot_Product ( E_P ( 2 * NNODE +1 : 3 * NNODE ) , DFX(:,3) )
!---------- ---------- ---------- ---------- ----------

         scalar_Lambda = ( u1_x + u2_y + u3_z ) * ( p1_x + p2_y + p3_z )
         scalar_Mu     = 2.0d0 * ( u1_x * p1_x + u2_y * p2_y + u3_z * p3_z ) &
                       + ( p1_y + p2_x ) * ( u1_y + u2_x ) &
                       + ( p1_z + p3_x ) * ( u1_z + u3_x ) &
                       + ( p3_y + p2_z ) * ( u3_y + u2_z )

         E_g_Mu     = E_g_Mu     + FN * scalar_Mu     * FAC0
         E_g_Lambda = E_g_Lambda + FN * scalar_Lambda * FAC0

      END DO
   END DO
END DO
!---------- ---------- ---------- ---------- ----------


RETURN
END


!************************************************
! Matrix free implementation of regularization. This is necessary for TV since it is nonlinear.
!************************************************
! REVISION : W, 30 July 2014: Be the change you want to see in the world.
! REVISION : T, 9 September 2014: Weighted-regularization scheme

SUBROUTINE ELEMENT_TV_TN_Reg_3D ( XYZ, INOD, NGP, MTEL, IEL, E_Mu, E_Lambda, E_Reg_Mu, E_Reg_Lambda, e_cost_reg_Lambda, e_cost_reg_Mu )
USE PARAMETERS
IMPLICIT DOUBLE PRECISION (A-H,O-Z)


! in
!---------- ---------- ---------- ---------- ----------
DIMENSION XYZ ( NJ, NDIM ), INOD ( NNODE , NEL ), NGP ( NEL ), MTEL ( NEL )
Dimension E_Mu ( NNODE ), E_Lambda ( NNODE )
DIMENSION ER   ( NNODE , NNODE )
!---------- ---------- ---------- ---------- ----------

! out
!---------- ---------- ---------- ---------- ----------
Dimension E_Reg_Mu     ( NNODE )
Dimension E_Reg_Lambda ( NNODE )
!---------- ---------- ---------- ---------- ----------


! local
!---------- ---------- ---------- ---------- ----------
DIMENSION XT ( NDIM, NNODE ), DJ ( NDIM, NDIM ), DJI ( NDIM, NDIM ), FN ( NNODE ), DFX ( NNODE, NDIM ), DFXI ( NNODE, NDIM )
DIMENSION XINT(9,9), WINT(9,9)

DIMENSION H1 ( NNODE, NNODE ), H2 ( NNODE, NNODE ), H3 ( NNODE, NNODE )
DIMENSION vec_help ( NNODE )
!---------- ---------- ---------- ---------- ----------


CALL GAUSS ( XINT, WINT )


! load material properties
!---------- ---------- ---------- ---------- ----------
NINT  = NGP(IEL)


! define factors and initialize
!---------- ---------- ---------- ---------- ---------- 
epsilon_TV        = 0.01d0   !*************************

E_Reg_Mu          = 0.0d0
E_Reg_Lambda      = 0.0d0
e_cost_reg_Mu     = 0.0d0
e_cost_reg_Lambda = 0.0d0


DO I = 1, NNODE
   K = INOD(I,IEL)
   DO J = 1, NDIM
      XT(J,I) = XYZ(K,J)
   END DO
END DO


! Gauss integration        
!---------- ---------- ---------- ---------- ----------
DO LZ = 1, NINT
   X3 = XINT(LZ,NINT)
   WZ = WINT(LZ,NINT)

   DO LY = 1, NINT
      X2 = XINT(LY,NINT)
      WY = WINT(LY,NINT)

      DO LX = 1, NINT
         X1 = XINT(LX,NINT)
         WX = WINT(LX,NINT)

         WSTAR = WX * WY * WZ


! Jacobian and shape functions
!---------- ---------- ---------- ---------- ----------
         SELECT CASE (NNODE)

            CASE (20)                            ! 20-noded elements (serendipity)
               CALL SHAPE_3_20 (X1, X2, X3, FN, DFXI)

            CASE (27)                            ! 27-noded elements (Lagrange)
               CALL SHAPE_3_27 (X1, X2, X3, FN, DFXI)

            CASE DEFAULT
               WRITE(*,'("FLAG ERROR IN EK_3D - NNODE")')
               STOP

         END SELECT


         DJ = MATMUL(XT,DFXI)

         DETJ = DJ(1,1) * DJ(2,2) * DJ(3,3) + DJ(2,1) * DJ(3,2) * DJ(1,3) + DJ(3,1) * DJ(1,2) * DJ(2,3) - &
                DJ(1,3) * DJ(2,2) * DJ(3,1) - DJ(2,3) * DJ(3,2) * DJ(1,1) - DJ(3,3) * DJ(1,2) * DJ(2,1)

         IF (DETJ.LE.0) THEN
            WRITE(*,*)'|J|<0 NEGATIVE DET. JACOBIAN'
            READ(*,*)
            STOP
         END IF

         DJI(1,1) = + DJ(3,3) * DJ(2,2) - DJ(3,2) * DJ(2,3)
         DJI(1,2) = - DJ(3,3) * DJ(1,2) + DJ(3,2) * DJ(1,3)
         DJI(1,3) = + DJ(2,3) * DJ(1,2) - DJ(2,2) * DJ(1,3)
 
         DJI(2,1) = - DJ(3,3) * DJ(2,1) + DJ(3,1) * DJ(2,3)
         DJI(2,2) = + DJ(3,3) * DJ(1,1) - DJ(3,1) * DJ(1,3)
         DJI(2,3) = - DJ(2,3) * DJ(1,1) + DJ(2,1) * DJ(1,3)

         DJI(3,1) = + DJ(3,2) * DJ(2,1) - DJ(3,1) * DJ(2,2)
         DJI(3,2) = - DJ(3,2) * DJ(1,1) + DJ(3,1) * DJ(1,2)
         DJI(3,3) = + DJ(2,2) * DJ(1,1) - DJ(2,1) * DJ(1,2)

         DJI = DJI / DETJ
              

         DFX = MATMUL(DFXI,DJI)
!---------- ---------- ---------- ---------- ---------- 
         FAC0 = WSTAR * DETJ


! Weighted-regularization scheme
!---------- ---------- ---------- ---------- ---------- $$$$$$$
         Z_Gauss = Dot_Product ( FN , XT ( 3 , : ) )
         f1      = 1.0d0 - Dabs ( Z_Gauss / w_Length )
         f2      =         Dabs ( Z_Gauss / w_Length )
         weight  = f1 * w_top + f2 * w_bot
         Select Case ( Weighted_Regularization )
            Case ('Yes')
               Fac1 = weight
            Case ('No')
               Fac1 = 1.0d0
         End Select
!---------- ---------- ---------- ---------- ---------- $$$$$$$


! form the matrix
!---------- ---------- ---------- ---------- ----------
         DO I = 1 , NNODE
            DO J = 1 , NNODE
               H1(I,J) = DFX(I,1) * DFX(J,1)
               H2(I,J) = DFX(I,2) * DFX(J,2)
               H3(I,J) = DFX(I,3) * DFX(J,3)
            END DO
         END DO

         ! regularization matrix
         ER = H1 + H2 + H3

         ! Mu
         vec_help     = MATMUL( ER , E_Mu )
         denominator  = DSQRT ( Dot_Product ( E_Mu , vec_help ) + epsilon_TV )
!         denominator  = 1.0d0      !(TN)
         E_Reg_Mu     = E_Reg_Mu + ( vec_help / denominator ) * FAC0 * Fac1
         ! cost of regularization
         numerator    = Dot_Product ( E_Mu , vec_help )
         e_cost_reg_Mu= e_cost_reg_Mu + ( numerator / denominator ) * FAC0 * Fac1
         !------------------------------------------------------------------------
         ! Lambda
         vec_help     = MATMUL( ER , E_Lambda )
         denominator  = DSQRT ( Dot_Product ( E_Lambda , vec_help ) + epsilon_TV )
!         denominator  = 1.0d0      !(TN)
         E_Reg_Lambda = E_Reg_Lambda + ( vec_help / denominator ) * FAC0 * Fac1
         ! cost of regularization
         numerator    = Dot_Product ( E_Lambda , vec_help )
         e_cost_reg_Lambda = e_cost_reg_Lambda + ( numerator / denominator ) * FAC0 * Fac1
!---------- ---------- ---------- ---------- ----------


      END DO
   END DO
END DO
!---------- ---------- ---------- ---------- ----------


RETURN
END

!************************************************
! gradient update (Mu & Lambda+2Mu)
!************************************************
! REVISION : F, 1 August 2014 : TBA

SUBROUTINE ELEMENT_g_3D_v2 ( XYZ, INOD, NGP, IEL, E_U, E_P, E_g_Mu, E_g_Lambda_2Mu )
USE PARAMETERS
IMPLICIT DOUBLE PRECISION (A-H,O-Z)


! in 
!---------- ---------- ---------- ---------- ----------
DIMENSION XYZ ( NJ, NDIM ), INOD ( NNODE, NEL ), NGP ( NEL )
Dimension E_U ( NNODE * NDIM )
Dimension E_P ( NNODE * NDIM )
!---------- ---------- ---------- ---------- ----------


! out
!---------- ---------- ---------- ---------- ----------
Dimension E_g_Mu         ( NNODE )
Dimension E_g_Lambda_2Mu ( NNODE )
!---------- ---------- ---------- ---------- ----------


! local
!---------- ---------- ---------- ---------- ----------
!Dimension E_core_Mu     ( NNODE * NDIM , NNODE * NDIM )
!Dimension E_core_Lambda ( NNODE * NDIM , NNODE * NDIM )
!---------- ---------- ---------- ---------- ----------


! local
!---------- ---------- ---------- ---------- ----------
DIMENSION XT(NDIM,NNODE), DJ(NDIM,NDIM), DJI(NDIM,NDIM), FN(NNODE), DFX(NNODE,NDIM), DFXI(NNODE,NDIM)
DIMENSION XINT(9,9), WINT(9,9)

!DIMENSION H1(NNODE, NNODE), H2(NNODE, NNODE), H3(NNODE, NNODE), H4(NNODE, NNODE), H5(NNODE, NNODE), H6(NNODE, NNODE), H7(NNODE, NNODE), H8(NNODE, NNODE), H9(NNODE, NNODE)
!DIMENSION G1(NNODE, NNODE), G2(NNODE, NNODE), G3(NNODE, NNODE), G4(NNODE, NNODE), G5(NNODE, NNODE), G6(NNODE, NNODE), G7(NNODE, NNODE), G8(NNODE, NNODE), G9(NNODE, NNODE)
!---------- ---------- ---------- ---------- ----------


CALL GAUSS ( XINT, WINT )


! load properties
!---------- ---------- ---------- ---------- ----------
NINT = NGP(IEL)


! define factors and initialize
!---------- ---------- ---------- ---------- ---------- 
E_core_Mu         = 0.0d0
E_core_Lambda     = 0.0d0
E_g_Mu            = 0.0d0
E_g_Lambda_2Mu    = 0.0d0


DO I = 1, NNODE
   K = INOD(I,IEL)
   DO J = 1, NDIM
      XT(J,I) = XYZ(K,J)
   END DO
END DO


! Gauss integration        
!---------- ---------- ---------- ---------- ----------
DO LZ = 1, NINT
   X3 = XINT(LZ,NINT)
   WZ = WINT(LZ,NINT)

   DO LY = 1, NINT
      X2 = XINT(LY,NINT)
      WY = WINT(LY,NINT)

      DO LX = 1, NINT
         X1 = XINT(LX,NINT)
         WX = WINT(LX,NINT)

         WSTAR = WX * WY * WZ


! Jacobian and shape functions    
!---------- ---------- ---------- ---------- ----------         
         SELECT CASE (NNODE)

            CASE (20)                            ! 20-noded elements (serendipity)
               CALL SHAPE_3_20 (X1, X2, X3, FN, DFXI)

            CASE (27)                            ! 27-noded elements (Lagrange)
               CALL SHAPE_3_27 (X1, X2, X3, FN, DFXI)

            CASE DEFAULT
               WRITE(*,'("FLAG ERROR IN EK_3D - NNODE")')
               STOP

         END SELECT


         DJ = MATMUL(XT,DFXI)

         DETJ = DJ(1,1) * DJ(2,2) * DJ(3,3) + DJ(2,1) * DJ(3,2) * DJ(1,3) + DJ(3,1) * DJ(1,2) * DJ(2,3) - &
                DJ(1,3) * DJ(2,2) * DJ(3,1) - DJ(2,3) * DJ(3,2) * DJ(1,1) - DJ(3,3) * DJ(1,2) * DJ(2,1)

         IF (DETJ.LE.0) THEN
            WRITE(*,*)'|J|<0 NEGATIVE DET. JACOBIAN'
            READ(*,*)
            STOP
         END IF

         DJI(1,1) = + DJ(3,3) * DJ(2,2) - DJ(3,2) * DJ(2,3)
         DJI(1,2) = - DJ(3,3) * DJ(1,2) + DJ(3,2) * DJ(1,3)
         DJI(1,3) = + DJ(2,3) * DJ(1,2) - DJ(2,2) * DJ(1,3)
 
         DJI(2,1) = - DJ(3,3) * DJ(2,1) + DJ(3,1) * DJ(2,3)
         DJI(2,2) = + DJ(3,3) * DJ(1,1) - DJ(3,1) * DJ(1,3)
         DJI(2,3) = - DJ(2,3) * DJ(1,1) + DJ(2,1) * DJ(1,3)

         DJI(3,1) = + DJ(3,2) * DJ(2,1) - DJ(3,1) * DJ(2,2)
         DJI(3,2) = - DJ(3,2) * DJ(1,1) + DJ(3,1) * DJ(1,2)
         DJI(3,3) = + DJ(2,2) * DJ(1,1) - DJ(2,1) * DJ(1,2)

         DJI = DJI / DETJ


         DFX = MATMUL(DFXI,DJI)
!---------- ---------- ---------- ---------- ---------- 
         FAC0 = WSTAR * DETJ


! form the components
!---------- ---------- ---------- ---------- ----------
         u1_x = Dot_Product ( E_U (            1 :     NNODE ) , DFX(:,1) )
         u2_x = Dot_Product ( E_U (     NNODE +1 : 2 * NNODE ) , DFX(:,1) )
         u3_x = Dot_Product ( E_U ( 2 * NNODE +1 : 3 * NNODE ) , DFX(:,1) )
         u1_y = Dot_Product ( E_U (            1 :     NNODE ) , DFX(:,2) )
         u2_y = Dot_Product ( E_U (     NNODE +1 : 2 * NNODE ) , DFX(:,2) )
         u3_y = Dot_Product ( E_U ( 2 * NNODE +1 : 3 * NNODE ) , DFX(:,2) )
         u1_z = Dot_Product ( E_U (            1 :     NNODE ) , DFX(:,3) )
         u2_z = Dot_Product ( E_U (     NNODE +1 : 2 * NNODE ) , DFX(:,3) )
         u3_z = Dot_Product ( E_U ( 2 * NNODE +1 : 3 * NNODE ) , DFX(:,3) )

         p1_x = Dot_Product ( E_P (            1 :     NNODE ) , DFX(:,1) )
         p2_x = Dot_Product ( E_P (     NNODE +1 : 2 * NNODE ) , DFX(:,1) )
         p3_x = Dot_Product ( E_P ( 2 * NNODE +1 : 3 * NNODE ) , DFX(:,1) )
         p1_y = Dot_Product ( E_P (            1 :     NNODE ) , DFX(:,2) )
         p2_y = Dot_Product ( E_P (     NNODE +1 : 2 * NNODE ) , DFX(:,2) )
         p3_y = Dot_Product ( E_P ( 2 * NNODE +1 : 3 * NNODE ) , DFX(:,2) )
         p1_z = Dot_Product ( E_P (            1 :     NNODE ) , DFX(:,3) )
         p2_z = Dot_Product ( E_P (     NNODE +1 : 2 * NNODE ) , DFX(:,3) )
         p3_z = Dot_Product ( E_P ( 2 * NNODE +1 : 3 * NNODE ) , DFX(:,3) )
!---------- ---------- ---------- ---------- ----------

         scalar_Lambda_2Mu = ( u1_x + u2_y + u3_z ) * ( p1_x + p2_y + p3_z )
         scalar_Mu     = -2.0d0 * ( p1_x * ( u2_y + u3_z ) + p2_y * ( u1_x + u3_z ) + p3_z * ( u1_x + u2_y ) ) &
                       + ( p1_y + p2_x ) * ( u1_y + u2_x ) &
                       + ( p1_z + p3_x ) * ( u1_z + u3_x ) &
                       + ( p3_y + p2_z ) * ( u3_y + u2_z )

         E_g_Mu         = E_g_Mu         + FN * scalar_Mu         * FAC0
         E_g_Lambda_2Mu = E_g_Lambda_2Mu + FN * scalar_Lambda_2Mu * FAC0

      END DO
   END DO
END DO
!---------- ---------- ---------- ---------- ----------


RETURN
END
