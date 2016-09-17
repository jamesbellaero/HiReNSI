!************************************************
! M-PML parameters at each Gauss point - 2D
!************************************************
! REVISION : M 4 Feb 2013

SUBROUTINE PML_FINDER ( XG, PML_PARAM, PML_ALPHA_BETA )
USE PARAMETERS
IMPLICIT DOUBLE PRECISION (A-H,O-Z)


! in
!---------- ---------- ---------- ---------- ----------
DIMENSION XG(NDIM), PML_PARAM(NDIM * 2 , 4)
!---------- ---------- ---------- ---------- ----------


! out
!---------- ---------- ---------- ---------- ----------	  
DIMENSION PML_ALPHA_BETA(19)
!---------- ---------- ---------- ---------- ----------


! initialize
!---------- ---------- ---------- ---------- ----------
PML_alpha_beta     = 0.0d0
PML_alpha_beta (1) = 1.0d0                       ! alpha_x (just part of the real part)
PML_alpha_beta (2) = 1.0d0                       ! alpha_y (just part of the real part)


! load PML parameters
!---------- ---------- ---------- ---------- ----------
! PML alpha_0
PML_alpha01 = PML_PARAM(1,1)
PML_alpha02 = PML_PARAM(2,1)
PML_alpha03 = PML_PARAM(3,1)
PML_alpha04 = PML_PARAM(4,1)
  
! PML beta_0
PML_beta01  = PML_PARAM(1,2)
PML_beta02  = PML_PARAM(2,2)
PML_beta03  = PML_PARAM(3,2)
PML_beta04  = PML_PARAM(4,2)

! PML length
PML_L1      = PML_PARAM(1,3) 
PML_L2      = PML_PARAM(2,3) 
PML_L3      = PML_PARAM(3,3) 
PML_L4      = PML_PARAM(4,3) 
  
! PML starting point location	  
PML_X01     = PML_PARAM(1,4)
PML_X02     = PML_PARAM(2,4)
PML_X03     = PML_PARAM(3,4)
PML_X04     = PML_PARAM(4,4)

! PML polynomial degree
I_PML_trig            = 0                                  ! (0) polynomial / (1) trigonometric functions, for PML profile
pi                    = 3.14159265358979323846264338327950288419716939937510d0
PML_polynomial_degree = 2.0d0
ratio                 = 0.00d0 !10.000d0 / 262.50d0        ! this is ratio/H in corrected M-PML; H characteristic length (e.g. largest dimension of problem)
!---------- ---------- ---------- ---------- ---------- 


! initialize diff = (s - s0), fac<1/2> = (s - s0) n / L, fac<3/4> = n / L
!---------- ---------- ---------- ---------- ----------
diffx       = 0.0d0
diffy       = 0.0d0

vnx         = 0.0d0
vny         = 0.0d0

FAC1        = 0.0d0
FAC2        = 0.0d0
FAC3        = 0.0d0
FAC4        = 0.0d0

PML_alpha0x = 0.0d0
PML_alpha0y = 0.0d0 
PML_beta0x  = 0.0d0
PML_beta0y  = 0.0d0

LOCx        = 0
LOCy        = 0
Lreg        = 0


! find the position of the Gauss point (right or left zones)
!---------- ---------- ---------- ---------- ----------
! the point XG(1) can be either in right or left zone (or neither)
IF (             XG(1) >= PML_X01 ) THEN         ! right PML zone
   diffx       = XG(1)  - PML_X01
   vnx         = 1.0d0
   PML_alpha0x = PML_alpha01
   PML_beta0x  = PML_beta01
   PML_Lx      = PML_L1
   FAC1        = diffx * vnx / PML_Lx
   FAC3        =         vnx / PML_Lx
   LOCx        = 1

ELSEIF (         XG(1) <= PML_X02 ) THEN         ! left PML zone
   diffx       = XG(1)  - PML_X02
   vnx         = -1.0d0
   PML_alpha0x = PML_alpha02
   PML_beta0x  = PML_beta02
   PML_Lx      = PML_L2
   FAC1        = diffx * vnx / PML_Lx
   FAC3        =         vnx / PML_Lx
   LOCx        = 2

END IF


! (top or bottom)
!---------- ---------- ---------- ---------- ----------
! the point XG(1) can be either in top or bottom zone (or neither)
! IF (              XG(2) > PML_X03 ) THEN       ! top PML zone
!     diffy       = XG(2) - PML_X03
!     ny          = + 1.0
!     PML_alpha0y = PML_alpha03
!     PML_beta0y  = PML_beta03
!     PML_Ly      = PML_L3                    
!     FAC2        = diffy * ny / PML_Ly

IF (             XG(2) <= PML_X04 ) THEN         ! bottom PML zone
   diffy       = XG(2)  - PML_X04
   vny         = -1.0d0
   PML_alpha0y = PML_alpha04
   PML_beta0y  = PML_beta04
   PML_Ly      = PML_L4
   FAC2        = diffy * vny / PML_Ly
   FAC4        =         vny / PML_Ly
   LOCy        = 4

END IF


! find the region that we are in it
!---------- ---------- ---------- ---------- ----------
! -------------------------
! |       |       |       |
! |   5   |   RD  |   1   |
! |       |       |       |
! -------------------------
! |       |       |       |
! |   4   |   3   |   2   |
! |       |       |       |
! -------------------------
IF     ( LOCx == 1 .AND. LOCy == 0 ) THEN
   Lreg = 1
ELSEIF ( LOCx == 1 .AND. LOCy == 4 ) THEN
   Lreg = 2
ELSEIF ( LOCx == 0 .AND. LOCy == 4 ) THEN
   Lreg = 3
ELSEIF ( LOCx == 2 .AND. LOCy == 4 ) THEN
   Lreg = 4
ELSEIF ( LOCx == 2 .AND. LOCy == 0 ) THEN
   Lreg = 5
ELSE
   write(*,*) 'Error in PML finder subroutine'
   write(*,'(A10,2F10.5)') 'XG =', XG(1), XG(2)
   write(*,*) 'LOC-x,-y =', LOCx, LOCy
   STOP
END IF


alpha_jump = 0.0d0
 beta_jump = 0.0d0


SELECT CASE ( I_PML_trig )

   CASE(0)

!write(*,*) 'classic pml'

   !SELECT CASE (Lreg) 

   !   CASE (1, 5)
   if ( Lreg == 1 .or. Lreg == 2 .or. Lreg == 4 .or. Lreg == 5) then
! compute PML parameters (honest: right-left)
!---------- ---------- ---------- ---------- ----------
         ! alpha_x
         PML_alpha_beta (1) = PML_alpha_beta (1) + alpha_jump + PML_alpha0x * FAC1 ** PML_polynomial_degree
         ! alpha_y
         PML_alpha_beta (2) = PML_alpha_beta (2) + (PML_alpha0x * FAC1 ** PML_polynomial_degree) * ratio
         ! beta_x
         PML_alpha_beta (3) = PML_alpha_beta (3) + beta_jump +  PML_beta0x * FAC1 ** PML_polynomial_degree
         ! beta_y
         PML_alpha_beta (4) = PML_alpha_beta (4) +  (PML_beta0x * FAC1 ** PML_polynomial_degree) * ratio
         ! d/dx alpha_y
         PML_alpha_beta (8) = (PML_alpha0x * PML_polynomial_degree * FAC3 * FAC1 ** (PML_polynomial_degree - 1)) * ratio
        ! d/dx beta_y
         PML_alpha_beta (10)=  (PML_beta0x * PML_polynomial_degree * FAC3 * FAC1 ** (PML_polynomial_degree - 1)) * ratio
         ! alpha_y_xt (y - x tilde)
         PML_alpha_beta (13)= - PML_alpha_beta (8)  * XG(2) / XG(1) * diffx
         ! d/dy alpha_y_xt
         PML_alpha_beta (17)= - PML_alpha_beta (8)
         ! beta_y_xt
         PML_alpha_beta (15)= - PML_alpha_beta (10) * XG(2) / XG(1) * diffx
         ! d/dy beta_y_xt
         PML_alpha_beta (19)= - PML_alpha_beta (10)

   end if

   !   CASE (3)
   if ( Lreg == 2 .or. Lreg == 3 .or. Lreg == 4 ) then
   !if (                Lreg == 3                ) then
! compute PML parameters (honest: bottom)
!---------- ---------- ---------- ---------- ----------
         ! alpha_x
         PML_alpha_beta (1) = PML_alpha_beta (1) + (PML_alpha0y * FAC2 ** PML_polynomial_degree) * ratio
         ! alpha_y
         PML_alpha_beta (2) = PML_alpha_beta (2) + alpha_jump + PML_alpha0y * FAC2 ** PML_polynomial_degree
         ! beta_x
         PML_alpha_beta (3) = PML_alpha_beta (3) +  (PML_beta0y * FAC2 ** PML_polynomial_degree) * ratio
         ! beta_y
         PML_alpha_beta (4) = PML_alpha_beta (4) + beta_jump + PML_beta0y * FAC2 ** PML_polynomial_degree
         ! d/dy alpha_x
         PML_alpha_beta (9) = (PML_alpha0y * PML_polynomial_degree * FAC4 * FAC2 ** (PML_polynomial_degree - 1)) * ratio
         ! d/dy beta_x
         PML_alpha_beta (11)=  (PML_beta0y * PML_polynomial_degree * FAC4 * FAC2 ** (PML_polynomial_degree - 1)) * ratio
         ! alpha_x_yt
         PML_alpha_beta (12)= - PML_alpha_beta (9)  * XG(1) / XG(2) * diffy
         ! d/dx  alpha_x_yt
         PML_alpha_beta (16)= - PML_alpha_beta (9)
         ! beta_x_yt
         PML_alpha_beta (14)= - PML_alpha_beta (11) * XG(1) / XG(2) * diffy
         ! d/dx  beta_x_yt
         PML_alpha_beta (18)= - PML_alpha_beta (11)

   end if

   !   CASE (2, 4)
! compute PML parameters (simply pml: corners)
!---------- ---------- ---------- ---------- ----------
         ! alpha_x
!         PML_alpha_beta (1) = PML_alpha_beta (1) +  PML_alpha0x * FAC1 ** PML_polynomial_degree
         ! alpha_y
!         PML_alpha_beta (2) = PML_alpha_beta (2) +  PML_alpha0y * FAC2 ** PML_polynomial_degree
         ! beta_x
!         PML_alpha_beta (3) = PML_alpha_beta (3) +   PML_beta0x * FAC1 ** PML_polynomial_degree
         ! beta_y
!         PML_alpha_beta (4) = PML_alpha_beta (4) +   PML_beta0y * FAC2 ** PML_polynomial_degree


!      CASE DEFAULT
!         WRITE(*,'("FLAG ERROR IN PML 2D")')
!         STOP


!   END SELECT

   CASE(1)                                       ! trigonometric PML profile

write(*,*) 'we should not be here: pml.f90' ; stop

   if ( Lreg == 1 .or. Lreg == 2 .or. Lreg == 4 .or. Lreg == 5) then
! compute PML parameters (honest: right-left)
!---------- ---------- ---------- ---------- ----------
         ! alpha_x
         PML_alpha_beta (1) = PML_alpha_beta (1) + alpha_jump + PML_alpha0x * 0.50d0 * ( 1.0d0 + dsin( pi * ( dabs(FAC1) - 0.50d0 ) ) ) 
         ! alpha_y
!         PML_alpha_beta (2) = PML_alpha_beta (2) + (PML_alpha0x * FAC1 ** PML_polynomial_degree) * ratio
         ! beta_x
         PML_alpha_beta (3) = PML_alpha_beta (3) + beta_jump +  PML_beta0x * 0.50d0 * ( 1.0d0 + dsin( pi * ( dabs(FAC1) - 0.50d0 ) ) )
         ! beta_y
!         PML_alpha_beta (4) = PML_alpha_beta (4) +  (PML_beta0x * FAC1 ** PML_polynomial_degree) * ratio
   end if

   if ( Lreg == 2 .or. Lreg == 3 .or. Lreg == 4 ) then
! compute PML parameters (honest: bottom)
!---------- ---------- ---------- ---------- ----------
         ! alpha_x
!         PML_alpha_beta (1) = PML_alpha_beta (1) + (PML_alpha0y * FAC2 ** PML_polynomial_degree) * ratio
         ! alpha_y
         PML_alpha_beta (2) = PML_alpha_beta (2) + alpha_jump + PML_alpha0y * 0.50d0 * ( 1.0d0 + dsin( pi * ( dabs(FAC2) - 0.50d0 ) ) )
         ! beta_x
!         PML_alpha_beta (3) = PML_alpha_beta (3) +  (PML_beta0y * FAC2 ** PML_polynomial_degree) * ratio
         ! beta_y
         PML_alpha_beta (4) = PML_alpha_beta (4) + beta_jump +  PML_beta0y * 0.50d0 * ( 1.0d0 + dsin( pi * ( dabs(FAC2) - 0.50d0 ) ) )
   end if
 
END SELECT


! compute PML parameters
!---------- ---------- ---------- ---------- ----------
! a = alpha_x * alpha_y
PML_alpha_beta (5) =  PML_alpha_beta (1) * PML_alpha_beta (2)
! b = alpha_x * beta_y + alpha_y * beta_x
PML_alpha_beta (6) =  PML_alpha_beta (1) * PML_alpha_beta (4) + PML_alpha_beta (2) * PML_alpha_beta (3)
! c = beta_x * beta_y
PML_alpha_beta (7) =  PML_alpha_beta (3) * PML_alpha_beta (4)  




!write(*,*) '------------------------------------------------'
!!   write(*,'(A10,2F10.5)') 'XG =', XG(1), XG(2)
!   write(*,*) 'LOC-x,-y =', LOCx, LOCy
!write(*,*) '------------------------------------------------'
!write(*,*) lreg
!write(*,*) '------------------------------------------------'
!!write(*,*) pml_param(:,2)
!write(*,*) '------------------------------------------------'
!write(*,*) pml_param(:,3)
!write(*,*) '------------------------------------------------'
!write(*,*) pml_param(:,4)
!write(*,*) '------------------------------------------------'
!write(*,*) PML_alpha_beta(1:7)
!write(*,*) '------------------------------------------------'
!pause


RETURN
END


!************************************************
! M-PML parameters at each Gauss point - 3D
!************************************************
! REVISION : M 18 March 2013

SUBROUTINE PML_FINDER_3D ( XG, PML_PARAM, PML_ALPHA_BETA )
USE PARAMETERS
IMPLICIT DOUBLE PRECISION (A-H,O-Z)


! in
!---------- ---------- ---------- ---------- ----------
DIMENSION XG(NDIM), PML_PARAM(NDIM * 2 , 4)
!---------- ---------- ---------- ---------- ----------


! out
!---------- ---------- ---------- ---------- ----------   
DIMENSION PML_ALPHA_BETA(19)
!---------- ---------- ---------- ---------- ----------


! initialize
!---------- ---------- ---------- ---------- ----------
PML_alpha_beta     = 0.0d0
PML_alpha_beta (1) = 1.0d0                       ! alpha_x (just part of the real part)
PML_alpha_beta (2) = 1.0d0                       ! alpha_y (just part of the real part)
PML_alpha_beta (3) = 1.0d0                       ! alpha_z (just part of the real part)


! load PML parameters
!---------- ---------- ---------- ---------- ----------
! PML alpha_0
PML_alpha01 = PML_PARAM(1,1)                     ! x = +1
PML_alpha02 = PML_PARAM(2,1)                     ! x = -1
PML_alpha03 = PML_PARAM(3,1)                     ! y = +1
PML_alpha04 = PML_PARAM(4,1)                     ! y = -1
PML_alpha05 = PML_PARAM(5,1)                     ! z = +1
PML_alpha06 = PML_PARAM(6,1)                     ! z = -1

! PML beta_0
PML_beta01  = PML_PARAM(1,2)
PML_beta02  = PML_PARAM(2,2)
PML_beta03  = PML_PARAM(3,2)
PML_beta04  = PML_PARAM(4,2)
PML_beta05  = PML_PARAM(5,2)
PML_beta06  = PML_PARAM(6,2)

! PML length
PML_L1      = PML_PARAM(1,3)
PML_L2      = PML_PARAM(2,3)
PML_L3      = PML_PARAM(3,3)
PML_L4      = PML_PARAM(4,3)
PML_L5      = PML_PARAM(5,3)
PML_L6      = PML_PARAM(6,3)

! PML starting point location         
PML_X01     = PML_PARAM(1,4)
PML_X02     = PML_PARAM(2,4)
PML_X03     = PML_PARAM(3,4)
PML_X04     = PML_PARAM(4,4)
PML_X05     = PML_PARAM(5,4)
PML_X06     = PML_PARAM(6,4)

! PML polynomial degree
PML_polynomial_degree = 2.00d0
ratio                 = 0.10d0
I_PML_trig            = 0                                  ! (0) polynomial / (1) trigonometric functions, for PML profile
pi                    = 3.14159265358979323846264338327950288419716939937510d0
!---------- ---------- ---------- ---------- ----------     


! initialize diff = (s - s0), fac<1,2,3> = (s - s0) n / L, fac<4,5,6> = n / L
!---------- ---------- ---------- ---------- ----------   
diffx          = 0.0d0
diffy          = 0.0d0
diffz          = 0.0d0

vnx            = 0.0d0
vny            = 0.0d0
vnz            = 0.0d0

FAC1           = 0.0d0
FAC2           = 0.0d0
FAC3           = 0.0d0
FAC4           = 0.0d0
FAC5           = 0.0d0
FAC6           = 0.0d0

PML_alpha0x    = 0.0d0
PML_alpha0y    = 0.0d0
PML_alpha0z    = 0.0d0
PML_beta0x     = 0.0d0
PML_beta0y     = 0.0d0
PML_beta0z     = 0.0d0

LOCx           = 0
LOCy           = 0
LOCz           = 0
Lreg           = 0

d_dy_alpha_x   = 0.0d0
d_dz_alpha_x   = 0.0d0
d_dx_alpha_y   = 0.0d0
d_dz_alpha_y   = 0.0d0
d_dx_alpha_z   = 0.0d0
d_dy_alpha_z   = 0.0d0

d_dy_beta_x    = 0.0d0
d_dz_beta_x    = 0.0d0
d_dx_beta_y    = 0.0d0
d_dz_beta_y    = 0.0d0
d_dx_beta_z    = 0.0d0
d_dy_beta_z    = 0.0d0
!goto 100
! FIND WHERE THE POINT IS LOCATED
! the point XG(1) can be either in the right or the left zone (or neither)
!---------- ---------- ---------- ---------- ----------
IF (             XG(1) >= PML_X01 ) THEN         ! right PML zone
   diffx       = XG(1)  - PML_X01
   vnx         = 1.0d0
   PML_alpha0x = PML_alpha01
   PML_beta0x  = PML_beta01
   PML_Lx      = PML_L1
   LOCx        = 1
   FAC1 = diffx * vnx / PML_Lx
   FAC4 =         vnx / PML_Lx

      ! alpha_x
      PML_alpha_beta (1) = PML_alpha_beta (1) +  PML_alpha0x * FAC1 ** PML_polynomial_degree
      ! beta_x
      PML_alpha_beta (4) = PML_alpha_beta (4) +  PML_beta0x  * FAC1 ** PML_polynomial_degree


ELSEIF (         XG(1) <= PML_X02 ) THEN         ! left PML zone
   diffx       = XG(1)  - PML_X02
   vnx         = -1.0d0
   PML_alpha0x = PML_alpha02
   PML_beta0x  = PML_beta02
   PML_Lx      = PML_L2
   LOCx        = -1
   FAC1 = diffx * vnx / PML_Lx
   FAC4 =         vnx / PML_Lx

      ! alpha_x
      PML_alpha_beta (1) = PML_alpha_beta (1) +  PML_alpha0x * FAC1 ** PML_polynomial_degree
      ! beta_x
      PML_alpha_beta (4) = PML_alpha_beta (4) +  PML_beta0x  * FAC1 ** PML_polynomial_degree


END IF

!goto 100
! the point XG(2) can be either at the front or the rear zone (or neither)
!---------- ---------- ---------- ---------- ----------
IF (             XG(2) >= PML_X03 ) THEN         ! front PML zone
   diffy       = XG(2)  - PML_X03
   vny         = 1.0d0
   PML_alpha0y = PML_alpha03
   PML_beta0y  = PML_beta03
   PML_Ly      = PML_L3                    
   LOCy        = 1
   FAC2 = diffy * vny / PML_Ly
   FAC5 =         vny / PML_Ly

      ! alpha_y
      PML_alpha_beta (2) = PML_alpha_beta (2) +  PML_alpha0y * FAC2 ** PML_polynomial_degree
      ! beta_y
      PML_alpha_beta (5) = PML_alpha_beta (5) +  PML_beta0y  * FAC2 ** PML_polynomial_degree


ELSEIF (         XG(2) <= PML_X04 ) THEN         ! rear PML zone
   diffy       = XG(2)  - PML_X04
   vny         = -1.0d0
   PML_alpha0y = PML_alpha04
   PML_beta0y  = PML_beta04
   PML_Ly      = PML_L4
   LOCy        = -1
   FAC2 = diffy * vny / PML_Ly
   FAC5 =         vny / PML_Ly

      ! alpha_y
      PML_alpha_beta (2) = PML_alpha_beta (2) +  PML_alpha0y * FAC2 ** PML_polynomial_degree
      ! beta_y
      PML_alpha_beta (5) = PML_alpha_beta (5) +  PML_beta0y  * FAC2 ** PML_polynomial_degree


END IF

100 continue
! the point XG(3) can be either on the top or at the bottom zone (or neither)
!---------- ---------- ---------- ---------- ----------
IF (             XG(3) >  PML_X05 ) THEN         ! top PML zone
   diffz       = XG(3)  - PML_X05
   vnz         = 1.0d0
   PML_alpha0z = PML_alpha05
   PML_beta0z  = PML_beta05
   PML_Lz      = PML_L5
   LOCz        = 1
   FAC3 = diffz * vnz / PML_Lz
   FAC6 =         vnz / PML_Lz
   write(*,*) 'Error in PML finder-3D subroutine: top PML zone detected'
   STOP

ELSEIF (         XG(3) <= PML_X06 ) THEN         ! bottom PML zone
   diffz       = XG(3)  - PML_X06
   vnz         = -1.0d0
   PML_alpha0z = PML_alpha06
   PML_beta0z  = PML_beta06
   PML_Lz      = PML_L6
   LOCz        = -1
   FAC3 = diffz * vnz / PML_Lz
   FAC6 =         vnz / PML_Lz


      ! alpha_z
      PML_alpha_beta (3) = PML_alpha_beta (3) +  PML_alpha0z * FAC3 ** PML_polynomial_degree
      ! beta_z
      PML_alpha_beta (6) = PML_alpha_beta (6) +  PML_beta0z  * FAC3 ** PML_polynomial_degree


END IF


! find the region that we are in it. each "plane" indicates one "volumetric" pml zone
!---------- ---------- ---------- ---------- ----------
!
!                     -------------------------
!                    /|       |       |       |
!                   / |   14  |   13  |   12  |
!                  /  |       |       |       |
!                 /   -------------------------
!                /    |       |       |       |
!               /     |   17  |   16  |   15  |
!              /      |       |       |       |
!             /       ------------------------- (rear zone: LOCx = -1)
!            /
!           -------------------------
!          /|       |       |       |
!         / |   8   |   RD  |   7   |
!        /  |       |       |       |
!       /   -------------------------
!      /    |       |       |       |
!     /     |   11  |   10  |   9   |
!    /      |       |       |       |
!   /       ------------------------- (middle zone: LOCx = 0)
!  /                       
! -------------------------
! |       |       |       |
! |   3   |   2   |   1   |
! |       |       |       |
! -------------------------
! |       |       |       |
! |   6   |   5   |   4   |
! |       |       |       |
! ------------------------- (front zone: LOCx = +1)

IF     ( LOCx ==  1 .AND. LOCy ==  1 .AND. LOCz ==  0 ) THEN
                                                             Lreg = 1
ELSEIF ( LOCx ==  1 .AND. LOCy ==  0 .AND. LOCz ==  0 ) THEN
                                                             Lreg = 2
ELSEIF ( LOCx ==  1 .AND. LOCy == -1 .AND. LOCz ==  0 ) THEN
                                                             Lreg = 3
ELSEIF ( LOCx ==  1 .AND. LOCy ==  1 .AND. LOCz == -1 ) THEN
                                                             Lreg = 4
ELSEIF ( LOCx ==  1 .AND. LOCy ==  0 .AND. LOCz == -1 ) THEN
                                                             Lreg = 5
ELSEIF ( LOCx ==  1 .AND. LOCy == -1 .AND. LOCz == -1 ) THEN
                                                             Lreg = 6
ELSEIF ( LOCx ==  0 .AND. LOCy ==  1 .AND. LOCz ==  0 ) THEN
                                                             Lreg = 7
ELSEIF ( LOCx ==  0 .AND. LOCy == -1 .AND. LOCz ==  0 ) THEN
                                                             Lreg = 8
ELSEIF ( LOCx ==  0 .AND. LOCy ==  1 .AND. LOCz == -1 ) THEN
                                                             Lreg = 9
ELSEIF ( LOCx ==  0 .AND. LOCy ==  0 .AND. LOCz == -1 ) THEN
                                                             Lreg = 10
ELSEIF ( LOCx ==  0 .AND. LOCy == -1 .AND. LOCz == -1 ) THEN
                                                             Lreg = 11
ELSEIF ( LOCx == -1 .AND. LOCy ==  1 .AND. LOCz ==  0 ) THEN
                                                             Lreg = 12
ELSEIF ( LOCx == -1 .AND. LOCy ==  0 .AND. LOCz ==  0 ) THEN
                                                             Lreg = 13
ELSEIF ( LOCx == -1 .AND. LOCy == -1 .AND. LOCz ==  0 ) THEN
                                                             Lreg = 14
ELSEIF ( LOCx == -1 .AND. LOCy ==  1 .AND. LOCz == -1 ) THEN
                                                             Lreg = 15
ELSEIF ( LOCx == -1 .AND. LOCy ==  0 .AND. LOCz == -1 ) THEN
                                                             Lreg = 16
ELSEIF ( LOCx == -1 .AND. LOCy == -1 .AND. LOCz == -1 ) THEN
                                                             Lreg = 17
ELSE
!   write(*,*) 'Error in PML finder subroutine'
!   write(*,'(A10,3F30.20)') 'XG =', XG(1), XG(2), XG(3)
!   write(*,*) 'LOC-x,-y,-z =', LOCx, LOCy, LOCz
!   STOP
!pause
END IF


! compute MPML alpha and beta parameters and store in an array.
!---------- ---------- ---------- ---------- ----------
! [ alpha_x, alpha_y, alpha_z, beta_x, beta_y, beta_z, a, b, c, d,                                              {1 2 3 - 4 5 6 - 7 8 9 10}
! d (alpha_x alpha_y), d (alpha_x alpha_z), d (alpha_y alpha_z),                                                {11 12 13}
! d (alpha_x beta_y + alpha_y beta_x), d (alpha_x beta_z + alpha_z beta_x),d (alpha_y beta_z + alpha_z beta_y), {14 15 16}
! d (beta_x beta_y), d (beta_x beta_z), d (beta_y beta_z) ]                                                     {17 18 19}


SELECT CASE ( I_PML_trig )

   CASE(0) ! polynomial PML
! ===================================================================================================================================================
!stop

if ( Lreg == 1 .or. Lreg == 2 .or. Lreg == 3 .or. Lreg == 4 .or. Lreg == 5 .or. Lreg == 6 .or. &
     Lreg == 12 .or. Lreg == 13 .or. Lreg == 14 .or. Lreg == 15 .or. Lreg == 16 .or. Lreg == 17) then
! compute PML parameters (X-strip)
!---------- ---------- ---------- ---------- ----------

end if


if ( Lreg == 1 .or. Lreg == 3 .or. Lreg == 4 .or. Lreg == 6 .or. Lreg == 7 .or. Lreg == 8 .or. &
     Lreg == 9 .or. Lreg == 11 .or. Lreg == 12 .or. Lreg == 14 .or. Lreg == 15 .or. Lreg == 17) then
! compute PML parameters (Y-strip)
!---------- ---------- ---------- ---------- ----------

end if


if ( Lreg == 4 .or. Lreg == 5 .or. Lreg == 6 .or. Lreg == 9 .or. Lreg == 10 .or. Lreg == 11 .or. Lreg == 15 .or. Lreg == 16 .or. Lreg == 17) then
! compute PML parameters (Z-strip)
!---------- ---------- ---------- ---------- ----------

end if


   CASE(1) ! trig PML <><><><><><><><><><><><><><><><><><><><><><><><>()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()
! ===================================================================================================================================================
!stop

if ( Lreg == 1 .or. Lreg == 2 .or. Lreg == 3 .or. Lreg == 4 .or. Lreg == 5 .or. Lreg == 6 .or. &
     Lreg == 12 .or. Lreg == 13 .or. Lreg == 14 .or. Lreg == 15 .or. Lreg == 16 .or. Lreg == 17) then
! compute PML parameters (X-strip)
!---------- ---------- ---------- ---------- ----------
      ! alpha_x
      PML_alpha_beta (1) = PML_alpha_beta (1) +  PML_alpha0x * 0.50d0 * ( 1.0d0 + dsin( pi * ( dabs(FAC1) - 0.50d0 ) ) )
      ! beta_x
      PML_alpha_beta (4) = PML_alpha_beta (4) +  PML_beta0x  * 0.50d0 * ( 1.0d0 + dsin( pi * ( dabs(FAC1) - 0.50d0 ) ) )
end if


if ( Lreg == 1 .or. Lreg == 3 .or. Lreg == 4 .or. Lreg == 6 .or. Lreg == 7 .or. Lreg == 8 .or. &
     Lreg == 9 .or. Lreg == 11 .or. Lreg == 12 .or. Lreg == 14 .or. Lreg == 15 .or. Lreg == 17) then
! compute PML parameters (Y-strip)
!---------- ---------- ---------- ---------- ----------
      ! alpha_y
      PML_alpha_beta (2) = PML_alpha_beta (2) +  PML_alpha0y * 0.50d0 * ( 1.0d0 + dsin( pi * ( dabs(FAC2) - 0.50d0 ) ) )
      ! beta_y
      PML_alpha_beta (5) = PML_alpha_beta (5) +  PML_beta0y  * 0.50d0 * ( 1.0d0 + dsin( pi * ( dabs(FAC2) - 0.50d0 ) ) )
end if


if ( Lreg == 4 .or. Lreg == 5 .or. Lreg == 6 .or. Lreg == 9 .or. Lreg == 10 .or. Lreg == 11 .or. Lreg == 15 .or. Lreg == 16 .or. Lreg == 17) then
! compute PML parameters (Z-strip)
!---------- ---------- ---------- ---------- ----------
      ! alpha_z
      PML_alpha_beta (3) = PML_alpha_beta (3) +  PML_alpha0z * 0.50d0 * ( 1.0d0 + dsin( pi * ( dabs(FAC3) - 0.50d0 ) ) )
      ! beta_z
      PML_alpha_beta (6) = PML_alpha_beta (6) +  PML_beta0z  * 0.50d0 * ( 1.0d0 + dsin( pi * ( dabs(FAC3) - 0.50d0 ) ) )
end if


! ===================================================================================================================================================
End Select

! a = alpha_x * alpha_y * alpha_z
!------------------------------------------------
PML_alpha_beta (7) =  PML_alpha_beta (1) * PML_alpha_beta (2) * PML_alpha_beta (3)


! b = alpha_x * alpha_y * beta_z + alpha_x * alpha_z * beta_y + alpha_y * alpha_z * beta_x
!------------------------------------------------
PML_alpha_beta (8) =  PML_alpha_beta (1) * PML_alpha_beta (2) * PML_alpha_beta (6) + &
                      PML_alpha_beta (1) * PML_alpha_beta (3) * PML_alpha_beta (5) + &
                      PML_alpha_beta (2) * PML_alpha_beta (3) * PML_alpha_beta (4)


! c = alpha_x * beta_y * beta_z + alpha_y * beta_x * beta_z + alpha_z * beta_x * beta_y
!------------------------------------------------
PML_alpha_beta (9) =  PML_alpha_beta (1) * PML_alpha_beta (5) * PML_alpha_beta (6) + &
                      PML_alpha_beta (2) * PML_alpha_beta (4) * PML_alpha_beta (6) + &
                      PML_alpha_beta (3) * PML_alpha_beta (4) * PML_alpha_beta (5)


! d = beta_x * beta_y * beta_z
!------------------------------------------------
PML_alpha_beta (10)=  PML_alpha_beta (4) * PML_alpha_beta (5) * PML_alpha_beta (6)


! components of div(Lambda)
!------------------------------------------------
PML_alpha_beta (11) = d_dz_alpha_x * PML_alpha_beta (2) + d_dz_alpha_y * PML_alpha_beta (1)
PML_alpha_beta (12) = d_dy_alpha_x * PML_alpha_beta (3) + d_dy_alpha_z * PML_alpha_beta (1)
PML_alpha_beta (13) = d_dx_alpha_y * PML_alpha_beta (3) + d_dx_alpha_z * PML_alpha_beta (2)

PML_alpha_beta (14) = d_dz_alpha_x * PML_alpha_beta (5) + PML_alpha_beta (1) * d_dz_beta_y + &
                      d_dz_alpha_y * PML_alpha_beta (4) + PML_alpha_beta (2) * d_dz_beta_x

PML_alpha_beta (15) = d_dy_alpha_x * PML_alpha_beta (6) + PML_alpha_beta (1) * d_dy_beta_z + &
                      d_dy_alpha_z * PML_alpha_beta (4) + PML_alpha_beta (3) * d_dy_beta_x

PML_alpha_beta (16) = d_dx_alpha_y * PML_alpha_beta (6) + PML_alpha_beta (2) * d_dx_beta_z + &
                      d_dx_alpha_z * PML_alpha_beta (5) + PML_alpha_beta (3) * d_dx_beta_y

PML_alpha_beta (17) = d_dz_beta_x  * PML_alpha_beta (5) + d_dz_beta_y  * PML_alpha_beta (4)
PML_alpha_beta (18) = d_dy_beta_x  * PML_alpha_beta (6) + d_dy_beta_z  * PML_alpha_beta (4)
PML_alpha_beta (19) = d_dx_beta_y  * PML_alpha_beta (6) + d_dx_beta_z  * PML_alpha_beta (5)


RETURN
END 
