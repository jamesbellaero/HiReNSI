! http://www.hsl.rl.ac.uk/archive/hslarchive/packages.new/ma21/ma21.pdf
! http://www.hsl.rl.ac.uk/cgi-bin/hsl2011/download.cgi?package=ma21&vers=arch
subroutine test
    implicit double precision (a-h, o-z)
    dimension a(2,2), b(2), w(5*(5+5))
     ia = 2
     n = 2
     e = 1
     
     a(1,1) = 10000
     a(2,2) = 1
     a(2,1) = 101
     a(1,2) = 99
     
     b(1) = 1
     b(2) = 2
     
     call ma21ad(a, ia, n, b, w, e)
     write(*,*) 'b = ',b 
     write(*,*) 'success', e     

     b(1) = 1
     b(2) = 2
              
     call ma21dd(a, ia, n, b, w, e)
     write(*,*) 'b = ',b 
     write(*,*) 'success', e       

return
end

!* COPYRIGHT (c) 1971 AEA Technology
!*######DATE 21 Dec 1992
!  Toolpack tool decs employed.
!  MA21OD reference removed.
!  Label 270 removed from MA21ID.
!  Arg DIAG to MA21ND made DOUBLE PRECISION.
!  Q reference removed from MA21CD.
!  DICT reference removed from MA21HD and MA21JD.
!  DFLOAT -> DBLE.
!  Arg dimensions set to *.
!*######DATE 12 Jan 1993. Bug corrected in MA21ND - full row or column
!  not always treated correctly. Opportunity taken to introduce an
!  iteration limit as in MC19, make SMIN a constant, and make the
!  inner loop more efficient and better structured.
!*######DATE 14 Jan 1993. DABS replaced by ABS globally.
!  REAL declarations in MA21ND replaced by DOUBLE PRECISION.
!*######DATE 18 Jan 1993. ALOG replaced by LOG globally.
!*######DATE 27 Jan 1993.
!  DATA replaced by PARAMETER.
!  SAVE statements inserted.
!  EPS4, EPS, IENT and K references removed - not used.
!C  EAT 21/6/93 EXTERNAL statement put in for block data so will work on VAXs.
!C  EAT 8/7/93 EXTERNAL statement removed
!C
      BLOCK DATA MA21OD
!   Save Statement ..
!.. Common blocks ..
      COMMON /MA21ED/LP,JSCALE,EA,EB
      DOUBLE PRECISION EA,EB
      INTEGER JSCALE,LP
!..
!..
!.. Save statement ..
      SAVE /MA21ED/
!..
!.. Data statements ..
      DATA LP/6/,JSCALE/1/,EA/0.0D0/,EB/0.0D0/
!..
!.. Executable Statements ..
      END
      SUBROUTINE MA21AD(A,IA,N,B,W,E)
!.. Parameters ..
      DOUBLE PRECISION ZERO,ONE,TWO,P5
      PARAMETER (ZERO=0.D0,ONE=1.D0,TWO=2.D0,P5=0.5D0)
!..
!.. Scalar Arguments ..
      DOUBLE PRECISION E
      INTEGER IA,N
!..
!.. Array Arguments ..
      DOUBLE PRECISION A(IA,*),B(*),W(*)
!..
!.. Local Scalars ..
      DOUBLE PRECISION DET,ENORM,ENORMA,EPSN
      INTEGER I,IENT,IP,IT
!..
!.. External Subroutines ..
      EXTERNAL MA21FD,MA21ID,MA21JD
!..
!.. Common blocks ..
      COMMON /MA21ED/LP,JSCALE,EA,EB
      DOUBLE PRECISION EA,EB
      INTEGER JSCALE,LP
!..
!.. Save statement ..
      SAVE /MA21ED/
!..
!.. Executable Statements ..
      IENT = 1
      IF (N.GT.0) GO TO 250
      IF (LP.GT.0) WRITE (LP,FMT=830) N

  830 FORMAT (' ERROR RETURN FROM MA21A/AD, VALUE OF N = ',I5)

      E = -TWO
      GO TO 850

  250 CALL MA21FD(A,IA,N,W,E,IP,DET,IENT,EPSN)
      IF (IP.EQ.0) GO TO 280
      IF (IP.GT.0) GO TO 290
      I = -IP
      IF (LP.GT.0) WRITE (LP,FMT=300) I

  300 FORMAT (' MA21A/AD HAS FOUND THAT PIVOT ',I5,' IS LESS THAN A SMALL NUMBER, RESULTS MAY BE UNRELIABLE')

      E = -ONE
      GO TO 280

  290 IF (LP.GT.0) WRITE (LP,FMT=310) IP

  310 FORMAT (' MA21A/AD HAS FOUND THAT PIVOT ',I5,' IS ZERO, RESULTS MAY BE UNRELIABLE')

      E = -ONE
  280 ENORM = ZERO
      IF (E.LE.ZERO) GO TO 270
      DO 260 I = 1,N
        W(4*N+I) = B(I)
        W(N+I) = ZERO
  260 CONTINUE
      IT = 0
!C
!FORWARD SUBSTITUTION
  270 CALL MA21ID(A,IA,N,W,B,ENORM,IP,ENORMA)
      IF (IP.EQ.0) GO TO 400
      IF (IP.GT.0) GO TO 410
      I = -IP
      IF (LP.GT.0) WRITE (LP,FMT=300) I
      GO TO 840

  410 IF (LP.GT.0) WRITE (LP,FMT=310) IP
      GO TO 840

  400 IF (E.LE.ZERO) GO TO 850
      CALL MA21JD(A,IA,W,N,B,ENORM,IT,ENORMA,E,IP)
      IF (IP.EQ.0) GO TO 850
      IF (IP.GT.0) GO TO 411
      I = -IP
      IF (LP.GT.0) WRITE (LP,FMT=300) I
      GO TO 840

  411 IF (LP.GT.0) WRITE (LP,FMT=310) IP
  840 E = -ONE
  850 RETURN

      END
      SUBROUTINE MA21BD(A,IA,N,W,E)
!.. Parameters ..
      DOUBLE PRECISION ZERO,ONE,TWO,P5
      PARAMETER (ZERO=0.D0,ONE=1.D0,TWO=2.D0,P5=0.5D0)
!..
!.. Scalar Arguments ..
      DOUBLE PRECISION E
      INTEGER IA,N
!..
!.. Array Arguments ..
      DOUBLE PRECISION A(IA,*),W(*)
!..
!.. Local Scalars ..
      DOUBLE PRECISION ANORM,DET,EPSN
      INTEGER I,IENT,IP,J
!..
!.. External Subroutines ..
      EXTERNAL MA21FD,MA21GD,MA21HD
!..
!.. Common blocks ..
      COMMON /MA21ED/LP,JSCALE,EA,EB
      DOUBLE PRECISION EA,EB
      INTEGER JSCALE,LP
!..
!.. Save statement ..
      SAVE /MA21ED/
!..
!.. Executable Statements ..
      IENT = 2
      IF (N.GT.0) GO TO 250
      IF (LP.GT.0) WRITE (LP,FMT=830) N

  830 FORMAT (' ERROR RETURN FROM MA21B/BD, VALUE OF N = ',I5)

      E = -TWO
      GO TO 850

  250 CALL MA21FD(A,IA,N,W,E,IP,DET,IENT,EPSN)
      IF (IP.EQ.0) GO TO 280
      IF (IP.GT.0) GO TO 290
      I = -IP
      IF (LP.GT.0) WRITE (LP,FMT=300) I

  300 FORMAT (' MA21B/BD HAS FOUND THAT PIVOT ',I5,' IS LESS THAN A SMALL NUMBER, RESULTS MAY BE UNRELIABLE')

      E = -ONE
      GO TO 280

  290 IF (LP.GT.0) WRITE (LP,FMT=310) IP

  310 FORMAT (' ERROR RETURN FROM MA21B/BD, PIVOT ',I5,' IS ZERO')

      E = -TWO
      GO TO 850
!C
!OVERWRITE LU FACTORISATION OF A BY INVERSE OF PERMUTED A.
  280 CALL MA21GD(A,IA,N,IP,EPSN,W,ANORM)
      IF (IP.EQ.0) GO TO 730
      IF (IP.GT.0) GO TO 740
      I = -IP
      IF (LP.GT.0) WRITE (LP,FMT=300) I
      E = -ONE
      GO TO 730

  740 IF (IP.LE.N) GO TO 750
      I = (IP-1)/N
      J = IP - N*I
      IF (LP.GT.0) WRITE (LP,FMT=760) I,J

  760 FORMAT (' MA21B/BD HAS FOUND THAT INVERSE ELEMENT (',I4,',',I4,') IS LARGE, RESULTS MAY BE UNRELIABLE')

      E = -ONE
      GO TO 730

  750 IF (LP.GT.0) WRITE (LP,FMT=310) IP
      E = -TWO
      GO TO 850

  730 IF (E.LE.ZERO) GO TO 850
      CALL MA21HD(A,IA,W,N,ANORM,E)
  850 RETURN

      END
      SUBROUTINE MA21CD(A,IA,N,DET,IDET,W)
!.. Parameters ..
      DOUBLE PRECISION BASE
      INTEGER LARGE
      PARAMETER (BASE=16.0D0,LARGE=16)
      DOUBLE PRECISION ZERO,ONE,TWO,P5
      PARAMETER (ZERO=0.D0,ONE=1.D0,TWO=2.D0,P5=0.5D0)
!..
!.. Scalar Arguments ..
      DOUBLE PRECISION DET
      INTEGER IA,IDET,N
!..
!.. Array Arguments ..
      DOUBLE PRECISION A(IA,*),W(*)
!..
!.. Local Scalars ..
      DOUBLE PRECISION C32,E,EPSN,VLARGE,VSMALL
      INTEGER I,IENT,IP
!..
!.. External Functions ..
      DOUBLE PRECISION FD05AD
      EXTERNAL FD05AD
!..
!.. External Subroutines ..
      EXTERNAL MA21FD
!..
!.. Intrinsic Functions ..
      INTRINSIC ABS
!..
!.. Common blocks ..
      COMMON /MA21ED/LP,JSCALE,EA,EB
      DOUBLE PRECISION EA,EB
      INTEGER JSCALE,LP
!..
!.. Save statement ..
      SAVE /MA21ED/
!..
!.. Executable Statements ..
      VLARGE = FD05AD(5)*.10D0
      VSMALL = FD05AD(3)*10.0D0
      IENT = 3
      C32 = BASE**LARGE
      DET = ONE
      E = ZERO
      IDET = 0
      IF (N.GT.0) GO TO 250
      IF (LP.GT.0) WRITE (LP,FMT=830) N

  830 FORMAT (' ERROR RETURN FROM MA21C/CD, VALUE OF N = ',I5)

      E = -TWO
      GO TO 850

  250 CALL MA21FD(A,IA,N,W,E,IP,DET,IENT,EPSN)
      IF (IP.EQ.0) GO TO 280
      IF (IP.GT.0) GO TO 291
      I = -IP
      IF (LP.GT.0) WRITE (LP,FMT=301) I

  301 FORMAT (' MA21C/CD HAS FOUND THAT PIVOT ',I5,' IS LESS THAN A SMALL NUMBER, RESULTS MAY BE UNRELIABLE')

      E = -ONE
      GO TO 280

  291 IF (LP.GT.0) WRITE (LP,FMT=310) LP

  310 FORMAT (' MA21C/CD HAS FOUND THAT PIVOT ',I5,' IS ZERO, RESULTS MAY BE UNRELIABLE')

      E = -ONE
  280 DO 300 I = 1,N
        IF (A(I,I).NE.ZERO) GO TO 200
        DET = ZERO
        IDET = 0
        GO TO 850

  200   IF (ABS(DET).GT.ONE) GO TO 270
        IF (ABS(A(I,I)).GE.VSMALL/ABS(DET)) GO TO 290
        DET = DET*C32
        IDET = IDET - LARGE
        GO TO 200

  270   IF (ABS(A(I,I)).LE.VLARGE/ABS(DET)) GO TO 290
        DET = DET/C32
        IDET = IDET + LARGE
        GO TO 200

  290   DET = DET*A(I,I)
  300 CONTINUE
  850 RETURN

      END
      SUBROUTINE MA21DD(A,IA,N,B,W,E)
!.. Parameters ..
      DOUBLE PRECISION ZERO,ONE,TWO,P5
      PARAMETER (ZERO=0.D0,ONE=1.D0,TWO=2.D0,P5=0.5D0)
!..
!.. Scalar Arguments ..
      DOUBLE PRECISION E
      INTEGER IA,N
!..
!.. Array Arguments ..
      DOUBLE PRECISION A(IA,*),B(*),W(*)
!..
!.. Local Scalars ..
      DOUBLE PRECISION ENORM,ENORMA
      REAL WSING
      INTEGER I,II,IP,IT
!..
!.. External Subroutines ..
      EXTERNAL MA21ID,MA21JD
!..
!.. Intrinsic Functions ..
      INTRINSIC DBLE,INT
!..
!.. Common blocks ..
      COMMON /MA21ED/LP,JSCALE,EA,EB
      DOUBLE PRECISION EA,EB
      INTEGER JSCALE,LP
!..
!.. Save statement ..
      SAVE /MA21ED/
!..
!.. Executable Statements ..
!C
!C
      ENORM = ZERO
      IP = 0
      IF (N.GT.0) GO TO 251
      IF (LP.GT.0) WRITE (LP,FMT=830) N

  830 FORMAT (' ERROR RETURN FROM MA21D/DD, VALUE OF N = ',I5)

      E = -TWO
      GO TO 850
!C
!CHECK WHETHER ON A PREVIOUS ENTRY A SMALL PIVOT WAS FOUND,IF
!SO PUT ON ERROR FLAG.
  251 IP = 0
      WSING = W(N)
      II = INT(WSING+.1)
      IF (II.GT.0) GO TO 250
      IP = II
      I = -II
      II = N
      W(N) = DBLE(N)
      GO TO 250

  250 IF (E.LE.ZERO) GO TO 270
      DO 260 I = 1,N
        W(4*N+I) = B(I)
        W(N+I) = ZERO
  260 CONTINUE
      IT = 0
!C
!FORWARD  SUBSTITUTION
  270 ENORMA = ENORM
      CALL MA21ID(A,IA,N,W,B,ENORM,IP,ENORMA)
      IF (IP.EQ.0) GO TO 400
      IF (IP.GT.0) GO TO 410
      I = -IP
      IF (LP.GT.0) WRITE (LP,FMT=300) I

  300 FORMAT (' MA21D/DD HAS FOUND THAT PIVOT ',I5,' IS LESS THAN A SMALL NUMBER, RESULTS MAY BE UNRELIABLE')

      GO TO 840

  410 IF (LP.GT.0) WRITE (LP,FMT=310) LP

  310 FORMAT (' MA21D/DD HAS FOUND THAT PIVOT ',I5,' IS ZERO, RESULTS MAY BE UNRELIABLE')

      GO TO 840

  400 IF (E.LE.ZERO) GO TO 850
      CALL MA21JD(A,IA,W,N,B,ENORM,IT,ENORMA,E,IP)
      IF (IP.EQ.0) GO TO 850
      IF (IP.GT.0) GO TO 411
      I = -IP
      IF (LP.GT.0) WRITE (LP,FMT=300) I
      GO TO 840

  411 IF (LP.GT.0) WRITE (LP,FMT=310) IP
  840 E = -ONE
  850 RETURN

      END
!C@PROCESS DIRECTIVE('IBMD')
      SUBROUTINE MA21FD(A,IA,N,W,E,IP,DET,IENT,EPSN)
!.. Parameters ..
      DOUBLE PRECISION ZERO,ONE,TWO,P5
      PARAMETER (ZERO=0.D0,ONE=1.D0,TWO=2.D0,P5=0.5D0)
!..
!.. Scalar Arguments ..
      DOUBLE PRECISION DET,E,EPSN
      INTEGER IA,IENT,IP,N
!..
!.. Array Arguments ..
      DOUBLE PRECISION A(IA,*),W(N,*)
!..
!.. Local Scalars ..
      DOUBLE PRECISION AA,AB,AC,AM,EPS,EPS4,SCPROD,WMAX,WW
      INTEGER I,II,IPROD,IT,J,L
!..
!.. External Functions ..
      DOUBLE PRECISION FD05AD
      EXTERNAL FD05AD
!..
!.. External Subroutines ..
      EXTERNAL MA21ND
!..
!.. Intrinsic Functions ..
      INTRINSIC ABS,DBLE,DMAX1
!..
!.. Common blocks ..
      COMMON /MA21ED/LP,JSCALE,EA,EB
      DOUBLE PRECISION EA,EB
      INTEGER JSCALE,LP
!..
!.. Save statement ..
      SAVE /MA21ED/
!..
!.. Executable Statements ..
      EPS = FD05AD(1)*10.0D0
      EPS4 = 4.0D0*EPS
      DET = ONE
      IP = 0
      IF (N.LT.1) GO TO 240
      IP = 0
      IF (JSCALE.GT.0) THEN
!   FIND SCALING FACTORS.
        CALL MA21ND(A,N,IA,W(1,5),W(1,1),IT)
        DO 110 I = 1,N
          W(I,3) = W(I,1)
          W(I,4) = W(I,2)
  110   CONTINUE
      END IF
!STORE IN W(J,1),J=1,N THE MAXIMAL ELEMENTS  OF COLUMNS AFTER
!APPLICATION OF SCALING.
      EPSN = ZERO
      IF (JSCALE.GT.0) THEN
        DO 140 J = 1,N
          DO 125 I = 1,N
            W(I,2) = W(I,3)*A(I,J)
  125     CONTINUE
          WMAX = -ONE
          DO 126 I = 1,N
            WW = ABS(W(I,2))
            IF (WMAX.LT.WW) THEN
              WMAX = WW
              II = I
            END IF

  126     CONTINUE
          W(J,1) = ABS(W(II,2))
          EPSN = DMAX1(EPSN,W(J,1)*W(J,4))
  140   CONTINUE

      ELSE
        DO 141 J = 1,N
          WMAX = -ONE
          DO 143 I = 1,N
            WW = ABS(A(I,J))
            IF (WMAX.LT.WW) THEN
              WMAX = WW
              II = I
            END IF

  143     CONTINUE
          W(J,1) = ABS(A(II,J))
          EPSN = DMAX1(EPSN,W(J,1))
  141   CONTINUE
      END IF

      EPSN = EPSN*EPS4
      IF ((E.LE.ZERO) .OR. (IENT.EQ.3)) GO TO 160
!MAKE COPY OF MATRIX
      DO 150 I = 1,N
        DO 150 J = 1,N
  150 W(I,J+5) = A(I,J)
!C
!FACTORISATION OF MATRIX INTO L*U, WHERE L IS A LOWER UNIT
!TRIANGLE AND U IS UPPER TRIANGLE
  160 DO 230 L = 1,N
        AM = ZERO
        II = L
!EVALUATE  ELEMENTS IN PIVOTAL COLUMN.
!CIBMD IGNORE RECRDEPS
        DO 170 I = L,N
          SCPROD = ZERO
!CIBMD PREFER SCALAR
          DO 165 IPROD = 1,L - 1
            SCPROD = SCPROD + A(IPROD,L)*A(I,IPROD)
  165     CONTINUE
          A(I,L) = A(I,L) - SCPROD
  170   CONTINUE
!LOOK FOR MAXIMUM ELEMENT IN PIVOTAL COLUMN.
        IF (JSCALE.GT.0) THEN
          DO 175 I = L,N
            AB = ABS(A(I,L))
            AB = AB*W(I,3)
            IF (AM.GE.AB) GO TO 175
            AM = AB
            II = I
  175     CONTINUE

        ELSE
          WMAX = -ONE
          DO 176 I = L,N
            WW = ABS(A(I,L))
            IF (WMAX.LT.WW) THEN
              WMAX = WW
              II = I
            END IF

  176     CONTINUE
          AM = ABS(A(II,L))
        END IF
!C
!TEST FOR SMALL OR ZERO PIVOT.
        AB = W(L,1)
        IF (AM.LE.EPS4*AB) IP = -L
        IF (AM.NE.ZERO) GO TO 190
        IP = L
!C
  190   IF (II.NE.L) THEN
!C
!   INTERCHANGE ROWS N AND II
!   INTERCHANGE EQUILIBRATION FACTORS
!   W(L,1)=II MEANS INTERCHANGE BETWEEN ROWS L AND II
          IF (IENT.EQ.3) DET = -DET
          IF (JSCALE.GT.0) THEN
            AA = W(L,3)
            W(L,3) = W(II,3)
            W(II,3) = AA
          END IF

          DO 210 I = 1,N
            AA = A(L,I)
            A(L,I) = A(II,I)
            A(II,I) = AA
  210     CONTINUE
        END IF

        W(L,1) = DBLE(II)
        IF (L.EQ.N) GO TO 240
        AA = ONE
        AC = A(L,L)
        IF (ABS(AC).NE.ZERO) AA = AA/AC
!UPPER TRIANGLE
        DO 222 I = L + 1,N
          A(I,L) = AA*A(I,L)
  222   CONTINUE
!CIBMD IGNORE RECRDEPS
        DO 227 I = L + 1,N
          SCPROD = ZERO
!CIBMD PREFER SCALAR
          DO 225 IPROD = 1,L - 1
            SCPROD = SCPROD + A(IPROD,I)*A(L,IPROD)
  225     CONTINUE
          A(L,I) = A(L,I) - SCPROD
  227   CONTINUE
  230 CONTINUE
!C
  240 RETURN

      END
      SUBROUTINE MA21GD(A,IA,N,IP,EPSN,W,ANORM)
!.. Parameters ..
      DOUBLE PRECISION ZERO,ONE,TWO,P5
      PARAMETER (ZERO=0.D0,ONE=1.D0,TWO=2.D0,P5=0.5D0)
!..
!.. Scalar Arguments ..
      DOUBLE PRECISION ANORM,EPSN
      INTEGER IA,IP,N
!..
!.. Array Arguments ..
      DOUBLE PRECISION A(IA,N),W(N,*)
!..
!.. Local Scalars ..
      DOUBLE PRECISION AA,AB,AC,ENORM,ERR,SCPROD
      REAL WSING
      INTEGER I,II,IPROD,J,K
!..
!.. Intrinsic Functions ..
      INTRINSIC ABS,DMAX1,INT
!..
!.. Common blocks ..
      COMMON /MA21ED/LP,JSCALE,EA,EB
      DOUBLE PRECISION EA,EB
      INTEGER JSCALE,LP
!..
!.. Save statement ..
      SAVE /MA21ED/
!..
!.. Executable Statements ..
!C
!OVERWRITE LU FACTORISATION OF A BY INVERSE OF PERMUTED A.
      IF (N.LT.2) GO TO 520
      DO 515 I = 2,N
        K = I - 1
        DO 510 J = 1,K
          SCPROD = ZERO
          DO 505 IPROD = J + 1,I - 1
            SCPROD = SCPROD + A(I,IPROD)*A(IPROD,J)
  505     CONTINUE
          A(I,J) = -A(I,J) - SCPROD
  510   CONTINUE
  515 CONTINUE
  520 DO 540 K = 1,N
        I = N + 1 - K
        ERR = ONE/EPSN
        IF (JSCALE.GT.0) ERR = ERR*W(I,4)
        DO 530 J = I,N
          W(J,2) = A(I,J)
  530   A(I,J) = ZERO
        A(I,I) = ONE
        AA = ONE/W(I,2)
        DO 540 J = 1,N
          AB = ONE
          IF (JSCALE.GT.0) AB = W(J,3)
          AC = A(I,J)
          IF (I.NE.N) THEN
            SCPROD = ZERO
            DO 535 IPROD = I + 1,N
              SCPROD = SCPROD + W(IPROD,2)*A(IPROD,J)
  535       CONTINUE
            AC = A(I,J) - SCPROD
          END IF

          A(I,J) = AC*AA
          IF (ABS(A(I,J)).GE.ERR*AB) IP = N*I + J
  540 CONTINUE
      ANORM = ZERO
      DO 590 K = 1,N
        I = N + 1 - K
        WSING = W(I,1)
        II = INT(WSING+.1)
        IF (II.EQ.I) GO TO 570
        IF (JSCALE.LE.0) GO TO 550
        AA = W(I,3)
        W(I,3) = W(II,3)
        W(II,3) = AA
  550   DO 560 J = 1,N
          AA = A(J,II)
          A(J,II) = A(J,I)
  560   A(J,I) = AA
        W(I,1) = W(II,1)
  570   ENORM = ZERO
        DO 580 J = 1,N
          AB = ABS(A(J,II))
          IF (JSCALE.GT.0) AB = AB/W(J,4)
  580   ENORM = DMAX1(ENORM,AB)
        W(II,1) = ENORM
        AB = ONE
        IF (JSCALE.GT.0) AB = W(II,3)
  590 ANORM = DMAX1(ANORM,ENORM/AB)
      RETURN

      END
      SUBROUTINE MA21HD(A,IA,W,N,ANORM,E)
!.. Parameters ..
      DOUBLE PRECISION ZERO,ONE,TWO,P5
      PARAMETER (ZERO=0.D0,ONE=1.D0,TWO=2.D0,P5=0.5D0)
!..
!.. Scalar Arguments ..
      DOUBLE PRECISION ANORM,E
      INTEGER IA,N
!..
!.. Array Arguments ..
      DOUBLE PRECISION A(IA,N),W(N,*)
!..
!.. Local Scalars ..
      DOUBLE PRECISION AA,AB,AC,AM,ANORMA,AXNORM,ENORM,EPS,ERR,SCPROD
      INTEGER I,IPROD,IRANDL,IRANDR,J,K,L
!..
!.. External Functions ..
      DOUBLE PRECISION FA01AD,FD05AD
      EXTERNAL FA01AD,FD05AD
!..
!.. External Subroutines ..
      EXTERNAL FA01CD,FA01DD
!..
!.. Intrinsic Functions ..
      INTRINSIC ABS,DMAX1
!..
!.. Common blocks ..
      COMMON /MA21ED/LP,JSCALE,EA,EB
      DOUBLE PRECISION EA,EB
      INTEGER JSCALE,LP
!..
!.. Save statement ..
      SAVE /MA21ED/
!..
!.. Executable Statements ..
      EPS = FD05AD(1)*10.0D0
      CALL FA01CD(IRANDL,IRANDR)
      CALL FA01DD(21845,21845)
      AXNORM = ANORM
  600 ANORMA = ANORM
      ANORM = ZERO
      L = 0
  610 L = L + 1
!C
!INVERSE OF A IS ITERATIVELY REFINED BY COLUMNS
!MAKE COPY OF APPROPIATE COLUMN.
      DO 620 I = 1,N
  620 W(I,2) = A(I,L)
!C
!COMPUTE RESIDUAL
      DO 450 J = 1,N
        AA = ZERO
        IF (L.EQ.J) AA = ONE
        SCPROD = ZERO
        DO 625 IPROD = 1,N
          SCPROD = SCPROD + W(J,IPROD+5)*W(IPROD,2)
  625   CONTINUE
        AC = AA - SCPROD
        AA = ZERO
        IF (EA) 400,440,380
!C
!MAKE PSEUDO RANDOM CHANGES TO ELEMENTS OF A AND B
  380   DO 390 K = 1,N
  390   AA = AA + FA01AD(-K)*W(J,K+5)*W(K,2)
        GO TO 420

  400   DO 410 K = 1,N
  410   AA = AA + FA01AD(-K)*W(K,2)
  420   AC = AC - ABS(EA)*AA
  440   W(J,5) = AC
  450 CONTINUE
!C
!CALCULATE THE CHANGE IN APPROPIATE COLUMN OF INVERSE OF A.
      ENORM = ZERO
      DO 640 I = 1,N
        W(I,2) = ZERO
        SCPROD = ZERO
        DO 455 IPROD = 1,N
          SCPROD = SCPROD + A(I,IPROD)*W(IPROD,5)
  455   CONTINUE
        W(I,2) = W(I,2) + SCPROD
        AM = ONE
        IF (JSCALE.GT.0) AM = W(I,4)
  640 ENORM = DMAX1(ENORM,ABS(W(I,2))/AM)
      AB = W(L,1)
      IF (ENORM.GT.P5*AB) GO TO 660
!C
!UPDATE APPROPIATE COLUMN OF INVERSE OF A.
      DO 650 J = 1,N
  650 A(J,L) = A(J,L) + W(J,2)
      W(L,1) = ENORM
      AB = ENORM
  660 ERR = ONE
      IF (JSCALE.GT.0) ERR = W(L,3)
      ANORM = DMAX1(ANORM,AB/ERR)
      IF (L.LT.N) GO TO 610
      IF (ANORM.GT.P5*ANORMA) GO TO 670
      IF (ANORM-EPS*AXNORM) 680,680,600
  670 ANORM = ANORMA
  680 IF (JSCALE.LE.0) GO TO 700
!C
!SET UP ACCURACY ESTIMATES FOR MATRIX INVERSE.
      ENORM = ZERO
      ERR = ZERO
      DO 690 I = 1,N
        AB = ANORM*W(I,3)
        ENORM = DMAX1(ENORM,AB)
        W(I,1) = AB
        AB = W(I,4)
        ERR = DMAX1(ERR,AB)
  690 W(I,2) = AB
      E = ENORM*ERR
      GO TO 710

  700 E = ANORM
  710 CALL FA01DD(IRANDL,IRANDR)
      RETURN

      END
      SUBROUTINE MA21ID(A,IA,N,W,B,ENORM,IP,ENORMA)
!.. Parameters ..
      DOUBLE PRECISION ZERO,ONE,TWO,P5
      PARAMETER (ZERO=0.D0,ONE=1.D0,TWO=2.D0,P5=0.5D0)
!..
!.. Scalar Arguments ..
      DOUBLE PRECISION ENORM,ENORMA
      INTEGER IA,IP,N
!..
!.. Array Arguments ..
      DOUBLE PRECISION A(IA,N),B(N),W(N,*)
!..
!.. Local Scalars ..
      DOUBLE PRECISION AA,AC,AM,SCPROD
      REAL WSING
      INTEGER I,II,IPROD,K
!..
!.. Intrinsic Functions ..
      INTRINSIC ABS,DMAX1,INT
!..
!.. Common blocks ..
      COMMON /MA21ED/LP,JSCALE,EA,EB
      DOUBLE PRECISION EA,EB
      INTEGER JSCALE,LP
!..
!.. Save statement ..
      SAVE /MA21ED/
!..
!.. Executable Statements ..
!C
!FORWARD SUBSTITUTION
      DO 280 I = 1,N
        WSING = W(I,1)
        II = INT(WSING+.1)
        AA = B(II)
        B(II) = B(I)
        SCPROD = ZERO
        DO 275 IPROD = 1,I - 1
          SCPROD = SCPROD + A(I,IPROD)*B(IPROD)
  275   CONTINUE
        B(I) = AA - SCPROD
  280 CONTINUE
      ENORMA = ENORM
      ENORM = ZERO
!C
!CALCULATE NORMS OF X,CHANGE IN SOLUTION
!BACKWARD SUBSTITUTION
      DO 310 K = 1,N
        I = N + 1 - K
        AA = A(I,I)
        IF (ABS(AA).EQ.ZERO) GO TO 290
        AC = B(I)
        IF (I.NE.N) THEN
          SCPROD = ZERO
          DO 285 IPROD = 1,N - I
            SCPROD = SCPROD + A(I,I+IPROD)*B(I+IPROD)
  285     CONTINUE
          AC = B(I) - SCPROD
        END IF

        B(I) = AC/AA
        GO TO 300

  290   IP = I
        B(I) = ZERO
  300   AM = ONE
        IF (JSCALE.GT.0) AM = W(I,4)
  310 ENORM = DMAX1(ENORM,ABS(B(I))/AM)
      RETURN

      END
      SUBROUTINE MA21JD(A,IA,W,N,B,ENORM,IT,ENORMA,E,IP)
!.. Parameters ..
      DOUBLE PRECISION ZERO,ONE,TWO,P5
      PARAMETER (ZERO=0.D0,ONE=1.D0,TWO=2.D0,P5=0.5D0)
!..
!.. Scalar Arguments ..
      DOUBLE PRECISION E,ENORM,ENORMA
      INTEGER IA,IP,IT,N
!..
!.. Array Arguments ..
      DOUBLE PRECISION A(IA,N),B(N),W(N,*)
!..
!.. Local Scalars ..
      DOUBLE PRECISION AA,AB,AC,EPS,ERR,SCPROD,XNORM
      INTEGER I,IPROD,IRANDL,IRANDR,J,K
!..
!.. External Functions ..
      DOUBLE PRECISION FA01AD,FD05AD
      EXTERNAL FA01AD,FD05AD
!..
!.. External Subroutines ..
      EXTERNAL FA01CD,FA01DD,MA21ID
!..
!.. Intrinsic Functions ..
      INTRINSIC ABS,DMAX1
!..
!.. Common blocks ..
      COMMON /MA21ED/LP,JSCALE,EA,EB
      DOUBLE PRECISION EA,EB
      INTEGER JSCALE,LP
!..
!.. Save statement ..
      SAVE /MA21ED/
!..
!.. Executable Statements ..
      EPS = FD05AD(1)*10.0D0
      XNORM = ENORM
      CALL FA01CD(IRANDL,IRANDR)
      CALL FA01DD(21845,21845)
!C
!UPDATE SOLUTION VECTOR X
  330 DO 340 I = 1,N
  340 W(I,2) = W(I,2) + B(I)
      IT = IT + 1
!C
!COMPUTE RESIDUAL
      DO 450 J = 1,N
        AA = W(J,5)
        SCPROD = ZERO
        DO 345 IPROD = 1,N
          SCPROD = SCPROD + W(J,IPROD+5)*W(IPROD,2)
  345   CONTINUE
        AC = AA - SCPROD
        AA = ZERO
        IF (EA) 400,430,380
!C
!MAKE PSEUDO RANDOM CHANGES TO ELEMENTS OF A AND B
  380   DO 390 K = 1,N
  390   AA = AA + FA01AD(-K)*W(J,K+5)*W(K,2)
        GO TO 420

  400   DO 410 K = 1,N
  410   AA = AA + FA01AD(-K)*W(K,2)
  420   AC = AC - ABS(EA)*AA
  430   AA = ABS(EB)*FA01AD(-J)
        IF (EB.GE.ZERO) AA = AA*W(J,5)
        AC = AA + AC
        B(J) = AC
  450 CONTINUE
      CALL MA21ID(A,IA,N,W,B,ENORM,IP,ENORMA)
      IF (IP.NE.0) GO TO 850
      IF (ENORM.GT.P5*ENORMA) GO TO 460
      IF (ENORM-EPS*XNORM) 470,470,330
  460 ENORM = ENORMA
  470 DO 480 I = 1,N
  480 B(I) = W(I,2)
      E = ENORM
      IF (JSCALE.LE.0) GO TO 710
!C
!SET UP ACCURACY ESTIMATES FOR SOLUTION VECTOR.
      ERR = ZERO
      DO 490 I = 1,N
        W(I,2) = ENORM*W(I,4)
        AB = W(I,2)
  490 ERR = DMAX1(ERR,AB)
      E = ERR
  710 CALL FA01DD(IRANDL,IRANDR)
!C
!C
  850 RETURN

      END
      SUBROUTINE MA21ND(A,N,NN,DIAG,RES,IS)
!.. Parameters ..
      DOUBLE PRECISION ZERO,ONE,P5
      INTEGER MAXIT
      DOUBLE PRECISION SMIN
      PARAMETER (ZERO=0.0D0,ONE=1.0D0,P5=0.5D0,MAXIT=100,SMIN=.01)
!..
!.. Scalar Arguments ..
      INTEGER IS,N,NN
!..
!.. Array Arguments ..
      DOUBLE PRECISION A(NN,*),DIAG(N),RES(N,4)
!..
!.. Local Scalars ..
      DOUBLE PRECISION E,E1,EM,Q,Q1,QM,R1,RSUM,S,S1,SM,SSUM,SUM,U,V
      REAL XLG
      INTEGER I,IDIAG,ITER,J,L,LL,LM,M,NP1
!..
!.. Intrinsic Functions ..
      INTRINSIC ABS,FLOAT,LOG,SIGN
!..
!.. Executable Statements ..
!RES IS USED TO RETURN SCALING FACTORS AS INTEGRAL
!     POWERS OF BASE, AND AS WORKSPACE
!IS IS SET TO 0 ON SUCCESSFUL COMPLETION, TO I IF ROW I HAS ONLY
!   ZERO ELEMENTS, TO -I IF COLUMN I HAS ONLY ZERO ELEMENTS
!DIAG IS USED TO HOLD COUNTS OF NON-ZEROS IN ROWS AND COLUMNS
! AND TO RETURN SCALING POWERS
!C
!SMIN IS USED IN A CONVERGENCE TEST ON (RESIDUAL NORM)**2
!SET UP CONSTANTS
      NP1 = N + 1
      XLG = 1./LOG(16.)
      IS = 0
!INITIALISE FOR ACCUMULATION OF SUMS AND PRODUCTS
      DO 2 L = 1,2
        DO 1 I = 1,N
          RES(I,L) = ZERO
          RES(I,L+2) = ZERO
          DIAG(I) = ZERO
    1   CONTINUE
    2 CONTINUE
      DO 4 I = 1,N
        DO 3 J = 1,N
          U = ABS(A(I,J))
          IF (U.EQ.ZERO) GO TO 3
          U = LOG(U)*XLG + 1.
!         COUNT NON-ZEROS IN ROW AND COLUMN
          DIAG(I) = DIAG(I) + ONE
          DIAG(J) = DIAG(J) + NP1
          RES(I,1) = RES(I,1) + U
          RES(J,3) = RES(J,3) + U
    3   CONTINUE
    4 CONTINUE
!COMPUTE RHS VECTORS TESTING FOR ZERO ROW OR COLUMN
      SSUM = ZERO
      J = 0
      DO 8 I = 1,N
        IDIAG = DIAG(I)
        LL = IDIAG/NP1
        LM = IDIAG - NP1*LL
        J = J + LM
        IF (LM.GT.0) GO TO 153
        LM = 1
        IS = I
  153   RES(I,1) = RES(I,1)/FLOAT(LM)
        IF (LL.GT.0) GO TO 154
        LL = 1
        IS = -I
  154   RES(I,3) = RES(I,3)/FLOAT(LL)
        DIAG(I) = LL*NP1 + LM
        SSUM = SSUM + RES(I,3)
    8 CONTINUE
      SM = SMIN*J
!SWEEP TO COMPUTE INITIAL RESIDUAL VECTOR
      RSUM = ZERO
      DO 110 I = 1,N
        SUM = SSUM
        IDIAG = DIAG(I)
        LM = IDIAG - NP1* (IDIAG/NP1)
        IF (LM.GE.N) GO TO 109
        SUM = ZERO
        DO 10 J = 1,N
          IF (A(I,J).EQ.ZERO) GO TO 10
          SUM = SUM + RES(J,3)
   10   CONTINUE
  109   RES(I,1) = RES(I,1) - SUM/FLOAT(LM)
        RSUM = RSUM + RES(I,1)
  110 CONTINUE
!INITIALISE ITERATION
      E = ZERO
      E1 = ZERO
      Q = 1.
      S = ZERO
      DO 11 I = 1,N
        IDIAG = DIAG(I)
        LM = IDIAG - (IDIAG/NP1)*NP1
        S = S + FLOAT(LM)*RES(I,1)**2
   11 CONTINUE
      L = 2
!ITERATION LOOP
      DO 20 ITER = 1,MAXIT
        IF (S.LE.SM) GO TO 100
        EM = E*E1
!C    SWEEP THROUGH MATRIX TO UPDATE RESIDUAL VECTOR
        DO 220 I = 1,N
          IDIAG = DIAG(I)
          LM = IDIAG/NP1
          SUM = RSUM
          IF (L.EQ.1) THEN
            LM = IDIAG - LM*NP1
            IF (LM.LT.N) THEN
              SUM = ZERO
              DO 21 J = 1,N
                IF (A(I,J).NE.ZERO) SUM = SUM + RES(J,2)
   21         CONTINUE
            END IF

          ELSE
            IF (LM.LT.N) THEN
              SUM = ZERO
              DO 22 J = 1,N
                IF (A(J,I).NE.ZERO) SUM = SUM + RES(J,1)
   22         CONTINUE
            END IF

          END IF

          RES(I,L) = RES(I,L) + SUM
  220   CONTINUE
        S1 = S
        S = ZERO
        RSUM = ZERO
        DO 23 I = 1,N
          V = -RES(I,L)/Q
          IDIAG = DIAG(I)
          LL = IDIAG/NP1
          IF (L.EQ.1) LL = IDIAG - LL*NP1
          RES(I,L) = V/FLOAT(LL)
          RSUM = RSUM + RES(I,L)
          S = S + V*RES(I,L)
   23   CONTINUE
        E1 = E
        E = Q*S/S1
        Q1 = Q
        Q = 1. - E
        M = 3 - L
        IF (S.GT.SM) GO TO 27
        E = M - 1
        M = 1
        Q = 1.
   27   IF (L.EQ.2) GO TO 25
        QM = Q*Q1
        DO 24 I = 1,N
          RES(I,4) = (EM*RES(I,4)+RES(I,2))/QM
   24   RES(I,3) = RES(I,3) + RES(I,4)
   25   L = M
        DO 26 I = 1,N
          IDIAG = DIAG(I)
          LL = IDIAG/NP1
          IF (L.EQ.1) LL = IDIAG - LL*NP1
          RES(I,L) = RES(I,L)*FLOAT(LL)*E
   26   CONTINUE
   20 CONTINUE
! SWEEP THROUGH MATRIX TO GET ROW SCALING POWERS
  100 DO 103 I = 1,N
        DO 103 J = 1,N
          U = ABS(A(I,J))
          IF (U.EQ.ZERO) GO TO 103
          U = LOG(U)*XLG - 1.
          RES(I,1) = RES(I,1) + RES(J,3) - U
  103 CONTINUE
! CONVERT POWERS TO INTEGERS
      DO 104 I = 1,N
        IDIAG = DIAG(I)
        LL = IDIAG - (IDIAG/NP1)*NP1
        V = RES(I,1)/FLOAT(LL)
!DIAG(I,1)=V+SIGN(P5,V)
        J = V + SIGN(P5,V)
        RES(I,1) = 16.**J
!DIAG(I,2)=-(RES(I,3)+SIGN(P5,RES(I,3)))
        R1 = RES(I,3)
        J = - (RES(I,3)+SIGN(P5,R1))
        RES(I,2) = 16.**J
  104 CONTINUE
      RETURN

      END





! COPYRIGHT (c) 1967 AEA Technology
!######DATE 4 Oct 1992
!       Toolpack tool decs employed.
!       SAVE statement for COMMON FA01ED added.
!  EAT 21/6/93 EXTERNAL statement put in for block data on VAXs.
!
!
      DOUBLE PRECISION FUNCTION FA01AD(I)
!     .. Scalar Arguments ..
      INTEGER I
!     ..
!     .. Local Scalars ..
      DOUBLE PRECISION R,S
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC DINT,MOD
!     ..
!     .. Common blocks ..
      COMMON /FA01ED/GL,GR
      DOUBLE PRECISION GL,GR
!     ..
!     .. Data block external statement
      EXTERNAL FA01FD
!     ..
!     .. Save statement ..
      SAVE /FA01ED/
!     ..
!     .. Executable Statements ..
      R = GR*9228907D0/65536D0
      S = DINT(R)
      GL = MOD(S+GL*9228907D0,65536D0)
      GR = R - S
      IF (I.GE.0) FA01AD = (GL+GR)/65536D0
      IF (I.LT.0) FA01AD = (GL+GR)/32768D0 - 1.D0
      GR = GR*65536D0
      RETURN

      END
      SUBROUTINE FA01BD(MAX,NRAND)
!     .. Scalar Arguments ..
      INTEGER MAX,NRAND
!     ..
!     .. External Functions ..
      DOUBLE PRECISION FA01AD
      EXTERNAL FA01AD
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC DBLE,INT
!     ..
!     .. Executable Statements ..
      NRAND = INT(FA01AD(1)*DBLE(MAX)) + 1
      RETURN

      END
      SUBROUTINE FA01CD(IL,IR)
!     .. Scalar Arguments ..
      INTEGER IL,IR
!     ..
!     .. Common blocks ..
      COMMON /FA01ED/GL,GR
      DOUBLE PRECISION GL,GR
!     ..
!     .. Save statement ..
      SAVE /FA01ED/
!     ..
!     .. Executable Statements ..
      IL = GL
      IR = GR
      RETURN

      END
      SUBROUTINE FA01DD(IL,IR)
!     .. Scalar Arguments ..
      INTEGER IL,IR
!     ..
!     .. Common blocks ..
      COMMON /FA01ED/GL,GR
      DOUBLE PRECISION GL,GR
!     ..
!     .. Save statement ..
      SAVE /FA01ED/
!     ..
!     .. Executable Statements ..
      GL = IL
      GR = IR
      RETURN

      END
      BLOCK DATA FA01FD
!     .. Common blocks ..
      COMMON /FA01ED/GL,GR
      DOUBLE PRECISION GL,GR
!     ..
!     .. Save statement ..
      SAVE /FA01ED/
!     ..
!     .. Data statements ..
      DATA GL/21845D0/
      DATA GR/21845D0/
!     ..
!     .. Executable Statements ..
      END


      DOUBLE PRECISION FUNCTION DDOT(N,DX,INCX,DY,INCY)
!
!     FORMS THE DOT PRODUCT OF TWO VECTORS.
!     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE.
!     JACK DONGARRA, LINPACK, 3/11/78.
!
      DOUBLE PRECISION DX(*),DY(*),DTEMP
      INTEGER I,INCX,INCY,IX,IY,M,MP1,N
!
      DDOT = 0.0D0
      DTEMP = 0.0D0
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20
!
!        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
!          NOT EQUAL TO 1
!
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DTEMP = DTEMP + DX(IX)*DY(IY)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      DDOT = DTEMP
      RETURN
!
!        CODE FOR BOTH INCREMENTS EQUAL TO 1
!
!
!        CLEAN-UP LOOP
!
   20 M = MOD(N,5)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        DTEMP = DTEMP + DX(I)*DY(I)
   30 CONTINUE
      IF( N .LT. 5 ) GO TO 60
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
        DTEMP = DTEMP + DX(I)*DY(I) + DX(I + 1)*DY(I + 1) +   DX(I + 2)*DY(I + 2) + DX(I + 3)*DY(I + 3) + DX(I + 4)*DY(I + 4)
   50 CONTINUE
   60 DDOT = DTEMP
      RETURN
      END


! COPYRIGHT (c) 1988 AEA Technology
!######DATE 21 Jan 1993
!       Toolpack tool decs employed.
!       SAVE statement added.
! 1/10/98 DC(3) not initialized to avoid SUN f90 failure
! 16 October 2001: STOP and WRITE statements removed.

      DOUBLE PRECISION FUNCTION FD05AD(INUM)
!----------------------------------------------------------------
!  Real constants for: IEEE double precision (8-byte arithmetic)
!
!  Obtained from H.S.L. subroutine ZE02AM.
!  Nick Gould and Sid Marlow, Harwell Laboratory, April 1988.
!----------------------------------------------------------------
!     .. Scalar Arguments ..
      INTEGER INUM
!     ..
!     .. Local Arrays ..
      DOUBLE PRECISION DC(5)
!     ..
!     .. Save statement ..
      SAVE DC
!     ..
!     .. Data statements ..
!
!  DC(1) THE SMALLEST POSITIVE NUMBER: 1.0 + DC(1) > 1.0.
!  DC(2) THE SMALLEST POSITIVE NUMBER: 1.0 - DC(2) < 1.0.
!  DC(3) THE SMALLEST NONZERO +VE REAL NUMBER.
!  DC(4) THE SMALLEST FULL PRECISION +VE REAL NUMBER.
!  DC(5) THE LARGEST FINITE +VE REAL NUMBER.
!
      DATA DC(1)/2.2204460492504D-16/
      DATA DC(2)/1.1102230246253D-16/
!     DATA DC(3)/4.9406564584126D-324/
      DATA DC(4)/2.2250738585073D-308/
      DATA DC(5)/1.7976931348622D+308/
!     ..
!     .. Executable Statements ..

      IF ( INUM .LE. 0 ) THEN
         FD05AD = DC( 1 )
      ELSE IF ( INUM .GE. 6 ) THEN
         FD05AD = DC( 5 )
      ELSE IF ( INUM .EQ. 3 ) THEN
         FD05AD = DC(4)/2.0D0**52
      ELSE
         FD05AD = DC( INUM )
      ENDIF
      RETURN
      END
