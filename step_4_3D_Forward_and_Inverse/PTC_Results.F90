
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                                                                                                                  ++
! Last Update:  06 Jan 2012                                                                                                                        ++
!                                                                                                                                                  ++
! Description: THIS Module prints out all the results produced by the codes                                                                        ++
!                                                                                                                                                  ++
!              Subroutines                     Description                                                                                         ++
!              ==============                  ==============                                                                                      ++
!              ResultsHis                                                                                                                          ++
!                                                                                                                                                  ++
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Module Results ;

Use Parameters ;

!Implicit None ;


  Interface ResDyn ;
    Module Procedure ResultsHeader, ResultsHis ;
  End Interface ;

Contains ;

! ***************************************************************************************************************************************************
!                                                                                                                                                  **
! Last Update:  06 Jan 2012                                                                                                                        **
! Description: THIS Subroutine writes down the header of all output files                                                                          **
!                                                                                                                                                  **
! ***************************************************************************************************************************************************

Subroutine ResultsHeader (                                                                                                      &
NNDH, NNVH, NNAH,                                                                                                               & ! Integer Variables
!                                                                                                                               & ! Real Variables
NDAN, NDPN, NVAN, NVPN, NAAN, NAPN                                                                                              & ! Integer Arrays
!                                                                                                                               & ! Real Arrays
!                                                                                                                               & ! Characters
!                                                                                                                               & ! Type 
) ;

!Implicit None ;

! =========================== Global Variables ======================================================================================================

! - Integer Variables -------------------------------------------------------------------------------------------------------------------------------
Integer, Intent(In)    :: NNDH, NNVH, NNAH ;

! - Real Variables ----------------------------------------------------------------------------------------------------------------------------------
!#Real (Kind=Dbl), Intent(In)    ::  ;
!#Real (Kind=Dbl), Intent(InOut) ::  ;
!#Real (Kind=Dbl), Intent(OUT)   ::  ;
! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
!#Complex, Intent(In)    ::  ;
!#Complex, Intent(InOut) ::  ;
!#Complex, Intent(OUT)   ::  ;
! - Integer Arrays ----------------------------------------------------------------------------------------------------------------------------------
Integer, Intent(In), Dimension (:  )  :: NDAN, NDPN, NVAN, NVPN, NAAN, NAPN ;

! - Real Arrays -------------------------------------------------------------------------------------------------------------------------------------
!#Real (Kind=Dbl), Intent(In), Dimension (:  )  :: 
!#Real (Kind=Dbl), Intent(InOut), Dimension (:  )  :: ;
!#Real (Kind=Dbl), Intent(OUT), Dimension (:  )  ::  ;
! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
!#Complex, Intent(In)    ::  ;
!#Complex, Intent(InOut) ::  ;
!#Complex, Intent(OUT)   ::  ;
! - Character Variables -----------------------------------------------------------------------------------------------------------------------------
!#Character   ::  ;
! - Logical Variables -------------------------------------------------------------------------------------------------------------------------------
!#Logical   ::  ;
! - Type DECLERATIONS -------------------------------------------------------------------------------------------------------------------------------
!#Type() :: ;
! =========================== Local Variables =======================================================================================================
! - Integer Variables -------------------------------------------------------------------------------------------------------------------------------
Integer  :: IO_Write ;               ! Holds error of write statements
Integer  :: UnFile ;                ! Holds unit of a file for error message

Integer   :: I ;                      ! Loop indeces.

! - Real Variables ----------------------------------------------------------------------------------------------------------------------------------
!#Real (Kind=Dbl)  ::  ;
! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
!#Complex  ::  ;
! - Integer Arrays ----------------------------------------------------------------------------------------------------------------------------------
!#Integer (Kind=Shrt)  ::  ;
! - Real Arrays -------------------------------------------------------------------------------------------------------------------------------------
!#Real (Kind=Dbl)  ::  ;
! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
!#Complex  ::  ;
! - Character Variables -----------------------------------------------------------------------------------------------------------------------------
!#Character   ::  ;
! - Logical Variables -------------------------------------------------------------------------------------------------------------------------------
!#Logical   ::  ;

! =========================== Subroutine CODE =======================================================================================================

UnFile = UN_HisD ;
Write (Unit = UnFile, FMT = "('History of Displacements and stresses ')"   , ADVANCE = 'YES', ASYNCHRONOUS = 'NO'                , IOSTAT = IO_Write, ERR = 1006 ) ;
Write (Unit = UnFile, FMT = "('Time',5x,<NNDH>('App#',I19,1x,'PETSc#',I19,1x))" , ADVANCE = 'YES', ASYNCHRONOUS = 'NO'                , IOSTAT = IO_Write, ERR = 1006 ) (NDAN(I), NDPN(I),I = 1, NNDH) ;

UnFile = UN_HisV ;
Write (Unit = UnFile, FMT = "('History of velociy ')"                        , ADVANCE = 'YES', ASYNCHRONOUS = 'NO'                , IOSTAT = IO_Write, ERR = 1006 ) ;
Write (Unit = UnFile, FMT = "('Time',5x,<NNVH>('App#',I19,1x,'PETSc#',I19,1x))" , ADVANCE = 'YES', ASYNCHRONOUS = 'NO'                , IOSTAT = IO_Write, ERR = 1006 ) (NVAN(I), NVPN(I),I = 1, NNVH) ;

UnFile = UN_HisA ;
Write (Unit = UnFile, FMT = "('History of Acceleration and stresses ')"   , ADVANCE = 'YES', ASYNCHRONOUS = 'NO'                , IOSTAT = IO_Write, ERR = 1006 ) ;
Write (Unit = UnFile, FMT = "('Time',5x,<NNAH>('App#',I19,1x,'PETSc#',I19,1x))" , ADVANCE = 'YES', ASYNCHRONOUS = 'NO'                , IOSTAT = IO_Write, ERR = 1006 ) (NAAN(I),NAPN(I),I = 1, NNAH) ;

UnFile = UN_Engy ;
Write (Unit = UnFile, FMT = "('Total Energy of the regular domain')"   , ADVANCE = 'YES', ASYNCHRONOUS = 'NO'                , IOSTAT = IO_Write, ERR = 1006 ) ;
Write (Unit = UnFile, FMT = "('Time',5x,'Energy')" , ADVANCE = 'YES', ASYNCHRONOUS = 'NO'                , IOSTAT = IO_Write, ERR = 1006 ) ;

Write(*     ,*) 'End Subroutine < ResultsHeader >' ;
Write(UnInf,*) 'End Subroutine < ResultsHeader >' ;
Return ;

! =============================================== ERROR IN Write STATEMENT ==========================================================================

1006    Write(*       , Fmt_Write1 ) UnFile, IO_Write ; Write( UnFile, Fmt_Write1 ) UnFile, IO_Write ;
        !#Call BEEP_FAIL ;
        Write(*, Fmt_FL) ;  Write(UnInf, Fmt_FL) ; Write(*, Fmt_End) ; Read(*,*) ;  Return ;


End Subroutine ResultsHeader ;


! ***************************************************************************************************************************************************
!                                                                                                                                                  **
! Last Update:  06 Jan 2012                                                                                                                        **
! Description: THIS Subroutine CALCULATES                                                                                                          **
!                                                                                                                                                  **
! ***************************************************************************************************************************************************

Subroutine ResultsHis    (                                                                                                      &
NDim,     NNDH, NNVH, NNAH, IStep,     IStart, IEnd,                                                                            & ! Integer Variables
Time, TlEnergy,                                                                                                                 & ! Real Variables
EqDis, EqVel, EqAcc,                                                                                                            & ! Integer Arrays
U, UD, UDD                                                                                                                      & ! Real Arrays
!                                                                                                                               & ! Characters
!                                                                                                                               & ! Type 
) ;

!Implicit None ;

! =========================== Global Variables ======================================================================================================

! - Integer Variables -------------------------------------------------------------------------------------------------------------------------------
Integer , Intent(In)    :: NDim ;
Integer , Intent(In)    :: NNDH, NNVH, NNAH ;
Integer , Intent(In)    :: IStart, IEnd, IStep ;

!PetscInt       :: IStart, IEnd ;           ! Holds ownership of U_PTC


! - Real Variables ----------------------------------------------------------------------------------------------------------------------------------
Real(8), Intent(In)    :: Time, TlEnergy ;

! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
!#Complex, Intent(In)    ::  ;
!#Complex, Intent(InOut) ::  ;
!#Complex, Intent(OUT)   ::  ;
! - Integer Arrays ----------------------------------------------------------------------------------------------------------------------------------
Integer, Intent(In), Dimension (:  )  :: EqDis, EqVel, EqAcc ;

! - Real Arrays -------------------------------------------------------------------------------------------------------------------------------------
Real(8), Intent(In), Dimension (:  )  :: U, UD, UDD ;

! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
!#Complex, Intent(In)    ::  ;
!#Complex, Intent(InOut) ::  ;
!#Complex, Intent(OUT)   ::  ;
! - Character Variables -----------------------------------------------------------------------------------------------------------------------------
!#Character   ::  ;
! - Logical Variables -------------------------------------------------------------------------------------------------------------------------------
!#Logical   ::  ;
! - Type DECLERATIONS -------------------------------------------------------------------------------------------------------------------------------
!#Type() :: ;

! =========================== Local Variables =======================================================================================================
! - Integer Variables -------------------------------------------------------------------------------------------------------------------------------
Integer   :: IO_Write ;               ! Holds error of write statements
Integer   :: UnFile ;                ! Holds unit of a file for error message

Integer   :: I ;                      ! Loop indeces.

! - Real Variables ----------------------------------------------------------------------------------------------------------------------------------
!#Real (Kind=Dbl)  ::  ;
! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
!#Complex  ::  ;
! - Integer Arrays ----------------------------------------------------------------------------------------------------------------------------------
!#Integer (Kind=Shrt)  ::  ;
! - Real Arrays -------------------------------------------------------------------------------------------------------------------------------------
Real(8)  :: Dis ( NNDH * NDim ) ;        ! Holds Displacemets of the node in which the history is requested in each step
Real(8)  :: Vel ( NNVH * NDim ) ;        ! Holds Velocity of the node in which the history is requested in each step
Real(8)  :: Acc ( NNAH * NDim ) ;        ! Holds Acceleration of the node in which the history is requested in each step

! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
!#Complex  ::  ;
! - Character Variables -----------------------------------------------------------------------------------------------------------------------------
!#Character   ::  ;
! - Logical Variables -------------------------------------------------------------------------------------------------------------------------------
!#Logical   ::  ;

! =========================== Subroutine CODE =======================================================================================================

! - Displacements and stresses ----------------------------------------------------------------------------------------------------------------------
!  If (  IStep == 1 ) Then ;  !?? check this and modify it@@@@@@
!    If ( Any ( EqDis  < IStart ) .Or. Any ( EqDis  > IEnd ) ) Then ;
!      Write (*            , FMT= "('Warning: Equation Number of node in which requested for history of displacement does not belong to this rank - check the node numbers and mesh partitioning','   IStart = ',I19,'   IEnd = ',I19,<NNDH*NDOF>(I19,2x))") IStart, IEnd, EqDis ;
!      Write (Unit = UnInf, FMT= "('Warning: Equation Number of node in which requested for history of displacement does not belong to this rank - check the node numbers and mesh partitioning','   IStart = ',I19,'   IEnd = ',I19,<NNDH*NDOF>(I19,2x))") IStart, IEnd, EqDis ;
!
!    End If ;
!  End If ;

!  Do I = 1, NNDH * NDOF ;
  Do I = 1, NNDH * NDim ;
    If ( EqDis ( I ) == -1 ) Then ; ! Degree of Freedom is a constraint
      Dis ( I ) = 0.0d0 ;
    Else ;
      Dis ( I ) = U ( EqDis ( I ) - IStart + 1 ) ;
    End If ;
  End Do ;

UnFile = UN_HisD ;
If ( NNDH /= 0 ) Write (Unit = UnFile, FMT = "(E14.6E3,2x, <NNDH>(<NDim>(E18.10E3,2x),3x) )" , ADVANCE = 'YES', ASYNCHRONOUS = 'NO'                , IOSTAT = IO_Write, ERR = 1006 ) Time, ( Dis ( I ), I = 1, NNDH * NDim );

! - Velocity ----------------------------------------------------------------------------------------------------------------------------------------
!  If ( IStep == 1 ) Then ;
!    If ( Any ( EqVel < IStart ) .Or. Any ( EqVel > IEnd ) ) Then ;
!      Write (*            ,FMT= "('Warning: Equation Number of node in which requested for history of velocity does not belong to this rank - check the node numbers and mesh partitioning','IStart = ',I19,'  IEnd = ',I19,<NNDH*NDOF>(I19,2x))") IStart, IEnd, EqVel ;
!      Write (Unit = UnInf,FMT= "('Warning: Equation Number of node in which requested for history of velocity does not belong to this rank - check the node numbers and mesh partitioning','IStart = ',I19,'  IEnd = ',I19,<NNDH*NDOF>(I19,2x))") IStart, IEnd, EqVel ;
!    End If ;
!  End If ;

  Do I = 1, NNVH * NDim ;
    If ( EqVel ( I ) == -1 ) Then ; ! Degree of Freedom is a constraint
      Vel ( I ) = 0.0d0 ;
    Else ;
      Vel ( I ) = UD ( EqVel ( I ) - IStart + 1 ) ;
    End If ;
  End Do ;

UnFile = UN_HisV ;
If ( NNVH /= 0 ) Write (Unit = UnFile, FMT = "(E14.6E3,2x, <NNVH>(<NDim>(E18.10E3,2x),3x) )" , ADVANCE = 'YES', ASYNCHRONOUS = 'NO'                , IOSTAT = IO_Write, ERR = 1006 ) Time, ( Vel ( I ), I = 1, NNVH * NDim );

! - Acceleration ------------------------------------------------------------------------------------------------------------------------------------
!  If ( IStep == 1 ) Then ;
!    If ( Any ( EqAcc < IStart ) .Or. Any ( EqVel > IEnd ) ) Then ;
!      Write (*            ,FMT= "('Warning: Equation Number of node in which requested for history of Acceleration does not belong to this rank - check the node numbers and mesh partitioning','IStart = ',I19,'  IEnd = ',I19,<NNDH*NDOF>(I19,2x))") IStart, IEnd, EqAcc ;
!      Write (Unit = UnInf,FMT= "('Warning: Equation Number of node in which requested for history of Acceleration does not belong to this rank - check the node numbers and mesh partitioning','IStart = ',I19,'  IEnd = ',I19,<NNDH*NDOF>(I19,2x))") IStart, IEnd, EqAcc ;
!    End If ;
!  End If ;

  Do I = 1, NNAH * NDim ;
    If ( EqAcc ( I ) == -1 ) Then ; ! Degree of Freedom is a constraint
      Acc ( I ) = 0.0d0 ;
    Else ;
      Acc ( I ) = UDD ( EqAcc ( I ) - IStart + 1 ) ;
    End If ;
  End Do ;

UnFile = UN_HisA ;
If ( NNAH /= 0  ) Write (Unit = UnFile, FMT = "(E14.6E3,2x, <NNAH>(<NDim>(E18.10E3,2x),3x) )" , ADVANCE = 'YES', ASYNCHRONOUS = 'NO'                , IOSTAT = IO_Write, ERR = 1006 ) Time, ( Acc ( I ), I = 1, NNAH * NDim );

! - Energy ------------------------------------------------------------------------------------------------------------------------------------------
UnFile = UN_Engy ;
Write (Unit = UnFile, FMT = "(E18.10E3,2x,E31.23E3)" , ADVANCE = 'YES', ASYNCHRONOUS = 'NO'                , IOSTAT = IO_Write, ERR = 1006 ) Time, TlEnergy ;

!Write(*     ,*) 'End Subroutine < ResultsHis >' ;
!Write(UnInf,*) 'End Subroutine < ResultsHis >' ;
Return ;

! =============================================== ERROR IN Write STATEMENT ==========================================================================

1006    Write(*       , Fmt_Write1 ) UnFile, IO_Write ; Write( UnFile, Fmt_Write1 ) UnFile, IO_Write ;
        !#Call BEEP_FAIL ;
        Write(*, Fmt_FL) ;  Write(UnInf, Fmt_FL) ; Write(*, Fmt_End) ; Read(*,*) ;  Return ;


End Subroutine ResultsHis ;




End Module Results ;
