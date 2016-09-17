
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                                                                                                                  ++
! Last Update:  10 Jan 2012                                                                                                                       ++
!                                                                                                                                                  ++
! Description: THIS Module WriteS DOWN ALL INFORMATION REGARDING THE CURRENT ANALYSIS.                                                             ++
!                                                                                                                                                  ++
!              Subroutines                     Description                                                                                         ++
!              ==============                  ==============                                                                                      ++
!              InfBasic                       WriteS DOWN TIME AND DATE                                                                           ++
!                                                                                                                                                  ++
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
Module Information ;

Use Parameters ;

Implicit None ;

  Interface Info ;
    Module Procedure InfBasic, InfTime ;
  End Interface  INFO;


Contains ;

! ***************************************************************************************************************************************************
!                                                                                                                                                  **
! Last Update:  07 JUNE 2011                                                                                                                       **
! Description: THIS Subroutine WriteS DOWN TIME AND DATE OF ANALYSIS.                                                                              **
!                                                                                                                                                  **
! ***************************************************************************************************************************************************

Subroutine InfBasic ( &
Iyr, Imon, Iday, Ih, Im, Is, I100th,            & ! Integer Variables
!                                               & ! Real Variables
!                                               & ! Integer Arrays
!                                               & ! Real Arrays
Name, Model_InDir, OutDir, InlDir                  & ! Characters
) ;

Implicit None ;

! =========================== Global Variables ======================================================================================================

! - Integer Arrays ------------------------------------------------------------------------------------------------------------------------------
Integer (Kind=Smll), Intent(In)    :: Iyr, Imon, Iday, Ih, Im, Is, I100th ;
!#Integer (Kind=Shrt), Intent(InOut) ::  ;
!#Integer (Kind=Shrt), Intent(OUT)   ::  ;
! - Real Variables ----------------------------------------------------------------------------------------------------------------------------------
!#Real (Kind=DBL), Intent(In)    ::  ;
!#Real (Kind=DBL), Intent(InOut) ::  ;
!#Real (Kind=DBL), Intent(OUT)   ::  ;
! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
!#Complex, Intent(In)    ::  ;
!#Complex, Intent(InOut) ::  ;
!#Complex, Intent(OUT)   ::  ;
! - Integer Arrays ----------------------------------------------------------------------------------------------------------------------------------
!#Integer (Kind=Shrt), Intent(In), Dimension (:  )  ::
!#Integer (Kind=Shrt), Intent(In), Dimension (:,:)  ::
!#Integer (Kind=Shrt), Intent(In)    ::  ;
!#Integer (Kind=Shrt), Intent(InOut) ::  ;
!#Integer (Kind=Shrt), Intent(OUT)   ::  ;
! - Real Arrays -------------------------------------------------------------------------------------------------------------------------------------
!#Real (Kind=DBL), Intent(In)    ::  ;
!#Real (Kind=DBL), Intent(InOut) ::  ;
!#Real (Kind=DBL), Intent(OUT)   ::  ;
! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
!#Complex, Intent(In)    ::  ;
!#Complex, Intent(InOut) ::  ;
!#Complex, Intent(OUT)   ::  ;
! - Character Variables -----------------------------------------------------------------------------------------------------------------------------
Character (Kind = 1, Len = 30 ) :: Name ;        ! Name of Input file
Character (Kind = 1, Len = 100) :: Model_InDir ;       ! Directory of input file.
Character (Kind = 1, Len = 100) :: InlDir ;      ! Directory of internal files.
Character (Kind = 1, Len = 100) :: OutDir ;      ! Directory of output files (Results)

! - Logical Variables -------------------------------------------------------------------------------------------------------------------------------
!#Logical   ::  ;
! =========================== LOCAL Variables =======================================================================================================
! - Integer Arrays ------------------------------------------------------------------------------------------------------------------------------
!#Integer (Kind=Shrt)  ::  ;
! - Real Variables ----------------------------------------------------------------------------------------------------------------------------------
!#Real (Kind=DBL)  ::  ;
! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
!#Complex  ::  ;
! - Integer Arrays ----------------------------------------------------------------------------------------------------------------------------------
!#Integer (Kind=Shrt)  ::  ;
! - Real Arrays -------------------------------------------------------------------------------------------------------------------------------------
!#Real (Kind=DBL)  ::  ;
! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
!#Complex  ::  ;
! - Character Variables -----------------------------------------------------------------------------------------------------------------------------
!#Character   ::  ;
! - Logical Variables -------------------------------------------------------------------------------------------------------------------------------
!#Logical   ::  ;
! - Type DECLERATIONS -------------------------------------------------------------------------------------------------------------------------------
!#Type() :: ;

! =========================== Subroutine CODE =======================================================================================================

! Write INFORMATION 
Write (UnInf, Fmt_DATE) Imon, Iday, Iyr,   Ih, Im, Is, I100th ;
Write (UnInf, Fmt_NM  ) Name, Model_InDir, OutDir, InlDir ;

Write (UnInf,*)" Analysis Type " ;
! Linear Analysis
! Time Domain
! Direct Approach
! Regular + Solid Elements
! PETSc library
! Direct Solver

Write(*     ,*) 'End Subroutine < InfBasic >' ;
!#Write(UnInf,*) 'End Subroutine < InfBasic >' ;
Return ;
End Subroutine InfBasic ;


! ***************************************************************************************************************************************************
!                                                                                                                                                  **
! Last Update:  10 Jan 2012                                                                                                                        **
! Description: THIS Subroutine writes down all timing information.                                                                                 **
!                                                                                                                                                  **
! ***************************************************************************************************************************************************

Subroutine InfTime  (                                                                                                           &
Nel,                                                                                                                            & ! Integer Variables
TimeE, TimeS, TimeInputE, TimeInputS, TimeIndexE, TimeIndexS, TimeAssemE, TimeAssemS, TimeEStiffE, TimeEStiffS, TimeSolveE, TimeSolveS & ! Real Variables
!                                                                                                                               & ! Integer Arrays
!                                                                                                                               & ! Real Arrays
!                                                                                                                               & ! Characters
!                                                                                                                               & ! Type 
) ;

Implicit None ;

! =========================== Global Variables ======================================================================================================

! - Integer Variables -------------------------------------------------------------------------------------------------------------------------------
Integer (Kind=Lng ), Intent(In)    :: NEl ;

! - Real Variables ----------------------------------------------------------------------------------------------------------------------------------
Real (Kind=Dbl), Intent(In)    :: TimeE, TimeS, TimeInputE, TimeInputS, TimeIndexE, TimeIndexS, TimeAssemE, TimeAssemS, TimeEStiffE, TimeEStiffS, TimeSolveE, TimeSolveS ;

! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
!#Complex, Intent(In)    ::  ;
!#Complex, Intent(InOut) ::  ;
!#Complex, Intent(OUT)   ::  ;
! - Integer Arrays ----------------------------------------------------------------------------------------------------------------------------------
!#Integer (Kind=Shrt), Intent(In), Dimension (:  )  :: ;
!#Integer (Kind=Shrt), Intent(In), Dimension (:,:)  :: ;
!#Integer (Kind=Shrt), Intent(In)    ::  ;
!#Integer (Kind=Shrt), Intent(InOut) ::  ;
!#Integer (Kind=Shrt), Intent(OUT)   ::  ;
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
!#Integer (Kind=Shrt)  ::  ;
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

Write(UnInf,*)

Write(*     ,Fmt_RUNTIME) "TOTAL"   , TimeE - TimeS ;

Write(UnInf,*)"---------- RUNNING TIME STATISTICS ----------" ;

Write(UnInf,Fmt_RUNTIME) "Reading Input files           ", TimeInputE  - TimeInputS ;
Write(UnInf,Fmt_RUNTIME) "Object Definition & Indexing  ", TimeIndexE  - TimeIndexS ;
Write(UnInf,Fmt_RUNTIME) "Global Matrix Assembly        ", TimeAssemE  - TimeAssemS ;
Write(UnInf,Fmt_RUNTIME) "Assemble per Element          ", (TimeAssemE - TimeAssemS) / NEl ;
Write(UnInf,Fmt_RUNTIME) "Effective Stiffness Matrix    ", TimeEStiffE - TimeEStiffS ;
Write(UnInf,Fmt_RUNTIME) "SOLVE                         ", TimeSolveE  - TimeSolveS ;
Write(UnInf,Fmt_RUNTIME) "TOTAL                         ", TimeE       - TimeS ;
Write(UnInf,*)

Write(*     ,*) 'End Subroutine < InfTime >' ;
!#Write(UnInf,*) 'End Subroutine < InfTime >' ;

Return ;
End Subroutine InfTime ;

End Module Information ;
