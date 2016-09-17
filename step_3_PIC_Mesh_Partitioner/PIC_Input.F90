
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                                                                                                                  ++
! Start:        01 Jul 2011                                                                                                                        ++
! Last Update:  2 April 2013                                                                                                                       ++
!                                                                                                                                                  ++
! Description: THIS Module Contains ALL Input Subroutines.                                                                                         ++
!                                                                                                                                                  ++
!              Subroutines                     ELEMENT NUMBER                                                                                      ++
!              ==============                  ==============                                                                                      ++
!              Input_BASIC                    COLLECTS ALL BASIC Variables.                                                                       ++
!              Input_Arrays                   COLLECTS ALL Arrays                                                                                 ++
!                                                                                                                                                  ++
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Module Input_Subroutines ;

Use Parameters ;

Implicit None ;

  Interface Input ;
    Module Procedure Input_BASIC, Input_Arrays ;
  End Interface ;

Contains ;

! ***************************************************************************************************************************************************
!                                                                                                                                                  **
! Start:        01 Oct 2011                                                                                                                        **
! Last Update:  2 April 2013                                                                                                                       **
! Description: THIS Subroutine COLLECTS ALL BASIC DATA FOR EACH ANALYSIS.                                                                          **
! Call by: Main.F90                                                                                                                                **
!                                                                                                                                                  **
! ***************************************************************************************************************************************************

Subroutine Input_BASIC (                                                                                                        &
NDOF, MaxNNode, NDim, Matrix_Type, MetisType,                                                                                   & ! Integer (1) Variables
NGroup, NPM, NMat, NLCase,                                                                                                      & ! Integer (2) Variables
NumFlag, NParts, KWay,                                                                                                          & ! Integer (4) Variables
NEl, NJ,                                                                                                                        & ! Integer (8) Variables
!                                                                                                                               & ! Real Variables
LoadC,                                                                                                                          & ! Integer Arrays
!                                                                                                                               & ! Real Arrays
!                                                                                                                               & ! Characters
Param                                                                                                                           & ! Type
) ;

Implicit None ;

! =========================== Global Variables ======================================================================================================

! - Integer Variables -------------------------------------------------------------------------------------------------------------------------------
Integer (Kind=Tiny), Intent(OUT)    :: NDOF, MaxNNode, NDim ;
Integer (Kind=Smll), Intent(OUT)    :: NGroup, NPM, NMat, NLCase ;
Integer (Kind=Lng ), Intent(OUT)    :: NEl, NJ ;

Integer (Kind=Tiny), Intent(OUT)    :: Matrix_Type, MetisType ;
Integer (Kind=Shrt), Intent(OUT)    :: NumFlag, NParts, KWay ;

! - Real Variables ----------------------------------------------------------------------------------------------------------------------------------
! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
! - Integer Arrays ----------------------------------------------------------------------------------------------------------------------------------
Integer (Kind=Tiny), Intent(OUT)  , Dimension (:  )  :: LoadC ;

! - Real Arrays -------------------------------------------------------------------------------------------------------------------------------------
! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
! - Character Variables -----------------------------------------------------------------------------------------------------------------------------
Character check*5 ;
! - Logical Variables -------------------------------------------------------------------------------------------------------------------------------
! - Type DECLERATIONS -------------------------------------------------------------------------------------------------------------------------------
Type ( BasicParam )    :: Param ;                ! Holds basic parameters of each load case

! =========================== LOCAL Variables =======================================================================================================

! - Integer Variables -------------------------------------------------------------------------------------------------------------------------------
Integer (Kind=Tiny)  :: I ;                     ! Loop Index

Integer (Kind=Smll)  :: IO_Read ;               ! Holds error of Read statements
Integer (Kind=Smll)  :: UnFile ;                ! Holds unit of a file for error message

! - Real Variables ----------------------------------------------------------------------------------------------------------------------------------
! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
! - Integer Arrays ----------------------------------------------------------------------------------------------------------------------------------
! - Real Arrays -------------------------------------------------------------------------------------------------------------------------------------
! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
! - Character Variables -----------------------------------------------------------------------------------------------------------------------------
! - Logical Variables -------------------------------------------------------------------------------------------------------------------------------
! - Type DECLERATIONS -------------------------------------------------------------------------------------------------------------------------------

! =========================== Subroutine CODE =======================================================================================================

Write(*,*)
Write(*,*)"<<<<<<<< Input Basic >>>>>>>>" ;

UnFile = UnInptMdl ;

!Read  (UnFile, *) check ; Write(*,*)'1', check ! 1
!Read  (UnFile, *) check ; Write(*,*)'2', check ! 2
!Read  (UnFile, *) check ; Write(*,*)'3', check ! 3
!Read  (UnFile, *) check ; Write(*,*)'4', check ! 4
!Read  (UnFile, *) check ; Write(*,*)'5', check ! 5
!Read  (UnFile, *) check ; Write(*,*)'6', check ! 6
!Read  (UnFile, *) check ; Write(*,*)'7', check ! 7

Read  (UnFile, *) ; ! 1
Read  (UnFile, *) ; ! 2
Read  (UnFile, *) ; ! 3
Read  (UnFile, *) ; ! 4
Read  (UnFile, *) ; ! 5
Read  (UnFile, *) ; ! 6
Read  (UnFile, *) ; ! 7
Read  (UnFile, *) NEl, NJ, NGroup, MaxNNode, NDOF, NDim ; ! 8
Read  (UnFile, *)

Write (*,*)"Basic Data" ;

! Indicating available load cases in the model and load case number
! Load Case. Indicates if we have static or dynamic loads or both. if 0 then the load is not in the model
! Static locations are 1: body force - 2: pressure - 3: concentrated force - 4: specified Displacements
! Dynamic load 5: 1: for dynamic pressure, 2: base acceleration 3:DRM
Read  (UnFile, *) ;
Read  (UnFile, *) ;
Read  (UnFile, *) ( LoadC (I), I = 1, 5 ) ;
Read  (UnFile, *) ;

NLCase = MAXVAL ( LoadC(1:4) ) ;
NPM    = 8 ;
NMat   = NGroup ;

Write(*,*)"Load Cases" ;

Allocate ( Param%IntM ( 7, 6),Param%RealM ( 7, 6)  ) ;

! - Load Parameters ---------------------------------------------------------------------------------------------------------------------------------
! initialize Load Parameters
Param%IntM  (:,:) = 0_Lng ;
Param%RealM (:,:) = 0._Dbl ;

  If (LoadC (1) /= 0_Tiny ) Then ;  ! Body Force
    Read  (UnInptMdl, *);
    Read  (UnInptMdl, *) Param%IntM( 1, 1), Param%IntM( 1, 2) ; ! NBLD   = NGroup, NPBL   = NDim ;
    Read  (UnInptMdl, *);
  End If ;

  If (LoadC (2) /= 0_Tiny ) Then ;  ! Pressure on elements
    Read  (UnInptMdl, *);
    !Read  (UnInptMdl, *) Param%IntL( 2, 1), Param%IntL( 2, 2), Param%IntL( 2, 3), Param%IntM( 2, 4),     Param%RealL( 2, 1), Param%RealL( 2, 2) ; ! PR_Type, LDIR, PRZERO, NIDBC        HREF, PR
    Read  (UnInptMdl, *) Param%IntM( 2, 4) ; ! NIDBC 
    Read  (UnInptMdl, *);
  End If ;

  If (LoadC (3) /= 0_Tiny ) Then ;  ! Concentrated force
    Read  (UnInptMdl, *) ;
    Read  (UnInptMdl, *) Param%IntM( 3, 1) ; ! NLN
    Read  (UnInptMdl, *) ;
  End If ;

  If (LoadC (4) /= 0_Tiny ) Then ;  ! Number of nodes (Supports) with predefined Displacements
    Read  (UnInptMdl, *) ;
    Read  (UnInptMdl, *) Param%IntM( 4, 1) ; ! NSDN
    Read  (UnInptMdl, *) ;
  End If ;

!  If ( LoadC (5) == 0_Tiny ) Then ;  ! Static Analysis
!    Read  (UnInptMdl, *) ;
!    Read  (UnInptMdl, *) Param%IntL( 5, 1 ) ; ! STStep
!    Read  (UnInptMdl, *) ;
!  End If ;

Write(*,*)"Static Loads" ;

  If ( LoadC (5) /= 0_Tiny ) Then ;  ! Dynamic Analysis
    Read  (UnInptMdl, *) ;
    !Read  (UnInptMdl, *) ( Param%IntP( 5, I), I = 1, 4 ), ( Param%RealP( 5, I), I = 1, 4 ) ; ! NSTEP, NTST, LoadType, NNodeDown          Delta, Gama, DT, t0,         :   Newmark Coefficients, time
    !Read  (UnInptMdl, *) ( Param%IntL( 5, I), I = 1, 2 ), ( Param%IntM( 5, I), I = 3, 3 ), ( Param%IntL( 5, I), I = 4, 6 ),        ( Param%RealL( 5, I), I = 3, 4 ) ; ! NSTEP, NTST, LoadType, NNodeDown, Solver(Implicit/Explicit 0/1), Analysis(Dynamics of structure/wave propagation 0/1),       DT, t0,            time
    Read  (UnInptMdl, *) Param%IntM( 5, 3), Param%IntM( 7, 6) ; ! LoadType,       PML existence 0: no PML in the domain 1: PML elements exist in the domain   ! pml analysis
    Read  (UnInptMdl, *) ;

    Read  (UnInptMdl, *) ;
    !Read  (UnInptMdl, *) ( Param%IntM( 6, I), I = 1, 3 ), ( Param%IntL( 6, I), I = 4, 6 ), Param%IntM( 7, 6) ; ! NNDH, NNVH, NNAH, NDamp, NEnergy, PARAM_Type, PML existence 0: no PML in the domain 1: PML elements exist in the domain   ! pml analysis
    Read  (UnInptMdl, *) ( Param%IntM( 6, I), I = 1, 6 ) ; ! NNDH, NNVH, NNAH, NDamp, NEnergy, PARAM_Type 
    Read  (UnInptMdl, *) ;

    Write(*,*)"Dynamic Loads" ;

      ! Loads
      If      ( Param%IntM( 5, 3) == 1_Tiny ) Then ;  ! Dynamic Pressure
        Read  (UnInptMdl, *) ;
        !Read  (UnInptMdl, *) Param%IntL( 2, 1), Param%IntL( 2, 2), Param%IntL( 2, 3), Param%IntM( 2, 4),     Param%RealL( 2, 1), Param%RealL( 2, 2) ; ! PR_Type, LDIR, PRZERO, NIDBC        HREF, PR
        Read  (UnInptMdl, *) Param%IntM( 2, 4) ; ! NIDBC
        Read  (UnInptMdl, *) ;

!        Read  (UnInptMdl, *) ;
!        Read  (UnInptMdl, *) Param%IntL( 7, 1) ; ! FUNC_Type  - SINE FUNCTION OR RICKER PULSE OR ... 
!        Read  (UnInptMdl, *) ;

        Write(*,*)"Dynamic pressure" ;

      Else If ( Param%IntM( 5, 3) == 2_Tiny ) Then ; ! Base Acceleration
!        Read  (UnInptMdl, *) ;
!        Read  (UnInptMdl, *) Param%IntL( 7, 1), Param%RealL( 7, 1) ; ! NASTEP, G   Number of base Accelerations to be read in the input file
!        Read  (UnInptMdl, *) ;

        Write(*,*)"Base Acceleration" ;

      Else If ( Param%IntM( 5, 3) == 3_Tiny ) Then ; ! Domain Reduction Method loading
        Read  (UnInptMdl, *) ;
        !Read  (UnInptMdl, *) ( Param%IntM( 7, I), I = 1, 2 ), ( Param%IntL( 7, I), I = 3, 4 ), Param%RealL( 7, 1) ; ! NNBndry_DRM, NNLayer_DRM, Wave_Type, wave_func     AngleInci: angle of incident wave plane with the system coordinate ;
        Read  (UnInptMdl, *) ( Param%IntM( 7, I), I = 1, 2 ) ; ! NNBndry_DRM, NNLayer_DRM 
        Read  (UnInptMdl, *) ;

        Write(*,*)"DRM" ;
      End If ;

      !  Implicit Solver
!      If      ( Param%IntL( 5, 5) == 0_Tiny ) Then ;
!
!        Read  (UnInptMdl, *);
!        Read  (UnInptMdl, *) ( Param%RealL( 5, I), I = 1, 2 ) ; ! Delta, Gama        Newmark Coefficients
!        Read  (UnInptMdl, *);
!
!      End If ;

  End If ;

Write(*,*)"Param Array" ;

! Partitioning Information regarding to the MeshPartitioner Program - ParMetis
READ  (UnInptMdl, *)
READ  (UnInptMdl, *)
READ  (UnInptMdl, *) NParts, &       ! NParts: Number of partitions
                  NumFlag, &      ! NumFlag Used to indicate which numbering scheme is Used for the element node array. NumFlag can take the following two values:  0 C-style numbering is assumed that starts from 0 - 1 Fortran-style numbering is assumed that starts from 1
                  Matrix_Type, &  ! Accuracy of calculations of Number of non-zero entries of the petsc objects. 0: approximate. 1: exact
                  MetisType, &    ! Metis Type 0: Metis. 1: ParMetis
                  KWay ;          ! Number of input partitions
READ  (UnInptMdl, *)


Write(*    ,*) 'End Subroutine < Input_BASIC >' ;
Write(UnInf,*) 'End Subroutine < Input_BASIC >' ;
Return ;
End Subroutine Input_BASIC ;


! ***************************************************************************************************************************************************
!                                                                                                                                                  **
! Last Update:  2 April 2013                                                                                                                       **
! Description: THIS Subroutine COLLECTS ALL Arrays.                                                                                                **
! Called by: Main - case ( ACN_LD_3 )                                                                                                              **
!                                                                                                                                                  **
! ***************************************************************************************************************************************************

Subroutine Input_Arrays  (                                                                                                      &
NDOF, NDim, MaxNNode,                                                                                                           & ! Integer (1) Variables
NPM, NMat,                                                                                                                      & ! Integer (2) Variables
!                                                                                                                               & ! Integer (4) Variables
NEl, NJ,                                                                                                                        & ! Integer (8) Variables
!                                                                                                                               & ! Real Variables
LoadC, IDBC,     MTEL, ELT, ELGR,         JLoad,  NDAN, NVAN, NAAN, NoBndry_DRM, NoLayer_DRM, INOD, ID,                         & ! Integer Arrays
Param,           PMat, PBLD, XYZ, UDis, PLoad,  PML_DIM                                                                         & ! Real Arrays
!                                                                                                                               & ! Characters
!                                                                                                                               & ! Type
) ;

Implicit None ;

! =========================== Global Variables ======================================================================================================

! - Integer Arrays ------------------------------------------------------------------------------------------------------------------------------
Integer (Kind=Tiny), Intent( IN )    :: NDOF, NDim, MaxNNode ;
Integer (Kind=Smll), Intent( In )    :: NPM, NMat ;
Integer (Kind=Lng ), Intent( IN )    :: NEl, NJ ;

! - Real Variables ----------------------------------------------------------------------------------------------------------------------------------
! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------

! - Integer Arrays ----------------------------------------------------------------------------------------------------------------------------------
Integer (Kind=Tiny), Intent(In )  , Dimension (:  )  :: LoadC ;

Integer (Kind=Smll), Intent(OUT)  , Dimension (:  )  :: MTEL, ELT ;

Integer (Kind=Shrt), Intent(OUT)  , Dimension (:  )  :: ELGR ;

Integer (Kind=Lng ), Intent(OUT)  , Dimension (:  )  :: JLoad, NDAN, NVAN, NAAN, NoBndry_DRM, NoLayer_DRM ;
Integer (Kind=Lng ), Intent(OUT)  , Dimension (:,:)  :: INOD, ID ;
Integer (Kind=Lng ), Intent(OUT)  , Dimension (:,:)  :: IDBC ;

! - Real Arrays -------------------------------------------------------------------------------------------------------------------------------------
Real (Kind=DBL), Intent(OUT)  , Dimension (:,:)  :: PMat, PBLD, XYZ, UDis, PLoad, PML_DIM ;

! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
! - Character Variables -----------------------------------------------------------------------------------------------------------------------------
! - Logical Variables -------------------------------------------------------------------------------------------------------------------------------
! - Type DECLERATIONS -------------------------------------------------------------------------------------------------------------------------------
Type ( BasicParam )    :: Param ;                ! Holds basic parameters of each load case

! =========================== LOCAL Variables =======================================================================================================
! - Integer Variables -------------------------------------------------------------------------------------------------------------------------------
Integer (Kind=Smll)  :: IO_Read ;                ! Holds error of Read statements.
Integer (Kind=Smll)  :: UnFile ;                ! Holds Unit of a file for error message.
Integer (Kind=Smll)  :: INode ;                  ! Loop index on NNode.
Integer (Kind=Smll)  :: NPBL ;                   ! Number of properties of body force load
Integer (Kind=Smll)  :: NBLD ;                   ! Number of body force loads

Integer (Kind=Shrt)  :: I1, I2 ;                 ! Loop indeces
Integer (Kind=Shrt)  :: NGR ;                    ! Number of Groups
Integer (Kind=Shrt)  :: NGEL ;                   ! Number of elements of each group
Integer (Kind=Shrt)  :: IElType ;                ! Load type of each element
Integer (Kind=Shrt)  :: LDSET ;                  ! Load Set
Integer (Kind=Shrt)  :: INC ;                    ! Increment of load type
Integer (Kind=Shrt)  :: LType ;                  ! Load type
Integer (Kind=Shrt)  :: NODE ;                   ! Node number for Reading constraints
Integer (Kind=Shrt)  :: NNBndry_DRM ;            ! Number of nodes on the DRM boundary for the Domain Reduction Method
Integer (Kind=Shrt)  :: NNLayer_DRM ;            ! Number of nodes on the DRM layer for the Domain Reduction Method
Integer (Kind=Shrt)  :: NLN ;                    ! Number of Loaded Nodes (Concentrated forces on nodes)
Integer (Kind=Shrt)  :: NSND ;                   ! Number of nodes (Supports) with predefined Displacements.
Integer (Kind=Shrt)  :: NASTEP ;                 ! Number of Steps of base acceleration in input file.
Integer (Kind=Shrt)  :: NTST ;                   ! Number of Steps (Times) which full informa
Integer (Kind=Shrt)  :: NNDH ;                   ! Number of Nodes which History of Displacement is required in dynamic analysis.
Integer (Kind=Shrt)  :: NNVH ;                   ! Number of Nodes which History of Velocity is required in dynamic analysis.
Integer (Kind=Shrt)  :: NNAH ;                   ! Number of Nodes which History of Acceleration is required in dynamic analysis.

Integer (Kind=Lng )  :: I, J, K ;                ! Loop indeces.
Integer (Kind=Lng )  :: NEQUATION ;              ! Counter on the equation number.
Integer (Kind=Lng )  :: IJ ;                     ! Loop index on NJ.
Integer (Kind=Lng )  :: IEl ;                    ! Loop index on NEl.

! - Real Variables ----------------------------------------------------------------------------------------------------------------------------------
! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
! - Integer Arrays ----------------------------------------------------------------------------------------------------------------------------------
! - Real Arrays -------------------------------------------------------------------------------------------------------------------------------------
! - Complex Variables -------------------------------------------------------------------------------------------------------------------------------
! - Character Variables -----------------------------------------------------------------------------------------------------------------------------
! - Logical Variables -------------------------------------------------------------------------------------------------------------------------------
! - Type DECLERATIONS -------------------------------------------------------------------------------------------------------------------------------

! =========================== Subroutine CODE =======================================================================================================

Write(*,*)
Write(*,*)"<<<<<<<< Input Arrays >>>>>>>>>" ;

! - Geometry ----------------------------------------------------------------------------------------------------------------------------------------

Write (*,*)"Reading Geometry ..." ;

Read  (UnInptMdl, *) ;
Read  (UnInptMdl, *) ;


!Write(*,*)"Group numbers ..." ;
!
!! Number of elements of each group
!Read  (UnInptMdl, *) ;
!Read  (UnInptMdl, *) ;
!  DO I = 1, NGroup ;
!    Read  (UnInptMdl, *) NGR, NGEL, IElType ; ! GROUP NUMBER, NUMBER OF ELEMENTS OF THIS GROUP , ELEMENT Type
!    NEG ( NGR, 1 ) = NGEL ; NEG ( NGR, 2 ) = IElType ;  !?? check this one, NOT REQUIRED
!  End Do ;
!Read  (UnInptMdl, *) ;

Write (*,*)"Material Properties ..." ;

! Material properties of each material
Read  (UnInptMdl, *) ;
Read  (UnInptMdl, *) ;
Read  (UnInptMdl, *) ;
  DO I = 1, NMat ;
    Read  (UnInptMdl, *) K, ( PMat ( K, J ), J = 1, NPM ) ; 
  End Do ;
Read  (UnInptMdl, *) ;


! PML TERRITORY
  If ( LoadC (5) /= 0_Tiny .AND. Param%IntM( 7, 6) == 1_Lng) Then ; ! Param%IntM( 7, 6) = pml exists
    Write (*,*)"PML teritory ..." ;
    Read  (UnInptMdl, *) ;
      DO J = 1, 2
        Read  (UnInptMdl, *) ;
        Read  (UnInptMdl, *) (PML_DIM ( I, J ), I = 1, 2 * NDim ) ; 
      End Do ;
    Read  (UnInptMdl, *) ;
  End If ;

! Coordinates of each node
Write (*,*)"Coordinates ..." ;
  DO K = 1, NJ ;
    Read  (UnInptXYZ, *) I, ( XYZ ( I, J ), J = 1, NDim ) ; 
  End Do ;

! Element connectivity
Write (*,*)"Element connectivities ..." ;
!INod = -1_Lng ;   ! We use this to detect those elements that have less nodes than MaxNNode
  ! Mind node's orientation of elements
  DO J = 1, NEl ;
    Read  (UnInptCnn, *) IEl, ( INOD ( I, IEl ), I = 1, MaxNNode ), MTEL ( IEl ), ELT ( IEl ), ELGR ( IEl ) ;
  End Do ;

! Constraints
Write (*,*)"Node constraints ..." ;
  DO J = 1, NJ ;
    Read  (UnInptCnt, *) NODE, ( ID ( NODE, I ), I = 1, NDOF ) ;
  End Do ;

! - Loads -------------------------------------------------------------------------------------------------------------------------------------------
Write (*,*)"Reading Load types ..." ;

! Body Force Load
  If ( LoadC (1) /= 0_Tiny ) Then ;
    ! Load types
!    LTEL=0
!    Read  (UnInptMdl, *) ;
!    Read  (UnInptMdl, *) LDSET ;
!      DO I = 1, LDSET ;
!        Read  (UnInptMdl, *) I1, I2, INC, LType ;
!          DO J = I1, I2, INC ;
!            LTEL ( J ) = LType ;
!          End Do ;
!      End Do ;
!    Read  (UnInptMdl, *) ;

    ! Body force load
    NBLD = Param%IntM( 1,1) ;
    NPBL = Param%IntM( 1,2) ;
    Read  (UnInptMdl, *) ;
    Read  (UnInptMdl, *) ;
      DO I = 1, NBLD ;
        Read  (UnInptMdl, *) ( PBLD ( I, J ) , J = 1, NPBL ) ;
      End Do
    Read  (UnInptMdl, *) ;
    Write(*,*)"Body Forces" ;
  End If ;

! Static Pressure
  If ( LoadC (2) /= 0_Tiny ) Then ;
    Read  (UnInptMdl, *) ;
      DO J = 1, NEl ;
        Read  (UnInptMdl, *) IEl, ( IDBC ( IEl, I ), I = 1, 2 * NDim ) ;
      End Do ;
    Read  (UnInptMdl, *) ;
    Write(*,*)"Pressure Loads ..." ;
  End If ;

! Concentrated Load
  If ( LoadC (3) /= 0_Tiny ) Then ;
    NLN = Param%IntM( 3,1) ; ! Number of Loaded Nodes
    Read  (UnInptMdl, *) ;
    Read  (UnInptMdl, *) ;
      IF ( LoadC ( 3 ) /= 0 ) Then ;
          DO I = 1, NLN ;
            Read  (UnInptMdl, *)  JLoad( I ), ( PLoad( J, I ), J = 1, NDim ) ;
          End Do
      EndIF
    Read  (UnInptMdl, *) ;
    Write(*,*)"Joint Loads" ;
  End If ;

! Supports displacements
  If ( LoadC (4) /= 0_Tiny ) Then ;
    NSND = Param%IntM( 4, 1 ) ; ! Number of Supports 
    Read  (UnInptMdl, *) ;
    Read  (UnInptMdl, *) ;
      IF ( LoadC ( 4 ) /= 0 ) Then ;
        DO I = 1, NSND ;
          Read  (UnInptMdl, *)  UDis( NDOF + 1, I ), ( UDis ( K, I ), K = 1, NDOF ) ;
        End Do ;
      EndIF ;
    Read  (UnInptMdl, *) ;
    Write(*,*)"Support Displacements" ;
  End If ;

Write(*,*)"Reading Dynamic Loads ..." ;

! Dynamic loads
  If ( LoadC (5) /= 0_Tiny ) Then ;

    If ( Param%IntM( 5, 3 ) == 1 ) Then ;    ! RICKER PULSE or sine function for dynamic pressure

      Read  (UnInptMdl, *) ;
        DO J = 1, Param%IntM ( 2, 4) ; ! Param%IntP ( 2, 4) = NIDBC
          Read  (UnInptMdl, *) IDBC ( J, 1 ), ( IDBC ( J, I ), I = 2, 2 * NDim + 1 ) ;
        End Do ;
      Read  (UnInptMdl, *) ;

!    Else If ( Param%IntM( 5, 3 ) == 2_Tiny ) Then ;   ! Base acceleration
!
!      NAStep = Param%IntL( 7, 1) ;
!      Read  (UnInptMdl, *) ;
!      Read  (UnInptMdl, *) ;
!        DO I = 1, NASTEP ;
!          Read  (UnInptMdl, *) ( BACL ( I, K ), K = 1, NDim ) ;
!        End Do ;

    Else If ( Param%IntM( 5, 3 ) == 3_Tiny ) Then ;    ! Domain Reduction Method (DRM)

!      Read  (UnInptMdl, *) ;
!      Read  (UnInptMdl, *) ;
!      Read  (UnInptMdl, *) InciWave(1), InciWave(2), InciWave(3), InciWave(4), InciWave(5) ;       ! Theta, Omega, amplitude, alpha1, alpha2 ( limits of the phase if wave_type = 0, i.e., sine function - alpha1 = central frequency if wave_type = 1, i.e., ricker pulse )
!      Read  (UnInptMdl, *) ;

      Write(*,*)"Reading nodes for DRM ..." ;

      ! Reading the node numbers on the DRM boundary
      NNBndry_DRM = Param%IntM( 7, 1) ;
      Read  (UnInptDRM, *) ;
        Do I = 1, NNBndry_DRM ;
          Read  (UnInptDRM, *) NoBndry_DRM ( I ) ;
        End Do ;
      !Read  (UnInptMdl, *) ;

      NNLayer_DRM = Param%IntM( 7, 2) ;
      ! Reading the node numbers on the DRM layer
      Read  (UnInptDRM, *) ;
        Do I = 1, NNLayer_DRM ;
          Read  (UnInptDRM, *) NoLayer_DRM ( I ) ;
        End Do ;
      !Read  (UnInptMdl, *) ;

    End If ;

  End If ;

  If ( LoadC (5) /= 0_Tiny ) Then ;
!    ! HISTORY OF NODE NUMBERS 
!    NTST = Param%IntL( 5, 2) ;
!    Read  (UnInptMdl, *) ;
!    Read  (UnInptMdl, *) ( STEP ( I ), I = 1, NTST ) ;
!    Read  (UnInptMdl, *) ;
!    Write(*,*)"Step Numbers" ;

    Write(*,*)"Reading node numbers for history ..." ;

    Write(*,*)"Displacements ..." ;
    NNDH = Param%IntM( 6, 1 ) ;
    Read  (UnInptMdl, *) ;
    Read  (UnInptMdl, *) ( NDAN ( I ), I = 1, NNDH ) ;
    Read  (UnInptMdl, *) ;
    Write(*,*)"Node Numbers for displacements" ;

    Write(*,*)"Displacements ..." ;
    NNVH = Param%IntM( 6, 2 ) ;
    Read  (UnInptMdl, *) ;
    Read  (UnInptMdl, *) ( NVAN ( I ), I = 1, NNVH ) ;
    Read  (UnInptMdl, *) ;
    Write(*,*)"Node Numbers for velocity" ;

    Write(*,*)"Displacements ..." ;
    NNAH = Param%IntM( 6, 3 ) ;
    Read  (UnInptMdl, *) ;
    Read  (UnInptMdl, *) ( NAAN ( I ), I = 1, NNAH ) ;
    Read  (UnInptMdl, *) ;
    Write(*,*)"Node Numbers for acceleration" ;
  End If ;

Write(*    ,*) 'End Subroutine < Input_Arrays >' ;
Write(UnInf,*) 'End Subroutine < Input_Arrays >' ;
Return ;

! - ERROR IN Read STATEMENT -------------------------------------------------------------------------------------------------------------------------
1003    Write(*       , Fmt_Read1 ) UnFile, IO_Read ; Write( UnFile, Fmt_Read1 ) UnFile, IO_Read ;
        !#Call BEEP_FAIL ;
        Write(*, Fmt_FL) ;  Write(UnInf, Fmt_FL) ; Write(*, Fmt_End) ; Read(*,*) ;  Return ;

! - End-OF-FILE IN Read STATEMENT -------------------------------------------------------------------------------------------------------------------
1004    Write(*       , Fmt_Read2 ) UnFile, IO_Read ; Write( UnFile, Fmt_Read2 ) UnFile, IO_Read ;
        !#Call BEEP_FAIL ;
        Write(*, Fmt_FL) ;  Write(UnInf, Fmt_FL) ; Write(*, Fmt_End) ; Read(*,*) ;  Return ;

! - End-OF-FILE IN Read STATEMENT -------------------------------------------------------------------------------------------------------------------
1005    Write(*       , Fmt_Read3 ) UnFile, IO_Read ; Write( UnFile, Fmt_Read3 ) UnFile, IO_Read ;
        !#Call BEEP_FAIL ;
        Write(*, Fmt_FL) ;  Write(UnInf, Fmt_FL) ; Write(*, Fmt_End) ; Read(*,*) ;  Return ;

End Subroutine Input_Arrays ;


End Module Input_Subroutines ;

