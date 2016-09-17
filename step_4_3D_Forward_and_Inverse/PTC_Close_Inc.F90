
! - Close ADDRESS FILE ------------------------------------------------------------------------------------------------------------------------------
 UnFile =  UN_ADR ;
 Close ( Unit = UnFile, STATUS = 'KEEP', ERR = 1002, IOSTAT = IO_File )

! Data file of model
 UnFile =  UnInptMdl ;
 Close ( Unit = UnFile, STATUS = 'KEEP', ERR = 1002, IOSTAT = IO_File ) ;

! Material Properties (.Mat)
 UnFile = UnInptMat ;
 Close ( Unit = UnFile, STATUS = 'KEEP', ERR = 1002, IOSTAT = IO_File ) ;

! Inversion Data Structure (.Inv_DS)
 UnFile = Un_Inversion_DS ;
 Close ( Unit = UnFile, STATUS = 'KEEP', ERR = 1002, IOSTAT = IO_File ) ;


! - Close Input FILE --------------------------------------------------------------------------------------------------------------------------------

  If ( Output_Type == 3 ) Then ;  ! HDF5 format



  Else ;   ! Binary, Formatted

    ! Coordinate file (.XYZ)
    UnFile = UnInptXYZ ;
    Close ( Unit = UnFile, STATUS = 'KEEP', ERR = 1002, IOSTAT = IO_File ) ;

    ! Connectivity file (.Cnn)
    UnFile = UnInptCnn ;
    Close ( Unit = UnFile, STATUS = 'KEEP', ERR = 1002, IOSTAT = IO_File ) ;

    ! Constraint file (.Cnt)
    UnFile = UnInptCnt ;
    Close ( Unit = UnFile, STATUS = 'KEEP', ERR = 1002, IOSTAT = IO_File ) ;

    ! Non-Zeros of Stiffness file (.NNZStiff)
    UnFile = UnInptStiff ;
    Close ( Unit = UnFile, STATUS = 'KEEP', ERR = 1002, IOSTAT = IO_File ) ;

    ! Non-Zeros of Damping file (.NNZDamp)
    UnFile = UnInptDamp ;
    Close ( Unit = UnFile, STATUS = 'KEEP', ERR = 1002, IOSTAT = IO_File ) ;

    ! Non-Zeros of Damping file (.NNZMass)
    UnFile = UnInptMass ;
    Close ( Unit = UnFile, STATUS = 'KEEP', ERR = 1002, IOSTAT = IO_File ) ;

    ! Application Ordering file (.App)
    UnFile = UnInptApp ;
    Close ( Unit = UnFile, STATUS = 'KEEP', ERR = 1002, IOSTAT = IO_File ) ;

  End If ;

