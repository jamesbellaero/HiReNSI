
! Simplicity is not a style. It's just the common sense which says yes, take that away, and it's a bit better.

!====================================================================================================================================================

! - End top file wraper for material visualization (only once at the end) ---------------------------------------------------------------------------
  If ( Rank == 0 ) Then ;

    UnFile = UN_OutallWrp_Mat_Vis ;
    Write (Unit = UnFile, FMT = "(A11)", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) '    </Grid>' ;
    Write (Unit = UnFile, FMT = "(A11)", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) '  </Domain>' ;
    Write (Unit = UnFile, FMT = "(A7 )", ADVANCE = 'YES', ASYNCHRONOUS = 'NO',                 IOSTAT = IO_Write, ERR = 1006 ) '</Xdmf>' ;

    Close ( Unit = UnFile, STATUS = 'KEEP', ERR = 1002, IOSTAT = IO_File ) ;

  End If ;

!====================================================================================================================================================
