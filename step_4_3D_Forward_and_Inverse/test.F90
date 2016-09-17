program main

implicit double precision (a-h, o-z)

do i = 1, 40
 
   Iter_File_10 = Mod ( i , 10 )
   Iter_File_20 = Mod ( i , 20 )

!      write(*,*) i, iter_file_10, iter_file_20

   If ( ( Iter_File_10 == 0 ) .AND. ( Iter_File_20 /= 0 ) ) Then
!      write(*,*) i
   Else If ( Iter_File_20 == 0 ) Then
!      write(*,*) i
   End If

end do

!--------------------------------------------------------------------------------------------------
K_LBFGS = 6
M_LBFGS = 5


Do I = K_LBFGS, 1, -1
!   Write(*,*) 'I = ', I
End Do
!
!
If ( (mod ( K_LBFGS , M_LBFGS )) /= 0 ) Then
   N_store = mod ( K_LBFGS , M_LBFGS )
Else
   N_store = M_LBFGS
End If
!
Write(*,*) 'N_store = ', N_store


end
