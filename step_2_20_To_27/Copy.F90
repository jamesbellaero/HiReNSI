Module Copyfile ;

!Implicit none ;

Contains;


Subroutine Copy_File (file_name, file_name_new) ;
! copies a file file_name to file_name_new
! file_name and file_name_new must include the path information and may include wildcard characters

USE ifport 


implicit character*300 (f)
character*300 fnam
logical*4 logical_result


len1 = len_trim(file_name); len2 = len_trim(file_name_new);
!fnam = 'copy/y ' //file_name(1:len1) //' '//file_name_new(1:len2)
fnam = 'cp ' //file_name(1:len1) //' '//file_name_new(1:len2)

l = len_trim(fnam) ;
logical_result = systemqq(fnam(1:l)) ;

Return ;
End Subroutine Copy_File ;

End Module Copyfile ;