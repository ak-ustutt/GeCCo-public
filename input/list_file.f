      subroutine list_file(lulog,luin)
 
      implicit none

      integer, intent(in) :: lulog, luin

      character(len=1024) :: line

      write(lulog,'(x,"+",77("-"),"+")')
      rewind luin
      file_loop: do
        read(luin,'(a)',end=100,err=200) line
        write(lulog,*) trim(line)
      end do file_loop

200   write(lulog,*) 'Error while reading ... try to continue ...'
100   continue

      write(lulog,'(x,"+",77("-"),"+")')

      rewind luin
    
      end

