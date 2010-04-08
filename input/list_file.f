      subroutine list_file(luout,luin)
 
      implicit none

      integer, intent(in) :: luout, luin

      character(len=1024) :: line

      write(luout,'(x,"+",77("-"),"+")')
      rewind luin
      file_loop: do
        read(luin,'(a)',end=100,err=200) line
        write(luout,*) trim(line)
      end do file_loop

200   write(luout,*) 'Error while reading ... try to continue ...'
100   continue

      write(luout,'(x,"+",77("-"),"+")')

      rewind luin
    
      end

