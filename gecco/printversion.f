*----------------------------------------------------------------------*
      subroutine printversion()
*----------------------------------------------------------------------*
      implicit none
*----------------------------------------------------------------------*
* date.h is automatically created by date.sh:      
*----------------------------------------------------------------------*
      include 'date.h'
      include 'stdunit.h'

      write (luout,'(/,x,2a,/)') vers
      write (luout,'(/,x,2a)') 'compiled with: ',cmp(1:len_trim(cmp))
      write (luout,'(x,2a)') 'level 1 optimization: ',trim(opt1)
      write (luout,'(x,2a)') 'level 2 optimization: ',trim(opt2)
      write (luout,'(x,2a/)')'level 3 optimization: ',trim(opt3)

      end
*----------------------------------------------------------------------*

