*----------------------------------------------------------------------*
      subroutine printversion(lu)
*----------------------------------------------------------------------*
      implicit none
*----------------------------------------------------------------------*
* date.h is automatically created by date.sh:      
*----------------------------------------------------------------------*
      include 'date.h'
      integer, intent(in) :: lu

      write (lu,'(/,x,2a,/)') vers
      write (lu,'(/,x,2a)') 'compiled with: ',cmp(1:len_trim(cmp))
      write (lu,'(x,2a)') 'level 1 optimization: ',trim(opt1)
      write (lu,'(x,2a)') 'level 2 optimization: ',trim(opt2)
      write (lu,'(x,2a/)')'level 3 optimization: ',trim(opt3)

      end
*----------------------------------------------------------------------*

