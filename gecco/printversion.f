*----------------------------------------------------------------------*
      subroutine printversion()
*----------------------------------------------------------------------*
      implicit none
*----------------------------------------------------------------------*
* date.h is automatically created by date.sh:      
*----------------------------------------------------------------------*
      include "date.h"

      write (6,'(/,x,2a,/)') vers
      write (6,'(/,x,2a)') 'compiled with: ',cmp(1:len_trim(cmp))
      write (6,'(x,2a)') 'level 1 optimization: ',opt1(1:len_trim(opt1))
      write (6,'(x,2a)') 'level 2 optimization: ',opt2(1:len_trim(opt2))
      write (6,'(x,2a/)')'level 3 optimization: ',opt3(1:len_trim(opt3))

      end
*----------------------------------------------------------------------*

