
*----------------------------------------------------------------------*
      subroutine wrt_occ(luout,iocc)
*----------------------------------------------------------------------*

      implicit none
      include 'opdim.h'

      integer, intent(in) ::
     &     luout, iocc(ngastp,2)
      character ::
     &     fmt*3

      write(fmt,'(i1,a)') ngastp,'i3'
      write(luout,'(x,"/",'//fmt//',"\")') iocc(1:ngastp,1)
      write(luout,'(x,"\",'//fmt//',"/")') iocc(1:ngastp,2)

      return
      end
