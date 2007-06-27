*----------------------------------------------------------------------*
      subroutine wrt_occ2(luout,iocc,jocc)
*----------------------------------------------------------------------*

      implicit none
      include 'opdim.h'

      integer, intent(in) ::
     &     luout,
     &     iocc(ngastp,2), jocc(ngastp,2)
      character ::
     &     fmt*3

      write(fmt,'(i1,a)') ngastp,'i3'
      write(luout,'(x,2("/",'//fmt//',"'//char(92)//'"))')
     &     iocc(1:ngastp,1),jocc(1:ngastp,1)
      write(luout,'(x,2("'//char(92)//'",'//fmt//',"/"))')
     &     iocc(1:ngastp,2),jocc(1:ngastp,2)

      return
      end
