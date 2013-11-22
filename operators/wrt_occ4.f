*----------------------------------------------------------------------*
      subroutine wrt_occ4(lulog,iocc,jocc,kocc,locc)
*----------------------------------------------------------------------*

      implicit none
      include 'opdim.h'

      integer, intent(in) ::
     &     lulog,
     &     iocc(ngastp,2), jocc(ngastp,2),
     &     kocc(ngastp,2), locc(ngastp,2)
      character ::
     &     fmt*3

      write(fmt,'(i1,a)') ngastp,'i3'
      write(lulog,'(x,4("/",'//fmt//',"'//char(92)//'"))')
     &     iocc(1:ngastp,1),jocc(1:ngastp,1),
     &     kocc(1:ngastp,1),locc(1:ngastp,1)
      write(lulog,'(x,4("'//char(92)//'",'//fmt//',"/"))')
     &     iocc(1:ngastp,2),jocc(1:ngastp,2),
     &     kocc(1:ngastp,2),locc(1:ngastp,2)

      return
      end
