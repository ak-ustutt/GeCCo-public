*----------------------------------------------------------------------*
      subroutine wrt_occ_n(lulog,iocc,nocc)
*----------------------------------------------------------------------*

      implicit none
      include 'opdim.h'

      integer, intent(in) ::
     &     lulog, nocc,
     &     iocc(ngastp,2,nocc)
      character ::
     &     fmt*3
      
      integer ::
     &     idxst, idxnd

      write(fmt,'(i1,a)') ngastp,'i3'
      idxnd = 0
      do while (idxnd.lt.nocc)
        idxst = idxnd+1
        idxnd = idxnd+min(4,nocc-idxnd)
        write(lulog,'(x,4("/",'//fmt//',"'//char(92)//'"))')
     &     iocc(1:ngastp,1,idxst:idxnd)
        write(lulog,'(x,4("'//char(92)//'",'//fmt//',"/"))')
     &     iocc(1:ngastp,2,idxst:idxnd)
      end do

      return
      end
