      integer function idx_int_graph(idxstr,nel,iy,igamorb,ngam)

      implicit none

      include 'multd2h.h'

      integer, intent(in) ::
     &     nel, ngam, idxstr(nel), iy(nel,ngam,*), igamorb(*)

      integer ::
     &     ii
      integer ::
     &     igamstr(nel)

c dbg
c      print *,'idxstr: ',idxstr
c dbg
      igamstr(1) = igamorb(idxstr(1))
      do ii = 2, nel
        igamstr(ii) = multd2h(igamstr(ii-1),igamorb(idxstr(ii)))
      end do

      idx_int_graph = 1
      do ii = 1, nel
        idx_int_graph = idx_int_graph + iy(ii,igamstr(ii),idxstr(ii))
c dbg
c        print *,'ii, g, o, iy, sum: ',ii, igamstr(ii),idxstr(ii),
c     &       iy(ii,igamstr(ii),idxstr(ii)), idx_int_graph
c dbg
      end do

      return
      end
