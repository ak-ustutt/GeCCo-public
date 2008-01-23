      subroutine arc_sort(arc,narc,nvtx)

      implicit none

      include 'opdim.h'
      include 'def_contraction.h'

      integer, intent(in) ::
     &     narc, nvtx
      type(cntr_arc), intent(inout) ::
     &     arc(narc)

      type(cntr_arc) ::
     &     arc_sv

      integer ::
     &     idx, ival, jdx, jval
      integer, external ::
     &     int_pack

      do idx = 2, narc
        ival = int_pack(arc(idx)%link,2,nvtx+1)
        arc_sv = arc(idx)
        jdx = idx-1
        do while(jdx.gt.0)
          jval = int_pack(arc(jdx)%link,2,nvtx+1)
          if (jval.le.ival) exit
          arc(jdx+1) = arc(jdx)
          jdx = jdx-1
        end do
        arc(jdx+1) = arc_sv
      end do

      return
      end
