      subroutine set_cnt_vtx_list(vtx_list,contr,arc_list,nlist)

      implicit none

      include 'opdim.h'
      include 'def_contraction.h'

      type(contraction), intent(in) ::
     &     contr
      integer, intent(in) ::
     &     nlist, arc_list(nlist)
      integer, intent(out) ::
     &     vtx_list(2*nlist)

      type(cntr_arc), pointer ::
     &     arc(:)
      integer ::
     &     idx, jdx

      arc => contr%arc
      
      jdx = 0
      do idx = 1, nlist
        vtx_list(jdx+1) = arc(arc_list(idx))%link(1)
        vtx_list(jdx+2) = arc(arc_list(idx))%link(2)
        jdx = jdx+2
      end do

      end
