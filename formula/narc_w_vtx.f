*----------------------------------------------------------------------*
      integer function narc_w_vtx(vtx_list,nlist,contr)
*----------------------------------------------------------------------*
*     return the number of arc which involve vertices in vtx_list(nlist)
*----------------------------------------------------------------------*
      include 'opdim.h'
      include 'def_contraction.h'

      type(contraction), intent(in) ::
     &     contr
      integer, intent(in) ::
     &     nlist, vtx_list(nlist)

      integer ::
     &     narc, iarc, il, ilist
      type(cntr_arc), pointer ::
     &     arc(:)

      narc = contr%narc
      arc => contr%arc
      narc_w_vtx = 0
      arc_loop: do iarc = 1, narc
        do il = 1, 2
          vtx = arc(iarc)%link(il)
          do ilist = 1, nlist
            if (vtx.eq.vtx_list(ilist)) then
              narc_w_vtx = narc_w_vtx+1
              cycle arc_loop
            end if
          end do
        end do
      end do arc_loop

      return
      end
