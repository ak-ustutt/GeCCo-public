      subroutine make_doubly_connected(flist,
     &     idxtbar,idxham,idxtop,op_info)
*----------------------------------------------------------------------*
*     ensure that in all terms TBAR is either connected to H or twice
*     to T
*----------------------------------------------------------------------*

      implicit none

      integer, parameter ::
     &     ntest =  000

      include 'stdunit.h'
      include 'opdim.h'
      include 'def_contraction.h'
      include 'mdef_operator_info.h'
      include 'ifc_operators.h'
      include 'def_orbinf.h'
      include 'def_formula_item.h'
      include 'def_formula.h'
      include 'def_del_list.h'
      include 'par_opnames_gen.h'

      type(formula_item), intent(inout), target::
     &     flist
      type(operator_info), intent(in) ::
     &     op_info
      integer, intent(in) ::
     &     idxtbar, idxtop, idxham

      logical ::
     &     ok, delete
      integer ::
     &     nvtx, narc, iarc, ivtx, jvtx, 
     &     idx_op, jdx_op,
     &     ntop, nt1, nham, ntbar
      character*64 ::
     &     op_name
      logical ::
     &     dagger

      type(cntr_vtx), pointer ::
     &     vertex(:)
      type(cntr_arc), pointer ::
     &     arc(:)
      type(formula_item), pointer ::
     &     form_pnt, form_pnt_next

c      integer, external ::
c     &     idx_oplist2

      if (ntest.ge.100) then
        call write_title(lulog,wst_dbg_subr,'make_doubly_connected')
      endif

      form_pnt => flist
      do 
        form_pnt_next => form_pnt%next
        ! Locate actual formula items.
        select case(form_pnt%command)
        case(command_end_of_formula)
          if(ntest.ge.1000) write(lulog,*) '[END]'
        case(command_set_target_init)
          if(ntest.ge.1000) write(lulog,*) '[INIT_TARGET]'
        case(command_add_contribution)

          nvtx = form_pnt%contr%nvtx
          narc = form_pnt%contr%narc
          vertex => form_pnt%contr%vertex
          arc => form_pnt%contr%arc

          ok = .true.
          ! look for tbar vertex
          do ivtx = 1, nvtx
            idx_op  = vertex(ivtx)%idx_op
c            iblk_op = vertex(ivtx)%iblk_op
            if (idx_op.eq.idxtbar) then
              ! look at arcs: find either two T's or
              ! the hamiltonian
              nham = 0
              ntop = 0
              do iarc = 1, narc
                jvtx = arc(iarc)%link(1)
                if (ivtx.ne.jvtx) cycle
                jvtx = arc(iarc)%link(2)
                jdx_op = vertex(jvtx)%idx_op
                if (jdx_op.eq.idxtop) ntop = ntop+1
                if (jdx_op.eq.idxham) nham = nham+1
              end do
              ok = ok.and.(ntop.ge.2.or.nham.eq.1)
            end if
            if (.not.ok) exit
          end do

          delete = .not.ok

          if (delete) then
            ! Print the deleted contraction.
            if(ntest.ge.1000)then
              write(lulog,*)'Deleted formula item:'
              call prt_contr2(lulog,form_pnt%contr,op_info)
            endif

            ! Delete the node.
            call delete_fl_node(form_pnt)
            deallocate(form_pnt)
          end if

        case default
          write(lulog,*)'command = ',form_pnt%command
          call quit(1,'make_doubly_connected','command undefined here')
        end select

        ! Exit or move to the next item.
        if(.not.associated(form_pnt_next))exit
        form_pnt => form_pnt_next

      enddo

      return
      end
