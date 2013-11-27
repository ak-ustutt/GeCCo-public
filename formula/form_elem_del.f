      subroutine form_elem_del(flist,del_list,op_info)
*----------------------------------------------------------------------*
*     Routine used to delete the elements described in del_list from the
*     formula, flist.
*     Usually called by truncate_form.f.
*     GWR February 2008
*----------------------------------------------------------------------*

      implicit none

      integer, parameter ::
     &     ntest = 000

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
      type(del_cond_list), intent(in) ::
     &     del_list
      type(operator_info), intent(in) ::
     &     op_info

      logical ::
     &     del_item
      integer ::
     &     or_dim, ordx, and_dim, anddx, nvtx, idxop, idx, idx_op,
     &     iblk_op, op_occur, nelec
      character*64 ::
     &     op_name
      logical ::
     &     dagger

      type(formula_item), pointer ::
     &     form_pnt
      type(del_cond), pointer ::
     &     temp_del

      integer, external ::
     &     idx_oplist2

      if(ntest.ge.100)then
        write(lulog,*) '------------------'
        write(lulog,*) ' Element deletion '
        write(lulog,*) '------------------'
      endif

      form_pnt => flist
      or_dim = del_list%or_dim
      do 
        ! Locate actual formula items.
        select case(form_pnt%command)
        case(command_end_of_formula)
          if(ntest.ge.1000) write(lulog,*) '[END]'
        case(command_set_target_init)
          if(ntest.ge.1000) write(lulog,*) '[INIT_TARGET]'
        case(command_add_contribution)
c          if(ntest.ge.1000)then
c            write(lulog,*)'[ADD]:'
c            call prt_contr2(lulog,form_pnt%contr,op_info)
c          endif

          nvtx = form_pnt%contr%nvtx
          or_loop: do ordx = 1, or_dim
            del_item = .true.
            and_dim = del_list%del_cond_item(ordx)%and_dim
            and_loop: do anddx = 1, and_dim

              temp_del =>
     &           del_list%del_cond_item(ordx)%del_cond_arr(anddx)

              op_name = temp_del%op_name
              idxop = idx_oplist2(trim(op_name),op_info)
              if (idxop.lt.0) then
                del_item = .false.
                exit and_loop
              end if

              dagger = temp_del%transposed

              op_occur = 0 

              vert_loop: do idx = 1, nvtx
                ! Is this vertex represented by the correct operator?
                if(form_pnt%contr%vertex(idx)%idx_op.eq.idxop .and.
     &             (form_pnt%contr%vertex(idx)%dagger.eqv.dagger))then
                  ! Properties of the vertex.
                  iblk_op = form_pnt%contr%vertex(idx)%iblk_op
                  nelec = op_info%op_arr(idxop)%op%ica_occ(1,iblk_op)

                  ! Check to see if this block has the correct occupancy.
                  if(nelec.lt.temp_del%part_num_restr(1).and.
     &                 temp_del%part_num_restr(1).ne.-1)cycle vert_loop
                  if(nelec.gt.temp_del%part_num_restr(2).and.
     &                 temp_del%part_num_restr(2).ne.-1)cycle vert_loop

                  ! Only increment the operator count if the occupancy 
                  ! conditions are met. (This allows the same operator
                  ! but with different occupancies to be treated)
                  op_occur = op_occur+1

                endif
              enddo vert_loop

              if(.not.(op_occur.ge.temp_del%num_op_restr(1)))then
                del_item = .false.
                exit and_loop
              endif
              if(.not.(op_occur.le.temp_del%num_op_restr(2).or.
     &             temp_del%num_op_restr(2).eq.-1))then
                del_item = .false.
                exit and_loop
              endif

            enddo and_loop

            if(del_item)then
              ! Print the deleted contraction.
              if(ntest.ge.1000)then
                write(lulog,*)'Deleted formula item:'
                call prt_contr2(lulog,form_pnt%contr,op_info)
              endif

              ! Delete the node.
              call delete_fl_node(form_pnt)
              if(associated(form_pnt%contr))then
                call dealloc_contr(form_pnt%contr)
                deallocate(form_pnt%contr)
              endif

              exit or_loop
            endif
          enddo or_loop

        case default
          write(lulog,*)'command = ',form_pnt%command
          call quit(1,'delete_non_fact','command undefined here')
        end select

        ! Exit or move to the next item.
        if(.not.associated(form_pnt%next))exit
        form_pnt => form_pnt%next

      enddo


      return
      end
      
      
