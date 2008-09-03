*----------------------------------------------------------------------*
      subroutine print_form_list_p(luout,form_head,op_info)
*----------------------------------------------------------------------*
*     print formula on linked list to unit luout
*     version that prints contractions in topo matrix form
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'mdef_operator_info.h'
      include 'def_contraction.h'
      include 'def_formula_item.h'

      integer, intent(in) ::
     &     luout
      type(formula_item), intent(in), target ::
     &     form_head
      type(operator_info), intent(in) ::
     &     op_info

      integer ::
     &     idx, nvtx, nj
      integer(8), pointer ::
     &     xlines(:,:), topo(:,:), vtx(:)

      type(formula_item), pointer ::
     &     form_ptr

      integer, external ::
     &     njres_contr

      idx = 0
      form_ptr => form_head
      do
        select case(form_ptr%command)
        case(command_end_of_formula)
          write(luout,*) '[END]'
        case(command_set_target_init)
          write(luout,*) '[INIT TARGET]',form_ptr%target
        case(command_set_target_update)
          write(luout,*) '[SET TARGET]',form_ptr%target
        case(command_new_intermediate)
          write(luout,*) '[NEW INTERMEDIATE]',form_ptr%target
        case(command_del_intermediate)
          write(luout,*) '[DELETE INTERMEDIATE]',form_ptr%target
        case(command_add_contribution)
          idx = idx+1
          write(luout,*) '[ADD]',form_ptr%target,'( term #',idx,')'
          write(luout,*) 'factor = ',form_ptr%contr%fac
          nj = njres_contr(form_ptr%contr)
          nvtx = form_ptr%contr%nvtx
          allocate(vtx(nvtx),xlines(nvtx,nj),topo(nvtx,nvtx))
          call pack_contr(form_ptr%contr%svertex,
     &                    vtx,topo,xlines,form_ptr%contr,nj)
          call prt_contr_p(luout,form_ptr%contr%svertex,
     &         vtx,topo,xlines,nvtx,nj)
          deallocate(vtx,xlines,topo)
        case(command_symmetrise)
          write(luout,*) '[SYMMETRISE TARGET]',form_ptr%target
        end select

        if (.not.associated(form_ptr%next)) exit
        form_ptr => form_ptr%next

      end do

      return
      end

