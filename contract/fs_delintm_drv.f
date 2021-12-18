*----------------------------------------------------------------------*
      subroutine fs_delintm_drv(fl_item,op_info)
*----------------------------------------------------------------------*
*     driver for [DEL INTERMEDIATE] operation
*----------------------------------------------------------------------*

      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'multd2h.h'
      include 'mdef_operator_info.h'
      include 'def_contraction.h'
      include 'def_formula_item.h'
      include 'def_orbinf.h'

      integer, parameter ::
     &     ntest = 00

      type(formula_item), intent(in) ::
     &     fl_item
      type(operator_info), intent(inout) ::
     &     op_info

      integer ::
     &     idx
      type(operator), pointer ::
     &     interm

      integer, external ::
     &     idx_oplist2

      if (ntest.ge.100) then
        call write_title(lulog,wst_dbg_subr,'fs_delintm_drv')
        write(lulog,'(2x,a)') trim(fl_item%interm%name)
        call print_op_occ(lulog,fl_item%interm)
      end if

      idx = idx_oplist2(trim(fl_item%label),op_info)
      interm => op_info%op_arr(idx)%op

      ! remove the associated list
      call del_me_list(interm%assoc_list,op_info)
      ! remove the operator
      call del_operator(interm%name,op_info)

      return
      end

