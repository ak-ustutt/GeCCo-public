      interface

c deactivated as currently problems with intel compiler occur
c      subroutine add_action(act_list,nactions,
c     &     action_type,nop_in,nop_out,
c     &     idxopdef_in,idxopdef_out,
c     &     idxopfile_in,idxopfile_out,
c     &     nform,
c     &     idx_formula)
c      implicit none
c      include 'stdunit.h'
c      include 'def_filinf.h'
c      include 'def_action.h'
c      include 'def_action_list.h'
c      type(action_list), intent(inout), target ::
c     &     act_list
c      integer, intent(inout) ::
c     &     nactions
c      integer, intent(in) ::
c     &     action_type
c      integer, intent(in), optional ::
c     &     nop_in,nop_out,
c     &     idxopdef_in(*),idxopdef_out(*),
c     &     idxopfile_in(2,*),idxopfile_out(2,*),
c     &     nform,
c     &     idx_formula(*)
c      end subroutine

      subroutine transpose_formula(form_head,op_info,multi)
      import :: formula_item
      import :: operator_info
      implicit none

      type(formula_item), intent(in), target ::
     &     form_head
      type(operator_info), intent(in) ::
     &     op_info
      logical, intent(in), optional ::
     &     multi

      end subroutine

      subroutine transpose_contr(contr,op_info,multi)
      import :: contraction
      import :: operator_info
      implicit none

      type(contraction) ::
     &     contr
      type(operator_info) ::
     &     op_info
      logical, intent(in), optional ::
     &     multi
      end subroutine

      end interface
