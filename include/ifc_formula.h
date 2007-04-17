      interface

      subroutine add_action(act_list,nactions,
     &     action_type,nop_in,nop_out,
     &     idxopdef_in,idxopdef_out,
     &     idxopfile_in,idxopfile_out,
     &     nform,
     &     idx_formula)
      implicit none
      include 'stdunit.h'
      include 'def_filinf.h'
      include 'def_action.h'
      include 'def_action_list.h'
      type(action_list), intent(inout), target ::
     &     act_list
      integer, intent(inout) ::
     &     nactions
      integer, intent(in) ::
     &     action_type
      integer, intent(in), optional ::
     &     nop_in,nop_out,
     &     idxopdef_in(*),idxopdef_out(*),
     &     idxopfile_in(2,*),idxopfile_out(2,*),
     &     nform,
     &     idx_formula(*)
      end subroutine


      end interface
