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


      end interface
