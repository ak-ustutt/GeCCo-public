*----------------------------------------------------------------------*
      subroutine set_multistate_operator(tgt_info,n_states,op,ini_state)
*----------------------------------------------------------------------*
*     clone operator op for each state, to op//state_label, using the
*     same target of op (assuming the same name)
*
*     yuri 2014
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'mdef_target_info.h'
      include 'ifc_targets.h'
      include 'par_actions.h'

      type(target_info), intent(inout) ::
     &     tgt_info
      integer, intent(in) ::
     &     n_states
      character(*), intent(in) ::
     &     op
      integer, intent(in) ::
     &     ini_state

      character(len=3) ::
     &     c_st
      integer ::
     &     i_state

      character(len_target_name), external ::
     &     state_label

      do i_state=ini_state,n_states,1
       c_st = state_label(i_state, .true.)
       call set_rule2(trim(op),CLONE_OP,tgt_info)
       call set_arg(trim(op),CLONE_OP,'LABEL',1,tgt_info,
     &      val_label=[trim(op)//trim(c_st)])
       call set_arg(trim(op),CLONE_OP,'TEMPLATE',1,tgt_info,
     &      val_label=[trim(op)])
      enddo
      
      end subroutine set_multistate_operator
