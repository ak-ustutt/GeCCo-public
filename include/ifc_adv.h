      interface

       subroutine mel_adv_state(mel,max_state,last_state)
       import :: me_list
       implicit none
       type(me_list), intent(inout) ::
     &      mel
       integer, intent(in) ::
     &      max_state
       logical, intent(out), optional ::
     &      last_state
       end subroutine

       subroutine mel_adv_state_wrap(label_mel,nmel,max_state,op_info,
     &     last_state)
       import :: operator_info
       implicit none
       integer, intent(in) ::
     &      max_state, nmel
       character(len=*), intent(in) ::
     &      label_mel(nmel)
       type(operator_info), intent(inout) ::
     &      op_info
       logical, intent(out), optional ::
     &      last_state
       end subroutine

       subroutine op_adv_state(label_op,nop,max_state,op_info,use_1,
     &     last_state)
       import :: operator_info
       implicit none
       integer, intent(in) ::
     &      max_state, nop
       character(*), intent(in) ::
     &      label_op(nop)
       type(operator_info), intent(inout) ::
     &      op_info
       logical, intent(in) ::
     &      use_1
       logical, intent(out), optional ::
     &      last_state
       end subroutine

      end interface
