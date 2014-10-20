*----------------------------------------------------------------------*
      subroutine mel_adv_state_wrap(label_mel,nmel,max_state,op_info,
     &     last_state)
*----------------------------------------------------------------------*
*     advance the state of the MELs in label_mel
*     wraper to mel_adv_state
*
*     yuri 2014
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'mdef_operator_info.h'

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
      end interface

      integer, parameter ::
     &     ntest = 2

      integer, intent(in) ::
     &     max_state, nmel
      character(len=*), intent(in) ::
     &     label_mel(nmel)
      type(operator_info), intent(inout) ::
     &     op_info
      logical, intent(out), optional ::
     &     last_state

      integer ::
     &     i_label, idx
      type(me_list), pointer ::
     &     mel_pnt

      integer, external ::
     &     idx_mel_list

      do i_label = 1,nmel

       ! the label must exist previously ...
       idx = idx_mel_list(label_mel(i_label),op_info)
       if (idx.lt.0)
     &      call quit(1,'mel_adv_state_wrap',
     &      'list does not already exist: "'//
     &      trim(label_mel(i_label)//'"'))
       
       if(ntest.GT.1) write(lulog,*) 'Advancing state in list "',
     &      trim(label_mel(i_label)),'".'
       
       mel_pnt => op_info%mel_arr(idx)%mel
              
       call mel_adv_state(mel_pnt,max_state,last_state)

      end do

      return
      end
