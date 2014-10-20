*----------------------------------------------------------------------*
      subroutine mel_adv_state(mel,max_state,last_state)
*----------------------------------------------------------------------*
*     advance the state of mel
*     if the state goes beyond max_state, restore to the first state
*     case that last_state is set to T.
*
*     yuri 2014
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'mdef_operator_info.h'

      integer, parameter ::
     &     ntest = 2

      type(me_list), intent(inout) ::
     &     mel
      integer, intent(in) ::
     &     max_state
      logical, intent(out), optional ::
     &     last_state

      integer ::
     &     i_state

      integer, external ::
     &     get_mel_record

      i_state = get_mel_record(mel)

      if (i_state.ge.max_state) then
       i_state = 1
       if(ntest.GT.1) write(lulog,FMT='(" All states processed.")')
       if(present(last_state)) last_state = .true.
      else
       i_state = i_state + 1
       if(present(last_state)) last_state = .false.
       if(ntest.GT.1) write(lulog,FMT='(" Next state: ",i0)') i_state
      end if
      
      call switch_mel_record(mel,i_state)

      return
      end
