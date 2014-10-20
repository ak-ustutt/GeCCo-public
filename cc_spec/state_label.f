*----------------------------------------------------------------------*
      function state_label(i_state, use_1)
*----------------------------------------------------------------------*
*     return a simple label to use as a sufix for state specific 
*     eneities: "_<i_state>"
*     if use_1 is false, don't use a label ("") for 1st state
*----------------------------------------------------------------------*

      implicit none
      include 'def_target.h'

      character(len_target_name) ::
     &     state_label

      integer, intent(in) ::
     &     i_state
      logical, intent(in) ::
     &     use_1
      
      if (i_state.eq.1.and..not.use_1) then
       state_label = ''
      else
       write(state_label, '("_",i0)') i_state
      end if
      
      return
      end
