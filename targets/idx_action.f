*----------------------------------------------------------------------*
      integer function idx_action(actname,tgt_info)
*----------------------------------------------------------------------*
*     given an action name, search tgt_info and return index of
*     of corresponding prototype action
*----------------------------------------------------------------------*

      implicit none

      include 'stdunit.h'
      include 'mdef_target_info.h'

      integer, parameter ::
     &     ntest = 00

      character, intent(in) ::
     &     actname*(*)
      type(target_info), intent(in) ::
     &     tgt_info

      integer ::
     &     iact

      if (ntest.ge.100) then
        call write_title(lulog,wst_dbg_subr,'this is idx_target')
        write(lulog,*) ' looking for: "',trim(actname),'"'
      end if

      idx_action = -1
      do iact = 1, tgt_info%nactions
        if (trim(actname).eq.
     &      trim(tgt_info%act_array(iact)%act%command)) then
          idx_action = iact
          exit
        end if
      end do

      if (ntest.ge.100) then
        write(lulog,*) 'result: ',idx_action
      end if

      return
      end
