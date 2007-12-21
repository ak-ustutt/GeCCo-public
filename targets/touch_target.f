*----------------------------------------------------------------------*
      subroutine touch_target(idx,touch_joined,tgt_info)
*----------------------------------------------------------------------*
*     mark target #idx as updated
*     i.e. increment counter in tgt_info and set the last_mod
*     information accordingly
*     touch_joined: if .true. touch the joined targets as well
*----------------------------------------------------------------------*
      implicit none

      include 'mdef_target_info.h'
      include 'event_counter.h'

      integer, intent(in) ::
     &     idx
      logical, intent(in) ::
     &     touch_joined
      type(target_info), intent(inout) ::
     &     tgt_info
      
      integer ::
     &     jdx, kdx

      event_time = event_time+1
      if (event_time.lt.0)
     &     call quit(1,'touch_target','event counter overflow')

      tgt_info%last_mod(idx) = event_time
      tgt_info%array(idx)%tgt%last_mod = event_time

      if (touch_joined) then
        do jdx = 1, tgt_info%array(idx)%tgt%n_joined_with
          kdx = tgt_info%array(idx)%tgt%idx_joined_with(jdx)
          
          tgt_info%last_mod(kdx) = event_time
          tgt_info%array(kdx)%tgt%last_mod = event_time

        end do
      end if

      return
      end
