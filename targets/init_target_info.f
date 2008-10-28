*----------------------------------------------------------------------*
      subroutine init_target_info(tgt_info)
*----------------------------------------------------------------------*

      implicit none

      include 'mdef_target_info.h'
      include 'event_counter.h'

      type(target_info), intent(inout) ::
     &     tgt_info

      ! no targets yet
      tgt_info%ntargets = 0
      ! ... and no action
      tgt_info%nactions = 0
      ! initialize target list
      allocate(tgt_info%list)
      nullify(tgt_info%list%tgt)
      nullify(tgt_info%list%prev)
      nullify(tgt_info%list%next)
      ! initialize action list
      allocate(tgt_info%act_list)
      nullify(tgt_info%act_list%act)
      nullify(tgt_info%act_list%prev)
      nullify(tgt_info%act_list%next)
      ! initialize pointer arrays
      nullify(tgt_info%array)
      nullify(tgt_info%act_array)

      return
      end
