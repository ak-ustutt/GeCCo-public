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
      ! initialize target list
      allocate(tgt_info%list)
      nullify(tgt_info%list%tgt)
      nullify(tgt_info%list%prev)
      nullify(tgt_info%list%next)
      ! initialize pointer array
      nullify(tgt_info%array)

      return
      end
