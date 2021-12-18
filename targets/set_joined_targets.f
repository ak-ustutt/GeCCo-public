*----------------------------------------------------------------------*
      subroutine set_joined_targets(name_target,name_joined,tgt_info)
*----------------------------------------------------------------------*
*     add entry "name_joined" to joined_with array of target 
*     referenced by "name_target" (must exist on tgt_info)
*----------------------------------------------------------------------*
      implicit none

      include 'mdef_target_info.h'

      type(target_info), intent(inout), target ::
     &     tgt_info
      character, intent(in) ::
     &     name_target*(*), name_joined*(*)

      integer ::
     &     idx
      type(target), pointer ::
     &     tgt
      character*(len_target_name), pointer ::
     &     new_joined_with(:)
      integer, external ::
     &     idx_target

      idx = idx_target(name_target,tgt_info)

      if (idx.le.0)
     &     call quit(1,'set_joined_targets',
     &     'target not (yet) defined: '//trim(name_target))
      if (len_trim(name_joined).gt.len_target_name)
     &     call quit(1,'set_joined_targets',
     &     'name for joined targets too long: '//trim(name_joined))

      tgt => tgt_info%array(idx)%tgt

      allocate(new_joined_with(tgt%n_joined_with+1))
      if (tgt%n_joined_with.gt.0) then
        new_joined_with(1:tgt%n_joined_with) =
     &       tgt%joined_with(1:tgt%n_joined_with)
        deallocate (tgt%joined_with)
      end if
      new_joined_with(tgt%n_joined_with+1)(1:len_target_name) = ' '
      new_joined_with(tgt%n_joined_with+1) = trim(name_joined)
      tgt%joined_with => new_joined_with
      tgt%n_joined_with = tgt%n_joined_with+1

      return
      end
