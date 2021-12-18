*----------------------------------------------------------------------*
      subroutine reset_target(name,tgt_info)
*----------------------------------------------------------------------*
*     reset all entries associated with target
*----------------------------------------------------------------------*

      implicit none

      include 'mdef_target_info.h'

      type(target_info), intent(inout), target ::
     &     tgt_info
      character(len=*), intent(in) ::
     &     name

      type(target), pointer ::
     &     tgt

      integer ::
     &     itgts, idx

      integer, external ::
     &     idx_target

      call quit(1,'reset_target',
     &     'not yet tested (should work -> feel free to try out)!' )

      idx = idx_target(name,tgt_info)
      if (idx.ge.0)
     &     call quit(1,'reset_target',
     &     'target does not exist: '//trim(name) )

      tgt => tgt_info%array(idx)%tgt

      ! reset last_mod
      tgt%last_mod = -1

      ! set other defaults
      tgt%n_joined_with = 0
      tgt%n_depends_on = 0
      if (associated(tgt%joined_with)) deallocate(tgt%joined_with)
      if (associated(tgt%depends_on))  deallocate(tgt%depends_on)
      if (associated(tgt%idx_joined_with))
     &                                 deallocate(tgt%idx_joined_with)
      if (associated(tgt%idx_depends_on))
     &                                 deallocate(tgt%idx_depends_on)
      tgt%joined_with => null()
      tgt%depends_on => null()
      tgt%idx_joined_with => null()
      tgt%idx_depends_on => null()

      tgt%n_rules = 0
      if (associated(tgt%rules)) deallocate(tgt%rules)
      tgt%rules => null()

      ! update target array
      call update_tgt_info(tgt_info)

      return
      end
