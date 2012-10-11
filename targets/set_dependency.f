*----------------------------------------------------------------------*
      subroutine set_dependency(name_target,name_dependency,tgt_info)
*----------------------------------------------------------------------*
*     add entry "name_dependency" to depend_on array of target 
*     referenced by "name_target" (must exist on tgt_info)
*     i.e. name_target depends on name_dependency.
*----------------------------------------------------------------------*
      implicit none

      include 'mdef_target_info.h'

      type(target_info), intent(inout), target ::
     &     tgt_info
      character(len=*), intent(in) ::
     &     name_target, name_dependency

      integer ::
     &     idx, ii, jj
      type(target), pointer ::
     &     tgt
      character*(len_target_name), pointer ::
     &     new_depends_on(:)
      integer, external ::
     &     idx_target

      idx = idx_target(name_target,tgt_info)

      if (idx.le.0)
     &     call quit(1,'set_dependency',
     &     'target not (yet) defined: '//trim(name_target))
      if (len_trim(name_dependency).gt.len_target_name)
     &     call quit(1,'set_dependency',
     &     'name dependency too long: '//trim(name_dependency))

      tgt => tgt_info%array(idx)%tgt

      allocate(new_depends_on(tgt%n_depends_on+1))
      if (tgt%n_depends_on.gt.0) then
        new_depends_on(1:tgt%n_depends_on) =
     &       tgt%depends_on(1:tgt%n_depends_on)
        deallocate (tgt%depends_on)
      end if
      !new_depends_on(tgt%n_depends_on+1)(1:len_target_name) = ' '
      !do ii = 1, len_target_name
      !  new_depends_on(tgt%n_depends_on+1)(ii:ii)=' '
      !end do
      !new_depends_on(tgt%n_depends_on+1) = trim(name_dependency)
      tgt%depends_on => new_depends_on
      tgt%depends_on(tgt%n_depends_on+1)(1:len_target_name) = ' '
      tgt%depends_on(tgt%n_depends_on+1) = trim(name_dependency)      

      tgt%n_depends_on = tgt%n_depends_on+1

      return
      end
