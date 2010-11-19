*----------------------------------------------------------------------*
      subroutine trunc_formula_dependencies(depend,idx_last,op_info)
*----------------------------------------------------------------------*
*     remove all dependencies for all items after the one with idx_last
*----------------------------------------------------------------------*

      implicit none

      include 'stdunit.h'
      include 'opdim.h'
      include 'def_dependency_info.h'
      include 'mdef_operator_info.h'

      type(dependency_info), intent(inout) ::
     &     depend
      type(operator_info), intent(in) ::
     &     op_info
      integer, intent(in) ::
     &     idx_last

      logical ::
     &     found
      integer ::
     &     idx

      found = .false.
      do idx = 1, depend%ntargets
        if (depend%idxlist(idx).eq.idx_last) then
          found = .true.
          cycle
        else if (.not.found) then
          cycle
        end if
        depend%depends_on_idxlist(1:depend%ndepend,idx) = 0
        ! we call it up to date
        call touch_file_rec(op_info%mel_arr(idx)%mel%fhand)
      end do

      return
      end
