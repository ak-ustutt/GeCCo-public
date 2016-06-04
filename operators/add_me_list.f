*----------------------------------------------------------------------*
      subroutine add_me_list(label,op_info)
*----------------------------------------------------------------------*
*     allocate a new slot for an me_list structure in op_info
*     the list is extended and the array is updated
*----------------------------------------------------------------------*

      implicit none

      include 'opdim.h'
      include 'mdef_operator_info.h'

      type(operator_info), intent(inout), target ::
     &     op_info
      character ::
     &     label*(*)

      type(me_list_list), pointer ::
     &     list_pnt

      integer ::
     &     imels

      list_pnt => op_info%mel_list
      ! advance to end of me_list list:
      do while (associated(list_pnt%next))
        list_pnt => list_pnt%next
      end do
      ! is last entry already in use?
      if (associated(list_pnt%mel)) then
        allocate(list_pnt%next)
        list_pnt%next%prev => list_pnt
        list_pnt => list_pnt%next
        list_pnt%next => null()
      end if
      allocate (list_pnt%mel)

      ! add user-supplied label (one free space needed)
      if (len_trim(label).ge.mxlen_melabel)
     &   call quit(1,'add_me_list','name too long: "'
     &     //trim(label)//'"')
      list_pnt%mel%label(1:mxlen_melabel) = ' '
      list_pnt%mel%label = label

      ! init all pointers
      list_pnt%mel%op => null()
      list_pnt%mel%fhand => null()
      list_pnt%mel%len_op_occ => null()
      list_pnt%mel%off_op_occ => null()
      list_pnt%mel%len_op_gmo => null()
      list_pnt%mel%off_op_gmo => null()
      list_pnt%mel%len_op_gmox => null()
      list_pnt%mel%off_op_gmox => null()
      list_pnt%mel%idx_graph => null()

      ! increment counter
      op_info%nmels = op_info%nmels+1

      ! update me_list array
      call update_mel_arr(op_info)

      return
      end
