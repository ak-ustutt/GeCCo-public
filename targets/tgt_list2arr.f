*----------------------------§------------------------------------------*
      subroutine tgt_list2arr(tgt_list,tgt_arr,ntgts)
*----------------------------------------------------------------------*
*     set up an array of pointers that point to the entries of the
*     list on tgt_list
*----------------------------------------------------------------------*

      implicit none
      include 'stdunit.h'
      include 'def_target.h'
      include 'def_target_list.h'
      include 'def_target_array.h'

      type(target_list), intent(in), target ::
     &     tgt_list
      integer, intent(in) ::
     &     ntgts
      type(target_array), intent(out) ::
     &     tgt_arr(ntgts)

      type(target_list), pointer ::
     &     current

      integer ::
     &     itgt

      current => tgt_list

      do itgt = 1, ntgts
        if (.not.associated(current%tgt))
     &       call quit(1,'tgt_list2arr','unallocated target on list')
        tgt_arr(itgt)%tgt => current%tgt
        if (itgt.lt.ntgts.and..not.(associated(current%next)))
     &       call quit(1,'tgt_list2arr','unexpected end of list')
        current => current%next
      end do

      return
      end
