*----------------------------------------------------------------------*
      subroutine del_operator(id_op,op_info)
*----------------------------------------------------------------------*
*     delete info for operator with ID id_op
*----------------------------------------------------------------------*

      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'mdef_operator_info.h'

      type(operator_info), intent(inout), target ::
     &     op_info
      integer, intent(in) ::
     &     id_op

      type(operator_list), pointer ::
     &     list_pnt
      type(operator), pointer ::
     &     op

      list_pnt => op_info%op_list
      ! advance to respective entry on operator list:
      do while (list_pnt%op%id.ne.id_op.and.associated(list_pnt%next))
        list_pnt => list_pnt%next
      end do
      if (list_pnt%op%id.ne.id_op) then
        write(luout,*) 'ID = ',id_op
        call quit(1,'del_operator','unknown ID')
      end if

      if (associated(list_pnt%prev)) list_pnt%prev%next => list_pnt%next
      if (associated(list_pnt%next)) list_pnt%next%prev => list_pnt%prev
      op => list_pnt%op
      if (associated(op%ihpvca_occ)) deallocate(op%ihpvca_occ)
      if (associated(op%ica_occ)) deallocate(op%ica_occ)
      if (associated(op%igasca_restr)) deallocate(op%igasca_restr)

      ! decrement counter
      op_info%nops = op_info%nops-1

      ! update operator array
      call update_op_arr(op_info)

      return
      end
