*----------------------------------------------------------------------*
      subroutine list_operators(lulog,op_info)
*----------------------------------------------------------------------*
*     print list of operators
*----------------------------------------------------------------------*

      implicit none

      include 'opdim.h'
      include 'mdef_operator_info.h'

      type(operator_info), intent(inout), target ::
     &     op_info
      integer, intent(in) ::
     &     lulog

      type(operator_list), pointer ::
     &     list_pnt

      integer ::
     &     idx

      idx = 0
      list_pnt => op_info%op_list
      ! advance to end of operator list:
      do 
        idx = idx+1
        write(lulog,'(3x,i4,2x,a8,2x,i4)') idx,list_pnt%op%name,
     &       list_pnt%op%n_occ_cls

        if (.not.associated(list_pnt%next)) exit
        list_pnt => list_pnt%next
      end do

      return
      end
