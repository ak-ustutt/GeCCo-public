*----------------------------------------------------------------------*
      subroutine dealloc_contr_list(contr_list)
*----------------------------------------------------------------------*
*     deallocate all subfields in contraction contr
*----------------------------------------------------------------------*
      implicit none
      
      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_contraction_list.h'

      type(contraction_list), intent(inout), target ::
     &     contr_list

      type(contraction_list), pointer ::
     &     current
      
      current => contr_list
      ! go to end of list
      do while(associated(current%next))
        current => current%next
      end do
     
      ! loop backward through list and deallocate
      do
        if (associated(current%contr)) then
          call dealloc_contr(current%contr)
          deallocate(current%contr)
        end if
        if (.not.associated(current%prev)) exit
        current => current%prev
        deallocate(current%next)
      end do

      return
      end
