*----------------------------------------------------------------------*
      subroutine new_formula_item(form_pnt,command,target)
*----------------------------------------------------------------------*
*     allocate new entry for linked list
*     on entry, form_pnt should point to last entry (containing [END])
*     on exit, form_pnt is ready to take the information coming
*     along with command "command" (i.e. contr or interm is alloc.d)
*----------------------------------------------------------------------*

      implicit none

      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_operator.h'
      include 'def_formula_item.h'
      
      type(formula_item), intent(inout), target ::
     &     form_pnt
      integer, intent(in) ::
     &     command, target

      ! set basic info and initialize arrays
      call init_fl_item(form_pnt,command,target)
      
      ! extend the list
      allocate(form_pnt%next)
      call init_fl_item_0(form_pnt%next) ! init to zero
      form_pnt%next%prev => form_pnt
      
      return
      end
