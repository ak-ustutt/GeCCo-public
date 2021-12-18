*----------------------------------------------------------------------*
      subroutine insert_fl_item(form_pnt,command,target)
*----------------------------------------------------------------------*
*     allocate new entry for linked list, insert it after the
*     entry that form_pnt is pointing to
*     on entry, form_pnt should point to last entry (containing [END])
*     on exit, form_pnt is ready to take the information coming
*     along with command "command" (i.e. contr or interm is alloc.d)
*----------------------------------------------------------------------*

      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'def_contraction.h'
      include 'def_operator.h'
      include 'def_formula_item.h'
      
      type(formula_item), intent(inout), target ::
     &     form_pnt
      integer, intent(in) ::
     &     command, target

      type(formula_item), pointer ::
     &     form_pnt_after

      if (.not.associated(form_pnt%next)) then
        write(lulog,*) 'illegal position for insert! command = ',
     &                  form_pnt%command
        call quit(1,'insert_fl_item','illegal position')
      end if
         

      ! remember the next entry in list
      form_pnt_after => form_pnt%next

      form_pnt%next => null()
      ! extend the list
      allocate(form_pnt%next)
      call init_fl_item_0(form_pnt%next)
      ! set all links
c      if (.not.associated(form_pnt%prev)) 
c     &     call quit(1,'insert_formula_item',
c     &                 'at present I cannot insert after 1st entry!')
c      form_pnt%next%prev => form_pnt%prev%next  ! better get the pointer from here
      form_pnt%next%prev => form_pnt
      form_pnt%next%next => form_pnt_after
      form_pnt%next%next%prev => form_pnt%next

      ! set basic info and initialize arrays
      call init_fl_item(form_pnt%next,command,target)
      
      return
      end

