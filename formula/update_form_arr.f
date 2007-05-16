*----------------------------------------------------------------------*
      subroutine update_form_arr(form_info)
*----------------------------------------------------------------------*
*     setup or update direct acces pointer array for linked list
*----------------------------------------------------------------------*

      implicit none

      include 'def_filinf.h'
      include 'mdef_formula_info.h'

      type(formula_info) ::
     &     form_info

      if (associated(form_info%form_arr)) deallocate(form_info%form_arr)

      if (form_info%nform.gt.0) then
        allocate(form_info%form_arr(form_info%nform))
        call form_list2arr(form_info%form_list,
     &       form_info%form_arr,form_info%nform)
      end if

      return
      end

