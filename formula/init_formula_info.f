      subroutine init_formula_info(form_info)

      implicit none

      include 'def_filinf.h'
      include 'mdef_formula_info.h'

      type(formula_info), intent(inout) ::
     &     form_info

      ! no formulae yet
      form_info%nform = 0
      ! initialize list
      allocate(form_info%form_list)
      nullify(form_info%form_list%form)
      nullify(form_info%form_list%prev)
      nullify(form_info%form_list%next)
      ! initialize pointer array
      nullify(form_info%form_arr)

      return
      end
