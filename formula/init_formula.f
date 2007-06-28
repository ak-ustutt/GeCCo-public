*----------------------------------------------------------------------*
      subroutine init_formula(form)
*----------------------------------------------------------------------*
*     initialize all pointers in first formula item and set to
*     [END]
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_operator.h'
      include 'def_formula_item.h'

      type(formula_item), intent(inout) ::
     &     form

      form%prev => null()
      form%next => null()
      form%contr => null()
      form%interm => null()
      form%command = command_end_of_formula
      form%target = 0

      return
      end
