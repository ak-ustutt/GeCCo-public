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
c      form%next => null()
      call init_fl_item_0(form)
c      form%contr => null()
c      form%interm => null()
c      form%command = command_end_of_formula
c      form%target = 0

      return
      end
