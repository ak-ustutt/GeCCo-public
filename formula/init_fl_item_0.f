      subroutine init_fl_item_0(fl_item)

      implicit none

      include 'opdim.h'
      include 'def_operator.h'
      include 'def_contraction.h'
      include 'def_formula_item.h'

      type(formula_item), intent(inout) ::
     &     fl_item

      fl_item%next => null()
      fl_item%contr => null()
      fl_item%interm => null()
      fl_item%incore = 0
      fl_item%parent1 => null()
      fl_item%parent2 => null()
      fl_item%label => null()
      fl_item%reo => null()
      fl_item%bcontr => null()
      fl_item%command = command_end_of_formula
      fl_item%target = 0

      return
      end
