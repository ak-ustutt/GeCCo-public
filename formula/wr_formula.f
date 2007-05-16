*----------------------------------------------------------------------*
      subroutine wr_formula(ffform,form)
*----------------------------------------------------------------------*
*     write next formula record to ffform
*----------------------------------------------------------------------*

      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'def_contraction.h'
      include 'def_operator.h'
      include 'def_filinf.h'
      include 'def_formula_item.h'

      type(filinf), intent(in) ::
     &     ffform
      type(formula_item), intent(in) ::
     &     form
      integer ::
     &     lu

      lu = ffform%unit

      ! write command record
      write(lu) form%command, form%target

      select case(form%command)
      case(command_add_contribution)
        call rw_contr_kernel(-1,lu,form%contr)
      case(command_new_intermediate)
        call quit(1,'wr_formula',
     &       'command "new intermediate" not yet implemented')
c        call rw_opdef_kernel(-1,lu,form%interm)
      end select

      return

      end
