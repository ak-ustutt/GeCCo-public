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
        call rw_opdef_kernel(-1,lu,form%interm,
     &                       form%parent1,form%parent2)
      case(command_reorder)
        call rw_reo_kernel(-1,lu,form%reo)
      case(command_add_bc_reo,command_bc_reo)
        call rw_bcontr_kernel(-1,lu,form%bcontr)
        call rw_reo_kernel(-1,lu,form%reo)
      case(command_add_intm,command_add_bc,command_bc)
        call rw_bcontr_kernel(-1,lu,form%bcontr)
      end select

      return

      end
