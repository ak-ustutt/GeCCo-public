*----------------------------------------------------------------------*
      logical function rd_formula(ffform,form)
*----------------------------------------------------------------------*
*     read next formula record from ffform
*----------------------------------------------------------------------*

      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'def_contraction.h'
      include 'def_filinf.h'
      include 'def_operator.h'
      include 'def_formula_item.h'

      type(filinf), intent(in) ::
     &     ffform
      type(formula_item), intent(inout) ::
     &     form
      integer ::
     &     lu, command, target

      rd_formula = .true.
      lu = ffform%unit

      ! read command record
      read(lu,end=100) command, target

c dbg
c      print *,'read: ',command,target
c dbg
      if (command.eq.command_end_of_formula) goto 100
      call new_formula_item(form,command,target)

      select case(form%command)
      case(command_add_contribution)
        form%contr%idx_res = form%target
        call rw_contr_kernel(+1,lu,form%contr)
      case(command_del_intermediate)
        read(lu) form%label
      case(command_new_intermediate)
        call rw_opdef_kernel(+1,lu,form%interm,form%incore,
     &                       form%parent1,form%parent2,
     &                       form%tra,form%tra1,form%tra2)
      case(command_reorder)
        call rw_reo_kernel(+1,lu,form%reo)
      case(command_add_bc_reo,command_bc_reo,command_add_reo)
        call rw_bcontr_kernel(+1,lu,form%bcontr)
        call rw_reo_kernel(+1,lu,form%reo)
      case(command_add_intm,command_cp_intm,command_add_bc,command_bc)
        call rw_bcontr_kernel(+1,lu,form%bcontr)
      end select

      return

      ! EOF encountered:
 100  rd_formula = .false.
      return

      end
