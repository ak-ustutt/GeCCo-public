      subroutine fl_item_rename_res(ok,fl_item,label_new,label_old)

      ! wrapper for renaming the result of a contraction/reordering

      implicit none

      include 'opdim.h'
      include 'def_orbinf.h'
      include 'def_operator.h'
      include 'def_contraction.h'
      include 'def_formula_item.h'


      character(len=19) ::
     &     i_am = 'fl_item_rename_intm'

      type(formula_item), intent(inout) ::
     &     fl_item
      character(len=*), intent(in) ::
     &     label_new, label_old
      logical ::
     &     ok


      select case (fl_item%command)
      case(command_bc,command_bc_reo,
     &     command_add_bc,command_add_bc_reo)
        if (trim(fl_item%bcontr%label_res).eq.trim(label_old)) then 
          ok = .true.
          fl_item%bcontr%label_res = label_new
        end if 
      case(command_reorder)
        if (trim(fl_item%reo%label_out).eq.trim(label_old)) then 
          ok = .true.
          fl_item%reo%label_out = label_new
        end if
      case default
          ok = .false.
      end select

      end

