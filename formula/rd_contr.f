*----------------------------------------------------------------------*
      logical function rd_contr(lu,contr,idx_res)
*----------------------------------------------------------------------*
*     read next contraction from lu and expand to long form on contr
*     return .false. if EOF
*     idx_res is obsolete now
*----------------------------------------------------------------------*

      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'def_contraction.h'
      include 'def_operator.h'
      include 'def_formula_item.h'

      integer, intent(in) ::
     &     lu, idx_res
      type(contraction), intent(out) ::
     &     contr
      integer ::
     &     command

      rd_contr = .false.

      ! read command record
      read(lu,end=100) command, contr%idx_res

      if (command.ne.command_add_contribution) then
        backspace lu
        return
      end if

      rd_contr = .true.

      call rw_contr_kernel(+1,lu,contr)
      return

      ! EOF encountered:
 100  rd_contr = .false.
      return

      end
