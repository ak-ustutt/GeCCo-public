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
      include 'def_formula.h'

      type(filinf), intent(in) ::
     &     ffform
      type(formula), intent(inout) ::
     &     form
      integer ::
     &     lu

      rd_formula = .true.
      lu = ffform%unit

      ! read command record
      read(lu,end=100) form%command, form%target

      select case(form%command)
      case(command_add_contribution)
        if (.not.associated(form%contr)) then
          allocate(form%contr)
          form%contr%mxvtx = 0
          form%contr%mxarc = 0
          form%contr%mxfac = 0
        end if
        form%contr%idx_res = form%target
        call rw_contr_kernel(+1,lu,form%contr)
      case(command_new_intermediate)
        call quit(1,'rd_formula',
     &       'command "new intermediate" not yet implemented')
c        call rw_opdef_kernel(+1,lu,form%interm)
      end select

      return

      ! EOF encountered:
 100  rd_formula = .false.
      return

      end
