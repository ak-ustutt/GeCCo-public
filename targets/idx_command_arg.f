*----------------------------------------------------------------------*
      integer function idx_command_arg(arg_label,rule)
*----------------------------------------------------------------------*
*     given an argument label, search rule and return index of
*     of corresponding argument
*----------------------------------------------------------------------*

      implicit none

      include 'stdunit.h'
      include 'mdef_target_info.h'

      integer, parameter ::
     &     ntest = 00

      character(len=*), intent(in) ::
     &     arg_label
      type(action), intent(in) ::
     &     rule

      integer ::
     &     iarg

      if (ntest.ge.100) then
        call write_title(lulog,wst_dbg_subr,
     &       'this is idx_command_arg')
        write(lulog,*) ' looking for: "',trim(arg_label),'"'
      end if

      idx_command_arg = -1
      do iarg = 1, rule%n_arguments
        if (trim(arg_label).eq.trim(rule%arg(iarg)%arg_label)) then
          idx_command_arg = iarg
          exit
        end if
      end do

      if (ntest.ge.100) then
        write(lulog,*) 'result: ',idx_command_arg
      end if

      return
      end
