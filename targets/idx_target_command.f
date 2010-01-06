*----------------------------------------------------------------------*
      integer function idx_target_command(mode,command,tgt)
*----------------------------------------------------------------------*
*     given a command name, search tgt and return index of
*     of corresponding rule
*     mode >= 0 first occurrence
*     mode <  0 last occurrence
*----------------------------------------------------------------------*

      implicit none

      include 'stdunit.h'
      include 'mdef_target_info.h'

      integer, parameter ::
     &     ntest = 00

      integer, intent(in) ::
     &     mode
      character(len=*), intent(in) ::
     &     command
      type(target), intent(in) ::
     &     tgt

      integer ::
     &     icmd

      if (ntest.ge.100) then
        call write_title(luout,wst_dbg_subr,
     &       'this is idx_target_command')
        write(luout,*) ' looking for: "',trim(command),'"'
      end if

      idx_target_command = -1

      if (mode.ge.0) then
        do icmd = 1, tgt%n_rules
          if (trim(command).eq.trim(tgt%rules(icmd)%command)) then
            idx_target_command = icmd
            exit
          end if
        end do
      else
        do icmd = tgt%n_rules, 1, -1
          if (trim(command).eq.trim(tgt%rules(icmd)%command)) then
            idx_target_command = icmd
            exit
          end if
        end do
      end if

      if (ntest.ge.100) then
        write(luout,*) 'result: ',idx_target_command
      end if

      return
      end
