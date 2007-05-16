*----------------------------------------------------------------------*
      integer function idx_oplist2(opname,op_info)
*----------------------------------------------------------------------*
*     given an operator name, search op_info and return index of
*     of corresponding operator
*----------------------------------------------------------------------*

      implicit none

      include 'stdunit.h'
      include 'mdef_operator_info.h'

      integer, parameter ::
     &     ntest = 00

      character, intent(in) ::
     &     opname*(*)
      type(operator_info), intent(in) ::
     &     op_info

      integer ::
     &     iop

      if (ntest.ge.100) then
        write(luout,*) '------------------'
        write(luout,*) 'this is idx_oplist'
        write(luout,*) '------------------'
        write(luout,*) ' looking for: "',trim(opname),'"'
      end if

      idx_oplist2 = -1
      do iop = 1, op_info%nops
        if (trim(opname).eq.trim(op_info%op_arr(iop)%op%name)) then
          idx_oplist2 = iop
          exit
        end if
      end do

      if (ntest.ge.100) then
        write(luout,*) 'result: ',idx_oplist2
      end if

      return
      end
