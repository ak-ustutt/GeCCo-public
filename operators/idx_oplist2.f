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
     &     iop, llabel

      if (ntest.ge.100) then
        write(lulog,*) '------------------'
        write(lulog,*) 'this is idx_oplist'
        write(lulog,*) '------------------'
        write(lulog,*) ' looking for: "',trim(opname),'"'
      end if

C      idx_oplist2 = -1
C      do iop = 1, op_info%nops
C        if (trim(opname).eq.trim(op_info%op_arr(iop)%op%name)) then
C          idx_oplist2 = iop
C          exit
C        end if
C      end do

      idx_oplist2 = -1
      llabel = len_trim(opname)
      do iop = 1, op_info%nops
        if (opname(1:llabel).eq.op_info%op_arr(iop)%op%name(1:llabel))
     &  then
          if (len_trim(op_info%op_arr(iop)%op%name).eq.llabel) then
            idx_oplist2 = iop
            exit
          end if
        end if
      end do

      if (ntest.ge.100) then
        write(lulog,*) 'result: ',idx_oplist2
      end if

      return
      end
