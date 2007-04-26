*----------------------------------------------------------------------*
      integer function idx_oplist(opname,ops,nops)
*----------------------------------------------------------------------*
*     given an operator name, search ops(nops) and return index of
*     of corresponding operator
*----------------------------------------------------------------------*

      implicit none

      include 'stdunit.h'
      include 'def_operator.h'

      integer, parameter ::
     &     ntest = 00

      character, intent(in) ::
     &     opname*(*)
      integer, intent(in) ::
     &     nops
      type(operator), intent(in) ::
     &     ops(nops)

      integer ::
     &     iop

      if (ntest.ge.100) then
        write(luout,*) '------------------'
        write(luout,*) 'this is idx_oplist'
        write(luout,*) '------------------'
        write(luout,*) ' looking for: "',trim(opname),'"'
      end if

      idx_oplist = -1
      do iop = 1, nops
        if (trim(opname).eq.trim(ops(iop)%name)) then
          idx_oplist = iop
          exit
        end if
      end do

      if (ntest.ge.100) then
        write(luout,*) 'result: ',idx_oplist
      end if

      return
      end
