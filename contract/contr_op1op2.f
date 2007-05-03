*----------------------------------------------------------------------*
      subroutine contr_op1op2(xfac,ffop1,ffop2,
     &     update,ffop1op2,xret,
     &     op1,op2,op1op2,
     &     iblkop1,iblkop2,iblkop1op2,
     &     iocc_ext1,iocc_ext2,iocc_cnt,
     &     irst_op1,irst_op2,irst_op1op2,
     &     mstop1,mstop2,mstop1op2,
     &     igamtop1,igamtop2,igamtop1op2,
     &     str_info,orb_info)
*----------------------------------------------------------------------*
*     wrapper for contraction routines
*----------------------------------------------------------------------*
      implicit none

      include 'routes.h'

      include 'opdim.h'
      include 'stdunit.h'
      include 'ioparam.h'
      include 'def_operator.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'def_orbinf.h'
      include 'def_filinf.h'

      logical, intent(in) ::
     &     update
      real(8), intent(in) ::
     &     xfac
      real(8), intent(inout) ::
     &     xret
      integer, intent(in) ::
     &     iblkop1, iblkop2, iblkop1op2,
     &     iocc_ext1(ngastp,2), iocc_ext2(ngastp,2), iocc_cnt(ngastp,2),
     &     irst_op1(*), irst_op2(*), irst_op1op2(*),
     &     mstop1,mstop2,mstop1op2,
     &     igamtop1,igamtop2,igamtop1op2
      type(filinf), intent(in) ::
     &     ffop1,ffop2,ffop1op2
      type(operator), intent(in) ::
     &     op1, op2, op1op2
      type(strinf), intent(in) ::
     &     str_info
      type(orbinf), intent(in) ::
     &     orb_info

      select case (irt_contr)
      case(0)
        call contr_op1op2_simple(xfac,ffop1,ffop2,
     &     update,ffop1op2,xret,
     &     op1,op2,op1op2,
     &     iblkop1,iblkop2,iblkop1op2,
     &     iocc_ext1,iocc_ext2,iocc_cnt,
     &     irst_op1,irst_op2,irst_op1op2,
     &     mstop1,mstop2,mstop1op2,
     &     igamtop1,igamtop2,igamtop1op2,
     &     str_info,orb_info)
      case default
        write(luout,*) 'contr_op1op2: route = ',irt_contr
        call quit(1,'contr_op1op2','route not implemented')
      end select

      return
      end
