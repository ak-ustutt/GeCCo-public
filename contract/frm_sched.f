*----------------------------------------------------------------------*
      subroutine frm_sched(xret,fffrm,
     &         op_info,str_info,strmap_info,orb_info)
*----------------------------------------------------------------------*
*     wrapper for actual scheduler routines
*----------------------------------------------------------------------*
      implicit none

      include 'routes.h'
      include 'contr_times.h'

      include 'stdunit.h'
      include 'opdim.h'
      include 'mdef_operator_info.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'def_strmapinf.h'
      include 'def_orbinf.h'

      real(8), intent(out) ::
     &     xret(*)
      type(filinf), intent(inout) ::
     &     fffrm
      type(operator_info) ::
     &     op_info
      type(strinf) ::
     &     str_info
      type(strmapinf) ::
     &     strmap_info
      type(orbinf) ::
     &     orb_info

      real(8) ::
     &     cpu, sys, wall, cpu0, sys0, wall0

      call atim_csw(cpu0,sys0,wall0)
      cnt_op1op2(1:3) = 0d0
      cnt_dloop(1:3) = 0d0
      cnt_kernel(1:3) = 0d0

      select case (irt_sched)
      case (0)
        call frm_sched0(xret,fffrm,
     &         op_info,str_info,strmap_info,orb_info)
      case default
        call quit(1,'frm_sched','illegal route')
      end select

      call atim_csw(cpu,sys,wall)
      if (iprlvl.ge.10) then
        call prtim(luout,'contraction kernel',
     &       cnt_kernel(1),cnt_kernel(2),-1d0)
c        call prtim(luout,'overhead contraction 1',
c     &       cnt_op1op2(1)-cnt_kernel(1),
c     &       cnt_op1op2(2)-cnt_kernel(2),
c     &       cnt_op1op2(3)-cnt_kernel(3))
        call prtim(luout,'overhead contraction 1',
     &       cnt_dloop(1)-cnt_kernel(1),
     &       cnt_dloop(2)-cnt_kernel(2),
     &       -1d0)
        call prtim(luout,'overhead contraction 2',
     &       cnt_op1op2(1)-cnt_dloop(1),
     &       cnt_op1op2(2)-cnt_dloop(2),
     &       -1d0)
        call prtim(luout,'overhead scheduler',
     &       cpu-cpu0-cnt_op1op2(1),
     &       sys-sys0-cnt_op1op2(2),
     &       -1d0)
        
      end if
      if (iprlvl.ge.5)
     &     call prtim(luout,'formula evaluation',
     &     cpu-cpu0,sys-sys0,wall-wall0)
      
      return
      end
