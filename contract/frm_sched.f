*----------------------------------------------------------------------*
      subroutine frm_sched(xret,fffrm,
     &         op_info,str_info,orb_info)
*----------------------------------------------------------------------*
*     wrapper for actual scheduler routines
*----------------------------------------------------------------------*
      implicit none

      include 'routes.h'

      include 'stdunit.h'
      include 'opdim.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'def_orbinf.h'
      include 'mdef_operator_info.h'

      real(8), intent(out) ::
     &     xret(*)
      type(filinf), intent(inout) ::
     &     fffrm
      type(operator_info) ::
     &     op_info
      type(strinf) ::
     &     str_info
      type(orbinf) ::
     &     orb_info

      real(8) ::
     &     cpu, sys, wall, cpu0, sys0, wall0

      call atim(cpu0,sys0,wall0)

      select case (irt_sched)
      case (0)
        call frm_sched0(xret,fffrm,
     &         op_info,str_info,orb_info)
      case default
        call quit(1,'frm_sched','illegal route')
      end select

      call atim(cpu,sys,wall)
      if (iprlvl.ge.5)
     &     call prtim(luout,'formula evaluation',
     &     cpu-cpu0,sys-sys0,wall-wall0)
      
      return
      end
