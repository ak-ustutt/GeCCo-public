*----------------------------------------------------------------------*
      subroutine frm_sched(xret,flist,depend_info,idxselect,nselect,
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
      include 'def_contraction.h'
      include 'def_formula_item.h'
      include 'def_dependency_info.h'

      real(8), intent(out) ::
     &     xret(*)
      integer, intent(in) ::
     &     nselect, idxselect
      type(formula_item), intent(inout) ::
     &     flist
      type(operator_info) ::
     &     op_info
      type(strinf) ::
     &     str_info
      type(strmapinf) ::
     &     strmap_info
      type(orbinf) ::
     &     orb_info
      type(dependency_info) ::
     &     depend_info

      real(8) ::
     &     cpu, sys, wall, cpu0, sys0, wall0

      call atim_csw(cpu0,sys0,wall0)
      cnt_op1op2(1:3) = 0d0
      cnt_dloop(1:3) = 0d0
      cnt_kernel(1:3) = 0d0
      cnt_coll1(1:2) = 0d0
      cnt_coll2(1:2) = 0d0
      cnt_dgemm(1:2) = 0d0
      cnt_scatt(1:2) = 0d0
      cnt_rd(1:3) = 0d0
      cnt_wr(1:3) = 0d0

      select case (irt_sched)
      case (0)
        call frm_sched1(xret,flist,depend_info,idxselect,nselect,
     &         op_info,str_info,strmap_info,orb_info)
      case default
        call quit(1,'frm_sched','illegal route')
      end select

      call atim_csw(cpu,sys,wall)
      if (iprlvl.ge.10) then
        call prtim(luout,'contraction kernel',
     &       cnt_kernel(1),cnt_kernel(2),-1d0)
        call prtim(luout,'     in collect 1',
     &       cnt_coll1(1),cnt_coll1(2),-1d0)
        call prtim(luout,'     in collect 2',
     &       cnt_coll2(1),cnt_coll2(2),-1d0)
        call prtim(luout,'         in dgemm',
     &       cnt_dgemm(1),cnt_dgemm(2),-1d0)
        call prtim(luout,'       in scatter',
     &       cnt_scatt(1),cnt_scatt(2),-1d0)
        call prtim(luout,'overhead contraction 1',
     &       cnt_dloop(1)-cnt_kernel(1),
     &       cnt_dloop(2)-cnt_kernel(2),
     &       -1d0)
        call prtim(luout,'IO read',
     &      cnt_rd(1),cnt_rd(2),-1d0)
        call prtim(luout,'IO write',
     &      cnt_wr(1),cnt_wr(2),-1d0)
        call prtim(luout,'overhead contraction 2',
     &       cnt_op1op2(1)-cnt_dloop(1)-cnt_rd(1)-cnt_wr(1),
     &       cnt_op1op2(2)-cnt_dloop(2)-cnt_rd(2)-cnt_wr(2),
     &       -1d0)
        call prtim(luout,'overhead scheduler',
     &       cpu-cpu0-cnt_op1op2(1),
     &       sys-sys0-cnt_op1op2(2),
     &       -1d0)
        
      end if
      if (iprlvl.ge.1)
     &     call prtim(luout,'formula evaluation',
     &     cpu-cpu0,sys-sys0,wall-wall0)
      
      return
      end
