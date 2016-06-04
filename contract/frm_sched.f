*----------------------------------------------------------------------*
      subroutine frm_sched(xret,flist,depend_info,idxselect,nselect,
     &         init,no_skip,op_info,str_info,strmap_info,orb_info)
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
      logical, intent(in) ::
     &     init, no_skip
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
c stat
      real(8) ::
     &     avg, sig
c stat

      call atim_csw(cpu0,sys0,wall0)
      cnt_op1op2(1:3) = 0d0
      cnt_dloop(1:3) = 0d0
      cnt_kernel(1:3) = 0d0
      cnt_coll1(1:2) = 0d0
      cnt_coll2(1:2) = 0d0
      cnt_dgemm(1:2) = 0d0
      cnt_scatt(1:2) = 0d0
      cnt_reo(1:2) = 0d0
      cnt_rd(1:3) = 0d0
      cnt_wr(1:3) = 0d0
      cnt_maxscr = 0
      cnt_used_reo = .false.
c dbg
      cnt_test(1:12) = 0d0
c dbg
c stat
      mm_call   = 0
      mm_dim1   = 0
      mm_dim1sq = 0
      mm_dim2   = 0
      mm_dim2sq = 0
      mm_cnt    = 0
      mm_cntsq  = 0
c stat

      select case (irt_sched)
      case (0)
        if (.not.init) call quit(1,'frm_sched',
     &                      'init option not (yet) for this scheduler')
        call frm_sched1(xret,flist,depend_info,idxselect,nselect,
     &         op_info,str_info,strmap_info,orb_info)
      case (1)
        call frm_sched2(xret,flist,depend_info,idxselect,nselect,
     &         init,no_skip,op_info,str_info,strmap_info,orb_info)
      case default
        call quit(1,'frm_sched','illegal route')
      end select

      call atim_csw(cpu,sys,wall)
      if (.false..and.iprlvl.ge.10) then
        call prtim(lulog,'contraction kernel',
     &       cnt_kernel(1),cnt_kernel(2),-1d0)
c dbg
        call prtim(lulog,'check kernel     ',
     &       cnt_test(1),cnt_test(2),-1d0)
        call prtim(lulog,'     in setup    ',
     &       cnt_test(3),cnt_test(4),-1d0)
        call prtim(lulog,'     inner loop ?',
     &       cnt_test(11),cnt_test(12),-1d0)
        call prtim(lulog,'     inner loop  ',
     &       cnt_test(9),cnt_test(10),-1d0)
        call prtim(lulog,'    middle loop  ',
     &       cnt_test(7)-cnt_test(9),cnt_test(8)-cnt_test(10),-1d0)
        call prtim(lulog,'     outer loop  ',
     &       cnt_test(5)-cnt_test(7),cnt_test(6)-cnt_test(8),-1d0)
c dbg
        call prtim(lulog,'     in collect 1',
     &       cnt_coll1(1),cnt_coll1(2),-1d0)
        call prtim(lulog,'     in collect 2',
     &       cnt_coll2(1),cnt_coll2(2),-1d0)
        call prtim(lulog,'         in dgemm',
     &       cnt_dgemm(1),cnt_dgemm(2),-1d0)
        call prtim(lulog,'       in scatter',
     &       cnt_scatt(1),cnt_scatt(2),-1d0)
        if (cnt_used_reo)
     &       call prtim(lulog,'additional reord.',
     &       cnt_reo(1),cnt_reo(2),-1d0)
        call prtim(lulog,'overhead contraction 1',
     &       cnt_dloop(1)-cnt_kernel(1)-cnt_reo(1),
     &       cnt_dloop(2)-cnt_kernel(2)-cnt_reo(2),
     &       -1d0)
        call prtim(lulog,'IO read',
     &      cnt_rd(1),cnt_rd(2),-1d0)
        call prtim(lulog,'IO write',
     &      cnt_wr(1),cnt_wr(2),-1d0)
        call prtim(lulog,'overhead contraction 2',
     &       cnt_op1op2(1)-cnt_dloop(1)-cnt_rd(1)-cnt_wr(1),
     &       cnt_op1op2(2)-cnt_dloop(2)-cnt_rd(2)-cnt_wr(2),
     &       -1d0)
        call prtim(lulog,'overhead scheduler',
     &       cpu-cpu0-cnt_op1op2(1),
     &       sys-sys0-cnt_op1op2(2),
     &       -1d0)
c stat
        write(lulog,*) 'calls to dgemm: ',mm_call
        avg = dble(mm_dim1)/dble(mm_call)
        sig = sqrt(abs(avg*avg - dble(mm_dim1sq)/dble(mm_call)))
        write(lulog,'(2x,a,2f12.2)') ' avg dim1, sigma:',avg,sig
        avg = dble(mm_dim2)/dble(mm_call)
        sig = sqrt(abs(avg*avg - dble(mm_dim2sq)/dble(mm_call)))
        write(lulog,'(2x,a,2f12.2)') ' avg dim2, sigma:',avg,sig
        avg = dble(mm_cnt)/dble(mm_call)
        sig = sqrt(abs(avg*avg - dble(mm_cntsq)/dble(mm_call)))
        write(lulog,'(2x,a,2f12.2)') ' avg cnt , sigma:',avg,sig
        write(lulog,'(/2x,a,i10,a,f12.3,a)')
     &       'max. scratch: ',cnt_maxscr,' = ',
     &       dble(cnt_maxscr)*8d0/1024d0**3d0,' Gbytes'
c stat
        
      end if
      if (iprlvl.ge.1)
     &     call prtim(lulog,'formula evaluation',
     &     cpu-cpu0,sys-sys0,wall-wall0)
      
      return
      end
