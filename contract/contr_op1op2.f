*----------------------------------------------------------------------*
      subroutine contr_op1op2(xfac,bc_sign,
     &     update,self,xret,type_xret,
     &     me_op1,me_op2,me_op1op2,me_op1op2tmp,
     &     tra_op1, tra_op2, tra_op1op2,
     &     iblkop1,iblkop2,iblkop1op2,iblkop1op2tmp,
     &     idoffop1,idoffop2,idoffop1op2,
     &     iocc_ex1,iocc_ex2,iocc_cnt,
     &     irst_ex1,irst_ex2,irst_cnt,mode_rst,
     &     iocc_op1, iocc_op2, iocc_op1op2,iocc_op1op2tmp,
     &     irst_op1,irst_op2,irst_op1op2,irst_op1op2tmp,
     &     merge_op1, merge_op2, merge_op1op2, merge_op2op1,
     &     njoined_op1, njoined_op2,njoined_op1op2, njoined_cnt,
     &     mstop1,mstop2,mstop1op2,
     &     igamtop1,igamtop2,igamtop1op2,
     &     reo_info,
     &     str_info,strmap_info,orb_info)
*----------------------------------------------------------------------*
*     wrapper for contraction routines
*     mode_rst = 0  ->  for use in frm_sched1
*     mode_rst = 2  ->  for use in frm_sched2
*----------------------------------------------------------------------*
      implicit none

      include 'routes.h'
      include 'contr_times.h'

      include 'opdim.h'
      include 'stdunit.h'
      include 'ioparam.h'
      include 'def_operator.h'
      include 'def_me_list.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'def_orbinf.h'
      include 'def_filinf.h'
      include 'def_strmapinf.h'
      include 'def_contraction_info.h'
      include 'def_reorder_info.h'

      logical, intent(in) ::
     &     update, self, tra_op1, tra_op2, tra_op1op2
      real(8), intent(in) ::
     &     xfac, bc_sign
      real(8), intent(inout) ::
     &     xret
      integer, intent(in) ::
     &     type_xret, mode_rst,
     &     iblkop1, iblkop2, iblkop1op2, iblkop1op2tmp,
     &     idoffop1, idoffop2, idoffop1op2,
     &     njoined_op1, njoined_op2, njoined_op1op2, njoined_cnt,
     &     iocc_ex1(ngastp,2,njoined_op1),
     &     iocc_ex2(ngastp,2,njoined_op2),
     &     iocc_cnt(ngastp,2,njoined_cnt),
     &     iocc_op1(ngastp,2,njoined_op1),
     &     iocc_op2(ngastp,2,njoined_op2),
     &     iocc_op1op2(ngastp,2,njoined_op1op2),
     &     iocc_op1op2tmp(ngastp,2,njoined_op1op2),
     &     merge_op1(*), merge_op2(*), merge_op1op2(*), merge_op2op1(*),
     &     irst_op1(*), irst_op2(*), irst_op1op2(*), irst_op1op2tmp(*),
     &     irst_ex1(*), irst_ex2(*), irst_cnt(*), 
     &     mstop1,mstop2,mstop1op2,
     &     igamtop1,igamtop2,igamtop1op2
      type(me_list), intent(in) ::
     &     me_op1, me_op2, me_op1op2, me_op1op2tmp
      type(strinf), intent(in) ::
     &     str_info
      type(strmapinf), intent(in) ::
     &     strmap_info
      type(orbinf), intent(in) ::
     &     orb_info
      type(reorder_info), intent(in) ::
     &     reo_info

      type(contraction_info) ::
     &     cnt_info
      
      integer ::
     &     idum
      real(8) ::
     &     cpu, sys, cpu0, sys0

c dbg
c      print *,'op1 ',trim(me_op1%op%name)
c      print *,'op2 ',trim(me_op2%op%name)
c      print *,'op1op2 ',trim(me_op1op2%op%name)
c      print *,'op1op2tmp ',trim(me_op1op2tmp%op%name)
c dbg

      if (mode_rst.ne.0.and.mode_rst.ne.2)
     &     call quit(1,'contr_op1op2',
     &     'illegal value for mode_rst (must be 0 or 2)')

c      call atim_cs(cpu0,sys0)
      select case (irt_contr)
      case(0)
        call quit(1,'contr_op1op2','route 0 is obsolete')
      case(1)
        call quit(1,'contr_op1op2','route 1 is obsolete')
c        call contr_op1op2_wmaps(xfac,ffop1,ffop2,
c     &     update,ffop1op2,xret,type_xret,
c     &     op1,op2,op1op2,
c     &     iblkop1,iblkop2,iblkop1op2,
c     &     idoffop1,idoffop2,idoffop1op2,
c     &     iocc_ex1,iocc_ex2,iocc_cnt,
c     &     irst_op1,irst_op2,irst_op1op2,
c     &     mstop1,mstop2,mstop1op2,
c     &     igamtop1,igamtop2,igamtop1op2,
c     &     str_info,strmap_info,orb_info)
      case(2,3)
        ! reform to "condensed" description of contraction
        call init_cnt_info(cnt_info,
     &     iocc_op1,iocc_ex1,njoined_op1,iocc_op2,iocc_ex2,njoined_op2,
     &     iocc_cnt,njoined_cnt,
     &     iocc_op1op2,njoined_op1op2,iocc_op1op2tmp,njoined_op1op2)

c dbg
c        print *,'bef. call to condens'
c        print *,'iblkop1op2,iblkop1op2tmp: ',iblkop1op2,iblkop1op2tmp
c        print *,'iocc_op1op2:'
c        call wrt_occ_n(6,iocc_op1op2,njoined_op1)
c dbg
        call condense_bc_info(
     &       cnt_info,self,
     &       iocc_op1, iocc_op2, iocc_op1op2, iocc_op1op2tmp,
     &       iocc_ex1,iocc_ex2,iocc_cnt,
     &       irst_op1, irst_op2, irst_op1op2, irst_op1op2tmp,
     &       irst_ex1, irst_ex2, irst_cnt, mode_rst,
     &       merge_op1, merge_op2, merge_op1op2, merge_op2op1,
     &       njoined_op1, njoined_op2,njoined_op1op2, njoined_cnt,
     &       str_info,orb_info)
c dbg
c        print *,'bef. call to contr'
c        print *,'iblkop1op2,iblkop1op2tmp: ',iblkop1op2,iblkop1op2tmp
c        print *,'iocc_op1op2:'
c        call wrt_occ_n(6,iocc_op1op2,njoined_op1op2)
c        call wrt_occ_n(6,iocc_op1op2tmp,njoined_op1op2)
c dbg
        if (.not.self) then
          call contr_op1op2_wmaps_c(xfac,bc_sign,
     &       update,xret,type_xret,
     &       me_op1,me_op2,me_op1op2,me_op1op2tmp,
     &       tra_op1, tra_op2, tra_op1op2,
     &       iblkop1,iblkop2,iblkop1op2,iblkop1op2tmp,
     &       idoffop1,idoffop2,idoffop1op2,
     &       cnt_info,reo_info,
     &       str_info,strmap_info,orb_info)
        else
          call trace_op(xfac,bc_sign,
     &     update,xret,type_xret,
     &     me_op1,me_op1op2,me_op1op2tmp,
     &     tra_op1, tra_op1op2,
     &     iblkop1,iblkop1op2,iblkop1op2tmp,
     &     idoffop1,idoffop1op2,
     &     cnt_info,reo_info,
     &     str_info,strmap_info,orb_info)

        end if

        call dealloc_cnt_info(cnt_info)

      case default
        write(lulog,*) 'contr_op1op2: route = ',irt_contr
        call quit(1,'contr_op1op2','route not implemented')
      end select
c      call atim_cs(cpu,sys)
c      cnt_op1op2(1) = cnt_op1op2(1)+cpu-cpu0
c      cnt_op1op2(2) = cnt_op1op2(2)+sys-sys0

      return
      end
