*----------------------------------------------------------------------*
      subroutine contr_op1op2(xfac,bc_sign,ffop1,ffop2,
     &     update,ffop1op2,xret,type_xret,
     &     op1,op2,op1op2,op1op2tmp,
     &     iblkop1,iblkop2,iblkop1op2,iblkop1op2tmp,
     &     idoffop1,idoffop2,idoffop1op2,
     &     iocc_ex1,iocc_ex2,iocc_cnt,
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
*----------------------------------------------------------------------*
      implicit none

      include 'routes.h'
      include 'contr_times.h'

      include 'opdim.h'
      include 'stdunit.h'
      include 'ioparam.h'
      include 'def_operator.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'def_orbinf.h'
      include 'def_filinf.h'
      include 'def_strmapinf.h'
      include 'def_contraction_info.h'
      include 'def_reorder_info.h'

      logical, intent(in) ::
     &     update
      real(8), intent(in) ::
     &     xfac, bc_sign
      real(8), intent(inout) ::
     &     xret
      integer, intent(in) ::
     &     type_xret,
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
     &     mstop1,mstop2,mstop1op2,
     &     igamtop1,igamtop2,igamtop1op2
      type(filinf), intent(in) ::
     &     ffop1,ffop2,ffop1op2
      type(operator), intent(in) ::
     &     op1, op2, op1op2, op1op2tmp
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
c      integer ::
c     &     nca_blk(2,7)
c      integer, pointer ::
c     &     cinfo_op1c(:,:),cinfo_op1a(:,:),
c     &     cinfo_op2c(:,:),cinfo_op2a(:,:),
c     &     cinfo_op1op2c(:,:),
c     &     cinfo_op1op2a(:,:),
c     &     cinfo_op1op2tmpc(:,:),
c     &     cinfo_op1op2tmpa(:,:),
c     &     cinfo_ex1c(:,:),cinfo_ex1a(:,:),
c     &     cinfo_ex2c(:,:),cinfo_ex2a(:,:),
c     &     cinfo_cntc(:,:),cinfo_cnta(:,:),
c     &     map_info1c(:),
c     &     map_info1a(:),
c     &     map_info2c(:),
c     &     map_info2a(:),
c     &     map_info12c(:),
c     &     map_info12a(:)
      
      real(8) ::
     &     cpu, sys, cpu0, sys0

      call atim_cs(cpu0,sys0)
      select case (irt_contr)
      case(0)
        call contr_op1op2_simple(xfac,ffop1,ffop2,
     &     update,ffop1op2,xret,type_xret,
     &     op1,op2,op1op2,
     &     iblkop1,iblkop2,iblkop1op2,
     &     idoffop1,idoffop2,idoffop1op2,
     &     iocc_ex1,iocc_ex2,iocc_cnt,
     &     irst_op1,irst_op2,irst_op1op2,
     &     mstop1,mstop2,mstop1op2,
     &     igamtop1,igamtop2,igamtop1op2,
     &     str_info,orb_info)
      case(1)
        call contr_op1op2_wmaps(xfac,ffop1,ffop2,
     &     update,ffop1op2,xret,type_xret,
     &     op1,op2,op1op2,
     &     iblkop1,iblkop2,iblkop1op2,
     &     idoffop1,idoffop2,idoffop1op2,
     &     iocc_ex1,iocc_ex2,iocc_cnt,
     &     irst_op1,irst_op2,irst_op1op2,
     &     mstop1,mstop2,mstop1op2,
     &     igamtop1,igamtop2,igamtop1op2,
     &     str_info,strmap_info,orb_info)
      case(2,3)
        ! reform to "condensed" description of contraction
        call init_cnt_info(cnt_info,
     &     iocc_op1,iocc_ex1,njoined_op1,iocc_op2,iocc_ex2,njoined_op2,
     &     iocc_cnt,njoined_cnt,
     &     iocc_op1op2,njoined_op1op2,iocc_op1op2tmp,njoined_op1op2)

        call condense_bc_info(
     &       cnt_info,
     &       iocc_op1, iocc_op2, iocc_op1op2, iocc_op1op2tmp,
     &       iocc_ex1,iocc_ex2,iocc_cnt,
     &       irst_op1, irst_op2, irst_op1op2, irst_op1op2tmp,
     &       merge_op1, merge_op2, merge_op1op2, merge_op2op1,
     &       njoined_op1, njoined_op2,njoined_op1op2, njoined_cnt,
     &       str_info,orb_info%ihpvgas,orb_info%ngas)
        call contr_op1op2_wmaps_c(xfac,bc_sign,ffop1,ffop2,
     &       update,ffop1op2,xret,type_xret,
     &       op1,op2,op1op2,op1op2tmp,
     &       iblkop1,iblkop2,iblkop1op2,iblkop1op2tmp,
     &       idoffop1,idoffop2,idoffop1op2,
     &       cnt_info,reo_info,
     &       str_info,strmap_info,orb_info)

        call dealloc_cnt_info(cnt_info)

      case default
        write(luout,*) 'contr_op1op2: route = ',irt_contr
        call quit(1,'contr_op1op2','route not implemented')
      end select
      call atim_cs(cpu,sys)
      cnt_op1op2(1) = cnt_op1op2(1)+cpu-cpu0
      cnt_op1op2(2) = cnt_op1op2(2)+sys-sys0

      return
      end
