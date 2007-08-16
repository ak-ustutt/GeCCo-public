*----------------------------------------------------------------------*
      subroutine contr_op1op2(xfac,bc_sign,ffop1,ffop2,
     &     update,ffop1op2,xret,type_xret,
     &     op1,op2,op1op2,
     &     iblkop1,iblkop2,iblkop1op2,
     &     idoffop1,idoffop2,idoffop1op2,
     &     iocc_ex1,iocc_ex2,iocc_cnt,
     &     iocc_op1, iocc_op2, iocc_op1op2,
     &     irst_op1,irst_op2,irst_op1op2,
     &     merge_op1, merge_op2, merge_op1op2, merge_op2op1,
     &     njoined_op1, njoined_op2,njoined_op1op2, njoined_cnt,
     &     mstop1,mstop2,mstop1op2,
     &     igamtop1,igamtop2,igamtop1op2,
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

      logical, intent(in) ::
     &     update
      real(8), intent(in) ::
     &     xfac, bc_sign
      real(8), intent(inout) ::
     &     xret
      integer, intent(in) ::
     &     type_xret,
     &     iblkop1, iblkop2, iblkop1op2,
     &     idoffop1, idoffop2, idoffop1op2,
     &     njoined_op1, njoined_op2, njoined_op1op2, njoined_cnt,
     &     iocc_ex1(ngastp,2,njoined_op1),
     &     iocc_ex2(ngastp,2,njoined_op2),
     &     iocc_cnt(ngastp,2,njoined_cnt),
     &     iocc_op1(ngastp,2,njoined_op1),
     &     iocc_op2(ngastp,2,njoined_op2),
     &     iocc_op1op2(ngastp,2,njoined_op1op2),
     &     merge_op1(*), merge_op2(*), merge_op1op2(*), merge_op2op1(*),
     &     irst_op1(*), irst_op2(*), irst_op1op2(*),
     &     mstop1,mstop2,mstop1op2,
     &     igamtop1,igamtop2,igamtop1op2
      type(filinf), intent(in) ::
     &     ffop1,ffop2,ffop1op2
      type(operator), intent(in) ::
     &     op1, op2, op1op2
      type(strinf), intent(in) ::
     &     str_info
      type(strmapinf), intent(in) ::
     &     strmap_info
      type(orbinf), intent(in) ::
     &     orb_info

      integer ::
     &     nca_blk(2,6)
      integer, pointer ::
     &     cinfo_op1c(:,:),cinfo_op1a(:,:),
     &     cinfo_op2c(:,:),cinfo_op2a(:,:),
     &     cinfo_op1op2c(:,:),
     &     cinfo_op1op2a(:,:),
     &     cinfo_ex1c(:,:),cinfo_ex1a(:,:),
     &     cinfo_ex2c(:,:),cinfo_ex2a(:,:),
     &     cinfo_cntc(:,:),cinfo_cnta(:,:),
     &     map_info1c(:),
     &     map_info1a(:),
     &     map_info2c(:),
     &     map_info2a(:),
     &     map_info12c(:),
     &     map_info12a(:)
      
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
      case(2)
        ! reform to "condensed" description of contraction
        call get_num_subblk(nca_blk(1,1),nca_blk(2,1),
     &       iocc_op1,njoined_op1)
        call get_num_subblk(nca_blk(1,2),nca_blk(2,2),
     &       iocc_op2,njoined_op2)
        call get_num_subblk(nca_blk(1,3),nca_blk(2,3),iocc_op1op2,
     &       njoined_op1op2)
        call get_num_subblk(nca_blk(1,4),nca_blk(2,4),
     &       iocc_ex1,njoined_op1)
        call get_num_subblk(nca_blk(1,5),nca_blk(2,5),
     &       iocc_ex2,njoined_op2)
        call get_num_subblk(nca_blk(1,6),nca_blk(2,6),
     &       iocc_cnt,njoined_cnt)
        allocate(
     &       cinfo_op1c(nca_blk(1,1),3),cinfo_op1a(nca_blk(2,1),3),
     &       cinfo_op2c(nca_blk(1,2),3),cinfo_op2a(nca_blk(2,2),3),
     &       cinfo_op1op2c(nca_blk(1,3),3),
     &       cinfo_op1op2a(nca_blk(2,3),3),
     &       cinfo_ex1c(nca_blk(1,4),3),cinfo_ex1a(nca_blk(2,4),3),
     &       cinfo_ex2c(nca_blk(1,5),3),cinfo_ex2a(nca_blk(2,5),3),
     &       cinfo_cntc(nca_blk(1,6),3),cinfo_cnta(nca_blk(2,6),3))
        allocate(
     &       map_info1c(max(1,nca_blk(1,1)*2*
     &       (njoined_op1+njoined_cnt))),
     &       map_info1a(max(1,nca_blk(2,1)*2*
     &       (njoined_op1+njoined_cnt))),
     &       map_info2c(max(1,nca_blk(1,2)*2*
     &       (njoined_op2+njoined_cnt))),
     &       map_info2a(max(1,nca_blk(2,2)*2*
     &       (njoined_op2+njoined_cnt))),
     &         map_info12c(max(1,nca_blk(1,3)*2*
     &       (njoined_op1+njoined_op2))),
     &       map_info12a(max(1,nca_blk(2,3)*2*
     &       (njoined_op1+njoined_op2))))
        call condense_bc_info(
     &       cinfo_op1c, cinfo_op1a, cinfo_op2c, cinfo_op2a,
     &       cinfo_op1op2c, cinfo_op1op2a,
     &       cinfo_ex1c, cinfo_ex1a, cinfo_ex2c, cinfo_ex2a,
     &       cinfo_cntc, cinfo_cnta,
     &       map_info1c, map_info1a,
     &       map_info2c, map_info2a,
     &       map_info12c, map_info12a,
     &       nca_blk,
     &       iocc_op1, iocc_op2, iocc_op1op2,
     &       iocc_ex1,iocc_ex2,iocc_cnt,
     &       irst_op1, irst_op2, irst_op1op2,
     &       merge_op1, merge_op2, merge_op1op2, merge_op2op1,
     &       njoined_op1, njoined_op2,njoined_op1op2, njoined_cnt,
     &       str_info,orb_info%ihpvgas,orb_info%ngas)
        call contr_op1op2_wmaps_c(xfac,bc_sign,ffop1,ffop2,
     &     update,ffop1op2,xret,type_xret,
     &     op1,op2,op1op2,
     &     iblkop1,iblkop2,iblkop1op2,
     &     idoffop1,idoffop2,idoffop1op2,
     &       nca_blk,
     &       cinfo_op1c, cinfo_op1a, cinfo_op2c, cinfo_op2a,
     &       cinfo_op1op2c, cinfo_op1op2a,
     &       cinfo_ex1c, cinfo_ex1a, cinfo_ex2c, cinfo_ex2a,
     &       cinfo_cntc, cinfo_cnta,
     &       map_info1c, map_info1a,
     &       map_info2c, map_info2a,
     &       map_info12c, map_info12a,
     &     mstop1,mstop2,mstop1op2,
     &     igamtop1,igamtop2,igamtop1op2,
     &     str_info,strmap_info,orb_info)
        deallocate(
     &       cinfo_op1c,cinfo_op1a,
     &       cinfo_op2c,cinfo_op2a,
     &       cinfo_op1op2c,cinfo_op1op2a,
     &       cinfo_ex1c,cinfo_ex1a,
     &       cinfo_ex2c,cinfo_ex2a,
     &       cinfo_cntc,cinfo_cnta)
        deallocate(
     &       map_info1c,
     &       map_info1a,
     &       map_info2c,
     &       map_info2a,
     &       map_info12c,
     &       map_info12a)

      case default
        write(luout,*) 'contr_op1op2: route = ',irt_contr
        call quit(1,'contr_op1op2','route not implemented')
      end select
      call atim_cs(cpu,sys)
      cnt_op1op2(1) = cnt_op1op2(1)+cpu-cpu0
      cnt_op1op2(2) = cnt_op1op2(2)+sys-sys0

      return
      end
