*----------------------------------------------------------------------*
      subroutine contr_op1op2_wmaps_c(xfac,casign,
     &     update,xret,type_xret,
     &     me_op1,me_op2,me_op1op2,me_op1op2tmp,
     &     tra_op1, tra_op2, tra_op1op2,
     &     iblkop1,iblkop2,iblkop1op2,iblkop1op2tmp,
     &     idoffop1,idoffop2,idoffop1op2,
     &     cnt_info,reo_info,
     &     str_info,strmap_info,orb_info)
*----------------------------------------------------------------------*
*     slightly improved version:
*       full blocks incore
*       addressing via precalculated maps
*       mainly written for testing maps
*
*      update: update matrix elements on ffop1op2
*              else overwrite
*
*     andreas, feb 2007
*
*     further improved version:
*      trying to reduce overhead due to set_len_str
*
*     andreas, june 2016
*
*----------------------------------------------------------------------*
      implicit none

      include 'routes.h'
      include 'contr_times.h'

      include 'opdim.h'
      include 'stdunit.h'
      include 'ioparam.h'
      include 'multd2h.h'
      include 'def_operator.h'
      include 'def_me_list.h'
      include 'def_orbinf.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'def_filinf.h'
      include 'def_strmapinf.h'
      include 'def_reorder_info.h'
      include 'def_contraction_info.h'
      include 'ifc_memman.h'
      include 'ifc_operators.h'
      include 'hpvxseq.h'

      integer, parameter ::
     &     ntest = 5000

      logical, intent(in) ::
     &     update
      real(8), intent(in) ::
     &     xfac, casign
      real(8), intent(inout), target ::
     &     xret(1)
      type(contraction_info), target ::
     &     cnt_info
      integer, intent(in) ::
     &     type_xret,
     &     iblkop1, iblkop2, iblkop1op2, iblkop1op2tmp,
     &     idoffop1,idoffop2,idoffop1op2
      logical, intent(in) ::
     &     tra_op1, tra_op2, tra_op1op2
      type(me_list), intent(in) ::
     &     me_op1, me_op2, me_op1op2, me_op1op2tmp
      type(strinf), intent(in) ::
     &     str_info
      type(strmapinf), intent(inout) ::
     &     strmap_info
      type(orbinf), intent(in) ::
     &     orb_info
      type(reorder_info), intent(in), target ::
     &     reo_info

      logical ::
     &     bufop1, bufop2, bufop1op2,
     &     first1, first2, first3, first4, first5,
     &     ms_fix1, ms_fix2, ms_fix12, ms_fix_tmp,
     &     reject, fix_success1, fix_success2, fix_success12,
     &     reo_op1op2, nonzero, use_tr_here,
     &     tra_op1_, tra_op2_
      integer ::
     &     mstop1,mstop2,mstop1op2,
     &     igamtop1,igamtop2,igamtop1op2,
     &     nc_op1, na_op1, nc_op2, na_op2,
     &     nc_ex1, na_ex1, nc_ex2, na_ex2,
     &     nc_op1op2, na_op1op2,
     &     nc_op1op2tmp, na_op1op2tmp,
     &     nc_cnt, na_cnt, idx_restr,
     &     nsym, isym, ifree, lenscr, lenblock, lenbuf,
     &     buftyp1, buftyp2, buftyp12,
     &     idxst_op1, idxst_op2, idxst_op1op2,
     &     ioff_op1, ioff_op2, ioff_op1op2,
     &     idxop1, idxop2, idxop1op2,
     &     lenop1, lenop2, lenop1op2,
     &     idxms_op1op2_last, igam_op1op2_last,
     &     mscmx_a, mscmx_c, msc_ac, msc_a, msc_c,
     &     msex1_a, msex1_c, msex2_a, msex2_c,
     &     igamc_ac, igamc_a, igamc_c,
     &     igamex1_a, igamex1_c, igamex2_a, igamex2_c,
     &     idxms, idxdis, lenmap, lbuf_op1op2, lblk_op1op2tmp,
     &     idxdis_op1op2, idx, ii, maxidxms
      integer ::
     &     ncblk_op1, nablk_op1, ncblk_ex1, nablk_ex1,
     &     ncblk_op2, nablk_op2, ncblk_ex2, nablk_ex2,
     &     ncblk_op1op2, nablk_op1op2, ncblk_op1op2tmp, nablk_op1op2tmp,
     &     ncblk_cnt, nablk_cnt,
     &     ncblk_op1op2_0, nablk_op1op2_0,
     &     ncblk_reo12,    nablk_reo12,
     &     iblkoff
      type(filinf), pointer ::
     &     ffop1,ffop2,ffop1op2
      type(operator), pointer ::
     &     op1, op2, op1op2, op1op2tmp
      integer, pointer ::
     &     cinfo_op1c(:,:),cinfo_op1a(:,:),
     &     cinfo_op2c(:,:),cinfo_op2a(:,:),
     &     cinfo_op1op2c(:,:),
     &     cinfo_op1op2a(:,:),
     &     cinfo_op1op2tmpc(:,:),
     &     cinfo_op1op2tmpa(:,:),
     &     cinfo_ex1c(:,:),cinfo_ex1a(:,:),
     &     cinfo_ex2c(:,:),cinfo_ex2a(:,:),
     &     cinfo_cntc(:,:),cinfo_cnta(:,:),
     &     map_info_1c(:),
     &     map_info_1a(:),
     &     map_info_2c(:),
     &     map_info_2a(:),
     &     map_info_12c(:),
     &     map_info_12a(:),
     &     dmap_op1c(:),dmap_op1a(:),
     &     dmap_op2c(:),dmap_op2a(:),
     &     dmap_op1op2tmpc(:),dmap_op1op2tmpa(:)


      real(8) ::
     &     xnrm, fac_scal, fac_scal0, fac_ab, xret_last
      real(8) ::
     &     cpu, sys, cpu0, sys0, cpu00, sys00

      real(8), pointer ::
     &     xop1(:), xop2(:), xop1op2(:), xscr(:)
      real(8), pointer ::
     &     xbf1(:), xbf2(:), xbf12(:), xbf12tmp(:), xop1op2blk(:)

      integer ::
     &     msbnd(2,3), igambnd(2,3),
     &     ms12i_a(3), ms12i_c(3), igam12i_a(3), igam12i_c(3),
     &     igam12i_raw(3)

      integer, pointer ::
     &     gmop1dis_c(:), gmop1dis_a(:),
     &     gmop2dis_c(:), gmop2dis_a(:),
     &     gmex1dis_c(:), gmex1dis_a(:),
     &     gmex2dis_c(:), gmex2dis_a(:),
     &     gmc_dis_c (:), gmc_dis_a (:),
     &     gmi_dis_c (:), gmi_dis_a (:),
     &     msop1dis_c(:), msop1dis_a(:),
     &     msop2dis_c(:), msop2dis_a(:),
     &     msex1dis_c(:), msex1dis_a(:),
     &     msex2dis_c(:), msex2dis_a(:),
     &     msc_dis_c (:), msc_dis_a (:),
     &     msi_dis_c (:), msi_dis_a (:),
     &     idxmsop1dis_c(:), idxmsop1dis_a(:),
     &     idxmsop2dis_c(:), idxmsop2dis_a(:),
     &     idxmsex1dis_c(:), idxmsex1dis_a(:),
     &     idxmsex2dis_c(:), idxmsex2dis_a(:),
     &     idxmsc_dis_c (:), idxmsc_dis_a (:),
     &     idxmsi_dis_c (:), idxmsi_dis_a (:),
     &     lstrex1(:),lstrex2(:),lstrcnt(:),
     &     lstrop1(:),lstrop2(:),lstrop1op2tmp(:)

      integer, pointer ::
     &     map_ex1ex2a(:), map_ex1ex2c(:),
     &     map_ex1cnta(:), map_ex1cntc(:),
     &     map_ex2cnta(:), map_ex2cntc(:)

      integer, pointer ::
     &     ndis_op1(:,:), d_gam_ms_op1(:,:,:), gam_ms_op1(:,:),
     &     len_gam_ms_op1(:,:),
     &     ndis_op2(:,:), d_gam_ms_op2(:,:,:), gam_ms_op2(:,:),
     &     len_gam_ms_op2(:,:),
     &     ndis_op1op2tmp(:,:), d_gam_ms_op1op2(:,:,:),
     &     gam_ms_op1op2(:,:),
     &     len_gam_ms_op1op2(:,:),
     &     len_gam_ms_op1op2tmp(:,:), len_d_gam_ms_op1op2tmp(:,:,:)

      integer, pointer ::
     &     cinfo_reo12c(:,:), cinfo_reo12a(:,:),
     &     cinfo_op1op2_0c(:,:), cinfo_op1op2_0a(:,:),
     &     map_info_reo1c(:), map_info_reo1a(:),
     &     map_info_reo2c(:), map_info_reo2a(:),
     &     mapca12(:), diag_idx12(:), diag_ca12(:),
     &     mapca1(:), diag_idx1(:), diag_ca1(:),
     &     mapca2(:), diag_idx2(:), diag_ca2(:),
     &     lenstr_array(:,:,:)

c dbg
      integer, pointer ::
     &     dum1_c(:), dum1_a(:), hpvx1_c(:), hpvx1_a(:),
     &     dum2_c(:), dum2_a(:), hpvx2_c(:), hpvx2_a(:)
      integer ::
     &     msd1(ngastp,2,me_op1%op%njoined),
     &     msd2(ngastp,2,me_op2%op%njoined),
     &     jdx, tot_c, tot_a
c dbg

      type(graph), pointer ::
     &     graphs(:)

      integer, external ::
     &     ielsum, ielprd, idx_msgmdst2, get_lenmap, idxlist,
     &     max_dis_blk, msa2idxms4op
      logical, external ::
     &     next_dist, next_msgamdist2, list_cmp,
     &     nondia_blk, nondia_distr
      real(8), external ::
     &     ddot

      if (ntest.gt.0) then
        call write_title(lulog,wst_dbg_subr,
     &       'contr_op1op2_wmaps_c at work')
      end if

      op1 => me_op1%op
      op2 => me_op2%op
      op1op2 => me_op1op2%op
      op1op2tmp => me_op1op2tmp%op


      ffop1 => me_op1%fhand
      ffop2 => me_op2%fhand
      ffop1op2 => me_op1op2%fhand

      if (ntest.ge.10) then
        write(lulog,*) 'list1:   ',trim(me_op1%label),' transp:',tra_op1
        write(lulog,*) 'list2:   ',trim(me_op2%label),' transp:',tra_op2
        write(lulog,*) 'list12:  ',trim(me_op1op2%label),
     &                                             ' transp:',tra_op1op2
        write(lulog,*) 'ffop1:   ',ffop1%name(1:len_trim(ffop1%name))
        write(lulog,*) 'ffop2:   ',ffop2%name(1:len_trim(ffop2%name))
        write(lulog,*) 'ffop1op2:',
     &       ffop1op2%name(1:len_trim(ffop1op2%name))
        write(lulog,*) 'xfac = ',xfac
        write(lulog,*) 'casign = ',casign
        if (type_xret.ne.0)
     &       write(lulog,*) 'xret on entry = ',xret(1)
        write(lulog,*) 'op1: ',trim(op1%name),
     &       ' block ',iblkop1
        write(lulog,*) 'op2: ',trim(op2%name),
     &       ' block ',iblkop2
        if (iblkop1op2.gt.0) then
          write(lulog,*) 'op1op2: ',trim(op1op2%name),
     &       ' block ',iblkop1op2
        else
          write(lulog,*) 'op1op2: scalar'
        end if
      end if

      ! flag whether non-zero contribution to op1op2 occurred
      nonzero = .false.

      reo_op1op2 = reo_info%n_op_reo.gt.0
      if (reo_op1op2.and..not.associated(reo_info%map_reo1c))
     &     call quit(1,'contr_op1op2_wmaps_c',
     &     'reo_info is not consistent')
      if (ntest.ge.10) write(lulog,*) 'reo_op1op2: ',reo_op1op2
      if (ntest.ge.10.and.reo_op1op2)
     &            write(lulog,*) 'reo_sign: ',reo_info%sign_reo
      if (ntest.ge.10 .and. reo_op1op2) then
        write(lulog,*) 'op1op2tmp: ',trim(op1op2tmp%name),
     &       ' block ',iblkop1op2tmp
      end if

      ! set pointers to cnt_info entries
      ncblk_op1 = cnt_info%ncblk_op1 ! 1,1
      nablk_op1 = cnt_info%nablk_op1 ! 2,1
      ncblk_ex1 = cnt_info%ncblk_ex1 ! 1,4
      nablk_ex1 = cnt_info%nablk_ex1 ! 2,4
      ncblk_op2 = cnt_info%ncblk_op2 ! 1,2
      nablk_op2 = cnt_info%nablk_op2 ! 2,2
      ncblk_ex2 = cnt_info%ncblk_ex2 ! 1,5
      nablk_ex2 = cnt_info%nablk_ex2 ! 2,5
      ncblk_cnt = cnt_info%ncblk_cnt ! 1,6
      nablk_cnt = cnt_info%nablk_cnt ! 2,6
      ncblk_op1op2 = cnt_info%ncblk_op1op2 ! 1,3
      nablk_op1op2 = cnt_info%nablk_op1op2 ! 2,3
      ncblk_op1op2tmp = cnt_info%ncblk_op1op2tmp ! 1,7
      nablk_op1op2tmp = cnt_info%nablk_op1op2tmp ! 2,7
c dbg - call statistics
      st_total = st_total+1
      if (ncblk_ex1.eq.0) st_ex1c_z = st_ex1c_z+1
      if (nablk_ex1.eq.0) st_ex1a_z = st_ex1a_z+1
      if (ncblk_ex2.eq.0) st_ex2c_z = st_ex2c_z+1
      if (nablk_ex2.eq.0) st_ex2a_z = st_ex2a_z+1
      if (ncblk_cnt.eq.0) st_cntc_z = st_cntc_z+1
      if (nablk_cnt.eq.0) st_cnta_z = st_cnta_z+1
      if (ncblk_ex1.eq.0.and.nablk_ex1.eq.0.or.
     &    ncblk_ex1.eq.0.and.nablk_cnt.eq.0.or.
     &    ncblk_cnt.eq.0.and.nablk_ex1.eq.0.or.
     &    ncblk_cnt.eq.0.and.nablk_cnt.eq.0)
     &         st_noop1_reo = st_noop1_reo+1
      if (ncblk_ex2.eq.0.and.nablk_ex2.eq.0.or.
     &    ncblk_ex2.eq.0.and.ncblk_cnt.eq.0.or.
     &    nablk_cnt.eq.0.and.nablk_ex2.eq.0.or.
     &    ncblk_cnt.eq.0.and.nablk_cnt.eq.0)
     &         st_noop2_reo = st_noop2_reo+1
      if (ncblk_ex1.eq.0.and.nablk_ex1.eq.0.or.
     &    ncblk_ex1.eq.0.and.nablk_ex2.eq.0.or.
     &    ncblk_ex2.eq.0.and.nablk_ex1.eq.0.or.
     &    ncblk_ex2.eq.0.and.nablk_ex2.eq.0)
     &         st_noop12_reo = st_noop12_reo+1
c dbg

      cinfo_op1c => cnt_info%cinfo_op1c
      cinfo_op1a => cnt_info%cinfo_op1a
      cinfo_ex1c => cnt_info%cinfo_ex1c
      cinfo_ex1a => cnt_info%cinfo_ex1a
      cinfo_op2c => cnt_info%cinfo_op2c
      cinfo_op2a => cnt_info%cinfo_op2a
      cinfo_ex2c => cnt_info%cinfo_ex2c
      cinfo_ex2a => cnt_info%cinfo_ex2a
      cinfo_op1op2c => cnt_info%cinfo_op1op2c
      cinfo_op1op2a => cnt_info%cinfo_op1op2a
      cinfo_op1op2tmpc => cnt_info%cinfo_op1op2tmpc
      cinfo_op1op2tmpa => cnt_info%cinfo_op1op2tmpa
      cinfo_cntc => cnt_info%cinfo_cntc
      cinfo_cnta => cnt_info%cinfo_cnta

c dbg
c      print *,'Print 1'
c      print *,'Create ',cinfo_op1op2tmpc
c      print *,'Annihilate ',cinfo_op1op2tmpa
c dbg

c dbg
c      print *,'-------------------------------------------------------'
c      print *,'       graph information table'
c      print *,'-------------------------------------------------------'
c      print *,'graphs for op1c: ',cinfo_op1c(1:ncblk_op1,2)
c      print *,'graphs for op1a: ',cinfo_op1a(1:nablk_op1,2)
c      print *
c      print *,'graphs for ex1c: ',cinfo_ex1c(1:ncblk_ex1,2)
c      print *,'graphs for ex1a: ',cinfo_ex1a(1:nablk_ex1,2)
c      print *,'graphs for cntc: ',cinfo_cntc(1:ncblk_cnt,2)
c      print *,'graphs for cnta: ',cinfo_cnta(1:nablk_cnt,2)
c      print *
c      print *
c      print *,'graphs for op2c: ',cinfo_op2c(1:ncblk_op2,2)
c      print *,'graphs for op2a: ',cinfo_op2a(1:nablk_op2,2)
c      print *
c      print *,'graphs for ex2c: ',cinfo_ex2c(1:ncblk_ex2,2)
c      print *,'graphs for ex2a: ',cinfo_ex2a(1:nablk_ex2,2)
c      print *
c      print *
c      print *,'graphs for op1op2c: ',cinfo_op1op2c(1:ncblk_op1op2,2)
c      print *,'graphs for op1op2a: ',cinfo_op1op2a(1:nablk_op1op2,2)
c      print *
c      print *,'-------------------------------------------------------'
c dbg

      map_info_1c => cnt_info%map_info_1c
      map_info_1a => cnt_info%map_info_1a
      map_info_2c => cnt_info%map_info_2c
      map_info_2a => cnt_info%map_info_2a
      map_info_12c => cnt_info%map_info_12c
      map_info_12a => cnt_info%map_info_12a

      allocate(
     &       gmop1dis_c(ncblk_op1), gmop1dis_a(nablk_op1),
     &       gmop2dis_c(ncblk_op2), gmop2dis_a(nablk_op2),
     &       gmex1dis_c(ncblk_ex1), gmex1dis_a(nablk_ex1),
     &       gmex2dis_c(ncblk_ex2), gmex2dis_a(nablk_ex2),
     &       gmc_dis_c(ncblk_cnt), gmc_dis_a(nablk_cnt),
     &       gmi_dis_c(ncblk_op1op2tmp), gmi_dis_a(nablk_op1op2tmp),
     &       msop1dis_c(ncblk_op1), msop1dis_a(nablk_op1),
     &       msop2dis_c(ncblk_op2), msop2dis_a(nablk_op2),
     &       msex1dis_c(ncblk_ex1), msex1dis_a(nablk_ex1),
     &       msex2dis_c(ncblk_ex2), msex2dis_a(nablk_ex2),
     &       msc_dis_c(ncblk_cnt), msc_dis_a(nablk_cnt),
     &       msi_dis_c(ncblk_op1op2tmp), msi_dis_a(nablk_op1op2tmp),
     &       idxmsop1dis_c(ncblk_op1), idxmsop1dis_a(nablk_op1),
     &       idxmsop2dis_c(ncblk_op2), idxmsop2dis_a(nablk_op2),
     &       idxmsex1dis_c(ncblk_ex1), idxmsex1dis_a(nablk_ex1),
     &       idxmsex2dis_c(ncblk_ex2), idxmsex2dis_a(nablk_ex2),
     &       idxmsc_dis_c(ncblk_cnt), idxmsc_dis_a(nablk_cnt),
     &       idxmsi_dis_c(ncblk_op1op2tmp),
     &       idxmsi_dis_a(nablk_op1op2tmp),
     &       lstrex1(ncblk_ex1+nablk_ex1),
     &       lstrex2(ncblk_ex2+nablk_ex2),
     &       lstrcnt(ncblk_cnt+nablk_cnt),
     &       lstrop1(ncblk_op1+nablk_op1),
     &       lstrop2(ncblk_op2+nablk_op2),
     &       lstrop1op2tmp(ncblk_op1op2tmp+nablk_op1op2tmp)
     &       )
      allocate(dmap_op1c(ncblk_op1),dmap_op1a(nablk_op1),
     &         dmap_op2c(ncblk_op2),dmap_op2a(nablk_op2),
     &         dmap_op1op2tmpc(ncblk_op1op2tmp),
     &         dmap_op1op2tmpa(nablk_op1op2tmp))

      call set_dis_tra_map(dmap_op1c,dmap_op1a,
     &     cinfo_op1c(1,3),cinfo_op1a(1,3),ncblk_op1,nablk_op1)
      call set_dis_tra_map(dmap_op2c,dmap_op2a,
     &     cinfo_op2c(1,3),cinfo_op2a(1,3),ncblk_op2,nablk_op2)
      call set_dis_tra_map(dmap_op1op2tmpc,dmap_op1op2tmpa,
     &     cinfo_op1op2tmpc(1,3),cinfo_op1op2tmpa(1,3),
     &     ncblk_op1op2tmp,nablk_op1op2tmp)

      idxst_op1 = me_op1%off_op_occ(iblkop1) + 1
      lenop1    = me_op1%len_op_occ(iblkop1)
      idxst_op2 = me_op2%off_op_occ(iblkop2) + 1
      lenop2    = me_op2%len_op_occ(iblkop2)
      if (iblkop1op2.gt.0) then
        ! refers to reordered op1op2 (if that makes a difference)
        idxst_op1op2 = me_op1op2%off_op_occ(iblkop1op2) + 1
        lenop1op2    = me_op1op2%len_op_occ(iblkop1op2)
      else
        idxst_op1op2 = 1
        lenop1op2 = 1
      end if

      mstop1 = me_op1%mst
      mstop2 = me_op2%mst
      mstop1op2 = me_op1op2%mst
      if (tra_op1) mstop1 = -mstop1
      if (tra_op2) mstop2 = -mstop2
      if (tra_op1op2) mstop1op2 = -mstop1op2
      igamtop1 = me_op1%gamt
      igamtop2 = me_op2%gamt
      igamtop1op2 = me_op1op2%gamt

      ! See which, if any, operators need to have their Ms values fixed.
      ms_fix1 = me_op1%fix_vertex_ms
      ms_fix2 = me_op2%fix_vertex_ms
      ms_fix12 = me_op1op2%fix_vertex_ms
      ms_fix_tmp = me_op1op2tmp%fix_vertex_ms

      if(ms_fix1) allocate(dum1_c(ncblk_op1),dum1_a(nablk_op1),
     &       hpvx1_c(ncblk_op1),hpvx1_a(nablk_op1))
      if(ms_fix2) allocate(dum2_c(ncblk_op2),dum2_a(nablk_op2),
     &       hpvx2_c(ncblk_op2),hpvx2_a(nablk_op2))

c dbg
c      print *,'op1 name, fix ', trim(me_op1%op%name),ms_fix1
c      print *,'op2 name, fix ', trim(me_op2%op%name),ms_fix2
c      print *,'op12 name, fix ', trim(me_op1op2%op%name),ms_fix12
c      print *,'optmp name, fix ', trim(me_op1op2tmp%op%name),ms_fix_tmp
c dbg

      ! use spin flip symmetry: only if all operators have a distinct
      ! alpha/beta symmetry
      use_tr_here = use_tr       .and.
     &     me_op1%absym   .ne.0  .and.
     &     me_op2%absym   .ne.0  .and.
     &     me_op1op2%absym.ne.0  .and.
     &     me_op1%absym*me_op2%absym.eq.me_op1op2%absym
      ! the last fix is to avoid a bug with spin densities, which are
      ! obtained by anti-symmetrizing the resulting ME list; in this
      ! case it can happen, that the symmetries do not match (can happen
      ! for all states with MS=0)

      fac_scal = 1d0

      if (multd2h(igamtop1,igamtop2).ne.igamtop1op2) then
        write(lulog,*) ' 1: gamma=',igamtop1,' op=',trim(op1%name),
     &                                       ' list=',trim(me_op1%label)
        write(lulog,*) ' 2: gamma=',igamtop2,' op=',trim(op2%name),
     &                                       ' list=',trim(me_op2%label)
        write(lulog,*)' R: gamma=',igamtop1op2,' op=',trim(op1op2%name),
     &                                    ' list=',trim(me_op1op2%label)
        call quit(1,'contr_op1op2_wmaps_c','inconsistent symmetries')
      end if

      if (me_op1%op%formal_blk(iblkop1).or.
c     &    me_op2%op%formal_blk(iblkop2).or.
c     &    me_op1op2%op%formal_blk(iblkop1op2)) then
     &    me_op2%op%formal_blk(iblkop2)) then
        write(lulog,*) me_op1%op%formal_blk(iblkop1),
     &                 me_op2%op%formal_blk(iblkop2),
     &                 me_op1op2%op%formal_blk(iblkop1op2)
        write(lulog,*) 'op1: ',trim(op1%name),
     &       ' block ',iblkop1
        write(lulog,*) 'op2: ',trim(op2%name),
     &       ' block ',iblkop2
        if (iblkop1op2.gt.0) then
          write(lulog,*) 'op1op2: ',trim(op1op2%name),
     &       ' block ',iblkop1op2
        else
          write(lulog,*) 'op1op2: scalar'
        end if

        call quit(1,'contr_op1op2_wmaps_c','called for formal block')
      else if (me_op1op2%op%formal_blk(iblkop1op2)) then
        ! do nothing for formal block
        if(ms_fix1) deallocate(dum1_c,dum1_a,hpvx1_c,hpvx1_a)
        if(ms_fix2) deallocate(dum2_c,dum2_a,hpvx2_c,hpvx2_a)

        deallocate(
     &       gmop1dis_c, gmop1dis_a,
     &       gmop2dis_c, gmop2dis_a,
     &       gmex1dis_c, gmex1dis_a,
     &       gmex2dis_c, gmex2dis_a,
     &       gmc_dis_c , gmc_dis_a ,
     &       gmi_dis_c , gmi_dis_a ,
     &       msop1dis_c, msop1dis_a,
     &       msop2dis_c, msop2dis_a,
     &       msex1dis_c, msex1dis_a,
     &       msex2dis_c, msex2dis_a,
     &       msc_dis_c , msc_dis_a ,
     &       msi_dis_c , msi_dis_a ,
     &       idxmsop1dis_c, idxmsop1dis_a,
     &       idxmsop2dis_c, idxmsop2dis_a,
     &       idxmsex1dis_c, idxmsex1dis_a,
     &       idxmsex2dis_c, idxmsex2dis_a,
     &       idxmsc_dis_c , idxmsc_dis_a ,
     &       idxmsi_dis_c , idxmsi_dis_a ,
     &       lstrex1,lstrex2,lstrcnt,
     &       lstrop1,lstrop2,lstrop1op2tmp)
        deallocate(dmap_op1c,dmap_op1a,
     &             dmap_op2c,dmap_op2a,
     &             dmap_op1op2tmpc,dmap_op1op2tmpa)
        return
      end if

      ! skip if result operator has zero length
      if (lenop1op2.eq.0) then
         if(ms_fix1) deallocate(dum1_c,dum1_a,hpvx1_c,hpvx1_a)
         if(ms_fix2) deallocate(dum2_c,dum2_a,hpvx2_c,hpvx2_a)

         deallocate(
     &        gmop1dis_c, gmop1dis_a,
     &        gmop2dis_c, gmop2dis_a,
     &        gmex1dis_c, gmex1dis_a,
     &        gmex2dis_c, gmex2dis_a,
     &        gmc_dis_c , gmc_dis_a ,
     &        gmi_dis_c , gmi_dis_a ,
     &        msop1dis_c, msop1dis_a,
     &        msop2dis_c, msop2dis_a,
     &        msex1dis_c, msex1dis_a,
     &        msex2dis_c, msex2dis_a,
     &        msc_dis_c , msc_dis_a ,
     &        msi_dis_c , msi_dis_a ,
     &        idxmsop1dis_c, idxmsop1dis_a,
     &        idxmsop2dis_c, idxmsop2dis_a,
     &        idxmsex1dis_c, idxmsex1dis_a,
     &        idxmsex2dis_c, idxmsex2dis_a,
     &        idxmsc_dis_c , idxmsc_dis_a ,
     &        idxmsi_dis_c , idxmsi_dis_a ,
     &        lstrex1,lstrex2,lstrcnt,
     &        lstrop1,lstrop2,lstrop1op2tmp)
         deallocate(dmap_op1c,dmap_op1a,
     &              dmap_op2c,dmap_op2a,
     &              dmap_op1op2tmpc,dmap_op1op2tmpa)
         return
      end if

c      ! we accept that certain non-totally symmetric operator blocks
c      ! may have zero length ...
c      if (lenop1.eq.0.and.igamtop1.ne.1 .or.
c     &    lenop2.eq.0.and.igamtop2.ne.1 .or.
c     &    lenop1op2.eq.0.and.igamtop1op2.ne.1) return
c
c      ! ... else we raise an error flag; check whether this is
c      ! useful or not ....
c      if (lenop1.le.0.or.lenop2.le.0.or.lenop1op2.le.0) then
c        write(lulog,*)
c     &       trim(op1%name),' ',
c     &       trim(op2%name),' ',
c     &       trim(op1op2%name)
c        write(lulog,*) 'lenop1, lenop2, lenop1op2: ',
c     &                  lenop1, lenop2, lenop1op2
c        call quit(1,'contr_op1op2_wmaps_c',
c     &     'zero length for operator?')
c      end if

      ifree = mem_setmark('contr1')

c      call atim_cs(cpu0,sys0)

      ! average scratch for each operator:
      !    about 1/6th of the available remaining core
      lenscr = ifree/6

      if (ffop1%buffered.and.ffop1%incore(iblkop1).ge.0) then
        buftyp1 = 0
        bufop1 = .true.
        xop1 => ffop1%buffer(idxst_op1:)
        ioff_op1 = idxst_op1-1
      else
        bufop1 = .false.
        ! LOWER incore requirements:
        ! check for length of operator
        call set_op_scratch(lenbuf,buftyp1,me_op1,iblkop1,
     &       lenscr,orb_info)
        ifree = mem_alloc_real(xbf1,lenbuf,'xbf1')
c        ifree = mem_alloc_real(xbf1,lenop1,'xbf1')
        xop1 => xbf1
        if (buftyp1.eq.0) then
          ioff_op1 = idxst_op1-1
          call get_vec(ffop1,xop1,idoffop1+idxst_op1,
     &                          idoffop1+idxst_op1-1+lenop1)
        end if
      end if
      if (ffop2%buffered.and.ffop2%incore(iblkop2).ge.0) then
        buftyp2 = 0
        bufop2 = .true.
        xop2 => ffop2%buffer(idxst_op2:)
        ioff_op2 = idxst_op2-1
      else
        bufop2 = .false.
        ! LOWER incore requirements:
        ! see above
        call set_op_scratch(lenbuf,buftyp2,me_op2,iblkop2,
     &       lenscr,orb_info)
        ifree = mem_alloc_real(xbf2,lenbuf,'xbf2')
c        ifree = mem_alloc_real(xbf2,lenop2,'xbf2')
        xop2 => xbf2
        if (buftyp2.eq.0) then
          ioff_op2 = idxst_op2-1
          call get_vec(ffop2,xop2,idoffop2+idxst_op2,
     &                          idoffop2+idxst_op2-1+lenop2)
        end if
      end if

      if (ntest.ge.100) write(lulog,*) ' bufop1/2: ',bufop1,bufop2

      ! get result vector as well (as we update)
      ! refers to reordered op1op2
      if (iblkop1op2.gt.0) then
        if (ffop1op2%buffered.and.ffop1op2%incore(iblkop1op2).ge.0) then
          buftyp12 = 0
          bufop1op2 = .true.
          xop1op2 => ffop1op2%buffer(idxst_op1op2:)
          ioff_op1op2 = idxst_op1op2-1
          if (.not.update) xop1op2(1:lenop1op2) = 0d0
        else
          bufop1op2 = .false.
          ! LOWER incore requirements:
          ! see above
          ! decide on buftyp12 here
          if (.not.reo_op1op2) then
            lenscr = 2*ifree/3
            call set_op_scratch(lenbuf,buftyp12,me_op1op2,iblkop1op2,
     &         lenscr,orb_info)
          else
            ! for reordered operators, we currently only can do:
            buftyp12 = 0
            lenbuf = lenop1op2
          end if
          ! presently: only buftyp12=0/1
c          if (buftyp12.gt.1) then
c            call warn('contr_op1op2_wmaps_c','setting buftyp12 to 1')
c            buftyp12 = 1
c            lenbuf = max_dis_blk(-1,me_op1op2,iblkop1op2,orb_info)
c          end if

          ! unset use_tr_here, if buftyp12!=0
          use_tr_here = use_tr_here.and.buftyp12.eq.0

          ifree = mem_alloc_real(xbf12,lenbuf,'xbf12')
          lbuf_op1op2 = lenbuf
          xop1op2 => xbf12
          if (buftyp12.eq.0) then
            ioff_op1op2 = idxst_op1op2-1
            if (update) then
              ! read from disc
              call get_vec(ffop1op2,xop1op2,idoffop1op2+idxst_op1op2,
     &                             idoffop1op2+idxst_op1op2-1+lenop1op2)
            else
              ! init with zero
              xop1op2(1:lenop1op2) = 0d0
            end if
          end if
        end if
        if (ntest.ge.100) write(lulog,*) ' bufop1op2: ',bufop1op2
      else
        bufop1op2 = .true.
        xop1op2 => xret
        if (ntest.ge.100) write(lulog,*) ' result is scalar '
      end if

      if (me_op1op2tmp%diag_type.ne.0) then
        allocate(mapca12(ncblk_op1op2tmp),diag_idx12(ncblk_op1op2tmp),
     &           diag_ca12(ncblk_op1op2tmp))
        if (nondia_blk(mapca12,diag_idx12,diag_ca12,
     &                 op1op2tmp%ihpvca_occ(1,1,
     &                        (iblkop1op2tmp-1)*op1op2tmp%njoined+1),
     &                 op1op2tmp%njoined,ncblk_op1op2tmp,
     &                 nablk_op1op2tmp,me_op1op2tmp%diag_type))
     &       call quit(1,'contr_op1op2_wmaps_c','non-diagonal block!')
      end if
      if (me_op1%diag_type.ne.0) then
        allocate(mapca1(ncblk_op1),diag_idx1(ncblk_op1),
     &           diag_ca1(ncblk_op1))
        if (nondia_blk(mapca1,diag_idx1,diag_ca1,
     &                 op1%ihpvca_occ(1,1,
     &                        (iblkop1-1)*op1%njoined+1),
     &                 op1%njoined,ncblk_op1,
     &                 nablk_op1,me_op1%diag_type))
     &       call quit(1,'contr_op1op2_wmaps_c','non-diagonal block!')
      end if
      if (me_op2%diag_type.ne.0) then
        allocate(mapca2(ncblk_op2),diag_idx2(ncblk_op2),
     &           diag_ca2(ncblk_op2))
        if (nondia_blk(mapca2,diag_idx2,diag_ca2,
     &                 op2%ihpvca_occ(1,1,
     &                        (iblkop2-1)*op2%njoined+1),
     &                 op2%njoined,ncblk_op2,
     &                 nablk_op2,me_op2%diag_type))
     &       call quit(1,'contr_op1op2_wmaps_c','non-diagonal block!')
      end if

c      call atim_cs(cpu,sys)
c      cnt_rd(1) = cnt_rd(1) + cpu-cpu0
c      cnt_rd(2) = cnt_rd(2) + sys-sys0

      if (ntest.ge.1000) then
        ! this will work if all blocks incore, only:
        write(lulog,*) 'operator 1 (',trim(op1%name),
     &                    ',list=',trim(me_op1%label),')'
        if (buftyp1.eq.0) then
          write(lulog,*) 'full list loaded:'
          call wrt_mel_buf(lulog,5,xop1,me_op1,iblkop1,iblkop1,
     &                  str_info,orb_info)
        else
          write(lulog,*) 'buftyp1 = ',buftyp1
          if (ntest.lt.10000) then
            write(lulog,*)
     &         'complete list not available for printout'
          else
            call wrt_mel_file(lulog,5,me_op1,iblkop1,iblkop1,
     &           str_info,orb_info)
          end if
        end if
        write(lulog,*) 'operator 2 (',trim(op2%name),
     &                    ',list=',trim(me_op2%label),')'
        if (buftyp2.eq.0) then
          write(lulog,*) 'full list loaded:'
          call wrt_mel_buf(lulog,5,xop2,me_op2,iblkop2,iblkop2,
     &                  str_info,orb_info)
        else
          write(lulog,*) 'buftyp2 = ',buftyp2
          if (ntest.lt.10000) then
            write(lulog,*)
     &         'complete list not available for printout'
          else
            call wrt_mel_file(lulog,5,me_op2,iblkop2,iblkop2,
     &           str_info,orb_info)
          end if
        end if

        if (iblkop1op2.gt.0) then
          write(lulog,*) 'operator 12 on entry (',trim(op1op2%name),
     &                                ',list=',trim(me_op1op2%label),')'

          if (buftyp12.eq.0) then
            write(lulog,*) 'full list loaded:'
            call wrt_mel_buf(lulog,5,xop1op2,me_op1op2,
     &                    iblkop1op2,iblkop1op2,
     &                    str_info,orb_info)
          else
            write(lulog,*) 'buftyp12 = ',buftyp12
            if (ntest.lt.10000) then
              write(lulog,*)
     &         'complete list not available for printout'
            else
              call wrt_mel_file(lulog,5,me_op1op2,iblkop1op2,iblkop1op2,
     &           str_info,orb_info)
            end if
          end if

        end if
      end if

      nsym = orb_info%nsym
      ! store info on graphs in more conventional array
      ! to allow calling of set_len_str2
      maxidxms =  str_info%max_idxms
      allocate(lenstr_array(nsym,str_info%max_idxms,str_info%ngraph))
      call set_lenstr_array(lenstr_array,nsym,
     &                      str_info%max_idxms,str_info)

      ! OLD graphs => str_info%g

      ndis_op1 => me_op1%off_op_gmox(iblkop1)%ndis
      gam_ms_op1 => me_op1%off_op_gmo(iblkop1)%gam_ms
      len_gam_ms_op1 => me_op1%len_op_gmo(iblkop1)%gam_ms
      d_gam_ms_op1 => me_op1%off_op_gmox(iblkop1)%d_gam_ms
      ndis_op2 => me_op2%off_op_gmox(iblkop2)%ndis
      gam_ms_op2 => me_op2%off_op_gmo(iblkop2)%gam_ms
      len_gam_ms_op2 => me_op2%len_op_gmo(iblkop2)%gam_ms
      d_gam_ms_op2 => me_op2%off_op_gmox(iblkop2)%d_gam_ms
      ndis_op1op2tmp => me_op1op2tmp%off_op_gmox(iblkop1op2tmp)%ndis
      gam_ms_op1op2 => me_op1op2%off_op_gmo(iblkop1op2)%gam_ms
      len_gam_ms_op1op2 => me_op1op2%len_op_gmo(iblkop1op2)%gam_ms
      len_gam_ms_op1op2tmp =>
     &                   me_op1op2tmp%len_op_gmo(iblkop1op2tmp)%gam_ms
      d_gam_ms_op1op2 => me_op1op2%off_op_gmox(iblkop1op2)%d_gam_ms
      len_d_gam_ms_op1op2tmp =>
     &                  me_op1op2tmp%len_op_gmox(iblkop1op2tmp)%d_gam_ms

      call sum_occ(nc_op1,cinfo_op1c,ncblk_op1)
      call sum_occ(na_op1,cinfo_op1a,nablk_op1)
      call sum_occ(nc_op2,cinfo_op2c,ncblk_op2)
      call sum_occ(na_op2,cinfo_op2a,nablk_op2)
      call sum_occ(nc_op1op2,cinfo_op1op2c,ncblk_op1op2)
      call sum_occ(na_op1op2,cinfo_op1op2a,nablk_op1op2)
      call sum_occ(nc_op1op2tmp,cinfo_op1op2tmpc,ncblk_op1op2tmp)
      call sum_occ(na_op1op2tmp,cinfo_op1op2tmpa,nablk_op1op2tmp)

c dbg
c      print *,'Print 2'
c      print *,'Create ',cinfo_op1op2tmpc
c      print *,'Annihilate ',cinfo_op1op2tmpa
c dbg

c dbg
      if (na_op1op2.ne.na_op1op2tmp)
     &     call quit(1,'contr_op1op2_wmaps_c','unexpected 1a')
      if (nc_op1op2.ne.nc_op1op2tmp) then
        write(lulog,*) 'OP1OP2 (C)   : ',nc_op1op2,
     &       ' <- ',cinfo_op1op2c(1:ncblk_op1op2,1)
        write(lulog,*) 'OP1OP2TMP (C): ',nc_op1op2tmp,
     &       ' <- ',cinfo_op1op2tmpc(1:ncblk_op1op2tmp,1)
        call quit(1,'contr_op1op2_wmaps_c','unexpected 1b')
      end if
c dbg
      call sum_occ(nc_ex1,cinfo_ex1c,ncblk_ex1)
      call sum_occ(na_ex1,cinfo_ex1a,nablk_ex1)
      call sum_occ(nc_ex2,cinfo_ex2c,ncblk_ex2)
      call sum_occ(na_ex2,cinfo_ex2a,nablk_ex2)
      call sum_occ(nc_cnt,cinfo_cntc,ncblk_cnt)
      call sum_occ(na_cnt,cinfo_cnta,nablk_cnt)

      tra_op1_ = tra_op1
c     a) incorrect for cases with X
c     b) no impact on run-time (needs to be re-checked)
c      ! decide, whether we might use the transposed operator
c      if (op1%hermitian.eq.1.and.ncblk_op1.eq.nablk_op1) then
c        tra_op1_ = tra_op1.xor.
c     &       (list_cmp(cinfo_op1c,cinfo_op1a,ncblk_op1).and.
c     &        na_cnt.gt.nc_cnt)
c      end if
      tra_op2_ = tra_op2
c      if (op2%hermitian.eq.1.and.ncblk_op2.eq.nablk_op2) then
c        tra_op2_ = tra_op2.xor.
c     &       (list_cmp(cinfo_op2c,cinfo_op2a,ncblk_op2).and.
c     &        nc_cnt.gt.na_cnt)
c      end if
cc dbg
c      if (tra_op1.neqv.tra_op1_) write(lustat,*) 'it happened 1'
c      if (tra_op2.neqv.tra_op2_) write(lustat,*) 'it happened 2'
cc dbg

      ! set up maps (if necessary)
      call strmap_man_c(1,lenmap,
     &     cinfo_ex1c(1,2),ncblk_ex1,
     &     cinfo_ex2c(1,2),ncblk_ex2,
     &     cinfo_op1op2tmpc(1,2),ncblk_op1op2tmp,map_info_12c,
     &     str_info,strmap_info,orb_info)
      ifree = mem_alloc_int(map_ex1ex2c,lenmap,'exmap_c')
      call strmap_man_c(1,lenmap,
     &     cinfo_ex2a(1,2),nablk_ex2,
     &     cinfo_ex1a(1,2),nablk_ex1,
     &     cinfo_op1op2tmpa(1,2),nablk_op1op2tmp,map_info_12a,
     &     str_info,strmap_info,orb_info)
      ifree = mem_alloc_int(map_ex1ex2a,lenmap,'exmap_a')
      call strmap_man_c(1,lenmap,
     &     cinfo_cntc(1,2),ncblk_cnt,
     &     cinfo_ex1c(1,2),ncblk_ex1,
     &     cinfo_op1c(1,2),ncblk_op1,map_info_1c,
     &     str_info,strmap_info,orb_info)
      ifree = mem_alloc_int(map_ex1cntc,lenmap,'op1map_c')
      call strmap_man_c(1,lenmap,
     &     cinfo_cnta(1,2),nablk_cnt,
     &     cinfo_ex1a(1,2),nablk_ex1,
     &     cinfo_op1a(1,2),nablk_op1,map_info_1a,
     &     str_info,strmap_info,orb_info)
      ifree = mem_alloc_int(map_ex1cnta,lenmap,'op1map_a')
      call strmap_man_c(1,lenmap,
     &     cinfo_cnta(1,2),nablk_cnt,
     &     cinfo_ex2c(1,2),ncblk_ex2,
     &     cinfo_op2c(1,2),ncblk_op2,map_info_2c,
     &     str_info,strmap_info,orb_info)
      ifree = mem_alloc_int(map_ex2cntc,lenmap,'op2map_c')
      call strmap_man_c(1,lenmap,
     &     cinfo_cntc(1,2),ncblk_cnt,
     &     cinfo_ex2a(1,2),nablk_ex2,
     &     cinfo_op2a(1,2),nablk_op2,map_info_2a,
     &     str_info,strmap_info,orb_info)
      ifree = mem_alloc_int(map_ex2cnta,lenmap,'op2map_a')

      if (reo_op1op2) then
        cinfo_reo12c => reo_info%cinfo_reo_c
        cinfo_reo12a => reo_info%cinfo_reo_a
        cinfo_op1op2_0c => reo_info%cinfo_opreo0c
        cinfo_op1op2_0a => reo_info%cinfo_opreo0a
        map_info_reo1c  => reo_info%map_reo1c
        map_info_reo1a  => reo_info%map_reo1a
        map_info_reo2c  => reo_info%map_reo2c
        map_info_reo2a  => reo_info%map_reo2a
        ncblk_op1op2_0  = reo_info%ncblk_reo0
        nablk_op1op2_0  = reo_info%nablk_reo0
        ncblk_reo12  = reo_info%ncblk_reo
        nablk_reo12  = reo_info%nablk_reo
        call strmap_man_c(2,lenmap,
     &     cinfo_op1op2_0c(1,2),ncblk_op1op2_0,
     &     cinfo_reo12c(1,2),ncblk_reo12,
     &     cinfo_op1op2tmpc(1,2),ncblk_op1op2tmp,map_info_reo1c,
     &     str_info,strmap_info,orb_info)
        call strmap_man_c(2,lenmap,
     &     cinfo_reo12a(1,2),nablk_reo12,
     &     cinfo_op1op2_0a(1,2),nablk_op1op2_0,
     &     cinfo_op1op2tmpa(1,2),nablk_op1op2tmp,map_info_reo1a,
     &     str_info,strmap_info,orb_info)
        call strmap_man_c(2,lenmap,
     &     cinfo_op1op2_0c(1,2),ncblk_op1op2_0,
     &     cinfo_reo12c(1,2),ncblk_reo12,
     &     cinfo_op1op2c(1,2),ncblk_op1op2,map_info_reo2c,
     &     str_info,strmap_info,orb_info)
        call strmap_man_c(2,lenmap,
     &     cinfo_reo12a(1,2),nablk_reo12,
     &     cinfo_op1op2_0a(1,2),nablk_op1op2_0,
     &     cinfo_op1op2a(1,2),nablk_op1op2,map_info_reo2a,
     &     str_info,strmap_info,orb_info)
      end if

      if (reo_op1op2) then
        lblk_op1op2tmp =
     &       max_dis_blk(0,me_op1op2tmp,iblkop1op2tmp,orb_info)
        ifree = mem_alloc_real(xbf12tmp,lblk_op1op2tmp,'xbf12tmp')
      end if

      if (irt_contr.gt.2) then
        ! preliminary: a wild guess
        lenblock = len_str_block
        lenscr = max(3*lenblock*lenblock,
     &        (max_dis_blk(0,me_op1op2tmp,iblkop1op2tmp,orb_info)
     &         + max_dis_blk(0,me_op1,iblkop1,orb_info)
     &         + max_dis_blk(0,me_op2,iblkop2,orb_info))*lenblock)
        ! but do at most request 90 per cent of remaining core
        lenscr = (ifree*9)/10
        ifree = mem_alloc_real(xscr,lenscr,'contr_scr')
      end if

      ! minimum Ms(A) for ...
      msbnd(1,1) = -na_op1 ! operator 1
      msbnd(1,2) = -na_op2 ! operator 2
      msbnd(1,3) = -na_op1op2 ! product
      ! if use_tr is active, restrict the operator with larger max MS:
      idx_restr = 1
      if (na_op2.gt.na_op1) idx_restr = 2
      ! maximum Ms(A) for ...
      msbnd(2,1) = -msbnd(1,1)
      msbnd(2,2) = -msbnd(1,2)
      msbnd(2,3) = -msbnd(1,3)
      ! max |Ms| for ...
      mscmx_a = na_cnt ! C(A)
      mscmx_c = nc_cnt ! C(C)
      ! minimum IRREP for operators
      igambnd(1,1) = 1
      igambnd(1,2) = 1
      igambnd(1,3) = 1
      ! maximum IRREP
      igambnd(2,1) = nsym
      igambnd(2,2) = nsym
      igambnd(2,3) = nsym
      ! loop Ms-cases of (Op1(A),Op2(A),Op1Op2(A))
      first1 = .true.
      idxms_op1op2_last = -1
      if (type_xret.ne.0.and.buftyp12.ne.0) then
        xret(1) = 0d0
        xret_last = 0d0
      end if
      ms_loop: do
        if (first1) then
          first1 = .false.
          ! initial Ms distribution
          ms12i_a(1:3) = msbnd(2,1:3)
        else
          ! next Ms distribution
          if (.not.next_dist(ms12i_a,3,msbnd,-2)) exit
        end if

        ms12i_c(1) = ms12i_a(1) + mstop1
        ms12i_c(2) = ms12i_a(2) + mstop2
        ms12i_c(3) = ms12i_a(3) + mstop1op2
        if (abs(ms12i_c(1)).gt.nc_op1) cycle ms_loop
        if (abs(ms12i_c(2)).gt.nc_op2) cycle ms_loop
        if (abs(ms12i_c(3)).gt.nc_op1op2) cycle ms_loop

        msc_ac = ms12i_a(1) + ms12i_a(2) - ms12i_a(3)

        if (mscmx_a+mscmx_c.lt.abs(msc_ac)) cycle ms_loop

        ! time-reversal symmetry (spin-flip) used:
        ! process only cases with MS(A) >= 0 for OP1, OP1OP2
        if (use_tr_here.and.
c     &      (ms12i_a(idx_restr).lt.0.or.ms12i_a(3).lt.0) then
     &      (ms12i_a(3).lt.0 .or.
     &       ms12i_a(3).eq.0.and.ms12i_a(idx_restr).lt.0)) then
c dbg
c          print *,'skipping: ',ms12i_a(1:3)
c dbg
          cycle ms_loop
        end if
        fac_scal0 = 1d0
        if (use_tr_here.and.ms12i_a(3).gt.0) fac_scal0 = fac_scal0*2d0
        if (use_tr_here.and.ms12i_a(3).eq.0.and.
     &              ms12i_a(idx_restr).gt.0) fac_scal0 = fac_scal0*2d0

        if (buftyp1.eq.1) then
          idxms = msa2idxms4op(ms12i_a(1),mstop1,na_op1,nc_op1)
          ioff_op1 = gam_ms_op1(1,idxms)
          lenblock = len_gam_ms_op1(1,idxms)
          do isym = 2, nsym
            lenblock = lenblock + len_gam_ms_op1(isym,idxms)
          end do
c dbg
c          print *,'loading MS blk for op1'
c dbg
          call get_vec(ffop1,xop1,idoffop1+ioff_op1+1,
     &                          idoffop1+ioff_op1+lenblock)
        end if
        if (buftyp2.eq.1) then
          idxms = msa2idxms4op(ms12i_a(2),mstop2,na_op2,nc_op2)
          ioff_op2 = gam_ms_op2(1,idxms)
          lenblock = len_gam_ms_op2(1,idxms)
          do isym = 2, nsym
            lenblock = lenblock + len_gam_ms_op2(isym,idxms)
          end do
c dbg
c          print *,'loading MS blk for op2'
c dbg
          call get_vec(ffop2,xop2,idoffop2+ioff_op2+1,
     &                          idoffop2+ioff_op2+lenblock)
        end if

        if (buftyp12.eq.1) then
          idxms = msa2idxms4op(ms12i_a(3),mstop1op2,na_op1op2,nc_op1op2)
          ioff_op1op2 = gam_ms_op1op2(1,idxms)
          lenblock = len_gam_ms_op1op2(1,idxms)
          do isym = 2, nsym
            lenblock = lenblock + len_gam_ms_op1op2(isym,idxms)
          end do
          if (update) then
c dbg
c          print *,'loading MS blk for op1op2 ',lenblock
c dbg
            call get_vec(ffop1op2,xop1op2,idoffop1op2+ioff_op1op2+1,
     &                                 idoffop1op2+ioff_op1op2+lenblock)
          else if (idxms_op1op2_last.ne.idxms) then
c dbg
c          print *,'zeroing MS blk for op1op2 ',lenblock
c dbg
            xop1op2(1:lenblock) = 0d0
          end if
          idxms_op1op2_last = idxms

        end if

        msc_loop: do msc_a = mscmx_a, -mscmx_a, -2
          msc_c = msc_ac - msc_a
          if (abs(msc_c).gt.mscmx_c) cycle

          ! Ms of string after lifting restrictions:
          !  Op1(C1,A1) -> Op1(C10,A10;CC,AC)
          msex1_c = ms12i_c(1) - msc_c
          if (abs(msex1_c).gt.nc_ex1)
     &         cycle msc_loop
          msex1_a = ms12i_a(1) - msc_a
          if (abs(msex1_a).gt.na_ex1)
     &         cycle msc_loop
          msex2_a = ms12i_a(2) - msc_c   ! other way
          msex2_c = ms12i_c(2) - msc_a   ! round (!!)
          if (abs(msex2_c).gt.nc_ex2)
     &         cycle msc_loop
          if (abs(msex2_a).gt.na_ex2)
     &         cycle msc_loop

          if (use_tr_here.and.
     &        ms12i_a(idx_restr).eq.0.and.ms12i_a(3).eq.0.and.
     &        msc_c-msc_a.gt.0) then
c dbg
c          print *,'skipping: ',ms12i_a(1:3),'msc_c,msc_a: ',msc_c,msc_a
c dbg
            cycle msc_loop
          end if
          fac_scal = fac_scal0
          if (use_tr_here.and.
     &         ms12i_a(idx_restr).eq.0.and.ms12i_a(3).eq.0.and.
     &         msc_c-msc_a.lt.0) fac_scal = fac_scal*2d0
c dbg
          print *,'processing: ',ms12i_a(1:3),
     &            'msc_c,msc_a: ',msc_c,msc_a,' fac = ',fac_scal
c dbg

          if (ntest.ge.100) then
            write(lulog,*) 'Current spin case:'
            write(lulog,*) ' OP1/OP2/INT (C) ->',ms12i_c(1:3)
            write(lulog,*) ' OP1/OP2/INT (A) ->',ms12i_a(1:3)
            write(lulog,*) ' CNT(C)/CNT(A)   ->',msc_c,msc_a
          end if

          ! loop IRREP cases of (Op1(A),Op2(A),Interm)
          first2 = .true.
          igam_op1op2_last = -1
          gam_loop: do
            if (first2) then
              first2 = .false.
              ! initial IRREP distribution
              igam12i_raw(1:3) = igambnd(1,1:3)
            else
              ! next IRREP distribution
              if (.not.next_dist(igam12i_raw,3,igambnd,+1)) exit
            end if

            if (.not.tra_op1_) then
              igam12i_a(1) = igam12i_raw(1)
              igam12i_c(1) = multd2h(igam12i_a(1),igamtop1)
            else
              igam12i_c(1) = igam12i_raw(1)
              igam12i_a(1) = multd2h(igam12i_c(1),igamtop1)
            end if
            if (.not.tra_op2_) then
              igam12i_a(2) = igam12i_raw(2)
              igam12i_c(2) = multd2h(igam12i_a(2),igamtop2)
            else
              igam12i_c(2) = igam12i_raw(2)
              igam12i_a(2) = multd2h(igam12i_c(2),igamtop2)
            end if
            if (.not.tra_op1op2) then
              igam12i_a(3) = igam12i_raw(3)
              igam12i_c(3) = multd2h(igam12i_a(3),igamtop1op2)
            else
              igam12i_c(3) = igam12i_raw(3)
              igam12i_a(3) = multd2h(igam12i_c(3),igamtop1op2)
            end if

            ! set up start addresses
            ! need to be modified, if more than one distribution
            ! exists, see below

            idxms = msa2idxms4op(ms12i_a(1),mstop1,
     &                           na_op1,nc_op1)
            if (len_gam_ms_op1(igam12i_raw(1),idxms).eq.0)
     &           cycle gam_loop
            idxms = msa2idxms4op(ms12i_a(2),mstop2,
     &                           na_op2,nc_op2)
            if (len_gam_ms_op2(igam12i_raw(2),idxms).eq.0)
     &           cycle gam_loop
            idxms = msa2idxms4op(ms12i_a(3),mstop1op2,
     &                           na_op1op2,nc_op1op2)
            if (len_gam_ms_op1op2(igam12i_raw(3),idxms).eq.0)
     &           cycle gam_loop

            idxms = msa2idxms4op(ms12i_a(1),mstop1,na_op1,nc_op1)
            if (buftyp1.eq.2) then
              ioff_op1 = gam_ms_op1(igam12i_raw(1),idxms)
              lenblock = len_gam_ms_op1(igam12i_raw(1),idxms)
c dbg
c          print *,'loading GAM blk for op1, ',igam12i_a(1),idxms
c dbg
              call get_vec(ffop1,xop1,idoffop1+ioff_op1+1,
     &             idoffop1+ioff_op1+lenblock)
            end if
            idxop1 = gam_ms_op1(igam12i_raw(1),idxms) + 1
     &             - ioff_op1

            idxms = msa2idxms4op(ms12i_a(2),mstop2,na_op2,nc_op2)
            if (buftyp2.eq.2) then
              ioff_op2 = gam_ms_op2(igam12i_raw(2),idxms)
              lenblock = len_gam_ms_op2(igam12i_raw(2),idxms)
c dbg
c          print *,'loading GAM blk for op2, ',igam12i_a(2),idxms
c dbg
              call get_vec(ffop2,xop2,idoffop2+ioff_op2+1,
     &             idoffop2+ioff_op2+lenblock)
            end if
            idxop2 = gam_ms_op2(igam12i_raw(2),idxms) + 1
     &             - ioff_op2
            idxms =
     &           msa2idxms4op(ms12i_a(3),mstop1op2,na_op1op2,nc_op1op2)

            ! inefficient, but it works ....
            if (buftyp12.eq.2) then
              ioff_op1op2 = gam_ms_op1op2(igam12i_raw(3),idxms)
              lenblock = len_gam_ms_op1op2(igam12i_raw(3),idxms)
              if (update) then
                call get_vec(ffop1op2,xop1op2,idoffop1op2+ioff_op1op2+1,
     &             idoffop1op2+ioff_op1op2+lenblock)
              else if (igam_op1op2_last.ne.igam12i_raw(3).and.
     &                 idxms_op1op2_last.ne.idxms) then
                xop1op2(1:lenblock) = 0d0
              else
                call get_vec(ffop1op2,xop1op2,idoffop1op2+ioff_op1op2+1,
     &             idoffop1op2+ioff_op1op2+lenblock)
              end if
              igam_op1op2_last = igam12i_raw(3)

            end if


c            idxms = (na_op1op2-ms12i_a(3))/2 + 1
            ! relevant for case where no reordering necessary
            ! then we have: op1op2tmp == op1op2
            if (iblkop1op2.gt.0)
     &           idxop1op2 = gam_ms_op1op2(igam12i_raw(3),idxms) + 1
c     &                - idxst_op1op2+1
     &                - ioff_op1op2
            if (reo_op1op2)
     &           lblk_op1op2tmp=len_gam_ms_op1op2tmp(igam12i_raw(3),
     &                                               idxms)
            if (iblkop1op2.eq.0) idxop1op2 = 1

            igamc_ac = multd2h(igam12i_a(1),igam12i_a(2))
            igamc_ac = multd2h(igamc_ac,igam12i_a(3))

            ! LOWER incore requirements:
            !  load here Op1, Op2, Op1Op2 block

            gamc_loop: do igamc_a = 1, nsym
              igamc_c = multd2h(igamc_a,igamc_ac)

              ! info on zero length CNT string?
              ! this we know for sure:
              if (ncblk_cnt.eq.0.and.igamc_c.gt.1) cycle gamc_loop
              if (nablk_cnt.eq.0.and.igamc_a.gt.1) cycle gamc_loop

              ! IRREPs after lifting restrictions (cf. above)
              igamex1_a = multd2h(igam12i_a(1),igamc_a)
              igamex1_c = multd2h(igam12i_c(1),igamc_c)
              igamex2_a = multd2h(igam12i_a(2),igamc_c) !  !!
              igamex2_c = multd2h(igam12i_c(2),igamc_a) !  !!

              if (ncblk_ex1.eq.0.and.igamex1_c.gt.1) cycle gamc_loop
              if (nablk_ex1.eq.0.and.igamex1_a.gt.1) cycle gamc_loop
              if (ncblk_ex2.eq.0.and.igamex2_c.gt.1) cycle gamc_loop
              if (nablk_ex2.eq.0.and.igamex2_a.gt.1) cycle gamc_loop

c              call atim_cs(cpu00,sys00)

              ! loop over distributions of current Ms and IRREP
              ! of Aex1 and Cex1 over ngastypes
              first3 = .true.
              caex1_loop: do
                if (.not.next_msgamdist2(first3,
     &             msex1dis_c,msex1dis_a,gmex1dis_c,gmex1dis_a,
     &             ncblk_ex1, nablk_ex1,
     &             cinfo_ex1c,cinfo_ex1a,
     &             msex1_c,msex1_a,igamex1_c,igamex1_a,nsym,
     &             ms_fix1,fix_success1))
     &             exit caex1_loop
                first3 = .false.
                if(ms_fix1.and..not.fix_success1)then
c dbg
c                  print *,'cycle caex1_loop'
c dbg
c                  cycle caex1_loop
                endif

c                call ms2idxms(idxmsex1dis_c,msex1dis_c,
c     &               cinfo_ex1c,ncblk_ex1)
                 do ii = 1, ncblk_ex1
                   idxmsex1dis_c(ii)
     &             =ishft(cinfo_ex1c(ii,1)-msex1dis_c(ii),-1)+1
                 end do
c                call ms2idxms(idxmsex1dis_a,msex1dis_a,
c     &               cinfo_ex1a,nablk_ex1)
                 do ii = 1, nablk_ex1
                   idxmsex1dis_a(ii)
     &             =ishft(cinfo_ex1a(ii,1)-msex1dis_a(ii),-1)+1
                 end do

                call set_len_str2(lstrex1,ncblk_ex1,nablk_ex1,
     &                  lenstr_array,nsym,maxidxms,
     &                  cinfo_ex1c(1,2),idxmsex1dis_c,
     &                                 gmex1dis_c,cinfo_ex1c(1,3),
     &                  cinfo_ex1a(1,2),idxmsex1dis_a,
     &                                 gmex1dis_a,cinfo_ex1a(1,3),
     &                  hpvxseq,.false.)

                ! test C and A separately to avoid overflow
                if ( ncblk_ex1+nablk_ex1.gt.0 .and.
     &               idxlist(0,lstrex1,
     &                          ncblk_ex1+nablk_ex1,1).gt.0)
     &               cycle

                ! loop over distributions of current Ms and IRREP
                ! of Aex2 and Cex2 over ngastypes
                first4 = .true.
                caex2_loop: do
                  if (.not.next_msgamdist2(first4,
     &               msex2dis_c,msex2dis_a,gmex2dis_c,gmex2dis_a,
     &               ncblk_ex2, nablk_ex2,
     &               cinfo_ex2c,cinfo_ex2a,
     &               msex2_c,msex2_a,igamex2_c,igamex2_a,nsym,
     &               ms_fix2,fix_success2))
     &               exit caex2_loop
                  first4 = .false.
                  if(ms_fix2.and..not.fix_success2)then
c dbg
c                    print *,'cycle caex2_loop'
c dbg
c                    cycle caex2_loop
                  endif

c                  call ms2idxms(idxmsex2dis_c,msex2dis_c,
c     &                 cinfo_ex2c,ncblk_ex2)
                  do ii = 1, ncblk_ex2
                     idxmsex2dis_c(ii)
     &               =ishft(cinfo_ex2c(ii,1)-msex2dis_c(ii),-1)+1
                  end do
c                  call ms2idxms(idxmsex2dis_a,msex2dis_a,
c     &                 cinfo_ex2a,nablk_ex2)
                  do ii = 1, nablk_ex2
                     idxmsex2dis_a(ii)
     &               =ishft(cinfo_ex2a(ii,1)-msex2dis_a(ii),-1)+1
                  end do

                  call set_len_str2(lstrex2,ncblk_ex2,nablk_ex2,
     &                 lenstr_array,nsym,maxidxms,
     &                 cinfo_ex2c(1,2),idxmsex2dis_c,
     &                                gmex2dis_c,cinfo_ex2c(1,3),
     &                 cinfo_ex2a(1,2),idxmsex2dis_a,
     &                                gmex2dis_a,cinfo_ex2a(1,3),
     &                 hpvxseq,.false.)

                  if ( ncblk_ex2+nablk_ex2.gt.0.and.
     &                 idxlist(0,lstrex2,
     &                 ncblk_ex2+nablk_ex2,1).gt.0)
     &                 cycle

                  ! get Ms and IRREP distribution of intermediate
                  call merge_msgmdis(msi_dis_c,gmi_dis_c,
     &                                 ncblk_op1op2tmp,
     &                                 msex1dis_c,gmex1dis_c,
     &                                 msex2dis_c,gmex2dis_c,
     &                                 map_info_12c)
                  call merge_msgmdis(msi_dis_a,gmi_dis_a,
     &                                 nablk_op1op2tmp,
     &                                 msex2dis_a,gmex2dis_a,
     &                                 msex1dis_a,gmex1dis_a,
     &                                 map_info_12a)

c                  call ms2idxms(idxmsi_dis_c,msi_dis_c,
c     &                   cinfo_op1op2tmpc,ncblk_op1op2tmp)
                  do ii = 1, ncblk_op1op2tmp
                     idxmsi_dis_c(ii)
     &               =ishft(cinfo_op1op2tmpc(ii,1)-msi_dis_c(ii),-1)+1
                  end do
c                  call ms2idxms(idxmsi_dis_a,msi_dis_a,
c     &                   cinfo_op1op2tmpa,nablk_op1op2tmp)
                  do ii = 1, nablk_op1op2tmp
                     idxmsi_dis_a(ii)
     &               =ishft(cinfo_op1op2tmpa(ii,1)-msi_dis_a(ii),-1)+1
                  end do

                  call set_len_str2(
     &                   lstrop1op2tmp,ncblk_op1op2tmp,nablk_op1op2tmp,
     &                   lenstr_array,nsym,maxidxms,
     &                   cinfo_op1op2tmpc(1,2),idxmsi_dis_c,
     &                                  gmi_dis_c,cinfo_op1op2tmpc(1,3),
     &                   cinfo_op1op2tmpa(1,2),idxmsi_dis_a,
     &                                  gmi_dis_a,cinfo_op1op2tmpa(1,3),
     &                   hpvxseq,.false.)

                  if ( ncblk_op1op2tmp+nablk_op1op2tmp.gt.0 .and.
     &                 idxlist(0,lstrop1op2tmp,
     &                          ncblk_op1op2tmp+nablk_op1op2tmp,1).gt.0)
     &                 cycle caex2_loop

                  ! get igrphext1,igrphext2->igrphop1op2 map
                  ! for given ms and irreps
c                  ifree = mem_setmark('ex_str')
                  ! special product with map_info ...
                  lenmap = get_lenmap(lstrex1,lstrex2,
     &                   map_info_12c,ncblk_op1op2tmp)
c dbg
c                  print *,'lenmap C: ',lenmap
c dbg

c                  ifree = mem_alloc_int(map_ex1ex2c,lenmap,'strmap_c')
                  ! for C: ex1,ex2 sequence !
                  call get_strmap_blk_c(map_ex1ex2c,
     &                 ncblk_ex1,ncblk_ex2,ncblk_op1op2tmp,
     &                 cinfo_ex1c,cinfo_ex2c,lstrex1,lstrex2,
     &                 cinfo_ex1c(1,2),cinfo_ex2c(1,2),
     &                                 cinfo_op1op2tmpc(1,2),
     &                 idxmsex1dis_c,idxmsex2dis_c,
     &                 gmex1dis_c,gmex2dis_c,map_info_12c,
     &                 strmap_info,nsym,str_info%ngraph)

                  lenmap = get_lenmap(lstrex2(ncblk_ex2+1),
     &                                lstrex1(ncblk_ex1+1),
     &                   map_info_12a,nablk_op1op2tmp)
c dbg
c                  print *,'lenmap A: ',lenmap
c dbg
c                  ifree = mem_alloc_int(map_ex1ex2a,lenmap,'strmap_a')
                  ! for A: ex2,ex1 sequence !
                  call get_strmap_blk_c(map_ex1ex2a,
     &                 nablk_ex2,nablk_ex1,nablk_op1op2tmp,
     &                 cinfo_ex2a,cinfo_ex1a,
     &                  lstrex2(ncblk_ex2+1),lstrex1(ncblk_ex1+1),
     &                 cinfo_ex2a(1,2),cinfo_ex1a(1,2),
     &                                 cinfo_op1op2tmpa(1,2),
     &                 idxmsex2dis_a,idxmsex1dis_a,
     &                 gmex2dis_a,gmex1dis_a,map_info_12a,
     &                 strmap_info,nsym,str_info%ngraph)

                  ! get distribution index
                  idxms = msa2idxms4op(ms12i_a(3),mstop1op2,
     &                                 na_op1op2,nc_op1op2)
c                  idxms = (na_op1op2-ms12i_a(3))/2 + 1
                  idxdis_op1op2 = 1
c dbg
c                  print *,'igam12i_a(3),idxms: ',igam12i_a(3),idxms
c                  print *,'ndis_op1op2tmp(igam12i_a(3),idxms):',
c     &                 ndis_op1op2tmp(igam12i_a(3),idxms)
c dbg

                  if (me_op1op2tmp%diag_type.ne.0) then
                    ! skip non-diagonal distributions ...
                    if (nondia_distr(mapca12,diag_idx12,diag_ca12,
     &                     msi_dis_c,msi_dis_a,gmi_dis_c,gmi_dis_a,
     &                     ncblk_op1op2tmp,me_op1op2tmp%msdiag,
     &                     me_op1op2tmp%gamdiag)) cycle caex2_loop
                  end if

                  if (iblkop1op2tmp.gt.0.and.
     &                 ndis_op1op2tmp(igam12i_raw(3),idxms).gt.1) then

c dbg
c                    print *,'here 1'
c dbg
c dbg
                    if(ms_fix12)then
c                      print *,'msop1dis_c',msop1dis_c
c                      print *,'msop1dis_a',msop1dis_a
                      do idx = 1, min(ncblk_op1op2tmp,nablk_op1op2tmp)
                        if(idxmsi_dis_c(idx).ne.idxmsi_dis_a(idx))
     &                       cycle caex2_loop
                      enddo
                    endif
c dbg

                    idxdis =
     &                  idx_msgmdst2(.true.,
     &                   iblkop1op2tmp,idxms,igam12i_raw(3),
     &                   cinfo_op1op2tmpc,idxmsi_dis_c,
     &                              gmi_dis_c,ncblk_op1op2tmp,
     &                   cinfo_op1op2tmpa,idxmsi_dis_a,
     &                              gmi_dis_a,nablk_op1op2tmp,
     &                   tra_op1op2,dmap_op1op2tmpc,dmap_op1op2tmpa,
     &                              me_op1op2tmp,nsym)

                    idxdis_op1op2 = idxdis

                    ! relevant for case w/o reordering
                    ! then we have op1op2tmp == op1op2
                    idxop1op2 =
     &                   d_gam_ms_op1op2(idxdis,igam12i_raw(3),idxms)+1
c     &                   - idxst_op1op2+1
     &                   - ioff_op1op2

                    if (reo_op1op2)
     &                 lblk_op1op2tmp =
     &                 len_d_gam_ms_op1op2tmp(idxdis,igam12i_raw(3),
     &                                        idxms)

                  end if

                  if (.not.reo_op1op2) then
                    ! direct update of result block
c                    xop1op2blk => xop1op2(idxop1op2:)
                    ! somehow, "=> xxx(i:)" may lead to seg. faults(?)
                    xop1op2blk => xop1op2(idxop1op2:lenop1op2)
                  else
                    ! put result to intermediate buffer
                    xop1op2blk => xbf12tmp
                    xbf12tmp(1:lblk_op1op2tmp) = 0d0
                  end if

                  ! loop over distributions of current Ms and IRREP
                  ! of AC and CC over ngastypes
                  first5 = .true.
                  cac_loop: do
                    if (.not.next_msgamdist2(first5,
     &                 msc_dis_c,msc_dis_a,gmc_dis_c,gmc_dis_a,
     &                 ncblk_cnt, nablk_cnt,
     &                 cinfo_cntc,cinfo_cnta,
     &                 msc_c,msc_a,igamc_c,igamc_a,nsym,
     &                 ms_fix12,fix_success12))
     &                 exit cac_loop
                    first5 = .false.
c                    if(ms_fix12.and..not.fix_success12)then
c dbg
c                      print *,'cycle cac_loop 1'
c dbg
c                      cycle cac_loop
c                    endif

                    ! length of contraction
c                    call ms2idxms(idxmsc_dis_c,msc_dis_c,
c     &                   cinfo_cntc,ncblk_cnt)
                    do ii = 1, ncblk_cnt
                      idxmsc_dis_c(ii)
     &                =ishft(cinfo_cntc(ii,1)-msc_dis_c(ii),-1)+1
                    end do
c                    call ms2idxms(idxmsc_dis_a,msc_dis_a,
c     &                   cinfo_cnta,nablk_cnt)
                    do ii = 1, nablk_cnt
                      idxmsc_dis_a(ii)
     &                =ishft(cinfo_cnta(ii,1)-msc_dis_a(ii),-1)+1
                    end do

                    ! length of contraction
                    call set_len_str2(lstrcnt,ncblk_cnt,nablk_cnt,
     &                  lenstr_array,nsym,maxidxms,
     &                  cinfo_cntc(1,2),idxmsc_dis_c,
     &                                  gmc_dis_c,cinfo_cntc(1,3),
     &                  cinfo_cnta(1,2),idxmsc_dis_a,
     &                                  gmc_dis_a,cinfo_cnta(1,3),
     &                  hpvxseq,.false.)

                    if ( ncblk_cnt+nablk_cnt.gt.0 .and.
     &                   idxlist(0,lstrcnt,
     &                   ncblk_cnt+nablk_cnt,1).gt.0)
     &                   cycle cac_loop

                    ! get Ms and IRREP distribution of op1
                    call merge_msgmdis(msop1dis_c,gmop1dis_c,
     &                                 ncblk_op1,
     &                                 msc_dis_c,gmc_dis_c,
     &                                 msex1dis_c,gmex1dis_c,
     &                                 map_info_1c)
                    call merge_msgmdis(msop1dis_a,gmop1dis_a,
     &                                 nablk_op1,
     &                                 msc_dis_a,gmc_dis_a,
     &                                 msex1dis_a,gmex1dis_a,
     &                                 map_info_1a)

c                    call ms2idxms(idxmsop1dis_c,msop1dis_c,
c     &                   cinfo_op1c,ncblk_op1)
                    do ii = 1, ncblk_op1
                      idxmsop1dis_c(ii)
     &                =ishft(cinfo_op1c(ii,1)-msop1dis_c(ii),-1)+1
                    end do
c                    call ms2idxms(idxmsop1dis_a,msop1dis_a,
c     &                   cinfo_op1a,nablk_op1)
                    do ii = 1, nablk_op1
                      idxmsop1dis_a(ii)
     &                =ishft(cinfo_op1a(ii,1)-msop1dis_a(ii),-1)+1
                    end do

                    call set_len_str2(
     &                   lstrop1,ncblk_op1,nablk_op1,
     &                   lenstr_array,nsym,maxidxms,
     &                   cinfo_op1c(1,2),idxmsop1dis_c,
     &                                    gmop1dis_c,cinfo_op1c(1,3),
     &                   cinfo_op1a(1,2),idxmsop1dis_a,
     &                                    gmop1dis_a,cinfo_op1a(1,3),
     &                   hpvxseq,.false.)

                    if ( ncblk_op1+nablk_op1.gt.0 .and.
     &                   idxlist(0,lstrop1,
     &                             ncblk_op1+nablk_op1,1).gt.0)
     &                   cycle cac_loop

c dbg
                    if(ms_fix1)then
                      iblkoff = (iblkop1-1)*op1%njoined
                      call condense_occ(dum1_c,dum1_a,
     &                     hpvx1_c,hpvx1_a,
     &                     op1%ihpvca_occ(1,1,iblkoff+1),
     &                     op1%njoined,hpvxblkseq)

                      call expand_occ(msd1,
     &                     me_op1%idx_graph(1,1,iblkoff+1),
     &                     ncblk_op1,nablk_op1,
     &                     msop1dis_c,msop1dis_a,
     &                     hpvx1_c,hpvx1_a,
     &                     op1%njoined)

                      do idx = 1, op1%njoined
                        tot_c = 0
                        tot_a = 0
                        do jdx = 1, ngastp
                          tot_c = tot_c + msd1(jdx,1,idx)
                          tot_a = tot_a + msd1(jdx,2,idx)
                        enddo
                        if(tot_c.ne.tot_a)then
c                          print *,'cycle cac_loop 2'
                          cycle cac_loop
                        endif
                      enddo

c                      do idx = 1, min(ncblk_op1,nablk_op1)
c                        if(idxmsop1dis_c(idx).ne.idxmsop1dis_a(idx))then
cc dbg
cc                          print *,'idxmsop1dis_c',idxmsop1dis_c
cc                          print *,'idxmsop1dis_a',idxmsop1dis_a
c                          print *,'cycle cac_loop 2'
cc dbg
c                          cycle cac_loop
c                        endif
c                      enddo
                    endif
c dbg

                    if (me_op1%diag_type.ne.0) then
                      ! skip non-diagonal distributions ...
                      if (nondia_distr(mapca1,diag_idx1,diag_ca1,
     &                      msop1dis_c,msop1dis_a,gmop1dis_c,gmop1dis_a,
     &                      ncblk_op1,me_op1%msdiag,
     &                      me_op1%gamdiag)) cycle cac_loop
                    end if

                    ! get distribution index
                    idxms =
     &                   msa2idxms4op(ms12i_a(1),mstop1,na_op1,nc_op1)
c                    idxms = (na_op1-ms12i_a(1))/2 + 1
                    if (ndis_op1(igam12i_raw(1),idxms).gt.1) then
c dbg
c                      print *,'here 2'
c dbg
                      idxdis =
     &                   idx_msgmdst2(.true.,
     &                     iblkop1,idxms,igam12i_raw(1),
     &                     cinfo_op1c,idxmsop1dis_c,
     &                              gmop1dis_c,ncblk_op1,
     &                     cinfo_op1a,idxmsop1dis_a,
     &                              gmop1dis_a,nablk_op1,
     &                     tra_op1_,dmap_op1c,dmap_op1a,me_op1,nsym)

                      idxop1 =
     &                     d_gam_ms_op1(idxdis,igam12i_raw(1),idxms) + 1
c     &                     - idxst_op1+1-ioff_op1
     &                     - ioff_op1
                    end if

                    xnrm = ddot(ielprd(lstrop1,
     &                   ncblk_op1+nablk_op1),
     &                   xop1(idxop1),1,xop1(idxop1),1)
                    if (xnrm.lt.1d-28) cycle cac_loop

                    ! get Ms and IRREP distribution of op2
                    ! remember: CNT^+ !
                    call merge_msgmdis(msop2dis_c,gmop2dis_c,
     &                                 ncblk_op2,
     &                                 msc_dis_a,gmc_dis_a,
     &                                 msex2dis_c,gmex2dis_c,
     &                                 map_info_2c)
                    call merge_msgmdis(msop2dis_a,gmop2dis_a,
     &                                 nablk_op2,
     &                                 msc_dis_c,gmc_dis_c,
     &                                 msex2dis_a,gmex2dis_a,
     &                                 map_info_2a)

c                    call ms2idxms(idxmsop2dis_c,msop2dis_c,
c     &                   cinfo_op2c,ncblk_op2)
                    do ii = 1, ncblk_op2
                      idxmsop2dis_c(ii)
     &                =ishft(cinfo_op2c(ii,1)-msop2dis_c(ii),-1)+1
                    end do
c                    call ms2idxms(idxmsop2dis_a,msop2dis_a,
c     &                   cinfo_op2a,nablk_op2)
                    do ii = 1, nablk_op2
                      idxmsop2dis_a(ii)
     &                =ishft(cinfo_op2a(ii,1)-msop2dis_a(ii),-1)+1
                    end do

                    call set_len_str2(
     &                   lstrop2,ncblk_op2,nablk_op2,
     &                   lenstr_array,nsym,maxidxms,
     &                   cinfo_op2c(1,2),idxmsop2dis_c,
     &                                    gmop2dis_c,cinfo_op2c(1,3),
     &                   cinfo_op2a(1,2),idxmsop2dis_a,
     &                                    gmop2dis_a,cinfo_op2a(1,3),
     &                   hpvxseq,.false.)

                    if ( ncblk_op2+nablk_op2.gt.0 .and.
     &                   idxlist(0,lstrop2,
     &                             ncblk_op2+nablk_op2,1).gt.0)
     &                   cycle cac_loop

                    if(ms_fix2)then
                      iblkoff = (iblkop2-1)*op2%njoined
                      call condense_occ(dum2_c,dum2_a,
     &                     hpvx2_c,hpvx2_a,
     &                     op2%ihpvca_occ(1,1,iblkoff+1),
     &                     op2%njoined,hpvxblkseq)

                      call expand_occ(msd2,
     &                     me_op2%idx_graph(1,1,iblkoff+1),
     &                     ncblk_op2,nablk_op2,
     &                     msop2dis_c,msop2dis_a,
     &                     hpvx2_c,hpvx2_a,
     &                     op2%njoined)

                      do idx = 1, op2%njoined
                        tot_c = 0
                        tot_a = 0
                        do jdx = 1, ngastp
                          tot_c = tot_c + msd2(jdx,1,idx)
                          tot_a = tot_a + msd2(jdx,2,idx)
                        enddo
                        if(tot_c.ne.tot_a)then
c                          print *,'cycle cac_loop 3'
                          cycle cac_loop
                        endif
                      enddo

c                      do idx = 1, min(ncblk_op2,nablk_op2)
c                        if(idxmsop2dis_c(idx).ne.idxmsop2dis_a(idx))then
cc dbg
c                          print *,'idxmsop2dis_c',idxmsop2dis_c
c                          print *,'idxmsop2dis_a',idxmsop2dis_a
c                          print *,'cycle cac_loop 3'
cc dbg
c                          cycle cac_loop
c                        endif
c                      enddo
                    endif

                    if (me_op2%diag_type.ne.0) then
                      ! skip non-diagonal distributions ...
                      if (nondia_distr(mapca2,diag_idx2,diag_ca2,
     &                      msop2dis_c,msop2dis_a,gmop2dis_c,gmop2dis_a,
     &                      ncblk_op2,me_op2%msdiag,
     &                      me_op2%gamdiag)) cycle cac_loop
                    end if

                    ! get distribution index
                    idxms =
     &                   msa2idxms4op(ms12i_a(2),mstop2,na_op2,nc_op2)
c                    idxms = (na_op2-ms12i_a(2))/2 + 1
                    if (ndis_op2(igam12i_raw(2),idxms).gt.1) then
c dbg
c                      print *,'here 3'
c dbg
                      idxdis =
     &                   idx_msgmdst2(.true.,
     &                     iblkop2,idxms,igam12i_raw(2),
     &                     cinfo_op2c,idxmsop2dis_c,
     &                              gmop2dis_c,ncblk_op2,
     &                     cinfo_op2a,idxmsop2dis_a,
     &                              gmop2dis_a,nablk_op2,
     &                     tra_op2_,dmap_op2c,dmap_op2a,me_op2,nsym)

                      idxop2 =
     &                     d_gam_ms_op2(idxdis,igam12i_raw(2),idxms)+1
c     &                     - idxst_op2+1
     &                     - ioff_op2
                    end if

                    xnrm = ddot(ielprd(lstrop2,
     &                   ncblk_op2+nablk_op2),
     &                   xop2(idxop2),1,xop2(idxop2),1)
                    if (xnrm.lt.1d-28) cycle cac_loop

                    ! if we get here, op1op2 will change:
                    nonzero = .true.

                    ! get igrphcnt,igrphext1->igrphop1 map
                    ! for given ms and irreps
c                    ifree = mem_setmark('cntstr')
                    lenmap = get_lenmap(lstrcnt,lstrex1,
     &                   map_info_1c,ncblk_op1)
c                    ifree = mem_alloc_int(map_ex1cntc,lenmap,'strmap_c')
                    call get_strmap_blk_c(map_ex1cntc,
     &                   ncblk_cnt,ncblk_ex1,ncblk_op1,
     &                   cinfo_cntc,cinfo_ex1c,lstrcnt,lstrex1,
     &                   cinfo_cntc(1,2),cinfo_ex1c(1,2),
     &                                   cinfo_op1c(1,2),
     &                   idxmsc_dis_c,idxmsex1dis_c,
     &                   gmc_dis_c,gmex1dis_c,map_info_1c,
     &                   strmap_info,nsym,str_info%ngraph)

                    lenmap = get_lenmap(lstrcnt(ncblk_cnt+1),
     &                                  lstrex1(ncblk_ex1+1),
     &                                  map_info_1a,nablk_op1)
                    call get_strmap_blk_c(map_ex1cnta,
     &                   nablk_op1op2,nablk_ex1,nablk_op1,
     &                   cinfo_cnta,cinfo_ex1a,
     &                     lstrcnt(ncblk_cnt+1),
     &                             lstrex1(ncblk_ex1+1),
     &                   cinfo_cnta(1,2),cinfo_ex1a(1,2),
     &                                   cinfo_op1a(1,2),
     &                   idxmsc_dis_a,idxmsex1dis_a,
     &                   gmc_dis_a,gmex1dis_a,map_info_1a,
     &                   strmap_info,nsym,str_info%ngraph)

                    ! get igrphcnt,igrphext2->igrphop2 map
                    ! for given ms and irreps
                    lenmap = get_lenmap(lstrcnt(ncblk_cnt+1),lstrex2,
     &                   map_info_2c,ncblk_op2)
c                    ifree = mem_alloc_int(map_ex2cntc,lenmap,'strmap_c')
c dbg
c                  print *,'getting map_ex2cntc:'
c dbg
                    call get_strmap_blk_c(map_ex2cntc,
     &                   nablk_op1op2,ncblk_ex2,ncblk_op2,
     &                   cinfo_cnta,cinfo_ex2c,
     &                     lstrcnt(ncblk_cnt+1),lstrex2,
     &                   cinfo_cnta(1,2),cinfo_ex2c(1,2),
     &                                   cinfo_op2c(1,2),
     &                   idxmsc_dis_a,idxmsex2dis_c,
     &                   gmc_dis_a,gmex2dis_c,map_info_2c,
     &                   strmap_info,nsym,str_info%ngraph)

                    lenmap = get_lenmap(lstrcnt,lstrex2(ncblk_ex2+1),
     &                   map_info_2a,nablk_op2)
c                    ifree = mem_alloc_int(map_ex2cnta,lenmap,'strmap_a')
                    call get_strmap_blk_c(map_ex2cnta,
     &                   ncblk_cnt,nablk_ex2,nablk_op2,
     &                   cinfo_cntc,cinfo_ex2a,
     &                       lstrcnt,lstrex2(ncblk_ex2+1),
     &                   cinfo_cntc(1,2),cinfo_ex2a(1,2),
     &                                   cinfo_op2a(1,2),
     &                   idxmsc_dis_c,idxmsex2dis_a,
     &                   gmc_dis_c,gmex2dis_a,map_info_2a,
     &                   strmap_info,nsym,str_info%ngraph)

c                    call atim_cs(cpu0,sys0)

                    ! make the contraction for this block
                    if (ntest.ge.100)
     &                   write(lulog,*) 'calling blk1blk2',
     &                   lenop1,idxop1,
     &                   lenop2,idxop2,
     &                   lenop1op2,idxop1op2
c                    if (ntest.ge.1000) then
c                      write(lulog,*) ' the maps:'
c                      write(lulog,*) ' X1X2(A): '
c                      call prt_strmap_c(map_ex1ex2a,
c     &                     cinfo_ex2,cinfo_ex2,
c     &                     cinfo_ex2(1,3),cinfo_ex2(1,3),
c     &                     lstrex2(1,2),lstrex1(1,2),
c     &                     nablk_ex2,nablk_ex1)
c                      write(lulog,*) ' X1X2(C): '
c                      call prt_strmap(map_ex1ex2c,
c     &                     iocc_ext1(1,1),iocc_ext2(1,1),
c     &                     lstrext1(1,1),lstrext2(1,1))
c                      write(lulog,*) ' X1C(A): '
c                      call prt_strmap(map_ex1cnta,
c     &                     iocc_cnt(1,2),iocc_ext1(1,2),
c     &                     lstrcnt(1,2),lstrext1(1,2))
c                      write(lulog,*) ' X1C(C): '
c                      call prt_strmap(map_ex1cntc,
c     &                     iocc_cnt(1,1),iocc_ext1(1,1),
c     &                     lstrcnt(1,1),lstrext1(1,1))
c                      write(lulog,*) ' X2C(A): '
c                      call prt_strmap(map_ex2cnta,
c     &                     iocc_cnt(1,1),iocc_ext2(1,2),
c     &                     lstrcnt(1,1),lstrext2(1,2))
c                      write(lulog,*) ' X2C(C): '
c                      call prt_strmap(map_ex2cntc,
c     &                     iocc_cnt(1,2),iocc_ext2(1,1),
c     &                     lstrcnt(1,2),lstrext2(1,1))
c                    end if
c dbg
c                    print *,'on call:'
c                    print *,'xop1: ',xop1(idxop1)
c                    print *,'xop2: ',xop2(idxop2)
c                    print *,'xop1op2: ',xop1op2(idxop1op2)
c dbg
c dbg
c                    if (lenop1op2.eq.20) then
c                      print *,' ::',idxop1op2,lenop1op2
c                      print *,gmi_dis_c
c                      print *,gmi_dis_a
c                      print *,lstrop1op2tmp
c                    end if
c                    call mem_check('before kernel')
c                    print *,'irt_contr',irt_contr
c                    print *,'lstrop1op2tmp: ',lstrop1op2tmp
c                    print *,'lstrex1: ',lstrex1
c                    print *,'lstrex2: ',lstrex2
c                    print *,'calling kernel'
c dbg
                    if (irt_contr.eq.2) then
c dbg
c                      if (lenop1op2.eq.1) then
c                        print *,'xop1op2blk before: ',xop1op2blk(1),
c     &                       xop1op2(1)
c                      end if
c dbg
                      call contr_blk1blk2_wmaps_c(
     &                     xfac*casign*fac_scal,
     &                   xop1op2blk,
     &                                 xop1(idxop1),xop2(idxop2),
     &                   tra_op1_, tra_op2_, tra_op1op2,
     &                   ncblk_op1,nablk_op1,ncblk_ex1,nablk_ex1,
     &                   ncblk_op2,nablk_op2,ncblk_ex2,nablk_ex2,
     &                   ncblk_cnt,nablk_cnt,
     &                           ncblk_op1op2tmp,nablk_op1op2tmp,
     &                   cinfo_op1c(1,3),cinfo_op1a(1,3),
     &                   cinfo_op2c(1,3),cinfo_op2a(1,3),
     &                   cinfo_op1op2tmpc(1,3),cinfo_op1op2tmpa(1,3),
     &                   lstrop1,lstrop2,lstrop1op2tmp,
     &                   lstrex1,lstrex2,lstrcnt,
     &                   map_info_12c, map_info_12a,
     &                   map_info_1c, map_info_1a,
     &                   map_info_2c, map_info_2a,
     &                   map_ex1ex2c, map_ex1ex2a,
     &                   map_ex1cntc, map_ex1cnta,
     &                   map_ex2cntc, map_ex2cnta
     &                   )
c dbg
c                      if (lenop1op2.eq.1) then
c                        print *,'xop1op2blk after: ',xop1op2blk(1),
c     &                       xop1op2(1)
c                      end if
c dbg
                    else
                      call contr_blk1blk2_blocked_mm(
     &                     xfac*casign*fac_scal,
     &                   xop1op2blk,
     &                                 xop1(idxop1),xop2(idxop2),
     &                   tra_op1_, tra_op2_, tra_op1op2,
     &                   xscr,lenscr,len_str_block,len_cnt_block,
     &                   ncblk_op1,nablk_op1,ncblk_ex1,nablk_ex1,
     &                   ncblk_op2,nablk_op2,ncblk_ex2,nablk_ex2,
     &                   ncblk_cnt,nablk_cnt,
     &                           ncblk_op1op2tmp,nablk_op1op2tmp,
     &                   cinfo_op1c(1,3),cinfo_op1a(1,3),
     &                   cinfo_op2c(1,3),cinfo_op2a(1,3),
     &                   cinfo_op1op2tmpc(1,3),cinfo_op1op2tmpa(1,3),
     &                   lstrop1,lstrop2,lstrop1op2tmp,
     &                   lstrex1,lstrex2,lstrcnt,
     &                   map_info_12c, map_info_12a,
     &                   map_info_1c, map_info_1a,
     &                   map_info_2c, map_info_2a,
     &                   map_ex1ex2c, map_ex1ex2a,
     &                   map_ex1cntc, map_ex1cnta,
     &                   map_ex2cntc, map_ex2cnta
     &                   )
                    end if
                    if (ntest.ge.100)
     &                   write(lulog,*) 'after blk1blk2'
c dbg
c                    call mem_check('after kernel')
c dbg

c                    call atim_cs(cpu,sys)
                    cnt_kernel(1) = cnt_kernel(1)+cpu-cpu0
                    cnt_kernel(2) = cnt_kernel(2)+sys-sys0

                  end do cac_loop

                  ! if necessary, reorder op1op2 block:
                  if (reo_op1op2.and.nonzero) then
c dbg
c          write(lulog,*) 'input block '
c          write(lulog,'(x,5g15.8)')    xbf12tmp(1:lblk_op1op2tmp)
c dbg
c          print *,'gamma c/a:',igam12i_c(3),igam12i_a(3)
c          print *,'msdis_c :',msi_dis_c
c          print *,'msdis_a :',msi_dis_a
c          print *,'gmdis_c :',gmi_dis_c
c          print *,'gmdis_a :',gmi_dis_a
c dbgend
cc          call wrt_mel_buf(lulog,5,xop1op2,me_op1op2,
cc     &         iblkop1op2,iblkop1op2,str_info,orb_info)
c          call wrt_mel_buf(lulog,5,xop1op2blk,me_op1op2tmp,
c     &         1,1,str_info,orb_info)
c dbg
c                    call atim_cs(cpu0,sys0)
                    cnt_used_reo = .true.
                    if (igam12i_a(3).ne.igam12i_raw(3))
     &                 call quit(1,'contr_op1op2_wmaps_c',
     &                           'check this case')
c dbg
c          print *,'lblk_op1op2tmp: ',lblk_op1op2tmp
c dbg
                    call reo_blk_wmaps_c(1d0,xop1op2,xop1op2blk,
     &                   lbuf_op1op2,lblk_op1op2tmp,
     &                   reo_info%sign_reo,
     &                   tra_op1op2, tra_op1op2,
     &                   ms12i_c(3),ms12i_a(3),
     &                                   igam12i_c(3),igam12i_a(3),
     &                   msi_dis_c,msi_dis_a,gmi_dis_c,gmi_dis_a,
     &                   ncblk_op1op2tmp,nablk_op1op2tmp,
     &                   cinfo_op1op2tmpc,cinfo_op1op2tmpa,
     &                   lstrop1op2tmp,
     &                   me_op1op2,iblkop1op2,
     &                   ncblk_op1op2,nablk_op1op2,
     &                   cinfo_op1op2c,cinfo_op1op2a,
     &                   reo_info%ncblk_reo,reo_info%nablk_reo,
     &                   reo_info%cinfo_reo_c,reo_info%cinfo_reo_a,
     &                   reo_info%ncblk_reo0,reo_info%nablk_reo0,
     &                   reo_info%cinfo_opreo0c,reo_info%cinfo_opreo0a,
     &                   reo_info%map_reo1c,reo_info%map_reo1a,
     &                   reo_info%map_reo2c,reo_info%map_reo2a,
     &                   nsym,str_info,strmap_info)
c dbg
c          write(lulog,*) 'reordered operator (',trim(op1op2%name),')'
c          call wrt_mel_buf(lulog,5,xop1op2,me_op1op2,
c     &         iblkop1op2,iblkop1op2,str_info,orb_info)
c dbg
c                    call atim_cs(cpu,sys)
                    cnt_reo(1) = cnt_reo(1)+cpu-cpu0
                    cnt_reo(2) = cnt_reo(2)+sys-sys0
                  end if

                end do caex2_loop
              end do caex1_loop

c              call atim_cs(cpu,sys)
              cnt_dloop(1) = cnt_dloop(1)+cpu-cpu00
              cnt_dloop(2) = cnt_dloop(2)+sys-sys00

            end do gamc_loop

            if (buftyp12.eq.2) then
              idxms =
     &           msa2idxms4op(ms12i_a(3),mstop1op2,na_op1op2,nc_op1op2)
              ioff_op1op2 = gam_ms_op1op2(igam12i_raw(3),idxms)
              lenblock = len_gam_ms_op1op2(igam12i_raw(3),idxms)
              if (type_xret.ne.0) then
                if (idxms_op1op2_last.eq.idxms.and.
     &              igam_op1op2_last.eq.igam12i_raw(3))
     &              xret(1) = xret(1) - xret_last ! no double counting
                xret_last = ddot(lenblock,xop1op2,1,xop1op2,1)
                xret(1) = xret(1) + xret_last
              end if
c dbg
c          print *,'punching GAM blk for op1op2, ',igam12i_a(3),idxms
c dbg
              call put_vec(ffop1op2,xop1op2,idoffop1op2+ioff_op1op2+1,
     &             idoffop1op2+ioff_op1op2+lenblock)

            end if

          end do gam_loop

        end do msc_loop

        if (buftyp12.eq.2)
     &       idxms_op1op2_last = idxms

        if (buftyp12.eq.1) then
          idxms = msa2idxms4op(ms12i_a(3),mstop1op2,na_op1op2,nc_op1op2)
          ioff_op1op2 = gam_ms_op1op2(1,idxms)
          lenblock = len_gam_ms_op1op2(1,idxms)
          do isym = 2, nsym
            lenblock = lenblock + len_gam_ms_op1op2(isym,idxms)
          end do
          if (type_xret.ne.0)
     &       xret(1) = xret(1) + ddot(lenblock,xop1op2,1,xop1op2,1)
c dbg
c          print *,'punching MS blk for op1op2 ',lenblock
c dbg
          if (.not.bufop1op2)
     &       call put_vec(ffop1op2,xop1op2,idoffop1op2+ioff_op1op2+1,
     &                                 idoffop1op2+ioff_op1op2+lenblock)
        end if

      end do ms_loop

      if (me_op1op2tmp%diag_type.ne.0)
     &         deallocate(mapca12,diag_idx12,diag_ca12)
      if (me_op1%diag_type.ne.0)
     &         deallocate(mapca1,diag_idx1,diag_ca1)
      if (me_op2%diag_type.ne.0)
     &         deallocate(mapca2,diag_idx2,diag_ca2)

c      if (use_tr_here.and..not.update) then
      if (use_tr_here) then
        if (buftyp12.ne.0)
     &       call quit(1,'contr_op1op2_wmaps_c','use_tr + buftyp!=0?')
        if (abs(me_op1op2%absym).ne.1) then
c          write(lulog,*) 'assuming AL-BE symmetry = +1'
c          fac_ab = +1
          call quit(1,'contr_op1op2_wmaps_c',
     &         'absym.ne.+/-1 for list: '//trim(me_op1op2%label))
        else
          fac_ab = dble(me_op1op2%absym)
        end if
        if (ntest.ge.1000) then
          if (iblkop1op2.gt.0
     &           ) then
           write(lulog,*) 'operator 12 bef. sym (',trim(op1op2%name),')'
            call wrt_mel_buf(lulog,5,xop1op2,me_op1op2,
     &             iblkop1op2,iblkop1op2,str_info,orb_info)
          end if
        end if
        call sym_ab_blk(xop1op2,
     &       xop1op2,0.5d0,fac_ab,
     &       me_op1op2,iblkop1op2,
     &       str_info,strmap_info,orb_info)
      end if

c      if (op1op2%name(1:3).eq.'_LT') then
      if (ntest.ge.1000) then
        if (iblkop1op2.gt.0
     &       ) then
          write(lulog,*) 'operator 12 on exit (',trim(op1op2%name),')'
          call wrt_mel_buf(lulog,5,xop1op2,me_op1op2,
     &         iblkop1op2,iblkop1op2,str_info,orb_info)
        end if
      end if

      if (type_xret.eq.2.and.buftyp12.eq.0) then
        xret(1) = xop1op2(1)
      else if (type_xret.eq.1.and.buftyp12.eq.0) then
        xret(1) = ddot(lenop1op2,xop1op2,1,xop1op2,1)
      end if

c      call atim_cs(cpu0,sys0)
      ! put result to disc
      if (.not.bufop1op2.and.buftyp12.eq.0) then
        call put_vec(ffop1op2,xop1op2,idoffop1op2+idxst_op1op2,
     &                    idoffop1op2+idxst_op1op2-1+lenop1op2)
      end if
c      call atim_cs(cpu,sys)
      cnt_wr(1) = cnt_wr(1)+cpu-cpu0
      cnt_wr(2) = cnt_wr(2)+sys-sys0

      if(ms_fix1) deallocate(dum1_c,dum1_a,hpvx1_c,hpvx1_a)
      if(ms_fix2) deallocate(dum2_c,dum2_a,hpvx2_c,hpvx2_a)

      deallocate(
     &     gmop1dis_c, gmop1dis_a,
     &     gmop2dis_c, gmop2dis_a,
     &     gmex1dis_c, gmex1dis_a,
     &     gmex2dis_c, gmex2dis_a,
     &     gmc_dis_c , gmc_dis_a ,
     &     gmi_dis_c , gmi_dis_a ,
     &     msop1dis_c, msop1dis_a,
     &     msop2dis_c, msop2dis_a,
     &     msex1dis_c, msex1dis_a,
     &     msex2dis_c, msex2dis_a,
     &     msc_dis_c , msc_dis_a ,
     &     msi_dis_c , msi_dis_a ,
     &     idxmsop1dis_c, idxmsop1dis_a,
     &     idxmsop2dis_c, idxmsop2dis_a,
     &     idxmsex1dis_c, idxmsex1dis_a,
     &     idxmsex2dis_c, idxmsex2dis_a,
     &     idxmsc_dis_c , idxmsc_dis_a ,
     &     idxmsi_dis_c , idxmsi_dis_a ,
     &     lstrex1,lstrex2,lstrcnt,
     &     lstrop1,lstrop2,lstrop1op2tmp
     &     )
      deallocate(dmap_op1c,dmap_op1a,
     &           dmap_op2c,dmap_op2a,
     &           dmap_op1op2tmpc,dmap_op1op2tmpa)
      deallocate(lenstr_array)

      ifree = mem_flushmark()

      if (ntest.ge.100) then
        if (type_xret.ne.0)
     &       write(lulog,*) 'xret on exit = ',xret(1)
      end if

      return
      end


