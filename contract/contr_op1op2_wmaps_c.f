*----------------------------------------------------------------------*
      subroutine contr_op1op2_wmaps_c(xfac,casign,
     &     update,xret,type_xret,
     &     me_op1,me_op2,me_op1op2,me_op1op2tmp,
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
     &     ntest = 000

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
     &     bufop1, bufop2, bufop1op2, tra_op1, tra_op2, tra_op1op2,
     &     first1, first2, first3, first4, first5,
     &     reo_op1op2, nonzero
      integer ::
     &     mstop1,mstop2,mstop1op2,
     &     igamtop1,igamtop2,igamtop1op2,
     &     nc_op1, na_op1, nc_op2, na_op2,
     &     nc_ex1, na_ex1, nc_ex2, na_ex2, 
     &     nc_op1op2, na_op1op2,
     &     nc_op1op2tmp, na_op1op2tmp,
     &     nc_cnt, na_cnt,
     &     nsym, ifree, lenscr, lenblock,
     &     idxst_op1, idxst_op2, idxst_op1op2,
     &     idxop1, idxop2, idxop1op2,
     &     lenop1, lenop2, lenop1op2,
     &     mscmx_a, mscmx_c, msc_ac, msc_a, msc_c,
     &     msex1_a, msex1_c, msex2_a, msex2_c,
     &     igamc_ac, igamc_a, igamc_c,
     &     igamex1_a, igamex1_c, igamex2_a, igamex2_c,
     &     idxms, idxdis, lenmap, lblk_op1op2tmp,
     &     idxdis_op1op2
      integer ::
     &     ncblk_op1, nablk_op1, ncblk_ex1, nablk_ex1, 
     &     ncblk_op2, nablk_op2, ncblk_ex2, nablk_ex2, 
     &     ncblk_op1op2, nablk_op1op2, ncblk_op1op2tmp, nablk_op1op2tmp, 
     &     ncblk_cnt, nablk_cnt,
     &     ncblk_op1op2_0, nablk_op1op2_0,
     &     ncblk_reo12,    nablk_reo12
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
     &     map_info_12a(:)

      real(8) ::
     &     xnrm
      real(8) ::
     &     cpu, sys, cpu0, sys0, cpu00, sys00

      real(8), pointer ::
     &     xop1(:), xop2(:), xop1op2(:), xscr(:)
      real(8), pointer ::
     &     xbf1(:), xbf2(:), xbf12(:), xbf12tmp(:), xop1op2blk(:)

      integer ::
     &     msbnd(2,3), igambnd(2,3),
     &     ms12i_a(3), ms12i_c(3), igam12i_a(3), igam12i_c(3)

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
     &     ndis_op2(:,:), d_gam_ms_op2(:,:,:), gam_ms_op2(:,:),
     &     ndis_op1op2tmp(:,:), d_gam_ms_op1op2(:,:,:),
     &     gam_ms_op1op2(:,:),
     &     len_gam_ms_op1op2tmp(:,:), len_d_gam_ms_op1op2tmp(:,:,:)

      integer, pointer ::
     &     cinfo_reo12c(:,:), cinfo_reo12a(:,:),
     &     cinfo_op1op2_0c(:,:), cinfo_op1op2_0a(:,:),
     &     map_info_reo1c(:), map_info_reo1a(:),
     &     map_info_reo2c(:), map_info_reo2a(:)

      type(graph), pointer ::
     &     graphs(:)

      integer, external ::
     &     ielsum, ielprd, idx_msgmdst2, get_lenmap, idxlist,
     &     max_dis_blk
      logical, external ::
     &     next_dist, next_msgamdist2, msa2idxms4op
      real(8), external ::
     &     ddot

      if (ntest.gt.0) then
        call write_title(luout,wst_dbg_subr,
     &       'contr_op1op2_wmaps_c at work')
      end if

      op1 => me_op1%op
      op2 => me_op2%op
      op1op2 => me_op1op2%op
      op1op2tmp => me_op1op2tmp%op

      ffop1 => me_op1%fhand
      ffop2 => me_op2%fhand
      ffop1op2 => me_op1op2%fhand

      ! any CA transposition necessary?
      tra_op1 = op1%dagger
      tra_op2 = op2%dagger
      tra_op1op2 = op1op2%dagger

      if (ntest.ge.10) then
        write(luout,*) 'list1:   ',trim(me_op1%label)
        write(luout,*) 'list2:   ',trim(me_op2%label)
        write(luout,*) 'list12:  ',trim(me_op1op2%label)
        write(luout,*) 'ffop1:   ',ffop1%name(1:len_trim(ffop1%name))
        write(luout,*) 'ffop2:   ',ffop2%name(1:len_trim(ffop2%name))
        write(luout,*) 'ffop1op2:',
     &       ffop1op2%name(1:len_trim(ffop1op2%name))
        write(luout,*) 'xfac = ',xfac
        write(luout,*) 'casign = ',casign
        if (type_xret.ne.0)
     &       write(luout,*) 'xret on entry = ',xret(1)
        write(luout,*) 'op1: ',trim(op1%name),
     &       ' block ',iblkop1
        write(luout,*) 'op2: ',trim(op2%name),
     &       ' block ',iblkop2
        if (iblkop1op2.gt.0) then
          write(luout,*) 'op1op2: ',trim(op1op2%name),
     &       ' block ',iblkop1op2
        else
          write(luout,*) 'op1op2: scalar'
        end if
      end if

      ! flag whether non-zero contribution to op1op2 occurred
      nonzero = .false.

      reo_op1op2 = reo_info%n_op_reo.gt.0
      if (reo_op1op2.and..not.associated(reo_info%map_reo1c))
     &     call quit(1,'contr_op1op2_wmaps_c',
     &     'reo_info is not consistent')
      if (ntest.ge.10) write(luout,*) 'reo_op1op2: ',reo_op1op2
      if (ntest.ge.10 .and. reo_op1op2) then
        write(luout,*) 'op1op2tmp: ',trim(op1op2tmp%name),
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
      igamtop1 = me_op1%gamt
      igamtop2 = me_op2%gamt
      igamtop1op2 = me_op1op2%gamt

      if (multd2h(igamtop1,igamtop2).ne.igamtop1op2)
     &     call quit(1,'contr_op1op2_wmaps_c','inconsistent symmetries')

      if (lenop1.le.0.or.lenop2.le.0.or.lenop1op2.le.0) then
        write(luout,*) 'lenop1, lenop2, lenop1op2: ',
     &                  lenop1, lenop2, lenop1op2
        call quit(1,'contr_op1op2_wmaps_c',
     &     'zero length for operator?')
      end if

      ifree = mem_setmark('contr1')

      if (ffop1%buffered.and.ffop1%incore(iblkop1).gt.0) then
        bufop1 = .true.
        xop1 => ffop1%buffer(idxst_op1:)
      else
        bufop1 = .false.
        ! LOWER incore requirements:
        ! check for length of operator
        ! if length < 2 * blocks * da_reclen, try to alloc
        !   full block
        ! else get largest symmetry block
        ifree = mem_alloc_real(xbf1,lenop1,'xbf1')
        xop1 => xbf1
        call get_vec(ffop1,xop1,idoffop1+idxst_op1,
     &                          idoffop1+idxst_op1-1+lenop1)
      end if
      if (ffop2%buffered.and.ffop2%incore(iblkop2).gt.0) then
        bufop2 = .true.
        xop2 => ffop2%buffer(idxst_op2:)
      else
        bufop2 = .false.
        ! LOWER incore requirements:
        ! see above
        ifree = mem_alloc_real(xbf2,lenop2,'xbf2')
        xop2 => xbf2
        call get_vec(ffop2,xop2,idoffop2+idxst_op2,
     &                          idoffop2+idxst_op2-1+lenop2)
      end if

      if (ntest.ge.100) write(luout,*) ' bufop1/2: ',bufop1,bufop2

      ! get result vector as well (as we update)
      ! refers to reordered op1op2
      if (iblkop1op2.gt.0) then
        if (ffop1op2%buffered.and.ffop1op2%incore(iblkop1op2).gt.0) then
          bufop1op2 = .true.
          xop1op2 => ffop1op2%buffer(idxst_op1op2:)
        else
          bufop1op2 = .false.
          ! LOWER incore requirements:
          ! see above
          ifree = mem_alloc_real(xbf12,lenop1op2,'xbf12')
          xop1op2 => xbf12
          if (update) then
            ! read from disc
            call get_vec(ffop1op2,xop1op2,idoffop1op2+idxst_op1op2,
     &                             idoffop1op2+idxst_op1op2-1+lenop1op2)
          else
            ! init with zero
            xop1op2(1:lenop1op2) = 0d0
          end if
        end if
        if (ntest.ge.100) write(luout,*) ' bufop1op2: ',bufop1op2
      else
        bufop1op2 = .true.
        xop1op2 => xret
        if (ntest.ge.100) write(luout,*) ' result is scalar '
      end if

      if (ntest.ge.1000) then
        ! this will work if all blocks incore, only:
        write(luout,*) 'operator 1 (',trim(op1%name),
     &                    ',list=',trim(me_op1%label),')'
        call wrt_mel_buf(luout,5,xop1,me_op1,iblkop1,iblkop1,
     &                  str_info,orb_info)
        write(luout,*) 'operator 2 (',trim(op2%name),
     &                    ',list=',trim(me_op2%label),')'
        call wrt_mel_buf(luout,5,xop2,me_op2,iblkop2,iblkop2,
     &                  str_info,orb_info)
        if (iblkop1op2.gt.0) then
          write(luout,*) 'operator 12 on entry (',trim(op1op2%name),
     &                                ',list=',trim(me_op1op2%label),')'

          call wrt_mel_buf(luout,5,xop1op2,me_op1op2,
     &                    iblkop1op2,iblkop1op2,
     &                    str_info,orb_info)
        end if
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
        ifree = mem_alloc_real(xscr,lenscr,'contr_scr')
      end if

      graphs => str_info%g

      ndis_op1 => me_op1%off_op_gmox(iblkop1)%ndis
      gam_ms_op1 => me_op1%off_op_gmo(iblkop1)%gam_ms
      d_gam_ms_op1 => me_op1%off_op_gmox(iblkop1)%d_gam_ms
      ndis_op2 => me_op2%off_op_gmox(iblkop2)%ndis
      gam_ms_op2 => me_op2%off_op_gmo(iblkop2)%gam_ms
      d_gam_ms_op2 => me_op2%off_op_gmox(iblkop2)%d_gam_ms
      ndis_op1op2tmp => me_op1op2tmp%off_op_gmox(iblkop1op2tmp)%ndis
      gam_ms_op1op2 => me_op1op2%off_op_gmo(iblkop1op2)%gam_ms
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
      if (na_op1op2.ne.na_op1op2tmp)
     &     call quit(1,'contr_op1op2_wmaps_c','unexpected 1a')
      if (nc_op1op2.ne.nc_op1op2tmp)
     &     call quit(1,'contr_op1op2_wmaps_c','unexpected 1b')
c dbg
      call sum_occ(nc_ex1,cinfo_ex1c,ncblk_ex1)
      call sum_occ(na_ex1,cinfo_ex1a,nablk_ex1)
      call sum_occ(nc_ex2,cinfo_ex2c,ncblk_ex2)
      call sum_occ(na_ex2,cinfo_ex2a,nablk_ex2)
      call sum_occ(nc_cnt,cinfo_cntc,ncblk_cnt)
      call sum_occ(na_cnt,cinfo_cnta,nablk_cnt)

      ! set up maps (if necessary)
      call strmap_man_c(lenmap,
     &     cinfo_ex1c(1,2),ncblk_ex1,
     &     cinfo_ex2c(1,2),ncblk_ex2,
     &     cinfo_op1op2tmpc(1,2),ncblk_op1op2tmp,map_info_12c,
     &     str_info,strmap_info,orb_info)
      ifree = mem_alloc_int(map_ex1ex2c,lenmap,'exmap_c')
      call strmap_man_c(lenmap,
     &     cinfo_ex2a(1,2),nablk_ex2,
     &     cinfo_ex1a(1,2),nablk_ex1,
     &     cinfo_op1op2tmpa(1,2),nablk_op1op2tmp,map_info_12a,
     &     str_info,strmap_info,orb_info)
      ifree = mem_alloc_int(map_ex1ex2a,lenmap,'exmap_a')
      call strmap_man_c(lenmap,
     &     cinfo_cntc(1,2),ncblk_cnt,
     &     cinfo_ex1c(1,2),ncblk_ex1,
     &     cinfo_op1c(1,2),ncblk_op1,map_info_1c,
     &     str_info,strmap_info,orb_info)
      ifree = mem_alloc_int(map_ex1cntc,lenmap,'op1map_c')
      call strmap_man_c(lenmap,
     &     cinfo_cnta(1,2),nablk_cnt,
     &     cinfo_ex1a(1,2),nablk_ex1,
     &     cinfo_op1a(1,2),nablk_op1,map_info_1a,
     &     str_info,strmap_info,orb_info)
      ifree = mem_alloc_int(map_ex1cnta,lenmap,'op1map_a')
      call strmap_man_c(lenmap,
     &     cinfo_cnta(1,2),nablk_cnt,
     &     cinfo_ex2c(1,2),ncblk_ex2,
     &     cinfo_op2c(1,2),ncblk_op2,map_info_2c,
     &     str_info,strmap_info,orb_info)
      ifree = mem_alloc_int(map_ex2cntc,lenmap,'op2map_c')
      call strmap_man_c(lenmap,
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
        call strmap_man_c(lenmap,
     &     cinfo_op1op2_0c(1,2),ncblk_op1op2_0,
     &     cinfo_reo12c(1,2),ncblk_reo12,
     &     cinfo_op1op2tmpc(1,2),ncblk_op1op2tmp,map_info_reo1c,
     &     str_info,strmap_info,orb_info)
        call strmap_man_c(lenmap,
     &     cinfo_reo12a(1,2),nablk_reo12,
     &     cinfo_op1op2_0a(1,2),nablk_op1op2_0,
     &     cinfo_op1op2tmpa(1,2),nablk_op1op2tmp,map_info_reo1a,
     &     str_info,strmap_info,orb_info)
        call strmap_man_c(lenmap,
     &     cinfo_op1op2_0c(1,2),ncblk_op1op2_0,
     &     cinfo_reo12c(1,2),ncblk_reo12,
     &     cinfo_op1op2c(1,2),ncblk_op1op2,map_info_reo2c,
     &     str_info,strmap_info,orb_info)
        call strmap_man_c(lenmap,
     &     cinfo_reo12a(1,2),nablk_reo12,
     &     cinfo_op1op2_0a(1,2),nablk_op1op2_0,
     &     cinfo_op1op2a(1,2),nablk_op1op2,map_info_reo2a,
     &     str_info,strmap_info,orb_info)
      end if

      ! minimum Ms(A) for ...
      msbnd(1,1) = -na_op1 ! operator 1
      msbnd(1,2) = -na_op2 ! operator 2        
      msbnd(1,3) = -na_op1op2 ! product
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
      nsym = orb_info%nsym
      igambnd(2,1) = nsym
      igambnd(2,2) = nsym
      igambnd(2,3) = nsym
      ! loop Ms-cases of (Op1(A),Op2(A),Op1Op2(A))
      first1 = .true.
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
          if (ntest.ge.100) then
            write(luout,*) 'Current spin case:'
            write(luout,*) ' OP1/OP2/INT (C) ->',ms12i_c(1:3)
            write(luout,*) ' OP1/OP2/INT (A) ->',ms12i_a(1:3)
            write(luout,*) ' CNT(C)/CNT(A)   ->',msc_c,msc_a
          end if

          ! loop IRREP cases of (Op1(A),Op2(A),Interm)
          first2 = .true.
          gam_loop: do
            if (first2) then
              first2 = .false.
              ! initial IRREP distribution
              igam12i_a(1:3) = igambnd(1,1:3)
            else
              ! next IRREP distribution
              if (.not.next_dist(igam12i_a,3,igambnd,+1)) exit
            end if
            ! set up start addresses
            ! need to be modified, if more than one distribution
            ! exists, see below
            idxms = msa2idxms4op(ms12i_a(1),mstop1,na_op1,nc_op1)
c            idxms = (na_op1-ms12i_a(1))/2 + 1
            idxop1 = gam_ms_op1(igam12i_a(1),idxms) + 1
     &             - idxst_op1+1
            idxms = msa2idxms4op(ms12i_a(2),mstop2,na_op2,nc_op2)
c            idxms = (na_op2-ms12i_a(2))/2 + 1
            idxop2 = gam_ms_op2(igam12i_a(2),idxms) + 1
     &             - idxst_op2+1
            idxms =
     &           msa2idxms4op(ms12i_a(3),mstop1op2,na_op1op2,nc_op1op2)
c            idxms = (na_op1op2-ms12i_a(3))/2 + 1
            ! relevant for case where no reordering necessary
            ! then we have: op1op2tmp == op1op2
            if (iblkop1op2.gt.0)
     &           idxop1op2 = gam_ms_op1op2(igam12i_a(3),idxms) + 1
     &                - idxst_op1op2+1
            if (reo_op1op2)
     &           lblk_op1op2tmp=len_gam_ms_op1op2tmp(igam12i_a(3),idxms)
            if (iblkop1op2.eq.0) idxop1op2 = 1

            igam12i_c(1) = multd2h(igam12i_a(1),igamtop1)
            igam12i_c(2) = multd2h(igam12i_a(2),igamtop2)
            igam12i_c(3) = multd2h(igam12i_a(3),igamtop1op2)

            igamc_ac = multd2h(igam12i_a(1),igam12i_a(2))
            igamc_ac = multd2h(igamc_ac,igam12i_a(3))

            ! LOWER incore requirements:
            !  load here Op1, Op2, Op1Op2 block

            gamc_loop: do igamc_a = 1, nsym
              igamc_c = multd2h(igamc_a,igamc_ac)

              ! IRREPs after lifting restrictions (cf. above)
              igamex1_a = multd2h(igam12i_a(1),igamc_a)
              igamex1_c = multd2h(igam12i_c(1),igamc_c)
              igamex2_a = multd2h(igam12i_a(2),igamc_c) !  !!
              igamex2_c = multd2h(igam12i_c(2),igamc_a) !  !!

              call atim_cs(cpu00,sys00)

              ! loop over distributions of current Ms and IRREP 
              ! of Aex1 and Cex1 over ngastypes
              first3 = .true.
              caex1_loop: do
                if (.not.next_msgamdist2(first3,
     &             msex1dis_c,msex1dis_a,gmex1dis_c,gmex1dis_a,
     &             ncblk_ex1, nablk_ex1,
     &             cinfo_ex1c,cinfo_ex1a,
     &             msex1_c,msex1_a,igamex1_c,igamex1_a,nsym))
     &             exit caex1_loop
                first3 = .false.

                call ms2idxms(idxmsex1dis_c,msex1dis_c,
     &               cinfo_ex1c,ncblk_ex1)
                call ms2idxms(idxmsex1dis_a,msex1dis_a,
     &               cinfo_ex1a,nablk_ex1)

                call set_len_str(lstrex1,ncblk_ex1,nablk_ex1,
     &                  graphs,
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
     &               msex2_c,msex2_a,igamex2_c,igamex2_a,nsym))
     &               exit caex2_loop
                  first4 = .false.

                  call ms2idxms(idxmsex2dis_c,msex2dis_c,
     &                 cinfo_ex2c,ncblk_ex2)
                  call ms2idxms(idxmsex2dis_a,msex2dis_a,
     &                 cinfo_ex2a,nablk_ex2)

                  call set_len_str(lstrex2,ncblk_ex2,nablk_ex2,
     &                 graphs,
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

                  call ms2idxms(idxmsi_dis_c,msi_dis_c,
     &                   cinfo_op1op2tmpc,ncblk_op1op2tmp)
                  call ms2idxms(idxmsi_dis_a,msi_dis_a,
     &                   cinfo_op1op2tmpa,nablk_op1op2tmp)

                  call set_len_str(
     &                   lstrop1op2tmp,ncblk_op1op2tmp,nablk_op1op2tmp,
     &                 graphs,
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
                  if (iblkop1op2tmp.gt.0.and.
     &                 ndis_op1op2tmp(igam12i_a(3),idxms).gt.1) then
                    idxdis =
     &                  idx_msgmdst2(
     &                   iblkop1op2tmp,idxms,igam12i_a(3),
     &                   cinfo_op1op2tmpc,idxmsi_dis_c,
     &                              gmi_dis_c,ncblk_op1op2tmp,
     &                   cinfo_op1op2tmpa,idxmsi_dis_a,
     &                              gmi_dis_a,nablk_op1op2tmp,
     &                   .false.,me_op1op2tmp,nsym)
                    idxdis_op1op2 = idxdis

                    ! relevant for case w/o reordering
                    ! then we have op1op2tmp == op1op2
                    idxop1op2 = 
     &                   d_gam_ms_op1op2(idxdis,igam12i_a(3),idxms)+1
     &                   - idxst_op1op2+1

                    if (reo_op1op2)
     &                 lblk_op1op2tmp =
     &                 len_d_gam_ms_op1op2tmp(idxdis,igam12i_a(3),idxms)

                  end if

                  if (.not.reo_op1op2) then
                    ! direct update of result block
                    xop1op2blk => xop1op2(idxop1op2:)
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
     &                 msc_c,msc_a,igamc_c,igamc_a,nsym))
     &                 exit cac_loop
                    first5 = .false.

                    ! length of contraction
                    call ms2idxms(idxmsc_dis_c,msc_dis_c,
     &                   cinfo_cntc,ncblk_cnt)
                    call ms2idxms(idxmsc_dis_a,msc_dis_a,
     &                   cinfo_cnta,nablk_cnt)

                    ! length of contraction
                    call set_len_str(lstrcnt,ncblk_cnt,nablk_cnt,
     &                  graphs,
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

                    call ms2idxms(idxmsop1dis_c,msop1dis_c,
     &                   cinfo_op1c,ncblk_op1)
                    call ms2idxms(idxmsop1dis_a,msop1dis_a,
     &                   cinfo_op1a,nablk_op1)
                    
                    call set_len_str(
     &                   lstrop1,ncblk_op1,nablk_op1,
     &                   graphs,
     &                   cinfo_op1c(1,2),idxmsop1dis_c,
     &                                    gmop1dis_c,cinfo_op1c(1,3),
     &                   cinfo_op1a(1,2),idxmsop1dis_a,
     &                                    gmop1dis_a,cinfo_op1a(1,3),
     &                   hpvxseq,.false.)

                    if ( ncblk_op1+nablk_op1.gt.0 .and.
     &                   idxlist(0,lstrop1,
     &                             ncblk_op1+nablk_op1,1).gt.0)
     &                   cycle cac_loop

                    ! get distribution index
                    idxms =
     &                   msa2idxms4op(ms12i_a(1),mstop1,na_op1,nc_op1)
c                    idxms = (na_op1-ms12i_a(1))/2 + 1
                    if (ndis_op1(igam12i_a(1),idxms).gt.1) then
                      idxdis =
     &                   idx_msgmdst2(
     &                     iblkop1,idxms,igam12i_a(1),
     &                     cinfo_op1c,idxmsop1dis_c,
     &                              gmop1dis_c,ncblk_op1,
     &                     cinfo_op1a,idxmsop1dis_a,
     &                              gmop1dis_a,nablk_op1,
     &                     .false.,me_op1,nsym)
                      idxop1 = 
     &                     d_gam_ms_op1(idxdis,igam12i_a(1),idxms) + 1
     &                     - idxst_op1+1
                    end if

c dbg
c                    print *,'call ddot(1)'
c dbg
                    xnrm = ddot(ielprd(lstrop1,
     &                   ncblk_op1+nablk_op1),
     &                   xop1(idxop1),1,xop1(idxop1),1)
                    if (xnrm.lt.1d-28) cycle cac_loop

                    ! get Ms and IRREP distribution of op2
                    ! remember: CNT^+ !
c dbg
c                    print *,'msc_dis_a:  ',msc_dis_a
c                    print *,'msex2dis_c: ',msex2dis_c
c dbg
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
c dbg
c                    print *,'msop2dis_c: ',msex2dis_c
c dbg

                    call ms2idxms(idxmsop2dis_c,msop2dis_c,
     &                   cinfo_op2c,ncblk_op2)
                    call ms2idxms(idxmsop2dis_a,msop2dis_a,
     &                   cinfo_op2a,nablk_op2)

                    call set_len_str(
     &                   lstrop2,ncblk_op2,nablk_op2,
     &                   graphs,
     &                   cinfo_op2c(1,2),idxmsop2dis_c,
     &                                    gmop2dis_c,cinfo_op2c(1,3),
     &                   cinfo_op2a(1,2),idxmsop2dis_a,
     &                                    gmop2dis_a,cinfo_op2a(1,3),
     &                   hpvxseq,.false.)

                    if ( ncblk_op2+nablk_op2.gt.0 .and.
     &                   idxlist(0,lstrop2,
     &                             ncblk_op2+nablk_op2,1).gt.0)
     &                   cycle cac_loop

                    ! get distribution index
                    idxms =
     &                   msa2idxms4op(ms12i_a(2),mstop2,na_op2,nc_op2)
c                    idxms = (na_op2-ms12i_a(2))/2 + 1
                    if (ndis_op2(igam12i_a(2),idxms).gt.1) then
                      idxdis =
     &                   idx_msgmdst2(
     &                     iblkop2,idxms,igam12i_a(2),
     &                     cinfo_op2c,idxmsop2dis_c,
     &                              gmop2dis_c,ncblk_op2,
     &                     cinfo_op2a,idxmsop2dis_a,
     &                              gmop2dis_a,nablk_op2,
     &                     .false.,me_op2,nsym)
                      idxop2 = 
     &                     d_gam_ms_op2(idxdis,igam12i_a(2),idxms) + 1
     &                     - idxst_op2+1
                    end if

c dbg
c                    print *,'call ddot(2)',idxop2,
c     &                   lstrop2(1:ncblk_op2+nablk_op2)
c dbg
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
     &                   idxmsc_dis_c,idxmsex2dis_a,
     &                   gmc_dis_c,gmex2dis_a,map_info_2a,
     &                   strmap_info,nsym,str_info%ngraph)                    

                    call atim_cs(cpu0,sys0)                    

                    ! make the contraction for this block
                    if (ntest.ge.100)
     &                   write(luout,*) 'calling blk1blk2',
     &                   lenop1,idxop1,
     &                   lenop2,idxop2,
     &                   lenop1op2,idxop1op2
c                    if (ntest.ge.1000) then
c                      write(luout,*) ' the maps:'
c                      write(luout,*) ' X1X2(A): '
c                      call prt_strmap_c(map_ex1ex2a,
c     &                     cinfo_ex2,cinfo_ex2,
c     &                     cinfo_ex2(1,3),cinfo_ex2(1,3),
c     &                     lstrex2(1,2),lstrex1(1,2),
c     &                     nablk_ex2,nablk_ex1)
c                      write(luout,*) ' X1X2(C): '
c                      call prt_strmap(map_ex1ex2c,
c     &                     iocc_ext1(1,1),iocc_ext2(1,1),
c     &                     lstrext1(1,1),lstrext2(1,1))
c                      write(luout,*) ' X1C(A): '
c                      call prt_strmap(map_ex1cnta,
c     &                     iocc_cnt(1,2),iocc_ext1(1,2),
c     &                     lstrcnt(1,2),lstrext1(1,2))
c                      write(luout,*) ' X1C(C): '
c                      call prt_strmap(map_ex1cntc,
c     &                     iocc_cnt(1,1),iocc_ext1(1,1),
c     &                     lstrcnt(1,1),lstrext1(1,1))
c                      write(luout,*) ' X2C(A): '
c                      call prt_strmap(map_ex2cnta,
c     &                     iocc_cnt(1,1),iocc_ext2(1,2),
c     &                     lstrcnt(1,1),lstrext2(1,2))
c                      write(luout,*) ' X2C(C): '
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
c dbg
                    if (irt_contr.eq.2) then
c dbg
c                      if (lenop1op2.eq.1) then
c                        print *,'xop1op2blk before: ',xop1op2blk(1),
c     &                       xop1op2(1)
c                      end if
c dbg
                      call contr_blk1blk2_wmaps_c(xfac*casign,
     &                   xop1op2blk,
     &                                 xop1(idxop1),xop2(idxop2),
     &                   tra_op1, tra_op2, tra_op1op2,
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
                      call contr_blk1blk2_blocked_mm(xfac*casign,
     &                   xop1op2blk,
     &                                 xop1(idxop1),xop2(idxop2),
     &                   tra_op1, tra_op2, tra_op1op2,
     &                   xscr,lenscr,len_str_block,
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
     &                   write(luout,*) 'after blk1blk2'
c dbg
c                    call mem_check('after kernel')
c dbg

                    call atim_cs(cpu,sys)
                    cnt_kernel(1) = cnt_kernel(1)+cpu-cpu0
                    cnt_kernel(2) = cnt_kernel(2)+sys-sys0

                  end do cac_loop
                  
                  ! if necessary, reorder op1op2 block:
                  if (reo_op1op2.and.nonzero) then
                    call reo_blk_wmaps_c(xop1op2,xop1op2blk,
     &                   reo_info%sign_reo,
     &                   tra_op1op2,
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
                  end if

                end do caex2_loop
              end do caex1_loop

              call atim_cs(cpu,sys)
              cnt_dloop(1) = cnt_dloop(1)+cpu-cpu00
              cnt_dloop(2) = cnt_dloop(2)+sys-sys00

            end do gamc_loop
          end do gam_loop

        end do msc_loop
      end do ms_loop

      if (ntest.ge.1000) then
        if (iblkop1op2.gt.0
     &       ) then
          write(luout,*) 'operator 12 on exit (',trim(op1op2%name),')'
          call wrt_mel_buf(luout,5,xop1op2,me_op1op2,
     &         iblkop1op2,iblkop1op2,str_info,orb_info)
        end if
      end if

c dbg
c      print *,'type_xret ',type_xret
c      print *,'xret' ,xret
c dbg
      if (type_xret.eq.2) then
        xret(1) = xop1op2(1)
      else if (type_xret.eq.1) then
        xret(1) = ddot(lenop1op2,xop1op2,1,xop1op2,1)
      end if

      ! put result to disc
      if (.not.bufop1op2) then
        call put_vec(ffop1op2,xop1op2,idoffop1op2+idxst_op1op2,
     &                    idoffop1op2+idxst_op1op2-1+lenop1op2)
      end if

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

      ifree = mem_flushmark()

      if (ntest.ge.100) then
        if (type_xret.ne.0)
     &       write(luout,*) 'xret on exit = ',xret(1)
      end if

      return
      end


