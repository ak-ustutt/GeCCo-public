      subroutine get_bc_info3(bc_sign,bc_sign_itf,possible,
     &     idx_op,iblk_op,
     &     iocc_ex1,iocc_ex2,iocc_cnt,
     &     iocc_op1,iocc_op2,iocc_op1op2,
     &     irestr_op1,irestr_op2,irestr_op1op2,
     &     tra_op1,tra_op2,tra_op1op2,
     &     mst_op,mst_op1op2,
     &     gamt_op,gamt_op1op2,
     &     njoined_op, njoined_op1op2, njoined_cnt,
     &     merge_op1, merge_op2, merge_op1op2, merge_op2op1,
     &     contr_in,occ_vtx_in,irestr_vtx_in,info_vtx,
     &     make_contr_red,
     &     contr_red,occ_vtx_red,irestr_vtx_red,info_vtx_red,
     &     set_reo, reo_info, reo_info_bef, !FIX
     &     iarc_contr,update_idxintm,idxintm,
     &     irestr_res,njoined_res,orb_info,op_info)

      implicit none

      include 'stdunit.h'
      include 'opdim.h'
      include 'ifc_operators.h'
      include 'def_contraction.h'
      include 'def_reorder_info.h'
      include 'mdef_operator_info.h'
      include 'def_orbinf.h'
      include 'multd2h.h'

      integer, parameter ::
     &     ntest = 00

      type(contraction), intent(in), target ::
     &     contr_in
      type(contraction), intent(out) ::
     &     contr_red
      type(reorder_info), intent(inout) ::
     &     reo_info, reo_info_bef
      type(orbinf), intent(in), target ::
     &     orb_info
      type(operator_info), intent(in) ::
     &     op_info
      logical, intent(in) ::
     &     set_reo, make_contr_red, update_idxintm
      logical, intent(out) ::
     &     possible
      integer, intent(inout) ::
     &     idxintm
      integer, intent(in) ::
     &     iarc_contr, njoined_res
      integer, intent(in) ::
     &     occ_vtx_in(ngastp,2,contr_in%nvtx+njoined_res),
     &     irestr_vtx_in(2,orb_info%ngas,2,2,contr_in%nvtx+njoined_res),
     &     info_vtx(2,contr_in%nvtx+njoined_res),
     &     irestr_res(2,orb_info%ngas,2,2,njoined_res)
      integer, intent(out) ::
     &     occ_vtx_red(ngastp,2,*),
     &     irestr_vtx_red(2,orb_info%ngas,2,2,*),
     &     info_vtx_red(2,*)
      real(8), intent(out) ::
     &     bc_sign, bc_sign_itf
      logical, intent(out) ::
     &     tra_op1, tra_op2, tra_op1op2
      integer, intent(out) ::
     &     idx_op(2), iblk_op(2),
     &     njoined_op(2), njoined_op1op2, njoined_cnt,
     &     iocc_ex1(ngastp,2,*), iocc_ex2(ngastp,2,*),
     &     iocc_cnt(ngastp,2,*),
     &     iocc_op1(ngastp,2,*), iocc_op2(ngastp,2,*),
     &     iocc_op1op2(ngastp,2,*),
     &     irestr_op1(2,orb_info%ngas,2,2,*),
     &     irestr_op2(2,orb_info%ngas,2,2,*),
     &     irestr_op1op2(2,orb_info%ngas,2,2,*),
     &     mst_op(2), mst_op1op2, gamt_op(2), gamt_op1op2,
     &     merge_op1(*), merge_op2(*), merge_op1op2(*), merge_op2op1(*)

      type(contraction), pointer ::
     &     contr, contr_pnt

      logical ::
     &     self, ok
      integer ::
     &     ld_mmap1, ld_mmap2, ld_mmap12, ngas,
     &     nvtx, ivtx, idx, iblk, idxnew_op1op2,
     &     ivtx1, ivtx2, isvtx1, isvtx2,
     &     len_list, nvtx_red, sh_sign, cnt_sign, itf_sign, ireo, jreo

      integer, pointer ::
     &     ireo_vtx_no(:), ireo_vtx_on(:),
     &     ireo_after_contr(:), ivtx_op1op2(:),
     &     arc_list(:),
     &     merge_map_op1(:,:,:), merge_map_op2(:,:,:),
     &     merge_map_op1op2(:,:,:),
     &     ihpvgas(:,:)
      integer ::
     &     occ_vtx(ngastp,2,contr_in%nvtx+njoined_res),
     &     irestr_vtx(2,orb_info%ngas,2,2,contr_in%nvtx+njoined_res)
      integer, allocatable ::
     &     svmap(:)


      integer, external ::
     &     imltlist
      logical, external ::
     &     allowed_contr

      if (ntest.ge.100) then
        call write_title(lulog,wst_dbg_subr,'get_bc_info3')
        call prt_contr3(lulog,contr_in,-1)
        write(lulog,*) 'info_vtx: '
        write(lulog,'(x,2i4)') info_vtx(1:2,1:contr_in%nvtx+njoined_res)
      end if

      ngas = orb_info%ngas
      ihpvgas => orb_info%ihpvgas

      nvtx = contr_in%nvtx

      ! test: check re-sort possibilities prior to contraction:
      allocate(contr)
      call init_contr(contr)
      call copy_contr(contr_in,contr)
      occ_vtx    = occ_vtx_in
      irestr_vtx = irestr_vtx_in

      if (.true..and.update_idxintm) then
c      if (.false..and.update_idxintm) then
        possible = .true.
        idxnew_op1op2 = idxintm
        call reorder_supvtx(possible,
     &       .true.,set_reo,.true.,reo_info_bef,
     &       contr,occ_vtx(1,1,njoined_res+1),
     &             irestr_vtx(1,1,1,1,njoined_res+1),idxnew_op1op2,
     &       orb_info)
        if (contr%nxarc.gt.0)
     &       call reorder_supvtx_x(possible,
     &         .true.,set_reo,.true.,reo_info_bef,
     &         contr,occ_vtx(1,1,njoined_res+1),
     &               irestr_vtx(1,1,1,1,njoined_res+1),idxnew_op1op2,
     &         orb_info)

        possible = .true.
        if (reo_info_bef%nreo.gt.0) idxintm = idxintm-1
        ! report that REO is before contraction!

c        ! to be commented out
c        if (reo_info_bef%nreo.gt.0) then
c          do ivtx = 1, nvtx
c            call dummy_restr(irestr_vtx(1,1,1,1,ivtx+njoined_res),
c     &           occ_vtx(1,1,ivtx+njoined_res),1,
c     &           orb_info)
c          end do
c        end if
        if (set_reo) call tidy_reo_info(reo_info_bef)

      end if

      ivtx1 = contr%arc(iarc_contr)%link(1)
      ivtx2 = contr%arc(iarc_contr)%link(2)

      isvtx1 = contr%svertex(ivtx1)
      njoined_op(1) = imltlist(isvtx1,contr%svertex,nvtx,1)
      if (ivtx1.le.nvtx) then
        idx_op(1) = contr%vertex(ivtx1)%idx_op
        ! iblkop fix:
        iblk_op(1) = (contr%vertex(ivtx1)%iblk_op-1)/njoined_op(1) + 1
        tra_op1   = contr%vertex(ivtx1)%dagger
      else
        write(lulog,*) 'ivtx1 = ',ivtx1
        write(lulog,*) 'nvtx = ',nvtx
        call prt_contr3(lulog,contr_in,-1)        
        call quit(1,'get_bc_info3','ivtx1>nvtx?')
      end if

      isvtx2 = contr%svertex(ivtx2)

      self = (isvtx1.eq.isvtx2)

      if (.not.self) then
        njoined_op(2) = imltlist(isvtx2,contr%svertex,nvtx,1)
        if (ivtx2.le.nvtx) then
          idx_op(2) = contr%vertex(ivtx2)%idx_op
          ! iblkop fix:
          iblk_op(2) = (contr%vertex(ivtx2)%iblk_op-1)/njoined_op(2) + 1
          tra_op2   = contr%vertex(ivtx2)%dagger
        else
          write(lulog,*) 'ivtx2 = ',ivtx2
          write(lulog,*) 'nvtx = ',nvtx
          call prt_contr3(lulog,contr_in,-1)        
          call quit(1,'get_bc_info3','ivtx2>nvtx?')
        end if
      else
        njoined_op(2) = 0
c        if (ivtx2.le.nvtx) then
          idx_op(2) = 0
          iblk_op(2) = 0
          tra_op2   = .false.
c        end if
      end if

c     &     call quit(1,'get_bc_info3','I am confused ....')

      mst_op(1) = info_vtx(1,ivtx1+njoined_res)
      gamt_op(1) = info_vtx(2,ivtx1+njoined_res)
      if (.not.self) then
        mst_op(2) = info_vtx(1,ivtx2+njoined_res)
        gamt_op(2) = info_vtx(2,ivtx2+njoined_res)
      else
        mst_op(2) = 0
        gamt_op(2) = 1
      end if

      ivtx1 = 0
      ivtx2 = 0
      do ivtx = 1, nvtx
        if (contr%svertex(ivtx).eq.isvtx1) then
          ivtx1 = ivtx1+1
          iocc_op1(1:ngastp,1:2,ivtx1) =
     &         occ_vtx(1:ngastp,1:2,ivtx+njoined_res)
          irestr_op1(1:2,1:ngas,1:2,1:2,ivtx1) =
     &         irestr_vtx(1:2,1:ngas,1:2,1:2,ivtx+njoined_res)
        end if
        if (contr%svertex(ivtx).eq.isvtx2.and..not.self) then
          ivtx2 = ivtx2+1
          iocc_op2(1:ngastp,1:2,ivtx2) =
     &         occ_vtx(1:ngastp,1:2,ivtx+njoined_res)
          irestr_op2(1:2,1:ngas,1:2,1:2,ivtx2) =
     &         irestr_vtx(1:2,1:ngas,1:2,1:2,ivtx+njoined_res)
        end if
      end do
      
      if (self) then
        iocc_op2(1:ngastp,1:2,1) = 0
        iocc_ex2(1:ngastp,1:2,1) = 0
        irestr_op2(1:2,1:ngas,1:2,1:2,1) = 0
      end if

      allocate(arc_list(contr%narc))
      ! set external and contraction indices
      ! get all involved arcs
      call get_associated_arcs(arc_list,len_list,iarc_contr,contr)

      ! set merging info for op1 and op2, i.e. how
      ! the blocks of op1 are related to the blocks of ex1 and cnt
      ld_mmap1  = max(njoined_op(1),len_list)
      if (self) ld_mmap1 = ld_mmap1*2
      ld_mmap2  = max(njoined_op(2),len_list)
      ld_mmap12 = max(njoined_op(1),njoined_op(2))
      allocate(merge_map_op1(ld_mmap1,2,njoined_op(1)),
     &         merge_map_op2(ld_mmap2,2,njoined_op(2)),
     &         merge_map_op1op2(ld_mmap12,2,nvtx))

      njoined_cnt = len_list
      if (self) njoined_cnt = 2*len_list
      call occ_op2ex(iocc_ex1,iocc_cnt,merge_map_op1,
     &               .true.,.true.,ld_mmap1,
     &               1,iocc_op1,njoined_op(1),isvtx1,
     &               contr,arc_list,len_list)
      if (.not.self) 
     &  call occ_op2ex(iocc_ex2,iocc_cnt,merge_map_op2,
     &               .false.,.true.,ld_mmap2,
     &               2,iocc_op2,njoined_op(2),isvtx2,
     &               contr,arc_list,len_list)

      call condense_merge_map(merge_op1,
     &                   merge_map_op1,ld_mmap1,njoined_op(1),.false.)
      if (.not.self)
     &     call condense_merge_map(merge_op2,
     &                   merge_map_op2,ld_mmap2,njoined_op(2),.false.)

      allocate(ireo_vtx_no(nvtx),ireo_vtx_on(nvtx),
     &         ireo_after_contr(nvtx),
     &         ivtx_op1op2(nvtx))

      if (update_idxintm) idxintm = idxintm-1

      call reduce_contr2(sh_sign,cnt_sign,itf_sign,
     &     iocc_op1op2,njoined_op1op2,
     &     ireo_vtx_no,ireo_vtx_on,ireo_after_contr,
     &     ivtx_op1op2,nvtx_red,
     &     merge_map_op1op2,ld_mmap12,
     &     make_contr_red,contr_red,idxintm,
     &     contr,isvtx1,isvtx2,arc_list,len_list,njoined_res,
     &     reo_info)

      call condense_merge_map(merge_op1op2,
     &                merge_map_op1op2,ld_mmap12,njoined_op1op2,.false.)
      ! the same for EX2/EX1 sequence
      call condense_merge_map(merge_op2op1,
     &                merge_map_op1op2,ld_mmap12,njoined_op1op2,.true.)

      ! final result? get restr from op_info
      if (contr%narc.eq.len_list) then
        idx = contr%idx_res
        iblk = contr%iblk_res        
        if (make_contr_red.and.reo_info%nreo.gt.0) then
          ! since (in the case of zero vtxs) the vertex order need not be
          ! the same as in the original operator, we need the svmap:
          allocate(svmap(nvtx_red))
          call svmap4contr2(svmap,contr_red,ok)
          if (.not.ok) call quit(1,'get_bc_info3',
     &        'final result should have unique svmap!')
c dbg
c        write(lulog,*) 'op1op2 incl. restrictions (1):'
c        do ivtx = 1, njoined_res
c          call wrt_occ_rstr(lulog,ivtx,iocc_op1op2(1,1,ivtx),
c     &          op_info%op_arr(idx)%op%igasca_restr(1,1,1,1,1,
c     &            (iblk-1)*njoined_res+ivtx),
c     &          orb_info%ngas,orb_info%nspin)
c        end do
c dbgend
          do ivtx = 1, nvtx_red
            if (svmap(ivtx).eq.0) then
              irestr_op1op2(:,:,:,:,ivtx) = 0
            else
              irestr_op1op2(:,:,:,:,ivtx)
     &         = op_info%op_arr(idx)%op%igasca_restr(:,:,:,:,1,
     &            (iblk-1)*njoined_res+svmap(ivtx))
            end if
          end do
          deallocate(svmap)
        else
          if (njoined_op1op2.gt.njoined_res)
     &         irestr_op1op2(:,:,:,:,njoined_res+1:njoined_op1op2) = 0
          irestr_op1op2(:,:,:,:,1:njoined_res)
     &         = op_info%op_arr(idx)%op%igasca_restr(:,:,:,:,1,
     &            (iblk-1)*njoined_res+1:
     &            (iblk-1)*njoined_res+njoined_res)
        end if
c dbg
c        print *,'fetching restr op1op2 from op_info!'
c        do idx = 1, njoined_op1op2
c          call wrt_occ_rstr(lulog,idx,iocc_op1op2(1,1,idx),
c     &                                irestr_op1op2(1,1,1,1,idx),
c     &          orb_info%ngas,orb_info%nspin)
c        end do
c dbg
      else

        call special_restr(irestr_op1op2,
     &     iocc_op1op2,njoined_op1op2,
     &     merge_op1op2,
     &     iocc_op1,iocc_ex1,irestr_op1,njoined_op(1),
     &     iocc_op2,iocc_ex2,irestr_op2,njoined_op(2),
     &     orb_info%ihpvgas,orb_info%nspin,orb_info%ngas)

      end if

      mst_op1op2 = mst_op(1) + mst_op(2)
      gamt_op1op2 = multd2h(gamt_op(1),gamt_op(2))

      ! check that we arrived at the correct result 
      ! (after last contraction):
      if (contr%narc.eq.len_list) then
c        if (njoined_op1op2.ne.njoined_res) then
c          write(lulog,*) 'njoined_op1op2, njoined_res: ',
c     &                    njoined_op1op2, njoined_res
c          call quit(1,'get_bc_info3','trap 1')
c        end if
        if (mst_op1op2.ne.info_vtx(1,1)) then
          write(lulog,*)'Inconsistency:'
          write(lulog,*)'mst_op(1:2):   ',mst_op(1:2),' -> ',mst_op1op2
          write(lulog,*)'info_vtx(1,1): ',info_vtx(1,1)
          write(lulog,*)'Possible reason: '//
     &         'Inconsistently defined spin states of ME lists!'
          call quit(1,'get_bc_info3','trap 2')
        end if
        if (gamt_op1op2.ne.info_vtx(2,1)) then
          write(lulog,*)'Inconsistency: '
          write(lulog,*)'gamt_op(1:2): ',gamt_op(1:2),' -> ',gamt_op1op2
          write(lulog,*)'info_vtx(2,1): ',info_vtx(2,1)
          write(lulog,*)'Possible reason: '//
     &         'Inconsistently defined IRREPs of ME lists!'
          call quit(1,'get_bc_info3','trap 3')
        end if
        tra_op1op2 = contr%dagger 
      else
        tra_op1op2 = .false.
      end if

      if (make_contr_red) then
        call occvtx_from_arcs(0,occ_vtx_red,contr_red,njoined_res)
        call reduce_contr_info(
     &       irestr_vtx_red,info_vtx_red,
     &       irestr_vtx,info_vtx,
     &       irestr_op1op2,mst_op1op2,gamt_op1op2,
     &       ireo_vtx_no,ivtx_op1op2,
     &       nvtx,nvtx_red,njoined_op1op2,njoined_res,ngas)
c        call reduce_fact_info(contr_red,contr,idx_contr+1,ireo_vtx_on)
      end if

      possible = .true.
      ! if this is not the last contraction and
      ! joined vertices exist: test for possible
      ! reorderings of vertices
      if (make_contr_red .and.
     &    contr%narc.ne.len_list .and.
     &    contr_red%nsupvtx.lt.nvtx
     &    ) then
        idxnew_op1op2 = idxintm
        call reorder_supvtx(possible,
     &       .true.,set_reo,.false.,reo_info,
     &       contr_red,occ_vtx_red(1,1,njoined_res+1),
     &              irestr_vtx_red(1,1,1,1,njoined_res+1),idxnew_op1op2,
     &       orb_info)
        if (contr%nxarc.gt.0)
     &       call reorder_supvtx_x(possible,
     &         .true.,set_reo,.false.,reo_info,
     &         contr_red,occ_vtx_red(1,1,njoined_res+1),
     &              irestr_vtx_red(1,1,1,1,njoined_res+1),idxnew_op1op2,
     &         orb_info)
        ! FIX - unclear, whether 2x reo to same vertex works
        ! also beware of consecutive reorderings (e.g. 1->2,2->3)
        ! other way round is o.k. (e.g. 1->2,3->1)
        do ireo = 1, reo_info%nreo
          do jreo = ireo+1, reo_info%nreo
            if (reo_info%reo(ireo)%from.eq.reo_info%reo(jreo)%from.or.
     &          reo_info%reo(ireo)%to  .eq.reo_info%reo(jreo)%to.or.
     &          reo_info%reo(ireo)%to  .eq.reo_info%reo(jreo)%from) then
              ! only forbid if non-zero overlap:
              if (iocc_nonzero(iocc_overlap(
     &            reo_info%reo(ireo)%occ_shift,.false.,
     &            reo_info%reo(jreo)%occ_shift,.false.)))
     &           possible = .false.
            end if
          end do
        end do

c        ! to be commented out
c        if (reo_info%nreo.gt.0) then
c          do ivtx = 1, nvtx
c            call dummy_restr(irestr_vtx_red(1,1,1,1,ivtx+njoined_res),
c     &           occ_vtx_red(1,1,ivtx+njoined_res),1,
c     &           orb_info)
c          end do
c        end if
        if (set_reo) call tidy_reo_info(reo_info)

      end if

      ! another FIX: skip (currently?) impossible contractions
      possible = possible.and.
     &           allowed_contr(contr,arc_list(1:len_list),len_list)

      bc_sign = dble(cnt_sign)
      bc_sign_itf = dble(itf_sign)

      deallocate(arc_list)
      deallocate(ireo_vtx_no,ireo_vtx_on,ivtx_op1op2)

      deallocate(ireo_after_contr)

      deallocate(merge_map_op1op2)
      deallocate(merge_map_op1,merge_map_op2)

      call dealloc_contr(contr)
      deallocate(contr)

      if (ntest.ge.100) then
        write(lulog,*) 'get_bc_info3 on exit:'
        write(lulog,*) 'idx_op/blk 1: ',idx_op(1), iblk_op(1)
        write(lulog,*) 'idx_op/blk 2: ',idx_op(2), iblk_op(2)
        write(lulog,*) 'transpose:    ',tra_op1, tra_op2, tra_op1op2
        write(lulog,*) 'MS:           ',mst_op(1), mst_op(2), mst_op1op2
        write(lulog,*) 'IRREP:        ',gamt_op(1), gamt_op(2),
     &                                                       gamt_op1op2
        write(lulog,*) 'sign: ',bc_sign
        write(lulog,*) 'sign (ITF): ',bc_sign_itf
        if (.not.self) write(lulog,*) 'op1, op2, op1op2:'
        if (     self) write(lulog,*) 'op1, tr(op1):'
        call wrt_occ_n(lulog,iocc_op1,njoined_op(1))
        call wrt_occ_n(lulog,iocc_op2,njoined_op(2))
        call wrt_occ_n(lulog,iocc_op1op2,njoined_op1op2)
        if (.not.self) write(lulog,*) 'ex1, ex2, cnt:'
        if (     self) write(lulog,*) 'ex1, cnt:'
        call wrt_occ_n(lulog,iocc_ex1,njoined_op(1))
        call wrt_occ_n(lulog,iocc_ex2,njoined_op(2))
        call wrt_occ_n(lulog,iocc_cnt,njoined_cnt)
        write(lulog,*) 'op1 incl. restrictions:'
        do idx = 1, njoined_op(1)
          call wrt_occ_rstr(lulog,idx,iocc_op1(1,1,idx),
     &                                irestr_op1(1,1,1,1,idx),
     &          orb_info%ngas,orb_info%nspin)
        end do
        write(lulog,*) 'op2 incl. restrictions:'
        do idx = 1, njoined_op(2)
          call wrt_occ_rstr(lulog,idx,iocc_op2(1,1,idx),
     &                                irestr_op2(1,1,1,1,idx),
     &          orb_info%ngas,orb_info%nspin)
        end do
        write(lulog,*) 'op1op2 incl. restrictions:'
        do idx = 1, njoined_op1op2
          call wrt_occ_rstr(lulog,idx,iocc_op1op2(1,1,idx),
     &                                irestr_op1op2(1,1,1,1,idx),
     &          orb_info%ngas,orb_info%nspin)
        end do
      end if

      return
      end
