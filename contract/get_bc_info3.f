      subroutine get_bc_info3(bc_sign,possible,
     &     idx_op,iblk_op,
     &     iocc_ex1,iocc_ex2,iocc_cnt,
     &     iocc_op1,iocc_op2,iocc_op1op2,
     &     irestr_op1,irestr_op2,irestr_op1op2,
     &     tra_op1,tra_op2,tra_op1op2,
     &     mst_op,mst_op1op2,
     &     gamt_op,gamt_op1op2,
     &     njoined_op, njoined_op1op2, njoined_cnt,
     &     merge_op1, merge_op2, merge_op1op2, merge_op2op1,
     &     contr,occ_vtx,irestr_vtx,info_vtx,
     &     make_contr_red,
     &     contr_red,occ_vtx_red,irestr_vtx_red,info_vtx_red,
     &     set_reo, reo_info,
     &     iarc_contr,idx_contr,idxnew_op1op2,
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

      type(contraction), intent(in) ::
     &     contr
      type(contraction), intent(out) ::
     &     contr_red
      type(reorder_info), intent(inout) ::
     &     reo_info
      type(orbinf), intent(in), target ::
     &     orb_info
      type(operator_info), intent(in) ::
     &     op_info
      logical, intent(in) ::
     &     set_reo, make_contr_red
      logical, intent(out) ::
     &     possible
      integer, intent(in) ::
     &     iarc_contr, idx_contr, idxnew_op1op2, njoined_res
      integer, intent(in) ::
     &     occ_vtx(ngastp,2,contr%nvtx+njoined_res),
     &     irestr_vtx(2,orb_info%ngas,2,2,contr%nvtx+njoined_res),
     &     info_vtx(2,contr%nvtx+njoined_res),
     &     irestr_res(2,orb_info%ngas,2,2,njoined_res)
      integer, intent(out) ::
     &     occ_vtx_red(ngastp,2,*),
     &     irestr_vtx_red(2,orb_info%ngas,2,2,*),
     &     info_vtx_red(2,*)
      real(8), intent(out) ::
     &     bc_sign
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

      integer ::
     &     ld_mmap1, ld_mmap2, ld_mmap12, ngas,
     &     nvtx, ivtx, idx,
     &     ivtx1, ivtx2, isvtx1, isvtx2,
     &     len_list, nvtx_red

      integer, pointer ::
     &     ireo_vtx_no(:), ireo_vtx_on(:), ivtx_op1op2(:),
     &     arc_list(:),
     &     merge_map_op1(:,:,:), merge_map_op2(:,:,:),
     &     merge_map_op1op2(:,:,:),
     &     ihpvgas(:,:)

      integer, external ::
     &     imltlist

      if (ntest.ge.100) then
        call write_title(luout,wst_dbg_subr,'get_bc_info3')
        call prt_contr3(luout,contr,-1)
      end if

      ngas = orb_info%ngas
      ihpvgas => orb_info%ihpvgas

      nvtx = contr%nvtx

      ivtx1 = contr%arc(iarc_contr)%link(1)
      ivtx2 = contr%arc(iarc_contr)%link(2)

      isvtx1 = contr%svertex(ivtx1)
      njoined_op(1) = imltlist(isvtx1,contr%svertex,nvtx,1)
      if (ivtx1.le.nvtx) then
        idx_op(1) = contr%vertex(ivtx1)%idx_op
        iblk_op(1) = (contr%vertex(ivtx1)%iblk_op-1)/njoined_op(1) + 1
        tra_op1   = contr%vertex(ivtx1)%dagger
      end if

      isvtx2 = contr%svertex(ivtx2)
      njoined_op(2) = imltlist(isvtx2,contr%svertex,nvtx,1)
      if (ivtx2.le.nvtx) then
        idx_op(2) = contr%vertex(ivtx2)%idx_op
        iblk_op(2) = (contr%vertex(ivtx2)%iblk_op-1)/njoined_op(2) + 1
        tra_op2   = contr%vertex(ivtx2)%dagger
      end if

      if (isvtx1.eq.isvtx2)
     &     call quit(1,'get_bc_info3','I am confused ....')

      mst_op(1) = info_vtx(1,ivtx1+njoined_res)
      mst_op(2) = info_vtx(1,ivtx2+njoined_res)
      gamt_op(1) = info_vtx(2,ivtx1+njoined_res)
      gamt_op(2) = info_vtx(2,ivtx2+njoined_res)

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
        if (contr%svertex(ivtx).eq.isvtx2) then
          ivtx2 = ivtx2+1
          iocc_op2(1:ngastp,1:2,ivtx2) =
     &         occ_vtx(1:ngastp,1:2,ivtx+njoined_res)
          irestr_op2(1:2,1:ngas,1:2,1:2,ivtx2) =
     &         irestr_vtx(1:2,1:ngas,1:2,1:2,ivtx+njoined_res)
        end if
      end do
      
      allocate(arc_list(contr%narc))
      ! set external and contraction indices
      ! get all involved arcs
      call get_associated_arcs(arc_list,len_list,iarc_contr,contr)

      ! set merging info for op1 and op2, i.e. how
      ! the blocks of op1 are related to the blocks of ex1 and cnt
      ld_mmap1  = max(njoined_op(1),len_list)
      ld_mmap2  = max(njoined_op(2),len_list)
      ld_mmap12 = max(njoined_op(1),njoined_op(2))
      allocate(merge_map_op1(ld_mmap1,2,njoined_op(1)),
     &         merge_map_op2(ld_mmap2,2,njoined_op(2)),
     &         merge_map_op1op2(ld_mmap12,2,nvtx))

      njoined_cnt = len_list
      call occ_op2ex(iocc_ex1,iocc_cnt,merge_map_op1,
     &               .true.,.true.,ld_mmap1,
     &               1,iocc_op1,njoined_op(1),isvtx1,
     &               contr,arc_list,len_list)
      call occ_op2ex(iocc_ex2,iocc_cnt,merge_map_op2,
     &               .false.,.true.,ld_mmap2,
     &               2,iocc_op2,njoined_op(2),isvtx2,
     &               contr,arc_list,len_list)

      call condense_merge_map(merge_op1,
     &                   merge_map_op1,ld_mmap1,njoined_op(1),.false.)
      call condense_merge_map(merge_op2,
     &                   merge_map_op2,ld_mmap2,njoined_op(2),.false.)

      allocate(ireo_vtx_no(nvtx),ireo_vtx_on(nvtx),ivtx_op1op2(nvtx))

      call reduce_contr2(iocc_op1op2,njoined_op1op2,
     &     ireo_vtx_no,ireo_vtx_on,ivtx_op1op2,nvtx_red,
     &     merge_map_op1op2,ld_mmap12,
     &     make_contr_red,contr_red,idxnew_op1op2,
     &     contr,isvtx1,isvtx2,arc_list,len_list,njoined_res)

      call condense_merge_map(merge_op1op2,
     &                merge_map_op1op2,ld_mmap12,njoined_op1op2,.false.)
      ! the same for EX2/EX1 sequence
      call condense_merge_map(merge_op2op1,
     &                merge_map_op1op2,ld_mmap12,njoined_op1op2,.true.)

      call dummy_restr(irestr_op1op2,
     &       iocc_op1op2,njoined_op1op2,orb_info)
      mst_op1op2 = mst_op(1) + mst_op(2)
      gamt_op1op2 = multd2h(gamt_op(1),gamt_op(2))

      ! check that we arrived at the correct result 
      ! (after last contraction):
      if (contr%narc.eq.len_list) then
        if (njoined_op1op2.ne.njoined_res) then
          write(luout,*) 'njoined_op1op2, njoined_res: ',
     &                    njoined_op1op2, njoined_res
          call quit(1,'get_bc_info3','trap 1')
        end if
        if (mst_op1op2.ne.info_vtx(1,1))
     &       call quit(1,'get_bc_info3','trap 2')
        if (gamt_op1op2.ne.info_vtx(2,1))
     &       call quit(1,'get_bc_info3','trap 3')
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
        call reduce_fact_info(contr_red,contr,idx_contr+1,ireo_vtx_on)
      end if

      deallocate(arc_list)
      deallocate(ireo_vtx_no,ireo_vtx_on,ivtx_op1op2)

      possible = .true.
      ! if this is not the last contraction and
      ! joined vertices exist: test for possible
      ! reorderings of vertices
c dbg
c      print *,'contr%narc,len_list: ',contr%narc,len_list
c      print *,'contr_red%nsupvtx,nvtx:  ',contr_red%nsupvtx,nvtx
c dbg
      if (make_contr_red .and.
     &    contr%narc.ne.len_list .and.
     &    contr_red%nsupvtx.lt.nvtx
     &    ) then
        call reorder_supvtx(possible,
     &       .true.,set_reo,reo_info,
     &       contr_red,occ_vtx_red(1,1,njoined_res+1),idxnew_op1op2)
        if (contr%nxarc.gt.0)
     &       call reorder_supvtx_x(possible,
     &         .true.,set_reo,reo_info,
     &         contr_red,occ_vtx_red(1,1,njoined_res+1),idxnew_op1op2)
        do ivtx = 1, nvtx
          call dummy_restr(irestr_vtx_red(1,1,1,1,ivtx+njoined_res),
     &         occ_vtx_red(1,1,ivtx+njoined_res),1,
     &         orb_info)
        end do

      end if

      ! calculate sign
      call sign_bc(bc_sign,
     &     isvtx1,isvtx2,contr%svertex,nvtx,
     &     iocc_op1,iocc_op2,iocc_cnt,
     &     njoined_op(1),njoined_op(2),njoined_op1op2,njoined_cnt,
     &     merge_map_op1,merge_map_op2,merge_map_op1op2,
     &     ld_mmap1,ld_mmap2,ld_mmap12)

      deallocate(merge_map_op1op2)
      deallocate(merge_map_op1,merge_map_op2)

      if (ntest.ge.100) then
        write(luout,*) 'get_bc_info3 on exit:'
        write(luout,*) 'idx_op/blk 1: ',idx_op(1), iblk_op(1)
        write(luout,*) 'idx_op/blk 2: ',idx_op(2), iblk_op(2)
        write(luout,*) 'transpose:    ',tra_op1, tra_op2, tra_op1op2
        write(luout,*) 'MS:           ',mst_op(1), mst_op(2), mst_op1op2
        write(luout,*) 'IRREP:        ',gamt_op(1), gamt_op(2),
     &                                                       gamt_op1op2
        write(luout,*) 'op1, op2, op1op2:'
        call wrt_occ_n(luout,iocc_op1,njoined_op(1))
        call wrt_occ_n(luout,iocc_op2,njoined_op(2))
        call wrt_occ_n(luout,iocc_op1op2,njoined_op1op2)
        write(luout,*) 'ex1, ex2, cnt:'
        call wrt_occ_n(luout,iocc_ex1,njoined_op(1))
        call wrt_occ_n(luout,iocc_ex2,njoined_op(2))
        call wrt_occ_n(luout,iocc_cnt,njoined_cnt)
      end if

      return
      end
