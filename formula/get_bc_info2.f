*----------------------------------------------------------------------*
      subroutine get_bc_info2(bc_sign,
     &     idx_op,iblk_op,
     &     iocc_ex1,iocc_ex2,iocc_cnt,
     &     iocc_op1,iocc_op2,iocc_op1op2,
     &     irestr_op1,irestr_op2,irestr_op1op2,
     &     mst_op,mst_op1op2,
     &     igamt_op,igamt_op1op2,
     &     njoined_op, njoined_op1op2, njoined_cnt,
     &     merge_op1, merge_op2, merge_op1op2, merge_op2op1,
     &     contr,njoined_res,occ_vtx,irestr_vtx,info_vtx,iarc_red,
     &     irestr_res,ihpvgas,ngas)
*----------------------------------------------------------------------*
*     set info for binary contraction
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'ifc_operators.h'
      include 'def_contraction.h'
      include 'multd2h.h'

      type(contraction), intent(in) ::
     &     contr
      integer, intent(in) ::
     &     ngas, iarc_red, njoined_res
      integer, intent(in) ::
     &     occ_vtx(ngastp,2,contr%nvtx+njoined_res),
     &     irestr_vtx(2,ngas,2,2,contr%nvtx+njoined_res),
     &     info_vtx(2,contr%nvtx+njoined_res),
     &     irestr_res(2,ngas,2,2,njoined_res),
     &     ihpvgas(ngas)
      real(8), intent(out) ::
     &     bc_sign
      integer, intent(out) ::
     &     idx_op(2), iblk_op(2),
     &     njoined_op(2), njoined_op1op2, njoined_cnt,
     &     iocc_ex1(ngastp,2,*), iocc_ex2(ngastp,2,*),
     &     iocc_cnt(ngastp,2,*),
     &     iocc_op1(ngastp,2,*), iocc_op2(ngastp,2,*),
     &     iocc_op1op2(ngastp,2,*),
     &     irestr_op1(2,ngas,2,2,*), irestr_op2(2,ngas,2,2,*),
     &     irestr_op1op2(2,ngas,2,2,*),
     &     mst_op(2), mst_op1op2, igamt_op(2), igamt_op1op2,
     &     merge_op1(*), merge_op2(*), merge_op1op2(*), merge_op2op1(*)

      logical ::
     &     merge
      integer ::
     &     ivtx1, ivtx2, ivtx, ivtxsuper1, ivtxsuper2,
     &     len_list, ilist, idum, ijoin,
     &     ld_mmap1, ld_mmap2, ld_mmap12
      integer ::
     &     arc_list(contr%narc), inum_ori(2,contr%nvtx),
     &     iocc_ex1ex2(ngastp,2,contr%nvtx)
      integer, pointer ::
     &     merge_map_op1(:,:,:), merge_map_op2(:,:,:),
     &     merge_map_op1op2(:,:,:)

      integer, external ::
     &     imltlist

      ! set up operator info
c dbg
c      print *,'current contraction ',iarc_red
c      call prt_contr3(luout,contr,occ_vtx(1,1,njoined_res+1))
c dbg
      ! set up operator 1 and 2
      ivtx1 = contr%arc(iarc_red)%link(1)
      ivtx2 = contr%arc(iarc_red)%link(2)

      ivtxsuper1 = contr%svertex(ivtx1)
      njoined_op(1) = imltlist(ivtxsuper1,contr%svertex,contr%nvtx,1)
      if (ivtx1.le.contr%nvtx) then
        idx_op(1) = contr%vertex(ivtx1)%idx_op
        iblk_op(1) = (contr%vertex(ivtx1)%iblk_op-1)/njoined_op(1) + 1
      end if

      ivtxsuper2 = contr%svertex(ivtx2)
      njoined_op(2) = imltlist(ivtxsuper2,contr%svertex,contr%nvtx,1)
      if (ivtx2.le.contr%nvtx) then
        idx_op(2) = contr%vertex(ivtx2)%idx_op
        iblk_op(2) = (contr%vertex(ivtx2)%iblk_op-1)/njoined_op(2) + 1
      end if

      if (ivtxsuper1.eq.ivtxsuper2)
     &     call quit(1,'get_bc_info','I am confused ....')

      mst_op(1) = info_vtx(1,ivtx1+njoined_res)
      mst_op(2) = info_vtx(1,ivtx2+njoined_res)
      igamt_op(1) = info_vtx(2,ivtx1+njoined_res)
      igamt_op(2) = info_vtx(2,ivtx2+njoined_res)

      ivtx1 = 0
      ivtx2 = 0
      do ivtx = 1, contr%nvtx
        if (contr%svertex(ivtx).eq.ivtxsuper1) then
          ivtx1 = ivtx1+1
          iocc_op1(1:ngastp,1:2,ivtx1) =
     &         occ_vtx(1:ngastp,1:2,ivtx+njoined_res)
          irestr_op1(1:2,1:ngas,1:2,1:2,ivtx1) =
     &         irestr_vtx(1:2,1:ngas,1:2,1:2,ivtx+njoined_res)
        end if
        if (contr%svertex(ivtx).eq.ivtxsuper2) then
          ivtx2 = ivtx2+1
          iocc_op2(1:ngastp,1:2,ivtx2) =
     &         occ_vtx(1:ngastp,1:2,ivtx+njoined_res)
          irestr_op2(1:2,1:ngas,1:2,1:2,ivtx2) =
     &         irestr_vtx(1:2,1:ngas,1:2,1:2,ivtx+njoined_res)
        end if
      end do
      
      ! set external and contraction indices
      ! get all involved arcs
      call get_associated_arcs(arc_list,len_list,iarc_red,contr)

      ! set merging info for op1 and op2, i.e. how
      ! the blocks of op1 are related to the blocks of ex1 and cnt
      ld_mmap1  = max(njoined_op(1),len_list)
      ld_mmap2  = max(njoined_op(2),len_list)
      ld_mmap12 = max(njoined_op(1),njoined_op(2))
      allocate(merge_map_op1(ld_mmap1,2,njoined_op(1)),
     &         merge_map_op2(ld_mmap2,2,njoined_op(2)))

      njoined_cnt = len_list
      call occ_op2ex(iocc_ex1,iocc_cnt,merge_map_op1,
     &               .true.,.true.,ld_mmap1,
     &               1,iocc_op1,njoined_op(1),ivtxsuper1,
     &               contr,arc_list,len_list)
      call occ_op2ex(iocc_ex2,iocc_cnt,merge_map_op2,
     &               .false.,.true.,ld_mmap2,
     &               2,iocc_op2,njoined_op(2),ivtxsuper2,
     &               contr,arc_list,len_list)

      call condense_merge_map(merge_op1,
     &                   merge_map_op1,ld_mmap1,njoined_op(1),.false.)
      call condense_merge_map(merge_op2,
     &                   merge_map_op2,ld_mmap2,njoined_op(2),.false.)

      ! set up EX1/EX2 in correct order
      call occ_ex1ex2(iocc_ex1ex2,inum_ori,
     &                iocc_ex1,iocc_ex2,
     &                njoined_op(1),njoined_op(2),
     &                ivtxsuper1,ivtxsuper2,
     &                contr%svertex,contr%nvtx)

      allocate(merge_map_op1op2(ld_mmap12,2,contr%nvtx))
      ! merge EX1/E2 + set up merging info
      call merge_ex1ex2(iocc_op1op2,njoined_op1op2,merge_map_op1op2,
     &                ld_mmap12,
     &                ivtxsuper1,ivtxsuper2,
     &                iocc_ex1ex2,inum_ori,njoined_op(1)+njoined_op(2),
     &                arc_list,len_list,
     &                contr,occ_vtx,njoined_res)

      call dummy_restr(irestr_op1op2,
     &       iocc_op1op2,njoined_op1op2,ihpvgas,ngas)
      mst_op1op2 = mst_op(1) + mst_op(2)
      igamt_op1op2 = multd2h(igamt_op(1),igamt_op(2))

      call condense_merge_map(merge_op1op2,
     &                merge_map_op1op2,ld_mmap12,njoined_op1op2,.false.)
        ! the same for EX2/EX1 sequence
      call condense_merge_map(merge_op2op1,
     &                merge_map_op1op2,ld_mmap12,njoined_op1op2,.true.)
 
      ! check that we arrived at the correct result 
      ! (after last contraction):
      if (contr%narc.eq.len_list) then
        if (njoined_op1op2.ne.njoined_res)
     &       call quit(1,'get_bc_info2','trap 1')
        if (mst_op1op2.ne.info_vtx(1,1))
     &       call quit(1,'get_bc_info2','trap 2')
        if (igamt_op1op2.ne.info_vtx(2,1))
     &       call quit(1,'get_bc_info2','trap 3')
      end if

      ! calculate sign
      call sign_bc(bc_sign,
     &     ivtxsuper1,ivtxsuper2,contr%svertex,contr%nvtx,
     &     iocc_op1,iocc_op2,iocc_cnt,
     &     njoined_op(1),njoined_op(2),njoined_op1op2,njoined_cnt,
     &     merge_map_op1,merge_map_op2,merge_map_op1op2,
     &     ld_mmap1,ld_mmap2,ld_mmap12)

      deallocate(merge_map_op1op2)
      deallocate(merge_map_op1,merge_map_op2)

      return
      end
