*----------------------------------------------------------------------*
      subroutine get_bc_info2(idx_op,iblk_op,
     &     iocc_ext1,iocc_ext2,iocc_cnt,
     &     iocc_op1,iocc_op2,iocc_op1op2,
     &     irestr_op1,irestr_op2,irestr_op1op2,
     &     mst_op,mst_op1op2,
     &     igamt_op,igamt_op1op2,
     &     njoined_op, njoined_op1op2, njoined_cnt,
     &     merge_op1, merge_op2, merge_op1op2,
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
      integer, intent(out) ::
     &     idx_op(2), iblk_op(2),
     &     njoined_op(2), njoined_op1op2, njoined_cnt,
     &     iocc_ext1(ngastp,2,*), iocc_ext2(ngastp,2,*),
     &     iocc_cnt(ngastp,2,*),
     &     iocc_op1(ngastp,2,*), iocc_op2(ngastp,2,*),
     &     iocc_op1op2(ngastp,2,*),
     &     irestr_op1(2,ngas,2,2,*), irestr_op2(2,ngas,2,2,*),
     &     irestr_op1op2(2,ngas,2,2,*),
     &     mst_op(2), mst_op1op2, igamt_op(2), igamt_op1op2,
     &     merge_op1(*), merge_op2(*), merge_op1op2(*)

      logical ::
     &     merge
      integer ::
     &     iarc, ivtx1, ivtx2, ivtx, ivtxsuper1, ivtxsuper2,
     &     idx, idx12, len_list, ilist, idum, ijoin,
     &     ld_mmap1, ld_mmap2, ld_mmap12, jdx 
      integer ::
     &     svmap(contr%nvtx), topomap(contr%nvtx,contr%nvtx),
     &     arc_list(contr%narc)
      integer, pointer ::
     &     merge_map_op1(:,:,:), merge_map_op2(:,:,:),
     &     merge_map_op1op2(:,:,:)

      integer, external ::
     &     imltlist
      logical, external ::
     &     merge_vtx1vtx2

      ! set up operator info
c dbg
c      print *,'current contraction'
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

      iocc_ext1(1:ngastp,1:2,1:njoined_op(1)) =
     &     iocc_op1(1:ngastp,1:2,1:njoined_op(1))
      iocc_ext2(1:ngastp,1:2,1:njoined_op(2)) =
     &     iocc_op2(1:ngastp,1:2,1:njoined_op(2))

      ! set merging info for op1 and op2, i.e. how
      ! the blocks of op1 are related to the blocks of ex1 and cnt
      ld_mmap1  = max(njoined_op(1),len_list)
      ld_mmap2  = max(njoined_op(2),len_list)
      ld_mmap12 = max(njoined_op(1),njoined_op(2))
      allocate(merge_map_op1(ld_mmap1,2,njoined_op(1)),
     &         merge_map_op2(ld_mmap2,2,njoined_op(2)))
      merge_map_op1 = 0
      merge_map_op2 = 0
      ! this is trivial: each block of ex1 corresponds to the
      ! block of op1 with the same vertex number:
      do ijoin = 1, njoined_op(1)
        merge_map_op1(1,1,ijoin) = ijoin 
      end do
      do ijoin = 1, njoined_op(2)
        merge_map_op2(1,1,ijoin) = ijoin 
      end do

      njoined_cnt = 0
      do ilist = 1, len_list
        iarc = arc_list(ilist)
c dbg
c        print *,'next iarc = ',iarc
c        print *,'current EX1/EX2:'
c        call wrt_occ_n(luout,iocc_ext1,njoined_op(1))
c        call wrt_occ_n(luout,iocc_ext2,njoined_op(2))
c        print *,'current CNT:'
c        call wrt_occ(luout,contr%arc(iarc)%occ_cnt)
c dbg

        ! current operator 1 and 2
        njoined_cnt = njoined_cnt+1
        if (contr%svertex(contr%arc(iarc)%link(1)).eq.ivtxsuper1) then
          ivtx1 = contr%arc(iarc)%link(1)
          ivtx2 = contr%arc(iarc)%link(2)
          iocc_cnt(1:ngastp,1:2,njoined_cnt) = contr%arc(iarc)%occ_cnt
        else
          ivtx2 = contr%arc(iarc)%link(1)
          ivtx1 = contr%arc(iarc)%link(2)
          iocc_cnt(1:ngastp,1:2,njoined_cnt) =
     &         iocc_dagger(contr%arc(iarc)%occ_cnt)
        end if
c dbg
c        print *,'ivtx1, ivtx2: ',ivtx1,ivtx2
c        print *,'-> super    : ',contr%svertex(ivtx1),
c     &       contr%svertex(ivtx2)
c        print *,'ivtxsuper1,ivtxsuper2: ',ivtxsuper1,ivtxsuper2 
c        print *,'svertex: ',contr%svertex(1:contr%nvtx)
c        print *,'updated CNT (last is current CNT)'
c        call wrt_occ_n(luout,iocc_cnt,njoined_cnt)
c dbg
      
        ! gives number of supervertices up to position ivtx1 == 
        !  # of primitive vertex in super-vertex
        idx = imltlist(ivtxsuper1,contr%svertex,ivtx1,1)
        iocc_ext1(1:ngastp,1:2,idx) = iocc_ext1(1:ngastp,1:2,idx)
     &                              - iocc_cnt(1:ngastp,1:2,njoined_cnt)
c dbg
c        print *,'1: idx = ',idx
c        print *,'resulting EXT1'
c        call wrt_occ(luout,iocc_ext1(1:ngastp,1:2,idx))
c dbg
        ! store merging info
        jdx = 1
        do while(merge_map_op1(jdx,2,idx).gt.0) ! get a free entry
          jdx = jdx+1
        end do
        merge_map_op1(jdx,2,idx) = njoined_cnt

        idx = imltlist(ivtxsuper2,contr%svertex,ivtx2,1)
        iocc_ext2(1:ngastp,1:2,idx) = iocc_ext2(1:ngastp,1:2,idx)
     &                 - iocc_dagger(iocc_cnt(1:ngastp,1:2,njoined_cnt))
c dbg
c        print *,'2: idx = ',idx
c        print *,'resulting EXT2'
c        call wrt_occ(luout,iocc_ext2(1:ngastp,1:2,idx))
c dbg

        ! store merging info
        jdx = 1
        do while(merge_map_op2(jdx,2,idx).gt.0)
          jdx = jdx+1
        end do
        merge_map_op2(jdx,2,idx) = njoined_cnt

      end do

      ! set intermediate/result
c      if (contr%narc.eq.len_list) then
c        njoined_op1op2 = njoined_res
c        iocc_op1op2(1:ngastp,1:2,1:njoined_res) =
c     &       occ_vtx(1:ngastp,1:2,1:njoined_res)
c        irestr_op1op2(1:2,1:ngas,1:2,1:2,1:njoined_res) =
c     &       irestr_vtx(1:2,1:ngas,1:2,1:2,1:njoined_res)
c        mst_op1op2 = info_vtx(1,1)
c        igamt_op1op2 = info_vtx(2,1)
c      else

        ! we have to find out whether we may merge the vertices
        ! after contraction ...
        ! get topology info
        if (njoined_res.eq.1) then
          svmap(1:contr%nvtx) = 1
        else
          call svmap4contr(svmap,contr,occ_vtx,njoined_res)
        end if

        call topomap4contr(1,topomap,idum,idum,idum,contr,
     &                       occ_vtx(1,1,njoined_res+1))

        allocate(merge_map_op1op2(ld_mmap12,2,contr%nvtx))
        merge_map_op1op2 = 0

        njoined_op1op2 = 
     &       imltlist(ivtxsuper1,contr%svertex,contr%nvtx,1)+
     &       imltlist(ivtxsuper2,contr%svertex,contr%nvtx,1)
        iocc_op1op2(1:ngastp,1:2,1:njoined_op1op2) = 0

        ! set up result occupations (merged, if possible)
        ! in case of merges, leave 0 occ for second vertex
        do ilist = 1, len_list
          iarc = arc_list(ilist)
          ivtx1 = contr%arc(iarc)%link(1)
          ivtx2 = contr%arc(iarc)%link(2)
c dbg
c          print *,'result vertices to be merged: ',ivtx1, ivtx2
c dbg
          
          merge = merge_vtx1vtx2(ivtx1,ivtx2,
     &                           contr%svertex,svmap,topomap,contr%nvtx)
c dbg
c          print *,'merge = ',merge
c          print *,'svmap: ',svmap(1:contr%nvtx)
c          print *,'topomap:'
c          call iwrtma(topomap,contr%nvtx,contr%nvtx,
c     &                        contr%nvtx,contr%nvtx)
c dbg

          idx = imltlist(ivtxsuper1,contr%svertex,ivtx1,1)
          idx12 = idx + imltlist(ivtxsuper2,contr%svertex,ivtx1,1)

          iocc_op1op2(1:ngastp,1:2,idx12) =
     &         iocc_ext1(1:ngastp,1:2,idx)

          jdx = 1
          do while(merge_map_op1op2(jdx,1,idx12).gt.0)
            jdx = jdx+1
          end do
          merge_map_op1op2(jdx,1,idx12) = idx

          if (merge) then
            idx = imltlist(ivtxsuper2,contr%svertex,ivtx2,1)
            iocc_op1op2(1:ngastp,1:2,idx12) =
     &         iocc_op1op2(1:ngastp,1:2,idx12) +
     &         iocc_ext2(1:ngastp,1:2,idx)

            jdx = 1
            do while(merge_map_op1op2(jdx,2,idx12).gt.0)
              jdx = jdx+1
            end do
            merge_map_op1op2(jdx,2,idx12) = idx

            idx12 = idx + imltlist(ivtxsuper2,contr%svertex,ivtx2,1)
            iocc_op1op2(1:ngastp,1:2,idx12) = 0

          else
            idx = imltlist(ivtxsuper2,contr%svertex,ivtx2,1)
            idx12 = idx + imltlist(ivtxsuper2,contr%svertex,ivtx2,1)
            iocc_op1op2(1:ngastp,1:2,idx12) =
     &         iocc_ext2(1:ngastp,1:2,idx)

            jdx = 1
            do while(merge_map_op1op2(jdx,2,idx12).gt.0)
              jdx = jdx+1
            end do
            merge_map_op1op2(jdx,2,idx12) = idx

          end if
        end do
          
c dbg
c        print *,'raw OP1OP2 occupation: '
c        call wrt_occ_n(luout,iocc_op1op2,
c     &       imltlist(ivtxsuper1,contr%svertex,contr%nvtx,1) +
c     &       imltlist(ivtxsuper2,contr%svertex,contr%nvtx,1))
c dbg
        ! remove zero occupations
        njoined_op1op2 = 0
        do ivtx = 1, imltlist(ivtxsuper1,contr%svertex,contr%nvtx,1) +
     &               imltlist(ivtxsuper2,contr%svertex,contr%nvtx,1)
          if (iocc_zero(iocc_op1op2(1:ngastp,1:2,ivtx))) cycle
          njoined_op1op2 = njoined_op1op2+1
          if (njoined_op1op2.lt.ivtx)
     &         iocc_op1op2(1:ngastp,1:2,njoined_op1op2) =
     &           iocc_op1op2(1:ngastp,1:2,ivtx)
        end do
        ! well, we keep at least one:
c dbg
c        if (njoined_op1op2.gt.1) print *,'HIER WIRD''S SPANNEND!'
c        print *,'final OP1OP2 occupation: '
c        call wrt_occ_n(luout,iocc_op1op2,njoined_op1op2)
c dbg
        njoined_op1op2 = max(1,njoined_op1op2)

        ! yes, the preliminary code
c        do ijoin = 1, njoined_op1op2
c          call fit_restr(irestr_op1op2(1,1,1,1,ijoin),
c     &         iocc_op1op2(1,1,ijoin),
c     &         irestr_res,ihpvgas,ngas)
c        end do
        call dummy_restr(irestr_op1op2,
     &       iocc_op1op2,njoined_op1op2,ihpvgas,ngas)
        mst_op1op2 = mst_op(1) + mst_op(2)
        igamt_op1op2 = multd2h(igamt_op(1),igamt_op(2))

        call condense_merge_map(merge_op1op2,
     &                        merge_map_op1op2,ld_mmap12,njoined_op1op2)
        deallocate(merge_map_op1op2)
 
c      end if

      call condense_merge_map(merge_op1,
     &                        merge_map_op1,ld_mmap1,njoined_op(1))
      call condense_merge_map(merge_op2,
     &                        merge_map_op2,ld_mmap2,njoined_op(2))

      deallocate(merge_map_op1,merge_map_op2)

      return
      end
