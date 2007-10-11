*----------------------------------------------------------------------*
      subroutine fact_cost2(possible,cost,iscale,
     &     contr,njoined_res,occ_vtx,irestr_vtx,info_vtx,iarc,
     &     op_info,str_info,ihpvgas,ngas,nsym)
*----------------------------------------------------------------------*
*     calculate the computational cost of a given contraction
*     for each contraction step, the number of flops is summed up
*     and the maximum size of the occurring intermediates is evaluated
*     returned on cost(1:3): cost(1)<-flops, 
*                            cost(2)<-memory(tot)
*                            cost(3)<-memory(block)
*     the maximum scaling of contraction and intermediates is returned
*     on iscale(ngastp,1:2): iscale(,1)<-contr, iscale(,2)<-mem
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'def_contraction.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'mdef_operator_info.h'
      include 'ifc_operators.h'
      include 'multd2h.h'
      
      integer, parameter ::
     &     ntest = 00
      
      integer, intent(in) ::
     &     ngas, ihpvgas(ngas), nsym, iarc, njoined_res
      type(contraction), intent(in) ::
     &     contr
      integer, intent(in) ::
     &     occ_vtx(ngastp,2,contr%nvtx),
     &     irestr_vtx(2,ngas,2,2,contr%nvtx),
     &     info_vtx(2,contr%nvtx)
      type(operator_info), intent(in) ::
     &     op_info
      type(strinf), intent(in) ::
     &     str_info
      real(8), intent(out) ::
     &     cost(3)
      integer, intent(out) ::
     &     iscale(ngastp,2)
      logical, intent(out) ::
     &     possible

      logical, parameter ::
     &     new_route = .true.
      integer ::
     &     nvtx, narc,
     &     np_op1op2, nh_op1op2, nx_op1op2, np_cnt, nh_cnt, nx_cnt
      integer ::
     &     iocc_cnt(ngastp,2,contr%nvtx),
     &     iocc_ex1(ngastp,2,contr%nvtx),
     &     iocc_ex2(ngastp,2,contr%nvtx),
     &     irst_res(2,ngas,2,2,contr%nvtx),
     &     iocc_op1(ngastp,2,contr%nvtx),
     &     iocc_op2(ngastp,2,contr%nvtx),
     &     iocc_op1op2(ngastp,2,contr%nvtx),
     &     irst_op1op2(2,ngas,2,2,contr%nvtx),
     &     igr_op1op2(ngastp,2,contr%nvtx),
     &     irst_op1(2,ngas,2,2,contr%nvtx),
     &     irst_op2(2,ngas,2,2,contr%nvtx), 
     &     mst_op(2), mst_op1op2,
     &     igamt_op(2), igamt_op1op2,
     &     njoined_op(2), njoined_op1op2, njoined_cnt,
     &     idar1(2), idar2(2),
     &     merge_op1(contr%nvtx*contr%nvtx), ! a bit too large, I guess ...
     &     merge_op2(contr%nvtx*contr%nvtx),
     &     merge_op1op2(contr%nvtx*contr%nvtx),
     &     merge_op2op1(contr%nvtx*contr%nvtx),
     &     nca_blk(2,7)
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
     &     flops, xmemtot, xmemblk, bc_sign

      integer, external ::
     &     idxlist, int_expand, int_pack, maxxlvl_op
      logical, external ::
     &     check_grph4occ
c dbg
c      integer idx1, idx2
c dbg
c dbg
c check whether
c      include 'par_opnames_gen.h'
c      integer idxtbar, ivtx
c      integer, external :: idx_oplist2
c      logical tbarda
c
c      idxtbar = idx_oplist2(op_tbar,op_info)
c      tbarda = .false.
c      do ivtx = 1, contr%nvtx
c        tbarda = tbarda.or.contr%vertex(ivtx)%idx_op.eq.idxtbar
c      end do
c      if (tbarda.and.contr%nsupvtx.gt.2) then
c        if (contr%vertex(contr%arc(iarc)%link(1))%idx_op.eq.idxtbar .or.
c     &      contr%vertex(contr%arc(iarc)%link(2))%idx_op.eq.idxtbar)then
c          
c          print *,'avoiding early Tbar contr ... arc = ',iarc
c          possible = .false.
c          return
c
c        end if
c      end if
c      if (tbarda) print *,'allowing arc = ',iarc
c dbg


      if (ntest.gt.0) then
        call write_title(luout,wst_dbg_subr,'this is fact_cost')
      end if

      nvtx = contr%nvtx
      narc = contr%narc

      ! preliminary version (valid for standard CC only):
      ! restriction on result:
c dbg
c      print *,'call to set_restr_prel'
c dbg
      call set_restr_prel(irst_res,contr,op_info,ihpvgas,ngas)
 
      call get_bc_info2(bc_sign,
     &     idar1,idar2,
     &     iocc_ex1,iocc_ex2,iocc_cnt,
     &     iocc_op1,iocc_op2,iocc_op1op2,
     &     irst_op1,irst_op2,irst_op1op2,
     &     mst_op,mst_op1op2,
     &     igamt_op,igamt_op1op2,
     &     njoined_op, njoined_op1op2, njoined_cnt,
     &     merge_op1,merge_op2,merge_op1op2,merge_op2op1,
     &     contr,njoined_res,occ_vtx,irestr_vtx,info_vtx,iarc,
     &     irst_res,ihpvgas,ngas)

      ! count particle, hole, (active) spaces involved:
      ! in intermediate
      nh_op1op2 = sum(iocc_ex1(ihole,1:2,1:njoined_op(1)))
     &          + sum(iocc_ex2(ihole,1:2,1:njoined_op(2)))
      np_op1op2 = sum(iocc_ex1(ipart,1:2,1:njoined_op(1)))
     &          + sum(iocc_ex2(ipart,1:2,1:njoined_op(2)))
      nx_op1op2 = sum(iocc_ex1(iextr,1:2,1:njoined_op(1)))
     &          + sum(iocc_ex2(iextr,1:2,1:njoined_op(2)))
      ! in contraction length
      nh_cnt = sum(iocc_cnt(ihole,1:2,1:njoined_cnt))
      np_cnt = sum(iocc_cnt(ipart,1:2,1:njoined_cnt))
      nx_cnt = sum(iocc_cnt(iextr,1:2,1:njoined_cnt))

      if (ntest.ge.100) then
        write(luout,*) 'op 1:'
        call wrt_occ_n(luout,iocc_op1,njoined_op(1))
        write(luout,*) 'op 2:'
        call wrt_occ_n(luout,iocc_op2,njoined_op(2))
        write(luout,*) 'externals 1:'
        call wrt_occ_n(luout,iocc_ex1,njoined_op(1))
        write(luout,*) 'externals 2:'
        call wrt_occ_n(luout,iocc_ex2,njoined_op(2))
        write(luout,*) 'contraction:'
        call wrt_occ_n(luout,iocc_cnt,njoined_cnt)
        write(luout,*) 'intermediate/result:'
        call wrt_occ_n(luout,iocc_op1op2,njoined_op1op2)
      end if

      if (ntest.ge.50) then
        write(luout,'(x,a,i2,a,i2,a)')
     &         'Contraction scales as  H^{',nh_op1op2+nh_cnt,
     &                                   '}P^{',np_op1op2+np_cnt,'}'
        write(luout,'(x,a,i2,a,i2,a)')
     &             'Intermediate scales as H^{',nh_op1op2,
     &                                   '}P^{',np_op1op2,'}'
      end if

      ! set iscale
      ! maximum is taken over total scaling, so H^2 .gt. P^1 
      if ( nh_op1op2+np_op1op2+nx_op1op2+nh_cnt+np_cnt+nx_cnt.gt.
     &           iscale(ihole,1)+iscale(ipart,1)+iscale(iextr,1) .or.
     &    (nh_op1op2+np_op1op2+nx_op1op2+nh_cnt+np_cnt+nx_cnt.eq.
     &           iscale(ihole,1)+iscale(ipart,1)+iscale(iextr,1).and.
     &     np_op1op2+np_cnt.gt.iscale(ipart,1) ) ) then
        iscale(ihole,1) = nh_op1op2+nh_cnt
        iscale(ipart,1) = np_op1op2+np_cnt
        iscale(iextr,1) = nx_op1op2+nx_cnt
      end if
      if ( nh_op1op2+np_op1op2+nx_op1op2.gt.
     &     iscale(ihole,2)+iscale(ipart,2)+iscale(iextr,2) .or.
     &    (nh_op1op2+np_op1op2+nx_op1op2.eq.
     &     iscale(ihole,2)+iscale(ipart,2)+iscale(iextr,2).and.
     &     np_op1op2.gt.iscale(ipart,2) ) ) then
        iscale(ihole,2) = nh_op1op2
        iscale(ipart,2) = np_op1op2
        iscale(iextr,2) = nx_op1op2
      end if

      possible = .true.     
      ! check whether intermediate can be addressed by
      ! the available graphs (preliminary fix)
      possible = possible.and.
     &     check_grph4occ(iocc_op1op2,irst_op1op2,
     &     str_info,ihpvgas,ngas,njoined_op1op2)
c dbg
      print *,'op1op2: ',possible
      call wrt_occ_n(6,iocc_op1op2,njoined_op1op2)
c dbg
      possible = possible.and.
     &     check_grph4occ(iocc_op1,irst_op1,
     &     str_info,ihpvgas,ngas,njoined_op(1))
c dbg
      print *,'op1: ',possible
      call wrt_occ_n(6,iocc_op1,njoined_op(1))
c dbg
      possible = possible.and.
     &     check_grph4occ(iocc_op2,irst_op2,
     &     str_info,ihpvgas,ngas,njoined_op(2))
c dbg
      print *,'op2: ',possible
      call wrt_occ_n(6,iocc_op2,njoined_op(2))
c dbg

      ! if not: do not allow this factorization
      if (.not.possible)
     &     cost(1:3) = huge(cost(1))

      if (possible) then
        if (new_route) then
          call get_num_subblk(nca_blk(1,1),nca_blk(2,1),
     &                        iocc_op1,njoined_op(1))
          call get_num_subblk(nca_blk(1,2),nca_blk(2,2),
     &                        iocc_op2,njoined_op(2))
          call get_num_subblk(nca_blk(1,3),nca_blk(2,3),iocc_op1op2,
     &                                                  njoined_op1op2)
          call get_num_subblk(nca_blk(1,4),nca_blk(2,4),
     &                        iocc_ex1,njoined_op(1))
          call get_num_subblk(nca_blk(1,5),nca_blk(2,5),
     &                        iocc_ex2,njoined_op(2))
          call get_num_subblk(nca_blk(1,6),nca_blk(2,6),
     &                        iocc_cnt,njoined_cnt)
          ! dummy setting for op1op2tmp (see contr_op1op2)
          nca_blk(1:2,7) = nca_blk(1:2,3)
          allocate(
     &         cinfo_op1c(nca_blk(1,1),3),cinfo_op1a(nca_blk(2,1),3),
     &         cinfo_op2c(nca_blk(1,2),3),cinfo_op2a(nca_blk(2,2),3),
     &         cinfo_op1op2c(nca_blk(1,3),3),
     &                                    cinfo_op1op2a(nca_blk(2,3),3),
     &         cinfo_ex1c(nca_blk(1,4),3),cinfo_ex1a(nca_blk(2,4),3),
     &         cinfo_ex2c(nca_blk(1,5),3),cinfo_ex2a(nca_blk(2,5),3),
     &         cinfo_cntc(nca_blk(1,6),3),cinfo_cnta(nca_blk(2,6),3))
          allocate(
     &         map_info1c(max(1,nca_blk(1,1)*2*
     &                                 (njoined_op(1)+njoined_cnt))),
     &         map_info1a(max(1,nca_blk(2,1)*2*
     &                                 (njoined_op(1)+njoined_cnt))),
     &         map_info2c(max(1,nca_blk(1,2)*2*
     &                                 (njoined_op(2)+njoined_cnt))),
     &         map_info2a(max(1,nca_blk(2,2)*2*
     &                                 (njoined_op(2)+njoined_cnt))),
     &         map_info12c(max(1,nca_blk(1,3)*2*
     &                               (njoined_op(1)+njoined_op(2)))),
     &         map_info12a(max(1,nca_blk(2,3)*2*
     &                               (njoined_op(1)+njoined_op(2)))))
          call condense_bc_info(
     &         cinfo_op1c, cinfo_op1a, cinfo_op2c, cinfo_op2a,
     &         cinfo_op1op2c, cinfo_op1op2a,
     &         cinfo_op1op2c, cinfo_op1op2a,
     &         cinfo_ex1c, cinfo_ex1a, cinfo_ex2c, cinfo_ex2a,
     &         cinfo_cntc, cinfo_cnta,
     &         map_info1c, map_info1a,
     &         map_info2c, map_info2a,
     &         map_info12c, map_info12a,
     &         nca_blk,
     &         iocc_op1, iocc_op2, iocc_op1op2, iocc_op1op2,
     &         iocc_ex1,iocc_ex2,iocc_cnt,
     &         irst_op1, irst_op2, irst_op1op2, irst_op1op2,
     &         merge_op1, merge_op2, merge_op1op2, merge_op2op1,
     &         njoined_op(1), njoined_op(2),njoined_op1op2, njoined_cnt,
     &         str_info,ihpvgas,ngas)
          call dummy_contr2(flops,xmemtot,xmemblk,
     &         nca_blk,
     &         cinfo_op1c, cinfo_op1a, cinfo_op2c, cinfo_op2a,
     &         cinfo_op1op2c, cinfo_op1op2a,
     &         cinfo_ex1c, cinfo_ex1a, cinfo_ex2c, cinfo_ex2a,
     &         cinfo_cntc, cinfo_cnta,
     &         map_info1c, map_info1a,
     &         map_info2c, map_info2a,
     &         map_info12c, map_info12a,
     &         mst_op(1),mst_op(2),mst_op1op2,
     &         igamt_op(1),igamt_op(2),igamt_op1op2,
     &         str_info,ngas,ihpvgas,nsym)          
          deallocate(
     &         cinfo_op1c,cinfo_op1a)
          deallocate(
     &         cinfo_op2c,cinfo_op2a)
          deallocate(
     &         cinfo_op1op2c)
          deallocate(
     &         cinfo_op1op2a)
          deallocate(
     &         cinfo_ex1c,cinfo_ex1a)
          deallocate(
     &         cinfo_ex2c,cinfo_ex2a)
          deallocate(
     &         cinfo_cntc,cinfo_cnta)
          deallocate(
     &         map_info1c)
          deallocate(
     &         map_info1a)
          deallocate(
     &         map_info2c)
          deallocate(
     &         map_info2a)
          deallocate(
     &         map_info12c)
          deallocate(
     &         map_info12a)
        else
          call dummy_contr(flops,xmemtot,xmemblk,
     &       iocc_op1,iocc_op2,iocc_ex1,iocc_ex2,
     &       iocc_op1op2,iocc_cnt,
     &       irst_op1,irst_op2,irst_op1op2,
     &       mst_op,mst_op1op2,igamt_op,igamt_op1op2,
     &       str_info,ngas,ihpvgas,nsym)
        end if

        cost(1) = cost(1)+flops
        cost(2) = max(cost(2),xmemtot)
        cost(3) = max(cost(3),xmemblk)

        if (ntest.ge.10)
     &       write(luout,'(x,a,3(g20.10,x))') '$$ cost: ',cost(1:3)
        if (ntest.ge.10)
     &       write(luout,'(x,a,2(2i4,2x))') '$$ scal: ',iscale(1:2,1:2)
      end if

      return
      end

