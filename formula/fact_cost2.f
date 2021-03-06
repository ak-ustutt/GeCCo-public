*----------------------------------------------------------------------*
      subroutine fact_cost2(possible,cost,iscale,
     &     iarc,njoined_res,nlevel,
     &     contr,occ_vtx,irestr_vtx,info_vtx,
     &     make_red,
     &     contr_red,occ_vtx_red,irestr_vtx_red,info_vtx_red,
     &     op_info,str_info,orb_info)
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
      include 'def_orbinf.h'
      include 'mdef_operator_info.h'
      include 'ifc_operators.h'
      include 'def_contraction_info.h'
      include 'def_reorder_info.h'
      include 'multd2h.h'
      
      integer, parameter ::
     &     ntest = 0
      logical, parameter ::
     &     formal = .false., exact = .false.
      
      logical, intent(in) ::
     &     make_red
      integer, intent(in) ::
     &     iarc, njoined_res, nlevel
      type(contraction), intent(in) ::
     &     contr
      type(contraction), intent(out) ::
     &     contr_red
      type(orbinf), intent(in), target ::
     &     orb_info
      integer, intent(in) ::
     &     occ_vtx(ngastp,2,contr%nvtx+njoined_res),
     &     irestr_vtx(2,orb_info%ngas,2,2,contr%nvtx+njoined_res),
     &     info_vtx(2,contr%nvtx+njoined_res)
      integer, intent(out) ::
     &     occ_vtx_red(ngastp,2,contr%nvtx+njoined_res),
     &     irestr_vtx_red(2,orb_info%ngas,2,2,contr%nvtx+njoined_res),
     &     info_vtx_red(2,contr%nvtx+njoined_res)
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

      integer ::
     &     nvtx, narc, ngas, nsym, idum,
     &     np_op1op2, nh_op1op2, nx_op1op2, np_cnt, nh_cnt, nx_cnt,
     &     nv_op1op2, nv_cnt
      integer ::
     &     iocc_cnt(ngastp,2,2*contr%nvtx),
     &     iocc_ex1(ngastp,2,contr%nvtx),
     &     iocc_ex2(ngastp,2,contr%nvtx),
     &     irst_res(2,orb_info%ngas,2,2,contr%nvtx),
     &     iocc_op1(ngastp,2,contr%nvtx),
     &     iocc_op2(ngastp,2,contr%nvtx),
     &     iocc_op1op2(ngastp,2,contr%nvtx),
     &     irst_op1op2(2,orb_info%ngas,2,2,contr%nvtx),
     &     igr_op1op2(ngastp,2,contr%nvtx),
     &     irst_op1(2,orb_info%ngas,2,2,contr%nvtx),
     &     irst_op2(2,orb_info%ngas,2,2,contr%nvtx), 
     &     mst_op(2), mst_op1op2,
     &     igamt_op(2), igamt_op1op2,
     &     njoined_op(2), njoined_op1op2, njoined_cnt,
     &     idxop(2), idar2(2),
     &     merge_op1(2*contr%nvtx*contr%nvtx), ! a bit too large, I guess ...
     &     merge_op2(2*contr%nvtx*contr%nvtx),
     &     merge_op1op2(2*contr%nvtx*contr%nvtx),
     &     merge_op2op1(2*contr%nvtx*contr%nvtx),
     &     iscale_new(ngastp)
      integer ::
     &     itf_index_info
      
      integer, pointer ::
     &     ihpvgas(:,:)
      
      type(contraction) ::
     &     contr_dum
      type(reorder_info) ::
     &     reo_dum
      type(contraction_info) ::
     &     cnt_info

      logical ::
     &     tra_op1, tra_op2, tra_op1op2
      real(8) ::
     &     flops, xmemtot, xmemblk, bc_sign

      integer, external ::
     &     idxlist, int_expand, int_pack, maxxlvl_op
      logical, external ::
     &     check_grph4occ
      real(8), external ::
     &     scale_rank

      if (ntest.gt.0) then
        call write_title(lulog,wst_dbg_subr,'this is fact_cost')
c dbg
c        call prt_contr3(lulog,contr,-1)
c dbg
      end if

      nsym = orb_info%nsym
      ngas = orb_info%ngas
      ihpvgas => orb_info%ihpvgas

      nvtx = contr%nvtx
      narc = contr%narc

      ! preliminary version (valid for standard CC only):
      ! restriction on result:
      call set_restr_prel(irst_res,contr,op_info,ihpvgas,ngas)

      itf_index_info = -1
      
      call get_bc_info3(bc_sign,idum,possible,
     &     idxop,idar2,
     &     iocc_ex1,iocc_ex2,iocc_cnt,
     &     iocc_op1,iocc_op2,iocc_op1op2,
     &     irst_op1,irst_op2,irst_op1op2,
     &     tra_op1,tra_op2,tra_op1op2,
     &     mst_op,mst_op1op2,
     &     igamt_op,igamt_op1op2,
     &     njoined_op, njoined_op1op2, njoined_cnt,
     &     merge_op1,merge_op2,merge_op1op2,merge_op2op1,
     &     itf_index_info,
     &     contr,occ_vtx,irestr_vtx,info_vtx,
     &     make_red,
     &     contr_red,occ_vtx_red,irestr_vtx_red,info_vtx_red,
     &     .false.,reo_dum,reo_dum,
     &     iarc,.false.,-nlevel,
     &     irst_res,njoined_res,orb_info,op_info,0)
c dbg
      print *,'iocc_op1op2 fresh from get_bc_info3'
      call wrt_occ_n(6,iocc_op1op2,njoined_op1op2)
c dbg

      ! count particle, hole, (active) spaces involved:
      ! in intermediate
      nh_op1op2 = sum(iocc_ex1(ihole,1:2,1:njoined_op(1)))
     &          + sum(iocc_ex2(ihole,1:2,1:njoined_op(2)))
      np_op1op2 = sum(iocc_ex1(ipart,1:2,1:njoined_op(1)))
     &          + sum(iocc_ex2(ipart,1:2,1:njoined_op(2)))
      nx_op1op2 = sum(iocc_ex1(iextr,1:2,1:njoined_op(1)))
     &          + sum(iocc_ex2(iextr,1:2,1:njoined_op(2)))
      nv_op1op2 = sum(iocc_ex1(ivale,1:2,1:njoined_op(1)))
     &          + sum(iocc_ex2(ivale,1:2,1:njoined_op(2)))
      ! in contraction length
      nh_cnt = sum(iocc_cnt(ihole,1:2,1:njoined_cnt))
      np_cnt = sum(iocc_cnt(ipart,1:2,1:njoined_cnt))
      nx_cnt = sum(iocc_cnt(iextr,1:2,1:njoined_cnt))
      nv_cnt = sum(iocc_cnt(ivale,1:2,1:njoined_cnt))

      if (ntest.ge.100) then
        write(lulog,*) 'op 1:'
        call wrt_occ_n(lulog,iocc_op1,njoined_op(1))
        write(lulog,*) 'op 2:'
        call wrt_occ_n(lulog,iocc_op2,njoined_op(2))
        write(lulog,*) 'externals 1:'
        call wrt_occ_n(lulog,iocc_ex1,njoined_op(1))
        write(lulog,*) 'externals 2:'
        call wrt_occ_n(lulog,iocc_ex2,njoined_op(2))
        write(lulog,*) 'contraction:'
        call wrt_occ_n(lulog,iocc_cnt,njoined_cnt)
        write(lulog,*) 'intermediate/result:'
        call wrt_occ_n(lulog,iocc_op1op2,njoined_op1op2)
      end if

      if (ntest.ge.50) then
        write(lulog,'(x,a,i2,a,i2,a)')
     &         'Contraction scales as  H^{',nh_op1op2+nh_cnt,
     &                                   '}P^{',np_op1op2+np_cnt,'}'
        write(lulog,'(x,a,i2,a,i2,a)')
     &             'Intermediate scales as H^{',nh_op1op2,
     &                                   '}P^{',np_op1op2,'}'
      end if

      ! set iscale
      ! maximum is taken over total scaling, so H^2 .gt. P^1 
      iscale_new = 0
      iscale_new(IHOLE) = nh_op1op2+nh_cnt
      iscale_new(IPART) = np_op1op2+np_cnt
      iscale_new(IEXTR) = nx_op1op2+nx_cnt
      iscale_new(IVALE) = nv_op1op2+nv_cnt

      if (scale_rank(iscale_new).gt.scale_rank(iscale(1,1)))
     &     iscale(1:ngastp,1) = iscale_new(1:ngastp)

      iscale_new = 0
      iscale_new(IHOLE) = nh_op1op2
      iscale_new(IPART) = np_op1op2
      iscale_new(IEXTR) = nx_op1op2
      iscale_new(IVALE) = nv_op1op2

      if (scale_rank(iscale_new).gt.scale_rank(iscale(1,2)))
     &     iscale(1:ngastp,2) = iscale_new(1:ngastp)

c      possible = .true.     
      ! check whether intermediate can be addressed by
      ! the available graphs (preliminary fix)
      possible = possible.and.
     &     check_grph4occ(iocc_op1op2,irst_op1op2,njoined_op1op2,
     &     str_info,orb_info)
      possible = possible.and.
     &     check_grph4occ(iocc_op1,irst_op1,njoined_op(1),
     &     str_info,orb_info)
      if (njoined_op(2).gt.0) possible = possible.and.
     &     check_grph4occ(iocc_op2,irst_op2,njoined_op(2),
     &     str_info,orb_info)

      ! if not: do not allow this factorization
      if (.not.possible)
     &     cost(1:3) = huge(cost(1))

      if (possible) then
        if (.not.formal) then
          call init_cnt_info(cnt_info,
     &         iocc_op1,iocc_ex1,njoined_op(1),
     &            iocc_op2,iocc_ex2,njoined_op(2),
     &         iocc_cnt,njoined_cnt,
     &         iocc_op1op2,njoined_op1op2,iocc_op1op2,njoined_op1op2)

          call condense_bc_info(
     &         cnt_info,idxop(2).eq.0,
     &         iocc_op1, iocc_op2, iocc_op1op2, iocc_op1op2,
     &         iocc_ex1,iocc_ex2,iocc_cnt,
     &         irst_op1, irst_op2, irst_op1op2, irst_op1op2,
     &         idum,idum,idum,0,
     &         merge_op1, merge_op2, merge_op1op2, merge_op2op1,
     &         njoined_op(1), njoined_op(2),njoined_op1op2, njoined_cnt,
     &         str_info,orb_info)
          if (exact) then
            call dummy_contr2(flops,xmemtot,xmemblk,
     &         cnt_info,
     &         mst_op(1),mst_op(2),mst_op1op2,
     &         igamt_op(1),igamt_op(2),igamt_op1op2,
     &         str_info,ngas,nsym)          
          else
            call fact_cost_estimate(flops,xmemtot,xmemblk,
     &           cnt_info,
     &           str_info,ngas,nsym)
          end if
          call dealloc_cnt_info(cnt_info)
          cost(1) = cost(1)+flops
          cost(2) = max(cost(2),xmemtot)
          cost(3) = max(cost(3),xmemblk)
        else
c          cost(1) = cost(1) + scale_rank(iscale(1,1))
          cost(1) = max(cost(1),scale_rank(iscale(1,1)))
          cost(2) = max(cost(2),scale_rank(iscale(1,2)))
          cost(3) = max(cost(3),scale_rank(iscale(1,2)))
        end if

        if (ntest.ge.10)
     &       write(lulog,'(x,a,3(g20.10,x))') '$$ cost: ',cost(1:3)
        if (ntest.ge.10)
     &       write(lulog,'(x,a,2(2i4,2x))') '$$ scal: ',iscale(1:2,1:2)
      end if

      return
      end

