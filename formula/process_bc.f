*----------------------------------------------------------------------*
      subroutine process_bc(mode,fl_fact,possible,cost,iscale,
     &     iarc,njoined_res,nlevel,idx_intm,iitem,
     &     contr,occ_vtx,irestr_vtx,info_vtx,
     &     contr_red,occ_vtx_red,irestr_vtx_red,info_vtx_red,
     &     op_info,str_info,orb_info)
*----------------------------------------------------------------------*
*
*     extract the next binary contraction and set up 
*     reduced contraction info
*
*     for the binary contraction:
*
*     mode == 'FIND'
*
*     calculate the computational cost of a given contraction
*     for each contraction step, the number of flops is summed up
*     and the maximum size of the occurring intermediates is evaluated
*     returned on cost(1:3): cost(1)<-flops, 
*                            cost(2)<-memory(tot)
*                            cost(3)<-memory(block)
*     the maximum scaling of contraction and intermediates is returned
*     on iscale(ngastp,1:2): iscale(,1)<-contr, iscale(,2)<-mem
*
*     mode == 'SET'
*
*     put info on contraction
*
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
      include 'def_formula_item.h'
      include 'multd2h.h'
      include 'routes.h'
      
      integer, parameter ::
     &     ntest_ = 00
      integer :: ntest
      integer :: ntest_bc
      logical, parameter ::
     &     formal = .false., exact = .false.

      character(len=*), intent(in) ::
     &     mode
      integer, intent(in) ::
     &     iarc, njoined_res, nlevel
      integer, intent(inout) ::
     &     idx_intm, iitem
      type(contraction), intent(in) ::
     &     contr
      type(contraction), intent(out) ::
     &     contr_red
      type(formula_item), intent(inout), target ::
     &     fl_fact
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
      type(strinf), intent(inout) ::
     &     str_info
      real(8), intent(out) ::
     &     cost(3)
      integer, intent(out) ::
     &     iscale(ngastp,2)
      logical, intent(out) ::
     &     possible

      character(len=len_opname) ::
     &     label, label1, label2, label_reo
      real(8) ::
     &     fact, fact_itf
      integer ::
     &     nvtx, narc, ngas, nsym, idum, target, command, nidx,
     &     np_op1op2, nh_op1op2, nx_op1op2, np_cnt, nh_cnt, nx_cnt,
     &     idxop_reo, idxop_ori, iblkop_reo, iblkop_ori,
     &     ireo, idxs_reo, mode_rst_cnt, nv_op1op2, nv_cnt, nreo_op1op2,
     &     iscal
      integer ::
     &     iocc_cnt(ngastp,2,contr%nvtx*(contr%nvtx+1)/2),
     &     iocc_ex1(ngastp,2,contr%nvtx),
     &     iocc_ex2(ngastp,2,contr%nvtx),
     &     irst_res(2,orb_info%ngas,2,2,contr%nvtx),
     &     iocc_op1(ngastp,2,contr%nvtx),
     &     iocc_op2(ngastp,2,contr%nvtx),
     &     iocc_ori(ngastp,2,contr%nvtx),
     &     iocc_reo(ngastp,2,contr%nvtx),
     &     iocc_op1op2(ngastp,2,contr%nvtx),
     &     irst_op1op2(2,orb_info%ngas,2,2,contr%nvtx),
     &     iocc_op1op2tmp(ngastp,2,contr%nvtx),
     &     irst_op1op2tmp(2,orb_info%ngas,2,2,contr%nvtx),
     &     irst_ex1(2,orb_info%ngas,2,2,contr%nvtx),
     &     irst_ex2(2,orb_info%ngas,2,2,contr%nvtx),
     &     irst_cnt(2,orb_info%ngas,2,2,contr%nvtx*(contr%nvtx+1)/2),
     &     igr_op1op2(ngastp,2,contr%nvtx),
     &     irst_ori(2,orb_info%ngas,2,2,contr%nvtx),
     &     irst_reo(2,orb_info%ngas,2,2,contr%nvtx), 
     &     irst_op1(2,orb_info%ngas,2,2,contr%nvtx),
     &     irst_op2(2,orb_info%ngas,2,2,contr%nvtx), 
     &     mst_op(2), mst_op1op2, mst_opreo,
     &     igamt_op(2), igamt_op1op2, igamt_opreo,
     &     njoined_op(2), njoined_op1op2, njoined_cnt, nj_ret,
     &     idxop(2), iblkop(2), iblkop1op2, iblkop1op2tmp,
     &     merge_op1(2*contr%nvtx*contr%nvtx), ! a bit too large, I guess ...
     &     merge_op2(2*contr%nvtx*contr%nvtx),
     &     merge_op1op2(2*contr%nvtx*contr%nvtx),
     &     merge_op2op1(2*contr%nvtx*contr%nvtx),
     &     merge_stp1(2*contr%nvtx*contr%nvtx),
     &     merge_stp1inv(2*contr%nvtx*contr%nvtx),
     &     merge_stp2(2*contr%nvtx*contr%nvtx),
     &     merge_stp2inv(2*contr%nvtx*contr%nvtx),
     &     merge_stp1_0(2*contr%nvtx*contr%nvtx),
     &     merge_stp1inv_0(2*contr%nvtx*contr%nvtx),
     &     merge_stp2_0(2*contr%nvtx*contr%nvtx),
     &     merge_stp2inv_0(2*contr%nvtx*contr%nvtx),
     &     iscale_new(ngastp)
      integer, pointer ::
     &     itf_index_info(:)
      
      integer, pointer ::
     &     ihpvgas(:,:)
      
      type(formula_item), pointer ::
     &     fl_pnt
      type(operator_array), pointer ::
     &     op_arr(:)

      type(contraction) ::
     &     contr_dum
      type(reorder_info) ::
     &     reo_info0, reo_info
      type(contraction_info) ::
     &     cnt_info

      logical ::
     &     tra_op1, tra_op2, tra_op1op2, tra_reo, tra_ori,
     &     reo_op1op2, reo_other, reo_before
      real(8) ::
     &     flops, xmemtot, xmemblk, bc_sign, bc_sign_itf, factor

      integer, external ::
     &     idxlist, int_expand, int_pack, maxxlvl_op,  get_nidx4contr
      logical, external ::
     &     check_grph4occ
      real(8), external ::
     &     scale_rank

      if (mode(1:3)=='SET') then
        ntest_bc = ntest_
        ntest = ntest_
      else
        ntest = 0
        ntest_bc = 0
      end if
      
      if (ntest.gt.0) then
        call write_title(lulog,wst_dbg_subr,'this is process_bc')
        write(lulog,*) 'mode = ',trim(mode)
      end if

      op_arr => op_info%op_arr

      fl_pnt => fl_fact

      nsym = orb_info%nsym
      ngas = orb_info%ngas
      ihpvgas => orb_info%ihpvgas

      nvtx = contr%nvtx
      narc = contr%narc

      ! reset reo_info
      call init_reo_info(reo_info0) ! FIX
      call init_reo_info(reo_info)

      nidx = get_nidx4contr(contr)
      allocate(itf_index_info(3+2*nidx))
      itf_index_info = 0

      ! extract BC
      call get_bc_info3(bc_sign,bc_sign_itf,possible,
     &     idxop,iblkop,
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
     &     .true.,
     &     contr_red,occ_vtx_red,irestr_vtx_red,info_vtx_red,
     &     .true.,reo_info,reo_info0, !FIX
     &     iarc,.true.,idx_intm,
     &     irst_res,njoined_res,orb_info,op_info,ntest_bc) ! irst_res is dummy

      reo_before = .false.
      reo_op1op2 = .false.
      reo_other  = .false.
      do ireo = 1, reo_info0%nreo
        reo_before = reo_before.or.
     &             reo_info0%reo(ireo)%reo_before
        reo_op1op2 = reo_op1op2.or.
     &       (     reo_info0%reo(ireo)%is_bc_result.and.
     &        .not.reo_info0%reo(ireo)%reo_before)
        reo_other  = reo_other.or.
     &       (.not.reo_info0%reo(ireo)%is_bc_result.and.
     &        .not.reo_info0%reo(ireo)%reo_before)
      end do
      if (reo_op1op2.or.reo_other)
     &     call quit(1,'process_bc','trap1')
      do ireo = 1, reo_info%nreo
        if (reo_info%reo(ireo)%reo_before)
     &     call quit(1,'process_bc','trap2')
        reo_op1op2 = reo_op1op2.or.
     &       (     reo_info%reo(ireo)%is_bc_result.and.
     &        .not.reo_info%reo(ireo)%reo_before)
        reo_other  = reo_other.or.
     &       (.not.reo_info%reo(ireo)%is_bc_result.and.
     &        .not.reo_info%reo(ireo)%reo_before)
      end do

      iocc_op1op2tmp = iocc_op1op2
      irst_op1op2tmp = irst_op1op2
      
      if (reo_op1op2.and.possible) then
        do ireo = 1, reo_info%nreo
          if (reo_info%reo(ireo)%is_bc_result.and.
     &         .not.reo_info%reo(ireo)%reo_before) then
            idxs_reo = reo_info%reo(ireo)%idxsuper
            exit
          end if
        end do
        call get_reo_info2(1,idxs_reo,
     &           iocc_op1op2,iocc_op1op2tmp,
     &           irst_op1op2,irst_op1op2tmp,
     &           njoined_op1op2,mst_op1op2,igamt_op1op2,
     &           merge_stp1,merge_stp1inv,merge_stp2,merge_stp2inv,
     &           occ_vtx_red,irestr_vtx_red,
     &                       contr_red%svertex,info_vtx_red,
     &                       njoined_res,contr_red%nvtx,
     &           reo_info,nreo_op1op2,str_info,orb_info)
      end if

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
        if (reo_op1op2) then
          write(lulog,*) 'intermediate/result bef. reo:'
          call wrt_occ_n(lulog,iocc_op1op2tmp,njoined_op1op2)
        end if
        write(lulog,*) 'intermediate/result:'
        call wrt_occ_n(lulog,iocc_op1op2,njoined_op1op2)
      end if

cmh: allow even if no graph exists --> will be created later (get_reo_info)
c      ! check whether intermediate can be addressed by
c      ! the available graphs (preliminary fix)
c      possible = possible.and.
c     &     check_grph4occ(iocc_op1op2,irst_op1op2,njoined_op1op2,
c     &     str_info,orb_info)
c      possible = possible.and.
c     &     check_grph4occ(iocc_op1,irst_op1,njoined_op(1),
c     &     str_info,orb_info)
c      if (njoined_op(2).gt.0) possible = possible.and.
c     &     check_grph4occ(iocc_op2,irst_op2,njoined_op(2),
c     &     str_info,orb_info)

      ! process restrictions:
      call ex1ex2cnt_restr(
     &         irst_ex1, irst_ex2, irst_cnt,
     &         -1,
     &         idxop(2).eq.0,
     &         iocc_op1, iocc_op2, iocc_op1op2, iocc_op1op2tmp,
     &         iocc_ex1,iocc_ex2,iocc_cnt,
     &         irst_op1, irst_op2, irst_op1op2, irst_op1op2tmp,
     &         merge_op1, merge_op2, merge_op1op2, merge_op2op1,
     &         njoined_op(1), njoined_op(2),njoined_op1op2, njoined_cnt,
     &         str_info,orb_info)

      if (mode.eq.'FIND'.or.lustat.gt.0) then
        ! count particle, hole, (active) spaces involved:
        ! in intermediate
        nh_op1op2 = sum(iocc_ex1(ihole,1:2,1:njoined_op(1)))
     &            + sum(iocc_ex2(ihole,1:2,1:njoined_op(2)))
        np_op1op2 = sum(iocc_ex1(ipart,1:2,1:njoined_op(1)))
     &            + sum(iocc_ex2(ipart,1:2,1:njoined_op(2)))
        nx_op1op2 = sum(iocc_ex1(iextr,1:2,1:njoined_op(1)))
     &            + sum(iocc_ex2(iextr,1:2,1:njoined_op(2)))
        nv_op1op2 = sum(iocc_ex1(ivale,1:2,1:njoined_op(1)))
     &            + sum(iocc_ex2(ivale,1:2,1:njoined_op(2)))
        ! in contraction length
        nh_cnt = sum(iocc_cnt(ihole,1:2,1:njoined_cnt))
        np_cnt = sum(iocc_cnt(ipart,1:2,1:njoined_cnt))
        nx_cnt = sum(iocc_cnt(iextr,1:2,1:njoined_cnt))
        nv_cnt = sum(iocc_cnt(ivale,1:2,1:njoined_cnt))

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
      end if
      if (mode.eq.'FIND') then

        if (scale_rank(iscale_new).gt.scale_rank(iscale(1,1)))
     &       iscale(1:ngastp,1) = iscale_new(1:ngastp)

        iscale_new = 0
        iscale_new(IHOLE) = nh_op1op2
        iscale_new(IPART) = np_op1op2
        iscale_new(IEXTR) = nx_op1op2
        iscale_new(IVALE) = nv_op1op2

        if (scale_rank(iscale_new).gt.scale_rank(iscale(1,2)))
     &       iscale(1:ngastp,2) = iscale_new(1:ngastp)

        ! if not: do not allow this factorization
        if (.not.possible)
     &     cost(1:3) = huge(cost(1))

        if (possible) then
          call init_cnt_info(cnt_info,
     &         iocc_op1,iocc_ex1,njoined_op(1),
     &            iocc_op2,iocc_ex2,njoined_op(2),
     &         iocc_cnt,njoined_cnt,
     &         iocc_op1op2,njoined_op1op2,iocc_op1op2tmp,njoined_op1op2)

c          mode_rst_cnt = 1 ! set and return irst_ex1/ex2/cnt
          mode_rst_cnt = 2 ! use as set before
          call condense_bc_info(
     &         cnt_info,idxop(2).eq.0,
     &         iocc_op1, iocc_op2, iocc_op1op2, iocc_op1op2tmp,
     &         iocc_ex1,iocc_ex2,iocc_cnt,
     &         irst_op1, irst_op2, irst_op1op2, irst_op1op2tmp,
     &         irst_ex1, irst_ex2, irst_cnt, mode_rst_cnt,
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
            call diag_type2cnt_info(cnt_info,idxop,contr_red,op_info)
            call fact_cost_estimate(flops,xmemtot,xmemblk,
     &           cnt_info,
     &           str_info,ngas,nsym)
          end if
          call dealloc_cnt_info(cnt_info)
c          if (flops.eq.0d0) call quit(1,'process_bc',
c     &        'operator with zero length?')
          ! reorderings are expensive:
          factor = 1d0
          if (reo_before) factor = factor+reo_factor
          if (reo_other)  factor = factor+reo_factor
          if (reo_op1op2) factor = factor+reo_factor

          ! get a penalty for large intermediates
          !iscal = max(1,sum(iscale(1:ngastp,2))-3)
          !factor = factor*2d0**dble(iscal)
          
          cost(1) = cost(1)+flops*factor
          cost(2) = max(cost(2),xmemtot)
          cost(3) = max(cost(3),xmemblk)

          if (ntest.ge.10)
     &        write(lulog,'(x,a,3(g20.10,x))') '$$ cost: ',cost(1:3)
          if (ntest.ge.10)
     &        write(lulog,'(x,a,2(2i4,2x))') '$$ scal: ',iscale(1:2,1:2)
        end if
      else if (mode.eq.'SET') then

        ! a reordering takes place before contraction?
        if (reo_before) then
          command = command_new_intermediate
          target = -1           ! unused here
          call new_formula_item(fl_pnt,command,target)

          do ireo = 1, reo_info0%nreo
            if (reo_info0%reo(ireo)%reo_before) then
              idxs_reo = reo_info0%reo(ireo)%idxsuper
              idxop_reo = reo_info0%reo(ireo)%idxop_new
              iblkop_reo = 1
              tra_reo = reo_info0%reo(ireo)%dagger_new
              idxop_ori = reo_info0%reo(ireo)%idxop_ori
              iblkop_ori = reo_info0%reo(ireo)%iblkop_ori
              tra_ori = reo_info0%reo(ireo)%dagger_ori
              exit
            end if
          end do

          if (idxop_reo.gt.0) then
            label_reo = op_arr(idxop_reo)%op%name
          else
            write(label_reo,'("_STIN",i4.4)') abs(idxop_reo)
          end if

          if (idxop_ori.gt.0) then
            label1 = op_arr(idxop_ori)%op%name
          else
            write(label1,'("_STIN",i4.4)') abs(idxop_ori)
          end if
          
          call get_reo_info2(-1,idxs_reo,
     &           iocc_reo,iocc_ori,
     &           irst_reo,irst_ori,
     &           nj_ret,mst_opreo,igamt_opreo,
     &           merge_stp1_0,merge_stp1inv_0,
     &                        merge_stp2_0,merge_stp2inv_0,
     &           occ_vtx,irestr_vtx,
     &                   contr%svertex,info_vtx,
     &                       njoined_res,contr%nvtx,
     &           reo_info0,idum,str_info,orb_info)

          ! FIX:
          if (nj_ret.gt.1.and..not.tra_ori) then
            iblkop_ori = (iblkop_ori-1)/nj_ret + 1
          else if (nj_ret.gt.1.and.tra_ori) then
            iblkop_ori = (iblkop_ori-1)/nj_ret + 1
          end if

          call store_def_intm(fl_pnt,
     &         label_reo,iocc_reo,irst_reo,nj_ret,1,
     &         label1,' ',tra_reo,tra_ori,.false.,
     &         orb_info)
          if (lustat.gt.0) then
            call print_form_item(lustat,iitem,fl_pnt,op_info)
            write(lustat,'(x,"Formal size of intermediate: '//
     &          'H^",i2," P^",i2," V^",i2," X^",i2)')
     &          nh_op1op2,np_op1op2,nv_op1op2,nx_op1op2
          end if
          
          fl_pnt => fl_pnt%next

          command = command_reorder
          target = -1           ! unused here
          call new_formula_item(fl_pnt,command,target)

          call store_reorder(fl_pnt,
     &       label_reo,label1,
     &       iblkop_reo,iblkop_ori,
     &       tra_reo, tra_ori,
     &       reo_info0%sign_reo,reo_info0%iocc_opreo0,
     &       reo_info0%from_to,reo_info0%iocc_reo,reo_info0%nreo,0,
     &       iocc_reo,irst_reo,nj_ret,
     &       iocc_ori,irst_ori,nj_ret,
     &       merge_stp1_0,merge_stp1inv_0,merge_stp2_0,merge_stp2inv_0,
     &       orb_info)
          if (lustat.gt.0)
     &       call print_form_item(lustat,iitem,fl_pnt,op_info)

          fl_pnt => fl_pnt%next

        end if

        target = -1 ! unused here
        if (idxop(1).gt.0) then
          label1 = op_arr(idxop(1))%op%name
        else
          write(label1,'("_STIN",i4.4)') abs(idxop(1))
        end if
        if (idxop(2).gt.0) then
          label2 = op_arr(idxop(2))%op%name
        else if (idxop(2).eq.0) then
          label2 = '---'
        else
          write(label2,'("_STIN",i4.4)') abs(idxop(2))
        end if

        ! not in final level? result is new intermediate
        if (contr_red%narc.gt.0) then
          command = command_new_intermediate
          target = -1           ! unused here
          write(label,'("_STIN",i4.4)') abs(idx_intm)
          call new_formula_item(fl_pnt,command,target)
          call store_def_intm(fl_pnt,
     &         label,iocc_op1op2,irst_op1op2,njoined_op1op2,1,
     &         label1,label2,tra_op1op2,tra_op1,tra_op2,
     &         orb_info)
          if (lustat.gt.0) then
            call print_form_item(lustat,iitem,fl_pnt,op_info)
            write(lustat,'(x,"Formal size of intermediate: '//
     &          'H^",i2," P^",i2," V^",i2," X^",i2)')
     &          nh_op1op2,np_op1op2,nv_op1op2,nx_op1op2
          end if
          fl_pnt => fl_pnt%next
          fact = bc_sign
          fact_itf = 1d0
          iblkop1op2 = 1
          command = command_bc
          if (reo_op1op2) command = command_bc_reo
        else
          label = op_arr(contr%idx_res)%op%name
          fact = bc_sign*contr%fac      
          fact_itf = dble(contr%total_sign)*contr%fac*contr%eqvl_fact
          iblkop1op2 = contr%iblk_res
          command = command_add_bc
          if (reo_op1op2) command = command_add_bc_reo
        end if
          
        call new_formula_item(fl_pnt,command,target)
        call store_bc(fl_pnt,
     &       fact,fact_itf,
     &       label,label1,label2,
     &       iblkop1op2,iblkop(1),iblkop(2),
     &       tra_op1op2,tra_op1,tra_op2,
     &       njoined_op1op2,njoined_op(1),njoined_op(2),
     &       iocc_op1op2tmp,iocc_op1,iocc_op2,
     &       irst_op1op2tmp,irst_op1,irst_op2,
     &       iocc_ex1,iocc_ex2,iocc_cnt,
     &       irst_ex1,irst_ex2,irst_cnt,njoined_cnt,
     &       merge_op1,merge_op2,
     &       merge_op1op2,merge_op2op1,
     &       itf_index_info,
     &       orb_info)
        if (reo_op1op2) then
          iblkop1op2tmp = 1
          call store_reorder(fl_pnt,
     &       label,label,
     &       iblkop1op2,iblkop1op2tmp,
     &       tra_op1op2,tra_op1op2,
     &       reo_info%sign_reo,reo_info%iocc_opreo0,
     &       reo_info%from_to,reo_info%iocc_reo,nreo_op1op2,
     &       reo_info%nreo-nreo_op1op2,
     &       iocc_op1op2,irst_op1op2,njoined_op1op2,
     &       iocc_op1op2tmp,irst_op1op2tmp,njoined_op1op2,
     &       merge_stp1,merge_stp1inv,merge_stp2,merge_stp2inv,
     &       orb_info)
        end if
        if (lustat.gt.0) then
          call print_form_item(lustat,iitem,fl_pnt,op_info)
          write(lustat,'(x,"Formal scaling of contraction: '//
     &          'H^",i2," P^",i2," V^",i2," X^",i2)')
     &          iscale_new(IHOLE),iscale_new(IPART),
     &          iscale_new(IVALE),iscale_new(IEXTR)
        end if
        if (reo_other) then
          fl_pnt => fl_pnt%next
          command = command_new_intermediate
          target = -1           ! unused here
          call new_formula_item(fl_pnt,command,target)

          do ireo = 1, reo_info%nreo
            if (.not.reo_info%reo(ireo)%is_bc_result.and.
     &          .not.reo_info%reo(ireo)%reo_before) then
              idxs_reo = reo_info%reo(ireo)%idxsuper
              idxop_reo = reo_info%reo(ireo)%idxop_new
              iblkop_reo = 1
              tra_reo = reo_info%reo(ireo)%dagger_new
              idxop_ori = reo_info%reo(ireo)%idxop_ori
              iblkop_ori = reo_info%reo(ireo)%iblkop_ori
              tra_ori = reo_info%reo(ireo)%dagger_ori
              exit
            end if
          end do

          if (idxop_reo.gt.0) then
            label_reo = op_arr(idxop_reo)%op%name
          else
            write(label_reo,'("_STIN",i4.4)') abs(idxop_reo)
          end if

          if (idxop_ori.gt.0) then
            label = op_arr(idxop_ori)%op%name
          else
            write(label,'("_STIN",i4.4)') abs(idxop_ori)
          end if
          
          call get_reo_info2(1,idxs_reo,
     &           iocc_reo,iocc_ori,
     &           irst_reo,irst_ori,
     &           nj_ret,mst_opreo,igamt_opreo,
     &           merge_stp1,merge_stp1inv,merge_stp2,merge_stp2inv,
     &           occ_vtx_red,irestr_vtx_red,
     &                       contr_red%svertex,info_vtx_red,
     &                       njoined_res,contr_red%nvtx,
     &           reo_info,idum,str_info,orb_info)

          ! FIX:
          if (nj_ret.gt.1.and..not.tra_ori) then
            iblkop_ori = (iblkop_ori-1)/nj_ret + 1
          else if (nj_ret.gt.1.and.tra_ori) then
            iblkop_ori = (iblkop_ori-1)/nj_ret + 1
          end if

          call store_def_intm(fl_pnt,
     &         label_reo,iocc_reo,irst_reo,nj_ret,1,
     &         label,' ',tra_reo,tra_ori,.false.,
     &         orb_info)
          if (lustat.gt.0) then
            call print_form_item(lustat,iitem,fl_pnt,op_info)
            write(lustat,'(x,"Formal size of intermediate: '//
     &          'H^",i2," P^",i2," V^",i2," X^",i2)')
     &          nh_op1op2,np_op1op2,nv_op1op2,nx_op1op2
          end if
          
          fl_pnt => fl_pnt%next

          command = command_reorder
          target = -1           ! unused here
          call new_formula_item(fl_pnt,command,target)

          call store_reorder(fl_pnt,
     &       label_reo,label,
     &       iblkop_reo,iblkop_ori,
     &       tra_reo, tra_ori,
     &       reo_info%sign_reo,reo_info%iocc_opreo0,
     &       reo_info%from_to,reo_info%iocc_reo,reo_info%nreo,0,
     &       iocc_reo,irst_reo,nj_ret,
     &       iocc_ori,irst_ori,nj_ret,
     &       merge_stp1,merge_stp1inv,merge_stp2,merge_stp2inv,
     &       orb_info)
          if (lustat.gt.0)
     &       call print_form_item(lustat,iitem,fl_pnt,op_info)
        end if
      else
        call quit(1,'process_bc','unknown mode: "'//trim(mode)//'"')
      end if

      deallocate(itf_index_info)

      call dealloc_reo_info(reo_info)
      call dealloc_reo_info(reo_info0) ! FIX

      return
      end

