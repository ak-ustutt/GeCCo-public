*----------------------------------------------------------------------*
      subroutine form_fact_new(fl_fact,contr,
     &     op_info,str_info,orb_info,iscale_stat,cost_stat,mem_stat,
     &     iitem)
*----------------------------------------------------------------------*
*     find optimum factorization of a given contraction
*     and write a sequence of unary/binary operations to fl_fact
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'routes.h'
      include 'def_contraction.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'def_orbinf.h'
      include 'def_reorder_info.h'
      include 'mdef_operator_info.h'
      include 'def_formula_item.h'
      
      integer, parameter ::
     &     maxcount = 1000000, ! at most 1000000 iterations
     &     ndisconn = 14,     ! at most 4 extra levels for disconnected
     &     ntest = 000                                   ! vertices

      type(contraction), intent(inout) ::
     &     contr
      type(formula_item), intent(inout) ::
     &     fl_fact

      type(operator_info), intent(inout) ::
     &     op_info
      type(strinf), intent(inout) ::
     &     str_info
      type(orbinf), intent(in), target ::
     &     orb_info
      integer, intent(inout) ::
     &     iscale_stat(ngastp,2), iitem
      real(8), intent(inout) ::
     &     cost_stat, mem_stat

      integer ::
     &     ngas, nsym, maxidx, iarc, ivtx, icount, njoined, ncost_eval,
     &     idx, narc_full, nvtx_full, nlevel, nlevel_best,
     &     command, target, idx_intm, idum, nvtx_max, nreo_op1op2
      integer, pointer ::
     &     ihpvgas(:,:),
     &     ifact(:,:), ifact_best(:,:), occ_vtx(:,:,:),
     &     irestr_vtx(:,:,:,:,:), info_vtx(:,:),
     &     irestr_res(:,:,:,:,:)

*----------------------------------------------------------------------*
*     ifact(4,*): describes factorization by giving
*            ( vertex-number 1, vertex-number 2, 
*                         result-vertex-number, arc-number )
*          all numbers might refer to fused vertices/arcs
*          the number is then an ordered list of all vertices/arcs fused
*          packed using the number of plain arcs (narc) or vertices (nvtx)
*          in the current contraction (as given in contr)
*----------------------------------------------------------------------*
      real(8) ::
     &     cost(3), costmin(3)
      integer ::
     &     iscale(ngastp,2), iscalemin(ngastp,2)
      type(operator), pointer ::
     &     op_res
      type(reorder_info) ::
     &     reo_info0

      logical ::
     &     possible, found, predef, reo_add
      real(8) ::
     &     cpu0, sys0, wall0, cpu, sys, wall

      integer(8), pointer ::
     &     vtx(:), topo(:,:), xlines(:,:), op_res_p(:), xlines_dum(:,:)
      integer, pointer ::
     &     svertex(:),
     &     iocc_ori(:,:,:), iocc_reo(:,:,:),
     &     irst_ori(:,:,:,:,:), irst_reo(:,:,:,:,:),
     &     merge_stp1(:), merge_stp1inv(:),
     &     merge_stp2(:), merge_stp2inv(:)

      logical, external ::
     &     next_fact
      integer, external ::
     &     int_pack, ifac
     
      if (ntest.ge.100) then
        call write_title(lulog,wst_dbg_subr,'form_fact_new at work!')
        call prt_contr2(lulog,contr,op_info)
      end if
c dbg
c      call check_xarcs(contr,op_info)
c dbg

      call atim_csw(cpu0,sys0,wall0)

      if (orb_info%ngas.eq.0)
     &     call quit(1,'form_fact_new','buggy ngas !')
      if (orb_info%nsym.eq.0)
     &     call quit(1,'form_fact_new','buggy nsym !')

      ngas = orb_info%ngas
      nsym = orb_info%nsym
      ihpvgas => orb_info%ihpvgas

      narc_full = contr%narc
      nvtx_full = contr%nvtx

      op_res => op_info%op_arr(contr%idx_res)%op
      njoined = op_res%njoined

      ! if only 1 vertex is present, we need not bother too much
      if (contr%nsupvtx.eq.1.and.narc_full.eq.0) then

        reo_add = .false.
        if (nvtx_full.gt.1) then
          ! we have to check whether reordering is necessary
          call init_reo_info(reo_info0)
          nvtx_max = max(nvtx_full,njoined)
          allocate(svertex(nvtx_full),vtx(nvtx_full),op_res_p(njoined),
     &             topo(nvtx_full,nvtx_full),xlines(nvtx_max,njoined),
     &             xlines_dum(nvtx_max,nvtx_max))
          xlines = 0
          xlines_dum = 0
          call pack_contr(svertex,vtx,topo,
     &                    xlines(1:nvtx_full,1:njoined),contr,njoined)
c dbg
c          write(lulog,*) 'checking for necessary reordering:'
c          call prt_contr_p(lulog,svertex,vtx,topo,
c     &            xlines(1:nvtx_full,1:njoined),nvtx_full,njoined)
c dbgend
          op_res_p = 0
          do idx = 1, njoined
            do ivtx = 1, nvtx_full
              op_res_p(idx) = op_res_p(idx) + xlines(ivtx,idx)
            end do
          end do

          xlines_dum = 0
          xlines_dum(1:nvtx_full,1:njoined)
     &          = xlines(1:nvtx_full,1:njoined)
          call set_final_reo(reo_info0,xlines,xlines_dum,op_res_p,
     &                     contr%vertex(1)%idx_op,nvtx_full,njoined)

          reo_add = (reo_info0%nreo.gt.0)
          deallocate(svertex,vtx,topo,xlines,op_res_p,xlines_dum)
          if (.not.reo_add) call dealloc_reo_info(reo_info0)
        end if

        if (reo_add) then
          allocate(iocc_reo(ngastp,2,nvtx_max),
     &             iocc_ori(ngastp,2,nvtx_max),
     &             irst_reo(2,ngas,2,2,nvtx_max),
     &             irst_ori(2,ngas,2,2,nvtx_max),
     &             merge_stp1(2*nvtx_max*nvtx_max),
     &             merge_stp1inv(2*nvtx_max*nvtx_max),
     &             merge_stp2(2*nvtx_max*nvtx_max),
     &             merge_stp2inv(2*nvtx_max*nvtx_max),
     &             occ_vtx(ngastp,2,nvtx_max+njoined),
     &             irestr_vtx(2,ngas,2,2,nvtx_max+njoined),
     &             info_vtx(2,nvtx_max+njoined),
     &             svertex(nvtx_max))
          iocc_reo = 0
          iocc_ori = 0
          irst_reo = 0
          irst_ori = 0
          occ_vtx = 0
          irestr_vtx = 0
          svertex(1:nvtx_max) = 1

          call occvtx4contr(0,occ_vtx,contr,op_info)

          call vtxinf4contr(irestr_vtx,info_vtx,contr,op_info,ngas)

          call get_reo_info2(-1,1,
     &           iocc_reo,iocc_ori,
     &           irst_reo,irst_ori,
     &           ivtx,idx,idum, ! <-- dummies
     &           merge_stp1,merge_stp1inv,merge_stp2,merge_stp2inv,
     &           occ_vtx,irestr_vtx,
     &                   svertex,info_vtx,
     &                       njoined,nvtx_max,
     &           reo_info0,nreo_op1op2,str_info,orb_info)

          command = command_add_reo
          target = contr%idx_res
          call new_formula_item(fl_fact,command,target)
          call store_add_intm(fl_fact,contr,op_info,orb_info)
          call store_reorder(fl_fact,
     &       op_res%name,op_info%op_arr(contr%vertex(1)%idx_op)%op%name,
     &       contr%iblk_res,(contr%vertex(1)%iblk_op-1)/nvtx_full+1,
     &       contr%dagger,contr%vertex(1)%dagger,
     &       reo_info0%sign_reo,reo_info0%iocc_opreo0,
     &       reo_info0%from_to,reo_info0%iocc_reo,nreo_op1op2,
     &       reo_info0%nreo-nreo_op1op2,
     &       iocc_reo,irst_reo,nvtx_max,
     &       iocc_ori,irst_ori,nvtx_max,
     &       merge_stp1,merge_stp1inv,merge_stp2,merge_stp2inv,
     &       orb_info)
          if (lustat.gt.0)
     &       call print_form_item(lustat,iitem,fl_fact,op_info)
          deallocate(iocc_reo,iocc_ori,irst_reo,irst_ori,
     &             merge_stp1,merge_stp1inv,merge_stp2,merge_stp2inv,
     &             occ_vtx,irestr_vtx,info_vtx,svertex)
          call dealloc_reo_info(reo_info0)
          return
        else
          command = command_add_intm
          target = contr%idx_res
          call new_formula_item(fl_fact,command,target)
          call store_add_intm(fl_fact,contr,op_info,orb_info)
          if (lustat.gt.0)
     &       call print_form_item(lustat,iitem,fl_fact,op_info)
          return
        end if
      else
        allocate(occ_vtx(ngastp,2,nvtx_full+njoined),
     &       irestr_vtx(2,orb_info%ngas,2,2,nvtx_full+njoined),
     &       info_vtx(2,nvtx_full+njoined))
      end if

      allocate(ifact(ld_inffac,narc_full+ndisconn),
     &     ifact_best(ld_inffac,narc_full+ndisconn))

      call occvtx4contr(0,occ_vtx,contr,op_info)

      call vtxinf4contr(irestr_vtx,info_vtx,contr,op_info,ngas)

      ! add 0-contractions, if necessary
cmh      if (nvtx_full.ne.njoined) call check_disconnected(contr)
      call check_disconnected(contr)
      
      found = .false.
      costmin = huge(costmin)
      nlevel = 1
      icount = 0
      ncost_eval = 0
      iscale = 0

      idx_intm = 0

      ! allow predefined contraction sequence:
      predef = contr%nfac.gt.0
      if (predef) ifact(4,1:contr%nfac) = contr%inffac(4,1:contr%nfac)

c dbg
c      print *,'now diving into the recursions, njoined = ',njoined
c dbg
      call form_fact_rec_new('FIND',predef,nlevel,ifact,fl_fact,
     &     cost,iscale,iitem,
     &     contr,occ_vtx,irestr_vtx,info_vtx)

      if (iprlvl.ge.10) then
        write(lulog,*) '# of dummy contractions: ',ncost_eval
      end if

      if (.not.found) then
        call prt_contr3(lulog,contr,-1)
        call prt_contr2(lulog,contr,op_info)
        call quit(1,'form_fact_new','Did not find any factorization!')
      end if

      if (ntest.ge.10) then
        write(lulog,*) 'optimal factorization: '
        write(lulog,'(x,a,g15.10)') 'flops: ',costmin(1)
        write(lulog,'(x,a,2g15.10)') 'mem:   ',costmin(2:3)
        write(lulog,'(x,a,"H^",i2," P^",i2,"V^",i2,"X^",i2)')
     &       'contraction scaling:  ',iscalemin(1:4,1)
        write(lulog,'(x,a,"H^",i2," P^",i2,"V^",i2,"X^",i2)')
     &       'intermediate scaling: ',iscalemin(1:4,2)
      end if

      ! statistics
      iscale_stat(1:4,1:2) = iscalemin(1:4,1:2)
      cost_stat = costmin(1)
      mem_stat  = costmin(2)

      idx_intm = 0

c dbg
c      print *,'now SETting'
c        call prt_contr2(lulog,contr,op_info)
c        call wrt_occ_n(lulog,occ_vtx,nvtx_full+njoined)
c dbg
      ! call kernel again for optimal sequence --> set fl_fact now
      call form_fact_rec_new('SET',.true.,nlevel,ifact_best,fl_fact,
     &     cost,iscale,iitem,
     &     contr,occ_vtx,irestr_vtx,info_vtx)

      deallocate(ifact,ifact_best,
     &     occ_vtx,irestr_vtx,info_vtx)

      call atim_csw(cpu,sys,wall)
      if (iprlvl.ge.5)
     &    call prtim(lulog,'time in form_fact_new',
     &     cpu-cpu0,sys-sys0,wall-wall0)

      if (ntest.ge.100) then
        write(lulog,*) 'generated formula'
        call print_form_list(lulog,fl_fact,op_info)
      end if

      return

      contains
*----------------------------------------------------------------------*
*     the recursive kernel:
*----------------------------------------------------------------------*
      recursive subroutine form_fact_rec_new(mode,predef,
     &     nlevel,ifact,fl_in,
     &     cost_in,iscale_in,iitem,
     &     contr,occ_vtx,irestr_vtx,info_vtx)
*----------------------------------------------------------------------*

      implicit none

      type(formula_item), intent(inout), target ::
     &     fl_in
      character(len=*), intent(in) ::
     &     mode
      integer, intent(inout) ::
     &     ifact(ld_inffac,*)
      type(contraction), intent(in) ::
     &     contr
      integer, intent(in) ::
     &     iscale_in(ngastp,2), nlevel,
     &     occ_vtx(ngastp,2,contr%nvtx+njoined),
     &     irestr_vtx(2,ngas,2,2,contr%nvtx+njoined),
     &     info_vtx(2,contr%nvtx+njoined)
      integer, intent(inout) ::
     &     iitem
      real(8), intent(in) ::
     &     cost_in(3)
      logical, intent(in) ::
     &     predef

      logical ::
     &     possible, new
      integer ::
     &     narc, ivtx_new, iarc, idx_op_new, ilist, len_list,
     &     iscale(ngastp,2), nlist
      real(8) ::
     &     cost(3)
      type(contraction) ::
     &     contr_red
      integer ::
     &     occ_vtx_red(ngastp,2,contr%nvtx+njoined),
     &     irestr_vtx_red(2,ngas,2,2,contr%nvtx+njoined),
     &     info_vtx_red(2,contr%nvtx+njoined),
     &     arc_list(contr%narc)
      type(formula_item), pointer ::
     &     fl_pnt

      integer, external ::
     &     joint_idx
     
      if (orb_info%nsym.eq.0)
     &     call quit(1,'form_fact_rec_new','buggy nsym !')
      if (ntest.ge.1000) then
        call write_title(lulog,wst_dbg_subr,
     &       'form_fact_rec_new rursively at work')
        write(lulog,*) 'nlevel = ',nlevel
        write(lulog,*) 'current cost: ',cost_in(1:3)
        if (nlevel.gt.1) then
          write(lulog,*) 'current ifact: '
          write(lulog,'(x,i5,"*",i5,"->",i5,"(",i5,")")')
     &         ifact(1:4,1:nlevel-1)
        end if
        write(lulog,*) 'result: njoined = ',njoined
        call wrt_occ_n(lulog,occ_vtx,njoined)
        write(lulog,*) 'current (reduced) contraction:'
        call prt_contr3(lulog,contr,occ_vtx(1,1,njoined+1))
      end if

      fl_pnt => fl_in

      call init_contr(contr_red)

      narc = contr%narc

      if (nlevel.gt.narc_full+ndisconn)
     &     call quit(1,'form_fact_rec_new','emergency exit')

      if (orb_info%nsym.eq.0) call quit(1,'form_fact_rec_new'
     &                                   ,'buggy nsym!')
      if (mode(1:4).eq.'FIND'.and..not.predef) then
        ! get list of (non-redundant) arcs, ordered according to
        ! contractraction strength (descending)
        call get_arc_list(arc_list,len_list,contr,orb_info)
      else if (mode(1:3).eq.'SET'
     &         .or.mode(1:4).eq.'FIND'.and.predef) then
        ! get arc form previously set ifact array
        len_list = 1
        arc_list(1) = ifact(4,nlevel)
      else
        call quit(1,'form_fact_rec_new'
     &             ,'unknown mode: "'//trim(MODE)//'"')
      end if

      ! only try first maxbranch elements of arc_list (if maxbranch is set)
      nlist = len_list
      if (maxbranch.gt.0) nlist = min(len_list,maxbranch)

      do ilist = 1, nlist

        iarc = arc_list(ilist)

        icount = icount+1
        if (ntest.ge.1000) then
          write(lulog,*) ' next arc, narc, icount: ',iarc,narc,icount
        end if

        ! reset
        cost = cost_in
        iscale = iscale_in

        if (icount.gt.maxcount)
     &       call quit(1,'form_fact_rec_new','too many cycles')

        new = .true.
        ! process binary contraction
        ! 'FIND' -- get cost of contraction
        ! 'SET'  -- set fl_fact
        call process_bc(mode,fl_pnt,possible,cost,iscale,
     &       iarc,njoined,nlevel,idx_intm,iitem,
     &       contr,occ_vtx,irestr_vtx,info_vtx,
     &       contr_red,occ_vtx_red,irestr_vtx_red,info_vtx_red,
     &       op_info,str_info,orb_info)
        ! advance pointer
        do while(fl_pnt%command.ne.command_end_of_formula)
          fl_pnt => fl_pnt%next
        end do

        if (ntest.ge.1000) then
          write(lulog,'(x,a,l1)')
     &         'possible: ',possible
          if (possible) then
            write(lulog,'(x,a,3f20.2)') 'cost:    ',cost
            write(lulog,'(x,a,3f20.2)') 'costmin: ',costmin
            if (cost(1).ge.costmin(1))
     &           write(lulog,*) 'too expensive !'
          end if
        end if

        ! the contraction is not allowed (for some reason)
        if (.not.possible) cycle

        ncost_eval = ncost_eval+1
        ! already more expensive than current best factorization?
        if (cost(1).ge.costmin(1)) cycle

        ! if the binary contraction is accepted, ...

        ! ... update ifact array ...
        ifact(1,nlevel) = contr%vertex(contr%arc(iarc)%link(1))%idx_op
        ifact(2,nlevel) = contr%vertex(contr%arc(iarc)%link(2))%idx_op
        ifact(3,nlevel) = -nlevel
        ifact(4,nlevel) = iarc
        ifact(5,nlevel) = iarc

        ! add 0-contractions, if necessary
        ! I think it is needed here as we else miss disconnected
        ! terms ....
        call check_disconnected(contr_red)

        ! any contraction left?
        if (contr_red%narc.gt.0) then

          ! add 0-contractions, if necessary          
          call form_fact_rec_new(mode,predef,nlevel+1,ifact,fl_pnt,
     &         cost,iscale,iitem,contr_red,occ_vtx_red,
     &                              irestr_vtx_red,info_vtx_red)

          if (ntest.ge.1000) then
            write(lulog,*) 'back in level ',nlevel
            write(lulog,*) ' current arc: ',iarc,narc
          end if

          ! reset intermediate counter
          idx_intm = 1 - nlevel

        else

          ! report, that we found at least one possible
          ! factorization
          found = .true.
          ! if we survived up to here, we are the current champion!
          costmin = cost
          iscalemin = iscale
          ifact_best(1:ld_inffac,1:nlevel) =
     &         ifact(1:ld_inffac,1:nlevel)
          nlevel_best = nlevel

          if (ntest.ge.1000) then
            write(lulog,*) 'currently best factorization:'
            write(lulog,'(x,i5,"*",i5,"->",i5,"(",i5,")")')
     &           ifact_best(1:4,1:nlevel_best)
            write(lulog,*) 'current costmin:',costmin
          end if

        end if

      end do

      if (ntest.ge.1000) then
        write(lulog,*) 'returning from level ',nlevel
      end if

      if (new) call dealloc_contr(contr_red)

      return

      end subroutine

      end
