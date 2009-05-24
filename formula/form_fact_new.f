*----------------------------------------------------------------------*
      subroutine form_fact_new(fl_fact,contr,
     &     op_info,str_info,orb_info,iscale_stat,cost_stat,mem_stat)
*----------------------------------------------------------------------*
*     find optimum factorization of a given contraction
*     and write a sequence of unary/binary operations to fl_fact
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'def_contraction.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'def_orbinf.h'
      include 'def_reorder_info.h'
      include 'mdef_operator_info.h'
      include 'def_formula_item.h'
      
      integer, parameter ::
     &     maxcount = 10000, ! at most 10000 iterations
     &     ndisconn = 3,     ! at most 3 extra levels for disconnected
     &     ntest = 000                                   ! vertices

      type(contraction), intent(inout) ::
     &     contr
      type(formula_item), intent(inout) ::
     &     fl_fact

      type(operator_info), intent(in) ::
     &     op_info
      type(strinf), intent(in) ::
     &     str_info
      type(orbinf), intent(in), target ::
     &     orb_info
      integer, intent(inout) ::
     &     iscale_stat(ngastp,2)
      real(8), intent(inout) ::
     &     cost_stat, mem_stat

      integer ::
     &     ngas, nsym, maxidx, iarc, ivtx, icount, njoined, ncost_eval,
     &     idx, narc_full, nvtx_full, nlevel, nlevel_best,
     &     command, target, idx_intm
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
     &     reo_dummy

      logical ::
     &     possible, found
      real(8) ::
     &     cpu0, sys0, wall0, cpu, sys, wall

      logical, external ::
     &     next_fact
      integer, external ::
     &     int_pack, ifac
     
      if (ntest.ge.100) then
        call write_title(luout,wst_dbg_subr,'form_fact_new at work!')
        call prt_contr2(luout,contr,op_info)
      end if
c dbg
c      call check_xarcs(contr,op_info)
c dbg

      call atim_csw(cpu0,sys0,wall0)

      if (orb_info%ngas.eq.0)
     &     call quit(1,'form_fact_new2','buggy ngas !')
      if (orb_info%nsym.eq.0)
     &     call quit(1,'form_fact_new2','buggy nsym !')

      ngas = orb_info%ngas
      nsym = orb_info%nsym
      ihpvgas => orb_info%ihpvgas

      narc_full = contr%narc
      nvtx_full = contr%nvtx

      ! if only 1 vertex is present, we need not bother too much
      if (contr%nsupvtx.eq.1.and.narc_full.eq.0) then
          command = command_add_intm
          target = contr%idx_res
          call new_formula_item(fl_fact,command,target)
          call store_add_intm(fl_fact,contr,op_info,orb_info)
        return
      end if

      op_res => op_info%op_arr(contr%idx_res)%op
      njoined = op_res%njoined

      allocate(ifact(ld_inffac,narc_full+ndisconn),
     &     ifact_best(ld_inffac,narc_full+ndisconn),
     &     occ_vtx(ngastp,2,nvtx_full+njoined),
     &     irestr_vtx(2,orb_info%ngas,2,2,nvtx_full+njoined),
     &     info_vtx(2,nvtx_full+njoined))

      call occvtx4contr(0,occ_vtx,contr,op_info)

      call vtxinf4contr(irestr_vtx,info_vtx,contr,op_info,ngas)

      ! add 0-contractions, if necessary
      if (nvtx_full.ne.njoined) call check_disconnected(contr)
      
      found = .false.
      costmin = huge(costmin)
      nlevel = 1
      icount = 0
      ncost_eval = 0

      idx_intm = 0
c dbg
c      print *,'now diving into the recursions, njoined = ',njoined
c dbg
      call form_fact_rec_new('FIND',nlevel,ifact,fl_fact,
     &     cost,iscale,
     &     contr,occ_vtx,irestr_vtx,info_vtx)

      if (iprlvl.ge.10) then
        write(luout,*) '# of dummy contractions: ',ncost_eval
      end if

      if (.not.found) then
        call prt_contr3(luout,contr,-1)
        call prt_contr2(luout,contr,op_info)
        call quit(1,'form_fact_new','Did not find any factorization!')
      end if

      if (ntest.ge.10) then
        write(luout,*) 'optimal factorization: '
        write(luout,'(x,a,g15.10)') 'flops: ',costmin(1)
        write(luout,'(x,a,2g15.10)') 'mem:   ',costmin(2:3)
        write(luout,'(x,a,"H^",i2," P^",i2,"V^",i2,"X^",i2)')
     &       'contraction scaling:  ',iscalemin(1:4,1)
        write(luout,'(x,a,"H^",i2," P^",i2,"V^",i2,"X^",i2)')
     &       'intermediate scaling: ',iscalemin(1:4,2)
      end if

      ! statistics
      iscale_stat(1:4,1:2) = iscalemin(1:4,1:2)
      cost_stat = costmin(1)
      mem_stat  = costmin(2)

      idx_intm = 0

c dbg
      print *,'now SETting'
        call prt_contr2(luout,contr,op_info)
        call wrt_occ_n(luout,occ_vtx,nvtx_full+njoined)
c dbg
      ! call kernel again for optimal sequence --> set fl_fact now
      call form_fact_rec_new('SET',nlevel,ifact_best,fl_fact,
     &     cost,iscale,
     &     contr,occ_vtx,irestr_vtx,info_vtx)

      deallocate(ifact,ifact_best,
     &     occ_vtx,irestr_vtx,info_vtx)

      call atim_csw(cpu,sys,wall)
      if (iprlvl.ge.5)
     &    call prtim(luout,'time in form_fact_new',
     &     cpu-cpu0,sys-sys0,wall-wall0)

      if (ntest.ge.100) then
        write(luout,*) 'generated formula'
        call print_form_list(luout,fl_fact,op_info)
      end if

      return

      contains
*----------------------------------------------------------------------*
*     the recursive kernel:
*----------------------------------------------------------------------*
      recursive subroutine form_fact_rec_new(mode,
     &     nlevel,ifact,fl_in,
     &     cost_in,iscale_in,
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
      real(8), intent(in) ::
     &     cost_in(3)

      logical ::
     &     possible, new
      integer ::
     &     narc, ivtx_new, iarc, idx_op_new, ilist, len_list,
     &     iscale(ngastp,2)
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
     &     call quit(1,'form_fact_new2_r1','buggy nsym !')
      if (ntest.ge.1000) then
        call write_title(luout,wst_dbg_subr,
     &       'form_fact_rec_new rursively at work')
        write(luout,*) 'nlevel = ',nlevel
        write(luout,*) 'current cost: ',cost_in(1:3)
        if (nlevel.gt.1) then
          write(luout,*) 'current ifact: '
          write(luout,'(x,i5,"*",i5,"->",i5,"(",i5,")")')
     &         ifact(1:4,1:nlevel-1)
        end if
        write(luout,*) 'result: njoined = ',njoined
        call wrt_occ_n(luout,occ_vtx,njoined)
        write(luout,*) 'current (reduced) contraction:'
        call prt_contr3(luout,contr,occ_vtx(1,1,njoined+1))
      end if

      fl_pnt => fl_in

      call init_contr(contr_red)

      narc = contr%narc

      if (nlevel.gt.narc_full+ndisconn)
     &     call quit(1,'form_fact_rec_new','emergency exit')

      if (orb_info%nsym.eq.0) call quit(1,'form_fact2_rc','buggy nsym!')
      if (mode(1:4).eq.'FIND') then
        ! get list of (non-redundant) arcs, ordered according to
        ! contractraction strength (descending)
        call get_arc_list(arc_list,len_list,contr,orb_info)
      else if (mode(1:3).eq.'SET') then
        ! get arc form previously set ifact array
        len_list = 1
        arc_list(1) = ifact(4,nlevel)
      else
        call quit(1,'form_fact_new','unknown mode: "'//trim(MODE)//'"')
      end if

      do ilist = 1, len_list

        iarc = arc_list(ilist)

        icount = icount+1
        if (ntest.ge.1000) then
          write(luout,*) ' next arc, narc, icount: ',iarc,narc,icount
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
     &       iarc,njoined,nlevel,idx_intm,
     &       contr,occ_vtx,irestr_vtx,info_vtx,
     &       contr_red,occ_vtx_red,irestr_vtx_red,info_vtx_red,
     &       op_info,str_info,orb_info)
        ! advance pointer
        do while(fl_pnt%command.ne.command_end_of_formula)
          fl_pnt => fl_pnt%next
        end do

        if (ntest.ge.1000) then
          write(luout,'(x,a,l1)')
     &         'possible: ',possible
          if (possible) then
            write(luout,'(x,a,3f20.2)') 'cost:    ',cost
            write(luout,'(x,a,3f20.2)') 'costmin: ',costmin
            if (cost(1).ge.costmin(1))
     &           write(luout,*) 'too expensive !'
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
          call form_fact_rec_new(mode,nlevel+1,ifact,fl_pnt,
     &         cost,iscale,contr_red,occ_vtx_red,
     &                              irestr_vtx_red,info_vtx_red)

          if (ntest.ge.1000) then
            write(luout,*) 'back in level ',nlevel
            write(luout,*) ' current arc: ',iarc,narc
          end if

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
            write(luout,*) 'currently best factorization:'
            write(luout,'(x,i5,"*",i5,"->",i5,"(",i5,")")')
     &           ifact_best(1:4,1:nlevel_best)
            write(luout,*) 'current costmin:',costmin
          end if

        end if

      end do

      if (ntest.ge.1000) then
        write(luout,*) 'returning from level ',nlevel
      end if

      if (new) call dealloc_contr(contr_red)

      return

      end subroutine

      end
