*----------------------------------------------------------------------*
      subroutine form_fact2(contr,
     &     op_info,str_info,orb_info,iscale_stat,cost_stat,mem_stat)
*----------------------------------------------------------------------*
*     find optimum factorization of a given contraction
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
      
      integer, parameter ::
     &     maxcount = 10000, ! at most 10000 iterations
     &     ndisconn = 3,     ! at most 3 extra levels for disconnected
     &     ntest = 00                                   ! vertices

      type(contraction), intent(inout) ::
     &     contr

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
     &     idx, narc_full, nvtx_full, nlevel, nlevel_best
      integer, pointer ::
     &     ihpvgas(:,:),
     &     ifact(:,:), ifact_best(:,:), occ_vtx(:,:,:),
     &     irestr_vtx(:,:,:,:,:), info_vtx(:,:),
     &     iarc_ori(:), ivtx_ori(:),
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
        call write_title(lulog,wst_dbg_subr,'form_fact at work!')
        call prt_contr2(lulog,contr,op_info)
      end if
c dbg
c      call check_xarcs(contr,op_info)
c dbg

      call atim_csw(cpu0,sys0,wall0)

      if (orb_info%ngas.eq.0) call quit(1,'form_fact2','buggy ngas!')
      if (orb_info%nsym.eq.0) call quit(1,'form_fact2','buggy nsym!')

      ngas = orb_info%ngas
      nsym = orb_info%nsym
      ihpvgas => orb_info%ihpvgas

      narc_full = contr%narc
      nvtx_full = contr%nvtx

      ! if only 1 vertex is present, we need not bother too much
      if (contr%nsupvtx.eq.1) then
        ! save factorization info
        if (contr%mxfac.gt.0) deallocate(contr%inffac)
        contr%nfac = narc_full
        contr%mxfac = narc_full
        if (narc_full.eq.1) then
          allocate(contr%inffac(ld_inffac,narc_full))
          contr%inffac(1,1) = 1
          contr%inffac(2,1) = 2
          ! naming convention, see function int_pack
          contr%inffac(3,1) = int_pack(contr%inffac,2,3)
          contr%inffac(4,1) = 1
          contr%inffac(5,1) = 1
        end if
        return
      end if

      op_res => op_info%op_arr(contr%idx_res)%op
      njoined = op_res%njoined

      allocate(ifact(ld_inffac,narc_full+ndisconn),
     &     ifact_best(ld_inffac,narc_full+ndisconn),
     &     occ_vtx(ngastp,2,nvtx_full+njoined),
     &     irestr_vtx(2,orb_info%ngas,2,2,nvtx_full+njoined),
     &     info_vtx(2,nvtx_full+njoined),
     &     iarc_ori(narc_full+ndisconn),ivtx_ori(nvtx_full),
     &     irestr_res(2,orb_info%ngas,2,2,njoined))

      do iarc = 1, narc_full+ndisconn
        iarc_ori(iarc) = iarc
      end do
      
      do ivtx = 1, nvtx_full
        ivtx_ori(ivtx) = ivtx
      end do

      call occvtx4contr(0,occ_vtx,contr,op_info)

      call vtxinf4contr(irestr_vtx,info_vtx,contr,op_info,ngas)

      if (njoined.eq.1) then
        call set_restr_prel(irestr_res,contr,op_info,
     &       ihpvgas,ngas)
      else
        call dummy_restr(irestr_res,occ_vtx,njoined,orb_info)
      end if

      ! add 0-contractions, if necessary
      if (nvtx_full.ne.njoined) call check_disconnected(contr)
      
      found = .false.
      costmin = huge(costmin)
      nlevel = 1
      icount = 0
      ncost_eval = 0
c dbg
c      print *,'now diving into the recursions, njoined = ',njoined
c dbg
      call form_fact_rec(nlevel,ifact,
     &     cost,iscale,
     &     contr,occ_vtx,irestr_vtx,info_vtx,
     &     ivtx_ori,iarc_ori)

      if (iprlvl.ge.10) then
        write(lulog,*) '# of dummy contractions: ',ncost_eval
      end if

      if (.not.found) then
        call prt_contr2(lulog,contr,op_info)
        call quit(1,'form_fact2','Did not find any factorization!')
      end if

      call resize_contr(contr,contr%nvtx,contr%narc,0,nlevel_best)
      contr%nfac = nlevel_best
      contr%inffac(1:ld_inffac,1:nlevel_best) =
     &     ifact_best(1:ld_inffac,1:nlevel_best)

      if (ntest.ge.10) then
        write(lulog,*) 'optimal factorization: '
        do idx = 1, contr%nfac
          write(lulog,*) contr%inffac(1:ld_inffac,idx)
        end do
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

      deallocate(ifact,ifact_best,
     &     occ_vtx,irestr_vtx,info_vtx,iarc_ori,ivtx_ori,irestr_res)

      call atim_csw(cpu,sys,wall)
      if (iprlvl.ge.5)
     &    call prtim(lulog,'time in form_fact',
     &     cpu-cpu0,sys-sys0,wall-wall0)

      return

      contains
*----------------------------------------------------------------------*
*     the recursive kernel:
*----------------------------------------------------------------------*
      recursive subroutine form_fact_rec(nlevel,ifact,
     &     cost_in,iscale_in,
     &     contr,occ_vtx,irestr_vtx,info_vtx,
     &     ivtx_ori,iarc_ori)
*----------------------------------------------------------------------*

      implicit none

      integer, intent(inout) ::
     &     ifact(ld_inffac,*)
      type(contraction), intent(in) ::
     &     contr
      integer, intent(in) ::
     &     iscale_in(ngastp,2), nlevel,
     &     ivtx_ori(contr%nvtx), iarc_ori(contr%nvtx),
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
     &     ivtx_ori_red(contr%nvtx),
     &     iarc_ori_red(contr%narc), arc_list(contr%narc)

      integer, external ::
     &     joint_idx
     
      if (orb_info%nsym.eq.0) call quit(1,'form_fact_rec','buggy nsym!')
      if (ntest.ge.1000) then
        call write_title(lulog,wst_dbg_subr,
     &       'form_fact_rec(ursively) at work')
        write(lulog,*) 'nlevel = ',nlevel
        write(lulog,*) 'current cost: ',cost_in(1:3)
        if (nlevel.gt.1) then
          write(lulog,*) 'current ifact: '
          write(lulog,'(x,i5,"*",i5,"->",i5,"(",i5,")")')
     &         ifact(1:4,1:nlevel-1)
        end if
        write(lulog,*) 'ivtx_ori: ',ivtx_ori(1:contr%nvtx)
        write(lulog,*) 'result: njoined = ',njoined
        call wrt_occ_n(lulog,occ_vtx,njoined)
        write(lulog,*) 'current (reduced) contraction:'
        call prt_contr3(lulog,contr,occ_vtx(1,1,njoined+1))
      end if

      call init_contr(contr_red)

      narc = contr%narc

      if (nlevel.gt.narc_full+ndisconn)
     &     call quit(1,'form_fact_rec','emergency exit')

      ! get list of (non-redundant) arcs, ordered according to
      ! contractraction strength (descending)
      if (orb_info%nsym.eq.0) call quit(1,'form_fact_rec','buggy nsym!')
      call get_arc_list(arc_list,len_list,contr,orb_info)

      do ilist = 1, len_list

        iarc = arc_list(ilist)

        icount = icount+1
        if (ntest.ge.1000) then
          write(lulog,*) ' next arc, narc, icount: ',iarc,narc,icount
        end if

        ! reset
        cost = cost_in
        iscale = iscale_in

        if (icount.gt.maxcount)
     &       call quit(1,'form_fact_rec','too many cycles')

        new = .true.
        ! evaluate cost of present binary contraction
        call fact_cost2(possible,cost,iscale,
     &       iarc,njoined,nlevel,
     &       contr,occ_vtx,irestr_vtx,info_vtx,
     &       new,
     &       contr_red,occ_vtx_red,irestr_vtx_red,info_vtx_red,
     &       op_info,str_info,orb_info)

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
        if (.not.new) then
        ifact(1,nlevel) = ivtx_ori(contr%arc(iarc)%link(1))
        ifact(2,nlevel) = ivtx_ori(contr%arc(iarc)%link(2))
        ivtx_new = joint_idx(ifact(1,nlevel),ifact(2,nlevel),
     &                              nvtx_full+1,.true.)
        ifact(3,nlevel) = ivtx_new
        ifact(4,nlevel) = iarc_ori(iarc)
        ifact(5,nlevel) = iarc
        else
        ifact(1,nlevel) = contr%vertex(contr%arc(iarc)%link(1))%idx_op
        ifact(2,nlevel) = contr%vertex(contr%arc(iarc)%link(2))%idx_op
        ifact(3,nlevel) = -nlevel
        ifact(4,nlevel) = iarc !iarc_ori(iarc)
        ifact(5,nlevel) = iarc
        end if

c        if (.not.new) then
c        ! ... and find out what the
c        ! contr structure looks like after this operation
c        call init_contr(contr_red)
c        call copy_contr(contr,contr_red)
c        occ_vtx_red = occ_vtx
c        ivtx_ori_red = ivtx_ori
c        iarc_ori_red = iarc_ori
c
c        idx_op_new = op_info%nops+nlevel
cc        call add_interm_info(irestr_interm,info_interm,
cc     &       idx_op_new,irestr_res,contr,occ_vtx)
c
c        irestr_vtx_red = irestr_vtx
c        info_vtx_red = info_vtx
c        
c        call reduce_contr(contr_red,occ_vtx_red,
c     &       possible,
c     &       iarc,idx_op_new,ivtx_new,
c     &       njoined,
c     &       .true.,ivtx_ori_red,iarc_ori_red,
c     &       .true.,irestr_vtx_red,info_vtx_red,irestr_res,
c     &       .false.,reo_dummy,orb_info)
c
c        if (.not.possible) cycle
c
c        end if ! old

        ! add 0-contractions, if necessary
        ! I think it is needed here as we else miss disconnected
        ! terms ....
        call check_disconnected(contr_red)
c dbg
c        print *,'calling check disc for'
c        call prt_contr3(lulog,contr_red,occ_vtx_red(1,1,njoined+1))
c dbg

        ! any contraction left?
        if (contr_red%narc.gt.0) then

          ! add 0-contractions, if necessary
c          ! why did we put it here?
c          call check_disconnected(contr_red)
          
          call form_fact_rec(nlevel+1,ifact,
     &         cost,iscale,contr_red,occ_vtx_red,
     &                              irestr_vtx_red,info_vtx_red,
     &         ivtx_ori_red,iarc_ori_red)

          if (ntest.ge.1000) then
            write(lulog,*) 'back in level ',nlevel
            write(lulog,*) ' current arc: ',iarc,narc
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
