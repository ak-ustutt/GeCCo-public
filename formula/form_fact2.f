*----------------------------------------------------------------------*
      subroutine form_fact2(contr,
     &     op_info,str_info,orb_info,iscale_stat)
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
      include 'mdef_operator_info.h'
      
      integer, parameter ::
     &     maxcount = 1000,  ! at most 1000 iterations
     &     ndisconn = 3,     ! at most 3 extral levels for disconnected
     &     ntest = 000                                   ! vertices

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

      integer ::
     &     ngas, nsym, maxidx, iarc, ivtx, icount,
     &     idx, narc_full, nvtx_full, nlevel, nlevel_best
c      integer, pointer ::
c     &     ihpvgas(:)
      integer, pointer ::
     &     ihpvgas(:),
     &     ifact(:,:), ifact_best(:,:), occ_vtx(:,:,:),
     &     irestr_vtx(:,:,:,:,:), info_vtx(:,:),
     &     iarc_ori(:), ivtx_ori(:)
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
     &     iscale(ngastp,2), iscalemin(ngastp,2),
     &     irestr_res(2,orb_info%ngas,2,2)

      logical ::
     &     possible, found
      real(8) ::
     &     cpu0, sys0, wall0, cpu, sys, wall

      logical, external ::
     &     next_fact
      integer, external ::
     &     int_pack, ifac
     
      if (ntest.ge.100) then
        call write_title(luout,wst_dbg_subr,'form_fact at work!')
        call prt_contr2(luout,contr,op_info)
      end if

      call atim_csw(cpu0,sys0,wall0)

      ngas = orb_info%ngas
      nsym = orb_info%nsym
      ihpvgas => orb_info%ihpvgas

      narc_full = contr%narc
      nvtx_full = contr%nvtx

      ! if no or only 1 arc is present, we need not bother too much
c      if (narc_full.le.0) then
      if (nvtx_full.eq.1) then
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

      allocate(ifact(ld_inffac,narc_full+ndisconn),
     &     ifact_best(ld_inffac,narc_full+ndisconn),
     &     occ_vtx(ngastp,2,nvtx_full+1),
     &     irestr_vtx(2,orb_info%ngas,2,2,nvtx_full+1),
     &     info_vtx(2,nvtx_full+1),
     &     iarc_ori(narc_full+ndisconn),ivtx_ori(nvtx_full))

      do iarc = 1, narc_full+ndisconn
        iarc_ori(iarc) = iarc
      end do
      
      do ivtx = 1, nvtx_full
        ivtx_ori(ivtx) = ivtx
      end do

      call set_restr_prel(irestr_res,contr,op_info,
     &     ihpvgas,ngas)

      call occvtx4contr(0,occ_vtx,contr,op_info)

      call vtxinf4contr(irestr_vtx,info_vtx,contr,op_info,ngas)

      ! add 0-contractions, if necessary
      call check_disconnected(contr)
      
      found = .false.
      costmin = huge(costmin)
      nlevel = 1
      icount = 0
c dbg
c      print *,'now diving into the recursions'
c dbg
      call form_fact_rec(nlevel,ifact,
     &     cost,iscale,
     &     contr,occ_vtx,irestr_vtx,info_vtx,
     &     ivtx_ori,iarc_ori)

      if (.not.found) then
        call prt_contr2(luout,contr,op_info)
        call quit(1,'form_fact','Did not find any factorization!')
      end if

      call resize_contr(contr,contr%nvtx,contr%narc,nlevel_best)
      contr%nfac = nlevel_best
      contr%inffac(1:ld_inffac,1:nlevel_best) =
     &     ifact_best(1:ld_inffac,1:nlevel_best)

      if (ntest.ge.100) then
        write(luout,*) 'optimal factorization: '
        do idx = 1, contr%nfac
          write(luout,*) contr%inffac(1:ld_inffac,idx)
        end do
        write(luout,'(x,a,g15.10)') 'flops: ',costmin(1)
        write(luout,'(x,a,2g15.10)') 'mem:   ',costmin(2:3)
        write(luout,'(x,a,"H^",i2," P^",i2,"V^",i2)')
     &       'contraction scaling:  ',iscalemin(1:3,1)
        write(luout,'(x,a,"H^",i2," P^",i2,"V^",i2)')
     &       'intermediate scaling: ',iscalemin(1:3,2)
      end if

      ! statistics
      if ( sum(iscalemin(1:3,1)).gt.sum(iscale_stat(1:3,1)).or.
     &    (sum(iscalemin(1:3,1)).eq.sum(iscale_stat(1:3,1)).and.
     &     iscalemin(2,1).gt.iscale_stat(2,1)) ) then
        iscale_stat(1:3,1) = iscalemin(1:3,1)
      end if

      if ( sum(iscalemin(1:3,2)).gt.sum(iscale_stat(1:3,2)).or.
     &    (sum(iscalemin(1:3,2)).eq.sum(iscale_stat(1:3,2)).and.
     &     iscalemin(2,2).gt.iscale_stat(2,2)) ) then
        iscale_stat(1:3,2) = iscalemin(1:3,2)
      end if

      deallocate(ifact,ifact_best)

      call atim_csw(cpu,sys,wall)
      if (iprlvl.ge.5)
     &    call prtim(luout,'time in form_fact',
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
     &     occ_vtx(ngastp,2,contr%nvtx+1),
     &     irestr_vtx(2,ngas,2,2,contr%nvtx+1),
     &     info_vtx(2,contr%nvtx+1)
      real(8), intent(in) ::
     &     cost_in(3)

      logical ::
     &     possible
      integer ::
     &     narc, ivtx_new, iarc, idx_op_new,
     &     iscale(ngastp,2)
      real(8) ::
     &     cost(3)
      type(contraction) ::
     &     contr_red
      integer ::
     &     occ_vtx_red(ngastp,2,contr%nvtx+1),
     &     irestr_vtx_red(2,ngas,2,2,contr%nvtx+1),
     &     info_vtx_red(2,contr%nvtx+1),
     &     ivtx_ori_red(contr%nvtx),
     &     iarc_ori_red(contr%narc)

      integer, external ::
     &     joint_idx
c dbg
     &    , idxlist
c dbg      
      
      if (ntest.ge.1000) then
        call write_title(luout,wst_dbg_subr,
     &       'form_fact_rec(ursively) at work')
        write(luout,*) 'nlevel = ',nlevel
        write(luout,*) 'current cost: ',cost_in(1:3)
        if (nlevel.gt.1) then
          write(luout,*) 'current ifact: '
          write(luout,'(x,i5,"*",i5,"->",i5,"(",i5,")")')
     &         ifact(1:4,1:nlevel-1)
        end if
        write(luout,*) 'current (reduced) contraction:'
        call prt_contr3(luout,contr,occ_vtx)
      end if

      narc = contr%narc

      if (nlevel.gt.narc_full+ndisconn)
     &     call quit(1,'form_fact_rec','emergency exit')

      do iarc = 1, narc

        icount = icount+1
        if (ntest.ge.1000) then
          write(luout,*) ' next arc, narc, icount: ',iarc,narc,icount
        end if

        ! reset
        cost = cost_in
        iscale = iscale_in

        if (icount.gt.maxcount)
     &       call quit(1,'form_fact_rec','too many cycles')

        ! evaluate cost of present binary contraction
        call fact_cost2(possible,cost,iscale,
     &       contr,occ_vtx,irestr_vtx,info_vtx,iarc,
     &       op_info,str_info,ihpvgas,ngas,nsym)

        if (ntest.ge.1000) then
          write(luout,'(x,a,l1)')
     &         'possible: ',possible
          if (possible) then
            write(luout,'(x,a,3f20.2)') 'cost:    ',cost
            write(luout,'(x,a,3f20.2)') 'costmin: ',costmin
          end if
        end if

        ! the contraction is not allowed (for some reason)
        if (.not.possible) cycle

        ! already more expensive than current best factorization?
        if (cost(1).ge.costmin(1)) cycle

        ! if the binary contraction is accepted, ...

        ! ... update ifact array ...
        ifact(1,nlevel) = ivtx_ori(contr%arc(iarc)%link(1))
        ifact(2,nlevel) = ivtx_ori(contr%arc(iarc)%link(2))
        ivtx_new = joint_idx(ifact(1,nlevel),ifact(2,nlevel),
     &                              nvtx_full+1,.true.)
        ifact(3,nlevel) = ivtx_new
        ifact(4,nlevel) = iarc_ori(iarc)
        ifact(5,nlevel) = iarc

c dbg
c        if (nlevel.gt.1.and.
c     &      (ifact(1,nlevel).gt.nvtx_full.and.
c     &       idxlist(ifact(1,nlevel),ifact,4*(nlevel-1),1).le.0 .or.
c     &       ifact(2,nlevel).gt.nvtx_full.and.
c     &       idxlist(ifact(2,nlevel),ifact,4*(nlevel-1),1).le.0) ) then
c          write(luout,*) 'something is wrong:'
c          write(luout,'(x,i5,"*",i5,"->",i5,"(",i5,")")')
c     &           ifact(1:4,1:nlevel)
c          write(luout,*) 'ivtx_ori = ',ivtx_ori(contr%nvtx)
c          write(luout,*) 'iarc_ori = ',iarc_ori(contr%narc)
c          stop 'WRONG'
c        end if
c dbg        

        ! ... and find out what the
        ! contr structure looks like after this operation
        call init_contr(contr_red)
        call copy_contr(contr,contr_red)
        occ_vtx_red = occ_vtx
        ivtx_ori_red = ivtx_ori
        iarc_ori_red = iarc_ori

        idx_op_new = op_info%nops+nlevel
c        call add_interm_info(irestr_interm,info_interm,
c     &       idx_op_new,irestr_res,contr,occ_vtx)
        
        call reduce_contr(contr_red,occ_vtx_red,
     &       iarc,idx_op_new,ivtx_new,
     &       .true.,ivtx_ori_red,iarc_ori_red)

        ! add 0-contractions, if necessary
        call check_disconnected(contr_red)
      
        ! any contraction left?
        if (contr_red%narc.gt.0) then
          
          irestr_vtx_red = irestr_vtx
          info_vtx_red = info_vtx

          call reduce_vtx_info(irestr_vtx_red,info_vtx_red,
     &                         contr,occ_vtx,iarc, ! pass the old contr!
     &                         irestr_res,orb_info)

          call form_fact_rec(nlevel+1,ifact,
     &         cost,iscale,contr_red,occ_vtx_red,
     &                              irestr_vtx_red,info_vtx_red,
     &         ivtx_ori_red,iarc_ori_red)

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

      end subroutine

      end
