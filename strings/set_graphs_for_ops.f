*----------------------------------------------------------------------*
      subroutine set_graphs_for_ops(str_info,op_list,nops,orb_info)
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
      implicit none
      include 'stdunit.h'
      include 'opdim.h'
      include 'def_operator.h'
      include 'def_operator_list.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'def_orbinf.h'
      include 'ifc_memman.h'

      integer, parameter ::
     &     ntest = 100

      type(strinf), intent(inout) ::
     &     str_info
      type(operator_list), intent(inout), target ::
     &     op_list
      integer, intent(in) ::
     &     nops
      type(orbinf), intent(in) ::
     &     orb_info

      type(operator_list), pointer ::
     &     current

      integer ::
     &     maxgraph, iop, ifree, ngas,
     &     max_igtyp, mxmx_igtyp, iprint

      iprint = max(ntest,iprlvl)

      if (iprint.gt.3) then
        write(luout,*) 'Setting up graphs'
      end if

      ! we allocate more space than required, which doesn't matter
      ! too much, as the arrays allocated here are still small
      maxgraph = 0
      current => op_list
      do iop = 1, nops
        if (.not.associated(current%op))
     &       call quit(0,'set_graphs_for_ops','buggy operator list (a)')
        maxgraph = maxgraph + ngastp*current%op%n_occ_cls
        if (iop.lt.nops.and..not.associated(current%next))
     &       call quit(0,'set_graphs_for_ops','buggy operator list (b)')
        if (iop.lt.nops) current => current%next
      end do
      if (iprint.gt.50) write(luout,*) 'maxgraph = ',maxgraph
      ngas = orb_info%ngas
      allocate(str_info%ispc_typ(maxgraph),
     &         str_info%ispc_occ(maxgraph),
     &         str_info%igas_restr(2,ngas,2,maxgraph))
      ifree = mem_register(maxgraph*(2+4*ngas),'strinfo 1')

      ! identify unique graphs
      str_info%ngraph = 0
      current => op_list
      mxmx_igtyp = 0
      do iop = 1, nops
        if (ntest.ge.100) write(luout,*) 'iop = ',iop
        ! allocate graph index
        allocate(current%op%idx_graph(ngastp,2,current%op%n_occ_cls))
        ifree = mem_register(maxgraph*(2+4*ngas),
     &       trim(current%op%name)//' idxg')
c dbg
        print *,'before unique_graph'
        call flush(6)
c dbg
        call unique_graph(str_info,max_igtyp,current%op,
     &                    orb_info%ihpvgas,orb_info%ngas)
c dbg
        print *,'after unique_graph'
        call flush(6)
c dbg
        mxmx_igtyp = max(mxmx_igtyp,max_igtyp)
        if (iop.lt.nops) current => current%next
      end do

      if (mxmx_igtyp.le.0)
     &     call quit(1,'set_graphs_for_ops',
     &     'strange value of mxmx_igtyp')
      ! paramter ld_gtab is defined in def_strinf.h
      allocate(str_info%gtab(ld_gtab,mxmx_igtyp))
      ifree = mem_register(ld_gtab*mxmx_igtyp,'strinfo #')

      ! set hash table for finding matching graphs
      call set_hash4gtyp(str_info,mxmx_igtyp)

c      allocate(cc_strinfo%g(cc_strinfo%ngraph))
c      allocate(leny(4,cc_strinfo%ngraph),lenwscr(3))
c
c      do igraph = 1, cc_strinfo%ngraph
c        ihpv = cc_strinfo%ispc_typ(igraph)
c        nexc = cc_strinfo%ispc_occ(igraph)
c        nspc = cc_orbinf%ngas_hpv(ihpv)
c        allocate(cc_strinfo%g(igraph)%yssg(nexc*nspc),
c     &           cc_strinfo%g(igraph)%wssg((nexc+1)*nspc))
c      end do
c
c      ! set up the graphs
c      ipass = 1
c      call set_graph(ipass,iprint,cc_strinfo,leny,idum,lenwscr,
c     &     cc_orbinf%ngas_hpv,cc_orbinf%nactt_hpv,
c     &     cc_orbinf%igamorb,cc_orbinf%mostnd,ngas,nsmob)
c
c      do igraph = 1, cc_strinfo%ngraph
c        allocate(cc_strinfo%g(igraph)%y4sg(leny(1,igraph)),
c     &           cc_strinfo%g(igraph)%yinf(leny(2,igraph)),
c     &           cc_strinfo%g(igraph)%idis_m(cc_strinfo%g(igraph)%ndis),
c     &           cc_strinfo%g(igraph)%lenstr_gm(nsmob,leny(4,igraph)),
c     &           cc_strinfo%g(igraph)%lenstr_dgm(leny(3,igraph)),
c     &           cc_strinfo%g(igraph)%ioffstr_dgm(leny(3,igraph)) )
c      end do
c      allocate(iwscr(lenwscr(1)+lenwscr(2)+lenwscr(3)))
c
c      ipass = 2
c      call set_graph(ipass,iprint,cc_strinfo,leny,iwscr,lenwscr,
c     &     cc_orbinf%ngas_hpv,cc_orbinf%nactt_hpv,
c     &     cc_orbinf%igamorb,cc_orbinf%mostnd,ngas,nsmob)
c
c      deallocate(leny,lenwscr,iwscr)      
c
      call quit(0,'test','test exit')
      end
