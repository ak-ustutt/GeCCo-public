*----------------------------------------------------------------------*
      subroutine set_graphs_for_ops(str_info,op_list,nops,orb_info)
*----------------------------------------------------------------------*
*     the possible index strings for a given H/P/V and C/A
*     can be described by graphs which are set up in this driver
*     routine
*
*     andreas, april 2007
*
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
      include 'par_globalmarks.h'

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
     &     maxgraph, iop, ifree, ngas, nsym, ipass, idum, mem,
     &     max_igtyp, mxmx_igtyp, iprint, ihpv, nocc, nspc, igraph

      integer, allocatable ::
     &     leny(:,:), lenwscr(:), iwscr(:)

      iprint = max(ntest,iprlvl)

      if (iprint.gt.3) then
        call write_title(luout,wst_section,'Setting up graphs')
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
      ! Allocate space for information on subspace type, occupancy, 
      ! and restrictions.
      allocate(str_info%ispc_typ(maxgraph),
     &         str_info%ispc_occ(maxgraph),
     &         str_info%igas_restr(2,ngas,2,maxgraph))
      ifree = mem_register(maxgraph*(2+4*ngas),'strinfo_base')

      ! identify unique graphs by checking occupancies of each space.
      str_info%ngraph = 0
      current => op_list
      mxmx_igtyp = 0

      do iop = 1, nops
        if(.not.current%op%formal)then
          if (ntest.ge.100) write(luout,*) 'iop = ',iop
          call unique_graph(str_info,max_igtyp,current%op,
     &         orb_info%ihpvgas,orb_info%ngas)
          mxmx_igtyp = max(mxmx_igtyp,max_igtyp)
        endif  
        if (iop.lt.nops) current => current%next
      end do

      if (mxmx_igtyp.le.0)
     &     call quit(1,'set_graphs_for_ops',
     &     'strange value of mxmx_igtyp')
      ! paramter ld_gtab is defined in def_strinf.h
      allocate(str_info%gtab(ld_gtab,mxmx_igtyp))
      ifree = mem_register(ld_gtab*mxmx_igtyp,'strinfo_hash')

      ! set hash table for finding matching graphs
      call set_hash4gtyp(str_info,mxmx_igtyp)

      ! -----------------------------
      ! and now to the actual graphs:
      ! -----------------------------
      allocate(str_info%g(str_info%ngraph))
      allocate(leny(4,str_info%ngraph),lenwscr(3))

      mem = 0
      do igraph = 1, str_info%ngraph
        ihpv = str_info%ispc_typ(igraph)
        nocc = str_info%ispc_occ(igraph)
        nspc = orb_info%ngas_hpv(ihpv)
        ! Allocate space for arc and vertex weights respectively.
        allocate(str_info%g(igraph)%yssg(nocc*nspc),
     &           str_info%g(igraph)%wssg((nocc+1)*nspc))
        mem = mem+nocc*nspc+(nocc+1)*nspc
      end do
      ifree = mem_register(mem,'strinfo_ssg')

      ngas = orb_info%ngas
      nsym = orb_info%nsym
      ! set up the graphs
      ipass = 1
      call set_graph(ipass,
     &     str_info,leny,idum,lenwscr,
     &     orb_info%ngas_hpv,orb_info%nactt_hpv,
     &     orb_info%igamorb,orb_info%mostnd,orb_info%idx_gas,ngas,nsym)

      mem = 0
      do igraph = 1, str_info%ngraph
        allocate(str_info%g(igraph)%y4sg(leny(1,igraph)),
     &           str_info%g(igraph)%yinf(leny(2,igraph)),
     &           str_info%g(igraph)%idis_m(str_info%g(igraph)%ndis),
     &           str_info%g(igraph)%lenstr_gm(nsym,leny(4,igraph)),
     &           str_info%g(igraph)%lenstr_dgm(leny(3,igraph)),
     &           str_info%g(igraph)%ioffstr_dgm(leny(3,igraph)) )
        mem = mem+leny(1,igraph)+leny(2,igraph)+str_info%g(igraph)%ndis
     &       +nsym*leny(4,igraph)+2*leny(3,igraph)
      end do
      ifree = mem_register(mem,'strinfo_graphs')

      allocate(iwscr(lenwscr(1)+lenwscr(2)+lenwscr(3)))

      ipass = 2
      call set_graph(ipass,
     &     str_info,leny,iwscr,lenwscr,
     &     orb_info%ngas_hpv,orb_info%nactt_hpv,
     &     orb_info%igamorb,orb_info%mostnd,orb_info%idx_gas,ngas,nsym)

      deallocate(leny,lenwscr,iwscr)      

      end
