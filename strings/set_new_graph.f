*----------------------------------------------------------------------*
      subroutine set_new_graph(igraph,str_info,orb_info)
*----------------------------------------------------------------------*
*     set up graph #igraph
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'def_orbinf.h'
      include 'ifc_memman.h'
      include 'par_globalmarks.h'

      integer, intent(in) ::
     &     igraph
      type(strinf), intent(inout) ::
     &     str_info
      type(orbinf), intent(inout) ::
     &     orb_info

      integer ::
     &     mem, ihpv, nocc, nspc, ngas, nsym, ifree, ipass, idum
      integer ::
     &     leny(4), lenwscr(3)
      character(32) ::
     &     label

      integer, pointer ::
     &     iwscr(:)

c dbg
      print *,'set_new_graph: ',igraph
c dbg

      call mem_pushmark()
      
      ifree = mem_gotomark(graph_def)

      mem = 0

      ihpv = str_info%ispc_typ(igraph)
      nocc = str_info%ispc_occ(igraph)
      nspc = orb_info%ngas_hpv(ihpv)
      ! Allocate space for arc and vertex weights respectively.
      allocate(str_info%g(igraph)%yssg(nocc*nspc),
     &         str_info%g(igraph)%wssg((nocc+1)*nspc))
      mem = mem+nocc*nspc+(nocc+1)*nspc
      write(label,'("str_ssg-",i2)'),igraph
      ifree = mem_register(mem,label)

      ngas = orb_info%ngas
      nsym = orb_info%nsym
      ! set up the graphs
      ipass = 1
      call set_graph2(ipass,igraph,
     &     str_info,leny,idum,lenwscr,
     &     orb_info%ngas_hpv,orb_info%nactt_hpv,
     &     orb_info%igamorb,orb_info%mostnd,orb_info%idx_gas,ngas,nsym)

      mem = 0
      allocate(str_info%g(igraph)%y4sg(leny(1)),
     &         str_info%g(igraph)%yinf(leny(2)),
     &         str_info%g(igraph)%idis_m(str_info%g(igraph)%ndis),
     &         str_info%g(igraph)%lenstr_gm(nsym,leny(4)),
     &         str_info%g(igraph)%lenstr_dgm(leny(3)),
     &         str_info%g(igraph)%ioffstr_dgm(leny(3)) )
      str_info%max_idxms = max(str_info%max_idxms,leny(4))
      mem = mem+leny(1)+leny(2)+str_info%g(igraph)%ndis
     &       +nsym*leny(4)+2*leny(3)
      write(label,'("str_graphs-",i2)'),igraph
      ifree = mem_register(mem,label)

      allocate(iwscr(lenwscr(1)+lenwscr(2)+lenwscr(3)))

      ipass = 2
      call set_graph2(ipass,igraph,
     &     str_info,leny,iwscr,lenwscr,
     &     orb_info%ngas_hpv,orb_info%nactt_hpv,
     &     orb_info%igamorb,orb_info%mostnd,orb_info%idx_gas,ngas,nsym)

      deallocate(iwscr)      

      call mem_popmark()

      return
      end
