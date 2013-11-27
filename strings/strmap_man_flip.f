*----------------------------------------------------------------------*
      subroutine strmap_man_flip(
     &     maxbuffer,
     &     igraph,nblk,
     &     str_info,strmap_info,orb_info)
*----------------------------------------------------------------------*
*     manage string maps:
*     generate new flip map, if necessary
*     maxbuffer returns the maximum number of integer words needed 
*     for an IRREP block
*----------------------------------------------------------------------*
      implicit none

      integer, parameter ::
     &     ntest = 00

      include 'opdim.h'
      include 'stdunit.h'
      include 'def_graph.h'
      include 'def_filinf.h'
      include 'def_strinf.h'
      include 'def_strmapinf.h'
      include 'def_orbinf.h'
      include 'ifc_memman.h'
      include 'par_globalmarks.h'

      integer, intent(out) ::
     &     maxbuffer
      integer, intent(in) ::
     &     nblk,
     &     igraph(nblk)
      type(strinf), intent(in) ::
     &     str_info
      type(strmapinf), intent(inout) ::
     &     strmap_info
      type(orbinf), intent(in) ::
     &     orb_info

      integer ::
     &     ngraph, ica, hpvx, idxgraph, idx,
     &     iocc, nsym, ifree, lenoff
      integer, pointer ::
     &     idx_flipmap(:)
      integer, external ::
     &     ifndmax

      if (ntest.ge.100) then
        call write_title(lulog,wst_dbg_subr,'This is strmap_man_flip')
        write(lulog,*) 'igraph: ',igraph
      end if

      maxbuffer = 0

      ngraph = strmap_info%mxgraph
      if (ngraph.lt.str_info%ngraph)
     &     call quit(1,'strmap_man_c',
     &     'you forgot to update the maps after adding a new graph')
      idx_flipmap => strmap_info%idx_flipmap

      if (strmap_info%ffstrmap%unit.le.0)
     &     call quit(1,'strmap_man',
     &     'ffstrmap must be open when calling strmap_man')

      do idx = 1, nblk
        iocc = str_info%ispc_occ(igraph(idx))

        if (ntest.ge.100) then
          write(lulog,*) 'need flip-map for: graph ',igraph(idx)
        end if

        ! does primitive map exist?
        idxgraph = igraph(idx)
        if (idx_flipmap(idxgraph).gt.-1) then
          ! check that offsets are set
          if (.not.associated(strmap_info%offsets_flip))
     &         call quit(1,'strmap_man_flip','offsets not initialized?')
          if (.not.associated(strmap_info%offsets_flip(idxgraph)%ms)
     &           .or.
     &         .not.associated(strmap_info%offsets_flip(idxgraph)%msgm))
     &         call quit(1,'strmap_man_flip','offsets not initialized?')
          if (ntest.ge.100) then
            write(lulog,*) 'I have this map already ...'
          end if
          maxbuffer = maxbuffer
     &             + strmap_info%maxlen_blk_flip(idxgraph)
          cycle
        end if

        if (ntest.ge.100) then
          write(lulog,*) 'I will generate this map ...'
        end if

        idx_flipmap(idxgraph) = strmap_info%idx_last+1
          
        ! get memory for offset arrays
        nsym = orb_info%nsym
        call mem_pushmark()
        ifree = mem_gotomark(strmaps)
        lenoff = (iocc+1)
        ifree = mem_alloc_int(strmap_info%offsets_flip(idxgraph)%ms,
     &       lenoff,'m_off')
        lenoff = lenoff*nsym
        ifree = mem_alloc_int(strmap_info%offsets_flip(idxgraph)%msgm,
     &                          lenoff,'mg_off')

        call mem_popmark()

        call set_flipmap(
     &         str_info%g(idxgraph),
     &         str_info%ispc_typ(idxgraph),
     &         str_info%ispc_occ(idxgraph),
     &         str_info%igas_restr(1,1,1,1,idxgraph),
C               ! ADAPT FOR OPEN SHELL  ^^^
     &         idxgraph,strmap_info,orb_info)

        maxbuffer = maxbuffer + strmap_info%maxlen_blk_flip(idxgraph)

      end do

      if (ntest.ge.100) then
        write(lulog,*) 'leaving strmap_man_flip() ...'
      end if

      return
      end
