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
     &     ngraph, ica, hpvx, idxmap, idx,
     &     iocc, nsym, ifree, lenoff
      integer, pointer ::
     &     idx_flipmap(:)
      integer, external ::
     &     ifndmax

      if (ntest.ge.100) then
        call write_title(luout,wst_dbg_subr,'This is strmap_man_flip')
        write(luout,*) 'igraph: ',igraph
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
          write(luout,*) 'need flip-map for: graph ',igraph(idx)
        end if

        ! does primitive map exist?
        idxmap = igraph(idx)
        if (idx_flipmap(idxmap).gt.-1) then
          ! check that offsets are set
          if (.not.associated(strmap_info%offsets_flip))
     &         call quit(1,'strmap_man_flip','offsets not initialized?')
          if (.not.associated(strmap_info%offsets_flip(idxmap)%ms)
     &           .or.
     &         .not.associated(strmap_info%offsets_flip(idxmap)%msgm))
     &         call quit(1,'strmap_man_flip','offsets not initialized?')
          if (ntest.ge.100) then
            write(luout,*) 'I have this map already ...'
          end if
          maxbuffer = maxbuffer
     &             + strmap_info%maxlen_blk_flip(idxmap)
          cycle
        end if

        if (ntest.ge.100) then
          write(luout,*) 'I will generate this map ...'
        end if

        idx_flipmap(idxmap) = strmap_info%idx_last+1
          
        ! get memory for offset arrays
        nsym = orb_info%nsym
        call mem_pushmark()
        ifree = mem_gotomark(strmaps)
        lenoff = (iocc+1)
        ifree = mem_alloc_int(strmap_info%offsets_flip(idxmap)%ms,
     &       lenoff,'m_off')
        lenoff = lenoff*nsym
        ifree = mem_alloc_int(strmap_info%offsets_flip(idxmap)%msgm,
     &                          lenoff,'mg_off')

        call mem_popmark()

        call set_flipmap(
     &         str_info%g(igraph(idx)),
     &         str_info%ispc_typ(igraph(idx)),
     &         str_info%ispc_occ(igraph(idx)),
     &         str_info%igas_restr(1,1,1,1,igraph(idx)),
C               ! ADAPT FOR OPEN SHELL  ^^^
     &         idxmap,strmap_info,orb_info)

        maxbuffer = maxbuffer + strmap_info%maxlen_blk_flip(idxmap)

      end do

      return
      end
