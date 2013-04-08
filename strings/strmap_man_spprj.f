*----------------------------------------------------------------------*
      subroutine strmap_man_spprj(
     &     maxbuffer,nmaps,
     &     igraph,msdis,nblk,
     &     str_info,strmap_info,orb_info)
*----------------------------------------------------------------------*
*     manage spin projection maps:
*     generate new spin projection map, if necessary
*     maxbuffer returns the maximum number of integer words needed 
*     for an IRREP block times combination cases
*----------------------------------------------------------------------*
      implicit none

      integer, parameter ::
     &     ntest = 100

      include 'opdim.h'
      include 'stdunit.h'
      include 'def_graph.h'
      include 'def_filinf.h'
      include 'def_strinf.h'
      include 'def_strmapinf.h'
      include 'def_orbinf.h'
      include 'ifc_memman.h'
      include 'par_globalmarks.h'

      integer, intent(in) ::
     &     nblk,
     &     igraph(nblk), msdis(nblk)
      integer, intent(out) ::
     &     maxbuffer, nmaps(nblk)
      type(strinf), intent(in) ::
     &     str_info
      type(strmapinf), intent(inout) ::
     &     strmap_info
      type(orbinf), intent(in) ::
     &     orb_info

      logical ::
     &     init
      integer ::
     &     ngraph, ica, hpvx, idxgraph, idx,
     &     iocc, nsym, ifree, lenoff, idxms
      integer, pointer ::
     &     idx_spprjmap(:)
      integer, external ::
     &     ifndmax

      if (ntest.ge.100) then
        call write_title(luout,wst_dbg_subr,'This is strmap_man_spprj')
        write(luout,*) 'igraph: ',igraph
        write(luout,*) 'msdis:  ',msdis
      end if

      maxbuffer = 0

      ngraph = strmap_info%mxgraph
      if (ngraph.lt.str_info%ngraph)
     &     call quit(1,'strmap_man_c',
     &     'you forgot to update the maps after adding a new graph')
      idx_spprjmap => strmap_info%idx_spprjmap

      if (strmap_info%ffstrmap%unit.le.0)
     &     call quit(1,'strmap_man',
     &     'ffstrmap must be open when calling strmap_man')

      blk_loop: do idx = 1, nblk
        iocc = str_info%ispc_occ(igraph(idx))

        ! number of elements to which element of sceleton is mapped (at most)
        nmaps(idx) = 2**iocc-1

        if (ntest.ge.100) then
          write(luout,*) 'need spin-projection map for: graph ',
     &                    igraph(idx), msdis(idx), iocc, nmaps(idx)
        end if

        idxms = (iocc-msdis(idx))/2 + 1

        ! does primitive map exist?
        idxgraph = igraph(idx)

        if (iocc.eq.1) then ! skip this trivial case
          maxbuffer = maxbuffer +
     &         ifndmax(
     &         str_info%g(idxgraph)%lenstr_gm,1,
     &         (str_info%ispc_occ(idxgraph)+1)*orb_info%nsym,1)

          cycle blk_loop ! go to next entry

        end if

        if (idx_spprjmap(idxgraph).gt.-1) then
          ! check that offsets are set
          if (.not.associated(strmap_info%offsets_spprj))
     &        call quit(1,'strmap_man_spprj','offsets not initialized?')
          if (.not.associated(strmap_info%offsets_spprj(idxgraph)%ms)
     &           .or.
     &        .not.associated(strmap_info%offsets_spprj(idxgraph)%msgm))
     &        call quit(1,'strmap_man_spprj','offsets not initialized?')
          if (strmap_info%offsets_spprj(idxgraph)%ms(idxms).gt.-1) then
            if (ntest.ge.100) then
              write(luout,*) 'I have this map already ...'
            end if
            maxbuffer = maxbuffer
     &             + strmap_info%maxlen_blk_spprj(idxgraph)

            ! go to next entry
            cycle blk_loop

          end if
        end if

        if (ntest.ge.100) then
          write(luout,*) 'I will generate this map ...'
        end if

        ! new occupation case or just new ms case?
        init = .false.
        if (idx_spprjmap(idxgraph).le.0) then
          init = .true.
          idx_spprjmap(idxgraph) = strmap_info%idx_last+1
          
          ! get memory for offset arrays
          nsym = orb_info%nsym
          call mem_pushmark()
          ifree = mem_gotomark(strmaps)
          lenoff = (iocc+1)
          ifree = mem_alloc_int(strmap_info%offsets_spprj(idxgraph)%ms,
     &       lenoff,'m_off')
          lenoff = lenoff*nsym
          ifree =mem_alloc_int(strmap_info%offsets_spprj(idxgraph)%msgm,
     &                          lenoff,'mg_off')

          call mem_popmark()

        end if

        call set_spprjmap(init,
     &         str_info%g(idxgraph),
     &         str_info%ispc_typ(idxgraph),
     &         str_info%ispc_occ(idxgraph),
     &         str_info%igas_restr(1,1,1,1,idxgraph),
C               ! ADAPT FOR OPEN SHELL  ^^^
     &         idxgraph,msdis(idx),strmap_info,orb_info)

        maxbuffer = maxbuffer + strmap_info%maxlen_blk_spprj(idxgraph)

      end do blk_loop

      if (ntest.ge.100) then
        write(luout,*) 'leaving strmap_man_spprj() ...'
      end if

      return
      end
