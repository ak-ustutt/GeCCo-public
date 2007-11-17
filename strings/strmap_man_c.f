*----------------------------------------------------------------------*
      subroutine strmap_man_c(
     &     maxbuffer,
     &     igraph1r,n1,
     &     igraph2r,n2,
     &     igraph12r,n12,map_info,
     &     str_info,strmap_info,orb_info)
*----------------------------------------------------------------------*
*     manage string maps:
*     generate new string map, if necessary
*     version for condensed contraction info
*     maxbuffer returns the maximum number of integer words needed 
*     for a MS/IRREP block
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
     &     n1, n2, n12,
     &     igraph1r(n1), igraph2r(n2), igraph12r(n12), map_info(*)
      type(strinf), intent(in) ::
     &     str_info
      type(strmapinf), intent(inout) ::
     &     strmap_info
      type(orbinf), intent(in) ::
     &     orb_info

      integer ::
     &     ngraph, ica, hpvx, idxmap, idx1, idx2, idx12,
     &     nsplit, idx_minf,
     &     iocc1, iocc2, nsym, ifree, lenoff, maxstr
      integer, pointer ::
     &     idx_strmap(:)
      integer, external ::
     &     ifndmax

      if (ntest.ge.100) then
        call write_title(luout,wst_dbg_subr,'This is strmap_man_c')
        write(luout,*) 'igraph1r: ',igraph1r
        write(luout,*) 'igraph2r: ',igraph2r
        write(luout,*) 'igraph12r: ',igraph12r
        write(luout,*) 'map_info: ',map_info(1:n12*2*(n1+n2))
      end if

      maxbuffer = 0

      ngraph = str_info%ngraph
      idx_strmap => strmap_info%idx_strmap

      if (strmap_info%ffstrmap%unit.le.0)
     &     call quit(1,'strmap_man',
     &     'ffstrmap must be open when calling strmap_man')

      idx_minf = 0
      idx12 = 0
      do while(idx12.lt.n12)
        idx12 = idx12+1
        idx1 = 0
        idx2 = 0
        idx_minf = idx_minf+1
        nsplit = map_info(idx_minf)
        if (nsplit.gt.1) call quit(1,'strmap_man_c','multi map needed')
        if (nsplit.eq.1) then
          idx_minf = idx_minf+1
          idx1 = map_info(idx_minf)
        end if
        idx_minf = idx_minf+1
        nsplit = map_info(idx_minf)
        if (nsplit.gt.1) call quit(1,'strmap_man_c','multi map needed')
        if (nsplit.eq.1) then
          idx_minf = idx_minf+1
          idx2 = map_info(idx_minf)
        end if
        ! do not consider trivial maps
        if (idx1.eq.0.or.idx2.eq.0) then
          ! set buffer space for trivial map
          if (idx1.gt.0) then
            maxstr = ifndmax(
     &           str_info%g(igraph1r(idx1))%lenstr_gm,1,
     &           (str_info%ispc_occ(igraph1r(idx1))+1)*orb_info%nsym,1)
          else if (idx2.gt.0) then
            maxstr = ifndmax(
     &           str_info%g(igraph2r(idx2))%lenstr_gm,1,
     &           (str_info%ispc_occ(igraph2r(idx2))+1)*orb_info%nsym,1)
          else
            maxstr = 1
          end if
          maxbuffer = maxbuffer + maxstr
          cycle
        end if
        if (ntest.ge.100) then
          write(luout,*) 'need map for: ',
     &         igraph1r(idx1),igraph2r(idx2)
        end if

        ! does primitive map exist?
        idxmap = ngraph*(igraph2r(idx2)-1)+igraph1r(idx1)
        if (idx_strmap(idxmap).gt.-1) then
          ! check that offsets are set
          if (.not.associated(strmap_info%offsets))
     &         call quit(1,'strmap_man','offsets not initialized?')
          if (.not.associated(strmap_info%offsets(idxmap)%msms).or.
     &          .not.associated(strmap_info%offsets(idxmap)%msmsgmgm))
     &           call quit(1,'strmap_man','offsets not initialized?')
          if (ntest.ge.100) then
            write(luout,*) 'I have this map already ...'
          end if
          maxbuffer = maxbuffer + strmap_info%maxlen_blk(idxmap)
          cycle
        end if

        if (ntest.ge.100) then
          write(luout,*) 'I will generate this map ...'
        end if

        idx_strmap(idxmap) = strmap_info%idx_last+1
          
        ! get memory for offset arrays
        iocc1 = str_info%ispc_occ(igraph1r(idx1))
        iocc2 = str_info%ispc_occ(igraph2r(idx2))
        nsym = orb_info%nsym
        call mem_pushmark()
        ifree = mem_gotomark(strmaps)
        lenoff = (iocc1+1)*(iocc2+1)
        ifree = mem_alloc_int(strmap_info%offsets(idxmap)%msms,
     &       lenoff,'mm_off')
        lenoff = lenoff*nsym*nsym
        ifree = mem_alloc_int(strmap_info%offsets(idxmap)%msmsgmgm,
     &                          lenoff,'mmgg_off')

        call mem_popmark()

        ! set primitive map
        call set_strmap(
     &         str_info%g(igraph1r(idx1)),
     &         str_info%ispc_typ(igraph1r(idx1)),
     &         str_info%ispc_occ(igraph1r(idx1)),
     &         str_info%igas_restr(1,1,1,igraph1r(idx1)),
     &         str_info%g(igraph2r(idx2)),
     &         str_info%ispc_typ(igraph2r(idx2)),
     &         str_info%ispc_occ(igraph2r(idx2)),
     &         str_info%igas_restr(1,1,1,igraph2r(idx2)),
     &         str_info%g(igraph12r(idx12)),
     &         str_info%ispc_typ(igraph12r(idx12)),
     &         str_info%ispc_occ(igraph12r(idx12)),
     &         idxmap,strmap_info,orb_info)

        maxbuffer = maxbuffer + strmap_info%maxlen_blk(idxmap)

      end do

      return
      end
