*----------------------------------------------------------------------*
      subroutine strmap_man_c(mode,
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
*     mode = 1: for contraction; mode = 2: for reordering
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
     &     n1, n2, n12, igraph12r(n12), map_info(*), mode
      integer, intent(in), target ::
     &     igraph1r(n1), igraph2r(n2)
      type(strinf), intent(in) ::
     &     str_info
      type(strmapinf), intent(inout) ::
     &     strmap_info
      type(orbinf), intent(in) ::
     &     orb_info

      integer ::
     &     mxgraph, ica, hpvx, idxmap, idxgr, idx1, idx2, idx12,
     &     idx1mx, idx2mx,
     &     nsplit, idx_minf,
     &     iocc1, iocc2, nsym, ifree, lenoff, maxstr
      logical ::
     &     error, is_one_map, is_fc_map
      integer, pointer ::
     &     idx_strmap(:), idx_fcmap(:), igraph1tmp(:), igraph2tmp(:)
      integer, external ::
     &     ifndmax

      if (ntest.ge.100) then
        call write_title(lulog,wst_dbg_subr,'This is strmap_man_c')
        write(lulog,*) 'igraph1r: ',igraph1r
        write(lulog,*) 'igraph2r: ',igraph2r
        write(lulog,*) 'igraph12r: ',igraph12r
        write(lulog,*) 'map_info: ',map_info(1:n12*2*(n1+n2))
      end if

      maxbuffer = 0

      mxgraph = strmap_info%mxgraph
      if (mxgraph.lt.str_info%ngraph)
     &     call quit(1,'strmap_man_c',
     &     'you forgot to update the maps after adding a new graph')
      idx_strmap => strmap_info%idx_strmap
      idx_fcmap => strmap_info%idx_fcmap

      if (strmap_info%ffstrmap%unit.le.0)
     &     call quit(1,'strmap_man_c',
     &     'ffstrmap must be open when calling strmap_man')

      idx_minf = 0
      idx12 = 0
      idx1mx = 0
      idx2mx = 0
      error = .false.
      do while(idx12.lt.n12)
        idx12 = idx12+1
        idx1 = 0
        idx2 = 0
        nsplit = 0
        igraph1tmp => igraph1r
        igraph2tmp => igraph2r
        idx_minf = idx_minf+1
        nsplit = nsplit + map_info(idx_minf)
        !  we allow nsplit=2 --> both indices from same occ.cls.
        !  but only if the other occ.cls. gives nsplit=0
        if (nsplit.gt.2.or.nsplit.eq.2.and.mode.eq.1) then
          error = .true.
          exit
        end if
        if (nsplit.eq.1) then
          idx_minf = idx_minf+1
          idx1 = map_info(idx_minf)
        else if (nsplit.eq.2) then
          idx_minf = idx_minf+2
          idx1 = map_info(idx_minf-1)
          idx2 = map_info(idx_minf)
          igraph2tmp => igraph1r
        end if
        idx_minf = idx_minf+1
        nsplit = nsplit + map_info(idx_minf)
        if (nsplit.gt.2.or.map_info(idx_minf).eq.2.and.mode.eq.1) then
          error = .true.
          exit
        end if
        if (map_info(idx_minf).eq.1) then
          idx_minf = idx_minf+1
          idx2 = map_info(idx_minf)
        else if (map_info(idx_minf).eq.2) then
          idx_minf = idx_minf+2
          idx1 = map_info(idx_minf-1)
          idx2 = map_info(idx_minf)
          igraph1tmp => igraph2r
        end if

c        if ((idx1.gt.0.and.idx1.lt.idx1mx) .or.
c     &       idx2.gt.0.and.idx2.lt.idx2mx) then
c          error = .true.
c          exit
c        end if
c        idx1mx = max(idx1,idx1mx)
c        idx2mx = max(idx2,idx2mx)

        is_one_map = idx1.eq.0.or.idx2.eq.0
        
        is_fc_map = .false.
        if (is_one_map) then
          if (idx1.ne.0) is_fc_map = 
     &          igraph1tmp(idx1).ne.igraph12r(idx12)
          if (idx2.ne.0) is_fc_map = 
     &          igraph2tmp(idx2).ne.igraph12r(idx12)
        end if

        ! do not consider trivial maps
        if (is_one_map.and..not.is_fc_map) then
          ! set buffer space for trivial map
          if (idx1.gt.0) then
            maxstr = ifndmax(
     &          str_info%g(igraph1tmp(idx1))%lenstr_gm,1,
     &          (str_info%ispc_occ(igraph1tmp(idx1))+1)*orb_info%nsym,1)
          else if (idx2.gt.0) then
            maxstr = ifndmax(
     &          str_info%g(igraph2tmp(idx2))%lenstr_gm,1,
     &          (str_info%ispc_occ(igraph2tmp(idx2))+1)*orb_info%nsym,1)
          else
            maxstr = 1
          end if
          maxbuffer = maxbuffer + maxstr
          cycle
        end if

        if (ntest.ge.100) then
          if (.not.is_fc_map) then
            write(lulog,*) 'need map for: ',
     &         igraph1tmp(idx1),igraph2tmp(idx2)
          else
            if (idx1.gt.0) then
              write(lulog,*) 'need fc map for: ',
     &             igraph1tmp(idx1),igraph12r(idx12)
            else
              write(lulog,*) 'need fc map for: ',
     &             igraph2tmp(idx2),igraph12r(idx12)
            end if
          end if
        end if

        ! process fc-map
        if (is_fc_map) then
          if (idx1.gt.0) then
            idxmap = (igraph1tmp(idx1)-1)*mxgraph + igraph12r(idx12)
            idxgr  = igraph1tmp(idx1)
          else
            idxmap = (igraph2tmp(idx2)-1)*mxgraph + igraph12r(idx12)
            idxgr  = igraph2tmp(idx2)
          end if
          if (idx_fcmap(idxmap).gt.-1) then
            ! check that offsets are set
            if (.not.associated(strmap_info%offsets_fc))
     &         call quit(1,'strmap_man_c'
     &                    ,'offsets (fc) not initialized?')
            if (.not.associated(strmap_info%offsets_fc(idxmap)%ms)
     &           .or.
     &         .not.associated(strmap_info%offsets_fc(idxmap)%msgm))
     &         call quit(1,'strmap_man_c'
     &                    ,'offsets (fc) not initialized?')
            if (ntest.ge.100) then
              write(lulog,*) 'I have this map already ...'
            end if
            maxbuffer = maxbuffer
     &           + strmap_info%maxlen_blk_fc(idxmap)
            cycle
          end if

          if (ntest.ge.100) then
            write(lulog,*) 'I will generate this fc map ...'
          end if

          idx_fcmap(idxmap) = strmap_info%idx_last + 1

          ! get memory for offset arrays
          iocc1 = str_info%ispc_occ(idxgr)
          nsym = orb_info%nsym
          call mem_pushmark()
          ifree = mem_gotomark(strmaps)
          lenoff = (iocc1+1)
          ifree = mem_alloc_int(strmap_info%offsets_fc(idxmap)%ms,
     &         lenoff,'m_off')
          lenoff = lenoff*nsym
          ifree = mem_alloc_int(strmap_info%offsets_fc(idxmap)%msgm,
     &                          lenoff,'mg_off')

          call mem_popmark()
           
          ! set primitive fc map
          call set_fcmap(
     &         str_info%g(idxgr),
     &         str_info%ispc_typ(idxgr),
     &         str_info%ispc_occ(idxgr),
     &         str_info%igas_restr(1,1,1,1,idxgr),
C               ! ADAPT FOR OPEN SHELL  ^^^
     &         str_info%g(igraph12r(idx12)),
     &         str_info%ispc_typ(igraph12r(idx12)),
     &         str_info%ispc_occ(igraph12r(idx12)),
     &         str_info%igas_restr(1,1,1,1,igraph12r(idx12)),
     &         idxmap,strmap_info,orb_info)

          maxbuffer = maxbuffer + strmap_info%maxlen_blk_fc(idxmap)

          cycle
        end if

        ! does primitive map exist?
        idxmap = mxgraph*(igraph2tmp(idx2)-1)+igraph1tmp(idx1)
        if (idx_strmap(idxmap).gt.-1) then
          ! check that offsets are set
          if (.not.associated(strmap_info%offsets))
     &         call quit(1,'strmap_man','offsets not initialized?')
          if (.not.associated(strmap_info%offsets(idxmap)%msms).or.
     &          .not.associated(strmap_info%offsets(idxmap)%msmsgmgm))
     &           call quit(1,'strmap_man','offsets not initialized?')
          if (ntest.ge.100) then
            write(lulog,*) 'I have this map already ...'
          end if
          maxbuffer = maxbuffer + strmap_info%maxlen_blk(idxmap)
          cycle
        end if

        if (ntest.ge.100) then
          write(lulog,*) 'I will generate this map ...'
        end if

        idx_strmap(idxmap) = strmap_info%idx_last+1
          
        ! get memory for offset arrays
        iocc1 = str_info%ispc_occ(igraph1tmp(idx1))
        iocc2 = str_info%ispc_occ(igraph2tmp(idx2))
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
     &         str_info%g(igraph1tmp(idx1)),
     &         str_info%ispc_typ(igraph1tmp(idx1)),
     &         str_info%ispc_occ(igraph1tmp(idx1)),
     &         str_info%igas_restr(1,1,1,1,igraph1tmp(idx1)),
C               ! ADAPT FOR OPEN SHELL  ^^^
     &         str_info%g(igraph2tmp(idx2)),
     &         str_info%ispc_typ(igraph2tmp(idx2)),
     &         str_info%ispc_occ(igraph2tmp(idx2)),
     &         str_info%igas_restr(1,1,1,1,igraph2tmp(idx2)),
C               ! ADAPT FOR OPEN SHELL  ^^^
     &         str_info%g(igraph12r(idx12)),
     &         str_info%ispc_typ(igraph12r(idx12)),
     &         str_info%ispc_occ(igraph12r(idx12)),
     &         str_info%igas_restr(1,1,1,1,igraph12r(idx12)),
     &         idxmap,strmap_info,orb_info)

        maxbuffer = maxbuffer + strmap_info%maxlen_blk(idxmap)

      end do

      if (error) then
        write(lulog,*) 'igraph1r: ',igraph1r
        write(lulog,*) 'igraph2r: ',igraph2r
        write(lulog,*) 'igraph12r: ',igraph12r
        write(lulog,*) 'map_info: ',map_info(1:n12*2*(n1+n2))
        write(lulog,*) 'mode: ',mode

        idx_minf = 0
        idx12 = 0
        do while(idx12.lt.n12)
          idx12 = idx12+1
          idx1 = 0
          idx2 = 0
          idx_minf = idx_minf+1
          nsplit = map_info(idx_minf)
          write(lulog,*) 'idx12, nsplit1: ',idx12,nsplit
          write(lulog,*) 'indices: ',
     &         map_info(idx_minf+1:idx_minf+nsplit)
          idx_minf = idx_minf+nsplit

          idx_minf = idx_minf+1
          nsplit = map_info(idx_minf)
          write(lulog,*) 'idx12, nsplit2: ',idx12,nsplit
          write(lulog,*) 'indices: ',
     &         map_info(idx_minf+1:idx_minf+nsplit)
          idx_minf = idx_minf+nsplit
        end do
        
        call quit(1,'strmap_man_c','Unsupported mapping type!')
      end if

      if (ntest.ge.100) then
        write(lulog,*) 'leaving strmap_man_c() ...'
      end if

      return
      end
