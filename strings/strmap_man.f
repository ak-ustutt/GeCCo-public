*----------------------------------------------------------------------*
      subroutine strmap_man(
     &     igraph1,dag1,
     &     igraph2,dag2,
     &     igraph12,dag12,
     &     str_info,strmap_info,orb_info)
*----------------------------------------------------------------------*
*     manage string maps:
*     generate new string map, if necessary
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

      logical, intent(in) ::
     &     dag1, dag2, dag12
      integer, intent(in) ::
     &     igraph1(2,ngastp), igraph2(2,ngastp), igraph12(2,ngastp)
      type(strinf), intent(in) ::
     &     str_info
      type(strmapinf), intent(inout) ::
     &     strmap_info
      type(orbinf), intent(in) ::
     &     orb_info

      integer ::
     &     ngraph, ica, hpvx, idxmap, ica1, ica2, ica12,
     &     iocc1, iocc2, nsym, ifree, lenoff
      integer, pointer ::
     &     idx_strmap(:)

      if (ntest.ge.100) then
        write(luout,*) '===================='
        write(luout,*) ' This is strmap_man'
        write(luout,*) '===================='
        write(luout,*) 'igraph1: '
        call wrt_occ(luout,igraph1)
        write(luout,*) 'igraph2: '
        call wrt_occ(luout,igraph2)
      end if

      ngraph = str_info%ngraph
      idx_strmap => strmap_info%idx_strmap

      if (strmap_info%ffstrmap%unit.le.0)
     &     call quit(1,'strmap_man',
     &     'ffstrmap must be open when calling strmap_man')

      ! loop over HPV/CA
      do ica = 1, 2
        ica1 = ica
        ica2 = ica
        ica12 = ica
        if (dag1) ica1 = 3-ica
        if (dag2) ica2 = 3-ica
        if (dag12) ica12 = 3-ica
        do hpvx = 1, ngastp
          ! trivial maps (where one graph corresponds to 
          ! a zero occupation) are set later (if necessary at all)
          if (igraph1(hpvx,ica1).eq.0.or.igraph2(hpvx,ica2).eq.0) cycle

          if (ntest.ge.100) then
            write(luout,*) 'need map for: ',
     &           igraph1(hpvx,ica1),igraph2(hpvx,ica2)
          end if

          ! does primitive map exist?
          idxmap = ngraph*(igraph2(hpvx,ica2)-1)+igraph1(hpvx,ica1)
          if (idx_strmap(idxmap).gt.-1) then
            ! check that offsets are set
            if (.not.associated(strmap_info%offsets))
     &           call quit(1,'strmap_man','offsets not initialized?')
            if (.not.associated(strmap_info%offsets(idxmap)%msms).or.
     &          .not.associated(strmap_info%offsets(idxmap)%msmsgmgm))
     &           call quit(1,'strmap_man','offsets not initialized?')
            if (ntest.ge.100) then
              write(luout,*) 'I have this map already ...'
            end if
            cycle
          end if

          if (ntest.ge.100) then
            write(luout,*) 'I will generate this map ...'
          end if

          idx_strmap(idxmap) = strmap_info%idx_last+1
          
          ! get memory for offset arrays
          iocc1 = str_info%ispc_occ(igraph1(hpvx,ica1))
          iocc2 = str_info%ispc_occ(igraph2(hpvx,ica2))
          nsym = orb_info%nsym
          call mem_pushmark()
          ifree = mem_gotomark(strmaps)
          lenoff = (iocc1+1)*(iocc2+1)
          ifree = mem_alloc_int(strmap_info%offsets(idxmap)%msms,
     &                          lenoff,'mm_off')
          lenoff = lenoff*nsym*nsym
          ifree = mem_alloc_int(strmap_info%offsets(idxmap)%msmsgmgm,
     &                          lenoff,'mmgg_off')

          call mem_popmark()

          ! set primitive map
          call set_strmap(
     &         str_info%g(igraph1(hpvx,ica1)),
     &         str_info%ispc_typ(igraph1(hpvx,ica1)),
     &         str_info%ispc_occ(igraph1(hpvx,ica1)),
     &         str_info%igas_restr(1,1,1,1,igraph1(hpvx,ica1)),
C               ! ADAPT FOR OPEN SHELL  ^^^
     &         str_info%g(igraph2(hpvx,ica2)),
     &         str_info%ispc_typ(igraph2(hpvx,ica2)),
     &         str_info%ispc_occ(igraph2(hpvx,ica2)),
     &         str_info%igas_restr(1,1,1,1,igraph2(hpvx,ica2)),
C               ! ADAPT FOR OPEN SHELL  ^^^
     &         str_info%g(igraph12(hpvx,ica12)),
     &         str_info%ispc_typ(igraph12(hpvx,ica12)),
     &         str_info%ispc_occ(igraph12(hpvx,ica12)),
     &         idxmap,strmap_info,orb_info)

        end do
      end do

      return
      end
