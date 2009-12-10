*----------------------------------------------------------------------*
      subroutine get_strmap_blk_c(strmap,
     &     n1,n2,n12,
     &     iocc1,iocc2,lstr1,lstr2,
     &     igrph1,igrph2,igrph12,
     &     idxms1,idxms2,igam1,igam2,map_info,
     &     strmap_info,nsym,ngraph)
*----------------------------------------------------------------------*
*     version for "condensed" quantities
*----------------------------------------------------------------------*

      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'hpvxseq.h'
      include 'def_filinf.h'
      include 'def_strmapinf.h'

      integer, intent(out), target ::
     &     strmap(*)
      integer, intent(in) ::
     &     n1, n2, n12, nsym, ngraph, map_info(*)
      integer, intent(in), target ::
     &     iocc1(n1), iocc2(n2), lstr1(n1), lstr2(n2),
     &     igrph1(n1), igrph2(n2), idxms1(n1), idxms2(n2),
     &     igam1(n1), igam2(n2), igrph12(*)
      type(strmapinf), intent(in), target ::
     &     strmap_info

      type(strmap_offsets), pointer ::
     &     offsets(:)
      type(flpmap_offsets), pointer ::
     &     offsets_fc(:)
      logical ::
     &     is_one_map, is_fc_map
      integer, pointer ::
     &     curmap(:)
      integer ::
     &     ioffmap, hpvx_idx, hpvx, idx, idx_minf, idxbuf,
     &     idxgrgr, idxmsms, idxgmgm, ioff, ilen, idx1, idx2, idx12,
     &     nsplit, mxgraph
      integer, pointer ::
     &     iocc1tmp(:), iocc2tmp(:), lstr1tmp(:), lstr2tmp(:),
     &     igrph1tmp(:), igrph2tmp(:), idxms1tmp(:), idxms2tmp(:),
     &     igam1tmp(:), igam2tmp(:)


      mxgraph = strmap_info%mxgraph
      if (mxgraph.lt.ngraph)
     &     call quit(1,'get_strmap_blk_c',
     &     'you forgot to update the maps after adding a graph')
      offsets => strmap_info%offsets
      offsets_fc => strmap_info%offsets_fc
      if (n1.lt.0.or.n2.lt.0.or.n12.lt.0) then
        write(luout,*) 'n1, n2, n12: ',n1, n2, n12
        call quit(1,'get_strmap_blk_c',
     &     'n1, n2, or n12: check call list')
      end if

      ioffmap = 0
      idx_minf = 0
      idx12 = 0
      do while(idx12.lt.n12)
        idx12 = idx12+1
c dbg
c        print *,'idx12 = ',idx12
c dbg
        idx1 = 0
        idx2 = 0
        nsplit = 0
        iocc1tmp => iocc1
        iocc2tmp => iocc2
        igrph1tmp => igrph1
        igrph2tmp => igrph2
        idxms1tmp => idxms1
        idxms2tmp => idxms2
        igam1tmp => igam1
        igam2tmp => igam2
        lstr1tmp => lstr1
        lstr2tmp => lstr2
        idx_minf = idx_minf+1
        nsplit = nsplit + map_info(idx_minf)
        !  we allow nsplit=2 --> both indices from same occ.cls.
        !  but only if the other occ.cls. gives nsplit=0
        if (nsplit.gt.2) call quit(1,'strmap_blk_c','multi map needed')
        if (nsplit.eq.1) then
          idx_minf = idx_minf+1
          idx1 = map_info(idx_minf)
        else if (nsplit.eq.2) then
          idx_minf = idx_minf+2
          idx1 = map_info(idx_minf-1)
          idx2 = map_info(idx_minf)
          iocc2tmp => iocc1
          igrph2tmp => igrph1
          idxms2tmp => idxms1
          igam2tmp => igam1
          lstr2tmp => lstr1
        end if
        idx_minf = idx_minf+1
        nsplit = nsplit + map_info(idx_minf)
        if (nsplit.gt.2) call quit(1,'strmap_blk_c','multi map needed')
        if (map_info(idx_minf).eq.1) then
          idx_minf = idx_minf+1
          idx2 = map_info(idx_minf)
        else if (map_info(idx_minf).eq.2) then
          idx_minf = idx_minf+2
          idx1 = map_info(idx_minf-1)
          idx2 = map_info(idx_minf)
          iocc1tmp => iocc2
          igrph1tmp => igrph2
          idxms1tmp => idxms2
          igam1tmp => igam2
          lstr1tmp => lstr2
        end if
c dbg
c        print *,'idx1/idx2:',idx1,idx2
c dbg
        if (idx1.eq.0.and.idx2.eq.0) then
          write(luout,*) 'n1,n2,n12: ',n1,n2,n12
          write(luout,*) 'map_info:  ',map_info(1:2*(n1+n2)*n12)
          write(luout,*) 'last index:',idx_minf
          write(luout,*) 'last idx12:',idx12
          call quit(1,'get_strmap_blk_c','buggy map_info?')
        end if

        is_one_map = idx1.eq.0.or.idx2.eq.0
        
        is_fc_map = .false.
        if (is_one_map) then
          if (idx1.ne.0) is_fc_map = igrph1tmp(idx1).ne.igrph12(idx12)
          if (idx2.ne.0) is_fc_map = igrph2tmp(idx2).ne.igrph12(idx12)
        end if

        if (.not.is_one_map) then
          ! get map from file
          idxgrgr = (igrph2tmp(idx2)-1)*mxgraph+igrph1tmp(idx1)
          idxmsms = (idxms2tmp(idx2)-1)*(iocc1tmp(idx1)+1)
     &              +idxms1tmp(idx1)
          idxgmgm = (idxmsms-1)*nsym*nsym+(igam2tmp(idx2)-1)*nsym
     &                                                 +igam1tmp(idx1)
          ! length of map
          ilen = lstr1tmp(idx1)*lstr2tmp(idx2)
c          ! have a look at the buffering table
c          idxmap = idxgrgr*nsym*nsym+idxgmgm
          ! offset of string map
          ioff = strmap_info%idx_strmap(idxgrgr)-1
          ! plus offset of current ms/ms block
          ioff = ioff + offsets(idxgrgr)%msms(idxmsms)
          ! plus offset of gam/gam block
          ioff = ioff + offsets(idxgrgr)%msmsgmgm(idxgmgm)
c dbg
c          print *,'fetching 1: ioff, len = ',ioff,ilen
c dbg
          call mem_iget(strmap_info%ffstrmap,
     &                  strmap(ioffmap+1),ioff+1,ioff+ilen)          

          ioffmap = ioffmap+ilen
        else if (is_fc_map) then
          ! get map from file
          if (idx1.gt.0) then
            idxgrgr = (igrph1tmp(idx1)-1)*mxgraph+igrph12(idx12)
            idxmsms =  idxms1tmp(idx1)
            idxgmgm =  igam1tmp(idx1)
            ilen = lstr1tmp(idx1)
          else
            idxgrgr = (igrph2tmp(idx2)-1)*mxgraph+igrph12(idx12)
            idxmsms =  idxms2tmp(idx2)
            idxgmgm =  igam2tmp(idx2)
            ilen = lstr2tmp(idx2)
          end if

          ! offset of string map
          ioff = strmap_info%idx_fcmap(idxgrgr)-1
          ! plus offset of ms/gam block
          ioff = ioff +
     &         offsets_fc(idxgrgr)%msgm((idxmsms-1)*nsym+idxgmgm)
c dbg
c          print *,'fetching 2: ioff, len = ',ioff,ilen
c dbg
          call mem_iget(strmap_info%ffstrmap,
     &                  strmap(ioffmap+1),ioff+1,ioff+ilen)          

          ioffmap = ioffmap+ilen
          
        else 
          ! set trivial map
          if (idx1.gt.0) ilen = lstr1tmp(idx1)
          if (idx2.gt.0) ilen = lstr2tmp(idx2)
c dbg
c          print *,'setting trivial map of length : ',ilen
c dbg
          curmap => strmap(ioffmap+1:ioffmap+ilen)
          do idx = 1, ilen
            curmap(idx) = idx
          end do
          ioffmap = ioffmap+ilen
        end if
      end do
c dbg
c      print *,'total len of map: ',ioffmap
c dbg

      return
      end
