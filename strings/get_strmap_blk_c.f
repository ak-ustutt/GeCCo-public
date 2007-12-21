*----------------------------------------------------------------------*
      subroutine get_strmap_blk_c(strmap,
     &     n1,n2,n12,
     &     iocc1,iocc2,lstr1,lstr2,
     &     igrph1,igrph2,
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
     &     n1, n2, n12,
     &     iocc1(n1), iocc2(n2), lstr1(n2), lstr2(n2),
     &     igrph1(n1), igrph2(n2), idxms1(n2), idxms2(n2),
     &     igam1(n1), igam2(n2), map_info(*),
     &     nsym, ngraph
      type(strmapinf), intent(in), target ::
     &     strmap_info

      type(strmap_offsets), pointer ::
     &     offsets(:)
      integer, pointer ::
     &     curmap(:)
      integer ::
     &     ioffmap, hpvx_idx, hpvx, idx, idxmap, idx_minf, idxbuf,
     &     idxgrgr, idxmsms, idxgmgm, ioff, ilen, idx1, idx2, idx12,
     &     nsplit, mxgraph


      mxgraph = strmap_info%mxgraph
      if (mxgraph.lt.ngraph)
     &     call quit(1,'get_strmap_blk_c',
     &     'you forgot to update the maps after adding a graph')
      offsets => strmap_info%offsets
c dbg
c      print *,'n1, n2, n12: ',n1, n2, n12
c dbg

      ioffmap = 0
      idx_minf = 0
      idx12 = 0
      do while(idx12.lt.n12)
        idx12 = idx12+1
        idx1 = 0
        idx2 = 0
        idx_minf = idx_minf+1
        nsplit = map_info(idx_minf)
        ! nsplit<=1 was tested in strmap_man_c
        if (nsplit.gt.1) call quit(1,'strmap_blk_c','multi map needed')
        if (nsplit.eq.1) then
          idx_minf = idx_minf+1
          idx1 = map_info(idx_minf)
        end if
        idx_minf = idx_minf+1
        nsplit = map_info(idx_minf)
        if (nsplit.gt.1) call quit(1,'strmap_blk_c','multi map needed')
        if (nsplit.eq.1) then
          idx_minf = idx_minf+1
          idx2 = map_info(idx_minf)
        end if
c dbg
c        print *,'idx1/idx2:',idx1,idx2
c dbg

        if (idx1.gt.0.and.idx2.gt.0) then
          ! get map from file
          idxgrgr = (igrph2(idx2)-1)*mxgraph+igrph1(idx1)
          idxmsms = (idxms2(idx2)-1)*(iocc1(idx1)+1)+idxms1(idx1)
          idxgmgm = (idxmsms-1)*nsym*nsym+(igam2(idx2)-1)*nsym
     &                                                 +igam1(idx1)
          ! length of map
          ilen = lstr1(idx1)*lstr2(idx2)
          ! have a look at the buffering table
          idxmap = idxgrgr*nsym*nsym+idxgmgm
          ! offset of string map
          ioff = strmap_info%idx_strmap(idxgrgr)-1
          ! plus offset of current ms/ms block
          ioff = ioff + offsets(idxgrgr)%msms(idxmsms)
          ! plus offset of gam/gam block
          ioff = ioff + offsets(idxgrgr)%msmsgmgm(idxgmgm)
          call mem_iget(strmap_info%ffstrmap,
     &                  strmap(ioffmap+1),ioff+1,ioff+ilen)          

          ioffmap = ioffmap+ilen
        else if (idx1.gt.0.or.idx2.gt.0) then
          ! set trivial map
          if (idx1.gt.0) ilen = lstr1(idx1)
          if (idx2.gt.0) ilen = lstr2(idx2)
c dbg
c          print *,'setting trivial map of length : ',ilen
c dbg
          curmap => strmap(ioffmap+1:ioffmap+ilen)
          do idx = 1, ilen
            curmap(idx) = idx
          end do
          ioffmap = ioffmap+ilen
        else
          write(luout,*) 'map_info:  ',map_info(1:2*(n1+n2)*n12)
          write(luout,*) 'last index:',idx_minf
          call quit(1,'get_strmap_blk_c','buggy map_info?')
        end if
      end do
c dbg
c      print *,'total len of map: ',ioffmap
c dbg

      return
      end
