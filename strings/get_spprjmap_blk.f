*----------------------------------------------------------------------*
      subroutine get_spprjmap_blk(strmap,
     &     nblk,iocc,lstr,igrph,idxms,igam,
     &     strmap_info,nsym,ngraph)
*----------------------------------------------------------------------*
*     obtain spin-projection string-map from buffer
*     call strmap_man_spprj before to make sure that it exists
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
     &     nblk,
     &     iocc(nblk), lstr(nblk), igrph(nblk), idxms(nblk),
     &     igam(nblk), nsym, ngraph
      type(strmapinf), intent(in), target ::
     &     strmap_info

      type(flpmap_offsets), pointer ::
     &     offsets(:)
      integer, pointer ::
     &     curmap(:)
      integer ::
     &     ioffmap, hpvx_idx, hpvx, idx, idxmap, idx_minf, idxbuf,
     &     idxgraph, ioff, ilen, iblk

      if (strmap_info%mxgraph.lt.ngraph)
     &     call quit(1,'get_flipmap_blk',
     &     'you forgot to update the maps after adding a graph')
      offsets => strmap_info%offsets_spprj

      ioffmap = 0

      do iblk = 1, nblk
        if (iocc(iblk).gt.1) then
c dbg
c          print *,'iblk = ',iblk
c      print *,'fetching G,MS/G: ',igrph(iblk),idxms(iblk),igam(iblk)
c dbg
          idxgraph = igrph(iblk)
          idxmap = igam(iblk)
          ! offset of string map (actually stored for ms, just in 
          ! case we need more maps, not only those starting from the
          ! canonical sceleton part of the ME list)
          ioff = offsets(idxgraph)%ms(idxms(iblk))
c dbg
c          print *,'idxgraph = ',idxgraph
c          print *,'ioff0: ',ioff
c dbg
          ! plus offset of current gm block
          ioff = ioff + offsets(idxgraph)%msgm(idxmap)
c dbg
c          print *,'idxgraph,idxmap: ',idxgraph,idxmap
c          print *,'idxmap,ioff1: ',idxmap,offsets(idxgraph)%msgm(idxmap)
c dbg
          ilen = lstr(iblk)*(2**iocc(iblk)-1)
c dbg          
c          print *,'ilen = ',ilen
c          print *,'fetching from: ',ioff+1,
c     &                              ioff+ilen
c          print *,'fetching map: ',ioffmap+1,ioffmap+ilen
c dbg
c dbg
c          print *,'call to mem_iget in get_flipmap_blk ',ioff,ilen
c dbg
          call mem_iget(strmap_info%ffstrmap,
     &         strmap(ioffmap+1),ioff+1,ioff+ilen)          
c dbg
c           print '(x,5i6,2x,5i6)',strmap(ioffmap+1:ioffmap+ilen)
c dbg
          ioffmap = ioffmap+ilen
        else
          ! set trivial map
          ilen = lstr(iblk)
c dbg
c          print *,'setting trivial map of length : ',ilen
c          print *,'settig map: ',ioffmap+1,ioffmap+ilen
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
