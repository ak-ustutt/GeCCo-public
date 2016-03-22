*----------------------------------------------------------------------*
      subroutine get_strmap_blk(strmap,
     &     iocc1,iocc2,lstr1,lstr2,
     &     igrph1,igrph2,ms1,ms2,igam1,igam2,
     &     strmap_info,nsym,ngraph)
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
     &     iocc1(ngastp), iocc2(ngastp), lstr1(ngastp), lstr2(ngastp),
     &     igrph1(ngastp), igrph2(ngastp), ms1(ngastp), ms2(ngastp),
     &     igam1(ngastp), igam2(ngastp), nsym, ngraph
      type(strmapinf), intent(in), target ::
     &     strmap_info

      type(strmap_offsets), pointer ::
     &     offsets(:)
      integer, pointer ::
     &     curmap(:)
      integer ::
     &     ioffmap, hpvx_idx, hpvx, idx, idxmap, idxbuf,
     &     idxgrgr, idxms1, idxms2, idxmsms, idxgmgm, ioff, ilen

c      ! does not work as long as ifort has problems with
c      ! user-defined types in explicit interfaces ...
c      ilen = 0
c      do hpvx = 1, ngastp
c        ilen = ilen+lstr1(hpvx)*lstr2(hpvx)
c      end do
c
c      ifree = mem_alloc_int(strmap,ilen,'strmap')

      offsets => strmap_info%offsets

      ioffmap = 0
      do hpvx_idx = 1, ngastp
        hpvx = hpvxseq(hpvx_idx)
c        hpvx = hpvx_idx
        if (igrph1(hpvx).gt.0.and.igrph2(hpvx).gt.0) then
          ! get map from file
          idxgrgr = (igrph2(hpvx)-1)*ngraph+igrph1(hpvx)
          idxms1 = (iocc1(hpvx)-ms1(hpvx))/2 + 1
          idxms2 = (iocc2(hpvx)-ms2(hpvx))/2 + 1
          idxmsms = (idxms2-1)*(iocc1(hpvx)+1)+idxms1
          idxgmgm = (idxmsms-1)*nsym*nsym+(igam2(hpvx)-1)*nsym
     &                                                   +igam1(hpvx)
          ! length of map
          ilen = lstr1(hpvx)*lstr2(hpvx)
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
        else if (igrph1(hpvx).gt.0.or.igrph2(hpvx).gt.0) then
          ! set trivial map
          ilen = max(lstr1(hpvx),lstr2(hpvx))
          curmap => strmap(ioffmap+1:ioffmap+ilen)
          do idx = 1, ilen
            curmap(idx) = idx
          end do
          ioffmap = ioffmap+ilen
        end if
      end do

      return
      end
