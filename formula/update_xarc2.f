*----------------------------------------------------------------------*
      subroutine update_xarc2(contr,xarcs_raw,nxarc_raw,ivtxder,nj_res)
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'def_contraction.h'

      integer, parameter ::
     &     ntest = 100

      type(contraction), intent(inout) ::
     &     contr
      integer, intent(in) ::
     &     nxarc_raw,ivtxder,nj_res
      type(cntr_arc), intent(in) ::
     &     xarcs_raw(nxarc_raw)

      integer ::
     &     nvtx, nxarc, idxsuper, ixarc, ivtx, iarc_raw, idx, iarc, idel
      integer, pointer ::
     &     del_list(:)

      if (ntest.ge.100) then
        call write_title(luout,wst_dbg_subr,'update_xarc2')
        write(luout,*) 'no. of external lines: ',contr%nxarc
        write(luout,*) 'raw xarc list: (',nxarc_raw,')'
        do ixarc = 1, nxarc_raw
          write(luout,'(3x,i3,x,i3,3x,4i4)')
     &         xarcs_raw(ixarc)%link,
     &         xarcs_raw(ixarc)%occ_cnt(1:ngastp,1)
          write(luout,'(3x,3x,x,3x,3x,4i4)')
     &         xarcs_raw(ixarc)%occ_cnt(1:ngastp,2)
        end do
      end if

      nvtx = contr%nvtx
      nxarc = contr%nxarc

      ! make sure that we have enough space on contr
      call resize_contr(contr,0,0,contr%nxarc+nxarc_raw,0)

      ! no external lines so far?
      if (nxarc.eq.0) then
        ! yes --> life is quite easy then
        idxsuper = 1 ! index of super-vertices of result
        ixarc    = 0 ! counter for xarcs
        do ivtx = 1, ivtxder-1
          do iarc_raw = 1, nxarc_raw
            if (xarcs_raw(iarc_raw)%link(1).eq.ivtx) then
              ixarc = ixarc+1
              contr%nxarc = ixarc
              contr%xarc(ixarc) = xarcs_raw(iarc_raw)
              contr%xarc(ixarc)%link(2) = idxsuper
            end if
          end do
        end do
        if (ixarc.gt.0) idxsuper = idxsuper+1
        ! quick fix for densities:
        if (ivtxder.eq.1.and.nj_res.eq.2) idxsuper=2
        do ivtx = ivtxder, nvtx
          do iarc_raw = 1, nxarc_raw
            if (xarcs_raw(iarc_raw)%link(1).eq.ivtx) then
              ixarc = ixarc+1
              contr%nxarc = ixarc
              contr%xarc(ixarc) = xarcs_raw(iarc_raw)
              contr%xarc(ixarc)%link(2) = idxsuper
            end if
          end do          
        end do
      else
        ! No --> must combine original and new xarcs referring to the 
        ! same contraction vertex.
        idxsuper = 1 ! index of super-vertices of result
        ixarc    = nxarc ! counter for xarcs
        do ivtx = 1, ivtxder-1
          do iarc_raw = 1, nxarc_raw
            if (xarcs_raw(iarc_raw)%link(1).eq.ivtx) then
              ixarc = ixarc+1
              contr%nxarc = ixarc
              contr%xarc(ixarc) = xarcs_raw(iarc_raw)
              contr%xarc(ixarc)%link(2) = idxsuper
            end if
          end do
        end do
        if (ixarc.gt.0) idxsuper = idxsuper+1
        ! quick fix for densities:
        if (ivtxder.eq.1.and.nj_res.eq.2) idxsuper=2
        do ivtx = ivtxder, nvtx
          do iarc_raw = 1, nxarc_raw
            if (xarcs_raw(iarc_raw)%link(1).eq.ivtx) then
              ixarc = ixarc+1
              contr%nxarc = ixarc
              contr%xarc(ixarc) = xarcs_raw(iarc_raw)
              contr%xarc(ixarc)%link(2) = idxsuper
            end if
          end do          
        end do

        ! Must now check if any of the new xarcs involve the same 
        ! vertices as the old xarcs. If so, they should be combined and
        ! the contraction updated
        if(ixarc.gt.0)then
          idel = 0
          allocate(del_list(ixarc))
          del_list(1:ixarc) = 0
          do idx = 1, nxarc
            do iarc = nxarc+1, ixarc
              if(contr%xarc(iarc)%link(1).eq.contr%xarc(idx)%link(1)
     &             .and.
     &             contr%xarc(iarc)%link(2).eq.contr%xarc(idx)%link(2))
     &             then
                contr%xarc(idx)%occ_cnt(1:ngastp,1:2) =
     &               contr%xarc(idx)%occ_cnt(1:ngastp,1:2) +
     &               contr%xarc(iarc)%occ_cnt(1:ngastp,1:2)
                idel = idel + 1
                contr%xarc(iarc)%link(1) = 0
                del_list(idel) = iarc
              endif
            enddo
          enddo

          do idx = 1, idel
            do iarc = del_list(idel), ixarc-1
              contr%xarc(iarc)%link(1:2) = contr%xarc(iarc+1)%link(1:2)
              contr%xarc(iarc)%occ_cnt(1:ngastp,1:2) =
     &             contr%xarc(iarc+1)%occ_cnt(1:ngastp,1:2)
              del_list(idx+1:idel) = del_list(idx+1:idel) - 1

            enddo
          enddo

          deallocate(del_list)

          contr%nxarc = contr%nxarc - idel

        endif

c        call quit(1,'update_xarc',
c     &       'did not have time for the complicated part yet')
      end if

      return
      end
