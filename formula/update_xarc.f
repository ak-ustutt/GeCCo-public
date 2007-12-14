*----------------------------------------------------------------------*
      subroutine update_xarc(contr,xarcs_raw,nxarc_raw,ivtxder,nj_res)
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'def_contraction.h'

      integer, parameter ::
     &     ntest = 00

      type(contraction), intent(inout) ::
     &     contr
      integer, intent(in) ::
     &     nxarc_raw,ivtxder,nj_res
      type(cntr_arc), intent(in) ::
     &     xarcs_raw(nxarc_raw)

      integer ::
     &     nvtx, idxsuper, ixarc, ivtx, iarc_raw

      if (ntest.ge.100) then
        call write_title(luout,wst_dbg_subr,'update_xarc')
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

      ! make sure that we have enough space on contr
      call resize_contr(contr,0,0,contr%nxarc+nxarc_raw,0)

      ! no external lines so far?
      if (contr%nxarc.eq.0) then
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
        call quit(1,'update_xarc',
     &       'did not have time for the complicated part yet')
      end if

      return
      end
