*----------------------------------------------------------------------*
      subroutine set_mapping_info(map_info_c,map_info_a,
     &                  nca_mode,
     &                  iocc1,njoined1,dagger1,
     &                  iocc2,njoined2,dagger2,
     &                  iocc12,merge_map,njoined12,hpvx_blk_seq)
*----------------------------------------------------------------------*
*     set for each non-zero index of iocc12:
*        ...,nidx1,idx1_1,idx1_2,...,idx1_nidx1,
*            nidx2,idx2_1,idx2_2,...,idx1_nidx2,....
*     nca_mode = 0: set both
*     nca_mode = 1: set C map only
*     nca_mode = 1: set A map only
*----------------------------------------------------------------------*
      implicit none
      
      include 'opdim.h'
      include 'stdunit.h'

      integer, parameter ::
     &     ntest = 00

      logical, intent(in) ::
     &     dagger1, dagger2
      integer, intent(in) ::
     &     nca_mode,
     &     njoined1, njoined2, njoined12,
     &     iocc1(ngastp,2,njoined1), iocc2(ngastp,2,njoined2),
     &     iocc12(ngastp,2,njoined12),
     &     merge_map(*), hpvx_blk_seq(ngastp)

      integer, intent(out) ::
     &     map_info_c(*),
     &     map_info_a(*)

      integer ::
     &     ica, icast, icand, ica1, ica2,
     &     hpvx, ijoin12, idx_hpvx,
     &     icount, nvtx1, nvtx2,
     &     idx_merge_map, ioff1_merge_map, ioff2_merge_map, 
     &     njoin1, idx, idx_base, ivtx, na12, nc12
      integer ::
     &     idxseq1(ngastp,2,njoined1), idxseq2(ngastp,2,njoined2)

      if (ntest.ge.100) then
        call write_title(luout,wst_dbg_subr,'set_mapping_info')
        write(luout,*) 'occ1 (dag = ',dagger1,')'
        call wrt_occ_n(luout,iocc1,njoined1)
        write(luout,*) 'occ2 (dag = ',dagger2,')'
        call wrt_occ_n(luout,iocc2,njoined2)
        write(luout,*) 'occ12'
        call wrt_occ_n(luout,iocc12,njoined12)
        write(luout,*) 'merge_map:'
        idx_base = 1
        do ijoin12 = 1, njoined12
          write(luout,'(i3,"->",i3,": ",10i3)')
     &                        ijoin12, merge_map(idx_base),
     &              merge_map(idx_base+1:idx_base+merge_map(idx_base))
          idx_base = idx_base + merge_map(idx_base) + 1
          write(luout,'(3x,"->",i3,": ",10i3)')
     &                             merge_map(idx_base),
     &              merge_map(idx_base+1:idx_base+merge_map(idx_base))
          idx_base = idx_base + merge_map(idx_base) + 1
        end do
      end if

      if (nca_mode.eq.0) then
        icast = 1
        icand = 2
      else if (nca_mode.eq.1) then
        icast = 1
        icand = 1
      else if (nca_mode.eq.2) then
        icast = 2
        icand = 2
      else
        call quit(1,'set_mapping_info','unknown nca_mode')
      end if
     
      ! fill map with -1
      call get_num_subblk(nc12,na12,iocc12,njoined12)
      if (nc12.gt.0.and.icast.eq.1)
     &     map_info_c(1:nc12*2*(njoined1+njoined2)) = -1
      if (na12.gt.0.and.icand.eq.2)
     &     map_info_a(1:na12*2*(njoined1+njoined2)) = -1
      ! set index sequence for iocc1/iocc2
      call set_idx_seq(idxseq1,iocc1,njoined1,.true.,hpvx_blk_seq)
      call set_idx_seq(idxseq2,iocc2,njoined2,.true.,hpvx_blk_seq)

      ! loop over iocc12 occupations and assemble map
      do ica = icast, icand
        ica1 = ica
        if (dagger1) ica1 = 3-ica
        ica2 = ica
        if (dagger2) ica2 = 3-ica
        idx_base = 1
        do idx_hpvx = 1, ngastp
          hpvx = hpvx_blk_seq(idx_hpvx)
          idx_merge_map = 1
          do ijoin12 = 1, njoined12
            ! read merge-map:
            nvtx1 = merge_map(idx_merge_map)
            ioff1_merge_map = idx_merge_map
            idx_merge_map = idx_merge_map+nvtx1+1
            nvtx2 = merge_map(idx_merge_map)
            ioff2_merge_map = idx_merge_map
            idx_merge_map = idx_merge_map+nvtx2+1

            if (iocc12(hpvx,ica,ijoin12).eq.0) cycle
c dbg
c            print *,'hpvx,ica,ijoin12,occ12:',hpvx,ica,ijoin12,
c     &           iocc12(hpvx,ica,ijoin12)
c dbg
            ! set counters for contributions to zero
            if (ica.eq.1) map_info_c(idx_base) = 0
            if (ica.eq.2) map_info_a(idx_base) = 0
            icount = 0
            ! number of iocc1 vertices contribution contribution to
            ! current iocc12 vertex
            ! loop over contribution vertices
            do idx = 1, nvtx1
              ivtx = merge_map(ioff1_merge_map+idx)
              ! nonzero?
              if (iocc1(hpvx,ica1,ivtx).eq.0) cycle
              ! increment counters
              icount = icount+1
              if (ica.eq.1) map_info_c(idx_base) = icount
              if (ica.eq.2) map_info_a(idx_base) = icount
              ! and store index of the contribution
              if (ica.eq.1) map_info_c(idx_base+icount)
     &             = idxseq1(hpvx,ica1,ivtx)
              if (ica.eq.2) map_info_a(idx_base+icount)
     &             = idxseq1(hpvx,ica1,ivtx)
            end do
            
            ! the same for iocc2 contributions
            idx_base = idx_base+icount+1

            ! set counters for contributions to zero
            if (ica.eq.1) map_info_c(idx_base) = 0
            if (ica.eq.2) map_info_a(idx_base) = 0
            icount = 0
            ! number of iocc1 vertices contribution contribution to
            ! current iocc12 vertex
            ! loop over contribution vertices
            do idx = 1, nvtx2
              ivtx = merge_map(ioff2_merge_map+idx)
              ! nonzero?
              if (iocc2(hpvx,ica2,ivtx).eq.0) cycle
              ! increment counters
              icount = icount+1
              if (ica.eq.1) map_info_c(idx_base) = icount
              if (ica.eq.2) map_info_a(idx_base) = icount
              ! and store index of the contribution
              if (ica.eq.1) map_info_c(idx_base+icount)
     &             = idxseq2(hpvx,ica2,ivtx)
              if (ica.eq.2) map_info_a(idx_base+icount)
     &             = idxseq2(hpvx,ica2,ivtx)
            end do

            ! increment base address of map
            idx_base = idx_base+icount+1

          end do

        end do
        
      end do ! ica

      if (ntest.ge.100) then
        write(luout,*) 'map_info_c: '
        if (icast.eq.1.and.nc12.gt.0)
     &      write(luout,'(10i5)')
     &       map_info_c(1:nc12*2*(njoined1+njoined2))
        write(luout,*) 'map_info_a: '
        if (icand.eq.2.and.na12.gt.0)
     &      write(luout,'(10i5)')
     &       map_info_a(1:na12*2*(njoined1+njoined2))
      end if

      return
      end
