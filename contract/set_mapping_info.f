*----------------------------------------------------------------------*
      subroutine set_mapping_info(map_info_c,map_info_a,
     &                  iocc1,njoined1,
     &                  iocc2,njoined2,
     &                  iocc12,merge_map,njoined12,hpvx_blk_seq)
*----------------------------------------------------------------------*
*     set for each non-zero index of iocc12:
*        ...,nidx1,idx1_1,idx1_2,...,idx1_nidx1,
*            nidx2,idx2_1,idx2_2,...,idx1_nidx2,....
*----------------------------------------------------------------------*
      implicit none
      
      include 'opdim.h'
      include 'stdunit.h'

      integer, parameter ::
     &     ntest = 00

      integer, intent(in) ::
     &     njoined1, njoined2, njoined12,
     &     iocc1(ngastp,2,njoined1), iocc2(ngastp,2,njoined2),
     &     iocc12(ngastp,2,njoined12),
     &     merge_map(*), hpvx_blk_seq(ngastp)

      integer, intent(out) ::
     &     map_info_c(*),
     &     map_info_a(*)

      integer ::
     &     ica, hpvx, idx_merge_map, ijoin12, idx_hpvx,
     &     icount, nvtx, njoin1, idx, idx_base, ivtx, na12, nc12
      integer ::
     &     idxseq1(ngastp,2,njoined1), idxseq2(ngastp,2,njoined2)

      if (ntest.ge.100) then
        call write_title(luout,wst_dbg_subr,'set_mapping_info')
        write(luout,*) 'occ1'
        call wrt_occ_n(luout,iocc1,njoined1)
        write(luout,*) 'occ2'
        call wrt_occ_n(luout,iocc1,njoined2)
        write(luout,*) 'occ12'
        call wrt_occ_n(luout,iocc12,njoined12)
        write(luout,*) 'merge_map:'
        idx_base = 1
        do ijoin12 = 1, njoined12
          write(luout,'(i3,"->",i3,": ",10i3)')
     &                        ijoin12, merge_map(idx_base),
     &              merge_map(idx_base+1:idx_base+merge_map(idx_base))
          idx_base = idx_base + merge_map(idx_base) + 1
        end do
      end if

      ! fill map with -1
c dbg
c      print *,'njoined12 -> ' , njoined12
c dbg
      call get_num_subblk(nc12,na12,iocc12,njoined12)
c      nc12 = sum(iocc12(1:ngastp,1,1:njoined12))
c      na12 = sum(iocc12(1:ngastp,1,1:njoined12))
c dbg
c      print *,'dims: ',nc12,njoined1,njoined2
c      print *,'dims: ',na12,njoined1,njoined2
c dbg
      if (nc12.gt.0) map_info_c(1:nc12*2*(njoined1+njoined2)) = -1
      if (na12.gt.0) map_info_a(1:na12*2*(njoined1+njoined2)) = -1
c dbg
c      print *,'1111'
c        if (nc12.gt.0)
c     &      write(luout,'(10i5)')
c     &       map_info_c(1:nc12*2*(njoined1+njoined2))
c
c dbg
      ! set index sequence for iocc1/iocc2
      call set_idx_seq(idxseq1,iocc1,njoined1,1,hpvx_blk_seq)
      call set_idx_seq(idxseq2,iocc2,njoined2,1,hpvx_blk_seq)
c dbg
c      print *,'2222'
c        if (nc12.gt.0)
c     &      write(luout,'(10i5)')
c     &       map_info_c(1:nc12*2*(njoined1+njoined2))
c
c dbg

      ! loop over iocc12 occupations and assemble map
      do ica = 1, 2
c dbg
c        print *,'ica = ',ica
c dbg
        idx_base = 1
        do idx_hpvx = 1, ngastp
          hpvx = hpvx_blk_seq(idx_hpvx)
          idx_merge_map = 1
          do ijoin12 = 1, njoined12
c dbg
c      print *,'top'
c        if (nc12.gt.0)
c     &      write(luout,'(10i5)')
c     &       map_info_c(1:nc12*2*(njoined1+njoined2))
c
c dbg
            if (iocc12(hpvx,ica,ijoin12).eq.0) cycle
            ! set counters for contributions to zero
            if (ica.eq.1) map_info_c(idx_base) = 0
            if (ica.eq.2) map_info_a(idx_base) = 0
            icount = 0
            ! number of iocc1 vertices contribution contribution to
            ! current iocc12 vertex
            nvtx = merge_map(idx_merge_map)
c dbg
c            print *,'1) idx_merge_map, nvtx: ',idx_merge_map, nvtx
c dbg
            ! loop over contribution vertices
            do idx = 1, nvtx
              idx_merge_map = idx_merge_map+1
              ivtx = merge_map(idx_merge_map)
c dbg
c            print *,'a) idx_merge_map, ivtx: ',idx_merge_map, ivtx
c dbg
              ! nonzero?
              if (iocc1(hpvx,ica,ivtx).eq.0) cycle
              ! increment counters
              icount = icount+1
              if (ica.eq.1) map_info_c(idx_base) = icount
              if (ica.eq.2) map_info_a(idx_base) = icount
c dbg
c              print *,'1) to ',idx_base,' -> ',icount
c              print *,'   to ',idx_base+icount,' -> ',
c     &             idxseq1(hpvx,ica,ivtx)
c dbg
              ! and store index of the contribution
              if (ica.eq.1) map_info_c(idx_base+icount)
     &             = idxseq1(hpvx,1,ivtx)
              if (ica.eq.2) map_info_a(idx_base+icount)
     &             = idxseq1(hpvx,2,ivtx)
            end do
            
c dbg
c      print *,'mid'
c        if (nc12.gt.0)
c     &      write(luout,'(10i5)')
c     &       map_info_c(1:nc12*2*(njoined1+njoined2))
c
c dbg
            ! the same for iocc2 contributions
c dbg
c            print *,'idx_base, icount: ',idx_base, icount
c dbg
            idx_base = idx_base+icount+1
            idx_merge_map = idx_merge_map+1
c dbg
c            print *,' idx_merge_map: ',idx_merge_map
c dbg

            ! set counters for contributions to zero
            if (ica.eq.1) map_info_c(idx_base) = 0
            if (ica.eq.2) map_info_a(idx_base) = 0
            icount = 0
            ! number of iocc1 vertices contribution contribution to
            ! current iocc12 vertex
            nvtx = merge_map(idx_merge_map)
c dbg
c            print *,'2) idx_merge_map, nvtx: ',idx_merge_map, nvtx
c dbg
            ! loop over contribution vertices
            do idx = 1, nvtx
              idx_merge_map = idx_merge_map+1
              ivtx = merge_map(idx_merge_map)
c dbg
c            print *,'a) idx_merge_map, ivtx: ',idx_merge_map, ivtx
c dbg
              ! nonzero?
              if (iocc2(hpvx,ica,ivtx).eq.0) cycle
              ! increment counters
              icount = icount+1
              if (ica.eq.1) map_info_c(idx_base) = icount
              if (ica.eq.2) map_info_a(idx_base) = icount
c dbg
c              call wrt_occ(luout,idxseq2)
c              print *,' -> ',hpvx,ica,ivtx
c              print *,'2) to ',idx_base,' -> ',icount
c              print *,'   to ',idx_base+icount,' -> ',
c     &             idxseq2(hpvx,ica,ivtx)
c dbg
              ! and store index of the contribution
              if (ica.eq.1) map_info_c(idx_base+icount)
     &             = idxseq2(hpvx,1,ivtx)
              if (ica.eq.2) map_info_a(idx_base+icount)
     &             = idxseq2(hpvx,2,ivtx)
            end do

            ! increment base address of map
            idx_base = idx_base+icount+1
c dbg
c            print *,'bot'
c        if (nc12.gt.0)
c     &      write(luout,'(10i5)')
c     &       map_info_c(1:nc12*2*(njoined1+njoined2))
c
c dbg

          end do

        end do
        
      end do ! ica

      if (ntest.ge.100) then
        write(luout,*) 'map_info_c: '
        if (nc12.gt.0)
     &      write(luout,'(10i5)')
     &       map_info_c(1:nc12*2*(njoined1+njoined2))
        write(luout,*) 'map_info_a: '
        if (na12.gt.0)
     &      write(luout,'(10i5)')
     &       map_info_a(1:na12*2*(njoined1+njoined2))
      end if

      return
      end
