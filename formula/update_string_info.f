      subroutine update_string_info(contr)

      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'def_contraction.h'

      integer, parameter ::
     &     ntest = 00
      
      type(contraction), intent(inout) ::
     &     contr

      integer ::
     &     iarc, ivtx1, ivtx2, vtx_last, ioff1, nidx1, ioff2, nidx2,
     &     idx, ica, icad, ihpvx, ncnt, icnt, ii, jj, sign
      integer, pointer ::
     &     idxpairs(:), ivtxoff(:), nidxvtx(:)
      integer, pointer ::
     &     occ_cnt(:,:)
      type(string_element), pointer ::
     &     string(:)

      if (ntest.ge.100) then
        write(lulog,*) 'update: string on entry'
        call print_string(contr%contr_string,contr%nidx,contr%nvtx)
        call print_string_cnt(contr%contr_string,contr%nidx,contr%nvtx)
        call print_string_idx(contr%contr_string,contr%nidx,contr%nvtx)
      end if
      
      string => contr%contr_string
      
      allocate(ivtxoff(contr%nvtx),nidxvtx(contr%nvtx),
     &     idxpairs(contr%nidx))
      ivtxoff = 0
      nidxvtx = 0

      vtx_last = 0
      do idx = 1, contr%nidx
        if (string(idx)%vtx.gt.vtx_last) then
          vtx_last = string(idx)%vtx
          ivtxoff(vtx_last) = idx-1
        end if
        nidxvtx(vtx_last) =  nidxvtx(vtx_last)+1
      end do

c      write(lulog,'(a,10i4)') 'off: ',ivtxoff(1:contr%nvtx)
c      write(lulog,'(a,10i4)') 'len: ',nidxvtx(1:contr%nvtx)
        
      do iarc = 1, contr%narc
        ivtx1 = contr%arc(iarc)%link(1)
        ivtx2 = contr%arc(iarc)%link(2)
        occ_cnt => contr%arc(iarc)%occ_cnt
        ioff1 = ivtxoff(ivtx1)
        nidx1 = nidxvtx(ivtx1)        
        ioff2 = ivtxoff(ivtx2)
        nidx2 = nidxvtx(ivtx2)
        do ica = 1, 2
          icad = 3-ica
          do ihpvx = 1, ngastp
            if (occ_cnt(ihpvx,ica)==0) cycle
            ncnt = occ_cnt(ihpvx,ica)

! find index pairs: there should be paired indices on both vertices
            icnt = 0
            search_loop: do ii = ioff1+1,ioff1+nidx1
              if (string(ii)%ca==ica.and.
     &            string(ii)%hpvx==ihpvx) then
                idx = string(ii)%idx
                do jj = ioff2+1,ioff2+nidx2
                  if (string(jj)%ca==icad.and.
     &                string(jj)%hpvx==ihpvx.and.
     &                string(jj)%idx==idx) then
                    icnt = icnt+1
                    idxpairs(icnt) = ii
                    icnt = icnt+1
                    idxpairs(icnt) = jj
                    if (icnt==2*ncnt) exit search_loop
                  end if
                end do
              end if
              if (ii==ioff1+nidx1.and.icnt<2*ncnt)
     &        call quit(1,'update_string_info','unsuccessful algorithm')
            end do search_loop
            do icnt = 1, 2*ncnt
              string(idxpairs(icnt))%cnt = iarc
            end do
          end do
        end do
      end do

      sign = 0
      do ivtx1 = 1, contr%nvtx
        ioff1 = ivtxoff(ivtx1)
        nidx1 = nidxvtx(ivtx1)
        call sort_string(sign,string(ioff1+1:ioff1+nidx1),nidx1)
      end do

      deallocate(ivtxoff,nidxvtx,idxpairs)
      
      if (ntest.ge.100) then
        write(lulog,*) 'update: string on exit'
        call print_string(contr%contr_string,contr%nidx,contr%nvtx)
        call print_string_cnt(contr%contr_string,contr%nidx,contr%nvtx)
        call print_string_idx(contr%contr_string,contr%nidx,contr%nvtx)
      end if
      
      end subroutine
