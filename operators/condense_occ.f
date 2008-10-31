*----------------------------------------------------------------------*
      subroutine condense_occ(occ_csub, occ_asub,
     &                        hpvx_csub,hpvx_asub,
     &                        occ,njoined,hpvx_blk_seq)
*----------------------------------------------------------------------*
*     put occupation on occ into condensed arrays occ_csub/occ_asub
*     hpvx info is then contained on hpvx_c/asub
*     use get_num_subblk() for getting dimensions of occ_csub etc.
*     we also reorder according to hpvx_blk_seq
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'stdunit.h'

      integer, parameter ::
     &     ntest = 00

      integer, intent(out) ::
     &     occ_csub(*), occ_asub(*),
     &     hpvx_csub(*), hpvx_asub(*)
      integer, intent(in) ::
     &     njoined, occ(ngastp,2,njoined), hpvx_blk_seq(ngastp)

      integer ::
     &     ijoin, hpvx, idx_hpvx, idx_c, idx_a

      if (ntest.ge.100) then
        call write_title(luout,wst_dbg_subr,'condense_occ')
        call wrt_occ_n(luout,occ,njoined)
      end if

      idx_c = 0
      idx_a = 0
      do idx_hpvx = 1, ngastp
        hpvx = hpvx_blk_seq(idx_hpvx)
        do ijoin = 1, njoined
          if (occ(hpvx,1,ijoin).gt.0) then 
            idx_c = idx_c + 1
            occ_csub(idx_c) = occ(hpvx,1,ijoin)
            hpvx_csub(idx_c) = hpvx
          end if
          if (occ(hpvx,2,ijoin).gt.0) then
            idx_a = idx_a + 1
            occ_asub(idx_a) = occ(hpvx,2,ijoin)
            hpvx_asub(idx_a) = hpvx
          end if
        end do
      end do
      
      if (ntest.ge.100) then
        write(luout,'(3x,"C: occ  ",10i5)') occ_csub(1:idx_c) 
        write(luout,'(3x,"   hpvx ",10i5)') hpvx_csub(1:idx_c) 
        write(luout,'(3x,"A: occ  ",10i5)') occ_asub(1:idx_a) 
        write(luout,'(3x,"   hpvx ",10i5)') hpvx_asub(1:idx_a) 
      end if

      return
      end
*----------------------------------------------------------------------*
