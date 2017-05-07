*----------------------------------------------------------------------*
      subroutine condense_occ(occ_csub, occ_asub,
     &                        hpvx_csub,hpvx_asub,
     &                        occ,njoined,hpvx_blk_seq)
*----------------------------------------------------------------------*
*     put occupation on occ into condensed arrays occ_csub/occ_asub
*     hpvx info is then contained on hpvx_c/asub
*     use get_num_subblk() for getting dimensions of occ_csub etc.
*     we also reorder according to hpvx_blk_seq
!  e.g.:
!    h p v x  h p v x  h p v x
!c  /1 0 1 0\/0 2 1 0\/0 0 0 0\
!a  \0 0 0 0/\0 2 0 0/\1 1 1 0/
!
! assuming hpvx_blk_seq = 4 2 3 1 -> x p v h :
!
!occ_csub =  2 1 1 1  ! number of creators
!hpvx_csub = 2 3 3 1  ! in space
!            p v v h
!           meaning the following creator string: p p ; v ; v ; h
!
!occ_asub  = 2 1 1 1
!hpvx_asub = 2 2 3 1
!            p p v h
!           meaning the following anihilator string: p p ; p ; v ; h
! NOTE: if called with idx_graph, collect the index of the associated graph to a subspace string on occ_csub, occ_asub
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
        call write_title(lulog,wst_dbg_subr,'condense_occ')
        call wrt_occ_n(lulog,occ,njoined)
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
        write(lulog,'(3x,"C: occ  ",10i5)') occ_csub(1:idx_c) 
        write(lulog,'(3x,"   hpvx ",10i5)') hpvx_csub(1:idx_c) 
        write(lulog,'(3x,"A: occ  ",10i5)') occ_asub(1:idx_a) 
        write(lulog,'(3x,"   hpvx ",10i5)') hpvx_asub(1:idx_a) 
      end if

      return
      end
*----------------------------------------------------------------------*
