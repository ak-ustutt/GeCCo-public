*----------------------------------------------------------------------*
      subroutine diag_condensed_occ(out_csub, out_asub,
     &                        in_csub,in_asub,
     &                        occ,njoined,hpvx_blk_seq)
*----------------------------------------------------------------------*
*     merge in_csub and in_asub to a condensed string of which
*     the input strings represent the diagonal elements
*
*     matthias, fall 2009 (adopted from condense_occ)
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'stdunit.h'

      integer, parameter ::
     &     ntest = 00

      integer, intent(out) ::
     &     out_csub(*), out_asub(*)
      integer, intent(in) ::
     &     in_csub(*), in_asub(*),
     &     njoined, occ(ngastp,2,njoined), hpvx_blk_seq(ngastp)

      integer ::
     &     ijoin, ihpv_dx, ihpv, ica, idxc, idxa, idx1, idx2

      if (ntest.ge.100) then
        call write_title(lulog,wst_dbg_subr,'diag_condensed_occ')
        call wrt_occ_n(lulog,occ,njoined)
      end if

      idxc = 0
      idxa = 0
      idx1 = 0
      idx2 = 0
      do ihpv_dx = 1, ngastp
        ihpv = hpvx_blk_seq(ihpv_dx)
        do ica = 1, 2
          do ijoin = 1, njoined
            if (occ(ihpv,ica,ijoin).eq.0) cycle
            idxc = idxc + 1
            if (ica.eq.1) then
              idx1 = idx1 + 1
              out_csub(idxc) = in_csub(idx1)
            else
              idx2 = idx2 + 1
              out_csub(idxc) = in_asub(idx2)
            end if
          end do
        end do
      end do
      idx1 = 0
      idx2 = 0
      do ihpv_dx = 1, ngastp
        ihpv = hpvx_blk_seq(ihpv_dx)
        do ica = 2, 1, -1
          do ijoin = 1, njoined
            if (occ(ihpv,ica,ijoin).eq.0) cycle
            idxa = idxa + 1
            if (ica.eq.1) then
              idx1 = idx1 + 1
              out_asub(idxa) = in_csub(idx1)
            else
              idx2 = idx2 + 1
              out_asub(idxa) = in_asub(idx2)
            end if
          end do
        end do
      end do      

      if (ntest.ge.100) then
        write(lulog,'(3x,"IN:  C ",10i5)') in_csub(1:idx1) 
        write(lulog,'(3x,"     A ",10i5)') in_asub(1:idx2) 
        write(lulog,'(3x,"OUT: C ",10i5)') out_csub(1:idxc) 
        write(lulog,'(3x,"     A ",10i5)') out_asub(1:idxa) 
      end if

      return
      end
*----------------------------------------------------------------------*
