*----------------------------------------------------------------------*
      subroutine get_diag_tuple1(idx,ca,
     &                        occ,njoined,hpvx_blk_seq)
*----------------------------------------------------------------------*
*     given a diagonal occupation, the indices of the c/a (sub)strings
*     which characterize the diagonal elements are returned.
*     There is no unique way of splitting the indices, here it is
*     chosen to split the occupation matrix into an "upper" and a
*     "lower" part, e.g. ("x" marks indices of lower part)
*            /0  0  0  0\
*            \0  0  1  0/
*            /1x 0  2  0\ <-- (only for occupied indices there must be a
*            \1  0  2x 0/      c/a flip for the middle vertex)
*            /0  0  1x 0\
*            \0  0  0  0/
*
*     matthias, feb 2010
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'stdunit.h'

      integer, parameter ::
     &     ntest = 00

      integer, intent(out) ::
     &     idx(*), ca(*)
      integer, intent(in) ::
     &     njoined, occ(ngastp,2*njoined), hpvx_blk_seq(ngastp)

      integer ::
     &     ijoin, hpvx, idx_hpvx, idx_c, idx_a, ica, ii, idx_dia,
     &     iflip

      if (ntest.ge.100) then
        call write_title(lulog,wst_dbg_subr,'get_diag_tuple1')
        call wrt_occ_n(lulog,occ,njoined)
      end if

      idx_c = 0
      idx_a = 0
      idx_dia = 0
      do idx_hpvx = 1, ngastp
        hpvx = hpvx_blk_seq(idx_hpvx)
        do ii = 1, 2*njoined
          ijoin = (ii-1)/njoined+1
          ica = 2-mod(ii,2)
          if (occ(hpvx,ii).gt.0) then
            if (hpvx.eq.1.and.mod(njoined,2).eq.1) then !occ: c/a flipped
              iflip = 1
            else
              iflip = 0
            end if
            if (ica.eq.1) then
              idx_c = idx_c + 1
              if (ii.gt.njoined-iflip) then
                idx_dia = idx_dia + 1
                idx(idx_dia) = idx_c
                ca(idx_dia) = 1
              end if
            else
              idx_a = idx_a + 1
              if (ii.gt.njoined+iflip) then
                idx_dia = idx_dia + 1
                idx(idx_dia) = idx_a
                ca(idx_dia) = 2
              end if
            end if
          end if
        end do
      end do
      
      if (ntest.ge.100) then
        write(lulog,'(3x,"idx: ",10i5)') idx(1:idx_dia) 
        write(lulog,'(3x,"c/a: ",10i5)') ca(1:idx_dia) 
      end if

      return
      end
*----------------------------------------------------------------------*
