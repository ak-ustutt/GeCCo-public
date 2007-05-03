*----------------------------------------------------------------------*
      integer function idx_str_blk(idxc,idxa,lenc,lena,iocc,ihpvsq)
*----------------------------------------------------------------------*
*
*     calculate address within block from given indices for 
*     C and A (each for H/P/V) with the convention:
*     HPV sequence as given by ihpvsq (innermost to outermost)
*     within each H/P/V: C is innermost loop
*      
*----------------------------------------------------------------------*

      implicit none

      include 'opdim.h'

      integer, intent(in) ::
     &     idxc(ngastp), idxa(ngastp),
     &     lenc(ngastp), lena(ngastp),
     &     iocc(ngastp,2), ihpvsq(ngastp)

      integer ::
     &     len_last, ihpvdx, ihpv

      ! we have to add a one to get the index
      idx_str_blk = 1
      len_last = 1
      do ihpvdx = 1, ngastp
        ihpv = ihpvsq(ihpvdx)
        ! the string indices are still idx-1, so the
        ! following is correct
        if (iocc(ihpv,1).gt.0)
     &       idx_str_blk = idx_str_blk+idxc(ihpv)*len_last
        if (iocc(ihpv,1).gt.0)
     &       len_last = len_last*lenc(ihpv)
        if (iocc(ihpv,2).gt.0)
     &       idx_str_blk = idx_str_blk+idxa(ihpv)*len_last
        if (iocc(ihpv,2).gt.0)
     &       len_last = len_last*lena(ihpv)
      end do
      
      return
      end
