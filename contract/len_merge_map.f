*----------------------------------------------------------------------*
      integer function len_merge_map(mmap,nvtx_tgt)
*----------------------------------------------------------------------*
*     return length of map
*----------------------------------------------------------------------*
      implicit none

      integer, intent(in) ::
     &     mmap(*), nvtx_tgt

      integer ::
     &     ivtx_tgt, ivtx, idx

      idx = 0
      do ivtx_tgt = 1, nvtx_tgt
        ivtx = mmap(idx+1)
        idx = idx + ivtx + 1
        ivtx = mmap(idx+1)
        idx = idx + ivtx + 1
      end do

      len_merge_map = idx

      return
      end
