      subroutine sum_occ(n,occr,nblk)
      
      implicit none

      integer, intent(out) ::
     &     n
      integer, intent(in) ::
     &     nblk, occr(nblk)

      n = 0
      if (nblk.eq.0) return
      n = sum(occr(1:nblk))
      return
      end subroutine
