      subroutine occ2caocc(ca_occ,occ,nj,nblk)

      implicit none
      include 'opdim.h'

      integer, intent(in) ::
     &     nj, nblk, occ(ngastp,2,nj,nblk)
      integer, intent(out) ::
     &     ca_occ(2,nblk)

      integer ::
     &     ij, iblk

      do iblk = 1, nblk
        ca_occ(1:2,iblk) = 0
        do ij = 1, nj
          ca_occ(1,iblk) = ca_occ(1,iblk) + sum(occ(1:ngastp,1,ij,iblk))
          ca_occ(2,iblk) = ca_occ(2,iblk) + sum(occ(1:ngastp,2,ij,iblk))
        end do
      end do

      return
      end
