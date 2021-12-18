      logical function occ_is_diag_blk(occ,njoined)

      ! TRUE if creation and annihilation occupation are the same

      implicit none

      include 'opdim.h'

      integer ::
     &     njoined,
     &     occ(ngastp,2,njoined)

      integer ::
     &     idx, ijoin_c, ijoin_a

      occ_is_diag_blk = .true.
      do ijoin_c = 1, njoined
        ijoin_a = njoined+1-ijoin_c
        do idx = 1, ngastp
          occ_is_diag_blk = occ_is_diag_blk.and.
     &                      occ(idx,1,ijoin_c).eq.occ(idx,2,ijoin_a)
        end do
      end do
      
      return
      end
