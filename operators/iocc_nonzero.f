      logical function iocc_nonzero(iocc)

      implicit none

      include 'opdim.h'

      integer ::
     &     iocc(ngastp,2)

      integer ::
     &     idx, jdx

      iocc_nonzero = .false.
      do jdx = 1, 2
        do idx = 1, ngastp
          iocc_nonzero = iocc_nonzero.or.iocc(idx,jdx).ne.0
        end do
      end do

      return
      end
