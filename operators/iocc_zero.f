      logical function iocc_zero(iocc)

      implicit none

      include 'opdim.h'

      integer ::
     &     iocc(ngastp,2)

      integer ::
     &     idx, jdx

      iocc_zero = .true.
      do jdx = 1, 2
        do idx = 1, ngastp
          iocc_zero = iocc_zero.and.iocc(idx,jdx).eq.0
        end do
      end do

      return
      end
