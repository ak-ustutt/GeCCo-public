
      integer function iocc_in_vec(inum,ivec,ndim)
*
*     count the number of times, elements of value inum occur in vector
*     ivec(ndim)
*
      implicit none

      integer, intent(in) ::
     &     inum, ndim, ivec(ndim)

      integer ::
     &     iocc, ii

      iocc = 0

      do ii = 1, ndim
        if (ivec(ii).eq.inum) iocc = iocc+1
      end do

      iocc_in_vec = iocc

      return
      end
