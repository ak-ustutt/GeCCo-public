      integer(4) function i4elsum(ivec,nelmnt)
*
* Sum of elements of integer vector IVEC
* Explicit integer(4) version
*
      implicit none

      integer(4), intent(in) ::      
     &   nelmnt, ivec(nelmnt)

      integer(4) ::
     &   isum, ielmnt
*
      isum = 0
      do ielmnt = 1, nelmnt
        isum = isum + ivec(ielmnt)
      end do
*
      i4elsum = isum
*
      return
      end
