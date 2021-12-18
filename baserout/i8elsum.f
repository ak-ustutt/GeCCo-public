      integer(8) function i8elsum(ivec,nelmnt)
*
* Sum of elements of integer vector IVEC
* Explicit integer(8) version
*
      implicit none

      integer(8), intent(in) ::      
     &   nelmnt, ivec(nelmnt)

      integer(8) ::
     &   isum, ielmnt
*
      isum = 0
      do ielmnt = 1, nelmnt
        isum = isum + ivec(ielmnt)
      end do
*
      i8elsum = isum
*
      return
      end
