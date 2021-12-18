*----------------------------------------------------------------------*
!>    finds the maximum value in an integer vector
!>    
!>    steps over an integer vector and returns the largest value encountered
!>    @param[in] ivec integer vector
!>    @param[in] idxoff offset where the function starts
!>    @param[in] lvec number of steps
!>    @param[in] inc stepsize
!>    @return 
*----------------------------------------------------------------------*
      integer function ifndmax(ivec,idxoff,lvec,inc)
*----------------------------------------------------------------------*

      implicit none

      integer, intent(in) ::
     &     ivec(*), lvec, inc, idxoff

      integer ::
     &     i, imx, idx

      imx = -huge(imx)
      idx = idxoff
      do i = 1, lvec
        imx = max(imx,ivec(idx))
        idx = idx + inc
      end do
      
      ifndmax = imx

      return
      end
