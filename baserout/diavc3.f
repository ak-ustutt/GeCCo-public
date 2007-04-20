      subroutine diavc3(vecout,vecin,diag,shift,ndim,vdsv)
*
* vecout(i)=vecin(i)/(diag(i)+shift)
*
* vdsv = sum(i) vecin(i) ** 2 /( diag(i) + shift )
*
*     adapted from jeppen olsen's original
* 
      implicit none

      include 'stdunit.h'

      integer, intent(in) ::
     &     ndim
      real(8), intent(in) ::
     &     vecin(ndim), diag(ndim), shift
      real(8), intent(out) ::
     &     vecout(ndim), vdsv

      real(8), parameter ::
     &     thres = 1d-10
      integer, parameter ::
     &     ntest = 00

      real(8) ::
     &     divide
      integer ::
     &     i

*
      vdsv = 0.0d0
      do i=1,ndim
*
        divide=diag(i)+shift
        if(abs(divide).le.thres) divide=thres
*
        vdsv = vdsv + vecin(i)*vecin(i) /divide
        vecout(i)=vecin(i)/divide
*
      end do
*
      if(ntest.ge.100) then
        write(luout,*) 'diavc3 : vecin, diag,vecout '
        do i = 1, ndim
          write(luout,'(3e15.8)') vecin(i),diag(i),vecout(i)
        end do
      end if
*
      return
      end
