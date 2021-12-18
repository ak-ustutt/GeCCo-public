      subroutine diavc(vecout,vecin,fac,diag,shift,ndim)
c
c     vecout(i)=fac*vecin(i)/(diag(i)+shift)
c
c     adapted from jeppe olsen's diavc2
c
      implicit none
      include 'stdunit.h'
      character(len=*),parameter::
     &     i_am="diavc"
      integer,parameter::
     &     ntest=00
      real(8), parameter ::
     &     thresh = 1d-10 

      integer, intent(in) ::
     &     ndim
      real(8), intent(in) ::
     &     fac, shift, diag(ndim)
      real(8), intent(inout) ::
     &     vecout(ndim), vecin(ndim)
c
      integer ::
     &     i
      real(8) ::
     &     divide

      if (ntest.ge.10)then
         write(lulog,*) "thresh, shift",thresh,shift
      end if
      do i=1,ndim
        divide=diag(i)+shift
        if(abs(divide).le.thresh) divide=thresh
        vecout(i)=fac*vecin(i)/divide
      end do

      return
      end
