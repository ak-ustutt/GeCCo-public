      real*8 function fndmnx(vector,ndim,minmax)
c
c     find smallest(minmax=1) or largest(minmax=2)
c     absolute value of elements in vector
c
c     adapted from Jeppe Olsen's version
c
      implicit none
      integer, intent(in) ::
     &     ndim, minmax
      real(8), intent(in) ::
     &     vector(ndim)
      integer ::
     &     i
c
      if(minmax.eq.1) then
        fndmnx=abs(vector(1))
        do i=2,ndim
          fndmnx=min(fndmnx,abs(vector(i)))
        end do
      end if
c
      if(minmax.eq.2) then
        fndmnx=abs(vector(1))
        do i=2,ndim
          fndmnx=max(fndmnx,abs(vector(i)))
        end do
      end if
c
      if(minmax.eq.-1) then
        fndmnx=vector(1)
        do i=2,ndim
          fndmnx=min(fndmnx,vector(i))
        end do
      end if
c
      if(minmax.eq.-2) then
        fndmnx=vector(1)
        do i=2,ndim
          fndmnx=max(fndmnx,vector(i))
        end do
      end if
      
      return
      end
