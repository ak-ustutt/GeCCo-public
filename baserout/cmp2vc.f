      subroutine cmp2vc(vec1,vec2,ndim,thres)
c
c compare two double precision  vectors vec1,and vec2
c
c only elements that differs by more than thre are printed
c
c     adapated from jeppe olsen's original
c
      implicit none

      include 'stdunit.h'

      integer, intent(in) ::
     &     ndim
      real(8), intent(in) ::
     &     thres,
     &     vec1(ndim),vec2(ndim)
c
      real(8) ::
     &     xmxdif, dif
      integer ::
     &     imxplc, i

      xmxdif = 0.0d0
      imxplc = 0
      write(lulog,*) ' comparison of two vectors '
      write(lulog,*) '      vector1      vector2        difference '
      do i = 1, ndim
        dif = vec1(i) - vec2 ( i )
        if( abs(dif ) .ge. xmxdif ) then
          xmxdif = abs(dif)
          imxplc = i
        end if
        if( abs ( dif ) .gt. thres ) then
          write(lulog,'(2x,i5,3e15.8)') i,vec1(i),vec2(i),dif
        end if
      end do
c
      if( xmxdif .eq. 0.0d0 ) then
        write(lulog,*) ' the two vectors are identical '
      else
        write(lulog,*) ' size and last place of largest deviation ',
     &  xmxdif,imxplc
      end if
c
      return
      end
