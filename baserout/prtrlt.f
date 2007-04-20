c-----------------------------------------------------------------------
      subroutine prtrlt(v,m)
      implicit real*8(a-h,o-z)
c
c     ----- print out the lower triangle of a symmetric matrix (stored
c           in packed canonical form (actually an upper triangle) !) -----
c
      include 'stdunit.h'
      dimension v(m*(m+1)/2)

      max=5
      imax = 0
      do while(imax.lt.m)
        imin = imax+1
        imax = min(imax+max,m)
        write(luout,'(/,5x,10(6x,i4,5x)/)') (i,i = imin,imax)
        do i=1,m
          ii = i*(i-1)/2
          mm = imin + ii
          kk = min(i,imax) + ii
          if(mm.le.kk) then
            write(luout,'(i4,1x,10e15.7)') i,(v(j),j=mm,kk)
          end if
        end do
      end do
      write(*,*)
      return
      end
