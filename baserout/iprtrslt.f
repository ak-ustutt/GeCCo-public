c-----------------------------------------------------------------------
      subroutine iprtrslt(iv,m)
      implicit real*8(a-h,o-z)
c
c     ----- print out the lower triangle of an (anti-symmetric) matrix 
c           stored as strict lower triangle (actually: upper triangle) -----
c
      include 'stdunit.h'
      dimension iv(m*(m-1)/2)

      max=5
      imax = 0
      do while(imax.lt.m-1)
        imin = imax+1
        imax = min(imax+max,m-1)
        write(lulog,'(/,5x,10(6x,i4,5x)/)') (i,i = imin,imax)
        do i=2,m
          ii = (i-1)*(i-2)/2
          mm = imin + ii
          kk = min(i-1,imax) + ii
          if(mm.le.kk) then
            write(lulog,'(i4,1x,10i15)') i,(iv(j),j=mm,kk)
          end if
        end do
      end do
      write(*,*)
      return
      end
