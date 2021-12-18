c-----------------------------------------------------------------------
      subroutine prtrslt(v,m)
      implicit real*8(a-h,o-z)
c
c     ----- print out the lower triangle of an (anti-symmetric) matrix 
c           stored as strict lower triangle (actually: upper triangle) -----
c
      include 'stdunit.h'
      dimension v(m*(m-1)/2)

      max=5
      imax = 0
      do while(imax.lt.m)
        imin = imax+1
        imax = min(imax+max,m)
        write(lulog,'(/,5x,10(6x,i4,5x)/)') (i,i = imin,imax)
        do i=2,m
          ii = i*(i-2)/2
          mm = imin + ii
          kk = min(i,imax) + ii
          if(mm.le.kk) then
            write(lulog,'(i4,1x,10e15.7)') i,(v(j),j=mm,kk)
          end if
        end do
      end do
      write(*,*)
      return
      end
