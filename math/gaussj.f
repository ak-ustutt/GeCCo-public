      subroutine gaussj(a,n,np)
C-----------------------------------------------------------------------
C     Gauss-Jordan elimination routine taken from Numerical Recipes
C     online book: www.nr.com.
C     Slightly modified to deal with double precision values.
C-----------------------------------------------------------------------
      implicit none
      integer n,nmax,np
      real*8 a(np,np)
      parameter (nmax=600)
      integer i,icol,irow,j,k,l,l1,indxc(nmax),indxr(nmax),ipiv(nmax)
      real*8 big,dum,pivinv
      do j=1,n
         ipiv(j)=0
      enddo
      do i=1,n
         big=0.d0
         do j=1,n
            if(ipiv(j).ne.1)then
               do k=1,n
                  if(ipiv(k).eq.0)then
                     if(dabs(a(j,k)).ge.big)then
                        big=dabs(a(j,k))
                        irow=j
                        icol=k
                     endif
                  endif
               enddo
            endif
         enddo
         ipiv(icol)=ipiv(icol)+1
         if(irow.ne.icol)then
            do l=1,n
               dum=a(irow,l)
               a(irow,l)=a(icol,l)
               a(icol,l)=dum
            enddo
         endif
         indxr(i)=irow
         indxc(i)=icol
         if(a(icol,icol).eq.0.d0)then
            write(6,'(''Singular matrix in gaussj'')')
            stop
         endif
         pivinv=1.d0/a(icol,icol)
         a(icol,icol)=1.d0
         do l=1,n
            a(icol,l)=a(icol,l)*pivinv
         enddo
         do l1=1,n
            if(l1.ne.icol)then
               dum=a(l1,icol)
               a(l1,icol)=0.d0
               do l=1,n
                  a(l1,l)=a(l1,l)-a(icol,l)*dum
               enddo
            endif
         enddo
      enddo
      do l=n,1,-1
         if(indxr(l).ne.indxc(l))then
            do k=1,n
               dum=a(k,indxr(l))
               a(k,indxr(l))=a(k,indxc(l))
               a(k,indxc(l))=dum
            enddo
         endif
      enddo
      return
      end
